"""Reproducer / evidence probe for ``TriangleGrid.nearest_surface`` query cost.

Background
----------
``TriangleGrid.nearest_surface`` is *exact*: it scans Chebyshev shells of grid
cells around the query point and stops only when no unexamined triangle can be
closer.  On a mesh whose bounding box is degenerate in one axis (a plate, a
shell patch, any open surface) the automatic resolution is driven by the
box *volume*, so the two in-plane axes get thousands of cells while the answer
distance stays O(1) metres: the shell search then has to clear a ball of
hundreds of cells per axis.  The same happens for a closed but thin body when
the query point sits a few centimetres outside it.

This script measures, per case, the grid geometry, the build time, and the wall
clock of a fixed batch of nearest-surface queries, each answer cross-checked
against ``brute_force_surface_distance``.  A per-query wall-clock ceiling is
enforced with ``signal.setitimer`` so a pathological case is *recorded as a
timeout* instead of hanging the run.

Usage::

    python scripts/probe_nearest_surface_cost.py --out artifacts/spatial_nearest_cost.json
    python scripts/probe_nearest_surface_cost.py --label after      # tag the run
"""
from __future__ import annotations

import argparse
import json
import math
import signal
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tensorcad.geometry import Mesh  # noqa: E402
from tensorcad.spatial_index import (  # noqa: E402
    TriangleGrid,
    brute_force_surface_distance,
)


class _Timeout(Exception):
    pass


def _alarm(signum, frame):  # noqa: ARG001
    raise _Timeout


def timed(fn, limit: float):
    """Run ``fn`` with a wall-clock ceiling; return (value, seconds, timed_out)."""
    old = signal.signal(signal.SIGALRM, _alarm)
    signal.setitimer(signal.ITIMER_REAL, limit)
    t0 = time.perf_counter()
    try:
        value = fn()
        timed_out = False
    except _Timeout:
        value = None
        timed_out = True
    finally:
        elapsed = time.perf_counter() - t0
        signal.setitimer(signal.ITIMER_REAL, 0.0)
        signal.signal(signal.SIGALRM, old)
    return value, elapsed, timed_out


def sheet_mesh(size: float = 1.0) -> Mesh:
    """Open 1x1 sheet at z=0: two triangles, bbox is degenerate in z."""
    s = size
    return Mesh(
        vertices=[(0.0, 0.0, 0.0), (s, 0.0, 0.0), (s, s, 0.0), (0.0, s, 0.0)],
        faces=[(0, 1, 2), (0, 2, 3)],
    )


def plate_box(thickness: float = 1e-3, size: float = 1.0) -> Mesh:
    """Closed thin box (a ship plating panel) with the same 1x1 footprint."""
    x0, y0, z0 = 0.0, 0.0, 0.0
    x1, y1, z1 = size, size, thickness
    v = [(x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
         (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1)]
    q = [(0, 2, 1), (0, 3, 2), (4, 5, 6), (4, 6, 7),
         (0, 1, 5), (0, 5, 4), (1, 2, 6), (1, 6, 5),
         (2, 3, 7), (2, 7, 6), (3, 0, 4), (3, 4, 7)]
    return Mesh(vertices=v, faces=q)


def unit_box() -> Mesh:
    return plate_box(thickness=1.0)


def hull_like(n_rings: int = 24, n_stations: int = 48) -> Mesh:
    """Closed, non-degenerate body: a wavy box-loft, N ~ 2*n_rings*n_stations."""
    verts: list[tuple[float, float, float]] = []
    faces: list[tuple[int, int, int]] = []
    for i in range(n_stations + 1):
        u = i / n_stations
        x = 10.0 * u
        beam = 2.0 * (1.0 - abs(2.0 * u - 1.0) ** 2)
        depth = 0.8 * (1.0 - abs(2.0 * u - 1.0) ** 2) + 1e-3
        for j in range(n_rings):
            a = 2.0 * math.pi * j / n_rings
            verts.append((x, 0.5 * beam * math.cos(a), 0.5 * depth * math.sin(a)))
    for i in range(n_stations):
        for j in range(n_rings):
            a = i * n_rings + j
            b = i * n_rings + (j + 1) % n_rings
            c = (i + 1) * n_rings + (j + 1) % n_rings
            d = (i + 1) * n_rings + j
            faces.append((a, b, c))
            faces.append((a, c, d))
    return Mesh(vertices=verts, faces=faces)


def query_points(lo, hi, n: int) -> list[tuple[float, float, float]]:
    """Deterministic satellite of probe points around the body's centre.

    Radii are fractions of the *body size* (bbox diagonal), not of each axis
    extent: on a plate that deliberately puts points many cell-widths away from
    the surface along the thin axis, which is the case that stresses the shell
    search.  ``n`` is rounded up to the next multiple of the direction count.
    """
    centre = tuple(0.5 * (lo[i] + hi[i]) for i in range(3))
    diag = math.sqrt(sum((hi[i] - lo[i]) ** 2 for i in range(3))) or 1.0
    # Fibonacci sphere: fixed, reproducible directions.
    dirs = []
    m = 8
    golden = math.pi * (3.0 - math.sqrt(5.0))
    for k in range(m):
        z = 1.0 - 2.0 * (k + 0.5) / m
        r = math.sqrt(max(0.0, 1.0 - z * z))
        theta = golden * k
        dirs.append((r * math.cos(theta), r * math.sin(theta), z))
    fractions = [0.05, 0.15, 0.3, 0.5, 0.75, 1.0]
    out = []
    for f in fractions:
        for d in dirs:
            out.append(tuple(centre[i] + 0.5 * diag * f * d[i] for i in range(3)))
    return out[: max(n, len(out))] if n <= len(out) else out


CASES = {
    "sheet_2tri": lambda: sheet_mesh(1.0),
    "plate_box_thin": lambda: plate_box(1e-3, 1.0),
    "unit_box": unit_box,
    "hull_like_2k": lambda: hull_like(24, 48),
}


def probe_case(name: str, mesh_fn, n_queries: int, build_limit: float, query_limit: float) -> dict:
    t0 = time.perf_counter()
    grid = TriangleGrid(mesh_fn())
    build_seconds = time.perf_counter() - t0
    mesh = grid.mesh
    xs = [v[0] for v in mesh.vertices]
    ys = [v[1] for v in mesh.vertices]
    zs = [v[2] for v in mesh.vertices]
    lo = (min(xs), min(ys), min(zs))
    hi = (max(xs), max(ys), max(zs))
    points = query_points(lo, hi, n_queries)

    total = 0.0
    timeouts = 0
    mismatches = 0
    worst_abs = 0.0
    checked = 0
    for p in points:
        hit, elapsed, timed_out = timed(lambda p=p: grid.nearest_surface(p), query_limit)
        total += elapsed
        if timed_out:
            timeouts += 1
            if timeouts >= 2:      # two strikes is enough evidence; stop burning wall clock
                break
            continue
        dist = math.inf if hit is None else hit.distance
        ref = brute_force_surface_distance(mesh, p)
        checked += 1
        if math.isinf(ref) != math.isinf(dist) or (
            not math.isinf(ref) and abs(dist - ref) > 1e-9 * max(1.0, ref)
        ):
            mismatches += 1
            worst_abs = max(worst_abs, abs(dist - ref))
    return {
        "n_triangles": len(mesh.faces),
        "resolution": list(grid.resolution),
        "n_cells": grid.n_cells,
        "n_oversized": grid.n_oversized,
        "build_seconds": round(build_seconds, 6),
        "queries_attempted": len(points) if not timeouts else checked + timeouts,
        "queries_checked_vs_bruteforce": checked,
        "seconds_total": round(total, 6),
        "seconds_per_query": round(total / max(1, checked + timeouts), 6),
        "timeouts": timeouts,
        "timeout_limit_seconds": query_limit,
        "bruteforce_mismatches": mismatches,
        "worst_abs_error": worst_abs,
    }


DOC_BEGIN = "<!-- nearest-cost:begin -->"
DOC_END = "<!-- nearest-cost:end -->"
METRIC_ORDER = (
    ("n_triangles", "三角面"),
    ("resolution", "自动分辨率"),
    ("n_cells", "网格单元"),
    ("build_seconds", "建索引 s"),
    ("seconds_per_query", "每查询 s"),
    ("timeouts", "超时数"),
    ("bruteforce_mismatches", "与全扫不符"),
)


def render_doc_block(runs: list[dict]) -> str:
    """Render the before/after table from the probe JSON (single source of truth)."""
    by_label: dict[str, dict] = {}
    for run in runs:
        by_label[str(run.get("label"))] = run
    cases: list[str] = []
    for run in runs:
        for name in run.get("cases", {}):
            if name not in cases:
                cases.append(name)
    # Fixed column order: the un-fixed code first, then the fix, so the table reads
    # left-to-right as the change.
    order = [tag for tag in ("before", "after") if tag in by_label]
    order += [tag for tag in sorted(by_label) if tag not in order]
    lines: list[str] = []
    lines.append("| 案例 | 指标 | " + " | ".join(order) + " |")
    lines.append("|---|---|" + "|".join("---" for _ in order) + "|")
    for case in cases:
        for key, label in METRIC_ORDER:
            cells = []
            for tag in order:
                rec = by_label[tag].get("cases", {}).get(case)
                value = "-" if rec is None else rec.get(key)
                if isinstance(value, list):
                    value = "×".join(str(v) for v in value)
                elif isinstance(value, float):
                    value = f"{value:.6g}"
                cells.append(str(value))
            lines.append(f"| `{case}` | {label} | " + " | ".join(cells) + " |")
    limits = order
    lines.append("")
    lines.append(
        "每查询上限 "
        + " s / ".join(
            f"{tag}: {by_label[tag]['cases'][next(iter(by_label[tag]['cases']))]['timeout_limit_seconds']:.0f}s"
            for tag in limits if by_label[tag].get("cases")
        )
        + "（超过记为“超时”，不是 0）。"
    )
    return "\n".join(lines) + "\n"


def check_doc(json_path: Path, doc_path: Path, write: bool) -> list[str]:
    """Compare (or rewrite) the marked table in ``doc_path`` with the JSON record."""
    payload = json.loads(json_path.read_text())
    runs = payload.get("runs", [])
    rendered = render_doc_block(runs)
    text = doc_path.read_text(encoding="utf-8")
    begin = text.index(DOC_BEGIN) + len(DOC_BEGIN)
    end = text.index(DOC_END)
    current = text[begin:end]
    if current.strip() == rendered.strip():
        return []
    if write:
        doc_path.write_text(text[:begin] + "\n" + rendered + "\n" + text[end:], encoding="utf-8")
        return []
    return [
        f"{doc_path.name}: nearest-cost table is stale -- re-run "
        f"{Path(__file__).name} --render-doc (table differs from {json_path.name})"
    ]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default="artifacts/spatial_nearest_cost.json")
    ap.add_argument("--label", default="run", help="tag recorded in the JSON")
    ap.add_argument("--queries", type=int, default=8)
    ap.add_argument("--query-limit", type=float, default=20.0,
                    help="per-query wall-clock ceiling (s); a hit is recorded as a timeout")
    ap.add_argument("--cases", default=",".join(CASES),
                    help=f"subset of {sorted(CASES)}")
    ap.add_argument("--doc", default="docs/spatial_index_nearest_cost.md",
                    help="markdown file whose nearest-cost table is rendered from the JSON")
    ap.add_argument("--render-doc", action="store_true",
                    help="rewrite the marked table in --doc from the recorded runs")
    ap.add_argument("--check-doc", action="store_true",
                    help="fail when the marked table in --doc is stale")
    args = ap.parse_args(argv)

    if args.render_doc or args.check_doc:
        out = Path(args.out)
        if not out.is_absolute():
            out = ROOT / out
        doc = Path(args.doc)
        if not doc.is_absolute():
            doc = ROOT / doc
        problems = check_doc(out, doc, write=args.render_doc)
        if args.check_doc and problems:
            for problem in problems:
                print(problem, file=sys.stderr)
            return 1
        print("ok" if not problems else "stale")
        return 0


    report: dict[str, object] = {
        "label": args.label,
        "python": sys.version.split()[0],
        "platform": sys.platform,
        "cases": {},
    }
    failed = 0
    for name in [c.strip() for c in args.cases.split(",") if c.strip()]:
        if name not in CASES:
            print(f"unknown case {name!r}", file=sys.stderr)
            return 2
        rec = probe_case(name, CASES[name], args.queries, 60.0, args.query_limit)
        report["cases"][name] = rec  # type: ignore[index]
        print(
            f"{name:16s} N={rec['n_triangles']:6d} res={rec['resolution']!s:18s} "
            f"cells={rec['n_cells']:9d} build={rec['build_seconds']:.4f}s "
            f"q={rec['seconds_per_query']:.5f}s/q timeouts={rec['timeouts']} "
            f"mismatch={rec['bruteforce_mismatches']}",
            flush=True,
        )
        if rec["bruteforce_mismatches"]:
            failed = 1

    out = Path(args.out)
    if not out.is_absolute():
        out = ROOT / out
    out.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if out.exists():
        try:
            previous = json.loads(out.read_text())
        except json.JSONDecodeError:
            previous = {}
        runs = previous.get("runs")
        if isinstance(runs, list):
            payload_runs = runs
        else:
            payload_runs = [previous] if previous else []
        payload_runs.append(report)
        out.write_text(json.dumps({"runs": payload_runs}, indent=2, sort_keys=True) + "\n")
    else:
        out.write_text(json.dumps({"runs": [report]}, indent=2, sort_keys=True) + "\n")
    print(f"wrote {out}")
    return failed


if __name__ == "__main__":
    sys.exit(main())