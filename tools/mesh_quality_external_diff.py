"""Differential cross-check of :mod:`tensorcad.meshquality` against external oracles.

The anchor this serves -- ``GAP-CAD-D1-MESH-QUALITY-LIB`` -- asks for
"与 meshmagick+Gmsh 指标对拍（≥3 网格、指标差异<1e-6 或差异解释）": compare the
library's per-element measures with **meshmagick** and **Gmsh** on at least three meshes,
and either reach 1e-6 or explain the difference.

Neither package imports on this host (no loong64 wheel for ``gmsh``/``vtk``, and
``meshmagick`` is not on the package index), so the differential is run against the code
those tools actually execute:

* ``verdict_oracle`` -- the **Verdict** library (Sandia, the definitions of VTK's
  ``vtkMeshQuality``). ``meshmagick.print_quality()`` delegates to
  ``vtk.vtkMeshQuality``, so Verdict *is* meshmagick's arithmetic.
* ``gmsh_oracle`` -- a harness over ``libgmsh`` calling
  ``gmsh::model::mesh::getElementQualities``, i.e. Gmsh's own
  ``MElement::gammaShapeMeasure`` / ``minSICNShapeMeasure`` / ...

Both oracles are built by ``tools/mesh_quality_oracles/build_oracles.sh`` and are *not*
vendored (third-party code, tens of MB). This driver therefore reports three states
honestly instead of silently skipping:

``RAN``            oracle executed, per-element differences below.
``NOT_BUILT``      the oracle binary was not supplied (see the build script).
``NOT_RUN``        oracle supplied but failed to execute (message recorded).

Usage::

    python tools/mesh_quality_external_diff.py --verdict-oracle <path> [--gmsh-oracle <path>]
    python tools/mesh_quality_external_diff.py --verdict-oracle <p> --require-oracles

Exit codes: 0 = every oracle that ran agrees within tolerance (or no oracle ran and
``--require-oracles`` was not given); 1 = an executed oracle disagreed; 2 = an oracle was
required but did not run.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from tensorcad.meshquality.metrics import (  # noqa: E402
    hex_metrics,
    quad_metrics,
    tet_metrics,
    tri_metrics,
)

ARTIFACT_DIR = REPO / "artifacts" / "mesh_quality_external_diff"
TOLERANCE = 1.0e-6

#: our metric -> external metric, per element kind and per oracle. Every pair here compares
#: the *same quantity*; differences beyond :data:`TOLERANCE` are reported, never hidden.
MAP: dict[str, dict[str, dict[str, str]]] = {
    "verdict": {
        "tri": {
            "area": "area",
            "edge_ratio": "edge_ratio",
            "min_angle": "min_angle",
            "max_angle": "max_angle",
            "scaled_jacobian": "scaled_jacobian",
            "skewness": "equiangle_skew",
        },
        "quad": {
            "aspect_ratio": "aspect_ratio",
            "edge_ratio": "edge_ratio",
            "min_angle": "min_angle",
            "max_angle": "max_angle",
            "scaled_jacobian": "scaled_jacobian",
            "skewness": "equiangle_skew",
        },
        "tet": {
            "volume": "volume",
            "edge_ratio": "edge_ratio",
            "scaled_jacobian": "scaled_jacobian",
        },
        "hex": {
            "scaled_jacobian": "scaled_jacobian",
            "skewness": "equiangle_skew",
        },
    },
    "gmsh": {
        "tri": {
            "gamma": "gamma",
            "eta": "eta",
            "min_edge": "minEdge",
            "max_edge": "maxEdge",
            "area": "volume",
        },
        "quad": {
            "eta": "eta",
            "min_edge": "minEdge",
            "max_edge": "maxEdge",
            "area": "volume",
        },
        "tet": {
            "gamma": "gamma",
            "eta": "eta",
            "min_edge": "minEdge",
            "max_edge": "maxEdge",
            "volume": "volume",
            "gmsh_icn": "minSICN",
        },
        "hex": {
            "eta": "eta",
            "min_edge": "minEdge",
            "max_edge": "maxEdge",
            "volume": "volume",
        },
    },
}

#: Unit changes: the oracle and we measure the same quantity in different units. The factor is
#: ``external / canonical``; the tool divides the oracle's values by it before comparing and
#: records the *measured* spread of the ratio, so a factor that is not really constant is
#: reported as a mismatch instead of being papered over.
SCALE: dict[tuple[str, str, str], tuple[float, str]] = {
    ("verdict", "tet", "scaled_jacobian"): (
        math.sqrt(2.0),
        "Verdict normalises the tet Jacobian by the *regular* tetrahedron (regular tet = 1); the "
        "canonical metric divides by the Jacobian's column norms, which give 1/sqrt(2) for the "
        "same regular tet. Exact unit change, ratio external/ours = sqrt(2).",
    ),
}

#: Pairs that are *known* to be different quantities and are reported as explained
#: differences instead of being compared (metric, external name, why).
EXPLAINED: dict[str, dict[str, list[tuple[str, str, str]]]] = {
    "verdict": {
        "tri": [
            ("gamma", "radius_ratio", "gamma = 1 / radius_ratio (Verdict radius_ratio = R/(2r)); "
                                      "compared after inversion"),
            ("aspect_ratio", "aspect_ratio",
             "different quantity, not a unit change: the canonical triangle aspect ratio is "
             "max edge / (2 * inradius-like edge form) scaled so an equilateral triangle is 1, "
             "Verdict's tri_aspect_ratio is R/(2r) with a different normalisation; the measured "
             "ratio external/ours is not constant across elements (spread recorded below), so no "
             "single factor reconciles them"),
        ],
        "tet": [
            ("skewness", "equiangle_skew",
             "different quantity: Verdict's tet_equiangle_skew (V_TetMetric.cpp) measures the six "
             "*interior dihedral* angles with equilibrium angle 70.5287793655093 deg; the "
             "canonical skewness measures the twelve *face* angles with equilibrium 60 deg. Both "
             "are equiangle skew, of different angle sets; ratio is not constant"),
            ("min_dihedral", "minimum_angle",
             "different quantity: Verdict's tet_minimum_angle uses the *unsigned* face-normal dot "
             "product, i.e. it reports the supplement of the interior dihedral for the same edge; "
             "the canonical min_dihedral is the interior dihedral, so the two agree only when the "
             "dihedral is 90 deg"),
        ],
    },
    "gmsh": {
        "tet": [
            ("min_sicn", "minSICN", "our min_sicn is normalised to the regular tetrahedron, Gmsh "
                                   "to its right-corner reference tetrahedron; gmsh_icn carries "
                                   "Gmsh's unit, so min_sicn is not the same number by construction"),
        ],
    },
}

METRIC_FN = {"tri": tri_metrics, "quad": quad_metrics, "tet": tet_metrics, "hex": hex_metrics}


# --------------------------------------------------------------------------------------
# meshes: >= 3 shapes, constant, reproducible, with analytic expectations
# --------------------------------------------------------------------------------------
def _grid_hexes(nx: int, ny: int, nz: int, scale=1.0, jitter=0.0):
    pts = []
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                x, y, z = i / nx, j / ny, k / nz
                if jitter:
                    h = jitter * (1.0 / max(nx, ny, nz))
                    x += h * math.sin(7.0 * i + 3.0 * j + k)
                    y += h * math.sin(5.0 * j + 2.0 * k + i)
                    z += h * math.sin(11.0 * k + j)
                pts.append((scale * x, scale * y, scale * z))
    idx = lambda i, j, k: i + (nx + 1) * (j + (ny + 1) * k)  # noqa: E731
    hexes = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                hexes.append(
                    (
                        idx(i, j, k),
                        idx(i + 1, j, k),
                        idx(i + 1, j + 1, k),
                        idx(i, j + 1, k),
                        idx(i, j, k + 1),
                        idx(i + 1, j, k + 1),
                        idx(i + 1, j + 1, k + 1),
                        idx(i, j + 1, k + 1),
                    )
                )
    return pts, hexes


def _kuhn_tets(n=1, scale=1.0):
    """Unit cube split into 6 tetrahedra (Kuhn), optionally n times per axis."""
    pts, hexes = _grid_hexes(n, n, n, scale=scale)
    # 6-tet split of a hex with node order (bottom 0..3 CCW, top 4..7)
    split = [
        (0, 1, 2, 6),
        (0, 2, 3, 6),
        (0, 3, 7, 6),
        (0, 7, 4, 6),
        (0, 4, 5, 6),
        (0, 5, 1, 6),
    ]
    tets = []
    for h in hexes:
        for s in split:
            tets.append(tuple(h[i] for i in s))
    return pts, tets


def _tri_surface(nx=4, ny=3, warp=0.0):
    pts = []
    for j in range(ny + 1):
        for i in range(nx + 1):
            x, y = i / nx, j / ny
            z = warp * math.sin(math.pi * x) * math.sin(math.pi * y)
            pts.append((x, y, z))
    idx = lambda i, j: i + (nx + 1) * j  # noqa: E731
    tris = []
    for j in range(ny):
        for i in range(nx):
            tris.append((idx(i, j), idx(i + 1, j), idx(i + 1, j + 1)))
            tris.append((idx(i, j), idx(i + 1, j + 1), idx(i, j + 1)))
    return pts, tris


def _quad_patch(nx=3, ny=2, skew=0.0):
    pts = []
    for j in range(ny + 1):
        for i in range(nx + 1):
            x, y = i / nx, j / ny
            pts.append((x + skew * y / (ny + 1), y, 0.0))
    idx = lambda i, j: i + (nx + 1) * j  # noqa: E731
    quads = []
    for j in range(ny):
        for i in range(nx):
            quads.append((idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)))
    return pts, quads


def _sliver_tets():
    """Degenerate-ish tets: the mesh-quality library must agree with the oracles *here*."""
    pts = [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
        (0.3, 0.3, 1.0e-4),
        (0.5, 0.5, 1.0),
        (0.9, 0.05, 2.0e-3),
    ]
    tets = [(0, 1, 2, 3), (0, 1, 2, 4), (1, 2, 3, 5), (0, 1, 4, 6)]
    return pts, tets


def meshes() -> dict[str, tuple[str, list, list]]:
    """The differential corpus: 6 meshes over the four element kinds."""
    out: dict[str, tuple[str, list, list]] = {}
    out["hex_unit"] = ("hex",) + _grid_hexes(2, 2, 2)
    out["hex_distorted"] = ("hex",) + _grid_hexes(3, 2, 2, jitter=0.35)
    out["tet_cube6"] = ("tet",) + _kuhn_tets(1)
    out["tet_sliver"] = ("tet",) + _sliver_tets()
    out["tri_warped"] = ("tri",) + _tri_surface(5, 4, warp=0.25)
    out["quad_skewed"] = ("quad",) + _quad_patch(4, 3, skew=0.4)
    return out


# --------------------------------------------------------------------------------------
# oracle plumbing
# --------------------------------------------------------------------------------------
def write_mesh_file(path: Path, kind: str, points, elems) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write("# mq-mesh v1\n")
        f.write(f"nodes {len(points)}\n")
        for p in points:
            f.write(f"{p[0]!r} {p[1]!r} {p[2]!r}\n")
        f.write(f"{kind} {len(elems)}\n")
        for e in elems:
            f.write(" ".join(str(int(i)) for i in e) + "\n")


def run_oracle(binary: str | None, path: Path, prefix: str) -> tuple[str, dict]:
    """Run one oracle. Returns ``(status, {kind: {metric: [values]}})``."""
    if not binary:
        return "NOT_BUILT", {}
    if not Path(binary).exists():
        return "NOT_BUILT", {}
    try:
        res = subprocess.run([binary, str(path)], capture_output=True, text=True, timeout=300)
    except Exception as exc:  # pragma: no cover - environment dependent
        return f"NOT_RUN: {exc}", {}
    if res.returncode != 0:
        return f"NOT_RUN: exit {res.returncode}: {res.stderr.strip()[:200]}", {}
    data: dict[str, dict[str, list]] = {}
    for line in res.stdout.splitlines():
        parts = line.split()
        if len(parts) != 5 or parts[0] != prefix:
            continue
        kind, metric, idx, value = parts[1], parts[2], int(parts[3]), float(parts[4])
        bucket = data.setdefault(kind, {}).setdefault(metric, [])
        while len(bucket) <= idx:
            bucket.append(float("nan"))
        bucket[idx] = value
    return "RAN", data


def compare(ours: dict, external: dict, mapping: dict[str, str]) -> dict:
    """Per-metric agreement statistics between two equal-length value lists."""
    out = {}
    for mine, theirs in mapping.items():
        a = ours.get(mine)
        b = external.get(theirs)
        if a is None or b is None:
            out[mine] = {"external": theirs, "status": "NOT_RETURNED"}
            continue
        if len(a) != len(b):
            out[mine] = {"external": theirs, "status": "LENGTH_MISMATCH",
                         "lengths": [len(a), len(b)]}
            continue
        diffs = []
        worst = (0.0, -1)
        for i, (x, y) in enumerate(zip(a, b)):
            if x != x and y != y:
                continue
            if x != x or y != y:
                diffs.append((float("inf"), i, x, y))
                continue
            scale = max(abs(x), abs(y), 1.0e-300)
            rel = abs(x - y) / scale
            diffs.append((rel, i, x, y))
            if rel > worst[0]:
                worst = (rel, i)
        if not diffs:
            out[mine] = {"external": theirs, "status": "EMPTY"}
            continue
        worst_rel = max(d[0] for d in diffs)
        i = max(diffs)[1]
        out[mine] = {
            "external": theirs,
            "status": "AGREE" if worst_rel <= TOLERANCE else "DIFFER",
            "max_rel_diff": worst_rel,
            "worst_index": i,
            "ours_at_worst": float(a[i]) if a[i] == a[i] else None,
            "external_at_worst": float(b[i]) if b[i] == b[i] else None,
            "n_elements": len(diffs),
            "n_beyond_tol": sum(1 for d in diffs if d[0] > TOLERANCE),
        }
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--verdict-oracle", default=os.environ.get("TP_VERDICT_ORACLE"))
    ap.add_argument("--gmsh-oracle", default=os.environ.get("TP_GMSH_ORACLE"))
    ap.add_argument("--verdict-version", default=os.environ.get("TP_VERDICT_VERSION", "unknown"))
    ap.add_argument("--gmsh-version", default=os.environ.get("TP_GMSH_VERSION", "unknown"))
    ap.add_argument("--require-oracles", action="store_true",
                    help="exit 2 unless both oracles actually ran (CI mode)")
    ap.add_argument("--artifact", default=str(ARTIFACT_DIR / "oracle_diff.json"))
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)

    corpus = meshes()
    report = {
        "tolerance": TOLERANCE,
        "verdict_oracle": {"path": args.verdict_oracle, "version": args.verdict_version},
        "gmsh_oracle": {"path": args.gmsh_oracle, "version": args.gmsh_version},
        "meshes": {},
    }
    statuses = {"verdict": None, "gmsh": None}
    failed = []
    with tempfile.TemporaryDirectory() as td:
        for name, (kind, points, elems) in corpus.items():
            path = Path(td) / f"{name}.mesh"
            write_mesh_file(path, kind, points, elems)
            ours = METRIC_FN[kind](points, elems, backend="scalar").metrics
            entry = {"kind": kind, "n_elements": len(elems), "oracles": {}}
            for oracle, binary, prefix, version in (
                ("verdict", args.verdict_oracle, "W", args.verdict_version),
                ("gmsh", args.gmsh_oracle, "V", args.gmsh_version),
            ):
                status, data = run_oracle(binary, path, prefix)
                if status != "RAN":
                    statuses[oracle] = statuses[oracle] or status
                    entry["oracles"][oracle] = {"status": status}
                    continue
                statuses[oracle] = "RAN"
                ext = data.get(kind, {})
                mapping = MAP[oracle].get(kind, {})
                ext_cmp = dict(ext)
                unit_changes = {}
                for mine, theirs in mapping.items():
                    key = (oracle, kind, mine)
                    vals = ext.get(theirs)
                    if vals is None or key not in SCALE:
                        continue
                    factor, why = SCALE[key]
                    a = list(ours.get(mine) or [])
                    ratios = [y / x for x, y in zip(a, vals) if x and abs(x) > 1.0e-9]
                    ext_cmp[theirs] = [v / factor for v in vals]
                    unit_changes[mine] = {
                        "factor": factor,
                        "why": why,
                        "measured_ratio_min": min(ratios) if ratios else None,
                        "measured_ratio_max": max(ratios) if ratios else None,
                        "constant": bool(ratios) and abs(max(ratios) - min(ratios)) <= 1.0e-12,
                    }
                    if not unit_changes[mine]["constant"]:
                        failed.append(
                            f"{name}/{oracle}/{mine}: declared unit factor {factor:.12g} is not "
                            f"constant (ratio {min(ratios):.12g}..{max(ratios):.12g})"
                        )
                cmp = compare({k: list(v) for k, v in ours.items()}, ext_cmp, mapping)
                explained = []
                for mine, theirs, why in EXPLAINED.get(oracle, {}).get(kind, []):
                    if mine == "gamma" and theirs == "radius_ratio" and "radius_ratio" in ext:
                        inv = [1.0 / v if v else float("inf") for v in ext["radius_ratio"]]
                        cmp["gamma_vs_1_over_radius_ratio"] = compare(
                            {k: list(v) for k, v in ours.items()}, {"x": inv}, {"gamma": "x"}
                        )["gamma"]
                        explained.append({"pair": [mine, theirs], "why": why})
                    else:
                        ev = {"pair": [mine, theirs], "why": why}
                        a = list(ours.get(mine) or [])
                        b = ext.get(theirs)
                        if a and b and len(a) == len(b):
                            rr = [y / x for x, y in zip(a, b) if x and abs(x) > 1.0e-9]
                            if rr:
                                ev["measured_ratio_min"] = min(rr)
                                ev["measured_ratio_max"] = max(rr)
                                ev["unit_change"] = abs(max(rr) - min(rr)) <= 1.0e-12
                        explained.append(ev)
                entry["oracles"][oracle] = {"status": "RAN", "version": version,
                                            "comparisons": cmp, "explained": explained,
                                            "unit_changes": unit_changes}
                for metric, c in cmp.items():
                    # A truncated oracle response (fewer values than elements) is a silent loss of
                    # evidence, not a pass: the CI gate must fail on it too.
                    if c.get("status") in ("DIFFER", "LENGTH_MISMATCH"):
                        rel = c.get("max_rel_diff")
                        detail = f"rel {rel:.3e}" if rel is not None else c["status"].lower()
                        failed.append(f"{name}/{oracle}/{metric}: {detail}")
            report["meshes"][name] = entry

    report["oracle_status"] = statuses
    report["mismatches"] = failed
    Path(args.artifact).parent.mkdir(parents=True, exist_ok=True)
    Path(args.artifact).write_text(json.dumps(report, indent=2, ensure_ascii=False))

    if not args.quiet:
        print(f"== mesh-quality external differential (tolerance {TOLERANCE:g}) ==")
        for oracle, status in statuses.items():
            print(f"oracle {oracle:8s}: {status}")
        for name, entry in report["meshes"].items():
            for oracle, res in entry["oracles"].items():
                if res["status"] != "RAN":
                    continue
                for metric, c in res["comparisons"].items():
                    if c.get("status") in ("AGREE", "DIFFER"):
                        print(f"  {name:14s} {oracle:8s} {metric:26s} "
                              f"max rel {c['max_rel_diff']:.3e} "
                              f"({c['n_beyond_tol']}/{c['n_elements']} beyond) {c['status']}")
        print(f"artifact: {args.artifact}")
        if failed:
            print("MISMATCHES:")
            for f in failed:
                print("  " + f)

    ran = [o for o, s in statuses.items() if s == "RAN"]
    if args.require_oracles and len(ran) < 2:
        print(f"REQUIRED oracle(s) did not run: {statuses}", file=sys.stderr)
        return 2
    if failed:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())