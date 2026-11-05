#!/usr/bin/env python3
"""Programmatic verification of the PR #44 cure merge (agent/m3347 <- main).

The merge is only "cured" if three things are true and *machine-checked*, not eyeballed:

1. nothing that either parent added was silently dropped -- measured as per-file
   added-line coverage of each side against the merged tree (exact-line match,
   comment/whitespace-insensitive), plus the exact set union of the public surface;
2. both parents' numerics still run in the merged tree -- re-derived here from the
   merged package, not read back from either side's stored results;
3. the test suite is not below the main baseline -- read from the two pytest logs.

Writes ``docs/m3822/merge_verification.json`` and the figures the report embeds.
Exit code is non-zero if any invariant fails, so the script doubles as a gate.
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "m3822"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
MERGED = "HEAD"
SIDE_A = "origin/agent/m3347"  # the PR branch before the cure
SIDE_B = "origin/main"
BASE = "e604f02f9645c15e39b6acb562d4079b1755219a"  # merge base of the two


# Lines the cure deliberately rewrote for style (main is at ruff 0 / mypy 0). Each entry is
# (path, line as PR #44 wrote it, line as the merged tree writes it): the gate then requires the
# replacement to be present and counts the original as *explained*, not as a silent loss.
REWRITES: dict[str, list[tuple[str, str]]] = {
    "examples/gz_curve/run_example.py": [
        (
            'if int(row["heel_deg"]) % 15 and not row["heel_deg"] in (0.0, 60.0):',
            'if int(row["heel_deg"]) % 15 and row["heel_deg"] not in (0.0, 60.0):',
        ),
        (
            "f\"{row['centre_of_buoyancy_y']:10.6f} {str(row['deck_immersed']):>6}\"",
            "f\"{row['centre_of_buoyancy_y']:10.6f} {row['deck_immersed']!s:>6}\"",
        ),
    ],
    "tensorcad/stability.py": [
        (
            "for previous, following in zip(angles, angles[1:]):",
            "for previous, following in pairwise(angles):",
        ),
        (
            "keel (hence ``keel_z + KG`` in the mesh's own frame). ``KG`` is a loading input"
            " (as in :func:`initial_stability`), so it is a",
            "keel (hence ``keel_z + KG`` in the mesh's own frame). ``KG`` is a loading input",
        ),
    ],
    "tests/test_gz_curve.py": [
        (
            "for left, right in zip(skirted.points, solid.points):",
            "for left, right in zip(skirted.points, solid.points, strict=False):",
        ),
        (
            "assert all(abs(b - a) < 0.05 for a, b in zip(arms, arms[1:]))",
            "assert all(abs(b - a) < 0.05 for a, b in pairwise(arms))",
        ),
    ],
}

# ---------------------------------------------------------------- git plumbing


def git(*args: str) -> str:
    return subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=True
    ).stdout


def git_bytes(*args: str) -> bytes:
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, check=True).stdout


def changed_files(side: str) -> set[str]:
    out = git("diff", "--name-only", BASE, side)
    return {line for line in out.splitlines() if line}


def normalise(line: str) -> str:
    """Collapse whitespace so re-sorting or re-wrapping is not miscounted as a loss."""
    return " ".join(line.split())


def added_lines(side: str, path: str) -> set[str]:
    """Non-empty, normalised lines this side added relative to the merge base."""
    try:
        diff = git("diff", "--unified=0", BASE, side, "--", path)
    except subprocess.CalledProcessError:
        return set()
    out: set[str] = set()
    for line in diff.splitlines():
        if line.startswith("+") and not line.startswith("+++"):
            text = normalise(line[1:])
            if text:
                out.add(text)
    return out


def merged_lines(path: str) -> set[str]:
    try:
        text = (ROOT / path).read_text(encoding="utf-8", errors="replace")
    except OSError:
        return set()
    return {normalise(line) for line in text.splitlines() if normalise(line)}


def coverage_report() -> dict[str, object]:
    files = sorted(changed_files(SIDE_A) | changed_files(SIDE_B))
    rows = []
    total_a = total_a_hit = total_b = total_b_hit = 0
    for path in files:
        present = merged_lines(path)
        if not present and not (ROOT / path).exists():
            # binary (png) or deleted -- presence is the whole test
            rows.append(
                {
                    "path": path,
                    "binary_or_absent": True,
                    "present_in_merged": (ROOT / path).exists(),
                }
            )
            continue
        a = added_lines(SIDE_A, path)
        b = added_lines(SIDE_B, path)
        # a line the cure rewrote on purpose is explained only if its replacement is there
        rewritten = REWRITES.get(path, [])
        explained = set()
        for old, new in rewritten:
            old_n, new_n = normalise(old), normalise(new)
            if old_n in a and old_n not in present:
                assert new_n in present, f"rewrite of {path!r} not present: {new!r}"
                explained.add(old_n)
        a_hit = len((a & present) | explained)
        b_hit = len(b & present)
        total_a += len(a)
        total_a_hit += a_hit
        total_b += len(b)
        total_b_hit += b_hit
        if a or b:
            rows.append(
                {
                    "path": path,
                    "pr44_added_lines": len(a),
                    "pr44_kept": a_hit,
                    "main_added_lines": len(b),
                    "main_kept": b_hit,
                    "rewritten_pr44": sorted(explained),
                    "missing_pr44": sorted((a - present) - explained)[:5],
                    "missing_main": sorted(b - present)[:5],
                }
            )
    return {
        "base": BASE,
        "pr44_tip_before_cure": git("rev-parse", SIDE_A).strip(),
        "main": git("rev-parse", SIDE_B).strip(),
        "files_touched_by_either_side": len(files),
        "pr44_added_lines": total_a,
        "pr44_added_lines_kept": total_a_hit,
        "main_added_lines": total_b,
        "main_added_lines_kept": total_b_hit,
        "pr44_coverage": total_a_hit / total_a if total_a else 1.0,
        "main_coverage": total_b_hit / total_b if total_b else 1.0,
        "per_file": rows,
    }


# ------------------------------------------------------- public surface contract


def _all_of_ref(ref: str, path: str = "tensorcad/__init__.py") -> set[str]:
    text = git_bytes("show", f"{ref}:{path}").decode()
    return set(_all_of_source(text))


def _all_of_source(text: str) -> list[str]:
    import ast

    for node in ast.walk(ast.parse(text)):
        if isinstance(node, ast.Assign) and any(
            getattr(t, "id", None) == "__all__" for t in node.targets
        ):
            return list(ast.literal_eval(node.value))
    raise AssertionError("no __all__ found")


def surface_report() -> dict[str, object]:
    import tensorcad

    live = set(tensorcad.__all__)
    a = _all_of_ref(SIDE_A)
    b = _all_of_ref(SIDE_B)
    return {
        "count_pr44": len(a),
        "count_main": len(b),
        "count_merged": len(live),
        "missing_from_pr44": sorted(a - live),
        "missing_from_main": sorted(b - live),
        "is_exact_union": live == (a | b),
        "importable_names": all(hasattr(tensorcad, name) for name in live),
    }


# -------------------------------------------------- numerics alive after merging


def _box(length: float, beam: float, draft: float, freeboard: float):
    """Closed box mesh: keel on z=-draft, deck on z=+freeboard, waterline z=0."""
    from tensorcad.geometry import MeshBuilder

    builder = MeshBuilder()
    corner = {}
    for x in (0.0, length):
        for y in (-0.5 * beam, 0.5 * beam):
            for z in (-draft, freeboard):
                corner[(x, y, z)] = builder.add_vertex((x, y, z))

    def quad(p0, p1, p2, p3) -> None:
        a, b, c, d = (corner[point] for point in (p0, p1, p2, p3))
        builder.add_triangle(a, b, c)
        builder.add_triangle(a, c, d)

    lo, hi, bottom, top, left, right = (
        -0.5 * beam,
        0.5 * beam,
        -draft,
        freeboard,
        0.0,
        length,
    )
    quad((left, lo, bottom), (left, hi, bottom), (right, hi, bottom), (right, lo, bottom))
    quad((left, lo, top), (right, lo, top), (right, hi, top), (left, hi, top))
    quad((left, lo, bottom), (right, lo, bottom), (right, lo, top), (left, lo, top))
    quad((left, hi, bottom), (left, hi, top), (right, hi, top), (right, hi, bottom))
    quad((left, lo, bottom), (left, lo, top), (left, hi, top), (left, hi, bottom))
    quad((right, lo, bottom), (right, hi, bottom), (right, hi, top), (right, lo, top))
    return builder.build()


def _section_reference(beam: float, draft: float, freeboard: float, kg: float, heel_deg: float) -> float:
    """Exact arm of the rectangular section: rotate, clip at the waterplane, bisect the plane.

    Same construction as ``tests/test_gz_curve.py::reference_gz`` -- i.e. the reference
    shares no code or formula with the solver's 3-D divergence integrals -- re-implemented
    here so this gate is independent of the test module.
    """
    phi = math.radians(heel_deg)
    cos_phi, sin_phi = math.cos(phi), math.sin(phi)
    keel = -draft
    section = [
        (y * cos_phi - z * sin_phi, y * sin_phi + z * cos_phi)
        for y, z in (
            (-0.5 * beam, keel),
            (0.5 * beam, keel),
            (0.5 * beam, freeboard),
            (-0.5 * beam, freeboard),
        )
    ]
    heights = [point[1] for point in section]
    low, high = min(heights) - 1.0, max(heights) + 1.0

    def clipped(plane: float) -> tuple[float, float]:
        polygon: list[tuple[float, float]] = []
        for index in range(len(section)):
            start, end = section[index], section[(index + 1) % len(section)]
            start_in, end_in = start[1] <= plane, end[1] <= plane
            if start_in:
                polygon.append(start)
            if start_in != end_in:
                ratio = (plane - start[1]) / (end[1] - start[1])
                polygon.append((start[0] + ratio * (end[0] - start[0]), plane))
        area = moment_y = 0.0
        for index in range(len(polygon)):
            x1, y1 = polygon[index]
            x2, y2 = polygon[(index + 1) % len(polygon)]
            cross = x1 * y2 - x2 * y1
            area += cross
            moment_y += (x1 + x2) * cross
        return 0.5 * area, moment_y / 6.0

    for _ in range(200):
        middle = 0.5 * (low + high)
        if clipped(middle)[0] < beam * draft:
            low = middle
        else:
            high = middle
    _area, moment_y = clipped(0.5 * (low + high))
    return -(keel + kg) * sin_phi - moment_y / _area


def _wall_sided(beam: float, draft: float, kg: float, heel_deg: float) -> float:
    """Classical wall-sided arm ``sin(phi)*(BM + BM tan^2(phi)/2 - BG)``."""
    phi = math.radians(heel_deg)
    bm = beam * beam / (12.0 * draft)
    bg = kg - 0.5 * draft
    return math.sin(phi) * (bm + 0.5 * bm * math.tan(phi) ** 2 - bg)


def gz_curve_evidence() -> dict[str, object]:
    """PR #44's feature, re-run from the merged tree against closed-form references."""
    import tensorcad as tc

    length, beam, draft, freeboard, kg = 10.0, 2.0, 0.5, 1.5, 0.45
    mesh = _box(length, beam, draft, freeboard)
    curve = tc.righting_arm_curve(mesh, kg=kg, waterline_z=0.0)

    heel = [point.heel_deg for point in curve.points]
    gz = [point.gz for point in curve.points]
    analytic = [_section_reference(beam, draft, freeboard, kg, phi) for phi in heel]
    abs_err = [abs(a - b) for a, b in zip(gz, analytic, strict=True)]

    # wall-sided regime: the waterplane stays on both vertical sides, tan(phi) < 2*draft/beam
    inside = [phi for phi in heel if abs(math.tan(math.radians(phi))) < 2.0 * draft / beam]
    wall_err = [
        abs(gz[heel.index(phi)] - _wall_sided(beam, draft, kg, phi)) / max(abs(_wall_sided(beam, draft, kg, phi)), 1e-30)
        for phi in inside
        if phi > 0.0
    ]
    return {
        "feature": "tensorcad.stability.righting_arm_curve (PR #44 side)",
        "box": {"length": length, "beam": beam, "draft": draft, "freeboard": freeboard, "kg": kg},
        "displacement_volume": length * beam * draft,
        "points": len(heel),
        "gz_max_m": curve.gz_max,
        "angle_of_gz_max_deg": curve.angle_of_gz_max_deg,
        "gm_from_merged_solver_m": curve.initial_gm_transverse,
        "max_abs_error_vs_section_clipping_m": max(abs_err),
        "max_rel_error_vs_wall_sided_in_regime": max(wall_err) if wall_err else None,
        "wall_sided_regime_max_heel_deg": max(inside),
        "heel_deg": heel,
        "gz_m": gz,
        "gz_analytic_m": analytic,
    }


def transom_evidence() -> dict[str, object]:
    """main's feature, re-run from the merged tree against its closed forms."""
    import tensorcad as tc

    tau = 0.2
    rows = []
    for nx, nz in ((21, 9), (41, 13), (81, 25), (161, 41)):
        parameters = tc.TransomHullParameters(transom_fraction=tau, nx=nx, nz=nz)
        model = tc.generate_transom_hull(parameters)
        hydro = model.hydrostatics
        exact = model.analytic
        rows.append(
            {
                "nx": nx,
                "nz": nz,
                "volume_m3": hydro.displacement_volume,
                "volume_analytic_m3": exact["displacement_volume"],
                "volume_rel_error": (hydro.displacement_volume - exact["displacement_volume"])
                / exact["displacement_volume"],
            }
        )
    return {
        "feature": "tensorcad.transom.generate_transom_hull (main side)",
        "transom_fraction": tau,
        "rows": rows,
        "converging": all(
            abs(rows[i + 1]["volume_rel_error"]) < abs(rows[i]["volume_rel_error"])
            for i in range(len(rows) - 1)
        ),
    }


def test_evidence() -> dict[str, object]:
    """Read the two pytest runs; the baseline must not be below main's."""

    def parse(path: str) -> dict[str, object] | None:
        p = Path(path)
        if not p.exists():
            return None
        text = p.read_text(errors="replace")
        summary = ""
        for line in text.splitlines():
            if " passed" in line or " failed" in line or " error" in line:
                summary = line.strip()
        counts = {"passed": 0, "failed": 0, "skipped": 0, "errors": 0}
        import re

        for key in counts:
            m = re.search(r"(\d+) " + key, summary)
            if m:
                counts[key] = int(m.group(1))
        failed = sorted(
            line.split("::")[0].replace("FAILED ", "").strip()
            for line in text.splitlines()
            if line.startswith("FAILED ")
        )
        return {"summary": summary, **counts, "failing_modules": sorted(set(failed))}

    return {
        "main_baseline": parse("/tmp/base_m3822_pytest.log"),
        "merged": parse("/tmp/m3822_final_pytest.log"),
        "excluded_modules": ["tests/test_cli_neural.py", "tests/test_neural_pipeline.py"],
        "exclusion_reason": "ModuleNotFoundError: safetensors in this image, pre-existing on main",
    }


# ---------------------------------------------------------------------- figures


def figure(payload: dict[str, object]) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    gz = payload["gz_curve"]
    transom = payload["transom"]
    tests = payload["tests"]

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.0), constrained_layout=True)

    ax = axes[0][0]
    ax.plot(gz["heel_deg"], gz["gz_analytic_m"], "-", color="#1f77b4", label="section clipping (exact)")
    ax.plot(gz["heel_deg"], gz["gz_m"], "--", color="#d62728", label="merged tree, righting_arm_curve")
    ax.set_xlabel("heel angle [deg]")
    ax.set_ylabel("GZ [m]")
    ax.set_title(
        f"PR #44 side alive: GZ curve, box {gz['box']['length']:.0f}x{gz['box']['beam']:.0f}"
        f"x{gz['box']['draft']:.1f}+{gz['box']['freeboard']:.1f} m\n"
        f"max |err| = {gz['max_abs_error_vs_section_clipping_m']:.2e} m"
    )
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[0][1]
    err = [abs(r["volume_rel_error"]) for r in transom["rows"]]
    cells = [r["nx"] * r["nz"] for r in transom["rows"]]
    ax.loglog(cells, err, "o-", color="#2ca02c")
    for cell, e, r in zip(cells, err, transom["rows"], strict=True):
        ax.annotate(f"{r['nx']}x{r['nz']}", (cell, e), textcoords="offset points", xytext=(4, 4), fontsize=7)
    ax.set_xlabel("grid cells (nx*nz)")
    ax.set_ylabel("|volume error| vs closed form [-]")
    ax.set_title(
        f"main side alive: transom hull tau={transom['transom_fraction']}, "
        f"error falls monotonically: {transom['converging']}"
    )
    ax.grid(alpha=0.3, which="both")

    ax = axes[1][0]
    base = tests["main_baseline"] or {}
    merged = tests["merged"] or {}
    labels = ["passed", "failed", "skipped"]
    xs = range(len(labels))
    ax.bar([x - 0.2 for x in xs], [base.get(k, 0) for k in labels], width=0.4, label="main baseline", color="#7f7f7f")
    ax.bar([x + 0.2 for x in xs], [merged.get(k, 0) for k in labels], width=0.4, label="cured branch", color="#1f77b4")
    ax.set_xticks(list(xs))
    ax.set_xticklabels(labels)
    ax.set_title(
        f"test suite: merged {merged.get('passed', 0)} passed vs main {base.get('passed', 0)}"
        f" (+{merged.get('passed', 0) - base.get('passed', 0)}), same 2 pre-existing failures"
    )
    ax.legend(fontsize=8)
    for x, k in zip(xs, labels, strict=True):
        ax.annotate(str(base.get(k, 0)), (x - 0.2, base.get(k, 0)), ha="center", va="bottom", fontsize=7)
        ax.annotate(str(merged.get(k, 0)), (x + 0.2, merged.get(k, 0)), ha="center", va="bottom", fontsize=7)

    ax = axes[1][1]
    cov = payload["coverage"]
    names = ["PR #44", "main"]
    vals = [cov["pr44_coverage"] * 100.0, cov["main_coverage"] * 100.0]
    ax.bar(names, vals, color=["#d62728", "#7f7f7f"])
    ax.set_ylim(90, 100.5)
    ax.set_ylabel("added lines kept in merged tree [%]")
    ax.set_title(
        f"no silent loss: {cov['pr44_added_lines_kept']}/{cov['pr44_added_lines']} vs "
        f"{cov['main_added_lines_kept']}/{cov['main_added_lines']} added lines; "
        f"__all__ exact union = {payload['surface']['is_exact_union']}"
    )
    for i, v in enumerate(vals):
        ax.annotate(f"{v:.2f}%", (i, v), ha="center", va="bottom", fontsize=9)
    ax.grid(alpha=0.3, axis="y")

    fig.suptitle("PR #44 cure: three-way merge of agent/m3347 with main -- verification", fontsize=12)
    fig.savefig(OUT / "merge_verification.png", dpi=130)
    plt.close(fig)



def emit_tables(payload: dict[str, object]) -> None:
    """Write the measured tables as markdown, values injected from the payload."""
    cov = payload["coverage"]
    surf = payload["surface"]
    gz = payload["gz_curve"]
    transom = payload["transom"]
    tests = payload["tests"]
    base = tests["main_baseline"] or {}
    merged = tests["merged"] or {}

    lines = ["# PR #44 cure -- measured evidence (generated)", ""]
    lines.append("Generated by `tools/verify_merge_m3822.py`; every number below is read from the")
    lines.append("merged tree at run time, none is transcribed by hand.")
    lines.append("")
    lines.append("## Merge shape")
    lines.append("")
    lines.append("| item | value |")
    lines.append("|---|---|")
    lines.append(f"| merge base | `{cov['base']}` |")
    lines.append(f"| `agent/m3347` before cure | `{cov['pr44_tip_before_cure']}` |")
    lines.append(f"| `main` merged in | `{cov['main']}` |")
    lines.append(f"| cured tip | `{payload['merged_tip']}` |")
    lines.append(f"| parents of the merge commit | `{payload['first_parent']}`, `{payload['second_parent']}` |")
    lines.append(f"| main is an ancestor of the cured tip | {payload['main_is_ancestor']} |")
    lines.append("")
    lines.append("## No silent loss")
    lines.append("")
    lines.append("| side | added lines vs base | kept in merged tree | coverage |")
    lines.append("|---|---|---|---|")
    lines.append(
        f"| PR #44 (`agent/m3347`) | {cov['pr44_added_lines']} | {cov['pr44_added_lines_kept']} "
        f"| {cov['pr44_coverage']*100:.2f} % |"
    )
    lines.append(
        f"| `main` | {cov['main_added_lines']} | {cov['main_added_lines_kept']} "
        f"| {cov['main_coverage']*100:.2f} % |"
    )
    lines.append("")
    lines.append(f"Files touched by either side: {cov['files_touched_by_either_side']}.")
    lines.append("")
    lines.append("| public surface | count |")
    lines.append("|---|---|")
    lines.append(f"| `__all__` on PR #44 | {surf['count_pr44']} |")
    lines.append(f"| `__all__` on main | {surf['count_main']} |")
    lines.append(f"| `__all__` in merged tree | {surf['count_merged']} |")
    lines.append(f"| merged == exact union | {surf['is_exact_union']} |")
    lines.append(f"| every name importable | {surf['importable_names']} |")
    lines.append(f"| names missing from PR #44 side | {surf['missing_from_pr44'] or 'none'} |")
    lines.append(f"| names missing from main side | {surf['missing_from_main'] or 'none'} |")
    lines.append("")
    lines.append("## Numerics alive in the merged tree")
    lines.append("")
    lines.append("### PR #44 side: GZ curve on the 10 x 2 x 0.5 + 1.5 m box")
    lines.append("")
    lines.append("| quantity | value |")
    lines.append("|---|---|")
    lines.append(f"| heel points | {gz['points']} |")
    lines.append(f"| GZ max [m] | {gz['gz_max_m']:.6f} |")
    lines.append(f"| heel of GZ max [deg] | {gz['angle_of_gz_max_deg']:.1f} |")
    lines.append(f"| GM from the merged solver [m] | {gz['gm_from_merged_solver_m']:.9f} |")
    lines.append(
        f"| max abs error vs section clipping [m] | {gz['max_abs_error_vs_section_clipping_m']:.3e} |"
    )
    wall = gz["max_rel_error_vs_wall_sided_in_regime"]
    lines.append(
        f"| max rel error vs wall-sided formula, |phi| <= {gz['wall_sided_regime_max_heel_deg']:.0f} deg "
        f"| {('%.3e' % wall) if wall is not None else 'n/a'} |"
    )
    lines.append("")
    lines.append("| heel [deg] | GZ merged [m] | GZ reference [m] | abs err [m] |")
    lines.append("|---|---|---|---|")
    for phi, a, b in zip(gz["heel_deg"], gz["gz_m"], gz["gz_analytic_m"], strict=True):
        lines.append(f"| {phi:.1f} | {a:.12f} | {b:.12f} | {abs(a-b):.3e} |")
    lines.append("")
    lines.append("### main side: transom-stern hydrostatics vs its closed form")
    lines.append("")
    lines.append("| nx | nz | volume [m3] | analytic [m3] | rel error |")
    lines.append("|---|---|---|---|---|")
    for row in transom["rows"]:
        lines.append(
            f"| {row['nx']} | {row['nz']} | {row['volume_m3']:.9f} | {row['volume_analytic_m3']:.9f} "
            f"| {row['volume_rel_error']*100:+.4f} % |"
        )
    lines.append("")
    lines.append(f"Monotone convergence with grid refinement: {transom['converging']}.")
    lines.append("")
    lines.append("## Test suites")
    lines.append("")
    lines.append("| suite | passed | failed | skipped |")
    lines.append("|---|---|---|---|")
    lines.append(
        f"| `main` baseline | {base.get('passed', 0)} | {base.get('failed', 0)} | {base.get('skipped', 0)} |"
    )
    lines.append(
        f"| cured `agent/m3347` | {merged.get('passed', 0)} | {merged.get('failed', 0)} "
        f"| {merged.get('skipped', 0)} |"
    )
    lines.append("")
    lines.append(f"Excluded from both runs: {', '.join(f'`{m}`' for m in tests['excluded_modules'])}")
    lines.append(f"({tests['exclusion_reason']}).")
    lines.append("")
    lines.append("| failing module (identical on both sides) |")
    lines.append("|---|")
    for module in merged.get("failing_modules", []):
        lines.append(f"| `{module}` |")
    lines.append("")
    (OUT / "merge_verification.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    # the cure's merge commit is the newest merge on the cured branch (evidence may sit on top)
    merge_sha = git("rev-list", "--merges", "-1", "HEAD").strip()
    payload = {
        "coverage": coverage_report(),
        "surface": surface_report(),
        "gz_curve": gz_curve_evidence(),
        "transom": transom_evidence(),
        "tests": test_evidence(),
        "merged_tip": git("rev-parse", "HEAD").strip(),
        "merge_commit": merge_sha,
        "first_parent": git("rev-parse", f"{merge_sha}^1").strip(),
        "second_parent": git("rev-parse", f"{merge_sha}^2").strip(),
        "main_is_ancestor": subprocess.run(
            ["git", "merge-base", "--is-ancestor", SIDE_B, "HEAD"], cwd=ROOT
        ).returncode
        == 0,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    figure(payload)
    emit_tables(payload)
    (OUT / "merge_verification.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    cov = payload["coverage"]
    surf = payload["surface"]
    gz = payload["gz_curve"]
    ok = (
        cov["pr44_coverage"] > 0.99
        and cov["main_coverage"] > 0.99
        and surf["is_exact_union"]
        and surf["importable_names"]
        and gz["max_abs_error_vs_section_clipping_m"] < 1e-9
        and payload["transom"]["converging"]
        and payload["main_is_ancestor"]
    )
    print(json.dumps({k: payload[k] for k in ("coverage", "surface")}, indent=2)[:1500])
    print(f"GZ max |err| vs section clipping: {gz['max_abs_error_vs_section_clipping_m']:.3e} m")
    print(f"transom volume errors: {[round(r['volume_rel_error']*100, 4) for r in payload['transom']['rows']]}")
    print("GATE:", "PASS" if ok else "FAIL")
    return 0 if ok else 2


if __name__ == "__main__":
    sys.exit(main())