"""Certify the IS Code areas of a barge past the wall-sided limit against the exact box solution.

    python3 tools/verify_box_gz_areas.py            # writes docs/box_gz_areas/ and syncs the doc

:mod:`tensorcad.statutory` already turns a ``GZ`` curve into the areas of IS Code 2008 part A
s. 2.2 -- 0-30, 0-40, 30-40 degrees -- and the weather criterion's gust area, and it does it with
an adaptive Simpson rule that reports its own uncertainty.  What it never had was a reference
*out there*: the closed-form antiderivative ``GM0 (1-cos) + BM0/2 (sec + cos - 2)`` stops at the
wall-sided limit, which for the doc barge is 21.8 degrees, i.e. before either regulated area is
even half finished.  The certification of that machinery therefore only ever covered the part of
the curve where the arithmetic is easiest.

``tensorcad.box_gz.box_gz_area`` supplies the missing reference: it integrates the exact
multi-regime curve, splitting the interval at the regime boundaries because ``KN'(phi)`` steps
there (the waterplane width jumps when a deck edge enters the water).  This harness checks, per
barge geometry:

1. **Inside the wall-sided window** the new integral must reproduce the closed-form
   antiderivative to round-off -- otherwise the new reference, not the code check, is suspect.
2. **The regulated areas** (0-30, 0-40, 30-40) and the wider windows (0-60, 0-89) must agree
   with :func:`tensorcad.statutory.area_under_gz_curve` on the same callable.  Two different
   quadratures, the same curve: this is what certifies the Simpson rule where no antiderivative
   exists.
3. **The mesh route** -- ``tensorcad.gz.righting_arm_curve`` off a closed-box mesh, integrated
   by the same statutory machinery from a table -- must agree to engineering tolerance, which is
   the claim a user of the mesh pipeline actually relies on.
4. **The split is load-bearing**: the same Gauss-Legendre rule without the split is reported at
   several panel counts, so the convergence order can be read off rather than asserted.
5. **Negative controls**: the wall-sided antiderivative evaluated past its limit, and a 1 mm
   ``KG`` shift, both of which *must* fail the thresholds above.  A comparison that cannot fail
   certifies nothing.

Everything comes from one run: ``docs/box_gz_areas/box_gz_areas.json``, the figure
``docs/box_gz_areas/box_gz_areas.png``, and -- through :func:`sync_doc_tables` -- the tables of
``docs/box_gz_areas.md``, replaced in place so a stale number is a diff rather than a typo.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib.pyplot as plt

from tensorcad import statutory
from tensorcad.box_gz import (
    box_gz,
    box_gz_area,
    box_regime_bounds_deg,
    box_section,
    box_vanishing_angle_deg,
)
from tensorcad.geometry import MeshBuilder
from tensorcad.gz import righting_arm_curve

ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS = ROOT / "docs" / "box_gz_areas"
DOCS = ROOT / "docs" / "box_gz_areas.md"

SPECS = [
    {"name": "doc barge", "length": 10.0, "breadth": 2.0, "draft": 0.8, "freeboard": 0.4, "kg": 0.55},
    {"name": "wide shallow", "length": 6.0, "breadth": 4.0, "draft": 0.5, "freeboard": 0.3, "kg": 0.9},
    {"name": "square deep", "length": 3.0, "breadth": 3.0, "draft": 1.0, "freeboard": 0.5, "kg": 1.2},
]

WINDOWS = ((0.0, 30.0), (0.0, 40.0), (30.0, 40.0), (0.0, 60.0), (0.0, 89.0))
ANTIDERIVATIVE_WINDOWS = (0.4, 0.8, 0.95)
UNSPLIT_PANELS = (8, 16, 32, 64, 128, 256)

# five-point Gauss-Legendre, repeated so the unsplit control is the *same* rule
NODES = (0.0, -0.5384693101056831, 0.5384693101056831, -0.9061798459386640, 0.9061798459386640)
WEIGHTS = (
    0.5688888888888889,
    0.4786286704993665,
    0.4786286704993665,
    0.2369268850561891,
    0.2369268850561891,
)

# claim thresholds, all relative
TOL_ANTIDERIVATIVE = 1e-11
TOL_STATUTORY = statutory.DEFAULT_AREA_TOLERANCE  # what the module promises
TOL_MESH = 1e-3
TOL_CONTROL = 1e-2  # a control that does not exceed this is not a control


def box_mesh(length: float, breadth: float, draft: float, freeboard: float):
    """Closed box with its bottom at ``-draft`` and its deck at ``+freeboard``."""
    builder = MeshBuilder()
    corner = {}
    for x in (0.0, length):
        for y in (-0.5 * breadth, 0.5 * breadth):
            for z in (-draft, freeboard):
                corner[(x, y, z)] = builder.add_vertex((x, y, z))

    def quad(p0, p1, p2, p3) -> None:
        a, b, c, d = (corner[point] for point in (p0, p1, p2, p3))
        builder.add_triangle(a, b, c)
        builder.add_triangle(a, c, d)

    lo, hi = -0.5 * breadth, 0.5 * breadth
    bottom, top = -draft, freeboard
    left, right = 0.0, length
    quad((left, lo, bottom), (left, hi, bottom), (right, hi, bottom), (right, lo, bottom))
    quad((left, lo, top), (right, lo, top), (right, hi, top), (left, hi, top))
    quad((left, lo, bottom), (right, lo, bottom), (right, lo, top), (left, lo, top))
    quad((left, hi, bottom), (left, hi, top), (right, hi, top), (right, hi, bottom))
    quad((left, lo, bottom), (left, lo, top), (left, hi, top), (left, hi, bottom))
    quad((right, lo, bottom), (right, hi, bottom), (right, hi, top), (right, lo, top))
    return builder.build()


def wall_sided_area(stop_deg: float, *, breadth: float, draft: float, kg: float) -> float:
    """``GM0 (1 - cos) + BM0/2 (sec + cos - 2)``, the antiderivative of the closed form."""
    bm = breadth**2 / (12.0 * draft)
    gm = 0.5 * draft + bm - kg
    phi = math.radians(stop_deg)
    return gm * (1.0 - math.cos(phi)) + 0.5 * bm * (1.0 / math.cos(phi) + math.cos(phi) - 2.0)


def unsplit_kn_area(spec: dict, depth: float, stop_deg: float, panels: int) -> float:
    """The same Gauss-Legendre rule over the whole interval: no split at the kinks."""
    low, high = 0.0, math.radians(stop_deg)
    width = (high - low) / panels
    half = 0.5 * width
    total = 0.0
    for index in range(panels):
        centre = low + (index + 0.5) * width
        for node, weight in zip(NODES, WEIGHTS, strict=True):
            total += weight * box_section(
                breadth=spec["breadth"],
                draft=spec["draft"],
                depth=depth,
                heel_deg=math.degrees(centre + half * node),
            )[2]
    return total * half


def _gz_callable(spec: dict, depth: float):
    def gz(angle: float) -> float:
        return box_gz(
            breadth=spec["breadth"],
            draft=spec["draft"],
            depth=depth,
            kg=spec["kg"],
            heel_deg=angle,
        ).righting_arm

    return gz


def run_spec(spec: dict) -> dict:
    breadth, draft = spec["breadth"], spec["draft"]
    depth = draft + spec["freeboard"]
    limits = {
        "breadth": breadth,
        "draft": draft,
        "depth": depth,
        "kg": spec["kg"],
    }
    bounds = box_regime_bounds_deg(breadth=breadth, draft=draft, depth=depth)
    wall_sided_limit = min(bounds["deck_immersion"], bounds["bilge_emergence"])

    exact = {
        f"{start:g}-{stop:g}": box_gz_area(**limits, start_deg=start, stop_deg=stop)
        for start, stop in WINDOWS
    }
    gz = _gz_callable(spec, depth)

    areas = []
    for (start, stop), result in zip(WINDOWS, (exact[f"{a:g}-{b:g}"] for a, b in WINDOWS), strict=True):
        theirs = statutory.area_under_gz_curve(gz, start, stop)
        areas.append(
            {
                "window_deg": [start, stop],
                "exact": result.area,
                "statutory": theirs,
                "relative_difference": abs(result.area - theirs) / abs(result.area),
                "quadrature_error": result.quadrature_error,
                "breakpoints_deg": list(result.breakpoints_deg),
                "heel_evaluations": result.heel_evaluations,
            }
        )

    antiderivative = []
    for fraction in ANTIDERIVATIVE_WINDOWS:
        stop = fraction * wall_sided_limit
        result = box_gz_area(**limits, stop_deg=stop)
        expected = wall_sided_area(stop, breadth=breadth, draft=draft, kg=spec["kg"])
        antiderivative.append(
            {
                "stop_deg": stop,
                "exact": result.area,
                "closed_form": expected,
                "relative_difference": abs(result.area - expected) / abs(expected),
            }
        )

    mesh = box_mesh(spec["length"], breadth, draft, spec["freeboard"])
    angles = [float(angle) for angle in range(0, 86)]
    curve = righting_arm_curve(mesh, spec["kg"], angles_deg=angles)
    table = (list(curve.angles_deg), list(curve.gz))
    mesh_areas = []
    for start, stop in ((0.0, 30.0), (0.0, 40.0), (0.0, 85.0)):
        reference = box_gz_area(**limits, start_deg=start, stop_deg=stop).area
        theirs = statutory.area_under_gz_curve(table, start, stop)
        mesh_areas.append(
            {
                "window_deg": [start, stop],
                "exact": reference,
                "mesh": theirs,
                "relative_difference": abs(reference - theirs) / abs(reference),
            }
        )

    reference_kn = box_gz_area(**limits, stop_deg=40.0, tolerance=1e-15).kn_area
    convergence = []
    for panels in UNSPLIT_PANELS:
        value = unsplit_kn_area(spec, depth, 40.0, panels)
        convergence.append(
            {
                "panels": panels,
                "heel_evaluations": panels * len(NODES),
                "absolute_error": abs(value - reference_kn),
            }
        )
    split = box_gz_area(**limits, stop_deg=40.0)

    # negative controls: both must break the thresholds above
    control_wall_sided = wall_sided_area(
        40.0, breadth=breadth, draft=draft, kg=spec["kg"]
    )
    control_kg = box_gz_area(**{**limits, "kg": spec["kg"] + 1e-3}, stop_deg=40.0)

    curve_angles = [float(angle) for angle in range(0, 90)]
    curve_points = [
        box_gz(breadth=breadth, draft=draft, depth=depth, kg=spec["kg"], heel_deg=angle)
        for angle in curve_angles
    ]

    return {
        "spec": spec,
        "depth": depth,
        "regime_bounds_deg": bounds,
        "wall_sided_limit_deg": wall_sided_limit,
        "vanishing_angle_deg": box_vanishing_angle_deg(**limits),
        "areas": areas,
        "antiderivative": antiderivative,
        "mesh_areas": mesh_areas,
        "convergence_unsplit": convergence,
        "split": {
            "area_0_40": split.area,
            "kn_area_0_40": split.kn_area,
            "panels": split.panels,
            "heel_evaluations": split.heel_evaluations,
            "quadrature_error": split.quadrature_error,
            "split_error": abs(split.kn_area - reference_kn),
            "breakpoints_deg": list(split.breakpoints_deg),
        },
        "controls": {
            "wall_sided_past_limit_area_0_40": control_wall_sided,
            "wall_sided_past_limit_relative_difference": abs(
                control_wall_sided - exact["0-40"].area
            )
            / abs(exact["0-40"].area),
            "kg_plus_1mm_absolute_change": abs(control_kg.area - exact["0-40"].area),
            "kg_plus_1mm_relative_change": abs(control_kg.area - exact["0-40"].area)
            / abs(exact["0-40"].area),
        },
        "curve": {
            "angles_deg": curve_angles,
            "gz": [point.righting_arm for point in curve_points],
            "regimes": [point.regime for point in curve_points],
            "cumulative_area": [
                box_gz_area(**limits, stop_deg=angle).area for angle in curve_angles[1:]
            ],
        },
    }


def plot(results: list[dict]) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))
    left, right = axes

    main = results[0]
    curve = main["curve"]
    angles = curve["angles_deg"]
    left.plot(angles, curve["gz"], "k-", lw=2.0, label="exact GZ, four regimes")
    twin = left.twinx()
    twin.plot(
        angles[1:],
        curve["cumulative_area"],
        "C3-",
        lw=1.8,
        label=r"$\int_0^{\phi} GZ\,d\phi$",
    )
    for name, bound in main["regime_bounds_deg"].items():
        if math.isfinite(bound) and bound < 89.0:
            left.axvline(bound, color="C0", ls=":", lw=1.0)
            left.annotate(
                name.replace("_", " "),
                xy=(bound, 0.02),
                xytext=(bound + 1.0, 0.02),
                rotation=90,
                fontsize=7,
                color="C0",
                va="bottom",
            )
    left.axhline(0.0, color="0.6", lw=0.8)
    left.set_xlabel("heel [deg]")
    left.set_ylabel("GZ [m]")
    twin.set_ylabel(r"area [m rad]")
    left.set_title(f"{main['spec']['name']}: the curve and the cumulative area")
    left.grid(alpha=0.3)

    markers = {"doc barge": "o", "wide shallow": "s", "square deep": "^"}
    for result in results:
        rows = result["convergence_unsplit"]
        right.loglog(
            [row["heel_evaluations"] for row in rows],
            [row["absolute_error"] for row in rows],
            marker=markers.get(result["spec"]["name"], "o"),
            ms=4,
            mfc="none",
            label=f"{result['spec']['name']}: unsplit rule",
        )
        right.axhline(
            max(result["split"]["split_error"], 1e-18),
            color="0.5",
            ls="--",
            lw=1.0,
        )
    right.set_xlabel("section evaluations")
    right.set_ylabel(r"error in $\int KN\,d\phi$ over 0-40 deg [m rad]")
    right.set_title("why the interval is split at the regime bounds")
    right.grid(alpha=0.3, which="both")
    right.legend(fontsize=7, loc="lower left")

    figure.tight_layout()
    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    figure.savefig(ARTIFACTS / "box_gz_areas.png", dpi=140)
    plt.close(figure)


def render_tables(results: list[dict]) -> dict[str, str]:
    """The markdown blocks, from the numbers this run produced."""
    benchmark = ["| barge | bounds, deg (deck / bilge) | A(0-30) m rad | A(0-40) m rad | A(30-40) m rad |"]
    benchmark.append("|---|---|---|---|---|")
    for result in results:
        by_window = {tuple(row["window_deg"]): row for row in result["areas"]}
        bounds = result["regime_bounds_deg"]
        benchmark.append(
            "| {name} | {deck:.4f} / {bilge:.4f} | {a30:.6f} | {a40:.6f} | {a3040:.6f} |".format(
                name=result["spec"]["name"],
                deck=bounds["deck_immersion"],
                bilge=bounds["bilge_emergence"],
                a30=by_window[(0.0, 30.0)]["exact"],
                a40=by_window[(0.0, 40.0)]["exact"],
                a3040=by_window[(30.0, 40.0)]["exact"],
            )
        )

    residuals = [
        "| barge | window, deg | exact m rad | statutory | relative | quadrature error | breakpoints |",
        "|---|---|---|---|---|---|---|",
    ]
    for result in results:
        for row in result["areas"]:
            start, stop = row["window_deg"]
            residuals.append(
                "| {name} | {start:g}-{stop:g} | {exact:.12f} | {theirs:.12f} | {rel:.2e} | "
                "{err:.1e} | {breaks} |".format(
                    name=result["spec"]["name"],
                    start=start,
                    stop=stop,
                    exact=row["exact"],
                    theirs=row["statutory"],
                    rel=row["relative_difference"],
                    err=row["quadrature_error"],
                    breaks=", ".join(f"{angle:.3f}" for angle in row["breakpoints_deg"]) or "--",
                )
            )

    interior = [
        "| barge | stop, deg | exact | closed-form antiderivative | relative |",
        "|---|---|---|---|---|",
    ]
    for result in results:
        for row in result["antiderivative"]:
            interior.append(
                "| {name} | {stop:.4f} | {exact:.12f} | {closed:.12f} | {rel:.2e} |".format(
                    name=result["spec"]["name"],
                    stop=row["stop_deg"],
                    exact=row["exact"],
                    closed=row["closed_form"],
                    rel=row["relative_difference"],
                )
            )

    mesh = [
        "| barge | window, deg | exact | mesh curve | relative |",
        "|---|---|---|---|---|",
    ]
    for result in results:
        for row in result["mesh_areas"]:
            start, stop = row["window_deg"]
            mesh.append(
                "| {name} | {start:g}-{stop:g} | {exact:.10f} | {theirs:.10f} | {rel:.2e} |".format(
                    name=result["spec"]["name"],
                    start=start,
                    stop=stop,
                    exact=row["exact"],
                    theirs=row["mesh"],
                    rel=row["relative_difference"],
                )
            )

    convergence = [
        "| barge | panels | section evaluations | error, unsplit | error, split at the bounds |",
        "|---|---|---|---|---|",
    ]
    for result in results:
        split_row = result["split"]
        for row in result["convergence_unsplit"]:
            convergence.append(
                "| {name} | {panels} | {evals} | {err:.2e} | -- |".format(
                    name=result["spec"]["name"],
                    panels=row["panels"],
                    evals=row["heel_evaluations"],
                    err=row["absolute_error"],
                )
            )
        convergence.append(
            "| {name} | {panels} | {evals} | -- | {err:.2e} |".format(
                name=result["spec"]["name"],
                panels=split_row["panels"],
                evals=split_row["heel_evaluations"],
                err=split_row["split_error"],
            )
        )

    controls = ["| control | value | against the threshold |", "|---|---|---|"]
    for result in results:
        row = result["controls"]
        controls.append(
            "| {name}: wall-sided antiderivative used past its limit, A(0-40) | {rel:.2e} | "
            "must exceed {tol:.0e} |".format(
                name=result["spec"]["name"],
                rel=row["wall_sided_past_limit_relative_difference"],
                tol=TOL_CONTROL,
            )
        )
        controls.append(
            "| {name}: KG + 1 mm, A(0-40) | {rel:.2e} relative ({abs_change:.3e} m rad) | "
            "must be non-zero |".format(
                name=result["spec"]["name"],
                rel=row["kg_plus_1mm_relative_change"],
                abs_change=row["kg_plus_1mm_absolute_change"],
            )
        )

    return {
        "benchmark": "\n".join(benchmark),
        "residuals": "\n".join(residuals),
        "interior": "\n".join(interior),
        "mesh": "\n".join(mesh),
        "convergence": "\n".join(convergence),
        "controls": "\n".join(controls),
    }


def sync_doc_tables(blocks: dict[str, str]) -> bool:
    """Replace each ``<!-- begin: name -->`` .. ``<!-- end: name -->`` block in place."""
    if not DOCS.exists():
        return False
    text = DOCS.read_text()
    changed = False
    for name, block in blocks.items():
        begin = f"<!-- begin: {name} -->"
        end = f"<!-- end: {name} -->"
        head = text.index(begin) + len(begin)
        tail = text.index(end)
        current = text[head:tail]
        replacement = "\n\n" + block + "\n\n"
        if current != replacement:
            text = text[:head] + replacement + text[tail:]
            changed = True
    DOCS.write_text(text)
    return changed


def claims(results: list[dict]) -> list[tuple[str, bool, str]]:
    """Every claim as (name, passed, detail); a failing one exits non-zero."""
    verdict: list[tuple[str, bool, str]] = []
    worst_antiderivative = max(
        row["relative_difference"] for result in results for row in result["antiderivative"]
    )
    verdict.append(
        (
            "inside the wall-sided window the integral is the closed-form antiderivative",
            worst_antiderivative <= TOL_ANTIDERIVATIVE,
            f"worst relative difference {worst_antiderivative:.2e} <= {TOL_ANTIDERIVATIVE:.0e}",
        )
    )
    worst_statutory = max(row["relative_difference"] for result in results for row in result["areas"])
    verdict.append(
        (
            "statutory's adaptive Simpson agrees past the wall-sided limit, within its own tolerance",
            worst_statutory <= TOL_STATUTORY,
            f"worst relative difference {worst_statutory:.2e} <= {TOL_STATUTORY:.0e}",
        )
    )
    worst_mesh = max(row["relative_difference"] for result in results for row in result["mesh_areas"])
    verdict.append(
        (
            "the heeled-mesh curve agrees to engineering tolerance",
            worst_mesh <= TOL_MESH,
            f"worst relative difference {worst_mesh:.2e} <= {TOL_MESH:.0e}",
        )
    )
    worst_control = min(
        result["controls"]["wall_sided_past_limit_relative_difference"] for result in results
    )
    verdict.append(
        (
            "the negative control breaks the threshold it is compared against",
            worst_control >= TOL_CONTROL,
            f"smallest control error {worst_control:.2e} >= {TOL_CONTROL:.0e}",
        )
    )
    smallest_kg = min(result["controls"]["kg_plus_1mm_relative_change"] for result in results)
    verdict.append(
        (
            "a millimetre of KG is visible in the area",
            smallest_kg > 0.0,
            f"smallest relative change {smallest_kg:.2e}",
        )
    )
    unsplit_gain = []
    for result in results:
        rows = result["convergence_unsplit"]
        unsplit_gain.append(rows[-1]["absolute_error"] / max(result["split"]["split_error"], 1e-18))
    verdict.append(
        (
            "the split is what makes the quadrature converge",
            min(unsplit_gain) > 10.0,
            "unsplit error / split error at the finest panel count: "
            + ", ".join(f"{gain:.0f}x" for gain in unsplit_gain),
        )
    )
    return verdict


def main() -> int:
    results = [run_spec(spec) for spec in SPECS]
    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    payload = {
        "harness": "tools/verify_box_gz_areas.py",
        "thresholds": {
            "antiderivative": TOL_ANTIDERIVATIVE,
            "statutory": TOL_STATUTORY,
            "mesh": TOL_MESH,
            "control": TOL_CONTROL,
        },
        "results": results,
    }
    (ARTIFACTS / "box_gz_areas.json").write_text(json.dumps(payload, indent=2) + "\n")
    plot(results)
    sync_doc_tables(render_tables(results))

    verdict = claims(results)
    width = max(len(name) for name, _, _ in verdict)
    for name, passed, detail in verdict:
        print(f"{'PASS' if passed else 'FAIL'}  {name:<{width}}  {detail}")
    print(f"\nwrote {ARTIFACTS.relative_to(ROOT)}/box_gz_areas.json and .png")
    return 0 if all(passed for _, passed, _ in verdict) else 1


if __name__ == "__main__":
    raise SystemExit(main())
