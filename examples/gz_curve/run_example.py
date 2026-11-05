"""Righting arms at finite heel: ``GZ(phi)`` checked against closed forms.

`tensorcad.heel` integrates the submerged body below a **tilted** water surface at
constant displacement. Nothing sampled, nothing assumed: the surface is clipped against
the half-space and the signed tetrahedron sum is exact for the polyhedron, so the only
error left is the one-dimensional root-find on the waterline offset. This script measures
where that claim holds and where the *references* stop holding, and exits non-zero if a
claim breaks.

1. **A rectangular barge.** The wall-sided formula
   `GZ = sin(phi) (GM + BM tan^2(phi) / 2)` is exact until the deck edge immerses; the
   mesh integral must reproduce it to round-off below that angle (`L = 20, B = 6,
   T = 2, KG = 2`, deck edge at 33.7 deg, wall-sided reference used to 30 deg).
2. **A shallow-freeboard barge**, where the deck edge immerses at 18.4 deg and the
   curve actually turns over: maximum righting arm and range of stability are reported,
   and the wall-sided formula is shown to *over-predict* past the deck edge instead of
   being quietly extended.
3. **A sphere**, where heel changes nothing (`GZ = sin(phi) (z_C - z_G)` at any angle,
   no deck to immerse). The residual there is the facet error of the mesh, and a
   refinement ladder shows it falling like `h^2` rather than asserting it is small.
4. **Consistency with `tensorcad.stability`**: `dGZ/dphi` at zero heel must equal the
   `GM` built from the waterplane moments, and `GZ` must be odd in `phi`.

Run it with ``python examples/gz_curve/run_example.py --output docs/gz_curve``.
matplotlib is optional: the JSON is written either way and the PNGs are skipped.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

from tensorcad import __version__
from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.heel import gz_curve, wall_sided_gz

BARGE = {"length": 20.0, "beam": 6.0, "draft": 2.0, "depth": 4.0, "kg": 2.0}
SHALLOW = {"length": 20.0, "beam": 6.0, "draft": 2.0, "depth": 3.0, "kg": 2.2}
SPHERE = {"radius": 1.0, "draft_fraction": 0.6, "kg": 0.3}


# ---------------------------------------------------------------------------
# reference bodies
# ---------------------------------------------------------------------------
def box(length: float, beam: float, draft: float, depth: float) -> Mesh:
    builder = MeshBuilder()
    corner = {}
    for x in (0.0, length):
        for y in (-0.5 * beam, 0.5 * beam):
            for z in (-draft, depth - draft):
                corner[(x, y, z)] = builder.add_vertex((x, y, z))

    def quad(*points: tuple[float, float, float]) -> None:
        a, b, c, d = (corner[point] for point in points)
        builder.add_triangle(a, b, c)
        builder.add_triangle(a, c, d)

    lo, hi = -0.5 * beam, 0.5 * beam
    keel, deck = -draft, depth - draft
    quad((0.0, lo, keel), (0.0, hi, keel), (length, hi, keel), (length, lo, keel))
    quad((0.0, lo, deck), (length, lo, deck), (length, hi, deck), (0.0, hi, deck))
    quad((0.0, lo, keel), (length, lo, keel), (length, lo, deck), (0.0, lo, deck))
    quad((0.0, hi, keel), (0.0, hi, deck), (length, hi, deck), (length, hi, keel))
    quad((0.0, lo, keel), (0.0, lo, deck), (0.0, hi, deck), (0.0, hi, keel))
    quad((length, lo, keel), (length, hi, keel), (length, hi, deck), (length, lo, deck))
    return builder.build()


def sphere(radius: float, centre_z: float, *, rings: int, segments: int) -> Mesh:
    builder = MeshBuilder()
    north = builder.add_vertex((0.0, 0.0, centre_z + radius))
    bands = []
    for ring in range(1, rings):
        theta = math.pi * ring / rings
        band = []
        for segment in range(segments):
            phi = 2.0 * math.pi * segment / segments
            band.append(
                builder.add_vertex(
                    (
                        radius * math.sin(theta) * math.cos(phi),
                        radius * math.sin(theta) * math.sin(phi),
                        centre_z + radius * math.cos(theta),
                    )
                )
            )
        bands.append(band)
    south = builder.add_vertex((0.0, 0.0, centre_z - radius))
    for segment in range(segments):
        nxt = (segment + 1) % segments
        builder.add_triangle(north, bands[0][segment], bands[0][nxt])
    for ring in range(len(bands) - 1):
        lower, upper = bands[ring], bands[ring + 1]
        for segment in range(segments):
            nxt = (segment + 1) % segments
            builder.add_triangle(lower[segment], upper[nxt], lower[nxt])
            builder.add_triangle(lower[segment], upper[segment], upper[nxt])
    for segment in range(segments):
        nxt = (segment + 1) % segments
        builder.add_triangle(south, bands[-1][nxt], bands[-1][segment])
    return builder.build()


def box_reference(case: dict[str, float]) -> tuple[float, float, float]:
    """``(KB, BM, deck_edge_deg)`` of the rectangular barge at zero heel."""
    beam, draft = case["beam"], case["draft"]
    return 0.5 * draft, beam**2 / (12.0 * draft), math.degrees(
        math.atan(2.0 * (case["depth"] - case["draft"]) / beam)
    )


# ---------------------------------------------------------------------------
# measurements
# ---------------------------------------------------------------------------
def measure_barge(case: dict[str, float], *, angles: tuple[float, ...]) -> dict[str, object]:
    mesh = box(case["length"], case["beam"], case["draft"], case["depth"])
    kb, bm, deck_edge = box_reference(case)
    curve = gz_curve(mesh, kg=case["kg"], angles_deg=angles)
    rows = []
    for section in curve.sections:
        reference = wall_sided_gz(kb=kb, bm=bm, kg=case["kg"], heel_deg=section.heel_deg)
        rows.append(
            {
                "heel_deg": section.heel_deg,
                "gz_m": section.gz,
                "wall_sided_m": reference,
                "wall_sided_valid": section.heel_deg <= deck_edge,
                "relative_error": (
                    (section.gz - reference) / reference if reference != 0.0 else 0.0
                ),
                "waterline_offset_m": section.waterline_offset,
                "displacement_volume_m3": section.displacement_volume,
            }
        )
    valid = [row for row in rows if row["wall_sided_valid"] and row["wall_sided_m"] != 0.0]
    worst = max((abs(row["relative_error"]) for row in valid), default=0.0)
    worst_offset = max((abs(row["waterline_offset_m"]) for row in valid), default=0.0)
    # the wall-sided formula has a tan^2 pole at 90 deg, so quote the over-prediction at
    # the largest angle where it is still a finite number
    capped = [row for row in rows if not row["wall_sided_valid"] and row["heel_deg"] <= 60.0]
    excess_row = capped[-1] if capped else None
    max_angle, max_gz = curve.max_gz
    return {
        "parameters": dict(case),
        "kb_m": kb,
        "bm_m": bm,
        "gm_m": curve.gm_transverse,
        "gm_closed_form": kb + bm - case["kg"],
        "displacement_volume_m3": curve.displacement_volume,
        "deck_edge_deg": deck_edge,
        "max_gz_m": max_gz,
        "angle_of_max_gz_deg": max_angle,
        "vanishing_angle_deg": curve.vanishing_angle_deg,
        "rows": rows,
        "worst_relative_error_vs_wall_sided": worst,
        "worst_waterline_offset_m": worst_offset,
        "excess_over_wall_sided_deg": None if excess_row is None else excess_row["heel_deg"],
        "excess_over_wall_sided_m": (
            0.0 if excess_row is None else excess_row["wall_sided_m"] - excess_row["gz_m"]
        ),
    }


def measure_sphere(*, rings: int = 24) -> dict[str, object]:
    radius = SPHERE["radius"]
    centre_z = radius - SPHERE["draft_fraction"] * radius
    keel_z = centre_z - radius
    kg = SPHERE["kg"]
    mesh = sphere(radius, centre_z, rings=rings, segments=2 * rings)
    curve = gz_curve(mesh, kg=kg, angles_deg=(0.0, 10.0, 20.0, 40.0, 60.0, 80.0))
    rows = []
    for section in curve.sections:
        reference = math.sin(math.radians(section.heel_deg)) * (centre_z - (keel_z + kg))
        rows.append(
            {
                "heel_deg": section.heel_deg,
                "gz_m": section.gz,
                "closed_form_m": reference,
                "error_m": section.gz - reference,
                "waterline_offset_m": section.waterline_offset,
                "displacement_volume_m3": section.displacement_volume,
            }
        )
    return {
        "parameters": dict(SPHERE),
        "rings": rings,
        "facets": 2 * rings * 2 * (rings - 1),
        "centre_z_m": centre_z,
        "keel_z_m": keel_z,
        "gm_closed_form_m": centre_z - (keel_z + kg),
        "gm_measured_m": curve.gm_transverse,
        "displacement_volume_m3": curve.displacement_volume,
        "rows": rows,
        "worst_absolute_error_m": max(abs(row["error_m"]) for row in rows),
    }


def measure_sphere_ladder(heel_deg: float = 60.0) -> list[dict[str, float]]:
    radius = SPHERE["radius"]
    centre_z = radius - SPHERE["draft_fraction"] * radius
    keel_z = centre_z - radius
    reference = math.sin(math.radians(heel_deg)) * (centre_z - (keel_z + SPHERE["kg"]))
    ladder = []
    for rings in (12, 24, 48):
        mesh = sphere(radius, centre_z, rings=rings, segments=2 * rings)
        curve = gz_curve(mesh, kg=SPHERE["kg"], angles_deg=(heel_deg,))
        error = curve.sections[0].gz - reference
        ladder.append(
            {
                "rings": rings,
                "segments": 2 * rings,
                "step_rad": math.pi / rings,
                "gz_m": curve.sections[0].gz,
                "reference_m": reference,
                "absolute_error_m": abs(error),
                "signed_error_m": error,
            }
        )
    for index in range(1, len(ladder)):
        previous = ladder[index - 1]["absolute_error_m"]
        ladder[index]["halving_factor"] = previous / ladder[index]["absolute_error_m"]
    return ladder


def measure_gm_consistency(case: dict[str, float]) -> dict[str, float]:
    """``dGZ/dphi`` at zero heel against the ``GM`` built from the waterplane moments."""
    mesh = box(case["length"], case["beam"], case["draft"], case["depth"])
    curve = gz_curve(mesh, kg=case["kg"], angles_deg=(0.0, 0.5, 1.0))
    slope = (curve.at(1.0) - curve.at(0.0)) / math.radians(1.0)
    kb, bm, _ = box_reference(case)
    return {
        "chord_slope_over_one_degree_m": slope,
        "gm_from_waterplane_m": curve.gm_transverse,
        "gm_closed_form_m": kb + bm - case["kg"],
        "relative_difference": (slope - curve.gm_transverse) / curve.gm_transverse,
    }


def measure_odd_symmetry(case: dict[str, float]) -> float:
    """Worst ``|GZ(phi) + GZ(-phi)|`` on a symmetric body."""
    mesh = box(case["length"], case["beam"], case["draft"], case["depth"])
    angles = (10.0, 20.0, 30.0)
    positive = gz_curve(mesh, kg=case["kg"], angles_deg=angles)
    negative = gz_curve(mesh, kg=case["kg"], angles_deg=tuple(-a for a in angles))
    return max(
        abs(plus.gz + minus.gz)
        for plus, minus in zip(positive.sections, reversed(negative.sections))
    )


# ---------------------------------------------------------------------------
# claims
# ---------------------------------------------------------------------------
def claim(text: str, value: float, threshold: float, *, below: bool = True) -> dict[str, object]:
    holds = abs(value) < threshold if below else abs(value) > threshold
    return {
        "claim": text,
        "value": value,
        "threshold": threshold,
        "direction": "<" if below else ">",
        "holds": bool(holds),
    }


def evaluate(results: dict[str, object]) -> list[dict[str, object]]:
    barge = results["barge"]
    shallow = results["shallow"]
    sphere = results["sphere"]
    ladder = results["sphere_ladder"]
    gm = results["gm_consistency"]
    symmetry = results["odd_symmetry"]
    claims = [
        claim(
            "barge GZ vs the wall-sided formula below the deck edge (worst relative)",
            barge["worst_relative_error_vs_wall_sided"],
            1e-9,
        ),
        claim(
            "barge waterline offset stays on z = 0 while the sides are wall-sided",
            barge["worst_waterline_offset_m"],
            1e-9,
        ),
        claim("barge GM from the waterplane moments vs KB + BM - KG", barge["gm_m"] - barge["gm_closed_form"], 1e-12),
        claim("GZ(0) = 0", abs(barge["rows"][0]["gz_m"]), 1e-15),
        claim("odd symmetry |GZ(phi) + GZ(-phi)|", symmetry, 1e-12),
        claim("small-angle chord slope vs GM (relative)", gm["relative_difference"], 2e-3),
        claim(
            "sphere GZ vs sin(phi)(z_C - z_G) at 24 rings (worst |error|, m)",
            sphere["worst_absolute_error_m"],
            2e-3,
        ),
        claim("sphere refinement: error falls at least 3x per doubling (48 vs 24 rings)", ladder[2]["halving_factor"], 3.0, below=False),
        claim("sphere 48-ring error (m)", ladder[2]["absolute_error_m"], 5e-4),
        claim(
            "constant displacement held at every sampled angle (relative, both bodies)",
            results["worst_volume_residual"],
            1e-9,
        ),
    ]
    claims.append(
        claim(
            "shallow barge has a finite range of stability (deg)",
            float(shallow["vanishing_angle_deg"] or 0.0),
            10.0,
            below=False,
        )
    )
    claims.append(
        claim(
            "wall-sided formula over-predicts past the deck edge (m, largest sample <= 60 deg)",
            shallow["excess_over_wall_sided_m"],
            0.5,
            below=False,
        )
    )
    return claims


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------
def figure_barge(results: dict[str, object], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axis = plt.subplots(figsize=(7.4, 4.4))
    for key, colour, label in (
        ("barge", "#1f77b4", "barge, freeboard 2 m, KG 2.0"),
        ("shallow", "#d62728", "barge, freeboard 1 m, KG 2.2"),
    ):
        case = results[key]
        heels = [row["heel_deg"] for row in case["rows"]]
        arms = [row["gz_m"] for row in case["rows"]]
        axis.plot(heels, arms, "-o", markersize=3, color=colour, label=f"{label}: GZ from heeled integration")
        valid = [(row["heel_deg"], row["wall_sided_m"]) for row in case["rows"] if row["wall_sided_valid"]]
        axis.plot(
            [point[0] for point in valid],
            [point[1] for point in valid],
            "--",
            linewidth=1.0,
            color=colour,
            label=f"{label}: wall-sided closed form (valid to {case['deck_edge_deg']:.1f}°)",
        )
        axis.axvline(case["deck_edge_deg"], color=colour, linewidth=0.6, linestyle=":", alpha=0.7)
        angle, arm = case["angle_of_max_gz_deg"], case["max_gz_m"]
        axis.plot([angle], [arm], "k*", markersize=8)
        if case["vanishing_angle_deg"] is not None:
            axis.plot([case["vanishing_angle_deg"]], [0.0], "kx", markersize=7)
            axis.annotate(
                f"vanishing {case['vanishing_angle_deg']:.1f}°",
                (case["vanishing_angle_deg"], 0.0),
                textcoords="offset points",
                xytext=(4, -12),
                fontsize=7,
                color=colour,
            )
    axis.axhline(0.0, color="k", linewidth=0.6)
    axis.set_xlabel("heel angle φ [deg]")
    axis.set_ylabel("righting arm GZ [m]")
    axis.set_title("Barge GZ: heeled integration vs the wall-sided formula, with deck-edge limits")
    axis.legend(fontsize=6.5, loc="upper right")
    axis.grid(alpha=0.3)
    figure.tight_layout()
    figure.savefig(path, dpi=140)
    plt.close(figure)


def figure_sphere(results: dict[str, object], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sphere = results["sphere"]
    heels = [row["heel_deg"] for row in sphere["rows"]]
    figure, (axis, residual) = plt.subplots(1, 2, figsize=(9.2, 3.8))
    axis.plot(heels, [row["gz_m"] for row in sphere["rows"]], "-o", markersize=3, label="heeled integration")
    axis.plot(
        [row["heel_deg"] for row in sphere["rows"]],
        [row["closed_form_m"] for row in sphere["rows"]],
        "k--",
        linewidth=1.0,
        label="sin φ (z$_C$ − z$_G$)",
    )
    axis.set_xlabel("heel angle φ [deg]")
    axis.set_ylabel("GZ [m]")
    axis.set_title(f"Sphere (R = {SPHERE['radius']:.1f} m, draft {SPHERE['draft_fraction']:.1f} R)")
    axis.legend(fontsize=8)
    axis.grid(alpha=0.3)
    residual.plot(heels, [1e3 * abs(row["error_m"]) for row in sphere["rows"]], "-s", markersize=3, color="#d62728")
    residual.set_xlabel("heel angle φ [deg]")
    residual.set_ylabel("|GZ − closed form| [mm]")
    residual.set_title("Residual is the facet error of this 24-ring mesh")
    residual.grid(alpha=0.3)
    figure.tight_layout()
    figure.savefig(path, dpi=140)
    plt.close(figure)


def figure_ladder(results: dict[str, object], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ladder = results["sphere_ladder"]
    steps = [row["step_rad"] for row in ladder]
    errors = [row["absolute_error_m"] for row in ladder]
    figure, axis = plt.subplots(figsize=(6.2, 4.0))
    axis.loglog(steps, errors, "-o", label="|GZ − sin φ (z$_C$ − z$_G$)| at 60°")
    axis.loglog(
        steps,
        [errors[0] * (step / steps[0]) ** 2 for step in steps],
        "k--",
        linewidth=0.9,
        label="slope 2 reference",
    )
    axis.set_xlabel("facet size h = π / rings [rad]")
    axis.set_ylabel("|absolute error| [m]")
    axis.invert_xaxis()
    axis.set_title("Sphere: the residual is mesh refinement, at second order")
    axis.legend(fontsize=8)
    axis.grid(alpha=0.3, which="both")
    figure.tight_layout()
    figure.savefig(path, dpi=140)
    plt.close(figure)


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
def render_report(results: dict[str, object]) -> str:
    barge = results["barge"]
    shallow = results["shallow"]
    sphere = results["sphere"]
    ladder = results["sphere_ladder"]
    gm = results["gm_consistency"]
    lines = [
        "# Righting arms at finite heel: GZ(φ) against closed forms",
        "",
        f"TensorCAD {results['version']} — `tensorcad.heel`. Generated by "
        "`python examples/gz_curve/run_example.py`; every number below is read from "
        "`results.json`.",
        "",
        "## 1. What is measured",
        "",
        "The submerged body below a tilted plane is built by clipping every triangle "
        "against `u·r ≤ c` and summing signed tetrahedra with the apex **on the water "
        "surface**, so the waterplane lid contributes exactly zero. Volume and first "
        "moment are therefore exact for the polyhedron; the only remaining error is the "
        "1-D solve for `c` that holds the displacement constant.",
        "",
        "- barge: "
        f"L = {BARGE['length']:.1f} m, B = {BARGE['beam']:.1f} m, T = {BARGE['draft']:.1f} m, "
        f"depth = {BARGE['depth']:.1f} m, KG = {BARGE['kg']:.2f} m; deck edge at "
        f"{barge['deck_edge_deg']:.2f}°; ∇ = {barge['displacement_volume_m3']:.4f} m³.",
        "- shallow barge: same L/B/T, depth "
        f"{SHALLOW['depth']:.1f} m, KG = {SHALLOW['kg']:.2f} m; deck edge at "
        f"{shallow['deck_edge_deg']:.2f}°.",
        "- sphere: R = "
        f"{SPHERE['radius']:.2f} m, draft = {SPHERE['draft_fraction']:.1f} R, KG = "
        f"{SPHERE['kg']:.2f} m, {sphere['rings']}-ring mesh ({sphere['facets']} facets).",
        "",
        "| body | GM from waterplane moments [m] | closed form [m] | max GZ [m] @ deg | vanishing angle [deg] |",
        "|---|---|---|---|---|",
    ]
    for label, case in (("barge", barge), ("shallow barge", shallow)):
        vanishing = case["vanishing_angle_deg"]
        lines.append(
            f"| {label} | {case['gm_m']:.6f} | {case['gm_closed_form']:.6f} | "
            f"{case['max_gz_m']:.4f} @ {case['angle_of_max_gz_deg']:.0f} | "
            f"{'—' if vanishing is None else f'{vanishing:.2f}'} |"
        )
    lines += [
        f"| sphere | {sphere['gm_measured_m']:.6f} | {sphere['gm_closed_form_m']:.6f} | "
        f"— | — |",
        "",
        "## 2. Barge against the wall-sided formula",
        "",
        "| φ [deg] | GZ mesh [m] | wall-sided [m] | reference valid | relative error | waterline offset c [m] |",
        "|---|---|---|---|---|---|",
    ]
    for row in barge["rows"]:
        lines.append(
            f"| {row['heel_deg']:.0f} | {row['gz_m']:+.9f} | {row['wall_sided_m']:+.9f} | "
            f"{'yes' if row['wall_sided_valid'] else 'no'} | {row['relative_error']:+.2e} | "
            f"{row['waterline_offset_m']:+.2e} |"
        )
    lines += [
        "",
        f"Worst relative deviation below the deck edge: "
        f"**{barge['worst_relative_error_vs_wall_sided']:.2e}**, and the solved waterline "
        f"offset never leaves the upright waterline by more than "
        f"{barge['worst_waterline_offset_m']:.2e} m — which is the wall-sided property "
        "(`V` linear in `c`) coming out of the solve rather than being assumed.",
        "",
        "## 3. Where the reference stops holding: the shallow barge",
        "",
        "| φ [deg] | GZ mesh [m] | wall-sided [m] | reference valid |",
        "|---|---|---|---|",
    ]
    for row in shallow["rows"]:
        lines.append(
            f"| {row['heel_deg']:.0f} | {row['gz_m']:+.6f} | {row['wall_sided_m']:+.6f} | "
            f"{'yes' if row['wall_sided_valid'] else 'no'} |"
        )
    lines += [
        "",
        f"Past the deck edge ({shallow['deck_edge_deg']:.2f}°) the unbounded-wall formula "
        f"over-predicts by {shallow['excess_over_wall_sided_m']:.4f} m at "
        f"{shallow['excess_over_wall_sided_deg']:.0f}°, and the measured curve turns over: "
        f"GZ_max = {shallow['max_gz_m']:.4f} m at {shallow['angle_of_max_gz_deg']:.0f}°, "
        f"vanishing at {shallow['vanishing_angle_deg']:.2f}°. The integration is the one "
        "to trust there; the formula is not a reference for an immersed deck.",
        "",
        "## 4. Sphere: an exact large-angle reference with no deck",
        "",
        "| φ [deg] | GZ mesh [m] | sin φ (z_C − z_G) [m] | error [mm] | waterline offset [m] |",
        "|---|---|---|---|---|",
    ]
    for row in sphere["rows"]:
        lines.append(
            f"| {row['heel_deg']:.0f} | {row['gz_m']:+.9f} | {row['closed_form_m']:+.9f} | "
            f"{1e3 * row['error_m']:+.3f} | {row['waterline_offset_m']:+.4f} |"
        )
    lines += [
        "",
        f"Worst |error| {1e3 * sphere['worst_absolute_error_m']:.3f} mm. This residual is "
        "**the faceting of the mesh**, not the integration — the waterline offset moves "
        "by "
        f"{abs(sphere['rows'][-1]['waterline_offset_m'] - sphere['rows'][0]['waterline_offset_m']):.3f} m "
        "between 0° and 80° while the displacement stays fixed, so the solve is doing real "
        "work, and the ladder below removes the doubt:",
        "",
        "| rings | facets | h = π/rings | GZ at 60° [m] | error [m] | factor per doubling |",
        "|---|---|---|---|---|---|",
    ]
    for row in ladder:
        factor = row.get("halving_factor")
        lines.append(
            f"| {row['rings']} | {2 * row['rings'] * 2 * (row['rings'] - 1)} | "
            f"{row['step_rad']:.4f} | {row['gz_m']:.9f} | {row['absolute_error_m']:.3e} | "
            f"{'—' if factor is None else f'{factor:.2f}×'} |"
        )
    lines += [
        "",
        "## 5. Consistency with `tensorcad.stability`",
        "",
        f"- 1° chord slope {gm['chord_slope_over_one_degree_m']:.6f} m/rad vs GM from the "
        f"waterplane moments {gm['gm_from_waterplane_m']:.6f} m "
        f"(closed form {gm['gm_closed_form_m']:.6f} m): relative difference "
        f"{gm['relative_difference']:+.2e}, the residual being the wall-sided tan² term "
        "that a one-degree chord already carries.",
        f"- worst odd-symmetry defect |GZ(φ) + GZ(−φ)| = {results['odd_symmetry']:.2e} m.",
        f"- worst relative displacement residual over all sampled angles: "
        f"{results['worst_volume_residual']:.2e}.",
        "",
        "## 6. Claims",
        "",
        "| claim | value | threshold | holds |",
        "|---|---|---|---|",
    ]
    for entry in results["claims"]:
        lines.append(
            f"| {entry['claim']} | {entry['value']:.3e} | {entry['direction']} "
            f"{entry['threshold']:.0e} | {'✅' if entry['holds'] else '❌'} |"
        )
    lines += [
        "",
        "## 7. What this does not do",
        "",
        "- No free-surface correction: tanks are a loading property, not a shape one.",
        "- No damage stability, no wind/launching heeling moment, no class rule.",
        "- `G` is on the centreline; an off-centre `KG` needs a transverse `y_G` input.",
        "- The closed body *is* the hull: a mesh capped at the design waterline has zero "
        "freeboard and raises `HeelError` instead of returning a number. Computing a "
        "ship's GZ needs the body closed at the deck.",
        "- Angles are sampled, so `vanishing_angle_deg` is a linear interpolation between "
        "two samples (and `None` when the curve has not crossed by the last sample).",
        "",
        "## Figures",
        "",
        "- `gz_barge.png` — barge GZ curves against the wall-sided formula, with the "
        "deck-edge angles marked and the shallow barge's maximum and vanishing angle.",
        "- `gz_sphere.png` — sphere GZ against `sin φ (z_C − z_G)`, and the millimetre residual.",
        "- `sphere_ladder.png` — that residual against facet size, at slope 2.",
    ]
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("docs/gz_curve"))
    arguments = parser.parse_args(argv)
    output = arguments.output
    output.mkdir(parents=True, exist_ok=True)

    barge_angles = tuple(float(angle) for angle in range(0, 65, 5))
    shallow_angles = tuple(float(angle) for angle in range(0, 95, 5))
    results: dict[str, object] = {
        "version": __version__,
        "barge": measure_barge(BARGE, angles=barge_angles),
        "shallow": measure_barge(SHALLOW, angles=shallow_angles),
        "sphere": measure_sphere(),
        "sphere_ladder": measure_sphere_ladder(),
        "gm_consistency": measure_gm_consistency(BARGE),
        "odd_symmetry": measure_odd_symmetry(BARGE),
    }
    worst_residual = 0.0
    for key in ("barge", "shallow"):
        for row in results[key]["rows"]:
            target = results[key]["displacement_volume_m3"]
            worst_residual = max(
                worst_residual, abs(row["displacement_volume_m3"] - target) / target
            )
    target = results["sphere"]["displacement_volume_m3"]
    for row in results["sphere"]["rows"]:
        worst_residual = max(
            worst_residual, abs(row["displacement_volume_m3"] - target) / target
        )
    results["worst_volume_residual"] = worst_residual

    claims = evaluate(results)
    results["claims"] = claims

    for name, builder in (
        ("gz_barge.png", figure_barge),
        ("gz_sphere.png", figure_sphere),
        ("sphere_ladder.png", figure_ladder),
    ):
        try:
            builder(results, output / name)
        except ImportError:  # matplotlib is optional
            print(f"[gz_curve] matplotlib missing, skipping {name}", file=sys.stderr)

    (output / "results.json").write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")
    (output / "report.md").write_text(render_report(results))

    barge = results["barge"]
    shallow = results["shallow"]
    sphere = results["sphere"]
    print(
        f"[gz_curve] barge worst relative error vs wall-sided "
        f"{barge['worst_relative_error_vs_wall_sided']:.3e} (deck edge {barge['deck_edge_deg']:.2f} deg)"
    )
    print(
        f"[gz_curve] shallow barge max GZ {shallow['max_gz_m']:.4f} m at "
        f"{shallow['angle_of_max_gz_deg']:.0f} deg, vanishing {shallow['vanishing_angle_deg']:.2f} deg"
    )
    print(
        f"[gz_curve] sphere worst |error| {sphere['worst_absolute_error_m']:.3e} m "
        f"at {sphere['rings']} rings"
    )
    for row in results["sphere_ladder"]:
        factor = row.get("halving_factor")
        print(
            f"[gz_curve] sphere {row['rings']:>2} rings  error {row['absolute_error_m']:.3e} m"
            + ("" if factor is None else f"  ({factor:.2f}x per doubling)")
        )
    failed = [entry for entry in claims if not entry["holds"]]
    for entry in claims:
        status = "PASS" if entry["holds"] else "FAIL"
        print(f"[gz_curve] {status}  {entry['claim']}  ({entry['value']:.3e} vs {entry['threshold']:.0e})")
    if failed:
        print(f"[gz_curve] {len(failed)} claim(s) do not hold", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())