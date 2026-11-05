#!/usr/bin/env python3
"""Verify the large-angle righting arm against references it cannot fake.

``examples/stability`` checked ``GM`` -- the slope of the righting-arm curve at zero
heel. This example checks the curve, which is the thing a designer actually reads,
and it does so against three independent references:

1. **A box section.** For a prismatic hull the heeled displacement integral separates
   into a 2D section clip, so ``analytic_box_righting_arm`` recomputes every point by
   a route that shares no code with the 3D tetrahedron integration: agreement is
   reported in metres and must be at round-off, at every heel angle up to 90 deg,
   including past deck-edge immersion.
2. **The wall-sided closed form** ``GZ = sin(phi)(GM + BM tan^2(phi)/2)``. It is exact
   only while both sides of the section are wet, i.e. below ``bilge_emergence_deg``;
   the example reports the agreement inside that window *and* the departure just
   outside it, so the validity limit is measured rather than asserted.
3. **The upright metacentric height**, ``GM = KB + BM - KG``, assembled from
   ``waterplane_moments`` (Green's theorem on the waterline) and the volume centroid
   measured here. ``dGZ/dphi`` at zero heel has to reproduce it -- a check that ties
   the heeled integral to a different integration entirely.

Then the same machinery is run on the parametric Wigley hull at three grid
resolutions, where no closed form for ``GZ`` exists: what is reported there is the
consistency with ``GM``, the odd symmetry in heel, the volume held to round-off, and
the mesh-convergence of the arm itself.

Writes ``results.json``, ``report.md`` and ``plot.png`` into ``--output`` and exits
non-zero if any claim stops holding.

Usage::

    python examples/righting/run_example.py [--output docs/righting]
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:  # allow running from a checkout without installing
    sys.path.insert(0, str(REPO_ROOT))

from tensorcad import mesh as meshlib
from tensorcad.righting import (
    analytic_box_righting_arm,
    bilge_emergence_deg,
    box_mesh,
    deck_edge_immersion_deg,
    righting_arm,
    righting_arm_curve,
    submerged_properties,
    upright_displacement,
    wall_sided_righting_arm,
)
from tensorcad.stability import waterplane_moments
from tensorcad.wigley import WigleyParameters, generate_wigley

# ---------------------------------------------------------------- reference cases

BOX_LENGTH = 10.0
BOX_BEAM = 2.0
BOX_DEPTH = 1.0
BOX_DRAFT = 0.5
BOX_KG = 0.55
BOX_BM = BOX_BEAM**2 / (12.0 * BOX_DRAFT)
BOX_KB = 0.5 * BOX_DRAFT
BOX_GM = BOX_KB + BOX_BM - BOX_KG

HULL_WATERLINE = -0.1  # the Wigley cap is at z = 0, so this leaves 0.1 m of freeboard
HULL_KG = 0.30
HULL_GRIDS = ((41, 13), (81, 25), (161, 41))
HULL_ANGLES = (0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0)


def box_study() -> dict[str, object]:
    mesh = box_mesh(BOX_LENGTH, BOX_BEAM, BOX_DEPTH)
    bilge = bilge_emergence_deg(beam=BOX_BEAM, draft=BOX_DRAFT)
    deck = deck_edge_immersion_deg(beam=BOX_BEAM, depth=BOX_DEPTH)
    angles = (0.0, 2.0, 5.0, 10.0, 15.0, 20.0, 25.0, 26.0, 30.0, 40.0, 45.0, 60.0, 75.0, 90.0)
    rows = []
    for heel_deg in angles:
        point = righting_arm(mesh, heel_deg, kg=BOX_KG, waterline_z=BOX_DRAFT)
        reference = analytic_box_righting_arm(
            BOX_LENGTH, BOX_BEAM, BOX_DEPTH, BOX_DRAFT, heel_deg, kg=BOX_KG
        )
        form = wall_sided_righting_arm(math.radians(heel_deg), gm=BOX_GM, bm=BOX_BM)
        row_residual = point.volume_residual
        rows.append(
            {
                "heel_deg": heel_deg,
                "gz_m": point.righting_arm,
                "gz_section_reference_m": reference,
                "difference_m": point.righting_arm - reference,
                "wall_sided_m": form,
                "wall_sided_difference_m": point.righting_arm - form,
                "waterplane_area_m2": point.waterplane_area,
                "volume_residual_m3": row_residual,
                "iterations": point.iterations,
                "centre_of_buoyancy_m": list(point.centre_of_buoyancy),
            }
        )
    curve = righting_arm_curve(mesh, angles, kg=BOX_KG, waterline_z=BOX_DRAFT)
    step = 1e-4
    slope = (
        righting_arm(mesh, step, kg=BOX_KG, waterline_z=BOX_DRAFT).righting_arm
        - righting_arm(mesh, -step, kg=BOX_KG, waterline_z=BOX_DRAFT).righting_arm
    ) / (2.0 * math.radians(step))
    inside = [row for row in rows if 0.0 < row["heel_deg"] < bilge]
    outside = [row for row in rows if row["heel_deg"] > bilge + 1.0]
    return {
        "parameters": {
            "length_m": BOX_LENGTH,
            "beam_m": BOX_BEAM,
            "depth_m": BOX_DEPTH,
            "draft_m": BOX_DRAFT,
            "kg_m": BOX_KG,
            "kb_m": BOX_KB,
            "bm_m": BOX_BM,
            "gm_m": BOX_GM,
            "displacement_volume_m3": BOX_LENGTH * BOX_BEAM * BOX_DRAFT,
        },
        "bilge_emergence_deg": bilge,
        "deck_edge_immersion_deg": deck,
        "rows": rows,
        "max_abs_difference_m": max(abs(row["difference_m"]) for row in rows),
        "max_abs_wall_sided_difference_inside_m": max(
            abs(row["wall_sided_difference_m"]) for row in inside
        ),
        "min_abs_wall_sided_difference_outside_m": min(
            abs(row["wall_sided_difference_m"]) for row in outside
        ),
        "max_relative_volume_residual": max(
            abs(row["volume_residual_m3"]) / (BOX_LENGTH * BOX_BEAM * BOX_DRAFT) for row in rows
        ),
        "max_iterations": max(row["iterations"] for row in rows),
        "slope_at_zero_m_per_rad": slope,
        "slope_relative_error": (slope - BOX_GM) / BOX_GM,
        "max_arm_m": curve.max_arm,
        "angle_of_max_arm_deg": curve.angle_of_max_arm_deg,
        "vanishing_stability_deg": curve.zero_crossing_deg,
        "area_to_30_deg_m_rad": curve.area(30.0),
        "area_to_40_deg_m_rad": curve.area(40.0),
        "watertight": meshlib.is_watertight(mesh),
        "enclosed_volume_m3": meshlib.total_signed_volume(mesh),
    }


def hull_study() -> dict[str, object]:
    case: dict[str, object] = {"waterline_z_m": HULL_WATERLINE, "kg_m": HULL_KG, "grids": []}
    for nx, nz in HULL_GRIDS:
        mesh = generate_wigley(WigleyParameters(nx=nx, nz=nz)).mesh
        displacement = upright_displacement(mesh, HULL_WATERLINE)
        moments = waterplane_moments(mesh, HULL_WATERLINE)
        keel_z = min(vertex[2] for vertex in mesh.vertices)
        kb = submerged_properties(mesh, (0.0, 0.0, 1.0), HULL_WATERLINE).centre_of_buoyancy[2] - keel_z
        bm = moments.transverse_second_moment / displacement
        gm = kb + bm - HULL_KG
        curve = righting_arm_curve(mesh, HULL_ANGLES, kg=HULL_KG, waterline_z=HULL_WATERLINE)
        step = 1.0
        slope = (
            righting_arm(mesh, step, kg=HULL_KG, waterline_z=HULL_WATERLINE).righting_arm
            - righting_arm(mesh, -step, kg=HULL_KG, waterline_z=HULL_WATERLINE).righting_arm
        ) / (2.0 * math.radians(step))
        symmetry = max(
            abs(
                righting_arm(mesh, angle, kg=HULL_KG, waterline_z=HULL_WATERLINE).righting_arm
                + righting_arm(mesh, -angle, kg=HULL_KG, waterline_z=HULL_WATERLINE).righting_arm
            )
            for angle in (5.0, 20.0, 40.0)
        )
        case["grids"].append(
            {
                "nx": nx,
                "nz": nz,
                "faces": len(mesh.faces),
                "displacement_volume_m3": displacement,
                "kb_m": kb,
                "bm_m": bm,
                "gm_m": gm,
                "slope_at_zero_m_per_rad": slope,
                "slope_relative_error": (slope - gm) / gm,
                "odd_symmetry_residual_m": symmetry,
                "max_relative_volume_residual": max(
                    point.relative_volume_residual for point in curve.points
                ),
                "max_iterations": max(point.iterations for point in curve.points),
                "angles_deg": list(curve.angles_deg),
                "arms_m": list(curve.arms),
                "max_arm_m": curve.max_arm,
                "angle_of_max_arm_deg": curve.angle_of_max_arm_deg,
                "area_to_30_deg_m_rad": curve.area(30.0),
                "area_to_40_deg_m_rad": curve.area(40.0),
                "vertex_offsets_volume_m3": meshlib.total_signed_volume(mesh),
            }
        )
    arms = [grid["arms_m"] for grid in case["grids"]]
    differences = [
        [abs(b - a) for a, b in zip(coarse, fine, strict=True)] for coarse, fine in itertools.pairwise(arms)
    ]
    case["successive_differences_m"] = differences
    case["index_20_deg"] = HULL_ANGLES.index(20.0)
    return case


def claims(box: dict[str, object], hull: dict[str, object]) -> list[tuple[str, bool, str]]:
    grids = hull["grids"]
    index = hull["index_20_deg"]
    coarse, medium, fine = (grid["arms_m"][index] for grid in grids)
    return [
        (
            "box GZ equals the 2D section reference at every heel angle",
            box["max_abs_difference_m"] < 1e-12,
            f"max |diff| = {box['max_abs_difference_m']:.2e} m over {len(box['rows'])} angles",
        ),
        (
            "wall-sided closed form is reproduced below bilge emergence",
            box["max_abs_wall_sided_difference_inside_m"] < 1e-12,
            f"max |diff| = {box['max_abs_wall_sided_difference_inside_m']:.2e} m "
            f"for 0 < phi < {box['bilge_emergence_deg']:.3f} deg",
        ),
        (
            "wall-sided closed form departs above bilge emergence",
            box["min_abs_wall_sided_difference_outside_m"] > 1e-3,
            f"min |diff| = {box['min_abs_wall_sided_difference_outside_m']:.2e} m above the limit",
        ),
        (
            "box GZ slope at zero heel is the analytic GM",
            abs(box["slope_relative_error"]) < 1e-6,
            f"dGZ/dphi = {box['slope_at_zero_m_per_rad']:.12f} vs GM = {box['parameters']['gm_m']:.12f}",
        ),
        (
            "box displacement is held to round-off at every heel angle",
            box["max_relative_volume_residual"] < 1e-12,
            f"max relative residual = {box['max_relative_volume_residual']:.2e}",
        ),
        (
            "hull GZ slope at zero heel matches GM from the waterplane moments",
            all(abs(grid["slope_relative_error"]) < 2e-3 for grid in grids),
            ", ".join(f"{grid['nx']}x{grid['nz']}: {grid['slope_relative_error']:+.2e}" for grid in grids),
        ),
        (
            "hull GZ is odd in heel",
            all(grid["odd_symmetry_residual_m"] < 1e-12 for grid in grids),
            f"max |GZ(phi) + GZ(-phi)| = {max(grid['odd_symmetry_residual_m'] for grid in grids):.2e} m",
        ),
        (
            "hull displacement is held to round-off at every heel angle",
            all(grid["max_relative_volume_residual"] < 1e-12 for grid in grids),
            f"max relative residual = {max(grid['max_relative_volume_residual'] for grid in grids):.2e}",
        ),
        (
            "hull GZ at 20 deg converges monotonically with grid refinement",
            coarse < medium < fine and abs(fine - medium) < abs(medium - coarse),
            f"{coarse:.6f} -> {medium:.6f} -> {fine:.6f} m",
        ),
        (
            "the offset solve stays cheap",
            max(grid["max_iterations"] for grid in grids) <= 10,
            f"max Newton iterations = {max(grid['max_iterations'] for grid in grids)}",
        ),
    ]


# ------------------------------------------------------------------------ plotting


def plot(box: dict[str, object], hull: dict[str, object], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(2, 2, figsize=(10.4, 7.6), dpi=90)
    colours = ("#8c8c8c", "#1f77b4", "#d62728")

    ax = axes[0][0]
    rows = box["rows"]
    angles = [row["heel_deg"] for row in rows]
    ax.plot(angles, [row["gz_m"] for row in rows], "o", ms=4.5, color="#1f77b4",
            label="mesh, 3D tetrahedron integration")
    ax.plot(angles, [row["gz_section_reference_m"] for row in rows], "-", lw=1.2,
            color="#d62728", label="2D section reference (independent)")
    ax.plot(angles, [row["wall_sided_m"] for row in rows], "--", lw=1.0, color="#2ca02c",
            label=r"wall-sided $\sin\phi\,(GM+\frac{1}{2}BM\tan^2\phi)$")
    ax.axvline(box["bilge_emergence_deg"], color="0.4", ls=":", lw=1.0)
    ax.axvline(box["deck_edge_immersion_deg"], color="0.4", ls="-.", lw=1.0)
    ax.annotate("bilge emerges", (box["bilge_emergence_deg"], 0.02), rotation=90,
                fontsize=7, ha="right", va="bottom", color="0.35")
    ax.annotate("deck edge immersed", (box["deck_edge_immersion_deg"], 0.02), rotation=90,
                fontsize=7, ha="right", va="bottom", color="0.35")
    ax.set_xlabel("heel angle [deg]")
    ax.set_ylabel("GZ [m]")
    ax.set_title("Box reference case: GZ(phi) at constant displacement")
    ax.legend(loc="upper right", fontsize=7)
    ax.grid(alpha=0.3)

    ax = axes[0][1]
    for colour, grid in zip(colours, hull["grids"], strict=True):
        ax.plot(grid["angles_deg"], grid["arms_m"], "o-", ms=3.5, lw=1.1, color=colour,
                label=f"{grid['nx']}x{grid['nz']} ({grid['faces']} faces)")
    ax.set_xlabel("heel angle [deg]")
    ax.set_ylabel("GZ [m]")
    ax.set_title(f"Wigley hull at z = {hull['waterline_z_m']} m, KG = {hull['kg_m']} m")
    ax.legend(loc="upper left", fontsize=7)
    ax.grid(alpha=0.3)

    ax = axes[1][0]
    differences = hull["successive_differences_m"]
    labels = ("81x25 - 41x13", "161x41 - 81x25")
    for colour, series, label in zip(("#1f77b4", "#d62728"), differences, labels, strict=True):
        series = [max(value, 1e-12) for value in series[1:]]
        ax.semilogy(hull["grids"][0]["angles_deg"][1:], series, "o-", ms=3.5, lw=1.1,
                    color=colour, label=label)
    ax.set_xlabel("heel angle [deg]")
    ax.set_ylabel("|change in GZ| between grids [m]")
    ax.set_title("Mesh convergence of the hull righting arm")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3, which="both")

    # Submerged sections of the box: what the 3D clip is actually integrating.
    ax = axes[1][1]
    half_beam = 0.5 * BOX_BEAM
    for heel_deg, colour in zip((0.0, 20.0, 45.0, 75.0), ("#1f77b4", "#2ca02c", "#ff7f0e", "#d62728"), strict=True):
        point = righting_arm(box_mesh(BOX_LENGTH, BOX_BEAM, BOX_DEPTH), heel_deg,
                             kg=BOX_KG, waterline_z=BOX_DRAFT)
        phi = math.radians(heel_deg)
        cosine, sine, offset = math.cos(phi), math.sin(phi), point.plane_offset
        rectangle = ((-half_beam, 0.0), (half_beam, 0.0), (half_beam, BOX_DEPTH), (-half_beam, BOX_DEPTH))
        heights = [-y * sine + z * cosine - offset for y, z in rectangle]
        wet = []
        for index in range(4):
            point_2d, next_2d = rectangle[index], rectangle[(index + 1) % 4]
            height, next_height = heights[index], heights[(index + 1) % 4]
            if height <= 0.0:
                wet.append(point_2d)
            if (height <= 0.0) != (next_height <= 0.0):
                fraction = height / (height - next_height)
                wet.append((point_2d[0] + fraction * (next_2d[0] - point_2d[0]),
                            point_2d[1] + fraction * (next_2d[1] - point_2d[1])))
        if len(wet) >= 3:
            ax.fill([p[0] for p in wet], [p[1] for p in wet], alpha=0.18, color=colour)
            ax.plot([p[0] for p in wet] + [wet[0][0]], [p[1] for p in wet] + [wet[0][1]],
                    lw=1.0, color=colour, label=f"{heel_deg:.0f} deg")
        ax.plot([point.centre_of_buoyancy[1]], [point.centre_of_buoyancy[2]], "x",
                ms=6, color=colour)
    ax.plot([], [], "x", color="k", ms=6, label="centre of buoyancy")
    ax.set_xlim(-half_beam - 0.15, half_beam + 0.15)
    ax.set_ylim(-0.05, BOX_DEPTH + 0.05)
    ax.set_xlabel("y [m]")
    ax.set_ylabel("z [m]")
    ax.set_title("Box: submerged section and B at four heel angles")
    ax.legend(fontsize=7, loc="upper left")
    ax.grid(alpha=0.3)
    ax.set_aspect("equal")

    figure.tight_layout()
    figure.savefig(path, pil_kwargs={"optimize": True})
    plt.close(figure)


# ------------------------------------------------------------------------- output


def markdown(box: dict[str, object], hull: dict[str, object], checks: list[tuple[str, bool, str]]) -> str:
    parameters = box["parameters"]
    lines = [
        "# Large-angle righting arms: measured evidence",
        "",
        "Generated by `examples/righting/run_example.py`; every number below is injected",
        "from the run, not transcribed.",
        "",
        "## Box reference case",
        "",
        f"L = {parameters['length_m']} m, B = {parameters['beam_m']} m, D = {parameters['depth_m']} m, "
        f"T = {parameters['draft_m']} m, KG = {parameters['kg_m']} m; KB = {parameters['kb_m']:.6f} m, "
        f"BM = {parameters['bm_m']:.6f} m, GM = {parameters['gm_m']:.6f} m.",
        "",
        "| heel [deg] | GZ mesh [m] | GZ 2D ref [m] | diff [m] | wall-sided [m] | diff [m] | "
        "waterplane [m^2] | Newton it. |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in box["rows"]:
        lines.append(
            f"| {row['heel_deg']:.1f} | {row['gz_m']:+.9f} | {row['gz_section_reference_m']:+.9f} | "
            f"{row['difference_m']:+.2e} | {row['wall_sided_m']:+.9f} | "
            f"{row['wall_sided_difference_m']:+.2e} | {row['waterplane_area_m2']:.6f} | "
            f"{row['iterations']} |"
        )
    lines += [
        "",
        f"Bilge emergence at {box['bilge_emergence_deg']:.4f} deg, deck edge immersed at "
        f"{box['deck_edge_immersion_deg']:.4f} deg.",
        f"vanishing-stability angle {box['vanishing_stability_deg']}; "
        f"GZ max {box['max_arm_m']:.6f} m at {box['angle_of_max_arm_deg']:.1f} deg; "
        f"area to 30 deg {box['area_to_30_deg_m_rad']:.6f} m rad, to 40 deg "
        f"{box['area_to_40_deg_m_rad']:.6f} m rad.",
        "",
        "## Parametric Wigley hull",
        "",
        "| grid | faces | displacement [m^3] | KB [m] | BM [m] | GM [m] | dGZ/dphi(0) | rel. | "
        "GZ(20 deg) [m] | max GZ [m] |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    index = hull["index_20_deg"]
    for grid in hull["grids"]:
        lines.append(
            f"| {grid['nx']}x{grid['nz']} | {grid['faces']} | {grid['displacement_volume_m3']:.6f} | "
            f"{grid['kb_m']:.6f} | {grid['bm_m']:.6f} | {grid['gm_m']:.6f} | "
            f"{grid['slope_at_zero_m_per_rad']:.6f} | {grid['slope_relative_error']:+.2e} | "
            f"{grid['arms_m'][index]:.6f} | {grid['max_arm_m']:.6f} |"
        )
    lines += [
        "",
        "| heel [deg] | " + " | ".join(f"GZ {g['nx']}x{g['nz']} [m]" for g in hull["grids"]) + " |",
        "|---:|" + "---:|" * len(hull["grids"]),
    ]
    for position, angle in enumerate(HULL_ANGLES):
        lines.append(
            f"| {angle:.0f} | "
            + " | ".join(f"{g['arms_m'][position]:+.6f}" for g in hull["grids"])
            + " |"
        )
    lines += [
        "",
        "## Claims",
        "",
        "| claim | result | evidence |",
        "|---|---|---|",
    ]
    for claim, ok, evidence in checks:
        lines.append(f"| {claim} | {'PASS' if ok else 'FAIL'} | {evidence} |")
    lines.append("")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("results"))
    args = parser.parse_args(argv)
    args.output.mkdir(parents=True, exist_ok=True)

    box = box_study()
    hull = hull_study()
    checks = claims(box, hull)

    plot(box, hull, args.output / "plot.png")
    (args.output / "results.json").write_text(
        json.dumps({"box": box, "hull": hull,
                    "checks": [{"claim": c, "passed": ok, "evidence": e} for c, ok, e in checks]},
                   indent=2, sort_keys=True)
        + "\n"
    )
    (args.output / "report.md").write_text(markdown(box, hull, checks))

    for claim, ok, evidence in checks:
        print(f"[{'PASS' if ok else 'FAIL'}] {claim}: {evidence}")
    print(f"\nwrote {args.output}/results.json, report.md, plot.png")
    return 0 if all(ok for _, ok, _ in checks) else 1


if __name__ == "__main__":
    raise SystemExit(main())
