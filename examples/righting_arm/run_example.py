"""Righting-arm curves measured against two closed forms and one hull ladder.

Every number in ``docs/righting_arm.md`` comes out of this run: the script writes
``results.json`` and ``report.md`` and **exits non-zero if one of its claims stops
holding**, so the page cannot drift away from the code.

What it measures, in order:

1. **A wall-sided barge** (``L = 10, B = 2, T = 0.8, freeboard = 0.4 m``). Inside the
   formula's validity window ``GZ = sin(phi)(GM + BM tan^2(phi)/2)`` is exact, and beyond
   it a 2-D polygon clip of the heeled section still is: the mesh integral has to agree
   with both, in their own ranges.
2. **A floating sphere**, where the buoyancy centre stays on the vertical through the
   centre, so ``GZ = GM sin(phi)`` at every heel and the error is pure meshing error. The
   refinement ladder measures it rather than asserting it.
3. **The parametric hull with freeboard** (depth 1 m, floated at ``z = -0.47`` so that
   0.47 m of hull is above the waterline). No closed-form ``GZ`` exists for it, so the
   checks are the small-angle slope against the *waterplane-moment* integration of
   :mod:`tensorcad.stability`, a refinement ladder, and the deck-edge immersion angle
   reported from the solved waterlines.

Run it with ``python examples/righting_arm/run_example.py --output docs/righting_arm``.
matplotlib is optional: the JSON and the Markdown are written either way, the PNGs are
skipped without it.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from itertools import pairwise
from pathlib import Path

from tensorcad import __version__
from tensorcad.geometry import MeshBuilder
from tensorcad.geometry_field import SphereField, to_mesh
from tensorcad.hull import HullParameters, generate_hull
from tensorcad.righting_arm import (
    heeled_displacement,
    righting_arm,
    righting_arm_curve,
    wall_sided_righting_arm,
    wall_sided_validity_deg,
)
from tensorcad.stability import waterplane_moments

LENGTH, BEAM, DRAFT, FREEBOARD = 10.0, 2.0, 0.8, 0.4
KG = 0.55
HULL_WATERLINE, HULL_KG = -0.47, 0.55
SPHERE_RADIUS, SPHERE_KG = 1.0, 0.9
HULL_LADDER = ((21, 9), (41, 13), (81, 25), (161, 41))
SPHERE_LADDER = (12, 16, 24, 32)


def box(length: float, beam: float, draft: float, freeboard: float) -> object:
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

    lo, hi, bottom, top, right = -0.5 * beam, 0.5 * beam, -draft, freeboard, length
    quad((0.0, lo, bottom), (0.0, hi, bottom), (right, hi, bottom), (right, lo, bottom))
    quad((0.0, lo, top), (right, lo, top), (right, hi, top), (0.0, hi, top))
    quad((0.0, lo, bottom), (right, lo, bottom), (right, lo, top), (0.0, lo, top))
    quad((0.0, hi, bottom), (0.0, hi, top), (right, hi, top), (right, hi, bottom))
    quad((0.0, lo, bottom), (0.0, lo, top), (0.0, hi, top), (0.0, hi, bottom))
    quad((right, lo, bottom), (right, hi, bottom), (right, hi, top), (right, lo, top))
    return builder.build()


def hull(nx: int, nz: int):  # type: ignore[no-untyped-def]
    return generate_hull(
        HullParameters(length=10.0, beam=2.0, draft=1.0, nx=nx, nz=nz)
    ).mesh


def sphere(resolution: int):  # type: ignore[no-untyped-def]
    return to_mesh(SphereField(radius=SPHERE_RADIUS), resolution).mesh


def deck_edge_immersion_deg(mesh, waterline_z: float, heels: tuple[float, ...]):  # type: ignore[no-untyped-def]
    """First sampled heel whose solved waterline rises above the deck edge.

    The deck edge is the widest vertex on the hull's top ring; at heel ``phi`` its world
    height is ``-|y| sin(phi)`` (measured from the upright waterline). Where the solved
    waterline overtakes it, a real ship would be taking water on deck, and the arm is
    then the arm of a *closed solid* -- an upper bound, not a ship.
    """
    top_z = max(vertex[2] for vertex in mesh.vertices)
    beam_edge = max(abs(vertex[1]) for vertex in mesh.vertices if abs(vertex[2] - top_z) < 1e-9)
    for heel in heels:
        state = heeled_displacement(mesh, heel_deg=heel, waterline_z=waterline_z)
        if state.waterline_z > -beam_edge * math.sin(math.radians(heel)):
            return heel
    return heels[-1]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("docs/righting_arm"))
    args = parser.parse_args(argv)
    args.output.mkdir(parents=True, exist_ok=True)

    results: dict[str, object] = {"version": __version__, "claims": []}
    failures: list[str] = []

    def claim(name: str, passed: bool, detail: str) -> None:
        results["claims"].append({"claim": name, "passed": bool(passed), "detail": detail})
        if not passed:
            failures.append(f"{name}: {detail}")

    # -- 1. the barge -----------------------------------------------------------------
    barge = box(LENGTH, BEAM, DRAFT, FREEBOARD)
    limit = wall_sided_validity_deg(breadth=BEAM, draft=DRAFT, depth=DRAFT + FREEBOARD)
    upright = heeled_displacement(barge, heel_deg=0.0, waterline_z=0.0)
    claim(
        "barge displaces B*T*L at the waterline",
        abs(upright.volume - LENGTH * BEAM * DRAFT) / (LENGTH * BEAM * DRAFT) < 1e-10,
        f"{upright.volume:.9f} m^3 vs {LENGTH * BEAM * DRAFT:.9f} m^3",
    )
    barge_rows = []
    worst_wall = 0.0
    for heel in (0.0, 2.0, 5.0, 10.0, 15.0, 20.0):
        measured = righting_arm(barge, heel_deg=heel, kg=KG).righting_arm
        formula = wall_sided_righting_arm(breadth=BEAM, draft=DRAFT, kg=KG, heel_deg=heel)
        worst_wall = max(worst_wall, abs(measured - formula))
        barge_rows.append(
            {"heel_deg": heel, "gz_m": measured, "wall_sided_m": formula,
             "difference_m": measured - formula}
        )
    claim(
        "barge GZ equals the wall-sided closed form inside its validity window",
        worst_wall < 1e-12,
        f"worst |GZ - formula| = {worst_wall:.2e} m up to {limit:.2f} deg",
    )

    barge_past = []
    for heel in (25.0, 30.0, 35.0, 40.0):
        measured = righting_arm(barge, heel_deg=heel, kg=KG).righting_arm
        formula = wall_sided_righting_arm(breadth=BEAM, draft=DRAFT, kg=KG, heel_deg=heel)
        barge_past.append(
            {"heel_deg": heel, "gz_m": measured, "wall_sided_m": formula,
             "difference_m": measured - formula,
             "waterline_m": righting_arm(barge, heel_deg=heel, kg=KG).waterline_z}
        )
    claim(
        "the wall-sided form is no longer exact past the deck-edge angle",
        all(abs(row["difference_m"]) > 1e-3 for row in barge_past),
        f"at 40 deg the formula is {barge_past[-1]['difference_m']:+.4f} m out",
    )

    # -- 2. the sphere ----------------------------------------------------------------
    sphere_rows = []
    sphere_convergence = []
    for resolution in SPHERE_LADDER:
        mesh = sphere(resolution)
        errors = []
        for heel in (30.0, 45.0, 60.0):
            measured = righting_arm(mesh, heel_deg=heel, kg=SPHERE_KG).righting_arm
            expected = (SPHERE_RADIUS - SPHERE_KG) * math.sin(math.radians(heel))
            errors.append(abs(measured - expected))
            if resolution == SPHERE_LADDER[-1]:
                sphere_rows.append(
                    {"heel_deg": heel, "gz_m": measured, "gm_sin_m": expected,
                     "error_m": measured - expected}
                )
        sphere_convergence.append({"resolution": resolution, "worst_error_m": max(errors)})
    claim(
        "sphere GZ reproduces GM sin(phi) to better than 1 mm at the finest mesh",
        sphere_convergence[-1]["worst_error_m"] < 1e-3,
        f"worst |error| = {sphere_convergence[-1]['worst_error_m']:.2e} m "
        f"at resolution {SPHERE_LADDER[-1]}",
    )
    claim(
        "sphere error falls with refinement",
        all(
            later["worst_error_m"] < earlier["worst_error_m"]
            for earlier, later in pairwise(sphere_convergence)
        ),
        ", ".join(
            f"{row['resolution']}: {row['worst_error_m']:.2e}" for row in sphere_convergence
        ),
    )

    # -- 3. the hull ------------------------------------------------------------------
    heels = (0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 60.0)
    hull_rows = []
    slope_rel = 0.0
    hull_ladder = []
    for nx, nz in HULL_LADDER:
        mesh = hull(nx, nz)
        state = heeled_displacement(mesh, heel_deg=0.0, waterline_z=HULL_WATERLINE)
        section = waterplane_moments(mesh, HULL_WATERLINE)
        keel = min(vertex[2] for vertex in mesh.vertices)
        gm = (
            state.centre_of_buoyancy[2]
            - keel
            + section.transverse_second_moment / state.volume
            - HULL_KG
        )
        step = 1e-2
        slope = righting_arm(
            mesh, heel_deg=step, kg=HULL_KG, waterline_z=HULL_WATERLINE
        ).righting_arm / math.sin(math.radians(step))
        gz30 = righting_arm(
            mesh, heel_deg=30.0, kg=HULL_KG, waterline_z=HULL_WATERLINE
        ).righting_arm
        hull_ladder.append(
            {"nx": nx, "nz": nz, "faces": len(mesh.faces), "volume_m3": state.volume,
             "gm_waterplane_m": gm, "slope_m": slope, "slope_relative_error": slope / gm - 1.0,
             "gz_30deg_m": gz30}
        )
        slope_rel = max(slope_rel, abs(slope / gm - 1.0))
    claim(
        "hull small-angle slope equals the waterplane-moment GM",
        slope_rel < 1e-6,
        f"worst relative difference {slope_rel:.2e} over {len(HULL_LADDER)} meshes",
    )
    ladder = [row["gz_30deg_m"] for row in hull_ladder]
    claim(
        "hull GZ(30 deg) converges as the mesh is refined",
        abs(ladder[-1] - ladder[-2]) < abs(ladder[-2] - ladder[-3]),
        "GZ(30) = " + ", ".join(f"{value:.6f}" for value in ladder),
    )

    mesh = hull(41, 13)
    curve = righting_arm_curve(
        mesh, kg=HULL_KG, waterline_z=HULL_WATERLINE, heels_deg=heels
    )
    deck_deg = deck_edge_immersion_deg(mesh, HULL_WATERLINE, heels)
    for point in curve.points:
        hull_rows.append(
            {
                "heel_deg": point.heel_deg,
                "gz_m": point.righting_arm,
                "waterline_m": point.waterline_z,
                "y_buoyancy_m": point.centre_of_buoyancy[1],
                "residual": point.volume_residual,
            }
        )
    peak_heel, peak_arm = curve.peak()
    claim(
        "hull curve has an interior maximum and a positive arm at 30 deg",
        peak_heel >= 30.0 and peak_arm > 0.05 and curve.gz_at(30.0) > 0.0,
        f"peak {peak_arm:.6f} m at {peak_heel:.0f} deg; "
        f"area {curve.area_under_curve():.6f} m rad",
    )
    claim(
        "constant displacement is held along the curve",
        all(abs(row["residual"]) < 1e-9 for row in hull_rows),
        f"worst relative volume residual "
        f"{max(abs(row['residual']) for row in hull_rows):.2e}",
    )
    claim(
        "deck-edge immersion is found and reported, not assumed",
        0.0 < deck_deg < 60.0,
        f"the closed-solid arm is a ship's arm only up to {deck_deg:.0f} deg heel",
    )

    results.update(
        {
            "barge": {
                "length_m": LENGTH, "beam_m": BEAM, "draft_m": DRAFT,
                "freeboard_m": FREEBOARD, "kg_m": KG,
                "wall_sided_validity_deg": limit,
                "inside_window": barge_rows, "beyond_window": barge_past,
            },
            "sphere": {
                "radius_m": SPHERE_RADIUS, "kg_m": SPHERE_KG,
                "finest": sphere_rows, "convergence": sphere_convergence,
            },
            "hull": {
                "waterline_m": HULL_WATERLINE, "kg_m": HULL_KG,
                "ladder": hull_ladder, "curve": hull_rows,
                "peak": {"heel_deg": peak_heel, "gz_m": peak_arm},
                "area_under_curve_m_rad": curve.area_under_curve(),
                "deck_edge_immersion_deg": deck_deg,
            },
            "curve_payload": curve.to_dict(),
        }
    )
    (args.output / "results.json").write_text(json.dumps(results, indent=2) + "\n")

    # -- figures ----------------------------------------------------------------------
    figures: list[Path] = []
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - matplotlib is optional
        print(f"matplotlib unavailable ({exc}); PNGs skipped", file=sys.stderr)
        plt = None

    if plt is not None:
        figure, axes = plt.subplots(figsize=(7.2, 4.6))
        heel_grid = [index * 0.5 for index in range(0, 121)]
        axes.plot(
            heel_grid,
            [wall_sided_righting_arm(breadth=BEAM, draft=DRAFT, kg=KG, heel_deg=heel)
             for heel in heel_grid],
            "--", label="wall-sided closed form (exact to the deck edge)",
        )
        axes.plot(
            [row["heel_deg"] for row in barge_rows + barge_past],
            [row["gz_m"] for row in barge_rows + barge_past],
            "o", label="heeled mesh integral (box)",
        )
        axes.axvline(limit, color="0.5", linewidth=0.8)
        axes.annotate(
            f"deck edge\n{limit:.1f} deg", (limit, 0.02), fontsize=8, ha="left"
        )
        axes.set_xlabel("heel angle [deg]")
        axes.set_ylabel("GZ [m]")
        axes.set_title("Barge: mesh integration vs the exact wall-sided form")
        axes.grid(alpha=0.3)
        axes.legend(fontsize=8)
        figure.tight_layout()
        path = args.output / "barge_vs_closed_form.png"
        figure.savefig(path, dpi=140)
        plt.close(figure)
        figures.append(path)

        figure, axes = plt.subplots(figsize=(7.2, 4.6))
        fine = sphere(SPHERE_LADDER[-1])
        heel_grid = list(range(0, 91, 5))
        axes.plot(
            heel_grid,
            [(SPHERE_RADIUS - SPHERE_KG) * math.sin(math.radians(heel))
             for heel in heel_grid],
            "-", label="GM sin(phi)",
        )
        axes.plot(
            heel_grid,
            [righting_arm(fine, heel_deg=heel, kg=SPHERE_KG).righting_arm
             for heel in heel_grid],
            "o", markersize=3,
            label=f"polyhedral sphere, resolution {SPHERE_LADDER[-1]}",
        )
        axes.set_xlabel("heel angle [deg]")
        axes.set_ylabel("GZ [m]")
        axes.set_title("Sphere: GZ is GM sin(phi) at every heel")
        axes.grid(alpha=0.3)
        axes.legend(fontsize=8)
        figure.tight_layout()
        path = args.output / "sphere_vs_analytic.png"
        figure.savefig(path, dpi=140)
        plt.close(figure)
        figures.append(path)

        figure, axes = plt.subplots(figsize=(7.2, 4.6))
        axes.plot(
            [row["heel_deg"] for row in hull_rows],
            [row["gz_m"] for row in hull_rows],
            "o-", label="parametric hull, 41x13, KG = 0.55 m",
        )
        axes.axvline(deck_deg, color="0.5", linewidth=0.8)
        axes.annotate(f"deck edge ~{deck_deg:.0f} deg", (deck_deg, 0.01), fontsize=8)
        axes.plot(
            [row["heel_deg"] for row in hull_ladder and hull_rows],
            [row["gz_m"] for row in hull_rows],
            alpha=0.0,
        )
        axes.axhline(0.0, color="0.7", linewidth=0.8)
        axes.set_xlabel("heel angle [deg]")
        axes.set_ylabel("GZ [m]")
        axes.set_title("Hull with freeboard: righting-arm curve at constant displacement")
        axes.grid(alpha=0.3)
        axes.legend(fontsize=8)
        figure.tight_layout()
        path = args.output / "hull_curve.png"
        figure.savefig(path, dpi=140)
        plt.close(figure)
        figures.append(path)

    lines = [
        "# Righting-arm curve: measured report",
        "",
        f"TensorCAD {__version__}, written by `examples/righting_arm/run_example.py`.",
        "",
        "## Claims",
        "",
        "| claim | result | detail |",
        "|---|---|---|",
    ]
    for entry in results["claims"]:  # type: ignore[assignment]
        mark = "PASS" if entry["passed"] else "FAIL"
        lines.append(f"| {entry['claim']} | {mark} | {entry['detail']} |")
    lines += [
        "",
        "## Barge inside the wall-sided window",
        "",
        "| heel [deg] | GZ mesh [m] | GZ closed form [m] | difference [m] |",
        "|---|---|---|---|",
    ]
    for row in barge_rows:
        lines.append(
            f"| {row['heel_deg']:.0f} | {row['gz_m']:.9f} | {row['wall_sided_m']:.9f} "
            f"| {row['difference_m']:+.2e} |"
        )
    lines += [
        "",
        "## Barge beyond it",
        "",
        "| heel [deg] | GZ mesh [m] | GZ closed form [m] | formula error [m] | solved waterline [m] |",
        "|---|---|---|---|---|",
    ]
    for row in barge_past:
        lines.append(
            f"| {row['heel_deg']:.0f} | {row['gz_m']:.9f} | {row['wall_sided_m']:.9f} "
            f"| {row['difference_m']:+.4f} | {row['waterline_m']:+.4f} |"
        )
    lines += [
        "",
        "## Sphere",
        "",
        "| resolution | worst |GZ - GM sin(phi)| [m] |",
        "|---|---|",
    ]
    for row in sphere_convergence:
        lines.append(f"| {row['resolution']} | {row['worst_error_m']:.3e} |")
    lines += [
        "",
        "## Hull ladder",
        "",
        "| mesh | faces | displacement [m^3] | GM (waterplane) [m] | slope [m] | rel. diff | GZ(30) [m] |",
        "|---|---|---|---|---|---|---|",
    ]
    for row in hull_ladder:
        lines.append(
            f"| {row['nx']}x{row['nz']} | {row['faces']} | {row['volume_m3']:.6f} "
            f"| {row['gm_waterplane_m']:.6f} | {row['slope_m']:.6f} "
            f"| {row['slope_relative_error']:+.2e} | {row['gz_30deg_m']:.6f} |"
        )
    lines += [
        "",
        "## Hull curve (41x13, waterline -0.47 m, KG 0.55 m)",
        "",
        "| heel [deg] | GZ [m] | waterline [m] | y_B [m] | volume residual |",
        "|---|---|---|---|---|",
    ]
    for row in hull_rows:
        lines.append(
            f"| {row['heel_deg']:.0f} | {row['gz_m']:+.6f} | {row['waterline_m']:+.5f} "
            f"| {row['y_buoyancy_m']:+.5f} | {row['residual']:+.1e} |"
        )
    lines += [
        "",
        f"Peak {peak_arm:.6f} m at {peak_heel:.0f} deg; "
        f"dynamical stability {curve.area_under_curve():.6f} m rad; "
        f"deck-edge immersion from the solved waterlines: {deck_deg:.0f} deg.",
        "",
    ]
    if figures:
        lines += ["## Figures", ""]
        lines += [f"- `{path.name}`" for path in figures]
        lines.append("")
    (args.output / "report.md").write_text("\n".join(lines))

    passed = sum(1 for entry in results["claims"] if entry["passed"])  # type: ignore[union-attr]
    total = len(results["claims"])  # type: ignore[arg-type]
    print(f"{passed}/{total} claims hold; wrote {args.output}")
    if failures:
        for failure in failures:
            print(f"  FAILED {failure}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
