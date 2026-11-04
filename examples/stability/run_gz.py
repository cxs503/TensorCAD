"""Righting arms against two closed forms, and a loading ladder that loses its range.

``examples/stability/run_example.py`` measures ``GM`` -- a waterplane moment and the
slope of ``GZ`` at zero heel. This script measures ``GZ`` itself, which is a different
integral: the hull is heeled, sunk until it displaces the same volume it displaces
upright, and the transverse offset between the centroid of the displaced volume and the
centre of gravity is read off. Three things are measured, in increasing order of how much
they can go wrong:

1. **A wall-sided box** against ``GZ = sin(phi) (GM + BM tan^2(phi) / 2)``. Up to the
   angle at which a deck corner enters the water this is exact, so the agreement has to
   be round-off -- and the script also measures where the closed form *stops* being
   valid, because agreeing past that point would mean the solver was not clipping.
2. **A regular prism** against the circular section it inscribes, whose transverse
   metacentre is the centre of the circle (``KM = R``), i.e. ``GZ = (R - KG) sin(phi)``
   at every heel angle. The polygon converges to it at second order, which is measured
   rather than assumed.
3. **A barge** where the geometry is the only input: the arm peaks, then the deck goes
   under and the range of stability ends. A ``KG`` ladder shows a curve disappearing
   into a negative one, which is what ``GM`` alone cannot tell a designer.

The script exits non-zero if a claim stops holding, so the numbers in
``docs/stability.md`` cannot drift away from the code.

Run it with ``python examples/stability/run_gz.py --output docs/stability``.
matplotlib is optional: the JSON and the markdown are written either way and the PNGs
are skipped.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any

from tensorcad import __version__
from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.stability import gz_curve
import itertools

SCHEMA = "tensorcad.examples.stability_gz/1"


# ---------------------------------------------------------------------------
# reference bodies
# ---------------------------------------------------------------------------
def box(length: float, beam: float, draft: float, freeboard: float) -> Mesh:
    """Closed box with the design waterline at ``z = 0``: draft below it, freeboard above."""
    builder = MeshBuilder()
    corner = {}
    for x in (0.0, length):
        for y in (-0.5 * beam, 0.5 * beam):
            for z in (-draft, freeboard):
                corner[(x, y, z)] = builder.add_vertex((x, y, z))

    def quad(p0: tuple[float, float, float], p1: tuple[float, float, float],
             p2: tuple[float, float, float], p3: tuple[float, float, float]) -> None:
        a, b, c, d = (corner[point] for point in (p0, p1, p2, p3))
        builder.add_triangle(a, b, c)
        builder.add_triangle(a, c, d)

    left, right = 0.0, length
    lo, hi = -0.5 * beam, 0.5 * beam
    quad((left, lo, -draft), (left, hi, -draft), (right, hi, -draft), (right, lo, -draft))
    quad((left, lo, freeboard), (right, lo, freeboard), (right, hi, freeboard), (left, hi, freeboard))
    quad((left, lo, -draft), (right, lo, -draft), (right, lo, freeboard), (left, lo, freeboard))
    quad((right, hi, -draft), (left, hi, -draft), (left, hi, freeboard), (right, hi, freeboard))
    quad((right, lo, -draft), (right, hi, -draft), (right, hi, freeboard), (right, lo, freeboard))
    quad((left, hi, -draft), (left, lo, -draft), (left, lo, freeboard), (left, hi, freeboard))
    return builder.build()


def prism(sides: int, radius: float, length: float) -> Mesh:
    """Regular ``sides``-gon prism about ``x``, offset so no vertex lies on the waterline."""
    builder = MeshBuilder()
    ring = [
        (radius * math.cos(2.0 * math.pi * (index + 0.5) / sides),
         radius * math.sin(2.0 * math.pi * (index + 0.5) / sides))
        for index in range(sides)
    ]
    near = [builder.add_vertex((0.0, y, z)) for y, z in ring]
    far = [builder.add_vertex((length, y, z)) for y, z in ring]
    for index in range(sides):
        following = (index + 1) % sides
        builder.add_triangle(near[index], near[following], far[following])
        builder.add_triangle(near[index], far[following], far[index])
    for index in range(1, sides - 1):
        builder.add_triangle(near[0], near[index + 1], near[index])
        builder.add_triangle(far[0], far[index], far[index + 1])
    return builder.build()


def wall_sided_arm(heel_degrees: float, *, gm: float, bm: float) -> float:
    phi = math.radians(heel_degrees)
    return math.sin(phi) * (gm + 0.5 * bm * math.tan(phi) ** 2)


# ---------------------------------------------------------------------------
# the three measurements
# ---------------------------------------------------------------------------
TALL = {"length": 6.0, "beam": 4.0, "draft": 4.0, "freeboard": 4.0}
TALL_KG = 1.0
BARGE = {"length": 4.0, "beam": 2.0, "draft": 1.0, "freeboard": 0.3}
CIRCULAR = {"radius": 1.0, "length": 6.0, "kg": 0.6}


def measure_wall_sided() -> dict[str, Any]:
    bm = TALL["beam"] ** 2 / (12.0 * TALL["draft"])
    gm = 0.5 * TALL["draft"] + bm - TALL_KG
    limit = math.degrees(math.atan(2.0 * TALL["freeboard"] / TALL["beam"]))
    angles = tuple(range(0, 61, 5))
    curve = gz_curve(box(**TALL), kg=TALL_KG, heel_angles_deg=angles)
    rows: list[dict[str, Any]] = []
    worst = 0.0
    for arm in curve.arms:
        reference = wall_sided_arm(arm.heel_angle_deg, gm=gm, bm=bm)
        difference = arm.righting_arm - reference
        worst = max(worst, abs(difference))
        rows.append({
            "heel_angle_deg": arm.heel_angle_deg,
            "righting_arm_m": arm.righting_arm,
            "closed_form_m": reference,
            "difference_m": difference,
            "sinkage_m": arm.sinkage,
        })
    # the closed form is invalid past the deck-corner immersion angle: measure the break
    outside = {
        angle: gz_curve(box(**TALL), kg=TALL_KG, heel_angles_deg=(angle,)).arms[0]
        for angle in (65.0, 70.0)
    }
    break_off = {
        str(angle): {
            "righting_arm_m": arm.righting_arm,
            "closed_form_m": wall_sided_arm(angle, gm=gm, bm=bm),
            "difference_m": arm.righting_arm - wall_sided_arm(angle, gm=gm, bm=bm),
        }
        for angle, arm in outside.items()
    }
    return {
        "geometry": dict(TALL, kg=TALL_KG),
        "bm_m": bm,
        "gm_m": gm,
        "deck_corner_immersion_deg": limit,
        "rows": rows,
        "max_abs_difference_m": worst,
        "max_volume_residual": curve.max_volume_residual,
        "beyond_the_limit": break_off,
        "claim_difference_below_1e_12": worst < 1e-12,
    }


def measure_circular() -> dict[str, Any]:
    radius, kg, expected = CIRCULAR["radius"], CIRCULAR["kg"], CIRCULAR["radius"] - CIRCULAR["kg"]
    angles = (10, 20, 30, 45, 60, 75)
    rows: list[dict[str, Any]] = []
    for sides in (24, 48, 96, 192):
        curve = gz_curve(prism(sides, radius, CIRCULAR["length"]), kg=kg, heel_angles_deg=angles)
        worst = 0.0
        arms: list[dict[str, Any]] = []
        for arm in curve.arms:
            sine = math.sin(math.radians(arm.heel_angle_deg))
            arms.append({"heel_angle_deg": arm.heel_angle_deg, "righting_arm_m": arm.righting_arm,
                         "gz_over_sine_m": arm.righting_arm / sine})
            worst = max(worst, abs(arm.righting_arm / sine - expected) / expected)
        rows.append({
            "sides": sides,
            "circumradius_m": radius,
            "relative_error": worst,
            "max_volume_residual": curve.max_volume_residual,
            "arms": arms,
        })
    ratios = [coarse["relative_error"] / fine["relative_error"] for coarse, fine in itertools.pairwise(rows)]
    return {
        "radius_m": radius,
        "kg_m": kg,
        "km_m": radius,
        "expected_gz_over_sine_m": expected,
        "heel_angles_deg": list(angles),
        "rows": rows,
        "error_ratios_per_doubling": ratios,
        "finest_relative_error": rows[-1]["relative_error"],
        "claim_second_order": all(ratio > 2.5 for ratio in ratios),
        "claim_finest_below_5e_4": rows[-1]["relative_error"] < 5e-4,
    }


def measure_barge() -> dict[str, Any]:
    angles = tuple(range(0, 141, 5))
    ladder: list[dict[str, Any]] = []
    for kg in (0.6, 0.65, 0.7):
        curve = gz_curve(box(**BARGE), kg=kg, heel_angles_deg=angles)
        ladder.append({
            "kg_m": kg,
            "km_m": 0.5 * BARGE["draft"] + BARGE["beam"] ** 2 / (12.0 * BARGE["draft"]),
            "gm_m": 0.5 * BARGE["draft"] + BARGE["beam"] ** 2 / (12.0 * BARGE["draft"]) - kg,
            "displacement_volume_m3": curve.displacement_volume,
            "hea_angles_deg": list(curve.heel_angles_deg),
            "righting_arms_m": list(curve.righting_arms),
            "max_arm_m": curve.max_arm.righting_arm,
            "angle_of_max_arm_deg": curve.angle_of_max_arm_deg,
            "vanishing_angle_deg": curve.vanishing_angle_deg,
            "max_volume_residual": curve.max_volume_residual,
        })
    vanishings = [row["vanishing_angle_deg"] for row in ladder]
    monotone = all(
        earlier is not None and later is not None and earlier > later
        for earlier, later in itertools.pairwise(vanishings)
    )
    return {
        "geometry": BARGE,
        "ladder": ladder,
        "claim_vanishing_angle_falls_with_kg": monotone,
    }


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------
def draw_figures(measurements: dict[str, Any], output: Path) -> list[str]:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:  # pragma: no cover - matplotlib is optional
        return []

    written: list[str] = []

    wall = measurements["wall_sided"]
    circular = measurements["circular"]
    barge = measurements["barge"]

    figure, (left, right) = plt.subplots(1, 2, figsize=(11.0, 4.2))
    heels = [row["heel_angle_deg"] for row in wall["rows"]]
    left.plot(heels, [row["righting_arm_m"] for row in wall["rows"]], "o",
              label="heeled displacement integration", markersize=5)
    left.plot(heels, [row["closed_form_m"] for row in wall["rows"]], "-",
              label="wall-sided closed form", linewidth=1.4)
    left.set_xlabel("heel angle [deg]")
    left.set_ylabel("righting arm GZ [m]")
    left.set_title(f"Wall-sided box: GZ to {max(heels):.0f} deg\n"
                   f"max |difference| = {wall['max_abs_difference_m']:.1e} m")
    left.legend(loc="lower right", fontsize=8)
    left.grid(alpha=0.3)

    for row in barge["ladder"]:
        right.plot(row["hea_angles_deg"], row["righting_arms_m"],
                   label=f"KG = {row['kg_m']:.2f} m (GM = {row['gm_m']:.3f} m)",
                   linewidth=1.4)
    right.axhline(0.0, color="black", linewidth=0.8)
    right.set_xlabel("heel angle [deg]")
    right.set_ylabel("righting arm GZ [m]")
    right.set_title("Barge: the range of stability shrinks with a rising KG")
    right.legend(fontsize=8)
    right.grid(alpha=0.3)
    figure.tight_layout()
    path = output / "gz_curves.png"
    figure.savefig(path, dpi=140)
    plt.close(figure)
    written.append(path.name)

    figure, axis = plt.subplots(figsize=(5.6, 4.0))
    sides = [row["sides"] for row in circular["rows"]]
    errors = [row["relative_error"] for row in circular["rows"]]
    axis.loglog(sides, errors, "o-", label="measured")
    guide = [errors[0] * (sides[0] / value) ** 2 for value in sides]
    axis.loglog(sides, guide, "--", label=r"$1/n^2$ guide", linewidth=1.0)
    axis.set_xlabel("polygon sides")
    axis.set_ylabel(r"max relative error in $GZ/\sin\phi$")
    axis.set_title("Inscribed prism -> circular section (KM = R)")
    axis.legend(fontsize=8)
    axis.grid(alpha=0.3, which="both")
    figure.tight_layout()
    path = output / "gz_convergence.png"
    figure.savefig(path, dpi=140)
    plt.close(figure)
    written.append(path.name)
    return written


# ---------------------------------------------------------------------------
# markdown report
# ---------------------------------------------------------------------------
def render_markdown(measurements: dict[str, Any], figures: list[str]) -> str:
    wall = measurements["wall_sided"]
    circular = measurements["circular"]
    barge = measurements["barge"]
    lines = [
        "# Righting arms, measured against two closed forms",
        "",
        f"TensorCAD {measurements['version']} -- generated by `examples/stability/run_gz.py`,",
        "which exits non-zero if a claim below stops holding.",
        "",
        "## 1. Wall-sided box",
        "",
        f"L = {wall['geometry']['length']} m, B = {wall['geometry']['beam']} m, "
        f"T = {wall['geometry']['draft']} m, freeboard = {wall['geometry']['freeboard']} m, "
        f"KG = {wall['geometry']['kg']} m.",
        f"BM = B^2/(12T) = {wall['bm_m']:.6f} m, GM = KB + BM - KG = {wall['gm_m']:.6f} m.",
        f"Closed form valid until a deck corner immerses at atan(2f/B) = "
        f"{wall['deck_corner_immersion_deg']:.3f} deg.",
        "",
        "| heel [deg] | GZ from the integration [m] | closed form [m] | difference [m] | sinkage [m] |",
        "|---|---|---|---|---|",
    ]
    for row in wall["rows"]:
        lines.append(
            f"| {row['heel_angle_deg']:.0f} | {row['righting_arm_m']: .9f} | {row['closed_form_m']: .9f} | "
            f"{row['difference_m']: .2e} | {row['sinkage_m']: .1e} |"
        )
    lines += [
        "",
        f"Max |difference| = {wall['max_abs_difference_m']:.2e} m over 0-60 deg, "
        f"max relative volume residual {wall['max_volume_residual']:.1e}.",
        f"Claim `< 1e-12 m`: **{'PASS' if wall['claim_difference_below_1e_12'] else 'FAIL'}**.",
        "",
        "Past the deck-corner angle the closed form is the wrong reference, and the two "
        "part company -- which is the check that the solver is clipping the body at the "
        "waterline rather than evaluating a formula:",
        "",
        "| heel [deg] | integration [m] | closed form [m] | difference [m] |",
        "|---|---|---|---|",
    ]
    for angle, row in wall["beyond_the_limit"].items():
        lines.append(
            f"| {float(angle):.0f} | {row['righting_arm_m']: .6f} | {row['closed_form_m']: .6f} | "
            f"{row['difference_m']: .2e} |"
        )
    lines += [
        "",
        "## 2. Inscribed prism against a circular section",
        "",
        f"R = {circular['radius_m']} m, KG = {circular['kg_m']} m, KM = R = {circular['km_m']} m, "
        f"so the reference is GZ/sin(phi) = {circular['expected_gz_over_sine_m']:.6f} m at every angle "
        f"({len(circular['heel_angles_deg'])} angles, {min(circular['heel_angles_deg']):.0f}-"
        f"{max(circular['heel_angles_deg']):.0f} deg).",
        "",
        "| sides | max relative error | volume residual |",
        "|---|---|---|",
    ]
    for row in circular["rows"]:
        lines.append(f"| {row['sides']} | {row['relative_error']:.3e} | {row['max_volume_residual']:.1e} |")
    lines += [
        "",
        "Error ratios per doubling of the side count: "
        + ", ".join(f"{ratio:.2f}" for ratio in circular["error_ratios_per_doubling"])
        + " (second order would be 4).",
        f"Claim `second order`: **{'PASS' if circular['claim_second_order'] else 'FAIL'}**; "
        f"claim `finest < 5e-4`: **{'PASS' if circular['claim_finest_below_5e_4'] else 'FAIL'}**.",
        "",
        "## 3. A barge, where the geometry is the only input",
        "",
        f"L = {barge['geometry']['length']} m, B = {barge['geometry']['beam']} m, "
        f"T = {barge['geometry']['draft']} m, freeboard = {barge['geometry']['freeboard']} m.",
        "",
        "| KG [m] | GM [m] | max arm [m] | at [deg] | vanishing angle [deg] |",
        "|---|---|---|---|---|",
    ]
    for row in barge["ladder"]:
        vanishing = "not reached" if row["vanishing_angle_deg"] is None else f"{row['vanishing_angle_deg']:.1f}"
        lines.append(
            f"| {row['kg_m']:.2f} | {row['gm_m']: .4f} | {row['max_arm_m']:.4f} | "
            f"{row['angle_of_max_arm_deg']:.0f} | {vanishing} |"
        )
    lines += [
        "",
        f"Claim `vanishing angle falls as KG rises`: "
        f"**{'PASS' if barge['claim_vanishing_angle_falls_with_kg'] else 'FAIL'}**.",
        "",
    ]
    if figures:
        lines += ["## Figures", ""]
        lines += [f"* `{name}`" for name in figures]
        lines.append("")
    lines += [
        "## Scope",
        "",
        "* No free-surface correction, no damage stability, no statutory criterion "
        "(`IS Code 2008`, weather criterion) and no downflooding opening: the arm is a "
        "solid-body, single-compartment quantity.",
        "* The body handed to the solver must be the *whole* closed body including "
        "freeboard. A wetted-surface mesh closed by a cap at the design waterline has "
        "zero freeboard, and heeling it drowns the deck: the arm is then negative and "
        "that is a statement about the body, not about the ship.",
        "* Trim is not solved: the hull is held at even keel, so the arm is transverse "
        "only and `x_G` is not used.",
        "",
    ]
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("docs/stability"))
    arguments = parser.parse_args(argv)
    arguments.output.mkdir(parents=True, exist_ok=True)

    measurements: dict[str, Any] = {
        "schema": SCHEMA,
        "version": __version__,
        "wall_sided": measure_wall_sided(),
        "circular": measure_circular(),
        "barge": measure_barge(),
    }
    figures = draw_figures(measurements, arguments.output)
    measurements["figures"] = figures

    (arguments.output / "gz_results.json").write_text(
        json.dumps(measurements, indent=2, sort_keys=True) + "\n"
    )
    (arguments.output / "gz_report.md").write_text(render_markdown(measurements, figures) + "\n")

    wall = measurements["wall_sided"]
    circular = measurements["circular"]
    barge = measurements["barge"]
    print(f"wall-sided box: max |GZ - closed form| = {wall['max_abs_difference_m']:.2e} m "
          f"over 0-60 deg (claim < 1e-12: {wall['claim_difference_below_1e_12']})")
    print(f"circular section: relative error at 192 sides = {circular['finest_relative_error']:.2e}, "
          f"ratios {[f'{value:.2f}' for value in circular['error_ratios_per_doubling']]}")
    for row in barge["ladder"]:
        vanishing = row["vanishing_angle_deg"]
        print(f"barge KG={row['kg_m']:.2f} GM={row['gm_m']: .4f}: "
              f"max arm {row['max_arm_m']:.4f} m at {row['angle_of_max_arm_deg']:.0f} deg, "
              f"vanishing {'n/a' if vanishing is None else f'{vanishing:.1f} deg'}")
    print(f"artifacts written to {arguments.output}")

    healthy = (
        wall["claim_difference_below_1e_12"]
        and circular["claim_second_order"]
        and circular["claim_finest_below_5e_4"]
        and barge["claim_vanishing_angle_falls_with_kg"]
    )
    return 0 if healthy else 1


if __name__ == "__main__":
    sys.exit(main())
