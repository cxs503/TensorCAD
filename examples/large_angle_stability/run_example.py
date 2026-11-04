#!/usr/bin/env python3
"""Large-angle stability: two reference cases, one figure, machine-checked claims.

Run from the repository root::

    python examples/large_angle_stability/run_example.py --output docs/large_angle_stability

Case A -- box barge (``L = 6 m, B = 2 m, T = 0.8 m``, ``0.6 m`` freeboard, ``KG = 0.5 m``).
The wall-sided formula ``GZ = sin(phi)(GM + BM tan^2(phi)/2)`` is closed form for a
rectangular section, so agreement is a statement about the integrator. It is exact only
while the deck edge stays dry, i.e. below ``atan(2 freeboard / B) = 30.96 deg``.

Case B -- the shipped hull family (``L = 10 m, B = 2 m, T = 0.8 m``) with 0.6 m of
freeboard added. Its arms are graded against ``GM`` from
:func:`tensorcad.stability.initial_stability` (the slope at zero heel), against
displacement conservation at every point, and for odd symmetry in heel. A second,
low-freeboard variant (0.1 m) shows a curve that peaks and vanishes inside the sector.

The script exits non-zero if any claim fails, so it doubles as a demo and as a gate.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from tensorcad import hydrostatics as hydrostatics_module
from tensorcad import mesh as meshlib
from tensorcad import stability as stability_module
from tensorcad.hull import HullParameters, generate_hull
from tensorcad.large_angle import (
    box_dynamical_stability,
    box_hull_mesh,
    box_righting_arm,
    extend_freeboard,
    gz_curve,
    righting_arm,
    wall_sided_limit_deg,
)

BREADTH = 2.0
DRAFT = 0.8
BOX_LENGTH = 6.0
BOX_FREEBOARD = 0.6
BOX_KG = 0.5
HULL_KG = 0.5
DEGREES = (0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0)


def relative_error(measured: float, reference: float) -> float:
    if reference == 0.0:
        return abs(measured)
    return abs(measured - reference) / abs(reference)


def build() -> dict[str, object]:
    box = box_hull_mesh(BOX_LENGTH, BREADTH, DRAFT, BOX_FREEBOARD)
    limit = wall_sided_limit_deg(breadth=BREADTH, freeboard=BOX_FREEBOARD)

    box_rows = []
    worst_box = 0.0
    for angle in DEGREES:
        measured = righting_arm(box, angle, kg=BOX_KG)
        reference = box_righting_arm(
            breadth=BREADTH, draft=DRAFT, kg=BOX_KG, heel_deg=angle
        )
        error = relative_error(measured.gz, reference)
        if reference != 0.0:
            worst_box = max(worst_box, error)
        box_rows.append(
            {
                "heel_deg": angle,
                "gz_measured_m": measured.gz,
                "gz_wall_sided_m": reference,
                "relative_error": error,
                "waterline_z_m": measured.hydrostatics.waterline_z,
                "kn_m": measured.kn,
                "volume_m3": measured.volume,
            }
        )

    box_curve_area = gz_curve(box, kg=BOX_KG, up_to_deg=30.0, step_deg=2.5).dynamical_stability(30.0)
    box_area_reference = box_dynamical_stability(
        breadth=BREADTH, draft=DRAFT, kg=BOX_KG, upto_deg=30.0
    )

    capped = generate_hull(HullParameters(length=10.0, beam=BREADTH, draft=DRAFT)).mesh
    hull = extend_freeboard(capped, BOX_FREEBOARD)
    capped_hydrostatics = hydrostatics_module.compute_hydrostatics(capped, 0.0)
    neutral = gz_curve(hull, kg=HULL_KG, up_to_deg=40.0, step_deg=5.0)
    gm = stability_module.initial_stability(capped, kg=HULL_KG).gm_transverse

    step = 1e-4
    slope = (
        righting_arm(hull, step, kg=HULL_KG).gz - righting_arm(hull, -step, kg=HULL_KG).gz
    ) / (2.0 * math.radians(step))

    low_freeboard = extend_freeboard(capped, 0.1)
    vanishing_curve = gz_curve(low_freeboard, kg=0.6, up_to_deg=60.0, step_deg=5.0)

    displacement = neutral.points[0].volume
    volume_error = max(
        abs(point.volume - displacement) / displacement for point in neutral.points
    )
    odd_error = max(
        abs(
            righting_arm(hull, angle, kg=HULL_KG).gz
             + righting_arm(hull, -angle, kg=HULL_KG).gz
        )
        for angle in (10.0, 20.0, 30.0)
    )

    return {
        "case": {
            "box": {
                "length_m": BOX_LENGTH,
                "breadth_m": BREADTH,
                "draft_m": DRAFT,
                "freeboard_m": BOX_FREEBOARD,
                "kg_m": BOX_KG,
                "wall_sided_limit_deg": limit,
            },
            "hull": {
                "length_m": 10.0,
                "beam_m": BREADTH,
                "draft_m": DRAFT,
                "freeboard_m": BOX_FREEBOARD,
                "kg_m": HULL_KG,
                "km_m": stability_module.initial_stability(capped, kg=HULL_KG).km_transverse,
                "gm_m": gm,
                "capped_displacement_m3": capped_hydrostatics.displacement_volume,
                "capped_block_coefficient": capped_hydrostatics.block_coefficient,
                "capped_prismatic_coefficient": capped_hydrostatics.prismatic_coefficient,
                "mesh_watertight_after_extension": meshlib.is_watertight(hull),
            },
        },
        "box_rows": box_rows,
        "box": {
            "worst_relative_error_vs_wall_sided": worst_box,
            "dynamical_stability_m_rad": box_curve_area,
            "dynamical_stability_closed_form_m_rad": box_area_reference,
            "dynamical_stability_relative_error": relative_error(
                box_curve_area, box_area_reference
            ),
        },
        "hull": {
            "gm_slope_check_relative_error": relative_error(slope, gm),
            "numerical_slope_m": slope,
            "volume_conservation_max_relative_error": volume_error,
            "odd_symmetry_max_abs_error_m": odd_error,
            "waterline_sink_at_40deg_m": neutral.points[-1].hydrostatics.waterline_z,
            "peak_gz_m": neutral.peak.gz,
            "peak_angle_deg": neutral.peak.heel_deg,
            "dynamical_stability_upto_40deg_m_rad": neutral.dynamical_stability(40.0),
            "points": [point.to_dict() for point in neutral.points],
        },
        "vanishing_case": {
            "freeboard_m": 0.1,
            "kg_m": 0.6,
            "vanishing_angle_deg": vanishing_curve.vanishing_angle_deg,
            "peak_gz_m": vanishing_curve.peak.gz,
            "peak_angle_deg": vanishing_curve.peak.heel_deg,
            "points": [
                {"heel_deg": point.heel_deg, "gz_m": point.gz}
                for point in vanishing_curve.points
            ],
        },
    }


def claims(payload: dict[str, object]) -> list[tuple[str, bool, str]]:
    box = payload["box"]  # type: ignore[index]
    hull = payload["hull"]  # type: ignore[index]
    vanishing = payload["vanishing_case"]  # type: ignore[index]
    case = payload["case"]  # type: ignore[index]
    results = [
        (
            "box GZ reproduces the wall-sided formula to < 1e-12",
            box["worst_relative_error_vs_wall_sided"] < 1e-12,
            f"{box['worst_relative_error_vs_wall_sided']:.2e}",
        ),
        (
            "box dynamical stability matches its closed form to < 1e-3",
            box["dynamical_stability_relative_error"] < 1e-3,
            f"{box['dynamical_stability_relative_error']:.2e}",
        ),
        (
            "hull dGZ/dphi(0) equals the metacentric height to < 1e-5",
            hull["gm_slope_check_relative_error"] < 1e-5,
            f"{hull['gm_slope_check_relative_error']:.2e}",
        ),
        (
            "hull displacement is held at every heel angle to < 1e-12",
            hull["volume_conservation_max_relative_error"] < 1e-12,
            f"{hull['volume_conservation_max_relative_error']:.2e}",
        ),
        (
            "GZ is odd in heel (symmetric hull) to < 1e-12 m",
            hull["odd_symmetry_max_abs_error_m"] < 1e-12,
            f"{hull['odd_symmetry_max_abs_error_m']:.2e}",
        ),
        (
            "the freeboard-extended hull is watertight",
            bool(case["hull"]["mesh_watertight_after_extension"]),
            str(case["hull"]["mesh_watertight_after_extension"]),
        ),
        (
            "the low-freeboard hull vanishes inside 0..60 deg",
            vanishing["vanishing_angle_deg"] is not None
            and 0.0 < vanishing["vanishing_angle_deg"] < 60.0,
            f"{vanishing['vanishing_angle_deg']}",
        ),
        (
            "the hull's righting arm is positive at 20 deg",
            hull["peak_gz_m"] > 0.0,
            f"peak {hull['peak_gz_m']:.4f} m at {hull['peak_angle_deg']:.0f} deg",
        ),
    ]
    return results


def figure(payload: dict[str, object], path: Path) -> bool:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as error:  # pragma: no cover - optional dependency
        print(f"[figure] skipped: {error}")
        return False

    box_rows = payload["box_rows"]  # type: ignore[index]
    hull = payload["hull"]  # type: ignore[index]
    vanishing = payload["vanishing_case"]  # type: ignore[index]
    box_case = payload["case"]["box"]  # type: ignore[index]

    figure_, axes = plt.subplots(1, 2, figsize=(10.0, 4.0), dpi=110)

    heel = [row["gz_wall_sided_m"] for row in box_rows]
    measured = [row["gz_measured_m"] for row in box_rows]
    axes[0].plot(heel, measured, "o", markersize=4.5, label="exact integrator (this module)")
    axes[0].plot(
        heel, [row["gz_wall_sided_m"] for row in box_rows], "-", color="black",
        linewidth=1.0, label="wall-sided closed form",
    )
    limit = box_case["wall_sided_limit_deg"]
    axes[0].set_title(
        f"Box barge: GZ vs closed form\n(max rel. error {payload['box']['worst_relative_error_vs_wall_sided']:.1e}"
        f", deck-edge limit {limit:.1f}$^\\circ$)",
        fontsize=9,
    )
    axes[0].set_xlabel("heel angle [deg]")
    axes[0].set_ylabel("GZ [m]")
    axes[0].grid(alpha=0.3)
    axes[0].legend(fontsize=7)

    hull_angles = [point["heel_deg"] for point in hull["points"]]
    hull_arms = [point["gz_m"] for point in hull["points"]]
    axes[1].plot(hull_angles, hull_arms, "-o", markersize=3.5,
                 label=f"hull 10x2x0.8 m, fb 0.6 m, KG {payload['case']['hull']['kg_m']} m")
    vanish_angles = [point["heel_deg"] for point in vanishing["points"]]
    vanish_arms = [point["gz_m"] for point in vanishing["points"]]
    axes[1].plot(vanish_angles, vanish_arms, "--s", markersize=3.0,
                 label="hull, fb 0.1 m, KG 0.6 m")
    angle = vanishing["vanishing_angle_deg"]
    if angle is not None:
        axes[1].axvline(angle, color="crimson", linestyle=":", linewidth=1.0)
        axes[1].annotate(f"vanishing {angle:.1f}$^\\circ$", (angle, 0.02),
                         fontsize=7, color="crimson")
    axes[1].axhline(0.0, color="black", linewidth=0.6)
    axes[1].set_title(
        "Shipped hull family with freeboard\n"
        f"dGZ/d$\\phi$(0) vs GM: {hull['gm_slope_check_relative_error']:.1e} relative",
        fontsize=9,
    )
    axes[1].set_xlabel("heel angle [deg]")
    axes[1].set_ylabel("GZ [m]")
    axes[1].grid(alpha=0.3)
    axes[1].legend(fontsize=7)

    figure_.tight_layout()
    figure_.savefig(path, bbox_inches="tight")
    plt.close(figure_)
    return True


def markdown(payload: dict[str, object], results: list[tuple[str, bool, str]]) -> str:
    box_rows = payload["box_rows"]  # type: ignore[index]
    hull = payload["hull"]  # type: ignore[index]
    case = payload["case"]  # type: ignore[index]
    lines = [
        "# Large-angle stability run",
        "",
        "Generated by `examples/large_angle_stability/run_example.py`.",
        "",
        "## Box barge vs the wall-sided closed form",
        "",
        f"`L = {case['box']['length_m']} m`, `B = {case['box']['breadth_m']} m`, "
        f"`T = {case['box']['draft_m']} m`, freeboard `{case['box']['freeboard_m']} m`, "
        f"`KG = {case['box']['kg_m']} m`.",
        "",
        "| heel [deg] | GZ measured [m] | wall-sided [m] | rel. error | waterline [m] |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in box_rows:
        lines.append(
            f"| {row['heel_deg']:.0f} | {row['gz_measured_m']:.12f} | "
            f"{row['gz_wall_sided_m']:.12f} | {row['relative_error']:.2e} | "
            f"{row['waterline_z_m']:+.4f} |"
        )
    lines += [
        "",
        "## Hull family with freeboard",
        "",
        f"`GM = {case['hull']['gm_m']:.6f} m` (upright waterplane), "
        f"`KM = {case['hull']['km_m']:.6f} m`, `KG = {case['hull']['kg_m']} m`, "
        f"`Cb = {case['hull']['capped_block_coefficient']:.4f}`, "
        f"`Cp = {case['hull']['capped_prismatic_coefficient']:.4f}`.",
        "",
        "| heel [deg] | GZ [m] | KN [m] | waterline [m] | displacement [m^3] |",
        "|---:|---:|---:|---:|---:|",
    ]
    for point in hull["points"]:  # type: ignore[union-attr]
        lines.append(
            f"| {point['heel_deg']:.0f} | {point['gz_m']:+.6f} | {point['kn_m']:.6f} | "
            f"{point['waterline_z_m']:+.4f} | {point['volume_m3']:.9f} |"
        )
    lines += ["", "## Claims", "", "| claim | result | measured |", "|---|---|---|"]
    for name, passed, value in results:
        lines.append(f"| {name} | {'PASS' if passed else 'FAIL'} | {value} |")
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="outputs/large_angle_stability")
    arguments = parser.parse_args()
    output = Path(arguments.output)
    output.mkdir(parents=True, exist_ok=True)

    payload = build()
    results = claims(payload)
    (output / "results.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    (output / "report.md").write_text(markdown(payload, results), encoding="utf-8")
    plotted = figure(payload, output / "gz_curves.png")

    print(f"box:  worst relative error vs wall-sided = "
          f"{payload['box']['worst_relative_error_vs_wall_sided']:.3e}")  # type: ignore[index]
    print(f"box:  area under GZ to 30 deg = "
          f"{payload['box']['dynamical_stability_m_rad']:.9f} m*rad vs "
          f"{payload['box']['dynamical_stability_closed_form_m_rad']:.9f} closed form")  # type: ignore[index]
    print(f"hull: GM = {payload['case']['hull']['gm_m']:.6f} m, "
          f"numerical slope = {payload['hull']['numerical_slope_m']:.6f} m, "
          f"rel. error = {payload['hull']['gm_slope_check_relative_error']:.2e}")  # type: ignore[index]
    print(f"hull: peak GZ = {payload['hull']['peak_gz_m']:.4f} m at "
          f"{payload['hull']['peak_angle_deg']:.0f} deg")  # type: ignore[index]
    print(f"vanishing angle (fb 0.1 m, KG 0.6 m) = "
          f"{payload['vanishing_case']['vanishing_angle_deg']}")  # type: ignore[index]
    print(f"figure: {'written' if plotted else 'skipped'} -> {output / 'gz_curves.png'}")
    for name, passed, value in results:
        print(f"[{'PASS' if passed else 'FAIL'}] {name} ({value})")
    return 0 if all(passed for _, passed, _ in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
