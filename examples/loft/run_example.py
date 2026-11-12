#!/usr/bin/env python3
"""Exact section skinning: five measured studies against closed-form references.

Every number in ``results.json`` is computed here; nothing is copied from a report.

1. **Cylinder.**  Five exact rational circles of radius ``R`` in five planes ``x = const``.
   The skin must *be* the cylinder: max radial deviation, and the area against ``2 pi R H``.
2. **Cone.**  Radii tapering linearly: the skin must be the cone ``r(x) = r0 + k x``, and
   the area must match the frustum ``pi (r0 + r1) sqrt(H^2 + (r0 - r1)^2)``.
3. **Parabola family.**  Polynomial sections: the skin reproduces ``z = y^2`` at dense
   ``(u, v)`` -- the exactness is geometric, not just at the interpolated rows.
4. **Continuity.**  Exact one-sided second-derivative jumps in the loft direction for
   ``degree_u = 3``, ``2`` and ``1``: the cubic skin is ``C2`` at round-off, the quadratic
   one jumps, the ruled one jumps in the first derivative.  A measurement that cannot
   fail proves nothing, so the reduced-degree runs are the control.
5. **Mixed degree + rotated sections.**  A degree-1 polyline next to a degree-2 ellipse,
   and ellipses rotating about the axis: both are reproduced exactly after elevation.

Writes ``results.json``, ``report.md`` and ``plot.png`` into ``--output`` and exits
non-zero if any claim stops holding.

Usage::

    python examples/loft/run_example.py [--output docs/loft]
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tensorcad.loft import loft_sections, section_deviation  # noqa: E402
from tensorcad.nurbs import NURBSCurve, circle_curve  # noqa: E402

TWO_PI = 2.0 * math.pi


def yz_circle(x0: float, radius: float, arcs: int = 3) -> NURBSCurve:
    base = circle_curve(radius, arcs=arcs)
    return NURBSCurve(
        degree=base.degree,
        knots=base.knots,
        control_points=tuple((x0, point[0], point[1]) for point in base.control_points),
        weights=base.weights,
    )


def parabola_section(x0: float) -> NURBSCurve:
    return NURBSCurve(
        degree=2,
        knots=(0.0, 0.0, 0.0, 1.0, 1.0, 1.0),
        control_points=((x0, -1.0, 1.0), (x0, 0.0, -1.0), (x0, 1.0, 1.0)),
    )


def polyline_section(points: tuple[tuple[float, float, float], ...]) -> NURBSCurve:
    count = len(points)
    interior = [index / (count - 1) for index in range(1, count - 1)]
    return NURBSCurve(
        degree=1, knots=(0.0, 0.0, *interior, 1.0, 1.0), control_points=points
    )


def ellipse_section(x0: float, semi: tuple[float, float], angle: float) -> NURBSCurve:
    base = circle_curve(1.0)
    cos_a, sin_a = math.cos(angle), math.sin(angle)
    points = []
    for point in base.control_points:
        y, z = semi[0] * point[0], semi[1] * point[1]
        points.append((x0, cos_a * y - sin_a * z, sin_a * y + cos_a * z))
    return NURBSCurve(
        degree=base.degree, knots=base.knots, control_points=tuple(points), weights=base.weights
    )


def tapered_family() -> list[NURBSCurve]:
    return [yz_circle(x, 0.5 + 0.15 * math.cos(1.3 * x)) for x in (0.0, 1.0, 2.0, 3.0, 4.0)]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="docs/loft", type=Path)
    args = parser.parse_args()
    output: Path = args.output
    output.mkdir(parents=True, exist_ok=True)

    claims: dict[str, bool] = {}
    sections_out: dict[str, object] = {}

    # -- 1. cylinder ---------------------------------------------------------
    radius, length = 0.5, 2.0
    cylinder = loft_sections([yz_circle(x, radius) for x in (0.0, 0.5, 1.0, 1.5, length)])
    radial = max(
        abs(math.hypot(*cylinder.surface.evaluate(i / 40, j / 180)[1:]) - radius)
        for i in range(41)
        for j in range(181)
    )
    area = cylinder.surface.area(order=24)
    area_analytic = TWO_PI * radius * length
    cylinder_measure = cylinder.measure(samples=101, deviation_samples=41)
    claims["cylinder_radial_deviation_lt_1e-13"] = radial < 1e-13
    claims["cylinder_area_lt_1e-12"] = abs(area - area_analytic) / area_analytic < 1e-12
    sections_out["cylinder"] = {
        "radius_m": radius,
        "length_m": length,
        "degree": list(cylinder.surface.degree),
        "control_points": [len(row) for row in cylinder.surface.control_net],
        "max_radial_deviation_m": radial,
        "area_m2": area,
        "area_analytic_m2": area_analytic,
        "area_relative_error": abs(area - area_analytic) / area_analytic,
        "section_defect_relative": cylinder_measure.relative_interpolation_defect,
        "second_derivative_jump_relative": cylinder_measure.relative_second_derivative_jump,
        "measured_continuous_order": cylinder_measure.measured_continuous_order,
    }

    # -- 2. cone -------------------------------------------------------------
    radii = {0.0: 0.6, 1.0: 0.4, 2.0: 0.2}
    cone = loft_sections([yz_circle(x, radii[x]) for x in sorted(radii)])
    cone_deviation = max(
        abs(math.hypot(*cone.surface.evaluate(i / 40, j / 180)[1:]) - (0.6 - 0.2 * cone.surface.evaluate(i / 40, j / 180)[0]))
        for i in range(41)
        for j in range(181)
    )
    cone_area = cone.surface.area(order=32)
    slant = math.sqrt(2.0 ** 2 + (0.6 - 0.2) ** 2)
    cone_area_analytic = math.pi * (0.6 + 0.2) * slant
    claims["cone_deviation_lt_1e-12"] = cone_deviation < 1e-12
    claims["cone_area_lt_1e-9"] = abs(cone_area - cone_area_analytic) / cone_area_analytic < 1e-9
    sections_out["cone"] = {
        "radii_m": [radii[x] for x in sorted(radii)],
        "stations_m": sorted(radii),
        "max_deviation_from_analytic_m": cone_deviation,
        "area_m2": cone_area,
        "area_analytic_m2": cone_area_analytic,
        "area_relative_error": abs(cone_area - cone_area_analytic) / cone_area_analytic,
    }

    # -- 3. polynomial family ------------------------------------------------
    parabola = loft_sections([parabola_section(x) for x in (0.0, 1.0, 2.0, 3.0, 4.0)])
    parabola_deviation = max(
        abs(parabola.surface.evaluate(i / 30, j / 60)[2] - parabola.surface.evaluate(i / 30, j / 60)[1] ** 2)
        for i in range(31)
        for j in range(61)
    )
    parabola_measure = parabola.measure(samples=101, deviation_samples=41)
    claims["parabola_deviation_lt_1e-13"] = parabola_deviation < 1e-13
    sections_out["parabola"] = {
        "degree": list(parabola.surface.degree),
        "max_abs_z_minus_y_squared_m": parabola_deviation,
        "section_defect_relative": parabola_measure.relative_interpolation_defect,
        "second_derivative_jump_relative": parabola_measure.relative_second_derivative_jump,
    }

    # -- 4. continuity control -----------------------------------------------
    continuity = {}
    for degree_u in (3, 2, 1):
        measurement = loft_sections(tapered_family(), degree_u=degree_u).measure()
        continuity[str(degree_u)] = {
            "measured_continuous_order": measurement.measured_continuous_order,
            "relative_first_derivative_jump": measurement.relative_first_derivative_jump,
            "relative_second_derivative_jump": measurement.relative_second_derivative_jump,
            "interior_knots": len(measurement.interior_knots_u),
        }
    claims["cubic_skin_is_c2"] = continuity["3"]["measured_continuous_order"] == 2
    claims["cubic_jump_lt_1e-12"] = continuity["3"]["relative_second_derivative_jump"] < 1e-12
    claims["quadratic_skin_jumps"] = continuity["2"]["measured_continuous_order"] == 1
    claims["ruled_skin_jumps_in_first_derivative"] = continuity["1"]["measured_continuous_order"] == 0
    sections_out["continuity"] = continuity

    # -- 5. mixed degree and rotating sections -------------------------------
    mixed = loft_sections(
        [
            polyline_section(((0.0, 0.0, 0.0), (0.0, 0.5, 0.1), (0.0, 1.0, 0.0))),
            ellipse_section(2.0, (0.5, 0.4), 0.0),
        ]
    )
    mixed_measure = mixed.measure(samples=81, deviation_samples=21)
    claims["mixed_degree_defect_lt_1e-13"] = mixed_measure.relative_interpolation_defect < 1e-13
    claims["mixed_degree_transverse_degree_2"] = mixed.surface.degree[1] == 2

    rotating = loft_sections(
        [ellipse_section(x, (0.5, 0.3), 0.35 * x) for x in (0.0, 0.7, 1.4, 2.1, 2.8)]
    )
    rotating_measure = rotating.measure(samples=91, deviation_samples=31)
    claims["rotating_sections_defect_lt_1e-13"] = (
        rotating_measure.relative_interpolation_defect < 1e-13
    )
    sections_out["mixed_and_rotating"] = {
        "mixed_degree": mixed_measure.to_dict(),
        "rotating": rotating_measure.to_dict(),
    }

    results = {
        "claims": claims,
        "all_claims_hold": all(claims.values()),
        "cylinder": sections_out["cylinder"],
        "cone": sections_out["cone"],
        "parabola": sections_out["parabola"],
        "continuity": continuity,
        "mixed_and_rotating": {
            "mixed_degree_defect_relative": mixed_measure.relative_interpolation_defect,
            "rotating_defect_relative": rotating_measure.relative_interpolation_defect,
        },
        "acceptance": {
            "surface_passes_through_sections_lt_1e-9_L": all(
                record.relative_interpolation_defect < 1e-9
                for record in (
                    cylinder_measure,
                    parabola_measure,
                    mixed_measure,
                    rotating_measure,
                )
            ),
            "G2_in_the_loft_direction": continuity["3"]["measured_continuous_order"] == 2,
        },
    }

    (output / "results.json").write_text(json.dumps(results, indent=2) + "\n")

    # -- figure -------------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

        figure = plt.figure(figsize=(13, 4.2))
        axis = figure.add_subplot(131, projection="3d")
        for i in range(0, 21):
            u = i / 20
            ring = [cylinder.surface.evaluate(u, j / 60) for j in range(61)]
            axis.plot([p[0] for p in ring], [p[1] for p in ring], [p[2] for p in ring],
                      color="tab:blue", lw=0.6)
        for x in (0.0, 0.5, 1.0, 1.5, 2.0):
            ring = [cylinder.surface.evaluate(x / 2.0, j / 120) for j in range(121)]
            axis.plot([p[0] for p in ring], [p[1] for p in ring], [p[2] for p in ring],
                      color="tab:red", lw=1.2)
        axis.set_title("exact skin of 5 circles = cylinder")
        axis.set_xlabel("x [m]"); axis.set_zlabel("y [m]")

        axis2 = figure.add_subplot(132)
        degrees = ["1", "2", "3"]
        first = [max(continuity[d]["relative_first_derivative_jump"], 1e-18) for d in degrees]
        second = [max(continuity[d]["relative_second_derivative_jump"], 1e-18) for d in degrees]
        width = 0.38
        axis2.bar([i - width / 2 for i in range(3)], first, width, label="|ΔS_u| / L")
        axis2.bar([i + width / 2 for i in range(3)], second, width, label="|ΔS_uu| / L")
        axis2.set_yscale("log")
        axis2.set_xticks(range(3))
        axis2.set_xticklabels([f"degree_u={d}" for d in degrees])
        axis2.set_title("one-sided derivative jumps at interior knots")
        axis2.axhline(1e-9, ls="--", color="k", lw=0.8)
        axis2.legend(fontsize=8)

        axis3 = figure.add_subplot(133)
        labels = ["cylinder", "parabola", "mixed deg", "rotating"]
        values = [
            max(cylinder_measure.relative_interpolation_defect, 1e-18),
            max(parabola_measure.relative_interpolation_defect, 1e-18),
            max(mixed_measure.relative_interpolation_defect, 1e-18),
            max(rotating_measure.relative_interpolation_defect, 1e-18),
        ]
        axis3.bar(range(len(values)), values, color="tab:green")
        axis3.set_yscale("log")
        axis3.set_xticks(range(len(labels)))
        axis3.set_xticklabels(labels, fontsize=8, rotation=15)
        axis3.set_title("section defect / L (relative)")
        axis3.axhline(1e-9, ls="--", color="k", lw=0.8)
        axis3.annotate("1e-9 L acceptance", xy=(0, 1e-9), xytext=(0.4, 3e-9), fontsize=7)
        figure.tight_layout()
        figure.savefig(output / "plot.png", dpi=150)
        plt.close(figure)
    except Exception as error:  # pragma: no cover - figure is best effort
        print(f"figure skipped: {error}")

    lines = [
        "# Exact section skinning -- measured results",
        "",
        f"| claim | holds |",
        f"|---|---|",
    ]
    for name, holds in claims.items():
        lines.append(f"| `{name}` | {'yes' if holds else '**NO**'} |")
    lines += [
        "",
        f"Cylinder: max radial deviation {radial:.3e} m, area rel. error "
        f"{results['cylinder']['area_relative_error']:.3e}.",
        f"Cone: max deviation {cone_deviation:.3e} m, area rel. error "
        f"{results['cone']['area_relative_error']:.3e}.",
        f"Parabola family: max |z - y^2| = {parabola_deviation:.3e} m.",
        "Continuity (relative, one-sided at interior knots): "
        + ", ".join(
            f"degree_u={d}: order {continuity[d]['measured_continuous_order']}"
            for d in ("1", "2", "3")
        ),
        "",
    ]
    (output / "report.md").write_text("\n".join(lines) + "\n")

    for name, holds in claims.items():
        print(f"{'PASS' if holds else 'FAIL'}  {name}")
    print(json.dumps({k: results[k] for k in ("acceptance",)}, indent=2))
    return 0 if all(claims.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())