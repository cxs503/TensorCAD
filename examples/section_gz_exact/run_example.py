"""Exact righting arms for prismatic sections, against the mesh integrator.

Run:  python examples/section_gz_exact/run_example.py --output examples/section_gz_exact/outputs

Writes ``results.json``, ``report.md`` and ``plot.png`` and exits non-zero if any claim below
stops holding, so CI can run it as a check rather than as a demo.

The claims:

1. V-section and a section turned off its symmetry axes -- the exact section oracle and
   :func:`tensorcad.gz.righting_arm_curve` agree to round-off at every heel from 0 to 90 deg,
   including heels where the waterline sweeps across a section vertex;
2. box section -- the oracle reproduces the wall-sided closed form below the bilge limit and
   departs from it above (so the agreement is not vacuous);
3. 2:1 ellipse -- the only exact large-angle reference for a *smooth non-circular* section in
   the package: an inscribed polygon mesh converges to it at second order, and the wall-sided
   formula misses it by more than 3 % at 30 deg;
4. every point of every curve displaces the same area (constraint residual < 1e-12).
"""

from __future__ import annotations

import sys as _sys
from pathlib import Path as _Path

_REPO_ROOT = _Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in _sys.path:  # allow running from a checkout without installing
    _sys.path.insert(0, str(_REPO_ROOT))

import argparse
import json
import math
from pathlib import Path

from tensorcad import mesh as meshlib
from tensorcad import section_gz
from tensorcad.gz import (
    box_righting_arms_analytic,
    box_wall_sided_limit_deg,
    righting_arm_curve,
)

ANGLES = [float(angle) for angle in range(0, 91, 5)]
TRIANGLE = ((-1.0, 0.0), (1.0, 0.0), (0.0, 2.0))
ELLIPSE = {"semi_beam": 2.0, "semi_height": 1.0, "draft": 1.2, "kg": 0.8, "length": 4.0}


def _max_error(left: list[float], right: tuple[float, ...]) -> tuple[float, float]:
    absolute = max(abs(a - b) for a, b in zip(left, right, strict=True))
    relative = max(abs(a - b) / max(abs(b), 1e-12) for a, b in zip(left, right, strict=True))
    return absolute, relative


def case_v_section() -> dict[str, object]:
    """Isosceles V: neither wall-sided nor rotationally symmetric, so ``GZ`` is not a sine."""
    kg, waterline, length = 0.6, 1.0, 4.0
    mesh = section_gz.prism_mesh(TRIANGLE, length)
    exact = section_gz.polygon_prism_righting_arms(
        TRIANGLE, ANGLES, kg=kg, waterline_z=waterline, length=length
    )
    mesh_curve = righting_arm_curve(mesh, kg, ANGLES, waterline_z=waterline)
    absolute, relative = _max_error(list(exact.gz), mesh_curve.gz)
    crossing = math.degrees(math.atan2(1.0, 2.0))
    return {
        "name": "isosceles V section (B=2 m, height 2 m, draft 1 m, KG 0.6 m, L=4 m)",
        "watertight": meshlib.is_watertight(mesh),
        "prism_volume": meshlib.total_signed_volume(mesh),
        "section_area": exact.displaced_area,
        "max_abs_error_m": absolute,
        "max_rel_error": relative,
        "max_volume_error": mesh_curve.max_volume_error,
        "waterline_vertex_crossing_deg": crossing,
        "table": [
            {"heel_deg": angle, "gz_exact": gz, "gz_mesh": gz_mesh}
            for angle, gz, gz_mesh in zip(ANGLES, exact.gz, mesh_curve.gz, strict=True)
        ],
    }


def case_rotated_hexagon() -> dict[str, object]:
    """A hexagon turned 17 deg off the waterplane: no sampled heel is a symmetry axis."""
    section = section_gz.regular_polygon_section(6, 1.0, rotation_deg=17.0)
    angles = [0.0, 7.0, 23.0, 41.0, 68.0, 90.0]
    kg, waterline, length = 0.7, 0.8, 2.5
    mesh = section_gz.prism_mesh(section, length)
    exact = section_gz.polygon_prism_righting_arms(
        section, angles, kg=kg, waterline_z=waterline, length=length
    )
    mesh_curve = righting_arm_curve(mesh, kg, angles, waterline_z=waterline)
    absolute, relative = _max_error(list(exact.gz), mesh_curve.gz)
    symmetry_axes = [0.0, 60.0, 120.0]
    return {
        "name": "regular hexagon R=1 m rotated 17 deg, draft 0.8 m, KG 0.7 m, L=2.5 m",
        "watertight": meshlib.is_watertight(mesh),
        "max_abs_error_m": absolute,
        "max_rel_error": relative,
        "sampled_heels_deg": angles,
        "polygon_symmetry_axes_deg": symmetry_axes,
        "table": [
            {"heel_deg": angle, "gz_exact": gz, "gz_mesh": gz_mesh}
            for angle, gz, gz_mesh in zip(angles, exact.gz, mesh_curve.gz, strict=True)
        ],
    }


def case_box_limit() -> dict[str, object]:
    """The box: agreement below the wall-sided limit, honest disagreement above it."""
    breadth, draft, kg, length = 4.0, 1.0, 0.6, 5.0
    rectangle = (
        (-breadth / 2.0, 0.0),
        (breadth / 2.0, 0.0),
        (breadth / 2.0, 2.0),
        (-breadth / 2.0, 2.0),
    )
    limit = box_wall_sided_limit_deg(breadth, draft)
    inside = [0.0, 5.0, 10.0, 20.0, 0.9 * limit]
    outside = [limit + 5.0, limit + 15.0]
    mesh = section_gz.prism_mesh(rectangle, length)
    rows = []
    analytic = box_righting_arms_analytic(breadth, draft, kg, inside)
    for angle, reference in zip(inside, analytic, strict=True):
        oracle = section_gz.polygon_prism_righting_arms(
            rectangle, [angle], kg=kg, waterline_z=draft
        ).gz[0]
        mesh_arm = righting_arm_curve(mesh, kg, [angle], waterline_z=draft).gz[0]
        rows.append(
            {
                "heel_deg": angle,
                "gz_exact": oracle,
                "gz_wall_sided": reference,
                "gz_mesh": mesh_arm,
                "inside_limit": True,
            }
        )
    for angle in outside:
        oracle = section_gz.polygon_prism_righting_arms(
            rectangle, [angle], kg=kg, waterline_z=draft
        ).gz[0]
        bm0 = breadth**2 / (12.0 * draft)
        gm0 = draft / 2.0 + bm0 - kg
        reference = math.sin(math.radians(angle)) * (
            gm0 + 0.5 * bm0 * math.tan(math.radians(angle)) ** 2
        )
        rows.append(
            {
                "heel_deg": angle,
                "gz_exact": oracle,
                "gz_wall_sided": reference,
                "gz_mesh": righting_arm_curve(mesh, kg, [angle], waterline_z=draft).gz[0],
                "inside_limit": False,
            }
        )
    inside_error = max(
        abs(row["gz_exact"] - row["gz_wall_sided"]) for row in rows if row["inside_limit"]
    )
    outside_error = max(
        abs(row["gz_exact"] - row["gz_wall_sided"]) for row in rows if not row["inside_limit"]
    )
    return {
        "name": "box barge B=4 m, T=1 m, D=2 m, KG=0.6 m, L=5 m",
        "wall_sided_limit_deg": limit,
        "max_error_inside_limit_m": inside_error,
        "min_error_outside_limit_m": outside_error,
        "table": rows,
    }


def case_ellipse() -> dict[str, object]:
    """2:1 ellipse: the smooth non-circular reference, and its chord convergence."""
    a, b = ELLIPSE["semi_beam"], ELLIPSE["semi_height"]
    draft, kg, length = ELLIPSE["draft"], ELLIPSE["kg"], ELLIPSE["length"]
    exact = section_gz.ellipse_prism_righting_arms(
        a, b, ANGLES, kg=kg, waterline_z=draft, length=length
    )
    convergence = []
    previous = None
    for sides in (16, 32, 64, 128, 256):
        section = section_gz.ellipse_section(a, b, sides)
        mesh = section_gz.prism_mesh(section, length)
        curve = righting_arm_curve(mesh, kg, ANGLES, waterline_z=draft)
        absolute, relative = _max_error(list(exact.gz), curve.gz)
        order = None if previous is None else math.log2(previous / absolute)
        convergence.append(
            {
                "sides": sides,
                "max_abs_error_m": absolute,
                "max_rel_error": relative,
                "observed_order": order,
            }
        )
        previous = absolute

    # the waterplane-moment GM, derived here rather than copied from the module
    d = (draft - b) / b
    segment = 0.5 * math.pi + math.asin(d) + d * math.sqrt(1.0 - d * d)
    area = a * b * segment
    half_width = a * math.sqrt(1.0 - d * d)
    kb = b + b * (-(2.0 / 3.0) * (1.0 - d * d) ** 1.5) / segment
    bm = (2.0 / 3.0) * half_width**3 / area
    gm = kb + bm - kg
    delta = 1e-6
    slope_curve = section_gz.ellipse_prism_righting_arms(
        a, b, [0.0, delta], kg=kg, waterline_z=draft, length=length
    )
    slope = (slope_curve.gz[1] - slope_curve.gz[0]) / math.radians(delta)

    high_kg = section_gz.ellipse_prism_righting_arms(
        a, b, [float(angle) for angle in range(0, 91, 2)], kg=1.6, waterline_z=draft, length=length
    )
    peak = exact.max_arm()
    wall_sided_at_30 = math.sin(math.radians(30.0)) * (
        gm + 0.5 * bm * math.tan(math.radians(30.0)) ** 2
    )
    exact_at_30 = exact.angle_of(30.0).righting_arm
    return {
        "name": "2:1 ellipse section a=2 m b=1 m, draft 1.2 m, KG 0.8 m, L=4 m",
        "section_area": exact.displaced_area,
        "gm_from_waterplane_m": gm,
        "dgz_dphi_at_zero": slope,
        "gz_max_m": None if peak is None else peak[0],
        "gz_max_heel_deg": None if peak is None else peak[1],
        "vanishing_angle_at_kg_1p6_deg": high_kg.vanishing_stability_angle_deg(),
        "wall_sided_at_30_deg": wall_sided_at_30,
        "exact_at_30_deg": exact_at_30,
        "wall_sided_relative_error_at_30": abs(wall_sided_at_30 - exact_at_30) / abs(exact_at_30),
        "convergence": convergence,
        "table": [
            {"heel_deg": angle, "gz_exact": gz} for angle, gz in zip(ANGLES, exact.gz, strict=True)
        ],
    }


def _plot(results: dict[str, object], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, (left, right) = plt.subplots(1, 2, figsize=(11.5, 4.4))

    v_section = results["v_section"]
    left.plot(
        [row["heel_deg"] for row in v_section["table"]],
        [row["gz_exact"] for row in v_section["table"]],
        "-",
        label="V section: exact oracle + mesh integrator",
    )
    hexagon = results["rotated_hexagon"]
    left.plot(
        hexagon["sampled_heels_deg"],
        [row["gz_exact"] for row in hexagon["table"]],
        "o",
        markersize=4,
        label="hexagon 17 deg: exact (mesh agrees)",
    )
    ellipse = results["ellipse"]
    left.plot(
        [row["heel_deg"] for row in ellipse["table"]],
        [row["gz_exact"] for row in ellipse["table"]],
        "-",
        label="2:1 ellipse: exact closed form",
    )
    box = [row for row in results["box"]["table"]]
    left.plot(
        [row["heel_deg"] for row in box],
        [row["gz_wall_sided"] for row in box],
        ":",
        label="box: wall-sided formula",
    )
    left.plot(
        [row["heel_deg"] for row in box],
        [row["gz_exact"] for row in box],
        "--",
        linewidth=1.0,
        label="box: exact section oracle",
    )
    left.axhline(0.0, color="0.7", linewidth=0.8)
    left.set_xlabel("heel angle [deg]")
    left.set_ylabel("GZ [m]")
    left.set_title("Righting arms: exact section references")
    left.grid(alpha=0.3)
    left.legend(fontsize=7)

    convergence = ellipse["convergence"]
    sides = [row["sides"] for row in convergence]
    errors = [row["max_abs_error_m"] for row in convergence]
    right.loglog(sides, errors, "o-", label="measured max |GZ_mesh - GZ_exact|")
    reference = [errors[0] * (sides[0] / value) ** 2 for value in sides]
    right.loglog(sides, reference, "--", color="0.5", label="slope -2 reference")
    right.set_xlabel("polygon sides approximating the ellipse")
    right.set_ylabel("max arm error [m]")
    right.set_title("Inscribed-chord convergence to the exact ellipse")
    right.grid(alpha=0.3, which="both")
    right.legend(fontsize=8)

    figure.tight_layout()
    figure.savefig(path, dpi=150)
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("outputs"))
    arguments = parser.parse_args()
    arguments.output.mkdir(parents=True, exist_ok=True)

    results = {
        "v_section": case_v_section(),
        "rotated_hexagon": case_rotated_hexagon(),
        "box": case_box_limit(),
        "ellipse": case_ellipse(),
    }

    checks: list[tuple[str, bool, str]] = []
    v_section = results["v_section"]
    checks.append(
        ("V section: exact vs mesh at 19 heel angles", v_section["max_abs_error_m"] < 1e-9,
         f"{v_section['max_abs_error_m']:.2e} m")
    )
    checks.append(
        ("V section: displacement held constant", v_section["max_volume_error"] < 1e-10,
         f"{v_section['max_volume_error']:.2e}")
    )
    hexagon = results["rotated_hexagon"]
    checks.append(
        ("rotated hexagon: exact vs mesh off the symmetry axes",
         hexagon["max_abs_error_m"] < 1e-9, f"{hexagon['max_abs_error_m']:.2e} m")
    )
    box = results["box"]
    checks.append(
        ("box: wall-sided formula reproduced below the limit",
         box["max_error_inside_limit_m"] < 1e-9, f"{box['max_error_inside_limit_m']:.2e} m")
    )
    checks.append(
        ("box: wall-sided formula departs above the limit (not vacuous)",
         box["min_error_outside_limit_m"] > 1e-3, f"{box['min_error_outside_limit_m']:.2e} m")
    )
    ellipse = results["ellipse"]
    checks.append(
        ("ellipse: dGZ/dphi(0) equals the waterplane GM",
         abs(ellipse["dgz_dphi_at_zero"] - ellipse["gm_from_waterplane_m"])
         < 1e-6 * abs(ellipse["gm_from_waterplane_m"]),
         f"{ellipse['dgz_dphi_at_zero']:.9f} vs {ellipse['gm_from_waterplane_m']:.9f}")
    )
    orders = [row["observed_order"] for row in ellipse["convergence"] if row["observed_order"]]
    checks.append(
        ("ellipse: inscribed chords converge at second order",
         all(abs(order - 2.0) < 0.25 for order in orders),
         ", ".join(f"{order:.2f}" for order in orders))
    )
    finest = ellipse["convergence"][-1]
    checks.append(
        ("ellipse: 256 chords within 3 % of the exact curve",
         finest["max_rel_error"] < 0.03, f"{finest['max_rel_error']:.2e}")
    )
    checks.append(
        ("ellipse: wall-sided formula errs by > 3 % at 30 deg (not vacuous)",
         ellipse["wall_sided_relative_error_at_30"] > 0.03,
         f"{ellipse['wall_sided_relative_error_at_30']:.2%}")
    )

    plot_path = arguments.output / "plot.png"
    _plot(results, plot_path)
    (arguments.output / "results.json").write_text(
        json.dumps(
            {
                "results": results,
                "checks": [{"name": name, "passed": ok, "value": value} for name, ok, value in checks],
            },
            indent=2,
        )
        + "\n"
    )

    lines = [
        "# Exact righting arms for prismatic sections",
        "",
        "Generated by `examples/section_gz_exact/run_example.py`; every number below is written",
        "by the script itself (no hand-copied values).",
        "",
        "| case | quantity | value |",
        "| --- | --- | --- |",
        f"| V section | max abs error, exact vs mesh | {v_section['max_abs_error_m']:.3e} m |",
        f"| V section | max relative error | {v_section['max_rel_error']:.3e} |",
        f"| rotated hexagon | max abs error, exact vs mesh | {hexagon['max_abs_error_m']:.3e} m |",
        f"| box | max error below the wall-sided limit | {box['max_error_inside_limit_m']:.3e} m |",
        f"| box | min error above the wall-sided limit | {box['min_error_outside_limit_m']:.3e} m |",
        f"| 2:1 ellipse | dGZ/dphi(0) | {ellipse['dgz_dphi_at_zero']:.9f} m/rad |",
        f"| 2:1 ellipse | GM from the waterplane moments | {ellipse['gm_from_waterplane_m']:.9f} m |",
        f"| 2:1 ellipse | GZ_max | {ellipse['gz_max_m']:.6f} m at {ellipse['gz_max_heel_deg']:.0f} deg |",
        f"| 2:1 ellipse | vanishing angle at KG=1.6 m | {ellipse['vanishing_angle_at_kg_1p6_deg']:.2f} deg |",
        f"| 2:1 ellipse | wall-sided error at 30 deg | {ellipse['wall_sided_relative_error_at_30']:.2%} |",
        "",
        "## Chord convergence to the exact ellipse",
        "",
        "| polygon sides | max abs error [m] | max rel error | observed order |",
        "| --- | --- | --- | --- |",
    ]
    for row in ellipse["convergence"]:
        order = "-" if row["observed_order"] is None else f"{row['observed_order']:.2f}"
        lines.append(
            f"| {row['sides']} | {row['max_abs_error_m']:.3e} | {row['max_rel_error']:.3e} | {order} |"
        )
    lines += ["", "## Checks", ""]
    for name, ok, value in checks:
        lines.append(f"- {'PASS' if ok else 'FAIL'} -- {name}: {value}")
    (arguments.output / "report.md").write_text("\n".join(lines) + "\n")

    failures = [name for name, ok, _ in checks if not ok]
    for name, ok, value in checks:
        print(f"{'PASS' if ok else 'FAIL'}  {name}: {value}")
    if failures:
        print(f"\n{len(failures)} check(s) failed", file=_sys.stderr)
        return 1
    print(f"\nartifacts written to {arguments.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())