"""Validation harness for :mod:`tensorcad.statutory`.

Every number in ``docs/statutory.md`` is printed by this script -- nothing in the
document is transcribed by hand. It also writes the figure the document embeds.

Run from the repository root::

    python tools/statutory_validation.py [--output docs/statutory]

The references are the two closed-form righting arms the package already derives:

* ``GZ = GM sin(phi)``                      -> ``int GZ dphi = GM (1 - cos phi)``
* ``GZ = sin(phi)(GM0 + BM0/2 tan^2 phi)``  -> ``GM0 (1 - cos phi) + BM0/2 (sec + cos - 2)``

and, for the mesh route, the same wall-sided barge built as a closed polyhedron and heeled
by :func:`tensorcad.gz.righting_arm_curve`.
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import tensorcad.intact_stability as iscode
from tensorcad import (
    INTACT_STABILITY_CRITERIA,
    area_under_gz_curve,
    box_righting_arms_analytic,
    check_intact_stability,
    gm0_from_gz_curve,
    righting_arm_curve,
)
from tensorcad import mesh as meshlib
from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.statutory import (
    AREA_30_TO_40_MIN_M_RAD,
    AREA_TO_30_MIN_M_RAD,
    AREA_TO_40_MIN_M_RAD,
    GZ_AT_30_MIN_M,
)

SINE_GM = 0.5
BREADTH, DRAFT, DEPTH, LENGTH = 6.0, 3.0, 8.0, 10.0
KG = 2.0
BM0 = BREADTH**2 / (12.0 * DRAFT)
GM0 = 0.5 * DRAFT + BM0 - KG


def sine_curve(angle: float) -> float:
    return SINE_GM * math.sin(math.radians(angle))


def sine_area(stop: float) -> float:
    return SINE_GM * (1.0 - math.cos(math.radians(stop)))


def wall_sided_curve(angle: float) -> float:
    tangent = math.tan(math.radians(angle))
    return math.sin(math.radians(angle)) * (GM0 + 0.5 * BM0 * tangent * tangent)


def wall_sided_area(stop: float) -> float:
    radians = math.radians(stop)
    return GM0 * (1.0 - math.cos(radians)) + 0.5 * BM0 * (
        (1.0 / math.cos(radians) + math.cos(radians)) - 2.0
    )


def tender_curve(angle: float) -> float:
    return 0.1 * math.sin(math.radians(angle))


def box() -> Mesh:
    builder = MeshBuilder()
    corners = [
        (0.0, -BREADTH / 2, 0.0),
        (LENGTH, -BREADTH / 2, 0.0),
        (LENGTH, BREADTH / 2, 0.0),
        (0.0, BREADTH / 2, 0.0),
        (0.0, -BREADTH / 2, DEPTH),
        (LENGTH, -BREADTH / 2, DEPTH),
        (LENGTH, BREADTH / 2, DEPTH),
        (0.0, BREADTH / 2, DEPTH),
    ]
    ids = [builder.add_vertex(point) for point in corners]
    for face in (
        (0, 3, 2),
        (0, 2, 1),
        (4, 5, 6),
        (4, 6, 7),
        (0, 1, 5),
        (0, 5, 4),
        (2, 3, 7),
        (2, 7, 6),
        (1, 2, 6),
        (1, 6, 5),
        (0, 4, 7),
        (0, 7, 3),
    ):
        builder.add_triangle(*[ids[index] for index in face])
    return meshlib.orient_outward(builder.build()).mesh


def relative_error(attained: float, reference: float) -> float:
    return abs(attained - reference) / abs(reference)


def area_table() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for name, curve, closed_form in (
        ("GZ = 0.5 sin(phi)", sine_curve, sine_area),
        ("wall-sided barge (B=6 m, T=3 m, KG=2 m)", wall_sided_curve, wall_sided_area),
    ):
        for start, stop in ((0.0, 30.0), (0.0, 40.0), (30.0, 40.0)):
            reference = closed_form(stop) - closed_form(start)
            attained = area_under_gz_curve(curve, start, stop)
            rows.append(
                {
                    "curve": name,
                    "interval_deg": f"{start:g}-{stop:g}",
                    "computed_m_rad": attained,
                    "closed_form_m_rad": reference,
                    "relative_error": relative_error(attained, reference),
                }
            )
    return rows


def gm0_table() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for name, curve, reference in (
        ("GZ = 0.5 sin(phi)", sine_curve, SINE_GM),
        ("wall-sided barge", wall_sided_curve, GM0),
    ):
        attained = gm0_from_gz_curve(curve)
        rows.append(
            {
                "curve": name,
                "computed_gm0_m": attained,
                "reference_gm0_m": reference,
                "relative_error": relative_error(attained, reference),
            }
        )
    return rows


def coarse_table_study() -> list[dict[str, object]]:
    """How the reported area uncertainty tracks the error of a coarsened table."""
    reference = wall_sided_area(30.0)
    rows: list[dict[str, object]] = []
    for step in (20.0, 10.0, 5.0, 2.5, 1.0):
        angles = [0.0]
        while angles[-1] < 50.0:
            angles.append(angles[-1] + step)
        levers = [wall_sided_curve(angle) for angle in angles]
        assessment = check_intact_stability((angles, levers))
        attained = assessment.area_to_30_m_rad or float("nan")
        rows.append(
            {
                "sample_step_deg": step,
                "samples": len(angles),
                "area_to_30_m_rad": attained,
                "trapezoid_error_m_rad": attained - reference,
                "reported_uncertainty_m_rad": assessment.check("2.2.1a").uncertainty,
                "covered": abs(attained - reference) <= max(assessment.check("2.2.1a").uncertainty, 1e-9),
            }
        )
    return rows


def verdict_table() -> list[dict[str, object]]:
    barge = box()
    angles = [0.0, 10.0, 20.0, 30.0, 40.0]
    mesh_curve = righting_arm_curve(barge, KG, angles, waterline_z=DRAFT)
    analytic = box_righting_arms_analytic(BREADTH, DRAFT, KG, angles, depth=DEPTH)
    cases = {
        "compliant analytic (GZ = 0.5 sin phi)": check_intact_stability(sine_curve),
        "wall-sided barge, analytic curve": check_intact_stability(wall_sided_curve),
        "wall-sided barge, mesh -> righting_arm_curve": check_intact_stability(mesh_curve),
        "tender control (GZ = 0.1 sin phi)": check_intact_stability(tender_curve),
    }
    rows: list[dict[str, object]] = []
    for name, assessment in cases.items():
        rows.append(
            {
                "case": name,
                "verdict": assessment.verdict(),
                "area_to_30_m_rad": assessment.area_to_30_m_rad,
                "area_to_40_m_rad": assessment.area_to_limit_m_rad,
                "area_30_to_40_m_rad": assessment.area_30_to_limit_m_rad,
                "gz_at_30_m": assessment.check("2.2.2").attained,
                "angle_of_max_gz_deg": assessment.angle_of_max_gz_deg,
                "gm0_m": assessment.gm0_m,
                "statuses": {item.code: item.status for item in assessment.checks},
            }
        )
    mesh_rows = [
        {
            "heel_deg": angle,
            "mesh_gz_m": gz,
            "closed_form_gz_m": expected,
            "absolute_error_m": abs(gz - expected),
        }
        for angle, gz, expected in zip(angles, mesh_curve.gz, analytic, strict=True)
    ]
    return rows, mesh_rows


def figure(sine_area_rows, coarse_rows, output_dir: pathlib.Path) -> pathlib.Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure_, axes = plt.subplots(1, 3, figsize=(15.0, 4.4))

    heels = [float(step) for step in range(0, 61, 1)]
    axes[0].plot(heels, [sine_curve(h) for h in heels], label="compliant: 0.5 sin(phi)")
    axes[0].plot(heels, [wall_sided_curve(h) for h in heels], label="wall-sided barge")
    axes[0].plot(heels, [tender_curve(h) for h in heels], label="tender control: 0.1 sin(phi)")
    axes[0].axhline(GZ_AT_30_MIN_M, color="k", linestyle=":", linewidth=1.0)
    axes[0].axvline(30.0, color="k", linestyle=":", linewidth=1.0)
    axes[0].annotate("2.2.2: GZ(30) >= 0.20 m", xy=(31.0, GZ_AT_30_MIN_M + 0.01), fontsize=8)
    axes[0].set_xlabel("heel angle [deg]")
    axes[0].set_ylabel("GZ [m]")
    axes[0].set_title("righting-arm curves")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.3)

    labels = ["2.2.1a\n0-30", "2.2.1b\n0-40", "2.2.1c\n30-40"]
    thresholds = [AREA_TO_30_MIN_M_RAD, AREA_TO_40_MIN_M_RAD, AREA_30_TO_40_MIN_M_RAD]
    cases = {
        "compliant 0.5 sin": check_intact_stability(sine_curve),
        "wall-sided barge": check_intact_stability(wall_sided_curve),
        "tender control": check_intact_stability(tender_curve),
    }
    width = 0.26
    for index, (name, assessment) in enumerate(cases.items()):
        values = [
            assessment.area_to_30_m_rad or 0.0,
            assessment.area_to_limit_m_rad or 0.0,
            assessment.area_30_to_limit_m_rad or 0.0,
        ]
        positions = [position + (index - 1) * width for position in range(len(labels))]
        axes[1].bar(positions, values, width=width, label=name)
    axes[1].plot(range(len(labels)), thresholds, "k_", markersize=26, label="threshold")
    axes[1].set_xticks(range(len(labels)))
    axes[1].set_xticklabels(labels, fontsize=8)
    axes[1].set_ylabel("area under GZ [m rad]")
    axes[1].set_title("section 2.2.1 areas against their thresholds")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.3, axis="y")

    steps = [row["sample_step_deg"] for row in coarse_rows]
    errors = [abs(row["trapezoid_error_m_rad"]) for row in coarse_rows]
    bands = [row["reported_uncertainty_m_rad"] for row in coarse_rows]
    axes[2].loglog(steps, errors, "o-", label="|area error| vs closed form")
    axes[2].loglog(steps, bands, "s--", label="reported uncertainty")
    axes[2].set_xlabel("table sample step [deg]")
    axes[2].set_ylabel("area uncertainty [m rad]")
    axes[2].set_title("the error bar tracks the sampling")
    axes[2].invert_xaxis()
    axes[2].legend(fontsize=8)
    axes[2].grid(alpha=0.3, which="both")

    figure_.tight_layout()
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "statutory_validation.png"
    figure_.savefig(path, dpi=150)
    plt.close(figure_)
    return path


#: statutory clause -> the sibling module's flattened label, matched on the quantity demanded.
CLAUSE_CROSSWALK = {
    "2.2.1a": "2.2.1",
    "2.2.1b": "2.2.2",
    "2.2.1c": "2.2.3",
    "2.2.2": "2.2.4",
    "2.2.3": "2.2.5",
    "2.2.4": "2.2.6",
}


def crosscheck_table() -> list[dict[str, object]]:
    """This module against :mod:`tensorcad.intact_stability` on the same curves.

    The sibling reports the same section of the same code by trapezoidal integration on the
    caller's grid and with ``gm0`` supplied; this module integrates a Hermite interpolant (or
    an adaptive quadrature) and derives ``GM0``. Two independent routes to one number are
    this package's usual way of validating a number, so the disagreement is tabulated rather
    than argued about.
    """
    angles = tuple(step * 0.25 for step in range(0, 4 * 41))
    levers = tuple(wall_sided_curve(angle) for angle in angles)
    rows: list[dict[str, object]] = []
    for stop in (30.0, 40.0):
        mine = area_under_gz_curve((angles, levers), 0.0, stop)
        theirs = iscode.gz_area_between(angles, levers, 0.0, stop)
        closed = wall_sided_area(stop)
        rows.append(
            {
                "quantity": f"area 0-{stop:g} deg [m rad]",
                "statutory": mine,
                "intact_stability": theirs,
                "closed_form": closed,
                "statutory_error": relative_error(mine, closed),
                "sibling_error": relative_error(theirs, closed),
                "difference": abs(mine - theirs),
            }
        )
    for name, curve, gm0 in (
        ("compliant 0.5 sin(phi)", sine_curve, SINE_GM),
        ("tender 0.1 sin(phi)", lambda angle: 0.1 * math.sin(math.radians(angle)), 0.1),
    ):
        table = (tuple(float(step) for step in range(41)), tuple(curve(a) for a in range(41)))
        mine = check_intact_stability(table, gm0=gm0)
        theirs = iscode.general_criteria(table, gm0_m=gm0)
        rows.append(
            {
                "quantity": f"verdict, {name}",
                "statutory": mine.verdict(),
                "intact_stability": theirs.verdict(),
                "closed_form": None,
                "statutory_error": None,
                "sibling_error": None,
                "difference": None,
                "failed_rules": sorted(item.code for item in mine.failed),
                "sibling_failed_rules": sorted(item.clause for item in theirs.failures),
                "crosswalked": sorted(CLAUSE_CROSSWALK[item.code] for item in mine.failed),
            }
        )
    return rows


def markdown(area_rows, gm0_rows, coarse_rows, verdict_rows, mesh_rows, crosscheck_rows) -> str:
    lines = ["### Area quadrature against the closed forms", ""]
    lines.append("| curve | interval [deg] | computed [m rad] | closed form [m rad] | rel. error |")
    lines.append("|---|---|---|---|---|")
    for row in area_rows:
        lines.append(
            f"| {row['curve']} | {row['interval_deg']} | {row['computed_m_rad']:.12f} | "
            f"{row['closed_form_m_rad']:.12f} | {row['relative_error']:.2e} |"
        )
    lines += ["", "### GM0 read from the slope at upright", ""]
    lines.append("| curve | computed GM0 [m] | reference [m] | rel. error |")
    lines.append("|---|---|---|---|")
    for row in gm0_rows:
        lines.append(
            f"| {row['curve']} | {row['computed_gm0_m']:.8f} | {row['reference_gm0_m']:.8f} | "
            f"{row['relative_error']:.2e} |"
        )
    lines += ["", "### Table coarsening: error and reported uncertainty", ""]
    lines.append(
        "| step [deg] | samples | area 0-30 [m rad] | error [m rad] "
        "| reported uncertainty [m rad] | covered |"
    )
    lines.append("|---|---|---|---|---|---|")
    for row in coarse_rows:
        lines.append(
            f"| {row['sample_step_deg']:g} | {row['samples']} | {row['area_to_30_m_rad']:.9f} | "
            f"{row['trapezoid_error_m_rad']:+.2e} | {row['reported_uncertainty_m_rad']:.2e} | "
            f"{row['covered']} |"
        )
    lines += ["", "### Verdicts", ""]
    for row in verdict_rows:
        statuses = ", ".join(f"{code}={status}" for code, status in row["statuses"].items())
        lines.append(f"* **{row['case']}** -> `{row['verdict']}` ({statuses})")
        lines.append(
            f"  areas {row['area_to_30_m_rad']} / {row['area_to_40_m_rad']} / "
            f"{row['area_30_to_40_m_rad']} m rad, GZ(30) = {row['gz_at_30_m']}, "
            f"GM0 = {row['gm0_m']}"
        )
    lines += ["", "### Cross-check against `tensorcad.intact_stability`", ""]
    lines.append("| quantity | `tensorcad.statutory` | `tensorcad.intact_stability` | closed form |")
    lines.append("|---|---|---|---|")
    for row in crosscheck_rows:
        if row["closed_form"] is None:
            mine_failed = ", ".join(row["failed_rules"]) or "none"
            theirs_failed = ", ".join(row["sibling_failed_rules"]) or "none"
            lines.append(
                f"| {row['quantity']} | `{row['statutory']}` (failed {mine_failed}) "
                f"| `{row['intact_stability']}` (failed {theirs_failed}) | - |"
            )
        else:
            lines.append(
                f"| {row['quantity']} | {row['statutory']:.12f} | {row['intact_stability']:.12f} | "
                f"{row['closed_form']:.12f} |"
            )
    lines += ["", "### Mesh route: GZ at the sampled heels against the closed form", ""]
    lines.append("| heel [deg] | mesh GZ [m] | closed form [m] | |error| [m] |")
    lines.append("|---|---|---|---|")
    for row in mesh_rows:
        lines.append(
            f"| {row['heel_deg']:g} | {row['mesh_gz_m']:.15f} | {row['closed_form_gz_m']:.15f} | "
            f"{row['absolute_error_m']:.2e} |"
        )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="docs/statutory")
    arguments = parser.parse_args()
    output_dir = pathlib.Path(arguments.output)

    area_rows = area_table()
    gm0_rows = gm0_table()
    coarse_rows = coarse_table_study()
    verdict_rows, mesh_rows = verdict_table()
    crosscheck_rows = crosscheck_table()
    report = markdown(area_rows, gm0_rows, coarse_rows, verdict_rows, mesh_rows, crosscheck_rows)
    path = figure(None, coarse_rows, output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "statutory_validation.json").write_text(
        json.dumps(
            {
                "criteria": {code: dict(spec) for code, spec in INTACT_STABILITY_CRITERIA.items()},
                "area_quadrature": area_rows,
                "gm0_from_slope": gm0_rows,
                "table_coarsening": coarse_rows,
                "verdicts": verdict_rows,
                "mesh_route": mesh_rows,
                "crosscheck_intact_stability": crosscheck_rows,
            },
            indent=2,
        )
        + "\n"
    )
    (output_dir / "report.md").write_text(report + "\n")
    print(report)
    print(f"\nfigure: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
