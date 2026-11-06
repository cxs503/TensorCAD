"""Intact-stability criteria on two hulls, with the area machinery checked first.

Run:  python examples/intact_stability/run_example.py --output docs/intact_stability

Writes ``results.json``, ``report.md`` and three figures into the output directory and
exits non-zero if any claim below stops holding, so CI can run it as a check rather than
as a demo.

The claims, in the order they are established:

1. **The area rule is exact for what it integrates.** On ``GZ = A sin(2 phi)``, whose
   integral is ``A/2 (cos 2a - cos 2b)``, the trapezoid rule on a 0.05 deg table is within
   1e-7 m*rad, and its error falls by 4 per halving of the step (second order).
2. **A real box barge reproduces the wall-sided closed form.**
   ``int_0^theta GZ dphi = GM (1-cos theta) + BM/2 (sec theta + cos theta - 2)`` holds to
   the sampling error below the bilge-emergence limit (45 deg for this barge, above the
   40 deg the Code quotes), and a ``KG`` shift moves the area by exactly
   ``-dKG (1 - cos theta)`` -- a negative control that would fail for any area that was
   not tracking the curve.
3. **The criteria answer the box barge's real numbers.** The barge passes the initial-GM
   clause (0.167 m > 0.15 m) and fails the 0-30 deg area clause (0.025 m*rad < 0.055), i.e.
   a low-freeboard box is not an IS-Code ship -- stated as a number, not a shrug.
4. **A hull with a real waterplane is assessed end to end**: a Wigley hull at a loading
   that puts ``GM`` at the Code's 0.15 m floor is heeled to 90 deg, and every clause is
   evaluated against its threshold, with the ``theta_f`` limitation written down rather
   than assumed away.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from tensorcad import mesh as meshlib
from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.gz import righting_arm_curve, wall_sided_righting_arm
from tensorcad.intact_stability import (
    IntactThresholds,
    area_between,
    assess_curve,
    assess_intact_stability,
    wall_sided_area,
)
from tensorcad.stability import initial_stability
from tensorcad.wigley import generate_wigley

LENGTH = 10.0
BREADTH = 1.6
DRAFT = 0.8
DEPTH = 2.2
KG = 0.5
KB0 = DRAFT / 2.0
BM0 = BREADTH**2 / (12.0 * DRAFT)
GM0 = KB0 + BM0 - KG


def box(minimum, maximum) -> Mesh:
    x0, y0, z0 = minimum
    x1, y1, z1 = maximum
    builder = MeshBuilder()
    corners = [
        (x0, y0, z0),
        (x1, y0, z0),
        (x1, y1, z0),
        (x0, y1, z0),
        (x0, y0, z1),
        (x1, y0, z1),
        (x1, y1, z1),
        (x0, y1, z1),
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


def angles_to(limit: float, step: float) -> tuple[float, ...]:
    return tuple(index * step for index in range(round(limit / step) + 1))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="docs/intact_stability", type=Path)
    parser.add_argument("--plots", action="store_true", help="write the PNG figures")
    args = parser.parse_args()
    out: Path = args.output
    out.mkdir(parents=True, exist_ok=True)

    claims: list[dict[str, object]] = []

    def claim(name: str, ok: bool, detail: str, numbers: dict[str, object] | None = None) -> None:
        claims.append({"claim": name, "holds": bool(ok), "detail": detail, "numbers": numbers or {}})

    # ------------------------------------------------------------------ 1. area rule
    amplitude = 0.5
    sine_results = []
    for step in (1.0, 0.5, 0.25, 0.1, 0.05):
        angles = angles_to(90.0, step)
        arms = tuple(amplitude * math.sin(2.0 * math.radians(a)) for a in angles)
        measured = area_between(angles, arms, 0.0, 40.0)
        exact = 0.5 * amplitude * (1.0 - math.cos(2.0 * math.radians(40.0)))
        sine_results.append({"step_deg": step, "area": measured, "exact": exact, "error": measured - exact})
    errors = [abs(row["error"]) for row in sine_results]
    # Least-squares slope of log(error) against log(step): 2.0 for a second-order rule.
    logs = [(math.log(row["step_deg"]), math.log(abs(row["error"]))) for row in sine_results]
    mean_x = sum(x for x, _ in logs) / len(logs)
    mean_y = sum(y for _, y in logs) / len(logs)
    slope = sum((x - mean_x) * (y - mean_y) for x, y in logs) / sum((x - mean_x) ** 2 for x, _ in logs)
    claim(
        "area rule is exact for A sin(2 phi) and second order in the step",
        errors[-1] < 1e-7 and 1.85 < slope < 2.15,
        f"error at 0.05 deg step {errors[-1]:.2e} m*rad; fitted log-log slope {slope:.4f} (2.0 = second order)",
        {"sine_table": sine_results, "fitted_slope": slope},
    )

    # ------------------------------------------------------- 2. box barge vs closed form
    barge = box((0.0, -BREADTH / 2, 0.0), (LENGTH, BREADTH / 2, DEPTH))
    barge_submerged = box((0.0, -BREADTH / 2, 0.0), (LENGTH, BREADTH / 2, DRAFT))
    gm = initial_stability(barge_submerged, kg=KG).gm_transverse

    barge_curves = {}
    barge_table = []
    for step in (5.0, 2.5, 1.25):
        angles = angles_to(40.0, step)
        curve = righting_arm_curve(barge, angles_deg=angles, kg=KG, waterline_z=DRAFT)
        barge_curves[step] = curve
        for name, upper in (("area_0_30", 30.0), ("area_0_40", 40.0)):
            measured = area_between(curve.angles_deg, curve.gz, 0.0, upper)
            exact = wall_sided_area(upper, gm=GM0, bm=BM0)
            barge_table.append(
                {
                    "step_deg": step,
                    "quantity": name,
                    "measured": measured,
                    "closed_form": exact,
                    "error": measured - exact,
                }
            )
    worst = max(abs(row["error"]) for row in barge_table)
    claim(
        "box barge areas reproduce the wall-sided closed form below the 45 deg limit",
        worst < 5e-4,
        f"worst |area error| over three steps and both intervals: {worst:.2e} m*rad "
        f"(clause thresholds are 0.055 and 0.090 m*rad)",
        {"barge_areas": barge_table, "gm_from_waterplane_m": gm, "gm_closed_form_m": GM0},
    )

    # negative control: KG shift
    base = righting_arm_curve(barge, angles_deg=angles_to(30.0, 0.5), kg=KG, waterline_z=DRAFT)
    shifted = righting_arm_curve(barge, angles_deg=angles_to(30.0, 0.5), kg=KG + 0.1, waterline_z=DRAFT)
    delta = area_between(shifted.angles_deg, shifted.gz, 0.0, 30.0) - area_between(
        base.angles_deg, base.gz, 0.0, 30.0
    )
    expected_delta = -0.1 * (1.0 - math.cos(math.radians(30.0)))
    claim(
        "negative control: +0.1 m KG moves the 0-30 deg area by -dKG(1-cos 30 deg)",
        abs(delta - expected_delta) < 2e-7,
        f"measured {delta:.9f} m*rad vs closed form {expected_delta:.9f} m*rad",
        {"measured": delta, "closed_form": expected_delta},
    )

    # ------------------------------------------------- 3. criteria on the box barge
    barge_assessment = assess_curve(barge_curves[2.5], initial_gm=gm)
    barge_payload = barge_assessment.to_dict()
    gm_clause = barge_assessment.criterion("IS-2008-PartA-2.2.6")
    area_clause = barge_assessment.criterion("IS-2008-PartA-2.2.1")
    claim(
        "box barge passes the initial-GM clause and fails the 0-30 deg area clause",
        gm_clause.satisfied is True and area_clause.satisfied is False,
        f"GM {gm_clause.achieved:.4f} m vs {gm_clause.required} m; area to 30 deg "
        f"{area_clause.achieved:.4f} m*rad vs {area_clause.required} m*rad",
        {"barge_assessment": barge_payload},
    )

    # ---------------------------------------------------------------- 4. a Wigley hull
    wigley = generate_wigley()
    wigley_mesh = wigley.mesh
    # Loading chosen so that GM sits on the Code's 0.15 m floor. KG is measured above the
    # keel and GM is linear in KG, so the loading is read off one upright solve: KG = KM - GM.
    t = wigley.parameters.draft
    km = initial_stability(wigley_mesh, kg=0.0).km_transverse  # G on the keel: GM = KM
    kg_floor = km - 0.15
    wigley_stability = initial_stability(wigley_mesh, kg=kg_floor)
    wigley_angles = angles_to(90.0, 5.0)
    wigley_curve = righting_arm_curve(wigley_mesh, angles_deg=wigley_angles, kg=kg_floor)
    wigley_assessment = assess_curve(wigley_curve, initial_gm=wigley_stability.gm_transverse)
    wigley_payload = wigley_assessment.to_dict()
    gm_wigley = wigley_assessment.criterion("IS-2008-PartA-2.2.6")
    claim(
        "Wigley hull: GM is put on the 0.15 m floor and every clause is evaluated",
        abs(gm_wigley.achieved - 0.15) < 1e-6 and len(wigley_assessment.criteria) == 6,
        "GM {:.4f} m; failed clauses: {}; undetermined: {}".format(
            gm_wigley.achieved,
            ", ".join(wigley_assessment.failed) or "none",
            ", ".join(wigley_assessment.undetermined) or "none",
        ),
        {"wigley_assessment": wigley_payload},
    )

    # A deck-edge (flooding-angle) variant, to show the 40 deg clauses being truncated.
    flooding = 25.0
    truncated = assess_intact_stability(
        wigley_curve.angles_deg,
        wigley_curve.gz,
        initial_gm=wigley_stability.gm_transverse,
        flooding_angle_deg=flooding,
        thresholds=IntactThresholds(),
    )
    claim(
        f"a {flooding:g} deg flooding angle shortens the area clauses and empties 2.2.3",
        truncated.criterion("IS-2008-PartA-2.2.3").status == "undetermined"
        and truncated.criterion("IS-2008-PartA-2.2.2").achieved is not None,
        "2.2.2 area measured to the {:.0f} deg flooding angle = {:.6f} m*rad; 2.2.3 reported as {}".format(
            flooding,
            truncated.criterion("IS-2008-PartA-2.2.2").achieved,
            truncated.criterion("IS-2008-PartA-2.2.3").status,
        ),
        {"flooding_angle_deg": flooding, "assessment": truncated.to_dict()},
    )

    all_hold = all(row["holds"] for row in claims)

    results = {
        "schema": "tensorcad.examples.intact_stability/1",
        "barge": {
            "length_m": LENGTH,
            "breadth_m": BREADTH,
            "draft_m": DRAFT,
            "depth_m": DEPTH,
            "kg_m": KG,
            "kb_m": KB0,
            "bm_m": BM0,
            "gm_closed_form_m": GM0,
            "gm_from_waterplane_m": gm,
        },
        "wigley": {
            "length_m": wigley.parameters.length,
            "beam_m": wigley.parameters.beam,
            "draft_m": t,
            "kg_m": kg_floor,
            "gm_m": wigley_stability.gm_transverse,
            "kb_m": wigley_stability.kb,
            "bm_m": wigley_stability.bm_transverse,
            "displacement_volume_m3": wigley_stability.displacement_volume,
            "max_volume_error": wigley_curve.max_volume_error,
            "gz_curve": {
                "heel_angles_deg": list(wigley_curve.angles_deg),
                "gz_m": list(wigley_curve.gz),
            },
        },
        "claims": claims,
        "all_claims_hold": all_hold,
    }
    (out / "results.json").write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")

    report_lines = ["# Intact-stability criteria: verification record", ""]
    report_lines.append(
        "Generated by `examples/intact_stability/run_example.py`; every number below is "
        "read from `results.json` by that script, none is typed in by hand."
    )
    report_lines.append("")
    report_lines.append("## Claims")
    report_lines.append("")
    report_lines.append("| claim | holds | detail |")
    report_lines.append("|---|---|---|")
    for row in claims:
        report_lines.append(f"| {row['claim']} | {'yes' if row['holds'] else 'NO'} | {row['detail']} |")
    report_lines.append("")
    report_lines.append("## Criteria, clause by clause")
    report_lines.append("")
    report_lines.append("| hull | clause | required | achieved | status | note |")
    report_lines.append("|---|---|---|---|---|---|")
    for label, assessment in (("box barge", barge_assessment), ("Wigley", wigley_assessment)):
        for result in assessment.criteria:
            achieved = "--" if result.achieved is None else f"{result.achieved:.6g}"
            report_lines.append(
                f"| {label} | {result.criterion} ({result.quantity}) | "
                f"{result.required:g} {result.units} | {achieved} | {result.status} | {result.note} |"
            )
    report_lines.append("")
    report_lines.append("## Box barge areas against the wall-sided closed form")
    report_lines.append("")
    report_lines.append("| step [deg] | quantity | measured [m*rad] | closed form [m*rad] | error [m*rad] |")
    report_lines.append("|---|---|---|---|---|")
    for row in barge_table:
        report_lines.append(
            f"| {row['step_deg']:g} | {row['quantity']} | {row['measured']:.9f} | "
            f"{row['closed_form']:.9f} | {row['error']:+.2e} |"
        )
    report_lines.append("")
    report_lines.append(
        "The comparison stops at 40 deg because the barge is wall-sided only up to "
        f"{math.degrees(math.atan(2.0 * DRAFT / BREADTH)):.2f} deg (bilge emergence), which is "
        "above the Code's 40 deg boundary -- deliberately: a comparison that stopped at the "
        "first interesting angle would prove less."
    )
    report_lines.append("")
    (out / "report.md").write_text("\n".join(report_lines) + "\n")

    if args.plots:
        write_figures(
            out,
            barge_curves[2.5],
            wigley_curve,
            wigley_assessment,
            barge_assessment,
            sine_results,
            barge_table,
            delta,
            expected_delta,
            GM0,
            BM0,
        )

    print(f"claims: {len(claims)}, all hold: {all_hold}")
    for row in claims:
        print(("  OK   " if row["holds"] else "  FAIL ") + row["claim"] + " -- " + str(row["detail"]))
    print(f"wrote {out / 'results.json'} and {out / 'report.md'}")
    return 0 if all_hold else 1


def write_figures(
    out,
    barge_curve,
    wigley_curve,
    wigley_assessment,
    barge_assessment,
    sine_results,
    barge_table,
    delta,
    expected_delta,
    gm0,
    bm0,
) -> None:
    """Three figures: the areas being measured, the verdicts, and the convergence."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Figure 1: the barge curve, the closed form, and the two shaded areas.
    figure, axes = plt.subplots(figsize=(7.2, 4.4))
    angles = list(barge_curve.angles_deg)
    arms = list(barge_curve.gz)
    closed = [wall_sided_righting_arm(a, gm=gm0, bm=bm0) for a in angles]
    axes.plot(angles, arms, "o-", label="integrated curve (box barge mesh)", markersize=3)
    axes.plot(angles, closed, "--", label="wall-sided closed form", linewidth=1)
    shade = [a for a in angles if a <= 40.0]
    shade_arms = arms[: len(shade)]
    axes.fill_between(
        [a for a in shade if a <= 30.0],
        [v for a, v in zip(shade, shade_arms, strict=True) if a <= 30.0],
        alpha=0.18,
        color="tab:blue",
        label="area 0-30 deg (clause 2.2.1)",
    )
    axes.fill_between(
        [a for a in shade if 30.0 <= a <= 40.0],
        [v for a, v in zip(shade, shade_arms, strict=True) if 30.0 <= a <= 40.0],
        alpha=0.28,
        color="tab:orange",
        label="area 30-40 deg (clause 2.2.3)",
    )
    axes.axhline(0.0, color="black", linewidth=0.6)
    axes.set_xlabel("heel angle [deg]")
    axes.set_ylabel("GZ [m]")
    axes.set_title("Box barge: GZ curve, closed form, and the areas the Code measures")
    axes.legend(fontsize=7)
    axes.grid(alpha=0.25)
    figure.tight_layout()
    figure.savefig(out / "figure1_barge_areas.png", dpi=140)
    plt.close(figure)

    # Figure 2: achieved vs required per clause, for both hulls.
    figure, axes = plt.subplots(figsize=(7.2, 4.4))
    labels = [result.criterion.replace("IS-2008-PartA-", "") for result in barge_assessment.criteria]
    width = 0.38
    positions = range(len(labels))
    for offset, (name, assessment, colour) in enumerate(
        (("box barge", barge_assessment, "tab:blue"), ("Wigley", wigley_assessment, "tab:green"))
    ):
        ratios = []
        for result in assessment.criteria:
            if result.achieved is None or result.required == 0.0:
                ratios.append(0.0)
            else:
                ratios.append(result.achieved / result.required)
        axes.bar(
            [p + (offset - 0.5) * width for p in positions],
            ratios,
            width=width,
            label=name,
            color=colour,
            alpha=0.85,
        )
    axes.axhline(1.0, color="black", linestyle="--", linewidth=1, label="requirement")
    axes.set_xticks(list(positions))
    axes.set_xticklabels(labels)
    axes.set_ylabel("achieved / required")
    axes.set_title("IS Code 2008 Part A 2.2: achieved over required (1.0 = exactly on the limit)")
    axes.legend(fontsize=7)
    axes.grid(axis="y", alpha=0.25)
    figure.tight_layout()
    figure.savefig(out / "figure2_criteria_ratios.png", dpi=140)
    plt.close(figure)

    # Figure 3: convergence of the area rule, plus the negative control.
    figure, axes = plt.subplots(1, 2, figsize=(9.0, 3.8))
    steps = [row["step_deg"] for row in sine_results]
    errors = [abs(row["error"]) for row in sine_results]
    axes[0].loglog(steps, errors, "o-", label="0-40 deg area error")
    reference = [errors[0] * (step / steps[0]) ** 2 for step in steps]
    axes[0].loglog(steps, reference, "k--", linewidth=1, label="slope -2 reference")
    axes[0].set_xlabel("heel step [deg]")
    axes[0].set_ylabel("|error| [m*rad]")
    axes[0].set_title("Second-order convergence on A sin(2 phi)")
    axes[0].legend(fontsize=7)
    axes[0].grid(alpha=0.25, which="both")
    axes[1].bar(["measured", "closed form"], [delta, expected_delta], color=["tab:red", "tab:gray"])
    axes[1].set_ylabel("change in 0-30 deg area [m*rad]")
    axes[1].set_title("Negative control: KG +0.1 m")
    axes[1].grid(axis="y", alpha=0.25)
    figure.tight_layout()
    figure.savefig(out / "figure3_convergence_control.png", dpi=140)
    plt.close(figure)


if __name__ == "__main__":
    raise SystemExit(main())
