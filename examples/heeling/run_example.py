"""Large-angle stability measured against a closed form that is exact for the body tested.

Everything on this page is *measured*: the script exits non-zero if one of its own claims
stops holding, so the numbers in ``docs/heeling.md`` cannot drift away from the code.

What it does, in order:

1. **The wall-sided barge.** A parallelepiped's cross-section cut by a heeled waterline
   gives ``y_B = tan(phi) b^2 / (3T)``, ``z_B = T/2 + b^2 tan^2(phi)/(6T)`` and therefore
   the textbook ``GZ = sin(phi)(GM + BM tan^2(phi)/2)`` *at every angle* while the
   waterline stays on both vertical sides. The mesh is the geometry, so the kernel's own
   error is round-off -- which is what the table reports.
2. **Where the closed form stops.** Deck-edge immersion at ``phi = atan(freeboard/b)``;
   past it the measured arm must fall *below* the wall-sided claim.
3. **Ninety degrees**, where the barge floats on its deck edge and ``KN = D/2`` exactly.
4. **The initial slope.** ``dGZ/dphi`` at zero, by central differences of the curve itself,
   against ``GM = KB + BM - KG`` assembled from an upright waterplane integral. The steps
   are refined to show the difference converging at second order.
5. **The Wagley/Wigley hull**, a curved body where the mesh is *not* the geometry: the same
   slope check, a GZ curve, and a refinement ladder that measures the chord error.

The curved hull is floated at 70 % of its design draught because the mesh is capped at the
design waterline: a deck at the waterline has no freeboard, and ``gz_curve`` refuses to
pretend otherwise. At 70 % draught the hull above the waterline is the freeboard.

Run it with ``python examples/heeling/run_example.py --output docs/heeling``.
matplotlib is optional: the JSON is written either way and the PNG is skipped.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

from tensorcad import __version__
from tensorcad.heeling import (
    barge_mesh,
    gz_curve,
    gz_initial_slope,
    submerged_properties,
    wall_sided_righting_arm,
    waterline_offset,
)
from tensorcad.mesh import total_signed_volume
from tensorcad.stability import waterplane_moments
from tensorcad.wigley import generate_wigley, reference_parameters

LENGTH, BEAM, DEPTH = 10.0, 2.0, 1.5
DRAFT = 0.8
KG = 0.5
KB = DRAFT / 2.0
BM = BEAM**2 / (12.0 * DRAFT)
GM = KB + BM - KG
DISPLACEMENT = LENGTH * BEAM * DRAFT
DECK_EDGE_DEGREES = math.degrees(math.atan((DEPTH - DRAFT) / (BEAM / 2.0)))
WALL_SIDED_HEELS = (2.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0)
SLOPE_STEPS = (0.8, 0.4, 0.2, 0.1)

claims: list[dict[str, object]] = []


def claim(name: str, ok: bool, detail: str, measured: object = None, reference: object = None) -> None:
    claims.append(
        {
            "claim": name,
            "ok": bool(ok),
            "detail": detail,
            "measured": measured,
            "reference": reference,
        }
    )


def relative_error(measured: float, reference: float) -> float:
    if reference == 0.0:
        return abs(measured)
    return abs(measured - reference) / abs(reference)


def wall_sided_table(barge) -> list[dict[str, float]]:
    curve = gz_curve(barge, WALL_SIDED_HEELS, kg=KG, displacement_volume=DISPLACEMENT)
    rows = []
    for point, arm in zip(curve.points, curve.gz):
        reference = wall_sided_righting_arm(point.heel_radians, kb=KB, bm=BM, kg=KG)
        rows.append(
            {
                "heel_degrees": point.heel_degrees,
                "gz_measured": arm,
                "gz_wall_sided": reference,
                "absolute_difference": arm - reference,
                "relative_difference": relative_error(arm, reference),
                "waterplane_area": point.waterplane_area,
                "centre_of_buoyancy": list(point.centre_of_buoyancy),
            }
        )
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("docs/heeling"))
    parser.add_argument("--no-figures", action="store_true")
    arguments = parser.parse_args()
    output = arguments.output
    output.mkdir(parents=True, exist_ok=True)

    results: dict[str, object] = {
        "tensorcad_version": __version__,
        "barge": {"length": LENGTH, "beam": BEAM, "depth": DEPTH, "draft": DRAFT, "kg": KG},
        "closed_form_barge": {"kb": KB, "bm": BM, "gm": GM, "displacement": DISPLACEMENT},
    }

    # 1. the barge: kernel error is round-off, so the table is a test of the machinery
    barge = barge_mesh(LENGTH, BEAM, DEPTH)
    results["barge_solid_volume"] = total_signed_volume(barge)
    table = wall_sided_table(barge)
    results["wall_sided_table"] = table
    worst = max(row["relative_difference"] for row in table)
    results["wall_sided_worst_relative_difference"] = worst
    claim(
        "barge GZ reproduces the wall-sided closed form at every angle below deck-edge immersion",
        worst < 1e-12,
        f"worst relative difference {worst:.3e} over {len(table)} angles",
        worst,
        1e-12,
    )

    # 2. deck-edge immersion ends the closed form's validity
    past = DECK_EDGE_DEGREES + 10.0
    past_curve = gz_curve(barge, (past,), kg=KG, displacement_volume=DISPLACEMENT)
    past_reference = wall_sided_righting_arm(math.radians(past), kb=KB, bm=BM, kg=KG)
    results["deck_edge_degrees"] = DECK_EDGE_DEGREES
    results["past_deck_edge"] = {
        "heel_degrees": past,
        "gz_measured": past_curve.gz[0],
        "gz_wall_sided": past_reference,
    }
    claim(
        "past deck-edge immersion the measured arm falls below the wall-sided claim",
        past_curve.gz[0] < past_reference,
        f"at {past:.2f} deg: measured {past_curve.gz[0]:.6f} < wall-sided {past_reference:.6f}",
        past_curve.gz[0],
        past_reference,
    )

    # 3. ninety degrees, exact: the waterplane is the L x D face
    upright = gz_curve(barge, (90.0,), kg=0.2, displacement_volume=DISPLACEMENT)
    results["ninety_degrees"] = {
        "kn_measured": upright.kn[0],
        "kn_expected": 0.5 * DEPTH,
        "gz_measured": upright.gz[0],
        "gz_expected": 0.5 * DEPTH - 0.2,
    }
    claim(
        "at 90 degrees KN equals half the depth, so GZ = D/2 - KG",
        relative_error(upright.kn[0], 0.5 * DEPTH) < 1e-12,
        f"KN {upright.kn[0]:.12f} vs D/2 {0.5 * DEPTH:.12f}",
        upright.kn[0],
        0.5 * DEPTH,
    )

    # 4. the initial slope against GM, refined in step size
    slope_rows = []
    for step in SLOPE_STEPS:
        result = gz_initial_slope(barge, kg=KG, waterline_z=DRAFT, step_degrees=step)
        slope_rows.append({"step_degrees": step, **result})
    results["barge_slope_ladder"] = slope_rows
    worst_slope = max(abs(row["slope"] - GM) for row in slope_rows)
    claim(
        "dGZ/dphi at zero reproduces GM on the barge",
        worst_slope < 1e-4,
        f"worst |slope - GM| = {worst_slope:.3e} over steps {SLOPE_STEPS}",
        worst_slope,
        1e-4,
    )
    orders = []
    for coarse, fine in zip(slope_rows, slope_rows[1:]):
        ratio = coarse["step_degrees"] / fine["step_degrees"]
        error_ratio = abs(coarse["difference"]) / max(abs(fine["difference"]), 1e-300)
        orders.append(
            {
                "step_ratio": ratio,
                "error_ratio": error_ratio,
                "observed_order": math.log(error_ratio) / math.log(ratio) if ratio > 1 else 0.0,
            }
        )
    results["slope_convergence"] = orders
    claim(
        "the slope difference is second order in the step (a central difference of a smooth curve)",
        all(row["observed_order"] > 1.8 for row in orders),
        "observed orders " + ", ".join(f"{row['observed_order']:.2f}" for row in orders),
    )

    # 5. the curved hull
    parameters = reference_parameters()
    waterline = -0.3 * parameters.draft
    mesh = generate_wigley(parameters).mesh
    target_volume = submerged_properties(mesh, (0.0, 0.0, 1.0), waterline)[0]
    wigley_slope = gz_initial_slope(mesh, kg=0.2, waterline_z=waterline, step_degrees=0.25)
    results["wigley"] = {
        "parameters": {
            "length": parameters.length,
            "beam": parameters.beam,
            "draft": parameters.draft,
            "nx": parameters.nx,
            "nz": parameters.nz,
        },
        "waterline_z": waterline,
        "displacement": target_volume,
        "solid_volume": total_signed_volume(mesh),
        "slope": wigley_slope,
        "waterplane_area_clip": submerged_properties(mesh, (0.0, 0.0, 1.0), waterline)[2],
        "waterplane_area_green": waterplane_moments(mesh, waterline).area,
    }
    claim(
        "the clipped waterplane area matches the Green-theorem waterplane area",
        relative_error(results["wigley"]["waterplane_area_clip"], results["wigley"]["waterplane_area_green"]) < 1e-12,
        f"{results['wigley']['waterplane_area_clip']:.9f} vs {results['wigley']['waterplane_area_green']:.9f}",
    )
    claim(
        "dGZ/dphi at zero reproduces GM on the curved hull",
        relative_error(wigley_slope["slope"], wigley_slope["gm"]) < 1e-4,
        f"slope {wigley_slope['slope']:.9f} vs GM {wigley_slope['gm']:.9f}",
        wigley_slope["slope"],
        wigley_slope["gm"],
    )

    heels = tuple(float(value) for value in range(0, 25, 5))
    wigley_curve = gz_curve(mesh, heels, kg=0.2, displacement_volume=target_volume)
    results["wigley_curve"] = {
        "angles_degrees": list(heels),
        "gz": list(wigley_curve.gz),
        "kn": list(wigley_curve.kn),
        "points": [point.to_dict() for point in wigley_curve.points],
    }
    conserved = max(
        abs(point.displacement_volume - target_volume) / target_volume for point in wigley_curve.points
    )
    claim(
        "every heel angle displaces the same volume (a GZ curve is constant displacement)",
        conserved < 1e-12,
        f"worst relative volume drift {conserved:.3e}",
        conserved,
        1e-12,
    )
    monotone = all(a < b for a, b in zip(wigley_curve.gz, wigley_curve.gz[1:]))
    claim(
        "the curved hull's righting arm grows monotonically over the sampled range",
        monotone,
        "GZ " + ", ".join(f"{value:.6f}" for value in wigley_curve.gz),
    )

    # the refinement ladder: the only error left is the chord error of the sides
    ladder = []
    fine = None
    for nx, nz in ((21, 9), (41, 13), (81, 25), (161, 41)):
        ladder_mesh = generate_wigley(reference_parameters(nx=nx, nz=nz)).mesh
        ladder_volume = submerged_properties(ladder_mesh, (0.0, 0.0, 1.0), waterline)[0]
        arm = gz_curve(ladder_mesh, (10.0,), kg=0.2, displacement_volume=ladder_volume).gz[0]
        ladder.append({"nx": nx, "nz": nz, "faces": len(ladder_mesh.faces), "gz_10deg": arm})
        fine = arm
    for row in ladder:
        row["difference_from_finest"] = row["gz_10deg"] - fine
        row["relative_difference"] = relative_error(row["gz_10deg"], fine)
    results["wigley_ladder"] = ladder
    shrinking = all(
        row["relative_difference"] >= next_row["relative_difference"]
        for row, next_row in zip(ladder, ladder[1:-1])
    )
    claim(
        "the curved-hull GZ converges as the mesh refines",
        shrinking and ladder[0]["relative_difference"] < 1e-2,
        "relative difference vs the 161x41 grid: "
        + ", ".join(f"{row['nx']}x{row['nz']} {row['relative_difference']:.2e}" for row in ladder[:-1]),
    )

    results["claims"] = claims
    results["all_claims_hold"] = all(bool(item["ok"]) for item in claims)
    (output / "results.json").write_text(json.dumps(results, indent=2, sort_keys=False) + "\n")

    write_report(output, results)
    if not arguments.no_figures:
        write_figure(output, table, past_curve.gz[0], past_reference, wigley_curve, ladder)

    failed = [item for item in claims if not item["ok"]]
    for item in claims:
        print(f"[{'PASS' if item['ok'] else 'FAIL'}] {item['claim']} -- {item['detail']}")
    if failed:
        print(f"{len(failed)} claim(s) failed", file=sys.stderr)
        return 1
    print(f"all {len(claims)} claims hold; evidence in {output}")
    return 0


def write_report(output: Path, results: dict) -> None:
    lines = ["# Large-angle stability (GZ) — measured evidence", ""]
    lines.append(f"TensorCAD {results['tensorcad_version']}; wall-sided barge "
                 f"L={LENGTH} m, B={BEAM} m, D={DEPTH} m, T={DRAFT} m, KG={KG} m.")
    lines.append("")
    lines.append("## Wall-sided barge vs the closed form")
    lines.append("")
    lines.append("| heel [deg] | GZ measured [m] | GZ closed form [m] | difference [m] | waterplane [m^2] |")
    lines.append("|---:|---:|---:|---:|---:|")
    for row in results["wall_sided_table"]:
        lines.append(
            f"| {row['heel_degrees']:.0f} | {row['gz_measured']:.12f} | {row['gz_wall_sided']:.12f} | "
            f"{row['absolute_difference']:.2e} | {row['waterplane_area']:.9f} |"
        )
    lines.append("")
    lines.append("## Initial slope vs GM")
    lines.append("")
    lines.append("| step [deg] | dGZ/dphi | GM [m] | difference [m] |")
    lines.append("|---:|---:|---:|---:|")
    for row in results["barge_slope_ladder"]:
        lines.append(
            f"| {row['step_degrees']} | {row['slope']:.12f} | {row['gm']:.12f} | {row['difference']:.3e} |"
        )
    lines.append("")
    lines.append("## Curved hull (Wigley at 70 % draught)")
    lines.append("")
    lines.append("| heel [deg] | KN [m] | GZ [m] | displacement [m^3] |")
    lines.append("|---:|---:|---:|---:|")
    for point, kn, arm in zip(
        results["wigley_curve"]["points"], results["wigley_curve"]["kn"], results["wigley_curve"]["gz"]
    ):
        lines.append(
            f"| {point['heel_degrees']:.0f} | {kn:.9f} | {arm:.9f} | {point['displacement_volume']:.12f} |"
        )
    lines.append("")
    lines.append("## Mesh refinement at 10 deg of heel")
    lines.append("")
    lines.append("| grid | faces | GZ(10 deg) [m] | difference from 161x41 [m] | relative |")
    lines.append("|---|---:|---:|---:|---:|")
    for row in results["wigley_ladder"]:
        lines.append(
            f"| {row['nx']}x{row['nz']} | {row['faces']} | {row['gz_10deg']:.9f} | "
            f"{row['difference_from_finest']:.2e} | {row['relative_difference']:.2e} |"
        )
    lines.append("")
    lines.append("## Claims")
    lines.append("")
    for item in results["claims"]:
        lines.append(f"- {'PASS' if item['ok'] else 'FAIL'} — {item['claim']}: {item['detail']}")
    lines.append("")
    (output / "report.md").write_text("\n".join(lines))


def write_figure(output: Path, table, past_gz: float, past_reference: float, wigley_curve, ladder) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:  # pragma: no cover - optional dependency
        print("matplotlib unavailable: skipping the figure")
        return
    figure, (left, right) = plt.subplots(1, 2, figsize=(11.0, 4.4))
    heels = [row["heel_degrees"] for row in table] + [DECK_EDGE_DEGREES + 10.0]
    measured = [row["gz_measured"] for row in table] + [past_gz]
    closed = [row["gz_wall_sided"] for row in table] + [past_reference]
    left.plot(heels, closed, "k-", label="wall-sided closed form")
    left.plot(heels, measured, "o", markerfacecolor="none", label="TensorCAD (mesh kernel)")
    left.axvline(DECK_EDGE_DEGREES, color="tab:red", linestyle=":", label="deck-edge immersion")
    left.set_xlabel("heel angle [deg]")
    left.set_ylabel("GZ [m]")
    left.set_title(f"Barge: kernel vs closed form (worst rel. diff {max(r['relative_difference'] for r in table):.1e})")
    left.grid(alpha=0.3)
    left.legend(fontsize=8)

    right.plot(wigley_curve.angles, wigley_curve.gz, "s-", color="tab:green", label="Wigley, GZ")
    right.plot(wigley_curve.angles, wigley_curve.kn, "^--", color="tab:blue", label="Wigley, KN")
    right.set_xlabel("heel angle [deg]")
    right.set_ylabel("righting arm [m]")
    right.set_title("Curved hull at 70 % draught (0.33582 m^3)")
    right.grid(alpha=0.3)
    right.legend(fontsize=8)
    figure.tight_layout()
    figure.savefig(output / "plot.png", dpi=150)
    plt.close(figure)


if __name__ == "__main__":
    raise SystemExit(main())