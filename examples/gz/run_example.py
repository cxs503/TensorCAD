"""Run: ``python examples/gz/run_example.py --output docs/gz``.

Writes ``results.json``, ``report.md`` and ``figures/gz-curve.png``, and **exits non-zero
if one of its claims stops holding**. Three artefacts, three purposes:

* the wall-sided barge is the *exactness* case: ``GZ = sin(phi)(GM + BM tan^2(phi)/2)`` is
  the textbook closed form for a hull whose waterline still cuts both side walls, so the
  measured curve must agree to rounding, not to "a few percent";
* the round bilge is the *independent code path*: the section of a circle is a circular
  segment, so its centre of buoyancy is on the vertical through the centre at every angle
  and ``KN = R sin(phi)`` exactly -- with the facetting as the only admitted error;
* the freeboarded barge carries the **IMO IS Code 2008 Part A 2.2** screening, which needs
  a curve from 0 to 40 degrees.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from tensorcad import HullParameters, MeshBuilder, generate_hull
from tensorcad.gz import gz_curve, heel_mesh, imo_intact_stability, submerged_solid
from tensorcad.stability import initial_stability

LENGTH, BEAM, DRAFT, DEPTH = 30.0, 6.0, 2.0, 10.0
KG_BARGE = 1.0
ANGLE_GRID = (0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 60.0)
CYLINDER_RADIUS, CYLINDER_LENGTH, CYLINDER_SEGMENTS = 3.0, 20.0, 360
KG_ROUND = 1.5


def box(length: float, beam: float, depth: float):
    builder = MeshBuilder()
    corner = {}
    for x in (0.0, length):
        for y in (-0.5 * beam, 0.5 * beam):
            for z in (0.0, depth):
                corner[(x, y, z)] = builder.add_vertex((x, y, z))

    def quad(p0, p1, p2, p3) -> None:
        a, b, c, d = (corner[point] for point in (p0, p1, p2, p3))
        builder.add_triangle(a, b, c)
        builder.add_triangle(a, c, d)

    lo, hi = -0.5 * beam, 0.5 * beam
    quad((0.0, lo, 0.0), (0.0, hi, 0.0), (length, hi, 0.0), (length, lo, 0.0))
    quad((0.0, lo, depth), (length, lo, depth), (length, hi, depth), (0.0, hi, depth))
    quad((0.0, lo, 0.0), (0.0, lo, depth), (0.0, hi, depth), (0.0, hi, 0.0))
    quad((length, lo, 0.0), (length, hi, 0.0), (length, hi, depth), (length, lo, depth))
    quad((0.0, lo, 0.0), (length, lo, 0.0), (length, lo, depth), (0.0, lo, depth))
    quad((0.0, hi, 0.0), (0.0, hi, depth), (length, hi, depth), (length, hi, 0.0))
    return builder.build()


def round_bilge(radius: float, length: float, segments: int):
    builder = MeshBuilder()
    left, right = [], []
    for index in range(segments):
        theta = 2.0 * math.pi * (index + 0.5) / segments
        y, z = radius * math.sin(theta), radius - radius * math.cos(theta)
        left.append(builder.add_vertex((0.0, y, z)))
        right.append(builder.add_vertex((length, y, z)))
    centre_left = builder.add_vertex((0.0, 0.0, radius))
    centre_right = builder.add_vertex((length, 0.0, radius))
    for index in range(segments):
        nxt = (index + 1) % segments
        builder.add_triangle(left[index], left[nxt], right[nxt])
        builder.add_triangle(left[index], right[nxt], right[index])
        builder.add_triangle(centre_left, left[nxt], left[index])
        builder.add_triangle(centre_right, right[index], right[nxt])
    return builder.build()


def wall_sided_gz(beam: float, draft: float, kg: float, heel_deg: float) -> float:
    kb = 0.5 * draft
    bm = beam * beam / (12.0 * draft)
    gm = kb + bm - kg
    tan_phi = math.tan(math.radians(heel_deg))
    return math.sin(math.radians(heel_deg)) * (gm + 0.5 * bm * tan_phi * tan_phi)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="docs/gz", type=Path)
    arguments = parser.parse_args()
    output: Path = arguments.output
    (output / "figures").mkdir(parents=True, exist_ok=True)

    claims: list[dict[str, object]] = []

    def claim(name: str, measured: float, reference: float, tolerance: float, note: str) -> None:
        error = abs(measured - reference)
        claims.append(
            {
                "claim": name,
                "measured": measured,
                "reference": reference,
                "abs_error": error,
                "relative_error": error / abs(reference) if reference else error,
                "tolerance": tolerance,
                "satisfied": error <= tolerance,
                "note": note,
            }
        )

    # ---- the wall-sided barge: exactness ------------------------------------------------
    barge = box(LENGTH, BEAM, DEPTH)
    curve = gz_curve(barge, kg=KG_BARGE, heel_angles_deg=ANGLE_GRID, waterline_z=DRAFT)
    # NB: ``compute_hydrostatics`` integrates the *whole* closed mesh (the package treats the
    # mesh as the wetted body), so it reports the 1800 m^3 of the barge's deck volume here.
    # The displaced volume of a waterline cut is ``submerged_solid``'s job, and that is what
    # ``gz_curve`` uses when ``displacement_volume`` is left to itself.
    displaced = curve.displacement_volume
    # wall-sided holds while tan(phi) < 2*T/B and tan(phi) < 2*(D-T)/B
    limit = min(2.0 * DRAFT / BEAM, 2.0 * (DEPTH - DRAFT) / BEAM)
    wall_sided_angles = [angle for angle in ANGLE_GRID if math.tan(math.radians(angle)) < limit - 1e-9]

    worst = 0.0
    rows = []
    for point in curve.points:
        reference = wall_sided_gz(BEAM, DRAFT, KG_BARGE, point.heel_angle_deg)
        in_range = point.heel_angle_deg in wall_sided_angles
        rows.append(
            {
                "heel_deg": point.heel_angle_deg,
                "kn_m": point.kn,
                "gz_m": point.gz,
                "wall_sided_gz_m": reference,
                "difference_m": point.gz - reference if in_range else None,
                "in_wall_sided_range": in_range,
                "waterline_z_m": point.waterline_z,
                "displacement_volume_m3": point.displacement_volume,
            }
        )
        if in_range:
            worst = max(worst, abs(point.gz - reference))
    claim(
        "wall-sided barge: max |GZ - closed form| over the wall-sided range",
        worst,
        0.0,
        1e-9,
        f"wall-sided holds for tan(phi) < {limit:.3f}, i.e. phi < {math.degrees(math.atan(limit)):.1f} deg",
    )
    for point in curve.points:
        if abs(point.displacement_volume - displaced) > 1e-9:
            raise SystemExit(
                f"displacement drifted at {point.heel_angle_deg} deg: {point.displacement_volume}"
            )

    # ---- the round bilge: an independent analytic path ----------------------------------
    round_mesh = round_bilge(CYLINDER_RADIUS, CYLINDER_LENGTH, CYLINDER_SEGMENTS)
    polygon_area = (
        0.5 * CYLINDER_SEGMENTS * CYLINDER_RADIUS**2 * math.sin(2.0 * math.pi / CYLINDER_SEGMENTS)
    )
    round_curve = gz_curve(
        round_mesh,
        kg=KG_ROUND,
        heel_angles_deg=ANGLE_GRID,
        displacement_volume=CYLINDER_LENGTH * polygon_area / 2.0,
        waterline_z=CYLINDER_RADIUS,
    )
    round_rows = []
    worst_round = 0.0
    worst_float = 0.0
    for point in round_curve.points:
        phi = math.radians(point.heel_angle_deg)
        reference = (CYLINDER_RADIUS - KG_ROUND) * math.sin(phi)
        floatation = CYLINDER_RADIUS * math.cos(phi)
        round_rows.append(
            {
                "heel_deg": point.heel_angle_deg,
                "kn_m": point.kn,
                "gz_m": point.gz,
                "analytic_gz_m": reference,
                "analytic_waterline_z_m": floatation,
                "waterline_z_m": point.waterline_z,
                "kn_error_m": point.kn - CYLINDER_RADIUS * math.sin(phi),
            }
        )
        worst_round = max(worst_round, abs(point.gz - reference))
        worst_float = max(worst_float, abs(point.waterline_z - floatation))
    claim(
        "round bilge: max |GZ - (R - KG) sin(phi)| at 360 facets",
        worst_round,
        0.0,
        2e-3,
        "the residual is the facetting of the circle, not the cut",
    )
    claim(
        "round bilge: max |floatation - R cos(phi)| (chord through the centre)",
        worst_float,
        0.0,
        2e-3,
        "a chord through the centre is the only area-bisecting chord of a circle",
    )

    # ---- IMO screening ------------------------------------------------------------------
    imo = imo_intact_stability(curve)
    if not imo.passed:
        raise SystemExit(f"the freeboarded barge fails IMO screening: {imo.to_dict()}")

    # ---- a real hull, small angles only -------------------------------------------------
    hull = generate_hull(HullParameters(length=12.0, beam=2.4, draft=0.9))
    hull_curve = gz_curve(hull.mesh, kg=0.40, heel_angles_deg=(0.0, 5.0, 10.0, 15.0, 20.0))
    hull_reference = initial_stability(hull.mesh, kg=0.40, waterline_z=0.0)
    claim(
        "hull: GM from the curve equal to initial_stability at the same KG",
        abs(hull_curve.gm_upright - hull_reference.gm_transverse),
        0.0,
        1e-12,
        "the two modules agree bit-for-bit where both are defined",
    )

    # ---- the heeled section, for the figure ---------------------------------------------
    heel_deg = 30.0
    heeled = heel_mesh(barge, heel_deg)
    solved_z, solid = (
        heeled,
        None,
    )[1] or (None, None)
    from tensorcad.gz import solve_waterline  # local: only the figure needs it

    solved_z, solid = solve_waterline(heeled, displaced)
    keel = 0.0
    angle = math.radians(heel_deg)
    corners = []
    for x in (0.0, LENGTH):
        for y in (-BEAM / 2.0, BEAM / 2.0):
            for z in (0.0, DEPTH):
                corners.append(
                    (y * math.cos(angle) + (z - keel) * math.sin(angle),
                     -y * math.sin(angle) + (z - keel) * math.cos(angle) + keel)
                )

    payload = {
        "case": {
            "barge": {"length_m": LENGTH, "beam_m": BEAM, "draft_m": DRAFT, "depth_m": DEPTH,
                      "kg_m": KG_BARGE},
            "round_bilge": {"radius_m": CYLINDER_RADIUS, "length_m": CYLINDER_LENGTH,
                            "segments": CYLINDER_SEGMENTS, "kg_m": KG_ROUND},
            "hull": {"length_m": 12.0, "beam_m": 2.4, "draft_m": 0.9, "kg_m": 0.40},
        },
        "barge": {
            "curve": curve.to_dict(),
            "table": rows,
            "gm_upright_m": curve.gm_upright,
            "displacement_volume_m3": displaced,
            "imo": imo.to_dict(),
        },
        "round_bilge": {"curve": round_curve.to_dict(), "table": round_rows},
        "hull": {"curve": hull_curve.to_dict(), "gm_upright_m": hull_curve.gm_upright},
        "heeled_section": {
            "heel_deg": heel_deg,
            "waterline_z_m": solved_z,
            "centre_buoyancy_m": solid.centroid,
            "corners_yz": corners,
        },
        "claims": claims,
        "all_claims_satisfied": all(entry["satisfied"] for entry in claims),
    }
    (output / "results.json").write_text(json.dumps(payload, indent=2) + "\n")

    # ---- figure -------------------------------------------------------------------------
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, (left_axis, right_axis) = plt.subplots(1, 2, figsize=(11.0, 4.4))
    barge_angles = [row["heel_deg"] for row in rows]
    left_axis.plot(barge_angles, [row["gz_m"] for row in rows], "o-", label="measured GZ")
    left_axis.plot(
        [row["heel_deg"] for row in rows if row["in_wall_sided_range"]],
        [row["wall_sided_gz_m"] for row in rows if row["in_wall_sided_range"]],
        "x--",
        label="wall-sided closed form",
    )
    left_axis.plot(barge_angles, [row["kn_m"] for row in rows], ":", label="KN (shape only)")
    left_axis.plot(
        [row["heel_deg"] for row in round_rows],
        [row["gz_m"] for row in round_rows],
        "-",
        label="round bilge GZ",
    )
    left_axis.plot(
        [row["heel_deg"] for row in round_rows],
        [row["analytic_gz_m"] for row in round_rows],
        "k.",
        markersize=3,
        label="round bilge analytic",
    )
    left_axis.axhline(0.0, color="0.7", linewidth=0.8)
    left_axis.set_xlabel("heel angle [deg]")
    left_axis.set_ylabel("righting arm [m]")
    left_axis.set_title("GZ curves against their closed forms")
    left_axis.grid(alpha=0.3)
    left_axis.legend(fontsize=7)

    right_axis.plot(
        [point[0] for point in corners],
        [point[1] for point in corners],
        ".",
        markersize=2,
        color="0.5",
    )
    order = [(0, 1), (2, 3), (4, 5), (6, 7), (0, 4), (1, 5), (2, 6), (3, 7), (0, 2), (1, 3), (4, 6), (5, 7)]
    for first, second in order:
        right_axis.plot(
            [corners[first][0], corners[second][0]],
            [corners[first][1], corners[second][1]],
            color="0.6",
            linewidth=0.6,
        )
    width = max(point[0] for point in corners) - min(point[0] for point in corners)
    right_axis.axhline(solved_z, color="tab:blue", linewidth=1.0)
    right_axis.plot([solid.centroid[1]], [solid.centroid[2]], "o", color="tab:red", label="B")
    right_axis.plot([0.0], [0.0], "s", color="black", label="K")
    right_axis.plot(
        [0.0, solid.centroid[1]],
        [0.0 + KG_BARGE, solved_z],
        "k:",
        linewidth=0.8,
    )
    right_axis.set_xlim(min(point[0] for point in corners) - 0.5, max(point[0] for point in corners) + 0.5)
    right_axis.set_ylim(-1.0, DEPTH + 1.0)
    right_axis.set_xlabel("y [m] (transverse)")
    right_axis.set_ylabel("z [m]")
    right_axis.set_title(f"barge at {heel_deg:.0f} deg heel: waterline z={solved_z:.3f} m")
    right_axis.grid(alpha=0.3)
    right_axis.legend(fontsize=7)
    right_axis.text(0.02, 0.02, f"KN = {solid.centroid[1]:.4f} m", transform=right_axis.transAxes, fontsize=7)
    figure.tight_layout()
    figure.savefig(output / "figures" / "gz-curve.png", dpi=140)
    plt.close(figure)

    # ---- report --------------------------------------------------------------------------
    lines = [
        "# Large-angle stability: the GZ curve from heeled plane cuts",
        "",
        "Produced by `python examples/gz/run_example.py`; every number below is read from",
        "`results.json` written by the same run, and the script exits non-zero if a claim breaks.",
        "",
        "## Case",
        "",
        f"* wall-sided barge: L={LENGTH} m, B={BEAM} m, T={DRAFT} m, depth {DEPTH} m, KG={KG_BARGE} m",
        f"* round bilge: R={CYLINDER_RADIUS} m, L={CYLINDER_LENGTH} m, {CYLINDER_SEGMENTS} facets, KG={KG_ROUND} m",
        "* hull: `generate_hull(12.0, 2.4, 0.9)`, KG=0.40 m (capped at its waterline: no freeboard)",
        "",
        "## Barge: measured against the wall-sided closed form",
        "",
        "| heel [deg] | KN [m] | GZ [m] | wall-sided GZ [m] | diff [m] | waterline z [m] |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        difference = "—" if row["difference_m"] is None else f"{row['difference_m']:.2e}"
        lines.append(
            f"| {row['heel_deg']:.0f} | {row['kn_m']:.6f} | {row['gz_m']:.6f} | "
            f"{row['wall_sided_gz_m']:.6f} | {difference} | {row['waterline_z_m']:.6f} |"
        )
    lines += [
        "",
        f"GM upright = {curve.gm_upright:.6f} m; displaced volume "
        f"{displaced:.3f} m^3 held constant at every angle.",
        "",
        "## Round bilge: measured against circular-segment geometry",
        "",
        "| heel [deg] | GZ [m] | analytic (R-KG) sin(phi) [m] | waterline z [m] | R cos(phi) [m] |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in round_rows:
        lines.append(
            f"| {row['heel_deg']:.0f} | {row['gz_m']:.6f} | {row['analytic_gz_m']:.6f} | "
            f"{row['waterline_z_m']:.6f} | {row['analytic_waterline_z_m']:.6f} |"
        )
    lines += ["", "## IMO IS Code 2008 Part A 2.2 screening", "", "| criterion | requirement | measured | met |", "|---|---|---:|:--:|"]
    for check in imo.checks:
        lines.append(
            f"| `{check.key}` | {check.requirement} | {check.measured:.4f} {check.unit} | "
            f"{'yes' if check.satisfied else 'NO'} |"
        )
    lines += ["", "## Claims", "", "| claim | measured | tolerance | met |", "|---|---:|---:|:--:|"]
    for entry in claims:
        lines.append(
            f"| {entry['claim']} | {entry['measured']:.3e} | {entry['tolerance']:.0e} | "
            f"{'yes' if entry['satisfied'] else 'NO'} |"
        )
    lines += [
        "",
        "## What this does not do",
        "",
        "* constant trim: the waterline is kept parallel to the still-water line, so this is the",
        "  classical cross curve, not a free-to-trim solve; a hull that trims as it heels will",
        "  differ, and the difference is second order in the trim angle;",
        "* no free-surface, grain or heeling-moment allowance, no damage stability, no class approval;",
        "* `generate_hull` caps the model at the design waterline, i.e. zero freeboard, so its deck",
        "  edge is under water from the first degree of heel: a GZ curve on such a model is the",
        "  no-flooding idealisation and is only quoted here over 0-20 deg. Large-angle work needs a",
        "  model with freeboard, which is what the barge fixture is;",
        "* a waterline that lies *along* mesh edges has no section at all; the module raises rather",
        "  than reporting half an area (see the closure guard in `tensorcad/gz.py`).",
        "",
    ]
    (output / "report.md").write_text("\n".join(lines))

    print(f"barge: worst wall-sided deviation {worst:.3e} m over {len(wall_sided_angles)} angles")
    print(f"round bilge: worst GZ deviation {worst_round:.3e} m, floatation {worst_float:.3e} m")
    print(f"IMO screening passed: {imo.passed}")
    print(f"hull GM: {hull_curve.gm_upright:.6f} m")
    failed = [entry["claim"] for entry in claims if not entry["satisfied"]]
    for name in failed:
        print(f"CLAIM FAILED: {name}")
    print(f"artefacts: {output}/results.json, {output}/report.md, {output}/figures/gz-curve.png")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())