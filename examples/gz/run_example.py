#!/usr/bin/env python3
"""Large-angle stability of a barge and a sphere, checked against closed forms.

Run it with no arguments to refresh ``docs/gz/``::

    python examples/gz/run_example.py

Every number quoted in ``docs/gz.md`` comes from ``results.json`` written here, and the
script exits 1 if any claim stops holding. Two geometries carry the argument:

* a **rectangular barge** with 0.4 m of freeboard -- its heeled displacement has a
  closed form (wall-sided) inside ``atan(2T/B)`` and an *independent* one everywhere
  else, because a prism's centre of buoyancy is its section centroid, which this script
  integrates by clipping a quadrilateral and applying the shoelace formula;
* a **sphere** from ``SphereField``, whose righting arm is ``(KM - KG) sin(phi)`` at
  every heel with ``KM = R``, so its error is pure discretisation.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from tensorcad import Mesh, MeshBuilder  # noqa: E402
from tensorcad import gz  # noqa: E402
from tensorcad.geometry_field import SphereField, to_mesh  # noqa: E402

LENGTH, BEAM, DRAFT, DEPTH = 10.0, 2.0, 0.8, 1.2
DEEP_DEPTH = 2.4
DEEP_KG = 0.3
LOADED_KG = 0.78
UNSTABLE_KG = 0.9
SPHERE_KG = 0.5
ANGLES = tuple(float(angle) for angle in range(0, 91, 5))


def barge(length=LENGTH, beam=BEAM, draft=DRAFT, depth=DEPTH) -> Mesh:
    """Closed box barge: waterline at z = 0, keel at -draft, deck at depth - draft."""
    builder = MeshBuilder()
    corner = {}
    for x in (0.0, length):
        for y in (-0.5 * beam, 0.5 * beam):
            for z in (-draft, -draft + depth):
                corner[(x, y, z)] = builder.add_vertex((x, y, z))

    def quad(p0, p1, p2, p3) -> None:
        a, b, c, d = (corner[point] for point in (p0, p1, p2, p3))
        builder.add_triangle(a, b, c)
        builder.add_triangle(a, c, d)

    lo, hi = -0.5 * beam, 0.5 * beam
    bottom, top = -draft, -draft + depth
    quad((0, lo, bottom), (0, hi, bottom), (length, hi, bottom), (length, lo, bottom))
    quad((0, lo, top), (length, lo, top), (length, hi, top), (0, hi, top))
    quad((0, lo, bottom), (length, lo, bottom), (length, lo, top), (0, lo, top))
    quad((0, hi, bottom), (0, hi, top), (length, hi, top), (length, hi, bottom))
    quad((0, lo, bottom), (0, lo, top), (0, hi, top), (0, hi, bottom))
    quad((length, lo, bottom), (length, hi, bottom), (length, hi, top), (length, lo, top))
    return builder.build()


def sphere(resolution=(20, 20, 20), radius=1.0) -> Mesh:
    return to_mesh(SphereField(radius=radius), resolution).mesh


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "docs" / "gz",
    )
    arguments = parser.parse_args()
    output = arguments.output
    output.mkdir(parents=True, exist_ok=True)

    claims: list[dict] = []

    def claim(name: str, measured: float, limit: str, passed: bool, unit: str = "") -> None:
        claims.append(
            {"claim": name, "measured": measured, "limit": limit, "passed": bool(passed),
             "unit": unit}
        )

    # ---- the loaded barge: peak, vanishing angle, range of stability ------------
    loaded = barge()
    loaded_curve = gz.gz_curve(loaded, ANGLES, LOADED_KG)
    deck_edge = math.degrees(gz.barge_deck_edge_angle(BEAM, DRAFT, DEPTH))

    # ---- the same barge against the independent two-dimensional reference -------
    differences = [
        abs(arm - gz.barge_righting_arm(BEAM, DRAFT, DEPTH, LOADED_KG, math.radians(angle)))
        for angle, arm in zip(loaded_curve.angles_deg, loaded_curve.arms)
    ]
    worst_section = max(differences)
    claim("barge: mesh arm vs independent 2-D section integration", worst_section,
          "< 1e-12 m", worst_section < 1e-12, "m")

    # ---- the deep barge against the wall-sided closed form ----------------------
    deep = barge(depth=DEEP_DEPTH)
    deep_curve = gz.gz_curve(deep, ANGLES, DEEP_KG)
    limit = gz.wall_sided_range(BEAM, DRAFT)
    wall_sided_rows = []
    worst_wall_sided = 0.0
    for angle in (5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0):
        radians = math.radians(angle)
        measured = gz.righting_arm(deep, radians, DEEP_KG).righting_arm
        closed = gz.wall_sided_righting_arm(BEAM, DRAFT, DEEP_KG, radians)
        error = abs(measured - closed) / max(abs(closed), 1e-30)
        worst_wall_sided = max(worst_wall_sided, error)
        wall_sided_rows.append({"heel_deg": angle, "mesh_m": measured, "closed_form_m": closed,
                                "relative_error": error})
    assert limit > math.radians(35.0)
    claim("deep barge: mesh arm vs wall-sided closed form", worst_wall_sided,
          "< 1e-12 relative", worst_wall_sided < 1e-12)

    # ---- displacement balance ---------------------------------------------------
    worst_balance = max(
        abs(point.relative_volume_error)
        for curve in (loaded_curve, deep_curve)
        for point in [
            gz.righting_arm(loaded if curve is loaded_curve else deep,
                            math.radians(angle), curve.kg)
            for angle in ANGLES
        ]
    )
    claim("heeled displacement balance (all angles, both barges)", worst_balance,
          "< 1e-13 relative", worst_balance < 1e-13)

    # ---- GM from the curve, as the heel shrinks ---------------------------------
    exact_gm = gz.wall_sided_metacentric_height(BEAM, DRAFT, DEEP_KG)
    bm = BEAM * BEAM / (12.0 * DRAFT)
    gm_ladder = []
    for angle in (4.0, 2.0, 1.0, 0.5, 0.25):
        arm = gz.righting_arm(deep, math.radians(angle), DEEP_KG).righting_arm
        estimate = arm / math.sin(math.radians(angle))
        gm_ladder.append({"heel_deg": angle, "gm_from_curve_m": estimate,
                          "error_m": estimate - exact_gm})
    # GZ/sin(phi) - GM ~ (BM/2) tan^2(phi) is the wall-sided second-order term, not a
    # quadrature error: remove it and the slope at 0.25 deg *is* GM.
    last = gm_ladder[-1]
    tangent = math.tan(math.radians(last["heel_deg"]))
    corrected = last["gm_from_curve_m"] - 0.5 * bm * tangent * tangent
    gm_error = abs(corrected - exact_gm)
    claim("GM from the curve at 0.25 deg, minus the BM/2*tan^2 term", gm_error, "< 1e-12 m",
          gm_error < 1e-12, "m")
    gm_order = [
        abs(gm_ladder[index]["error_m"]) / abs(gm_ladder[index + 1]["error_m"])
        for index in range(len(gm_ladder) - 1)
    ]
    claim("halving the heel quarters the GM deficit (second order)", min(gm_order), "> 3.9",
          min(gm_order) > 3.9)

    # ---- angle of loll against the closed form ----------------------------------
    unstable = gz.gz_curve(deep, ANGLES, UNSTABLE_KG)
    gm_unstable = gz.wall_sided_metacentric_height(BEAM, DRAFT, UNSTABLE_KG)
    expected_loll = math.degrees(math.atan(math.sqrt(-2.0 * gm_unstable / bm)))
    loll_error = abs((unstable.loll_angle_deg or float("nan")) - expected_loll)
    claim("angle of loll vs tan^2(phi) = -2 GM / BM", loll_error, "< 1e-6 deg",
          loll_error < 1e-6, "deg")

    # ---- the vanishing angle is a zero of the curve -----------------------------
    vanishing = loaded_curve.vanishing_angle_deg
    assert vanishing is not None
    vanishing_arm = gz.righting_arm(loaded, math.radians(vanishing), LOADED_KG).righting_arm
    claim("GZ at the bisected angle of vanishing stability", abs(vanishing_arm), "< 1e-12 m",
          abs(vanishing_arm) < 1e-12, "m")

    # ---- dynamical stability against the exact wall-sided area ------------------
    exact_area = gz.wall_sided_dynamical_stability(BEAM, DRAFT, DEEP_KG, math.radians(30.0))
    area_curve = gz.gz_curve(deep, [float(angle) for angle in range(0, 31, 5)], DEEP_KG)
    area_error = abs(area_curve.dynamical_stability_m_rad - exact_area) / exact_area
    claim("dynamical stability 0-30 deg vs closed-form integral", area_error,
          "< 1e-3 relative", area_error < 1e-3)
    loaded_area = loaded_curve.dynamical_stability_m_rad

    # ---- the sphere: pure discretisation error ---------------------------------
    sphere_rows = []
    for resolution in ((12, 12, 12), (16, 16, 16), (20, 20, 20)):
        mesh = sphere(resolution)
        worst = 0.0
        balance = 0.0
        for angle in (15.0, 30.0, 45.0, 60.0, 75.0, 90.0):
            radians = math.radians(angle)
            point = gz.righting_arm(mesh, radians, SPHERE_KG)
            exact = (1.0 - SPHERE_KG) * math.sin(radians)
            worst = max(worst, abs(point.righting_arm - exact) / exact)
            balance = max(balance, abs(point.relative_volume_error))
        sphere_rows.append({"resolution": list(resolution), "triangles": len(mesh.faces),
                            "worst_relative_error": worst,
                            "worst_volume_error": balance})
    claim("sphere (20^3) worst GZ error vs (KM - KG) sin(phi)",
          sphere_rows[-1]["worst_relative_error"], "< 5e-3",
          sphere_rows[-1]["worst_relative_error"] < 5e-3)
    ratio = sphere_rows[0]["worst_relative_error"] / sphere_rows[-1]["worst_relative_error"]
    claim("sphere: error ratio between 12^3 and 20^3", ratio, "> 2", ratio > 2.0)

    # ---- linearity in the loading ----------------------------------------------
    angle = math.radians(40.0)
    base = gz.righting_arm(loaded, angle, LOADED_KG).righting_arm
    shifted = gz.righting_arm(loaded, angle, LOADED_KG + 0.2).righting_arm
    linearity = abs((shifted - base) + 0.2 * math.sin(angle))
    claim("dGZ/dKG = -sin(phi)", linearity, "< 1e-13 m", linearity < 1e-13, "m")

    # ---- the trap: a mesh capped at the waterline ------------------------------
    refused = False
    try:
        gz.righting_arm(barge(depth=DRAFT), math.radians(5.0), DEEP_KG)
    except gz.HeelError:
        refused = True
    claim("a mesh capped at the design waterline is refused", float(refused), "== 1",
          refused)

    # ------------------------------------------------------------------ figures
    figure, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    right = axes[0]
    right.plot(loaded_curve.angles_deg, loaded_curve.arms, "o-", label=f"barge, KG={LOADED_KG} m")
    right.plot(deep_curve.angles_deg, deep_curve.arms, "s-", label=f"deep barge, KG={DEEP_KG} m")
    right.plot(unstable.angles_deg, unstable.arms, "^-", label=f"deep barge, KG={UNSTABLE_KG} m")
    reference_angles = [angle for angle in ANGLES if angle <= math.degrees(limit)]
    right.plot(
        reference_angles,
        [gz.wall_sided_righting_arm(BEAM, DRAFT, DEEP_KG, math.radians(angle))
         for angle in reference_angles],
        "k--", lw=1.0, label="wall-sided closed form (deep barge)",
    )
    right.axhline(0.0, color="0.6", lw=0.8)
    right.axvline(deck_edge, color="0.6", lw=0.8, ls=":")
    right.annotate(f"deck edge {deck_edge:.1f} deg", xy=(deck_edge, 0.005),
                   xytext=(deck_edge + 3.0, 0.10), fontsize=8,
                   arrowprops={"arrowstyle": "->", "lw": 0.6})
    if unstable.loll_angle_deg:
        right.plot([unstable.loll_angle_deg], [0.0], "r*", ms=9, label="angle of loll (closed form)")
    right.set_xlabel("heel angle [deg]")
    right.set_ylabel("GZ [m]")
    right.set_title("righting arms: mesh vs closed form")
    right.grid(alpha=0.3)
    right.legend(fontsize=7, loc="upper left")

    left = axes[1]
    resolutions = [row["triangles"] for row in sphere_rows]
    left.loglog(resolutions, [row["worst_relative_error"] for row in sphere_rows], "o-",
                label="sphere: worst GZ error vs analytic")
    left.loglog([6, 8, 10], [1.0 / (value * value) for value in (6, 8, 10)], "k--", lw=0.8,
                label="second order guide")
    left.set_xlabel("triangles")
    left.set_ylabel("relative error")
    left.set_title("sphere: discretisation error, and the GM ladder")
    left.grid(alpha=0.3, which="both")
    left.legend(fontsize=7, loc="lower left")
    inset = left.inset_axes([0.52, 0.52, 0.44, 0.4])
    inset.semilogy([row["heel_deg"] for row in gm_ladder],
                   [abs(row["error_m"]) for row in gm_ladder], "s-", ms=3)
    inset.set_title("|GZ/sin(phi) - GM| [m]", fontsize=6)
    inset.tick_params(labelsize=6)
    inset.grid(alpha=0.3)
    figure.tight_layout()
    figure.savefig(output / "gz_curves.png", dpi=150)
    plt.close(figure)

    # ------------------------------------------------------------------ results
    payload = {
        "geometry": {
            "barge": {"length_m": LENGTH, "beam_m": BEAM, "draft_m": DRAFT,
                      "depth_m": DEPTH, "deep_depth_m": DEEP_DEPTH,
                      "deck_edge_deg": deck_edge,
                      "wall_sided_limit_deg": math.degrees(limit)},
            "sphere": {"radius_m": 1.0, "kg_m": SPHERE_KG},
        },
        "loading": {"loaded_kg_m": LOADED_KG, "deep_kg_m": DEEP_KG,
                    "unstable_kg_m": UNSTABLE_KG},
        "hydrostatics": {
            "displaced_volume_m3": LENGTH * BEAM * DRAFT,
            "km_m": gz.wall_sided_metacentric_height(BEAM, DRAFT, 0.0),
            "gm_loaded_m": gz.wall_sided_metacentric_height(BEAM, DRAFT, LOADED_KG),
            "gm_deep_m": exact_gm,
            "gm_unstable_m": gm_unstable,
        },
        "loaded_curve": loaded_curve.to_dict(),
        "deep_curve": deep_curve.to_dict(),
        "unstable_curve": unstable.to_dict(),
        "curve_table": [
            {"heel_deg": angle,
             "mesh_m": arm,
             "section_reference_m": gz.barge_righting_arm(BEAM, DRAFT, DEPTH, LOADED_KG,
                                                          math.radians(angle)),
             "difference_m": difference}
            for angle, arm, difference in zip(loaded_curve.angles_deg, loaded_curve.arms,
                                              differences)
        ],
        "wall_sided_table": wall_sided_rows,
        "gm_ladder": gm_ladder,
        "sphere_ladder": sphere_rows,
        "angles_of_interest": {
            "deck_edge_deg": deck_edge,
            "wall_sided_limit_deg": math.degrees(limit),
            "angle_of_maximum_deg": loaded_curve.angle_of_maximum_deg,
            "maximum_arm_m": loaded_curve.maximum_arm,
            "vanishing_angle_deg": vanishing,
            "dynamical_stability_m_rad": loaded_area,
            "loll_angle_deg": unstable.loll_angle_deg,
            "expected_loll_deg": expected_loll,
        },
        "claims": claims,
        "all_claims_hold": all(item["passed"] for item in claims),
    }
    (output / "results.json").write_text(json.dumps(payload, indent=2) + "\n")

    # ------------------------------------------------------------------ report
    lines = [
        "# Large-angle stability: GZ from the heeled displacement",
        "",
        f"Barge {LENGTH} x {BEAM} x {DEPTH} m, draft {DRAFT} m (`KG = {LOADED_KG} m`); "
        f"deep barge depth {DEEP_DEPTH} m (`KG = {DEEP_KG}` / `{UNSTABLE_KG} m`); "
        f"sphere `R = 1 m`, `KG = {SPHERE_KG} m`.",
        "",
        "## Numbers",
        "",
        "| item | value |",
        "|---|---:|",
    ]
    for key, value in payload["angles_of_interest"].items():
        lines.append(f"| {key} | {value} |")
    lines += ["", "## Claims", "", "| claim | measured | limit | holds |", "|---|---:|---|:--:|"]
    for item in claims:
        lines.append(
            f"| {item['claim']} | {item['measured']:.4e} {item['unit']} | {item['limit']} | "
            f"{'yes' if item['passed'] else 'NO'} |"
        )
    lines += [
        "",
        "## Curve, sampled",
        "",
        "| heel [deg] | mesh GZ [m] | section reference [m] | difference [m] |",
        "|---:|---:|---:|---:|",
    ]
    for row in payload["curve_table"]:
        lines.append(
            f"| {row['heel_deg']:.0f} | {row['mesh_m']:.9f} | {row['section_reference_m']:.9f} | "
            f"{row['difference_m']:.2e} |"
        )
    (output / "report.md").write_text("\n".join(lines) + "\n")

    for item in claims:
        state = "ok  " if item["passed"] else "FAIL"
        print(f"{state} {item['claim']}: {item['measured']:.4e} {item['unit']} ({item['limit']})")
    print(f"\nGZ max {loaded_curve.maximum_arm:.6f} m at {loaded_curve.angle_of_maximum_deg} deg, "
          f"vanishing {vanishing:.3f} deg, dynamical stability {loaded_area:.6f} m rad")
    print(f"artefacts in {output}")
    return 0 if payload["all_claims_hold"] else 1


if __name__ == "__main__":
    sys.exit(main())