#!/usr/bin/env python3
"""Verify the righting-arm integration against the one closed form that exists for it.

A ``GZ`` curve has no general analytic solution, which is exactly why a curve that
looks smooth is not evidence of anything. This example therefore checks it against
two references that were computed by different code:

1. **A rectangular barge against the wall-sided closed form**
   ``GZ = sin(phi)(GM + BM tan^2(phi)/2)``, which is *exact* while the waterplane
   stays on the vertical sides. The barge is a polygon and the module clips
   polygons exactly, so the agreement has to be round-off; anything else is a
   defect, not a grid effect.
2. **The same barge's initial slope against ``GM = KB + BM - KG``** rebuilt from
   the geometry (``KB = T/2``, ``BM = B^2/(12T)``), by Richardson extrapolation of
   ``GZ/phi`` from 1 and 0.5 degrees.
3. **A hull's initial slope against ``KB + I_T/V - KG``**, where ``I_T`` comes from
   :func:`tensorcad.stability.waterplane_moments` -- the surface-cut integration,
   not the station areas this module uses.
4. **The station-area volume against the closed-surface tetrahedron volume** of the
   same mesh (:func:`tensorcad.hydrostatics.compute_hydrostatics`), which must be
   the same number up to round-off.
5. **Where the model stops describing a ship**: the deck-edge immersion angle,
   compared with the classic ``(B/2) tan(phi) = freeboard`` limit, and the
   convergence of the curve with the station grid.

Writes ``results.json``, ``report.md`` and ``righting_arms.png`` into ``--output``
and exits non-zero if any claim stops holding.

Usage::

    python examples/righting/run_example.py [--output docs/righting]
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:  # allow running from a checkout without installing
    sys.path.insert(0, str(REPO_ROOT))

import tensorcad as tc
from tensorcad.righting import (
    righting_arm_curve,
    submerged_body,
    wall_sided_righting_arm,
)

BARGE = {"length": 10.0, "beam": 2.0, "draft": 1.0, "freeboard": 1.0, "kg": 0.8}
HULL = {"length": 10.0, "beam": 2.0, "draft": 1.0, "nx": 41, "nz": 13}
HULL_WATERLINE = -0.3   # 0.3 m of freeboard above the floating waterline
HULL_KG = 0.75

#: Angles the barge is still wall-sided at: a 1.0 m freeboard over a 1.0 m half
#: breadth means the deck edge goes under at 45 degrees.
BARGE_ANGLES = tuple(float(a) for a in range(0, 46, 5))
CURVE_ANGLES = tuple(float(a) for a in range(0, 91, 5))


def barge_mesh(length: float, beam: float, draft: float, freeboard: float) -> tc.Mesh:
    """Watertight rectangular barge: ``draft`` below ``z = 0``, deck at ``z = freeboard``."""
    half = beam / 2.0
    vertices = (
        (0.0, -half, -draft),
        (length, -half, -draft),
        (length, half, -draft),
        (0.0, half, -draft),
        (0.0, -half, freeboard),
        (length, -half, freeboard),
        (length, half, freeboard),
        (0.0, half, freeboard),
    )
    faces = (
        (0, 2, 1), (0, 3, 2),
        (4, 5, 6), (4, 6, 7),
        (0, 1, 5), (0, 5, 4),
        (1, 2, 6), (1, 6, 5),
        (2, 3, 7), (2, 7, 6),
        (3, 0, 4), (3, 4, 7),
    )
    return tc.Mesh(vertices=vertices, faces=faces)


def richardson_slope(curve: tc.RightingArmCurve) -> float:
    """``dGZ/dphi`` at zero heel from the 0.5 and 1 degree points."""
    one = curve.points[[p.heel_deg for p in curve.points].index(1.0)]
    half = curve.points[[p.heel_deg for p in curve.points].index(0.5)]
    coarse = one.gz / math.radians(1.0)
    fine = half.gz / math.radians(0.5)
    return (4.0 * fine - coarse) / 3.0


def build_results() -> dict:
    results: dict = {"barge": {}, "hull": {}, "claims": []}

    # -------------------------------------------------------------- barge: exact
    barge_geometry = barge_mesh(BARGE["length"], BARGE["beam"], BARGE["draft"], BARGE["freeboard"])
    barge = results["barge"]
    kb = BARGE["draft"] / 2.0
    bm = BARGE["beam"] ** 2 / (12.0 * BARGE["draft"])
    gm = kb + bm - BARGE["kg"]
    barge_curve = righting_arm_curve(
        barge_geometry, kg=BARGE["kg"], heel_angles_deg=BARGE_ANGLES, station_intervals=8
    )
    rows = []
    worst = 0.0
    for point in barge_curve.points:
        reference = wall_sided_righting_arm(point.heel_deg, gm=gm, bm=bm)
        deviation = abs(point.gz - reference) / reference if reference else 0.0
        worst = max(worst, deviation)
        rows.append(
            {
                "heel_deg": point.heel_deg,
                "gz_m": point.gz,
                "gz_wall_sided_m": reference,
                "relative_deviation": deviation,
                "submerged_volume_m3": point.submerged_volume,
            }
        )
    slope_curve = righting_arm_curve(
        barge_geometry, kg=BARGE["kg"], heel_angles_deg=(0.0, 0.5, 1.0), station_intervals=8
    )
    barge.update(
        {
            "parameters": BARGE,
            "kb_m": kb,
            "bm_m": bm,
            "gm_m": gm,
            "volume_m3": BARGE["length"] * BARGE["beam"] * BARGE["draft"],
            "rows": rows,
            "worst_relative_deviation": worst,
            "slope_from_curve": richardson_slope(slope_curve),
            "slope_relative_error": abs(richardson_slope(slope_curve) - gm) / gm,
            "deck_immersion_onset_deg": barge_curve.deck_immersion_onset_deg,
            "classic_deck_edge_deg": math.degrees(math.atan2(BARGE["freeboard"], BARGE["beam"] / 2.0)),
            "full_curve": [point.to_dict() for point in barge_curve.points],
            "gz_max_m": barge_curve.gz_max,
            "angle_gz_max_deg": barge_curve.angle_gz_max_deg,
        }
    )

    # --------------------------------------------------------------- hull: checks
    hull_geometry = tc.generate_hull(tc.HullParameters(**HULL)).mesh
    hull = results["hull"]
    volume_check = submerged_body(hull_geometry, heel_deg=0.0, waterline_z=0.0, station_intervals=80)
    hydrostatics = tc.compute_hydrostatics(hull_geometry, 0.0)
    hull_curve = righting_arm_curve(
        hull_geometry, kg=HULL_KG, waterline_z=HULL_WATERLINE, heel_angles_deg=CURVE_ANGLES, station_intervals=80
    )
    slope = righting_arm_curve(
        hull_geometry, kg=HULL_KG, waterline_z=HULL_WATERLINE, heel_angles_deg=(0.0, 0.5, 1.0), station_intervals=80
    )
    curve_slope = richardson_slope(slope)
    kb_hull = hull_curve.points[0].vertical_centre_buoyancy - hull_curve.keel_z
    bm_hull = tc.waterplane_moments(hull_geometry, HULL_WATERLINE).transverse_second_moment / hull_curve.displacement_volume
    gm_hull = kb_hull + bm_hull - HULL_KG

    coarse = righting_arm_curve(
        hull_geometry, kg=HULL_KG, waterline_z=HULL_WATERLINE, heel_angles_deg=CURVE_ANGLES, station_intervals=40
    )
    finer = righting_arm_curve(
        hull_geometry, kg=HULL_KG, waterline_z=HULL_WATERLINE, heel_angles_deg=CURVE_ANGLES, station_intervals=160
    )

    zero_freeboard_rejected = False
    try:
        righting_arm_curve(hull_geometry, kg=HULL_KG, waterline_z=0.0, heel_angles_deg=(0.0, 5.0))
    except tc.RightingError:
        zero_freeboard_rejected = True

    hull.update(
        {
            "parameters": {**HULL, "kg": HULL_KG, "waterline_z": HULL_WATERLINE},
            "freeboard_m": hull_curve.freeboard,
            "displacement_volume_m3": hull_curve.displacement_volume,
            "surface_integral_volume_m3": hydrostatics.displacement_volume,
            "volume_relative_difference": abs(
                volume_check.volume - hydrostatics.displacement_volume
            )
            / hydrostatics.displacement_volume,
            "vcb_relative_difference": abs(
                volume_check.vertical_centre_buoyancy - hydrostatics.vertical_centre_buoyancy
            )
            / abs(hydrostatics.vertical_centre_buoyancy),
            "kb_m": kb_hull,
            "bm_m": bm_hull,
            "gm_m": gm_hull,
            "slope_from_curve": curve_slope,
            "slope_relative_error": abs(curve_slope - gm_hull) / gm_hull,
            "gz_max_m": hull_curve.gz_max,
            "angle_gz_max_deg": hull_curve.angle_gz_max_deg,
            "angle_vanishing_stability_deg": hull_curve.angle_vanishing_stability_deg,
            "deck_immersion_onset_deg": hull_curve.deck_immersion_onset_deg,
            "zero_freeboard_rejected": zero_freeboard_rejected,
            "grid_convergence": {
                "intervals": [40, 80, 160],
                "volume_m3": [
                    coarse.displacement_volume,
                    hull_curve.displacement_volume,
                    finer.displacement_volume,
                ],
                "gz_max_m": [coarse.gz_max, hull_curve.gz_max, finer.gz_max],
                "angle_vanishing_deg": [
                    coarse.angle_vanishing_stability_deg,
                    hull_curve.angle_vanishing_stability_deg,
                    finer.angle_vanishing_stability_deg,
                ],
            },
            "curve": [point.to_dict() for point in hull_curve.points],
        }
    )

    # ------------------------------------------------------------------- claims
    claims = []

    def claim(name: str, measured: float, limit: float, unit: str, detail: str) -> None:
        claims.append(
            {
                "claim": name,
                "measured": measured,
                "limit": limit,
                "unit": unit,
                "passed": bool(measured <= limit),
                "detail": detail,
            }
        )

    claim(
        "barge GZ matches the wall-sided closed form (0-45 deg)",
        worst,
        1e-9,
        "relative",
        "exact for a polygonal barge below deck-edge immersion",
    )
    claim(
        "barge dGZ/dphi at 0 recovers KB + BM - KG",
        barge["slope_relative_error"],
        1e-6,
        "relative",
        "Richardson extrapolation from 1 and 0.5 deg",
    )
    claim(
        "barge deck-edge immersion occurs at the classic angle",
        abs(barge["deck_immersion_onset_deg"] - barge["classic_deck_edge_deg"]),
        1e-9,
        "deg",
        "(B/2) tan(phi) = freeboard",
    )
    claim(
        "station-area volume equals the closed-surface volume",
        hull["volume_relative_difference"],
        1e-9,
        "relative",
        "Green-on-sections vs tetrahedra on the same mesh",
    )
    claim(
        "hull dGZ/dphi at 0 recovers KB + I_T/V - KG",
        hull["slope_relative_error"],
        2e-3,
        "relative",
        "I_T from tensorcad.stability.waterplane_moments (independent integration)",
    )
    claim(
        "hull curve is grid-converged (40 vs 80 intervals)",
        abs(coarse.gz_max - hull_curve.gz_max) / hull_curve.gz_max,
        5e-3,
        "relative",
        "GZ_max over the 0-90 deg curve",
    )
    claim(
        "hull curve vanishing angle is grid-converged (40 vs 80 intervals)",
        abs(coarse.angle_vanishing_stability_deg - hull_curve.angle_vanishing_stability_deg)
        / hull_curve.angle_vanishing_stability_deg,
        1e-2,
        "relative",
        "interpolated zero crossing after GZ_max",
    )
    claim(
        "zero-freeboard hull is rejected rather than answered",
        0.0 if zero_freeboard_rejected else 1.0,
        0.0,
        "bool",
        "a body capped at its waterline has no freeboard to heel into",
    )
    results["claims"] = claims
    return results


def render(results: dict, output: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, (left, right) = plt.subplots(1, 2, figsize=(12.0, 4.8))

    barge = results["barge"]
    angles = [row["heel_deg"] for row in barge["rows"]]
    left.plot(angles, [row["gz_m"] for row in barge["rows"]], "o-", label="integration")
    left.plot(
        angles,
        [row["gz_wall_sided_m"] for row in barge["rows"]],
        "k--",
        label="wall-sided closed form",
    )
    left.axvline(barge["deck_immersion_onset_deg"], color="tab:red", ls=":", lw=1.0)
    left.annotate(
        "deck edge immersed\n(formula no longer valid)",
        xy=(barge["deck_immersion_onset_deg"], barge["rows"][-1]["gz_m"] * 0.55),
        fontsize=8,
        color="tab:red",
    )
    left.set_xlabel("heel angle [deg]")
    left.set_ylabel("GZ [m]")
    left.set_title(
        "Barge L=10 B=2 T=1 m, KG=0.8 m\n"
        f"max deviation {barge['worst_relative_deviation']:.1e} (round-off)"
    )
    left.legend(fontsize=8)
    left.grid(alpha=0.3)

    hull = results["hull"]
    heeled = [point["heel_deg"] for point in hull["curve"] if point["gz_m"] is not None]
    arms = [point["gz_m"] for point in hull["curve"] if point["gz_m"] is not None]
    right.plot(heeled, arms, "o-", color="tab:blue", label="GZ (station integration)")
    onset = hull["deck_immersion_onset_deg"]
    if onset is not None and not (isinstance(onset, float) and math.isnan(onset)):
        right.axvspan(onset, 90.0, color="tab:red", alpha=0.12)
        right.annotate(
            "deck edge immersed:\nwatertight-deck extrapolation",
            xy=((onset + 90.0) / 2.0, 0.02),
            fontsize=8,
            ha="center",
            color="tab:red",
        )
    right.axhline(0.0, color="k", lw=0.8)
    right.axvline(hull["angle_gz_max_deg"], color="tab:green", ls=":", lw=1.0)
    right.annotate(
        f"GZ_max {hull['gz_max_m']:.4f} m @ {hull['angle_gz_max_deg']:.1f} deg",
        xy=(hull["angle_gz_max_deg"], hull["gz_max_m"]),
        xytext=(hull["angle_gz_max_deg"] + 8.0, hull["gz_max_m"] * 0.85),
        fontsize=8,
        arrowprops={"arrowstyle": "->", "lw": 0.8},
    )
    vanishing = hull["angle_vanishing_stability_deg"]
    if vanishing is not None:
        right.axvline(vanishing, color="tab:orange", ls="--", lw=1.0)
        right.annotate(
            f"vanishing angle {vanishing:.1f} deg",
            xy=(vanishing, hull["gz_max_m"] * 0.4),
            xytext=(vanishing + 6.0, hull["gz_max_m"] * 0.45),
            fontsize=8,
            color="tab:orange",
            arrowprops={"arrowstyle": "->", "lw": 0.8, "color": "tab:orange"},
        )
    right.set_xlabel("heel angle [deg]")
    right.set_ylabel("GZ [m]")
    right.set_title(
        "Hull (sin^0.7 sqrt family), 0.3 m freeboard, KG=0.75 m\n"
        f"slope at 0 vs KB+I_T/V-KG: {hull['slope_relative_error']:.1e} relative"
    )
    right.grid(alpha=0.3)
    right.legend(fontsize=8)

    figure.tight_layout()
    figure.savefig(output / "righting_arms.png", dpi=160)
    plt.close(figure)


def render_markdown(results: dict) -> str:
    barge, hull, claims = results["barge"], results["hull"], results["claims"]
    lines: list[str] = []
    add = lines.append
    add("# Righting-arm curves: evidence")
    add("")
    add("Generated by `examples/righting/run_example.py`; every number in this file was")
    add("written by the program, none was typed.")
    add("")
    add("## 1. Case")
    add("")
    add("| body | L [m] | B [m] | T [m] | freeboard [m] | KG [m] |")
    add("|---|--:|--:|--:|--:|--:|")
    add(
        f"| barge | {barge['parameters']['length']} | {barge['parameters']['beam']} | "
        f"{barge['parameters']['draft']} | {barge['parameters']['freeboard']} | "
        f"{barge['parameters']['kg']} |"
    )
    add(
        f"| hull (`sin^0.7 sqrt` family) | {hull['parameters']['length']} | "
        f"{hull['parameters']['beam']} | {hull['parameters']['draft']} | "
        f"{hull['freeboard_m']:.3f} | {hull['parameters']['kg']} |"
    )
    add("")
    add("## 2. Barge against the wall-sided closed form")
    add("")
    add("| heel [deg] | GZ integration [m] | GZ closed form [m] | relative | volume [m^3] |")
    add("|--:|--:|--:|--:|--:|")
    for row in barge["rows"]:
        add(
            f"| {row['heel_deg']:.1f} | {row['gz_m']:.9f} | {row['gz_wall_sided_m']:.9f} | "
            f"{row['relative_deviation']:.2e} | {row['submerged_volume_m3']:.6f} |"
        )
    add("")
    add(f"GM = {barge['gm_m']:.6f} m; slope from the curve = {barge['slope_from_curve']:.9f} m "
        f"(relative error {barge['slope_relative_error']:.2e}).")
    add("")
    add("## 3. Hull curve")
    add("")
    add("| heel [deg] | GZ [m] | y_B [m] | z_B [m] | deck immersion |")
    add("|--:|--:|--:|--:|--:|")
    for point in hull["curve"]:
        add(
            f"| {point['heel_deg']:.1f} | {point['gz_m']:.6f} | "
            f"{point['transverse_centre_buoyancy_m']:.6f} | {point['vertical_centre_buoyancy_m']:.6f} | "
            f"{point['deck_immersion_fraction']:.2f} |"
        )
    add("")
    add("## 4. Convergence with the station grid")
    add("")
    add("| intervals | volume [m^3] | GZ_max [m] | vanishing angle [deg] |")
    add("|--:|--:|--:|--:|")
    convergence = hull["grid_convergence"]
    for index, intervals in enumerate(convergence["intervals"]):
        add(
            f"| {intervals} | {convergence['volume_m3'][index]:.9f} | "
            f"{convergence['gz_max_m'][index]:.6f} | "
            f"{convergence['angle_vanishing_deg'][index]:.3f} |"
        )
    add("")
    add("## 5. Claims")
    add("")
    add("| claim | measured | limit | unit | passed |")
    add("|---|--:|--:|---|:--:|")
    for item in claims:
        add(
            f"| {item['claim']} | {item['measured']:.3e} | {item['limit']:.0e} | "
            f"{item['unit']} | {'yes' if item['passed'] else 'NO'} |"
        )
    add("")
    add("![righting arms](righting_arms.png)")
    add("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("docs/righting"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    results = build_results()
    (args.output / "results.json").write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")
    (args.output / "report.md").write_text(render_markdown(results))
    render(results, args.output)

    failed = [item for item in results["claims"] if not item["passed"]]
    for item in results["claims"]:
        status = "PASS" if item["passed"] else "FAIL"
        print(f"[{status}] {item['claim']}: {item['measured']:.3e} {item['unit']} (limit {item['limit']:.0e})")
    print(f"\nwrote {args.output}/results.json, report.md, righting_arms.png")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
