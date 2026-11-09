#!/usr/bin/env python3
"""Verify the heeled GZ chain on the transom-stern hull family.

Produces ``artifacts/transom_gz/evidence.json`` (every number in the doc and the
report is read from this file, none retyped) and three figures.

Run:  python3 tools/verify_transom_gz.py
"""

from __future__ import annotations

import json
import math
from pathlib import Path

from tensorcad.gz import (
    deck_edge_immersion_deg,
    heeled_hydrostatics,
    righting_arm_curve as gz_curve,
    with_freeboard,
)
from tensorcad.righting_arm import righting_arm_curve as bisection_curve
from tensorcad.stability import initial_stability, waterplane_moments
from tensorcad.transom import generate_transom_hull

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "artifacts" / "transom_gz"
FREEBOARD = 0.15
KG = 0.3
HEELS = (0.0, 2.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0)


def main() -> dict:
    OUT.mkdir(parents=True, exist_ok=True)
    model = generate_transom_hull()
    mesh = model.mesh
    banded = with_freeboard(mesh, FREEBOARD, waterline_z=0.0)
    keel_z = mesh.bounding_box().minimum[2]

    # --- the floating condition on the *unbanded* hull: GM from waterplane moments
    state = initial_stability(mesh, kg=KG, waterline_z=0.0)
    section = waterplane_moments(mesh, 0.0)
    upright_volume = model.hydrostatics.displacement_volume

    # --- two independent heeled-displacement integrators
    fast = gz_curve(banded, kg=KG, angles_deg=HEELS, waterline_z=0.0)
    slow = bisection_curve(banded, kg=KG, heels_deg=HEELS, waterline_z=0.0)
    fast_gz = list(fast.gz)
    slow_gz = [point.righting_arm for point in slow.points]
    cross = [abs(a - b) for a, b in zip(fast_gz, slow_gz)]

    rows = []
    for index, heel in enumerate(HEELS):
        point = heeled_hydrostatics(banded, heel, waterline_z=0.0)
        rows.append(
            {
                "heel_deg": heel,
                "gz_fast": fast_gz[index],
                "gz_bisection": slow_gz[index],
                "cross_route_abs_diff": cross[index],
                "kn": point.kn,
                "displacement_volume": point.displacement_volume,
                "volume_residual": (point.displacement_volume - upright_volume)
                / upright_volume,
                "waterplane_area": point.waterplane_area,
                "wetted_surface_area": point.wetted_surface_area,
                "waterline_rise": point.waterline_offset,
                "kb_z_from_keel": point.centre_of_buoyancy_z,
                "centre_of_buoyancy_y": point.centre_of_buoyancy_y,
            }
        )

    trailing = [
        {"heel_deg": p.heel_deg, "displacement_volume": p.displacement_volume,
         "waterline_z": p.waterline_z, "volume_residual": p.volume_residual}
        for p in slow.points
    ]

    small = 2.0
    slope = fast_gz[HEELS.index(small)] / math.sin(math.radians(small))
    gm = state.gm_transverse
    # deck-edge immersion: the heel at which the *low* deck edge goes under
    half_beam = model.hydrostatics.waterplane_area  # placeholder, replaced below
    half_breadth = max(v[1] for v in mesh.vertices)
    deck_edge = deck_edge_immersion_deg(half_breadth, FREEBOARD)

    evidence = {
        "hull": {
            "family": "transom-stern (tensorcad.transom), default parameters",
            "length": model.parameters.length,
            "beam": model.parameters.beam,
            "draft": model.parameters.draft,
            "transom_fraction": model.parameters.transom_fraction,
            "nx": model.parameters.nx,
            "nz": model.parameters.nz,
            "keel_z": keel_z,
            "vertices": len(mesh.vertices),
            "faces": len(mesh.faces),
        },
        "freeboard": FREEBOARD,
        "kg": KG,
        "upright": {
            "displacement_volume_mesh": upright_volume,
            "displacement_volume_analytic": model.analytic["displacement_volume"],
            "displacement_gap_percent": 100.0
            * (upright_volume / model.analytic["displacement_volume"] - 1.0),
            "kb": state.kb,
            "bm_transverse": state.bm_transverse,
            "gm_transverse": gm,
            "waterplane_area_mesh": section.area,
            "transverse_second_moment": section.transverse_second_moment,
            "upright_gz_fast": fast_gz[0],
            "upright_gz_bisection": slow_gz[0],
        },
        "deck_edge_immersion_deg": deck_edge,
        "small_angle": {
            "heel_deg": small,
            "gz_over_sin": slope,
            "gm": gm,
            "slope_minus_gm": slope - gm,
            "slope_gap_percent": 100.0 * (slope / gm - 1.0),
        },
        "heel_table": rows,
        "bisection_volume_residual_max": max(
            abs(item["volume_residual"]) for item in trailing
        ),
        "cross_route_max_abs_diff": max(cross),
        "gz_max": max(fast_gz), "gz_max_heel_deg": HEELS[fast_gz.index(max(fast_gz))],
        "surface_area_ratio_40_over_0": rows[-1]["wetted_surface_area"]
        / rows[0]["wetted_surface_area"],
    }
    (OUT / "evidence.json").write_text(json.dumps(evidence, indent=2) + "\n")

    figures(rows, evidence, OUT)
    print(json.dumps({k: v for k, v in evidence.items() if k != "heel_table"}, indent=2))
    for row in rows:
        print(
            f"{row['heel_deg']:5.1f} GZ={row['gz_fast']:.6f} bis={row['gz_bisection']:.6f} "
            f"d={row['cross_route_abs_diff']:.2e} Vres={row['volume_residual']:+.2e} "
            f"Aw={row['waterplane_area']:.4f} rise={row['waterline_rise']:+.5f}"
        )
    return evidence


def figures(rows, evidence, out: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    heels = [row["heel_deg"] for row in rows]
    plt.figure(figsize=(6.4, 4.0))
    plt.plot(heels, [row["gz_fast"] for row in rows], "o-", label="plane-sweep integrator (gz)")
    plt.plot(heels, [row["gz_bisection"] for row in rows], "x--",
             label="bisection integrator (righting_arm)")
    gm = evidence["upright"]["gm_transverse"]
    plt.plot(heels, [gm * math.sin(math.radians(h)) for h in heels], ":",
             label=r"$GM\sin\varphi$ (small-angle limit)")
    plt.axvline(evidence["deck_edge_immersion_deg"], color="grey", lw=1,
                label="deck-edge immersion")
    plt.xlabel("heel [deg]")
    plt.ylabel("GZ [m]")
    plt.title("Transom-stern hull: righting arm, two independent integrators")
    plt.grid(alpha=0.3)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(out / "fig1_gz_routes.png", dpi=140)
    plt.close()

    plt.figure(figsize=(6.4, 4.0))
    plt.semilogy(heels[1:], [row["cross_route_abs_diff"] for row in rows[1:]], "o-",
                 label="|GZ_sweep - GZ_bisection|")
    plt.semilogy(heels[1:], [abs(row["volume_residual"]) for row in rows[1:]], "s--",
                 label="|V(phi) - V(0)| / V(0)")
    plt.xlabel("heel [deg]")
    plt.ylabel("relative / absolute residual")
    plt.title("Consistency of the heeled equilibrium")
    plt.grid(alpha=0.3, which="both")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(out / "fig2_consistency.png", dpi=140)
    plt.close()

    plt.figure(figsize=(6.4, 4.0))
    plt.plot(heels, [row["waterline_rise"] for row in rows], "o-")
    plt.xlabel("heel [deg]")
    plt.ylabel("waterline rise of the floating attitude [m]")
    plt.title("Waterline solved at constant displacement (immersed transom)")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(out / "fig3_waterline.png", dpi=140)
    plt.close()


if __name__ == "__main__":
    main()