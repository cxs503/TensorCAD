"""The damaged ship's righting arm, measured against an exact oracle.

Run: ``python examples/damage_gz/run_example.py``

Produces ``outputs/metrics.json`` (every number in ``docs/damage_gz.md`` is read from it,
never typed by hand), ``outputs/table.md`` (the document's tables) and the figures in
``docs/damage_gz/``.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "tests"))
sys.path.insert(0, str(REPO))

from damage_gz_oracle import exact_damaged_arm, wall_sided_arm  # noqa: E402
from test_damage_gz import (  # noqa: E402
    BM0,
    DISPLACEMENT,
    DECK_Z,
    GM0,
    KEEL_Z,
    KG,
    L,
    B,
    T,
    box_mesh,
)

from tensorcad.damage import (  # noqa: E402
    Compartment,
    damage_gm,
    flooded_equilibrium,
    flooded_righting_arm,
    flooded_righting_arm_curve,
)

HEELS = tuple(float(value) for value in range(0, 61, 5))


def main() -> dict:
    hull = box_mesh()
    midship = Compartment(
        4.0,
        6.0,
        y0=-1.0,
        y1=1.0,
        z0=KEEL_Z,
        z1=0.0,
        permeability=0.6,
        name="midship",
    )
    intact = flooded_righting_arm_curve(hull, (), DISPLACEMENT, kg=KG, heels_deg=HEELS)
    damaged = flooded_righting_arm_curve(hull, midship, DISPLACEMENT, kg=KG, heels_deg=HEELS)
    state = flooded_equilibrium(hull, midship, DISPLACEMENT, kg=KG)
    gm = damage_gm(state, KG)

    rows = []
    for point in damaged.points:
        exact = exact_damaged_arm(
            heel_deg=point.heel_deg,
            displacement_volume=DISPLACEMENT,
            length=L,
            beam=B,
            keel_z=KEEL_Z,
            deck_z=DECK_Z,
            kg=KG,
            compartment=(4.0, 6.0, KEEL_Z, 0.0, 0.6),
        )
        wall = wall_sided_arm(point.heel_deg, gm=GM0, bm=BM0)
        rows.append(
            {
                "heel_deg": point.heel_deg,
                "library_gz_m": point.righting_arm,
                "oracle_gz_m": exact.righting_arm,
                "residual_m": point.righting_arm - exact.righting_arm,
                "intact_library_gz_m": intact.points[int(point.heel_deg / 5.0)].righting_arm,
                "wall_sided_gz_m": wall,
                "waterline_z_m": point.waterline.base_z,
                "kb_m": point.kb,
                "tcb_m": point.tcb,
                "lcb_m": point.lcb,
                "volume_m3": point.volume,
            }
        )

    slopes = {}
    for heel in (1e-2, 1e-3, 1e-4):
        point = flooded_righting_arm(hull, midship, DISPLACEMENT, heel_deg=heel, kg=KG)
        slopes[heel] = point.righting_arm / math.sin(math.radians(heel))

    metrics = {
        "case": {
            "length_m": L,
            "beam_m": B,
            "draft_m": T,
            "freeboard_m": DECK_Z,
            "displacement_volume_m3": DISPLACEMENT,
            "kg_m": KG,
            "compartment": midship.to_dict(),
        },
        "intact": {
            "gm_m": GM0,
            "bm_m": BM0,
            "wall_sided_limit_deg": math.degrees(math.atan(2.0 * T / B)),
            "max_arm_m": intact.max_arm,
            "heel_at_max_arm_deg": intact.heel_at_max_arm_deg,
        },
        "damaged": {
            "kb_m": state.kb,
            "gm_m": gm,
            "sinkage_m": state.sinkage,
            "max_arm_m": damaged.max_arm,
            "heel_at_max_arm_deg": damaged.heel_at_max_arm_deg,
            "vanishing_angle_deg": damaged.vanishing_angle_deg,
            "positive_range_deg": damaged.positive_range_deg,
            "area_to_max_arm_m_rad": damaged.area_to_max_arm,
            "small_angle_slopes": {str(k): v for k, v in slopes.items()},
            "small_angle_slope_error": {str(k): v - gm for k, v in slopes.items()},
        },
        "rows": rows,
        "worst_residual_m": max(abs(row["residual_m"]) for row in rows),
        "wall_sided_error_at_60_deg_m": rows[-1]["wall_sided_gz_m"] - rows[-1]["library_gz_m"],
    }

    outputs = HERE / "outputs"
    outputs.mkdir(exist_ok=True)
    (outputs / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")

    table = [
        "| heel [deg] | library GZ [m] | exact oracle GZ [m] | residual [m] | intact GZ [m] | wall-sided GZ [m] | waterline z [m] | KB [m] |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        table.append(
            f"| {row['heel_deg']:.0f} | {row['library_gz_m']:.9f} | {row['oracle_gz_m']:.9f} | "
            f"{row['residual_m']:.2e} | {row['intact_library_gz_m']:.9f} | "
            f"{row['wall_sided_gz_m']:.9f} | {row['waterline_z_m']:.6f} | {row['kb_m']:.6f} |"
        )
    (outputs / "table.md").write_text("\n".join(table) + "\n")

    _figures(metrics, intact, damaged)
    print(json.dumps({k: v for k, v in metrics.items() if k != "rows"}, indent=2))
    return metrics


def _figures(metrics: dict, intact, damaged) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    target = REPO / "docs" / "damage_gz"
    target.mkdir(parents=True, exist_ok=True)
    heels = [row["heel_deg"] for row in metrics["rows"]]

    fig, axis = plt.subplots(figsize=(7.0, 4.6))
    axis.plot(heels, [row["library_gz_m"] for row in metrics["rows"]], "o-", label="damaged GZ (section integration)")
    axis.plot(heels, [row["oracle_gz_m"] for row in metrics["rows"]], "x", ms=8, label="damaged GZ (exact oracle)")
    axis.plot(heels, [row["intact_library_gz_m"] for row in metrics["rows"]], "s--", label="intact GZ")
    axis.plot(heels, [row["wall_sided_gz_m"] for row in metrics["rows"]], ":", label="wall-sided closed form (intact)")
    axis.axvline(metrics["intact"]["wall_sided_limit_deg"], color="0.5", lw=0.8)
    axis.text(metrics["intact"]["wall_sided_limit_deg"] + 0.6, 0.05, "wall-sided form\nstops here", fontsize=7, color="0.35")
    axis.set_xlabel("heel angle [deg]")
    axis.set_ylabel("GZ [m]")
    axis.set_title("Box barge, midships compartment flooded ($\\mu$=0.6), KG=0.5 m")
    axis.grid(alpha=0.3)
    axis.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(target / "arm_curves.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.8))
    residuals = [abs(row["residual_m"]) for row in metrics["rows"]]
    axes[0].semilogy(heels, [max(value, 1e-18) for value in residuals], "o-")
    axes[0].set_xlabel("heel angle [deg]")
    axes[0].set_ylabel("|library - oracle| [m]")
    axes[0].set_title("the integrator versus the exact answer")
    axes[0].grid(alpha=0.3, which="both")
    axes[0].axhline(1e-12, color="0.4", lw=0.8, ls="--")
    axes[0].text(2.0, 2e-12, "1e-12 m", fontsize=7, color="0.35")
    axes[1].plot(heels, [row["waterline_z_m"] for row in metrics["rows"]], "o-", label="waterline z")
    axes[1].plot(heels, [row["kb_m"] for row in metrics["rows"]], "s-", label="KB above keel")
    axes[1].set_xlabel("heel angle [deg]")
    axes[1].set_ylabel("[m]")
    axes[1].set_title("the ship sinks as it heels")
    axes[1].grid(alpha=0.3)
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(target / "residuals_and_sinkage.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()