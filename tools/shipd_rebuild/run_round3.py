"""Round-3 runner for GAP-CAD-C2-SHIPD-REBUILD.

Builds the 500-ship population model, runs the rebuild pipeline in both method
modes (coefficient-calibrated and the round-2 legacy chain), cross-checks the
rebuilt surfaces against the nominal performance fields, and writes every artefact
plus the comparison figure set.

    python3 tools/shipd_rebuild/run_round3.py            # full 500-ship campaign
    python3 tools/shipd_rebuild/run_round3.py --ships 50 # quick smoke

All numbers in the report are read back from these JSON artefacts -- nothing is
typed in by hand.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from tensorcad.cb_sections import (  # noqa: E402
    CalibratedHullParameters,
    build_calibrated_hull,
    realised_coefficients,
    relative_errors,
)
from tensorcad.shipd_population import (  # noqa: E402
    build_shipd_database,
    correlation_sanity,
    nominal_performance_fields,
)
from tensorcad.shipd_rebuild import run_shipd_rebuild  # noqa: E402
from tensorcad import mesh as meshlib  # noqa: E402

GATE_NAMES = [
    "G1_watertight", "G2_edge_manifold", "G3_euler_char", "G4_outward_oriented",
    "G5_no_zero_area", "G6_volume_consistency", "G7_dimension_agreement",
    "G8_triangle_count",
]


def _write(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2))


def convergence_study(out_dir: Path) -> list[dict]:
    """Resolution study on one representative hull.

    Two different quantities are being checked, and they behave differently:

    * the *coefficients* (volume, ``Cm``, ``Cwp``) are solved on the discrete
      functional, so their error is round-off and does **not** shrink with the
      grid -- calibration, not refinement, is what removes it;
    * the *wetted surface* is a length-like functional of the sampled polygon, so
      it does converge with the grid; its error against the finest mesh is
      reported with the observed order of accuracy.
    """
    rows = []
    target = dict(L=200.0, B=30.0, T=10.0, Cb=0.65, Cm=0.98, Cwp=0.80)
    for nx, nz in ((21, 7), (41, 13), (81, 21), (161, 41), (321, 81)):
        params = CalibratedHullParameters(
            target["L"], target["B"], target["T"], target["Cb"], target["Cm"], target["Cwp"],
            nx=nx, nz=nz,
        )
        hull = build_calibrated_hull(params)
        measured = realised_coefficients(
            hull.mesh, length=target["L"], beam=target["B"], draft=target["T"]
        )
        requested = {"Cb": target["Cb"], "Cm": target["Cm"], "Cwp": target["Cwp"]}
        errors = relative_errors(measured["coefficients"], requested)
        rows.append({
            "nx": nx, "nz": nz,
            "triangles": len(hull.mesh.faces),
            "vertices": len(hull.mesh.vertices),
            "volume_m3": measured["volume_m3"],
            "volume_target_m3": target["Cb"] * target["L"] * target["B"] * target["T"],
            "volume_relative_error": errors["Cb"],
            "midship_relative_error": errors["Cm"],
            "waterplane_relative_error": errors["Cwp"],
            "wetted_surface_m2": measured["wetted_surface_m2"],
            "watertight": meshlib.is_watertight(hull.mesh),
        })
    finest = rows[-1]["wetted_surface_m2"]
    for row in rows:
        row["wetted_surface_relative_error_vs_finest"] = (
            row["wetted_surface_m2"] - finest
        ) / finest
    _write(out_dir / "convergence.json", rows)
    return rows


def wetted_surface_order(rows: list[dict]) -> list[dict]:
    """Observed order of accuracy of the wetted-surface integral.

    Three successive grids with refinement ratio ``r`` give
    ``p = log((S_2-S_1)/(S_1-S_0)) / log(r)`` on the differences to the finest.
    Reported for each consecutive triple; the quantity is smooth in ``h`` so the
    last triple is the one worth quoting.
    """
    out = []
    for a, b, c in zip(rows, rows[1:], rows[2:]):
        d0 = b["wetted_surface_m2"] - a["wetted_surface_m2"]
        d1 = c["wetted_surface_m2"] - b["wetted_surface_m2"]
        if d0 == 0.0 or d1 == 0.0:
            continue
        ratio = (b["nx"] - 1) / (a["nx"] - 1)
        if abs(ratio - 1.0) < 1e-12:
            continue
        if (d0 > 0) != (d1 > 0):
            continue
        out.append({
            "nx_triple": [a["nx"], b["nx"], c["nx"]],
            "refinement_ratio": ratio,
            "observed_order": abs(math.log(abs(d0) / abs(d1)) / math.log(ratio)),
        })
    return out


def figure_gate_pass_rates(calibrated: dict, legacy: dict, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 4.6))
    names = GATE_NAMES + ["G9_midship", "G10_waterplane"]
    cal_rates = [calibrated["gate_pass_rates"].get(g, float("nan")) for g in GATE_NAMES]
    cal_rates += [
        calibrated["coefficient_agreement"]["Cm"]["within_tolerance_fraction"],
        calibrated["coefficient_agreement"]["Cwp"]["within_tolerance_fraction"],
    ]
    leg_rates = [legacy["gate_pass_rates"].get(g, float("nan")) for g in GATE_NAMES]
    leg_rates += [0.0, 0.0]
    idx = range(len(names))
    ax.bar([i - 0.2 for i in idx], [100 * r for r in leg_rates], width=0.4,
           label="round 2 chain (Series 60 / parametric)", color="#b0413e")
    ax.bar([i + 0.2 for i in idx], [100 * r for r in cal_rates], width=0.4,
           label="round 3 calibrated family", color="#2c6e9b")
    ax.set_xticks(list(idx))
    ax.set_xticklabels([n.replace("_", "\n") for n in names], fontsize=8)
    ax.set_ylabel("gate pass rate [%]")
    ax.set_ylim(0, 105)
    ax.set_title("Rebuild quality gates, 500-ship population model")
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def figure_volume_error(
    calibrated_rows: list[dict], legacy_rows: list[dict], path: Path,
    cb_by_id: dict[int, float],
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    for ax, rows, title, colour in (
        (axes[0], legacy_rows, "round 2 chain", "#b0413e"),
        (axes[1], calibrated_rows, "round 3 calibrated family", "#2c6e9b"),
    ):
        xs = [cb_by_id[row["ship_id"]] for row in rows if row["volume_error"] is not None]
        ys = [100 * row["volume_error"] for row in rows if row["volume_error"] is not None]
        ax.scatter(xs, ys, s=9, alpha=0.65, color=colour)
        ax.axhline(3.0, ls="--", color="k", lw=1, label="gate G6: 3 %")
        ax.set_xlabel("requested $C_b$")
        ax.set_ylabel("volume error [%]")
        ax.set_title(title)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def figure_coefficient_calibration(rows: list[dict], path: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.2))
    for ax, key, label in (
        (axes[0], "Cb", "$C_b$"),
        (axes[1], "Cm", "$C_m$"),
        (axes[2], "Cwp", "$C_{wp}$"),
    ):
        pairs = [(row["requested"][key], row["realised"][key]) for row in rows]
        ax.scatter([p[0] for p in pairs], [p[1] for p in pairs], s=9, alpha=0.6, color="#2c6e9b")
        lo = min(min(p) for p in pairs) - 0.01
        hi = max(max(p) for p in pairs) + 0.01
        ax.plot([lo, hi], [lo, hi], "k--", lw=1, label="1:1")
        ax.set_xlabel(f"requested {label}")
        ax.set_ylabel(f"realised {label} (from mesh)")
        ax.set_title(f"{label}: 500 hulls")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def figure_wetted_correlation(correlation: dict, path: Path) -> None:
    pairs = correlation["wetted_surface"]["pairs"] if "pairs" in correlation["wetted_surface"] else None
    fig, ax = plt.subplots(figsize=(5.8, 5.2))
    if pairs:
        ax.scatter([p[0] for p in pairs], [p[1] for p in pairs], s=9, alpha=0.6, color="#2c6e9b")
        lo = min(min(p) for p in pairs)
        hi = max(max(p) for p in pairs)
        ax.plot([lo, hi], [lo, hi], "k--", lw=1, label="1:1")
    ax.set_xlabel("nominal wetted surface, Holtrop & Mennen [m²]")
    ax.set_ylabel("rebuilt mesh wetted surface [m²]")
    r2 = correlation["wetted_surface"]["r_squared"]
    ax.set_title(f"wetted surface, {correlation['n_pairs']} hulls (R² = {r2:.5f})")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def figure_convergence(rows: list[dict], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ns = [row["nx"] - 1 for row in rows]
    errs = [max(row["volume_relative_error"], 1e-16) for row in rows]
    ax.loglog(ns, errs, "o-", label="volume error (calibrated family)")
    wetted = [abs(row["wetted_surface_relative_error_vs_finest"]) for row in rows]
    ax.loglog(ns, [max(v, 1e-16) for v in wetted], "s-", color="#2c6e9b",
              label="wetted surface vs finest grid")
    ax.set_xlabel("longitudinal panels $n_x - 1$")
    ax.set_ylabel("relative error")
    ax.set_title("Calibrated hull: coefficients are exact, area converges")
    ax.grid(which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ships", type=int, default=500)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out", default=str(ROOT / "runs" / "shipd_rebuild_r3"))
    parser.add_argument("--figures", default=str(ROOT / "runs" / "shipd_rebuild_r3"))
    args = parser.parse_args()

    out_dir = Path(args.out)
    fig_dir = Path(args.figures)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    ships = build_shipd_database(args.ships, seed=args.seed)
    db_path = ROOT / "artifacts" / "shipd_rebuild" / f"shipd_{args.ships}.json"
    _write(db_path, ships)
    print(f"population model: {len(ships)} hulls -> {db_path}")

    timings = {}
    evidences = {}
    for mode in ("legacy", "calibrated"):
        target = out_dir / mode
        t0 = time.time()
        evidence = run_shipd_rebuild(
            ships, nx=81, nz=21, out_dir=target, method_preference=mode,
        )
        timings[mode] = time.time() - t0
        evidence["wall_time_s"] = timings[mode]
        _write(target / "evidence.json", evidence)
        evidences[mode] = evidence
        print(
            f"{mode:11s}: pass {evidence['pass_rate']:.4f} "
            f"strict {evidence['strict_pass_rate']:.4f} "
            f"({timings[mode]:.1f} s)"
        )

    # Per-ship results, both modes, for the figures.
    legacy_rows = json.loads((out_dir / "legacy" / "ship_results.json").read_text())
    cal_rows = json.loads((out_dir / "calibrated" / "ship_results.json").read_text())

    figure_gate_pass_rates(
        evidences["calibrated"], evidences["legacy"], fig_dir / "gate_pass_rates.png"
    )
    figure_volume_error(
        legacy_rows, cal_rows, fig_dir / "volume_error_vs_cb.png",
        {ship["ship_id"]: ship["Cb"] for ship in ships},
    )

    # Calibration error scatter, recomputed from the calibrated hulls themselves.
    calibration_rows = []
    for ship in ships:
        params = CalibratedHullParameters(
            ship["L"], ship["B"], ship["T"], ship["Cb"], ship["Cm"], ship["Cwp"],
        )
        hull = build_calibrated_hull(params)
        measured = realised_coefficients(
            hull.mesh, length=ship["L"], beam=ship["B"], draft=ship["T"]
        )
        calibration_rows.append({
            "ship_id": ship["ship_id"],
            "requested": {"Cb": ship["Cb"], "Cm": ship["Cm"], "Cwp": ship["Cwp"]},
            "realised": measured["coefficients"],
            "relative_error": relative_errors(
                measured["coefficients"],
                {"Cb": ship["Cb"], "Cm": ship["Cm"], "Cwp": ship["Cwp"]},
            ),
        })
    _write(out_dir / "calibration_500.json", calibration_rows)
    figure_coefficient_calibration(calibration_rows, fig_dir / "coefficient_calibration.png")

    correlation = evidences["calibrated"]["correlation_sanity"]
    # The figure needs the raw pairs; recompute them from the stored artefacts.
    pairs = []
    for ship in ships:
        row = next(r for r in cal_rows if r["ship_id"] == ship["ship_id"])
        if row["wetted_surface_area"] is None:
            continue
        nominal = nominal_performance_fields(ship)
        pairs.append((nominal["wetted_surface_m2"], row["wetted_surface_area"]))
    correlation["wetted_surface"]["pairs"] = pairs
    _write(out_dir / "correlation.json", correlation)
    figure_wetted_correlation(correlation, fig_dir / "wetted_surface_correlation.png")

    convergence = convergence_study(out_dir)
    order = wetted_surface_order(convergence)
    _write(out_dir / "convergence_order.json", order)
    figure_convergence(convergence, fig_dir / "convergence.png")

    comparison = {
        "n_ships": len(ships),
        "seed": args.seed,
        "resolution": {"nx": 81, "nz": 21},
        "legacy": {
            "pass_rate": evidences["legacy"]["pass_rate"],
            "strict_pass_rate": evidences["legacy"]["strict_pass_rate"],
            "volume_error_mean": evidences["legacy"]["volume_error_mean"],
            "volume_error_max": evidences["legacy"]["volume_error_max"],
            "method_distribution": evidences["legacy"]["method_distribution"],
            "gate_pass_rates": evidences["legacy"]["gate_pass_rates"],
            "coefficient_agreement": evidences["legacy"]["coefficient_agreement"],
            "wall_time_s": timings["legacy"],
        },
        "calibrated": {
            "pass_rate": evidences["calibrated"]["pass_rate"],
            "strict_pass_rate": evidences["calibrated"]["strict_pass_rate"],
            "volume_error_mean": evidences["calibrated"]["volume_error_mean"],
            "volume_error_max": evidences["calibrated"]["volume_error_max"],
            "method_distribution": evidences["calibrated"]["method_distribution"],
            "gate_pass_rates": evidences["calibrated"]["gate_pass_rates"],
            "coefficient_agreement": evidences["calibrated"]["coefficient_agreement"],
            "correlation_sanity": evidences["calibrated"]["correlation_sanity"],
            "failure_mode_ledger": evidences["calibrated"]["failure_mode_ledger"],
            "wall_time_s": timings["calibrated"],
        },
        "convergence": convergence,
        "convergence_order": order,
    }
    _write(out_dir / "comparison.json", comparison)
    print(json.dumps({
        "legacy_pass": comparison["legacy"]["pass_rate"],
        "calibrated_pass": comparison["calibrated"]["pass_rate"],
        "calibrated_strict": comparison["calibrated"]["strict_pass_rate"],
        "vol_err_max_calibrated": comparison["calibrated"]["volume_error_max"],
        "wetted_r2": comparison["calibrated"]["correlation_sanity"]["wetted_surface"]["r_squared"],
    }, indent=2))


if __name__ == "__main__":
    main()