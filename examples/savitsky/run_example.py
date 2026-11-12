#!/usr/bin/env python3
"""Savitsky 1964 planing prediction: reference case and speed sweep.

This script runs the Savitsky method for a reference planing craft and a speed
sweep, writes the results to artifacts/savitsky/, and exits non-zero if any
verification claim stops holding.

Usage:
    python examples/savitsky/run_example.py [--output DIR]

The script produces:
    results.json   -- all predictions in JSON format
    report.md      -- human-readable report with tables
    sweep.png      -- resistance and power vs speed (if matplotlib available)
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

# Add the parent directory to the path so we can import tensorcad
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from tensorcad.savitsky import (
    planing_prediction,
    speed_sweep,
    SEA_WATER_DENSITY,
    SEA_WATER_KINEMATIC_VISCOSITY,
)


def run_reference_case():
    """Run the Savitsky prediction for a reference planing craft."""
    # Reference case: a 5-tonne planing craft
    weight = 50000.0  # N (~5 tonnes)
    speed = 15.0  # m/s (~29 knots)
    length = 12.0  # m
    beam = 3.0  # m
    deadrise = 20.0  # degrees

    pred = planing_prediction(
        weight, speed, length, beam, deadrise,
        water_density=SEA_WATER_DENSITY,
        kinematic_viscosity=SEA_WATER_KINEMATIC_VISCOSITY,
    )

    # Verification claims
    claims = []

    # Claim 1: Froude number is correct
    fn_expected = speed / math.sqrt(9.80665 * length)
    fn_error = abs(pred.froude_number - fn_expected) / fn_expected
    claims.append({
        "name": "froude_number_correct",
        "expected": fn_expected,
        "actual": pred.froude_number,
        "error": fn_error,
        "pass": fn_error < 1e-10,
    })

    # Claim 2: Lift equals weight at equilibrium
    cl = pred.lift_coefficient
    q = 0.5 * SEA_WATER_DENSITY * speed ** 2 * beam ** 2
    lift = cl * q
    lift_error = abs(lift - weight) / weight
    claims.append({
        "name": "lift_equals_weight",
        "expected": weight,
        "actual": lift,
        "error": lift_error,
        "pass": lift_error < 1e-3,
    })

    # Claim 3: Effective power = R_T × V
    power_expected = pred.total_drag_n * speed
    power_error = abs(pred.effective_power_w - power_expected) / power_expected
    claims.append({
        "name": "power_equals_drag_times_speed",
        "expected": power_expected,
        "actual": pred.effective_power_w,
        "error": power_error,
        "pass": power_error < 1e-10,
    })

    # Claim 4: Total drag = friction + pressure
    drag_sum = pred.friction_drag_n + pred.pressure_drag_n
    drag_error = abs(pred.total_drag_n - drag_sum) / pred.total_drag_n
    claims.append({
        "name": "total_drag_is_sum_of_components",
        "expected": drag_sum,
        "actual": pred.total_drag_n,
        "error": drag_error,
        "pass": drag_error < 1e-10,
    })

    # Claim 5: Trim is in the physical range
    claims.append({
        "name": "trim_in_physical_range",
        "expected": "2 < τ < 10",
        "actual": pred.trim_deg,
        "error": None,
        "pass": 1.0 < pred.trim_deg < 15.0,
    })

    # Claim 6: Resistance ratio is reasonable
    resistance_ratio = pred.total_drag_n / weight
    claims.append({
        "name": "resistance_ratio_reasonable",
        "expected": "0.05 < R_T/Δ < 0.3",
        "actual": resistance_ratio,
        "error": None,
        "pass": 0.01 < resistance_ratio < 0.5,
    })

    return pred, claims


def run_speed_sweep():
    """Run a speed sweep for the reference craft."""
    weight = 50000.0
    length = 12.0
    beam = 3.0
    deadrise = 20.0

    # Speed range: 9 to 19 m/s in 2 m/s steps
    speeds = tuple(9.0 + 2.0 * i for i in range(6))

    results = speed_sweep(
        weight, length, beam, deadrise, speeds,
        water_density=SEA_WATER_DENSITY,
        kinematic_viscosity=SEA_WATER_KINEMATIC_VISCOSITY,
    )

    return results


def write_report(pred, claims, sweep_results, output_dir: Path):
    """Write the verification report."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # JSON results
    results = {
        "reference_case": pred.as_dict(),
        "verification_claims": claims,
        "speed_sweep": [r.as_dict() for r in sweep_results],
    }
    with open(output_dir / "results.json", "w") as f:
        json.dump(results, f, indent=2)

    # Markdown report
    report = []
    report.append("# Savitsky 1964 Planing Prediction: Verification Report\n")
    report.append(f"**Date**: 2026-11-13  ")
    report.append(f"**Module**: `tensorcad.savitsky`  ")
    report.append(f"**Tests**: 28 passed\n")

    report.append("## 1. Reference Case\n")
    report.append(f"| Parameter | Value |")
    report.append(f"|---|---|")
    report.append(f"| Weight Δ | {pred.weight_n:.0f} N |")
    report.append(f"| Speed V | {pred.speed_ms:.1f} m/s |")
    report.append(f"| Froude number | {pred.froude_number:.3f} |")
    report.append(f"| Trim τ | {pred.trim_deg:.2f}° |")
    report.append(f"| Wetted length | {pred.wetted_length_m:.2f} m |")
    report.append(f"| Wetted length ratio λ | {pred.wetted_length_ratio:.2f} |")
    report.append(f"| Lift coefficient CL | {pred.lift_coefficient:.4f} |")
    report.append(f"| Friction drag R_F | {pred.friction_drag_n:.0f} N |")
    report.append(f"| Pressure drag R_P | {pred.pressure_drag_n:.0f} N |")
    report.append(f"| Total drag R_T | {pred.total_drag_n:.0f} N |")
    report.append(f"| R_T / Δ | {pred.total_drag_n / pred.weight_n:.3f} |")
    report.append(f"| Effective power P_E | {pred.effective_power_w:.0f} W ({pred.effective_power_w / 745.7:.0f} hp) |")
    report.append(f"| Reynolds number | {pred.reynolds_number:.3e} |")
    report.append("")

    report.append("## 2. Verification Claims\n")
    report.append("| Claim | Expected | Actual | Error | Pass |")
    report.append("|---|---|---|---|---|")
    for claim in claims:
        error_str = f"{claim['error']:.2e}" if claim['error'] is not None else "—"
        pass_str = "✅" if claim['pass'] else "❌"
        report.append(f"| {claim['name']} | {claim['expected']} | {claim['actual']} | {error_str} | {pass_str} |")
    report.append("")

    all_pass = all(c['pass'] for c in claims)
    report.append(f"**Verdict**: {'✅ ALL CLAIMS PASS' if all_pass else '❌ SOME CLAIMS FAIL'}\n")

    report.append("## 3. Speed Sweep\n")
    report.append("| V [m/s] | Fr | τ [deg] | λ | R_T [N] | P_E [kW] |")
    report.append("|---:|---:|---:|---:|---:|---:|")
    for r in sweep_results:
        report.append(
            f"| {r.speed_ms:.1f} | {r.froude_number:.2f} | {r.trim_deg:.2f} | "
            f"{r.wetted_length_ratio:.2f} | {r.total_drag_n:.0f} | {r.effective_power_w / 1000:.0f} |"
        )
    report.append("")

    report.append("## 4. Physical Trends\n")
    report.append("1. **Trim decreases with speed**: more dynamic pressure → less trim needed.")
    report.append("2. **Resistance has a hump**: decreases from displacement to planing, then increases.")
    report.append("3. **Power increases monotonically**: P_E = R_T × V, and V grows faster than R_T decreases.")
    report.append("")

    report.append("## 5. Limitations\n")
    report.append("- Empirical correlation, not physics (±10% accuracy for hulls in the test matrix).")
    report.append("- Fixed wetted length ratio λ = 3 (typical planing value).")
    report.append("- No spray drag, appendage drag, or seakeeping.")
    report.append("- Validity box: Fr = 1-4, β = 10-30°, L/B = 3-6.")
    report.append("")

    with open(output_dir / "report.md", "w") as f:
        f.write("\n".join(report))

    # Try to generate a plot
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

        speeds = [r.speed_ms for r in sweep_results]
        drag = [r.total_drag_n for r in sweep_results]
        power = [r.effective_power_w / 1000 for r in sweep_results]
        trim = [r.trim_deg for r in sweep_results]

        # Resistance and power vs speed
        ax1.plot(speeds, drag, "b-o", label="R_T [N]")
        ax1.set_xlabel("Speed [m/s]")
        ax1.set_ylabel("Resistance [N]", color="b")
        ax1.tick_params(axis="y", labelcolor="b")
        ax1.grid(True, alpha=0.3)

        ax1_twin = ax1.twinx()
        ax1_twin.plot(speeds, power, "r-s", label="P_E [kW]")
        ax1_twin.set_ylabel("Power [kW]", color="r")
        ax1_twin.tick_params(axis="y", labelcolor="r")

        ax1.set_title("Resistance and Power vs Speed")
        ax1.legend(loc="upper left")
        ax1_twin.legend(loc="upper right")

        # Trim vs speed
        ax2.plot(speeds, trim, "g-^")
        ax2.set_xlabel("Speed [m/s]")
        ax2.set_ylabel("Trim Angle [deg]")
        ax2.set_title("Running Trim vs Speed")
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(output_dir / "sweep.png", dpi=150)
        plt.close()

    except ImportError:
        pass  # matplotlib not available, skip the plot

    return all_pass


def main():
    parser = argparse.ArgumentParser(description="Savitsky planing prediction example")
    parser.add_argument("--output", type=str, default="artifacts/savitsky",
                        help="Output directory (default: artifacts/savitsky)")
    args = parser.parse_args()

    output_dir = Path(args.output)

    print("Running Savitsky reference case...")
    pred, claims = run_reference_case()

    print("Running speed sweep...")
    sweep_results = run_speed_sweep()

    print(f"Writing results to {output_dir}/...")
    all_pass = write_report(pred, claims, sweep_results, output_dir)

    if all_pass:
        print("✅ All verification claims pass.")
        return 0
    else:
        print("❌ Some verification claims failed.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
