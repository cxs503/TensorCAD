"""Generate figures for CASE-CAD-CURVATURE-SPHERE benchmark report."""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[3]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def make_convergence_figure(result: dict, out: Path) -> None:
    """Figure 1: Convergence of area-weighted mean curvature vs mesh resolution."""
    g2 = result["gates"]["G2_integral_marching_tet"]
    conv = g2["convergence"]

    fig, ax = plt.subplots(figsize=(8, 5))
    levels = [c["level"] for c in conv]
    errors = [c["error_pct"] for c in conv]
    faces = [c["faces"] for c in conv]

    ax.semilogy(levels, errors, "o-", color="#1f77b4", linewidth=2, markersize=8)
    ax.axhline(y=2.0, color="red", linestyle="--", linewidth=1.5, label="2% tolerance")
    ax.set_xlabel("Marching-tet resolution level", fontsize=12)
    ax.set_ylabel("Relative error (%)", fontsize=12)
    ax.set_title("Sphere mean curvature convergence (area-weighted, S=2/R convention)", fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    # Add face count annotations
    for i, (lev, err, fc) in enumerate(zip(levels, errors, faces)):
        ax.annotate(f"{fc} faces", (lev, err), textcoords="offset points",
                    xytext=(10, 5), fontsize=8, color="gray")

    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def make_error_cloud_figure(result: dict, out: Path) -> None:
    """Figure 2: Pointwise error distribution (raw vs conditioned)."""
    g3 = result["gates"]["G3_pointwise_conditioned"]
    g7 = result["gates"]["G7_raw_diagnostic"]

    fig, ax = plt.subplots(figsize=(8, 5))

    categories = ["Raw\nmarching-tet", "Taubin+Project\nconditioned"]
    medians = [g3["raw"]["median"], g3["conditioned"]["median"]]
    p95s = [g3["raw"]["p95"], g3["conditioned"]["p95"]]
    errors_med = [abs(m - 2.0) / 2.0 * 100 for m in medians]
    errors_p95 = [abs(p - 2.0) / 2.0 * 100 for p in p95s]

    x = np.arange(len(categories))
    width = 0.35

    bars1 = ax.bar(x - width/2, errors_med, width, label="Median error", color="#1f77b4")
    bars2 = ax.bar(x + width/2, errors_p95, width, label="P95 error", color="#ff7f0e")

    ax.axhline(y=2.0, color="red", linestyle="--", linewidth=1.5, label="2% tolerance")
    ax.set_ylabel("Relative error (%)", fontsize=12)
    ax.set_title("Pointwise mean curvature error: raw vs conditioned mesh", fontsize=13)
    ax.set_xticks(x)
    ax.set_xticklabels(categories, fontsize=10)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, axis="y")

    # Add value labels
    for bar in bars1:
        height = bar.get_height()
        ax.annotate(f"{height:.1f}%", xy=(bar.get_x() + bar.get_width()/2, height),
                    xytext=(0, 5), textcoords="offset points", ha="center", fontsize=9)
    for bar in bars2:
        height = bar.get_height()
        ax.annotate(f"{height:.1f}%", xy=(bar.get_x() + bar.get_width()/2, height),
                    xytext=(0, 5), textcoords="offset points", ha="center", fontsize=9)

    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def make_readings_figure(result: dict, out: Path) -> None:
    """Figure 3: Multiple radii — measured vs reference."""
    g5 = result["gates"]["G5_multiple_radii"]

    fig, ax = plt.subplots(figsize=(8, 5))
    radii = []
    measured = []
    reference = []
    for row in g5["rows"]:
        # Extract R from quantity string like "sphere R=0.5, area-weighted S"
        qty = row["quantity"]
        R = float(qty.split("R=")[1].split(",")[0])
        radii.append(R)
        measured.append(row["measured"])
        reference.append(row["reference"])

    ax.plot(radii, reference, "s-", color="red", linewidth=2, label="Reference 2/R", markersize=8)
    ax.plot(radii, measured, "o-", color="#1f77b4", linewidth=2, label="Measured (area-weighted)", markersize=8)
    ax.set_xlabel("Sphere radius R (m)", fontsize=12)
    ax.set_ylabel("Mean curvature sum S = 2/R (1/m)", fontsize=12)
    ax.set_title("Sphere mean curvature vs radius: reference vs measured", fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    # Add error annotations
    for R, m, ref in zip(radii, measured, reference):
        err = abs(m - ref) / ref * 100
        ax.annotate(f"{err:.3f}%", (R, m), textcoords="offset points",
                    xytext=(5, -12), fontsize=8, color="green" if err < 2 else "red")

    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def main() -> None:
    result_path = Path("runs/case_cad_curvature_sphere/result.json")
    if not result_path.exists():
        print("Run the benchmark first: python benchmarks/verified/sphere_mean_curvature/run_benchmark.py --out runs/case_cad_curvature_sphere")
        sys.exit(1)

    with open(result_path) as f:
        result = json.load(f)

    fig_dir = Path("benchmarks/verified/sphere_mean_curvature/figures")
    fig_dir.mkdir(parents=True, exist_ok=True)

    make_convergence_figure(result, fig_dir / "convergence.png")
    make_error_cloud_figure(result, fig_dir / "error_cloud.png")
    make_readings_figure(result, fig_dir / "readings.png")

    print(f"Figures written to {fig_dir}")


if __name__ == "__main__":
    main()
