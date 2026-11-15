#!/usr/bin/env python
"""Generate publication-quality figures for the baseline fleet verification report.

Reads evidence.json and produces:
1. hull_profiles.png — waterline half-breadth profiles (top view)
2. body_plans.png — midship section comparison
3. convergence.png — volume convergence on three grids
4. volume_errors.png — volume error vs 1% tolerance bar chart
5. fleet_summary.png — comprehensive overview table
"""

import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

EVIDENCE_PATH = Path(__file__).parent.parent.parent / "runs" / "baseline_fleet_m21286" / "evidence.json"
MANIFEST_PATH = Path(__file__).parent.parent.parent / "tensorcad" / "fleet" / "data" / "fleet_manifest.json"
OUTPUT_DIR = Path(__file__).parent.parent.parent / "artifacts" / "baseline_fleet_m21286"


def load_data():
    evidence = json.loads(EVIDENCE_PATH.read_text())
    manifest = json.loads(MANIFEST_PATH.read_text())
    return evidence, manifest


def fig_hull_profiles(evidence, manifest):
    """Waterline half-breadth profiles (top view) for all 6 ships."""
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()

    colors = {"KCS": "#1f77b4", "KVLCC2": "#ff7f0e", "ONRT": "#2ca02c",
              "JBC": "#d62728", "DTC": "#9467bd", "Series60": "#8c564b"}

    for idx, ship_data in enumerate(evidence["ships"]):
        name = ship_data["ship"]
        record = None
        for s in manifest["ships"]:
            if s["name"] == name:
                record = s
                break
        if record is None:
            continue

        ax = axes[idx]
        L = record["length"]
        B = record["beam"]
        T = record["draft"]
        Cb = record["block_coefficient"]

        # Generate profile from the section law
        nx, nz = 241, 49
        xs = np.linspace(0, L, nx)
        zs = np.linspace(-T, 0, nz)

        if ship_data["route"] == "series60_form":
            # Series 60 section law
            Cm = record.get("midship_coefficient", 0.98)
            m = Cm / (1.0 - Cm) if Cm < 1.0 else 1.0
            # Simplified profile at mid-draft
            z_mid = -T / 2
            xi = xs / L
            # Deck taper
            taper = np.where(xi <= 0.5, (2 * xi) ** 0.5, (2 * (1 - xi)) ** 0.5)
            # Vertical distribution
            z_norm = abs(z_mid) / T
            breadth = 0.5 * B * taper * (z_norm ** (1.0 / m))
        else:
            # Parametric sections - simplified
            Cm = record.get("midship_coefficient", 0.98)
            m = Cm / (1.0 - Cm) if Cm < 1.0 else 1.0
            xi = xs / L
            # Beta-like distribution
            a, b = 0.5, 0.5  # simplified
            g = np.where((xi > 0) & (xi < 1),
                        (xi ** a) * ((1 - xi) ** b), 0)
            g_max = g.max()
            if g_max > 0:
                g /= g_max
            z_norm = 0.5
            breadth = 0.5 * B * g * (z_norm ** (1.0 / m))

        ax.plot(xs, breadth, color=colors.get(name, "black"), linewidth=2, label="Starboard")
        ax.plot(xs, -breadth, color=colors.get(name, "black"), linewidth=2, label="Port")
        ax.set_title(f"{name}  (Cb={Cb:.3f})", fontsize=11, fontweight="bold")
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_aspect("equal")
        ax.axhline(0, color="gray", linewidth=0.5)
        ax.axvline(L / 2, color="gray", linewidth=0.5, linestyle="--", alpha=0.5)
        ax.grid(True, alpha=0.3)

    plt.suptitle("Baseline Fleet — Waterline Half-Breadth Profiles (Top View, z = T/2)",
                 fontsize=14, fontweight="bold", y=1.01)
    plt.tight_layout()
    out = OUTPUT_DIR / "hull_profiles.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {out}")


def fig_body_plans(evidence, manifest):
    """Midship section comparison."""
    fig, ax = plt.subplots(figsize=(10, 8))

    colors = {"KCS": "#1f77b4", "KVLCC2": "#ff7f0e", "ONRT": "#2ca02c",
              "JBC": "#d62728", "DTC": "#9467bd", "Series60": "#8c564b"}

    for ship_data in evidence["ships"]:
        name = ship_data["ship"]
        record = None
        for s in manifest["ships"]:
            if s["name"] == name:
                record = s
                break
        if record is None:
            continue

        B = record["beam"]
        T = record["draft"]
        Cm = record.get("midship_coefficient", 0.98)
        m = Cm / (1.0 - Cm) if Cm < 1.0 else 1.0

        nz = 50
        zs = np.linspace(-T, 0, nz)
        z_norm = abs(zs) / T
        breadth = 0.5 * B * (z_norm ** (1.0 / m))

        ax.plot(breadth, zs, color=colors.get(name, "black"), linewidth=2.5,
                label=f"{name} (Cm={Cm:.4f})")
        ax.plot(-breadth, zs, color=colors.get(name, "black"), linewidth=2.5)

    ax.set_xlabel("Half-breadth y [m]", fontsize=12)
    ax.set_ylabel("z [m] (keel at bottom)", fontsize=12)
    ax.set_title("Baseline Fleet — Midship Section Comparison", fontsize=14, fontweight="bold")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_aspect("equal")

    plt.tight_layout()
    out = OUTPUT_DIR / "body_plans.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {out}")


def fig_convergence(evidence, manifest):
    """Volume convergence on three grids."""
    fig, ax = plt.subplots(figsize=(10, 6))

    colors = {"KCS": "#1f77b4", "KVLCC2": "#ff7f0e", "ONRT": "#2ca02c",
              "JBC": "#d62728", "DTC": "#9467bd", "Series60": "#8c564b"}

    for ship_data in evidence["ships"]:
        name = ship_data["ship"]
        conv = ship_data.get("convergence", {})
        if not conv or "volumes" not in conv:
            continue

        volumes = conv["volumes"]
        grids = conv["grids"]
        ref_vol = ship_data["record"]["reference_volume"]
        errors = [100.0 * abs(v - ref_vol) / ref_vol for v in volumes]

        grid_labels = [f"{g[0]}×{g[1]}" for g in grids]
        ax.semilogy(grid_labels, errors, "o-", color=colors.get(name, "black"),
                    linewidth=2, markersize=8, label=f"{name} (p={conv.get('observed_order', '?'):.2f})")

    ax.axhline(1.0, color="red", linestyle="--", linewidth=2, alpha=0.7, label="1% tolerance")
    ax.axhline(0.1, color="green", linestyle="--", linewidth=1.5, alpha=0.5, label="0.1%")
    ax.set_xlabel("Grid Resolution (nx × nz)", fontsize=12)
    ax.set_ylabel("Volume Error [%]", fontsize=12)
    ax.set_title("Baseline Fleet — Volume Convergence (log scale)", fontsize=14, fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, which="both")

    plt.tight_layout()
    out = OUTPUT_DIR / "convergence.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {out}")


def fig_volume_errors(evidence, manifest):
    """Volume error bar chart vs 1% tolerance."""
    fig, ax = plt.subplots(figsize=(10, 5))

    names = [s["ship"] for s in evidence["ships"]]
    errors = [s.get("volume_relative_error", 0) * 100 for s in evidence["ships"]]
    ref_bases = [s.get("reference_basis", "") for s in evidence["ships"]]
    passed = [s.get("passed", False) for s in evidence["ships"]]

    colors = []
    for p, base in zip(passed, ref_bases):
        if p:
            colors.append("#2ca02c" if base == "literature_value_recalled" else "#1f77b4")
        else:
            colors.append("#d62728")

    bars = ax.bar(range(len(names)), errors, color=colors, edgecolor="black", linewidth=0.5)
    ax.axhline(1.0, color="red", linestyle="--", linewidth=2, alpha=0.7, label="1% tolerance")

    for i, (name, err, p) in enumerate(zip(names, errors, passed)):
        ax.text(i, err + 0.05, f"{err:.3f}%", ha="center", va="bottom",
                fontsize=9, fontweight="bold" if p else "normal")

    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, fontsize=10)
    ax.set_ylabel("Volume Error [%]", fontsize=12)
    ax.set_title("Baseline Fleet — Volume Error vs 1% Tolerance (241×49)",
                 fontsize=14, fontweight="bold")
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, axis="y")
    ax.set_ylim(0, max(errors) * 1.3)

    plt.tight_layout()
    out = OUTPUT_DIR / "volume_errors.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {out}")


def fig_fleet_summary(evidence, manifest):
    """Comprehensive fleet overview."""
    fig, ax = plt.subplots(figsize=(16, 8))
    ax.axis("off")

    # Build table data
    headers = ["Ship", "Class", "L [m]", "B [m]", "T [m]", "Cb", "Cm",
               "V_ref [m³]", "V_mesh [m³]", "Error [%]", "Triangles", "Euler", "Status"]

    rows = []
    for ship_data in evidence["ships"]:
        name = ship_data["ship"]
        record = ship_data.get("record", {})
        topo = ship_data.get("topology", {})
        vol_err = ship_data.get("volume_relative_error")
        err_str = f"{vol_err*100:.3f}" if vol_err is not None else "N/A"
        status = "✅ PASS" if ship_data.get("passed") else "❌ FAIL"

        rows.append([
            name,
            record.get("ship_class", "")[:20],
            f"{record.get('length', 0):.1f}",
            f"{record.get('beam', 0):.1f}",
            f"{record.get('draft', 0):.1f}",
            f"{record.get('block_coefficient', 0):.4f}",
            f"{record.get('midship_coefficient', 0):.4f}",
            f"{record.get('reference_volume', 0):,.0f}",
            f"{topo.get('signed_volume', 0):,.1f}",
            err_str,
            f"{topo.get('triangle_count', 0):,}",
            topo.get("euler_characteristic", "?"),
            status,
        ])

    table = ax.table(cellText=rows, colLabels=headers, loc="center",
                     cellLoc="center", bbox=[0, 0, 1, 1])
    table.auto_set_font_size(False)
    table.set_fontsize(9)

    # Style header
    for j in range(len(headers)):
        table[0, j].set_facecolor("#1f77b4")
        table[0, j].set_text_props(color="white", fontweight="bold")

    # Color status column
    for i in range(1, len(rows) + 1):
        if "PASS" in rows[i - 1][-1]:
            table[i, -1].set_facecolor("#d4edda")
        else:
            table[i, -1].set_facecolor("#f8d7da")

    ax.set_title("Baseline Fleet — Comprehensive Verification Summary (241×49)",
                 fontsize=14, fontweight="bold", pad=20)

    plt.tight_layout()
    out = OUTPUT_DIR / "fleet_summary.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {out}")


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    evidence, manifest = load_data()

    print("Generating fleet figures...")
    fig_hull_profiles(evidence, manifest)
    fig_body_plans(evidence, manifest)
    fig_convergence(evidence, manifest)
    fig_volume_errors(evidence, manifest)
    fig_fleet_summary(evidence, manifest)
    print("Done.")


if __name__ == "__main__":
    main()
