#!/usr/bin/env python3
"""
Deformation Family Visualization Generator for m585.
Generates ≥4 publication-quality figures from m542 (DTMB 5415) and m546 (KCS)
bow deformation families.

Figures:
  1. Baseline vs representative variants overlay – side view (with deformation zone)
  2. Baseline vs representative variants overlay – top view (with deformation zone)
  3. Non-deformation zone deviation cloud map (should be zero)
  4. Displacement & parameter distribution plots
  5. Block coefficient vs variant index
  6. Gate assertion summary heatmap
"""
from __future__ import annotations

import json
import math
import os
import struct
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

# ---------------------------------------------------------------------------
# STL loader (binary)
# ---------------------------------------------------------------------------
def load_stl_binary(path: str):
    """Load binary STL, return (vertices, faces) as numpy arrays."""
    with open(path, "rb") as f:
        header = f.read(80)
        n_tri = struct.unpack("<I", f.read(4))[0]
        dtype = np.dtype([
            ("normal", "<f4", (3,)),
            ("v1", "<f4", (3,)),
            ("v2", "<f4", (3,)),
            ("v3", "<f4", (3,)),
            ("attr", "<u2"),
        ])
        data = np.fromfile(f, dtype=dtype, count=n_tri)
    vertices = np.vstack([data["v1"], data["v2"], data["v3"]])
    faces = np.arange(len(data) * 3).reshape(-1, 3)
    return vertices, faces


def load_stl_to_profile(path: str, n_samples: int = 200):
    """Extract a 2D profile (x, half-breadth) from STL for side/top view plotting."""
    verts, _ = load_stl_binary(path)
    x = verts[:, 0]
    y = np.abs(verts[:, 1])  # half-breadth (symmetric hull)
    z = verts[:, 2]

    # Side view: max half-breadth at each x (waterline profile)
    # Top view: max z at each x (deck/profile)
    x_unique = np.linspace(x.min(), x.max(), n_samples)
    y_profile = np.zeros(n_samples)
    z_profile = np.zeros(n_samples)

    for i, xv in enumerate(x_unique):
        mask = np.abs(x - xv) < (x.max() - x.min()) / n_samples * 1.5
        if mask.any():
            y_profile[i] = y[mask].max()
            z_profile[i] = z[mask].max()

    return x_unique, y_profile, z_profile, x.min(), x.max()


def load_stl_full(path: str):
    """Load full STL and return vertices for deviation analysis."""
    verts, faces = load_stl_binary(path)
    return verts, faces


# ---------------------------------------------------------------------------
# Data loaders
# ---------------------------------------------------------------------------
def load_5415_family(base: str):
    """Load DTMB 5415 deformation family from m542."""
    base = Path(base)
    summary_path = base / "family_summary.json"
    with open(summary_path) as f:
        summary = json.load(f)

    variants = {}
    for vid, vdata in summary["variants"].items():
        vdir = base / vid
        stl_path = vdir / "hull.stl"
        params_path = vdir / "params.json"
        gates_path = vdir / "gates.json"

        params = {}
        if params_path.exists():
            with open(params_path) as f:
                params = json.load(f)

        gates = {}
        if gates_path.exists():
            with open(gates_path) as f:
                gates = json.load(f)

        variants[vid] = {
            "params": params,
            "gates": gates,
            "volume": vdata.get("volume_m3", 0),
            "all_passed": vdata.get("all_passed", False),
            "stl_path": str(stl_path) if stl_path.exists() else None,
        }

    return summary, variants


def load_kcs_family(base: str):
    """Load KCS deformation family from m546."""
    base = Path(base)
    summary_path = base / "variant_summary.json"
    with open(summary_path) as f:
        summary = json.load(f)

    variants = {}
    for vdata in summary["variants"]:
        vid = vdata["id"]
        vdir = base / vid
        stl_files = list(vdir.glob("*.stl"))
        stl_path = str(stl_files[0]) if stl_files else None

        params_files = list(vdir.glob("*_params.json"))
        params = {}
        if params_files:
            with open(params_files[0]) as f:
                params = json.load(f)

        variants[vid] = {
            "params": params,
            "displacement": vdata.get("displacement", 0),
            "Cb": vdata.get("Cb", 0),
            "all_passed": vdata.get("assertions", {}).get("all_passed", False),
            "stl_path": stl_path,
        }

    return summary, variants


# ---------------------------------------------------------------------------
# Figure 1 & 2: Baseline vs representative variants overlay
# ---------------------------------------------------------------------------
def plot_overlay_side(variants_5415, variants_kcs, outpath: str):
    """Side view overlay: baseline + representative deformed variants."""
    fig, axes = plt.subplots(1, 2, figsize=(16, 5))

    # --- 5415 family ---
    ax = axes[0]
    # Baseline (v000)
    v000 = variants_5415["v000"]
    if v000["stl_path"]:
        x, y, z, xmin, xmax = load_stl_to_profile(v000["stl_path"])
        ax.plot(x, z, "k-", linewidth=2.5, label="Baseline (v000)", zorder=10)

    # Representative variants: v005, v010, v015, v020, v025, v029
    reps = ["v005", "v010", "v015", "v020", "v025", "v029"]
    cmap = plt.cm.viridis
    for i, vid in enumerate(reps):
        v = variants_5415.get(vid)
        if v and v["stl_path"]:
            x, y, z, _, _ = load_stl_to_profile(v["stl_path"])
            color = cmap(0.2 + 0.7 * i / len(reps))
            ax.plot(x, z, color=color, linewidth=1.2, alpha=0.7, label=vid)

    # Deformation zone shading
    deform_limit = 3.5  # m from AP
    ax.axvline(x=deform_limit, color="red", linestyle="--", linewidth=1.5, alpha=0.6)
    ax.axvspan(deform_limit, xmax, alpha=0.1, color="red", label="Deformation zone")

    ax.set_xlabel("x [m] (from AP)", fontsize=11)
    ax.set_ylabel("z [m] (vertical, +up)", fontsize=11)
    ax.set_title("DTMB 5415 Bow Deformation Family – Side View", fontsize=13, fontweight="bold")
    ax.legend(fontsize=8, loc="upper left", ncol=2)
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)

    # --- KCS family ---
    ax = axes[1]
    v000_kcs = variants_kcs["kcs_v000"]
    if v000_kcs["stl_path"]:
        x, y, z, xmin, xmax = load_stl_to_profile(v000_kcs["stl_path"])
        ax.plot(x, z, "k-", linewidth=2.5, label="Baseline (kcs_v000)", zorder=10)

    reps_kcs = ["kcs_v005", "kcs_v010", "kcs_v015", "kcs_v020", "kcs_v025", "kcs_v031"]
    for i, vid in enumerate(reps_kcs):
        v = variants_kcs.get(vid)
        if v and v["stl_path"]:
            x, y, z, _, _ = load_stl_to_profile(v["stl_path"])
            color = cmap(0.2 + 0.7 * i / len(reps_kcs))
            ax.plot(x, z, color=color, linewidth=1.2, alpha=0.7, label=vid)

    # KCS deformation zone: x/Lpp > 0.68, Lpp=230m → x > 156.4m
    # But model scale: check actual x range
    deform_x_kcs = xmin + 0.68 * (xmax - xmin)
    ax.axvline(x=deform_x_kcs, color="red", linestyle="--", linewidth=1.5, alpha=0.6)
    ax.axvspan(deform_x_kcs, xmax, alpha=0.1, color="red", label="Deformation zone (x/L>0.68)")

    ax.set_xlabel("x [m] (from AP)", fontsize=11)
    ax.set_ylabel("z [m] (vertical, +up)", fontsize=11)
    ax.set_title("KCS Bow Deformation Family – Side View", fontsize=13, fontweight="bold")
    ax.legend(fontsize=8, loc="upper left", ncol=2)
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {outpath}")


def plot_overlay_top(variants_5415, variants_kcs, outpath: str):
    """Top view overlay: baseline + representative deformed variants."""
    fig, axes = plt.subplots(1, 2, figsize=(16, 5))

    # --- 5415 family ---
    ax = axes[0]
    v000 = variants_5415["v000"]
    if v000["stl_path"]:
        x, y, z, xmin, xmax = load_stl_to_profile(v000["stl_path"])
        ax.plot(x, y, "k-", linewidth=2.5, label="Baseline (v000)", zorder=10)
        ax.plot(x, -y, "k-", linewidth=2.5, zorder=10)

    reps = ["v005", "v010", "v015", "v020", "v025", "v029"]
    cmap = plt.cm.plasma
    for i, vid in enumerate(reps):
        v = variants_5415.get(vid)
        if v and v["stl_path"]:
            x, y, z, _, _ = load_stl_to_profile(v["stl_path"])
            color = cmap(0.2 + 0.7 * i / len(reps))
            ax.plot(x, y, color=color, linewidth=1.2, alpha=0.7)
            ax.plot(x, -y, color=color, linewidth=1.2, alpha=0.7)

    deform_limit = 3.5
    ax.axvline(x=deform_limit, color="red", linestyle="--", linewidth=1.5, alpha=0.6)
    ax.axvspan(deform_limit, xmax, alpha=0.1, color="red", label="Deformation zone")

    ax.set_xlabel("x [m] (from AP)", fontsize=11)
    ax.set_ylabel("y [m] (half-breadth, ±stbd/port)", fontsize=11)
    ax.set_title("DTMB 5415 – Top View (Waterline Profile)", fontsize=13, fontweight="bold")
    ax.legend(fontsize=8, loc="upper left")
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)

    # --- KCS family ---
    ax = axes[1]
    v000_kcs = variants_kcs["kcs_v000"]
    if v000_kcs["stl_path"]:
        x, y, z, xmin, xmax = load_stl_to_profile(v000_kcs["stl_path"])
        ax.plot(x, y, "k-", linewidth=2.5, label="Baseline (kcs_v000)", zorder=10)
        ax.plot(x, -y, "k-", linewidth=2.5, zorder=10)

    reps_kcs = ["kcs_v005", "kcs_v010", "kcs_v015", "kcs_v020", "kcs_v025", "kcs_v031"]
    for i, vid in enumerate(reps_kcs):
        v = variants_kcs.get(vid)
        if v and v["stl_path"]:
            x, y, z, _, _ = load_stl_to_profile(v["stl_path"])
            color = cmap(0.2 + 0.7 * i / len(reps_kcs))
            ax.plot(x, y, color=color, linewidth=1.2, alpha=0.7)
            ax.plot(x, -y, color=color, linewidth=1.2, alpha=0.7)

    deform_x_kcs = xmin + 0.68 * (xmax - xmin)
    ax.axvline(x=deform_x_kcs, color="red", linestyle="--", linewidth=1.5, alpha=0.6)
    ax.axvspan(deform_x_kcs, xmax, alpha=0.1, color="red", label="Deformation zone (x/L>0.68)")

    ax.set_xlabel("x [m] (from AP)", fontsize=11)
    ax.set_ylabel("y [m] (half-breadth, ±stbd/port)", fontsize=11)
    ax.set_title("KCS – Top View (Waterline Profile)", fontsize=13, fontweight="bold")
    ax.legend(fontsize=8, loc="upper left")
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {outpath}")


# ---------------------------------------------------------------------------
# Figure 3: Non-deformation zone deviation cloud map
# ---------------------------------------------------------------------------
def plot_non_deform_deviation(variants_5415, variants_kcs, outpath: str):
    """
    Non-deformation zone deviation: for each variant, compute max deviation
    in the non-deform region. Should be ~0.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # --- 5415 ---
    ax = axes[0]
    deform_limit = 3.5
    vids = sorted(variants_5415.keys())
    max_devs = []
    for vid in vids:
        v = variants_5415[vid]
        if not v["stl_path"]:
            max_devs.append(0)
            continue
        verts, _ = load_stl_full(v["stl_path"])
        # Non-deform region: x <= deform_limit
        mask = verts[:, 0] <= deform_limit
        if mask.sum() == 0:
            max_devs.append(0)
            continue
        # Compare z deviation from baseline
        v0_verts, _ = load_stl_full(variants_5415["v000"]["stl_path"])
        v0_mask = v0_verts[:, 0] <= deform_limit
        # Simple approach: max |z| difference in non-deform region
        z_deform = verts[mask, 2]
        z_base = v0_verts[v0_mask, 2]
        # Resample to same grid for comparison
        x_deform = verts[mask, 0]
        x_base = v0_verts[v0_mask, 0]
        # Use nearest neighbor for quick comparison
        if len(z_base) > 0 and len(z_deform) > 0:
            # Bin by x
            x_bins = np.linspace(x_deform.min(), x_deform.max(), 50)
            dev_vals = []
            for xb in x_bins:
                mask_d = np.abs(x_deform - xb) < (x_bins[1] - x_bins[0]) * 1.5
                mask_b = np.abs(x_base - xb) < (x_bins[1] - x_bins[0]) * 1.5
                if mask_d.any() and mask_b.any():
                    dev_vals.append(abs(z_deform[mask_d].mean() - z_base[mask_b].mean()))
            max_devs.append(max(dev_vals) if dev_vals else 0)
        else:
            max_devs.append(0)

    x_idx = np.arange(len(vids))
    colors = ["green" if d < 1e-6 else "red" for d in max_devs]
    ax.bar(x_idx, max_devs, color=colors, alpha=0.8, width=0.8)
    ax.axhline(y=1e-6, color="red", linestyle="--", linewidth=1, alpha=0.5, label="Threshold (1e-6 m)")
    ax.set_xlabel("Variant index", fontsize=11)
    ax.set_ylabel("Max deviation in non-deform zone [m]", fontsize=10)
    ax.set_title("DTMB 5415 – Non-Deformation Zone Deviation\n(should be ≈ 0)", fontsize=12, fontweight="bold")
    ax.set_xticks(x_idx[::5])
    ax.set_xticklabels([vids[i] for i in x_idx[::5]], rotation=45, fontsize=7)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, axis="y")

    # --- KCS ---
    ax = axes[1]
    vids_kcs = sorted(variants_kcs.keys())
    max_devs_kcs = []
    for vid in vids_kcs:
        v = variants_kcs[vid]
        if not v["stl_path"]:
            max_devs_kcs.append(0)
            continue
        verts, _ = load_stl_full(v["stl_path"])
        # KCS non-deform: x/Lpp <= 0.68
        xmin, xmax = verts[:, 0].min(), verts[:, 0].max()
        deform_x = xmin + 0.68 * (xmax - xmin)
        mask = verts[:, 0] <= deform_x
        if mask.sum() == 0:
            max_devs_kcs.append(0)
            continue
        v0_verts, _ = load_stl_full(variants_kcs["kcs_v000"]["stl_path"])
        v0_mask = v0_verts[:, 0] <= deform_x
        z_deform = verts[mask, 2]
        z_base = v0_verts[v0_mask, 2]
        x_deform = verts[mask, 0]
        x_base = v0_verts[v0_mask, 0]
        if len(z_base) > 0 and len(z_deform) > 0:
            x_bins = np.linspace(x_deform.min(), x_deform.max(), 50)
            dev_vals = []
            for xb in x_bins:
                mask_d = np.abs(x_deform - xb) < (x_bins[1] - x_bins[0]) * 1.5
                mask_b = np.abs(x_base - xb) < (x_bins[1] - x_bins[0]) * 1.5
                if mask_d.any() and mask_b.any():
                    dev_vals.append(abs(z_deform[mask_d].mean() - z_base[mask_b].mean()))
            max_devs_kcs.append(max(dev_vals) if dev_vals else 0)
        else:
            max_devs_kcs.append(0)

    x_idx_kcs = np.arange(len(vids_kcs))
    colors_kcs = ["green" if d < 1e-6 else "red" for d in max_devs_kcs]
    ax.bar(x_idx_kcs, max_devs_kcs, color=colors_kcs, alpha=0.8, width=0.8)
    ax.axhline(y=1e-6, color="red", linestyle="--", linewidth=1, alpha=0.5, label="Threshold (1e-6 m)")
    ax.set_xlabel("Variant index", fontsize=11)
    ax.set_ylabel("Max deviation in non-deform zone [m]", fontsize=10)
    ax.set_title("KCS – Non-Deformation Zone Deviation\n(should be ≈ 0)", fontsize=12, fontweight="bold")
    ax.set_xticks(x_idx_kcs[::5])
    ax.set_xticklabels([vids_kcs[i] for i in x_idx_kcs[::5]], rotation=45, fontsize=7)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {outpath}")

    return max_devs, max_devs_kcs


# ---------------------------------------------------------------------------
# Figure 4: Displacement & parameter distribution
# ---------------------------------------------------------------------------
def plot_param_distribution(variants_5415, variants_kcs, outpath: str):
    """Parameter distribution and displacement scatter plots."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # --- 5415: Volume distribution ---
    ax = axes[0, 0]
    vids = sorted(variants_5415.keys())
    volumes = [variants_5415[v]["volume"] for v in vids]
    ax.bar(np.arange(len(vids)), volumes, color="steelblue", alpha=0.8, width=0.8)
    ax.axhline(y=volumes[0], color="red", linestyle="--", linewidth=1.5, label=f"Baseline: {volumes[0]:.4f} m³")
    ax.set_xlabel("Variant", fontsize=10)
    ax.set_ylabel("Volume [m³]", fontsize=10)
    ax.set_title("DTMB 5415 – Displacement Volume Distribution", fontsize=12, fontweight="bold")
    ax.set_xticks(np.arange(len(vids))[::5])
    ax.set_xticklabels([vids[i] for i in np.arange(len(vids))[::5]], rotation=45, fontsize=7)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, axis="y")

    # --- 5415: Parameter scatter (bulb_volume vs entrance_fullness) ---
    ax = axes[0, 1]
    bulb_vols = [variants_5415[v]["params"].get("bulb_volume", 1.0) for v in vids]
    entrance_f = [variants_5415[v]["params"].get("entrance_fullness", 1.0) for v in vids]
    scatter = ax.scatter(bulb_vols, entrance_f, c=np.arange(len(vids)), cmap="viridis", s=60, alpha=0.8, edgecolors="black", linewidth=0.5)
    ax.scatter([1.0], [1.0], c="red", s=120, marker="*", label="Baseline", zorder=10)
    ax.set_xlabel("Bulb Volume Factor", fontsize=10)
    ax.set_ylabel("Entrance Fullness Factor", fontsize=10)
    ax.set_title("DTMB 5415 – Parameter Space (bulb_volume vs entrance_fullness)", fontsize=11, fontweight="bold")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    plt.colorbar(scatter, ax=ax, label="Variant index")

    # --- KCS: Displacement distribution ---
    ax = axes[1, 0]
    vids_kcs = sorted(variants_kcs.keys())
    displacements = [variants_kcs[v]["displacement"] for v in vids_kcs]
    ax.bar(np.arange(len(vids_kcs)), displacements, color="coral", alpha=0.8, width=0.8)
    ax.axhline(y=displacements[0], color="red", linestyle="--", linewidth=1.5, label=f"Baseline: {displacements[0]:.1f} t")
    ax.set_xlabel("Variant", fontsize=10)
    ax.set_ylabel("Displacement [t]", fontsize=10)
    ax.set_title("KCS – Displacement Distribution", fontsize=12, fontweight="bold")
    ax.set_xticks(np.arange(len(vids_kcs))[::5])
    ax.set_xticklabels([vids_kcs[i] for i in np.arange(len(vids_kcs))[::5]], rotation=45, fontsize=7)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, axis="y")

    # --- KCS: Cb vs displacement ---
    ax = axes[1, 1]
    cbs = [variants_kcs[v]["Cb"] for v in vids_kcs]
    ax.scatter(cbs, displacements, c=np.arange(len(vids_kcs)), cmap="plasma", s=60, alpha=0.8, edgecolors="black", linewidth=0.5)
    ax.scatter([cbs[0]], [displacements[0]], c="red", s=120, marker="*", label="Baseline", zorder=10)
    ax.set_xlabel("Block Coefficient Cb", fontsize=10)
    ax.set_ylabel("Displacement [t]", fontsize=10)
    ax.set_title("KCS – Cb vs Displacement", fontsize=12, fontweight="bold")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    plt.colorbar(scatter, ax=ax, label="Variant index")

    plt.tight_layout()
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {outpath}")


# ---------------------------------------------------------------------------
# Figure 5: Gate assertion summary heatmap
# ---------------------------------------------------------------------------
def plot_gate_heatmap(variants_5415, variants_kcs, outpath: str):
    """Heatmap of gate assertion results for all variants."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # --- 5415 ---
    ax = axes[0]
    vids = sorted(variants_5415.keys())
    gate_names = ["watertight_manifold", "displacement_constraint", "non_deform_region"]
    gate_matrix = np.zeros((len(gate_names), len(vids)))
    for j, vid in enumerate(vids):
        v = variants_5415[vid]
        gates = v.get("gates", {})
        if isinstance(gates, dict) and "gates" in gates:
            gates = gates["gates"]
        for i, gn in enumerate(gate_names):
            g = gates.get(gn, {})
            passed = g.get("passed", False) if isinstance(g, dict) else False
            gate_matrix[i, j] = 1.0 if passed else 0.0

    im = ax.imshow(gate_matrix, cmap="RdYlGn", aspect="auto", vmin=0, vmax=1)
    ax.set_yticks(range(len(gate_names)))
    ax.set_yticklabels(["Watertight\nManifold", "Displacement\nConstraint", "Non-Deform\nRegion"], fontsize=9)
    ax.set_xticks(range(0, len(vids), 5))
    ax.set_xticklabels([vids[i] for i in range(0, len(vids), 5)], rotation=45, fontsize=7)
    ax.set_title("DTMB 5415 – Gate Assertion Heatmap", fontsize=12, fontweight="bold")
    ax.set_xlabel("Variant", fontsize=10)
    plt.colorbar(im, ax=ax, ticks=[0, 1], label="Pass/Fail")

    # --- KCS ---
    ax = axes[1]
    vids_kcs = sorted(variants_kcs.keys())
    gate_names_kcs = ["watertight_manifold", "displacement_constraint", "non_deformed_region"]
    gate_matrix_kcs = np.zeros((len(gate_names_kcs), len(vids_kcs)))
    for j, vid in enumerate(vids_kcs):
        v = variants_kcs[vid]
        assertions = v.get("assertions", {})
        # Check if assertions is nested
        if not assertions:
            # Try to get from the variant data directly
            assertions = {}
        for i, gn in enumerate(gate_names_kcs):
            passed = assertions.get(gn, False)
            gate_matrix_kcs[i, j] = 1.0 if passed else 0.0

    im = ax.imshow(gate_matrix_kcs, cmap="RdYlGn", aspect="auto", vmin=0, vmax=1)
    ax.set_yticks(range(len(gate_names_kcs)))
    ax.set_yticklabels(["Watertight\nManifold", "Displacement\nConstraint", "Non-Deformed\nRegion"], fontsize=9)
    ax.set_xticks(range(0, len(vids_kcs), 5))
    ax.set_xticklabels([vids_kcs[i] for i in range(0, len(vids_kcs), 5)], rotation=45, fontsize=7)
    ax.set_title("KCS – Gate Assertion Heatmap", fontsize=12, fontweight="bold")
    ax.set_xlabel("Variant", fontsize=10)
    plt.colorbar(im, ax=ax, ticks=[0, 1], label="Pass/Fail")

    plt.tight_layout()
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {outpath}")


# ---------------------------------------------------------------------------
# Figure 6: 3D hull comparison (wireframe)
# ---------------------------------------------------------------------------
def plot_3d_comparison(variants_5415, variants_kcs, outpath: str):
    """3D wireframe comparison of baseline vs deformed variant."""
    fig = plt.figure(figsize=(16, 6))

    # --- 5415: 3D view ---
    ax1 = fig.add_subplot(121, projection="3d")
    v000 = variants_5415["v000"]
    if v000["stl_path"]:
        verts, faces = load_stl_full(v000["stl_path"])
        # Subsample for performance
        step = max(1, len(faces) // 2000)
        faces_sub = faces[::step]
        ax1.plot_trisurf(verts[:, 0], verts[:, 1], verts[:, 2],
                         triangles=faces_sub, color="lightblue", alpha=0.3, edgecolor="none")

    v015 = variants_5415.get("v015")
    if v015 and v015["stl_path"]:
        verts, faces = load_stl_full(v015["stl_path"])
        step = max(1, len(faces) // 2000)
        faces_sub = faces[::step]
        ax1.plot_trisurf(verts[:, 0], verts[:, 1], verts[:, 2],
                         triangles=faces_sub, color="coral", alpha=0.5, edgecolor="none")

    ax1.set_xlabel("x [m]", fontsize=9)
    ax1.set_ylabel("y [m]", fontsize=9)
    ax1.set_zlabel("z [m]", fontsize=9)
    ax1.set_title("DTMB 5415 – 3D Comparison\nBaseline (blue) vs v015 (coral)", fontsize=11, fontweight="bold")
    ax1.view_init(elev=20, azim=-60)

    # --- KCS: 3D view ---
    ax2 = fig.add_subplot(122, projection="3d")
    v000_kcs = variants_kcs["kcs_v000"]
    if v000_kcs["stl_path"]:
        verts, faces = load_stl_full(v000_kcs["stl_path"])
        step = max(1, len(faces) // 2000)
        faces_sub = faces[::step]
        ax2.plot_trisurf(verts[:, 0], verts[:, 1], verts[:, 2],
                         triangles=faces_sub, color="lightblue", alpha=0.3, edgecolor="none")

    v015_kcs = variants_kcs.get("kcs_v015")
    if v015_kcs and v015_kcs["stl_path"]:
        verts, faces = load_stl_full(v015_kcs["stl_path"])
        step = max(1, len(faces) // 2000)
        faces_sub = faces[::step]
        ax2.plot_trisurf(verts[:, 0], verts[:, 1], verts[:, 2],
                         triangles=faces_sub, color="coral", alpha=0.5, edgecolor="none")

    ax2.set_xlabel("x [m]", fontsize=9)
    ax2.set_ylabel("y [m]", fontsize=9)
    ax2.set_zlabel("z [m]", fontsize=9)
    ax2.set_title("KCS – 3D Comparison\nBaseline (blue) vs kcs_v015 (coral)", fontsize=11, fontweight="bold")
    ax2.view_init(elev=20, azim=-60)

    plt.tight_layout()
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {outpath}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    outdir = Path(__file__).parent / "artifacts"
    outdir.mkdir(parents=True, exist_ok=True)

    m542_base = "/root/tp-agent/work/m542/deform5415"
    m546_base = "/root/tp-agent/work/m546/artifacts/kcs_variants"

    print("Loading DTMB 5415 family...")
    summary_5415, variants_5415 = load_5415_family(m542_base)
    print(f"  Loaded {len(variants_5415)} variants")

    print("Loading KCS family...")
    summary_kcs, variants_kcs = load_kcs_family(m546_base)
    print(f"  Loaded {len(variants_kcs)} variants")

    print("\nGenerating figures...")

    print("  [1/6] Side view overlay...")
    plot_overlay_side(variants_5415, variants_kcs, str(outdir / "fig01_side_overlay.png"))

    print("  [2/6] Top view overlay...")
    plot_overlay_top(variants_5415, variants_kcs, str(outdir / "fig02_top_overlay.png"))

    print("  [3/6] Non-deformation zone deviation...")
    max_devs_5415, max_devs_kcs = plot_non_deform_deviation(
        variants_5415, variants_kcs, str(outdir / "fig03_non_deform_deviation.png")
    )

    print("  [4/6] Parameter distribution...")
    plot_param_distribution(variants_5415, variants_kcs, str(outdir / "fig04_param_distribution.png"))

    print("  [5/6] Gate assertion heatmap...")
    plot_gate_heatmap(variants_5415, variants_kcs, str(outdir / "fig05_gate_heatmap.png"))

    print("  [6/6] 3D comparison...")
    plot_3d_comparison(variants_5415, variants_kcs, str(outdir / "fig06_3d_comparison.png"))

    # Summary stats
    print("\n=== Summary ===")
    print(f"DTMB 5415: {summary_5415['total_variants']} variants, "
          f"{summary_5415['passed']} passed, pass_rate={summary_5415['pass_rate']}")
    print(f"KCS: {summary_kcs['total_variants']} variants, "
          f"{summary_kcs['passed_all_assertions']} passed all assertions")

    # Non-deform zone verification
    n_zero_5415 = sum(1 for d in max_devs_5415 if d < 1e-6)
    n_zero_kcs = sum(1 for d in max_devs_kcs if d < 1e-6)
    print(f"Non-deform zone: 5415 {n_zero_5415}/{len(max_devs_5415)} at zero, "
          f"KCS {n_zero_kcs}/{len(max_devs_kcs)} at zero")

    print(f"\nAll figures saved to: {outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
