#!/usr/bin/env python3
"""DTMB 5415 parametric reconstruction example.

Generates:
  - examples/dtmb5415/results.json   (numeric results, programmatic)
  - examples/dtmb5415/fig_reference_hull.png   (3-view of reference hull)
  - examples/dtmb5415/fig_gradient_parity.png  (analytic vs FD gradient)
  - examples/dtmb5415/fig_quality_summary.png  (quality metrics across configs)

Run from the repo root:
    PYTHONPATH=. python examples/dtmb5415/run_example.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent


def main() -> int:
    from tensorcad.dtmb5415 import (DATASET, DRAFT_M, DTMB5415Space,
                                    half_breadth_mask, project_half_breadth)

    space = DTMB5415Space()
    results: dict = {}

    # ------------------------------------------------------------------
    # 1. Reference hull (v = 0)
    # ------------------------------------------------------------------
    v0 = space.zero()
    verts = space.vertices(v0)
    faces = space.faces()
    q = space.quality(v0)
    hydro = space.hydrostatics(v0)
    vol = space.volume(v0)
    grad = space.volume_gradient(v0)

    results["reference_hull"] = {
        "design_vector": "zeros (reference)",
        "n_vertices": int(q["n_vertices"]),
        "n_faces": int(q["n_faces"]),
        "shell_volume_m3": round(float(vol), 4),
        "displacement_volume_m3": round(float(hydro["displacement_volume_m3"]), 4),
        "waterplane_area_m2": round(float(hydro["waterplane_area_m2"]), 4),
        "lcb_m": round(float(hydro["lcb_m"]), 4),
        "cb": round(float(hydro["cb"]), 4),
        "cwp": round(float(hydro["cwp"]), 4),
        "min_angle_deg": round(float(q["min_angle_deg"]), 1),
        "edge_median_m": round(float(q["edge_median_m"]), 4),
        "edge_max_m": round(float(q["edge_max_m"]), 4),
        "half_breadth_raw_min_m": round(float(q["half_breadth_raw_min_m"]), 6),
        "half_breadth_min_m": round(float(q["half_breadth_min_m"]), 6),
        "half_breadth_clamped_vertices": int(q["half_breadth_clamped_vertices"]),
        "volume_gradient_norm": round(float(np.linalg.norm(grad)), 4),
    }

    # ------------------------------------------------------------------
    # 2. Gradient parity check (all 27 modes)
    # ------------------------------------------------------------------
    h = 1e-5
    fd_grads = []
    ana_grads = space.volume_gradient(v0)
    for mode in range(27):
        e = np.zeros(27)
        e[mode] = h
        fd = (space.volume(e) - space.volume(-e)) / (2 * h)
        fd_grads.append(float(fd))

    rel_errors = []
    for i in range(27):
        fd = fd_grads[i]
        ana = ana_grads[i]
        if abs(fd) > 1e-10:
            rel_errors.append(abs(ana - fd) / abs(fd))
        else:
            rel_errors.append(abs(ana - fd))

    results["gradient_parity"] = {
        "max_relative_error": round(float(max(rel_errors)), 10),
        "mean_relative_error": round(float(np.mean(rel_errors)), 10),
        "all_modes_pass_1pct": all(e < 0.01 for e in rel_errors),
    }

    # ------------------------------------------------------------------
    # 3. Quality across sampled configurations
    # ------------------------------------------------------------------
    fixture = np.load(Path(__file__).resolve().parents[2] / "docs" / "data" / "dtmb5415" / "holdout_fixture.npz")
    fixture_vectors = [fixture["design_vector"][i].astype(np.float64) for i in range(len(fixture["design_vector"]))]

    rng = np.random.default_rng(20700)
    n_random = max(0, 20 - len(fixture_vectors))
    random_vectors = [rng.uniform(-1.0, 1.0, size=27) for _ in range(n_random)]
    all_vectors = fixture_vectors + random_vectors

    quality_summary = []
    for i, v in enumerate(all_vectors):
        q = space.quality(v)
        quality_summary.append({
            "config": i,
            "n_vertices": int(q["n_vertices"]),
            "n_faces": int(q["n_faces"]),
            "min_angle_deg": round(float(q["min_angle_deg"]), 1),
            "half_breadth_raw_min_m": round(float(q["half_breadth_raw_min_m"]), 6),
            "half_breadth_min_m": round(float(q["half_breadth_min_m"]), 6),
            "half_breadth_clamped_vertices": int(q["half_breadth_clamped_vertices"]),
            "edge_max_m": round(float(q["edge_max_m"]), 4),
        })

    results["quality_summary"] = {
        "n_configs": len(quality_summary),
        "all_pass": all(
            q["min_angle_deg"] >= 20.0 and
            q["half_breadth_min_m"] >= -1e-6 and
            q["edge_max_m"] <= 8.0
            for q in quality_summary
        ),
        "min_angle_range": [
            round(min(q["min_angle_deg"] for q in quality_summary), 1),
            round(max(q["min_angle_deg"] for q in quality_summary), 1),
        ],
        "max_clamped_vertices": max(q["half_breadth_clamped_vertices"] for q in quality_summary),
    }

    # ------------------------------------------------------------------
    # 4. Model accuracy (holdout)
    # ------------------------------------------------------------------
    rms_errors = []
    for v, true in zip(fixture["design_vector"], fixture["true_cloud"]):
        pred = space.cloud(v.astype(np.float64))
        e = np.linalg.norm(pred - true, axis=1)
        rms_errors.append(float(np.sqrt((e ** 2).mean())))

    results["model_accuracy"] = {
        "n_holdout": len(rms_errors),
        "rms_mean_m": round(float(np.mean(rms_errors)), 4),
        "rms_median_m": round(float(np.median(rms_errors)), 4),
        "rms_max_m": round(float(max(rms_errors)), 4),
        "rms_pct_beam_mean": round(float(np.mean(rms_errors) / space.beam_full_m * 100), 2),
    }

    # ------------------------------------------------------------------
    # Write results
    # ------------------------------------------------------------------
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "results.json").write_text(json.dumps(results, indent=2))
    print(f"Results written to {OUT / 'results.json'}")

    # ------------------------------------------------------------------
    # Generate figures
    # ------------------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    except ImportError:
        print("matplotlib not available, skipping figures")
        return 0

    # Figure 1: Reference hull 3-view
    fig = plt.figure(figsize=(14, 10))

    ax1 = fig.add_subplot(221, projection="3d")
    ax1.plot_trisurf(verts[:, 0], verts[:, 1], verts[:, 2],
                     triangles=faces, cmap="Blues", alpha=0.8, edgecolor="none")
    ax1.set_xlabel("x [m]")
    ax1.set_ylabel("y [m]")
    ax1.set_zlabel("z [m]")
    ax1.set_title("Reference Hull (v=0) — 3D")
    ax1.set_box_aspect([1, 0.15, 0.1])

    ax2 = fig.add_subplot(222)
    # Body plan: y-z at several x stations
    x_stations = np.linspace(verts[:, 0].min(), verts[:, 0].max(), 9)
    for xs in x_stations:
        mask = np.abs(verts[:, 0] - xs) < 1.5
        if mask.sum() > 10:
            ax2.plot(verts[mask, 1], verts[mask, 2], "o", markersize=1, alpha=0.5)
    ax2.set_xlabel("y [m]")
    ax2.set_ylabel("z [m]")
    ax2.set_title("Body Plan")
    ax2.set_aspect("equal")
    ax2.axhline(y=0, color="k", linewidth=0.5)

    ax3 = fig.add_subplot(223)
    # Profile view: x-z at y=0 (centreline) and max beam
    centreline = np.abs(verts[:, 1]) < 0.1
    ax3.plot(verts[centreline, 0], verts[centreline, 2], "b-", linewidth=0.5, label="Centreline")
    ax3.set_xlabel("x [m]")
    ax3.set_ylabel("z [m]")
    ax3.set_title("Profile View")
    ax3.legend()

    ax4 = fig.add_subplot(224)
    # Plan view: x-y at waterline
    wl_mask = np.abs(verts[:, 2] - DRAFT_M) < 0.5
    ax4.plot(verts[wl_mask, 0], verts[wl_mask, 1], "b.", markersize=1, alpha=0.3, label="Starboard")
    ax4.plot(verts[wl_mask, 0], -verts[wl_mask, 1], "r.", markersize=1, alpha=0.3, label="Port")
    ax4.set_xlabel("x [m]")
    ax4.set_ylabel("y [m]")
    ax4.set_title("Plan View (near waterline)")
    ax4.set_aspect("equal")
    ax4.legend()

    plt.tight_layout()
    fig.savefig(OUT / "fig_reference_hull.png", dpi=150)
    plt.close(fig)
    print(f"Figure saved: {OUT / 'fig_reference_hull.png'}")

    # Figure 2: Gradient parity
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(range(27), ana_grads, "bo-", label="Analytic", markersize=4)
    ax.plot(range(27), fd_grads, "r^--", label="Finite Difference", markersize=4)
    ax.set_xlabel("Design variable index")
    ax.set_ylabel("dV/dv [m³]")
    ax.set_title("Volume Gradient: Analytic vs Finite Difference")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fig.savefig(OUT / "fig_gradient_parity.png", dpi=150)
    plt.close(fig)
    print(f"Figure saved: {OUT / 'fig_gradient_parity.png'}")

    # Figure 3: Quality summary
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))

    angles = [q["min_angle_deg"] for q in quality_summary]
    axes[0].bar(range(len(angles)), angles, color="steelblue")
    axes[0].axhline(y=20, color="r", linestyle="--", label="Gate (20°)")
    axes[0].set_xlabel("Config")
    axes[0].set_ylabel("Min angle [°]")
    axes[0].set_title("Minimum Angle")
    axes[0].legend()

    clamped = [q["half_breadth_clamped_vertices"] for q in quality_summary]
    axes[1].bar(range(len(clamped)), clamped, color="coral")
    axes[1].set_xlabel("Config")
    axes[1].set_ylabel("Clamped vertices")
    axes[1].set_title("Symmetry-Plane Clamped Vertices")

    edges = [q["edge_max_m"] for q in quality_summary]
    axes[2].bar(range(len(edges)), edges, color="seagreen")
    axes[2].axhline(y=8.0, color="r", linestyle="--", label="Gate (8 m)")
    axes[2].set_xlabel("Config")
    axes[2].set_ylabel("Max edge [m]")
    axes[2].set_title("Maximum Edge Length")
    axes[2].legend()

    plt.tight_layout()
    fig.savefig(OUT / "fig_quality_summary.png", dpi=150)
    plt.close(fig)
    print(f"Figure saved: {OUT / 'fig_quality_summary.png'}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
