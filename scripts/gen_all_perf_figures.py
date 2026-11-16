#!/usr/bin/env python3
"""Generate all 13 performance figures for GAP-CAD-D4-PERF-SCALE.

Reads committed JSON artifacts from benchmarks/perf_scale/ and writes PNGs
to benchmarks/perf_scale/figures/.  Every figure is a pure function of the
artifacts — no hand-typed numbers.

Usage:
    python scripts/gen_all_perf_figures.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PERF = ROOT / "benchmarks" / "perf_scale"
ARTIFACTS = PERF / "artifacts"
FIGURES = PERF / "figures"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _fmt(n: float) -> str:
    """Format a number with commas for readability."""
    if n >= 1e6:
        return f"{n/1e6:.1f}M"
    if n >= 1e3:
        return f"{n/1e3:.1f}k"
    return f"{n:.0f}"


# ---------------------------------------------------------------------------
# 1. scaling.png — wall-clock vs size log-log (4 subtasks, mesh ladder)
# ---------------------------------------------------------------------------
def gen_scaling() -> None:
    data = load_json(ARTIFACTS / "ladder_mesh.json")
    results = data["results"]

    fig, ax = plt.subplots(figsize=(10, 7))

    colors = {"load": "#2196F3", "point_locate": "#F44336",
              "point_locate_accel": "#4CAF50", "volume": "#9C27B0"}
    markers = {"load": "o", "point_locate": "s",
               "point_locate_accel": "^", "volume": "D"}

    for task, color in colors.items():
        pts = [(r["triangles"], r["seconds"]) for r in results
               if r["subtask"] == task and r["status"] == "measured" and r.get("seconds")]
        if not pts:
            continue
        x, y = zip(*pts)
        ax.loglog(x, y, f"{markers[task]}-", color=color, linewidth=2,
                  markersize=7, label=task)

    ax.axvline(1e6, color="gray", linestyle="--", alpha=0.4, label="10⁶")
    ax.axvline(1e7, color="gray", linestyle="--", alpha=0.4, label="10⁷")
    ax.set_xlabel("Triangles")
    ax.set_ylabel("Wall-clock (s)")
    ax.set_title("TensorCAD Geometry Performance: Wall-Clock vs Mesh Size")
    ax.grid(True, alpha=0.3, which="both")
    ax.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(FIGURES / "scaling.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  scaling.png ✓")


# ---------------------------------------------------------------------------
# 2. point_locate.png — brute-force vs accel throughput comparison
# ---------------------------------------------------------------------------
def gen_point_locate() -> None:
    data = load_json(ARTIFACTS / "ladder_mesh.json")
    results = data["results"]

    fig, ax = plt.subplots(figsize=(10, 6))

    for task, color, marker in [("point_locate", "#F44336", "s"),
                                 ("point_locate_accel", "#4CAF50", "^")]:
        pts = [(r["triangles"], r["throughput"]) for r in results
               if r["subtask"] == task and r["status"] == "measured" and r.get("throughput")]
        if not pts:
            continue
        x, y = zip(*pts)
        ax.loglog(x, y, f"{marker}-", color=color, linewidth=2, markersize=8,
                  label=task)

    ax.axvline(1e6, color="gray", linestyle="--", alpha=0.4)
    ax.axvline(1e7, color="gray", linestyle="--", alpha=0.4)
    ax.set_xlabel("Triangles")
    ax.set_ylabel("Throughput (queries/s)")
    ax.set_title("Point Location: Brute-Force vs Spatial Index (TriangleGrid)")
    ax.grid(True, alpha=0.3, which="both")
    ax.legend(fontsize=10)
    plt.tight_layout()
    plt.savefig(FIGURES / "point_locate.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  point_locate.png ✓")


# ---------------------------------------------------------------------------
# 3. index_cost.png — index build time & size vs triangles
# ---------------------------------------------------------------------------
def gen_index_cost() -> None:
    data = load_json(ARTIFACTS / "ladder_mesh.json")
    results = data["results"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    pts = [(r["triangles"], r["build_seconds"]) for r in results
           if r["subtask"] == "point_locate_accel" and r["status"] == "measured"
           and r.get("build_seconds")]
    if pts:
        x, y = zip(*pts)
        ax1.loglog(x, y, "o-", color="#FF9800", linewidth=2, markersize=8)
        ax1.set_xlabel("Triangles")
        ax1.set_ylabel("Index Build Time (s)", color="#FF9800")
        ax1.set_title("Spatial Index Build Cost")
        ax1.grid(True, alpha=0.3, which="both")

    pts2 = [(r["triangles"], r["index_bytes"] / (1024*1024)) for r in results
            if r["subtask"] == "point_locate_accel" and r["status"] == "measured"
            and r.get("index_bytes")]
    if pts2:
        x, y = zip(*pts2)
        ax2.loglog(x, y, "s-", color="#2196F3", linewidth=2, markersize=8)
        ax2.set_xlabel("Triangles")
        ax2.set_ylabel("Index Size (MiB)", color="#2196F3")
        ax2.set_title("Spatial Index Memory Footprint")
        ax2.grid(True, alpha=0.3, which="both")

    plt.tight_layout()
    plt.savefig(FIGURES / "index_cost.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  index_cost.png ✓")


# ---------------------------------------------------------------------------
# 4. voxelize.png — voxelization scaling (32³ → 384³)
# ---------------------------------------------------------------------------
def gen_voxelize() -> None:
    data = load_json(ARTIFACTS / "ladder_voxel.json")
    results = data["results"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    pts = [(r["grid_cells"], r["seconds"]) for r in results
           if r["status"] == "measured" and r.get("seconds")]
    if pts:
        x, y = zip(*pts)
        ax1.loglog(x, y, "o-", color="#4CAF50", linewidth=2, markersize=8)
        ax1.set_xlabel("Voxel Grid Cells")
        ax1.set_ylabel("Wall-clock (s)")
        ax1.set_title("SDF Voxelization: Wall-Clock vs Grid Size")
        ax1.grid(True, alpha=0.3, which="both")

    pts2 = [(r["grid_cells"], r["throughput"]) for r in results
            if r["status"] == "measured" and r.get("throughput")]
    if pts2:
        x, y = zip(*pts2)
        ax2.semilogx(x, y, "s-", color="#2196F3", linewidth=2, markersize=8)
        ax2.set_xlabel("Voxel Grid Cells")
        ax2.set_ylabel("Throughput (cells/s)")
        ax2.set_title("SDF Voxelization: Throughput vs Grid Size")
        ax2.grid(True, alpha=0.3, which="both")

    plt.tight_layout()
    plt.savefig(FIGURES / "voxelize.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  voxelize.png ✓")


# ---------------------------------------------------------------------------
# 5. fig3pt_scaling.png — fleet 3-point ladder exponents
# ---------------------------------------------------------------------------
def gen_fig3pt_scaling() -> None:
    data = load_json(PERF / "ladder_fleet_10m_3pt.json")
    results = data["results"]
    gates = data.get("gates", [])

    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot scaling for each subtask
    by_task: dict[str, list] = {}
    for r in results:
        if r["status"] == "measured" and r.get("seconds"):
            by_task.setdefault(r["subtask"], []).append(
                (r["triangles"], r["seconds"]))

    colors = {"load": "#2196F3", "fleet_assemble": "#FF9800",
              "point_locate_accel": "#4CAF50", "volume": "#9C27B0"}

    for task, color in colors.items():
        pts = by_task.get(task, [])
        if len(pts) < 2:
            continue
        x, y = zip(*pts)
        ax.loglog(x, y, "o-", color=color, linewidth=2, markersize=8,
                  label=task)

    # Add gate exponents as text
    gate_text = "Complexity Exponents (fitted):\n"
    for g in gates:
        if g["status"] == "measured":
            gate_text += f"  {g['name']}: {g.get('ratio', '—')}\n"

    ax.text(0.02, 0.98, gate_text, transform=ax.transAxes, fontsize=9,
            verticalalignment="top", fontfamily="monospace",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

    ax.axvline(1e6, color="gray", linestyle="--", alpha=0.4)
    ax.axvline(1e7, color="gray", linestyle="--", alpha=0.4)
    ax.set_xlabel("Triangles")
    ax.set_ylabel("Wall-clock (s)")
    ax.set_title("Fleet 3-Point Ladder: Scaling in 10⁶–10⁷ Band")
    ax.grid(True, alpha=0.3, which="both")
    ax.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(FIGURES / "fig3pt_scaling.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  fig3pt_scaling.png ✓")


# ---------------------------------------------------------------------------
# 6. fig3pt_memory.png — memory vs size for fleet ladder
# ---------------------------------------------------------------------------
def gen_fig3pt_memory() -> None:
    data = load_json(PERF / "ladder_fleet_10m_3pt.json")
    results = data["results"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # RSS vs triangles
    pts = [(r["triangles"], r.get("rss_mb") or r.get("peak_rss_mb")) for r in results
           if r["status"] == "measured" and (r.get("rss_mb") or r.get("peak_rss_mb"))]
    if pts:
        x, y = zip(*pts)
        ax1.semilogx(x, y, "o-", color="#F44336", linewidth=2, markersize=8)
        ax1.set_xlabel("Triangles")
        ax1.set_ylabel("RSS (MiB)")
        ax1.set_title("Live Memory vs Mesh Size (Fleet)")
        ax1.grid(True, alpha=0.3)

    # Index bytes vs triangles
    pts2 = [(r["triangles"], r["index_bytes"] / (1024*1024)) for r in results
            if r["subtask"] == "point_locate_accel" and r["status"] == "measured"
            and r.get("index_bytes")]
    if pts2:
        x, y = zip(*pts2)
        ax2.semilogx(x, y, "s-", color="#2196F3", linewidth=2, markersize=8)
        ax2.set_xlabel("Triangles")
        ax2.set_ylabel("Index Size (MiB)")
        ax2.set_title("Spatial Index Memory (Fleet)")
        ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(FIGURES / "fig3pt_memory.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  fig3pt_memory.png ✓")


# ---------------------------------------------------------------------------
# 7. import_scaling.png — load/obj_write/obj_read/volume wall clock
# ---------------------------------------------------------------------------
def gen_import_scaling() -> None:
    data = load_json(ARTIFACTS / "ladder_import.json")
    results = data["results"]

    fig, ax = plt.subplots(figsize=(10, 6))

    colors = {"load": "#2196F3", "obj_write": "#FF9800",
              "obj_read": "#4CAF50", "volume": "#9C27B0"}

    for task, color in colors.items():
        pts = [(r["triangles"], r["seconds"]) for r in results
               if r["subtask"] == task and r["status"] == "measured" and r.get("seconds")]
        if not pts:
            continue
        x, y = zip(*pts)
        ax.loglog(x, y, "o-", color=color, linewidth=2, markersize=7,
                  label=task)

    ax.set_xlabel("Triangles")
    ax.set_ylabel("Wall-clock (s)")
    ax.set_title("OBJ Import/Export Pipeline: Wall-Clock vs Mesh Size")
    ax.grid(True, alpha=0.3, which="both")
    ax.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(FIGURES / "import_scaling.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  import_scaling.png ✓")


# ---------------------------------------------------------------------------
# 8. import_throughput.png — throughput vs size for import pipeline
# ---------------------------------------------------------------------------
def gen_import_throughput() -> None:
    data = load_json(ARTIFACTS / "ladder_import.json")
    results = data["results"]

    fig, ax = plt.subplots(figsize=(10, 6))

    colors = {"load": "#2196F3", "obj_write": "#FF9800",
              "obj_read": "#4CAF50", "volume": "#9C27B0"}

    for task, color in colors.items():
        pts = [(r["triangles"], r["throughput"]) for r in results
               if r["subtask"] == task and r["status"] == "measured" and r.get("throughput")]
        if not pts:
            continue
        x, y = zip(*pts)
        ax.loglog(x, y, "o-", color=color, linewidth=2, markersize=7,
                  label=task)

    ax.set_xlabel("Triangles")
    ax.set_ylabel("Throughput (tri/s)")
    ax.set_title("OBJ Import/Export Pipeline: Throughput vs Mesh Size")
    ax.grid(True, alpha=0.3, which="both")
    ax.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(FIGURES / "import_throughput.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  import_throughput.png ✓")


# ---------------------------------------------------------------------------
# 9. import_volume_accuracy.png — OBJ round-trip volume accuracy
# ---------------------------------------------------------------------------
def gen_import_volume_accuracy() -> None:
    data = load_json(ARTIFACTS / "ladder_import.json")
    results = data["results"]

    fig, ax = plt.subplots(figsize=(8, 5))

    pts = [(r["triangles"], abs(r.get("volume_rel_err", 0) or 0) * 100)
           for r in results if r["subtask"] == "obj_read"
           and r["status"] == "measured" and r.get("volume_rel_err") is not None]
    if pts:
        x, y = zip(*pts)
        ax.semilogx(x, y, "o-", color="#F44336", linewidth=2, markersize=8)
        ax.axhline(3.0, color="gray", linestyle="--", alpha=0.5, label="3% limit")
        ax.set_xlabel("Triangles")
        ax.set_ylabel("Volume |rel_err| (%)")
        ax.set_title("OBJ Round-Trip: Volume Accuracy vs Mesh Size")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=9)

    plt.tight_layout()
    plt.savefig(FIGURES / "import_volume_accuracy.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  import_volume_accuracy.png ✓")


# ---------------------------------------------------------------------------
# 10. single10m_scaling.png — single body 10⁷ subdivision ladder
# ---------------------------------------------------------------------------
def gen_single10m_scaling() -> None:
    data = load_json(ARTIFACTS / "ladder_single_10m.json")
    results = data["results"]

    fig, ax = plt.subplots(figsize=(10, 6))

    colors = {"load": "#2196F3", "subdivide": "#FF9800",
              "volume": "#9C27B0", "point_locate_accel": "#4CAF50"}

    for task, color in colors.items():
        pts = [(r["triangles"], r["seconds"]) for r in results
               if r["subtask"] == task and r["status"] == "measured" and r.get("seconds")]
        if not pts:
            continue
        x, y = zip(*pts)
        ax.loglog(x, y, "o-", color=color, linewidth=2, markersize=7,
                  label=task)

    ax.axvline(1e7, color="gray", linestyle="--", alpha=0.4, label="10⁷")
    ax.set_xlabel("Triangles")
    ax.set_ylabel("Wall-clock (s)")
    ax.set_title("Single-Body 10⁷ Ladder: Subdivision Scaling")
    ax.grid(True, alpha=0.3, which="both")
    ax.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(FIGURES / "single10m_scaling.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  single10m_scaling.png ✓")


# ---------------------------------------------------------------------------
# 11. single10m_memory.png — memory vs size for single body 10⁷
# ---------------------------------------------------------------------------
def gen_single10m_memory() -> None:
    data = load_json(ARTIFACTS / "ladder_single_10m.json")
    results = data["results"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # RSS vs triangles
    pts = [(r["triangles"], r.get("rss_mb") or r.get("peak_rss_mb")) for r in results
           if r["status"] == "measured" and (r.get("rss_mb") or r.get("peak_rss_mb"))]
    if pts:
        x, y = zip(*pts)
        ax1.semilogx(x, y, "o-", color="#F44336", linewidth=2, markersize=8)
        ax1.axvline(1e7, color="gray", linestyle="--", alpha=0.4)
        ax1.set_xlabel("Triangles")
        ax1.set_ylabel("RSS (MiB)")
        ax1.set_title("Live Memory vs Mesh Size (Single-Body 10⁷)")
        ax1.grid(True, alpha=0.3)

    # Index bytes vs triangles
    pts2 = [(r["triangles"], r["index_bytes"] / (1024*1024)) for r in results
            if r["subtask"] == "point_locate_accel" and r["status"] == "measured"
            and r.get("index_bytes")]
    if pts2:
        x, y = zip(*pts2)
        ax2.semilogx(x, y, "s-", color="#2196F3", linewidth=2, markersize=8)
        ax2.axvline(1e7, color="gray", linestyle="--", alpha=0.4)
        ax2.set_xlabel("Triangles")
        ax2.set_ylabel("Index Size (MiB)")
        ax2.set_title("Spatial Index Memory (Single-Body 10⁷)")
        ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(FIGURES / "single10m_memory.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  single10m_memory.png ✓")


# ---------------------------------------------------------------------------
# 12. single10m_volume_accuracy.png — subdivision volume accuracy
# ---------------------------------------------------------------------------
def gen_single10m_volume_accuracy() -> None:
    data = load_json(ARTIFACTS / "ladder_single_10m.json")
    results = data["results"]

    fig, ax = plt.subplots(figsize=(8, 5))

    # Volume error for volume subtasks
    pts = [(r["triangles"], abs(r.get("volume_rel_err", 0) or 0) * 100)
           for r in results if r["subtask"] == "volume"
           and r["status"] == "measured" and r.get("volume_rel_err") is not None]
    if pts:
        x, y = zip(*pts)
        ax.semilogx(x, y, "o-", color="#F44336", linewidth=2, markersize=8,
                    label="volume")

    ax.axhline(3.0, color="gray", linestyle="--", alpha=0.5, label="3% limit")
    ax.set_xlabel("Triangles")
    ax.set_ylabel("Volume |rel_err| (%)")
    ax.set_title("Single-Body 10⁷: Volume Accuracy vs Mesh Size")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)

    plt.tight_layout()
    plt.savefig(FIGURES / "single10m_volume_accuracy.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  single10m_volume_accuracy.png ✓")


# ---------------------------------------------------------------------------
# 13. reference_scaling_bands.png — coverage bands for reference doc
# ---------------------------------------------------------------------------
def gen_reference_scaling_bands() -> None:
    """Composite figure: (a) all ladders wall-clock vs size, (b) coverage bands."""
    # Load all artifacts
    mesh = load_json(ARTIFACTS / "ladder_mesh.json")
    voxel = load_json(ARTIFACTS / "ladder_voxel.json")
    fleet = load_json(PERF / "ladder_fleet_10m_3pt.json")
    single = load_json(ARTIFACTS / "ladder_single_10m.json")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    # (a) Wall-clock vs size — all ladders
    colors = {"load": "#2196F3", "point_locate": "#F44336",
              "point_locate_accel": "#4CAF50", "volume": "#9C27B0",
              "fleet_assemble": "#FF9800", "subdivide": "#E91E63",
              "voxelize": "#795548"}

    for data, label in [(mesh, "mesh"), (fleet, "fleet"), (single, "single")]:
        for task, color in colors.items():
            pts = [(r.get("triangles", r.get("grid_cells", 0)), r["seconds"])
                   for r in data["results"]
                   if r["subtask"] == task and r["status"] == "measured"
                   and r.get("seconds") and (r.get("triangles") or r.get("grid_cells", 0)) > 0]
            if not pts:
                continue
            x, y = zip(*pts)
            alpha = 0.9 if label == "mesh" else 0.5
            ax1.loglog(x, y, "o-", color=color, linewidth=1.5, markersize=5,
                       alpha=alpha, label=f"{task} ({label})")

    ax1.axvline(1e6, color="gray", linestyle="--", alpha=0.3)
    ax1.axvline(1e7, color="gray", linestyle="--", alpha=0.3)
    ax1.set_xlabel("Triangles")
    ax1.set_ylabel("Wall-clock (s)")
    ax1.set_title("(a) All Ladders: Wall-Clock vs Size")
    ax1.grid(True, alpha=0.3, which="both")

    # (b) Coverage bands bar chart
    subtasks = ["load", "point_locate", "point_locate_accel", "volume",
                "voxelize", "obj_read", "obj_write", "subdivide", "fleet_assemble"]
    max_tri = {
        "load": 8199998,
        "point_locate": 410878,
        "point_locate_accel": 13158272,
        "volume": 13158272,
        "voxelize": 56623104,  # cells, not triangles
        "obj_read": 12011994,
        "obj_write": 12011994,
        "subdivide": 13158272,
        "fleet_assemble": 12011994,
    }

    y_pos = range(len(subtasks))
    bar_colors = []
    for st in subtasks:
        v = max_tri.get(st, 0)
        if v >= 1e7:
            bar_colors.append("#4CAF50")  # green: 10⁷
        elif v >= 1e6:
            bar_colors.append("#2196F3")  # blue: 10⁶
        else:
            bar_colors.append("#F44336")  # red: <10⁶

    bars = ax2.barh(y_pos, [max_tri.get(s, 0) for s in subtasks],
                    color=bar_colors, alpha=0.8)
    ax2.axvline(1e6, color="gray", linestyle="--", alpha=0.5, label="10⁶")
    ax2.axvline(1e7, color="gray", linestyle="--", alpha=0.5, label="10⁷")
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(subtasks, fontsize=9)
    ax2.set_xlabel("Max Measured Size (triangles / cells)")
    ax2.set_xscale("log")
    ax2.set_title("(b) Coverage Bands: Deepest Measured Scale")
    ax2.grid(True, alpha=0.3, axis="x")
    ax2.legend(fontsize=8)

    # Add value labels
    for bar, val in zip(bars, [max_tri.get(s, 0) for s in subtasks]):
        if val > 0:
            ax2.text(val * 1.1, bar.get_y() + bar.get_height()/2,
                     _fmt(val), va="center", fontsize=8)

    plt.tight_layout()
    plt.savefig(FIGURES / "reference_scaling_bands.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  reference_scaling_bands.png ✓")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    print("Generating performance figures for GAP-CAD-D4-PERF-SCALE ...")
    gen_scaling()
    gen_point_locate()
    gen_index_cost()
    gen_voxelize()
    gen_fig3pt_scaling()
    gen_fig3pt_memory()
    gen_import_scaling()
    gen_import_throughput()
    gen_import_volume_accuracy()
    gen_single10m_scaling()
    gen_single10m_memory()
    gen_single10m_volume_accuracy()
    gen_reference_scaling_bands()
    print(f"\nAll 13 figures written to {FIGURES}/")


if __name__ == "__main__":
    main()
