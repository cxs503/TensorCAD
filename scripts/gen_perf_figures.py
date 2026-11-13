"""Generate performance scaling figures for the GAP-CAD-D4-PERF-SCALE report."""
from __future__ import annotations

import base64
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load_baseline(path: str | Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def generate_scaling_figure(baseline_path: str | Path, output_png: str | Path) -> str:
    """Generate a log-log scaling figure from baseline data and return base64 PNG."""
    data = load_baseline(baseline_path)
    results = data.get("results", [])

    # Group by subtask
    by_task: dict[str, list[dict]] = {}
    for r in results:
        by_task.setdefault(r["subtask"], []).append(r)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("TensorCAD Geometry Performance Scaling\n(GAP-CAD-D4-PERF-SCALE)", fontsize=14, fontweight="bold")

    task_labels = {
        "load": ("Hull Generation", "Triangles", "triangles/s"),
        "point_locate": ("Ray Casting (Point Location)", "Triangles", "queries/s"),
        "voxelize": ("SDF Voxelization", "Grid Cells", "cells/s"),
        "volume": ("Volume Integration", "Triangles", "triangles/s"),
    }

    for idx, (task, (title, x_label, y_label)) in enumerate(task_labels.items()):
        ax = axes[idx // 2, idx % 2]
        entries = by_task.get(task, [])
        if not entries:
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(task)
            continue

        # X: size metric (triangles or cells)
        if task == "voxelize":
            x = np.array([e["grid_cells"] for e in entries if e.get("grid_cells", 0) > 0])
        else:
            x = np.array([e["triangles"] for e in entries if e["triangles"] > 0])

        # Y: throughput
        y = np.array([e["throughput"] for e in entries if e.get("throughput") is not None])
        t = np.array([e["seconds"] for e in entries if e.get("seconds") is not None])

        if len(x) < 2:
            ax.text(0.5, 0.5, "Insufficient data", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(title)
            continue

        # Plot throughput vs size (log-log)
        ax.loglog(x, y, "o-", linewidth=2, markersize=8, label="Throughput")

        # Plot time vs size (secondary axis)
        ax2 = ax.twinx()
        ax2.loglog(x, t, "s--", color="red", linewidth=1.5, markersize=6, alpha=0.7, label="Wall time")

        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label, color="blue")
        ax2.set_ylabel("Seconds", color="red")
        ax.set_title(title)
        ax.grid(True, alpha=0.3)

        # Combine legends
        lines1, labels1 = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax.legend(lines1 + lines2, labels1 + labels2, loc="best", fontsize=8)

    plt.tight_layout()
    plt.savefig(output_png, dpi=150, bbox_inches="tight")
    plt.close()

    # Return base64
    png_bytes = Path(output_png).read_bytes()
    return base64.b64encode(png_bytes).decode("ascii")


def generate_throughput_comparison(baseline_path: str | Path, output_png: str | Path) -> str:
    """Generate a bar chart comparing throughput across sub-tasks."""
    data = load_baseline(baseline_path)
    results = data.get("results", [])

    by_task: dict[str, list[dict]] = {}
    for r in results:
        by_task.setdefault(r["subtask"], []).append(r)

    fig, ax = plt.subplots(figsize=(12, 6))

    tasks = ["load", "point_locate", "voxelize", "volume"]
    colors = ["#2196F3", "#FF9800", "#4CAF50", "#9C27B0"]

    for task, color in zip(tasks, colors, strict=True):
        entries = by_task.get(task, [])
        if not entries:
            continue
        # Use the largest size entry for comparison
        if task == "voxelize":
            entry = max(entries, key=lambda e: e.get("grid_cells", 0))
            size = entry.get("grid_cells", 0)
        else:
            entry = max(entries, key=lambda e: e["triangles"])
            size = entry["triangles"]
        tp = entry.get("throughput")
        if tp is None:
            continue
        ax.bar(task, tp, color=color, alpha=0.8, label=f"{task} ({size:,})")
        ax.text(task, tp * 1.05, f"{tp:,.0f}", ha="center", va="bottom", fontsize=9, fontweight="bold")

    ax.set_ylabel("Throughput (ops/s, log scale)")
    ax.set_yscale("log")
    ax.set_title("Peak Throughput by Sub-Task (Largest Size Measured)")
    ax.grid(True, alpha=0.3, axis="y")
    ax.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(output_png, dpi=150, bbox_inches="tight")
    plt.close()

    png_bytes = Path(output_png).read_bytes()
    return base64.b64encode(png_bytes).decode("ascii")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", required=True, help="Path to baseline JSON")
    parser.add_argument("--output-dir", default="artifacts/perf_benchmark", help="Output directory for PNGs")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    b64_1 = generate_scaling_figure(args.baseline, out_dir / "scaling.png")
    b64_2 = generate_throughput_comparison(args.baseline, out_dir / "throughput.png")

    print(f"Scaling figure: {out_dir / 'scaling.png'} ({len(b64_1)} bytes base64)")
    print(f"Throughput figure: {out_dir / 'throughput.png'} ({len(b64_2)} bytes base64)")
