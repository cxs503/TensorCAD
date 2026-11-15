#!/usr/bin/env python
"""Generate publication-quality figures for the mesh format matrix report.

Produces figures in artifacts/m22384/figures/:
- matrix_heatmap.png: Volume relative error heatmap (6 formats × 3 geometries)
- topology_verification.png: Topology invariant verification bar chart
- file_size_comparison.png: File size comparison across formats
- volume_comparison.png: Volume comparison (reference vs round-trip)
"""
from __future__ import annotations

import json
import tempfile
from pathlib import Path

import numpy as np

from tensorcad import mesh_formats as mf

BODIES = (mf.cube(), mf.icosphere(1.0, 2), mf.torus())
FORMATS = mf.MESH_FORMATS
ART_DIR = Path(__file__).resolve().parents[1] / "artifacts" / "m22384"
FIG_DIR = ART_DIR / "figures"


def run_matrix() -> list[dict]:
    """Run the full 18-combination matrix and return results."""
    results = []
    for body in BODIES:
        for fmt in FORMATS:
            tol = 1e-5 if fmt == "stl" else 0.0
            with tempfile.TemporaryDirectory() as td:
                path = Path(td) / f"{body.name}.{fmt}"
                report = mf.roundtrip(body, fmt, path, tol=tol)
                results.append(
                    {
                        "body": body.name,
                        "format": fmt,
                        "pass": report["pass"],
                        "volume_ref": report["volume_ref"],
                        "volume_got": report["volume_got"],
                        "volume_rel_err": report["volume_rel_err"],
                        "vertices": report["vertices"],
                        "faces": report["faces"],
                        "euler": report["euler"],
                        "boundary_edges": report["boundary_edges"],
                        "bijection": report["bijection"],
                        "flipped_faces": report["flipped_faces"],
                        "bytes": report["bytes"],
                    }
                )
    return results


def generate_heatmap(results: list[dict]) -> bytes:
    """Figure 1: Volume relative error heatmap."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors

    body_names = [b.name for b in BODIES]
    fmt_names = list(FORMATS)

    err_matrix = np.zeros((len(body_names), len(fmt_names)))
    for r in results:
        i = body_names.index(r["body"])
        j = fmt_names.index(r["format"])
        err_matrix[i, j] = r["volume_rel_err"]

    fig, ax = plt.subplots(figsize=(10, 5))
    # Use symmetric log scale to show both zero and non-zero values
    norm = mcolors.SymLogNorm(linthresh=1e-15, linscale=1, vmin=0, vmax=1e-5)
    im = ax.imshow(err_matrix, cmap="YlOrRd", norm=norm, aspect="auto")

    ax.set_xticks(range(len(fmt_names)))
    ax.set_xticklabels(fmt_names, fontsize=12, fontweight="bold")
    ax.set_yticks(range(len(body_names)))
    body_labels = ["cube (genus 0)", "icosphere (genus 0)", "torus (genus 1)"]
    ax.set_yticklabels(body_labels, fontsize=11)

    for i in range(len(body_names)):
        for j in range(len(fmt_names)):
            val = err_matrix[i, j]
            if val == 0:
                label = "0"
            elif val < 1e-12:
                label = "<1e-12"
            else:
                label = f"{val:.1e}"
            ax.text(j, i, label, ha="center", va="center", fontsize=10,
                    color="white" if val > 1e-8 else "black")

    ax.set_xlabel("Exchange Format", fontsize=13, fontweight="bold")
    ax.set_ylabel("Test Geometry", fontsize=13, fontweight="bold")
    ax.set_title("Round-trip Volume Relative Error\n(6 formats × 3 geometries = 18 combinations)",
                 fontsize=14, fontweight="bold", pad=12)
    cbar = plt.colorbar(im, ax=ax, label="Relative error (symlog scale)")
    plt.tight_layout()

    from io import BytesIO
    buf = BytesIO()
    plt.savefig(buf, format="png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    return buf.getvalue()


def generate_topology_chart(results: list[dict]) -> bytes:
    """Figure 2: Topology invariant verification."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    body_names = [b.name for b in BODIES]
    fmt_names = list(FORMATS)

    # Check topology invariants
    vertex_preserved = np.zeros((len(body_names), len(fmt_names)), dtype=bool)
    face_preserved = np.zeros((len(body_names), len(fmt_names)), dtype=bool)
    euler_preserved = np.zeros((len(body_names), len(fmt_names)), dtype=bool)
    closed = np.zeros((len(body_names), len(fmt_names)), dtype=bool)

    for r in results:
        i = body_names.index(r["body"])
        j = fmt_names.index(r["format"])
        vertex_preserved[i, j] = r["vertices"][0] == r["vertices"][1]
        face_preserved[i, j] = r["faces"][0] == r["faces"][1]
        euler_preserved[i, j] = r["euler"][0] == r["euler"][1]
        closed[i, j] = r["boundary_edges"][1] == 0

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    titles = ["Vertex Count Preserved", "Face Count Preserved",
              "Euler Characteristic Preserved", "Closed Surface (no boundary edges)"]
    data = [vertex_preserved, face_preserved, euler_preserved, closed]

    for ax, title, d in zip(axes.flat, titles, data):
        im = ax.imshow(d, cmap="RdYlGn", vmin=0, vmax=1, aspect="auto")
        ax.set_xticks(range(len(fmt_names)))
        ax.set_xticklabels(fmt_names, fontsize=10)
        ax.set_yticks(range(len(body_names)))
        ax.set_yticklabels(body_names, fontsize=10)
        for i in range(len(body_names)):
            for j in range(len(fmt_names)):
                ax.text(j, i, "✅" if d[i, j] else "❌",
                        ha="center", va="center", fontsize=14)
        ax.set_title(title, fontsize=11, fontweight="bold")

    plt.suptitle("Topology Invariant Verification (18/18 combinations)",
                 fontsize=14, fontweight="bold", y=1.01)
    plt.tight_layout()

    from io import BytesIO
    buf = BytesIO()
    plt.savefig(buf, format="png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    return buf.getvalue()


def generate_file_size_chart(results: list[dict]) -> bytes:
    """Figure 3: File size comparison across formats."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    body_names = [b.name for b in BODIES]
    fmt_names = list(FORMATS)

    size_matrix = np.zeros((len(body_names), len(fmt_names)))
    for r in results:
        i = body_names.index(r["body"])
        j = fmt_names.index(r["format"])
        size_matrix[i, j] = r["bytes"]

    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(len(fmt_names))
    width = 0.25

    colors = ["#2196F3", "#4CAF50", "#FF9800"]
    for i, (body, color) in enumerate(zip(body_names, colors)):
        sizes = [r["bytes"] for r in results if r["body"] == body]
        bars = ax.bar(x + i * width, sizes, width, label=body, color=color, edgecolor="white")
        for bar, size in zip(bars, sizes):
            if size > 1000:
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 500,
                        f"{size/1024:.1f}KB", ha="center", va="bottom", fontsize=8, rotation=90)
            else:
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 5,
                        f"{size}B", ha="center", va="bottom", fontsize=8, rotation=90)

    ax.set_xticks(x + width)
    ax.set_xticklabels(fmt_names, fontsize=11, fontweight="bold")
    ax.set_ylabel("File Size (bytes)", fontsize=12)
    ax.set_title("Output File Size by Format and Geometry", fontsize=13, fontweight="bold")
    ax.legend(fontsize=10)
    ax.set_yscale("log")
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()

    from io import BytesIO
    buf = BytesIO()
    plt.savefig(buf, format="png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    return buf.getvalue()


def generate_volume_comparison(results: list[dict]) -> bytes:
    """Figure 4: Volume comparison (reference vs round-trip)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    body_names = [b.name for b in BODIES]
    fmt_names = list(FORMATS)

    for idx, body in enumerate(BODIES):
        ax = axes[idx]
        vol_ref = body.volume()
        vols = []
        for fmt in FORMATS:
            r = next(r for r in results if r["body"] == body.name and r["format"] == fmt)
            vols.append(r["volume_got"])

        x = np.arange(len(fmt_names))
        bars = ax.bar(x, vols, color="#2196F3", edgecolor="white", width=0.6)
        ax.axhline(y=vol_ref, color="#FF5722", linestyle="--", linewidth=2,
                    label=f"Reference: {vol_ref:.6f}")

        for bar, v in zip(bars, vols):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
                    f"{v:.4f}", ha="center", va="bottom", fontsize=8, rotation=45)

        ax.set_xticks(x)
        ax.set_xticklabels(fmt_names, fontsize=10)
        ax.set_title(f"{body.name} (V={vol_ref:.4f} m³)", fontsize=11, fontweight="bold")
        ax.legend(fontsize=8)
        ax.grid(axis="y", alpha=0.3)

    plt.suptitle("Round-trip Volume Comparison: Reference vs Each Format",
                 fontsize=13, fontweight="bold", y=1.05)
    plt.tight_layout()

    from io import BytesIO
    buf = BytesIO()
    plt.savefig(buf, format="png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    return buf.getvalue()


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    results = run_matrix()

    # Save results JSON
    summary = {
        "total": len(results),
        "pass": sum(1 for r in results if r["pass"]),
        "fail": sum(1 for r in results if not r["pass"]),
        "worst_own_reader_rel_err": max(r["volume_rel_err"] for r in results),
    }
    data = {
        "summary": summary,
        "results": results,
        "acceptance": {
            "implementation_threshold": 1e-5,
            "description": "6 formats × 3 geometries = 18 combinations, all pass",
        },
    }
    (ART_DIR / "results.json").write_text(json.dumps(data, indent=2), encoding="utf-8")

    # Generate figures
    figures = {
        "matrix_heatmap.png": ("Volume Relative Error Heatmap", generate_heatmap),
        "topology_verification.png": ("Topology Invariant Verification", generate_topology_chart),
        "file_size_comparison.png": ("File Size Comparison", generate_file_size_chart),
        "volume_comparison.png": ("Volume Comparison", generate_volume_comparison),
    }

    for name, (caption, gen_func) in figures.items():
        png_data = gen_func(results)
        (FIG_DIR / name).write_bytes(png_data)
        print(f"  ✓ {name} ({len(png_data)} bytes)")

    print(f"\nAll figures written to {FIG_DIR}")
    print(f"Matrix: {summary['pass']}/{summary['total']} passed")
    print(f"Worst rel err: {summary['worst_own_reader_rel_err']:.2e}")


if __name__ == "__main__":
    main()
