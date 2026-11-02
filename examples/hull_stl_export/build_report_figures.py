#!/usr/bin/env python3
"""Regenerate the two report figures from the example artifacts.

Reads ``artifacts/hull_stl_export/{resolution_study,summary,coarse_gate}.json``
(produced by ``run_example.py``) and writes ``convergence.png`` and
``field_resolution.png`` next to them, plus ``report_numbers.json`` with every
number the figures use. Nothing is hand-copied: the tables in the report are read
from these files.

Needs matplotlib (an example-only dependency; the core package stays standard
library). Usage::

    python examples/hull_stl_export/run_example.py --output artifacts/hull_stl_export
    python examples/hull_stl_export/build_report_figures.py --artifacts artifacts/hull_stl_export
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:  # allow running from a checkout without installing
    sys.path.insert(0, str(REPO_ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from tensorcad import HullParameters, generate_hull, read_stl
from tensorcad import mesh as meshlib


def _half_breadth(x: np.ndarray, z: np.ndarray, length: float, beam: float, draft: float) -> np.ndarray:
    """The analytic breadth distribution the hull samples: sin(pi x/L)^0.7 * sqrt((z+T)/T)."""
    taper = np.sin(np.pi * np.clip(x, 0.0, length) / length) ** 0.7
    profile = np.sqrt(np.clip((z + draft) / draft, 0.0, 1.0))
    return 0.5 * beam * taper * profile


def _facet_cloud(mesh) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    xs: list[float] = []
    zs: list[float] = []
    areas: list[float] = []
    for face in mesh.faces:
        xs.append(sum(mesh.vertices[v][0] for v in face) / 3.0)
        zs.append(sum(mesh.vertices[v][2] for v in face) / 3.0)
        areas.append(meshlib.triangle_area(mesh, face))
    return np.array(xs), np.array(zs), np.array(areas)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifacts", default="artifacts/hull_stl_export")
    args = parser.parse_args(argv)
    root = Path(args.artifacts)
    if not root.is_absolute():
        root = REPO_ROOT / root

    rows = json.loads((root / "resolution_study.json").read_text())
    summary = json.loads((root / "summary.json").read_text())
    coarse_gate = json.loads((root / "coarse_gate.json").read_text())
    plan = summary["plan"]
    stl = summary["stl"]

    # ---- Figure 1: measured error against the analytic displacement ----
    figure, axes = plt.subplots(figsize=(7.2, 4.6))
    triangles = [row["triangles"] for row in rows]
    errors = [abs(row["relative_error"]) for row in rows]
    axes.loglog(triangles, errors, "o-", color="#1f4e79", lw=1.6, ms=7, label="ladder grids")
    for row, count, error in zip(rows, triangles, errors, strict=True):
        axes.annotate(f"{row['nx']}x{row['nz']}", (count, error), textcoords="offset points",
                      xytext=(6, 6), fontsize=8, color="#1f4e79")
    axes.loglog([plan["triangles"]], [abs(plan["relative_error"])], "*", ms=16, color="#c00000",
                label=f"resolved {plan['parameters']['nx']}x{plan['parameters']['nz']} "
                      f"(target {100 * summary['target_volume_error']:g} %)")
    axes.loglog([stl["triangles"]], [abs(stl["relative_volume_error_vs_analytic"])], "s", ms=8,
                mfc="none", mec="#2e7d32", mew=1.8, label="binary STL read back")
    axes.axhline(0.02, color="#7f7f7f", ls="--", lw=1.0)
    axes.text(triangles[0], 0.024, "gate default tolerance 2 %", fontsize=8, color="#555555")
    axes.axhline(summary["target_volume_error"], color="#c00000", ls=":", lw=1.2)
    axes.text(triangles[0], summary["target_volume_error"] * 1.08,
              f"requested tolerance {100 * summary['target_volume_error']:g} %",
              fontsize=8, color="#c00000")
    anchor_x, anchor_y = triangles[2], errors[2]
    guide = np.array([triangles[1], triangles[-1]])
    axes.loglog(guide, anchor_y * (guide / anchor_x) ** -2.0, color="#999999", ls="-.", lw=0.9)
    axes.set_xlabel("triangle count")
    axes.set_ylabel("|relative volume error| vs analytic integral")
    axes.set_title("Hull displacement error against the closed-form volume\n"
                   "(analytic reference, independent of the mesh)")
    axes.grid(True, which="both", ls=":", lw=0.5)
    axes.legend(fontsize=8, loc="lower left")
    figure.tight_layout()
    figure.savefig(root / "convergence.png", dpi=160)
    plt.close(figure)

    # ---- Figure 2: the analytic field and the sampled resolution ----
    length = float(plan["parameters"]["length"])
    beam = float(plan["parameters"]["beam"])
    draft = float(plan["parameters"]["draft"])
    coarse_mesh = generate_hull(HullParameters(length, beam, draft, 11, 5)).mesh
    resolved_mesh = read_stl(summary["artifacts"]["stl"])
    grid_x = np.linspace(0.0, length, 240)
    grid_z = np.linspace(-draft, 0.0, 120)
    mesh_x, mesh_z = np.meshgrid(grid_x, grid_z)
    breadth = _half_breadth(mesh_x, mesh_z, length, beam, draft)

    figure, axes = plt.subplots(1, 2, figsize=(11.0, 4.3), sharey=True)
    for axis, mesh, label in (
        (axes[0], coarse_mesh, f"coarse 11x5, {len(coarse_mesh.faces)} facets"),
        (axes[1], resolved_mesh, f"resolved, {len(resolved_mesh.faces)} facets (STL read back)"),
    ):
        contour = axis.contourf(mesh_x, mesh_z, breadth, levels=18, cmap="viridis")
        xs, zs, areas = _facet_cloud(mesh)
        axis.scatter(xs, zs, s=np.clip(np.sqrt(areas) * 55.0, 1.5, 45.0), c=areas,
                     cmap="autumn_r", edgecolors="none", alpha=0.85)
        axis.set_title(label, fontsize=10)
        axis.set_xlabel("x [m] (aft 0 -> forward L)")
        axis.set_xlim(0.0, length)
        axis.set_ylim(-draft, 0.0)
        colorbar = figure.colorbar(contour, ax=axis, pad=0.02)
        colorbar.set_label("half-breadth b(x,z) [m]", fontsize=8)
    axes[0].set_ylabel("z [m] (keel -T -> design waterline 0)")
    figure.suptitle("Analytic breadth field with facet centroids (marker size ~ local facet area)",
                    fontsize=11)
    figure.tight_layout()
    figure.savefig(root / "field_resolution.png", dpi=160)
    plt.close(figure)

    payload = {
        "ladder": [
            {"grid": f"{row['nx']}x{row['nz']}", "triangles": row["triangles"],
             "volume_m3": row["volume_m3"], "error_pct": 100.0 * row["relative_error"]}
            for row in rows
        ],
        "plan": plan,
        "stl": stl,
        "coarse_gate_passed": coarse_gate["passed"],
        "coarse_gate_failures": coarse_gate["failures"],
        "coarse_volume_error_pct": 100.0 * coarse_gate["measurements"]["relative_volume_error"],
        "analytic_volume_m3": summary["analytic_volume_m3"],
        "resolved_mesh_error_pct": 100.0 * summary["resolved_mesh_relative_error"],
    }
    (root / "report_numbers.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    for name in ("convergence.png", "field_resolution.png"):
        print(f"{name}: {(root / name).stat().st_size} bytes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
