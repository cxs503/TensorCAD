#!/usr/bin/env python
"""Generate figures for IN31-C01 mission m5479 report."""
from __future__ import annotations

import base64
import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tensorcad.geometry_audit import (
    audit_fairness,
    audit_watertightness,
    gaussian_curvature,
    inject_bump,
    inject_seam_offset,
    offsets_from_function,
)
from tensorcad.hydrostatics import compute_hydrostatics
from tensorcad.wigley import (
    WigleyParameters,
    analytic_block_coefficient,
    analytic_volume,
    analytic_waterplane_area,
    analytic_wetted_surface_area,
    build_mesh,
    half_breadth,
    reference_parameters,
    wigley_analytic_hydrostatics,
)

ARTIFACTS_DIR = REPO_ROOT / "artifacts" / "m5479_geometry_audit"
ARTIFACTS_DIR.mkdir(parents=True, exist_ok=True)

WIGLEY_PARAMS = WigleyParameters(length=10.0, beam=1.0, draft=0.625, nx=81, nz=21)


def encode_png(path: Path) -> str:
    """Read a PNG and return base64 string."""
    return base64.b64encode(path.read_bytes()).decode("ascii")


def fig_c1_watertightness():
    """Figure C1: Seam gap detection recall."""
    rng = __import__("random").Random(20261106)
    wigley_mesh, _ = build_mesh(WIGLEY_PARAMS)

    fractions = []
    detected_flags = []
    for _ in range(24):
        fraction = 0.001 + 0.019 * rng.random()
        delta = fraction * WIGLEY_PARAMS.length
        x_split = WIGLEY_PARAMS.length * (0.3 + 0.4 * rng.random())
        defected, _ = inject_seam_offset(wigley_mesh, delta=delta, x_split=x_split)
        audit = audit_watertightness(defected, length=WIGLEY_PARAMS.length)
        fractions.append(fraction)
        detected_flags.append(any(d.kind in ("seam_gap", "open_boundary") for d in audit.defects))

    fig, ax = plt.subplots(figsize=(8, 5))
    colors = ["green" if d else "red" for d in detected_flags]
    ax.scatter([f * 100 for f in fractions], [1 if d else 0 for d in detected_flags],
               c=colors, s=60, edgecolors="black", linewidth=0.5, zorder=3)
    ax.axvline(x=0.1, color="orange", linestyle="--", linewidth=1.5, label="Tolerance (0.1 % L)")
    ax.axhline(y=0.5, color="gray", linestyle=":", linewidth=1)
    ax.set_xlabel("Injected gap (% of L)")
    ax.set_ylabel("Detected (1=yes, 0=no)")
    ax.set_title("C1: Watertightness — Seam Gap Detection Recall")
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["No", "Yes"])
    ax.legend()
    ax.grid(True, alpha=0.3)
    path = ARTIFACTS_DIR / "fig_c1_watertightness.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def fig_c2_fairness():
    """Figure C2: Curvature field on clean vs bumped Wigley."""
    wigley_mesh, _ = build_mesh(WIGLEY_PARAMS)
    base_field = gaussian_curvature(wigley_mesh)

    sigma = 3.5 * WIGLEY_PARAMS.length / 80.0
    centre = (7.0, 0.35, -0.35)
    closest = min(wigley_mesh.vertices, key=lambda p: math.dist(p, centre))
    amplitude = 12.0 * base_field.mean_abs_curvature * sigma * sigma
    bumped = inject_bump(wigley_mesh, centre=closest, amplitude=amplitude, sigma=sigma, field=base_field)
    bump_field = gaussian_curvature(bumped)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for ax, field, title in [
        (axes[0], base_field, "Clean Wigley — Gaussian Curvature"),
        (axes[1], bump_field, "Bumped Wigley — Gaussian Curvature"),
    ]:
        values = [field.curvature[i] for i in range(len(field.curvature)) if field.interior[i]]
        sc = ax.scatter(
            [wigley_mesh.vertices[i][0] for i in range(len(wigley_mesh.vertices)) if field.interior[i]],
            [wigley_mesh.vertices[i][2] for i in range(len(wigley_mesh.vertices)) if field.interior[i]],
            c=values, cmap="RdBu_r", s=3, vmin=-5, vmax=5
        )
        ax.set_xlabel("x (m)")
        ax.set_ylabel("z (m)")
        ax.set_title(title)
        ax.set_aspect("equal")
        fig.colorbar(sc, ax=ax, label="K (1/m²)")

    fig.suptitle("C2: Fairness — Gaussian Curvature Field", fontsize=14, y=1.02)
    path = ARTIFACTS_DIR / "fig_c2_fairness.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def fig_c3_c4_hydrostatics():
    """Figure C3/C4: Hydrostatics comparison — mesh vs analytic."""
    wigley_mesh, _ = build_mesh(WIGLEY_PARAMS)
    hydro = compute_hydrostatics(wigley_mesh, waterline_z=0.0)
    ref = wigley_analytic_hydrostatics(WIGLEY_PARAMS)

    names = ["Volume\n(m³)", "Waterplane\n(m²)", "Wetted Area\n(m²)", "Cb", "Cp", "Am\n(m²)"]
    keys = ["displacement_volume", "waterplane_area", "wetted_surface_area",
            "block_coefficient", "prismatic_coefficient", "maximum_section_area"]

    mesh_vals = [getattr(hydro, k) for k in keys]
    ref_vals = [getattr(ref, k) for k in keys]

    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(len(names))
    width = 0.35

    bars1 = ax.bar(x - width/2, mesh_vals, width, label="Mesh (divergence theorem)", color="steelblue")
    bars2 = ax.bar(x + width/2, ref_vals, width, label="Analytic (closed form)", color="coral")

    ax.set_ylabel("Value")
    ax.set_title("C3/C4: Hydrostatics — Mesh vs Analytic Closed Forms")
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.legend()
    ax.grid(True, alpha=0.3, axis="y")

    # Add error percentages on top
    for i, (m, r) in enumerate(zip(mesh_vals, ref_vals)):
        if abs(r) > 1e-12:
            err = abs(m - r) / abs(r) * 100
            ax.text(i, max(m, r) * 1.05, f"{err:.2f}%", ha="center", fontsize=8, color="darkred")

    path = ARTIFACTS_DIR / "fig_c3_c4_hydrostatics.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def main():
    print("Generating figures...")
    paths = {
        "fig_c1_watertightness.png": fig_c1_watertightness(),
        "fig_c2_fairness.png": fig_c2_fairness(),
        "fig_c3_c4_hydrostatics.png": fig_c3_c4_hydrostatics(),
    }
    for name, path in paths.items():
        print(f"  {name}: {path} ({path.stat().st_size} bytes)")
    print("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
