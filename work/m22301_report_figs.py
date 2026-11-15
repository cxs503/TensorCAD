#!/usr/bin/env python3
"""Figures for the m22301 report: prism/pyramid interior-face non-orthogonality.

Every number is read from artifacts/mesh_quality_prism_pyramid/audit.json, which is written
by tools/prism_pyramid_audit.py -- nothing is transcribed by hand.

    python work/m22301_report_figs.py        # writes /tmp/figs/*.png
"""
import base64
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from tensorcad.meshquality.topology import CELL_FACES, non_orthogonality  # noqa: E402

AUDIT = json.loads((ROOT / "artifacts/mesh_quality_prism_pyramid/audit.json").read_text())
OUT = Path("/tmp/figs")
OUT.mkdir(parents=True, exist_ok=True)


def records(kind):
    return sorted((r for r in AUDIT["non_orthogonality"] if r["kind"] == kind),
                  key=lambda r: r["shear"])


def fig_convergence():
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.6))
    markers = {"prism": ("o", "#1f77b4"), "pyramid": ("s", "#d62728")}

    ax = axes[0]
    for kind in ("prism", "pyramid"):
        rows = records(kind)
        shear = np.array([r["shear"] for r in rows])
        closed = np.array([r["closed_form_deg"] for r in rows])
        lib = np.array([r["library_deg"] for r in rows])
        second = np.array([r["second_route_deg"] for r in rows])
        m, c = markers[kind]
        fine = np.linspace(0.0, 3.2, 200)
        ax.plot(fine, np.degrees(np.arctan(fine)), "-", color=c, lw=1.0,
                label=f"{kind}: closed form atan(shear)")
        ax.plot(shear, closed, "^", ms=11, mfc="none", mec=c)
        ax.plot(shear, lib, m, ms=6, color=c, label=f"{kind}: library")
        ax.plot(shear, second, "_", ms=14, mew=1.6, color=c, label=f"{kind}: atan2 route")
    ax.axhline(70.0, ls="--", color="0.4", lw=1.0)
    ax.text(0.05, 71.5, "shipped gate: max non-orthogonality <= 70 deg", fontsize=8, color="0.3")
    ax.set_xlabel("shear s (x += s*z)")
    ax.set_ylabel("max internal-face non-orthogonality, deg")
    ax.set_title("(a) library vs closed form vs atan2 route")
    ax.legend(fontsize=6.5, loc="lower right")
    ax.grid(alpha=0.3)

    ax = axes[1]
    for kind in ("prism", "pyramid"):
        rows = records(kind)
        shear = np.array([r["shear"] for r in rows])
        err_closed = np.maximum([r["library_minus_closed_form"] for r in rows], 1e-18)
        err_second = np.maximum([r["library_minus_second_route"] for r in rows], 1e-18)
        m, c = markers[kind]
        ax.semilogy(shear + 0.01, err_closed, m, color=c, ms=6,
                    label=f"{kind}: |library - closed form|")
        ax.semilogy(shear + 0.03, err_second, m, color=c, ms=6, mfc="none",
                    label=f"{kind}: |library - atan2 route|")
    ax.axhline(1e-6, ls="--", color="0.4", lw=1.0)
    ax.text(0.05, 2e-6, "mission acceptance 1e-6 deg", fontsize=8, color="0.3")
    ax.axhline(1e-15, ls=":", color="0.6", lw=1.0)
    ax.text(0.05, 1.3e-15, "double-precision round-off 1e-15", fontsize=8, color="0.5")
    ax.set_xlabel("shear s (x += s*z)")
    ax.set_ylabel("absolute difference, deg")
    ax.set_title("(b) residual against both independent routes")
    ax.legend(fontsize=6.5, loc="lower right")
    ax.grid(alpha=0.3, which="both")
    ax.set_ylim(1e-18, 1e-2)

    ax = axes[2]
    labels, values, colors = [], [], []
    for kind in ("prism", "pyramid"):
        rows = records(kind)
        labels.append(f"{kind}\nlibrary vs atan2")
        values.append(max(r["library_minus_second_route"] for r in rows))
        colors.append(markers[kind][1])
    bars = ax.bar(labels, values, color=colors, width=0.5)
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, value * 1.6, f"{value:.2e}",
                ha="center", fontsize=9)
    ax.set_yscale("log")
    ax.set_ylabel("worst |library - atan2|, deg")
    ax.set_title("(c) worst residual per kind over the sweep")
    ax.set_ylim(1e-16, 1e-12)
    ax.grid(alpha=0.3, axis="y", which="both")

    fig.suptitle("Prism/pyramid interior-face non-orthogonality: two table rows + a closed "
                 "form, not new arithmetic", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    path = OUT / "fig1_non_orthogonality.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def prism_stack(layers=3, shear=0.0):
    points = []
    for level in range(layers + 1):
        z = float(level)
        points += [(shear * z, 0.0, z), (shear * z + 1.0, 0.0, z), (shear * z, 1.0, z)]
    cells = [(3 * k, 3 * k + 1, 3 * k + 2, 3 * k + 3, 3 * k + 4, 3 * k + 5) for k in range(layers)]
    return points, cells


def pyramid_pair(shear=0.0):
    base = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0)]
    points = [(x + shear * z, y, z) for x, y, z in base + [(0.5, 0.5, 1.0), (0.5, 0.5, -1.0)]]
    return points, [(0, 1, 2, 3, 4), (0, 1, 2, 3, 5)]


def _project(points, cell, face):
    return [points[cell[i]] for i in face]


def _draw_cell(ax, points, cell, face_table, color, alpha=0.28):
    for face in face_table:
        loop = _project(points, cell, face)
        xs = [p[0] for p in loop] + [loop[0][0]]
        ys = [p[2] for p in loop] + [loop[0][2]]
        ax.fill(xs, ys, color=color, alpha=alpha, lw=0.8, ec="0.25")


def fig_geometry():
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.0))
    shear = 1.0

    ax = axes[0]
    points, cells = prism_stack(3, shear)
    for cell in cells:
        _draw_cell(ax, points, cell, CELL_FACES["prism"], "#1f77b4")
    centres = [(sum(points[i][0] for i in cell) / 6.0, sum(points[i][2] for i in cell) / 6.0)
               for cell in cells]
    ax.plot([c[0] for c in centres], [c[1] for c in centres], "k.-", lw=1.0, ms=4,
            label="node centroids (cell centres)")
    for k in range(len(cells) - 1):
        z_face = float(k + 1)
        x0 = shear * z_face
        ax.plot([x0, x0 + 1.0], [z_face, z_face], "-", color="crimson", lw=2.2)
    ax.plot([], [], "-", color="crimson", lw=2.2, label="internal face (measured)")
    faces = non_orthogonality(points, cells, "prism")
    ax.set_title(f"(a) 3 stacked prisms, shear s={shear:g}\n"
                 f"internal={faces.internal_faces}, boundary={faces.boundary_faces}, "
                 f"max={faces.maximum:.6f} deg")
    ax.set_xlabel("x")
    ax.set_ylabel("z")
    ax.set_aspect("equal")
    ax.legend(fontsize=7, loc="upper left")
    ax.grid(alpha=0.3)

    ax = axes[1]
    points, cells = pyramid_pair(shear)
    for index, cell in enumerate(cells):
        _draw_cell(ax, points, cell, CELL_FACES["pyramid"],
                   "#d62728" if index == 0 else "#2ca02c")
    centres = [(sum(points[i][0] for i in cell) / 5.0, sum(points[i][2] for i in cell) / 5.0)
               for cell in cells]
    ax.plot([c[0] for c in centres], [c[1] for c in centres], "k.-", lw=1.0, ms=5,
            label="node centroids (cell centres)")
    ax.plot([0.0, 1.0], [0.0, 0.0], "-", color="crimson", lw=2.2,
            label="internal face (measured)")
    faces = non_orthogonality(points, cells, "pyramid")
    ax.set_title(f"(b) 2 pyramids sharing a base quad, shear s={shear:g}\n"
                 f"internal={faces.internal_faces}, boundary={faces.boundary_faces}, "
                 f"max={faces.maximum:.6f} deg")
    ax.set_xlabel("x")
    ax.set_ylabel("z")
    ax.set_aspect("equal")
    ax.legend(fontsize=7, loc="upper left")
    ax.grid(alpha=0.3)

    fig.suptitle("The two conforming stacks: an apex-up/apex-down pyramid pair is the "
                 "hex-to-tet transition; a prism column is the boundary layer", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    path = OUT / "fig2_stacks.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def fig_context():
    rows = AUDIT["meshes"]
    names = [r["name"] for r in rows]
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 4.4))
    ax = axes[0]
    values = [max(r["max_abs_difference_second_route"], 1e-19) for r in rows]
    colors = ["#2ca02c" if all(r["criteria"].values()) else "#d62728" for r in rows]
    ax.barh(names, values, color=colors)
    for index, value in enumerate(values):
        ax.text(value * 1.6, index, f"{value:.2e}", va="center", fontsize=8)
    ax.set_xscale("log")
    ax.axvline(1e-6, ls="--", color="0.4", lw=1.0)
    ax.text(1.2e-6, -0.45, "acceptance 1e-6", fontsize=8, color="0.3")
    ax.set_xlabel("max |library - second route| over all metrics")
    ax.set_title("(a) metric set vs finite-difference / Cayley-Menger route")
    ax.set_xlim(1e-13, 1e-4)
    ax.grid(alpha=0.3, axis="x", which="both")

    ax = axes[1]
    nonzero = [r for r in rows if r["volume_error"] is not None]
    values = []
    for row in nonzero:
        values.append(max(row["volume_error"], 1e-19))
    ax.barh([r["name"] for r in nonzero], values, color="#1f77b4")
    for index, value in enumerate(values):
        ax.text(value * 1.6, index, f"{value:.2e}", va="center", fontsize=8)
    ax.set_xscale("log")
    ax.set_xlabel("|summed element volume - exact geometric volume|")
    ax.set_title("(b) volume vs closed-form geometric volume")
    ax.set_xlim(1e-19, 1e-12)
    ax.grid(alpha=0.3, axis="x", which="both")
    fig.suptitle("Context: the same audit's per-element criteria, all satisfied "
                 "(the two FAIL gates are deliberate defect meshes)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    path = OUT / "fig3_context.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


if __name__ == "__main__":
    for path in (fig_convergence(), fig_geometry(), fig_context()):
        print(path, path.stat().st_size, "bytes")
        print("  b64 head:", base64.b64encode(path.read_bytes())[:40].decode())