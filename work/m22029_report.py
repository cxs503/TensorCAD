#!/usr/bin/env python3
"""m22029 report generator: named-boundary SU2 export of the body grid.

Reads ``runs/m19950_volume_bridge/bridge_matrix.json`` (written by
``examples/meshio_volume_bridge/run_example.py``) plus a grid-refinement study computed here,
and emits the figures and the markdown table for the mission report.  No number in the
report is typed by hand: everything is read back from the JSON or computed in this script.

    python work/m22029_report.py            # writes runs/m22029_report/{figures,report.json}
"""
from __future__ import annotations

import base64
import json
import math
import sys
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from tensorcad import volume_mesh as vmlib
from tensorcad.meshio_exchange import volume_mesh_patches

OUT = REPO / "runs" / "m22029_report"
BRIDGE = REPO / "runs" / "m19950_volume_bridge" / "bridge_matrix.json"


def hull_section_grid(length: float, beam: float, draft: float, nx: int, nz: int):
    xs = tuple(-length / 2 + i * (length / (nx - 1)) for i in range(nx))
    zs = tuple(j * (draft / (nz - 1)) for j in range(nz))

    def breadth(i: int, j: int) -> float:
        x, z = xs[i], zs[j]
        return max(0.0, (beam / 2) * (1 - (2 * x / length) ** 2) * (1 - (z / draft) ** 2))

    return vmlib.SectionGrid(xs=xs, zs=zs, breadth=breadth)


def wetted_area(length: float, beam: float, draft: float,
                nx: int = 400, nz: int = 200) -> float:
    """``2 \\iint sqrt(1 + b_x^2 + b_z^2) dx dz`` from the breadth function (midpoint rule)."""
    total = 0.0
    dx, dz = length / nx, draft / nz
    for i in range(nx):
        x = -length / 2 + (i + 0.5) * dx
        for j in range(nz):
            z = (j + 0.5) * dz
            bx = (beam / 2) * (-8.0 * x / length**2) * (1 - (z / draft) ** 2)
            bz = (beam / 2) * (1 - (2 * x / length) ** 2) * (-2.0 * z / draft**2)
            total += float(np.sqrt(1.0 + bx * bx + bz * bz))
    return 2.0 * total * dx * dz


def refinement_study(length: float, beam: float, draft: float) -> list[dict[str, Any]]:
    """Mesh ``hull`` patch area against the quadrature, at three grid resolutions."""
    reference = wetted_area(length, beam, draft)
    rows = []
    for nx, nz in ((25, 9), (49, 17), (97, 25)):
        mesh = vmlib.build_volume_mesh(hull_section_grid(length, beam, draft, nx, nz),
                                       first_layer=0.004, growth=1.25, layers=2)
        patches = volume_mesh_patches(mesh)
        area = patches.area("hull")
        rows.append({
            "nx": nx, "nz": nz,
            "hull_area": area,
            "reference": reference,
            "rel_error": abs(area - reference) / reference,
            "watertight": patches.is_watertight(),
            "slit_cancelled": patches.cancelled_faces,
        })
    return rows


def figure_skin(patches, path: Path) -> None:
    """The skin that goes into the file: triangles coloured by patch, in plan and profile."""
    fig = plt.figure(figsize=(13.0, 4.6))
    colours = {"hull": "#1f5fa8", "bottom": "#d35400", "top": "#27ae60",
               "fore": "#8e44ad", "aft": "#c0392b", "symmetry": "#7f8c8d"}
    for panel, (a, b, label) in enumerate(((0, 1, "y [m] (port $\\rightarrow$ starboard)"),
                                          (0, 2, "z [m]")), start=1):
        ax = fig.add_subplot(1, 3, panel)
        for name, faces in patches.patches:
            segments = []
            for face in faces:
                for i, j in ((0, 1), (1, 2), (2, 0)):
                    segments.append(((patches.vertices[face[i]][a], patches.vertices[face[i]][b]),
                                     (patches.vertices[face[j]][a], patches.vertices[face[j]][b])))
            for (x0, y0), (x1, y1) in segments:
                ax.plot([x0, x1], [y0, y1], color=colours.get(name, "k"), linewidth=0.25,
                        alpha=0.9)
        ax.set_xlabel("x [m]")
        ax.set_ylabel(label)
        ax.set_title(f"({'ab'[panel - 1]}) named skin, {panel} view "
                     f"({len(patches.skin())} triangles)")
        ax.grid(alpha=0.2)
        if panel == 2:
            ax.set_aspect("equal", adjustable="datalim")

    ax = fig.add_subplot(1, 3, 3)
    names = list(patches.counts())
    areas = [patches.area(name) for name in names]
    ax.bar(range(len(names)), areas, color=[colours.get(n, "k") for n in names])
    for index, (name, area) in enumerate(zip(names, areas)):
        ax.text(index, area, f"{area:.4f}", ha="center", va="bottom", fontsize=7)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names)
    ax.set_ylabel("area [m$^2$]")
    ax.set_title("(c) patch areas of the exported skin")
    ax.grid(alpha=0.2, axis="y")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def figure_checks(markers: list[dict[str, Any]], study: list[dict[str, Any]],
                  path: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.2))

    ax = axes[0]
    labels, values = [], []
    for row in markers:
        ref = row["hull_area_reference"]
        labels += [f"{row['body'].split(' ')[0]}\nhull"]
        values.append(abs(row["areas"]["hull"] - ref) / ref)
        for name, reference in row.get("surface_area_reference", {}).items():
            labels.append(f"{row['body'].split(' ')[0]}\n{name}")
            values.append(abs(row["areas"][name] - reference) / reference)
    floor = 1e-18
    ax.bar(range(len(values)), [max(v, floor) for v in values], color="#1f5fa8")
    ax.set_yscale("log")
    ax.set_xticks(range(len(values)))
    ax.set_xticklabels(labels, fontsize=6, rotation=45, ha="right")
    ax.axhline(1e-12, color="#c0392b", ls="--", lw=1, label="round-off floor")
    ax.set_ylabel("relative area error")
    ax.set_title("(a) patch areas vs closed-form anchors")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.2, axis="y")

    ax = axes[1]
    resolution = [row["nx"] for row in study]
    errors = [row["rel_error"] for row in study]
    ax.loglog(resolution, errors, "o-", color="#1f5fa8", label="mesh hull area vs quadrature")
    ax.loglog(resolution, [errors[0] * (resolution[0] / n) ** 1.0 for n in resolution],
              "k--", lw=1, label="O(h) reference")
    ax.set_xlabel("stations nx (levels nz = (nx-1)/3 + 1)")
    ax.set_ylabel("relative error")
    ax.set_title("(b) discretisation error closes under refinement")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.25, which="both")

    ax = axes[2]
    rows = [row for row in markers]
    for index, row in enumerate(rows):
        ax.bar(index, row["slit_triangles_cancelled"], color="#c0392b",
               label="slit pairs cancelled" if index == 0 else None)
        ax.bar(index, row["skin_triangles"], bottom=row["slit_triangles_cancelled"],
               color="#1f5fa8", label="skin (exported)" if index == 0 else None)
        ax.text(index, row["raw_boundary_triangles"], f"{row['raw_boundary_triangles']}",
                ha="center", va="bottom", fontsize=7)
    ax.set_xticks(range(len(rows)))
    ax.set_xticklabels([row["body"] for row in rows], fontsize=7)
    ax.set_ylabel("triangles")
    ax.set_title("(c) raw 'boundary' faces -> exported skin")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.2, axis="y")

    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    bridge = json.loads(BRIDGE.read_text())
    markers = bridge["markers"]
    study = refinement_study(10.0, 1.0, 0.625)

    hull_mesh = vmlib.build_volume_mesh(hull_section_grid(10.0, 1.0, 0.625, 25, 9),
                                        first_layer=0.004, growth=1.25, layers=4)
    patches = volume_mesh_patches(hull_mesh)
    fig_a = OUT / "skin_patches.png"
    fig_b = OUT / "verification.png"
    figure_skin(patches, fig_a)
    figure_checks(markers, study, fig_b)

    table = ["| body | patch | triangles | area [m²] |", "|---|---|---|---|"]
    for row in markers:
        for name, count in row["patches"].items():
            table.append(f"| {row['body']} | `{name}` | {count} | {row['areas'][name]:.9f} |")

    payload = {"markers": markers, "refinement": study, "patch_table": "\n".join(table),
               "figures": [str(fig_a.relative_to(REPO)), str(fig_b.relative_to(REPO))]}
    (OUT / "report.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload["refinement"], indent=1))
    print(payload["patch_table"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())