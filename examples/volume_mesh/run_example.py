#!/usr/bin/env python3
"""Measure the volume mesher: exact box, hull convergence, wall layer, VTK, refusals.

Writes ``artifacts/volume_mesh/results.json`` (every number the document quotes),
``artifacts/volume_mesh/plot.png`` (the four panels) and ``artifacts/volume_mesh/hull.vtk``
(the mesh itself), and exits non-zero if any claim stops holding.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from tensorcad.hull import HullParameters, analytic_volume, build_mesh  # noqa: E402
from tensorcad.mesh import total_signed_volume  # noqa: E402
from tensorcad.volume_mesh import (  # noqa: E402
    SectionGrid,
    VolumeMeshError,
    box_section_grid,
    build_volume_mesh,
    volume_mesh_from_hull,
)

ARTIFACT_DIR = REPO_ROOT / "artifacts" / "volume_mesh"
LENGTH, BEAM, DRAFT = 10.0, 2.0, 0.8
ANALYTIC = analytic_volume(LENGTH, BEAM, DRAFT)


def stat(section: dict) -> dict:
    """Trim a Jacobian section to the numbers the document quotes (None-safe)."""
    return {
        "count": section["count"],
        "min": section["min"],
        "median": section["median"],
        "mean": section["mean"],
        "below_threshold": section["below_threshold"],
        "below_threshold_fraction": section["below_threshold_fraction"],
        "non_positive": section["non_positive"],
    }


def wall_normal(mesh, i: int, j: int, k: int) -> tuple[float, float, float]:
    """Outward normal at a wall node, from the mesh's own wall nodes (cross of tangents)."""
    nx, nz, breadth_cells = mesh.shape
    i0, i1 = max(0, i - 1), min(nx - 1, i + 1)
    j0, j1 = max(0, j - 1), min(nz - 1, j + 1)
    p = mesh.vertices[mesh.node_id(i, j, k)]
    a = mesh.vertices[mesh.node_id(i1, j, k)]
    b = mesh.vertices[mesh.node_id(i0, j, k)]
    c = mesh.vertices[mesh.node_id(i, j1, k)]
    d = mesh.vertices[mesh.node_id(i, j0, k)]
    t1 = tuple(a[m] - b[m] for m in range(3))
    t2 = tuple(c[m] - d[m] for m in range(3))
    normal = (
        t1[1] * t2[2] - t1[2] * t2[1],
        t1[2] * t2[0] - t1[0] * t2[2],
        t1[0] * t2[1] - t1[1] * t2[0],
    )
    length = math.sqrt(sum(component * component for component in normal))
    return tuple(component / length for component in normal)  # type: ignore[return-value]


def normal_clearance(mesh, i: int, j: int, k_wall: int, inner_node) -> float:
    """Distance from ``inner_node`` to the wall *plane* at wall node ``(i, j, k_wall)``.

    Measured independently of how the mesher placed the node: the normal comes from the
    cross product of the wall's own node tangents, the origin is the wall node, and the
    quantity is the perpendicular projection of the gap.
    """
    wall_node = mesh.vertices[mesh.node_id(i, j, k_wall)]
    normal = wall_normal(mesh, i, j, k_wall)
    gap = tuple(inner_node[m] - wall_node[m] for m in range(3))
    return abs(sum(gap[m] * normal[m] for m in range(3)))


def rows_graded(column) -> bool:
    """True when the column carries at least two genuinely growing layers."""
    return len(column) > 1 and column[1] > column[0] * 1.01


def main() -> int:
    ARTIFACT_DIR.mkdir(parents=True, exist_ok=True)
    results: dict = {"case": {"length": LENGTH, "beam": BEAM, "draft": DRAFT, "analytic_volume_m3": ANALYTIC}}
    failures: list[str] = []

    # ---------------------------------------------------------------- 1. box control
    box = build_volume_mesh(box_section_grid(length=LENGTH, beam=BEAM, draft=DRAFT, nx=11, nz=7),
                            first_layer=0.01, growth=1.2, layers=4)
    box_analytic = LENGTH * BEAM * DRAFT
    box_volume_error = box.volume() / box_analytic - 1.0
    box_summary = box.summary()
    results["box"] = {
        "analytic_volume_m3": box_analytic,
        "mesh_volume_m3": box.volume(),
        "relative_error": box_volume_error,
        "first_layer_requested_m": 0.01,
        "first_layer_measured_m": box.columns[0][0],
        "growth_requested": 1.2,
        "growth_measured": box.columns[0][1] / box.columns[0][0],
        "counts": box_summary["counts"],
        "conforming": box_summary["conforming"],
        "quality_hex": {k: stat(v) for k, v in box_summary["quality"]["hex"].items()},
        "quality_tet": {k: stat(v) for k, v in box_summary["quality"]["tet"].items()},
        "collapsed_cells": box_summary["quality"]["collapsed_cells"],
    }
    if abs(box_volume_error) > 1e-12:
        failures.append(f"box volume error {box_volume_error:.3e} is not round-off")
    if abs(box.columns[0][0] - 0.01) > 1e-15:
        failures.append("box first layer is not the requested thickness")

    # ---------------------------------------------------------------- 2. hull convergence
    resolutions = [(41, 13), (81, 21), (161, 41)]
    convergence = []
    meshes = {}
    for nx, nz in resolutions:
        parameters = HullParameters(length=LENGTH, beam=BEAM, draft=DRAFT, nx=nx, nz=nz)
        mesh = volume_mesh_from_hull(parameters, first_layer=0.005, growth=1.2, layers=4)
        summary = mesh.summary()
        surface, _offsets = build_mesh(parameters)
        convergence.append(
            {
                "nx": nx,
                "nz": nz,
                "cells": summary["counts"],
                "volume_m3": summary["volume_m3"],
                "volume_relative_error": summary["volume_m3"] / ANALYTIC - 1.0,
                "surface_mesh_volume_m3": total_signed_volume(surface),
                "surface_mesh_relative_error": total_signed_volume(surface) / ANALYTIC - 1.0,
                "volume_split_relative_difference": summary["volume_split_relative_difference"],
                "conforming": summary["conforming"],
                "boundary_triangles": summary["boundary_triangles"],
                "collapsed_columns": summary["collapsed_columns"],
                "columns": summary["columns"],
                "graded_columns": summary["graded_columns"],
                "collapsed_cells": summary["quality"]["collapsed_cells"],
                "quality_hex": {k: stat(v) for k, v in summary["quality"]["hex"].items()},
                "quality_tet": {k: stat(v) for k, v in summary["quality"]["tet"].items()},
            }
        )
        meshes[(nx, nz)] = mesh
        if not summary["conforming"]:
            failures.append(f"hull {nx}x{nz} split is not conforming")
        if summary["volume_split_relative_difference"] > 1e-12:
            failures.append(f"hull {nx}x{nz} hex and tet volumes disagree")
        if abs(summary["volume_m3"] - total_signed_volume(surface)) > 1e-9 * ANALYTIC:
            failures.append(f"hull {nx}x{nz} volume mesh does not reproduce the surface volume")
    order = []
    for coarse, fine in zip(convergence, convergence[1:]):
        e_coarse, e_fine = abs(coarse["volume_relative_error"]), abs(fine["volume_relative_error"])
        if e_coarse > 0.0 and e_fine > 0.0:
            order.append(
                {
                    "coarse_stations": coarse["nx"],
                    "fine_stations": fine["nx"],
                    "observed_order": math.log(e_coarse / e_fine)
                    / math.log((fine["nx"] - 1) / (coarse["nx"] - 1)),
                }
            )
    results["convergence"] = {
        "rows": convergence,
        "observed_order_station_refinement": order,
        "analytic_volume_m3": ANALYTIC,
    }

    # ---------------------------------------------------------------- 3. wall layer
    parameters = HullParameters(length=LENGTH, beam=BEAM, draft=DRAFT, nx=81, nz=21)
    layer_rows = []
    for first_layer in (0.002, 0.005, 0.02):
        for growth, layers in ((1.0, 2), (1.2, 4), (1.5, 4), (1.2, 8)):
            mesh = volume_mesh_from_hull(
                parameters, first_layer=first_layer, growth=growth, layers=layers, split_tets=False
            )
            summary = mesh.summary()
            # node one step inside the starboard wall, at midship and mid-draught
            i = (parameters.nx - 1) // 2
            j = (parameters.nz - 1) // 2
            column = mesh.columns[i * parameters.nz + j]
            wall_node = mesh.vertices[mesh.node_id(i, j, 2 * layers)]
            inner_node = mesh.vertices[mesh.node_id(i, j, 2 * layers - 1)]
            measured = normal_clearance(mesh, i, j, 2 * layers, inner_node)
            ratio_measured = column[1] / column[0] if rows_graded(column) else None
            layer_rows.append(
                {
                    "first_layer_requested_m": first_layer,
                    "growth_requested": growth,
                    "layers": layers,
                    "first_layer_y_step_m": column[0],
                    "growth_measured": ratio_measured,
                    "normal_clearance_m": measured,
                    "normal_clearance_relative_error": measured / first_layer - 1.0,
                    "wall_node_gap_m": math.dist(wall_node, inner_node),
                    "tangent_normal": [round(v, 12) for v in wall_normal(mesh, i, j, 2 * layers)],
                    "graded_columns": summary["graded_columns"],
                    "columns": summary["columns"],
                    "quality_hex_all_min": summary["quality"]["hex"]["all"]["min"],
                    "volume_relative_error": summary["volume_m3"] / ANALYTIC - 1.0,
                }
            )
    worst = max(abs(row["normal_clearance_relative_error"]) for row in layer_rows)
    results["wall_layer"] = {"rows": layer_rows, "worst_normal_clearance_relative_error": worst}
    if worst > 0.05:
        failures.append(f"wall-normal first layer missed the request by {worst:.2%}")

    # ---------------------------------------------------------------- 4. VTK round trip
    mesh = meshes[(41, 13)]
    vtk_path = mesh.to_vtk(ARTIFACT_DIR / "hull.vtk")
    signed = 0.0
    points, cells = [], []
    with vtk_path.open() as handle:
        lines = handle.read().splitlines()
    n_points = int(lines[4].split()[1])
    for line in lines[5 : 5 + n_points]:
        points.append(tuple(float(v) for v in line.split()))
    cells_header = lines[5 + n_points].split()
    n_cells = int(cells_header[1])
    for line in lines[6 + n_points : 6 + n_points + n_cells]:
        ids = [int(v) for v in line.split()[1:]]
        cells.append(ids)
    for ids in cells:
        p0, p1, p2, p3, p4, p5, p6, p7 = (points[i] for i in ids)
        for tet in ((p0, p1, p2, p6), (p0, p2, p3, p6), (p0, p3, p7, p6), (p0, p7, p4, p6), (p0, p4, p5, p6), (p0, p5, p1, p6)):
            a, b, c, d = tet
            ab = [b[k] - a[k] for k in range(3)]
            ac = [c[k] - a[k] for k in range(3)]
            ad = [d[k] - a[k] for k in range(3)]
            det = (
                ab[0] * (ac[1] * ad[2] - ac[2] * ad[1])
                - ab[1] * (ac[0] * ad[2] - ac[2] * ad[0])
                + ab[2] * (ac[0] * ad[1] - ac[1] * ad[0])
            )
            signed += det / 6.0
    results["vtk"] = {
        "path": str(vtk_path.relative_to(REPO_ROOT)),
        "points": n_points,
        "cells": n_cells,
        "volume_from_cells_m3": signed,
        "mesh_volume_m3": mesh.volume(),
        "positive_orientation": signed > 0,
    }
    if abs(signed - mesh.volume()) > 1e-9 * ANALYTIC:
        failures.append("VTK round trip does not reproduce the mesh volume")
    if signed <= 0:
        failures.append("VTK cells are not positively oriented")

    # ---------------------------------------------------------------- 5. refusals
    refusals = []
    attempts = {
        "negative first layer": lambda: build_volume_mesh(box_section_grid(length=1, beam=1, draft=1), first_layer=-0.1),
        "growth below one": lambda: build_volume_mesh(box_section_grid(length=1, beam=1, draft=1), first_layer=0.1, growth=0.5),
        "zero layers": lambda: build_volume_mesh(box_section_grid(length=1, beam=1, draft=1), first_layer=0.1, layers=0),
        "non-finite breadth": lambda: SectionGrid(xs=(0.0, 1.0), zs=(0.0, 1.0), breadth=lambda i, j: math.nan),
        "negative breadth": lambda: SectionGrid(xs=(0.0, 1.0), zs=(0.0, 1.0), breadth=lambda i, j: -1.0),
        "non-monotone stations": lambda: SectionGrid(xs=(0.0, 0.0), zs=(0.0, 1.0), breadth=lambda i, j: 1.0),
        "one station": lambda: SectionGrid(xs=(0.0,), zs=(0.0, 1.0), breadth=lambda i, j: 1.0),
        "unknown vtk cell type": lambda: box.to_vtk(ARTIFACT_DIR / "x.vtk", cells="wedge"),
    }
    for name, attempt in attempts.items():
        try:
            attempt()
        except VolumeMeshError:
            refusals.append({"case": name, "refused": True})
        else:
            refusals.append({"case": name, "refused": False})
            failures.append(f"{name} was not refused")
    results["refusals"] = refusals

    # ---------------------------------------------------------------- figure
    figure, axes = plt.subplots(2, 2, figsize=(11.5, 8.0))
    ax = axes[0][0]
    mid = meshes[(81, 21)]
    i = (81 - 1) // 2
    for j in range(21 - 1):
        for k in range(2 * 4):
            ids = mid.cell_nodes(i, j, k)
            face = [mid.vertices[ids[a]] for a in (4, 5, 6, 7)]
            ax.plot([p[1] for p in face + [face[0]]], [p[2] for p in face + [face[0]]], "-", lw=0.4, color="k")
    ax.set_title("midship section: graded wall layers")
    ax.set_xlabel("y [m]")
    ax.set_ylabel("z [m]")
    ax.set_aspect("equal")

    ax = axes[0][1]
    ns = [(row["nx"] - 1) for row in convergence]
    ax.loglog(ns, [abs(row["volume_relative_error"]) for row in convergence], "o-", label="volume mesh vs analytic")
    ax.loglog(ns, [abs(row["surface_mesh_relative_error"]) for row in convergence], "s--", label="hull surface mesh vs analytic")
    ax.set_xlabel("station intervals")
    ax.set_ylabel("|relative volume error|")
    ax.set_title("volume conservation follows the surface")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)

    ax = axes[1][0]
    requested = [row["first_layer_requested_m"] for row in layer_rows]
    measured = [row["normal_clearance_m"] for row in layer_rows]
    ax.plot([min(requested) * 0.8, max(requested) * 1.2], [min(requested) * 0.8, max(requested) * 1.2], "k--", lw=1, label="requested")
    ax.plot(requested, measured, "o", label="measured normal clearance")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("requested first layer [m]")
    ax.set_ylabel("measured wall-normal clearance [m]")
    ax.set_title("the first layer is a normal clearance, not a y step")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)

    ax = axes[1][1]
    hex_all = [row["quality_hex"]["all"]["median"] for row in convergence]
    tet_all = [row["quality_tet"]["all"]["median"] for row in convergence]
    ax.plot(ns, hex_all, "o-", label="hex (median)")
    ax.plot(ns, tet_all, "s-", label="tet split (median)")
    ax.axhline(0.1, color="r", ls=":", label="distortion threshold")
    ax.set_xlabel("station intervals")
    ax.set_ylabel("scaled Jacobian")
    ax.set_title("cell quality by cell type")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    figure.tight_layout()
    figure.savefig(ARTIFACT_DIR / "plot.png", dpi=140)

    results["failures"] = failures
    (ARTIFACT_DIR / "results.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(json.dumps({"failures": failures, "worst_clearance_error": worst}, indent=2))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())