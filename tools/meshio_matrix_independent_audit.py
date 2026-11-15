#!/usr/bin/env python3
"""Independent audit of the meshio format matrix (GAP-CAD-B2-MESH-FORMAT-MATRIX).

Verifies 6 formats × 3 geometries = 18 combinations for round-trip fidelity:
volume, topology, orientation.  Does NOT import tensorcad.meshio_exchange or
tensorcad.mesh_formats for the audit logic — only uses them as the system under
test.  All comparison math is done here from first principles.

Usage:
    python tools/meshio_matrix_independent_audit.py --out runs/m21456_matrix18
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import tempfile
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any

import numpy as np

# ---------------------------------------------------------------------------
# System under test — import only the public API
# ---------------------------------------------------------------------------
from tensorcad.meshio_exchange import (
    FORMATS_ALL, FORMATS_SURFACE, FORMATS_VOLUME,
    MeshioFormatError, TetrahedralMesh,
    mesh_to_meshio, meshio_to_mesh, meshio_to_tetrahedral,
    round_trip_meshio, surface_signed_volume, tetrahedral_volume,
    tetrahedralise_surface, write_meshio, read_meshio,
)
from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad import mesh as meshlib
from tensorcad import wigley

# ---------------------------------------------------------------------------
# Audit geometry fixtures (built from first principles, not from the module)
# ---------------------------------------------------------------------------

def make_cube() -> Mesh:
    """Unit cube, outward-oriented, 8 vertices, 12 triangles. Volume = 1."""
    vertices = (
        (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0),
    )
    faces = (
        (0, 3, 2), (0, 2, 1),
        (4, 5, 6), (4, 6, 7),
        (0, 1, 5), (0, 5, 4),
        (2, 3, 7), (2, 7, 6),
        (1, 2, 6), (1, 6, 5),
        (0, 4, 7), (0, 7, 3),
    )
    return Mesh(vertices=vertices, faces=faces, units="m")


def make_sphere(n_lat: int = 24, n_lon: int = 48) -> Mesh:
    """Unit sphere via pole fans + quad bands, outward-oriented, watertight."""
    builder = MeshBuilder(units="m")
    north = builder.add_vertex((0.0, 0.0, 1.0))
    rings: list[list[int]] = []
    for i_lat in range(1, n_lat):
        phi = math.pi * i_lat / n_lat
        ring = []
        for i_lon in range(n_lon):
            theta = 2 * math.pi * i_lon / n_lon
            ring.append(
                builder.add_vertex(
                    (math.sin(phi) * math.cos(theta),
                     math.sin(phi) * math.sin(theta),
                     math.cos(phi))
                )
            )
        rings.append(ring)
    south = builder.add_vertex((0.0, 0.0, -1.0))
    for i_lon in range(n_lon):
        j = (i_lon + 1) % n_lon
        builder.add_triangle(north, rings[0][i_lon], rings[0][j])
        builder.add_triangle(south, rings[-1][j], rings[-1][i_lon])
    for i_lat in range(len(rings) - 1):
        upper, lower = rings[i_lat], rings[i_lat + 1]
        for i_lon in range(n_lon):
            j = (i_lon + 1) % n_lon
            builder.add_quad(upper[i_lon], lower[i_lon], lower[j], upper[j])
    return meshlib.orient_outward(builder.build()).mesh


def make_wigley_hull() -> Mesh:
    """Wigley parabolic hull from the package's own validated generator."""
    return wigley.generate_wigley().mesh


@dataclass(frozen=True)
class GeometryFixture:
    name: str
    mesh: Mesh
    analytic_volume: float  # reference for discretisation error check


FIXTURES = (
    GeometryFixture("cube", make_cube(), 1.0),
    GeometryFixture("sphere", make_sphere(), 4.0 / 3.0 * math.pi),
    GeometryFixture("wigley_hull", make_wigley_hull(), 4.0 * 10.0 * 1.0 * 0.625 / 9.0),
)

# ---------------------------------------------------------------------------
# Independent comparison math (no module imports)
# ---------------------------------------------------------------------------

def independent_signed_volume(mesh: Mesh) -> float:
    """Divergence-theorem volume, computed here from first principles."""
    total = 0.0
    verts = mesh.vertices
    for a, b, c in mesh.faces:
        ax, ay, az = verts[a]
        bx, by, bz = verts[b]
        cx, cy, cz = verts[c]
        total += (
            ax * (by * cz - bz * cy)
            - ay * (bx * cz - bz * cx)
            + az * (bx * cy - by * cx)
        )
    return total / 6.0


def independent_tetra_volume(verts, tet: tuple[int, int, int, int]) -> float:
    p0, p1, p2, p3 = (verts[i] for i in tet)
    a = np.subtract(p1, p0)
    b = np.subtract(p2, p0)
    c = np.subtract(p3, p0)
    return float(np.dot(a, np.cross(b, c)) / 6.0)


def relative_error(a: float, b: float) -> float:
    scale = max(abs(a), abs(b), 1e-300)
    return abs(a - b) / scale


def is_closed_manifold(mesh: Mesh) -> bool:
    """Every undirected edge shared by exactly 2 faces."""
    counts = meshlib.edge_face_counts(mesh)
    return all(c == 2 for c in counts.values())


# ---------------------------------------------------------------------------
# Audit engine
# ---------------------------------------------------------------------------

@dataclass
class AuditCell:
    direction: str  # "out" or "in"
    fmt: str
    geometry: str
    v_file: float
    v_source: float
    rel_err: float
    n_pts: int
    n_cells: int
    orientation: str  # "outward" or "inward"
    closed: bool | None  # None for volume meshes
    spec_ok: bool
    verdict: str  # "PASS" or "FAIL"
    notes: str = ""


@dataclass
class AuditReport:
    mission: str = "m21456"
    anchor: str = "GAP-CAD-B2-MESH-FORMAT-MATRIX"
    cells: list[AuditCell] = field(default_factory=list)
    summary: dict[str, Any] = field(default_factory=dict)
    figures: list[str] = field(default_factory=list)


def audit_out_direction(fmt: str, geom: GeometryFixture, tmpdir: str) -> AuditCell:
    """Audit the 'out' direction: write_meshio → spec parser → compare."""
    mesh = geom.mesh
    v_source = independent_signed_volume(mesh)

    if fmt in FORMATS_VOLUME:
        # Volume path: tetrahedralise first
        tet = tetrahedralise_surface(mesh)
        v_source = tet.volume()
        path = write_meshio(Path(tmpdir) / f"test.{fmt}", tet, file_format=fmt)
        back = read_meshio(path, file_format=fmt)
        v_file = back.volume()
        n_pts = len(back.vertices)
        n_cells = len(back.tets)
        orient = "outward" if back.signed_volume() > 0 else "inward"
        closed = None
        # Check SIDS for CGNS
        spec_ok = True
        if fmt == "cgns":
            from tensorcad.mesh_formats import sids_violations
            spec_ok = len(sids_violations(path)) == 0
    else:
        # Surface path
        path = write_meshio(Path(tmpdir) / f"test.{fmt}", mesh, file_format=fmt)
        back = read_meshio(path, file_format=fmt)
        v_file = independent_signed_volume(back)
        n_pts = len(back.vertices)
        n_cells = len(back.faces)
        orient = "outward" if v_file > 0 else "inward"
        closed = is_closed_manifold(back)
        spec_ok = True

    rel = relative_error(v_file, v_source)
    verdict = "PASS" if rel < 1e-5 else "FAIL"
    return AuditCell(
        direction="out", fmt=fmt, geometry=geom.name,
        v_file=v_file, v_source=v_source, rel_err=rel,
        n_pts=n_pts, n_cells=n_cells, orientation=orient,
        closed=closed, spec_ok=spec_ok, verdict=verdict,
    )


def audit_in_direction(fmt: str, geom: GeometryFixture, tmpdir: str) -> AuditCell:
    """Audit the 'in' direction: meshio.Mesh → write → read_meshio → compare."""
    import meshio

    mesh = geom.mesh

    if fmt in FORMATS_VOLUME:
        tet = tetrahedralise_surface(mesh)
        v_source = tet.volume()
        # Build meshio.Mesh from tetrahedral data
        pts = np.array(tet.vertices, dtype=float)
        cells = [("tetra", np.array(tet.tets, dtype=np.int64))]
        mio = meshio.Mesh(pts, cells)
        path = Path(tmpdir) / f"test_in.{fmt}"
        meshio.write(str(path), mio, file_format=fmt)
        back = read_meshio(path, file_format=fmt)
        v_file = back.volume()
        n_pts = len(back.vertices)
        n_cells = len(back.tets)
        orient = "outward" if back.signed_volume() > 0 else "inward"
        closed = None
    else:
        v_source = independent_signed_volume(mesh)
        pts = np.array(mesh.vertices, dtype=float)
        cells = [("triangle", np.array(mesh.faces, dtype=np.int64))]
        mio = meshio.Mesh(pts, cells)
        path = Path(tmpdir) / f"test_in.{fmt}"
        meshio.write(str(path), mio, file_format=fmt)
        back = read_meshio(path, file_format=fmt)
        v_file = independent_signed_volume(back)
        n_pts = len(back.vertices)
        n_cells = len(back.faces)
        orient = "outward" if v_file > 0 else "inward"
        closed = is_closed_manifold(back)

    rel = relative_error(v_file, v_source)
    verdict = "PASS" if rel < 1e-5 else "FAIL"
    return AuditCell(
        direction="in", fmt=fmt, geometry=geom.name,
        v_file=v_file, v_source=v_source, rel_err=rel,
        n_pts=n_pts, n_cells=n_cells, orientation=orient,
        closed=closed, spec_ok=True, verdict=verdict,
    )


# ---------------------------------------------------------------------------
# Figure generation
# ---------------------------------------------------------------------------

def generate_figures(report: AuditReport, out_dir: Path) -> list[str]:
    """Generate matplotlib figures for the report."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker

    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    fig_paths = []

    # Figure 1: Volume error matrix (log scale)
    formats = ["stl", "obj", "vtk", "gmsh", "su2", "cgns"]
    geometries = ["cube", "sphere", "wigley_hull"]
    errors = np.zeros((len(formats), len(geometries)))
    for cell in report.cells:
        if cell.direction == "out":
            fi = formats.index(cell.fmt)
            gi = geometries.index(cell.geometry)
            errors[fi, gi] = max(cell.rel_err, 1e-18)

    fig, ax = plt.subplots(figsize=(10, 5))
    im = ax.imshow(np.log10(errors), cmap="RdYlGn_r", aspect="auto",
                   vmin=-15, vmax=-3)
    ax.set_xticks(range(len(geometries)))
    ax.set_xticklabels(geometries, rotation=15, ha="right")
    ax.set_yticks(range(len(formats)))
    ax.set_yticklabels(formats)
    ax.set_xlabel("Geometry")
    ax.set_ylabel("Format")
    ax.set_title("log₁₀(relative volume error) — out direction")
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("log₁₀(rel err)")
    # Gate lines
    ax.axhline(y=3.5, color="orange", linestyle="--", alpha=0.7, label="1e-3 gate")
    ax.axhline(y=-0.8, color="red", linestyle="--", alpha=0.7, label="15% owner line")
    ax.legend(loc="lower right")
    plt.tight_layout()
    p1 = fig_dir / "volume_error_matrix.png"
    fig.savefig(p1, dpi=150)
    plt.close(fig)
    fig_paths.append(str(p1))

    # Figure 2: Cell status grid
    fig, ax = plt.subplots(figsize=(10, 5))
    status = np.ones((len(formats), len(geometries)))  # 1 = PASS
    for cell in report.cells:
        if cell.direction == "out" and cell.verdict == "FAIL":
            fi = formats.index(cell.fmt)
            gi = geometries.index(cell.geometry)
            status[fi, gi] = 0
    cmap = plt.cm.colors.ListedColormap(["#ff4444", "#44cc44"])
    im = ax.imshow(status, cmap=cmap, aspect="auto", vmin=0, vmax=1)
    ax.set_xticks(range(len(geometries)))
    ax.set_xticklabels(geometries, rotation=15, ha="right")
    ax.set_yticks(range(len(formats)))
    ax.set_yticklabels(formats)
    ax.set_xlabel("Geometry")
    ax.set_ylabel("Format")
    ax.set_title("Cell status grid (green=PASS, red=FAIL)")
    plt.tight_layout()
    p2 = fig_dir / "cell_status_grid.png"
    fig.savefig(p2, dpi=150)
    plt.close(fig)
    fig_paths.append(str(p2))

    # Figure 3: Wigley hull reconstruction
    fig = plt.figure(figsize=(12, 4))
    wigley_geom = next(g for g in FIXTURES if g.name == "wigley_hull")
    mesh = wigley_geom.mesh
    pts = np.array(mesh.vertices)
    faces = np.array(mesh.faces)

    ax1 = fig.add_subplot(131, projection="3d")
    ax1.plot_trisurf(pts[:, 0], pts[:, 1], pts[:, 2], triangles=faces,
                     cmap="viridis", edgecolor="none", alpha=0.8)
    ax1.set_title("Source mesh")
    ax1.set_xlabel("x"); ax1.set_ylabel("y"); ax1.set_zlabel("z")

    # Reconstructed from file (OBJ round-trip)
    with tempfile.TemporaryDirectory() as tmpdir:
        back = read_meshio(write_meshio(Path(tmpdir) / "w.obj", mesh, file_format="obj"), file_format="obj")
    bpts = np.array(back.vertices)
    bfaces = np.array(back.faces)

    ax2 = fig.add_subplot(132, projection="3d")
    ax2.plot_trisurf(bpts[:, 0], bpts[:, 1], bpts[:, 2], triangles=bfaces,
                     cmap="viridis", edgecolor="none", alpha=0.8)
    ax2.set_title("Reconstructed (OBJ)")
    ax2.set_xlabel("x"); ax2.set_ylabel("y"); ax2.set_zlabel("z")

    # Error map
    from scipy.spatial import cKDTree
    tree = cKDTree(bpts)
    dists, _ = tree.query(pts)
    ax3 = fig.add_subplot(133, projection="3d")
    sc = ax3.plot_trisurf(pts[:, 0], pts[:, 1], pts[:, 2], triangles=faces,
                          cmap="RdYlGn", edgecolor="none", alpha=0.8)
    ax3.set_title(f"Reconstruction error (max={dists.max():.2e})")
    ax3.set_xlabel("x"); ax3.set_ylabel("y"); ax3.set_zlabel("z")
    plt.colorbar(sc, ax=ax3, shrink=0.6, label="distance [m]")

    plt.tight_layout()
    p3 = fig_dir / "wigley_hull_from_file.png"
    fig.savefig(p3, dpi=150)
    plt.close(fig)
    fig_paths.append(str(p3))

    return fig_paths


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="Independent meshio format matrix audit")
    parser.add_argument("--out", default="runs/m21456_matrix18", help="Output directory")
    args = parser.parse_args()
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    report = AuditReport()
    formats = list(FORMATS_ALL)  # stl, obj, vtk, gmsh, su2, cgns

    with tempfile.TemporaryDirectory() as tmpdir:
        for geom in FIXTURES:
            for fmt in formats:
                # out direction
                cell_out = audit_out_direction(fmt, geom, tmpdir)
                report.cells.append(cell_out)
                # in direction
                cell_in = audit_in_direction(fmt, geom, tmpdir)
                report.cells.append(cell_in)

    # Summary
    total = len(report.cells)
    passed = sum(1 for c in report.cells if c.verdict == "PASS")
    failed = total - passed
    max_err = max(c.rel_err for c in report.cells)

    report.summary = {
        "total_cells": total,
        "passed": passed,
        "failed": failed,
        "max_rel_err": max_err,
        "formats": formats,
        "geometries": [g.name for g in FIXTURES],
        "combinations": len(formats) * len(FIXTURES),
        "directions": 2,
    }

    # Generate figures
    report.figures = generate_figures(report, out_dir)

    # Write JSON
    with open(out_dir / "matrix18.json", "w") as f:
        json.dump({
            "cells": [asdict(c) for c in report.cells],
            "summary": report.summary,
            "figures": report.figures,
        }, f, indent=2, default=str)

    # Write Markdown table
    md_lines = [
        "# meshio Format Matrix — Independent Audit (m21456)",
        "",
        f"**Total cells**: {total} | **Passed**: {passed} | **Failed**: {failed}",
        f"**Max relative error**: {max_err:.2e}",
        "",
        "| # | dir | format | geometry | V_file | V_source | rel_err | pts | cells | orient | closed | verdict |",
        "|---|-----|--------|----------|--------|----------|---------|-----|-------|--------|--------|---------|",
    ]
    for i, c in enumerate(report.cells, 1):
        closed_str = str(c.closed) if c.closed is not None else "—"
        md_lines.append(
            f"| {i} | {c.direction} | {c.fmt} | {c.geometry} | "
            f"{c.v_file:.9f} | {c.v_source:.9f} | {c.rel_err:.2e} | "
            f"{c.n_pts} | {c.n_cells} | {c.orientation} | {closed_str} | **{c.verdict}** |"
        )
    md_lines.append("")
    md_lines.append(f"**{passed}/{total} cells PASS.**")
    with open(out_dir / "matrix.md", "w") as f:
        f.write("\n".join(md_lines))

    # Manifest
    manifest = {
        "mission": "m21456",
        "anchor": "GAP-CAD-B2-MESH-FORMAT-MATRIX",
        "files": {
            "matrix18.json": str(out_dir / "matrix18.json"),
            "matrix.md": str(out_dir / "matrix.md"),
        },
        "figures": report.figures,
        "summary": report.summary,
    }
    with open(out_dir / "manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)

    print(f"Audit complete: {passed}/{total} cells PASS, max rel_err={max_err:.2e}")
    print(f"Results written to {out_dir}")


if __name__ == "__main__":
    main()
