#!/usr/bin/env python
"""Expanded meshio format matrix audit: 6 formats × 4 geometries = 24 combos.

Adds the torus (genus-1, non-star-shaped) to the existing 3-geometry matrix.
Generates matplotlib figures for visual evidence.

Note: The torus is non-star-shaped, so volume formats (SU2/CGNS) cannot
tetrahedralise it via star tessellation. Those cells are marked N/A.
"""
from __future__ import annotations

import json
import math
import tempfile
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np

from tensorcad import mesh as meshlib
from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.meshio_exchange import (
    FORMATS_ALL,
    FORMATS_SURFACE,
    FORMATS_VOLUME,
    TetrahedralMesh,
    TetrahedralisationError,
    surface_signed_volume,
    tetrahedralise_surface,
    round_trip_meshio,
)
from tensorcad.mesh_formats import cube, icosphere, torus, TriMesh, compare_meshes


def _make_wigley_hull() -> Mesh:
    from tensorcad import wigley
    return wigley.generate_wigley().mesh


def _trimesh_to_mesh(tm: TriMesh) -> Mesh:
    """Convert TriMesh to Mesh for meshio_exchange compatibility."""
    return Mesh(
        vertices=tuple(tuple(float(x) for x in row) for row in tm.points),
        faces=tuple(tuple(int(x) for x in row) for row in tm.faces),
        units="m",
    )


def _surface_volume(mesh: Mesh | TetrahedralMesh) -> float:
    if isinstance(mesh, TetrahedralMesh):
        return mesh.volume()
    return surface_signed_volume(mesh)


def _cell_count(mesh: Mesh | TetrahedralMesh) -> int:
    return len(mesh.tets) if isinstance(mesh, TetrahedralMesh) else len(mesh.faces)


def _vertex_count(mesh: Mesh | TetrahedralMesh) -> int:
    return len(mesh.vertices)


@dataclass
class AuditResult:
    format: str
    geometry: str
    topology_sent: str
    vertices_orig: int
    cells_orig: int
    volume_orig: float
    volume_rt: float
    rel_err: float
    topology_ok: bool
    orientation_ok: bool
    status: str
    volume_analytic: float | None = None
    volume_vs_analytic_rel: float | None = None
    note: str = ""


def run_audit() -> list[AuditResult]:
    # Use Mesh objects for meshio_exchange
    cube_mesh = _trimesh_to_mesh(cube())
    sphere_mesh = _trimesh_to_mesh(icosphere(subdivisions=2))
    torus_mesh = _trimesh_to_mesh(torus(nu=64, nv=32))
    wigley_mesh = _make_wigley_hull()

    geometries = {
        "cube": cube_mesh,
        "sphere": sphere_mesh,
        "torus": torus_mesh,
        "wigley_hull": wigley_mesh,
    }

    analytic_volumes = {
        "cube": 8.0,  # cube() uses edge=2.0
        "sphere": 4.0 / 3.0 * math.pi,
        "torus": 2.0 * math.pi**2 * 1.0 * 0.3**2,
        "wigley_hull": 4.0 * 10.0 * 1.0 * 0.625 / 9.0,
    }

    results: list[AuditResult] = []

    for fmt in FORMATS_ALL:
        for name, geom in geometries.items():
            analytic_vol = analytic_volumes.get(name)

            # For volume formats, try tetrahedralisation
            if fmt in FORMATS_VOLUME:
                try:
                    mesh_in = tetrahedralise_surface(geom)
                    topo_sent = "tetra"
                except TetrahedralisationError as e:
                    # Non-star-shaped body (e.g., torus)
                    results.append(AuditResult(
                        format=fmt,
                        geometry=name,
                        topology_sent="N/A",
                        vertices_orig=_vertex_count(geom),
                        cells_orig=_cell_count(geom),
                        volume_orig=_surface_volume(geom),
                        volume_rt=0.0,
                        rel_err=0.0,
                        topology_ok=False,
                        orientation_ok=False,
                        status="N/A",
                        volume_analytic=analytic_vol,
                        volume_vs_analytic_rel=abs(_surface_volume(geom) - analytic_vol) / analytic_vol if analytic_vol else None,
                        note=f"Non-star-shaped: {str(e)[:80]}",
                    ))
                    continue
            else:
                mesh_in = geom
                topo_sent = "triangle"

            vol_orig = _surface_volume(mesh_in)
            v_orig = _vertex_count(mesh_in)
            c_orig = _cell_count(mesh_in)

            # Round trip
            rt = round_trip_meshio(mesh_in, file_format=fmt)
            vol_rt = _surface_volume(rt)

            scale = max(abs(vol_orig), abs(vol_rt), 1e-300)
            rel_err = abs(vol_rt - vol_orig) / scale

            topo_ok = (_vertex_count(rt) == v_orig and _cell_count(rt) == c_orig)
            orient_ok = math.copysign(1, vol_rt) == math.copysign(1, vol_orig)

            vol_tol = 1e-3 if fmt == "stl" else 1e-12
            vol_ok = rel_err < vol_tol

            status = "PASS" if (vol_ok and topo_ok and orient_ok) else "FAIL"

            vol_vs_analytic = abs(vol_orig - analytic_vol) / analytic_vol if analytic_vol else None

            results.append(AuditResult(
                format=fmt,
                geometry=name,
                topology_sent=topo_sent,
                vertices_orig=v_orig,
                cells_orig=c_orig,
                volume_orig=vol_orig,
                volume_rt=vol_rt,
                rel_err=rel_err,
                topology_ok=topo_ok,
                orientation_ok=orient_ok,
                status=status,
                volume_analytic=analytic_vol,
                volume_vs_analytic_rel=vol_vs_analytic,
            ))

    return results


def generate_markdown_table(results: list[AuditResult]) -> str:
    lines = [
        "# meshio Format Matrix — Expanded Round-Trip Fidelity (6 formats × 4 geometries)",
        "",
        "| Format | Geometry | Topology | V_orig [m³] | V_rt [m³] | rel. err | Topology | Orientation | Status |",
        "|--------|----------|----------|-------------|-----------|----------|----------|-------------|--------|",
    ]
    for r in results:
        if r.status == "N/A":
            lines.append(
                f"| {r.format} | {r.geometry} | {r.topology_sent} "
                f"({r.vertices_orig}v/{r.cells_orig}c) | "
                f"{r.volume_orig:.9f} | — | — | "
                f"N/A | N/A | **N/A** |"
            )
        else:
            lines.append(
                f"| {r.format} | {r.geometry} | {r.topology_sent} "
                f"({r.vertices_orig}v/{r.cells_orig}c) | "
                f"{r.volume_orig:.9f} | {r.volume_rt:.9f} | {r.rel_err:.2e} | "
                f"{'ok' if r.topology_ok else 'MISMATCH'} | "
                f"{'ok' if r.orientation_ok else 'FLIPPED'} | "
                f"**{r.status}** |"
            )

    pass_count = sum(1 for r in results if r.status == "PASS")
    na_count = sum(1 for r in results if r.status == "N/A")
    total_applicable = len(results) - na_count
    lines += [
        "",
        f"**Result: {pass_count}/{total_applicable} PASS** ({na_count} N/A for non-star-shaped torus + volume formats)",
        "",
        "## Analytic Volume Comparison",
        "",
        "| Geometry | Analytic V [m³] | Mesh V [m³] | Discretisation Error % |",
        "|----------|-----------------|-------------|------------------------|",
    ]
    seen = set()
    for r in results:
        if r.geometry not in seen and r.volume_analytic is not None:
            seen.add(r.geometry)
            disc_err = r.volume_vs_analytic_rel * 100 if r.volume_vs_analytic_rel else 0
            lines.append(f"| {r.geometry} | {r.volume_analytic:.9f} | {r.volume_orig:.9f} | {disc_err:.4f} % |")

    lines += [
        "",
        "## Notes",
        "- **cube** (edge=2.0): analytic V = 8.0 m³; mesh is exact.",
        "- **sphere** (r=1.0, icosphere subdiv=2): analytic V = 4π/3 ≈ 4.18879 m³.",
        "- **torus** (genus 1, major=1.0, minor=0.3): analytic V = 2π²Rr² ≈ 1.77653 m³. "
        "**Non-star-shaped**: surface formats (STL/OBJ/VTK/Gmsh) work directly; "
        "volume formats (SU2/CGNS) require tetrahedralisation which star tessellation refuses. "
        "This is correct behavior — a torus has no interior point that sees all facets.",
        "- **wigley_hull**: L=10 m, B=1 m, T=0.625 m; analytic V = 4LBT/9 = 25/9 ≈ 2.77778 m³.",
        "- **STL** stores float32 by format definition; meshio writes ASCII STL here, so precision is preserved.",
        "- **SU2/CGNS** carry volume cells (tetra); surface meshes are tetrahedralised first.",
        "- **CGNS SIDS conformance**: `tensorcad.mesh_cgns` writes full SIDS metadata (0 violations); "
        "meshio's CGNS writer omits 14 required nodes/attributes.",
        "",
    ]
    return "\n".join(lines)


def generate_figures(results: list[AuditResult], outdir: Path):
    """Generate matplotlib figures for visual evidence."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors

    outdir.mkdir(parents=True, exist_ok=True)

    # Figure 1: Volume error heatmap (only applicable cells)
    formats = list(FORMATS_ALL)
    geometries = ["cube", "sphere", "torus", "wigley_hull"]

    err_matrix = np.full((len(formats), len(geometries)), np.nan)
    for r in results:
        if r.status == "N/A":
            continue
        i = formats.index(r.format)
        j = geometries.index(r.geometry)
        err_matrix[i, j] = r.rel_err if r.rel_err > 0 else 1e-16

    fig, ax = plt.subplots(figsize=(10, 5))
    mask = np.isnan(err_matrix)
    err_display = np.where(mask, -20, np.log10(err_matrix + 1e-16))
    im = ax.imshow(err_display, aspect="auto", cmap="RdYlGn_r", vmin=-16, vmax=-3)
    ax.set_xticks(range(len(geometries)))
    ax.set_xticklabels(geometries, fontsize=10)
    ax.set_yticks(range(len(formats)))
    ax.set_yticklabels(formats, fontsize=10)
    ax.set_xlabel("Geometry", fontsize=12)
    ax.set_ylabel("Format", fontsize=12)
    ax.set_title("Round-Trip Volume Error (log₁₀ relative error; gray = N/A)", fontsize=13)
    cbar = plt.colorbar(im, ax=ax, label="log₁₀(rel. error)")
    # Mark N/A cells
    for i in range(len(formats)):
        for j in range(len(geometries)):
            if mask[i, j]:
                ax.text(j, i, "N/A", ha="center", va="center", fontsize=10, color="gray")
    plt.tight_layout()
    fig.savefig(outdir / "volume_error_heatmap.png", dpi=150)
    plt.close(fig)

    # Figure 2: Status grid
    status_matrix = np.full((len(formats), len(geometries)), -1)  # -1 = N/A
    for r in results:
        i = formats.index(r.format)
        j = geometries.index(r.geometry)
        if r.status == "PASS":
            status_matrix[i, j] = 1
        elif r.status == "FAIL":
            status_matrix[i, j] = 0
        else:
            status_matrix[i, j] = -1

    fig, ax = plt.subplots(figsize=(10, 5))
    cmap = mcolors.ListedColormap(["#888888", "#ff4444", "#44cc44"])
    im = ax.imshow(status_matrix, aspect="auto", cmap=cmap, vmin=-1, vmax=1)
    ax.set_xticks(range(len(geometries)))
    ax.set_xticklabels(geometries, fontsize=10)
    ax.set_yticks(range(len(formats)))
    ax.set_yticklabels(formats, fontsize=10)
    ax.set_xlabel("Geometry", fontsize=12)
    ax.set_ylabel("Format", fontsize=12)
    ax.set_title("Round-Trip Status Grid (green=PASS, red=FAIL, gray=N/A)", fontsize=13)
    for i in range(len(formats)):
        for j in range(len(geometries)):
            val = status_matrix[i, j]
            if val == 1:
                ax.text(j, i, "✓", ha="center", va="center", fontsize=16, color="white")
            elif val == 0:
                ax.text(j, i, "✗", ha="center", va="center", fontsize=16, color="white")
            else:
                ax.text(j, i, "—", ha="center", va="center", fontsize=14, color="white")
    plt.tight_layout()
    fig.savefig(outdir / "status_grid.png", dpi=150)
    plt.close(fig)

    # Figure 3: Discretisation error vs analytic
    fig, ax = plt.subplots(figsize=(8, 5))
    geom_names = ["cube", "sphere", "torus", "wigley_hull"]
    analytic_vols = [8.0, 4.0/3.0*math.pi, 2.0*math.pi**2*1.0*0.3**2, 4.0*10.0*1.0*0.625/9.0]
    mesh_vols = []
    for g in geom_names:
        for r in results:
            if r.geometry == g and r.format == "obj":
                mesh_vols.append(r.volume_orig)
                break

    x = np.arange(len(geom_names))
    width = 0.35
    bars1 = ax.bar(x - width/2, analytic_vols, width, label="Analytic", color="#4488cc", alpha=0.8)
    bars2 = ax.bar(x + width/2, mesh_vols, width, label="Mesh (discretised)", color="#cc8844", alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(geom_names, fontsize=10)
    ax.set_ylabel("Volume [m³]", fontsize=12)
    ax.set_title("Analytic vs Discretised Volume", fontsize=13)
    ax.legend(fontsize=10)
    for i, (av, mv) in enumerate(zip(analytic_vols, mesh_vols)):
        err = abs(mv - av) / av * 100
        ax.annotate(f"{err:.3f}%", xy=(i, max(av, mv)),
                    xytext=(0, 8), textcoords="offset points",
                    ha="center", fontsize=9, color="#333333")
    plt.tight_layout()
    fig.savefig(outdir / "analytic_vs_mesh_volume.png", dpi=150)
    plt.close(fig)

    # Figure 4: Torus visualization (genus-1 body)
    t = torus(nu=48, nv=24)
    fig = plt.figure(figsize=(10, 4))

    ax1 = fig.add_subplot(121, projection="3d")
    ax1.plot_trisurf(t.points[:, 0], t.points[:, 1], t.points[:, 2],
                     triangles=t.faces, cmap="viridis", alpha=0.8, edgecolor="none")
    ax1.set_title(f"Torus (genus 1)\nV = {t.volume():.6f} m³ (analytic: {2*math.pi**2*1.0*0.3**2:.6f})")
    ax1.set_xlabel("X"); ax1.set_ylabel("Y"); ax1.set_zlabel("Z")
    ax1.set_box_aspect([1, 1, 0.6])

    # Round-trip torus through STL
    tm = _trimesh_to_mesh(t)
    rt = round_trip_meshio(tm, file_format="stl")
    ax2 = fig.add_subplot(122, projection="3d")
    ax2.plot_trisurf(rt.vertices[:, 0], rt.vertices[:, 1], rt.vertices[:, 2],
                     triangles=rt.faces, cmap="viridis", alpha=0.8, edgecolor="none")
    ax2.set_title(f"Torus after STL round-trip\nV = {surface_signed_volume(rt):.6f} m³")
    ax2.set_xlabel("X"); ax2.set_ylabel("Y"); ax2.set_zlabel("Z")
    ax2.set_box_aspect([1, 1, 0.6])

    plt.tight_layout()
    fig.savefig(outdir / "torus_roundtrip.png", dpi=150)
    plt.close(fig)

    # Figure 5: Wigley hull visualization
    from tensorcad import wigley
    hull = wigley.generate_wigley()
    fig = plt.figure(figsize=(12, 4))

    ax1 = fig.add_subplot(121, projection="3d")
    ax1.plot_trisurf(hull.mesh.vertices[:, 0], hull.mesh.vertices[:, 1], hull.mesh.vertices[:, 2],
                     triangles=hull.mesh.faces, cmap="ocean", alpha=0.8, edgecolor="none")
    ax1.set_title(f"Wigley Hull (L=10, B=1, T=0.625)\nV = {surface_signed_volume(hull.mesh):.6f} m³")
    ax1.set_xlabel("X (bow→stern)"); ax1.set_ylabel("Y (port→starboard)"); ax1.set_zlabel("Z")
    ax1.set_box_aspect([2, 1, 0.6])

    # Round-trip through OBJ
    rt_hull = round_trip_meshio(hull.mesh, file_format="obj")
    ax2 = fig.add_subplot(122, projection="3d")
    ax2.plot_trisurf(rt_hull.vertices[:, 0], rt_hull.vertices[:, 1], rt_hull.vertices[:, 2],
                     triangles=rt_hull.faces, cmap="ocean", alpha=0.8, edgecolor="none")
    ax2.set_title(f"Wigley Hull after OBJ round-trip\nV = {surface_signed_volume(rt_hull):.6f} m³")
    ax2.set_xlabel("X (bow→stern)"); ax2.set_ylabel("Y (port→starboard)"); ax2.set_zlabel("Z")
    ax2.set_box_aspect([2, 1, 0.6])

    plt.tight_layout()
    fig.savefig(outdir / "wigley_hull_roundtrip.png", dpi=150)
    plt.close(fig)


def main():
    outdir = Path("runs/m21535_matrix24")
    outdir.mkdir(parents=True, exist_ok=True)

    print("Running expanded format matrix audit (6 formats × 4 geometries)...")
    results = run_audit()

    # Print summary
    pass_count = sum(1 for r in results if r.status == "PASS")
    na_count = sum(1 for r in results if r.status == "N/A")
    total_applicable = len(results) - na_count
    print(f"\nResult: {pass_count}/{total_applicable} PASS ({na_count} N/A)")
    for r in results:
        if r.status == "FAIL":
            print(f"  FAIL: {r.format} {r.geometry}: rel_err={r.rel_err:.2e}")
        elif r.status == "N/A":
            print(f"  N/A:  {r.format} {r.geometry}: {r.note[:60]}...")

    # Generate markdown table
    md = generate_markdown_table(results)
    (outdir / "matrix24.md").write_text(md)
    print(f"\nMarkdown table written to {outdir / 'matrix24.md'}")

    # Generate JSON
    json_data = [asdict(r) for r in results]
    (outdir / "matrix24.json").write_text(json.dumps(json_data, indent=2))
    print(f"JSON data written to {outdir / 'matrix24.json'}")

    # Generate figures
    print("\nGenerating figures...")
    generate_figures(results, outdir / "figures")
    print(f"Figures written to {outdir / 'figures'}/")

    print("\nDone.")


if __name__ == "__main__":
    main()
