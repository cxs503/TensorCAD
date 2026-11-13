"""meshio format matrix: round-trip fidelity for 6 formats × 3 geometries.

Verifies that exporting a surface mesh through meshio and reading it back
preserves:
- **Volume** (divergence-theorem signed volume, relative error < 0.1%)
- **Topology** (vertex count, face count; exact match for OBJ/VTK/Gmsh)
- **Orientation** (signed volume sign unchanged → outward normals preserved)

Format capability
=================
STL / OBJ / VTK / Gmsh → REAL (surface triangles supported)
SU2 / CGNS             → NOT_SUPPORTED for surface meshes (volume-only writers)

The 12 supported combinations (4 formats × 3 geometries) are tested;
the 4 unsupported ones (SU2/CGNS × 3 geometries) are verified to raise
``MeshioFormatError`` with a clear message.
"""
from __future__ import annotations

import math
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pytest

from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad import mesh as meshlib
from tensorcad.meshio_exchange import (
    FORMATS_SURFACE,
    FORMATS_VOLUME,
    MeshioFormatError,
    mesh_to_meshio,
    meshio_to_mesh,
    round_trip_meshio,
)

# ---------------------------------------------------------------------------
# Test geometries
# ---------------------------------------------------------------------------

def _make_cube() -> Mesh:
    """Unit cube, outward-oriented, 8 vertices, 12 triangles."""
    vertices = (
        (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0),
    )
    faces = (
        (0, 3, 2), (0, 2, 1),  # bottom z=0
        (4, 5, 6), (4, 6, 7),  # top    z=1
        (0, 1, 5), (0, 5, 4),  # front  y=0
        (2, 3, 7), (2, 7, 6),  # back   y=1
        (1, 2, 6), (1, 6, 5),  # right  x=1
        (0, 4, 7), (0, 7, 3),  # left   x=0
    )
    return Mesh(vertices=vertices, faces=faces, units="m")


def _make_sphere(n_lat: int = 16, n_lon: int = 32) -> Mesh:
    """Approximate unit sphere via latitude-longitude tessellation.

    Analytic volume = 4/3·π ≈ 4.18879.
    """
    builder = MeshBuilder(units="m")
    # Build vertices
    vertex_map: list[list[int]] = []
    for i_lat in range(n_lat + 1):
        phi = math.pi * i_lat / n_lat  # 0..pi
        row: list[int] = []
        for i_lon in range(n_lon):
            theta = 2 * math.pi * i_lon / n_lon  # 0..2pi
            x = math.sin(phi) * math.cos(theta)
            y = math.sin(phi) * math.sin(theta)
            z = math.cos(phi)
            row.append(builder.add_vertex((x, y, z)))
        vertex_map.append(row)

    # Build triangles
    for i_lat in range(n_lat):
        for i_lon in range(n_lon):
            i_lon_next = (i_lon + 1) % n_lon
            a = vertex_map[i_lat][i_lon]
            b = vertex_map[i_lat][i_lon_next]
            c = vertex_map[i_lat + 1][i_lon_next]
            d = vertex_map[i_lat + 1][i_lon]
            builder.add_quad(a, b, c, d)

    return builder.build()


def _make_wigley() -> Mesh:
    """Wigley III hull form (parabolic waterlines and sections).

    L=2, B=0.3, D=0.125.  Analytic volume ≈ 0.075 (Cb≈0.5 for Wigley).
    """
    builder = MeshBuilder(units="m")
    nx, ny, nz = 20, 10, 8
    L, B, D = 2.0, 0.3, 0.125

    def half_beam(x: float, z: float) -> float:
        """Wigley III half-beam at station x and height z."""
        xi = x / (L / 2)  # -1..1
        zeta = z / D  # 0..1
        return (B / 2) * (1 - xi**2) * (1 - zeta**2) * (1 + 0.2 * xi**2)

    # Build grid vertices
    v: list[list[list[int]]] = []
    for iz in range(nz + 1):
        z = D * iz / nz
        row: list[list[int]] = []
        for ix in range(nx + 1):
            x = -L / 2 + L * ix / nx
            yb = half_beam(x, z)
            # port (-y) and starboard (+y)
            vp = builder.add_vertex((x, -yb, z))
            vs = builder.add_vertex((x, yb, z))
            row.append((vp, vs))
        v.append(row)

    # Build quads on each side
    for iz in range(nz):
        for ix in range(nx):
            # starboard side
            a = v[iz][ix][1]
            b = v[iz][ix + 1][1]
            c = v[iz + 1][ix + 1][1]
            d = v[iz + 1][ix][1]
            builder.add_quad(a, b, c, d)
            # port side (reverse order for outward normal)
            a = v[iz][ix][0]
            b = v[iz + 1][ix][0]
            c = v[iz + 1][ix + 1][0]
            d = v[iz][ix + 1][0]
            builder.add_quad(a, b, c, d)

    # Bottom (z=0)
    for ix in range(nx):
        a = v[0][ix][0]
        b = v[0][ix][1]
        c = v[0][ix + 1][1]
        d = v[0][ix + 1][0]
        builder.add_quad(a, b, c, d)

    return builder.build()


@dataclass(frozen=True)
class GeometryCase:
    name: str
    mesh: Mesh
    expected_volume: float  # approximate, for sanity check


GEOMETRIES: tuple[GeometryCase, ...] = (
    GeometryCase("cube", _make_cube(), 1.0),
    GeometryCase("sphere", _make_sphere(), 4.0 / 3.0 * math.pi),
    GeometryCase("wigley", _make_wigley(), 0.075),
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _signed_volume(mesh: Mesh) -> float:
    return meshlib.total_signed_volume(mesh)


def _relative_error(a: float, b: float) -> float:
    if b == 0:
        return float("inf") if a != 0 else 0.0
    return abs(a - b) / abs(b)


# ---------------------------------------------------------------------------
# Tests: supported formats (STL, OBJ, VTK, Gmsh)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("geom", GEOMETRIES, ids=lambda g: g.name)
@pytest.mark.parametrize("fmt", FORMATS_SURFACE)
def test_round_trip_volume_preserved(fmt: str, geom: GeometryCase) -> None:
    """Round-trip volume must match original within 0.1%."""
    original = geom.mesh
    rt = round_trip_meshio(original, file_format=fmt)
    vol_orig = _signed_volume(original)
    vol_rt = _signed_volume(rt)
    rel_err = _relative_error(vol_rt, vol_orig)
    assert rel_err < 1e-3, (
        f"{fmt}/{geom.name}: volume {vol_rt} vs {vol_orig} "
        f"(rel_err={rel_err:.2e}, limit=1e-3)"
    )


@pytest.mark.parametrize("geom", GEOMETRIES, ids=lambda g: g.name)
@pytest.mark.parametrize("fmt", FORMATS_SURFACE)
def test_round_trip_topology_preserved(fmt: str, geom: GeometryCase) -> None:
    """Vertex count and face count must be preserved."""
    original = geom.mesh
    rt = round_trip_meshio(original, file_format=fmt)
    assert len(rt.vertices) == len(original.vertices), (
        f"{fmt}/{geom.name}: vertex count {len(rt.vertices)} vs {len(original.vertices)}"
    )
    assert len(rt.faces) == len(original.faces), (
        f"{fmt}/{geom.name}: face count {len(rt.faces)} vs {len(original.faces)}"
    )


@pytest.mark.parametrize("geom", GEOMETRIES, ids=lambda g: g.name)
@pytest.mark.parametrize("fmt", FORMATS_SURFACE)
def test_round_trip_orientation_preserved(fmt: str, geom: GeometryCase) -> None:
    """Signed volume sign must be unchanged (outward normals preserved)."""
    original = geom.mesh
    rt = round_trip_meshio(original, file_format=fmt)
    vol_orig = _signed_volume(original)
    vol_rt = _signed_volume(rt)
    assert math.copysign(1, vol_rt) == math.copysign(1, vol_orig), (
        f"{fmt}/{geom.name}: orientation flipped "
        f"(vol_orig={vol_orig}, vol_rt={vol_rt})"
    )


@pytest.mark.parametrize("geom", GEOMETRIES, ids=lambda g: g.name)
@pytest.mark.parametrize("fmt", ("obj", "vtk", "gmsh"))
def test_exact_topology_for_lossless_formats(fmt: str, geom: GeometryCase) -> None:
    """OBJ, VTK, and Gmsh must preserve exact vertex positions and face indices."""
    original = geom.mesh
    rt = round_trip_meshio(original, file_format=fmt)
    # Vertex positions must match exactly
    for i, (v_orig, v_rt) in enumerate(zip(original.vertices, rt.vertices)):
        assert v_orig == pytest.approx(v_rt, abs=1e-12), (
            f"{fmt}/{geom.name}: vertex {i} mismatch {v_orig} vs {v_rt}"
        )
    # Face indices must match exactly
    for i, (f_orig, f_rt) in enumerate(zip(original.faces, rt.faces)):
        assert f_orig == f_rt, (
            f"{fmt}/{geom.name}: face {i} mismatch {f_orig} vs {f_rt}"
        )


# ---------------------------------------------------------------------------
# Tests: unsupported formats (SU2, CGNS)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("geom", GEOMETRIES, ids=lambda g: g.name)
@pytest.mark.parametrize("fmt", FORMATS_VOLUME)
def test_volume_formats_reject_surface_mesh(fmt: str, geom: GeometryCase) -> None:
    """SU2 and CGNS must raise MeshioFormatError for surface triangle meshes."""
    with pytest.raises(MeshioFormatError, match="volume cells"):
        round_trip_meshio(geom.mesh, file_format=fmt)


# ---------------------------------------------------------------------------
# Tests: meshio conversion functions
# ---------------------------------------------------------------------------

def test_mesh_to_meshio_preserves_data() -> None:
    mesh = _make_cube()
    mio = mesh_to_meshio(mesh)
    assert mio.points.shape == (8, 3)
    assert len(mio.cells) == 1
    assert mio.cells[0].type == "triangle"
    assert mio.cells[0].data.shape == (12, 3)


def test_meshio_to_mesh_roundtrip() -> None:
    mesh = _make_cube()
    mio = mesh_to_meshio(mesh)
    recovered = meshio_to_mesh(mio, units=mesh.units, coordinate_frame=mesh.coordinate_frame)
    assert recovered.vertices == mesh.vertices
    assert recovered.faces == mesh.faces
    assert recovered.units == mesh.units


# ---------------------------------------------------------------------------
# Summary table generation (for artifacts)
# ---------------------------------------------------------------------------

def generate_format_matrix_table() -> str:
    """Generate a markdown table of all 18 combinations."""
    lines = [
        "# meshio Format Matrix — Round-Trip Fidelity Report",
        "",
        "| Format | Geometry | Volume Preserved | Topology Preserved | Orientation Preserved | Status |",
        "|--------|----------|-----------------|-------------------|----------------------|--------|",
    ]

    for geom in GEOMETRIES:
        for fmt in FORMATS_SURFACE:
            original = geom.mesh
            rt = round_trip_meshio(original, file_format=fmt)
            vol_orig = _signed_volume(original)
            vol_rt = _signed_volume(rt)
            vol_ok = _relative_error(vol_rt, vol_orig) < 1e-3
            topo_ok = len(rt.vertices) == len(original.vertices) and len(rt.faces) == len(original.faces)
            orient_ok = math.copysign(1, vol_rt) == math.copysign(1, vol_orig)
            status = "✅ PASS" if (vol_ok and topo_ok and orient_ok) else "❌ FAIL"
            lines.append(
                f"| {fmt} | {geom.name} | {'✅' if vol_ok else '❌'} "
                f"({vol_rt:.6f} vs {vol_orig:.6f}) | "
                f"{'✅' if topo_ok else '❌'} "
                f"({len(rt.vertices)}v/{len(rt.faces)}f) | "
                f"{'✅' if orient_ok else '❌'} | "
                f"{status} |"
            )

    for geom in GEOMETRIES:
        for fmt in FORMATS_VOLUME:
            lines.append(
                f"| {fmt} | {geom.name} | N/A | N/A | N/A | "
                f"⚠️ NOT_SUPPORTED (volume-only format) |"
            )

    lines.append("")
    lines.append("## Notes")
    lines.append("- STL reorders vertices and faces but preserves volume and orientation.")
    lines.append("- OBJ, VTK, Gmsh preserve exact vertex positions and face indices.")
    lines.append("- SU2 and CGNS only support volume cells (tetra/hex/wedge/pyramid); surface meshes require tetrahedralisation first.")
    lines.append("")
    return "\n".join(lines)


if __name__ == "__main__":
    table = generate_format_matrix_table()
    print(table)
