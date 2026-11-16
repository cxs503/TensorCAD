"""Unified mesh import: auto-detect + read all 9 dialects.

Tests the new :mod:`tensorcad.mesh_import` module against the 6-format matrix
geometries, verifying that:

1. Auto-detection works for every supported extension
2. Every format round-trips with volume/topology/orientation preserved
3. Foreign dialects (binary STL, VTK XML, Gmsh 2.2) are handled correctly
4. The import result metadata is accurate
5. Unknown extensions raise ``UnknownFormatError``
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.mesh_import import (
    ImportResult,
    UnknownFormatError,
    detect_format,
    read_mesh,
    SUPPORTED_FORMATS,
    FORMAT_EXTENSIONS,
    FOREIGN_EXTENSIONS,
)
from tensorcad.meshio_exchange import (
    FORMATS_SURFACE,
    FORMATS_VOLUME,
    TetrahedralMesh,
    surface_signed_volume,
    tetrahedralise_surface,
    write_meshio,
)
from tensorcad.mesh_foreign_dialects import (
    write_stl_binary,
    write_gmsh22,
    write_vtu_ascii,
)

# ---------------------------------------------------------------------------
# Test fixtures
# ---------------------------------------------------------------------------

def _make_cube() -> Mesh:
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


def _make_sphere(n_lat: int = 24, n_lon: int = 48) -> Mesh:
    """Unit sphere: pole fans + quad bands, oriented outward, watertight."""
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
                    (math.sin(phi) * math.cos(theta), math.sin(phi) * math.sin(theta), math.cos(phi))
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

    from tensorcad import mesh as meshlib
    return meshlib.orient_outward(builder.build()).mesh


def _make_wigley_hull():
    """Real benchmark hull: tensorcad.wigley's closed, deck-capped Wigley."""
    from tensorcad import wigley
    return wigley.generate_wigley().mesh


SURFACE_GEOMETRIES = {
    "cube": _make_cube(),
    "sphere": _make_sphere(),
    "wigley_hull": _make_wigley_hull(),
}

ANALYTIC_VOLUME = {
    "cube": 1.0,
    "sphere": 4.0 / 3.0 * math.pi,
    "wigley_hull": 4.0 * 10.0 * 1.0 * 0.625 / 9.0,
}


# ---------------------------------------------------------------------------
# Auto-detection tests
# ---------------------------------------------------------------------------

class TestDetectFormat:
    """Format auto-detection from file extension."""

    def test_standard_extensions(self, tmp_path: Path) -> None:
        for ext, fmt in FORMAT_EXTENSIONS.items():
            p = tmp_path / f"test{ext}"
            p.write_bytes(b"dummy")
            if ext == ".stl":
                # STL needs content to distinguish ASCII vs binary
                p.write_bytes(b"solid test\n")
                assert detect_format(p) == "stl"
            else:
                assert detect_format(p) == fmt

    def test_foreign_extensions(self, tmp_path: Path) -> None:
        for ext, fmt in FOREIGN_EXTENSIONS.items():
            if ext == ".stl":
                continue  # handled above
            p = tmp_path / f"test{ext}"
            p.write_bytes(b"dummy")
            assert detect_format(p) == fmt

    def test_binary_stl_detection(self, tmp_path: Path) -> None:
        """Binary STL does not start with 'solid'."""
        p = tmp_path / "test.stl"
        # 80-byte header + uint32 count + one facet
        header = b"\x00" * 80
        count = (1).to_bytes(4, "little")
        facet = b"\x00" * 50
        p.write_bytes(header + count + facet)
        assert detect_format(p) == "stl_binary"

    def test_ascii_stl_detection(self, tmp_path: Path) -> None:
        p = tmp_path / "test.stl"
        p.write_text("solid test\nendsolid test\n")
        assert detect_format(p) == "stl"

    def test_unknown_extension(self, tmp_path: Path) -> None:
        p = tmp_path / "test.xyz"
        p.write_bytes(b"dummy")
        with pytest.raises(UnknownFormatError, match="unknown extension"):
            detect_format(p)

    def test_case_insensitive(self, tmp_path: Path) -> None:
        p = tmp_path / "test.OBJ"
        p.write_bytes(b"dummy")
        assert detect_format(p) == "obj"


# ---------------------------------------------------------------------------
# Surface format round-trip tests
# ---------------------------------------------------------------------------

SURFACE_FORMATS = ("stl", "obj", "vtk", "gmsh")


@pytest.mark.parametrize("fmt", SURFACE_FORMATS)
@pytest.mark.parametrize("geom_name", SURFACE_GEOMETRIES.keys())
def test_surface_format_round_trip(fmt: str, geom_name: str, tmp_path: Path) -> None:
    """Write a surface mesh, read it back via read_mesh, verify invariants."""
    geom = SURFACE_GEOMETRIES[geom_name]
    ext = {"stl": ".stl", "obj": ".obj", "vtk": ".vtk", "gmsh": ".msh"}[fmt]
    path = tmp_path / f"{geom_name}{ext}"

    write_meshio(path, geom, file_format=fmt)
    result = read_mesh(path)

    assert isinstance(result.mesh, Mesh)
    assert result.cell_kind == "tri"
    assert result.file_format == fmt
    assert result.orientation == "outward"
    assert result.cell_count == len(geom.faces)

    # Volume preserved
    ref_vol = surface_signed_volume(geom)
    assert abs(result.volume - ref_vol) / max(abs(ref_vol), 1e-300) < 1e-9


# ---------------------------------------------------------------------------
# Volume format round-trip tests
# ---------------------------------------------------------------------------

VOLUME_FORMATS = ("su2", "cgns")


@pytest.mark.parametrize("fmt", VOLUME_FORMATS)
@pytest.mark.parametrize("geom_name", ("cube", "sphere"))
def test_volume_format_round_trip(fmt: str, geom_name: str, tmp_path: Path) -> None:
    """Write a tetrahedral mesh, read it back via read_mesh, verify invariants."""
    geom = SURFACE_GEOMETRIES[geom_name]
    tet = tetrahedralise_surface(geom)
    ext = {"su2": ".su2", "cgns": ".cgns"}[fmt]
    path = tmp_path / f"{geom_name}{ext}"

    write_meshio(path, tet, file_format=fmt)
    result = read_mesh(path)

    assert isinstance(result.mesh, TetrahedralMesh)
    assert result.cell_kind == "tetra"
    assert result.file_format == fmt
    assert result.cell_count == len(tet.tets)

    # Volume preserved
    ref_vol = tet.volume()
    assert abs(result.volume - ref_vol) / max(abs(ref_vol), 1e-300) < 1e-9


# ---------------------------------------------------------------------------
# Foreign dialect tests
# ---------------------------------------------------------------------------

FOREIGN_SURFACE = ("stl_binary", "vtu", "gmsh22")


@pytest.mark.parametrize("dialect", FOREIGN_SURFACE)
@pytest.mark.parametrize("geom_name", ("cube", "sphere"))
def test_foreign_surface_round_trip(dialect: str, geom_name: str, tmp_path: Path) -> None:
    """Write with spec writer, read via read_mesh (auto-detect), verify."""
    geom = SURFACE_GEOMETRIES[geom_name]
    ext = {"stl_binary": ".stl", "vtu": ".vtu", "gmsh22": ".msh2"}[dialect]
    path = tmp_path / f"{geom_name}{ext}"

    writers = {
        "stl_binary": write_stl_binary,
        "vtu": write_vtu_ascii,
        "gmsh22": write_gmsh22,
    }
    writers[dialect](path, geom)
    result = read_mesh(path)

    assert isinstance(result.mesh, Mesh)
    assert result.cell_kind == "tri"
    assert result.dialect == dialect
    assert result.orientation == "outward"

    ref_vol = surface_signed_volume(geom)
    tol = 1e-6 if dialect == "stl_binary" else 1e-12
    assert abs(result.volume - ref_vol) / max(abs(ref_vol), 1e-300) < tol


def test_foreign_tetra_round_trip(tmp_path: Path) -> None:
    """Gmsh 2.2 can carry tetrahedra too."""
    geom = SURFACE_GEOMETRIES["cube"]
    tet = tetrahedralise_surface(geom)
    path = tmp_path / "cube_tetra.msh2"

    write_gmsh22(path, tet)
    result = read_mesh(path)

    assert isinstance(result.mesh, TetrahedralMesh)
    assert result.cell_kind == "tetra"
    assert result.dialect == "gmsh22"
    assert abs(result.volume - tet.volume()) / tet.volume() < 1e-12


# ---------------------------------------------------------------------------
# Explicit format override
# ---------------------------------------------------------------------------

def test_explicit_format_override(tmp_path: Path) -> None:
    """Caller can override auto-detection."""
    geom = SURFACE_GEOMETRIES["cube"]
    path = tmp_path / "cube.obj"
    write_meshio(path, geom, file_format="obj")

    # Auto-detect
    result_auto = read_mesh(path)
    assert result_auto.file_format == "obj"

    # Explicit override
    result_explicit = read_mesh(path, file_format="obj")
    assert result_explicit.file_format == "obj"
    assert result_explicit.cell_count == result_auto.cell_count


# ---------------------------------------------------------------------------
# Error handling
# ---------------------------------------------------------------------------

def test_file_not_found() -> None:
    with pytest.raises(FileNotFoundError):
        read_mesh("/nonexistent/path/cube.obj")


def test_unsupported_format(tmp_path: Path) -> None:
    p = tmp_path / "test.xyz"
    p.write_bytes(b"dummy")
    with pytest.raises(UnknownFormatError):
        read_mesh(p)


# ---------------------------------------------------------------------------
# ImportResult metadata
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("geom_name", SURFACE_GEOMETRIES.keys())
def test_import_result_metadata(geom_name: str, tmp_path: Path) -> None:
    """ImportResult carries accurate metadata."""
    geom = SURFACE_GEOMETRIES[geom_name]
    path = tmp_path / f"{geom_name}.obj"
    write_meshio(path, geom, file_format="obj")

    result = read_mesh(path)

    assert result.source_path == str(path)
    assert result.vertex_count == len(geom.vertices)
    assert result.cell_count == len(geom.faces)
    assert result.cell_kind == "tri"
    assert result.orientation == "outward"

    # Volume matches analytic (within discretisation for curved bodies)
    analytic = ANALYTIC_VOLUME.get(geom_name)
    if analytic is not None:
        rel_err = abs(result.volume - analytic) / analytic
        limit = 0.0 if geom_name == "cube" else 0.01
        assert rel_err <= limit, f"{geom_name}: volume error {rel_err:.3e} > {limit}"


# ---------------------------------------------------------------------------
# Cross-format volume agreement
# ---------------------------------------------------------------------------

def test_all_formats_agree_on_cube_volume(tmp_path: Path) -> None:
    """All 9 dialects must agree on the cube's volume (exactly 1.0)."""
    geom = SURFACE_GEOMETRIES["cube"]
    tet = tetrahedralise_surface(geom)

    volumes: dict[str, float] = {}

    # Surface formats
    for fmt in SURFACE_FORMATS:
        ext = {"stl": ".stl", "obj": ".obj", "vtk": ".vtk", "gmsh": ".msh"}[fmt]
        path = tmp_path / f"cube{ext}"
        write_meshio(path, geom, file_format=fmt)
        result = read_mesh(path)
        volumes[fmt] = result.volume

    # Volume formats
    for fmt in VOLUME_FORMATS:
        ext = {"su2": ".su2", "cgns": ".cgns"}[fmt]
        path = tmp_path / f"cube{ext}"
        write_meshio(path, tet, file_format=fmt)
        result = read_mesh(path)
        volumes[fmt] = result.volume

    # Foreign dialects
    write_stl_binary(tmp_path / "cube.stl", geom)
    volumes["stl_binary"] = read_mesh(tmp_path / "cube.stl").volume

    write_vtu_ascii(tmp_path / "cube.vtu", geom)
    volumes["vtu"] = read_mesh(tmp_path / "cube.vtu").volume

    write_gmsh22(tmp_path / "cube.msh2", geom)
    volumes["gmsh22"] = read_mesh(tmp_path / "cube.msh2").volume

    # All must be exactly 1.0 (within float64 round-off)
    for fmt, vol in volumes.items():
        assert abs(vol - 1.0) < 1e-12, f"{fmt}: volume {vol} != 1.0"


# ---------------------------------------------------------------------------
# Supported formats list
# ---------------------------------------------------------------------------

def test_supported_formats_complete() -> None:
    """SUPPORTED_FORMATS includes all standard + foreign formats."""
    expected = set(FORMAT_EXTENSIONS.values()) | set(FOREIGN_EXTENSIONS.values())
    # Remove duplicate .stl entry (both ASCII and binary map to .stl)
    expected.discard("stl_binary")  # .stl maps to both
    assert set(SUPPORTED_FORMATS) >= {"stl", "obj", "vtk", "gmsh", "su2", "cgns", "stl_binary", "vtu", "gmsh22"}
