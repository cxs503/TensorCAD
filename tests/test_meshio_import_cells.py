"""Import side of the meshio format matrix: files TensorCAD did not write.

The round-trip matrix (``tests/test_meshio_format_matrix.py``) only ever reads
back what TensorCAD itself wrote, so it can only prove self-consistency.  Real
imports look different:

* surface cages arrive as **quads** (structured hull panels, Salome/Pointwise
  Gmsh exports, quad OBJ) or **polygons**, not as triangles;
* the same geometry arrives in **several cell blocks** (Gmsh regions, CGNS
  zones, a boundary patch stored next to the volume cells).

Both cases were broken before this module existed:

* a six-quad unit cube -- the simplest closed body -- raised
  ``MeshioFormatError: meshio mesh contains no triangle cells`` for OBJ, VTK
  and Gmsh;
* a file with a quad cage plus one stray triangle block silently returned **one
  triangle** (volume wrong, no error); a volume file with two tetra blocks
  silently returned only the first zone.

The tests below assert the invariants that matter for a geometry foundation,
not merely that the call returns: signed volume against the analytic value,
closure (watertight), manifoldness, and positive orientation.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

meshio = pytest.importorskip("meshio")

from tensorcad import mesh as meshlib
from tensorcad.geometry import Mesh
from tensorcad.meshio_exchange import (
    SURFACE_CELL_TYPES,
    VOLUME_CELL_TYPES,
    MeshioFormatError,
    TetrahedralMesh,
    meshio_to_mesh,
    meshio_to_tetrahedral,
    read_meshio,
    surface_signed_volume,
    tetrahedral_volume,
    write_meshio,
)

# ---------------------------------------------------------------------------
# Quad / polygon fixtures (exact analytic volumes)
# ---------------------------------------------------------------------------

#: Outward-oriented quad cage of the unit cube (volume 1, area 6).
CUBE_POINTS = np.array(
    [
        [0.0, 0.0, 0.0],  # 0
        [1.0, 0.0, 0.0],  # 1
        [1.0, 1.0, 0.0],  # 2
        [0.0, 1.0, 0.0],  # 3
        [0.0, 0.0, 1.0],  # 4
        [1.0, 0.0, 1.0],  # 5
        [1.0, 1.0, 1.0],  # 6
        [0.0, 1.0, 1.0],  # 7
    ],
    dtype=float,
)
CUBE_QUADS = np.array(
    [
        [0, 3, 2, 1],  # bottom z=0
        [4, 5, 6, 7],  # top    z=1
        [0, 1, 5, 4],  # front  y=0
        [2, 3, 7, 6],  # back   y=1
        [1, 2, 6, 5],  # right  x=1
        [0, 4, 7, 3],  # left   x=0
    ],
    dtype=np.int64,
)


def _pentagonal_prism() -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Closed pentagonal prism: two pentagon caps + five quad sides.

    Regular pentagon of circumradius R=1, height H=1, so the exact volume is
    ``(5/2) R^2 sin(72 deg) * H``.  Mixed cell types in separate blocks is
    exactly what CAD tessellations and Gmsh region exports look like, and the
    n-gon caps force the polygon fan path.
    """
    radius, height = 1.0, 1.0
    angles = [2.0 * math.pi * k / 5 for k in range(5)]
    points = np.array(
        [(radius * math.cos(a), radius * math.sin(a), 0.0) for a in angles]
        + [(radius * math.cos(a), radius * math.sin(a), height) for a in angles]
    )
    pentagons = np.array(
        [
            list(range(4, -1, -1)),  # bottom cap, reversed -> outward (-z)
            list(range(5, 10)),      # top cap, CCW from above -> outward (+z)
        ]
    )
    quads = np.array([(i, (i + 1) % 5, 5 + (i + 1) % 5, 5 + i) for i in range(5)])
    exact = 2.5 * radius**2 * math.sin(2.0 * math.pi / 5) * height
    return points, pentagons, quads, exact


def _assert_closed_outward(mesh: Mesh, expected_volume: float, tol: float = 1e-12) -> None:
    assert meshlib.is_closed(mesh) if hasattr(meshlib, "is_closed") else meshlib.is_watertight(mesh)
    assert meshlib.is_edge_manifold(mesh)
    assert meshlib.zero_area_faces(mesh) == ()
    volume = surface_signed_volume(mesh)
    assert volume > 0.0, "imported cage must stay outward-oriented"
    assert abs(volume - expected_volume) <= tol, (volume, expected_volume)


# ---------------------------------------------------------------------------
# 1. Quad cages: the classic "structured panel" import
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("file_format,ext", [("obj", ".obj"), ("vtk", ".vtk"), ("gmsh", ".msh")])
def test_quad_cube_imports_from_file(tmp_path, file_format: str, ext: str) -> None:
    """A quad cage written by meshio (i.e. by another tool) reads back as triangles."""
    path = write_meshio_quad(tmp_path / f"cube{ext}", file_format)
    mesh = read_meshio(path, file_format=file_format)
    assert isinstance(mesh, Mesh)
    assert len(mesh.faces) == 12, "six quads must split into 12 triangles"
    assert len(mesh.vertices) == 8
    _assert_closed_outward(mesh, 1.0, tol=1e-9)


def write_meshio_quad(path, file_format: str):
    mio = meshio.Mesh(CUBE_POINTS, [("quad", CUBE_QUADS)])
    meshio.write(str(path), mio, file_format=file_format)
    return path


def test_quad_cube_in_memory() -> None:
    mio = meshio.Mesh(CUBE_POINTS, [("quad", CUBE_QUADS)])
    mesh = meshio_to_mesh(mio)
    assert len(mesh.faces) == 12
    _assert_closed_outward(mesh, 1.0)


def test_polygon_and_quad_cage_in_memory() -> None:
    """n-gon caps plus quad sides (CAD tessellation topology) are triangulated."""
    points, pentagons, quads, exact = _pentagonal_prism()
    mio = meshio.Mesh(points, [("polygon", pentagons), ("quad", quads)])
    mesh = meshio_to_mesh(mio)
    assert all(len(face) == 3 for face in mesh.faces)
    assert len(mesh.faces) == 2 * 3 + 5 * 2, "2 pentagons fan to 3 each, 5 quads to 2 each"
    _assert_closed_outward(mesh, exact, tol=1e-15)


def test_quad_winding_preserved() -> None:
    """Triangulation must not flip a facet: face normals stay outward."""
    mio = meshio.Mesh(CUBE_POINTS, [("quad", CUBE_QUADS)])
    mesh = meshio_to_mesh(mio)
    for face in mesh.faces:
        normal = meshlib.face_normal(mesh, face)
        centroid = meshlib.face_centroid(mesh, face)
        centre = (0.5, 0.5, 0.5)
        outward = tuple(centroid[k] - centre[k] for k in range(3))
        assert sum(normal[k] * outward[k] for k in range(3)) > 0.0


def test_degenerate_quad_triangle_dropped() -> None:
    """A quad collapsed onto a repeated index contributes no zero-area triangle."""
    quads = np.array([[0, 0, 1, 2]])
    mio = meshio.Mesh(
        np.array([[0.0, 0, 0], [1.0, 0, 0], [1.0, 1.0, 0]]), [("quad", quads)]
    )
    mesh = meshio_to_mesh(mio)
    assert mesh.faces == ((0, 1, 2),)


# ---------------------------------------------------------------------------
# 2. Multi-block files: no silent geometry loss
# ---------------------------------------------------------------------------

def test_mixed_blocks_all_consumed() -> None:
    """Quad cage + separate triangle block: every facet must survive."""
    extra = np.array([[0, 1, 2]])
    mio = meshio.Mesh(CUBE_POINTS, [("quad", CUBE_QUADS), ("triangle", extra)])
    mesh = meshio_to_mesh(mio)
    assert len(mesh.faces) == 13, "the stray triangle block must not be dropped"
    # the 12 cube triangles are still there and still close the cube
    cube_only = Mesh(vertices=mesh.vertices, faces=mesh.faces[:12])
    _assert_closed_outward(cube_only, 1.0)


def test_two_triangle_blocks_concatenated() -> None:
    """meshio may split one triangle set into several blocks (e.g. per entity)."""
    first = np.array([[0, 3, 2], [0, 2, 1], [4, 5, 6]])
    second = np.array([[4, 6, 7], [0, 1, 5], [0, 5, 4]])
    third = np.array([[2, 3, 7], [2, 7, 6], [1, 2, 6], [1, 6, 5], [0, 4, 7], [0, 7, 3]])
    mio = meshio.Mesh(CUBE_POINTS, [("triangle", first), ("triangle", second), ("triangle", third)])
    mesh = meshio_to_mesh(mio)
    assert len(mesh.faces) == 12
    _assert_closed_outward(mesh, 1.0)


def test_two_tetra_blocks_concatenated() -> None:
    """Multi-zone volume files: both zones must survive."""
    points = np.array(
        [
            [0.0, 0, 0], [1.0, 0, 0], [0.0, 1.0, 0], [0.0, 0, 1.0],
            [2.0, 0, 0], [3.0, 0, 0], [2.0, 1.0, 0], [2.0, 0, 1.0],
        ]
    )
    zone1 = np.array([[0, 1, 2, 3]])
    zone2 = np.array([[4, 5, 6, 7]])
    mio = meshio.Mesh(points, [("tetra", zone1), ("tetra", zone2)])
    tet_mesh = meshio_to_tetrahedral(mio)
    assert isinstance(tet_mesh, TetrahedralMesh)
    assert len(tet_mesh.tets) == 2
    total = sum(tetrahedral_volume(tet_mesh.vertices, tet) for tet in tet_mesh.tets)
    assert math.isclose(total, 1.0 / 6.0 + 1.0 / 6.0, rel_tol=1e-12)


def test_boundary_line_block_ignored() -> None:
    """Boundary ``line`` sets are not geometry: they must not corrupt the surface."""
    lines = np.array([[0, 1], [1, 2]])
    mio = meshio.Mesh(
        CUBE_POINTS,
        [("quad", CUBE_QUADS), ("line", lines)],
    )
    mesh = meshio_to_mesh(mio)
    assert len(mesh.faces) == 12
    _assert_closed_outward(mesh, 1.0)


# ---------------------------------------------------------------------------
# 3. Refusals: wrong topology must be diagnosed, never silently accepted
# ---------------------------------------------------------------------------

def test_volume_only_file_imported_as_surface_is_diagnosed() -> None:
    points = np.array([[0.0, 0, 0], [1.0, 0, 0], [0.0, 1.0, 0], [0.0, 0, 1.0]])
    mio = meshio.Mesh(points, [("tetra", np.array([[0, 1, 2, 3]]))])
    with pytest.raises(MeshioFormatError) as excinfo:
        meshio_to_mesh(mio)
    message = str(excinfo.value)
    assert "tetra" in message and "meshio_to_tetrahedral" in message


def test_no_cells_at_all_is_diagnosed() -> None:
    mio = meshio.Mesh(CUBE_POINTS.copy(), [("vertex", np.arange(8).reshape(-1, 1))])
    with pytest.raises(MeshioFormatError) as excinfo:
        meshio_to_mesh(mio)
    assert "vertex" in str(excinfo.value)


def test_hexahedra_are_refused_with_a_reason() -> None:
    points = CUBE_POINTS.copy()
    hexes = np.array([[0, 1, 2, 3, 4, 5, 6, 7]])
    mio = meshio.Mesh(points, [("hexahedron", hexes)])
    with pytest.raises(MeshioFormatError) as excinfo:
        meshio_to_tetrahedral(mio)
    message = str(excinfo.value)
    assert "no tetra cells" in message and "hexahedron" in message


def test_out_of_range_index_rejected() -> None:
    mio = meshio.Mesh(CUBE_POINTS[:3].copy(), [("triangle", np.array([[0, 1, 7]]))])
    with pytest.raises(MeshioFormatError) as excinfo:
        meshio_to_mesh(mio)
    assert "outside" in str(excinfo.value) or "references vertex" in str(excinfo.value)


# ---------------------------------------------------------------------------
# 4. Contract surface
# ---------------------------------------------------------------------------

def test_surface_cell_types_are_the_documented_set() -> None:
    assert SURFACE_CELL_TYPES == ("triangle", "quad", "polygon")
    assert "tetra" in VOLUME_CELL_TYPES


def test_quad_volume_matches_triangle_volume_of_same_body() -> None:
    """Both representations of the same closed body must agree to machine precision."""
    tri_mesh = meshio_to_mesh(meshio.Mesh(CUBE_POINTS, [("quad", CUBE_QUADS)]))
    tri_mio = meshio.Mesh(
        CUBE_POINTS,
        [
            ("triangle", np.array([[0, 3, 2], [0, 2, 1], [4, 5, 6], [4, 6, 7],
                                   [0, 1, 5], [0, 5, 4], [2, 3, 7], [2, 7, 6],
                                   [1, 2, 6], [1, 6, 5], [0, 4, 7], [0, 7, 3]]))
        ],
    )
    ref = meshio_to_mesh(tri_mio)
    assert surface_signed_volume(tri_mesh) == pytest.approx(surface_signed_volume(ref), abs=1e-15)