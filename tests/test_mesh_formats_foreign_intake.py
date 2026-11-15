"""Foreign-intake parity for the native six-format matrix.

Anchor: ``GAP-CAD-B2-MESH-FORMAT-MATRIX`` — "meshio 全格式矩阵进出 (STL/OBJ/VTK/
SU2/CGNS/Gmsh)", acceptance "往返保真：体积/拓扑/取向不变".

The 18-cell round-trip matrix in :mod:`tests.test_mesh_formats` writes every file
with *this package's own writer*, so it proves writer and reader agree with each
other — the "出" half.  The "进" half is a different object: what actually arrives
from a CAD kernel, Pointwise or Salome is a **quad-dominant** cage, a planar
``NDIME=2`` SU2 file, a CGNS zone split over several ``Elements_t`` blocks, or a
file carrying stray line elements.  That axis used to be asymmetric: the
meshio-facing importer triangulated quads (``test_meshio_import_cells``) while
the native reader set — the one this module calls the product reader set —
crashed on a Gmsh quad cage with ``faces must be (m,3), got (0,)``, refused a
legal SU2 quad, and rejected planar SU2 outright.  These tests pin the fix, and
they pin the *policy*, not just the exit code:

* quads are fan-triangulated on the ``0-2`` diagonal, the same choice the meshio
  path makes, so one physical file reads the same through either entrance;
* nothing is dropped without trace: non-surface cells (lines, points) must be
  acknowledged through ``notes=[]`` or the read raises;
* content that is genuinely not a surface (SU2 volume cells, a file with no
  surface element at all) raises with the offending type/count quoted.
"""

from __future__ import annotations

import numpy as np
import pytest

from tensorcad import mesh_formats as mf

#: Unit cube, outward winding, as a 6-quad cage — the simplest closed foreign body.
QUAD_CAGE_POINTS = np.array(
    [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0],
        [0.0, 1.0, 1.0],
    ]
)

QUAD_CAGE_QUADS = [
    [0, 3, 2, 1],  # bottom, normal -z
    [4, 5, 6, 7],  # top, normal +z
    [0, 1, 5, 4],  # -y
    [1, 2, 6, 5],  # +x
    [2, 3, 7, 6],  # +y
    [3, 0, 4, 7],  # -x
]

#: The reference triangulation: fan from node 0 of each quad (diagonal 0-2).
_REF_FACES = np.array([face for q in QUAD_CAGE_QUADS for face in mf._fan(q)], dtype=np.int64)
_REF = mf.TriMesh(QUAD_CAGE_POINTS, _REF_FACES, "quad_cube_reference")

QUAD_FORMATS = ("obj", "vtk", "gmsh", "su2", "cgns")


def _write_obj(path, quads) -> None:
    body = "".join(f"v {x} {y} {z}\n" for x, y, z in QUAD_CAGE_POINTS)
    body += "".join("f " + " ".join(str(i + 1) for i in quad) + "\n" for quad in quads)
    path.write_text(body, encoding="utf-8")


def _write_vtk(path, quads) -> None:
    head = "# vtk DataFile Version 3.0\nquad cage\nASCII\nDATASET UNSTRUCTURED_GRID\n"
    head += f"POINTS {len(QUAD_CAGE_POINTS)} float\n"
    head += "".join(f"{x} {y} {z}\n" for x, y, z in QUAD_CAGE_POINTS)
    head += f"CELLS {len(quads)} {5 * len(quads)}\n"
    head += "".join("4 " + " ".join(map(str, quad)) + "\n" for quad in quads)
    head += f"CELL_TYPES {len(quads)}\n" + "9\n" * len(quads)
    path.write_text(head, encoding="utf-8")


def _write_gmsh(path, quads) -> None:
    body = "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n"
    body += f"$Nodes\n{len(QUAD_CAGE_POINTS)}\n"
    body += "".join(f"{i + 1} {x} {y} {z}\n" for i, (x, y, z) in enumerate(QUAD_CAGE_POINTS))
    body += "$EndNodes\n"
    body += f"$Elements\n{len(quads)}\n"
    body += "".join(
        f"{i + 1} 3 0 " + " ".join(str(t + 1) for t in quad) + "\n"
        for i, quad in enumerate(quads)
    )
    body += "$EndElements\n"
    path.write_text(body, encoding="utf-8")


def _write_su2(path, quads) -> None:
    body = "NDIME= 3\n"
    body += f"NPOIN= {len(QUAD_CAGE_POINTS)}\n"
    body += "".join(f"{x} {y} {z}\n" for x, y, z in QUAD_CAGE_POINTS)
    body += f"NELEM= {len(quads)}\n"
    body += "".join("9 " + " ".join(map(str, quad)) + " 0\n" for quad in quads)
    body += "NMARK= 0\n"
    path.write_text(body, encoding="utf-8")


def _write_cgns(path, quads) -> None:
    """Hand-built SIDS/HDF5 zone with quad blocks — deliberately not ``write_cgns``.

    Using our own writer would test our writer, not the *foreign* grammar: the
    point of a ``QUAD_4`` block is that no TensorCAD writer emits one.
    """
    h5py = pytest.importorskip("h5py")
    with h5py.File(path, "w") as f:
        mf._sids_attrs(f, {"version": "3.4", "format": "0.9"})
        mf._sids_node(
            f, "CGNSLibraryVersion", "CGNSLibraryVersion_t", "R4", 1, np.array([4.3], dtype=np.float32)
        )
        base = mf._sids_node(f, "Base", "CGNSBase_t", "I4", 1, np.array([3, 3], dtype=np.int32))
        zone = mf._sids_node(
            base, "Zone", "Zone_t", "I4", 1,
            np.array([len(QUAD_CAGE_POINTS), len(quads), 0], dtype=np.int32),
        )
        zt = mf._sids_node(zone, "ZoneType", "ZoneType_t", "C1", 1)
        zt.create_dataset(" data", data=np.bytes_("Unstructured")).attrs[" type"] = np.bytes_("C1")
        gc = mf._sids_node(zone, "GridCoordinates", "GridCoordinates_t", "MT", 1)
        for k, axis in enumerate("XYZ"):
            mf._sids_node(
                gc, f"Coordinate{axis}", "DataArray_t", "R8", k + 1,
                np.ascontiguousarray(QUAD_CAGE_POINTS[:, k]),
            )
        el = mf._sids_node(zone, "Elements", "Elements_t", "I4", 2, np.array([mf.CGNS_QUAD_4, 0], dtype=np.int32))
        mf._sids_node(
            el, "ElementRange", "IndexRange_t", "I4", 1,
            np.array([1, len(quads)], dtype=np.int32),
        )
        mf._sids_node(
            el, "ElementConnectivity", "DataArray_t", "I4", 2,
            (np.array(quads) + 1).astype(np.int32).ravel(),
        )


_WRITERS = {"obj": _write_obj, "vtk": _write_vtk, "gmsh": _write_gmsh, "su2": _write_su2, "cgns": _write_cgns}


def _ext(fmt: str) -> str:
    return "msh" if fmt == "gmsh" else fmt


@pytest.fixture
def quad_cage_file(tmp_path, fmt):
    path = tmp_path / f"quad_cage.{_ext(fmt)}"
    _WRITERS[fmt](path, QUAD_CAGE_QUADS)
    return path


# ---------------------------------------------------------------------------
# the intake matrix: a quad cage through every native reader
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("fmt", QUAD_FORMATS)
def test_quad_cage_intake_is_geometrically_exact(fmt, quad_cage_file):
    """6 quads in, 12 triangles out, volume exactly 1.0 — same as the tri cage."""
    got = mf.read(quad_cage_file, fmt)
    assert len(got.faces) == 2 * len(QUAD_CAGE_QUADS), "every quad must yield two triangles"
    report = mf.compare_meshes(_REF, got, tol=0.0)
    assert report["pass"], report
    assert report["volume_rel_err"] == 0.0
    assert got.euler() == 2 and got.boundary_edges() == 0, "closed genus-0 body must stay closed"
    assert abs(got.volume() - 1.0) < 1e-12


@pytest.mark.parametrize("fmt", QUAD_FORMATS)
def test_quad_intake_does_not_launder_orientation(fmt, tmp_path):
    """An inward-wound cage must read back with *negative* volume, not be fixed up."""
    flipped = tmp_path / f"flipped.{_ext(fmt)}"
    _WRITERS[fmt](flipped, [list(reversed(quad)) for quad in QUAD_CAGE_QUADS])
    got = mf.read(flipped, fmt)
    assert len(got.faces) == 2 * len(QUAD_CAGE_QUADS)
    assert got.volume() == pytest.approx(-1.0, abs=1e-12), "orientation flip must survive intake"


@pytest.mark.parametrize("fmt", ("obj", "vtk", "gmsh", "su2"))
def test_quad_triangulation_agrees_with_the_meshio_path(fmt, quad_cage_file):
    """Two independently written importers must triangulate the same cage alike."""
    meshio = pytest.importorskip("meshio")
    from tensorcad.meshio_exchange import meshio_to_mesh

    native = mf.read(quad_cage_file, fmt)
    foreign = meshio_to_mesh(meshio.read(str(quad_cage_file), file_format=fmt))
    via_meshio = mf.TriMesh(
        np.array(foreign.vertices, dtype=float),
        np.array(foreign.faces, dtype=np.int64),
        "meshio_path",
    )
    assert len(via_meshio.faces) == len(native.faces)
    assert via_meshio.volume() == pytest.approx(native.volume(), abs=1e-12)
    assert via_meshio.euler() == native.euler() == 2
    forward = {tuple(sorted(f)) for f in native.faces.tolist()}
    reverse = {tuple(sorted(f)) for f in via_meshio.faces.tolist()}
    assert forward == reverse, "the two paths must triangulate the same quads"


# ---------------------------------------------------------------------------
# planar SU2 (NDIME=2) — legal SU2 that used to raise "points must be (n,3)"
# ---------------------------------------------------------------------------
def test_su2_ndime2_planar_surface_is_read_with_zero_fill(tmp_path):
    path = tmp_path / "planar.su2"
    path.write_text(
        "NDIME= 2\nNPOIN= 4\n0.0 0.0\n1.0 0.0\n1.0 1.0\n0.0 1.0\n"
        "NELEM= 2\n5 0 1 2 0\n5 0 2 3 0\nNMARK= 0\n",
        encoding="utf-8",
    )
    got = mf.read_su2(path)
    assert len(got.faces) == 2
    assert got.area() == pytest.approx(1.0, abs=1e-15)
    assert np.all(got.points[:, 2] == 0.0), "a 2-D file must land in the z=0 plane"


def test_su2_ndime2_quad_marker_block_only(tmp_path):
    """A marker-only geometry file (NELEM=0) with quad markers is a surface patch."""
    path = tmp_path / "marker_quad.su2"
    path.write_text(
        "NDIME= 2\nNPOIN= 4\n0.0 0.0\n1.0 0.0\n1.0 1.0\n0.0 1.0\n"
        "NELEM= 0\nNMARK= 1\nMARKER_TAG= wall\nMARKER_ELEMS= 1\n9 0 1 2 3\n",
        encoding="utf-8",
    )
    got = mf.read_su2(path)
    assert len(got.faces) == 2
    assert got.area() == pytest.approx(1.0, abs=1e-15)


def test_su2_ndime_outside_2_3_is_refused(tmp_path):
    path = tmp_path / "bad.su2"
    path.write_text("NDIME= 9\nNPOIN= 3\n0 0 0\n1 0 0\n0 1 0\nNELEM= 1\n5 0 1 2 0\n", encoding="utf-8")
    with pytest.raises(ValueError, match="NDIME=9"):
        mf.read_su2(path)


# ---------------------------------------------------------------------------
# multi-block CGNS: a hull patch exported per-patch arrives as several blocks
# ---------------------------------------------------------------------------
def test_cgns_multiblock_tri_and_quad_blocks_are_both_consumed(tmp_path):
    h5py = pytest.importorskip("h5py")
    path = tmp_path / "multiblock.cgns"
    tri_block = [[0, 1, 2], [0, 2, 3]]
    quad_block = [[4, 5, 6, 7], [0, 1, 5, 4]]
    with h5py.File(path, "w") as f:
        mf._sids_attrs(f, {"version": "3.4", "format": "0.9"})
        mf._sids_node(f, "CGNSLibraryVersion", "CGNSLibraryVersion_t", "R4", 1, np.array([4.3], np.float32))
        base = mf._sids_node(f, "Base", "CGNSBase_t", "I4", 1, np.array([3, 3], np.int32))
        zone = mf._sids_node(base, "Zone", "Zone_t", "I4", 1, np.array([8, 4, 0], np.int32))
        zt = mf._sids_node(zone, "ZoneType", "ZoneType_t", "C1", 1)
        zt.create_dataset(" data", data=np.bytes_("Unstructured")).attrs[" type"] = np.bytes_("C1")
        gc = mf._sids_node(zone, "GridCoordinates", "GridCoordinates_t", "MT", 1)
        for k, axis in enumerate("XYZ"):
            mf._sids_node(
                gc, f"Coordinate{axis}", "DataArray_t", "R8", k + 1,
                np.ascontiguousarray(QUAD_CAGE_POINTS[:, k]),
            )
        for i, (etype, conn) in enumerate(((mf.CGNS_TRI_3, tri_block), (mf.CGNS_QUAD_4, quad_block))):
            el = mf._sids_node(zone, f"Elements{i + 1}", "Elements_t", "I4", 2, np.array([etype, 0], np.int32))
            first = 1 + sum(len(c) for _e, c in ((mf.CGNS_TRI_3, tri_block), (mf.CGNS_QUAD_4, quad_block))[:i])
            mf._sids_node(
                el, "ElementRange", "IndexRange_t", "I4", 1,
                np.array([first, first + len(conn) - 1], np.int32),
            )
            mf._sids_node(
                el, "ElementConnectivity", "DataArray_t", "I4", 2,
                (np.array(conn) + 1).astype(np.int32).ravel(),
            )
    assert mf.sids_violations(path) == [], "the extended checker must accept well-formed multi-block files"
    notes: list[str] = []
    got = mf.read_cgns(path, notes=notes)
    assert len(got.faces) == len(tri_block) + 2 * len(quad_block), "both blocks, quads triangulated"
    assert len(notes) == 2 and all("type-" in n for n in notes)


def test_cgns_disjoint_element_ranges_are_refused(tmp_path):
    """Two blocks claiming [1,2] and [5,6] leave cells 3-4 undefined: not conformant."""
    h5py = pytest.importorskip("h5py")
    path = tmp_path / "gap.cgns"
    points = QUAD_CAGE_POINTS
    with h5py.File(path, "w") as f:
        mf._sids_attrs(f, {"version": "3.4", "format": "0.9"})
        mf._sids_node(f, "CGNSLibraryVersion", "CGNSLibraryVersion_t", "R4", 1, np.array([4.3], np.float32))
        base = mf._sids_node(f, "Base", "CGNSBase_t", "I4", 1, np.array([3, 3], np.int32))
        zone = mf._sids_node(base, "Zone", "Zone_t", "I4", 1, np.array([len(points), 4, 0], np.int32))
        zt = mf._sids_node(zone, "ZoneType", "ZoneType_t", "C1", 1)
        zt.create_dataset(" data", data=np.bytes_("Unstructured")).attrs[" type"] = np.bytes_("C1")
        gc = mf._sids_node(zone, "GridCoordinates", "GridCoordinates_t", "MT", 1)
        for k, axis in enumerate("XYZ"):
            mf._sids_node(
                gc, f"Coordinate{axis}", "DataArray_t", "R8", k + 1,
                np.ascontiguousarray(points[:, k]),
            )
        for i, (first, conn) in enumerate(((1, [[0, 1, 2]]), (5, [[0, 1, 2]]))):
            el = mf._sids_node(zone, f"Elements{i + 1}", "Elements_t", "I4", 2, np.array([mf.CGNS_TRI_3, 0], np.int32))
            mf._sids_node(el, "ElementRange", "IndexRange_t", "I4", 1, np.array([first, first], np.int32))
            mf._sids_node(
                el, "ElementConnectivity", "DataArray_t", "I4", 2,
                (np.array(conn) + 1).astype(np.int32).ravel(),
            )
    violations = mf.sids_violations(path)
    assert any("does not continue" in v for v in violations), violations
    with pytest.raises(ValueError, match="non-conformant"):
        mf.read_cgns(path)


# ---------------------------------------------------------------------------
# no silent loss, and loud refusal of non-surface content
# ---------------------------------------------------------------------------
def test_su2_volume_cells_are_refused_by_the_surface_reader(tmp_path):
    path = tmp_path / "tet.su2"
    path.write_text(
        "NDIME= 3\nNPOIN= 4\n0 0 0\n1 0 0\n0 1 0\n0 0 1\nNELEM= 1\n10 0 1 2 3 0\nNMARK= 0\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="volume cell"):
        mf.read_su2(path)


@pytest.mark.parametrize("fmt", ("su2", "gmsh"))
def test_non_surface_elements_must_be_acknowledged_not_dropped(fmt, tmp_path):
    """A line element beside a triangle: raise by default, counted when notes=[] is passed."""
    if fmt == "su2":
        path = tmp_path / "mixed.su2"
        path.write_text(
            "NDIME= 3\nNPOIN= 3\n0 0 0\n1 0 0\n0 1 0\n"
            "NELEM= 1\n5 0 1 2 0\nNMARK= 1\nMARKER_TAG= wall\nMARKER_ELEMS= 1\n3 0 1\n",
            encoding="utf-8",
        )
        reader = mf.read_su2
    else:
        path = tmp_path / "mixed.msh"
        path.write_text(
            "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n3\n"
            "1 0 0 0\n2 1 0 0\n3 0 1 0\n$EndNodes\n$Elements\n2\n"
            "1 2 0 1 2 3\n2 1 0 1 2\n$EndElements\n",
            encoding="utf-8",
        )
        reader = mf.read_gmsh
    with pytest.raises(ValueError, match="pass notes="):
        reader(path)
    notes: list[str] = []
    got = reader(path, notes=notes)
    assert len(got.faces) == 1 and len(notes) == 1
    assert ("1" in notes[0]) and ("line" in notes[0])


def test_gmsh_without_surface_elements_says_so(tmp_path):
    path = tmp_path / "points_only.msh"
    path.write_text(
        "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n2\n1 0 0 0\n2 1 0 0\n$EndNodes\n"
        "$Elements\n0\n$EndElements\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="no triangles or quads"):
        mf.read_gmsh(path)


def test_su2_without_surface_cells_says_so(tmp_path):
    path = tmp_path / "empty.su2"
    path.write_text("NDIME= 3\nNPOIN= 3\n0 0 0\n1 0 0\n0 1 0\nNELEM= 0\nNMARK= 0\n", encoding="utf-8")
    with pytest.raises(ValueError, match="no surface triangles or quads"):
        mf.read_su2(path)

# ---------------------------------------------------------------------------
# Polygon (n>4) intake — CAD kernels often produce pentagons/hexagons
# ---------------------------------------------------------------------------

#: A single pentagonal face on the XY plane (5 vertices, 1 pentagon → 3 triangles).
PENTAGON_POINTS = np.array([
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.2, 0.8, 0.0],
    [0.5, 1.3, 0.0],
    [-0.2, 0.6, 0.0],
])
PENTAGON = [[0, 1, 2, 3, 4]]
# Fan triangulation from vertex 0: (0,1,2), (0,2,3), (0,3,4) → 3 triangles (n-2)
_PENTAGON_TRIS = np.array([[0, 1, 2], [0, 2, 3], [0, 3, 4]], dtype=np.int64)
_PENTAGON_REF = mf.TriMesh(PENTAGON_POINTS, _PENTAGON_TRIS, "pentagon_reference")


def _write_pentagon_obj(path) -> None:
    body = "".join(f"v {x} {y} {z}\n" for x, y, z in PENTAGON_POINTS)
    body += "f 1 2 3 4 5\n"
    path.write_text(body, encoding="utf-8")


def _write_pentagon_vtk(path) -> None:
    head = "# vtk DataFile Version 3.0\npentagon\nASCII\nDATASET UNSTRUCTURED_GRID\n"
    head += f"POINTS {len(PENTAGON_POINTS)} float\n"
    head += "".join(f"{x} {y} {z}\n" for x, y, z in PENTAGON_POINTS)
    head += f"CELLS 1 6\n5 0 1 2 3 4\n"
    head += "CELL_TYPES 1\n7\n"  # VTK_POLYGON = 7
    path.write_text(head, encoding="utf-8")


_PENTAGON_WRITERS = {"obj": _write_pentagon_obj, "vtk": _write_pentagon_vtk}


@pytest.mark.parametrize("fmt", ("obj", "vtk"))
def test_pentagon_intake_fan_triangulates_correctly(fmt, tmp_path):
    """A 5-node polygon must yield 3 triangles (n-2) via fan triangulation from vertex 0.

    A single planar face is not a closed surface, so we check triangulation
    directly rather than using compare_meshes (which requires closedness).
    """
    path = tmp_path / f"pentagon.{_ext(fmt)}"
    _PENTAGON_WRITERS[fmt](path)
    got = mf.read(path, fmt)
    assert len(got.faces) == 3, f"pentagon must yield 3 triangles, got {len(got.faces)}"
    # Check the triangulation matches the reference (same face set)
    got_faces = {tuple(sorted(f)) for f in got.faces.tolist()}
    ref_faces = {tuple(sorted(f)) for f in _PENTAGON_REF.faces.tolist()}
    assert got_faces == ref_faces, f"triangulation mismatch: {got_faces} vs {ref_faces}"
    assert got.area() == pytest.approx(_PENTAGON_REF.area(), abs=1e-15)


def test_hexagon_intake_via_obj(tmp_path):
    """A 6-node polygon must yield 4 triangles (n-2) via fan triangulation."""
    hex_pts = np.array([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.5, 0.8, 0.0],
        [1.0, 1.6, 0.0], [0.0, 1.6, 0.0], [-0.5, 0.8, 0.0],
    ])
    hex_faces = np.array([
        [0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 5]
    ], dtype=np.int64)
    hex_ref = mf.TriMesh(hex_pts, hex_faces, "hexagon_reference")

    path = tmp_path / "hexagon.obj"
    body = "".join(f"v {x} {y} {z}\n" for x, y, z in hex_pts)
    body += "f 1 2 3 4 5 6\n"
    path.write_text(body, encoding="utf-8")

    got = mf.read(path, "obj")
    assert len(got.faces) == 4, f"hexagon must yield 4 triangles, got {len(got.faces)}"
    got_faces = {tuple(sorted(f)) for f in got.faces.tolist()}
    ref_faces = {tuple(sorted(f)) for f in hex_ref.faces.tolist()}
    assert got_faces == ref_faces, f"triangulation mismatch: {got_faces} vs {ref_faces}"
    assert got.area() == pytest.approx(hex_ref.area(), abs=1e-15)
