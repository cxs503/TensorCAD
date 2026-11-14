"""Independent audit of the 6-format exchange matrix (anchor GAP-CAD-B2).

The matrix in ``tests/test_meshio_format_matrix.py`` writes with the package and
reads with the package.  These tests verify the *files* from the other side:

* ``out``: ``write_meshio`` -> parse the bytes with
  :mod:`tensorcad.format_conformance` (spec readers, no meshio, no TensorCAD
  reader) -> volume / topology / orientation / closedness must match the source;
* ``in``: meshio writes from a plain ``meshio.Mesh`` -> ``read_meshio`` must
  recover the same geometry.

CGNS gets an extra gate: the file we export must satisfy the SIDS node/attribute
table (``cgns_issues`` returns nothing), while a meshio-written CGNS file must
fail it -- that defect is what this anchor was missing.
"""
from __future__ import annotations

import numpy as np
import pytest

from test_meshio_format_matrix import ANALYTIC_VOLUME, GEOMETRIES

from tensorcad import format_conformance as fc
from tensorcad.mesh_cgns import read_cgns, write_cgns
from tensorcad.meshio_exchange import (
    FORMATS_ALL,
    FORMATS_VOLUME,
    read_meshio,
    surface_signed_volume,
    write_meshio,
)

GEOM_NAMES = ("cube", "sphere", "wigley_hull")
VOLUME_TOL = {"stl": 1e-5}
DEFAULT_TOL = 1e-12


def _geom(name):
    return next(g for g in GEOMETRIES if g.name == name)


def _volume(mesh):
    return mesh.volume() if hasattr(mesh, "tets") else surface_signed_volume(mesh)


# ---------------------------------------------------------------------------
# "out": our writer, verified by an independent parser of the file bytes
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("fmt", FORMATS_ALL)
@pytest.mark.parametrize("geom_name", GEOM_NAMES)
def test_export_survives_spec_parser(fmt, geom_name, tmp_path):
    geom = _geom(geom_name)
    source = geom.tetra if fmt in FORMATS_VOLUME else geom.mesh
    path = write_meshio(tmp_path / f"{fmt}_{geom_name}", source, file_format=fmt)

    parsed, meta = fc.parse_file(path, fmt)

    source_volume = _volume(source)
    assert parsed.volume() == pytest.approx(source_volume, rel=VOLUME_TOL.get(fmt, DEFAULT_TOL))
    assert parsed.orientation() == "outward"
    assert len(parsed.points) == len(source.vertices)
    expected_cells = len(source.tets) if fmt in FORMATS_VOLUME else len(source.faces)
    assert len(parsed.cells) == expected_cells
    if fmt not in FORMATS_VOLUME:
        assert parsed.open_edges() == 0, "exported surface is not closed"
    # No format-spec violation in any encoding we emit.
    assert meta.get("issues") == []


def test_stl_spec_parser_checks_binary_file_size(tmp_path):
    """A binary STL whose declared triangle count overruns the file is invalid.

    (meshio writes ASCII STL, so the binary layout is exercised with a
    hand-built file: header count 2 but only one 50-byte facet present.)
    """
    path = tmp_path / "truncated.stl"
    path.write_bytes(b"\0" * 80 + (2).to_bytes(4, "little") + b"\0" * 50)
    with pytest.raises(fc.FormatParseError, match="84\+50n"):
        fc.parse_stl(path)


def test_stl_spec_parser_reports_what_the_ascii_file_declares(tmp_path):
    """Dropping a facet from an ASCII STL must show up as one facet fewer."""
    path = write_meshio(tmp_path / "cube.stl", _geom("cube").mesh, file_format="stl")
    data = path.read_bytes()
    assert data.lstrip().startswith(b"solid"), "meshio writes ASCII STL"
    path.write_bytes(data[:-50])
    parsed, meta = fc.parse_stl(path)
    assert meta["declared_triangles"] == len(_geom("cube").mesh.faces) - 1
    assert len(parsed.cells) == meta["declared_triangles"]


def test_obj_spec_parser_rejects_zero_index(tmp_path):
    """OBJ indices are 1-based; index 0 is a spec violation, not '-1'."""
    path = tmp_path / "bad.obj"
    path.write_text("v 0 0 0\nv 1 0 0\nv 0 1 0\nf 0 1 2\n")
    with pytest.raises(fc.FormatParseError, match="1-based"):
        fc.parse_obj(path)


def test_su2_spec_parser_rejects_wrong_node_count(tmp_path):
    """SU2 element type 10 (tetra) must carry exactly 4 nodes."""
    path = tmp_path / "bad.su2"
    path.write_text(
        "NDIME= 3\nNPOIN= 3\n0 0 0\n1 0 0\n0 1 0\nNELEM= 1\n10 0 1 2\nNMARK= 0\n"
    )
    with pytest.raises(fc.FormatParseError, match="expected 4"):
        fc.parse_su2(path)


def test_gmsh_spec_parser_only_claims_msh41(tmp_path):
    path = tmp_path / "old.msh"
    path.write_text("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
    with pytest.raises(fc.FormatParseError, match="4.1"):
        fc.parse_gmsh(path)


# ---------------------------------------------------------------------------
# "in": meshio writes, our reader must recover the geometry
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("fmt", FORMATS_ALL)
@pytest.mark.parametrize("geom_name", GEOM_NAMES)
def test_import_from_foreign_writer(fmt, geom_name, tmp_path):
    meshio = pytest.importorskip("meshio")
    geom = _geom(geom_name)
    if fmt in FORMATS_VOLUME:
        points = np.array(geom.tetra.vertices, dtype=np.float64)
        cells = [("tetra", np.array(geom.tetra.tets, dtype=np.int64))]
        source_volume = geom.tetra.volume()
    else:
        points = np.array(geom.mesh.vertices, dtype=np.float64)
        cells = [("triangle", np.array(geom.mesh.faces, dtype=np.int64))]
        source_volume = surface_signed_volume(geom.mesh)
    path = tmp_path / f"foreign_{fmt}"
    meshio.write(str(path), meshio.Mesh(points, cells), file_format=fmt)

    back = read_meshio(path, file_format=fmt)

    assert _volume(back) == pytest.approx(source_volume, rel=VOLUME_TOL.get(fmt, DEFAULT_TOL))
    assert len(back.vertices) == len(points)
    assert (len(back.tets) if hasattr(back, "tets") else len(back.faces)) == len(cells[0][1])


# ---------------------------------------------------------------------------
# CGNS/SIDS: the defect this anchor was carrying
# ---------------------------------------------------------------------------


def test_meshio_cgns_writer_is_not_sids_conforming(tmp_path):
    """Pin the finding: meshio's CGNS output violates the SIDS node table."""
    meshio = pytest.importorskip("meshio")
    geom = _geom("cube")
    path = tmp_path / "meshio.cgns"
    meshio.write(
        str(path),
        meshio.Mesh(
            np.array(geom.tetra.vertices, dtype=np.float64),
            [("tetra", np.array(geom.tetra.tets, dtype=np.int64))],
        ),
        file_format="cgns",
    )
    problems = fc.cgns_issues(path)
    assert problems, "meshio's CGNS writer is expected to be non-conforming"
    joined = " ".join(problems)
    assert "CGNSLibraryVersion" in joined
    assert "ElementType" in joined
    assert any("label" in p for p in problems)


@pytest.mark.parametrize("geom_name", GEOM_NAMES)
def test_cgns_sids_export_is_conforming_and_readable_both_ways(geom_name, tmp_path):
    geom = _geom(geom_name)
    path = write_cgns(tmp_path / f"{geom_name}.cgns", geom.tetra)

    assert fc.cgns_issues(path) == []

    our_back = read_cgns(path)
    assert _volume(our_back) == pytest.approx(geom.tetra.volume(), rel=1e-12)
    assert len(our_back.tets) == len(geom.tetra.tets)

    meshio = pytest.importorskip("meshio")
    foreign = meshio.read(str(path), file_format="cgns")
    assert foreign.cells[0].type == "tetra"
    assert len(foreign.cells[0].data) == len(geom.tetra.tets)


def test_cgns_surface_export_round_trips_through_spec_reader(tmp_path):
    """A TRI_3 CGNS file is legal SIDS but meshio cannot read it -- we can."""
    geom = _geom("cube")
    path = write_cgns(tmp_path / "surface.cgns", geom.mesh)
    assert fc.cgns_issues(path) == []
    back = read_cgns(path)
    assert surface_signed_volume(back) == pytest.approx(1.0, rel=1e-12)
    assert len(back.faces) == len(geom.mesh.faces)


def test_meshio_round_trip_uses_the_sids_writer(tmp_path):
    """``round_trip_meshio`` must go through the conforming CGNS exporter."""
    from tensorcad.meshio_exchange import round_trip_meshio

    geom = _geom("cube")
    back = round_trip_meshio(geom.tetra, file_format="cgns")
    assert back.volume() == pytest.approx(geom.tetra.volume(), rel=1e-12)
    assert len(back.tets) == len(geom.tetra.tets)


def test_hull_volume_from_file_matches_analytic_within_discretisation():
    """Context for the hull cell: the facet volume is below the analytic body by
    the polygonal discretisation, not by a format loss (measured 4.3e-3)."""
    geom = _geom("wigley_hull")
    source = surface_signed_volume(geom.mesh)
    analytic = ANALYTIC_VOLUME["wigley_hull"]
    assert abs(source - analytic) / analytic < 1e-2