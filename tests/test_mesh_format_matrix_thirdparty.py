"""Third-party verification of the 6-format x 3-geometry mesh exchange matrix.

``tests/test_meshio_format_matrix.py`` proves round-trip fidelity by writing
with meshio and reading back with meshio.  That is a *self-consistency* check:
a private convention shared by the writer and the reader (a wrong cell-type
code, a byte order, an offset convention) passes it while another tool would
read the file differently or not at all.

This module re-derives the same 18 combinations from the **file bytes** with
the independent parsers in :mod:`tensorcad.mesh_interop_independent` -- no
meshio, no ``tensorcad.mesh_cgns`` -- and asserts the ledger acceptance:

    往返保真：体积/拓扑/取向不变（6 格式 x 3 几何 = 18 组合全绿表）

with a stricter orientation invariant than the in-package matrix.  That matrix
checks ``sign(volume)``, which two opposite facet flips can keep while the
surface is broken; here every cell carries its own winding sign and the two
multisets must match exactly.

The negative controls at the bottom are the reason the table means anything:
they mutate a file that the table calls green and require the harness to catch
it, so a PASS here cannot come from a check that merely executes.
"""
from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pytest

from tensorcad.mesh_interop_independent import (
    ANALYTIC_VOLUME_3,
    IndependentReadError,
    oriented_signature,
    read_independent,
    read_obj,
    read_stl_ascii,
    run_format_matrix,
    surface_signed_volume,
)

GEOMETRIES = ("cube", "capsule", "wigley_hull")
FORMATS = ("stl", "obj", "vtk", "gmsh", "su2", "cgns")
VOLUME_FORMATS = ("su2", "cgns")

#: The in-package matrix accepts 1e-3 for STL; that tolerance dates from when
#: STL was assumed to store float32 positions.  ASCII STL is exact here, so the
#: independent check holds every format to float64 round-off.
VOLUME_TOL = 1e-12


@pytest.fixture(scope="module")
def matrix(tmp_path_factory) -> list[dict]:
    workdir = tmp_path_factory.mktemp("thirdparty-matrix")
    rows = run_format_matrix(workdir=workdir)
    assert len(rows) == len(GEOMETRIES) * len(FORMATS) == 18
    return rows


@pytest.fixture(scope="module")
def matrix_by_key(matrix: list[dict]) -> dict[tuple[str, str], dict]:
    return {(r["geometry"], r["format"]): r for r in matrix}


def _row(matrix_by_key, geometry: str, file_format: str) -> dict:
    return matrix_by_key[(geometry, file_format)]


@pytest.mark.parametrize("file_format", FORMATS)
@pytest.mark.parametrize("geometry", GEOMETRIES)
def test_matrix_is_green(matrix_by_key, geometry: str, file_format: str) -> None:
    """Every one of the 18 combinations survives a third-party read."""
    row = _row(matrix_by_key, geometry, file_format)
    assert row["ok"], f"{geometry}/{file_format}: {row}"
    assert row["cell_count_preserved"], row
    assert row["winding_preserved"], (
        f"{geometry}/{file_format}: per-cell winding changed when the file was "
        f"re-read by an independent parser"
    )
    assert row["volume_relative_error"] <= VOLUME_TOL, row


@pytest.mark.parametrize("file_format", FORMATS)
@pytest.mark.parametrize("geometry", GEOMETRIES)
def test_volume_matches_across_all_formats(matrix_by_key, geometry: str, file_format: str) -> None:
    """The six formats must agree with each other to the last bit.

    A format that silently dropped a facet would land on its own volume; this
    cross-format identity is what makes the table an equivalence proof rather
    than six independent round-trips.
    """
    reference = _row(matrix_by_key, geometry, "obj")["volume_independent"]
    mine = _row(matrix_by_key, geometry, file_format)["volume_independent"]
    assert abs(mine - reference) / max(abs(reference), 1e-300) <= 1e-12


@pytest.mark.parametrize("geometry", GEOMETRIES)
def test_discrete_volume_agrees_with_analytic(matrix_by_key, geometry: str) -> None:
    """The exchanged volume also matches the *analytic* volume of the fixture.

    Discretisation error is expected for the curved fixtures (faceted surface),
    so the bound is 1 % for the capsule and the Wigley hull and exact for the
    cube; the point is that no format adds to it.
    """
    analytic = ANALYTIC_VOLUME_3[geometry]
    for fmt in FORMATS:
        v = _row(matrix_by_key, geometry, fmt)["volume_independent"]
        rel = abs(v - analytic) / analytic
        limit = 0.0 if geometry == "cube" else 0.01
        assert rel <= limit, f"{geometry}/{fmt}: {rel:.3e} > {limit:g} vs analytic {analytic}"


@pytest.mark.parametrize("geometry", GEOMETRIES)
def test_format_identity_records(matrix_by_key, geometry: str) -> None:
    """The bytes carry the records their format promises (not just a file name)."""
    stl = _row(matrix_by_key, geometry, "stl")
    assert stl["identity"]["starts_solid"] and stl["identity"]["ends_endsolid"]
    # STL stores a redundant normal per facet; it must agree with the winding.
    assert stl["stl_normals"]["min_cos"] >= 1.0 - 1e-9

    obj = _row(matrix_by_key, geometry, "obj")
    assert obj["identity"]["has_v"] and obj["identity"]["has_f"]

    vtk = _row(matrix_by_key, geometry, "vtk")
    assert vtk["identity"]["is_legacy_vtk"] and vtk["identity"]["has_unstructured_grid"]
    assert vtk["declared_meta"]["encoding"] == "BINARY"
    assert vtk["declared_meta"]["cell_types"] == [5]  # VTK_TRIANGLE

    gmsh = _row(matrix_by_key, geometry, "gmsh")
    assert gmsh["identity"]["version_line"].startswith("4.1")
    assert gmsh["identity"]["has_end_elements"]
    assert gmsh["declared_meta"]["element_types"] == [2]  # gmsh triangle
    tags = gmsh["declared_meta"]["element_tag_range"]
    assert tags == [1, gmsh["cells"]], "MSH 4.1 element tags must run 1..N"

    su2 = _row(matrix_by_key, geometry, "su2")
    assert su2["identity"]["has_ndime"] and su2["identity"]["has_nelem"]
    assert su2["declared_meta"]["ndime"] == 3
    assert su2["declared_meta"]["element_types"] == [10]  # SU2 tetra

    cgns = _row(matrix_by_key, geometry, "cgns")
    assert cgns["identity"]["hdf5_magic"] and cgns["identity"]["has_cgns_library_version"]
    assert cgns["declared_meta"]["cell_dimension"] == 3
    assert cgns["declared_meta"]["physical_dimension"] == 3
    assert cgns["declared_meta"]["cgns_library_version"] >= 3.0


@pytest.mark.parametrize("geometry", GEOMETRIES)
def test_volume_formats_are_tetrahedra(matrix_by_key, geometry: str) -> None:
    """SU2/CGNS rows must be tetra meshes, not the boundary surface renamed."""
    for fmt in VOLUME_FORMATS:
        row = _row(matrix_by_key, geometry, fmt)
        assert row["cell_kind"] == "tetra"
        assert row["cells"] == row["cells_reference"]


def test_stl_vertex_count_is_a_soup(matrix_by_key) -> None:
    """Documented limitation: STL duplicates every vertex per facet.

    The acceptance claim is "volume/topology/orientation preserved", and for STL
    topology is preserved at the *facet* level only -- vertex identity is lost by
    the format, so 3 x facets vertices come back.  Stating it keeps the green
    table honest.
    """
    row = _row(matrix_by_key, "capsule", "stl")
    assert row["vertex_count_preserved"] is False
    assert row["vertices"] == 3 * row["cells"] == 3 * row["cells_reference"]


# ---------------------------------------------------------------------------
# Negative controls: the harness must catch a corrupted file
# ---------------------------------------------------------------------------

def _write_cube_obj(tmp_path: Path) -> Path:
    from tensorcad.mesh_interop_independent import campaign_mesh

    mesh = campaign_mesh("cube")
    path = tmp_path / "cube.obj"
    edges = set()
    lines = ["# hand-written obj (no meshio)"]
    for v in mesh.vertices:
        lines.append(f"v {v[0]!r} {v[1]!r} {v[2]!r}")
    for f in mesh.faces:
        lines.append("f " + " ".join(str(i + 1) for i in f))
        for a, b in ((f[0], f[1]), (f[1], f[2]), (f[2], f[0])):
            edges.add((min(a, b), max(a, b)))
    assert len(edges) == 18, "cube should be a closed 12-triangle surface"
    path.write_text("\n".join(lines) + "\n")
    return path


def test_harness_detects_a_flipped_facet(tmp_path: Path) -> None:
    """Re-winding one facet keeps sign(volume) ~ but must not pass the check."""
    from tensorcad.mesh_interop_independent import campaign_mesh

    good = read_obj(_write_cube_obj(tmp_path))
    reference = oriented_signature(
        np.asarray(campaign_mesh("cube").vertices), campaign_mesh("cube").faces, cell_kind="tri"
    )
    assert oriented_signature(good.vertices, good.cells, cell_kind="tri") == reference

    text = _write_cube_obj(tmp_path).read_text().splitlines()
    for i, line in enumerate(text):
        if line.startswith("f "):
            a, b, c = line.split()[1:]
            text[i] = f"f {a} {c} {b}"  # flip exactly one facet
            break
    flipped_path = tmp_path / "flipped.obj"
    flipped_path.write_text("\n".join(text) + "\n")
    flipped = read_obj(flipped_path)
    assert oriented_signature(flipped.vertices, flipped.cells, cell_kind="tri") != reference

    # The flipped facet is the bottom face, whose plane passes through the
    # origin: it contributes exactly zero to the divergence-theorem sum in
    # either winding, so the volume does not move by even one ulp.  This is the
    # concrete reason sign(volume) -- what the in-package matrix checks -- is
    # not a sufficient orientation invariant, and why the per-cell signature
    # above exists.
    assert abs(surface_signed_volume(flipped.vertices, flipped.cells) - 1.0) < 1e-15

    # Flip an off-origin facet as well: now the volume does move, by 2/3 of the
    # facet's tetrahedron (a = (0,0,0), b, c all on the far side).
    text = _write_cube_obj(tmp_path).read_text().splitlines()
    for i, line in enumerate(text):
        if line.startswith("f 5 6 7"):  # face (4, 5, 6), plane z = 1 (misses the origin)
            a, b, c = line.split()[1:]
            text[i] = f"f {a} {c} {b}"
            break
    (tmp_path / "flipped2.obj").write_text("\n".join(text) + "\n")
    flipped2 = read_obj(tmp_path / "flipped2.obj")
    assert oriented_signature(flipped2.vertices, flipped2.cells, cell_kind="tri") != reference
    assert abs(surface_signed_volume(flipped2.vertices, flipped2.cells) - 1.0) > 0.1


def test_harness_detects_a_moved_vertex(tmp_path: Path) -> None:
    good = read_obj(_write_cube_obj(tmp_path))
    lines = (tmp_path / "cube.obj").read_text().splitlines()
    for i, line in enumerate(lines):
        if line.startswith("v "):
            lines[i] = "v 0.0 0.0 0.5"  # move the first corner by half an edge
            break
    moved_path = tmp_path / "moved.obj"
    moved_path.write_text("\n".join(lines) + "\n")
    moved = read_obj(moved_path)
    assert abs(surface_signed_volume(moved.vertices, moved.cells) - 1.0) > 1e-3
    assert oriented_signature(moved.vertices, moved.cells, cell_kind="tri") != oriented_signature(
        good.vertices, good.cells, cell_kind="tri"
    )


def test_reader_rejects_truncated_gmsh(tmp_path: Path) -> None:
    from tensorcad.meshio_exchange import write_meshio
    from tensorcad.mesh_interop_independent import campaign_mesh

    path = write_meshio(tmp_path / "cut.msh", campaign_mesh("cube"), file_format="gmsh")
    raw = path.read_bytes()
    (tmp_path / "cut.msh").write_bytes(raw[: len(raw) // 2])
    with pytest.raises((IndependentReadError, ValueError, KeyError, IndexError)):
        read_independent(tmp_path / "cut.msh", file_format="gmsh")


def test_reader_rejects_non_cgns_hdf5(tmp_path: Path) -> None:
    from tensorcad.meshio_exchange import tetrahedralise_surface, write_meshio
    from tensorcad.mesh_interop_independent import campaign_mesh

    path = write_meshio(
        tmp_path / "x.cgns", tetrahedralise_surface(campaign_mesh("cube")), file_format="cgns"
    )
    raw = bytearray(path.read_bytes())
    raw[0:8] = b"NOTHDF5!"
    (tmp_path / "x.cgns").write_bytes(bytes(raw))
    with pytest.raises(IndependentReadError, match="magic"):
        read_independent(tmp_path / "x.cgns", file_format="cgns")


def test_reader_rejects_su2_point_count_lie(tmp_path: Path) -> None:
    from tensorcad.meshio_exchange import tetrahedralise_surface, write_meshio
    from tensorcad.mesh_interop_independent import campaign_mesh

    path = write_meshio(
        tmp_path / "x.su2", tetrahedralise_surface(campaign_mesh("cube")), file_format="su2"
    )
    lines = path.read_text().splitlines()
    for i, line in enumerate(lines):
        if line.startswith("NPOIN="):
            declared = int(line.split("=")[1])
            lines[i] = f"NPOIN= {declared + 1}"  # claim one point too many
            break
    (tmp_path / "lie.su2").write_text("\n".join(lines) + "\n")
    with pytest.raises((IndependentReadError, ValueError)):
        read_independent(tmp_path / "lie.su2", file_format="su2")


def test_stl_parser_rejects_binary_stl(tmp_path: Path) -> None:
    """A binary STL must be refused loudly, never half-parsed as ASCII."""
    header = b"\x00\x00\x00\x00" * 20
    path = tmp_path / "binary.stl"
    path.write_bytes(b"binary stl" + header)
    with pytest.raises(IndependentReadError):
        read_stl_ascii(path)


def test_unknown_format_is_refused(tmp_path: Path) -> None:
    path = tmp_path / "x.nas"
    path.write_text("GRID 1 0 0 0\n")
    with pytest.raises(IndependentReadError, match="no independent reader"):
        read_independent(path, file_format="nastran")


def test_verifier_module_is_independent_of_meshio() -> None:
    """Importing the verifier must not pull meshio in (subprocess, so clean)."""
    import subprocess
    import sys

    code = (
        "import sys, tensorcad.mesh_interop_independent as m; "
        "print('meshio' in sys.modules); print(m.read_independent.__name__)"
    )
    out = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    ).stdout.split()
    assert out[0] == "False", "the independent readers must not depend on meshio"
