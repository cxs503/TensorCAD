"""Foreign-dialect exchange matrix + real third-party tool execution.

Anchor: ``GAP-CAD-B2-MESH-FORMAT-MATRIX`` (meshio 全格式矩阵进出).  The
18-combination matrix (``docs/mesh_format_matrix_thirdparty/``) is green, and
its verification report names its own ceiling:

    No external CAD/solver tool was executed ... "another tool can read it" is
    evidenced by spec-conformance of the bytes plus independent parsing, not by
    running that tool.

and lists the dialects it cannot ingest (binary STL, MSH 2.2, VTK XML).  This
module closes those two gaps and nothing else:

1. **foreign dialects produced by other tools are read** -- binary STL (what
   every CAD system writes), ASCII MSH 2.2 (what most solver decks ship) and
   VTK XML ``.vtu`` (ParaView) are parsed by
   :mod:`tensorcad.mesh_foreign_dialects` and must reproduce volume, cell count
   and per-cell winding;
2. **a genuinely foreign implementation is executed** -- trimesh (its own
   STL/OBJ loaders and its own mass properties) reads files written by
   TensorCAD, and trimesh-written files are read back by the independent
   readers, so the "another tool agrees" claim rests on running that tool
   rather than on our reading of a specification.

The tolerance is per dialect and *justified*: text/float64 dialects must be
exact (``1e-12``), binary STL stores float32 positions and is bounded by
``STL_BINARY_RELATIVE_TOLERANCE``.  A binary-STL row is required to show a
non-zero error on a geometry that float32 cannot represent exactly -- otherwise
the tolerance would be hiding the format's real precision.
"""
from __future__ import annotations

import base64
import json
import struct
import subprocess
import sys
import zlib
from pathlib import Path

import numpy as np
import pytest

from tensorcad.mesh_foreign_dialects import (
    STL_BINARY_RELATIVE_TOLERANCE,
    _STL_RECORD,
    read_gmsh22,
    read_stl_binary,
    read_vtk_xml,
    run_foreign_matrix,
    trimesh_consumer_rows,
    write_gmsh22,
    write_stl_binary,
    write_vtu_ascii,
)
from tensorcad.mesh_interop_independent import (
    IndependentReadError,
    campaign_mesh,
    file_identity,
    read_independent,
)

ROOT = Path(__file__).resolve().parents[1]
GEOMETRIES = ("cube", "capsule", "wigley_hull")
DIALECTS = ("stl_binary", "gmsh22", "vtu")


@pytest.fixture(scope="module")
def foreign(tmp_path_factory) -> dict:
    workdir = tmp_path_factory.mktemp("foreign_dialects")
    rows = run_foreign_matrix(GEOMETRIES, workdir=workdir / "files")
    consumer = trimesh_consumer_rows(GEOMETRIES, workdir=workdir / "consumer")
    return {
        "rows": rows,
        "consumer": consumer,
        "index": {(r["geometry"], r["dialect"], r["producer"], r["cell_kind"]): r for r in rows},
        "workdir": workdir,
    }


SURFACE_CASES = [
    (g, d, p) for g in GEOMETRIES for d in DIALECTS for p in ("tensorcad", "meshio", "trimesh")
]


# ---------------------------------------------------------------------------
# 1. the matrix itself
# ---------------------------------------------------------------------------

def test_matrix_has_expected_number_of_rows(foreign: dict) -> None:
    # 3 geometries x 6 (dialect, producer) surface rows + 3 x 2 tetra rows
    assert len(foreign["rows"]) == 24
    assert {r["cell_kind"] for r in foreign["rows"]} == {"tri", "tetra"}


def test_every_row_is_green(foreign: dict) -> None:
    failed = [
        (r["geometry"], r["dialect"], r["producer"], r["cell_kind"])
        for r in foreign["rows"]
        if not r["ok"]
    ]
    assert failed == [], f"foreign-dialect matrix has failing rows: {failed}"


@pytest.mark.parametrize("geometry,dialect,producer", SURFACE_CASES)
def test_surface_row_present_and_green(foreign: dict, geometry, dialect, producer) -> None:
    if producer == "trimesh" and dialect != "stl_binary":
        pytest.skip("trimesh only produces STL/OBJ; it is the binary STL producer here")
    if producer == "meshio" and dialect == "stl_binary":
        pytest.skip("meshio has no binary STL writer (this is why trimesh produces it)")
    row = foreign["index"][(geometry, dialect, producer, "tri")]
    assert row["ok"], row
    assert row["cell_count_preserved"]
    assert row["winding_preserved"]
    assert row["volume_relative_error"] <= row["volume_tolerance"]


@pytest.mark.parametrize("geometry", GEOMETRIES)
@pytest.mark.parametrize("dialect", ("gmsh22", "vtu"))
def test_tetra_rows_are_green(foreign: dict, geometry, dialect) -> None:
    row = foreign["index"][(geometry, dialect, "tensorcad", "tetra")]
    assert row["cell_kind"] == "tetra"
    assert row["ok"], row
    assert row["volume_relative_error"] <= 1e-12


@pytest.mark.parametrize("geometry", GEOMETRIES)
def test_producers_agree_with_each_other(foreign: dict, geometry) -> None:
    """TensorCAD's spec writers and the third-party writers must agree exactly.

    Two independent writers of the same dialect producing different volumes
    would mean at least one of them writes a file that is not what it claims.
    """
    for dialect in DIALECTS:
        volumes = [
            r["volume_independent"]
            for r in foreign["rows"]
            if r["geometry"] == geometry and r["dialect"] == dialect and r["cell_kind"] == "tri"
        ]
        assert len(volumes) >= 2, (geometry, dialect)
        # binary STL truncates to float32, so compare it at that precision only
        tolerance = STL_BINARY_RELATIVE_TOLERANCE if dialect == "stl_binary" else 1e-12
        spread = (max(volumes) - min(volumes)) / max(abs(max(volumes)), 1e-300)
        assert spread <= tolerance, f"{geometry}/{dialect}: producers disagree by {spread:.3e}"


@pytest.mark.parametrize("geometry", GEOMETRIES)
def test_float64_dialects_are_exact(foreign: dict, geometry) -> None:
    for dialect in ("gmsh22", "vtu"):
        for row in foreign["rows"]:
            if row["geometry"] == geometry and row["dialect"] == dialect:
                assert row["volume_relative_error"] <= 1e-12, row


@pytest.mark.parametrize("geometry", ("capsule", "wigley_hull"))
def test_binary_stl_error_is_float32_bounded_but_nonzero(foreign: dict, geometry) -> None:
    """The binary STL tolerance must be *used*: float32 really does cost precision.

    If this error were exactly zero the 1e-6 tolerance would be decoration and a
    genuine loss could hide under it.
    """
    row = foreign["index"][(geometry, "stl_binary", "tensorcad", "tri")]
    assert 0.0 < row["volume_relative_error"] < STL_BINARY_RELATIVE_TOLERANCE
    assert row["declared_meta"]["position_dtype"] == "float32"


def test_byte_identity_facts_are_recorded(foreign: dict) -> None:
    for row in foreign["rows"]:
        identity = row["identity"]
        assert identity, row
        if row["dialect"] == "stl_binary":
            assert identity["position_width_bits"] == 32
            assert identity["declared_facets"] == row["cells"]
            assert identity["bytes_expected"] == identity["bytes_actual"]
        elif row["dialect"] == "gmsh22":
            assert identity["has_meshformat"] and identity["version_line"].startswith("2.2")
            assert identity["has_end_elements"]
        else:
            assert identity["xml_declaration"] and identity["vtkfile_tag"]
            assert identity["type_declared"]


def test_stl_binary_normals_agree_with_winding(foreign: dict) -> None:
    """Binary STL also stores a normal: it must match the facet's own winding."""
    for row in foreign["rows"]:
        if row["dialect"] != "stl_binary":
            continue
        normals = row["stl_normals"]
        assert normals is not None
        assert normals["facet_count"] == row["cells"]
        assert normals["min_cos"] >= 1.0 - 1e-6, (row["geometry"], row["producer"], normals)


# ---------------------------------------------------------------------------
# 2. the third-party tool (trimesh) is actually executed
# ---------------------------------------------------------------------------

def test_third_party_rows_were_produced(foreign: dict) -> None:
    if not foreign["consumer"]:
        pytest.skip("trimesh not available in this environment")
    assert len(foreign["consumer"]) == 9
    assert {r["consumer"] for r in foreign["consumer"]} == {"trimesh"}


def test_third_party_reader_agrees(foreign: dict) -> None:
    if not foreign["consumer"]:
        pytest.skip("trimesh not available in this environment")
    for row in foreign["consumer"]:
        assert row["ok"], row
        assert row["thirdparty_watertight"], row
        assert row["thirdparty_winding_consistent"], row
        assert row["volume_relative_error"] <= row["volume_tolerance"], row


def test_third_party_tool_identity_recorded(foreign: dict) -> None:
    if not foreign["consumer"]:
        pytest.skip("trimesh not available in this environment")
    versions = {r["consumer_version"] for r in foreign["consumer"]}
    assert all(isinstance(v, str) and v for v in versions), versions


def test_third_party_written_file_is_read_by_us(foreign: dict) -> None:
    """trimesh-produced binary STL parsed by our reader (bytes we did not write)."""
    row = foreign["index"][("wigley_hull", "stl_binary", "trimesh", "tri")]
    assert row["ok"] and row["cell_count_preserved"]
    assert row["identity"]["header_len"] == 80


# ---------------------------------------------------------------------------
# 3. the harness must be able to fail (mutation checks)
# ---------------------------------------------------------------------------

def _cube_surface():
    return campaign_mesh("cube")


def _flip_facet_in_binary_stl(path: Path, facet: int = 2) -> Path:
    """Reverse the winding of one facet, keeping everything else byte-identical."""
    raw = bytearray(path.read_bytes())
    count = struct.unpack("<I", raw[80:84])[0]
    assert count > facet
    offset = 84 + 50 * facet
    record = np.frombuffer(bytes(raw[offset : offset + 50]), dtype=_STL_RECORD)[0]
    vertices = record["vertices"].copy()
    vertices[[1, 2]] = vertices[[2, 1]]
    raw[offset : offset + 50] = (
        record["normal"].tobytes() + vertices.tobytes() + record["attribute"].tobytes()
    )
    assert len(raw) == 84 + 50 * count
    flipped = path.with_name("flipped_" + path.name)
    flipped.write_bytes(bytes(raw))
    return flipped


def test_harness_detects_flipped_facet_in_binary_stl(tmp_path: Path) -> None:
    from tensorcad.mesh_interop_independent import oriented_signature

    written = write_stl_binary(tmp_path / "cube.stl", _cube_surface())
    reference = read_stl_binary(written)
    flipped = read_stl_binary(_flip_facet_in_binary_stl(written, facet=2))
    assert flipped.volume() != reference.volume()
    assert oriented_signature(
        reference.vertices, reference.cells, cell_kind="tri"
    ) != oriented_signature(flipped.vertices, flipped.cells, cell_kind="tri")


def test_volume_alone_is_blind_to_a_facet_through_the_origin(tmp_path: Path) -> None:
    """Why the matrix compares per-cell winding and not just sign(volume).

    Cube facet 0 lies in the plane ``z = 0``, which passes through the origin:
    reversing it changes the divergence-theorem sum by exactly zero, so a
    volume-only gate reports "preserved" for a mesh that is no longer closed.
    The winding signature still catches it.
    """
    from tensorcad.mesh_interop_independent import oriented_signature

    written = write_stl_binary(tmp_path / "cube.stl", _cube_surface())
    reference = read_stl_binary(written)
    flipped = read_stl_binary(_flip_facet_in_binary_stl(written, facet=0))
    assert flipped.volume() == reference.volume()
    assert oriented_signature(
        reference.vertices, reference.cells, cell_kind="tri"
    ) != oriented_signature(flipped.vertices, flipped.cells, cell_kind="tri")


def test_harness_detects_flipped_facet_in_vtu(tmp_path: Path) -> None:
    from tensorcad.mesh_interop_independent import oriented_signature

    written = write_vtu_ascii(tmp_path / "cube.vtu", _cube_surface())
    text = written.read_text()
    parsed = read_vtk_xml(written)
    conn = np.concatenate([np.asarray(c) for c in parsed.cells])
    first = conn[:3].copy()
    conn[:3] = first[[0, 2, 1]]  # facet 0 lies in z = 0
    mutated = written.with_name("flipped.vtu")
    start = text.index('Name="connectivity" format="ascii">') + len('Name="connectivity" format="ascii">')
    end = text.index("<", start)
    mutated.write_text(text[:start] + " ".join(str(int(i)) for i in conn) + text[end:])
    flipped = read_vtk_xml(mutated)
    # facet 0 is in the plane z = 0: the volume is *unchanged* by construction,
    # so this test also documents the blind spot of a volume-only gate.
    assert flipped.volume() == parsed.volume()
    assert oriented_signature(
        parsed.vertices, parsed.cells, cell_kind="tri"
    ) != oriented_signature(flipped.vertices, flipped.cells, cell_kind="tri")


def test_harness_detects_moved_vertex_in_gmsh22(tmp_path: Path) -> None:
    written = write_gmsh22(tmp_path / "cube.msh", _cube_surface())
    reference = read_gmsh22(written)
    text = written.read_text().replace("\n1 0 0 0\n", "\n1 0 0 0.5\n", 1)
    assert text != written.read_text()
    moved = written.with_name("moved.msh")
    moved.write_text(text)
    assert read_gmsh22(moved).volume() != reference.volume()


# ---------------------------------------------------------------------------
# 4. malformed input must raise, never silently produce a smaller mesh
# ---------------------------------------------------------------------------

def test_binary_stl_truncated_raises(tmp_path: Path) -> None:
    written = write_stl_binary(tmp_path / "cube.stl", _cube_surface())
    raw = written.read_bytes()
    truncated = written.with_name("truncated.stl")
    truncated.write_bytes(raw[: 84 + 50 * 3])  # declares 12 facets, carries 3
    with pytest.raises(IndependentReadError, match="declares 12 facets"):
        read_stl_binary(truncated)


def test_binary_reader_rejects_ascii_stl(tmp_path: Path) -> None:
    ascii_stl = tmp_path / "ascii.stl"
    ascii_stl.write_text(
        "solid x\n facet normal 0 0 1\n  outer loop\n   vertex 0 0 0\n"
        "   vertex 1 0 0\n   vertex 0 1 0\n  endloop\n endfacet\nendsolid x\n"
    )
    with pytest.raises(IndependentReadError, match="ASCII STL"):
        read_stl_binary(ascii_stl)


def test_gmsh22_rejects_msh41_header(tmp_path: Path) -> None:
    path = tmp_path / "v41.msh"
    path.write_text("$MeshFormat\n4.1 0 8\n$EndMeshFormat\n")
    with pytest.raises(IndependentReadError, match="not 2.x"):
        read_gmsh22(path)


def test_gmsh22_rejects_binary_flag(tmp_path: Path) -> None:
    path = tmp_path / "bin.msh"
    path.write_text("$MeshFormat\n2.2 1 8\n$EndMeshFormat\n$Nodes\n0\n$EndNodes\n")
    with pytest.raises(IndependentReadError, match="binary is not parsed"):
        read_gmsh22(path)


def test_gmsh22_rejects_element_tag_gap(tmp_path: Path) -> None:
    written = write_gmsh22(tmp_path / "cube.msh", _cube_surface())
    text = written.read_text().replace("$Elements\n12\n1 2", "$Elements\n12\n7 2", 1)
    broken = written.with_name("gap.msh")
    broken.write_text(text)
    with pytest.raises(IndependentReadError, match="1..N"):
        read_gmsh22(broken)


def test_gmsh22_rejects_node_count_lie(tmp_path: Path) -> None:
    written = write_gmsh22(tmp_path / "cube.msh", _cube_surface())
    text = written.read_text().replace("$Nodes\n8\n", "$Nodes\n9\n", 1)
    broken = written.with_name("nodelie.msh")
    broken.write_text(text)
    with pytest.raises(IndependentReadError, match="declared 9 nodes"):
        read_gmsh22(broken)


def test_vtu_rejects_point_count_lie(tmp_path: Path) -> None:
    written = write_vtu_ascii(tmp_path / "cube.vtu", _cube_surface())
    text = written.read_text().replace('NumberOfPoints="8"', 'NumberOfPoints="9"', 1)
    broken = written.with_name("lie.vtu")
    broken.write_text(text)
    with pytest.raises(IndependentReadError, match="point coordinates decoded"):
        read_vtk_xml(broken)


def test_vtu_rejects_appended_data(tmp_path: Path) -> None:
    path = tmp_path / "appended.vtu"
    path.write_text(
        '<?xml version="1.0"?>\n<VTKFile type="UnstructuredGrid" version="1.0">\n'
        '<UnstructuredGrid><Piece NumberOfPoints="3" NumberOfCells="1">\n'
        '<Points><DataArray type="Float64" Name="Points" NumberOfComponents="3" '
        'format="appended" offset="0"/></Points></Piece></UnstructuredGrid></VTKFile>\n'
    )
    with pytest.raises(IndependentReadError, match="appended"):
        read_vtk_xml(path)


def test_vtu_rejects_non_grid_cell_type(tmp_path: Path) -> None:
    """A quad (VTK type 9) must not be silently read as something else."""
    written = write_vtu_ascii(tmp_path / "cube.vtu", _cube_surface())
    text = written.read_text()
    start = text.index('Name="types" format="ascii">') + len('Name="types" format="ascii">')
    end = text.index("<", start)
    types = text[start:end].split()
    types[0] = "9"
    broken = written.with_name("quad.vtu")
    broken.write_text(text[:start] + " ".join(types) + text[end:])
    with pytest.raises(IndependentReadError, match="not one of tri/tetra"):
        read_vtk_xml(broken)


def test_vtu_rejects_short_compressed_payload(tmp_path: Path) -> None:
    """Regression: a wrong block layout used to decode to 2 values, not an error."""
    payload = np.arange(24, dtype=np.float64).tobytes()
    compressed = zlib.compress(payload)
    header = struct.pack("<4I", 1, 32768, len(payload), len(compressed))
    text = base64.b64encode(header).decode() + base64.b64encode(compressed).decode()
    good = (tmp_path / "compressed.vtu")
    good.write_text(
        '<?xml version="1.0"?>\n<VTKFile type="UnstructuredGrid" version="1.0" '
        'byte_order="LittleEndian" compressor="vtkZLibDataCompressor">\n'
        '<UnstructuredGrid><Piece NumberOfPoints="8" NumberOfCells="1">\n'
        '<Points><DataArray type="Float64" Name="Points" NumberOfComponents="3" '
        f'format="binary">{text}</DataArray></Points>\n'
        '<Cells><DataArray type="Int64" Name="connectivity" format="ascii">0 1 2 3</DataArray>'
        '<DataArray type="Int64" Name="offsets" format="ascii">4</DataArray>'
        '<DataArray type="UInt8" Name="types" format="ascii">10</DataArray></Cells>\n'
        "</Piece></UnstructuredGrid></VTKFile>\n"
    )
    parsed = read_vtk_xml(good)
    assert parsed.vertices.shape == (8, 3)
    assert np.allclose(parsed.vertices.reshape(-1), np.arange(24, dtype=np.float64))

    # now lie in the header about the block sizes: must raise, not truncate
    bad_header = struct.pack("<4I", 2, 8, len(payload), len(compressed))
    broken = tmp_path / "broken.vtu"
    broken.write_text(good.read_text().replace(
        base64.b64encode(header).decode(), base64.b64encode(bad_header).decode(), 1
    ))
    with pytest.raises(IndependentReadError):
        read_vtk_xml(broken)


def test_unknown_dialect_is_refused() -> None:
    with pytest.raises(IndependentReadError, match="no independent reader"):
        read_independent("whatever.xyz", file_format="xyz")


def test_foreign_identity_refuses_unknown_dialect(tmp_path: Path) -> None:
    from tensorcad.mesh_foreign_dialects import foreign_identity

    path = tmp_path / "x.bin"
    path.write_bytes(b"0000")
    with pytest.raises(IndependentReadError, match="no foreign identity"):
        foreign_identity(path, dialect="ply")


# ---------------------------------------------------------------------------
# 5. independence of the new modules
# ---------------------------------------------------------------------------

def test_foreign_module_does_not_import_meshio() -> None:
    code = (
        "import sys; import tensorcad.mesh_foreign_dialects, tensorcad.mesh_interop_independent; "
        "print('meshio' in sys.modules)"
    )
    out = subprocess.run(
        [sys.executable, "-c", code], cwd=ROOT, capture_output=True, text=True
    )
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip() == "False", out.stdout


def test_thirdparty_oracle_does_not_import_meshio() -> None:
    code = (
        "import sys; import tensorcad.mesh_interop_thirdparty; "
        "print('meshio' in sys.modules)"
    )
    out = subprocess.run(
        [sys.executable, "-c", code], cwd=ROOT, capture_output=True, text=True
    )
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip() == "False", out.stdout


def test_matrix_json_is_serialisable(foreign: dict) -> None:
    payload = json.dumps(foreign["rows"])
    assert len(payload) > 1000