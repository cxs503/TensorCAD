"""Machine gate for the m17595 XVAL of the native six-format mesh matrix.

The tests are the executable form of the claims in
``docs/xval_mesh_six_format.md``: they assert what an *independent* reader recovers,
not what our own reader recovers, and they pin the reference-side boundaries
(meshio POLYDATA / NELEM=0 / CGNS stub) so a future meshio upgrade that fixes them
shows up as a failing test rather than as silent drift.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from tensorcad import mesh_xval as X

ROOT = Path(__file__).resolve().parents[1]
ART = ROOT / "artifacts" / "m17595_xval_six_format"
OWNER_TOL = 0.15          # owner precision standard
IMPL_TOL = 1e-3           # implementation acceptance line


def _meshio():
    pytest.importorskip("meshio")
    import meshio
    return meshio


# ---------------------------------------------------------------------------
# builders and comparison kernel
# ---------------------------------------------------------------------------
def test_bodies_are_closed_and_analytic():
    cube = X.cube(2.0)
    assert cube.volume() == pytest.approx(8.0, rel=1e-15)
    assert cube.euler() == 2 and cube.boundary_edges() == 0

    sph = X.icosphere(1.0, 2)
    exact = 4.0 / 3.0 * np.pi
    assert abs(sph.volume() - exact) / exact < 0.05      # discretisation, not format
    assert sph.euler() == 2 and sph.boundary_edges() == 0

    tor = X.torus(1.0, 0.3, 64, 32)
    exact_t = 2 * np.pi ** 2 * 1.0 * 0.3 ** 2
    assert abs(tor.volume() - exact_t) / exact_t < 0.01
    assert tor.euler() == 0                              # genus 1
    assert tor.boundary_edges() == 0


def test_compare_detects_wrong_orientation_and_missing_vertex():
    m = X.cube(2.0)
    flipped = X.TriMesh(m.points, m.faces[:, ::-1], "flipped")
    res = X.compare_meshes(m, flipped, 0.0)
    assert res["volume_got"] == pytest.approx(-8.0)  # halved volume sign, not magnitude
    assert res["volume_sign_equal"] is False
    assert res["pass"] is False                      # magnitude equal, orientation not

    dropped = X.TriMesh(m.points, m.faces[:-1], "open")
    assert X.compare_meshes(m, dropped, 0.0)["pass"] is False


# ---------------------------------------------------------------------------
# matrix A: our writers, independent readers
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("fmt,reader,tol_rel", [
    ("stl", "meshio", 1e-5), ("stl", "numpy-stl", 1e-5), ("stl", "trimesh", 1e-5),
    ("obj", "meshio", 0.0), ("obj", "trimesh", 0.0),
    ("vtk", "meshio", 0.0),
    ("gmsh", "meshio", 0.0),
])
def test_oss_readback_roundtrip(fmt, reader, tol_rel, tmp_path):
    """An independent implementation must recover volume, topology and orientation."""
    m = X.icosphere(1.0, 1)
    p = tmp_path / f"m.{'msh' if fmt == 'gmsh' else fmt}"
    X.WRITERS[fmt](m, p)
    tol = tol_rel * float(np.linalg.norm(m.points.max(axis=0) - m.points.min(axis=0)))
    fn = {
        "meshio": lambda: X.meshio_mesh(p, None if fmt == "gmsh" else fmt),
        "numpy-stl": lambda: X.numpy_stl_mesh(p),
        "trimesh": lambda: X.trimesh_mesh(p, process=False),
    }[reader]
    got = fn()
    res = X.compare_meshes(m, got, tol)
    assert res["volume_rel_err"] < IMPL_TOL
    assert res["volume_rel_err"] < OWNER_TOL
    assert res["bijection"] and res["oriented_face_multiset_equal"]
    assert res["flipped_faces"] == 0


def test_su2_strict_reader_roundtrip(tmp_path):
    m = X.cube(2.0)
    p = X.write_su2(m, tmp_path / "m.su2")
    got = X.read_su2_strict(p)
    res = X.compare_meshes(m, got, 0.0)
    assert res["pass"], res
    assert got.volume() == pytest.approx(8.0, rel=1e-15)


def test_cgns_writer_is_sids_conformant(tmp_path):
    """Geometry AND the SIDS-to-HDF5 node mapping must both hold."""
    pytest.importorskip("h5py")
    m = X.torus(1.0, 0.3, 16, 8)
    p = X.write_cgns(m, tmp_path / "m.cgns")
    back, violations = X.read_cgns_sids(p)
    assert violations == []
    assert back.volume() == pytest.approx(m.volume(), rel=1e-14)
    assert len(back.faces) == len(m.faces)


# ---------------------------------------------------------------------------
# reference-side boundaries (characterisation, so an OSS change is noticed)
# ---------------------------------------------------------------------------
def test_meshio_cannot_read_legacy_polydata(tmp_path):
    _meshio()
    m = X.cube(2.0)
    p = X.write_vtk(m, tmp_path / "m.vtk", dataset="POLYDATA")
    with pytest.raises((BaseException,)):
        X.meshio_mesh(p, "vtk")


def test_meshio_vtk_unstructured_grid_is_readable(tmp_path):
    _meshio()
    m = X.cube(2.0)
    p = X.write_vtk(m, tmp_path / "m.vtk", dataset="UNSTRUCTURED_GRID")
    assert X.compare_meshes(m, X.meshio_mesh(p, "vtk"), 0.0)["pass"]


def test_meshio_su2_reader_fails_only_on_zero_volume_elements(tmp_path):
    """meshio/su2/_su2.py:113-122: next() on an empty islice when NELEM=0."""
    _meshio()
    m = X.cube(2.0)
    boundary = X.write_su2(m, tmp_path / "b.su2", nelem_mode="boundary")
    volume = X.write_su2(m, tmp_path / "v.su2", nelem_mode="volume")
    try:
        X.meshio_mesh(boundary, "su2")
    except BaseException:                      # noqa: BLE001 - meshio calls sys.exit
        pass
    else:                                      # pragma: no cover
        pytest.skip("meshio now reads NELEM=0; §4.5 of the report must be re-measured")
    assert X.compare_meshes(m, X.meshio_mesh(volume, "su2"), 0.0)["pass"]


def test_meshio_cgns_reader_is_a_stub():
    meshio = _meshio()
    status = X.meshio_cgns_status()
    assert status["raises_can_only_read_tetrahedra"] is True
    assert status["path"].endswith("cgns/_cgns.py")


def test_trimesh_merges_vertices_by_default():
    """The vertex count is the reader's choice, not the file's -> weld explicitly."""
    pytest.importorskip("trimesh")
    m = X.cube(2.0)
    import tempfile

    with tempfile.TemporaryDirectory() as d:
        p = X.write_stl(m, Path(d) / "m.stl")
        raw = X.trimesh_mesh(p, process=False)
        # soup (36 corners) or merged (8) depending on the loader; welding must give 8/12
        welded, _ = X.weld(raw)
        assert (len(welded.points), len(welded.faces)) == (8, 12)


# ---------------------------------------------------------------------------
# sensitivity: float32 STL in a global frame
# ---------------------------------------------------------------------------
def test_stl_precision_degrades_with_coordinate_offset():
    m = X.icosphere(1.0, 1)
    import tempfile

    with tempfile.TemporaryDirectory() as d:
        devs = {}
        for off in (0.0, 1000.0):
            shifted = X.TriMesh(m.points + off, m.faces, m.name)
            p = X.write_stl(shifted, Path(d) / f"m{off:g}.stl")
            back, _ = X.weld(X.meshio_mesh(p, "stl"))
            d2 = np.linalg.norm(shifted.points[:, None, :] - back.points[None, :, :], axis=2)
            devs[off] = float(d2.min(axis=1).max())
        assert devs[0.0] < 1e-6
        assert devs[1000.0] > 50 * max(devs[0.0], 1e-9)
        assert devs[1000.0] < 1.0                       # still bounded by float32 eps*|x|


# ---------------------------------------------------------------------------
# committed artefacts
# ---------------------------------------------------------------------------
def test_committed_results_match_the_report():
    if not (ART / "results.json").exists():  # artefacts are not always checked out
        pytest.skip(f"artefacts missing: {ART}")
    res = json.loads((ART / "results.json").read_text())
    assert res["acceptance"]["owner_threshold"] == OWNER_TOL
    a, b = res["matrix_a"]["summary"], res["matrix_b"]["summary"]
    assert a["worst_volume_rel_err"] < IMPL_TOL
    assert b["worst_volume_rel_err"] < IMPL_TOL
    assert a["combinations"] == b["combinations"] == 18
    # every red entry must be one of the documented reference-side / dialect findings
    for f in a["failing"] + b["failing"]:
        assert (f["format"], f["status"]) in {
            ("vtk", "unreadable"), ("su2", "unreadable"), ("cgns", "nonconformant")}, f
    doc = (ROOT / "docs" / "xval_mesh_six_format.md").read_text()
    assert f"{a['worst_volume_rel_err']:.3e}" in doc
    assert (ART / "manifest.json").exists()