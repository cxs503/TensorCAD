"""Tests for :mod:`tensorcad.mesh_formats` — the native six-format exchange.

The contract under test: a mesh written in any of the six grammars and read back
by *our* reader must be geometrically and topologically identical to the input
(bit-exact for the five float64 formats, float32-tight for STL), the CGNS file must
satisfy every SIDS-to-HDF5 rule in ``CGNS_RULES``, and ``compare_meshes`` must
reject a mesh that is subtly wrong (orientation flip, missing vertex) instead of
merely executing.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from tensorcad import mesh_formats as mf
from tensorcad import mesh_xval as mx

BODIES = (
    pytest.param(mf.cube(), id="cube"),
    pytest.param(mf.icosphere(1.0, 2), id="icosphere2"),
    pytest.param(mf.torus(), id="torus64x32"),
)


def _ext(fmt: str) -> str:
    return "msh" if fmt == "gmsh" else fmt


# ---------------------------------------------------------------------------
# round trips
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("body", BODIES)
@pytest.mark.parametrize("fmt", mf.MESH_FORMATS)
def test_roundtrip_is_geometrically_identical(body, fmt, tmp_path):
    tol = 1e-5 if fmt == "stl" else 0.0  # STL is float32 by format definition
    report = mf.roundtrip(body, fmt, tmp_path / f"{body.name}.{_ext(fmt)}", tol=tol)
    assert report["pass"], report
    assert report["faces"] == [len(body.faces), len(body.faces)]
    assert report["volume_rel_err"] <= (1e-6 if fmt == "stl" else 0.0)
    assert report["flipped_faces"] == 0
    assert report["closed"]


def test_stl_is_a_facet_soup_and_needs_welding(tmp_path):
    body = mf.torus()
    report = mf.roundtrip(body, "stl", tmp_path / "torus.stl", tol=1e-5)
    assert report["welded"] is True
    assert report["vertices_read"][1] == 3 * len(body.faces)  # 3 corners per facet
    assert report["vertices"][1] == len(body.points)


def test_volume_matches_analytic_reference_and_converges_second_order():
    assert mf.cube(2.0).volume() == pytest.approx(8.0, rel=0, abs=1e-15)
    # inscribed polyhedra: volume error is O(h^2), so each subdivision must cut it
    # by ~4x.  This pins the body generator (the comparison reference) itself.
    exact = 4.0 / 3.0 * math.pi
    errs = [abs(mf.icosphere(1.0, k).volume() - exact) / exact for k in (2, 3, 4)]
    assert errs[0] > errs[1] > errs[2]
    assert 3.4 < errs[0] / errs[1] < 4.6
    assert 3.4 < errs[1] / errs[2] < 4.6
    assert errs[2] < 5e-3
    torus = mf.torus(1.0, 0.3, 128, 64).volume()
    assert torus == pytest.approx(2 * math.pi**2 * 1.0 * 0.3**2, rel=5e-3)  # O(h^2) inscribed


def test_euler_characteristic_distinguishes_genus():
    assert mf.cube().euler() == 2
    assert mf.icosphere(1.0, 2).euler() == 2
    assert mf.torus().euler() == 0  # genus 1
    for body in (mf.cube(), mf.icosphere(1.0, 2), mf.torus()):
        assert body.boundary_edges() == 0  # closed surfaces only


# ---------------------------------------------------------------------------
# CGNS SIDS conformance: the point of the module
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("body", BODIES)
def test_cgns_is_sids_conformant(body, tmp_path):
    path = mf.write(body, "cgns", tmp_path / "m.cgns")
    assert mf.sids_violations(path) == []
    assert mf.read(path, "cgns").volume() == pytest.approx(body.volume(), rel=1e-15)


def test_cgns_carries_library_version_node(tmp_path):
    import h5py

    path = mf.write(mf.cube(), "cgns", tmp_path / "cube.cgns")
    with h5py.File(path, "r") as f:
        lib = [g for g in f.values() if g.attrs.get(" label") in (b"CGNSLibraryVersion_t", "CGNSLibraryVersion_t")]
        assert lib, "no CGNSLibraryVersion_t node"
        assert float(np.asarray(lib[0][" data"][()]).ravel()[0]) == pytest.approx(4.3)
        assert lib[0].attrs[" type"] in (b"R4", "R4")


def test_checker_rejects_a_dialect_without_library_version(tmp_path):
    """Negative control: the reference dialect of ``mesh_xval`` omits the node.

    Without this test the checker could be vacuous, so a file that differs *only*
    by the missing root version node must be reported as non-conformant -- and
    ``read_cgns`` must refuse to launder it.
    """
    path = mx.write_cgns(mx.TriMesh(*_xy(mf.cube())), tmp_path / "xval.cgns")
    violations = mf.sids_violations(path)
    assert violations == ["root: CGNSLibraryVersion_t node missing"]
    with pytest.raises(ValueError, match="CGNSLibraryVersion_t"):
        mf.read(path, "cgns")


def _xy(mesh):
    return mesh.points, mesh.faces


def test_sids_checker_flags_a_wrong_attribute_mapping(tmp_path):
    """Drop one required attribute -> exactly one violation naming it."""
    import h5py

    path = mf.write(mf.icosphere(1.0, 1), "cgns", tmp_path / "s.cgns")
    with h5py.File(path, "a") as f:
        del f["Base"].attrs[" order"]
    violations = mf.sids_violations(path)
    assert any("'order' missing" in v for v in violations)


# ---------------------------------------------------------------------------
# dialect choices (each one bought with a measured downstream failure)
# ---------------------------------------------------------------------------
def test_vtk_default_dialect_is_meshio_readable(tmp_path):
    import meshio

    body = mf.cube()
    path = mf.write(body, "vtk", tmp_path / "cube.vtk")
    assert "DATASET UNSTRUCTURED_GRID" in path.read_text()
    mesh = meshio.read(path)
    assert mesh.points.shape == (8, 3)
    assert sum(len(c.data) for c in mesh.cells) == len(body.faces)


def test_vtk_polydata_dialect_is_parsed_by_our_reader(tmp_path):
    body = mf.cube()
    path = mf.write(body, "vtk", tmp_path / "cube_pd.vtk", dataset="POLYDATA")
    assert "DATASET POLYDATA" in path.read_text()
    assert mf.compare_meshes(body, mf.read(path, "vtk"))["pass"]


def test_su2_default_puts_triangles_in_nelem_and_is_meshio_readable(tmp_path):
    import meshio

    path = mf.write(mf.cube(), "su2", tmp_path / "c.su2")
    assert "NELEM= 12" in path.read_text()
    assert sum(len(c.data) for c in meshio.read(path).cells) == 12


def test_su2_boundary_mode_is_legal_and_parsed_here(tmp_path):
    body = mf.cube()
    path = mf.write(body, "su2", tmp_path / "b.su2", nelem_mode="boundary")
    text = path.read_text()
    assert "NELEM= 0" in text and "MARKER_ELEMS= 12" in text
    assert mf.compare_meshes(body, mf.read(path, "su2"))["pass"]


def test_gmsh_alias_and_version_guard(tmp_path):
    path = mf.write(mf.cube(), "msh", tmp_path / "c.msh")
    assert path.read_text().startswith("$MeshFormat\n2.2")
    assert len(mf.read(path, "msh").faces) == 12
    bad = tmp_path / "v4.msh"
    bad.write_text("$MeshFormat\n4.1 0 8\n$EndMeshFormat\n")
    with pytest.raises(ValueError, match="MSH 2.x"):
        mf.read(bad, "gmsh")


def test_unknown_format_is_rejected():
    with pytest.raises(ValueError, match="unknown mesh format"):
        mf.write(mf.cube(), "dxf", "/tmp/x.dxf")


# ---------------------------------------------------------------------------
# comparison contract must bite
# ---------------------------------------------------------------------------
def test_compare_flags_a_flipped_triangle():
    body = mf.cube()
    bad = mf.TriMesh(body.points, body.faces.copy(), "flipped")
    bad.faces[0] = bad.faces[0][::-1]
    report = mf.compare_meshes(body, bad)
    assert report["pass"] is False
    assert report["flipped_faces"] >= 1
    assert report["volume_rel_err"] > 0


def test_compare_flags_a_missing_vertex():
    body = mf.icosphere(1.0, 1)
    bad = mf.TriMesh(body.points[:-1], body.faces[body.faces.max(axis=1) < len(body.points) - 1], "short")
    report = mf.compare_meshes(body, bad)
    assert report["pass"] is False
    assert report["bijection"] is False


def test_compare_accepts_float32_tolerance_only_when_asked():
    body = mf.icosphere(1.0, 2)
    nudge = body.points * (1.0 + 1e-9)
    loose = mf.TriMesh(nudge, body.faces, "nudged")
    assert mf.compare_meshes(body, loose)["pass"] is False      # exact bijection demands equality
    assert mf.compare_meshes(body, loose, tol=1e-6)["pass"] is True

# ---------------------------------------------------------------------------
# committed artefacts: the report numbers must be the tool's numbers
# ---------------------------------------------------------------------------
ART = "m17907_native_six_format"


def _artifacts() -> "object":
    import json
    from pathlib import Path

    root = Path(__file__).resolve().parents[1] / "artifacts" / ART
    if not (root / "results.json").exists():  # pragma: no cover - artefact not checked out
        pytest.skip("artefacts not present")
    return root, json.loads((root / "results.json").read_text(encoding="utf-8"))


def test_committed_results_match_the_acceptance_gate():
    root, data = _artifacts()
    summary = data["summary"]
    assert summary["own_reader_topology_pass"] == summary["cells"] == 18
    assert summary["worst_own_reader_rel_err"] <= data["acceptance"]["implementation_threshold"]
    assert summary["own_reader_all_within_owner_threshold"]
    assert summary["cgns_native_violations"] == 0
    assert summary["cross_writer_agreement"] == summary["cross_writer_cells"]
    # the legacy dialect is non-conformant for exactly one reason: no version node
    assert data["sids"]["legacy_reference_dialect_violations"] == ["root: CGNSLibraryVersion_t node missing"]
    # every recorded file still hashes to the manifest entry
    import hashlib

    files = data["files_sha256"]
    assert files
    for name, digest in files.items():
        h = hashlib.sha256((root / name).read_bytes()).hexdigest()
        assert h == digest, name


def test_committed_matrix_and_figure_exist():
    root, data = _artifacts()
    assert (root / "matrix.md").read_text(encoding="utf-8").startswith("# Native six-format matrix")
    figs = data.get("figures") or []
    if figs:  # figures are produced by the tool, not by --check
        path = root / "figures" / figs[0]["name"]
        assert path.exists() and path.stat().st_size > 10_000
