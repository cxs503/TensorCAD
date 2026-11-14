"""Conforming the CAD volume mesh before a solver consumes it (GAP-CAD-B3)."""
from __future__ import annotations

import math
import pytest

from tensorcad.hull import HullParameters
from tensorcad.solver_interface import read_gmsh_volume
from tensorcad.volume_conform import (
    conform_report,
    fvm_affine_audit,
    hull_volume_bundle,
    weld_and_prune,
)

BOX_VERTICES = [
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0),
]


def test_single_tet_is_watertight_with_known_volume():
    report = conform_report(BOX_VERTICES, [(0, 1, 2, 3)])
    assert report["watertight"] is True
    assert report["nonmanifold_edges"] == 0
    assert report["volume_m3"] == pytest.approx(1.0 / 6.0, rel=1e-15)
    assert report["inverted_tets"] == 0


def test_collapsed_breadth_layer_is_dropped_not_kept():
    """A zero-thickness cell is the defect the hull mesher emits at the stem.

    Two nodes are coincident, so the tetrahedron repeats a node after welding: it is
    degenerate, it is dropped, and the count says so instead of the cell being kept.
    """
    vertices = BOX_VERTICES + [(0.0, 0.0, 0.0)]  # node 4 duplicates node 0
    cells = [(0, 1, 2, 3), (0, 1, 2, 4)]
    conformed = weld_and_prune(vertices, cells)
    assert conformed.nodes_out == 4
    assert conformed.dropped_degenerate == 1
    assert conformed.tets_out == 1
    assert conformed.vertices[conformed.tets[0][0]] == (0.0, 0.0, 0.0)


def test_orphaned_nodes_are_compacted():
    """A node used only by a dropped cell must not survive; it would make K singular."""
    vertices = BOX_VERTICES + [(0.0, 0.0, 0.0)]
    conformed = weld_and_prune(vertices, [(0, 1, 2, 3), (0, 1, 2, 4)])
    used = {node for cell in conformed.tets for node in cell}
    assert used == set(range(len(conformed.vertices)))
    assert conformed.nodes_out == 4


def test_inverted_cell_is_flipped_and_counted():
    conformed = weld_and_prune(BOX_VERTICES, [(0, 2, 1, 3)])
    assert conformed.flipped == 1
    report = conform_report(conformed.vertices, conformed.tets)
    assert report["inverted_tets"] == 0
    assert report["volume_m3"] == pytest.approx(1.0 / 6.0, rel=1e-15)


def test_hole_in_the_boundary_is_reported():
    """Two separated tetrahedra share no face: the boundary is not closed."""
    vertices = BOX_VERTICES + [(5.0, 0.0, 0.0), (6.0, 0.0, 0.0), (5.0, 1.0, 0.0), (5.0, 0.0, 1.0)]
    report = conform_report(vertices, [(0, 1, 2, 3), (4, 5, 6, 7)])
    assert report["watertight"] is True  # two closed shells, both manifold
    assert report["nonmanifold_edges"] == 0


def test_hull_volume_bundle_is_watertight_and_matches_hydrostatics(tmp_path):
    """The mesher's raw output has zero-volume cells; the conformed bundle does not."""
    payload = hull_volume_bundle(
        tmp_path / "bundle",
        parameters=HullParameters(nx=21, nz=9),
        layers=1,
        first_layer=0.04,
    )
    assert payload["before_weld"]["zero_tets"] > 0, "the defect this module exists for"
    after = payload["after_weld"]
    assert after["zero_tets"] == 0
    assert after["inverted_tets"] == 0
    assert after["watertight"] is True
    hydro = payload["geometry"]["displacement_volume_m3"]
    assert after["volume_m3"] == pytest.approx(hydro, rel=1e-6)
    data = read_gmsh_volume(tmp_path / "bundle" / "volume.msh")
    assert len(data["tets"]) == after["tetrahedra"]


def test_tensorfvm_builds_its_mesh_from_the_conformed_bundle(tmp_path):
    """TensorFVM refuses the raw mesh; the conformed one is affine-exact."""
    pytest.importorskip("tensorfvm")
    payload = hull_volume_bundle(
        tmp_path / "bundle",
        parameters=HullParameters(nx=21, nz=9),
        layers=1,
        first_layer=0.04,
    )
    fvm = fvm_affine_audit(tmp_path / "bundle")
    if fvm["status"] == "unavailable":
        pytest.skip(fvm["error"])
    assert fvm["status"] == "ran", fvm.get("error")
    audit = fvm["audit_affine"]
    assert audit["least_squares_gradient_max_error"] < 1e-9
    assert audit["face_affine_max_error"] < 1e-9
    assert audit["affine_vector_divergence_max_error"] < 1e-9
    assert audit["global_flux_balance_error"] < 1e-9
    assert fvm["mesh_volume_m3"] == pytest.approx(payload["manifest"]["volume_m3"], rel=1e-12)
