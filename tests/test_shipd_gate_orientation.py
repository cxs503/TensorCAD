"""Negative controls for the rebuild gate set: does G4 actually see an inverted patch?

The round-3/round-4 acceptance statistic is "N gates passed on N ships". A gate set
is only an instrument if each gate fires on the defect it claims to detect, so these
tests inject defects into an otherwise clean mesh and check two things at once:

* the defective mesh *is* caught by the gate that claims that property, and
* the other gates are documented as blind to it (that is the reason the check had to
  be added -- a solely volume-based "outward oriented" gate passed a mesh with one
  inverted triangle).

Nothing here is a tolerance tweak: ``G4_outward_oriented`` now requires both a
positive divergence-theorem volume and agreement of every shared edge's direction
between its two faces.
"""
import math

import pytest

from tensorcad import mesh as meshlib
from tensorcad.cb_sections import CalibratedHullParameters, build_calibrated_hull
from tensorcad.geometry import Mesh
from tensorcad.shipd_rebuild import _evaluate_gates


def _hull_ship():
    return {
        "ship_id": 0,
        "ship_type": "container",
        "L": 200.0, "B": 30.0, "T": 10.0,
        "Cb": 0.65, "Cm": 0.98, "Cwp": 0.80,
        "LCB": 0.5,
    }


def _clean_mesh(nx=41, nz=13):
    ship = _hull_ship()
    params = CalibratedHullParameters(
        length=ship["L"], beam=ship["B"], draft=ship["T"],
        block_coefficient=ship["Cb"], midship_coefficient=ship["Cm"],
        waterplane_coefficient=ship["Cwp"], nx=nx, nz=nz,
    )
    return ship, build_calibrated_hull(params).mesh


def _flip(mesh: Mesh, indices) -> Mesh:
    faces = list(mesh.faces)
    for i in indices:
        a, b, c = faces[i]
        faces[i] = (a, c, b)
    return Mesh(vertices=mesh.vertices, faces=faces, units=mesh.units,
                coordinate_frame=mesh.coordinate_frame)


def _gates(ship, mesh):
    return {g.name: g for g in _evaluate_gates(ship, mesh, "cb_calibrated_sections", "-", 0.0).gates}


class TestCleanMesh:
    def test_generator_is_consistently_wound(self):
        _, mesh = _clean_mesh()
        assert meshlib.edge_winding_conflicts(mesh) == ()
        assert meshlib.is_consistently_oriented(mesh)

    def test_all_eight_gates_pass(self):
        ship, mesh = _clean_mesh()
        gates = _gates(ship, mesh)
        failed = {name: g.detail for name, g in gates.items() if not g.passed}
        assert failed == {}
        assert len(gates) == 8


class TestSingleInvertedTriangle:
    """One flipped face: the documented blind spot of a volume-sign-only G4."""

    def test_winding_test_sees_it(self):
        _, mesh = _clean_mesh()
        bad = _flip(mesh, [1000])
        # The three edges of the flipped triangle now disagree with their
        # neighbours; the interior edges of the patch are unaffected.
        assert len(meshlib.edge_winding_conflicts(bad)) == 3

    def test_g4_fails(self):
        ship, mesh = _clean_mesh()
        gates = _gates(ship, _flip(mesh, [1000]))
        assert not gates["G4_outward_oriented"].passed
        assert "winding-inconsistent" in gates["G4_outward_oriented"].detail

    def test_every_other_gate_is_blind_to_it(self):
        """This is why the directed-edge test was added, not a tolerance change."""
        ship, mesh = _clean_mesh()
        clean = _gates(ship, mesh)
        bad = _gates(ship, _flip(mesh, [1000]))
        for name in ("G1_watertight", "G2_edge_manifold", "G3_euler_char",
                     "G5_no_zero_area", "G6_volume_consistency",
                     "G7_dimension_agreement", "G8_triangle_count"):
            assert bad[name].passed, f"{name} unexpectedly fired on one flipped face"
        # The volume gate stays inside its 3 % band: flipping one face perturbs the
        # divergence-theorem sum by one tetrahedron (measured 1.05e-3 relative on
        # this 41x13 hull, 3000x the clean mesh's 2.2e-15 but 28x inside the gate).
        assert clean["G6_volume_consistency"].value < 1e-12
        assert bad["G6_volume_consistency"].value < 0.03
        assert bad["G6_volume_consistency"].value > 1e-6
        assert bad["G4_outward_oriented"].value > 0

    def test_repair_restores_agreement(self):
        """``orient_outward`` is the repair path and closes the defect."""
        ship, mesh = _clean_mesh()
        bad = _flip(mesh, [1000])
        repaired = meshlib.orient_outward(bad)
        assert repaired.flipped_faces >= 1
        assert meshlib.edge_winding_conflicts(repaired.mesh) == ()
        assert _gates(ship, repaired.mesh)["G4_outward_oriented"].passed


class TestLargerDefects:
    def test_inverted_patch(self):
        ship, mesh = _clean_mesh()
        patch = range(500, 540)
        bad = _flip(mesh, patch)
        gates = _gates(ship, bad)
        assert len(meshlib.edge_winding_conflicts(bad)) > 0
        assert not gates["G4_outward_oriented"].passed
        assert gates["G1_watertight"].passed  # topology survives, orientation does not

    def test_globally_inverted_mesh(self):
        """The other half of G4: consistent but inward -> caught by the sign."""
        ship, mesh = _clean_mesh()
        inward = _flip(mesh, range(len(mesh.faces)))
        assert meshlib.edge_winding_conflicts(inward) == ()  # globally consistent
        gates = _gates(ship, inward)
        assert not gates["G4_outward_oriented"].passed
        assert gates["G4_outward_oriented"].value < 0

    def test_both_halves_are_needed(self):
        """Neither test alone covers both defects; the gate requires both."""
        ship, mesh = _clean_mesh()
        inward = _flip(mesh, range(len(mesh.faces)))
        # sign test alone would pass the single-face defect ...
        assert meshlib.total_signed_volume(_flip(mesh, [1000])) > 0
        # ... and the winding test alone would pass the global inversion.
        assert meshlib.edge_winding_conflicts(inward) == ()


class TestPopulationMeshes:
    """The production generator must not need the new test to pass by luck."""

    @pytest.mark.parametrize("cb", [0.45, 0.60, 0.72, 0.86])
    def test_type_representative_hulls_are_consistent(self, cb):
        params = CalibratedHullParameters(
            length=150.0, beam=25.0, draft=8.0,
            block_coefficient=cb,
            midship_coefficient=1.0 - 0.20 * (1.0 - cb) ** 1.8,
            waterplane_coefficient=(1.0 + 2.0 * cb) / 3.0,
            nx=41, nz=13,
        )
        mesh = build_calibrated_hull(params).mesh
        assert meshlib.edge_winding_conflicts(mesh) == ()
        assert meshlib.total_signed_volume(mesh) > 0
        assert math.isfinite(meshlib.total_signed_volume(mesh))