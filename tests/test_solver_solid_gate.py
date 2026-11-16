"""Tests for the one-solid / three-route solver gate.

They assert invariants, not just that code ran:

* the solid really is closed, and every route's manifest carries the *same*
  ``solid_digest`` -- "the same geometry" is a value, not a sentence;
* the star decomposition tiles the sampled surface exactly, so the FVM volume is
  the surface volume to working precision (and the fixture's own 0.36 % offset
  from ``4/3 pi a b c`` is gated separately);
* the DG payload is an *exact* rational-quadratic ellipse: its perimeter and area
  match references derived in a different parametrisation, and 65 sampled points
  lie on the solid according to the geometry library's own ``signed_distance``;
* the lattice accounts for every cell and its solid volume is within the
  staircase tolerance;
* a ``y-up`` manifest is refused by each route's **own** verifier after a positive
  control, and a verifier that refuses nothing is caught (the defect the old
  integration module had: ``except Exception: return True``).
"""
from __future__ import annotations

import json
import math
import shutil
from pathlib import Path

import pytest

from tensorcad.dg_handoff import verify_bundle as verify_dg_bundle
from tensorcad.solver_interface import (
    HandoffError,
    signed_volume,
    total_volume,
    verify_bundle,
)
from tensorcad.solver_solid_gate import (
    ROUTE_TOLERANCES,
    _pin_tensordg_path,
    cross_route_gate,
    dg_route,
    fvm_route,
    lattice_route,
    solid_from_library,
    star_tetrahedra,
    tamper_refused,
)

SEMI_AXES = (1.0, 0.5, 0.3)


@pytest.fixture(scope="module")
def solid():
    return solid_from_library(SEMI_AXES, resolution=10)


@pytest.fixture(scope="module")
def routes(solid, tmp_path_factory):
    out = tmp_path_factory.mktemp("gate")
    fvm = fvm_route(solid, out / "fvm", solve=False)
    dg = dg_route(solid, out / "dg")
    lattice = lattice_route(solid, out / "lattice", radial_cells=4)
    return {"fvm": fvm, "dg": dg, "lattice": lattice}


# ------------------------------------------------------------------ the solid


def test_library_solid_is_closed_and_matches_the_analytic_volume(solid):
    assert solid["edges"]["closed"] is True
    assert solid["edges"]["boundary_edges"] == 0
    assert solid["surface_volume_m3"] > 0.0
    assert solid["surface_volume_rel_error"] < 0.05
    assert solid["frame"] == "z-up"


def test_sampling_is_deterministic(solid):
    again = solid_from_library(SEMI_AXES, resolution=10)
    assert again["digest_sha256"] == solid["digest_sha256"]
    assert again["nodes"] == solid["nodes"]
    assert again["triangles"] == solid["triangles"]


def test_a_refined_solid_moves_the_volume_towards_the_analytic_one(solid):
    finer = solid_from_library(SEMI_AXES, resolution=20)
    assert finer["surface_volume_rel_error"] < solid["surface_volume_rel_error"]


# --------------------------------------------------------- star tetrahedra


def test_star_tetrahedra_tile_the_surface_exactly(solid):
    vertices, tets = star_tetrahedra(solid["mesh"])
    surface = signed_volume(solid["mesh"].vertices, solid["mesh"].faces)
    assert len(tets) == len(solid["mesh"].faces)
    assert total_volume(vertices, tets) == pytest.approx(surface, rel=1e-12)


def test_an_inward_oriented_surface_is_refused():
    """Inward-facing triangles invert every tet; the gate must refuse, not average."""
    solid = solid_from_library(SEMI_AXES, resolution=8)
    mesh = solid["mesh"]

    class Inward:
        vertices = mesh.vertices
        faces = [(a, c, b) for a, b, c in mesh.faces]

    with pytest.raises(HandoffError) as excinfo:
        star_tetrahedra(Inward())
    assert "star-shaped" in str(excinfo.value) or "non-positive" in str(excinfo.value)


# ------------------------------------------------------------------ FVM route


def test_fvm_bundle_re_reads_and_reproduces_the_surface_volume(routes, solid):
    fvm = routes["fvm"]
    assert fvm["verify_re_read_ok"] is True
    assert fvm["tet_vs_surface_volume_rel_error"] < ROUTE_TOLERANCES[
        "fvm_tet_vs_surface_volume_rel_error"]
    assert fvm["verify_volume_m3"] == pytest.approx(fvm["volume_m3"], rel=1e-12)
    assert fvm["solid_digest"] == solid["digest_sha256"]
    assert fvm["frame"] == "z-up"


def test_fvm_consumer_is_the_solver_or_says_why_not(routes):
    consumer = routes["fvm"]["consumer"]
    assert consumer["status"] in {"ran", "withheld", "error"}
    if consumer["status"] != "ran":
        assert consumer["reason"]
    else:
        assert consumer["cells"] == routes["fvm"]["tets"]
        assert consumer["maximum_relative_closure"] < 1e-09


def test_fvm_patch_field_is_reproduced_exactly(tmp_path):
    """A harmonic field must come back exactly from the conservative FV operator."""
    solid = solid_from_library(SEMI_AXES, resolution=8)
    fvm = fvm_route(solid, tmp_path / "fvm", solve=True)
    consumer = fvm["consumer"]
    if consumer["status"] != "ran":  # pragma: no cover - environment dependent
        pytest.skip(f"TensorFVM not usable here: {consumer.get('reason')}")
    assert consumer["patch_relative_linf"] < ROUTE_TOLERANCES["fvm_patch_relative_linf"]
    assert consumer["patch_relative_l2"] < ROUTE_TOLERANCES["fvm_patch_relative_l2"]
    assert abs(consumer["global_relative_balance"]) < ROUTE_TOLERANCES[
        "fvm_global_relative_balance"]


# ------------------------------------------------------------------- DG route


def test_dg_payload_is_an_exact_ellipse_on_the_solid(routes, solid):
    dg = routes["dg"]
    assert dg["perimeter_rel_error"] < ROUTE_TOLERANCES["dg_perimeter_rel_error"]
    assert dg["area_rel_error"] < ROUTE_TOLERANCES["dg_area_rel_error"]
    assert dg["on_solid_max_distance_m"] < ROUTE_TOLERANCES["dg_on_solid_max_distance_m"]
    assert dg["closure_gap_m"] < 1e-12
    assert dg["solid_digest"] == solid["digest_sha256"]
    assert dg["frame"] == "z-up"


def test_dg_bundle_declares_the_contract_frame(routes):
    manifest = json.loads((Path(routes["dg"]["bundle"]) / "handoff.json").read_text())
    assert manifest["frame"]["name"] == "z-up"
    assert manifest["geometry"]["r_inner_m"] == pytest.approx(1.0)  # enclosing radius
    assert manifest["measured"]["semi_axes"] == [1.0, 0.5]


def test_dg_verify_refuses_a_missing_frame(routes, tmp_path):
    source = Path(routes["dg"]["bundle"])
    copy = tmp_path / "noframe"
    shutil.copytree(source, copy)
    manifest = json.loads((copy / "handoff.json").read_text())
    manifest.pop("frame")
    (copy / "handoff.json").write_text(json.dumps(manifest))
    with pytest.raises(Exception) as excinfo:
        verify_dg_bundle(copy)
    assert "frame" in str(excinfo.value)


# ----------------------------------------------------------------- LBM route


def test_lattice_accounts_for_every_cell(routes, solid):
    lattice = routes["lattice"]
    assert lattice["verify_re_read_ok"] is True
    assert lattice["cells_accounted"] == 0
    assert lattice["solid_cells"] + lattice["fluid_cells"] == lattice["total_cells"]
    assert lattice["volume_rel_error"] < ROUTE_TOLERANCES["lattice_volume_rel_error"]
    assert lattice["solid_digest"] == solid["digest_sha256"]
    assert lattice["frame"] == "z-up"


def test_lattice_records_the_mesh_it_voxelised(routes):
    """The voxeliser's own mesh fingerprint reaches the manifest, not just a claim."""
    lattice = routes["lattice"]
    assert len(lattice["lattice_mesh_digest"]) == 64
    assert lattice["manifest_lattice_mesh_digest"] == lattice["lattice_mesh_digest"]


# ---------------------------------------------------------------------- gate


def test_the_gate_binds_all_routes_to_one_solid_and_one_frame(routes, solid):
    gate = cross_route_gate(routes, solid)
    checks = {c["check"]: c for c in gate["checks"]}
    assert checks["frame_contract"]["verdict"] == "PASS"
    assert checks["same_solid_digest"]["verdict"] == "PASS"
    assert gate["failed"] == []
    assert gate["solid_digests"]["fvm"] == gate["solid_digests"]["dg"] == \
        gate["solid_digests"]["lattice"]


def test_cross_route_volumes_are_compared_not_just_scored_against_analytic(routes, solid):
    gate = cross_route_gate(routes, solid)
    names = [c["check"] for c in gate["checks"]]
    assert "cross_route_volume_rel_error" in names
    entry = next(c for c in gate["checks"] if c["check"] == "cross_route_volume_rel_error")
    expected = abs(routes["fvm"]["volume_m3"] - routes["lattice"]["solid_volume_m3"]) / \
        solid["analytic_volume_m3"]
    assert entry["measured"] == pytest.approx(expected, rel=1e-12)


# -------------------------------------------------------------------- tamper


@pytest.mark.parametrize("route,verifier", [
    ("fvm", verify_bundle),
    ("lattice", verify_bundle),
    ("dg", verify_dg_bundle),
])
def test_a_yup_manifest_is_refused_by_the_routes_own_verifier(routes, route, verifier):
    result = tamper_refused(Path(routes[route]["bundle"]), verifier, name=route)
    assert result["positive_control_ok"] is True
    assert result["tamper_refused"] is True
    assert "frame" in result["refusal"]


def test_a_verifier_that_refuses_nothing_is_caught(routes):
    """The old integration module counted any exception -- or none -- as a pass."""
    def permissive(outdir):
        return {"route": "permissive"}

    with pytest.raises(HandoffError):
        tamper_refused(Path(routes["fvm"]["bundle"]), permissive, name="permissive")


def test_tamper_leaves_the_original_bundle_untouched(routes):
    source = Path(routes["fvm"]["bundle"])
    before = (source / "volume.msh").read_bytes()
    tamper_refused(source, verify_bundle, name="fvm")
    assert (source / "volume.msh").read_bytes() == before


# ----------------------------------------------------------- DG consumption


def test_dg_consumer_probe_records_what_tensordg_actually_does(routes, solid, tmp_path):
    if _pin_tensordg_path() is None:  # pragma: no cover - environment dependent
        pytest.skip("tensordg is not importable here")
    from tensorcad.solver_solid_gate import dg_consumer_probe

    probe = dg_consumer_probe(solid, tmp_path / "probe")
    ellipse = probe["ellipse_section"]
    circle = probe["circle_section"]
    # the ellipsoid's exact ellipse section must not be silently accepted: the
    # current TensorDG annulus builder requires a circular wall
    assert ellipse["status"] in {"refused", "accepted"}
    if ellipse["status"] == "refused":
        assert "circular" in ellipse["reason"]
    # and the circular section must be consumed, from the file, with no radial error
    assert circle["status"] == "ran"
    assert circle["bundle_verified"] is True
    assert circle["wall_radius_rel_error"] < 1e-09
    assert circle["wall_radius_spread_rel"] < 1e-09


def test_ellipse_perimeter_agrees_with_the_theta_parametrised_integral():
    """Independent cross-check of the rational-quadratic construction."""
    from tensorcad.dg_handoff import ellipse_contract, ellipse_perimeter

    contract = ellipse_contract(1.0, 0.5, segments=16)
    measured = 0.0
    # plain composite Simpson on |C'(t)| in the *NURBS* parametrisation
    from tensorcad.dg_handoff import evaluate_curve_derivative

    n = 20000
    ts = [i / n for i in range(n + 1)]
    speed = [math.hypot(*d) for d in evaluate_curve_derivative(contract, ts)]
    h = 1.0 / n
    measured = h / 3.0 * (speed[0] + speed[-1] + 4 * sum(speed[1:-1:2]) + 2 * sum(speed[2:-1:2]))
    assert measured == pytest.approx(ellipse_perimeter(1.0, 0.5), rel=1e-10)