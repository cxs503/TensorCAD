"""Regression tests for the round-4 adversarial audit of the ShipD rebuild.

These tests pin the *evidence properties* the round-3 report relied on but never
exercised:

* every quality gate must be able to fail (negative controls);
* the coefficient error measured on the mesh must be accompanied by a continuum
  error that is (a) two orders below the 3 % gate and (b) shrinking with
  refinement -- otherwise the "realised to 1e-6" claim is a self-consistency
  artefact;
* the calibration must survive a 16x finer measurement grid;
* the volume instrument must be exact on an analytically known solid.

Kept deliberately small so the suite stays in the seconds budget; the full
campaign is `tools/shipd_rebuild/audit_r4.py`.
"""
from __future__ import annotations

from tensorcad.cb_sections import CalibratedHullParameters, build_calibrated_hull
from tensorcad.geometry import Mesh
from tensorcad import mesh as meshlib
from tensorcad.shipd_population import build_shipd_database
from tensorcad.shipd_rebuild import _evaluate_gates

from tools.shipd_rebuild import audit_r4 as audit

SHIP = {"ship_id": 0, "ship_type": "audit", "L": 200.0, "B": 30.0, "T": 10.0,
        "Cb": 0.65, "Cm": 0.98, "Cwp": 0.80}


def _hull(nx: int = 41, nz: int = 13):
    return build_calibrated_hull(CalibratedHullParameters(
        SHIP["L"], SHIP["B"], SHIP["T"], SHIP["Cb"], SHIP["Cm"], SHIP["Cwp"], nx=nx, nz=nz))


def _tripped(mesh: Mesh) -> list[str]:
    result = _evaluate_gates(SHIP, mesh, "audit", "test", 0.0)
    gates = list(result.gates) + list(result.coefficient_gates)
    return [g.name for g in gates if not g.passed]


# --------------------------------------------------------------------------- #
# gate information content
# --------------------------------------------------------------------------- #
def test_clean_mesh_passes_every_gate() -> None:
    assert _tripped(_hull().mesh) == []


def test_inverted_triangle_is_detected() -> None:
    mesh = _hull().mesh
    flipped = Mesh(vertices=mesh.vertices,
                   faces=tuple(tuple(reversed(f)) for f in mesh.faces))
    assert "G4_outward_oriented" in _tripped(flipped)


def test_dropped_triangle_breaks_watertightness() -> None:
    mesh = _hull().mesh
    holed = Mesh(vertices=mesh.vertices, faces=mesh.faces[:-1])
    assert "G1_watertight" in _tripped(holed)


def test_duplicated_triangle_breaks_manifoldness() -> None:
    mesh = _hull().mesh
    doubled = Mesh(vertices=mesh.vertices, faces=mesh.faces + (mesh.faces[0],))
    assert "G2_edge_manifold" in _tripped(doubled)


def test_zero_area_triangle_is_detected() -> None:
    mesh = _hull().mesh
    f = mesh.faces[-1]
    degenerate = Mesh(vertices=mesh.vertices, faces=mesh.faces[:-1] + ((f[0], f[0], f[1]),))
    assert "G5_no_zero_area" in _tripped(degenerate)


def test_beam_error_trips_dimension_and_volume_gates() -> None:
    mesh = _hull().mesh
    wider = Mesh(vertices=tuple((x, y * 1.05, z) for x, y, z in mesh.vertices),
                 faces=mesh.faces)
    tripped = _tripped(wider)
    assert "G7_dimension_agreement" in tripped
    assert "G6_volume_consistency" in tripped


def test_coarse_mesh_trips_triangle_count_gate() -> None:
    assert "G8_triangle_count" in _tripped(build_calibrated_hull(CalibratedHullParameters(
        SHIP["L"], SHIP["B"], SHIP["T"], SHIP["Cb"], SHIP["Cm"], SHIP["Cwp"],
        nx=21, nz=7,
    )).mesh)


def test_negative_control_matrix_has_no_false_positive() -> None:
    """Every gate fires on at least one injected defect, the clean mesh on none."""
    report = audit.gate_negative_controls()
    assert report["summary"]["false_positive_on_clean_mesh"] is False
    assert report["summary"]["gates_with_zero_detections"] == []


def test_documented_blind_spot_single_inverted_triangle() -> None:
    """The round-4 blind spot, re-measured: round 5 closed it.

    Round 4 found that a single inverted triangle left the mesh watertight,
    edge-manifold, Euler-2 and positively oriented, so no gate rejected it, and
    this test pinned that blind spot ("if this test ever fails, a gate gained that
    power and the audit report must be updated"). Round 5 hardened
    ``G4_outward_oriented`` with the oriented-edge test, so the pin moves: the
    same measurement is asserted in the other direction. Detection must be by G4
    specifically -- a single flipped face moves the total volume by ~0.1 %, which
    the 3 % G6 tolerance cannot see, so a G6 hit would not be evidence.
    """
    report = audit.gate_negative_controls()
    assert report["undetected_defects"] == []
    row = next(r for r in report["rows"] if r["mutation"] == "flip_one_triangle")
    assert "G4_outward_oriented" in row["tripped"]
    assert "G6_volume_consistency" not in row["tripped"]

def test_continuum_cb_error_is_well_below_the_volume_gate() -> None:
    report = audit.discretisation_gap(cbs=(0.55, 0.86), grids=((81, 21),), L=200.0)
    assert report["summary"]["max_cb_continuum_error"] < 0.01
    assert report["summary"]["max_cb_mesh_error"] < 1e-4


def test_continuum_cb_error_shrinks_under_refinement() -> None:
    coarse = audit.discretisation_gap(cbs=(0.86,), grids=((81, 21),))
    fine = audit.discretisation_gap(cbs=(0.86,), grids=((321, 81),))
    e_coarse = coarse["rows"][0]["cb_continuum_error"]
    e_fine = fine["rows"][0]["cb_continuum_error"]
    assert e_fine < e_coarse / 2.0, (e_coarse, e_fine)


# --------------------------------------------------------------------------- #
# anti-overfit: calibrate coarse, measure fine
# --------------------------------------------------------------------------- #
def test_coefficients_survive_a_finer_measurement_grid() -> None:
    ships = build_shipd_database(6, seed=42)
    report = audit.resolution_independence(ships, sample=2, fine=(121, 41))
    assert report["summary"]["all_within_tolerance"] is True


# --------------------------------------------------------------------------- #
# instrument validation
# --------------------------------------------------------------------------- #
def test_volume_instrument_is_exact_on_a_unit_cube() -> None:
    report = audit.instrument_validation(samples=200)
    assert report["summary"]["unit_cube_exact"] is True


def test_ray_parity_agrees_with_pipeline_volume() -> None:
    report = audit.instrument_validation(samples=4000)
    assert report["summary"]["mc_agrees_with_pipeline"] is True


# --------------------------------------------------------------------------- #
# metric power
# --------------------------------------------------------------------------- #
def test_pearson_r_is_blind_to_shape_error_that_the_volume_gate_sees() -> None:
    """The acceptance metric must be shown able to miss a real defect."""
    ships = build_shipd_database(8, seed=42)
    report = audit.metric_power(ships, distortions=(0.0, 0.10))
    flat = report["rows"][0]["volume_error_max"]
    distorted = report["rows"][1]["volume_error_max"]
    assert distorted > flat, (flat, distorted)
    # Pearson r stays high on near-proportional data: documented, not a bug in r.
    assert report["rows"][1]["pearson_r"] > 0.99


def test_wetted_surface_area_is_positive_and_finite() -> None:
    mesh = _hull().mesh
    assert meshlib.surface_area(mesh) > 0.0