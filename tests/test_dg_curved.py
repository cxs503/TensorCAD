"""Invariants of the CAD -> DG curved-element handoff.

Every test asserts a *property* (analytic value, refusal of a malformed
contract, convergence of a measured quantity), never "the code ran".
"""
from __future__ import annotations

import math
import subprocess
import sys
from pathlib import Path

import pytest

from tensorcad import dg_curved as dgc
from tensorcad.nurbs import cylinder_solid_mesh
from tensorcad.solver_interface import ellipsoid_mesh

ROOT = Path(__file__).resolve().parents[1]


def _cylinder_contract(u_samples: int = 96, v_samples: int = 6, radius: float = 1.0):
    mesh = cylinder_solid_mesh(radius, 2.0, u_samples=u_samples, v_samples=v_samples)
    contour = dgc.section_contour(mesh, 0.86)
    fit = dgc.fit_circle(contour)
    contract = dgc.rational_circle_contract(
        (fit["centre_x"], fit["centre_y"]), fit["radius"])
    return contour, fit, contract


# ---------------------------------------------------------------------------
# sectioning
# ---------------------------------------------------------------------------
def test_section_of_cylinder_is_a_closed_ring_on_the_surface():
    mesh = cylinder_solid_mesh(1.0, 2.0, u_samples=48, v_samples=6)
    contour = dgc.section_contour(mesh, 0.86)
    assert len(contour) >= 48
    radii = [math.hypot(x, y) for x, y in contour]
    # the slice cuts the triangulation, so every point sits on a facet: the
    # deviation from the exact radius is the O(h^2) sagitta, not zero.
    h = 2.0 * math.pi / 48
    assert max(abs(r - 1.0) for r in radii) < 0.25 * h * h


def test_section_refuses_ambiguous_plane_through_vertices():
    mesh = cylinder_solid_mesh(1.0, 2.0, u_samples=16, v_samples=4)
    with pytest.raises(dgc.ContractError, match="passes through a mesh vertex"):
        dgc.section_contour(mesh, 1.0)   # a lateral vertex ring sits at z = 1.0


def test_section_refuses_plane_outside_the_solid():
    mesh = cylinder_solid_mesh(1.0, 2.0, u_samples=16, v_samples=4)
    with pytest.raises(dgc.ContractError, match="outside the mesh extent"):
        dgc.section_contour(mesh, 5.0)


# ---------------------------------------------------------------------------
# circle fit + contract representation
# ---------------------------------------------------------------------------
def test_fit_recovers_an_exact_circle_to_round_off():
    pts = [(1.7 + 0.8 * math.cos(t), -0.3 + 0.8 * math.sin(t))
           for t in [2.0 * math.pi * i / 64 for i in range(64)]]
    fit = dgc.fit_circle(pts)
    assert abs(fit["centre_x"] - 1.7) < 1e-12
    assert abs(fit["centre_y"] + 0.3) < 1e-12
    assert abs(fit["radius"] - 0.8) < 1e-12
    assert fit["max_residual"] < 1e-12


def test_emitted_nurbs_is_the_declared_circle():
    _, _, contract = _cylinder_contract()
    report = dgc.curve_report(contract, samples=2048,
                              analytic_centre=(0.0, 0.0), analytic_radius=1.0)
    # representation invariant: rational quadratic NURBS circle == the circle
    assert report["self_consistency_max_m"] / report["declared_radius_m"] < 1e-15
    assert report["closure_gap"] < 1e-15
    # the *fitted* radius carries the section faceting (~1e-4 at 96 segments),
    # so the area error against the analytic circle is that, not round-off
    assert abs(report["area_rel_error"]) < 1e-3
    assert abs(report["radial_deviation_max_m"]) < 1e-3


def test_area_and_perimeter_are_sampling_limited_not_geometry_limited():
    _, _, contract = _cylinder_contract()
    coarse = dgc.curve_report(contract, samples=512, analytic_centre=(0.0, 0.0),
                              analytic_radius=1.0)
    fine = dgc.curve_report(contract, samples=2048, analytic_centre=(0.0, 0.0),
                            analytic_radius=1.0)
    assert abs(fine["area_rel_error"]) < abs(coarse["area_rel_error"])
    assert abs(fine["perimeter_rel_error"]) < abs(coarse["perimeter_rel_error"])


def test_contract_inherits_the_faceting_error_of_the_cad_section():
    fine_contour, fine_fit, fine_contract = _cylinder_contract(u_samples=192)
    coarse_contour, coarse_fit, coarse_contract = _cylinder_contract(u_samples=48)
    fine = dgc.curve_report(fine_contract, samples=2048,
                            analytic_centre=(0.0, 0.0), analytic_radius=1.0,
                            section_points=fine_contour)
    coarse = dgc.curve_report(coarse_contract, samples=2048,
                              analytic_centre=(0.0, 0.0), analytic_radius=1.0,
                              section_points=coarse_contour)
    assert fine["radial_deviation_max_m"] < coarse["radial_deviation_max_m"]
    ratio = (coarse["section_point_circle_residual_max_m"]
             / fine["section_point_circle_residual_max_m"])
    assert ratio > 10.0, f"faceting error did not converge at O(h^2), ratio={ratio}"


# ---------------------------------------------------------------------------
# validation: every refusal has a reason
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("mutate, message", [
    (lambda c: c.update(schema="tensordg.curved_geometry/2"), "unsupported schema"),
    (lambda c: c.update(kind="nurbs_surface"), "unsupported geometry kind"),
    (lambda c: c.pop("units"), "declare metres"),
    (lambda c: c.update(units="mm"), "declare metres"),
    (lambda c: c.update(frame="y-up"), "z-up frame"),
    (lambda c: c.pop("knots"), "missing field"),
    (lambda c: c.update(weights=[-1.0] * 9), "strictly positive"),
    (lambda c: c.update(knots=c["knots"][:-1]), "knot count"),
    (lambda c: c.update(knots=[1.0, 0.0] + c["knots"][2:]), "non-decreasing"),
    (lambda c: c.update(control_points=c["control_points"][:-1] + [[9.0, 9.0]]),
     "identical first/last control point"),
    (lambda c: c.update(degree=0), "degree must be >= 1"),
])
def test_validate_refuses_malformed_contracts(mutate, message):
    _, _, contract = _cylinder_contract(u_samples=32)
    broken = dict(contract)
    mutate(broken)
    with pytest.raises(dgc.ContractError, match=message):
        dgc.validate_contract(broken)


def test_handoff_writes_contract_and_record(tmp_path):
    _, _, contract = _cylinder_contract()
    report = dgc.curve_report(contract, samples=512)
    record = dgc.write_handoff(tmp_path, contract, report, extra={"k": "v"})
    assert (tmp_path / "contract.json").exists()
    assert record["frame"] == "z-up" and record["units"] == "m"
    assert record["k"] == "v"
    assert record["geometry"]["contract_digest"] == report["contract_digest"]


# ---------------------------------------------------------------------------
# non-circular section: general path + honest consumer boundary
# ---------------------------------------------------------------------------
def test_closed_bspline_contract_tracks_an_elliptical_section():
    ellipsoid = ellipsoid_mesh((1.2, 0.9, 0.8), n_theta=48, n_phi=96)
    contour = dgc.section_contour(ellipsoid, 0.08)
    contract = dgc.bspline_closed_contract(contour, degree=3, n_control=16)
    dgc.validate_contract(contract)          # structurally a valid CAD contract
    assert contract["control_points"][0] == contract["control_points"][-1]
    report = dgc.curve_report(contract, samples=2048, section_points=contour)
    # the least-squares fit tracks the section to well under a percent of a
    a = 1.2 * math.sqrt(1.0 - (0.08 / 0.8) ** 2)
    # measured: 0.83 % of the semi-axis for a 16-control-point fit
    assert report["section_point_deviation_max_m"] < 2e-2 * a
    spread = dgc.fit_circle(contour)["max_residual"] / dgc.fit_circle(contour)["radius"]
    assert spread > 0.1, "an ellipse must not pass as a circle"


# ---------------------------------------------------------------------------
# consumer side (skipped when TensorDG is not importable in this environment)
# ---------------------------------------------------------------------------
def test_dg_consumes_the_cad_circle_contract():
    dgc_mod = pytest.importorskip("tensordg.curved_geometry")
    _, fit, contract = _cylinder_contract()
    validated = dgc_mod.validate_geometry_contract(contract)
    mesh = dgc_mod.mesh_from_contract(validated, 2, 16, 2.0, 2)
    radii = []
    for mapping in mesh.flat_maps():
        x, y = mapping.eval(__import__("numpy").array([-1.0]),
                            __import__("numpy").array([0.0]))
        radii.append(float(__import__("numpy").mean(
            __import__("numpy").hypot(__import__("numpy").ravel(x),
                                      __import__("numpy").ravel(y)))))
    assert min(radii) > fit["radius"] - 1e-9


def test_dg_refuses_the_non_circular_contract_with_a_reason():
    dgc_mod = pytest.importorskip("tensordg.curved_geometry")
    ellipsoid = ellipsoid_mesh((1.2, 0.9, 0.8), n_theta=32, n_phi=64)
    contour = dgc.section_contour(ellipsoid, 0.08)
    contract = dgc.bspline_closed_contract(contour, degree=3, n_control=16)
    validated = dgc_mod.validate_geometry_contract(contract)
    with pytest.raises(ValueError, match="not circular"):
        dgc_mod.mesh_from_contract(validated, 2, 16, 2.0, 2)


# ---------------------------------------------------------------------------
# the core lane stays standard-library only
# ---------------------------------------------------------------------------
def test_module_imports_no_third_party_code():
    """``dg_curved`` lives in the standard-library lane of the package.

    The aggregate ``import tensorcad`` pulls in NumPy-backed modules
    (``mesh_formats``), so the property asserted here is the one that actually
    belongs to this module: its own imports are stdlib plus
    :mod:`tensorcad.geometry`.
    """
    source = (ROOT / "tensorcad" / "dg_curved.py").read_text()
    tree = __import__("ast").parse(source)
    imported = set()
    for node in __import__("ast").walk(tree):
        if isinstance(node, __import__("ast").Import):
            imported.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, __import__("ast").ImportFrom) and node.level == 0:
            imported.add((node.module or "").split(".")[0])
    assert imported <= {"__future__", "hashlib", "json", "math",
                        "collections", "pathlib", "typing"}
