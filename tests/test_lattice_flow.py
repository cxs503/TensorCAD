"""Tests for the CAD -> lattice -> 3-D LBM duct route (:mod:`tensorcad.lattice_flow`).

What these pin, and why each one is here rather than "it ran":

* the closed form is arithmetic, so a wrong formula cannot hide behind a solver;
* the lattice that leaves the CAD side is the contract's lattice (route, frame,
  cell ordering) and its staircase area is a measured, converging quantity;
* the flow solver returns the *shape* of Poiseuille flow (centre maximum, monotone
  towards the wall, no transverse motion) and the right magnitude;
* the mask arithmetic that makes bounce-back periodic is exercised on a wrap case
  that a wrong `%` would silently pass;
* a bundle of the wrong route is refused before any solver state is allocated.
"""
from __future__ import annotations

import json
import math

import pytest

torch = pytest.importorskip("torch")
pytest.importorskip("tensorlbm")

from tensorcad.lattice_flow import (  # noqa: E402
    duct_bundle,
    duct_lattice,
    gate_order,
    lattice_flow_from_bundle,
    poiseuille_reference,
    stream3d_probe,
    _shifted_wall_mask,
)
from tensorcad.solver_interface import (  # noqa: E402
    COORDINATE_CONTRACT,
    HandoffError,
    read_lattice_table,
    verify_bundle,
)

RADIUS = 0.05
BODY_FORCE = 2.0e-5
TAU = 0.8


def _small_case(tmp_path, radial_cells=5, axial_cells=2):
    bundle = tmp_path / f"duct_r{radial_cells}"
    duct_bundle(bundle, radius=RADIUS, radial_cells=radial_cells, axial_cells=axial_cells)
    return bundle, lattice_flow_from_bundle(
        bundle, tau=TAU, body_force=BODY_FORCE, max_steps=6000, tol=1e-13
    )


def test_closed_form_is_poiseuille():
    ref = poiseuille_reference(2.0e-5, 8.0, 0.1)
    assert ref["u_max"] == pytest.approx(2.0e-5 * 64 / (4 * 0.1))
    assert ref["u_mean"] == pytest.approx(ref["u_max"] / 2.0)
    assert ref["flow_rate"] == pytest.approx(math.pi * ref["u_max"] * 64 / 2.0)


def test_reference_refuses_degenerate_input():
    with pytest.raises(ValueError):
        poiseuille_reference(1.0, 0.0, 0.1)
    with pytest.raises(ValueError):
        poiseuille_reference(1.0, 1.0, 0.0)


def test_gate_order_recovers_first_order():
    # errors halving with the cell size is order 1
    result = gate_order([4.0e-3, 2.0e-3, 1.0e-3], [2.0e-2, 1.0e-2, 5.0e-3])
    assert result["orders"] == pytest.approx([1.0, 1.0])
    # a doubled error is detected as negative order, not silently as "converged"
    inverted = gate_order([1.0e-3, 2.0e-3], [1.0e-2, 5.0e-3])
    assert inverted["orders"][0] == pytest.approx(-1.0)


@pytest.mark.parametrize("radial_cells", [5, 10])
def test_staircase_area_is_close_and_stable(radial_cells):
    built = duct_lattice(radius=RADIUS, length=0.02, radial_cells=radial_cells)
    analytic = math.pi * RADIUS**2
    rel = built["cross_section_area_m2"] / analytic - 1.0
    assert abs(rel) < 0.08, rel
    # every periodic plane carries the same cross-section: the duct is straight
    nx, ny, nz = built["shape"]
    per_plane = [0] * nz
    for _, _, k in built["fluid_cells"]:
        per_plane[k] += 1
    assert len(set(per_plane)) == 1
    # and the fluid reaches both ends of the periodic axis, or the duct is blocked
    assert all(per_plane[k] > 0 for k in range(nz))


def test_staircase_area_converges_with_resolution():
    coarse = duct_lattice(radius=RADIUS, length=0.02, radial_cells=5)
    fine = duct_lattice(radius=RADIUS, length=0.02, radial_cells=20)
    analytic = math.pi * RADIUS**2
    assert abs(fine["cross_section_area_m2"] / analytic - 1.0) < abs(
        coarse["cross_section_area_m2"] / analytic - 1.0
    )


def test_bundle_is_the_contract_bundle(tmp_path):
    bundle = tmp_path / "duct"
    duct_bundle(bundle, radius=RADIUS, radial_cells=6, axial_cells=2)
    manifest = verify_bundle(bundle)
    assert manifest["route"] == "lattice"
    assert manifest["frame"] == COORDINATE_CONTRACT["name"]
    table = read_lattice_table(bundle / "lattice.json")
    nx, ny, nz = table["shape"]
    assert nx == ny > nz  # the box is the cross-section; the axis is short and periodic
    assert table["axis"] == {"i": "+x", "j": "+y", "k": "+z"}
    assert table["walls"]["axial_axis"] == "z"
    handoff = json.loads((bundle / "handoff.json").read_text(encoding="utf-8"))
    assert handoff["tables"]["lattice_json"] == "lattice.json"
    assert "cross_section" in json.dumps(handoff)  # the writer records the voxel area


def test_duct_flow_matches_poiseuille_and_is_axial(tmp_path):
    bundle, result = _small_case(tmp_path, radial_cells=5, axial_cells=2)
    metrics = result["metrics"]
    assert metrics["converged"]
    assert metrics["steps"] < 6000
    # no motion across the duct: a body force along +z must not create a swirl
    assert metrics["transverse_speed_max"] < 1e-6
    # mass is conserved by the plain source to round-off (sum of S_i = 0)
    assert metrics["mass_drift"] == pytest.approx(0.0, abs=1e-9)
    # magnitude against the closed form at the area-equivalent radius
    ratio = metrics["flow_rate_lattice"] / metrics["reference_flow_rate_lattice"]
    assert 0.80 < ratio < 1.20, ratio
    ratio_max = metrics["u_max_lattice"] / metrics["reference_u_max_lattice"]
    assert 0.85 < ratio_max < 1.15, ratio_max
    # mean is half the maximum for a developed duct: the profile is not a plug
    assert metrics["u_mean_lattice"] / metrics["u_max_lattice"] == pytest.approx(0.5, abs=0.12)


def test_flow_rate_scales_linearly_with_body_force(tmp_path):
    """Stokes regime: doubling the force must double the flux, not 1.9x or 2.1x."""
    bundle = tmp_path / "duct"
    duct_bundle(bundle, radius=RADIUS, radial_cells=5, axial_cells=2)
    fast = dict(tau=TAU, max_steps=6000, tol=1e-13)
    one = lattice_flow_from_bundle(bundle, body_force=BODY_FORCE, **fast)["metrics"]
    two = lattice_flow_from_bundle(bundle, body_force=2.0 * BODY_FORCE, **fast)["metrics"]
    assert two["flow_rate_lattice"] / one["flow_rate_lattice"] == pytest.approx(2.0, rel=0.01)


def test_wall_mask_wraps_on_the_periodic_axis():
    import numpy as np

    fluid = np.ones((4, 3, 3), dtype=bool)
    fluid[0, 1, 1] = False  # a single wall cell on the first periodic plane
    seen = _shifted_wall_mask(~fluid, (1, 0, 0))  # wall at k+1
    assert seen[3, 1, 1]  # the last plane sees the wall across the periodic seam
    assert not seen[0, 1, 1]
    assert not seen[2, 1, 1]


def test_wrong_frame_is_refused(tmp_path):
    """A bundle whose declared frame is not the contract's is refused before solving."""
    bundle = tmp_path / "duct"
    duct_bundle(bundle, radius=RADIUS, radial_cells=5, axial_cells=2)
    handoff_path = bundle / "handoff.json"
    handoff = json.loads(handoff_path.read_text(encoding="utf-8"))
    handoff["frame"]["name"] = "y-up"
    handoff_path.write_text(json.dumps(handoff), encoding="utf-8")
    with pytest.raises(HandoffError):
        lattice_flow_from_bundle(bundle, max_steps=10)


def test_stream3d_probe_reports_a_measurement():
    probe = stream3d_probe()
    assert probe["status"] in ("measured", "unavailable")
    if probe["status"] == "measured":
        assert probe["expected_dx_lattice"] == 1
        assert probe["agrees_with_pull_convention"] == (probe["dx_lattice"] == 1)