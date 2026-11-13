"""Verification of the Rayleigh-Benard onset reference bundle (:mod:`tensorcad.rb_onset`).

The tests are ordered the way the evidence chain is built:

1. the numerical machinery (Chebyshev operators, interpolation) is exercised on
   problems with known answers, so the discretisation cannot silently be wrong;
2. the marginal-state solver is validated against the **exact** stress-free
   closed form ``Ra = (pi^2 + a^2)^3 / a^2`` at three wavenumbers -- an
   independent, analytic benchmark of the whole scan (Chebyshev order, Rayleigh
   bisection and the ``min_a`` step);
3. only then is the rigid-rigid reference (1707.762, Chandrasekhar 1961)
   asserted, and the case-pool verdict checked;
4. the solver-side harness (growth rate, Ra bisection) is self-tested against an
   injected probe with a known onset, and shown to reject a non-bracketing probe
   instead of fabricating a number;
5. the cell specification and bundle digest are pinned, so a consumer (TensorLBM /
   TensorFVM / TensorDG) reads the documented geometry.
"""
from __future__ import annotations

import hashlib
import json
import math

import pytest

from tensorcad.rb_onset import (
    CASE_SPEC,
    CHANDRASEKHAR_1961,
    FREE_FREE_ANALYTIC,
    NONDIMENSIONALISATION,
    RB_ONSET_ACCEPTANCE,
    bisect_onset,
    chebyshev_differentiation_matrix,
    chebyshev_interpolate,
    chebyshev_lobatto,
    critical_rayleigh,
    determinant_sign_log,
    growth_rate,
    linear_stability_matrix,
    neutral_curve,
    nusselt_number,
    onset,
    onset_mode,
    rb_cell,
    write_bundle,
)

EXACT_FREE_FREE = 27.0 * math.pi**4 / 4.0  # 657.5113644795163
RIGID_REFERENCE = 1707.762


# --------------------------------------------------------------------------- #
# (1) numerical machinery on known answers
# --------------------------------------------------------------------------- #


def test_chebyshev_lobatto_spans_the_layer():
    nodes = chebyshev_lobatto(8)
    assert len(nodes) == 9
    assert nodes[0] == pytest.approx(0.0, abs=1e-15)
    assert nodes[-1] == pytest.approx(1.0, abs=1e-15)
    assert nodes == sorted(nodes), "nodes must run from the hot plate to the cold plate"
    with pytest.raises(ValueError):
        chebyshev_lobatto(0)


def test_differentiation_matrix_is_exact_for_polynomials():
    order = 10
    d = chebyshev_differentiation_matrix(order)
    nodes = chebyshev_lobatto(order)
    for degree in range(1, order):
        values = [z**degree for z in nodes]
        derivative = [sum(d[i][j] * values[j] for j in range(order + 1)) for i in range(order + 1)]
        exact = [degree * z ** (degree - 1) for z in nodes]
        assert max(abs(a - b) for a, b in zip(derivative, exact, strict=True)) < 1e-9 * max(1.0, degree)


def test_barycentric_interpolation_is_exact_for_polynomials():
    nodes = chebyshev_lobatto(8)
    values = [3.0 * z**3 - 2.0 * z + 0.5 for z in nodes]
    for z in (0.0, 0.13, 0.5, 0.77, 1.0):
        assert chebyshev_interpolate(nodes, values, z) == pytest.approx(3.0 * z**3 - 2.0 * z + 0.5, abs=1e-12)
    with pytest.raises(ValueError):
        chebyshev_interpolate(nodes, values[:-1], 0.5)


def test_determinant_sign_and_log_agree_with_direct_determinant():
    # det [[2,1],[1,3]] = 5 ; det [[1,2],[2,4]] = 0
    sign, log_abs = determinant_sign_log([[2.0, 1.0], [1.0, 3.0]])
    assert sign == 1
    assert math.exp(log_abs) == pytest.approx(5.0)
    sign, log_abs = determinant_sign_log([[1.0, 2.0], [2.0, 4.0]])
    assert sign == 0 and log_abs == -math.inf
    sign, _ = determinant_sign_log([[-1.0, 0.0], [0.0, 2.0]])
    assert sign == -1


def test_stability_matrix_is_square_and_rejects_bad_input():
    matrix = linear_stability_matrix(3.0, 1700.0, order=10, walls="rigid")
    assert len(matrix) == 22
    assert all(len(row) == 22 for row in matrix)
    with pytest.raises(ValueError):
        linear_stability_matrix(0.0, 1700.0, order=10)
    with pytest.raises(ValueError):
        linear_stability_matrix(3.0, 1700.0, order=4)
    with pytest.raises(ValueError):
        linear_stability_matrix(3.0, 1700.0, order=10, walls="sticky")


# --------------------------------------------------------------------------- #
# (2) exact free-free benchmark of the whole scan
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("wavenumber", [1.5, 2.221441469079183, 3.0])
def test_free_free_neutral_curve_matches_closed_form(wavenumber):
    """Stress-free plates admit the exact marginal curve ``(pi^2 + a^2)^3 / a^2``."""
    exact = (math.pi**2 + wavenumber**2) ** 3 / wavenumber**2
    value = critical_rayleigh(wavenumber, order=16, walls="free", tol=1e-7)
    assert value == pytest.approx(exact, rel=1e-6)


def test_free_free_minimum_is_the_sine_mode_value():
    scan = onset(
        order=12,
        walls="free",
        wavenumber_range=(1.5, 3.5),
        wavenumber_points=7,
        tol=1e-5,
    )
    assert scan.rayleigh_critical == pytest.approx(EXACT_FREE_FREE, rel=1e-4)
    assert scan.wavenumber_critical == pytest.approx(math.pi / math.sqrt(2.0), rel=2e-3)
    assert scan.verdict == "PASS"
    assert scan.reference == pytest.approx(EXACT_FREE_FREE)
    # the scan must sit above the analytic minimum everywhere
    assert all(ra >= EXACT_FREE_FREE * (1.0 - 1e-6) for _, ra in scan.wavenumber_scan)


# --------------------------------------------------------------------------- #
# (3) rigid-rigid reference and the case-pool verdict
# --------------------------------------------------------------------------- #


def test_rigid_rigid_onset_matches_chandrasekhar():
    scan = onset(
        order=16,
        walls="rigid",
        wavenumber_range=(2.0, 4.5),
        wavenumber_points=7,
        tol=1e-5,
    )
    assert scan.reference == RIGID_REFERENCE
    assert scan.rayleigh_critical == pytest.approx(RIGID_REFERENCE, rel=1e-4)
    assert scan.wavenumber_critical == pytest.approx(CHANDRASEKHAR_1961["wavenumber_critical"], rel=5e-3)
    assert scan.verdict == "PASS"
    assert scan.accepted
    assert abs(scan.relative_error) <= RB_ONSET_ACCEPTANCE["tolerance_relative"]
    payload = scan.as_dict()
    assert payload["walls"] == "rigid"
    assert payload["reference_source"].startswith("Chandrasekhar")


def test_stability_criterion_flips_across_the_marginal_state():
    """Below Ra_c only the trivial stationary perturbation exists, above it does not."""
    a_c = CHANDRASEKHAR_1961["wavenumber_critical"]
    order = 16
    ra_c = critical_rayleigh(a_c, order=order)
    below = determinant_sign_log(linear_stability_matrix(a_c, 0.9 * ra_c, order))[0]
    above = determinant_sign_log(linear_stability_matrix(a_c, 1.1 * ra_c, order))[0]
    assert below != 0 and above != 0
    assert below != above, "det M must change sign when the conduction state destabilises"


def test_critical_rayleigh_reports_no_onset_when_scan_is_too_short():
    assert critical_rayleigh(3.0, order=12, walls="rigid", rayleigh_max=100.0) is None


def test_neutral_curve_drops_unreachable_wavenumbers():
    curve = neutral_curve([3.0, 3.1], order=12, rayleigh_max=100.0)
    assert curve == []


# --------------------------------------------------------------------------- #
# (4) solver-side harness
# --------------------------------------------------------------------------- #


def test_growth_rate_recovers_a_known_exponential():
    sigma, amplitude = 0.37, 2.5
    times = [0.1 * i for i in range(40)]
    history = [amplitude * math.exp(sigma * t) for t in times]
    assert growth_rate(times, history) == pytest.approx(sigma, rel=1e-9)
    # discarding the startup window must not bias an exactly exponential history
    assert growth_rate(times, history, discard=5) == pytest.approx(sigma, rel=1e-9)
    with pytest.raises(ValueError):
        growth_rate(times, history[:-1])
    with pytest.raises(ValueError):
        growth_rate([1.0, 1.0], [1.0, 2.0])


def test_bisect_onset_recovers_an_injected_onset():
    """The harness must find a known onset from a probe it has never seen before."""
    planted = 1550.0

    def probe(ra: float) -> float:
        return (ra - planted) / 250.0 + 0.37 * (ra - planted) ** 3 / 250.0**3

    verdict = bisect_onset(
        probe,
        rayleigh_stable=800.0,
        rayleigh_unstable=3000.0,
        tol=1e-6,
        reference=RIGID_REFERENCE,
    )
    assert verdict.rayleigh_critical == pytest.approx(planted, rel=1e-5)
    assert verdict.relative_error is not None and not verdict.accepted  # 1550 is 9% off
    assert verdict.stable_side < verdict.rayleigh_critical < verdict.unstable_side
    assert verdict.history[0][1] < 0.0 < verdict.history[1][1]
    assert verdict.as_dict()["iterations"] == verdict.iterations


def test_bisect_onset_accepts_a_probe_that_matches_the_reference():
    def probe(ra: float) -> float:
        return (ra - RIGID_REFERENCE) / 400.0

    verdict = bisect_onset(
        probe, rayleigh_stable=1000.0, rayleigh_unstable=2400.0, tol=1e-6, reference=RIGID_REFERENCE
    )
    assert verdict.accepted is True


def test_bisect_onset_rejects_a_probe_without_a_sign_change():
    def always_stable(ra: float) -> float:
        return -1.0

    with pytest.raises(ValueError, match="does not bracket"):
        bisect_onset(always_stable, rayleigh_stable=1000.0, rayleigh_unstable=2000.0)
    with pytest.raises(ValueError):
        bisect_onset(always_stable, rayleigh_stable=2000.0, rayleigh_unstable=1000.0)


def test_nusselt_number_is_one_for_conduction_and_above_one_for_convection():
    nz = 20
    dz = 1.0 / nz
    conduction = [1.0 - k * dz for k in range(nz + 1)]  # nodal, starts on the hot wall
    assert nusselt_number(conduction, dz) == pytest.approx(1.0, abs=1e-12)
    # a hotter wall implies a steeper near-wall gradient: convective enhancement
    steeper = [conduction[0] + 0.01, *conduction[1:]]
    assert nusselt_number(steeper, dz) == pytest.approx(1.3, rel=1e-9)
    colder_wall = [conduction[0] - 0.01, *conduction[1:]]
    assert nusselt_number(colder_wall, dz) < 1.0
    with pytest.raises(ValueError):
        nusselt_number([1.0, 0.0], 0.5)
    with pytest.raises(ValueError):
        nusselt_number([1.0, 0.5, 0.0], 0.0)


# --------------------------------------------------------------------------- #
# (5) critical mode shape and cell specification
# --------------------------------------------------------------------------- #


def test_critical_mode_satisfies_the_boundary_conditions():
    a_c = CHANDRASEKHAR_1961["wavenumber_critical"]
    mode = onset_mode(a_c, RIGID_REFERENCE, order=16, iterations=30)
    assert mode.residual < 1e-10, f"mode is not a marginal solution: {mode.residual}"
    # Psi = D Psi = 0 on both rigid plates
    assert mode.psi[0] == pytest.approx(0.0, abs=1e-9)
    assert mode.psi[-1] == pytest.approx(0.0, abs=1e-9)
    assert abs(mode.psi[1]) < 0.05 * max(abs(v) for v in mode.psi)
    # temperature perturbation is normalised to unit peak and vanishes on the plates
    assert max(abs(v) for v in mode.theta) == pytest.approx(1.0, abs=1e-9)
    assert mode.theta[0] == pytest.approx(0.0, abs=1e-9)
    assert mode.theta[-1] == pytest.approx(0.0, abs=1e-9)
    # the roll core sits inside the layer
    peak_z = mode.nodes[max(range(len(mode.psi)), key=lambda i: abs(mode.psi[i]))]
    assert 0.2 < peak_z < 0.8
    # the 2-D reconstruction is courtesy of the same profile
    assert mode.stream_function(0.0, 0.5) == pytest.approx(
        chebyshev_interpolate(mode.nodes, mode.psi, 0.5), abs=1e-12
    )
    assert mode.temperature_perturbation(0.0, 0.5) == pytest.approx(
        chebyshev_interpolate(mode.nodes, mode.theta, 0.5), abs=1e-12
    )
    assert json.loads(json.dumps(mode.as_dict()))["walls"] == "rigid"


def test_rb_cell_specification_is_consumable():
    cell = rb_cell()
    assert cell.walls == "rigid"
    assert cell.width == pytest.approx(2.0 * math.pi / CHANDRASEKHAR_1961["wavenumber_critical"])
    assert cell.roll_wavelength == pytest.approx(cell.width)
    assert cell.height == 1.0
    assert cell.boundary["bottom"]["thermal"].startswith("T = 1")
    assert cell.boundary["top"]["thermal"].startswith("T = 0")
    assert cell.boundary["bottom"]["kinematic"] == "no-slip"
    assert cell.boundary["left"]["type"] == "periodic"
    assert cell.boundary["left"]["pairs_with"] == "right"
    assert cell.monitor["onset_reference"] == RIGID_REFERENCE
    assert cell.monitor["conduction_limit"] == 1.0
    assert cell.base_state(0.0) == pytest.approx(1.0)
    assert cell.base_state(0.75) == pytest.approx(0.25)
    xs, zs = cell.cell_centres(8, 4)
    assert len(xs) == 8 and len(zs) == 4
    assert xs[0] > 0.0 and xs[-1] < cell.width
    assert zs[0] == pytest.approx(0.125)
    assert [list(r) for r in cell.as_dict()["resolutions"]] == [[64, 32], [128, 64], [256, 128]]
    assert cell.as_dict()["nondimensionalisation"] == NONDIMENSIONALISATION
    with pytest.raises(ValueError):
        cell.cell_centres(0, 4)


def test_rb_cell_free_slip_variant_uses_the_exact_reference():
    cell = rb_cell(walls="free")
    assert cell.boundary["bottom"]["kinematic"] == "free-slip"
    assert cell.monitor["onset_reference"] == pytest.approx(EXACT_FREE_FREE)
    with pytest.raises(ValueError):
        rb_cell(wavenumber=0.0)
    with pytest.raises(ValueError):
        rb_cell(walls="slippery")


def test_case_spec_matches_the_ledger_entry():
    assert CASE_SPEC["case"] == "LBM-RB-ONSET-RA1708"
    assert CASE_SPEC["reference"] == RIGID_REFERENCE
    assert CASE_SPEC["tolerance_relative"] == RB_ONSET_ACCEPTANCE["tolerance_relative"] == 0.02
    assert CASE_SPEC["analytic_sanity_check"]["rayleigh_critical"] == pytest.approx(EXACT_FREE_FREE)
    assert FREE_FREE_ANALYTIC["rayleigh_critical"] == pytest.approx(EXACT_FREE_FREE)
    assert CASE_SPEC["method"].startswith("Ra bisection scan")


def test_write_bundle_is_digest_protected(tmp_path):
    scan = onset(order=12, walls="free", wavenumber_range=(1.5, 3.5), wavenumber_points=7, tol=1e-5)
    cell = rb_cell(walls="free")
    target = write_bundle(
        tmp_path / "rb_onset.json", scan=scan, cell=cell, extra={"note": "unit test"}
    )
    payload = json.loads(target.read_text())
    digest = payload.pop("sha256")
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    assert digest == hashlib.sha256(canonical.encode("utf-8")).hexdigest()
    assert payload["onset"]["verdict"] == "PASS"
    assert payload["cell"]["walls"] == "free"
    assert payload["note"] == "unit test"
