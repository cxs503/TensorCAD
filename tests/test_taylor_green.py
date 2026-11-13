"""Tests for :mod:`tensorcad.taylor_green` (Taylor-Green reference solutions).

Verification philosophy: "exact" must be a *measured* statement.  Every test
below either evaluates an analytic residual, compares two independent
derivations of the same quantity, or checks series consistency across orders.
"""
from __future__ import annotations

import json
import math

import pytest

from tensorcad.taylor_green import (
    BeltramiABC,
    TaylorGreen3DSeries,
    TaylorGreenMode,
    case_spec,
    compare_metrics,
    fit_decay_rate,
    metrics_table,
    mode_table,
    sample_grid,
    write_case_bundle,
)

POINTS = [(0.31, 1.17, 2.05), (2.7, 0.4, 5.1), (4.2, 5.6, 0.9), (1.0, 3.0, 6.0)]


# ---------------------------------------------------------------------------
# closed forms
# ---------------------------------------------------------------------------
def test_invariant_mode_is_an_exact_solution():
    mode = TaylorGreenMode(amplitude=0.7, viscosity=0.03)
    for x, y, z in POINTS:
        for t in (0.0, 0.5, 2.0):
            res = mode.ns_residual(x, y, z, t)
            assert max(abs(c) for c in res) < 1e-14


def test_invariant_mode_divergence_free():
    mode = TaylorGreenMode(amplitude=0.7, viscosity=0.03)
    for x, y, z in POINTS:
        g = mode.velocity_gradient(x, y, z, 0.4)
        assert abs(g[0][0] + g[1][1] + g[2][2]) < 1e-15


def test_invariant_mode_energy_identity():
    """``dE/dt = -2 nu Z`` for the closed form, checked analytically."""
    nu, a = 0.03, 0.7
    mode = TaylorGreenMode(amplitude=a, viscosity=nu)
    for t in (0.0, 0.3, 1.5):
        d_e = -4.0 * nu * mode.energy(t)
        assert abs(d_e + 2.0 * nu * mode.enstrophy(t)) < 1e-15


def test_beltrami_is_exact_and_z_dependent():
    abc = BeltramiABC(amplitude=0.4, wavenumber=1.0, viscosity=0.02)
    for x, y, z in POINTS:
        for t in (0.0, 0.7, 2.5):
            assert abc.curl_eigen_residual(x, y, z, t) < 1e-14
            assert max(abs(c) for c in abc.ns_residual(x, y, z, t)) < 1e-13
    # genuinely three-dimensional: the velocity depends on z
    assert abs(abc.velocity(0.0, 0.0, 0.0, 0.0)[0]
               - abc.velocity(0.0, 0.0, 0.9, 0.0)[0]) > 1e-3


def test_fit_decay_rate_recovers_analytic_rate():
    mode = TaylorGreenMode(amplitude=1.0, viscosity=0.05)
    times = [0.0, 0.5, 1.0, 1.5, 2.0]
    energies = [mode.energy(t) for t in times]
    fit = fit_decay_rate(times, energies)
    assert abs(fit["rate"] - mode.decay_rate()) < 1e-12
    assert fit["max_log_residual"] < 1e-12


# ---------------------------------------------------------------------------
# series reference
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def series10() -> TaylorGreen3DSeries:
    return TaylorGreen3DSeries(order=8, amplitude=0.25, viscosity=0.01)


def test_series_initial_condition(series10: TaylorGreen3DSeries):
    for x, y, z in POINTS:
        got = series10.velocity(x, y, z, 0.0)
        want = series10.initial_velocity(x, y, z)
        assert max(abs(got[i] - want[i]) for i in range(3)) < 1e-15


def test_series_divergence_free(series10: TaylorGreen3DSeries):
    for x, y, z in POINTS:
        for t in (0.0, 1.0, 3.0):
            assert abs(series10.divergence(x, y, z, t)) < 1e-13


def test_series_defect_spectrum_is_order_complete(series10: TaylorGreen3DSeries):
    spec = series10.defect_spectrum()
    assert len(spec) == series10.order + 1
    assert max(spec[:series10.order]) < 1e-9  # every order balances
    assert spec[series10.order] > 1e-12  # the first unfixed (truncation) term


def test_series_dissipation_identity(series10: TaylorGreen3DSeries):
    # dE/dt = -2 nu Z holds exactly for the NS solution; for the truncated
    # polynomial it holds up to truncation error, so the window is t <= 2.
    for t in (0.0, 1.0, 2.0):
        lhs = series10.dissipation_rate(t)
        rhs = -2.0 * series10.viscosity * series10.enstrophy(t)
        assert abs(lhs - rhs) < 1e-8 * max(1.0, abs(rhs))


def test_series_converges_with_order():
    coarse = TaylorGreen3DSeries(order=6, amplitude=0.25, viscosity=0.01)
    fine = TaylorGreen3DSeries(order=10, amplitude=0.25, viscosity=0.01)
    for t, tol in ((1.0, 1e-7), (2.0, 1e-6), (3.0, 2e-5), (4.0, 5e-4)):
        e_fine = fine.energy(t)
        rel = abs(coarse.energy(t) - e_fine) / e_fine
        assert rel < tol, (t, rel)


def test_series_energy_decays_monotonically():
    ref = TaylorGreen3DSeries(order=10, amplitude=0.25, viscosity=0.01)
    energies = [ref.energy(t) for t in (0.0, 1.0, 2.0, 3.0, 4.0)]
    assert all(b < a for a, b in zip(energies, energies[1:]))
    assert abs(energies[0] - 0.25 ** 2 / 8.0) < 1e-15  # A^2/8


# ---------------------------------------------------------------------------
# consumption helpers
# ---------------------------------------------------------------------------
def test_compare_metrics_flags_tolerance():
    ref = {"energy": 1.0, "enstrophy": 2.0}
    good = compare_metrics(ref, {"energy": 1.05, "enstrophy": 1.9})
    assert good["all_pass"] and abs(good["max_rel_error"] - 0.05) < 1e-12
    bad = compare_metrics(ref, {"energy": 1.2, "enstrophy": 2.0})
    assert not bad["all_pass"]
    assert bad["rows"][0]["pass"] is False
    assert compare_metrics(ref, {"energy": 1.09}, tolerance=0.1)["all_pass"]


def test_mode_table_and_metrics_table(series10: TaylorGreen3DSeries):
    modes = mode_table(series10, 2.0, max_modes=8)
    assert modes and all(len(m["mode"]) == 3 for m in modes)
    amps = [m["amplitude"] for m in modes]
    assert amps == sorted(amps, reverse=True)
    table = metrics_table(series10, (1.0, 2.0))
    assert table["case"] == "CASE-TAYLOR-GREEN-3D"
    for row in table["rows"]:
        assert abs(row["dissipation"] - row["two_nu_Z"]) < 1e-8 * max(
            1.0, abs(row["two_nu_Z"]))


def test_sample_grid_shape(series10: TaylorGreen3DSeries):
    grid = sample_grid(series10, 4, 1.0)
    assert grid["level"] == 4
    assert len(grid["u"]) == 4 and len(grid["u"][0][0]) == 4
    assert abs(grid["spacing"] - 2.0 * math.pi / 4) < 1e-15


def test_case_spec_reports_provenance(series10: TaylorGreen3DSeries):
    spec = case_spec(series10, (1.0, 2.0))
    assert spec["case"] == "CASE-TAYLOR-GREEN-3D"
    assert spec["relative_tolerance"] == pytest.approx(0.08)
    assert spec["parameters"]["amplitude_A"] == pytest.approx(0.25)
    # dE/dt at t=0 must equal -2 nu Z0 exactly
    par = spec["parameters"]
    assert par["initial_dissipation_dE_dt"] == pytest.approx(
        -2.0 * par["kinematic_viscosity_nu"] * par["enstrophy_Z0"], abs=1e-15)


def test_write_case_bundle_is_deterministic(tmp_path):
    series = TaylorGreen3DSeries(order=6, amplitude=0.25, viscosity=0.01)
    man1 = write_case_bundle(tmp_path / "a", series=series, times=(1.0, 2.0),
                             grid_levels=(4,))
    man2 = write_case_bundle(tmp_path / "b", series=series, times=(1.0, 2.0),
                             grid_levels=(4,))
    assert set(man1["files"]) == set(man2["files"])
    assert man1["files"] == man2["files"]  # sha256 of every artefact
    spec = json.loads((tmp_path / "a" / "case_spec.json").read_text())
    assert spec["case"] == "CASE-TAYLOR-GREEN-3D"
    for name in man1["files"]:
        assert (tmp_path / "a" / name).exists()