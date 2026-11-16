"""Tests: manufactured-solution verification of the hull CAD -> FVM route.

The references are checked before they are used to judge anything (a reference that
silently stops being harmonic would turn the whole study into a tautology), then the
cheapest geometry level of each family is actually solved through the handoff.
"""
from __future__ import annotations

import math

import pytest

from tensorcad import fvm_manufactured as fm
from tensorcad.fvm_handoff import barge_bundle, hull_volume_bundle

BOUNDS = ((0.0, -0.5, -0.6), (10.0, 0.5, 0.0))


# --------------------------------------------------------------------------- references


@pytest.mark.parametrize("reference", [fm.QUADRATIC_A, fm.QUADRATIC_B, fm.wave_reference(2.0)])
def test_each_reference_is_harmonic_under_numerical_differentiation(reference):
    probe = fm.harmonic_probe(reference, BOUNDS)
    assert probe["passed"], probe
    assert probe["harmonic_defect"] < probe["tolerance"]


def test_a_non_harmonic_field_is_rejected_by_the_same_probe():
    """Without this, ``harmonic_probe`` could be passing everything."""
    broken = fm.Reference(name="x_squared", field=lambda p: p[:, 0] ** 2, closed_form="phi = x^2")
    probe = fm.harmonic_probe(broken, BOUNDS)
    assert not probe["passed"]
    assert probe["harmonic_defect"] > 1.0


def test_wave_reference_requires_a_positive_wavenumber():
    with pytest.raises(ValueError):
        fm.wave_reference(0.0)


def test_wave_reference_closed_form_states_its_wavenumber():
    reference = fm.wave_reference(3.5, y0=0.5)
    assert reference.wavenumber == 3.5
    assert "3.5" in reference.closed_form
    assert reference.kh(0.02) == pytest.approx(0.07)


# --------------------------------------------------------------------------- skew census


def test_skew_census_agrees_with_the_solver_audit(tmp_path):
    bundle = tmp_path / "box"
    barge_bundle(bundle, length=2.0, beam=0.4, draft=0.3, nx=9, nz=5, layers=2)
    stats = fm.mesh_skew_stats(bundle)
    assert 0.0 <= stats["mean_nonorthogonality_degrees"] <= stats["p95_nonorthogonality_degrees"]
    assert stats["p95_nonorthogonality_degrees"] <= stats["maximum_nonorthogonality_degrees"] + 1e-9
    assert stats["internal_faces"] > 0
    # The census is measured on the file, so it must reproduce the solver's own
    # internal-face angle exactly.  The solver's audit additionally folds boundary
    # faces into its maximum (its ``d`` for a boundary face is centre -> face centre),
    # so its headline number can only be larger -- never equal by construction.
    import numpy as np

    from tensorcad.fvm_handoff import _open_bundle

    opened = _open_bundle(bundle)
    mesh = opened.mesh
    vectors = mesh.area_vectors
    delta = mesh.cell_centers[mesh.neighbor] - mesh.cell_centers[mesh.owner]
    inside = mesh.neighbor >= 0
    projection = np.einsum("ij,ij->i", vectors, delta)
    angle = np.degrees(
        np.arccos(
            np.clip(
                projection / (np.linalg.norm(vectors, axis=1) * np.linalg.norm(delta, axis=1)), -1.0, 1.0
            )
        )
    )
    internal_max = float(angle[inside].max())
    assert stats["maximum_nonorthogonality_degrees"] == pytest.approx(internal_max, abs=1e-9)
    audit = mesh.audit()
    assert float(audit["maximum_nonorthogonality_degrees"]) >= internal_max - 1e-9


# --------------------------------------------------------------------------- solves


def test_quadratic_reference_is_reproduced_on_the_curved_hull(tmp_path):
    """The hull leg previously carried only the linear patch test."""
    bundle = tmp_path / "hull"
    hull_volume_bundle(bundle, length=10.0, beam=1.0, draft=0.6, nx=13, nz=7, layers=2)
    result = fm.manufactured_from_bundle(bundle, fm.QUADRATIC_A)
    assert result.status == "ran", result.error
    assert not result.gates_failed, result.gates_failed
    metrics = result.metrics
    assert metrics["frame"] == "z-up"
    assert metrics["relative_linf"] < 3e-3
    assert metrics["relative_l2"] < 3e-3
    assert metrics["mean_nonorthogonality_degrees"] > 30.0  # this really is a skewed mesh
    assert metrics["declared_volume_m3"] > 0.0


def test_two_quadratics_are_independent_references(tmp_path):
    """A second, differently oriented quadratic -- not the same field rescaled."""
    bundle = tmp_path / "hull_b"
    hull_volume_bundle(bundle, length=10.0, beam=1.0, draft=0.6, nx=13, nz=7, layers=2)
    a = fm.manufactured_from_bundle(bundle, fm.QUADRATIC_A)
    b = fm.manufactured_from_bundle(bundle, fm.QUADRATIC_B)
    assert a.status == "ran" and b.status == "ran"
    assert a.metrics["relative_linf"] != pytest.approx(b.metrics["relative_linf"], rel=1e-6)


def test_wave_error_grows_with_curvature_resolution_on_one_mesh(tmp_path):
    """Same mesh, two frequencies: the error follows ``k h``, not the mesh alone."""
    bundle = tmp_path / "box_wave"
    barge_bundle(bundle, length=2.0, beam=0.4, draft=0.3, nx=9, nz=5, layers=2)
    low = fm.manufactured_from_bundle(bundle, fm.wave_reference(0.5, y0=0.2), tolerance=0.05)
    high = fm.manufactured_from_bundle(bundle, fm.wave_reference(4.0, y0=0.2), tolerance=0.50)
    assert low.status == "ran" and high.status == "ran"
    assert high.metrics["relative_linf"] > 5.0 * low.metrics["relative_linf"]
    assert low.metrics["k_h"] < high.metrics["k_h"]


# --------------------------------------------------------------------------- order bookkeeping


def test_observed_order_is_arithmetic_and_flags_the_asymptotic_range():
    coarse = {"cell_size_m": 0.2, "relative_linf": 4e-2, "k_h": 0.4}
    fine = {"cell_size_m": 0.1, "relative_linf": 1e-2, "k_h": 0.2}
    row = fm.observed_order(coarse, fine)
    assert row["order"] == pytest.approx(2.0)
    assert row["cell_size_ratio"] == pytest.approx(2.0)
    assert row["asymptotic"] is False  # k h = 0.4 / 0.2 is outside the band

    near = fm.observed_order(
        {"cell_size_m": 0.002, "relative_linf": 4e-6, "k_h": 0.004},
        {"cell_size_m": 0.001, "relative_linf": 1e-6, "k_h": 0.002},
    )
    assert near["asymptotic"] is True
    assert near["order"] == pytest.approx(2.0)


def test_convergence_rows_walks_consecutive_levels():
    rows = [
        {"cell_size_m": 0.4, "relative_linf": 16e-3, "k_h": 0.8},
        {"cell_size_m": 0.2, "relative_linf": 4e-3, "k_h": 0.4},
        {"cell_size_m": 0.1, "relative_linf": 1e-3, "k_h": 0.2},
    ]
    orders = [row["order"] for row in fm.convergence_rows(rows)]
    assert orders == pytest.approx([2.0, 2.0])


def test_preasymptotic_countercheck_states_the_curvature_resolution():
    check = fm.preasymptotic_countercheck(2.0, 0.4, 0.3)
    gamma = math.pi * math.sqrt(1.0 / 2.0**2 + 1.0 / 0.4**2)
    assert check["gamma_1_per_m"] == pytest.approx(gamma)
    assert gamma * 0.0679 > 0.15  # the earlier round's coarsest box level
    assert check["asymptotic_threshold_k_h"] == 0.15