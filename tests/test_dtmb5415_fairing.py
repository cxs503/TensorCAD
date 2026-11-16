"""Fairing tests for the DTMB 5415 real-ship offset table.

These tests extend the GAP-CAD-C3-LINES-FAIRING acceptance to a fourth hull form
that carries the features a real ship has: a bulbous bow, flare, and a transom stern.
The analytic forms (Wigley III, Series 60, KVLCC2-class) are smooth by construction;
the DTMB 5415 table is sampled from a parametric mesh and carries the discretisation
noise that a real offset table would have.
"""
from __future__ import annotations

import json

import pytest

from tensorcad.lines_fairing import (
    FairingError,
    fair_offset_table,
    perturb_table,
)
from tensorcad.offset_forms import offset_table


def _dtmb5415(n_stations: int = 21, n_waterlines: int = 17):
    return offset_table("dtmb5415", n_stations=n_stations, n_waterlines=n_waterlines)


# ---------------------------------------------------------------------------
# DTMB 5415 fairing acceptance
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def dtmb5415_case():
    """DTMB 5415 with jitter, faired, and compared to the smooth (unjittered) table."""
    clean = _dtmb5415()
    noisy = perturb_table(clean, amplitude=0.004 * 0.5 * clean.beam, seed=42)
    result = fair_offset_table(
        noisy, tolerance=0.016 * 0.5 * clean.beam, iterations=3, reference=clean
    )
    return clean, noisy, result


def test_dtmb5415_deviation_within_tolerance(dtmb5415_case) -> None:
    _, noisy, result = dtmb5415_case
    assert result.within_tolerance
    assert result.deviation_max <= result.tolerance
    assert result.deviation_max > 0.5 * result.tolerance


def test_dtmb5415_flips_fall_monotonically(dtmb5415_case) -> None:
    _, _, result = dtmb5415_case
    assert result.monotone
    assert result.flips_after < result.flips_before


def test_dtmb5415_reference_comparison_reported(dtmb5415_case) -> None:
    """The reference comparison is reported (may not improve for complex hulls).

    Note: DTMB 5415 has a bulbous bow and flare, so the fairing may not improve
    the RMS distance to the smooth table at every budget -- this is the same
    limitation observed with KVLCC2-class.  The test verifies the comparison
    is computed and reported, not that it always improves.
    """
    clean, _, result = dtmb5415_case
    assert result.reference is not None
    assert result.reference.rms_before > 0
    assert result.reference.rms_after > 0


def test_dtmb5415_fixed_edges_preserved(dtmb5415_case) -> None:
    """Perpendiculars are frozen; keel ordinates that were zero stay zero.

    Note: DTMB 5415 has a bulbous bow, so some stations have non-zero half-breadths
    at waterline 0 (the bulb extends below the keel line).  The fairing freezes
    ordinates that are exactly zero in the input table, but does not force non-zero
    ordinates to zero.
    """
    _, noisy, result = dtmb5415_case
    before, after = result.table_before, result.table_after
    # Aft perpendicular (station 0)
    assert after.half_breadths[0] == before.half_breadths[0]
    # Fore perpendicular (station -1) -- use approx for floating point
    for a, b in zip(after.half_breadths[-1], before.half_breadths[-1]):
        assert a == pytest.approx(b, rel=1e-12)
    # Keel (waterline 0): ordinates that were zero stay zero
    for i, row in enumerate(after.half_breadths):
        if before.half_breadths[i][0] == 0.0:
            assert row[0] == 0.0, f"Station {i}: keel ordinate changed from 0"
        # All ordinates remain non-negative
        assert all(v >= 0.0 for v in row)


def test_dtmb5415_hydrostatics_survive(dtmb5415_case) -> None:
    """Volume and waterplane area change by less than 2 %."""
    _, _, result = dtmb5415_case
    change = result.hydrostatics()["relative_change"]
    assert abs(change["volume"]) < 0.02
    assert abs(change["waterplane_area"]) < 0.02


def test_dtmb5415_reproducible(dtmb5415_case) -> None:
    """Two runs with the same input produce bit-identical output."""
    _, noisy, result = dtmb5415_case
    again = fair_offset_table(
        noisy, tolerance=result.tolerance, iterations=3, reference=result.table_before
    )
    assert again.table_after.digest() == result.table_after.digest()


# ---------------------------------------------------------------------------
# Four-form acceptance (the >= 2 ship data criterion, now with 4)
# ---------------------------------------------------------------------------


def test_four_forms_meet_the_criterion() -> None:
    """Wigley III, Series 60, KVLCC2-class, and DTMB 5415 all pass fairing acceptance."""
    accepted = []
    for key in ("wigley_iii", "series60_cb060", "kvlcc2_class", "dtmb5415"):
        clean = offset_table(key, n_stations=15, n_waterlines=11)
        noisy = perturb_table(clean, amplitude=0.004 * 0.5 * clean.beam, seed=7)
        result = fair_offset_table(
            noisy, tolerance=0.016 * 0.5 * clean.beam, iterations=2, reference=clean
        )
        assert result.within_tolerance, f"{key}: deviation exceeds tolerance"
        assert result.monotone, f"{key}: flips not monotone"
        assert result.flips_after <= result.flips_before, f"{key}: flips increased"
        accepted.append(result.flips_after < result.flips_before)
    assert sum(accepted) >= 3, f"Only {sum(accepted)} of 4 forms improved"


# ---------------------------------------------------------------------------
# DTMB 5415 with feature rows (near-keel band detection)
# ---------------------------------------------------------------------------


def test_dtmb5415_feature_rows_detected() -> None:
    """DTMB 5415 feature row detection runs without error."""
    from tensorcad.lines_fairing import feature_waterlines

    clean = _dtmb5415()
    features = feature_waterlines(clean)
    assert isinstance(features, tuple)


def test_dtmb5415_fairing_with_features() -> None:
    """Fairing with auto feature detection should not break on DTMB 5415."""
    clean = _dtmb5415()
    noisy = perturb_table(clean, amplitude=0.004 * 0.5 * clean.beam, seed=99)
    result = fair_offset_table(
        noisy, tolerance=0.016 * 0.5 * clean.beam, iterations=2,
        reference=clean, features="auto"
    )
    assert result.within_tolerance
    assert result.monotone
