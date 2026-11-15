"""Tests for the hull-form offsets added in round 4: parallel body and keel rise.

The round-3 family is a flat-keel sine planform. Round 4 generalises it with two real
hull-form features and must not move the round-3 results, so the first test here is a
bit-for-bit regression of the old family.
"""
from __future__ import annotations

import math

import pytest

from tensorcad import mesh as meshlib
from tensorcad.cb_sections import (
    KEEL_RISE_EXPONENT,
    CalibratedHullParameters,
    build_calibrated_hull,
    keel_ratio_profile,
    power_mean,
    realised_coefficients,
    section_area_ratio,
    sine_power_mean,
    taper_profile,
    weighted_power_mean,
)


def _params(**kwargs) -> CalibratedHullParameters:
    base = dict(
        length=200.0, beam=30.0, draft=10.0,
        block_coefficient=0.65, midship_coefficient=0.95,
        waterplane_coefficient=0.78, nx=81, nz=21,
    )
    base.update(kwargs)
    return CalibratedHullParameters(**base)


# --------------------------------------------------------------------------- #
# Taper profile
# --------------------------------------------------------------------------- #
def test_zero_parallel_body_reproduces_the_sine_law():
    taper = taper_profile(81, 0.0)
    last = 80
    for i, value in enumerate(taper):
        expected = 0.0 if i in (0, last) else math.sin(math.pi * i / last)
        assert value == pytest.approx(expected, abs=1e-15)
    # ... and therefore the round-3 helper is unchanged.
    assert power_mean(taper, 0.92) == pytest.approx(sine_power_mean(0.92, 81), abs=1e-15)


def test_end_stations_are_exactly_zero():
    """Non-zero ends leave the loft open and the volume short by a whole end cap."""
    for fraction in (0.0, 0.15, 0.3, 0.45):
        taper = taper_profile(81, fraction)
        assert taper[0] == 0.0 and taper[-1] == 0.0
        assert 0.0 <= min(taper) and max(taper) <= 1.0


def test_parallel_body_is_flat_topped():
    taper = taper_profile(81, 0.3)
    middle = [value for i, value in enumerate(taper) if abs(i / 80 - 0.5) <= 0.15]
    assert middle and all(value == 1.0 for value in middle)
    assert taper[1] < taper[20] < 1.0  # the entrance is finer than the body


def test_taper_profile_rejects_impossible_fractions():
    for fraction in (-0.1, 1.0, 1.5):
        with pytest.raises(ValueError):
            taper_profile(21, fraction)


# --------------------------------------------------------------------------- #
# Keel ratio profile
# --------------------------------------------------------------------------- #
def test_keel_ratio_is_one_where_the_taper_is_one():
    taper = taper_profile(81, 0.3)
    keel = keel_ratio_profile(taper, 0.4)
    for g, h in zip(taper, keel):
        assert h == pytest.approx(0.4 + 0.6 * g**KEEL_RISE_EXPONENT, rel=1e-12)
    assert max(keel) == 1.0


def test_flat_keel_ratio_is_one_everywhere():
    assert set(keel_ratio_profile(taper_profile(41, 0.2), 1.0)) == {1.0}


def test_keel_ratio_rejects_out_of_range():
    for ratio in (0.0, -0.2, 1.2):
        with pytest.raises(ValueError):
            keel_ratio_profile((0.0, 0.5, 1.0, 0.5, 0.0), ratio)


# --------------------------------------------------------------------------- #
# Round-3 regression: the flat-keel family must be untouched
# --------------------------------------------------------------------------- #
def test_flat_keel_family_still_realises_its_coefficients():
    params = _params()
    hull = build_calibrated_hull(params)
    measured = realised_coefficients(
        hull.mesh, length=params.length, beam=params.beam, draft=params.draft
    )
    for name, want in (("Cb", 0.65), ("Cm", 0.95), ("Cwp", 0.78)):
        assert measured["coefficients"][name] == pytest.approx(want, rel=1e-9)
    assert meshlib.is_watertight(hull.mesh)


def test_flat_keel_breadth_exponent_matches_the_round3_value():
    """The round-3 solve for Cwp = 0.78 on 81 stations is a fixed number."""
    hull = build_calibrated_hull(_params())
    assert hull.breadth_exponent == pytest.approx(0.4451, abs=2e-3)


# --------------------------------------------------------------------------- #
# Keel-rise family
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("end_draft_ratio", [0.9, 0.7, 0.5])
def test_keel_rise_still_realises_its_coefficients(end_draft_ratio):
    params = _params(block_coefficient=0.55, waterplane_coefficient=0.70,
                     end_draft_ratio=end_draft_ratio)
    if params.problems():
        pytest.skip("coefficient vector cannot be hollowed out that far")
    hull = build_calibrated_hull(params)
    measured = realised_coefficients(
        hull.mesh, length=params.length, beam=params.beam, draft=params.draft
    )
    # The flat-keel family is exact to round-off; a rising keel makes the sectional
    # area vary quadratically between stations, so the loft's trapezoid rule leaves
    # ~2e-5 -- still four orders inside the 3 % gate (documented in the report).
    assert measured["coefficients"]["Cb"] == pytest.approx(0.55, rel=1e-4)
    assert measured["coefficients"]["Cm"] == pytest.approx(0.95, rel=1e-8)
    assert measured["coefficients"]["Cwp"] == pytest.approx(0.70, rel=1e-8)
    assert meshlib.is_watertight(hull.mesh)


def test_keel_rise_reduces_wetted_surface_at_fixed_coefficients():
    """The whole point of the feature: same Cb/Cm/Cwp, less wetted surface."""
    flat = build_calibrated_hull(
        _params(block_coefficient=0.55, waterplane_coefficient=0.70)
    )
    risen = build_calibrated_hull(
        _params(block_coefficient=0.55, waterplane_coefficient=0.70, end_draft_ratio=0.8)
    )
    def _wetted(hull):
        return realised_coefficients(
            hull.mesh, length=200.0, beam=30.0, draft=10.0
        )["wetted_surface_m2"]

    assert _wetted(risen) < _wetted(flat)


def test_keel_rise_keeps_the_keel_at_the_requested_draft():
    params = _params(block_coefficient=0.55, waterplane_coefficient=0.70,
                     end_draft_ratio=0.8)
    hull = build_calibrated_hull(params)
    keel_z = min(vertex[2] for vertex in hull.mesh.vertices)
    # The deepest point is amidships at the full design draft, never deeper.
    assert keel_z == pytest.approx(-params.draft, rel=1e-9)


def test_full_hull_cannot_be_hollowed_out_indefinitely():
    """A full tanker can be hollowed out only so far: the ends still owe volume.

    The refusal has to be monotone in the cut-away, and it has to be *reported*: a
    hull that quietly misses Cb would poison every downstream subject.
    """
    def problems_at(ratio: float) -> list[str]:
        return _params(block_coefficient=0.86, midship_coefficient=0.99,
                       waterplane_coefficient=0.90, end_draft_ratio=ratio).problems()

    assert not problems_at(1.0)
    refused = [ratio for ratio in (0.6, 0.5, 0.4, 0.3, 0.2, 0.15) if problems_at(ratio)]
    assert refused, "a deep cut-away must eventually be refused for a full hull"
    assert "end_draft_ratio" in " ".join(problems_at(refused[0]))
    assert max(refused) < 1.0


def test_keel_ratio_out_of_range_is_reported_by_problems():
    assert _params(end_draft_ratio=0.0).problems()
    assert _params(end_draft_ratio=1.5).problems()
    assert _params(parallel_body_fraction=-0.1).problems()


# --------------------------------------------------------------------------- #
# Discrete functionals
# --------------------------------------------------------------------------- #
def test_weighted_power_mean_equals_unweighted_with_unit_weights():
    taper = taper_profile(41, 0.2)
    for beta in (0.3, 1.0, 2.5):
        assert weighted_power_mean(taper, beta, [1.0] * len(taper)) == pytest.approx(
            power_mean(taper, beta), rel=1e-15
        )


def test_weighted_power_mean_requires_matching_lengths():
    with pytest.raises(ValueError):
        weighted_power_mean((0.1, 0.2), 1.0, (1.0,))


def test_section_area_curve_integrates_to_the_mesh_volume():
    from tensorcad.cb_sections import sectional_area_curve

    params = _params(block_coefficient=0.55, waterplane_coefficient=0.70,
                     end_draft_ratio=0.8)
    hull = build_calibrated_hull(params)
    xs, areas = sectional_area_curve(hull)
    last = len(xs) - 1
    integral = (math.fsum(areas) - 0.5 * (areas[0] + areas[-1])) / last * params.length
    volume = meshlib.total_signed_volume(hull.mesh)
    assert integral == pytest.approx(volume, rel=1e-4)


def test_section_area_ratio_is_the_polygon_area_not_its_limit():
    """R(p) at the flat-bottom limit is slightly below 1 on a finite level grid."""
    assert section_area_ratio(1.0e-4, 21) < 1.0
    assert section_area_ratio(1.0e-4, 21) > 0.99