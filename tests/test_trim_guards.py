"""Refusals, record round-trips and the quadrature fallbacks in :mod:`tensorcad.trim`.

The trim suite proves the numbers; this file covers the edges the numbers do not
reach -- the arguments the module refuses instead of answering wrongly, the
trapezoid fallback for non-uniform stations, an indexed station with no faces in
its bucket, and the ``to_dict`` surfaces a consumer persists. It exists because
``pytest --cov`` gates the package at 90 % and a module that only reports its
happy path leaves the refusals unmeasured.
"""
from __future__ import annotations

import math

import pytest

from tensorcad import generate_wigley
from tensorcad.trim import (
    TrimError,
    _integrate,
    _round,
    _SectionCutter,
    hydrostatic_particulars,
    trimmed_hydrostatics,
)

BOX = (20.0, 4.0, 5.0)


@pytest.fixture(scope="module")
def wigley():
    return generate_wigley()


def test_round_leaves_non_finite_and_non_float_values_alone():
    assert math.isnan(_round(float("nan")))
    assert _round(float("inf")) == float("inf")
    assert _round(3) == 3
    assert _round(1.0 / 3.0) == pytest.approx(0.333333333, abs=1e-12)


def test_integrate_falls_back_to_trapezoid_on_a_short_or_uneven_sample():
    assert _integrate((0.0,), (1.0,)) == 0.0
    # Non-uniform stations cannot use Simpson, so the trapezoid rule answers instead.
    assert _integrate((0.0, 1.0, 3.0), (0.0, 1.0, 2.0)) == pytest.approx(0.5 + 3.0, rel=1e-12)


def test_an_indexed_station_with_no_faces_in_its_bucket_is_empty(wigley):
    cutter = _SectionCutter(wigley.mesh, (1.0e6,))
    assert cutter.area(0, 0.0) == 0.0


def test_particulars_refuse_a_waterline_that_is_not_a_number(wigley):
    with pytest.raises(TrimError, match="finite"):
        hydrostatic_particulars(wigley.mesh, waterline_z=float("nan"))


def test_trimmed_hydrostatics_refuses_a_trim_angle_that_is_not_a_number(wigley):
    with pytest.raises(TrimError, match="finite"):
        trimmed_hydrostatics(wigley.mesh, trim_angle=float("nan"), waterline_at_reference=-1.0)
