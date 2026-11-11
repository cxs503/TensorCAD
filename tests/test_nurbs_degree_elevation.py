"""NURBS curve degree elevation (absorbed from ISO 10303-42 / Piegl & Tiller A5.9).

Degree elevation is a standard CAD-kernel operation (OCCT ``Geom_BSplineCurve::
IncreaseDegree``, the Bezier form in every textbook): it rewrites a curve of degree
``p`` as one of degree ``q > p`` *without moving a single point of the curve*.
TensorCAD already offers the inverse operation family (knot insertion for surfaces),
but a curve could not be raised to a common degree, which is a prerequisite for
stitching patches of different degree into one assembly.

The tests below are written against the operation's contract, not its implementation:

* geometry and parametrisation are preserved to round-off (sampled at knots and in
  between, rational circles stay exactly on the circle);
* the knot vector grows by the textbook rule (every non-zero knot multiplicity
  ``m`` becomes ``m + t``; a clamped end ``p + 1`` becomes ``q + 1``);
* a degree-1 Bezier elevates to the analytic midpoint control point;
* the affine/Bezier identity holds for a single-span polynomial curve against
  Bernstein-basis evaluation written independently in the test;
* invalid ``times`` values are rejected instead of silently returning something.
"""
from __future__ import annotations

import math
from itertools import pairwise

import pytest

from tensorcad.nurbs import NURBSCurve, NURBSError, circle_curve

RADIUS = 2.0


def _sample_parameters(curve: NURBSCurve, per_span: int = 17) -> list[float]:
    """Every distinct knot, plus ``per_span`` interior samples of every span."""
    low, high = curve.domain
    distinct = sorted({value for value in curve.knots if low <= value <= high})
    parameters = list(distinct)
    for left, right in pairwise(distinct):
        for index in range(1, per_span):
            parameters.append(left + (right - left) * index / per_span)
    return parameters


def _max_deviation(first: NURBSCurve, second: NURBSCurve, per_span: int = 17) -> float:
    return max(
        math.dist(first.evaluate(u), second.evaluate(u))
        for u in _sample_parameters(first, per_span)
    )


def _bernstein(degree: int, index: int, value: float) -> float:
    return math.comb(degree, index) * value**index * (1.0 - value) ** (degree - index)


def _bezier_evaluate(points: list[tuple[float, float, float]], value: float) -> tuple[float, float, float]:
    degree = len(points) - 1
    return tuple(
        sum(_bernstein(degree, i, value) * point[axis] for i, point in enumerate(points))
        for axis in range(3)
    )


def test_degree_elevation_preserves_the_circle_exactly() -> None:
    """A rational circle raised to degree 3 and 5 is still the same circle."""
    original = circle_curve(RADIUS, arcs=3)
    for times in (1, 2, 3):
        elevated = original.elevate_degree(times)
        assert elevated.degree == original.degree + times
        assert elevated.rational
        assert all(weight > 0.0 for weight in elevated.weights)
        deviation = max(
            abs(math.hypot(*elevated.evaluate(u)[:2]) - RADIUS)
            for u in _sample_parameters(elevated)
        )
        assert deviation < 1e-12, (times, deviation)
        assert _max_deviation(original, elevated) < 1e-12, times
        # Parametrisation, not just the point set: the domain and the sample
        # correspondence are unchanged as well.
        assert elevated.domain == original.domain
        assert abs(elevated.arc_length() - original.arc_length()) < 1e-12


def test_elevation_follows_the_knot_multiplicity_rule() -> None:
    """Every knot multiplicity ``m`` becomes ``m + t``; clamped ends stay clamped."""
    original = circle_curve(RADIUS, arcs=3)
    degree, times = original.degree, 2
    elevated = original.elevate_degree(times)

    def multiplicities(knots: tuple[float, ...]) -> dict[float, int]:
        counts: dict[float, int] = {}
        for value in knots:
            counts[value] = counts.get(value, 0) + 1
        return counts

    before = multiplicities(original.knots)
    after = multiplicities(elevated.knots)
    assert set(before) == set(after)
    for value, count in before.items():
        assert after[value] == count + times, value
    assert elevated.knots[0] == elevated.knots[degree + times] == original.domain[0]
    assert elevated.knots[-1] == elevated.knots[-degree - times - 1] == original.domain[1]
    # Three spans sharing two joins: span_count * q + 1 control points.
    assert len(elevated.control_points) == 3 * (degree + times) + 1
    assert len(elevated.knots) == len(elevated.control_points) + elevated.degree + 1


def test_single_span_polynomial_matches_independent_bezier_evaluation() -> None:
    """A degree-3 Bezier elevated to 6 equals Bernstein evaluation of its control net."""
    points = [(0.0, 0.0, 0.0), (1.0, 2.0, 0.0), (2.0, 0.0, 0.0), (3.0, 3.0, 0.0)]
    curve = NURBSCurve(
        degree=3,
        knots=(0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0),
        control_points=tuple(points),
    )
    elevated = curve.elevate_degree(3)
    assert elevated.degree == 6
    assert len(elevated.control_points) == 7
    for index in range(101):
        u = index / 100.0
        assert math.dist(elevated.evaluate(u), _bezier_evaluate(points, u)) < 1e-12, u


def test_degree_one_line_elevates_to_the_analytic_midpoint() -> None:
    line = NURBSCurve(
        degree=1,
        knots=(0.0, 0.0, 1.0, 1.0),
        control_points=((0.0, 0.0, 0.0), (2.0, 4.0, -6.0)),
    )
    elevated = line.elevate_degree(1)
    assert elevated.degree == 2
    assert elevated.control_points[0] == pytest.approx((0.0, 0.0, 0.0))
    assert elevated.control_points[1] == pytest.approx((1.0, 2.0, -3.0))
    assert elevated.control_points[2] == pytest.approx((2.0, 4.0, -6.0))
    assert _max_deviation(line, elevated) < 1e-15


def test_elevating_twice_equals_elevating_by_two() -> None:
    """The operation is additive in degree, to round-off, on a non-rational curve."""
    original = NURBSCurve(
        degree=3,
        knots=(0.0, 0.0, 0.0, 0.0, 0.4, 1.0, 1.0, 1.0, 1.0),
        control_points=((0, 0, 0), (1, 2, 0), (2, 0, 0), (3, 3, 0), (4, 1, 0)),
    )
    twice = original.elevate_degree(1).elevate_degree(1)
    direct = original.elevate_degree(2)
    assert twice.degree == direct.degree == 5
    assert len(twice.control_points) == len(direct.control_points)
    assert _max_deviation(twice, direct) < 1e-12


def test_zero_times_is_identity_and_negative_times_is_rejected() -> None:
    curve = circle_curve(RADIUS, arcs=3)
    assert curve.elevate_degree(0) is curve
    with pytest.raises(NURBSError):
        curve.elevate_degree(-1)
    with pytest.raises(NURBSError):
        curve.elevate_degree(1.5)  # type: ignore[arg-type]

def test_continuity_is_preserved_and_the_representation_grows_minimally() -> None:
    """A C2 join must stay C2: the naive per-span Bezier route would demote it to C0.

    The cubic below has one interior knot of multiplicity 1 (declared C2). Elevating
    to degree 4 must give that knot multiplicity 2 and 7 control points; a per-span
    Bezier construction would instead return multiplicity 4 / 9 control points.
    """
    original = NURBSCurve(
        degree=3,
        knots=(0.0, 0.0, 0.0, 0.0, 0.4, 1.0, 1.0, 1.0, 1.0),
        control_points=((0, 0, 0), (1, 2, 0), (2, 0, 0), (3, 3, 0), (4, 1, 0)),
    )
    elevated = original.elevate_degree(1)
    assert elevated.degree == 4
    assert elevated.knots.count(0.4) == 2
    assert len(elevated.control_points) == 7
    assert elevated.continuity(order=2)["declared_minimum_order"] == 2
    assert _max_deviation(original, elevated) < 1e-12


def test_discontinuous_curve_is_rejected_with_a_named_limitation() -> None:
    """An interior multiplicity of ``degree + 1`` is a break, not a C0 join.

    The two sides share no control point, so the single-curve collocation matrix is
    singular. The operation refuses with the remedy (split first) instead of returning
    a plausible-looking wrong curve.
    """
    broken = NURBSCurve(
        degree=3,
        knots=(0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0),
        control_points=(
            (0, 0, 0), (1, 1, 0), (2, 0, 0), (3, 1, 0),
            (4, 2, 0), (5, 3, 0), (6, 1, 0), (7, 0, 0),
        ),
    )
    with pytest.raises(NURBSError, match="discontinuous"):
        broken.elevate_degree(1)
