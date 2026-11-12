"""Piecewise (break-aware) degree elevation: splitting, exactness, analytic anchors.

The reference object is a degree-2 curve with a real break in it: an exact rational
quarter circle of radius 2 (control polygon ``(2,0), (2,2), (0,2)``, weights
``1, sqrt(2)/2, 1``) followed by a straight segment to ``(-2, 0)``. The join carries
multiplicity ``degree + 1 = 3``, so it is not a ``C0`` join but a break: the two sides
share no control point, the tangent turns by 45 degrees across it, and every
single-curve operation that differentiates (``arc_length``, ``curvature``,
``elevate_degree``) is singular there.

Analytic anchors (no reference implementation involved):

* the arc piece is a circle of radius 2, so ``|P(u)| == 2`` for every ``u``;
* the arc piece has length ``pi R / 2 = pi`` and the straight piece length ``2 sqrt 2``;
* the straight piece is the linear map ``u -> (-2u, 2 - 2u)`` on its own parameter
  interval, derived here by hand from the Bernstein form;
* the kink angle between the two one-sided tangents is 45 degrees.
"""
from __future__ import annotations

import math

import pytest

from tensorcad.nurbs import NURBSCurve, NURBSError, circle_curve

RADIUS = 2.0
ARC_LENGTH = math.pi * RADIUS / 2.0          # pi
LINE_LENGTH = math.hypot(2.0, 2.0)           # 2 sqrt 2
TOTAL_LENGTH = ARC_LENGTH + LINE_LENGTH      # pi + 2 sqrt 2 = 5.970019778335983


def composite_curve() -> NURBSCurve:
    """Quarter circle followed by a straight segment, joined at a genuine break."""
    half_angle_weight = math.cos(math.pi / 4.0)  # sqrt(2)/2: exact rational quarter circle
    return NURBSCurve(
        degree=2,
        knots=(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0),
        control_points=(
            (2.0, 0.0, 0.0),
            (2.0, 2.0, 0.0),
            (0.0, 2.0, 0.0),
            (0.0, 2.0, 0.0),
            (-1.0, 1.0, 0.0),
            (-2.0, 0.0, 0.0),
        ),
        weights=(1.0, half_angle_weight, 1.0, 1.0, 1.0, 1.0),
    )


def jump_curve() -> NURBSCurve:
    """Degree 3, three pieces, two breaks: the joins do not even share a point."""
    knots = (0.0, 0.0, 0.0, 0.0, 0.4, 1.0, 1.0, 1.0, 1.0, 1.5, 1.5, 1.5, 1.5, 2.0, 2.0, 2.0, 2.0)
    points = tuple((float(index), float((index * index) % 5), 0.0) for index in range(len(knots) - 4))
    return NURBSCurve(degree=3, knots=knots, control_points=points)


def interior_parameters(curve: NURBSCurve, samples: int = 97) -> list[float]:
    """Parameters inside the domain, never on a break (evaluation there is one-sided)."""
    low, high = curve.domain
    breaks = {
        value
        for value in curve.knots
        if low < value < high
        and sum(1 for knot in curve.knots if knot == value) > curve.degree
    }
    return [
        low + (high - low) * (index + 0.5) / samples
        for index in range(samples)
        if all(abs(low + (high - low) * (index + 0.5) / samples - b) > 1e-12 for b in breaks)
    ]


# ---------------------------------------------------------------------------
# splitting
# ---------------------------------------------------------------------------
def test_split_returns_one_piece_per_continuous_piece() -> None:
    pieces = composite_curve().split_at_discontinuities()
    assert len(pieces) == 2
    assert [piece.degree for piece in pieces] == [2, 2]
    assert [piece.domain for piece in pieces] == [(0.0, 1.0), (1.0, 2.0)]
    # a pure slice: three control points per single-span quadratic piece
    assert [len(piece.control_points) for piece in pieces] == [3, 3]
    # weights travel untouched, so the arc piece is still the exact rational arc
    assert pieces[0].weights == (1.0, math.cos(math.pi / 4.0), 1.0)
    assert pieces[1].weights == (1.0, 1.0, 1.0)
    assert pieces[0].rational and not pieces[1].rational


def test_split_pieces_reproduce_every_interior_point() -> None:
    curve = composite_curve()
    pieces = curve.split_at_discontinuities()
    deviation = 0.0
    for piece in pieces:
        for u in interior_parameters(curve):
            if piece.domain[0] <= u <= piece.domain[1]:
                deviation = max(deviation, math.dist(curve.evaluate(u), piece.evaluate(u)))
    assert deviation < 1e-14


def test_split_is_a_pure_slice_on_a_three_piece_curve() -> None:
    """Interior knots *inside* a piece are reproduced verbatim -- nothing is inserted."""
    curve = jump_curve()
    pieces = curve.split_at_discontinuities()
    assert len(pieces) == 3
    assert pieces[0].knots == (0.0, 0.0, 0.0, 0.0, 0.4, 1.0, 1.0, 1.0, 1.0)
    assert pieces[1].knots == (1.0, 1.0, 1.0, 1.0, 1.5, 1.5, 1.5, 1.5)
    assert pieces[2].knots == (1.5, 1.5, 1.5, 1.5, 2.0, 2.0, 2.0, 2.0)
    assert [len(piece.control_points) for piece in pieces] == [5, 4, 4]
    # the C2 interior knot of the first piece survives the split
    assert pieces[0].continuity()["passed"] is True
    deviation = max(
        math.dist(curve.evaluate(u), piece.evaluate(u))
        for piece in pieces
        for u in [
            piece.domain[0] + (piece.domain[1] - piece.domain[0]) * (index + 0.5) / 61.0
            for index in range(61)
        ]
    )
    assert deviation == 0.0


def test_a_polyline_with_duplicated_vertices_splits_into_its_segments() -> None:
    """The realistic case: a chine profile stored as one degree-1 polygon."""
    curve = NURBSCurve(
        degree=1,
        knots=(0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0),
        control_points=(
            (0.0, 0.0, 0.0),
            (2.0, 0.0, 0.0),
            (2.0, 0.0, 0.0),
            (2.0, 2.0, 0.0),
            (2.0, 2.0, 0.0),
            (0.0, 2.0, 0.0),
        ),
    )
    pieces = curve.split_at_discontinuities()
    assert [piece.domain for piece in pieces] == [(0.0, 1.0), (1.0, 2.0), (2.0, 3.0)]
    assert sum(piece.arc_length() for piece in pieces) == pytest.approx(6.0, abs=1e-14)
    elevated = curve.elevate_degree_piecewise(1)
    assert [piece.degree for piece in elevated] == [2, 2, 2]
    # elevating a straight degree-1 segment puts its middle control point exactly on
    # the analytic midpoint of that segment
    expected_first = ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0))
    expected_last = ((2.0, 2.0, 0.0), (1.0, 2.0, 0.0), (0.0, 2.0, 0.0))
    for raised, expected in ((elevated[0], expected_first), (elevated[-1], expected_last)):
        for point, target in zip(raised.control_points, expected, strict=True):
            assert math.dist(point, target) < 1e-15


def test_a_smooth_curve_is_returned_as_itself() -> None:
    circle = circle_curve(RADIUS, arcs=3)
    pieces = circle.split_at_discontinuities()
    assert len(pieces) == 1
    assert pieces[0] is circle  # identity, not a copy


def test_the_broken_curve_cannot_be_differentiated_across_the_break() -> None:
    """Why the split exists: the whole curve's derivative net has a zero-length interval."""
    with pytest.raises(NURBSError, match="zero-length knot interval"):
        composite_curve().arc_length()


def test_elevation_of_a_broken_curve_still_refuses_and_names_the_remedy() -> None:
    with pytest.raises(NURBSError, match="discontinuous") as error:
        composite_curve().elevate_degree(1)
    assert "elevate_degree_piecewise" in str(error.value)


# ---------------------------------------------------------------------------
# analytic anchors on the pieces (before and after elevation)
# ---------------------------------------------------------------------------
def test_the_arc_piece_is_the_analytic_circle() -> None:
    arc = composite_curve().split_at_discontinuities()[0]
    radii = [math.hypot(*arc.evaluate(u)[:2]) for u in (index / 200.0 for index in range(201))]
    assert max(abs(radius - RADIUS) for radius in radii) < 1e-14
    assert arc.evaluate(0.0) == pytest.approx((RADIUS, 0.0, 0.0), abs=1e-15)
    assert arc.evaluate(1.0) == pytest.approx((0.0, RADIUS, 0.0), abs=1e-15)
    # the 45-degree point of the rational quadratic is (sqrt 2, sqrt 2)
    assert arc.evaluate(0.5) == pytest.approx((math.sqrt(2.0), math.sqrt(2.0), 0.0), abs=1e-15)


def test_the_straight_piece_is_the_analytic_line() -> None:
    """Independent reference: the Bernstein form of this piece is ``u -> (-2u, 2 - 2u)``."""
    line = composite_curve().split_at_discontinuities()[1]
    assert line.domain == (1.0, 2.0)  # parameters stay on the original ones
    worst = 0.0
    for index in range(101):
        u = index / 100.0
        expected = (-2.0 * u, 2.0 - 2.0 * u, 0.0)  # local parameter u = knot - 1
        worst = max(worst, math.dist(line.evaluate(1.0 + u), expected))
    assert worst < 1e-15


def test_piecewise_arc_length_matches_the_analytic_length() -> None:
    pieces = composite_curve().split_at_discontinuities()
    lengths = [piece.arc_length() for piece in pieces]
    assert lengths[0] == pytest.approx(ARC_LENGTH, abs=1e-13)
    assert lengths[1] == pytest.approx(LINE_LENGTH, abs=1e-13)
    assert sum(lengths) == pytest.approx(TOTAL_LENGTH, abs=1e-13)


def test_the_kink_angle_is_45_degrees_and_survives_elevation() -> None:
    pieces = composite_curve().split_at_discontinuities()
    elevated = composite_curve().elevate_degree_piecewise(3)
    for left, right in (
        (pieces[0].tangent(1.0), pieces[1].tangent(1.0)),
        (elevated[0].tangent(1.0), elevated[1].tangent(1.0)),
    ):
        cosine = max(-1.0, min(1.0, sum(a * b for a, b in zip(left, right, strict=True))))
        assert math.degrees(math.acos(cosine)) == pytest.approx(45.0, abs=1e-12)


# ---------------------------------------------------------------------------
# piecewise elevation
# ---------------------------------------------------------------------------
def test_piecewise_elevation_reproduces_the_geometry_exactly() -> None:
    curve = composite_curve()
    pieces = curve.split_at_discontinuities()
    elevated = curve.elevate_degree_piecewise(2)
    assert [piece.degree for piece in elevated] == [4, 4]
    deviation = 0.0
    for original, raised in zip(pieces, elevated, strict=True):
        low, high = original.domain
        for index in range(201):
            u = low + (high - low) * index / 200.0
            deviation = max(deviation, math.dist(original.evaluate(u), raised.evaluate(u)))
    assert deviation < 1e-14
    # the raised arc is still exactly the rational circle: elevation is exact for
    # rational geometry because the solve happens in homogeneous coordinates
    raised_arc = elevated[0]
    radii = [abs(math.hypot(*raised_arc.evaluate(index / 200.0)[:2]) - RADIUS) for index in range(201)]
    assert max(radii) < 1e-13


def test_piecewise_elevation_keeps_the_arc_length() -> None:
    elevated = composite_curve().elevate_degree_piecewise(3)
    assert sum(piece.arc_length() for piece in elevated) == pytest.approx(TOTAL_LENGTH, abs=1e-12)


def test_growth_is_the_minimum_times_spans_per_piece() -> None:
    """Each piece grows by exactly ``times * (#spans of the piece)`` control points."""
    curve = jump_curve()
    pieces = curve.split_at_discontinuities()
    elevated = curve.elevate_degree_piecewise(1)
    for original, raised in zip(pieces, elevated, strict=True):
        low, high = original.domain
        spans = len({value for value in original.knots if low <= value <= high}) - 1
        assert len(raised.control_points) - len(original.control_points) == spans
        assert raised.degree == original.degree + 1


def test_piecewise_elevation_does_not_close_the_breaks() -> None:
    """A break stays a break: the two sides of each jump stay exactly where they were."""
    curve = jump_curve()
    pieces = curve.split_at_discontinuities()
    elevated = curve.elevate_degree_piecewise(2)
    for index in (0, 1):
        left, right = pieces[index], pieces[index + 1]
        raised_left, raised_right = elevated[index], elevated[index + 1]
        before = math.dist(left.evaluate(left.domain[1]), right.evaluate(right.domain[0]))
        after = math.dist(
            raised_left.evaluate(raised_left.domain[1]),
            raised_right.evaluate(raised_right.domain[0]),
        )
        assert after == pytest.approx(before, abs=1e-14)
        assert before > 1.0  # it really is a jump, not a shared point


def test_piecewise_elevation_of_a_smooth_curve_equals_direct_elevation() -> None:
    circle = circle_curve(RADIUS, arcs=3)
    direct = circle.elevate_degree(2)
    piecewise = circle.elevate_degree_piecewise(2)
    assert len(piecewise) == 1
    assert piecewise[0].degree == direct.degree
    deviation = max(
        math.dist(direct.evaluate(index / 300.0), piecewise[0].evaluate(index / 300.0))
        for index in range(301)
    )
    assert deviation < 1e-13


def test_zero_times_is_the_identity_and_negative_times_is_rejected() -> None:
    curve = composite_curve()
    pieces = curve.split_at_discontinuities()
    unchanged = curve.elevate_degree_piecewise(0)
    assert [piece.degree for piece in unchanged] == [2, 2]
    assert [piece.knots for piece in unchanged] == [piece.knots for piece in pieces]
    for bad in (-1, 1.5, True):
        with pytest.raises(NURBSError, match="non-negative integer"):
            curve.elevate_degree_piecewise(bad)
