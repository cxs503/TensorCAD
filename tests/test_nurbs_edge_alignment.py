"""Tests for compatible edge representations (knot insertion + degree alignment).

The claim under test is narrow and checkable: knot insertion changes the
representation and not the curve, and ``align_edges`` leaves two edges of the *same*
parametrised curve with one representation (equal degree, equal knots, equal control
points and weights) while *reporting* -- not hiding -- the gap when the two edges are
not the same curve.
"""
from __future__ import annotations

import math

import pytest

from tensorcad.edge_alignment import (
    SIDES,
    align_edges,
    boundary_curve,
    edge_report,
    insert_knot,
    reverse_curve,
)
from tensorcad.highorder.nurbs_geom import insert_knot as insert_surface_knot
from tensorcad.highorder.nurbs_geom import split_patch
from tensorcad.nurbs import NURBSCurve, NURBSError, circle_curve, sphere_surface

RADIUS = 1.0


def _parameters(curve: NURBSCurve, samples: int = 257) -> list[float]:
    low, high = curve.domain
    return [low + (high - low) * index / (samples - 1) for index in range(samples)]


def _max_deviation(first: NURBSCurve, second: NURBSCurve) -> float:
    return max(math.dist(first.evaluate(u), second.evaluate(u)) for u in _parameters(first))


def _seam_pair() -> tuple[NURBSCurve, NURBSCurve]:
    """Two patch edges that describe one curve but disagree in degree and knots.

    A sphere is split at ``u = 0.5`` so the two halves share an interface *exactly*.
    The neighbour's representation is then changed the way a real assembly does: one
    extra knot in ``v`` and one degree elevation in ``v``. The shared edge is a curve
    in ``v``, so it now has degree 3 and nine control points instead of degree 2 and
    five -- the same geometry, no shared data.
    """
    surface = sphere_surface(RADIUS)
    left, right = split_patch(surface, 0, 0.5)
    neighbour = insert_surface_knot(right, 1, 0.37, 1).elevate_degree(0, 1)
    return boundary_curve(left, "u1"), boundary_curve(neighbour, "u0")


# -- knot insertion ---------------------------------------------------------
def test_insertion_preserves_the_curve_and_the_domain() -> None:
    curve = circle_curve(2.0, arcs=3)
    inserted = insert_knot(insert_knot(insert_knot(curve, 0.25, 2), 0.7), 0.9)
    assert inserted.degree == curve.degree
    assert inserted.domain == curve.domain
    assert len(inserted.control_points) == len(curve.control_points) + 4
    assert len(inserted.knots) == len(curve.knots) + 4
    assert _max_deviation(curve, inserted) < 1e-12


def test_insertion_grows_the_multiplicity_once_per_insertion() -> None:
    curve = circle_curve(RADIUS, arcs=3)
    before = curve.knots.count(0.5) if 0.5 in curve.knots else 0
    inserted = insert_knot(curve, 0.5, 2)
    assert inserted.knots.count(0.5) == before + 2
    assert len(inserted.knots) == len(curve.knots) + 2


def test_insertion_is_exact_for_a_rational_circle() -> None:
    curve = circle_curve(2.0, arcs=3)
    inserted = insert_knot(curve, 0.25, 2)
    assert inserted.rational
    radius = [math.hypot(*inserted.evaluate(u)[:2]) for u in _parameters(inserted)]
    assert max(abs(value - 2.0) for value in radius) < 1e-12
    assert abs(inserted.arc_length() - 2.0 * math.pi * 2.0) / (2.0 * math.pi * 2.0) < 1e-12


def test_insertion_keeps_a_unit_weight_net_non_rational() -> None:
    """An all-equal weight net must not drift onto the rational path.

    The blend ``alpha * w + (1 - alpha) * w`` returns ``w +/- 1e-16``, which would
    export a "rational" curve whose radii differ in the last bits and whose
    ``rational`` flag changes. The weight is restored exactly instead.
    """
    curve = NURBSCurve(
        degree=2,
        knots=(0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0),
        control_points=((0, 0, 0), (1, 1, 0), (2, 0, 0), (3, 1, 0)),
    )
    inserted = insert_knot(curve, 0.25, 1)
    assert set(inserted.weights) == {1.0}
    assert not inserted.rational
    assert _max_deviation(curve, inserted) < 1e-12


def test_insertion_past_the_degree_is_refused_with_the_remedy() -> None:
    curve = circle_curve(RADIUS, arcs=3)  # interior multiplicity 2 == degree
    assert curve.knots.count(1.0 / 3.0) == 2
    with pytest.raises(NURBSError, match="exceed degree"):
        insert_knot(curve, 1.0 / 3.0, 1)


def test_insertion_outside_the_domain_is_refused() -> None:
    curve = circle_curve(RADIUS, arcs=3)
    with pytest.raises(NURBSError, match="outside the domain"):
        insert_knot(curve, 1.5)


def test_zero_and_negative_insertions() -> None:
    curve = circle_curve(RADIUS, arcs=3)
    assert insert_knot(curve, 0.5, 0) is curve
    with pytest.raises(NURBSError):
        insert_knot(curve, 0.5, -1)
    with pytest.raises(NURBSError):
        insert_knot(curve, 0.5, 1.5)  # type: ignore[arg-type]


# -- boundary curves --------------------------------------------------------
def test_every_boundary_curve_lies_on_the_surface() -> None:
    surface = sphere_surface(RADIUS)
    (u_low, u_high), (v_low, v_high) = surface.domain
    for side in SIDES:
        curve = boundary_curve(surface, side)
        if side in ("v0", "v1"):
            v = v_low if side == "v0" else v_high
            worst = max(
                math.dist(curve.evaluate(u), surface.evaluate(u, v)) for u in _parameters(curve)
            )
        else:
            u = u_low if side == "u0" else u_high
            worst = max(
                math.dist(curve.evaluate(v), surface.evaluate(u, v)) for v in _parameters(curve)
            )
        assert worst < 1e-12, f"{side}: {worst!r}"


def test_boundary_sides_follow_their_parameter_direction() -> None:
    surface = sphere_surface(RADIUS)
    degree_u, degree_v = surface.degree
    assert boundary_curve(surface, "v0").degree == degree_u
    assert boundary_curve(surface, "v1").degree == degree_u
    assert boundary_curve(surface, "u0").degree == degree_v
    assert boundary_curve(surface, "u1").degree == degree_v
    with pytest.raises(ValueError, match="side must be one of"):
        boundary_curve(surface, "u2")


def test_reverse_curve_traces_the_same_points_backwards() -> None:
    curve = circle_curve(RADIUS, arcs=3)
    reversed_curve = reverse_curve(curve)
    low, high = curve.domain
    worst = max(
        math.dist(curve.evaluate(low + high - u), reversed_curve.evaluate(u))
        for u in _parameters(curve)
    )
    assert worst < 1e-12
    assert reversed_curve.domain == curve.domain


# -- alignment --------------------------------------------------------------
def test_alignment_gives_a_real_patch_seam_one_representation() -> None:
    first, second = _seam_pair()
    before = edge_report(first, second)
    assert before["degree"] == [2, 3]
    assert before["control_points"] == [5, 9]
    assert before["knots_equal"] is False
    assert before["watertight"] is False
    # The geometry was always coincident -- only the representation disagreed.
    assert before["max_gap"] < 1e-12

    aligned_first, aligned_second = align_edges(first, second)
    after = edge_report(aligned_first, aligned_second)
    assert after["degree"] == [3, 3]
    assert after["knots_equal"] is True
    assert after["control_points"][0] == after["control_points"][1]
    assert after["max_control_difference"] < 1e-12
    assert after["max_weight_difference"] < 1e-12
    assert after["watertight"] is True
    assert after["max_gap"] < 1e-12


def test_alignment_keeps_each_edge_on_its_own_curve() -> None:
    first, second = _seam_pair()
    aligned_first, aligned_second = align_edges(first, second)
    assert _max_deviation(first, aligned_first) < 1e-12
    assert _max_deviation(second, aligned_second) < 1e-12


def test_alignment_is_minimal() -> None:
    """The aligned edges carry the union of the breakpoints and nothing more."""
    first, second = _seam_pair()
    aligned_first, aligned_second = align_edges(first, second)
    for curve in (aligned_first, aligned_second):
        assert len(curve.control_points) == len(curve.knots) - curve.degree - 1
    union = {knot for knot in first.knots if first.domain[0] < knot < first.domain[1]}
    union |= {knot for knot in second.knots if second.domain[0] < knot < second.domain[1]}
    for value in union:
        assert aligned_first.knots.count(value) == max(
            first.knots.count(value), second.knots.count(value)
        )


def test_alignment_does_not_hide_a_geometric_gap() -> None:
    first = circle_curve(1.0, arcs=3)
    second = circle_curve(1.002, arcs=3)
    aligned_first, aligned_second = align_edges(first, second)
    report = edge_report(aligned_first, aligned_second)
    assert report["same_degree"] is True
    assert report["knots_equal"] is True
    assert report["watertight"] is False
    assert report["max_gap"] == pytest.approx(0.002, rel=1e-6)


def test_alignment_reports_a_weight_mismatch() -> None:
    """Same point set, different parametrisation: the weights are what disagree.

    Three 120-degree arcs and six 60-degree arcs are both the unit circle, and after
    alignment their dehomogenised control polygons agree to round-off -- but the
    weights are 0.75 against ``cos(30 deg)``, so the two curves visit the circle at
    different parameters. Reporting the gap (~3.9% of the radius) is the honest
    outcome; returning ``watertight`` because the polygons matched would not be.
    """
    first = circle_curve(1.0, arcs=3)
    second = circle_curve(1.0, arcs=6)
    aligned_first, aligned_second = align_edges(first, second)
    report = edge_report(aligned_first, aligned_second)
    assert report["knots_equal"] is True
    assert report["max_control_difference"] < 1e-12
    assert report["max_weight_difference"] > 0.1
    assert report["max_gap"] == pytest.approx(0.0385, abs=5e-4)
    assert report["watertight"] is False


def test_alignment_requires_a_common_parameter_domain() -> None:
    first = circle_curve(1.0, arcs=3)
    second = insert_knot(first, 0.5, 1)
    scaled = NURBSCurve(
        degree=second.degree,
        knots=tuple(2.0 * knot for knot in second.knots),
        control_points=second.control_points,
        weights=second.weights,
    )
    with pytest.raises(NURBSError, match="common parameter domain"):
        align_edges(first, scaled)
