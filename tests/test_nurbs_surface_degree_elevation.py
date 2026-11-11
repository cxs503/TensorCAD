"""Surface degree elevation: exactness on the rational anchors, and the traps.

The curve-level operation is verified in ``tests/test_nurbs_degree_elevation.py``.
Here the tensor-product extension is held to the same standard: the elevated patch
must trace the *same* surface, its area must not move, the knot multiplicities must
follow the ISO 10303-42 rule, and the two directions must commute.
"""

from __future__ import annotations

import math
import random

import pytest

from tensorcad.nurbs import (
    NURBSCurve,
    NURBSError,
    NURBSSurface,
    circle_curve,
    cylinder_surface,
    sphere_surface,
)

GRID = [(u / 12.0, v / 7.0) for u in range(13) for v in range(8)]


def _max_deviation(left: NURBSSurface, right: NURBSSurface) -> float:
    worst = 0.0
    for u, v in GRID:
        a = left.evaluate(u, v)
        b = right.evaluate(u, v)
        worst = max(worst, math.dist(a, b))
    return worst


def _elliptic_patch(
    degree=(3, 3),
    knots=((0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0), (0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0)),
    seed=20261111,
) -> NURBSSurface:
    """A non-rational patch with no symmetries, so a wrong net cannot look right."""
    rng = random.Random(seed)
    net = tuple(
        tuple(
            (rng.uniform(-1.0, 1.0), rng.uniform(-1.0, 1.0), 0.15 * (i + j) + rng.uniform(-0.2, 0.2))
            for j in range(len(knots[1]) - degree[1] - 1)
        )
        for i in range(len(knots[0]) - degree[0] - 1)
    )
    return NURBSSurface(degree=degree, knots=knots, control_net=net)


class TestRationalAnchors:
    """The cylinder and the sphere are exact rational surfaces, so they are decisive."""

    def test_cylinder_area_and_points_are_unchanged(self) -> None:
        radius, height = 1.37, 5.2
        cylinder = cylinder_surface(radius, height)
        elevated = cylinder.elevate_degree(2, 1)

        assert elevated.degree == (4, 2)
        assert elevated.rational
        assert _max_deviation(elevated, cylinder) < 1e-13

        analytic = 2.0 * math.pi * radius * height
        assert cylinder.area() == pytest.approx(analytic, rel=1e-12)
        assert elevated.area() == pytest.approx(analytic, rel=1e-12)

    def test_cylinder_radius_is_still_exactly_radius(self) -> None:
        radius = 1.37
        elevated = cylinder_surface(radius, 5.2).elevate_degree(3, 2)
        for u, v in GRID:
            x, y, _ = elevated.evaluate(u, v)
            assert math.hypot(x, y) == pytest.approx(radius, abs=1e-13)

    def test_sphere_area_and_radius_are_unchanged(self) -> None:
        radius = 2.2
        sphere = sphere_surface(radius)
        elevated = sphere.elevate_degree(1, 2)

        assert elevated.degree == (3, 4)
        assert _max_deviation(elevated, sphere) < 1e-13
        for u, v in GRID:
            assert math.sqrt(sum(c * c for c in elevated.evaluate(u, v))) == pytest.approx(
                radius, abs=1e-13
            )
        analytic = 4.0 * math.pi * radius * radius
        assert elevated.area() == pytest.approx(analytic, rel=1e-12)

    def test_normals_keep_their_direction(self) -> None:
        cylinder = cylinder_surface(2.0, 3.0)
        elevated = cylinder.elevate_degree(1, 1)
        for u, v in GRID:
            a = cylinder.normal(u, v)
            b = elevated.normal(u, v)
            assert sum(x * y for x, y in zip(a, b, strict=True)) > 1.0 - 1e-12


class TestRepresentationRules:
    def test_knot_multiplicities_grow_by_times(self) -> None:
        cylinder = cylinder_surface(1.0, 2.0)  # three 120-degree arcs: interior mult 2
        third = 1.0 / 3.0
        assert cylinder.knots[0].count(third) == 2
        elevated = cylinder.elevate_degree(2, 1)
        u_knots = elevated.knots[0]
        assert elevated.degree[0] == 4
        assert u_knots.count(third) == 4
        assert u_knots.count(2.0 * third) == 4
        assert elevated.knots[1] == (0.0, 0.0, 0.0, 1.0, 1.0, 1.0)

    def test_growth_is_the_minimum_times_spans(self) -> None:
        cylinder = cylinder_surface(1.0, 2.0)
        u_count = len(cylinder.control_net)
        v_count = len(cylinder.control_net[0])
        u_spans = len(cylinder.interior_knots[0]) + 1
        elevated = cylinder.elevate_degree(3, 2)
        assert len(elevated.control_net) == u_count + 3 * u_spans
        assert len(elevated.control_net[0]) == v_count + 2 * 1

    def test_a_c1_join_stays_c1(self) -> None:
        patch = NURBSSurface(
            degree=(3, 1),
            knots=((0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0), (0.0, 0.0, 1.0, 1.0)),
            control_net=tuple(
                ((float(i), float(j), float(i * i - j)) for j in range(2)) for i in range(6)
            ),
        )
        assert patch.continuity(order=1)["passed"]
        elevated = patch.elevate_degree(2, 1)
        assert elevated.degree == (5, 2)
        assert elevated.interior_knots[0] == (1.0,)
        assert elevated.knots[0].count(1.0) == 4
        assert elevated.continuity(order=1)["passed"]

    def test_the_two_directions_commute(self) -> None:
        patch = _elliptic_patch()
        u_first = patch.elevate_degree(2, 1)
        v_first = patch.elevate_degree(0, 1).elevate_degree(2, 0)
        assert u_first.degree == v_first.degree
        assert u_first.knots == v_first.knots
        assert _max_deviation(u_first, v_first) < 1e-12

    def test_elevating_twice_in_a_direction_equals_elevating_by_two(self) -> None:
        patch = _elliptic_patch()
        twice = patch.elevate_degree(1, 0).elevate_degree(1, 0)
        direct = patch.elevate_degree(2, 0)
        assert twice.knots == direct.knots
        assert _max_deviation(twice, direct) < 1e-12

    def test_zero_elevation_returns_the_same_object(self) -> None:
        patch = _elliptic_patch()
        assert patch.elevate_degree(0, 0) is patch


class TestNonRationalPatch:
    def test_the_patch_is_reproduced_at_the_same_parameters(self) -> None:
        patch = _elliptic_patch()
        for times_u, times_v in ((1, 1), (2, 3)):
            elevated = patch.elevate_degree(times_u, times_v)
            assert elevated.degree == (3 + times_u, 3 + times_v)
            assert elevated.domain == patch.domain
            assert _max_deviation(elevated, patch) < 1e-12
            assert not elevated.rational

    def test_boundary_curves_are_unchanged(self) -> None:
        patch = _elliptic_patch()
        elevated = patch.elevate_degree(1, 2)
        (u0, u1), (v0, v1) = patch.domain
        for u in (u0, 0.25 * u1 + 0.75 * u0, 0.5 * (u0 + u1), u1):
            for v in (v0, v1):
                assert math.dist(elevated.evaluate(u, v), patch.evaluate(u, v)) < 1e-13
        for v in (v0, 0.5 * (v0 + v1), v1):
            for u in (u0, u1):
                assert math.dist(elevated.evaluate(u, v), patch.evaluate(u, v)) < 1e-13

    def test_derivatives_agree(self) -> None:
        patch = _elliptic_patch()
        elevated = patch.elevate_degree(1, 1)
        for u, v in GRID[::5]:
            left = patch.derivatives(u, v, order=1)
            right = elevated.derivatives(u, v, order=1)
            assert set(left) == set(right)
            for key in left:
                if key == (0, 0):
                    continue
                assert math.dist(left[key], right[key]) < 1e-11


class TestRationalTrap:
    """Elevating the 3-D net and forgetting the weights is a *different* surface."""

    def test_a_unit_weight_circle_is_not_the_circle(self) -> None:
        radius = 2.0
        circle = circle_curve(radius)
        blind = NURBSCurve(degree=2, knots=circle.knots, control_points=circle.control_points)
        elevated = blind.elevate_degree(2)
        worst = max(
            abs(math.dist(elevated.evaluate(i / 41.0), (0.0, 0.0, 0.0)) - radius)
            for i in range(42)
        )
        assert worst > 1e-2  # 500x above the homogeneous route's round-off

    def test_a_unit_weight_cylinder_loses_its_radius(self) -> None:
        radius = 1.0
        cylinder = cylinder_surface(radius, 3.0)
        blind = NURBSSurface(
            degree=cylinder.degree, knots=cylinder.knots, control_net=cylinder.control_net
        )
        elevated = blind.elevate_degree(2, 0)
        worst = max(
            abs(math.hypot(*elevated.evaluate(u, v)[:2]) - radius) for u, v in GRID
        )
        assert worst > 1e-3
        # ... while the weight-carrying route is exact on the same columns
        assert _max_deviation(cylinder.elevate_degree(2, 0), cylinder) < 1e-13


class TestValidation:
    def test_zero_times_are_rejected(self) -> None:
        cylinder = cylinder_surface(1.0, 1.0)
        with pytest.raises(NURBSError, match="times_u"):
            cylinder.elevate_degree(-1, 0)
        with pytest.raises(NURBSError, match="times_v"):
            cylinder.elevate_degree(0, 1.5)  # type: ignore[arg-type]
        with pytest.raises(NURBSError, match="times_u"):
            cylinder.elevate_degree(True, 0)  # type: ignore[arg-type]

    def test_a_discontinuous_direction_is_refused(self) -> None:
        # degree 2 with an interior multiplicity 3 is a break, not a C0 join.
        broken = NURBSSurface(
            degree=(2, 1),
            knots=((0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0), (0.0, 0.0, 1.0, 1.0)),
            control_net=tuple(((float(i), float(j), 0.0) for j in range(2)) for i in range(6)),
        )
        with pytest.raises(NURBSError, match="discontinuous"):
            broken.elevate_degree(1, 0)
        # the other direction is fine, so the error names the culprit
        assert broken.elevate_degree(0, 1).degree == (2, 2)
