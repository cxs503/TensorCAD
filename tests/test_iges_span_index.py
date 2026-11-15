"""Regression tests for the knot-span index used by the type-128/126 evaluators.

The defect (fixed 2026-11-15, anchor GAP-CAD-T1-BASELINE-FLEET): ``_find_span``
early-outed at ``knots[n]`` and returned ``n-1``, so from the second-to-last support
knot onwards the span was one too small.  ``_basis`` then indexed
``span - degree + i`` from ``-1``, Python wrapped that to the end of the knot vector,
every basis function came back ``0.0``, and ``evaluate`` returned ``(0, 0, 0)`` --
a silent wrong answer, not an exception.  A direction with a single knot span (a
Bezier patch, ``count == degree + 1``) hits it over its whole domain: measured on the
official KVLCC2 IGES, 5 of 19 surfaces evaluated to the origin and 6 more were
distorted at their boundary rings.
"""

from __future__ import annotations

import math
import os

import pytest

from tensorcad.iges import BSplineCurve, BSplineSurface, _basis, _find_span

#: official archives are cached outside the repository (no redistribution licence);
#: the tests below run only when the bytes happen to be there.
KVLCC2_CANDIDATES = (
    "/root/sim-assets/m20163-fedab06121c6/artifacts/baseline_fleet/raw/KVLCC2/kvlcc2.igs",
    "/root/sim-assets/m20045-285e6d47d6f1/.cache/baseline_fleet/iges/kvlcc2.igs",
)


def _bezier_patch() -> BSplineSurface:
    """A flat patch whose v direction is a single (Bezier) span, degree 2."""
    count_u, count_v, degree = 3, 3, 2
    knots = (0.0, 0.0, 0.0, 2.0, 2.0, 2.0)  # one span in both directions
    control = []
    for j in range(count_v):
        for i in range(count_u):
            control.append((float(i), float(j), 0.5))
    return BSplineSurface(
        degree_u=degree, degree_v=degree, knots_u=knots, knots_v=knots,
        control=tuple(control), weights=(1.0,) * (count_u * count_v),
        count_u=count_u, count_v=count_v,
    )


class TestFindSpan:
    def test_single_span_returns_the_last_index(self):
        """n = count-1 is the span of a Bezier direction, not n-1."""
        knots = (0.0, 0.0, 0.0, 2.0, 2.0, 2.0)
        assert _find_span(2, 2, knots, 1.0) == 2
        assert _find_span(2, 2, knots, 0.0) == 2 or _find_span(2, 2, knots, 0.0) == 2

    def test_interior_span_is_unchanged(self):
        knots = (0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0, 3.0)
        assert _find_span(5, 3, knots, 0.5) == 3
        assert _find_span(5, 3, knots, 1.5) == 4
        # the value that used to collapse the span: the second-to-last support knot
        assert _find_span(5, 3, knots, 2.5) == 5
        assert _find_span(5, 3, knots, 3.0) == 5

    def test_basis_partition_of_unity_on_the_last_span(self):
        knots = (0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0, 3.0)
        values, _ = _basis(5, 3, knots, 2.5)
        assert math.isclose(sum(values), 1.0, abs_tol=1e-12)
        assert max(values) > 0.0


class TestSurfaceEvaluation:
    def test_bezier_direction_is_not_the_origin(self):
        s = _bezier_patch()
        p = s.evaluate(1.0, 1.0)
        assert p != (0.0, 0.0, 0.0)
        assert math.isclose(p[0], 1.0, abs_tol=1e-12)
        assert math.isclose(p[1], 1.0, abs_tol=1e-12)
        assert math.isclose(p[2], 0.5, abs_tol=1e-12)

    def test_bezier_domain_endpoints_interpolate_the_corners(self):
        s = _bezier_patch()
        assert s.evaluate(0.0, 0.0) == (0.0, 0.0, 0.5)
        assert s.evaluate(2.0, 2.0) == (2.0, 2.0, 0.5)

    def test_every_sample_stays_inside_the_control_net(self):
        s = _bezier_patch()
        (lo), (hi) = s.bounding_box()
        for a in range(9):
            for b in range(9):
                p = s.evaluate(2.0 * a / 8, 2.0 * b / 8)
                for k in range(3):
                    assert lo[k] - 1e-12 <= p[k] <= hi[k] + 1e-12

    def test_curve_evaluates_on_a_bezier_span(self):
        curve = BSplineCurve(
            degree=2, knots=(0.0, 0.0, 0.0, 2.0, 2.0, 2.0),
            control=((0.0, 0.0), (1.0, 2.0), (2.0, 0.0)),
            weights=(1.0, 1.0, 1.0), dim=2,
        )
        mid = curve.evaluate(1.0)
        assert mid != (0.0, 0.0)
        assert math.isclose(mid[0], 1.0, abs_tol=1e-12)
        assert math.isclose(mid[1], 1.0, abs_tol=1e-12)


@pytest.mark.parametrize("path", KVLCC2_CANDIDATES)
def test_official_kvlcc2_has_no_degenerate_surface(path):
    if not os.path.exists(path):
        pytest.skip("official KVLCC2 IGES not cached in this environment")
    from tensorcad.iges import read_iges

    iges = read_iges(path)
    surfaces = iges.surfaces()
    assert len(surfaces) == 19
    for index, s in enumerate(surfaces):
        (x0, y0, z0), (x1, y1, z1) = s.bounding_box()
        u0, u1 = s.domain_u
        v0, v1 = s.domain_v
        for a in range(13):
            for b in range(13):
                p = s.evaluate(u0 + (u1 - u0) * a / 12, v0 + (v1 - v0) * b / 12)
                assert p != (0.0, 0.0, 0.0), f"surface {index} evaluates to the origin"
                assert x0 - 1e-6 <= p[0] <= x1 + 1e-6, f"surface {index} outside its control net"
                assert y0 - 1e-6 <= p[1] <= y1 + 1e-6, f"surface {index} outside its control net"
                assert z0 - 1e-6 <= p[2] <= z1 + 1e-6, f"surface {index} outside its control net"
