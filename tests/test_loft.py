"""Tests for :mod:`tensorcad.loft` -- exact skinning through arbitrary section curves.

Every numerical claim is measured against a reference computed here (a parabola, a
circle/cylinder/cone closed form, or a knot-normalised copy of the same curve), never
against a constant copied out of a document.
"""
from __future__ import annotations

import math

import pytest

from tensorcad.loft import (
    Loft,
    LoftError,
    knot_union,
    loft_sections,
    loft_surface,
    section_deviation,
)
from tensorcad.nurbs import NURBSCurve, circle_curve

TWO_PI = 2.0 * math.pi


# ---------------------------------------------------------------------------
# Fixtures: section families built from exact representations
# ---------------------------------------------------------------------------
def parabola_section(x0: float) -> NURBSCurve:
    """Exact parabola ``(x0, 2t - 1, (2t - 1)^2)`` as a degree-2 Bezier."""
    return NURBSCurve(
        degree=2,
        knots=(0.0, 0.0, 0.0, 1.0, 1.0, 1.0),
        control_points=((x0, -1.0, 1.0), (x0, 0.0, -1.0), (x0, 1.0, 1.0)),
    )


def yz_circle(x0: float, radius: float = 0.5, *, arcs: int = 3) -> NURBSCurve:
    """Exact rational circle in the plane ``x = x0`` (a rotation of ``circle_curve``)."""
    base = circle_curve(radius, arcs=arcs)
    return NURBSCurve(
        degree=base.degree,
        knots=base.knots,
        control_points=tuple((x0, point[0], point[1]) for point in base.control_points),
        weights=base.weights,
    )


def yz_ellipse(x0: float, semi: tuple[float, float], *, angle: float = 0.0) -> NURBSCurve:
    """Exact rational ellipse in ``x = x0``, rotated about the ``x`` axis by ``angle``."""
    base = circle_curve(1.0)
    cos_a, sin_a = math.cos(angle), math.sin(angle)
    points = []
    for point in base.control_points:
        y = semi[0] * point[0]
        z = semi[1] * point[1]
        points.append((x0, cos_a * y - sin_a * z, sin_a * y + cos_a * z))
    return NURBSCurve(
        degree=base.degree, knots=base.knots, control_points=tuple(points), weights=base.weights
    )


def line_section(p0: tuple[float, float, float], p1: tuple[float, float, float]) -> NURBSCurve:
    return NURBSCurve(
        degree=1, knots=(0.0, 0.0, 1.0, 1.0), control_points=(p0, p1)
    )


def polyline_section(points: tuple[tuple[float, float, float], ...]) -> NURBSCurve:
    count = len(points)
    interior = [index / (count - 1) for index in range(1, count - 1)]
    knots = (0.0, 0.0, *interior, 1.0, 1.0)
    return NURBSCurve(degree=1, knots=knots, control_points=points)


# ---------------------------------------------------------------------------
# Basic structure
# ---------------------------------------------------------------------------
def test_loft_returns_a_tensor_product_patch_with_the_expected_degree():
    loft = loft_sections([parabola_section(x) for x in (0.0, 1.0, 2.0, 3.0, 4.0)])
    assert isinstance(loft, Loft)
    # five sections -> cubic in the loft direction; parabola is degree 2 transversally
    assert loft.surface.degree == (3, 2)
    assert loft.n_sections == 5
    assert loft.n_v == 3
    assert loft.surface.rational is False


def test_degree_v_raises_the_transverse_degree():
    loft = loft_sections(
        [parabola_section(x) for x in (0.0, 1.0, 2.0)], degree_v=4
    )
    assert loft.surface.degree[1] == 4
    # elevation is exact, so the section is still reproduced everywhere
    measurement = loft.measure()
    assert measurement.relative_interpolation_defect < 1e-14


def test_two_sections_give_a_ruled_patch():
    loft = loft_sections([line_section((0.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
                          line_section((2.0, 0.0, 1.0), (2.0, 1.0, 1.0))])
    assert loft.surface.degree[0] == 1, "two sections cannot support a cubic in u"
    # the ruled surface is the linear blend, exactly
    for v in (0.0, 0.25, 1.0):
        midpoint = loft.surface.evaluate(0.5, v)
        expected = (
            1.0,
            v,
            0.5 * (0.0 + 1.0),
        )
        assert max(abs(midpoint[axis] - expected[axis]) for axis in range(3)) < 1e-15


def test_loft_surface_wrapper_returns_the_same_patch():
    sections = [parabola_section(x) for x in (0.0, 1.0, 2.0)]
    loft = loft_sections(sections)
    surface = loft_surface(sections)
    assert surface.degree == loft.surface.degree
    assert surface.knots == loft.surface.knots
    assert surface.control_net == loft.surface.control_net


# ---------------------------------------------------------------------------
# Exact interpolation of every section
# ---------------------------------------------------------------------------
def test_surface_restricts_to_each_section_curve():
    sections = [parabola_section(x) for x in (0.0, 1.0, 2.0, 3.0, 4.0)]
    loft = loft_sections(sections)
    measurement = loft.measure(samples=101, deviation_samples=41)
    assert measurement.relative_interpolation_defect < 1e-14
    assert measurement.relative_deviation < 1e-14
    for defect in measurement.interpolation_defect:
        assert defect < 1e-14 * measurement.scale


def test_weights_are_carried_row_by_row():
    sections = [yz_circle(x) for x in (0.0, 0.5, 1.0, 1.5)]
    loft = loft_sections(sections)
    for i, section in enumerate(sections):
        for j, weight in enumerate(section.weights):
            assert loft.surface.weights[i][j] == pytest.approx(weight, rel=0.0, abs=1e-15)


def test_polynomial_skin_reproduces_the_analytic_surface_everywhere():
    """The exact skin of parabolas is the parabola surface, not just its sample rows."""
    loft = loft_sections([parabola_section(x) for x in (0.0, 1.0, 2.0, 3.0, 4.0)])
    worst = 0.0
    for i in range(31):
        u = i / 30
        for j in range(61):
            v = j / 60
            point = loft.surface.evaluate(u, v)
            worst = max(worst, abs(point[2] - point[1] * point[1]))
    assert worst < 1e-13


def test_exact_circle_skin_is_the_cylinder():
    radius, length = 0.5, 2.0
    sections = [yz_circle(x, radius) for x in (0.0, 0.5, 1.0, 1.5, length)]
    loft = loft_sections(sections)
    worst = 0.0
    for i in range(21):
        u = i / 20
        for j in range(91):
            v = j / 90
            point = loft.surface.evaluate(u, v)
            worst = max(worst, abs(math.hypot(point[1], point[2]) - radius))
    assert worst < 1e-14


def test_cylinder_area_matches_two_pi_r_h():
    radius, length = 0.5, 2.0
    loft = loft_sections([yz_circle(x, radius) for x in (0.0, 0.5, 1.0, 1.5, length)])
    area = loft.surface.area(order=24)
    analytic = TWO_PI * radius * length
    assert abs(area - analytic) / analytic < 1e-12


def test_exact_circle_skin_is_the_cone():
    radii = {0.0: 0.6, 0.5: 0.5, 1.0: 0.4, 1.5: 0.3, 2.0: 0.2}
    sections = [yz_circle(x, radii[x]) for x in sorted(radii)]
    loft = loft_sections(sections)
    worst = 0.0
    for i in range(21):
        u = i / 20
        for j in range(91):
            point = loft.surface.evaluate(u, j / 90)
            expected = 0.6 - 0.2 * point[0]  # linear taper over x in [0, 2]
            worst = max(worst, abs(math.hypot(point[1], point[2]) - expected))
    assert worst < 1e-13


def test_cone_area_matches_the_frustum_closed_form():
    radii = {0.0: 0.6, 1.0: 0.4, 2.0: 0.2}
    loft = loft_sections([yz_circle(x, radii[x]) for x in sorted(radii)])
    area = loft.surface.area(order=32)
    slant = math.sqrt(2.0 ** 2 + (0.6 - 0.2) ** 2)
    analytic = math.pi * (0.6 + 0.2) * slant
    assert abs(area - analytic) / analytic < 1e-9


# ---------------------------------------------------------------------------
# Exact re-representation: knot normalisation, insertion, degree elevation
# ---------------------------------------------------------------------------
def test_knot_normalisation_does_not_move_the_curve():
    circle = circle_curve(0.75, arcs=4)
    scaled = NURBSCurve(
        degree=circle.degree,
        knots=tuple(3.0 + 5.0 * value for value in circle.knots),
        control_points=circle.control_points,
        weights=circle.weights,
    )
    loft = loft_sections([scaled, yz_circle(1.0, 0.75, arcs=3)])
    measurement = loft.measure(samples=65, deviation_samples=17)
    assert measurement.relative_interpolation_defect < 1e-14


def test_sections_with_different_knot_vectors_are_both_reproduced():
    """arcs=3 and arcs=4 circles have different interior knots; the union must be exact."""
    sections = [yz_circle(0.0, 0.5, arcs=3), yz_circle(1.0, 0.5, arcs=4)]
    loft = loft_sections(sections)
    measurement = loft.measure(samples=91, deviation_samples=41)
    assert measurement.relative_interpolation_defect < 1e-14
    # ... but the *interior* of the strip is not a circle: the two sections are different
    # parametrisations of the same circle, and the skin blends them point-to-point.  This
    # is the documented limitation (no section reparametrisation), pinned here so it
    # cannot be mistaken for exactness.
    interior = loft.surface.evaluate(0.5, 0.37)
    assert abs(math.hypot(interior[1], interior[2]) - 0.5) > 1e-9


def test_identically_parametrised_sections_blend_to_the_same_circle():
    sections = [yz_circle(x, 0.5, arcs=3) for x in (0.0, 0.5, 1.0)]
    loft = loft_sections(sections)
    for i in range(9):
        point = loft.surface.evaluate(i / 8, 0.37)
        assert abs(math.hypot(point[1], point[2]) - 0.5) < 1e-14


def test_mixed_degrees_are_unified_by_elevation():
    polyline = polyline_section(((0.0, 0.0, 0.0), (0.0, 0.5, 0.1), (0.0, 1.0, 0.0)))
    ellipse = yz_ellipse(2.0, (0.5, 0.4))
    loft = loft_sections([polyline, ellipse])
    assert loft.surface.degree[1] == 2, "the degree-1 section must be elevated to degree 2"
    measurement = loft.measure(samples=81, deviation_samples=17)
    assert measurement.relative_interpolation_defect < 1e-14


def test_polyline_sections_keep_their_corners():
    left = polyline_section(((0.0, 0.0, 0.0), (0.0, 0.5, 0.5), (0.0, 1.0, 0.0)))
    right = polyline_section(((2.0, 0.0, 0.0), (2.0, 0.5, 0.5), (2.0, 1.0, 0.0)))
    loft = loft_sections([left, right])
    # the ridge at v = 0.5 is a corner: evaluate on both sides of it
    for u in (0.0, 0.3, 1.0):
        left_side = loft.surface.evaluate(u, 0.4)
        right_side = loft.surface.evaluate(u, 0.6)
        assert abs(left_side[2] - 0.5) > 1e-3
        assert left_side[2] == pytest.approx(right_side[2], rel=0.0, abs=1e-15)


# ---------------------------------------------------------------------------
# Continuity in the loft direction -- and the same measurement can fail
# ---------------------------------------------------------------------------
def _tapered_family() -> list[NURBSCurve]:
    """Sections whose radius varies non-polynomially along the loft direction."""
    return [
        yz_circle(x, 0.5 + 0.15 * math.cos(1.3 * x)) for x in (0.0, 1.0, 2.0, 3.0, 4.0)
    ]


def test_cubic_skin_is_c2_in_the_loft_direction():
    measurement = loft_sections(_tapered_family(), degree_u=3).measure()
    assert measurement.interior_knots_u, "the fixture must have interior knots to test"
    assert measurement.relative_second_derivative_jump < 1e-9
    assert measurement.relative_first_derivative_jump < 1e-9
    assert measurement.measured_continuous_order == 2


def test_quadratic_skin_reports_a_second_derivative_jump():
    measurement = loft_sections(_tapered_family(), degree_u=2).measure()
    assert measurement.relative_second_derivative_jump > 1e-3
    assert measurement.relative_first_derivative_jump < 1e-9
    assert measurement.measured_continuous_order == 1


def test_piecewise_linear_skin_reports_a_first_derivative_jump():
    measurement = loft_sections(_tapered_family(), degree_u=1).measure()
    assert measurement.relative_first_derivative_jump > 1e-3
    assert measurement.measured_continuous_order == 0


# ---------------------------------------------------------------------------
# Sections that are not planar-parallel
# ---------------------------------------------------------------------------
def test_rotated_sections_are_interpolated_exactly():
    sections = [
        yz_ellipse(x, (0.5, 0.3), angle=0.35 * x) for x in (0.0, 0.7, 1.4, 2.1, 2.8)
    ]
    loft = loft_sections(sections)
    measurement = loft.measure(samples=91, deviation_samples=33)
    assert measurement.relative_interpolation_defect < 1e-14
    assert measurement.measured_continuous_order == 2


def test_skewed_sections_crossing_in_space_are_not_silently_rejected():
    """Sections need not be nested: a sheared stack lifts the transverse rows."""
    sections = [
        line_section((x, 0.0, 0.0), (x, 1.0, 0.3 * (x - 2.0)))
        for x in (0.0, 1.0, 2.0, 3.0)
    ]
    loft = loft_sections(sections)
    measurement = loft.measure(samples=41, deviation_samples=17)
    assert measurement.relative_interpolation_defect < 1e-14
    assert all(math.isfinite(value) for value in loft.surface.evaluate(0.5, 0.5))


def test_section_deviation_is_a_geometric_measurement():
    sections = [yz_circle(x, 0.4) for x in (0.0, 1.0, 2.0)]
    loft = loft_sections(sections)
    # the middle section equals the surface row, so the distance must be round-off
    worst, mean = section_deviation(loft.surface, loft.u_parameters[1], sections[1], samples=33)
    assert worst < 1e-13
    assert mean < 1e-13


# ---------------------------------------------------------------------------
# Parametrisation and reporting
# ---------------------------------------------------------------------------
def test_loft_parameters_are_anchored_and_increasing():
    loft = loft_sections([parabola_section(x) for x in (0.0, 1.0, 2.0, 3.0)])
    params = loft.u_parameters
    assert params[0] == 0.0 and params[-1] == 1.0
    assert all(params[i + 1] > params[i] for i in range(len(params) - 1))


def test_measurement_serialises_to_plain_json_types():
    measurement = loft_sections(_tapered_family()).measure()
    payload = measurement.to_dict()
    assert payload["degree"] == [3, 2]
    assert payload["n_sections"] == 5
    assert isinstance(payload["interpolation_defect"], list)
    assert payload["relative_interpolation_defect"] < 1e-14
    assert payload["measured_continuous_order"] == 2


def test_knot_union_is_the_maximum_degree_and_multiplicity():
    degree, knots = knot_union(
        [yz_circle(0.0, 0.5, arcs=3), yz_circle(1.0, 0.5, arcs=4)], degree=3
    )
    assert degree == 3
    assert knots[0] == 0.0 and knots[-1] == 1.0
    for value in (0.25, 1.0 / 3.0, 0.5, 2.0 / 3.0, 0.75):
        assert value in knots, f"union knot vector is missing {value}"


# ---------------------------------------------------------------------------
# Failures are first-class
# ---------------------------------------------------------------------------
def test_single_section_is_rejected():
    with pytest.raises(LoftError, match="at least 2 sections"):
        loft_sections([parabola_section(0.0)])


def test_non_curve_section_is_rejected():
    with pytest.raises(LoftError, match="NURBSCurve"):
        loft_sections([parabola_section(0.0), ((0.0, 0.0, 0.0), (1.0, 1.0, 1.0))])  # type: ignore[list-item]


def test_zero_degree_is_rejected():
    with pytest.raises(LoftError, match="degree_u"):
        loft_sections([parabola_section(0.0), parabola_section(1.0)], degree_u=0)


def test_coincident_sections_are_rejected():
    section = yz_circle(0.0, 0.5)
    with pytest.raises(LoftError, match="coincide"):
        loft_sections([section, section])


def test_a_break_in_a_section_names_the_split_remedy():
    # degree 1 with an interior knot of multiplicity 2 = degree + 1 = a genuine break
    broken = NURBSCurve(
        degree=1,
        knots=(0.0, 0.0, 0.5, 0.5, 1.0, 1.0),
        control_points=((0.0, 0.0, 0.0), (0.0, 0.3, 0.0), (0.0, 0.6, 0.2), (0.0, 1.0, 0.0)),
    )
    with pytest.raises(LoftError, match="split_at_discontinuities"):
        loft_sections([broken, polyline_section(((1.0, 0.0, 0.0), (1.0, 1.0, 0.0)))])


def test_error_messages_name_the_failing_section_index():
    broken = NURBSCurve(
        degree=1,
        knots=(0.0, 0.0, 0.5, 0.5, 1.0, 1.0),
        control_points=((0.0, 0.0, 0.0), (0.0, 0.3, 0.0), (0.0, 0.6, 0.2), (0.0, 1.0, 0.0)),
    )
    with pytest.raises(LoftError, match="section 2"):
        loft_sections(
            [
                polyline_section(((0.0, 0.0, 0.0), (0.0, 1.0, 0.0))),
                polyline_section(((1.0, 0.0, 0.0), (1.0, 1.0, 0.0))),
                broken,
            ]
        )