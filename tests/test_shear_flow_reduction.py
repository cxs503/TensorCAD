"""The reduction of the two shear-flow solvers to a *third*, hand-checkable answer.

``tests/test_shear_flow_multicell.py`` compares the graph solver with the merged one. That is
a consistency check between two implementations of the same idealisation, and it cannot catch
a shared misreading of the physics. This module adds the two routes that can:

1. **Closed form** (:func:`closed_form_box`) -- for the symmetric rectangular box under a
   vertical shear through the centre, Bredt's circulation vanishes by symmetry and the answer
   is the first moment measured from the *symmetry cut*, which is arithmetic a reader can
   check on paper. The numbers are written out in ``docs/shear_flow_reduction.md``.
2. **The differential relation** (:class:`PerimeterOracle`) -- ``dq/ds = -(V/I) t (v - v_c)``
   integrated around the perimeter with the no-twist condition as the only closure. It forms
   no first moment at all, so it shares no algebra with either solver.

Both solvers are checked against both routes, and the oracle against itself (a cut may be
moved, and a loop may be walked either way, without changing the field).
"""
from __future__ import annotations

import math

import pytest
from shear_flow_reduction import (
    PerimeterOracle,
    canonical_q,
    canonical_span,
    difference_on_points,
    field_at,
    flow_scale,
    points_along,
    relative_error,
)

from tensorcad.shear_flow import ThinWalledSection, Wall, shear_flow
from tensorcad.shear_flow_multicell import (
    MultiCellSection,
    multi_cell_box_section,
    multi_cell_shear_flow,
)

WIDTH = 4.0
HEIGHT = 3.0
THICKNESS = 0.02
FORCE = 8000.0

#: The rectangular box the whole file is about, walked counter-clockwise from the origin.
BOX_CHAIN = (
    Wall(((0.0, 0.0), (WIDTH, 0.0)), THICKNESS, "bottom"),
    Wall(((WIDTH, 0.0), (WIDTH, HEIGHT)), THICKNESS, "side_s"),
    Wall(((WIDTH, HEIGHT), (0.0, HEIGHT)), THICKNESS, "top"),
    Wall(((0.0, HEIGHT), (0.0, 0.0)), THICKNESS, "side_p"),
)


def closed_form_box(site: str, u: float) -> float:
    """The righting flow of the box, from the symmetry cut, in closed form.

    ``I = 2 (t * b * (h/2)^2) + 2 (t * h^3 / 12) = 0.45 m^4`` and ``z_c = h/2`` for this
    body. Cutting the perimeter where it must carry no flow -- the mid-point of either
    flange, by symmetry -- leaves no circulation to determine, so for a point on the
    starboard half of the top flange at distance ``s`` from the mid-point

        ``q = (V / I) t (h/2) s``

    (the first moment of the flange area between the cut and the point) and on the starboard
    web at distance ``u`` below the top corner

        ``q = (V / I) [ t (b/2) (h/2) + t ((h/2) u - u^2 / 2) ]``.

    Both are three lines of arithmetic: at the corner ``V * 0.06 / 0.45 = 1066.67 N/m`` and
    at mid-height ``V * 0.0825 / 0.45 = 1466.67 N/m``.
    """
    b, h, t, v = WIDTH, HEIGHT, THICKNESS, FORCE
    second_moment = 2.0 * t * b * (h / 2.0) ** 2 + 2.0 * t * h**3 / 12.0
    if site == "flange":
        return v * t * (h / 2.0) * u / second_moment
    if site == "web":
        first_moment = t * (b / 2.0) * (h / 2.0) + t * ((h / 2.0) * u - u * u / 2.0)
        return v * first_moment / second_moment
    raise ValueError(f"unknown site {site!r}")


def the_two_solvers():
    merged = shear_flow(ThinWalledSection(BOX_CHAIN), FORCE)
    graph = multi_cell_shear_flow(multi_cell_box_section([WIDTH], HEIGHT, THICKNESS), FORCE)
    return merged, graph


class TestClosedFormAnchor:
    """What a reader can check on paper, checked against both solvers."""

    def test_the_corner_and_mid_height_values_are_the_paper_ones(self):
        """``V * 0.06 / 0.45 = 1066.67`` at the corner, ``V * 0.0825 / 0.45 = 1466.67`` mid-web."""
        assert closed_form_box("web", 0.0) == pytest.approx(1066.6666666666667, rel=1e-12)
        assert closed_form_box("web", HEIGHT / 2.0) == pytest.approx(1466.6666666666667,
                                                                   rel=1e-12)

    @pytest.mark.parametrize("u", [0.2, 0.5, 1.5, 2.5])
    def test_the_web_flow_matches_the_closed_form(self, u):
        merged, graph = the_two_solvers()
        expected = closed_form_box("web", u)
        point = (WIDTH, HEIGHT - u)
        assert field_at(merged, point) == pytest.approx(expected, rel=1e-12)
        assert field_at(graph, point) == pytest.approx(expected, rel=1e-12)

    def test_the_flange_flow_is_linear_from_the_symmetry_cut(self):
        merged, graph = the_two_solvers()
        for s in (0.5, 1.0, 1.5):
            expected = closed_form_box("flange", s)
            point = (WIDTH / 2.0 + s, HEIGHT)  # starboard half of the top flange
            # The sign is the wall's canonical direction (the flow runs -y on that half);
            # the magnitude is the closed form.
            assert abs(field_at(merged, point)) == pytest.approx(expected, rel=1e-12)
            assert abs(field_at(graph, point)) == pytest.approx(expected, rel=1e-12)
        # Both halves: equal and opposite about the symmetry cut, which is why q_0 = 0.
        assert field_at(merged, (1.0, HEIGHT)) == pytest.approx(
            -field_at(merged, (3.0, HEIGHT)), rel=1e-12
        )

    def test_the_circulation_vanishes_by_symmetry(self):
        """The symmetry cut is not a choice of convenience: ``q_0`` is exactly zero here."""
        oracle = PerimeterOracle(BOX_CHAIN, shear_force=FORCE, coordinate=lambda p: p[1],
                                 cut_wall=2, cut_fraction=0.5)
        assert oracle.centroid == pytest.approx(HEIGHT / 2.0, abs=1e-12)
        assert oracle.second_moment == pytest.approx(0.45, rel=1e-12)
        # Not bit-exact zero: the panel sums are floating point, so "zero" is a few 1e-11
        # of the 1467 N/m the section carries.
        assert abs(oracle.circulation) < 1e-9 * FORCE / HEIGHT
        # Cut anywhere else and a circulation is needed to make the twist vanish -- the
        # no-twist condition is doing work, not being satisfied by accident.
        elsewhere = PerimeterOracle(BOX_CHAIN, shear_force=FORCE, coordinate=lambda p: p[1],
                                    cut_wall=0, cut_fraction=0.25)
        assert abs(elsewhere.circulation) > 100.0


class TestDifferentialOracle:
    """The third route: the differential relation, no first moment anywhere."""

    def test_both_solvers_reproduce_the_oracle_on_the_whole_perimeter(self):
        merged, graph = the_two_solvers()
        oracle = PerimeterOracle(BOX_CHAIN, shear_force=FORCE, coordinate=lambda p: p[1],
                                 cut_wall=2, cut_fraction=0.5)
        scale = flow_scale(merged)
        points = points_along(BOX_CHAIN, count=17)
        assert difference_on_points(merged, oracle, points) < 1e-10 * scale
        assert difference_on_points(graph, oracle, points) < 1e-10 * scale

    def test_the_oracle_reproduces_the_applied_shear_force(self):
        oracle = PerimeterOracle(BOX_CHAIN, shear_force=FORCE, coordinate=lambda p: p[1],
                                 cut_wall=1, cut_fraction=0.25)
        assert oracle.resultant[1] == pytest.approx(FORCE, rel=1e-13)
        assert oracle.resultant[0] == pytest.approx(0.0, abs=1e-9)

    @pytest.mark.parametrize("cut_wall,cut_fraction", [(0, 0.25), (1, 0.5), (2, 0.5), (3, 0.75)])
    @pytest.mark.parametrize("direction", [1, -1])
    def test_the_oracle_does_not_depend_on_the_cut_or_the_walk(self, cut_wall, cut_fraction,
                                                               direction):
        """A cut is where you start reading, and a loop can be walked either way."""
        reference = PerimeterOracle(BOX_CHAIN, shear_force=FORCE, coordinate=lambda p: p[1],
                                    cut_wall=2, cut_fraction=0.5)
        other = PerimeterOracle(BOX_CHAIN, shear_force=FORCE, coordinate=lambda p: p[1],
                                cut_wall=cut_wall, cut_fraction=cut_fraction,
                                direction=direction)
        scale = flow_scale(the_two_solvers()[0])
        points = points_along(BOX_CHAIN, count=13)
        assert difference_on_points(reference, other, points) < 1e-10 * scale
        # The circulation is a property of the *cut*, so it must move with it -- and the
        # field must not. Both statements are asserted, not just the convenient one.
        if (cut_wall, cut_fraction) == (2, 0.5):
            assert abs(other.circulation) == pytest.approx(abs(reference.circulation), abs=1e-9)
        else:
            assert abs(other.circulation) > 100.0

    def test_the_flange_flow_is_antisymmetric_about_the_shear_centre(self):
        """Symmetric box: the flange flow about ``y = b/2`` is equal and opposite, both routes.

        That antisymmetry is exactly what makes the applied shear pass through the centroid
        with no torque, i.e. what puts the shear centre at ``y = b/2`` -- and it is not an
        input to either solver.
        """
        merged, graph = the_two_solvers()
        for offset in (0.5, 1.0, 1.75):
            left, right = (WIDTH / 2.0 - offset, HEIGHT), (WIDTH / 2.0 + offset, HEIGHT)
            assert field_at(merged, left) == pytest.approx(-field_at(merged, right), rel=1e-12)
            assert field_at(graph, left) == pytest.approx(-field_at(graph, right), rel=1e-12)


class TestComparatorTraps:
    """The three traps that made the reduction look red (see ``tests/shear_flow_reduction.py``)."""

    def test_a_reversed_wall_reports_the_same_field(self):
        """Trap 1: ``q_at(f)`` means a different point when the arc sense is flipped."""
        merged, _graph = the_two_solvers()
        side = next(edge for edge in merged.edges if edge.label == "side_s")
        flipped = type(side)(
            label=side.label,
            thickness=side.thickness,
            start=side.end,
            end=side.start,
            length=side.length,
            flow_coefficients=(-side.flow_coefficients[0], -side.flow_coefficients[1],
                               -side.flow_coefficients[2]),
            s_values=side.s_values,
            q_values=side.q_values,
            stress_values=side.stress_values,
        )
        for fraction in (0.0, 0.25, 0.5, 0.75, 1.0):
            assert canonical_q(side, fraction) == pytest.approx(
                canonical_q(flipped, fraction), abs=1e-9
            )
        assert canonical_span(side) == canonical_span(flipped)

    def test_the_error_metric_is_scale_relative(self):
        """Trap 2: at a zero crossing the relative error is 100 % and means nothing."""
        measured, reference, scale = 9.09e-13, 4.55e-13, 1466.6666666666667
        naive = abs(measured - reference) / reference
        assert naive > 0.9  # the relative metric calls round-off a total failure
        assert relative_error(measured, reference, scale=scale) < 1e-15

    def test_a_spliced_flange_and_an_unspliced_one_agree(self):
        """Trap 3: two 2 m walls and one 4 m wall are the same flange.

        The *same* single cell, described twice: once with a 4 m top flange and once with it
        split at the mid-point (which is what a body has to look like when an internal web
        meets the flange there). The wall-by-wall comparison cannot see this equality -- the
        wall sets differ -- so it is compared as a field.
        """
        merged, _ = the_two_solvers()
        spliced_walls = (
            Wall(((0.0, 0.0), (WIDTH, 0.0)), THICKNESS, "bottom"),
            Wall(((WIDTH, 0.0), (WIDTH, HEIGHT)), THICKNESS, "side_s"),
            Wall(((WIDTH, HEIGHT), (WIDTH / 2.0, HEIGHT)), THICKNESS, "top_s"),
            Wall(((WIDTH / 2.0, HEIGHT), (0.0, HEIGHT)), THICKNESS, "top_p"),
            Wall(((0.0, HEIGHT), (0.0, 0.0)), THICKNESS, "side_p"),
        )
        spliced = shear_flow(ThinWalledSection(spliced_walls), FORCE)
        spans = {canonical_span(edge) for edge in spliced.edges}
        assert canonical_span(merged.edges[2]) not in spans  # the 4 m top wall is gone
        scale = flow_scale(merged)
        points = points_along(BOX_CHAIN, count=9)
        assert difference_on_points(spliced, merged, points) < 1e-10 * scale
        # And the graph solver on the same spliced body agrees too -- it *must* split walls
        # at junctions, so every multi-cell body it sees is spliced somewhere.
        graph = multi_cell_shear_flow(MultiCellSection(spliced_walls), FORCE)
        assert difference_on_points(graph, merged, points) < 1e-10 * scale
        assert math.isclose(field_at(spliced, (1.0, HEIGHT)), field_at(merged, (1.0, HEIGHT)))
