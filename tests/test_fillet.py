"""Fillet and chamfer: exact values first, then the discretisation the mesh actually sees.

The tests are ordered by what they can prove.  The `SectionMoments` block checks the
Green's-theorem machinery against hand-integrable shapes (a square, a triangle, a translated
square through the parallel-axis theorem).  The fillet block checks the *operation* against
closed forms that are independent of that machinery -- the rounded-square area
``a^2 - (4 - pi) r^2``, and the circle limit ``r = a/2``, where a filleted square must be a
circle to round-off in area, perimeter, centroid and second moment.  The failure block pins the
refusal messages and the clamp semantics, including the short-edge limit where two fillet arcs
meet tangentially on a vanished edge.  The solid block is the mesh side: watertightness, Euler
characteristic, and a *measured* second-order convergence of the mesh volume to the closed form.
"""

from __future__ import annotations

import json
import math

import pytest

from tensorcad.fillet import (
    SCHEMA_FILLET,
    SCHEMA_FILLET_SOLID,
    FilletError,
    chamfer_section,
    chamfer_solid,
    fillet_section,
    fillet_solid,
    max_chamfer_setback,
    max_fillet_radius,
    polygon_moments,
    rectangle_section,
    section_turns,
)

SQUARE = rectangle_section(2.0, 2.0, centred=True)
SQUARE_AREA = 4.0
SQUARE_PERIMETER = 8.0
L_SHAPE = ((0.0, 0.0), (2.0, 0.0), (2.0, 1.0), (1.0, 1.0), (1.0, 2.0), (0.0, 2.0))
L_AREA = 3.0


def _fillet_area_change(radius: float) -> float:
    """Closed form ``(4 - pi) r^2`` for four right-angle corners."""
    return (4.0 - math.pi) * radius * radius


def _order(errors: list[float]) -> list[float]:
    return [
        math.log(errors[index] / errors[index + 1]) / math.log(2.0)
        for index in range(len(errors) - 1)
    ]


# ------------------------------------------------------------------ polygon moments


class TestPolygonMoments:
    def test_unit_square(self) -> None:
        moments = polygon_moments(((0, 0), (1, 0), (1, 1), (0, 1)))
        assert moments.area == pytest.approx(1.0, rel=0.0, abs=1e-15)
        assert moments.perimeter == pytest.approx(4.0, rel=0.0, abs=1e-15)
        assert moments.centroid == pytest.approx((0.5, 0.5), abs=1e-15)
        assert moments.second_u == pytest.approx(1.0 / 3.0, rel=1e-15)
        assert moments.second_v == pytest.approx(1.0 / 3.0, rel=1e-15)
        assert moments.polar_second_moment == pytest.approx(2.0 / 3.0, rel=1e-15)

    def test_right_triangle(self) -> None:
        moments = polygon_moments(((0, 0), (1, 0), (0, 1)))
        assert moments.area == pytest.approx(0.5, rel=0.0, abs=1e-15)
        assert moments.perimeter == pytest.approx(2.0 + math.sqrt(2.0), rel=1e-15)
        assert moments.centroid == pytest.approx((1.0 / 3.0, 1.0 / 3.0), abs=1e-15)
        assert moments.second_u == pytest.approx(1.0 / 12.0, rel=1e-14)
        assert moments.second_v == pytest.approx(1.0 / 12.0, rel=1e-14)

    def test_clockwise_input_is_normalised(self) -> None:
        clockwise = tuple(reversed(rectangle_section(3.0, 1.0)))
        counter = rectangle_section(3.0, 1.0)
        a = polygon_moments(clockwise)
        b = polygon_moments(counter)
        assert a.area == pytest.approx(b.area, rel=1e-15)
        assert a.centroid == pytest.approx(b.centroid, abs=1e-15)
        assert a.second_u == pytest.approx(b.second_u, rel=1e-14)

    def test_parallel_axis_shift(self) -> None:
        base = polygon_moments(SQUARE)
        shift = (5.0, -3.0)
        moved = polygon_moments(tuple((u + shift[0], v + shift[1]) for u, v in SQUARE))
        assert moved.area == pytest.approx(base.area, rel=1e-15)
        assert moved.centroid[0] == pytest.approx(base.centroid[0] + shift[0], abs=1e-14)
        assert moved.centroid[1] == pytest.approx(base.centroid[1] + shift[1], abs=1e-14)
        assert moved.second_u == pytest.approx(
            base.second_u + 2.0 * shift[0] * base.centroid[0] * base.area + shift[0] ** 2 * base.area,
            rel=1e-12,
        )

    def test_far_from_the_origin(self) -> None:
        # The moments are evaluated about a local origin and shifted back, so a section a
        # kilometre away keeps its own second moment instead of losing it to cancellation.
        base = polygon_moments(SQUARE)
        shift = 1e3
        far = polygon_moments(tuple((u + shift, v + shift) for u, v in SQUARE))
        expected = (
            base.second_u
            + 2.0 * shift * base.centroid[0] * base.area
            + shift**2 * base.area
        )
        assert far.second_u == pytest.approx(expected, abs=1e-6)
        assert far.second_u - shift**2 * base.area == pytest.approx(base.second_u, abs=1e-7)


# ------------------------------------------------------------------ fillet exact values


class TestFilletExactArea:
    @pytest.mark.parametrize("radius", [0.0, 0.05, 0.2, 0.5, 0.9, 1.0])
    def test_rounded_square_area_and_perimeter(self, radius: float) -> None:
        shape = fillet_section(SQUARE, radius, arc_points=16)
        assert shape.moments.area == pytest.approx(
            SQUARE_AREA - _fillet_area_change(radius), rel=1e-14, abs=1e-15
        )
        assert shape.moments.perimeter == pytest.approx(
            SQUARE_PERIMETER - 8.0 * radius + 2.0 * math.pi * radius, rel=1e-14
        )

    def test_the_circle_limit(self) -> None:
        # A square of side 2r filleted with radius r *is* the circle of radius r.
        shape = fillet_section(SQUARE, 1.0, arc_points=32)
        assert shape.moments.area == pytest.approx(math.pi, rel=1e-15)
        assert shape.moments.perimeter == pytest.approx(2.0 * math.pi, rel=1e-15)
        assert shape.moments.centroid == pytest.approx((0.0, 0.0), abs=1e-15)
        assert shape.moments.second_u == pytest.approx(math.pi / 4.0, rel=1e-15)
        assert shape.moments.second_v == pytest.approx(math.pi / 4.0, rel=1e-15)

    @pytest.mark.parametrize("sides", [3, 4, 5, 6, 8, 12])
    def test_regular_polygon_filleted_by_its_inradius_is_a_circle(self, sides: int) -> None:
        circumradius = 1.0
        inradius = circumradius * math.cos(math.pi / sides)
        ring = tuple(
            (
                circumradius * math.cos(2.0 * math.pi * index / sides),
                circumradius * math.sin(2.0 * math.pi * index / sides),
            )
            for index in range(sides)
        )
        shape = fillet_section(ring, inradius, arc_points=16)
        assert shape.moments.area == pytest.approx(math.pi * inradius**2, rel=1e-13)
        assert shape.moments.perimeter == pytest.approx(2.0 * math.pi * inradius, rel=1e-13)

    def test_reentrant_corner_adds_material(self) -> None:
        radius = 0.2
        shape = fillet_section(L_SHAPE, radius, arc_points=16)
        # Five convex corners lose, the one re-entrant corner gains the same amount.
        assert shape.moments.area == pytest.approx(
            L_AREA - 4.0 * (1.0 - math.pi / 4.0) * radius**2, rel=1e-14
        )
        kinds = [corner.kind for corner in shape.corners]
        assert kinds.count("reflex") == 1

    def test_a_collinear_vertex_changes_nothing(self) -> None:
        with_midpoint = ((0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0))
        plain = ((0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0))
        a = fillet_section(with_midpoint, 0.2, arc_points=16)
        b = fillet_section(plain, 0.2, arc_points=16)
        assert a.moments.area == pytest.approx(b.moments.area, rel=1e-15)
        assert a.moments.perimeter == pytest.approx(b.moments.perimeter, rel=1e-15)

    def test_zero_radius_is_the_polygon(self) -> None:
        shape = fillet_section(L_SHAPE, 0.0)
        reference = polygon_moments(L_SHAPE)
        assert shape.moments.area == pytest.approx(reference.area, rel=1e-15)
        assert shape.moments.second_u == pytest.approx(reference.second_u, rel=1e-14)
        assert shape.clamped is False

    @pytest.mark.parametrize("radius", [0.1, 0.5, 1.0])
    def test_the_two_area_derivations_agree(self, radius: float) -> None:
        shape = fillet_section(SQUARE, radius)
        assert shape.area_crosscheck_error < 1e-14

    @pytest.mark.parametrize("radius", [0.1, 0.5, 1.0])
    def test_g1_continuity_of_the_arc_tangents(self, radius: float) -> None:
        shape = fillet_section(SQUARE, radius)
        assert shape.g1_error_deg < 1e-9
        assert shape.g1_error_deg < 1.0  # the acceptance bound, in degrees

    def test_turns(self) -> None:
        turns = section_turns(L_SHAPE)
        assert [round(math.degrees(turn), 9) for turn in turns] == [90.0] * 3 + [-90.0] + [90.0] * 2

    def test_report_is_json_serialisable(self) -> None:
        payload = fillet_section(SQUARE, 0.2).to_dict()
        assert payload["schema"] == SCHEMA_FILLET
        assert json.loads(json.dumps(payload))["achieved"] == pytest.approx(0.2)


# ------------------------------------------------------------------ chamfer


class TestChamfer:
    def test_square_chamfer_area_and_perimeter(self) -> None:
        shape = chamfer_section(SQUARE, 0.2)
        assert shape.moments.area == pytest.approx(SQUARE_AREA - 4.0 * 0.5 * 0.2**2, rel=1e-15)
        assert shape.moments.perimeter == pytest.approx(
            SQUARE_PERIMETER - 8.0 * 0.2 + 4.0 * math.sqrt(2.0) * 0.2, rel=1e-14
        )
        assert shape.area_crosscheck_error < 1e-15

    def test_maximum_setback_turns_the_square_into_its_inscribed_diamond(self) -> None:
        shape = chamfer_section(SQUARE, 1.0)
        assert shape.moments.area == pytest.approx(2.0, rel=1e-14)
        assert shape.moments.perimeter == pytest.approx(4.0 * math.sqrt(2.0), rel=1e-14)

    def test_a_chamfer_is_explicitly_not_g1(self) -> None:
        shape = chamfer_section(SQUARE, 0.2)
        assert shape.g1_error_deg == pytest.approx(45.0, abs=1e-9)

    def test_max_chamfer_setback(self) -> None:
        assert max_chamfer_setback(rectangle_section(2.0, 0.4)) == pytest.approx(0.2, rel=1e-15)


# ------------------------------------------------------------------ feasibility and refusals


class TestFeasibility:
    def test_max_radius_of_a_square(self) -> None:
        assert max_fillet_radius(SQUARE) == pytest.approx(1.0, rel=1e-15)

    def test_a_radius_above_the_limit_is_refused_with_the_limit_named(self) -> None:
        with pytest.raises(FilletError) as info:
            fillet_section(SQUARE, 1.5)
        message = str(info.value)
        assert "overruns an edge" in message
        assert "1.0" in message

    def test_clamping_gives_the_circle_at_the_short_edge_limit(self) -> None:
        shape = fillet_section(SQUARE, 5.0, allow_clamp=True, arc_points=16)
        assert shape.clamped is True
        assert shape.achieved == pytest.approx(1.0, rel=1e-15)
        assert shape.moments.area == pytest.approx(math.pi, rel=1e-14)
        assert shape.moments.perimeter == pytest.approx(2.0 * math.pi, rel=1e-14)

    def test_short_edge_limit(self) -> None:
        strip = rectangle_section(2.0, 0.2)
        assert max_fillet_radius(strip) == pytest.approx(0.1, rel=1e-12)
        with pytest.raises(FilletError):
            fillet_section(strip, 0.5)
        clamped = fillet_section(strip, 0.5, allow_clamp=True)
        assert clamped.achieved == pytest.approx(0.1, rel=1e-12)

    def test_a_needle_vertex_cannot_be_filleted_and_falls_back(self) -> None:
        # Base 2.0 and height 1e-6: the three corners are within 1e-6 of a straight line, so
        # the largest radius that fits is of that order, not of the size of the triangle.
        needle = ((0.0, 0.0), (2.0, 0.0), (1.0, 1e-6))
        assert max_fillet_radius(needle) < 1e-5
        with pytest.raises(FilletError):
            fillet_section(needle, 0.1)
        clamped = fillet_section(needle, 0.1, allow_clamp=True)
        assert clamped.clamped is True
        assert clamped.achieved < 1e-5
        # At the limit the sliver is filleted away to its own incircle, which is exactly the
        # geometry of "the radius that consumes every edge".
        inradius = 2.0 * 1e-6 / (2.0 + 2.0 * math.hypot(1.0, 1e-6))
        assert inradius == pytest.approx(clamped.achieved, rel=1e-9)
        assert clamped.moments.area == pytest.approx(math.pi * clamped.achieved**2, rel=1e-6)

    @pytest.mark.parametrize(
        "bad, message",
        [
            (((0.0, 0.0), (1.0, 1.0)), "at least 3 points"),
            (((0.0, 0.0), (1.0, 0.0), (1.0, 0.0), (2.0, 0.0)), "degenerate"),
            (((0.0, 0.0), (1.0, 0.0), (float("nan"), 1.0)), "not finite"),
            (((0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (1.0, -1.0), (0.0, 2.0)), "not a simple polygon"),
        ],
    )
    def test_bad_sections_are_refused(self, bad: tuple, message: str) -> None:
        with pytest.raises(FilletError) as info:
            fillet_section(bad, 0.1)
        assert message in str(info.value)

    @pytest.mark.parametrize("radius", [-0.1, float("nan"), float("inf")])
    def test_bad_radii_are_refused(self, radius: float) -> None:
        with pytest.raises(FilletError):
            fillet_section(SQUARE, radius)

    def test_a_reverse_given_polygon_gets_the_same_answer(self) -> None:
        forward = fillet_section(SQUARE, 0.2)
        backward = fillet_section(tuple(reversed(SQUARE)), 0.2)
        assert backward.moments.area == pytest.approx(forward.moments.area, rel=1e-14)


# ------------------------------------------------------------------ discretisation


class TestSampledOutline:
    def test_area_converges_at_second_order(self) -> None:
        exact = fillet_section(SQUARE, 1.0, arc_points=512).moments.area
        errors = [
            abs(fillet_section(SQUARE, 1.0, arc_points=points).sampled_moments().area - exact)
            for points in (4, 8, 16, 32)
        ]
        assert min(_order(errors)) > 1.9

    def test_second_moment_converges_at_second_order(self) -> None:
        exact = fillet_section(SQUARE, 1.0, arc_points=512).moments.second_u
        errors = [
            abs(fillet_section(SQUARE, 1.0, arc_points=points).sampled_moments().second_u - exact)
            for points in (4, 8, 16, 32)
        ]
        assert min(_order(errors)) > 1.85


# ------------------------------------------------------------------ solids


class TestFilletSolid:
    def test_a_filleted_prism_is_a_closed_solid(self) -> None:
        solid = fillet_solid(SQUARE, 3.0, 0.2, arc_points=8)
        assert solid.watertight is True
        assert solid.edge_manifold is True
        assert solid.euler_characteristic == 2
        assert solid.zero_area_faces == 0
        assert solid.volume_exact == pytest.approx(
            fillet_section(SQUARE, 0.2).moments.area * 3.0, rel=1e-14
        )

    def test_mesh_volume_converges_at_second_order(self) -> None:
        errors = []
        for points in (4, 8, 16, 32):
            solid = fillet_solid(SQUARE, 3.0, 0.5, arc_points=points)
            errors.append(abs(solid.volume_mesh - solid.volume_exact) / solid.volume_exact)
        assert errors[-1] < 1e-4
        assert min(_order(errors)) > 1.9

    def test_variable_radius_exact_volume(self) -> None:
        solid = fillet_solid(SQUARE, 2.0, 0.1, radius_end=0.4, stations=4, arc_points=16)
        start, end = solid.size_start, solid.size_end
        # A(r) = A0 + K r**2, so the volume is A0 L + K L (r0**2 + r0 r1 + r1**2) / 3 with A0
        # the *unfilleted* area -- the closed form the module reports, recomputed by hand.
        expected = (
            polygon_moments(SQUARE).area * 2.0
            + solid.area_coefficient * 2.0 * (start**2 + start * end + end**2) / 3.0
        )
        assert solid.volume_exact == pytest.approx(expected, rel=1e-13)
        assert solid.volume_mesh == pytest.approx(solid.volume_exact, rel=1e-4)
        assert solid.watertight is True

    def test_a_tapered_fillet_stays_tangent_to_the_flat_faces(self) -> None:
        # Tangency is a property of the surface, not of the chord count: the ruling direction
        # lies in the flat face, so the deviation stays at round-off for any radius profile.
        solid = fillet_solid(SQUARE, 2.0, 0.05, radius_end=0.9, stations=6, arc_points=4)
        assert solid.g1_error_deg < 1e-9
        assert solid.g1_error_deg < 1.0

    def test_the_chord_crease_shrinks_with_the_arc_resolution(self) -> None:
        coarse = fillet_solid(SQUARE, 1.0, 0.5, arc_points=4)
        fine = fillet_solid(SQUARE, 1.0, 0.5, arc_points=64)
        assert coarse.sampled_crease_deg > 5.0
        assert fine.sampled_crease_deg < 1.5
        assert fine.sampled_crease_deg == pytest.approx(coarse.sampled_crease_deg / 16.0, rel=1e-6)

    def test_a_reflex_section_can_be_swept(self) -> None:
        solid = fillet_solid(L_SHAPE, 1.0, 0.1, arc_points=8)
        assert solid.watertight is True
        assert solid.euler_characteristic == 2
        assert solid.volume_mesh == pytest.approx(solid.volume_exact, rel=1e-4)

    def test_clamping_preserves_the_radius_ratio(self) -> None:
        solid = fillet_solid(SQUARE, 2.0, 0.5, radius_end=2.0, allow_clamp=True)
        assert solid.clamped is True
        assert solid.size_end / solid.size_start == pytest.approx(4.0, rel=1e-12)
        assert solid.size_end == pytest.approx(max_fillet_radius(SQUARE), rel=1e-12)

    def test_an_overrunning_sweep_is_refused(self) -> None:
        with pytest.raises(FilletError):
            fillet_solid(SQUARE, 2.0, 0.5, radius_end=5.0)

    def test_a_zero_radius_sweep_is_refused(self) -> None:
        with pytest.raises(FilletError):
            fillet_solid(SQUARE, 2.0, 0.0)

    def test_solid_report_is_json_serialisable(self) -> None:
        payload = fillet_solid(SQUARE, 1.0, 0.1).to_dict()
        assert payload["schema"] == SCHEMA_FILLET_SOLID
        assert json.loads(json.dumps(payload))["watertight"] is True


class TestChamferSolid:
    def test_a_chamfered_prism_is_closed_and_exact(self) -> None:
        solid = chamfer_solid(SQUARE, 3.0, 0.2)
        assert solid.watertight is True
        assert solid.euler_characteristic == 2
        assert solid.volume_mesh == pytest.approx(solid.volume_exact, rel=1e-14)
        assert solid.volume_exact == pytest.approx(
            chamfer_section(SQUARE, 0.2).moments.area * 3.0, rel=1e-14
        )

    def test_variable_chamfer(self) -> None:
        solid = chamfer_solid(SQUARE, 2.0, 0.1, setback_end=0.5, stations=3)
        assert solid.watertight is True
        assert solid.volume_mesh == pytest.approx(solid.volume_exact, rel=1e-12)


# ------------------------------------------------------------------ ship application


def _rounded_rectangle_inertia(width: float, height: float, radius: float) -> float:
    """`int u**2 dA` of a centred rounded rectangle, by strip integration.

    Independent of the module's Green's-theorem route: for each `u` the vertical extent is
    analytic (flat top in the middle, a circular cap over the fillet), and Simpson integrates
    the two smooth pieces.
    """
    half_width = width / 2.0
    half_height = height / 2.0
    flat = half_width - radius
    steps = 20000

    def top(u: float) -> float:
        if abs(u) <= flat:
            return half_height
        return (half_height - radius) + math.sqrt(max(0.0, radius**2 - (abs(u) - flat) ** 2))

    total = 0.0
    for lower, upper in ((0.0, flat), (flat, half_width)):
        if upper <= lower:
            continue
        h = (upper - lower) / steps
        piece = 0.0
        for step in range(steps + 1):
            u = lower + step * h
            weight = 1.0 if step in (0, steps) else (4.0 if step % 2 else 2.0)
            piece += weight * u * u * 2.0 * top(u)
        total += piece * h / 3.0
    return 2.0 * total  # the outline is symmetric in u, so the left half mirrors the right


class TestBilgeRadiusApplication:
    @pytest.mark.parametrize("radius", [0.0, 0.1, 0.25, 0.4])
    def test_waterplane_inertia_matches_strip_integration(self, radius: float) -> None:
        beam, draft = 2.0, 1.0
        shape = fillet_section(rectangle_section(beam, draft, centred=True), radius, arc_points=16)
        reference = _rounded_rectangle_inertia(beam, draft, radius)
        # The fillet radius of a *waterplane* rounds the corners inward, exactly as the section
        # builder does, so the two routes must agree.
        assert shape.moments.second_u == pytest.approx(reference, rel=1e-6)

    def test_the_box_barge_limit(self) -> None:
        beam, draft = 2.0, 1.0
        shape = fillet_section(rectangle_section(beam, draft, centred=True), 0.0)
        # BM = I_T / volume: a box barge recovers B**2 / (12 T) exactly.
        assert shape.moments.second_u / shape.moments.area == pytest.approx(
            beam**2 / (12.0 * draft), rel=1e-14
        )

    def test_a_bilge_radius_reduces_bm(self) -> None:
        beam, draft = 2.0, 1.0
        values = []
        for radius in (0.0, 0.1, 0.2, 0.3):
            shape = fillet_section(rectangle_section(beam, draft, centred=True), radius)
            values.append(shape.moments.second_u / shape.moments.area)
        assert values == sorted(values, reverse=True)
        # The loss leads with r**2 (the corner inertia carries an r**4 term too), so the slope
        # (base - BM) / r**2 falls slowly with r rather than staying constant.
        base = values[0]
        slopes = [(base - value) / radius**2 for radius, value in zip((0.1, 0.2, 0.3), values[1:])]
        assert slopes == sorted(slopes, reverse=True)
        assert slopes[0] == pytest.approx(slopes[1], rel=0.1)
        assert slopes[0] == pytest.approx(0.27, rel=0.05)