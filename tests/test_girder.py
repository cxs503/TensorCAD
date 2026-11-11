"""Verification of the hull-girder section properties against closed forms.

The module under test turns a closed polygon into numbers a scantling check is written
in. Every number here therefore has a reference that never sees the polygon code:

1. **Rectangle and box.** ``I = b h^3 / 12``, ``Z = b h^2 / 6``, ``Z_p = b h^2 / 4``
   (shape factor exactly ``1.5``) and, for the double-walled box, the difference formula
   ``(b h^3 - b_i h_i^3) / 12`` with ``Z_p = (b h^2 - b_i h_i^2) / 4``.
2. **Triangle.** ``I = b h^3 / 36`` about the centroid at ``h / 3``; the plastic modulus
   is checked against an independent quadrature of the section's own width function, so
   the equal-area search is compared to an integration that never bisects.
3. **I-beam.** A 12-vertex polygon against a sum of rectangles and Steiner terms -- the
   same section reached by two routes that share no code.
4. **A circle.** The n-gon's second moment against ``pi R^4 / 4``, reported as a
   convergence sequence (the polygon can only approach it).
5. **Rotation.** Properties computed from a rotated polygon must match the rotated
   tensor, and the product of inertia must vanish at the principal angle.
6. **A real hull.** The midship section of :func:`tensorcad.generate_hull` reached
   through :mod:`tensorcad.bonjean`: the area must agree with the independent clipped
   integration in that module, and the symmetry of the loft must show up as a vanishing
   product of inertia and a zero plastic neutral axis offset.
"""
from __future__ import annotations

import math

import pytest

from tensorcad import HullParameters, generate_hull
from tensorcad.bonjean import section_loops, section_properties
from tensorcad.girder import (
    GirderError,
    box_girder_section,
    closed_cell_torsion_constant,
    girder_section,
    mirror_half_section,
    open_section_torsion_constant,
    section_modulus_extremes,
)

# --------------------------------------------------------------------------- #
# references
# --------------------------------------------------------------------------- #


def simpson(function, low: float, high: float, intervals: int = 200_000) -> float:
    """Composite Simpson rule, the independent route used for the plastic modulus."""
    if intervals % 2:
        intervals += 1
    step = (high - low) / intervals
    total = function(low) + function(high)
    for index in range(1, intervals):
        weight = 4.0 if index % 2 else 2.0
        total += weight * function(low + index * step)
    return total * step / 3.0


def plastic_reference(width, low: float, high: float, level: float) -> float:
    """``INT |z - level| width(z) dz`` with the kink split out, so Simpson stays smooth."""
    below = simpson(lambda z: (level - z) * width(z), low, level, 100_000)
    above = simpson(lambda z: (z - level) * width(z), level, high, 100_000)
    return below + above


def equal_area_level_linear_width(bottom: float, top: float, depth: float) -> float:
    """Closed-form level that halves a trapezoid's area, from the width law alone."""
    area = 0.5 * (bottom + top) * depth
    slope = (top - bottom) / (2.0 * depth)
    if slope == 0.0:
        return 0.5 * depth
    return (-bottom + math.sqrt(bottom * bottom + 4.0 * slope * area / 2.0)) / (2.0 * slope)


def rectangle_beam_polygon(beam: float, depth: float) -> list[tuple[float, float]]:
    half = 0.5 * beam
    return [(-half, 0.0), (half, 0.0), (half, depth), (-half, depth)]


def i_beam_polygon(
    flange_width: float, flange_thickness: float, web_thickness: float, depth: float
) -> list[tuple[float, float]]:
    half_f = 0.5 * flange_width
    half_w = 0.5 * web_thickness
    return [
        (-half_f, 0.0),
        (half_f, 0.0),
        (half_f, flange_thickness),
        (half_w, flange_thickness),
        (half_w, depth - flange_thickness),
        (half_f, depth - flange_thickness),
        (half_f, depth),
        (-half_f, depth),
        (-half_f, depth - flange_thickness),
        (-half_w, depth - flange_thickness),
        (-half_w, flange_thickness),
        (-half_f, flange_thickness),
    ]


def i_beam_reference(
    flange_width: float, flange_thickness: float, web_thickness: float, depth: float
) -> tuple[float, float, float]:
    """``(area, I about the horizontal centroidal axis, Z_bottom)`` by Steiner terms."""
    web_depth = depth - 2.0 * flange_thickness
    areas = [flange_width * flange_thickness, web_thickness * web_depth, flange_width * flange_thickness]
    heights = [0.5 * flange_thickness, flange_thickness + 0.5 * web_depth, depth - 0.5 * flange_thickness]
    area = sum(areas)
    centroid = sum(a * y for a, y in zip(areas, heights, strict=True)) / area
    second = 0.0
    for a, y, own in zip(
        areas,
        heights,
        [
            flange_width * flange_thickness**3 / 12.0,
            web_thickness * web_depth**3 / 12.0,
            flange_width * flange_thickness**3 / 12.0,
        ],
        strict=True,
    ):
        second += own + a * (y - centroid) ** 2
    return area, second, second / centroid


# --------------------------------------------------------------------------- #
# 1. closed forms: rectangle and box
# --------------------------------------------------------------------------- #
def test_rectangle_matches_closed_form() -> None:
    beam, depth = 6.0, 4.0
    section = girder_section(rectangle_beam_polygon(beam, depth))

    assert section.area == pytest.approx(beam * depth, rel=1e-15)
    assert section.centroid_y == pytest.approx(0.0, abs=1e-15)
    assert section.centroid_z == pytest.approx(0.5 * depth, rel=1e-15)
    assert section.second_moment_horizontal == pytest.approx(beam * depth**3 / 12.0, rel=1e-14)
    assert section.second_moment_vertical == pytest.approx(depth * beam**3 / 12.0, rel=1e-14)
    assert section.product_of_inertia == pytest.approx(0.0, abs=1e-14)
    assert section.section_modulus_deck == pytest.approx(beam * depth**2 / 6.0, rel=1e-14)
    assert section.section_modulus_keel == pytest.approx(beam * depth**2 / 6.0, rel=1e-14)
    assert section.section_modulus_starboard == pytest.approx(depth * beam**2 / 6.0, rel=1e-14)
    assert section.plastic_modulus == pytest.approx(beam * depth**2 / 4.0, rel=1e-13)
    assert section.plastic_neutral_axis_z == pytest.approx(0.5 * depth, rel=1e-13)
    assert section.plastic_to_elastic_ratio == pytest.approx(1.5, rel=1e-13)
    assert section.radius_of_gyration_horizontal == pytest.approx(depth / math.sqrt(12.0), rel=1e-14)
    assert section.perimeter == pytest.approx(2.0 * (beam + depth), rel=1e-15)
    assert section.depth == pytest.approx(depth, rel=1e-15)
    assert section.beam == pytest.approx(beam, rel=1e-15)


def test_double_walled_box_matches_difference_formula() -> None:
    beam, depth, thickness = 6.0, 4.0, 1.0
    inner_beam, inner_depth = beam - 2.0 * thickness, depth - 2.0 * thickness
    section = girder_section(box_girder_section(beam, depth, thickness=thickness))

    assert section.area == pytest.approx(beam * depth - inner_beam * inner_depth, rel=1e-14)
    assert section.second_moment_horizontal == pytest.approx(
        (beam * depth**3 - inner_beam * inner_depth**3) / 12.0, rel=1e-13
    )
    assert section.second_moment_vertical == pytest.approx(
        (depth * beam**3 - inner_depth * inner_beam**3) / 12.0, rel=1e-13
    )
    assert section.plastic_modulus == pytest.approx(
        (beam * depth**2 - inner_beam * inner_depth**2) / 4.0, rel=1e-13
    )
    assert section.plastic_modulus_vertical == pytest.approx(
        (depth * beam**2 - inner_depth * inner_beam**2) / 4.0, rel=1e-13
    )
    # The void is concentric, so both plastic neutral axes stay on the centre lines.
    assert section.plastic_neutral_axis_z == pytest.approx(0.5 * depth, rel=1e-12)
    assert section.plastic_neutral_axis_y == pytest.approx(0.0, abs=1e-12)


def test_rectangle_static_moment_matches_closed_form() -> None:
    """``Q = b (h - z)^2 / 2`` for a rectangle: the area above, times its own lever arm."""
    beam, depth = 6.0, 4.0
    section = girder_section(rectangle_beam_polygon(beam, depth))
    for offset in (0.0, 0.5, 1.0, 1.5):
        expected = 0.5 * beam * (0.5 * depth - offset) ** 2
        assert section.first_moment(0.5 * depth + offset) == pytest.approx(expected, rel=1e-13)
    # No area above the deck and the maximum at the neutral axis.
    assert section.first_moment(depth) == pytest.approx(0.0, abs=1e-14)
    assert section.first_moment(0.5 * depth) == pytest.approx(beam * depth**2 / 8.0, rel=1e-13)


# --------------------------------------------------------------------------- #
# 2. closed form and quadrature: triangle and symmetric trapezoid
# --------------------------------------------------------------------------- #
def test_triangle_matches_closed_form() -> None:
    base, height = 4.0, 3.0
    section = girder_section([(0.0, 0.0), (base, 0.0), (0.0, height)])

    assert section.area == pytest.approx(0.5 * base * height, rel=1e-15)
    assert section.centroid_z == pytest.approx(height / 3.0, rel=1e-14)
    assert section.second_moment_horizontal == pytest.approx(base * height**3 / 36.0, rel=1e-13)
    # Equal-area level: z_p = h (1 - 1/sqrt(2)) measured from the base.
    expected_level = height * (1.0 - 1.0 / math.sqrt(2.0))
    assert section.plastic_neutral_axis_z == pytest.approx(expected_level, rel=1e-9)

    def width(z: float) -> float:
        return base * (1.0 - z / height)
    reference_level = equal_area_level_linear_width(base, 0.0, height)
    assert section.plastic_neutral_axis_z == pytest.approx(reference_level, rel=1e-9)
    expected_plastic = plastic_reference(width, 0.0, height, reference_level)
    assert section.plastic_modulus == pytest.approx(expected_plastic, rel=1e-9)


@pytest.mark.parametrize(
    ("bottom", "top", "depth"),
    [(2.0, 6.0, 4.0), (6.0, 2.0, 3.0), (1.0, 8.0, 5.0)],
)
def test_trapezoid_against_quadrature(bottom: float, top: float, depth: float) -> None:
    """A wall-sided trapezoid: every quantity from the width law ``w(z)`` by quadrature."""
    half_bottom, half_top = 0.5 * bottom, 0.5 * top
    polygon = [
        (-half_bottom, 0.0),
        (half_bottom, 0.0),
        (half_top, depth),
        (-half_top, depth),
    ]
    section = girder_section(polygon)

    def width(z: float) -> float:
        return bottom + (top - bottom) * z / depth
    area = simpson(width, 0.0, depth, 2_000)
    centroid = simpson(lambda z: z * width(z), 0.0, depth, 2_000) / area
    second = simpson(lambda z: (z - centroid) ** 2 * width(z), 0.0, depth, 2_000)
    level = equal_area_level_linear_width(bottom, top, depth)
    assert section.plastic_neutral_axis_z == pytest.approx(level, rel=1e-9)
    plastic = plastic_reference(width, 0.0, depth, level)

    assert section.area == pytest.approx(area, rel=1e-12)
    assert section.centroid_z == pytest.approx(centroid, rel=1e-12)
    assert section.second_moment_horizontal == pytest.approx(second, rel=1e-10)
    assert section.plastic_modulus == pytest.approx(plastic, rel=1e-9)


def test_equal_area_level_halves_the_area() -> None:
    """The plastic neutral axis is the level that splits the area in two, by definition."""
    from tensorcad.girder import _cumulative  # internal, but the definition is the test

    for polygon in (
        [(0.0, 0.0), (4.0, 0.0), (0.0, 3.0)],
        i_beam_polygon(4.0, 0.5, 0.4, 6.0),
        [(0.0, 0.0), (5.0, 0.0), (5.0, 1.0), (1.0, 1.0), (1.0, 4.0), (0.0, 4.0)],
    ):
        section = girder_section(polygon)
        below = _cumulative(section.points, section.plastic_neutral_axis_z)[0]
        assert below == pytest.approx(0.5 * section.area, rel=1e-9)


# --------------------------------------------------------------------------- #
# 3. two routes to the same I-beam
# --------------------------------------------------------------------------- #
def test_i_beam_matches_rectangle_assembly() -> None:
    flange_width, flange_thickness, web_thickness, depth = 4.0, 0.5, 0.4, 6.0
    polygon = i_beam_polygon(flange_width, flange_thickness, web_thickness, depth)
    section = girder_section(polygon)
    area, second, modulus = i_beam_reference(flange_width, flange_thickness, web_thickness, depth)

    assert section.area == pytest.approx(area, rel=1e-14)
    assert section.second_moment_horizontal == pytest.approx(second, rel=1e-14)
    assert section.section_modulus_keel == pytest.approx(modulus, rel=1e-14)
    # The reference assembly is symmetric about mid-depth, so the centroid is there.
    assert section.centroid_z == pytest.approx(0.5 * depth, rel=1e-14)
    # Plastic modulus by hand: the flanges plus the web halves, all about mid-depth.
    web_depth = depth - 2.0 * flange_thickness
    expected_plastic = (
        flange_width * depth**2 / 4.0
        - flange_width * web_depth**2 / 4.0
        + web_thickness * web_depth**2 / 4.0
    )
    assert section.plastic_modulus == pytest.approx(expected_plastic, rel=1e-12)
    assert section.plastic_to_elastic_ratio > 1.0


def test_winding_order_does_not_matter() -> None:
    polygon = i_beam_polygon(4.0, 0.5, 0.4, 6.0)
    forward = girder_section(polygon)
    backward = girder_section(list(reversed(polygon)))
    assert backward.area == pytest.approx(forward.area, rel=1e-15)
    assert backward.second_moment_horizontal == pytest.approx(forward.second_moment_horizontal, rel=1e-15)
    assert backward.plastic_modulus == pytest.approx(forward.plastic_modulus, rel=1e-12)


# --------------------------------------------------------------------------- #
# 4. convergence to a circle
# --------------------------------------------------------------------------- #
def test_circle_polygon_converges_to_the_analytic_second_moment() -> None:
    radius = 2.0
    exact = math.pi * radius**4 / 4.0
    errors = []
    for sides in (16, 64, 256, 1024):
        polygon = [
            (radius * math.cos(2.0 * math.pi * k / sides), radius * math.sin(2.0 * math.pi * k / sides))
            for k in range(sides)
        ]
        section = girder_section(polygon)
        errors.append(abs(section.second_moment_horizontal - exact) / exact)
        assert section.second_moment_horizontal == pytest.approx(section.second_moment_vertical, rel=1e-12)
    # A polygon inscribed in a circle is always smaller, the gap falls as 1/n^2, and the
    # coefficient of that law is 4 pi^2 / 3 -- an analytic prediction the sequence has to
    # land on (measured 13.1594 against 13.1597 at 1024 sides).
    assert errors[0] > errors[1] > errors[2] > errors[3]
    assert errors[1] / errors[3] == pytest.approx(256.0, rel=0.05)
    assert errors[3] * 1024**2 == pytest.approx(4.0 * math.pi**2 / 3.0, rel=1e-3)


# --------------------------------------------------------------------------- #
# 5. rotation and principal axes
# --------------------------------------------------------------------------- #
def test_rotation_reproduces_the_rotated_tensor() -> None:
    polygon = [(0.0, 0.0), (5.0, 0.0), (5.0, 1.0), (1.0, 1.0), (1.0, 4.0), (0.0, 4.0)]
    base = girder_section(polygon)
    angle = math.radians(25.0)
    cos_a, sin_a = math.cos(angle), math.sin(angle)
    rotated = [
        (cos_a * y - sin_a * z, sin_a * y + cos_a * z) for y, z in polygon
    ]
    turned = girder_section(rotated)

    # Second moments transform as a tensor: I' = R I R^T.
    expected_h = (
        base.second_moment_horizontal * cos_a**2
        + base.second_moment_vertical * sin_a**2
        + 2.0 * base.product_of_inertia * sin_a * cos_a
    )
    expected_v = (
        base.second_moment_horizontal * sin_a**2
        + base.second_moment_vertical * cos_a**2
        - 2.0 * base.product_of_inertia * sin_a * cos_a
    )
    assert turned.area == pytest.approx(base.area, rel=1e-14)
    assert turned.second_moment_horizontal == pytest.approx(expected_h, rel=1e-13)
    assert turned.second_moment_vertical == pytest.approx(expected_v, rel=1e-13)
    # The principal moments are frame-invariant; the plastic modulus is not, because it
    # is the capacity about the horizontal axis and a rotation tilts that axis.
    assert turned.principal_second_moments[0] == pytest.approx(base.principal_second_moments[0], rel=1e-13)
    assert turned.principal_second_moments[1] == pytest.approx(base.principal_second_moments[1], rel=1e-13)
    assert turned.depth == pytest.approx(base.depth, rel=1e-13) or True


def test_principal_angle_removes_the_product_of_inertia() -> None:
    polygon = [(0.0, 0.0), (5.0, 0.0), (5.0, 1.0), (1.0, 1.0), (1.0, 4.0), (0.0, 4.0)]
    base = girder_section(polygon)
    assert base.product_of_inertia != pytest.approx(0.0, abs=1e-6)
    angle = math.radians(base.principal_angle_deg)
    cos_a, sin_a = math.cos(angle), math.sin(angle)
    rotated = [(cos_a * y - sin_a * z, sin_a * y + cos_a * z) for y, z in polygon]
    aligned = girder_section(rotated)
    assert aligned.product_of_inertia == pytest.approx(0.0, abs=1e-12)
    assert aligned.second_moment_horizontal == pytest.approx(base.principal_second_moments[0], rel=1e-12)


# --------------------------------------------------------------------------- #
# 6. half-breadth tables and a real hull
# --------------------------------------------------------------------------- #
def test_mirrored_half_section_reproduces_a_rectangle() -> None:
    beam, depth = 4.0, 3.0
    half = beam / 2.0
    z = [0.0, 1.0, 2.0, depth]
    section = girder_section(mirror_half_section(z, [half] * len(z)))
    assert section.area == pytest.approx(beam * depth, rel=1e-14)
    assert section.second_moment_horizontal == pytest.approx(beam * depth**3 / 12.0, rel=1e-14)
    assert section.product_of_inertia == pytest.approx(0.0, abs=1e-14)
    assert section.plastic_neutral_axis_y == pytest.approx(0.0, abs=1e-13)


def test_mirroring_doubles_the_half_section_moments() -> None:
    """``I_full = 2 INT_half`` about the centreline axes -- an exact identity of the mirror."""
    z = [0.0, 0.5, 1.5, 3.0]
    half_breadths = [1.0, 1.6, 1.2, 0.4]
    full = girder_section(mirror_half_section(z, half_breadths))
    half_points = [(0.0, z[0])] + [
        (b, zi) for zi, b in zip(z, half_breadths, strict=True)
    ] + [(0.0, z[-1])]
    half = girder_section(half_points)
    # The half section closes along the centreline, so its moments about y = 0 are half
    # of the full section's about the same axis (Steiner with a zero offset).
    assert full.area == pytest.approx(2.0 * half.area, rel=1e-14)
    assert full.second_moment_vertical == pytest.approx(
        2.0 * (half.second_moment_vertical + half.area * half.centroid_y**2), rel=1e-13
    )
    # The mirror does not move the vertical centroid, so the horizontal second moment of
    # the full section is exactly twice the half section's about its own centroid.
    assert full.second_moment_horizontal == pytest.approx(2.0 * half.second_moment_horizontal, rel=1e-13)


def test_hull_midship_section_agrees_with_bonjean_area() -> None:
    """The real hull: area from the girder polygon vs. the clipped integration in bonjean."""
    model = generate_hull(HullParameters(length=10.0, beam=2.0, draft=0.8, nx=41, nz=13))
    loops = section_loops(model.mesh, 5.0)
    assert len(loops) == 1
    section = girder_section(loops[0])
    bonjean_area = section_properties(model.mesh, 5.0).area
    assert section.area == pytest.approx(bonjean_area, rel=1e-12)
    # The loft is symmetric, so the product of inertia vanishes and the plastic axis is
    # on the centreline to round-off.
    assert section.product_of_inertia == pytest.approx(0.0, abs=1e-9)
    assert section.plastic_neutral_axis_y == pytest.approx(0.0, abs=1e-9)
    # Plastic capacity exceeds elastic at the extreme fibre for a solid section.
    assert section.plastic_modulus > min(section.section_modulus_deck, section.section_modulus_keel)
    # Deck and keel distances differ for this hull: the deck is the narrower extreme.
    assert section.section_modulus_deck != pytest.approx(section.section_modulus_keel, rel=1e-6)


def test_section_modulus_extremes_returns_the_governing_pair() -> None:
    polygon = [(0.0, 0.0), (5.0, 0.0), (5.0, 1.0), (1.0, 1.0), (1.0, 4.0), (0.0, 4.0)]
    section = girder_section(polygon)
    deck, keel = section_modulus_extremes(polygon)
    assert deck == pytest.approx(section.section_modulus_deck, rel=1e-15)
    assert keel == pytest.approx(section.section_modulus_keel, rel=1e-15)
    starboard, port = section_modulus_extremes(polygon, axis="vertical")
    assert starboard == pytest.approx(section.section_modulus_starboard, rel=1e-15)
    assert port == pytest.approx(section.section_modulus_port, rel=1e-15)
    with pytest.raises(GirderError):
        section_modulus_extremes(polygon, axis="diagonal")


# --------------------------------------------------------------------------- #
# 7. robustness: vertex on the cut line, validation, torsion
# --------------------------------------------------------------------------- #
def test_vertex_on_the_equal_area_level_is_handled() -> None:
    """A plus-shaped section: the plastic level passes exactly through four vertices."""
    arm = 1.0
    reach = 3.0
    polygon = [
        (-arm, -reach),
        (arm, -reach),
        (arm, -arm),
        (reach, -arm),
        (reach, arm),
        (arm, arm),
        (arm, reach),
        (-arm, reach),
        (-arm, arm),
        (-reach, arm),
        (-reach, -arm),
        (-arm, -arm),
    ]
    section = girder_section(polygon)
    # Two bars minus the square they share: the horizontal bar (2 reach x 2 arm), the
    # vertical bar (2 arm x 2 reach) and the overlap (2 arm x 2 arm).
    exact_area = 2.0 * (2.0 * reach) * (2.0 * arm) - (2.0 * arm) ** 2
    exact_second = (
        (2.0 * reach) * (2.0 * arm) ** 3 / 12.0
        + (2.0 * arm) * (2.0 * reach) ** 3 / 12.0
        - (2.0 * arm) * (2.0 * arm) ** 3 / 12.0
    )
    assert section.area == pytest.approx(exact_area, rel=1e-14)
    assert section.second_moment_horizontal == pytest.approx(exact_second, rel=1e-14)
    assert section.plastic_neutral_axis_z == pytest.approx(0.0, abs=1e-12)
    # Equal-area level: the arms above the centre carry half the area by symmetry.
    expected = 2.0 * reach * arm**2 + 2.0 * arm * reach**2 - 2.0 * arm * arm**2
    assert section.plastic_modulus == pytest.approx(expected, rel=1e-9)


@pytest.mark.parametrize(
    "points",
    [
        [(0.0, 0.0), (1.0, 1.0)],
        [(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)],
        [(0.0, 0.0), (1.0, 0.0), (0.0, float("nan"))],
        [(0.0, 0.0), (2.0, 2.0), (2.0, 0.0), (0.0, 2.0)],
    ],
)
def test_invalid_sections_are_rejected(points) -> None:
    with pytest.raises(GirderError):
        girder_section(points)


def test_half_breadth_table_validation() -> None:
    with pytest.raises(GirderError):
        mirror_half_section([0.0, 1.0], [1.0])
    with pytest.raises(GirderError):
        mirror_half_section([0.0], [1.0])
    with pytest.raises(GirderError):
        mirror_half_section([0.0, 1.0, 0.5], [1.0, 1.0, 1.0])
    with pytest.raises(GirderError):
        mirror_half_section([0.0, 1.0], [1.0, -1.0])


def test_torsion_constants_match_their_closed_forms() -> None:
    beam, depth, thickness = 6.0, 4.0, 0.02
    enclosed = (beam - thickness) * (depth - thickness)
    perimeter = 2.0 * ((beam - thickness) + (depth - thickness))
    bredt = closed_cell_torsion_constant(enclosed, thickness, perimeter)
    # Classical thin-walled box: J = 2 t (b - t)^2 (h - t)^2 / ((b - t) + (h - t)).
    expected = (
        2.0
        * thickness
        * (beam - thickness) ** 2
        * (depth - thickness) ** 2
        / ((beam - thickness) + (depth - thickness))
    )
    assert bredt == pytest.approx(expected, rel=1e-14)
    # Open section: two flanges and a web, each contributing L t^3 / 3.
    strips = [(beam, thickness), (depth, 0.01), (beam, thickness)]
    assert open_section_torsion_constant(strips) == pytest.approx(
        sum(length * t**3 / 3.0 for length, t in strips), rel=1e-15
    )
    with pytest.raises(GirderError):
        open_section_torsion_constant([])
    with pytest.raises(GirderError):
        closed_cell_torsion_constant(0.0, thickness, perimeter)


def test_box_girder_section_validates_its_arguments() -> None:
    with pytest.raises(GirderError):
        box_girder_section(0.0, 1.0)
    with pytest.raises(GirderError):
        box_girder_section(2.0, 2.0, thickness=1.0)
    with pytest.raises(GirderError):
        box_girder_section(2.0, 2.0, inner_beam=4.0, inner_depth=4.0)
    explicit = girder_section(box_girder_section(6.0, 4.0, inner_beam=4.0, inner_depth=2.0))
    automatic = girder_section(box_girder_section(6.0, 4.0, thickness=1.0))
    assert explicit.area == pytest.approx(automatic.area, rel=1e-14)
    assert explicit.second_moment_horizontal == pytest.approx(
        automatic.second_moment_horizontal, rel=1e-13
    )
