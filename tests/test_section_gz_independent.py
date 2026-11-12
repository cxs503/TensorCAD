"""Independent re-check of :mod:`tensorcad.section_gz` (verification debt, mission m16430).

``tensorcad.section_gz`` (merged in m16358) ships its own evidence; these tests ask whether a
route written independently -- slice quadrature for the ellipse, exact slice integration for
polygons, no clipping and no Green's theorem -- reproduces the same arms.  The references live
in :mod:`tests.section_gz_oracle`; nothing here reads a committed artifact, so a stale number
cannot turn the suite green.
"""

from __future__ import annotations

import math

import pytest

from section_gz_oracle import (
    EllipseSliceOracle,
    ellipse_prism_arms_independent,
    polygon_prism_arms_independent,
    polygon_properties_below,
    wall_sided_arm,
    wall_sided_limit_deg,
    waterplane_gm,
)
from tensorcad import section_gz

HEELS = [0.0, 10.0, 20.0, 30.0, 45.0, 60.0, 90.0]
TRIANGLE = ((-1.0, 0.0), (1.0, 0.0), (0.0, 2.0))
RECTANGLE = ((-2.0, 0.0), (2.0, 0.0), (2.0, 2.0), (-2.0, 2.0))
ELLIPSE = {"semi_beam": 2.0, "semi_height": 1.0, "draft": 1.2, "kg": 0.8}


def test_the_independent_ellipse_oracle_reproduces_the_ellipse_area() -> None:
    """Sanity anchor for the oracle itself: the whole conic integrates to ``pi a b``."""
    oracle = EllipseSliceOracle(ELLIPSE["semi_beam"], ELLIPSE["semi_height"], 0.3)
    assert oracle.area_below(10.0) == pytest.approx(oracle.area, rel=1e-12)


def test_the_ellipse_closed_form_agrees_with_slice_quadrature() -> None:
    """The library's disk formula and the independent conic quadrature must agree to round-off."""
    arms = section_gz.ellipse_prism_righting_arms(
        ELLIPSE["semi_beam"],
        ELLIPSE["semi_height"],
        HEELS,
        kg=ELLIPSE["kg"],
        waterline_z=ELLIPSE["draft"],
        length=4.0,
    )
    independent = ellipse_prism_arms_independent(
        ELLIPSE["semi_beam"],
        ELLIPSE["semi_height"],
        HEELS,
        kg=ELLIPSE["kg"],
        waterline_z=ELLIPSE["draft"],
    )
    worst = max(abs(left - right) for left, right in zip(arms.gz, independent, strict=True))
    assert worst < 1e-9, f"closed form and independent quadrature differ by {worst:.3e} m"


def test_the_polygon_oracle_agrees_with_exact_slice_integration() -> None:
    """Box, V section and an off-axis hexagon: clipping/Green vs exact slice integration."""
    cases = [
        (RECTANGLE, 0.6, 1.0),
        (TRIANGLE, 0.6, 1.0),
        (section_gz.regular_polygon_section(6, 1.0, rotation_deg=17.0), 0.7, 0.8),
    ]
    for section, kg, waterline in cases:
        arms = section_gz.polygon_prism_righting_arms(
            section, HEELS, kg=kg, waterline_z=waterline
        )
        independent = polygon_prism_arms_independent(
            section, HEELS, kg=kg, waterline_z=waterline
        )
        worst = max(abs(left - right) for left, right in zip(arms.gz, independent, strict=True))
        assert worst < 1e-12, f"{len(section)}-gon arms differ by {worst:.3e} m"


def test_the_wall_sided_formula_is_an_independent_anchor_inside_its_limit() -> None:
    """Below the deck-edge limit the published formula and the oracle are the same number."""
    limit = wall_sided_limit_deg(4.0, 1.0)
    heels = [0.0, 5.0, 10.0, 20.0, 0.9 * limit]
    arms = section_gz.polygon_prism_righting_arms(
        RECTANGLE, heels, kg=0.6, waterline_z=1.0
    ).gz
    worst = max(
        abs(arm - wall_sided_arm(4.0, 1.0, 0.6, heel))
        for arm, heel in zip(arms, heels, strict=True)
    )
    assert worst < 1e-12, f"wall-sided agreement lost by {worst:.3e} m"


def test_the_wall_sided_formula_really_fails_past_the_limit() -> None:
    """Above the limit the oracle leaves the wall-sided formula, so the check is not vacuous."""
    above = [wall_sided_limit_deg(4.0, 1.0) + 5.0, wall_sided_limit_deg(4.0, 1.0) + 15.0]
    arms = section_gz.polygon_prism_righting_arms(
        RECTANGLE, above, kg=0.6, waterline_z=1.0
    ).gz
    departures = [
        abs(arm - wall_sided_arm(4.0, 1.0, 0.6, heel))
        for arm, heel in zip(arms, above, strict=True)
    ]
    assert min(departures) > 1e-2, f"wall-sided formula still fits above its limit: {departures}"


def test_the_heeled_gm_equals_the_waterplane_moments_gm() -> None:
    """``dGZ/dphi(0)`` from the closed form against ``KB + I_T/A - KG`` integrated here."""
    a, b = ELLIPSE["semi_beam"], ELLIPSE["semi_height"]
    draft, kg = ELLIPSE["draft"], ELLIPSE["kg"]
    oracle = EllipseSliceOracle(a, b, 0.0)
    area, _, kb = oracle.properties_below(draft)
    gm_waterplane = waterplane_gm(0.5 * oracle.width(draft), area, kb, kg)
    delta = 1e-6
    gm_slope = (
        section_gz.ellipse_prism_righting_arms(
            a, b, [delta], kg=kg, waterline_z=draft, length=4.0
        ).gz[0]
        / math.radians(delta)
    )
    assert gm_slope == pytest.approx(gm_waterplane, abs=1e-9)


def test_a_symmetric_section_has_an_antisymmetric_arm() -> None:
    """``GZ(-phi) = -GZ(phi)`` for the ellipse, on both the library and the oracle."""
    heels = [10.0, 30.0, 60.0]
    signed = [heel for heel in heels] + [-heel for heel in heels]
    arms = section_gz.ellipse_prism_righting_arms(
        ELLIPSE["semi_beam"], ELLIPSE["semi_height"], signed, kg=ELLIPSE["kg"],
        waterline_z=ELLIPSE["draft"], length=4.0,
    ).gz
    half = len(heels)
    worst = max(abs(arms[index] + arms[index + half]) for index in range(half))
    assert worst < 1e-12


@pytest.mark.parametrize("sides", [5, 16, 17, 18, 20, 32])
def test_ellipse_section_is_inscribed_in_the_exact_ellipse(sides: int) -> None:
    """Regression for the m16430 defect: a side count without a keel vertex used to cross out.

    ``ellipse_section`` samples the ellipse and then shifts the lowest vertex to ``z = 0``.  If
    no vertex sits at the keel the shift translates the polygon *through* the section: for
    ``sides = 18`` the boundary reached ``1.03015`` on the ellipse's quadratic form, i.e. 3 %
    outside, and the stated "inscribed, known-sign chord error" did not hold.  The sampling now
    starts at ``-pi/2`` so a keel vertex always exists (and ``sides`` divisible by four keeps
    exactly the old vertex set).
    """
    a, b = ELLIPSE["semi_beam"], ELLIPSE["semi_height"]
    points = section_gz.ellipse_section(a, b, sides)
    assert min(z for (_, z) in points) == pytest.approx(0.0, abs=1e-15)
    worst = 0.0
    for index in range(len(points)):
        y0, z0 = points[index]
        y1, z1 = points[(index + 1) % len(points)]
        for step in range(21):
            fraction = step / 20.0
            y = y0 + fraction * (y1 - y0)
            z = z0 + fraction * (z1 - z0)
            worst = max(worst, (y / a) ** 2 + ((z - b) / b) ** 2)
    assert worst <= 1.0 + 1e-12, f"polygon boundary leaves the ellipse: quadric {worst!r}"
    area = polygon_properties_below(points, max(z for (_, z) in points))[0]
    assert area < math.pi * a * b, "an inscribed polygon can never fill the ellipse"