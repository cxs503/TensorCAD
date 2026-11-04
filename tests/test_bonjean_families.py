"""The two Bonjean references that did not exist, and the meshes they check.

Covered here:

* :func:`tensorcad.bonjean.analytic_section_area` is bit-for-bit the shipped
  expression when no exponents are given, and splits at midship when they are
  (:func:`tensorcad.bonjean.taper_exponent_at`);
* :func:`tensorcad.wigley.analytic_section_area` against a Gauss-Legendre
  quadrature of the surface itself (:func:`tensorcad.wigley.half_breadth`),
  which is an independent path, and against the three limits the family has
  (no area below the keel, the maximum section at the design waterline, a full
  section above it);
* both references against the meshes the generators loft -- including the
  negative control: the asymmetric hull must *fail* a symmetric reference, or
  the agreement proves nothing about the fore/aft split.
"""
from __future__ import annotations

import math

import pytest

from tensorcad.bonjean import analytic_section_area, section_properties, taper_exponent_at
from tensorcad.hull import BREADTH_EXPONENT, HullParameters, build_mesh, taper_integral
from tensorcad.wigley import (
    WigleyParameters,
    analytic_maximum_section_area,
    analytic_section_area as wigley_section_area,
    build_mesh as build_wigley_mesh,
    half_breadth,
)

LENGTH, BEAM, DRAFT = 10.0, 2.0, 0.8
FOR, AFT = 1.4, 0.7
WIGLEY = {"length": 10.0, "beam": 1.0, "draft": 0.625}


def hull_mesh(nx: int = 41, nz: int = 25, fore: float = FOR, aft: float = AFT):
    parameters = HullParameters(
        length=LENGTH, beam=BEAM, draft=DRAFT, nx=nx, nz=nz,
        fore_exponent=fore, aft_exponent=aft,
    )
    mesh, _ = build_mesh(parameters)
    return parameters, mesh


def hull_reference(station: float, waterline: float, fore: float = FOR, aft: float = AFT) -> float:
    return analytic_section_area(
        LENGTH, BEAM, DRAFT, station, waterline, fore_exponent=fore, aft_exponent=aft
    )


def wigley_reference(station: float, waterline: float, parameters: WigleyParameters | None = None) -> float:
    p = parameters or WigleyParameters(**WIGLEY)
    return wigley_section_area(
        p.length, p.beam, p.draft, station, waterline, p.longitudinal_exponent, p.vertical_exponent
    )


def gauss_legendre(function, lower: float, upper: float, order: int = 400) -> float:
    """Composite Gauss-Legendre: an integration path independent of the algebra."""
    total = 0.0
    for node in range(order):
        a = lower + (upper - lower) * node / order
        b = lower + (upper - lower) * (node + 1) / order
        mid, half = 0.5 * (a + b), 0.5 * (b - a)
        for point, weight in ((0.0, 0.8888888888888888),
                              (-0.7745966692414834, 0.5555555555555556),
                              (0.7745966692414834, 0.5555555555555556)):
            total += weight * function(mid + half * point) * half
    return total


# ----------------------------------------------------- the shipped expression --
def test_defaults_reproduce_the_shipped_formula_bit_for_bit():
    for station in (0.0, 1.0, 2.5, 5.0, 7.5, 9.0, 10.0):
        for waterline in (-DRAFT, -0.6, -0.2, 0.0):
            zeta = waterline + DRAFT
            expected = 0.0
            if zeta > 0.0:
                expected = (
                    (2.0 / 3.0)
                    * BEAM
                    * math.pow(math.sin(math.pi * station / LENGTH), BREADTH_EXPONENT)
                    * math.pow(zeta, 1.5)
                    / math.sqrt(DRAFT)
                )
            assert analytic_section_area(LENGTH, BEAM, DRAFT, station, waterline) == expected


def test_explicit_shipped_exponents_equal_the_defaults():
    station, waterline = 3.0, -0.3
    default = analytic_section_area(LENGTH, BEAM, DRAFT, station, waterline)
    explicit = analytic_section_area(
        LENGTH, BEAM, DRAFT, station, waterline,
        fore_exponent=BREADTH_EXPONENT, aft_exponent=BREADTH_EXPONENT,
    )
    assert default == explicit


def test_taper_exponent_splits_at_midship():
    assert taper_exponent_at(LENGTH, 0.0, FOR, AFT) == AFT
    assert taper_exponent_at(LENGTH, LENGTH / 2 - 1e-9, FOR, AFT) == AFT
    # Midship belongs to the fore half, matching hull._breadth_exponent, where
    # sin(pi/2) = 1 makes the choice immaterial.
    assert taper_exponent_at(LENGTH, LENGTH / 2, FOR, AFT) == FOR
    assert taper_exponent_at(LENGTH, LENGTH, FOR, AFT) == FOR


def test_taper_exponent_rejects_bad_arguments():
    with pytest.raises(ValueError):
        taper_exponent_at(0.0, 0.0)
    with pytest.raises(ValueError):
        taper_exponent_at(LENGTH, 5.0, fore_exponent=0.0)
    with pytest.raises(ValueError):
        taper_exponent_at(LENGTH, 5.0, aft_exponent=float("nan"))


# ------------------------------------------------------- the split closed form --
@pytest.mark.parametrize("station", [1.0, 2.5, 4.9, 5.0, 6.0, 8.5, 9.9])
@pytest.mark.parametrize("waterline", [-0.7, -0.4, -0.1, 0.0])
def test_split_closed_form_matches_quadrature_of_the_surface(station, waterline):
    """Integrate the half-breadth the generator lofts, independently of the algebra."""
    exponent = taper_exponent_at(LENGTH, station, FOR, AFT)
    taper = math.sin(math.pi * station / LENGTH) ** exponent

    # ``zeta = H s^2``, ``dzeta = 2 H s ds``: the sqrt taper becomes analytic in s,
    # so a three-point rule is exact to rounding instead of converging like h^1.5.
    def integrand(s: float) -> float:
        height = waterline + DRAFT
        return 2.0 * BEAM * taper * s * math.sqrt(height * s * s / DRAFT)

    height = waterline + DRAFT
    quadrature = height * gauss_legendre(integrand, 0.0, 1.0, order=20) if height > 0.0 else 0.0
    closed = hull_reference(station, waterline)
    assert closed == pytest.approx(quadrature, rel=1e-12, abs=1e-15)


def test_split_closed_form_is_the_shipped_one_scaled_by_the_taper_ratio():
    station = 8.5
    shipped = analytic_section_area(LENGTH, BEAM, DRAFT, station, -0.3)
    split = hull_reference(station, -0.3)
    ratio = math.sin(math.pi * station / LENGTH) ** (FOR - BREADTH_EXPONENT)
    assert split == pytest.approx(shipped * ratio, rel=1e-12)


def test_split_closed_form_is_zero_below_the_keel_and_positive_above():
    assert hull_reference(3.0, -DRAFT) == 0.0
    assert hull_reference(3.0, -DRAFT - 0.5) == 0.0
    assert hull_reference(3.0, -0.4) > 0.0
    # The closed form keeps growing above the design waterline (the family has no
    # deck); the *mesh* stops there because the loft is capped at z = 0. Both
    # sides of that are asserted, the reference here and the cap in
    # ``test_mesh_area_stops_at_the_deck``.
    assert hull_reference(3.0, 0.4) > hull_reference(3.0, 0.0)


def test_hull_volume_is_the_integral_of_the_section_areas():
    """The taper pair and the section reference agree on the displacement.

    ``analytic_volume`` integrates the taper through :func:`taper_integral_pair`;
    integrating the section reference along ``x`` must return the same number,
    which ties the new station-resolved form to the published whole-hull one.
    """
    volume = gauss_legendre(lambda x: hull_reference(x, 0.0), 0.0, LENGTH, order=400)
    from tensorcad.hull import analytic_volume

    expected = analytic_volume(LENGTH, BEAM, DRAFT, fore_exponent=FOR, aft_exponent=AFT)
    # ``sin(pi x/L)^p`` has a cusp at both ends, so the composite rule converges
    # slowly there (1.4e-7 at order 400, 9e-9 at 2000): the tolerance is the
    # quadrature's, not the identity's.
    assert volume == pytest.approx(expected, rel=1e-6)
    assert expected > 0.0
    assert taper_integral(FOR) != taper_integral(AFT)


# -------------------------------------------------------- the Wigley reference --
@pytest.mark.parametrize("station", [0.0, 1.0, 3.0, 5.0, 7.0, 9.5, 10.0])
@pytest.mark.parametrize("waterline", [-0.625, -0.5, -0.25, -0.05, 0.0])
def test_wigley_closed_form_matches_quadrature_of_half_breadth(station, waterline):
    p = WigleyParameters(**WIGLEY)

    def integrand(z: float) -> float:
        # ``half_breadth`` takes depth below the design waterline: at height z
        # above the keel that depth is ``T + z``.
        return 2.0 * half_breadth(p, station, p.draft + z)

    height = min(max(waterline + p.draft, 0.0), p.draft)
    quadrature = gauss_legendre(integrand, -height, 0.0) if height > 0.0 else 0.0
    assert wigley_reference(station, waterline) == pytest.approx(quadrature, rel=1e-11, abs=1e-15)


def test_wigley_limits():
    p = WigleyParameters(**WIGLEY)
    assert wigley_section_area(p.length, p.beam, p.draft, 5.0, -p.draft) == 0.0
    assert wigley_section_area(p.length, p.beam, p.draft, 5.0, -2 * p.draft) == 0.0
    # Midship at the design waterline is the maximum section the hydrostatics
    # module already published.
    assert wigley_section_area(p.length, p.beam, p.draft, 5.0, 0.0) == pytest.approx(
        analytic_maximum_section_area(p.beam, p.draft), rel=1e-15
    )
    # The deck cap is flat: nothing above the design waterline.
    assert wigley_section_area(p.length, p.beam, p.draft, 5.0, 0.2) == wigley_section_area(
        p.length, p.beam, p.draft, 5.0, 0.0
    )
    # Ends close on the centreline.
    assert wigley_section_area(p.length, p.beam, p.draft, 0.0, 0.0) == 0.0


def test_wigley_closed_form_is_monotone_in_the_waterline():
    p = WigleyParameters(**WIGLEY)
    levels = [-p.draft + p.draft * index / 20 for index in range(21)]
    areas = [wigley_section_area(p.length, p.beam, p.draft, 3.0, z) for z in levels]
    assert all(later > earlier for earlier, later in zip(areas, areas[1:]))


def test_wigley_rejects_bad_arguments():
    with pytest.raises(ValueError):
        wigley_section_area(0.0, 1.0, 0.625, 0.0, 0.0)
    with pytest.raises(ValueError):
        wigley_section_area(10.0, 1.0, 0.625, 11.0, 0.0)
    with pytest.raises(ValueError):
        wigley_section_area(10.0, 1.0, 0.625, 5.0, 0.0, vertical_exponent=0.0)
    with pytest.raises(ValueError):
        wigley_section_area(10.0, 1.0, 0.625, 5.0, float("inf"))


# ------------------------------------------------------------- against meshes --
def test_wigley_mesh_reproduces_the_closed_form():
    parameters = WigleyParameters(nx=81, nz=13, **WIGLEY)
    mesh, _ = build_wigley_mesh(parameters)
    worst = 0.0
    for station in (1.0, 2.5, 5.0, 7.5, 9.0):
        for waterline in (-0.5, -0.3125, -0.15625, -0.05):
            measured = section_properties(mesh, station, waterline).area
            worst = max(worst, abs(measured / wigley_reference(station, waterline, parameters) - 1.0))
    assert worst < 0.01


def test_asymmetric_mesh_reproduces_the_split_closed_form():
    _, mesh = hull_mesh(81, 25)
    stations = tuple(fraction * LENGTH for fraction in (0.1, 0.25, 0.4, 0.6, 0.75, 0.9))
    worst = 0.0
    for station in stations:
        for waterline in (-0.4, -0.3, -0.2, -0.1, -0.02):
            measured = section_properties(mesh, station, waterline).area
            worst = max(worst, abs(measured / hull_reference(station, waterline) - 1.0))
    assert worst < 0.01, worst


def test_symmetric_reference_fails_on_the_asymmetric_hull():
    """Negative control: ignoring the fore/aft split must be visible on the mesh."""
    _, mesh = hull_mesh(81, 25)
    station = 0.9 * LENGTH
    measured = section_properties(mesh, station, -0.2).area
    symmetric = analytic_section_area(LENGTH, BEAM, DRAFT, station, -0.2)
    split = hull_reference(station, -0.2)
    assert abs(measured / symmetric - 1.0) > 0.1
    assert abs(measured / split - 1.0) < 0.02


def test_swapping_the_exponents_mirrors_the_curves():
    _, forward = hull_mesh(41, 13, fore=FOR, aft=AFT)
    _, swapped = hull_mesh(41, 13, fore=AFT, aft=FOR)
    for fraction in (0.15, 0.3, 0.45, 0.7, 0.85):
        a = section_properties(forward, fraction * LENGTH, -0.2).area
        b = section_properties(swapped, (1.0 - fraction) * LENGTH, -0.2).area
        assert a == pytest.approx(b, rel=1e-6)


def test_refinement_reduces_the_error_quadratically_in_the_waterline_count():
    errors = []
    for nz in (13, 25):
        _, mesh = hull_mesh(41, nz)
        worst = max(
            abs(section_properties(mesh, fraction * LENGTH, -0.2).area / hull_reference(fraction * LENGTH, -0.2) - 1.0)
            for fraction in (0.1, 0.3, 0.5, 0.7, 0.9)
        )
        errors.append(worst)
    assert errors[1] < errors[0]
    assert errors[0] / errors[1] > 2.0


def test_midship_ring_is_shared_by_both_halves():
    """The two exponents meet at midship: the ring there is the same geometry.

    The *cut* just forward of midship is not compared: ``section_properties``
    nudges the plane by ``1e-9`` and that side is genuinely the fore half, so the
    lofted vertices a hair forward of midship differ once the exponents differ.
    The ring itself has ``sin(pi/2) = 1`` under any exponent and must be identical.
    """
    _, symmetric = hull_mesh(41, 13, fore=BREADTH_EXPONENT, aft=BREADTH_EXPONENT)
    _, split = hull_mesh(41, 13)
    a = [y for x, y, z in symmetric.vertices if x == LENGTH / 2]
    b = [y for x, y, z in split.vertices if x == LENGTH / 2]
    assert a and a == b
    assert max(abs(y) for y in a) == pytest.approx(BEAM / 2, rel=1e-15)