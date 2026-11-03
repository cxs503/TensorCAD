"""Tests for the closed-form hydrostatics of the parametric hull family.

The mesh-integrated hydrostatics in :mod:`tensorcad.hydrostatics` were only ever
pinned against a box, whose tetrahedra integrate exactly. The generated hull is
the first case with a curved surface, a keel singularity and a waterline cap, and
it had no reference value beyond its own refined mesh. These tests supply one: the
analytic surface the generator samples is separable, so its volume, waterplane
area, midship section area, LCB, LCF, VCB and form coefficients all have closed
forms, and the wetted area has a quadrature reference.

Two directions are tested:

* the closed forms themselves, against an *independent* composite-Simpson
  quadrature of the same integrands (a formula checked only against itself proves
  nothing);
* the mesh integrals, against the closed forms, at four resolutions, with the sign
  of the error and the convergence direction pinned as well as its magnitude.
"""
from __future__ import annotations

import math

import pytest

from tensorcad import HullParameters, compute_hydrostatics, generate_hull
from tensorcad.hull import (
    BREADTH_EXPONENT,
    TAPER_INTEGRAL,
    analytic_hydrostatics,
    analytic_longitudinal_centre,
    analytic_maximum_section_area,
    analytic_relative_errors,
    analytic_vertical_centre_buoyancy,
    analytic_volume,
    analytic_waterplane_area,
    analytic_wetted_surface_area,
)


def simpson(function, lower: float, upper: float, panels: int) -> float:
    """Composite Simpson on ``panels`` (even) sub-intervals."""
    if panels % 2:
        panels += 1
    step = (upper - lower) / panels
    total = function(lower) + function(upper)
    for index in range(1, panels):
        weight = 4.0 if index % 2 else 2.0
        total += weight * function(lower + index * step)
    return total * step / 3.0


def simpson_2d(function, panels: int) -> float:
    """Composite Simpson product rule on the unit square."""
    nodes, weights = _simpson_nodes(panels)
    total = 0.0
    for tau, w_tau in zip(nodes, weights, strict=False):
        for u, w_u in zip(nodes, weights, strict=False):
            total += w_tau * w_u * function(tau, u)
    return total


def _simpson_nodes(panels: int) -> tuple[list[float], list[float]]:
    if panels % 2:
        panels += 1
    step = 1.0 / panels
    nodes = [index * step for index in range(panels + 1)]
    weights = [step / 3.0 * (1.0 if index in (0, panels) else (4.0 if index % 2 else 2.0))
               for index in range(panels + 1)]
    return nodes, weights


# --------------------------------------------------------------------------- #
# The closed forms, against an independent quadrature.
# --------------------------------------------------------------------------- #

def taper_integral_by_simpson(panels: int = 4000) -> float:
    """``int_0^1 sin(pi*x)^0.7 dx`` by Simpson after removing the tip singularity.

    ``sin(pi*x)^0.7`` behaves like ``x^0.7`` at the ends, so a Simpson rule in
    ``x`` itself converges only algebraically and 20 000 panels still leave an
    ``1e-8`` error. Substituting ``x = (1 - cos(theta))/2`` makes the integrand
    smooth in ``theta`` and the same rule then reaches ``1e-13``.
    """
    def integrand(theta: float) -> float:
        return 0.5 * math.sin(theta) * math.sin(0.5 * math.pi * (1.0 - math.cos(theta))) ** BREADTH_EXPONENT
    return simpson(integrand, 0.0, math.pi, panels)


def sqrt_profile_by_simpson(panels: int = 2000) -> float:
    """``int_0^1 sqrt(w) dw`` with ``w = u^2``: the integrand becomes ``2u^2``."""
    return simpson(lambda u: 2.0 * u * u, 0.0, 1.0, panels)


def test_taper_integral_matches_simpson() -> None:
    """``TAPER_INTEGRAL`` is a Beta function; Simpson on sin^0.7 must find it."""
    assert pytest.approx(taper_integral_by_simpson(), rel=1e-11) == TAPER_INTEGRAL


def test_closed_form_volume_matches_independent_quadrature() -> None:
    length, beam, draft = 10.0, 2.0, 0.8
    taper = taper_integral_by_simpson()
    profile = sqrt_profile_by_simpson()
    reference = length * beam * draft * taper * profile
    assert analytic_volume(length, beam, draft) == pytest.approx(reference, rel=1e-11)


def test_closed_form_waterplane_and_section_areas_match_quadrature() -> None:
    length, beam, draft = 10.0, 2.0, 0.8
    taper = taper_integral_by_simpson()
    assert analytic_waterplane_area(length, beam) == pytest.approx(
        beam * length * taper, rel=1e-11
    )
    # Section area: full breadth = beam*sqrt(w), integrated over the draft.
    assert analytic_maximum_section_area(beam, draft) == pytest.approx(
        beam * draft * sqrt_profile_by_simpson(), rel=1e-11
    )


def test_closed_form_vertical_centre_matches_quadrature() -> None:
    """VCB is a ratio of two z-integrals; check both, not just the ratio."""
    draft = 0.8
    # w = u^2 removes the 1/sqrt(w) keel term, exactly as in the module: the
    # breadth is beam*u and dz = 2*draft*u du, so area = int 2*u^2 and
    # moment = int 2*u^2*(u^2 - 1)*draft.
    area = simpson(lambda u: 2.0 * u * u, 0.0, 1.0, 2000)
    moment = simpson(lambda u: 2.0 * u * u * (u * u - 1.0) * draft, 0.0, 1.0, 2000)
    reference = moment / area
    assert analytic_vertical_centre_buoyancy(draft) == pytest.approx(reference, rel=1e-13)
    assert analytic_vertical_centre_buoyancy(draft) == pytest.approx(-0.4 * draft)


def test_longitudinal_centre_is_the_midpoint() -> None:
    assert analytic_longitudinal_centre(10.0) == pytest.approx(5.0)
    # Sanity: the taper is symmetric, so its first moment about the midpoint is zero.
    half = simpson(lambda x: (x - 0.5) * math.sin(math.pi * x) ** BREADTH_EXPONENT, 0.0, 1.0, 20000)
    assert abs(half) < 1e-12


# --------------------------------------------------------------------------- #
# The wetted area: no closed form, so two quadrature families must agree.
# --------------------------------------------------------------------------- #

def _wetted_area_by_simpson(length: float, beam: float, draft: float, panels: int) -> float:
    """The same surface element as the module, integrated by Simpson instead of Gauss.

    Surface element from ``(x, (beam/2)*sin(pi*x/length)^p*sqrt(w), draft*(w-1))``
    with ``w = u^2`` and ``x = length*psi/pi``, ``psi = (pi/2)*tau^3``.
    """
    p = BREADTH_EXPONENT

    def integrand(tau: float, u: float) -> float:
        psi = 0.5 * math.pi * tau ** 3
        dpsi_dtau = 1.5 * math.pi * tau * tau
        sine = math.sin(psi)
        if sine <= 0.0:
            return 0.0
        s = sine ** p
        ds_dx = p * sine ** (p - 1.0) * math.cos(psi) * (math.pi / length)
        element = math.sqrt(
            4.0 * u * u * draft * draft
            + (draft * beam * ds_dx * u * u) ** 2
            + (0.5 * beam * s) ** 2
        )
        return dpsi_dtau * element

    return (4.0 * length / math.pi) * simpson_2d(integrand, panels)


def test_wetted_area_quadrature_is_converged() -> None:
    length, beam, draft = 10.0, 2.0, 0.8
    coarse = analytic_wetted_surface_area(length, beam, draft, nodes=60)
    shipped = analytic_wetted_surface_area(length, beam, draft, nodes=80)
    fine = analytic_wetted_surface_area(length, beam, draft, nodes=120)
    assert shipped == pytest.approx(fine, rel=1e-8), "Gauss rule is not converged"
    assert coarse == pytest.approx(fine, rel=1e-7)


def test_wetted_area_matches_an_independent_simpson_product_rule() -> None:
    length, beam, draft = 10.0, 2.0, 0.8
    gauss = analytic_wetted_surface_area(length, beam, draft)
    simpson_value = _wetted_area_by_simpson(length, beam, draft, panels=400)
    # Simpson at 400 panels is an O(h^4) rule on a smooth integrand: agreement to
    # 1e-6 relative means neither family is hiding a systematic term.
    assert gauss == pytest.approx(simpson_value, rel=1e-6)


def test_wetted_area_is_positive_and_scales_quadratically() -> None:
    base = analytic_wetted_surface_area(10.0, 2.0, 0.8)
    scaled = analytic_wetted_surface_area(20.0, 4.0, 1.6)
    assert base > 0.0
    assert scaled == pytest.approx(4.0 * base, rel=1e-12), "area is homogeneous of degree 2"


# --------------------------------------------------------------------------- #
# The mesh path, against the closed forms.
# --------------------------------------------------------------------------- #

RESOLUTIONS = ((11, 5), (21, 9), (41, 13), (161, 41))


def _errors(parameters: HullParameters) -> dict[str, float]:
    return analytic_relative_errors(parameters)["relative_error"]


def test_mesh_hydrostatics_converge_to_the_closed_forms() -> None:
    volume, waterplane, wetted = [], [], []
    for nx, nz in RESOLUTIONS:
        errors = _errors(HullParameters(nx=nx, nz=nz))
        volume.append(abs(errors["displacement_volume"]))
        waterplane.append(abs(errors["waterplane_area"]))
        wetted.append(abs(errors["wetted_surface_area"]))
    for series, label in ((volume, "volume"), (waterplane, "waterplane"), (wetted, "wetted area")):
        assert series == sorted(series, reverse=True), f"{label} error is not monotone in h"
    assert volume[-1] < 2.0e-3
    assert waterplane[-1] < 3.0e-4
    assert wetted[-1] < 4.0e-3


def test_default_resolution_is_inside_three_percent_on_every_quantity() -> None:
    report = analytic_relative_errors(HullParameters(nx=41, nz=13))
    for name, error in report["relative_error"].items():
        assert abs(error) < 0.03, f"{name} is {100 * error:+.3f}% off its closed form"


def test_sampled_surface_is_inscribed_so_it_under_reports() -> None:
    """The chords between stations cut inside the surface: volume and area come out low."""
    errors = _errors(HullParameters(nx=41, nz=13))
    for name in (
        "displacement_volume",
        "waterplane_area",
        "wetted_surface_area",
        "maximum_section_area",
        "block_coefficient",
        "prismatic_coefficient",
    ):
        assert errors[name] < 0.0, f"{name} should be under-reported, got {errors[name]:+.3e}"


def test_coefficients_are_the_taper_integral() -> None:
    reference = analytic_hydrostatics(10.0, 2.0, 0.8)
    assert reference.block_coefficient == pytest.approx((2.0 / 3.0) * TAPER_INTEGRAL, rel=1e-14)
    assert reference.prismatic_coefficient == pytest.approx(TAPER_INTEGRAL, rel=1e-14)
    assert reference.prismatic_coefficient / reference.block_coefficient == pytest.approx(1.5)


def test_reference_agrees_with_the_shipped_volume_formula() -> None:
    reference = analytic_hydrostatics(10.0, 2.0, 0.8)
    assert reference.displacement_volume == pytest.approx(analytic_volume(10.0, 2.0, 0.8))
    assert reference.displacement_mass == pytest.approx(1025.0 * reference.displacement_volume)
    assert "closed form" in reference.basis
    assert reference.waterline_z == 0.0


def test_vertical_centre_error_is_bounded_at_the_default_grid() -> None:
    parameters = HullParameters(nx=41, nz=13)
    mesh, _ = generate_hull(parameters).mesh, None
    measured = compute_hydrostatics(mesh, mesh.bounding_box().maximum[2], closed=True)
    exact = analytic_vertical_centre_buoyancy(parameters.draft)
    assert measured.vertical_centre_buoyancy == pytest.approx(exact, rel=1e-2)
    assert measured.vertical_centre_buoyancy > exact, "the chord rule sits above the keel"


def test_longitudinal_centres_land_on_the_midpoint_at_every_resolution() -> None:
    for nx, nz in RESOLUTIONS:
        hydro = generate_hull(HullParameters(nx=nx, nz=nz)).hydrostatics
        assert hydro.longitudinal_centre_buoyancy == pytest.approx(5.0, abs=1e-9)
        assert hydro.longitudinal_centre_flotation == pytest.approx(5.0, abs=1e-9)


def test_refinement_estimate_is_optimistic_about_the_true_error() -> None:
    """The coarse/refined relative change under-states the error against the closed form.

    ``estimate_discretisation_error`` differences two approximations; the closed
    form says how far the coarse one really is. Recording the gap keeps anyone
    from quoting the estimator as an error bound.
    """
    report = analytic_relative_errors(HullParameters(nx=41, nz=13))
    outside = report["outside_discretisation_estimate"]
    assert "displacement_volume" in outside
    measured = abs(outside["displacement_volume"]["measured_relative_error"])
    reported = outside["displacement_volume"]["reported_relative_change"]
    assert measured > reported
    assert measured / reported < 5.0, "the estimator would then be useless, not merely optimistic"


@pytest.mark.parametrize("bad", [0.0, -1.0, float("nan"), float("inf")])
def test_invalid_dimensions_are_rejected(bad: float) -> None:
    with pytest.raises(ValueError):
        analytic_hydrostatics(bad, 2.0, 0.8)
    with pytest.raises(ValueError):
        analytic_waterplane_area(10.0, bad)
    with pytest.raises(ValueError):
        analytic_maximum_section_area(2.0, bad)
    with pytest.raises(ValueError):
        analytic_vertical_centre_buoyancy(bad)
    with pytest.raises(ValueError):
        analytic_wetted_surface_area(10.0, 2.0, bad)
    with pytest.raises(ValueError):
        analytic_longitudinal_centre(bad)


def test_density_is_an_explicit_input() -> None:
    fresh = analytic_hydrostatics(10.0, 2.0, 0.8, water_density=1000.0)
    sea = analytic_hydrostatics(10.0, 2.0, 0.8, water_density=1025.0)
    assert fresh.displacement_volume == pytest.approx(sea.displacement_volume)
    assert fresh.displacement_mass == pytest.approx(sea.displacement_volume * 1000.0)
    with pytest.raises(ValueError):
        analytic_hydrostatics(10.0, 2.0, 0.8, water_density=0.0)
