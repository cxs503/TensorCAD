"""Verification of :mod:`tensorcad.analytic_solids` — the golden suite's area reference.

The point of these tests is *independence*: the ellipsoid area reference is what the
golden suite compares a marching-tetrahedra mesh against, so it must not be able to
inherit a mesher defect. Three routes to the same number are checked here:

* the quadrature of the exact area element (what the suite uses),
* the elementary spheroid closed forms (oblate ``atanh``, prolate ``asin``) — no shared
  code with the quadrature,
* the sphere, where both must hit ``4 pi r^2``.

A reference that agrees with itself but not with these is worthless, and a reference that
returns a number when it has not converged is worse than none.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from tensorcad import analytic_solids as solids


# --------------------------------------------------------------------------------------
# the sphere: the one case with an elementary answer
# --------------------------------------------------------------------------------------
def test_quadrature_reproduces_the_sphere_area() -> None:
    for radius in (0.25, 1.0, 7.5):
        area = solids.ellipsoid_surface_area([radius] * 3)
        expected = 4.0 * math.pi * radius**2
        assert abs(area - expected) / expected < 1e-12


def test_spheroid_closed_form_reproduces_the_sphere() -> None:
    assert solids.spheroid_surface_area(2.0, 2.0) == pytest.approx(4.0 * math.pi * 4.0, rel=1e-15)


# --------------------------------------------------------------------------------------
# quadrature vs the closed form, over the interesting range of aspect ratios
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize("ratio", [1.05, 1.5, 2.0, 4.0, 10.0])
def test_quadrature_agrees_with_the_oblate_closed_form(ratio: float) -> None:
    """Oblate: two equal equatorial semi-axes larger than the polar one."""
    area = solids.ellipsoid_surface_area([ratio, ratio, 1.0])
    expected = solids.spheroid_surface_area(ratio, 1.0)
    assert abs(area - expected) / expected < 1e-10


@pytest.mark.parametrize("ratio", [1.05, 1.5, 2.0, 4.0, 10.0])
def test_quadrature_agrees_with_the_prolate_closed_form(ratio: float) -> None:
    """Prolate: elongated along the axis of revolution (the ``asin`` branch)."""
    area = solids.ellipsoid_surface_area([1.0, 1.0, ratio])
    expected = solids.spheroid_surface_area(1.0, ratio)
    assert abs(area - expected) / expected < 1e-10


def test_closed_form_branches_meet_at_the_sphere() -> None:
    """Oblate and prolate must converge to the same limit as the eccentricity vanishes."""
    sphere = 4.0 * math.pi
    for eps in (1e-2, 1e-3, 1e-4):
        oblate = solids.spheroid_surface_area(1.0, 1.0 - eps)
        prolate = solids.spheroid_surface_area(1.0, 1.0 + eps)
        assert abs(oblate - sphere) / sphere < 10 * eps
        assert abs(prolate - sphere) / sphere < 10 * eps


# --------------------------------------------------------------------------------------
# properties the reference must have
# --------------------------------------------------------------------------------------
def test_area_is_invariant_under_a_permutation_of_the_semi_axes() -> None:
    """The area is a property of the body, not of which axis the parametrisation calls z."""
    reference = solids.ellipsoid_surface_area([2.0, 1.5, 0.5])
    for permutation in ([1.5, 0.5, 2.0], [0.5, 2.0, 1.5], [0.5, 1.5, 2.0]):
        assert solids.ellipsoid_surface_area(permutation) == pytest.approx(reference, rel=1e-12)


def test_a_coarse_start_is_escalated_until_converged() -> None:
    """Starting far too coarse must not decide the answer: the order is doubled instead."""
    reference = solids.ellipsoid_surface_area([2.0, 1.5, 0.5])
    for start in (2, 4, 8, 16):
        assert solids.ellipsoid_surface_area([2.0, 1.5, 0.5], order=start) == pytest.approx(
            reference, rel=1e-11
        )


def test_a_non_converged_integral_raises_instead_of_returning() -> None:
    """Fail closed: an impossible convergence demand must not yield a number."""
    with pytest.raises(ArithmeticError, match="not converged"):
        solids.ellipsoid_surface_area([2.0, 1.5, 0.5], tolerance=1e-20, max_order=128)


def test_a_body_beyond_the_cap_raises_rather_than_returning_a_truncated_area() -> None:
    """A 200:1 disc needs more than the panel ladder is allowed to spend; it refuses, loudly.

    This is the documented applicability limit, asserted so it stays a known refusal rather
    than becoming a silent loss of digits on a slender body.
    """
    with pytest.raises(ArithmeticError, match="not converged"):
        solids.ellipsoid_surface_area([10.0, 10.0, 0.05], max_order=128, max_panels=2)


def test_the_unsplit_rule_alone_still_refuses_a_slender_body() -> None:
    """The historical limit is real and stays visible: one tensor block cannot do 20:1."""
    with pytest.raises(ArithmeticError, match="not converged"):
        solids.ellipsoid_surface_area([20.0, 20.0, 1.0], max_panels=1)


# --------------------------------------------------------------------------------------
# panel splitting: the corner in the area element, and what removing it buys
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize("ratio", [20.0, 100.0])
def test_panel_splitting_reaches_the_slender_discs_the_unsplit_rule_refused(ratio: float) -> None:
    """Splitting ``u`` at its corner restores geometric convergence where order alone stalls.

    The check is still the elementary oblate closed form, i.e. the same independent route
    the rest of this file uses; the panel count only changes how the integral is evaluated.
    """
    area = solids.ellipsoid_surface_area([ratio, ratio, 1.0])
    expected = solids.spheroid_surface_area(ratio, 1.0)
    assert abs(area - expected) / expected < 1e-12


def test_panel_splitting_handles_a_body_that_corners_in_both_directions() -> None:
    """(8, 1, 0.1) has a corner in ``cos u`` and another in ``sin v``: no closed form exists,
    so it is checked against a different quadrature family (cell midpoint rule) at 1e-5."""
    area = solids.ellipsoid_surface_area([8.0, 1.0, 0.1])
    assert abs(_midpoint_area([8.0, 1.0, 0.1]) - area) / area < 1e-5


def test_the_committed_reference_value_is_unmoved_by_the_panel_refactor() -> None:
    """The suite's blessed number must not drift: ``panels=1`` is the historical rule and it
    still returns the bit-identical area that ``docs/golden`` was generated from."""
    assert solids.ellipsoid_surface_area([2.0, 1.5, 0.5]) == 22.076658276751232


def _midpoint_area(semi_axes: tuple[float, float, float], cells: int = 800) -> float:
    """Area element summed on a cell-centre grid: shares no nodes with the Gauss rule."""
    a, b, c = semi_axes
    u = (np.arange(cells) + 0.5) * math.pi / cells
    v = (np.arange(cells) + 0.5) * 2.0 * math.pi / cells
    sin_u, cos_u = np.sin(u)[:, None], np.cos(u)[:, None]
    sin_v, cos_v = np.sin(v)[None, :], np.cos(v)[None, :]
    element = sin_u * np.sqrt(
        a * a * b * b * cos_u**2 + c * c * sin_u**2 * (b * b * cos_v**2 + a * a * sin_v**2)
    )
    return float(element.sum()) * (math.pi / cells) * (2.0 * math.pi / cells)


# --------------------------------------------------------------------------------------
# the elementary spheroid forms: the flat limit and the near-spherical digits
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("polar", "tolerance"),
    # the residual is k^2 log(2/k) with k = polar / equatorial: 1.4e-11 at k = 1e-6 and
    # 2e-17 at k = 1e-9, so the tolerance is the size of the next term, not a fudge
    [(1e-6, 1e-10), (1e-9, 1e-15), (1e-12, 1e-15), (1e-300, 1e-15)],
)
def test_the_flat_spheroid_limit_is_the_disc_area(polar: float, tolerance: float) -> None:
    """A spheroid flattened without limit tends to a disc of area ``2 pi a^2``.

    Regression: the old ``atanh(e)`` branch was handed ``atanh(1.0)`` once ``e`` rounded to
    one, and raised a domain error on a perfectly well-posed limit.
    """
    assert solids.spheroid_surface_area(1.0, polar) == pytest.approx(
        2.0 * math.pi, rel=tolerance
    )


@pytest.mark.parametrize("delta", [1e-3, 1e-6, 1e-9, 1e-12])
def test_the_closed_forms_agree_with_the_quadrature_near_the_sphere(delta: float) -> None:
    """Digits near ``a = c`` are where a cancelling eccentricity shows up: both branches must
    still meet the quadrature to 1e-13 relative."""
    for polar in (1.0 - delta, 1.0 + delta):
        quadrature = solids.ellipsoid_surface_area([1.0, 1.0, polar])
        closed_form = solids.spheroid_surface_area(1.0, polar)
        assert abs(closed_form - quadrature) / quadrature < 1e-13


@pytest.mark.parametrize("ratio", [1e3, 1e6, 1e9])
def test_the_needle_limit_keeps_its_precision(ratio: float) -> None:
    """``asin(e)`` is a subtraction from ``pi/2`` for a needle, and the old branch paid for
    it: 2.8e-11 at 1e6:1 and 6.4e-10 at 1e9:1 against the quadrature. Reflecting the arcsine
    keeps it at 1e-14, checked against the quadrature, which converges here."""
    quadrature = solids.ellipsoid_surface_area([1.0, 1.0, ratio])
    closed_form = solids.spheroid_surface_area(1.0, ratio)
    assert abs(closed_form - quadrature) / quadrature < 1e-13


@pytest.mark.parametrize(
    "semi_axes",
    [(1.0, 1.0), (1.0, 1.0, 1.0, 1.0), (1.0, 0.0, 1.0), (-1.0, 1.0, 1.0), (math.nan, 1.0, 1.0)],
)
def test_invalid_semi_axes_are_rejected(semi_axes: tuple[float, ...]) -> None:
    with pytest.raises(ValueError):
        solids.ellipsoid_surface_area(semi_axes)


def test_volume_is_the_elementary_product() -> None:
    assert solids.ellipsoid_volume([2.0, 1.5, 0.5]) == pytest.approx(2.0 * math.pi, rel=1e-15)
    assert solids.ellipsoid_volume([1.0, 1.0, 1.0]) == pytest.approx(4.0 / 3.0 * math.pi, rel=1e-15)


def test_area_bounds_the_volume_equivalent_sphere() -> None:
    """Isoperimetric sanity: an ellipsoid is not smaller than the sphere of equal volume."""
    axes = [2.0, 1.5, 0.5]
    volume = solids.ellipsoid_volume(axes)
    equal_volume = 4.0 * math.pi * (3.0 * volume / (4.0 * math.pi)) ** (2.0 / 3.0)
    assert solids.ellipsoid_surface_area(axes) > equal_volume
