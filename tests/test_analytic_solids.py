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
    """A 200:1 disc needs a far higher order than the cap; the reference refuses, loudly.

    This is the documented applicability limit of the tensor rule, asserted so it stays a
    known refusal rather than becoming a silent loss of digits on a slender body.
    """
    with pytest.raises(ArithmeticError, match="not converged"):
        solids.ellipsoid_surface_area([10.0, 10.0, 0.05], max_order=128)


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
