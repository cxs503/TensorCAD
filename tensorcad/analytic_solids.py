"""Exact geometry of quadric primitives, computed without any mesher.

Why this module exists
----------------------
A regression suite is only as good as the references it compares against. A reference
that is produced by the same mesher the suite is checking proves nothing: a defect in
the mesher moves the measurement and the reference together, and the check stays green.
So references live here, derived from closed forms or from quadrature of an exact
parametric area element, and they never touch :mod:`tensorcad.geometry_field` or any
marching-tetrahedra code.

The triaxial ellipsoid surface area is the interesting case. It has no elementary closed
form -- it is an elliptic integral -- which is why it was for a while only drift-locked in
the golden suite rather than asserted. It does, however, have an exact *integral*
representation. With ``r(u, v) = (a sin u cos v, b sin u sin v, c cos u)``, the area
element is

.. math::

    |r_u \\times r_v| = \\sin u \\sqrt{a^2 b^2 \\cos^2 u
                      + c^2 \\sin^2 u \\,(b^2 \\cos^2 v + a^2 \\sin^2 v)}

which is smooth for positive semi-axes, so a tensor-product Gauss-Legendre rule on
``u in [0, pi]``, ``v in [0, 2 pi)`` converges geometrically. The order is not a magic
number: :func:`ellipsoid_surface_area` doubles it and refuses to return a value that has
not stopped moving, the same discipline the golden suite applies to its tolerances. Nor is
the parametrisation above innocuous on a flat or slender body, where the element picks up
a corner from ``|cos u|`` or ``|sin v|`` and stops converging at any order; the reference
then splits the interval at that corner, which is a change of integrand rather than a
request for more nodes.

For spheroids the area is elementary, and :func:`spheroid_surface_area` exposes it so the
quadrature can be checked against something that shares no code with it.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np

__all__ = [
    "ellipsoid_surface_area",
    "ellipsoid_volume",
    "spheroid_surface_area",
]

#: Quadrature order of the coarse pass; the returned value comes from twice this.
DEFAULT_ORDER = 64
#: The two passes must agree to this relative level before a value is returned.
CONVERGENCE_TOLERANCE = 1e-11
#: Quadrature order is doubled up to this cap before the unsplit rule is declared
#: insufficient. The tensor rule needs an order proportional to the equatorial-to-polar
#: axis ratio, and cost grows as order**2, so the cap bounds a call at about a second.
MAX_ORDER = 512
#: Sub-interval cap for the panel-split rule the reference falls back to. On a flat or
#: slender body the area element has a corner where ``cos u`` or ``cos v`` changes sign,
#: and no order of a single tensor block removes it; splitting the intervals at those
#: corners restores geometric convergence. The ladder is 1, 2, 4, ... panels per direction,
#: and a split value is returned only when the next panel count corroborates it.
MAX_PANELS = 8


def _semi_axes(semi_axes: Sequence[float]) -> tuple[float, float, float]:
    axes = tuple(float(value) for value in semi_axes)
    if len(axes) != 3:
        raise ValueError(f"semi_axes must have three entries, got {len(axes)}")
    if any(not math.isfinite(value) or value <= 0.0 for value in axes):
        raise ValueError(f"semi_axes must be finite and positive, got {axes}")
    return axes  # type: ignore[return-value]


def ellipsoid_volume(semi_axes: Sequence[float]) -> float:
    """Volume of an ellipsoid with the given semi-axes: ``4/3 pi a b c``."""
    a, b, c = _semi_axes(semi_axes)
    return 4.0 / 3.0 * math.pi * a * b * c


def _area_quadrature(semi_axes: Sequence[float], order: int, panels: int = 1) -> float:
    """Tensor Gauss-Legendre integral of the exact area element.

    ``panels`` splits ``u in [0, pi]`` and ``v in [0, 2 pi)`` into that many equal
    sub-intervals each. That is what removes the corner the single block stalls on: for a
    flat or slender body the element behaves like ``|cos u|`` (or ``|sin v|`` in ``v``),
    which is smooth only one sub-interval at a time. ``panels=1`` is the unsplit rule and
    returns bit-identical numbers to the single-block implementation it replaces.
    """
    a, b, c = _semi_axes(semi_axes)
    if order < 2:
        raise ValueError(f"order must be at least 2, got {order}")
    if panels < 1:
        raise ValueError(f"panels must be at least 1, got {panels}")
    nodes, weights = np.polynomial.legendre.leggauss(order)
    total = 0.0
    for i in range(panels):
        # shift the [-1, 1] rule onto u in [i, i+1] * pi / panels, carrying the Jacobian
        u_lo, u_hi = math.pi * i / panels, math.pi * (i + 1) / panels
        u = 0.5 * (u_hi - u_lo) * (nodes + 1.0) + u_lo
        w_u = 0.5 * (u_hi - u_lo) * weights
        sin_u, cos_u = np.sin(u)[:, None], np.cos(u)[:, None]
        for j in range(panels):
            v_lo, v_hi = 2.0 * math.pi * j / panels, 2.0 * math.pi * (j + 1) / panels
            v = 0.5 * (v_hi - v_lo) * (nodes + 1.0) + v_lo
            w_v = 0.5 * (v_hi - v_lo) * weights
            sin_v, cos_v = np.sin(v)[None, :], np.cos(v)[None, :]
            element = sin_u * np.sqrt(
                a * a * b * b * cos_u**2
                + c * c * sin_u**2 * (b * b * cos_v**2 + a * a * sin_v**2)
            )
            total += float(w_u @ element @ w_v)
    return total


def _converged_area(
    semi_axes: Sequence[float],
    order: int,
    tolerance: float,
    max_order: int,
    panels: int,
) -> float | None:
    """Area at one panel count, or ``None`` if the order cap is exhausted first."""
    coarse = _area_quadrature(semi_axes, order, panels)
    while 2 * order <= max_order:
        fine = _area_quadrature(semi_axes, 2 * order, panels)
        if fine <= 0.0:
            raise ArithmeticError(f"non-positive area {fine!r} for semi-axes {tuple(semi_axes)}")
        if abs(fine - coarse) / fine <= tolerance:
            return fine
        coarse, order = fine, 2 * order
    return None


def ellipsoid_surface_area(
    semi_axes: Sequence[float],
    *,
    order: int = DEFAULT_ORDER,
    tolerance: float = CONVERGENCE_TOLERANCE,
    max_order: int = MAX_ORDER,
    max_panels: int = MAX_PANELS,
) -> float:
    """Surface area of a (triaxial) ellipsoid, by converged Gauss-Legendre quadrature.

    The area element above is integrated at ``order`` and ``2 * order``; a value is returned
    only when the two passes agree to ``tolerance`` relative, so the number is a converged
    integral rather than one quadrature's opinion. Slender ellipsoids need a higher order for
    the same accuracy, so the order is doubled (up to ``max_order``) before anything else is
    tried; if the order cap is exhausted, the interval is split into panels (up to
    ``max_panels``) and the search restarts, because the obstacle on a flat or slender body is
    a corner in the integrand rather than a resolution shortfall. A split value is returned
    only when the next panel count corroborates it. Raises :class:`ArithmeticError` if
    neither route attests a value -- it never returns a silently truncated result.

    Applicability (measured, not assumed: see ``artifacts/golden_ellipsoid_area``): the
    unsplit rule converges for axis ratios up to about 15:1, and the panel-split rule for
    ratios up to about 1000:1, at 0.1-2 s. Beyond that the call raises; raising ``max_order``
    or ``max_panels`` buys accuracy at quadratic and linear cost respectively.
    """
    if max_panels < 1:
        raise ValueError(f"max_panels must be at least 1, got {max_panels}")
    previous: float | None = None
    panels = 1
    while panels <= max_panels:
        value = _converged_area(semi_axes, order, tolerance, max_order, panels)
        if value is not None:
            if panels == 1:
                # the unsplit rule attests on order refinement alone, as it always has
                return value
            if previous is not None and abs(value - previous) <= tolerance * value:
                return value
            previous = value
        panels *= 2
    raise ArithmeticError(
        f"ellipsoid area quadrature has not converged: order up to {max_order} per "
        f"direction and panels up to {max_panels} did not attest a value to "
        f"{tolerance:.1e} relative for semi-axes {tuple(semi_axes)}"
    )


def spheroid_surface_area(equatorial: float, polar: float) -> float:
    """Surface area of a spheroid (two equal semi-axes) in elementary closed form.

    ``equatorial`` is the radius in the plane of symmetry, ``polar`` the radius along the
    axis of revolution. Oblate (``polar < equatorial``) uses ``atanh``, prolate
    (``polar > equatorial``) uses ``asin`` on the eccentricity; the sphere is the limit of
    either. Kept separate from the quadrature so the two can be checked against each other.
    """
    a = float(equatorial)
    c = float(polar)
    if not (math.isfinite(a) and math.isfinite(c)) or a <= 0.0 or c <= 0.0:
        raise ValueError(f"spheroid radii must be finite and positive, got {(a, c)}")
    if math.isclose(a, c, rel_tol=1e-15):
        return 4.0 * math.pi * a * a
    if c < a:  # oblate: flattened along the axis of revolution
        # k = c / a in (0, 1). The eccentricity comes from e^2 = 1 - k^2 = (1 - k)(1 + k)
        # rather than a literal 1 - k^2, so it stays accurate when the spheroid is nearly
        # spherical: 1 - k is exact for k in [0.5, 1], while 1 - k^2 cancels.
        k = c / a
        kk = k * k
        if kk == 0.0:  # flatter than the square of the smallest double: the disc limit
            return 2.0 * math.pi * a * a
        flat = 1.0 - k
        e = math.sqrt(flat * (1.0 + k))
        # atanh(e) = log((1 + e) / k) = log1p(e) - log1p(-(1 - k)). The logarithmic form is
        # finite as k -> 0 -- a flat spheroid's area tends to 2 pi a^2 -- whereas the
        # literal atanh(e) would be asked for atanh(1.0) there and raise a domain error.
        atanh_e = math.log1p(e) - math.log1p(-flat)
        return 2.0 * math.pi * a * a * (1.0 + kk / e * atanh_e)
    # prolate: elongated along the axis of revolution
    k = a / c
    e = math.sqrt((1.0 - k) * (1.0 + k))
    # asin(e) = pi/2 - asin(k): take whichever of the two is not close to pi/2, so the
    # needle limit (k -> 0) keeps full precision instead of subtracting from pi/2.
    asin_e = math.asin(e) if k >= 0.5 else 0.5 * math.pi - math.asin(k)
    return 2.0 * math.pi * a * a * (1.0 + asin_e / (k * e))
