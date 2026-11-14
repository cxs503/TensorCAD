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
not stopped moving, the same discipline the golden suite applies to its tolerances.

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
#: Quadrature order is doubled up to this cap before convergence is declared impossible.
#: The tensor rule needs an order proportional to the equatorial-to-polar axis ratio (the
#: area element develops a boundary layer of width ~c/a), and cost grows as order**2, so
#: the cap bounds a call at about a second. A caller with a more slender body passes a
#: larger ``max_order`` rather than getting a silently truncated area.
MAX_ORDER = 512


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


def _area_quadrature(semi_axes: Sequence[float], order: int) -> float:
    a, b, c = _semi_axes(semi_axes)
    if order < 2:
        raise ValueError(f"order must be at least 2, got {order}")
    nodes, weights = np.polynomial.legendre.leggauss(order)
    # u in [0, pi], v in [0, 2 pi): shift the [-1, 1] rule, carrying the Jacobian
    u = 0.5 * math.pi * (nodes + 1.0)
    w_u = 0.5 * math.pi * weights
    v = math.pi * (nodes + 1.0)
    w_v = math.pi * weights
    sin_u, cos_u = np.sin(u)[:, None], np.cos(u)[:, None]
    sin_v, cos_v = np.sin(v)[None, :], np.cos(v)[None, :]
    element = sin_u * np.sqrt(
        a * a * b * b * cos_u**2 + c * c * sin_u**2 * (b * b * cos_v**2 + a * a * sin_v**2)
    )
    return float(w_u @ element @ w_v)


def ellipsoid_surface_area(
    semi_axes: Sequence[float],
    *,
    order: int = DEFAULT_ORDER,
    tolerance: float = CONVERGENCE_TOLERANCE,
    max_order: int = MAX_ORDER,
) -> float:
    """Surface area of a (triaxial) ellipsoid, by converged Gauss-Legendre quadrature.

    The area element above is integrated at ``order`` and ``2 * order``; a value is returned
    only when the two passes agree to ``tolerance`` relative, so the number is a converged
    integral rather than one quadrature's opinion. Slender ellipsoids need a higher order for
    the same accuracy, so the order is doubled (up to ``max_order``) before convergence is
    declared impossible. Raises :class:`ArithmeticError` in that case -- it never returns a
    silently truncated result.

    Applicability: with the default cap this is verified for axis ratios up to about 5:1
    (the golden suite's 4:3:1 body converges at order 128). More slender bodies may exhaust
    the cap and raise; raising ``max_order`` buys accuracy at quadratic cost in runtime.
    """
    coarse = _area_quadrature(semi_axes, order)
    while True:
        fine = _area_quadrature(semi_axes, 2 * order)
        if fine <= 0.0:
            raise ArithmeticError(f"non-positive area {fine!r} for semi-axes {tuple(semi_axes)}")
        change = abs(fine - coarse) / fine
        if change <= tolerance:
            return fine
        if 2 * order >= max_order:
            raise ArithmeticError(
                f"ellipsoid area quadrature has not converged: order {order} -> {coarse!r}, "
                f"order {2 * order} -> {fine!r} (relative change {change:.3e} > "
                f"{tolerance:.1e}) for semi-axes {tuple(semi_axes)}"
            )
        coarse, order = fine, 2 * order


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
        e = math.sqrt(1.0 - (c * c) / (a * a))
        return 2.0 * math.pi * a * a * (1.0 + (1.0 - e * e) / e * math.atanh(e))
    # prolate: elongated along the axis of revolution
    e = math.sqrt(1.0 - (a * a) / (c * c))
    return 2.0 * math.pi * a * a * (1.0 + c / (a * e) * math.asin(e))
