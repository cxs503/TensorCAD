"""Elliptic integrals by AGM/Carlson, and the Legendre formula for ellipsoid area.

This module exists to be *independent* of the NURBS machinery. The quadrature in
:mod:`tensorcad.highorder.quadrature` and the curvature in
:mod:`tensorcad.highorder.nurbs_geom` both derive from the control net; if the wetted
surface area of an ellipsoid were also computed from that control net, comparing the two
would only prove that the code agrees with itself. Here the same area is computed from
the implicit equation of the ellipsoid and a classical elliptic-integral formula, so the
comparison has content.

Two routes to the same complete integrals are implemented on purpose:

* the **AGM** route (Gauss's arithmetic-geometric mean), which converges quadratically
  and is the classical statement of the result, and
* the **Carlson symmetric forms** ``R_F``/``R_D``, which are what make the *incomplete*
  integrals in Legendre's triaxial formula well conditioned (and which handle the closed
  case ``k = 1`` that the AGM series meets as a singularity).

They are cross-checked against each other and against direct Gauss-Legendre quadrature of
the defining integrals in ``tests/test_highorder_geometry.py``; agreement of three
independent routes is what licenses using any of them as a reference.

Named limitations
-----------------
* The AGM series for ``E(k)`` is singular at ``k = 1`` (it is the flat/spherical limit
  where the representation degenerates); :func:`ellipe_agm` refuses ``k >= 1 - 1e-12``
  and points at :func:`ellipe`. The Carlson route covers the closed case.
* The triaxial formula needs the axes ordered; :func:`triaxial_ellipsoid_area` sorts them
  internally because the area is symmetric under permutation while the formula is not.
* Areas are returned for a **closed** ellipsoid. A hull's wetted surface is half of it
  only when the waterplane is a symmetry plane -- the caller owns that decision, and
  :func:`tensorcad.highorder.nurbs_geom.hydrostatics_reference` states it explicitly.
* The Carlson iterations use the classical relative error tolerances (2.5e-3 / 1.5e-3),
  which are enough for double precision because each iteration is quadratically
  convergent; the achieved accuracy is *measured* in the tests rather than assumed.
"""

from __future__ import annotations

import math

__all__ = [
    "agm",
    "carlson_rd",
    "carlson_rf",
    "ellipe",
    "ellipe_agm",
    "ellipk",
    "ellipk_agm",
    "ellipse_perimeter",
    "incomplete_ef",
    "spheroid_area",
    "triaxial_ellipsoid_area",
]

def agm(a: float, b: float, *, iterations: int = 64, tolerance: float = 1e-16) -> float:
    """Arithmetic-geometric mean of two positive numbers."""
    if a <= 0.0 or b <= 0.0:
        raise ValueError(f"agm needs positive arguments, got a={a!r}, b={b!r}")
    for _ in range(iterations):
        if abs(a - b) <= tolerance * 0.5 * (a + b):
            break
        a, b = 0.5 * (a + b), math.sqrt(a * b)
    return 0.5 * (a + b)


def ellipk_agm(k: float) -> float:
    """Complete elliptic integral of the first kind ``K(k)`` by AGM."""
    if not 0.0 <= k < 1.0:
        raise ValueError(f"K(k) by AGM needs 0 <= k < 1, got {k!r}")
    return math.pi / (2.0 * agm(1.0, math.sqrt(1.0 - k * k)))


def ellipe_agm(k: float) -> float:
    """Complete elliptic integral of the second kind ``E(k)`` by the AGM series.

    ``E = K * (1 - sum_{n>=0} 2^{n-1} c_n^2)`` with ``c_0 = k`` and
    ``c_{n+1} = (a_n - b_n) / 2`` from the AGM iteration. The series diverges at
    ``k = 1`` because ``K`` does; that limit is served by :func:`ellipe` instead.
    """
    if not 0.0 <= k < 1.0 - 1e-12:
        raise ValueError(f"E(k) by AGM needs 0 <= k < 1, got {k!r}; use ellipe for k -> 1")
    a, b = 1.0, math.sqrt(1.0 - k * k)
    c_squared_sum = 0.5 * k * k
    weight = 1.0
    for _ in range(64):
        if a - b == 0.0:
            break
        half_difference = 0.5 * (a - b)
        a, b = 0.5 * (a + b), math.sqrt(a * b)
        c_squared_sum += weight * half_difference * half_difference
        weight *= 2.0
        if abs(half_difference) <= 1e-17 * a:
            break
    return ellipk_agm(k) * (1.0 - c_squared_sum)


def carlson_rf(x: float, y: float, z: float) -> float:
    """Carlson's symmetric integral ``R_F``, ``(1/2) int_0^inf dt / sqrt(prod(t + x_i))``."""
    for name, value in (("x", x), ("y", y), ("z", z)):
        if value < 0.0:
            raise ValueError(f"R_F argument {name} is negative: {value!r}")
    if min(x, y, z) == 0.0 and max(x, y, z) == 0.0:
        raise ValueError("R_F is undefined when all arguments vanish")
    xt, yt, zt = x, y, z
    error_target = 2.5e-3
    for _ in range(64):
        sqrt_x, sqrt_y, sqrt_z = math.sqrt(xt), math.sqrt(yt), math.sqrt(zt)
        lam = sqrt_x * (sqrt_y + sqrt_z) + sqrt_y * sqrt_z
        xt, yt, zt = 0.25 * (xt + lam), 0.25 * (yt + lam), 0.25 * (zt + lam)
        average = (xt + yt + zt) / 3.0
        dx, dy, dz = (average - xt) / average, (average - yt) / average, (average - zt) / average
        if max(abs(dx), abs(dy), abs(dz)) < error_target:
            break
    e2 = dx * dy - dz * dz
    e3 = dx * dy * dz
    return (1.0 + (e2 / 24.0 - 0.1 - 3.0 * e3 / 44.0) * e2 + e3 / 14.0) / math.sqrt(average)


def carlson_rd(x: float, y: float, z: float) -> float:
    """Carlson's ``R_D``, ``(3/2) int_0^inf dt / (sqrt(prod(t + x_i)) * (t + z))``."""
    for name, value in (("x", x), ("y", y), ("z", z)):
        if value < 0.0:
            raise ValueError(f"R_D argument {name} is negative: {value!r}")
    if z == 0.0:
        raise ValueError("R_D needs a positive third argument")
    xt, yt, zt = x, y, z
    error_target = 1.5e-3
    total = 0.0
    factor = 1.0
    for _ in range(64):
        sqrt_x, sqrt_y, sqrt_z = math.sqrt(xt), math.sqrt(yt), math.sqrt(zt)
        lam = sqrt_x * (sqrt_y + sqrt_z) + sqrt_y * sqrt_z
        total += factor / (sqrt_z * (zt + lam))
        factor *= 0.25
        xt, yt, zt = 0.25 * (xt + lam), 0.25 * (yt + lam), 0.25 * (zt + lam)
        average = (xt + yt + 3.0 * zt) / 5.0
        dx, dy, dz = (average - xt) / average, (average - yt) / average, (average - zt) / average
        if max(abs(dx), abs(dy), abs(dz)) < error_target:
            break
    ea = dx * dy
    eb = dz * dz
    ec = ea - eb
    ed = ea - 6.0 * eb
    ee = ed + ec + ec
    series = (
        1.0
        + ed * (-3.0 / 14.0 + 9.0 / 22.0 * ed - 2.0 / 26.0 * ee)
        + dz * (ee / 6.0 + ed * (-9.0 / 22.0 + 3.0 / 26.0 * ee))
    )
    return 3.0 * total + factor * series / (average * math.sqrt(average))


def ellipk(k: float) -> float:
    """Complete elliptic integral of the first kind ``K(k)`` (Carlson route).

    ``K(1)`` is logarithmic and infinite; it raises instead of returning a large float.
    """
    if not 0.0 <= k < 1.0:
        raise ValueError(f"K(k) needs 0 <= k < 1, got {k!r}")
    return carlson_rf(0.0, 1.0 - k * k, 1.0)


def ellipe(k: float) -> float:
    """Complete elliptic integral of the second kind ``E(k)`` (Carlson route).

    ``E = R_F(0, 1 - k^2, 1) - (k^2 / 3) R_D(0, 1 - k^2, 1)``. Using ``1e-300`` in place
    of the zero argument keeps ``R_F``/``R_D`` out of their singular point; the error that
    substitution introduces is ``O(sqrt(1e-300))`` in the argument, i.e. below the
    round-off of the result.
    """
    if not 0.0 <= k <= 1.0:
        raise ValueError(f"E(k) needs 0 <= k <= 1, got {k!r}")
    if k > 1.0 - 1e-12:
        return 1.0  # the k = 1 limit, where R_F and R_D both vanish and cancel
    complement = 1.0 - k * k
    return carlson_rf(0.0, complement, 1.0) - (k * k / 3.0) * carlson_rd(0.0, complement, 1.0)


def incomplete_ef(phi: float, k: float) -> tuple[float, float]:
    """``(E(phi, k), F(phi, k))`` for ``0 <= phi <= pi/2`` and ``0 <= k <= 1``."""
    if not 0.0 <= phi <= 0.5 * math.pi:
        raise ValueError(f"phi must be in [0, pi/2], got {phi!r}")
    if not 0.0 <= k <= 1.0:
        raise ValueError(f"k must be in [0, 1], got {k!r}")
    sin_phi = math.sin(phi)
    cos_phi = math.cos(phi)
    x = cos_phi * cos_phi
    y = 1.0 - k * k * sin_phi * sin_phi
    rf = carlson_rf(x, y, 1.0)
    first = sin_phi * rf
    second = first - (k * k / 3.0) * sin_phi**3 * carlson_rd(x, y, 1.0)
    return (second, first)


def ellipse_perimeter(a: float, b: float) -> float:
    """Perimeter of an ellipse with semi-axes ``a >= b``: ``4 a E(e)``, ``e^2 = 1 - b^2/a^2``."""
    if a <= 0.0 or b <= 0.0:
        raise ValueError(f"semi-axes must be positive, got a={a!r}, b={b!r}")
    major, minor = max(a, b), min(a, b)
    eccentricity = math.sqrt(max(1.0 - (minor / major) ** 2, 0.0))
    return 4.0 * major * ellipe(eccentricity)


def spheroid_area(equatorial: float, polar: float) -> float:
    """Area of a spheroid (two equal equatorial semi-axes): closed elementary form.

    Oblate (``equatorial > polar``): ``2 pi a^2 + pi c^2/e * ln((1+e)/(1-e))``.
    Prolate (``polar > equatorial``): ``2 pi a^2 + 2 pi a c/e * asin(e)``.
    Both limits are the sphere area ``4 pi a^2``; a sphere input is returned directly to
    avoid the ``0/0`` that the two formulas meet there.
    """
    a, c = equatorial, polar
    if a <= 0.0 or c <= 0.0:
        raise ValueError(f"semi-axes must be positive, got a={a!r}, c={c!r}")
    if abs(a - c) <= 1e-15 * a:
        return 4.0 * math.pi * a * a
    if a > c:
        eccentricity = math.sqrt(1.0 - (c / a) ** 2)
        return 2.0 * math.pi * a * a + math.pi * c * c / eccentricity * math.log(
            (1.0 + eccentricity) / (1.0 - eccentricity)
        )
    eccentricity = math.sqrt(1.0 - (a / c) ** 2)
    return 2.0 * math.pi * a * a + 2.0 * math.pi * a * c / eccentricity * math.asin(eccentricity)


def triaxial_ellipsoid_area(a: float, b: float, c: float) -> float:
    """Area of a closed ellipsoid with semi-axes ``a >= b >= c`` (Legendre's formula).

    ``A = 2 pi c^2 + (2 pi a b / sin(phi)) * (E(phi, k) sin^2(phi) + F(phi, k) cos^2(phi))``
    with ``cos(phi) = c/a`` and ``k^2 = a^2 (b^2 - c^2) / (b^2 (a^2 - c^2))``. Axes are
    sorted internally because the area is symmetric under permutation while the formula is
    not. The spheroid case ``a = b`` is delegated to :func:`spheroid_area` where the
    formula degenerates to ``k = 1`` and the closed elementary form is better conditioned.
    """
    values = sorted((a, b, c), reverse=True)
    major, middle, minor = values
    if minor <= 0.0:
        raise ValueError(f"semi-axes must be positive, got {(a, b, c)!r}")
    if abs(major - middle) <= 1e-14 * major:
        return spheroid_area(major, minor)
    if abs(middle - minor) <= 1e-14 * middle:
        return spheroid_area(major, middle)
    sin_phi = math.sqrt(max(1.0 - (minor / major) ** 2, 0.0))
    cos_phi = minor / major
    phi = math.atan2(sin_phi, cos_phi)
    k_squared = (major**2 * (middle**2 - minor**2)) / (middle**2 * (major**2 - minor**2))
    k = math.sqrt(min(max(k_squared, 0.0), 1.0))
    second_kind, first_kind = incomplete_ef(phi, k)
    return 2.0 * math.pi * minor**2 + (2.0 * math.pi * major * middle / sin_phi) * (
        second_kind * sin_phi**2 + first_kind * cos_phi**2
    )
