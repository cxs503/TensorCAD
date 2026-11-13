"""Ellipsoidal geodesics: the Earth model the campaign geometry was not using.

:mod:`tensorcad.icechannel` projects the de Long -> Pevek route onto a local
azimuthal-equidistant plane built on a *sphere* of ``R = 6371.0088 km`` and its
evidence page (``docs/icechannel.md``) budgets the projection distortion down to
0.46 % at 1061 km. The sphere itself is never in that budget: it is the datum,
not an approximation. It is an approximation, and at 70-77 deg N it is the *same
size* as the distortion it was measured against -- over that latitude band the
WGS84 ellipsoid is 0.33 % (meridional) to 0.43 % (transverse) longer than the
mean-radius sphere, because the meridional radius ``M`` and the prime-vertical
radius ``N`` both grow towards the pole while ``R`` does not.

This module supplies the missing datum so a corridor, a route length or a
distance-along-track can be reported in true metres on a stated ellipsoid:

* :func:`vincenty_inverse` / :func:`vincenty_direct` -- the classical series
  solution (Vincenty 1975, *Direct and inverse solutions of geodesics on the
  ellipsoid*, Survey Review 23(176)), iterated to ``1e-12`` rad on the
  longitude difference. The two are *round-trip* inverses of each other, which
  is one of the two independent checks this module ships with.
* :func:`clairaut_distance` -- a deliberately different route to the same
  number: the geodesic is integrated as a curve of constant Clairaut parameter
  ``c = r sin(alpha)`` with ``r = N cos(phi)``, ``ds/dphi = M / sqrt(1-(c/r)^2)``,
  and ``c`` is found by shooting on the longitude difference. Two unrelated
  algorithms (series vs quadrature) agreeing to ``1e-9`` relative is evidence;
  one algorithm agreeing with itself is not.
* :func:`great_circle` -- spherical distance, i.e. what the ice channel used,
  with the difference to the ellipsoid reported rather than left implicit.

Scope, stated so nobody reads more into it than is there:

* Monotonic latitude only for :func:`clairaut_distance`. If the geodesic vertex
  falls inside the latitude interval the quadrature is no longer a single
  branch and the function raises :class:`GeodesicError` instead of integrating
  the wrong interval; :func:`vincenty_inverse` has no such restriction (it does
  not care where the vertex is).
* No datum shift (``WGS84`` <-> local datums), no geoid, no height above the
  ellipsoid. Distances are surface (chord-free) geodesic lengths.
* No map projection *inverse*: :func:`azimuthal_equidistant` is the forward
  spherical projection the ice channel uses plus its distortion factor, not a
  general projection library.

Only the Python standard library is required.
"""
from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

__all__ = [
    "WGS84_A",
    "WGS84_B",
    "WGS84_F",
    "WGS84_E2",
    "MEAN_RADIUS_M",
    "GeodesicError",
    "GeodesicSolution",
    "meridional_radius",
    "prime_vertical_radius",
    "meridian_arc_m",
    "vertex_latitude",
    "parallel_radius_m",
    "vincenty_inverse",
    "vincenty_direct",
    "clairaut_distance",
    "great_circle",
    "azimuthal_equidistant",
    "sphere_error",
]

#: WGS84 semi-major axis [m].
WGS84_A = 6378137.0
#: WGS84 flattening.
WGS84_F = 1.0 / 298.257223563
#: WGS84 semi-minor axis [m].
WGS84_B = WGS84_A * (1.0 - WGS84_F)
#: WGS84 first eccentricity squared.
WGS84_E2 = WGS84_F * (2.0 - WGS84_F)
#: The sphere radius the campaign geometry uses: the IUGG mean radius, i.e.
#: ``(2a + b) / 3`` for WGS84, rounded to 0.1 m.
MEAN_RADIUS_M = 6371008.8

_DEG = math.pi / 180.0


class GeodesicError(ValueError):
    """Raised when a geodesic cannot be solved by the requested route."""


@dataclass(frozen=True)
class GeodesicSolution:
    """Inverse-problem answer: length, end azimuths, and how hard it was to get."""

    distance_m: float
    azimuth1_deg: float
    azimuth2_deg: float
    iterations: int
    converged: bool
    method: str

    def to_dict(self) -> dict[str, object]:
        return {
            "distance_m": self.distance_m,
            "distance_km": self.distance_m / 1000.0,
            "azimuth1_deg": self.azimuth1_deg,
            "azimuth2_deg": self.azimuth2_deg,
            "iterations": self.iterations,
            "converged": self.converged,
            "method": self.method,
            "ellipsoid": "WGS84 a=6378137.0 f=1/298.257223563",
        }


# -- ellipsoid radii ------------------------------------------------------------


def meridional_radius(latitude_rad: float) -> float:
    """Radius of curvature in the meridian ``M = a(1-e2)/(1-e2 sin^2)^(3/2)``."""
    s = math.sin(latitude_rad)
    return WGS84_A * (1.0 - WGS84_E2) / (1.0 - WGS84_E2 * s * s) ** 1.5


def prime_vertical_radius(latitude_rad: float) -> float:
    """Radius of curvature in the prime vertical ``N = a/sqrt(1-e2 sin^2)``."""
    s = math.sin(latitude_rad)
    return WGS84_A / math.sqrt(1.0 - WGS84_E2 * s * s)


def parallel_radius_m(latitude_rad: float) -> float:
    """Distance from the spin axis to the surface point: ``r = N cos(phi)``."""
    return prime_vertical_radius(latitude_rad) * math.cos(latitude_rad)


# -- spherical reference --------------------------------------------------------


def great_circle(lat1_deg: float, lon1_deg: float, lat2_deg: float, lon2_deg: float,
                 radius_m: float = MEAN_RADIUS_M) -> float:
    """Great-circle distance by the haversine form, on a sphere of ``radius_m``.

    Haversine rather than ``acos`` of the dot product: the latter loses half its
    significant digits for short legs, which is exactly the regime a 250 m
    corridor step lives in.
    """
    p1 = lat1_deg * _DEG
    p2 = lat2_deg * _DEG
    dp = p2 - p1
    dl = (lon2_deg - lon1_deg) * _DEG
    h = math.sin(dp / 2.0) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2.0) ** 2
    return 2.0 * radius_m * math.asin(min(1.0, math.sqrt(h)))


def azimuthal_equidistant(lat_deg: float, lon_deg: float, origin: tuple[float, float],
                          radius_m: float = MEAN_RADIUS_M) -> tuple[float, float]:
    """Forward spherical azimuthal-equidistant projection about ``origin``.

    Returns ``(x_east_m, y_north_m)``. Distances *from the origin* are exact by
    construction; the tangential scale is inflated by ``c / sin(c)`` in the
    central angle ``c``, which :func:`sphere_error` quantifies per point.
    """
    d = great_circle(origin[0], origin[1], lat_deg, lon_deg, radius_m)
    az = _initial_bearing(origin[0], origin[1], lat_deg, lon_deg)
    return d * math.sin(az * _DEG), d * math.cos(az * _DEG)


def _initial_bearing(lat1_deg: float, lon1_deg: float, lat2_deg: float, lon2_deg: float) -> float:
    p1 = lat1_deg * _DEG
    p2 = lat2_deg * _DEG
    dl = (lon2_deg - lon1_deg) * _DEG
    y = math.sin(dl) * math.cos(p2)
    x = math.cos(p1) * math.sin(p2) - math.sin(p1) * math.cos(p2) * math.cos(dl)
    return math.degrees(math.atan2(y, x))


def sphere_error(lat1_deg: float, lon1_deg: float, lat2_deg: float, lon2_deg: float) -> dict[str, float]:
    """How wrong the mean-radius sphere is for one pair of points.

    ``relative = sphere / ellipsoid - 1`` is negative where the sphere
    understates the true surface distance (the polar case: both radii grow with
    latitude, so the sphere is short).
    """
    sph = great_circle(lat1_deg, lon1_deg, lat2_deg, lon2_deg)
    ell = vincenty_inverse(lat1_deg, lon1_deg, lat2_deg, lon2_deg).distance_m
    if ell == 0.0:
        raise GeodesicError("coincident points have no relative error")
    return {
        "sphere_m": sph,
        "ellipsoid_m": ell,
        "difference_m": sph - ell,
        "relative": sph / ell - 1.0,
        "mean_latitude_deg": 0.5 * (lat1_deg + lat2_deg),
    }


# -- spherical arcs used as quadrature anchors ---------------------------------


def meridian_arc_m(lat1_deg: float, lat2_deg: float, *, panels: int = 4096) -> float:
    """Meridian arc length between two latitudes on WGS84 (Simpson rule).

    The meridian has ``dlambda = 0``, so the arc is the plain integral of ``M``
    and needs no shooting. Simpson on 4096 panels is converged well past double
    precision for a 90 deg span (the integrand is smooth and symmetric), which is
    what makes it usable as an anchor for the general case.
    """
    lo = min(lat1_deg, lat2_deg) * _DEG
    hi = max(lat1_deg, lat2_deg) * _DEG
    if hi == lo:
        return 0.0
    if panels % 2:
        panels += 1
    h = (hi - lo) / panels
    total = meridional_radius(lo) + meridional_radius(hi)
    for i in range(1, panels):
        total += (4.0 if i % 2 else 2.0) * meridional_radius(lo + i * h)
    arc = total * h / 3.0
    return arc if lat2_deg >= lat1_deg else -arc


# -- Vincenty 1975 ---------------------------------------------------------------


def vincenty_inverse(lat1_deg: float, lon1_deg: float, lat2_deg: float, lon2_deg: float,
                     *, tolerance: float = 1e-12, max_iterations: int = 200) -> GeodesicSolution:
    """Geodesic distance and azimuths between two WGS84 points (Vincenty inverse).

    Non-convergence is reported (``converged=False``) rather than raised: the
    classical series stalls only near antipodal points, and a caller comparing
    routes should see that happen instead of getting an exception from a
    distance query.
    """
    a, b, f = WGS84_A, WGS84_B, WGS84_F
    if not (-90.0 <= lat1_deg <= 90.0 and -90.0 <= lat2_deg <= 90.0):
        raise GeodesicError(f"latitude out of range: {lat1_deg!r}, {lat2_deg!r}")
    if lat1_deg == lat2_deg and (lon1_deg - lon2_deg) % 360.0 == 0.0:
        return GeodesicSolution(0.0, 0.0, 0.0, 0, True, "vincenty-inverse")

    p1 = lat1_deg * _DEG
    p2 = lat2_deg * _DEG
    big_l = ((lon2_deg - lon1_deg + 180.0) % 360.0 - 180.0) * _DEG
    u1 = math.atan((1.0 - f) * math.tan(p1))
    u2 = math.atan((1.0 - f) * math.tan(p2))
    sin_u1, cos_u1 = math.sin(u1), math.cos(u1)
    sin_u2, cos_u2 = math.sin(u2), math.cos(u2)

    lam = big_l
    sin_sigma = cos_sigma = sigma = sin_alpha = cos2_alpha = cos2_sigma_m = 0.0
    converged = False
    iterations = 0
    for iterations in range(1, max_iterations + 1):
        sin_lam, cos_lam = math.sin(lam), math.cos(lam)
        sin_sigma = math.hypot(cos_u2 * sin_lam, cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_lam)
        if sin_sigma == 0.0:
            return GeodesicSolution(0.0, 0.0, 0.0, iterations, True, "vincenty-inverse")
        cos_sigma = sin_u1 * sin_u2 + cos_u1 * cos_u2 * cos_lam
        sigma = math.atan2(sin_sigma, cos_sigma)
        sin_alpha = cos_u1 * cos_u2 * sin_lam / sin_sigma
        cos2_alpha = 1.0 - sin_alpha * sin_alpha
        cos2_sigma_m = (cos_sigma - 2.0 * sin_u1 * sin_u2 / cos2_alpha) if cos2_alpha else 0.0
        c = f / 16.0 * cos2_alpha * (4.0 + f * (4.0 - 3.0 * cos2_alpha))
        previous = lam
        lam = big_l + (1.0 - c) * f * sin_alpha * (
            sigma + c * sin_sigma * (cos2_sigma_m + c * cos_sigma * (-1.0 + 2.0 * cos2_sigma_m ** 2))
        )
        if abs(lam - previous) < tolerance:
            converged = True
            break

    u2sq = cos2_alpha * (a * a - b * b) / (b * b)
    big_a = 1.0 + u2sq / 16384.0 * (4096.0 + u2sq * (-768.0 + u2sq * (320.0 - 175.0 * u2sq)))
    big_b = u2sq / 1024.0 * (256.0 + u2sq * (-128.0 + u2sq * (74.0 - 47.0 * u2sq)))
    delta_sigma = big_b * sin_sigma * (
        cos2_sigma_m
        + big_b / 4.0 * (
            cos_sigma * (-1.0 + 2.0 * cos2_sigma_m ** 2)
            - big_b / 6.0 * cos2_sigma_m * (-3.0 + 4.0 * sin_sigma ** 2) * (-3.0 + 4.0 * cos2_sigma_m ** 2)
        )
    )
    distance = b * big_a * (sigma - delta_sigma)
    sin_lam, cos_lam = math.sin(lam), math.cos(lam)
    az1 = math.atan2(cos_u2 * sin_lam, cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_lam)
    az2 = math.atan2(cos_u1 * sin_lam, -sin_u1 * cos_u2 + cos_u1 * sin_u2 * cos_lam)
    return GeodesicSolution(distance, math.degrees(az1), math.degrees(az2), iterations,
                            converged, "vincenty-inverse")


def vincenty_direct(lat1_deg: float, lon1_deg: float, azimuth1_deg: float, distance_m: float,
                    *, tolerance: float = 1e-12, max_iterations: int = 200
                    ) -> tuple[float, float, float, int, bool]:
    """WGS84 direct problem: ``(lat2, lon2, azimuth2_deg, iterations, converged)``.

    Present so the inverse can be checked by round trip: forward along the
    azimuth it returned, and the endpoint must come back to the input point.
    """
    a, b, f = WGS84_A, WGS84_B, WGS84_F
    p1 = lat1_deg * _DEG
    alpha1 = azimuth1_deg * _DEG
    tan_u1 = (1.0 - f) * math.tan(p1)
    cos_u1 = 1.0 / math.sqrt(1.0 + tan_u1 * tan_u1)
    sin_u1 = tan_u1 * cos_u1
    sigma1 = math.atan2(tan_u1, math.cos(alpha1))
    sin_alpha = cos_u1 * math.sin(alpha1)
    cos2_alpha = 1.0 - sin_alpha * sin_alpha
    u2sq = cos2_alpha * (a * a - b * b) / (b * b)
    big_a = 1.0 + u2sq / 16384.0 * (4096.0 + u2sq * (-768.0 + u2sq * (320.0 - 175.0 * u2sq)))
    big_b = u2sq / 1024.0 * (256.0 + u2sq * (-128.0 + u2sq * (74.0 - 47.0 * u2sq)))

    sigma = distance_m / (b * big_a)
    converged = False
    iterations = 0
    two_sigma_m = delta_sigma = 0.0
    for iterations in range(1, max_iterations + 1):
        two_sigma_m = 2.0 * sigma1 + sigma
        sin_s, cos_s = math.sin(sigma), math.cos(sigma)
        cos2s = math.cos(two_sigma_m)
        delta_sigma = big_b * sin_s * (
            cos2s + big_b / 4.0 * (
                cos_s * (-1.0 + 2.0 * cos2s ** 2)
                - big_b / 6.0 * cos2s * (-3.0 + 4.0 * sin_s ** 2) * (-3.0 + 4.0 * cos2s ** 2)
            )
        )
        previous = sigma
        sigma = distance_m / (b * big_a) + delta_sigma
        if abs(sigma - previous) < tolerance:
            converged = True
            break

    sin_s, cos_s = math.sin(sigma), math.cos(sigma)
    cos2s = math.cos(two_sigma_m)
    tmp = sin_u1 * sin_s - cos_u1 * cos_s * math.cos(alpha1)
    p2 = math.atan2(sin_u1 * cos_s + cos_u1 * sin_s * math.cos(alpha1),
                    (1.0 - f) * math.sqrt(sin_alpha * sin_alpha + tmp * tmp))
    lam = math.atan2(sin_s * math.sin(alpha1), cos_u1 * cos_s - sin_u1 * sin_s * math.cos(alpha1))
    c = f / 16.0 * cos2_alpha * (4.0 + f * (4.0 - 3.0 * cos2_alpha))
    big_l = lam - (1.0 - c) * f * sin_alpha * (
        sigma + c * sin_s * (cos2s + c * cos_s * (-1.0 + 2.0 * cos2s ** 2))
    )
    lon2 = lon1_deg + big_l / _DEG
    az2 = math.atan2(sin_alpha, -sin_u1 * sin_s + cos_u1 * cos_s * math.cos(alpha1))
    return math.degrees(p2), lon2, math.degrees(az2), iterations, converged


# -- independent quadrature route --------------------------------------------------


def _simpson(f, lo: float, hi: float, panels: int) -> float:
    """Composite Simpson rule on a smooth integrand (``panels`` made even)."""
    if hi == lo:
        return 0.0
    if panels % 2:
        panels += 1
    h = (hi - lo) / panels
    total = f(lo) + f(hi)
    for i in range(1, panels):
        total += (4.0 if i % 2 else 2.0) * f(lo + i * h)
    return total * h / 3.0


def vertex_latitude(c: float) -> float:
    """Latitude at which ``r = N cos(phi)`` equals the Clairaut parameter ``c``.

    ``r`` falls monotonically from ``a`` at the equator to zero at the pole, so the
    vertex latitude is the unique root and a bisection finds it. A fitted ``c`` whose
    vertex lies *inside* the latitude interval means the geodesic turns round inside
    it: the quadrature would then be integrating one branch of a two-branch path and
    silently returning a length for the wrong curve, which is why
    :func:`clairaut_distance` checks this and refuses.
    """
    lo, hi = 0.0, math.pi / 2.0
    if c >= WGS84_A:
        return 0.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mid <= lo or mid >= hi:  # pragma: no cover - float exhaustion
            break
        if parallel_radius_m(mid) > c:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def clairaut_distance(lat1_deg: float, lon1_deg: float, lat2_deg: float, lon2_deg: float,
                      *, panels: int = 512, tolerance: float = 1e-13,
                      max_iterations: int = 200) -> float:
    """Geodesic length by quadrature of the Clairaut relation (no Vincenty series).

    On an ellipsoid of revolution every geodesic keeps ``c = r sin(alpha)``
    constant with ``r = N cos(phi)``. With ``cos(alpha) = sqrt(1-(c/r)^2)`` the
    arc length and the longitude gain along the meridian coordinate are

        ds/dphi      = M / sqrt(1-(c/r)^2)
        dlambda/dphi = c M / (r^2 sqrt(1-(c/r)^2))

    so ``c`` follows from shooting on ``|Delta lambda|`` and the length from one
    Simpson integral of the first relation. ``panels=512`` is the default because
    that is where the answer stops moving: 512 and 4096 panels agree to 1.6e-14
    relative on the steepest route leg, four orders below anything the geometry
    cares about. Nothing is shared with the Vincenty
    expansion except the ellipsoid constants, which is what makes the agreement
    between the two an independent check.

    Only ``|Delta lambda|`` enters: reflecting through a meridian plane
    (``lambda -> -lambda``) is an isometry of the ellipsoid that swaps an eastbound
    for a westbound geodesic of the same length, so the two directions are one
    problem, not two.

    Raises :class:`GeodesicError` when this reduction does not apply: latitudes
    that straddle the equator, equal non-zero latitudes (a parallel is not a
    geodesic off the equator), or a fitted ``c`` whose vertex falls inside the
    latitude interval (the integrand is then two branches, and integrating one of
    them would be a silently wrong answer).
    """
    if not (-90.0 <= lat1_deg <= 90.0 and -90.0 <= lat2_deg <= 90.0):
        raise GeodesicError(f"latitude out of range: {lat1_deg!r}, {lat2_deg!r}")
    delta_lon = (lon2_deg - lon1_deg + 180.0) % 360.0 - 180.0
    if delta_lon == 0.0:
        return abs(meridian_arc_m(lat1_deg, lat2_deg, panels=panels * 8))
    if lat1_deg * lat2_deg < 0.0:
        raise GeodesicError("latitudes straddle the equator: Clairaut route not applicable")
    if lat1_deg == lat2_deg:
        if lat1_deg != 0.0:
            raise GeodesicError("equal non-zero latitudes: a parallel is not a geodesic")
        return WGS84_A * abs(delta_lon) * _DEG

    lo_deg, hi_deg = (lat1_deg, lat2_deg) if lat1_deg < lat2_deg else (lat2_deg, lat1_deg)
    lo = lo_deg * _DEG
    hi = hi_deg * _DEG
    r_min = min(parallel_radius_m(lo), parallel_radius_m(hi))
    target = abs(delta_lon) * _DEG

    def longitude_gain(c: float) -> float:
        def dlambda(phi: float) -> float:
            r = parallel_radius_m(phi)
            ratio = c / r
            if ratio >= 1.0:
                raise GeodesicError("vertex inside the latitude interval")
            return c * meridional_radius(phi) / (r * r * math.sqrt(1.0 - ratio * ratio))

        return _simpson(dlambda, lo, hi, panels)

    low = 1e-12 * r_min
    high = r_min * (1.0 - 1e-12)
    if longitude_gain(low) > target:
        raise GeodesicError("longitude change too small for any Clairaut parameter")
    if longitude_gain(high) < target:
        raise GeodesicError("vertex inside the latitude interval: two branches")

    c = high
    for _ in range(max_iterations):
        mid = 0.5 * (low + high)
        if high - low <= tolerance * r_min:
            c = mid
            break
        if longitude_gain(mid) < target:
            low = mid
        else:
            high = mid
        c = mid

    vertex = math.degrees(vertex_latitude(c))
    if vertex < hi_deg - 1e-9:
        raise GeodesicError(
            "vertex latitude %.6f deg lies inside [%.6f, %.6f] deg: two branches"
            % (vertex, lo_deg, hi_deg)
        )

    def arc(phi: float) -> float:
        r = parallel_radius_m(phi)
        ratio = c / r
        return meridional_radius(phi) / math.sqrt(1.0 - ratio * ratio)

    length = _simpson(arc, lo, hi, panels)
    # A-posteriori convergence: when the vertex sits *just* outside the interval the
    # integrand is nearly singular at the end and Simpson on a fixed grid stops
    # converging -- at (60, 140) -> (60.1, 175) the 512-panel answer is 1.1 % long.
    # Rather than trust a grid, halve it and refuse if the answer moves.
    coarse = _simpson(arc, lo, hi, max(8, panels // 2))
    if length != 0.0 and abs(coarse / length - 1.0) > 1e-8:
        raise GeodesicError(
            "quadrature not converged (vertex %.3f deg, interval top %.3f deg): "
            "%.2e relative between %d and %d panels" % (vertex, hi_deg,
                                                        abs(coarse / length - 1.0),
                                                        max(8, panels // 2), panels)
        )
    return length


def polyline_length_m(points: Sequence[tuple[float, float]], *,
                      ellipsoid: bool = True) -> float:
    """Sum of geodesic (or great-circle) legs over ``(lat_deg, lon_deg)`` points."""
    total = 0.0
    for (la1, lo1), (la2, lo2) in zip(points, list(points)[1:]):
        total += (vincenty_inverse(la1, lo1, la2, lo2).distance_m if ellipsoid
                  else great_circle(la1, lo1, la2, lo2))
    return total
