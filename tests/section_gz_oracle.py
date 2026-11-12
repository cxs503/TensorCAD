"""Independent references for :mod:`tensorcad.section_gz`.

``tensorcad.section_gz`` computes exact righting arms of a prismatic body two ways: a convex
polygon is rotated and clipped by the waterline (Sutherland-Hodgman) and the clipped area and
centroid come from Green's theorem; an ellipse uses the closed form of the unit disk cut by a
line.  This module re-derives the same quantities by a *different* method so the two can be
compared without sharing a code path:

* **slices, not clipping.**  For a fixed heel the section is an ellipse (or a convex polygon)
  in the world frame; the horizontal chord at height ``zeta`` has a half width ``w(zeta)/2``
  and a chord centre ``y_c(zeta)``.  The submerged area is ``int w dzeta`` and the first
  moment is ``int w y_c dzeta``.  Nothing is clipped, no polygon vertex list is rebuilt, and
  Green's theorem is never used.
* **the conic is found numerically, not analytically.**  For the ellipse the quadric is
  recovered from the quadratic form ``(y/a)^2 + ((z-b)/b)^2`` evaluated at six points, then
  the chord geometry is read off the conic coefficients.  The disk formula
  ``pi/2 + asin(d) + d sqrt(1-d^2)`` used by the library never appears here.
* **the integrals are quadrature, not closed form.**  ``int w dzeta`` is evaluated by
  composite Simpson after the substitution ``zeta = z_c + H sin(theta)``, which removes the
  square-root end behaviour of the chord.  1024 panels put the quadrature error around
  ``1e-11`` relative on an integral of order one -- three orders below the agreement the
  checks below need, and far below the platform's 3 % criterion.
* **polygons are integrated exactly, not sampled.**  Between consecutive vertex heights the
  chord width and the chord centre are both linear in ``zeta``, so the area is a trapezoid and
  the first moments are exact with a three-point Simpson rule (exact for a quadratic).

The heel convention is the library's (see :func:`tensorcad.section_gz`): the section point
``(y, z)`` maps to ``(y cos(phi) + z sin(phi), -y sin(phi) + z cos(phi))``, ``KN`` is the
submerged centroid's ``y``, and ``GZ = KN - KG sin(phi)`` with ``KG`` above the keel.

Only the Python standard library is used.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

Point2 = tuple[float, float]

#: Bisection budget for the waterline level.  Independent of the library's, generous enough
#: that the level is solved to the neighbouring float.
_LEVEL_ITERATIONS = 120

#: Simpson panels for the conic integrals.
_SIMPSON_PANELS = 1024


def _simpson(function, low: float, high: float, panels: int = _SIMPSON_PANELS) -> float:
    """Composite Simpson on ``[low, high]``; ``high == low`` returns zero."""
    span = high - low
    if span == 0.0:
        return 0.0
    if panels % 2:
        panels += 1
    step = span / panels
    total = function(low) + function(high)
    for index in range(1, panels, 2):
        total += 4.0 * function(low + index * step)
    for index in range(2, panels, 2):
        total += 2.0 * function(low + index * step)
    return total * step / 3.0


def _bisect_level(area_of_level, target: float, low: float, high: float) -> float:
    """Solve ``area_of_level(level) = target`` for an area that increases with the level."""
    if not (area_of_level(low) <= target <= area_of_level(high)):
        raise ValueError(
            f"displaced area {target!r} is outside [{area_of_level(low)!r}, "
            f"{area_of_level(high)!r}] on {low!r}..{high!r}"
        )
    for _ in range(_LEVEL_ITERATIONS):
        middle = 0.5 * (low + high)
        if middle <= low or middle >= high:
            break
        if area_of_level(middle) < target:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


class EllipseSliceOracle:
    """One heeled elliptical section, described by its horizontal chords.

    The section is the ellipse of semi-axes ``semi_beam`` (y) and ``semi_height`` (z) whose
    centre sits at ``(0, semi_height)`` (keel at ``z = 0``), heeled by ``heel_rad``.
    """

    def __init__(self, semi_beam: float, semi_height: float, heel_rad: float) -> None:
        a, b = float(semi_beam), float(semi_height)
        if not (a > 0.0 and b > 0.0):
            raise ValueError("semi-axes must be positive")
        cosine, sine = math.cos(heel_rad), math.sin(heel_rad)

        def body(y: float, z: float) -> tuple[float, float]:
            """World -> body, the inverse of the library's heel rotation."""
            return (cosine * y - sine * z, sine * y + cosine * z)

        def quadric(y: float, z: float) -> float:
            body_y, body_z = body(y, z)
            return (body_y / a) ** 2 + ((body_z - b) / b) ** 2

        at_origin = quadric(0.0, 0.0)
        at_y = quadric(1.0, 0.0)
        at_minus_y = quadric(-1.0, 0.0)
        at_z = quadric(0.0, 1.0)
        at_minus_z = quadric(0.0, -1.0)
        at_yyz = quadric(1.0, 1.0)
        m11 = 0.5 * (at_y + at_minus_y) - at_origin
        m22 = 0.5 * (at_z + at_minus_z) - at_origin
        l1 = 0.5 * (at_y - at_minus_y)
        l2 = 0.5 * (at_z - at_minus_z)
        m12 = 0.5 * (at_yyz - m11 - m22 - l1 - l2 - at_origin)
        if not (m11 > 0.0):
            raise ValueError("the conic is not an ellipse with a horizontal chord")
        det = m11 * m22 - m12 * m12
        if not (det > 0.0):
            raise ValueError("the conic is degenerate")
        # Q(y, z) <= 1 with a fixed z is m11 y^2 + b_q y + c_q = 0;  the discriminant is a
        # downward parabola in z whose roots bracket the section and whose width gives the
        # chord:  D(z) = -4 det z^2 + 4 P z + C.
        p_coeff = m12 * l1 - m11 * l2
        c_coeff = l1 * l1 - 4.0 * m11 * (at_origin - 1.0)
        root = math.sqrt(max(0.0, p_coeff * p_coeff + det * c_coeff))
        self.semi_beam = a
        self.semi_height = b
        self.heel_rad = heel_rad
        self.m11, self.m12 = m11, m12
        self.l1 = l1
        self.z_centre = p_coeff / (2.0 * det)
        self.half_height = root / (2.0 * det)
        # w(zeta) = _scale * sqrt(1 - u^2) and w dzeta = _scale * H cos^2(theta) dtheta
        self._scale = 2.0 * self.half_height * math.sqrt(det) / m11

    @property
    def area(self) -> float:
        return math.pi * self.semi_beam * self.semi_height

    def _theta(self, level: float) -> float:
        u = (level - self.z_centre) / self.half_height
        return math.asin(max(-1.0, min(1.0, u)))

    def width(self, level: float) -> float:
        """Horizontal chord length at world height ``level``."""
        u = (level - self.z_centre) / self.half_height
        if abs(u) >= 1.0:
            return 0.0
        return self._scale * math.sqrt(max(0.0, 1.0 - u * u))

    def chord_centre(self, level: float) -> float:
        """World ``y`` of the midpoint of that chord."""
        return -(2.0 * self.m12 * level + self.l1) / (2.0 * self.m11)

    def area_below(self, level: float) -> float:
        theta = self._theta(level)
        integral = _simpson(lambda t: math.cos(t) ** 2, -0.5 * math.pi, theta)
        return self._scale * self.half_height * integral

    def moment_y_below(self, level: float) -> float:
        """``int y dA`` below ``level`` (first moment about the keel centreline plane)."""
        theta = self._theta(level)
        half, centre, scale = self.half_height, self.z_centre, self._scale

        def integrand(t: float) -> float:
            zeta = centre + half * math.sin(t)
            return scale * half * math.cos(t) ** 2 * self.chord_centre(zeta)

        return _simpson(integrand, -0.5 * math.pi, theta)

    def moment_z_below(self, level: float) -> float:
        """``int z dA`` below ``level``."""
        theta = self._theta(level)
        half, centre, scale = self.half_height, self.z_centre, self._scale

        def integrand(t: float) -> float:
            zeta = centre + half * math.sin(t)
            return scale * half * math.cos(t) ** 2 * zeta

        return _simpson(integrand, -0.5 * math.pi, theta)

    def properties_below(self, level: float) -> tuple[float, float, float]:
        """``(area, centroid_y, centroid_z)`` below ``level``."""
        area = self.area_below(level)
        if area <= 0.0:
            return 0.0, 0.0, 0.0
        return (
            area,
            self.moment_y_below(level) / area,
            self.moment_z_below(level) / area,
        )


def ellipse_prism_arms_independent(
    semi_beam: float,
    semi_height: float,
    angles_deg: Sequence[float],
    *,
    kg: float,
    waterline_z: float,
) -> list[float]:
    """Righting arms of an elliptical prism, by slice quadrature instead of the disk formula."""
    target = EllipseSliceOracle(semi_beam, semi_height, 0.0).area_below(waterline_z)
    if not (0.0 < target <= math.pi * semi_beam * semi_height):
        raise ValueError("the upright waterline does not cut the section")
    arms: list[float] = []
    for angle in angles_deg:
        heel_rad = math.radians(float(angle))
        oracle = EllipseSliceOracle(semi_beam, semi_height, heel_rad)
        level = _bisect_level(
            oracle.area_below,
            target,
            oracle.z_centre - oracle.half_height,
            oracle.z_centre + oracle.half_height,
        )
        area, centre_y, _ = oracle.properties_below(level)
        arms.append(centre_y - kg * math.sin(heel_rad))
    return arms


def _chord(section: Sequence[Point2], level: float) -> tuple[float, float]:
    """``(width, centre_y)`` of the horizontal chord of a convex polygon at ``level``."""
    crossings: list[float] = []
    count = len(section)
    for index in range(count):
        y0, z0 = section[index]
        y1, z1 = section[(index + 1) % count]
        if z0 == z1:
            if z0 == level:
                crossings.append(y0)
                crossings.append(y1)
            continue
        low, high = (z0, z1) if z0 < z1 else (z1, z0)
        if low <= level <= high:
            fraction = (level - z0) / (z1 - z0)
            crossings.append(y0 + fraction * (y1 - y0))
    if len(crossings) < 2:
        return 0.0, 0.0
    return max(crossings) - min(crossings), 0.5 * (max(crossings) + min(crossings))


def polygon_properties_below(
    section: Sequence[Point2], level: float
) -> tuple[float, float, float]:
    """``(area, centroid_y, centroid_z)`` of a convex polygon below ``level``, exactly.

    The section is integrated slice by slice between consecutive vertex heights: the chord
    width and the chord centre are both linear in ``zeta`` there, so a trapezoid gives the
    area and a three-point Simpson rule (exact for a quadratic) gives the two first moments.
    """
    heights = sorted({z for (_, z) in section})
    bottom, top = heights[0], heights[-1]
    ceiling = min(level, top)
    if ceiling <= bottom:
        return 0.0, 0.0, 0.0
    cuts = sorted({z for z in heights if bottom < z < ceiling} | {ceiling})
    area = 0.0
    moment_y = 0.0
    moment_z = 0.0
    lower = bottom
    for upper in cuts:
        if upper <= lower:
            continue
        width_low, centre_low = _chord(section, lower)
        width_high, centre_high = _chord(section, upper)
        middle = 0.5 * (lower + upper)
        width_mid, centre_mid = _chord(section, middle)
        span = upper - lower
        area += 0.5 * (width_low + width_high) * span
        weight = span / 6.0
        moment_y += weight * (
            width_low * centre_low + 4.0 * width_mid * centre_mid + width_high * centre_high
        )
        moment_z += weight * (
            width_low * lower + 4.0 * width_mid * middle + width_high * upper
        )
        lower = upper
    if area <= 0.0:
        return 0.0, 0.0, 0.0
    return area, moment_y / area, moment_z / area


def polygon_prism_arms_independent(
    section: Sequence[Point2],
    angles_deg: Sequence[float],
    *,
    kg: float,
    waterline_z: float,
) -> list[float]:
    """Righting arms of a convex-polygon prism, by slice integration instead of clipping."""
    points = tuple((float(y), float(z)) for (y, z) in section)
    keel = min(z for (_, z) in points)
    shifted = tuple((y, z - keel) for (y, z) in points)
    base_level = waterline_z - keel
    target = polygon_properties_below(shifted, base_level)[0]
    if target <= 0.0:
        raise ValueError("the upright waterline does not cut the section")
    arms: list[float] = []
    for angle in angles_deg:
        heel_rad = math.radians(float(angle))
        cosine, sine = math.cos(heel_rad), math.sin(heel_rad)
        heeled = tuple((y * cosine + z * sine, -y * sine + z * cosine) for (y, z) in shifted)
        heights = [z for (_, z) in heeled]
        level = _bisect_level(
            lambda value: polygon_properties_below(heeled, value)[0],
            target,
            min(heights),
            max(heights),
        )
        _, centre_y, _ = polygon_properties_below(heeled, level)
        arms.append(centre_y - kg * math.sin(heel_rad))
    return arms


def wall_sided_arm(breadth: float, draft: float, kg: float, heel_deg: float) -> float:
    """The published wall-sided righting arm ``sin(phi) (GM + BM tan^2(phi) / 2)``.

    ``KB = draft/2``, ``BM = breadth^2 / (12 draft)`` and ``GM = KB + BM - KG``: the standard
    box-barge result (Biran, *Ship Hydrostatics and Stability*, eq. 4.62; Rawson & Tupper,
    *Basic Ship Theory*), transcribed here rather than read from the library.
    """
    heel_rad = math.radians(float(heel_deg))
    kb = 0.5 * draft
    bm = breadth * breadth / (12.0 * draft)
    gm = kb + bm - kg
    return math.sin(heel_rad) * (gm + 0.5 * bm * math.tan(heel_rad) ** 2)


def wall_sided_limit_deg(breadth: float, draft: float) -> float:
    """Heel at which the waterline leaves the vertical sides of a box of the given draft.

    The freeboard is taken equal to the draft (the box used in the checks is twice as deep as
    it is immersed), so the last wall-sided heel solves ``tan(phi) = draft / (breadth/2)``.
    """
    return math.degrees(math.atan(draft / (0.5 * breadth)))


def waterplane_gm(semi_beam_waterline: float, submerged_area: float, kb: float, kg: float) -> float:
    """``GM = KB + BM - KG`` with ``BM = I_T / A`` and ``I_T = 2 y_w^3 / 3`` per unit length."""
    inertia = 2.0 * semi_beam_waterline**3 / 3.0
    return kb + inertia / submerged_area - kg