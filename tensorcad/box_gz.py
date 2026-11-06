"""Exact righting arms for a rectangular barge at *every* heel, not just the wall-sided window.

Why this module exists
----------------------
:func:`tensorcad.gz.box_righting_arms_analytic` is exact and one line, but it only holds while
the heeled waterline still cuts both side walls. Past that heel the waterline reaches a bilge or
a deck edge and the closed form stops representing the ship -- the module says so and refuses:

    >>> box_righting_arms_analytic(breadth=2.0, draft=0.8, kg=0.55, angles_deg=(40.0,),
    ...                            depth=1.2)
    Traceback (most recent call last):
    ...
    tensorcad.gz.GzError: heel 40.0 exceeds the wall-sided limit 21.801 degrees: ...

That refusal is right, and it is also where large-angle stability starts to matter: a barge's
range of stability, its angle of vanishing stability and the area under ``GZ`` all live past
21.8 degrees. What was missing was not a solver -- :mod:`tensorcad.gz` and
:mod:`tensorcad.righting_arm` both integrate a heeled mesh -- but an *independent* answer out
there, because two integrators run on the same mesh agree even when both are modelling the same
wrong thing.

The closed form generalises.  A rectangular barge's section is a rectangle
``|y| <= B/2, 0 <= z <= D`` and, at heel ``phi``, the immersed part is that rectangle cut by a
straight waterline.  A straight line cuts a rectangle into a polygon whose area, and whose first
moments, are elementary: the *shoelace* sums over the clipped polygon are exact in exact
arithmetic, and the only unknown is where the waterline sits.  Holding the displaced volume at
its upright value fixes it: ``A(c) = B * T`` where ``A(c)`` is that exact area and ``c`` the
waterline height.  ``A`` is strictly increasing, so the inverse is unique, and

    dA/dc = the length of the waterline chord inside the section

gives Newton's method an exact derivative.  Nothing here is discretised: what
:mod:`tensorcad.gz` calls "skipped discretisation error" is, for this geometry, zero.

The four regimes are worth naming because the *shape* of ``GZ(phi)`` changes at each boundary,
and because a reference is only trustworthy if it is continuous across them:

===============  ===========================================  ==============================
regime           waterline cuts                               entered when
===============  ===========================================  ==============================
``wall_sided``   both side walls                              ``phi = 0``
``deck_immersed``  the deck and the low side wall             ``tan(phi) -> 2(D-T)/B``
``bilge_emerged``  the bottom and the high side wall         ``tan(phi) -> 2T/B``
``bilge_and_deck`` the bottom and the deck                   both of the above
===============  ===========================================  ==============================

Frames
------
The section is written in the body frame (keel at ``z = 0``, centreline at ``y = 0``) and
rotated into the world frame by the same convention as :mod:`tensorcad.righting_arm`,
``Y = y cos(phi) + z sin(phi)``, ``Z = -y sin(phi) + z cos(phi)``, so the arithmetic below is
``KN(phi) = Y_B(phi)`` and ``GZ = KN - KG sin(phi)``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

__all__ = [
    "BoxGzError",
    "BoxGzPoint",
    "DEFAULT_ANGLES_DEG",
    "MAX_HEEL_DEG",
    "box_gz",
    "box_gz_curve",
    "box_regime_bounds_deg",
    "box_section",
    "box_vanishing_angle_deg",
]

#: Above this heel the section is nearly a line and the "waterline" language stops paying; the
#: clip itself stays exact, but a user asking for 90 degrees wants a different model.
MAX_HEEL_DEG = 89.0

DEFAULT_ANGLES_DEG = tuple(float(angle) for angle in range(0, 90, 5))


class BoxGzError(ValueError):
    """Raised for a barge that is not a rectangle, or for a heel outside the model."""


@dataclass(frozen=True)
class BoxGzPoint:
    """One point of the exact curve, with the geometry that produced it."""

    heel_deg: float
    regime: str
    righting_arm: float
    kn: float
    waterline_z: float
    submerged_area: float
    centre_of_buoyancy: tuple[float, float]


def _require_positive(name: str, value: float) -> float:
    if not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0.0:
        raise BoxGzError(f"{name} must be a finite positive number: {value!r}")
    return float(value)


def _geometry(breadth: float, draft: float, depth: float) -> tuple[float, float, float]:
    breadth = _require_positive("breadth", breadth)
    draft = _require_positive("draft", draft)
    depth = _require_positive("depth", depth)
    if depth <= draft:
        raise BoxGzError(f"depth {depth!r} must exceed draft {draft!r}")
    return breadth, draft, depth


def _heel_radians(heel_deg: float) -> float:
    if not isinstance(heel_deg, (int, float)) or not math.isfinite(heel_deg):
        raise BoxGzError(f"heel must be a finite number of degrees: {heel_deg!r}")
    if abs(heel_deg) >= MAX_HEEL_DEG + 1.0:
        raise BoxGzError(
            f"heel {heel_deg!r} is at or beyond 90 degrees: the section degenerates to a line"
        )
    return math.radians(heel_deg)


def _rotated_section(breadth: float, depth: float, phi: float) -> list[tuple[float, float]]:
    """Section rectangle corners in the world frame, counter-clockwise."""
    half = 0.5 * breadth
    cos_heel, sin_heel = math.cos(phi), math.sin(phi)
    corners = ((-half, 0.0), (half, 0.0), (half, depth), (-half, depth))
    return [(y * cos_heel + z * sin_heel, -y * sin_heel + z * cos_heel) for (y, z) in corners]


def _clip_below(
    polygon: Sequence[tuple[float, float]], height: float
) -> list[tuple[float, float]]:
    """Sutherland-Hodgman clip of a convex polygon against ``Z <= height``."""
    out: list[tuple[float, float]] = []
    count = len(polygon)
    for index in range(count):
        current = polygon[index]
        following = polygon[(index + 1) % count]
        delta_current = current[1] - height
        delta_following = following[1] - height
        if delta_current <= 0.0:
            out.append(current)
        if (delta_current < 0.0 < delta_following) or (delta_following < 0.0 < delta_current):
            fraction = delta_current / (delta_current - delta_following)
            out.append(
                (
                    current[0] + fraction * (following[0] - current[0]),
                    height,
                )
            )
    return out


def _polygon_moments(
    polygon: Sequence[tuple[float, float]]
) -> tuple[float, float, float]:
    """``(area, int Y dA, int Z dA)`` of an implicitly closed polygon, by shoelace sums."""
    area2 = 0.0
    first_y = 0.0
    first_z = 0.0
    count = len(polygon)
    for index in range(count):
        y1, z1 = polygon[index]
        y2, z2 = polygon[(index + 1) % count]
        cross = y1 * z2 - y2 * z1
        area2 += cross
        first_y += (y1 + y2) * cross
        first_z += (z1 + z2) * cross
    return (0.5 * area2, first_y / 6.0, first_z / 6.0)


def _chord_at(polygon: Sequence[tuple[float, float]], height: float) -> float:
    """Length of the horizontal chord of a convex polygon at ``Z = height`` (>= 0)."""
    ys: list[float] = []
    count = len(polygon)
    for index in range(count):
        y1, z1 = polygon[index]
        y2, z2 = polygon[(index + 1) % count]
        if z1 == z2:
            if z1 == height:
                ys.extend((y1, y2))
            continue
        if min(z1, z2) <= height <= max(z1, z2):
            fraction = (height - z1) / (z2 - z1)
            ys.append(y1 + fraction * (y2 - y1))
    if not ys:
        return 0.0
    return max(ys) - min(ys)


def _solve_waterline(
    polygon: Sequence[tuple[float, float]], target_area: float, *, tolerance: float = 0.0
) -> tuple[float, float, float, float]:
    """Vertical position of the waterline for a given immersed area.

    Returns ``(height, area, int Y dA, int Z dA)``.  Bisection brackets the answer on the
    polygon's vertical extent (``A`` is strictly increasing there), then Newton steps on the
    exact derivative ``dA/dh = chord(h)`` polish it.  The bracket loop exits when the interval
    no longer resolves in double precision, so the returned area equals ``target_area`` to
    round-off rather than to a tolerance anyone had to choose.
    """
    low = min(z for _, z in polygon)
    high = max(z for _, z in polygon)
    total_area = _polygon_moments(polygon)[0]
    if target_area >= total_area:
        raise BoxGzError(
            f"the section cannot hold {target_area!r}: the whole rectangle is {total_area!r}"
        )
    for _ in range(200):
        middle = 0.5 * (low + high)
        if middle <= low or middle >= high:
            break
        if _polygon_moments(_clip_below(polygon, middle))[0] < target_area:
            low = middle
        else:
            high = middle
    height = 0.5 * (low + high)
    area, first_y, first_z = _polygon_moments(_clip_below(polygon, height))
    for _ in range(3):  # Newton polish; the derivative is the exact chord length
        chord = _chord_at(polygon, height)
        if chord <= 1e-15:
            break
        height += (target_area - area) / chord
        area, first_y, first_z = _polygon_moments(_clip_below(polygon, height))
        if abs(area - target_area) <= tolerance:
            break
    return height, area, first_y, first_z


def _waterline_limits(
    phi: float, breadth: float, height: float
) -> tuple[float, float]:
    """Waterline heights at the two side walls, lower first.

    Writing the water plane in the body frame, ``z = height/cos(phi) + y tan(phi)``, the two
    walls sit at ``intercept -/+ (B/2) tan(phi)``.  Sorting them matters: at negative heel
    ``tan(phi)`` changes sign and the walls swap roles, so using the signed expression without
    sorting labels a port heel with the starboard regime (this was a real bug here).
    """
    cos_heel = math.cos(phi)
    if abs(cos_heel) < 1e-12:
        raise BoxGzError("the section degenerates to a line at 90 degrees of heel")
    offset = 0.5 * breadth * abs(math.tan(phi))
    intercept = height / cos_heel
    return intercept - offset, intercept + offset


def _regime(phi: float, breadth: float, depth: float, height: float) -> str:
    """Which pair of section edges the waterline cuts, from the body-frame intercept."""
    lower, upper = _waterline_limits(phi, breadth, height)
    scale = 1e-12 * max(abs(height), breadth, depth)
    bilge_emerged = lower < -scale
    deck_immersed = upper > depth + scale
    if bilge_emerged and deck_immersed:
        return "bilge_and_deck"
    if bilge_emerged:
        return "bilge_emerged"
    if deck_immersed:
        return "deck_immersed"
    return "wall_sided"


def box_section(
    *,
    breadth: float,
    draft: float,
    depth: float,
    heel_deg: float,
) -> tuple[float, float, float, float]:
    """``(height, area, Y_B, Z_B)`` of the heeled rectangle holding ``breadth * draft``.

    The reference itself: no mesh, no discretisation, one clip and one waterline solve.
    """
    breadth, draft, depth = _geometry(breadth, draft, depth)
    phi = _heel_radians(heel_deg)
    polygon = _rotated_section(breadth, depth, phi)
    target = breadth * draft
    height, area, first_y, first_z = _solve_waterline(polygon, target)
    if area <= 0.0:
        raise BoxGzError("the clip produced an empty section")
    return height, area, first_y / area, first_z / area


def box_gz(
    *,
    breadth: float,
    draft: float,
    depth: float,
    kg: float,
    heel_deg: float,
) -> BoxGzPoint:
    """Exact ``GZ`` for one heel, valid across all four section regimes."""
    breadth, draft, depth = _geometry(breadth, draft, depth)
    phi = _heel_radians(heel_deg)
    if not isinstance(kg, (int, float)) or not math.isfinite(kg):
        raise BoxGzError(f"KG must be finite: {kg!r}")
    height, area, y_b, z_b = box_section(
        breadth=breadth, draft=draft, depth=depth, heel_deg=heel_deg
    )
    # ``Y_B`` is already the world-frame abscissa of the buoyancy centre, so it *is* KN:
    # rotating it again would be the classic double-rotation error (it used to be here).
    kn = y_b
    return BoxGzPoint(
        heel_deg=float(heel_deg),
        regime=_regime(phi, breadth, depth, height),
        righting_arm=kn - float(kg) * math.sin(phi),
        kn=kn,
        waterline_z=height,
        submerged_area=area,
        centre_of_buoyancy=(y_b, z_b),
    )


def box_gz_curve(
    *,
    breadth: float,
    draft: float,
    depth: float,
    kg: float,
    angles_deg: Sequence[float] = DEFAULT_ANGLES_DEG,
) -> tuple[BoxGzPoint, ...]:
    """The whole curve, one :class:`BoxGzPoint` per heel."""
    heels = tuple(float(angle) for angle in angles_deg)
    if not heels:
        raise BoxGzError("at least one heel angle is needed")
    return tuple(
        box_gz(breadth=breadth, draft=draft, depth=depth, kg=kg, heel_deg=heel)
        for heel in heels
    )


def box_regime_bounds_deg(
    *, breadth: float, draft: float, depth: float
) -> dict[str, float]:
    """Heels at which a box leaves the wall-sided regime, and where it leaves both walls.

    ``deck_immersion`` is where the low-side deck edge touches the waterline,
    ``bilge_emergence`` where the high-side bilge leaves it, and ``both`` the heel at which the
    section is cut by the bottom and the deck at once.  Each is a root of a monotone residual in
    ``phi`` and is bracketed to machine precision.
    """
    breadth, draft, depth = _geometry(breadth, draft, depth)

    def residual(heel_deg: float, which: str) -> float:
        height, _, _, _ = box_section(
            breadth=breadth, draft=draft, depth=depth, heel_deg=heel_deg
        )
        lower, upper = _waterline_limits(math.radians(heel_deg), breadth, height)
        if which == "deck_immersion":
            return upper - depth
        if which == "bilge_emergence":
            return -lower
        return min(upper - depth, -lower)

    def bisect(which: str) -> float:
        low, high = 0.0, MAX_HEEL_DEG
        if residual(low, which) >= 0.0:
            return 0.0
        if residual(high, which) < 0.0:
            return math.inf
        for _ in range(120):
            middle = 0.5 * (low + high)
            if middle <= low or middle >= high:
                break
            if residual(middle, which) < 0.0:
                low = middle
            else:
                high = middle
        return 0.5 * (low + high)

    return {
        "deck_immersion": bisect("deck_immersion"),
        "bilge_emergence": bisect("bilge_emergence"),
        "both_walls_lost": bisect("both"),
    }


def box_vanishing_angle_deg(
    *, breadth: float, draft: float, depth: float, kg: float
) -> float:
    """Heel at which ``GZ`` returns to zero (the range of stability), or ``inf``.

    A barge that is too stiff or too heavy never vanishes below ``89`` degrees; that is
    reported as ``math.inf`` rather than as a fabricated angle.
    """
    breadth, draft, depth = _geometry(breadth, draft, depth)

    def gz(heel_deg: float) -> float:
        return box_gz(
            breadth=breadth, draft=draft, depth=depth, kg=kg, heel_deg=heel_deg
        ).righting_arm

    low, high = 0.0, MAX_HEEL_DEG
    if gz(1e-9) <= 0.0:  # GZ(0) is identically zero; the sign that matters is the slope
        return 0.0
    if gz(high) > 0.0:
        return math.inf
    for _ in range(120):
        middle = 0.5 * (low + high)
        if middle <= low or middle >= high:
            break
        if gz(middle) > 0.0:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)