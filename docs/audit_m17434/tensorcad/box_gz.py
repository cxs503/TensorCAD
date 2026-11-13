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
from collections.abc import Sequence
from dataclasses import dataclass
from itertools import pairwise

__all__ = [
    "DEFAULT_ANGLES_DEG",
    "MAX_HEEL_DEG",
    "BoxGzArea",
    "BoxGzError",
    "BoxGzPoint",
    "box_gz",
    "box_gz_area",
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


# --------------------------------------------------------------------------------------
# areas under the curve: the quantity the stability codes actually regulate
# --------------------------------------------------------------------------------------

#: Five-point Gauss-Legendre abscissae and weights on ``[-1, 1]``.  The rule integrates
#: polynomials up to degree nine exactly; on a smooth piece that is a ``2**10`` error ratio
#: against the same rule on doubled panels, which is what makes the doubling estimate below
#: a bound rather than a wish.
_GAUSS5_NODES = (
    0.0,
    -0.5384693101056831,
    0.5384693101056831,
    -0.9061798459386640,
    0.9061798459386640,
)
_GAUSS5_WEIGHTS = (
    0.5688888888888889,
    0.4786286704993665,
    0.4786286704993665,
    0.2369268850561891,
    0.2369268850561891,
)
_GAUSS5_GAIN = 2**10  # error ratio between the n-panel and the 2n-panel rule
_MIN_PANELS = 4
_MAX_PANELS = 512


def _area_heel(name: str, value: float) -> float:
    """A heel the area integral may touch -- the model's span, not the clip's."""
    if not isinstance(value, (int, float)) or not math.isfinite(value):
        raise BoxGzError(f"{name} must be a finite number of degrees: {value!r}")
    if abs(value) > MAX_HEEL_DEG:
        raise BoxGzError(
            f"{name} {value!r} is outside the model span "
            f"-{MAX_HEEL_DEG} .. {MAX_HEEL_DEG} degrees"
        )
    return float(value)


def _kn(
    *, breadth: float, draft: float, depth: float, heel_deg: float
) -> float:
    """``KN(phi) = Y_B(phi)`` from the exact section, without the ``KG`` term."""
    return box_section(breadth=breadth, draft=draft, depth=depth, heel_deg=heel_deg)[2]


def _kn_area_piece(
    *,
    breadth: float,
    draft: float,
    depth: float,
    start_rad: float,
    stop_rad: float,
    tolerance: float,
) -> tuple[float, float, int, int]:
    """``(integral, error estimate, panels, evaluations)`` over one smooth piece.

    Composite five-point Gauss-Legendre with panel doubling.  The estimate is the classic
    Richardson difference divided by the rule's known gain, so it is only as honest as the
    smoothness of the piece it is given -- which is why the caller splits at the regime
    boundaries before calling.
    """
    panels = _MIN_PANELS
    evaluations = 0
    previous: float | None = None
    value = 0.0
    error = math.inf
    while True:
        width = (stop_rad - start_rad) / panels
        half = 0.5 * width
        total = 0.0
        for index in range(panels):
            centre = start_rad + (index + 0.5) * width
            for node, weight in zip(_GAUSS5_NODES, _GAUSS5_WEIGHTS, strict=True):
                total += weight * _kn(
                    breadth=breadth,
                    draft=draft,
                    depth=depth,
                    heel_deg=math.degrees(centre + half * node),
                )
        value = total * half
        evaluations += len(_GAUSS5_NODES) * panels
        if previous is not None:
            error = abs(value - previous) / (_GAUSS5_GAIN - 1)
            if error <= max(tolerance * abs(value), 1e-16 * max(abs(value), 1.0)):
                break
            if panels >= _MAX_PANELS:
                break
        previous = value
        panels *= 2
    return value, error, panels, evaluations


@dataclass(frozen=True)
class BoxGzArea:
    """``int GZ dphi`` between two heels, with the split that produced it.

    ``area = kn_area - kg_term``, where ``kg_term = KG (cos(start) - cos(stop))`` is exact.
    ``quadrature_error`` is the summed panel-doubling estimate over the smooth pieces, and
    ``breakpoints_deg`` are the regime boundaries the interval was split at, because
    ``KN'(phi)`` steps there: the waterplane width jumps when a deck edge enters the water.
    """

    start_deg: float
    stop_deg: float
    area: float
    kn_area: float
    kg_term: float
    quadrature_error: float
    breakpoints_deg: tuple[float, ...]
    panels: int
    heel_evaluations: int


def box_gz_area(
    *,
    breadth: float,
    draft: float,
    depth: float,
    kg: float,
    start_deg: float = 0.0,
    stop_deg: float = 30.0,
    tolerance: float = 1e-12,
) -> BoxGzArea:
    """Area under the exact righting-arm curve, in m rad, over ``[start_deg, stop_deg]``.

    This is the number the stability codes regulate -- the 0-30, 0-40 and 30-40 degree areas of
    IS Code 2008 part A s. 2.2, and the gust area of the weather criterion -- and for a barge
    it is the one that used to stop at the wall-sided limit, where the closed-form
    antiderivative ``GM0 (1 - cos) + BM0/2 (sec + cos - 2)`` ceases to represent the ship.
    Nothing here is bounded by that window.

    The interval is cut at the regime boundaries of :func:`box_regime_bounds_deg` and each
    smooth piece is integrated separately.  That is not cosmetic: ``d(KN)/dphi`` is the
    waterplane second moment over the displaced volume, and the waterplane width *steps* when a
    deck edge enters the water, so a single high-order rule spanning that heel falls like
    ``h**2.5`` instead of ``h**10`` (measured on the doc barge, ``int KN dphi`` over 0-40 deg:
    4.0e-8 with 8 panels and 6.6e-12 with 256, against 8 panels *per piece* and round-off with
    the split).  ``tolerance`` is the relative target of the panel-doubling estimate, never a
    claim about the geometry -- the section is exact, the quadrature is not.

    >>> result = box_gz_area(breadth=2.0, draft=0.8, depth=1.2, kg=0.55, stop_deg=40.0)
    >>> round(result.area, 8)
    0.06703591
    >>> result.quadrature_error < 1e-14
    True
    >>> [round(angle, 4) for angle in result.breakpoints_deg]
    [21.8014]
    """
    breadth, draft, depth = _geometry(breadth, draft, depth)
    if not isinstance(kg, (int, float)) or not math.isfinite(kg):
        raise BoxGzError(f"KG must be finite: {kg!r}")
    start = _area_heel("start heel", start_deg)
    stop = _area_heel("stop heel", stop_deg)
    if stop <= start:
        raise BoxGzError(f"empty heel interval: {start!r} to {stop!r} degrees")
    if not isinstance(tolerance, (int, float)) or not math.isfinite(tolerance) or tolerance <= 0:
        raise BoxGzError(f"tolerance must be a finite positive number: {tolerance!r}")

    bounds = box_regime_bounds_deg(breadth=breadth, draft=draft, depth=depth)
    # the geometry is symmetric, so every boundary has a mirror on the other side of upright;
    # missing the mirror is what makes an integral over a port interval creep instead of converge
    candidates = {
        side * angle
        for angle in bounds.values()
        if math.isfinite(angle)
        for side in (1.0, -1.0)
    }
    breakpoints = tuple(sorted(angle for angle in candidates if start < angle < stop))
    edges = (start, *breakpoints, stop)

    kn_area = 0.0
    error = 0.0
    panels = 0
    evaluations = 0
    for low, high in pairwise(edges):
        piece, piece_error, piece_panels, piece_evaluations = _kn_area_piece(
            breadth=breadth,
            draft=draft,
            depth=depth,
            start_rad=math.radians(low),
            stop_rad=math.radians(high),
            tolerance=tolerance,
        )
        kn_area += piece
        error += piece_error
        panels = max(panels, piece_panels)
        evaluations += piece_evaluations

    kg_term = float(kg) * (math.cos(math.radians(start)) - math.cos(math.radians(stop)))
    return BoxGzArea(
        start_deg=start,
        stop_deg=stop,
        area=kn_area - kg_term,
        kn_area=kn_area,
        kg_term=kg_term,
        quadrature_error=error,
        breakpoints_deg=breakpoints,
        panels=panels,
        heel_evaluations=evaluations,
    )
