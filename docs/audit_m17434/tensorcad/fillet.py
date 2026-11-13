"""Fillet and chamfer of polygonal sections, with exact arc-aware moments.

``GAP-CAD-012`` in ``docs/oss_inventory_cad/GAP-ANCHORS.md`` records that the package could
build a hull from parametric sections but had **no way to round a corner**: every edge in
:mod:`tensorcad.hull`, :mod:`tensorcad.parametric_hull`, :mod:`tensorcad.series60` and
:mod:`tensorcad.section_gz` is hard, and a zero-width stem is a modelling decision rather than
a radius.  A designer who wants a bilge radius, a rounded deck edge or a bracket toe gets no
geometry, no failure mode and no way to check the result.

This module adds the operation and -- more importantly for a verification-first package -- an
**exact** account of what it did to the section.  Rounding a corner of radius ``r`` and turn
angle ``phi`` (the exterior angle, positive on a convex corner) replaces two straight lengths
``t = r tan(|phi|/2)`` by a circular arc of sweep ``|phi|``.  Everything then follows in closed
form:

* area change ``dA = -sign(phi) r^2 (tan(|phi|/2) - |phi|/2)``
  (a convex corner loses material, a re-entrant corner gains it),
* perimeter change ``dP = r |phi| - 2 r tan(|phi|/2)`` (negative: a fillet shortens the
  outline),
* the arc's contribution to area, centroid and second moments comes from Green's theorem with
  the arc parameterised in closed form, so the exact values carry no discretisation at all.

:class:`SectionShape` therefore carries **two independent values of the area**: the
Green's-theorem sum over the exact path (lines plus arcs) and the per-corner closed form above.
They are derived differently and both are stored; ``area_crosscheck_error`` is their relative
difference, which is round-off (measured ~1e-16, see ``docs/fillet_chamfer.md``).

Feasibility is a result, not an exception path bolted on.  An edge of length ``L`` between
corners with turn angles ``phi1``, ``phi2`` can host at most
``L / (tan(|phi1|/2) + tan(|phi2|/2))``.  :func:`fillet_section` refuses a larger radius and
names the binding edge; ``allow_clamp=True`` clamps to the exact limit, which is *not* a
degenerate geometry -- at ``r = r_max`` the two fillet arcs meet tangentially on the vanished
edge.  A square of side ``2 r`` is the bubble test: clamping yields the exact circle.

Sections may be convex or re-entrant (chamfers likewise).  Solids are swept from section to
section and the size may vary linearly along the sweep axis.  The area of a filleted section is
*exactly* quadratic in the size, so the exact volume is
``A0 L + K L (r0^2 + r0 r1 + r1^2) / 3``.  A swept fillet surface stays tangent to the adjacent
flat faces for any radius profile, because the tangency direction lies in that flat face; the
module checks this and reports ``g1_error_deg``.

**What this is not.**  No B-rep, no radius that varies *around* one section, no blend of two
fillet surfaces where three faces meet (no corner sphere), no rolling-ball fillet over a
triangle mesh, no OCP backend.  It is a 2.5-D kernel with an exact oracle, which is what the
section and hull builders in this package actually consume.

Only the Python standard library is used.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

from .geometry import Mesh, MeshBuilder
from .mesh import (
    boundary_edges,
    edge_face_counts,
    remove_degenerate_faces,
    total_signed_volume,
    weld,
    zero_area_faces,
)

SCHEMA_FILLET = "tensorcad.geometry.fillet_section/1"
SCHEMA_FILLET_SOLID = "tensorcad.geometry.fillet_solid/1"

Point2 = tuple[float, float]

#: Below this turn a vertex is collinear: no arc, no length consumed, no area change.
_COLLINEAR_TURN = 1e-12


class FilletError(ValueError):
    """Raised when a section, a size or a sweep cannot be built."""


# ------------------------------------------------------------------ exact moments


@dataclass(frozen=True)
class SectionMoments:
    """Exact area properties of a closed section outline.

    ``second_u`` is ``int u**2 dA`` (second moment about the ``v`` axis) and ``second_v`` is
    ``int v**2 dA``; together they give the polar moment.  For a waterplane, ``int u**2 dA`` is
    the transverse inertia ``I_T`` behind ``BM = I_T / volume``.
    """

    area: float
    perimeter: float
    centroid: Point2
    second_u: float
    second_v: float

    @property
    def polar_second_moment(self) -> float:
        return self.second_u + self.second_v

    def to_dict(self) -> dict[str, float]:
        return {
            "area": self.area,
            "perimeter": self.perimeter,
            "centroid_u": self.centroid[0],
            "centroid_v": self.centroid[1],
            "second_u": self.second_u,
            "second_v": self.second_v,
            "polar_second_moment": self.polar_second_moment,
        }


def _line_terms(u0: float, v0: float, u1: float, v1: float) -> tuple[float, ...]:
    """Green's-theorem contributions of the straight segment ``(u0,v0) -> (u1,v1)``."""
    du = u1 - u0
    dv = v1 - v0
    area = dv * (u0 + du / 2.0)
    mu = dv / 2.0 * (u0 * u0 + u0 * du + du * du / 3.0)
    mv = -du / 2.0 * (v0 * v0 + v0 * dv + dv * dv / 3.0)
    iu = dv / 3.0 * (u0**3 + 1.5 * u0 * u0 * du + u0 * du * du + du**3 / 4.0)
    iv = -du / 3.0 * (v0**3 + 1.5 * v0 * v0 * dv + v0 * dv * dv + dv**3 / 4.0)
    return area, mu, mv, iu, iv


def _arc_terms(cu: float, cv: float, r: float, a0: float, sweep: float) -> tuple[float, ...]:
    """Green's-theorem contributions of an arc, integrated in closed form.

    The arc is ``(cu + r cos a, cv + r sin a)`` for ``a`` in ``[a0, a0 + sweep]``.  Every
    integrand is a polynomial of degree at most four in ``cos``/``sin``, so the contribution is
    elementary and exact.
    """
    a1 = a0 + sweep
    sin0, sin1 = math.sin(a0), math.sin(a1)
    cos0, cos1 = math.cos(a0), math.cos(a1)
    s1 = sin1 - sin0
    c1 = cos1 - cos0
    s2 = math.sin(2.0 * a1) - math.sin(2.0 * a0)
    s4 = math.sin(4.0 * a1) - math.sin(4.0 * a0)
    i_cos = s1
    i_cos2 = sweep / 2.0 + s2 / 4.0
    i_cos3 = s1 - (sin1**3 - sin0**3) / 3.0
    i_cos4 = 3.0 * sweep / 8.0 + s2 / 4.0 + s4 / 32.0
    i_sin = -c1
    i_sin2 = sweep / 2.0 - s2 / 4.0
    i_sin3 = -c1 + (cos1**3 - cos0**3) / 3.0
    i_sin4 = 3.0 * sweep / 8.0 - s2 / 4.0 + s4 / 32.0
    area = cu * r * i_cos + r * r * i_cos2
    mu = r / 2.0 * (cu * cu * i_cos + 2.0 * cu * r * i_cos2 + r * r * i_cos3)
    mv = r / 2.0 * (cv * cv * i_sin + 2.0 * cv * r * i_sin2 + r * r * i_sin3)
    iu = r / 3.0 * (
        cu**3 * i_cos + 3.0 * cu * cu * r * i_cos2 + 3.0 * cu * r * r * i_cos3 + r**3 * i_cos4
    )
    iv = r / 3.0 * (
        cv**3 * i_sin + 3.0 * cv * cv * r * i_sin2 + 3.0 * cv * r * r * i_sin3 + r**3 * i_sin4
    )
    return area, mu, mv, iu, iv


#: A path element: ``("line", u0, v0, u1, v1)`` or ``("arc", cu, cv, r, a0, sweep)``.
_Element = tuple


def _path_moments(elements: Sequence[_Element], perimeter: float) -> SectionMoments:
    """Accumulate exact moments over a closed path, evaluated about a local origin.

    Shifting the origin to the mean of the path's vertices keeps ``int u**2 dA`` away from the
    cancellation that ``u**2 - 2 u c + c**2`` suffers when the section sits far from the
    coordinate origin; the shift is undone exactly at the end.
    """
    us: list[float] = []
    vs: list[float] = []
    for element in elements:
        if element[0] == "line":
            us.extend((element[1], element[3]))
            vs.extend((element[2], element[4]))
        else:
            cu, cv, r = element[1], element[2], element[3]
            for angle in (element[4], element[4] + element[5]):
                us.append(cu + r * math.cos(angle))
                vs.append(cv + r * math.sin(angle))
    if not us:
        raise FilletError("the section path is empty")
    su = sum(us) / len(us)
    sv = sum(vs) / len(vs)
    area = mu = mv = iu = iv = 0.0
    for element in elements:
        if element[0] == "line":
            da, dmu, dmv, diu, div = _line_terms(
                element[1] - su, element[2] - sv, element[3] - su, element[4] - sv
            )
        else:
            da, dmu, dmv, diu, div = _arc_terms(
                element[1] - su, element[2] - sv, element[3], element[4], element[5]
            )
        area += da
        mu += dmu
        mv += dmv
        iu += diu
        iv += div
    if not area > 0.0:
        raise FilletError(f"the section path is not positively oriented: area={area!r}")
    mu_local = mu
    mv_local = mv
    mu = mu_local + su * area
    mv = mv_local + sv * area
    iu += 2.0 * su * mu_local + su * su * area
    iv += 2.0 * sv * mv_local + sv * sv * area
    return SectionMoments(
        area=area,
        perimeter=perimeter,
        centroid=(mu / area, mv / area),
        second_u=iu,
        second_v=iv,
    )


# ------------------------------------------------------------------ validation


def _finite(value: object) -> bool:
    return isinstance(value, (int, float)) and math.isfinite(float(value))


def _clean_ring(points: Sequence[Sequence[float]], *, what: str = "section") -> tuple[Point2, ...]:
    """Validate a closed outline and return it counter-clockwise, without duplicate points."""
    if len(points) < 3:
        raise FilletError(f"a {what} needs at least 3 points, got {len(points)}")
    cleaned: list[Point2] = []
    for index, point in enumerate(points):
        try:
            u, v = point
        except (TypeError, ValueError) as exc:
            raise FilletError(f"point {index} of the {what} is not a pair: {point!r}") from exc
        if not (_finite(u) and _finite(v)):
            raise FilletError(f"point {index} of the {what} is not finite: {point!r}")
        candidate = (float(u), float(v))
        if cleaned and candidate == cleaned[-1]:
            continue
        cleaned.append(candidate)
    if len(cleaned) > 1 and cleaned[0] == cleaned[-1]:
        cleaned.pop()
    if len(cleaned) < 3:
        raise FilletError(f"the {what} collapses to fewer than 3 distinct points")
    area = _signed_area(cleaned)
    if area == 0.0:
        raise FilletError(f"the {what} is degenerate (zero area)")
    if area < 0.0:
        cleaned.reverse()
    _reject_self_intersection(cleaned, what)
    return tuple(cleaned)


def _signed_area(points: Sequence[Point2]) -> float:
    total = 0.0
    count = len(points)
    for index in range(count):
        u0, v0 = points[index]
        u1, v1 = points[(index + 1) % count]
        total += u0 * v1 - u1 * v0
    return total / 2.0


def _orientation(a: Point2, b: Point2, c: Point2) -> float:
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def _on_segment(p: Point2, q: Point2, r: Point2) -> bool:
    """True when ``p`` lies on the closed segment ``qr`` (the three are collinear)."""
    return min(q[0], r[0]) <= p[0] <= max(q[0], r[0]) and min(q[1], r[1]) <= p[1] <= max(q[1], r[1])


def _segments_cross(a: Point2, b: Point2, c: Point2, d: Point2) -> bool:
    """True when the two closed segments share a point that is not a shared endpoint."""
    d1 = _orientation(c, d, a)
    d2 = _orientation(c, d, b)
    d3 = _orientation(a, b, c)
    d4 = _orientation(a, b, d)
    if ((d1 > 0.0 > d2) or (d1 < 0.0 < d2)) and ((d3 > 0.0 > d4) or (d3 < 0.0 < d4)):
        return True
    for p, q, r in ((a, c, d), (b, c, d), (c, a, b), (d, a, b)):
        if _orientation(q, r, p) == 0.0 and _on_segment(p, q, r) and p not in (q, r):
            return True
    return False


def _reject_self_intersection(points: Sequence[Point2], what: str) -> None:
    count = len(points)
    for i in range(count):
        a, b = points[i], points[(i + 1) % count]
        for j in range(i + 1, count):
            if (i + 1) % count == j or (j + 1) % count == i:
                continue
            c, d = points[j], points[(j + 1) % count]
            if _segments_cross(a, b, c, d):
                raise FilletError(
                    f"the {what} is not a simple polygon: edges {i} and {j} intersect"
                )


def section_turns(points: Sequence[Sequence[float]]) -> tuple[float, ...]:
    """Exterior turn angle at each vertex of a counter-clockwise simple polygon.

    Positive on a convex corner (the outline turns left), negative on a re-entrant one and
    ``0.0`` on a collinear vertex.  The turn angle is ``pi - interior angle``.
    """
    ring = _clean_ring(points)
    return _turns_of(ring)


def _turns_of(ring: Sequence[Point2]) -> tuple[float, ...]:
    count = len(ring)
    turns: list[float] = []
    for index in range(count):
        previous = ring[index - 1]
        current = ring[index]
        following = ring[(index + 1) % count]
        in_u, in_v = current[0] - previous[0], current[1] - previous[1]
        out_u, out_v = following[0] - current[0], following[1] - current[1]
        cross = in_u * out_v - in_v * out_u
        dot = in_u * out_u + in_v * out_v
        turn = math.atan2(cross, dot)
        turns.append(0.0 if abs(turn) < _COLLINEAR_TURN else turn)
    return tuple(turns)


def _tangent_factor(turn: float) -> float:
    """``tan(|turn| / 2)``; unbounded as the interior angle approaches zero."""
    half = abs(turn) / 2.0
    if half >= math.pi / 2.0 - 1e-15:
        return math.inf
    return math.tan(half)


def max_fillet_radius(points: Sequence[Sequence[float]]) -> float:
    """Largest radius every corner can take without the fillets overrunning an edge."""
    ring = _clean_ring(points)
    return _max_radius_of(ring, _turns_of(ring))


def _max_radius_of(ring: Sequence[Point2], turns: Sequence[float]) -> float:
    count = len(ring)
    limit = math.inf
    for index in range(count):
        u0, v0 = ring[index]
        u1, v1 = ring[(index + 1) % count]
        length = math.hypot(u1 - u0, v1 - v0)
        demand = _tangent_factor(turns[index]) + _tangent_factor(turns[(index + 1) % count])
        if demand == 0.0:
            continue
        limit = min(limit, length / demand)
    return limit


def _min_edge_half(ring: Sequence[Point2]) -> float:
    count = len(ring)
    limit = math.inf
    for index in range(count):
        u0, v0 = ring[index]
        u1, v1 = ring[(index + 1) % count]
        limit = min(limit, math.hypot(u1 - u0, v1 - v0) / 2.0)
    return limit


def max_chamfer_setback(points: Sequence[Sequence[float]]) -> float:
    """Largest setback every corner can take without the chamfers overrunning an edge."""
    return _min_edge_half(_clean_ring(points))


def polygon_moments(points: Sequence[Sequence[float]]) -> SectionMoments:
    """Exact area, perimeter, centroid and second moments of a simple polygon."""
    ring = _clean_ring(points)
    perimeter = 0.0
    elements: list[_Element] = []
    count = len(ring)
    for index in range(count):
        u0, v0 = ring[index]
        u1, v1 = ring[(index + 1) % count]
        elements.append(("line", u0, v0, u1, v1))
        perimeter += math.hypot(u1 - u0, v1 - v0)
    return _path_moments(elements, perimeter)


# ------------------------------------------------------------------ shapes


@dataclass(frozen=True)
class CornerReport:
    """What happened at one vertex."""

    index: int
    vertex: Point2
    turn_deg: float
    set_back: float
    radius: float
    kind: str

    def to_dict(self) -> dict[str, object]:
        return {
            "index": self.index,
            "vertex": [self.vertex[0], self.vertex[1]],
            "turn_deg": self.turn_deg,
            "set_back": self.set_back,
            "radius": self.radius,
            "kind": self.kind,
        }


@dataclass(frozen=True)
class SectionShape:
    """A filleted (or chamfered) section: exact moments plus a mesh-ready outline.

    ``moments`` is exact -- Green's theorem over the real path, arcs included -- and
    ``area_crosscheck_error`` is its relative difference from the per-corner closed form, i.e.
    two independent derivations of the same area.  ``points`` is the same outline with each arc
    chorded into ``arc_points`` segments, which is what a mesh can use.

    ``g1_error_deg`` is the largest angle between a fillet's end tangent and the edge it joins;
    for a chamfer it is the crease angle the chamfer introduces (a chamfer is *not* G1, and
    saying so is the point).
    """

    kind: str
    points: tuple[Point2, ...]
    moments: SectionMoments
    corners: tuple[CornerReport, ...]
    requested: float
    achieved: float
    clamped: bool
    arc_points: int
    g1_error_deg: float
    area_crosscheck_error: float
    base_area: float
    area_coefficient: float

    def area_of(self, size: float) -> float:
        """Closed-form area of the same section built with another uniform size."""
        return self.base_area + self.area_coefficient * size * size

    def sampled_moments(self) -> SectionMoments:
        """Moments of the chorded outline -- the discretisation a mesh actually sees."""
        return polygon_moments(self.points)

    def to_dict(self) -> dict[str, object]:
        return {
            "schema": SCHEMA_FILLET,
            "kind": self.kind,
            "outline_points": len(self.points),
            "moments": self.moments.to_dict(),
            "sampled_moments": self.sampled_moments().to_dict(),
            "base_area": self.base_area,
            "area_coefficient": self.area_coefficient,
            "requested": self.requested,
            "achieved": self.achieved,
            "clamped": self.clamped,
            "arc_points": self.arc_points,
            "g1_error_deg": self.g1_error_deg,
            "area_crosscheck_error": self.area_crosscheck_error,
            "corners": [corner.to_dict() for corner in self.corners],
        }


def _resolve_size(
    turns: Sequence[float], limit: float, requested: float, allow_clamp: bool, kind: str
) -> tuple[float, bool]:
    if not _finite(requested) or requested < 0.0:
        raise FilletError(f"the {kind} size must be finite and non-negative: {requested!r}")
    if requested == 0.0:
        return 0.0, False
    if not any(abs(turn) > _COLLINEAR_TURN for turn in turns):
        raise FilletError(f"the section has no corner to {kind}")
    # A request that sits *on* the limit is what a designer writes when they mean the limit
    # (a polygon filleted by its own inradius, say), and the limit itself carries the rounding
    # of ``tan`` and ``hypot``; treat anything within 1e-9 relative as on the limit.
    if requested <= limit * (1.0 + 1e-9):
        return float(min(requested, limit)), requested > limit
    if not allow_clamp:
        raise FilletError(
            f"the requested {kind} size {requested!r} overruns an edge: the largest feasible "
            f"size is {limit!r} (pass allow_clamp=True to clamp to it)"
        )
    return float(limit), True


def _corner_geometry(
    previous: Point2, current: Point2, following: Point2, turn: float, size: float
) -> tuple[Point2, Point2, float, float, float, Point2]:
    """Tangent points, arc centre, start angle and outgoing edge direction at one corner."""
    in_u = current[0] - previous[0]
    in_v = current[1] - previous[1]
    in_len = math.hypot(in_u, in_v)
    out_u = following[0] - current[0]
    out_v = following[1] - current[1]
    out_len = math.hypot(out_u, out_v)
    in_u, in_v = in_u / in_len, in_v / in_len
    out_u, out_v = out_u / out_len, out_v / out_len
    set_back = size * _tangent_factor(turn)
    start = (current[0] - set_back * in_u, current[1] - set_back * in_v)
    end = (current[0] + set_back * out_u, current[1] + set_back * out_v)
    # The centre is one radius from the incoming edge, on the interior side for a convex
    # corner and on the void side for a re-entrant one.
    normal_u, normal_v = -in_v, in_u
    sign = 1.0 if turn > 0.0 else -1.0
    centre_u = start[0] + sign * size * normal_u
    centre_v = start[1] + sign * size * normal_v
    a0 = math.atan2(start[1] - centre_v, start[0] - centre_u)
    return start, end, centre_u, centre_v, a0, (out_u, out_v)


def _wrap(angle: float) -> float:
    return (angle + math.pi) % (2.0 * math.pi) - math.pi


def _corner_path(
    ring: Sequence[Point2],
    turns: Sequence[float],
    size: float,
    arc_points: int,
    *,
    kind: str,
    dense: bool,
) -> tuple[list[_Element], list[CornerReport], list[Point2], float, float]:
    """Path elements, corner reports, chorded outline, perimeter and G1 error of one section."""
    count = len(ring)
    starts: list[Point2] = []
    ends: list[Point2] = []
    arcs: list[_Element | None] = []
    corners: list[CornerReport] = []
    sampled: list[Point2] = []
    g1_error = 0.0
    chord_total = 0.0
    arc_total = 0.0
    for index in range(count):
        current = ring[index]
        turn = turns[index]
        collinear = abs(turn) < _COLLINEAR_TURN
        if size == 0.0 or collinear:
            starts.append(current)
            ends.append(current)
            arcs.append(None)
            if dense and not collinear:
                sampled.extend([current] * (arc_points + 1))
            else:
                sampled.append(current)
            corners.append(
                CornerReport(
                    index=index,
                    vertex=current,
                    turn_deg=math.degrees(turn),
                    set_back=0.0,
                    radius=0.0,
                    kind="convex" if turn > 0.0 else ("reflex" if turn < 0.0 else "collinear"),
                )
            )
            continue
        start, end, centre_u, centre_v, a0, _out = _corner_geometry(
            ring[index - 1], current, ring[(index + 1) % count], turn, size
        )
        starts.append(start)
        ends.append(end)
        in_u = current[0] - ring[index - 1][0]
        in_v = current[1] - ring[index - 1][1]
        in_len = math.hypot(in_u, in_v)
        if kind == "fillet":
            sweep = turn
            a1 = a0 + sweep
            end_angle = math.atan2(end[1] - centre_v, end[0] - centre_u)
            drift = max(
                abs(_wrap(a1 - end_angle)),
                math.hypot(
                    centre_u + size * math.cos(a1) - end[0],
                    centre_v + size * math.sin(a1) - end[1],
                ),
            )
            if drift > 1e-9 * max(1.0, size):
                raise FilletError(f"fillet construction drifted at vertex {index}: {drift!r}")
            arcs.append(("arc", centre_u, centre_v, size, a0, sweep))
            arc_total += size * abs(sweep)
            for step in range(arc_points + 1):
                angle = a0 + sweep * step / arc_points
                sampled.append(
                    (centre_u + size * math.cos(angle), centre_v + size * math.sin(angle))
                )
            # The arc's tangent at its start, for a sweep of either sign.
            direction = 1.0 if sweep > 0.0 else -1.0
            tangent = (-direction * math.sin(a0), direction * math.cos(a0))
            cos_error = min(1.0, abs((tangent[0] * in_u + tangent[1] * in_v) / in_len))
            g1_error = max(g1_error, math.degrees(math.acos(cos_error)))
            corner_radius = size
        else:
            arcs.append(("line", start[0], start[1], end[0], end[1]))
            sampled.extend((start, end))
            chord_u = end[0] - start[0]
            chord_v = end[1] - start[1]
            chord_len = math.hypot(chord_u, chord_v)
            chord_total += chord_len
            cos_crease = min(1.0, abs((chord_u * in_u + chord_v * in_v) / (chord_len * in_len)))
            g1_error = max(g1_error, math.degrees(math.acos(cos_crease)))
            corner_radius = 0.0
        corners.append(
            CornerReport(
                index=index,
                vertex=current,
                turn_deg=math.degrees(turn),
                set_back=size * _tangent_factor(turn),
                radius=corner_radius,
                kind="convex" if turn > 0.0 else "reflex",
            )
        )
    elements: list[_Element] = [arc for arc in arcs if arc is not None]
    perimeter = 0.0
    for index in range(count):
        u0, v0 = ends[index]
        u1, v1 = starts[(index + 1) % count]
        elements.append(("line", u0, v0, u1, v1))
        perimeter += math.hypot(u1 - u0, v1 - v0)
    perimeter += chord_total + arc_total
    return elements, corners, sampled, perimeter, g1_error


def _dedupe_ring(points: Sequence[Point2], scale: float) -> tuple[Point2, ...]:
    tolerance = 1e-12 * max(1.0, scale)
    cleaned: list[Point2] = []
    for point in points:
        if cleaned and math.hypot(
            point[0] - cleaned[-1][0], point[1] - cleaned[-1][1]
        ) <= tolerance:
            continue
        cleaned.append(point)
    while len(cleaned) > 2 and math.hypot(
        cleaned[0][0] - cleaned[-1][0], cleaned[0][1] - cleaned[-1][1]
    ) <= tolerance:
        cleaned.pop()
    if len(cleaned) < 3:
        raise FilletError("the filleted outline collapses to fewer than 3 points")
    return tuple(cleaned)


def _area_coefficient(turns: Sequence[float], kind: str) -> float:
    """``A(size) = base_area + coefficient * size**2`` -- exact for both operations."""
    total = 0.0
    for turn in turns:
        if abs(turn) < _COLLINEAR_TURN:
            continue
        sign = 1.0 if turn > 0.0 else -1.0
        if kind == "fillet":
            total -= sign * (_tangent_factor(turn) - abs(turn) / 2.0)
        else:
            total -= sign * math.sin(abs(turn)) / 2.0
    return total


def _build_shape(
    points: Sequence[Sequence[float]],
    size: float,
    *,
    kind: str,
    arc_points: int,
    allow_clamp: bool,
) -> SectionShape:
    if arc_points < 1:
        raise FilletError(f"arc_points must be at least 1, got {arc_points}")
    ring = _clean_ring(points)
    turns = _turns_of(ring)
    limit = _max_radius_of(ring, turns) if kind == "fillet" else _min_edge_half(ring)
    achieved, clamped = _resolve_size(turns, limit, size, allow_clamp, kind)
    elements, corners, sampled, perimeter, g1_error = _corner_path(
        ring, turns, achieved, arc_points, kind=kind, dense=False
    )
    moments = _path_moments(elements, perimeter)
    base_area = _signed_area(ring)
    coefficient = _area_coefficient(turns, kind)
    analytic = base_area + coefficient * achieved * achieved
    scale = max(1.0, max(abs(u) + abs(v) for u, v in ring))
    return SectionShape(
        kind=kind,
        points=_dedupe_ring(sampled, scale),
        moments=moments,
        corners=tuple(corners),
        requested=float(size),
        achieved=achieved,
        clamped=clamped,
        arc_points=arc_points,
        g1_error_deg=g1_error,
        area_crosscheck_error=abs(moments.area - analytic) / abs(analytic),
        base_area=base_area,
        area_coefficient=coefficient,
    )


def fillet_section(
    points: Sequence[Sequence[float]],
    radius: float,
    *,
    arc_points: int = 8,
    allow_clamp: bool = False,
) -> SectionShape:
    """Round every corner of a simple polygon with a common radius.

    ``radius=0`` returns the polygon itself (moments from the same Green's-theorem machinery).
    A radius that does not fit the shortest edge raises :class:`FilletError` naming the largest
    feasible radius, unless ``allow_clamp=True``, in which case the radius is clamped to that
    limit and ``clamped`` is set.
    """
    return _build_shape(
        points, radius, kind="fillet", arc_points=arc_points, allow_clamp=allow_clamp
    )


def chamfer_section(
    points: Sequence[Sequence[float]],
    setback: float,
    *,
    allow_clamp: bool = False,
) -> SectionShape:
    """Cut every corner of a simple polygon with a straight chamfer (``arc_points`` is 1)."""
    return _build_shape(
        points, setback, kind="chamfer", arc_points=1, allow_clamp=allow_clamp
    )


def rectangle_section(width: float, height: float, *, centred: bool = False) -> tuple[Point2, ...]:
    """Axis-aligned rectangle as a counter-clockwise outline (``u`` right, ``v`` up)."""
    if not (_finite(width) and _finite(height)) or width <= 0.0 or height <= 0.0:
        raise FilletError(f"a rectangle needs positive finite sides: {width!r} x {height!r}")
    if centred:
        hu, hv = width / 2.0, height / 2.0
        return ((-hu, -hv), (hu, -hv), (hu, hv), (-hu, hv))
    return ((0.0, 0.0), (width, 0.0), (width, height), (0.0, height))


# ------------------------------------------------------------------ solids


@dataclass(frozen=True)
class FilletSolid:
    """A section swept along ``+x`` with a size that may vary linearly, as a closed mesh.

    ``volume_exact`` uses the closed form ``A0 L + K L (r0^2 + r0 r1 + r1^2) / 3`` and
    ``volume_mesh`` integrates the actual triangles, so the ``mesh_volume_error`` in
    :meth:`to_dict` is the chord discretisation and nothing else.
    """

    kind: str
    mesh: Mesh
    length: float
    size_start: float
    size_end: float
    stations: int
    arc_points: int
    shape_start: SectionShape
    shape_end: SectionShape
    base_area: float
    area_coefficient: float
    volume_exact: float
    volume_mesh: float
    g1_error_deg: float
    sampled_crease_deg: float
    clamped: bool
    watertight: bool
    edge_manifold: bool
    euler_characteristic: int
    zero_area_faces: int

    @property
    def mesh_volume_error(self) -> float:
        return abs(self.volume_mesh - self.volume_exact) / abs(self.volume_exact)

    def to_dict(self) -> dict[str, object]:
        return {
            "schema": SCHEMA_FILLET_SOLID,
            "kind": self.kind,
            "length": self.length,
            "size_start": self.size_start,
            "size_end": self.size_end,
            "stations": self.stations,
            "arc_points": self.arc_points,
            "vertices": len(self.mesh.vertices),
            "faces": len(self.mesh.faces),
            "base_area": self.base_area,
            "area_coefficient": self.area_coefficient,
            "volume_exact": self.volume_exact,
            "volume_mesh": self.volume_mesh,
            "mesh_volume_error": self.mesh_volume_error,
            "g1_error_deg": self.g1_error_deg,
            "sampled_crease_deg": self.sampled_crease_deg,
            "clamped": self.clamped,
            "watertight": self.watertight,
            "edge_manifold": self.edge_manifold,
            "euler_characteristic": self.euler_characteristic,
            "zero_area_faces": self.zero_area_faces,
            "start_section": self.shape_start.to_dict(),
            "end_section": self.shape_end.to_dict(),
        }


def _point_in_triangle(p: Point2, a: Point2, b: Point2, c: Point2) -> bool:
    d1 = _orientation(a, b, p)
    d2 = _orientation(b, c, p)
    d3 = _orientation(c, a, p)
    has_negative = d1 < 0.0 or d2 < 0.0 or d3 < 0.0
    has_positive = d1 > 0.0 or d2 > 0.0 or d3 > 0.0
    return not (has_negative and has_positive)


def _triangulate(ring: Sequence[Point2]) -> list[tuple[int, int, int]]:
    """Ear-clipping triangulation of a simple counter-clockwise polygon."""
    remaining = list(range(len(ring)))
    triangles: list[tuple[int, int, int]] = []
    guard = 0
    while len(remaining) > 3:
        guard += 1
        if guard > 4 * len(ring) + 16:
            raise FilletError("cap triangulation stalled; the section is not a simple polygon")
        clipped = False
        for position in range(len(remaining)):
            i0 = remaining[position - 1]
            i1 = remaining[position]
            i2 = remaining[(position + 1) % len(remaining)]
            a, b, c = ring[i0], ring[i1], ring[i2]
            if _orientation(a, b, c) <= 0.0:
                continue
            if any(
                _point_in_triangle(ring[j], a, b, c)
                for j in remaining
                if j not in (i0, i1, i2)
            ):
                continue
            triangles.append((i0, i1, i2))
            remaining.pop(position)
            clipped = True
            break
        if not clipped:
            raise FilletError("cap triangulation found no ear; the section is not simple")
    if len(remaining) == 3:
        triangles.append((remaining[0], remaining[1], remaining[2]))
    return triangles


def _sweep(
    points: Sequence[Sequence[float]],
    length: float,
    size_start: float,
    size_end: float | None,
    *,
    kind: str,
    arc_points: int,
    stations: int,
    allow_clamp: bool,
) -> FilletSolid:
    if not _finite(length) or length <= 0.0:
        raise FilletError(f"the sweep length must be positive and finite: {length!r}")
    if stations < 1:
        raise FilletError(f"stations must be at least 1, got {stations}")
    ring = _clean_ring(points)
    turns = _turns_of(ring)
    limit = _max_radius_of(ring, turns) if kind == "fillet" else _min_edge_half(ring)
    end_request = size_start if size_end is None else size_end
    for request in (size_start, end_request):
        if not _finite(request) or request < 0.0:
            raise FilletError(f"the {kind} size must be finite and non-negative: {request!r}")
    start, end = float(size_start), float(end_request)
    clamped = False
    biggest = max(start, end)
    if biggest <= limit * (1.0 + 1e-9):
        start, end = min(start, limit), min(end, limit)
    else:
        if not allow_clamp:
            raise FilletError(
                f"the requested {kind} size {biggest!r} overruns an edge: the largest feasible "
                f"size is {limit!r} (pass allow_clamp=True to clamp to it)"
            )
        scale = limit / biggest if biggest > 0.0 else 1.0
        start *= scale
        end *= scale
        clamped = True
    if max(start, end) == 0.0:
        raise FilletError(f"a {kind} solid needs a positive {kind} size on at least one end")
    base_area = _signed_area(ring)
    coefficient = _area_coefficient(turns, kind)
    volume_exact = (
        base_area * length
        + coefficient * length * (start * start + start * end + end * end) / 3.0
    )

    rings: list[tuple[Point2, ...]] = []
    g1_error = 0.0
    crease = 0.0
    for station in range(stations + 1):
        fraction = station / stations
        size = start + (end - start) * fraction
        _elements, corners, sampled, _perimeter, station_g1 = _corner_path(
            ring, turns, size, arc_points, kind=kind, dense=True
        )
        g1_error = max(g1_error, station_g1)
        rings.append(tuple(sampled))
        # Chord crease at this station: adjacent chords differ in normal by the angle they
        # subtend, |turn| / arc_points (the flat-to-first-chord crease is half of that, so the
        # adjacent-chord pair is the worst one).
        for corner in corners:
            if corner.kind == "collinear" or abs(corner.turn_deg) < 1e-12:
                continue
            if kind == "fillet":
                crease = max(crease, abs(corner.turn_deg) / arc_points)
    builder = MeshBuilder()
    indices = [
        [builder.add_vertex((length * station / stations, u, v)) for (u, v) in ring_points]
        for station, ring_points in enumerate(rings)
    ]
    count = len(rings[0])
    for station in range(stations):
        for index in range(count):
            nxt = (index + 1) % count
            a, b = indices[station][index], indices[station][nxt]
            c, d = indices[station + 1][nxt], indices[station + 1][index]
            builder.add_triangle(a, b, c)
            builder.add_triangle(a, c, d)
    # Caps: the sampled ring may carry coincident consecutive points (the short-edge limit),
    # so triangulate a de-duplicated copy and keep the original indices.
    unique: list[int] = []
    for index, point in enumerate(rings[0]):
        if unique:
            previous = rings[0][unique[-1]]
            if math.hypot(point[0] - previous[0], point[1] - previous[1]) <= 1e-12:
                continue
        unique.append(index)
    cap_points = [rings[0][i] for i in unique]
    if len(cap_points) < 3:
        raise FilletError("the swept section collapses")
    for i0, i1, i2 in _triangulate(cap_points):
        i0, i1, i2 = unique[i0], unique[i1], unique[i2]
        builder.add_triangle(indices[0][i0], indices[0][i2], indices[0][i1])
        builder.add_triangle(indices[-1][i0], indices[-1][i1], indices[-1][i2])
    mesh = remove_degenerate_faces(weld(builder.build(), 1e-12))
    volume_mesh = total_signed_volume(mesh)
    counts = edge_face_counts(mesh)
    edges = len(counts)
    faces = len(mesh.faces)
    vertices = len(mesh.vertices)
    shape_start = _build_shape(
        ring, start, kind=kind, arc_points=arc_points, allow_clamp=False
    )
    shape_end = _build_shape(ring, end, kind=kind, arc_points=arc_points, allow_clamp=False)
    return FilletSolid(
        kind=kind,
        mesh=mesh,
        length=length,
        size_start=start,
        size_end=end,
        stations=stations,
        arc_points=arc_points,
        shape_start=shape_start,
        shape_end=shape_end,
        base_area=base_area,
        area_coefficient=coefficient,
        volume_exact=volume_exact,
        volume_mesh=volume_mesh,
        g1_error_deg=g1_error,
        sampled_crease_deg=crease if kind == "fillet" else g1_error,
        clamped=clamped,
        watertight=not boundary_edges(mesh),
        edge_manifold=all(value == 2 for value in counts.values()),
        euler_characteristic=vertices - edges + faces,
        zero_area_faces=len(zero_area_faces(mesh)),
    )


def fillet_solid(
    points: Sequence[Sequence[float]],
    length: float,
    radius: float,
    *,
    radius_end: float | None = None,
    arc_points: int = 8,
    stations: int = 1,
    allow_clamp: bool = False,
) -> FilletSolid:
    """Sweep a filleted section along ``+x`` into a closed mesh with an exact volume.

    ``radius_end`` (optional) varies the radius linearly along the sweep; the fillet surface
    stays tangent to the adjacent flat faces at every station, which is what ``g1_error_deg``
    measures.  Sizes above the edge limit are refused, or clamped with ``allow_clamp=True``.
    """
    return _sweep(
        points,
        length,
        radius,
        radius_end,
        kind="fillet",
        arc_points=arc_points,
        stations=stations,
        allow_clamp=allow_clamp,
    )


def chamfer_solid(
    points: Sequence[Sequence[float]],
    length: float,
    setback: float,
    *,
    setback_end: float | None = None,
    stations: int = 1,
    allow_clamp: bool = False,
) -> FilletSolid:
    """Sweep a chamfered section along ``+x`` into a closed mesh with an exact volume."""
    return _sweep(
        points,
        length,
        setback,
        setback_end,
        kind="chamfer",
        arc_points=1,
        stations=stations,
        allow_clamp=allow_clamp,
    )


__all__ = [
    "CornerReport",
    "FilletError",
    "FilletSolid",
    "SCHEMA_FILLET",
    "SCHEMA_FILLET_SOLID",
    "SectionMoments",
    "SectionShape",
    "chamfer_section",
    "chamfer_solid",
    "fillet_section",
    "fillet_solid",
    "max_chamfer_setback",
    "max_fillet_radius",
    "polygon_moments",
    "rectangle_section",
    "section_turns",
]