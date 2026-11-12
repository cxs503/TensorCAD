"""Independent oracles for the mid-surface tests: distances and section offsets.

Nothing here imports :mod:`tensorcad.midsurface`. The two questions the module answers are
answered a second time, by different arithmetic:

* *is the sheet where it claims to be* -- a point-to-triangle distance against the source
  surface (brute force over every triangle, no grid, no rays), and
* *how big should the sheet be* -- a 2-D half-plane clipping of the section, extruded, which
  reproduces the "a corner is the intersection of the two mid-planes" convention the module
  uses without sharing a line of code with it.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

from tensorcad.geometry import Mesh, Vec3


def point_triangle_distance(
    point: Sequence[float], a: Sequence[float], b: Sequence[float], c: Sequence[float]
) -> float:
    """Distance from ``point`` to the closed triangle ``abc`` (Ericson, *Real-Time Collision*)."""

    def sub(u, v):
        return (u[0] - v[0], u[1] - v[1], u[2] - v[2])

    def dot(u, v):
        return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]

    ab = sub(b, a)
    ac = sub(c, a)
    ap = sub(point, a)
    d1 = dot(ab, ap)
    d2 = dot(ac, ap)
    if d1 <= 0.0 and d2 <= 0.0:
        closest = a
    else:
        bp = sub(point, b)
        d3 = dot(ab, bp)
        d4 = dot(ac, bp)
        cp = sub(point, c)
        d5 = dot(ab, cp)
        d6 = dot(ac, cp)
        vc = d1 * d4 - d3 * d2
        vb = d5 * d2 - d1 * d6
        va = d3 * d6 - d5 * d4
        if d3 >= 0.0 and d4 <= d3:
            closest = b
        elif d6 >= 0.0 and d5 <= d6:
            closest = c
        elif vc <= 0.0 and d1 >= 0.0 and d3 <= 0.0:
            denom = d1 - d3
            closest = tuple(a[axis] + ab[axis] * (d1 / denom) for axis in range(3))
        elif vb <= 0.0 and d2 >= 0.0 and d6 <= 0.0:
            denom = d2 - d6
            closest = tuple(a[axis] + ac[axis] * (d2 / denom) for axis in range(3))
        elif va <= 0.0 and (d4 - d3) >= 0.0 and (d5 - d6) >= 0.0:
            denom = (d4 - d3) + (d5 - d6)
            edge = sub(c, b)
            closest = tuple(b[axis] + edge[axis] * ((d4 - d3) / denom) for axis in range(3))
        else:
            total = va + vb + vc
            v = vb / total
            w = vc / total
            closest = tuple(a[axis] + ab[axis] * v + ac[axis] * w for axis in range(3))
    diff = sub(point, closest)
    return math.sqrt(dot(diff, diff))


def distance_to_surface(point: Sequence[float], mesh: Mesh) -> float:
    """Smallest distance from ``point`` to any triangle of ``mesh`` -- brute force."""
    return min(
        point_triangle_distance(
            point, mesh.vertices[face[0]], mesh.vertices[face[1]], mesh.vertices[face[2]]
        )
        for face in mesh.faces
    )


def clip_polygon_2d(
    polygon: Sequence[tuple[float, float]], a: float, b: float, c: float
) -> list[tuple[float, float]]:
    """Sutherland-Hodgman against ``a x + b y + c >= 0`` (convex subject, convex clip)."""
    result: list[tuple[float, float]] = []
    if not polygon:
        return result
    previous = polygon[-1]
    previous_value = a * previous[0] + b * previous[1] + c
    for current in polygon:
        current_value = a * current[0] + b * current[1] + c
        if current_value >= 0.0:
            if previous_value < 0.0:
                result.append(_cross_point(previous, current, previous_value, current_value))
            result.append(current)
        elif previous_value >= 0.0:
            result.append(_cross_point(previous, current, previous_value, current_value))
        previous, previous_value = current, current_value
    return result


def _cross_point(a, b, va, vb):
    denominator = va - vb
    s = va / denominator if abs(denominator) > 1e-300 else 0.0
    return (a[0] + (b[0] - a[0]) * s, a[1] + (b[1] - a[1]) * s)


def offset_section(
    section: Sequence[tuple[float, float]], distance: float
) -> list[tuple[float, float]]:
    """Erode a CCW convex section by ``distance`` with one half-plane per edge.

    For a CCW polygon the interior lies to the left of every edge, so the inward constraint of
    edge ``p -> q`` is ``n . (x - p) >= distance`` with ``n`` the left normal. Two offset
    half-planes therefore meet at the offset corner -- the same convention the module's
    mid-plane clipping realises in 3-D, reached here by 2-D arithmetic alone.
    """
    polygon = list(section)
    for index, point in enumerate(section):
        nxt = section[(index + 1) % len(section)]
        dx = nxt[0] - point[0]
        dy = nxt[1] - point[1]
        length = math.hypot(dx, dy)
        # Left normal of the direction (dx, dy) is (-dy, dx); it points into a CCW polygon.
        a = -dy / length
        b = dx / length
        c = -(a * point[0] + b * point[1]) - distance
        polygon = clip_polygon_2d(polygon, a, b, c)
    return polygon


def polygon_area(polygon: Sequence[tuple[float, float]]) -> float:
    total = 0.0
    for index, point in enumerate(polygon):
        nxt = polygon[(index + 1) % len(polygon)]
        total += point[0] * nxt[1] - nxt[0] * point[1]
    return 0.5 * total


def polygon_perimeter(polygon: Sequence[tuple[float, float]]) -> float:
    total = 0.0
    for index, point in enumerate(polygon):
        nxt = polygon[(index + 1) % len(polygon)]
        total += math.hypot(nxt[0] - point[0], nxt[1] - point[1])
    return total


def prism_mid_sheet(section, length: float, thickness: float) -> tuple[float, float]:
    """``(area, volume)`` of the mid-sheet of a prismatic thin-walled solid.

    The mid-sheet of a prism of convex section ``S``, length ``H`` and wall ``t`` is the eroded
    section extruded over ``H - t``: a closed sheet, so it has both an area and a volume.
    """
    eroded = offset_section(section, 0.5 * thickness)
    area_section = polygon_area(eroded)
    perimeter = polygon_perimeter(eroded)
    inner_length = length - thickness
    return (perimeter * inner_length + 2.0 * area_section, area_section * inner_length)


def regular_section(inradius: float, sides: int) -> tuple[tuple[float, float], ...]:
    """Independent regular polygon on the origin, inradius to facet (CCW)."""
    angle = math.pi / sides
    circumradius = inradius / math.cos(angle)
    return tuple(
        (
            circumradius * math.cos(2.0 * angle * index + angle),
            circumradius * math.sin(2.0 * angle * index + angle),
        )
        for index in range(sides)
    )


def regular_polygon_area(inradius: float, sides: int) -> float:
    """Closed form ``n r^2 tan(pi/n)`` -- a second opinion on :func:`polygon_area`."""
    return sides * inradius * inradius * math.tan(math.pi / sides)


def regular_polygon_perimeter(inradius: float, sides: int) -> float:
    return 2.0 * sides * inradius * math.tan(math.pi / sides)