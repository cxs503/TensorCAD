"""Triangle/triangle intersection and mesh self-intersection detection.

A closed manifold mesh can still be geometrically useless: if the surface folds
back through itself the enclosed volume is not a solid, and a downstream volume
mesher either refuses it or silently produces an inverted cell layer. Watertight
topology (:func:`tensorcad.mesh.is_watertight`) says nothing about that, so the
PC6 geometry service runs an explicit geometric test.

The primitive is the standard edge/face formulation:

* two non-coplanar triangles intersect **iff** at least one closed edge of one
  triangle pierces the interior of the other (the intersection segment of two
  convex triangles always begins or ends on a boundary edge of one of them);
* coplanar pairs are decided by point-in-triangle containment in a projected
  plane, because both triangles then share the same plane and no edge lifts out
  of it.

Everything is pure standard library. The mesh-level search is pruned with a
uniform hash grid over triangle bounding boxes, so the cost is close to linear
in the face count instead of the quadratic all-pairs scan.
"""
from __future__ import annotations

import math
from collections.abc import Sequence

from .geometry import Face, Mesh, Vec3

#: Relative tolerance on plane distances and barycentric coordinates. Every
#: quantity is scaled by the triangle's own edge length before comparison, so the
#: same tolerance is meaningful for a 0.01 m part and a 200 m hull.
RELATIVE_TOLERANCE = 1e-10


def _sub(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _cross(a: Vec3, b: Vec3) -> Vec3:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _dot(a: Vec3, b: Vec3) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _norm(a: Vec3) -> float:
    return math.sqrt(_dot(a, a))


def _scale(triangle: Sequence[Vec3]) -> float:
    """Characteristic size of a triangle: its longest edge."""
    return max(
        _norm(_sub(triangle[1], triangle[0])),
        _norm(_sub(triangle[2], triangle[1])),
        _norm(_sub(triangle[0], triangle[2])),
    )


def _plane(triangle: Sequence[Vec3]) -> tuple[Vec3, float, float]:
    """Unit normal, plane offset and characteristic size of a triangle.

    The normal is canonicalised (the largest-magnitude component is made positive)
    so that two triangles on the same plane always report the same offset, even
    when their winding makes the raw normals antiparallel. Without that, a
    back-to-back pair 250 m apart looked coplanar.
    """
    normal = _cross(_sub(triangle[1], triangle[0]), _sub(triangle[2], triangle[0]))
    length = _norm(normal)
    if length == 0.0:
        return (0.0, 0.0, 0.0), 0.0, _scale(triangle)
    unit = (normal[0] / length, normal[1] / length, normal[2] / length)
    dominant = max(range(3), key=lambda index: abs(unit[index]))
    if unit[dominant] < 0.0:
        unit = (-unit[0], -unit[1], -unit[2])
    return unit, _dot(unit, triangle[0]), _scale(triangle)


def _point_in_triangle(point: Vec3, triangle: Sequence[Vec3], tolerance: float) -> bool:
    """Containment test for a point already known to lie in the triangle's plane."""
    normal, _, _ = _plane(triangle)
    if normal == (0.0, 0.0, 0.0):
        return False
    # Project away the dominant normal axis: the 2-D test is then exact up to the
    # rounding of the two remaining coordinates.
    axis = max(range(3), key=lambda index: abs(normal[index]))
    u_axis, v_axis = [(1, 2), (2, 0), (0, 1)][axis]
    a = (triangle[0][u_axis], triangle[0][v_axis])
    b = (triangle[1][u_axis], triangle[1][v_axis])
    c = (triangle[2][u_axis], triangle[2][v_axis])
    p = (point[u_axis], point[v_axis])
    denominator = (b[1] - c[1]) * (a[0] - c[0]) + (c[0] - b[0]) * (a[1] - c[1])
    if denominator == 0.0:
        return False
    w0 = ((b[1] - c[1]) * (p[0] - c[0]) + (c[0] - b[0]) * (p[1] - c[1])) / denominator
    w1 = ((c[1] - a[1]) * (p[0] - c[0]) + (a[0] - c[0]) * (p[1] - c[1])) / denominator
    w2 = 1.0 - w0 - w1
    return w0 >= -tolerance and w1 >= -tolerance and w2 >= -tolerance


def _segment_pierces(
    start: Vec3, end: Vec3, triangle: Sequence[Vec3], tolerance: float
) -> bool:
    """True when the closed segment ``start``-``end`` crosses the triangle interior.

    Coplanar segments are reported as *not* piercing; the coplanar case is decided
    separately by :func:`triangles_intersect`.
    """
    normal, offset, size = _plane(triangle)
    if normal == (0.0, 0.0, 0.0):
        return False
    limit = tolerance * max(size, 1.0)
    d0 = _dot(normal, start) - offset
    d1 = _dot(normal, end) - offset
    if abs(d0) <= limit and abs(d1) <= limit:
        return False  # coplanar with the triangle's plane
    if d0 * d1 > 0.0:
        return False  # both ends on the same side
    if d0 == d1:
        return False  # parallel to the plane, never crosses
    t = d0 / (d0 - d1)
    point = (
        start[0] + t * (end[0] - start[0]),
        start[1] + t * (end[1] - start[1]),
        start[2] + t * (end[2] - start[2]),
    )
    return _point_in_triangle(point, triangle, 1e-9)


def triangles_intersect(
    first: Sequence[Vec3], second: Sequence[Vec3], *, tolerance: float = RELATIVE_TOLERANCE
) -> bool:
    """True when two closed triangles share at least one point.

    Touching at a shared vertex or edge counts as intersecting; mesh callers skip
    adjacency explicitly (see :func:`find_self_intersections`), because two faces
    that share a vertex are *supposed* to meet there.
    """
    if len(first) != 3 or len(second) != 3:  # pragma: no cover - caller contract
        raise ValueError("triangles_intersect takes two triangles")
    a_min = tuple(min(point[axis] for point in first) for axis in range(3))
    a_max = tuple(max(point[axis] for point in first) for axis in range(3))
    b_min = tuple(min(point[axis] for point in second) for axis in range(3))
    b_max = tuple(max(point[axis] for point in second) for axis in range(3))
    pad = tolerance * max(_scale(first), _scale(second), 1.0)
    for axis in range(3):
        if a_max[axis] + pad < b_min[axis] or b_max[axis] + pad < a_min[axis]:
            return False

    n1, d1, _ = _plane(first)
    n2, d2, _ = _plane(second)
    if n1 == (0.0, 0.0, 0.0) or n2 == (0.0, 0.0, 0.0):
        return False  # degenerate triangle carries no area, nothing to intersect
    parallel = _norm(_cross(n1, n2)) <= tolerance
    size = max(_scale(first), _scale(second), 1.0)
    if parallel and abs(d1 - d2) <= tolerance * size:
        # Coplanar: containment of any vertex in the other triangle is the test.
        return (
            any(_point_in_triangle(point, second, 1e-9) for point in first)
            or any(_point_in_triangle(point, first, 1e-9) for point in second)
        )
    for i in range(3):
        if _segment_pierces(first[i], first[(i + 1) % 3], second, tolerance):
            return True
        if _segment_pierces(second[i], second[(i + 1) % 3], first, tolerance):
            return True
    return False


def _grid_cell(
    bbox_min: Vec3, bbox_max: Vec3, origin: Vec3, size: float
) -> tuple[int, int, int, int, int, int]:
    low = tuple(max(0, math.floor((bbox_min[axis] - origin[axis]) / size)) for axis in range(3))
    high = tuple(max(0, math.floor((bbox_max[axis] - origin[axis]) / size)) for axis in range(3))
    return (low[0], low[1], low[2], high[0], high[1], high[2])


def find_self_intersections(
    mesh: Mesh, *, limit: int = 64, cells_per_axis: int = 24
) -> tuple[tuple[int, int], ...]:
    """Pairs of non-adjacent faces that intersect each other.

    Faces sharing a vertex are skipped: a closed manifold mesh has them meeting at
    that vertex by construction, and reporting those would drown the real folds.
    The returned list is sorted and truncated to ``limit`` pairs -- a non-empty
    result is the finding, the exact count beyond the limit is not.
    """
    if limit < 1:
        raise ValueError("limit must be at least 1")
    if len(mesh.faces) < 2:
        return ()
    points = mesh.vertices
    boxes: list[tuple[Vec3, Vec3]] = []
    origin: Vec3 = (
        min(point[0] for point in points),
        min(point[1] for point in points),
        min(point[2] for point in points),
    )
    diagonal = max(
        max(point[axis] for point in points) - origin[axis] for axis in range(3)
    )
    if diagonal <= 0.0:
        return ()
    size = diagonal / max(cells_per_axis, 1)
    grid: dict[tuple[int, int, int], list[int]] = {}
    for index, face in enumerate(mesh.faces):
        corners = [points[vertex] for vertex in face]
        low: Vec3 = (
            min(point[0] for point in corners),
            min(point[1] for point in corners),
            min(point[2] for point in corners),
        )
        high: Vec3 = (
            max(point[0] for point in corners),
            max(point[1] for point in corners),
            max(point[2] for point in corners),
        )
        boxes.append((low, high))
        x0, y0, z0, x1, y1, z1 = _grid_cell(low, high, origin, size)
        for ix in range(x0, x1 + 1):
            for iy in range(y0, y1 + 1):
                for iz in range(z0, z1 + 1):
                    grid.setdefault((ix, iy, iz), []).append(index)

    found: set[tuple[int, int]] = set()
    for bucket in grid.values():
        count = len(bucket)
        for i in range(count):
            first = bucket[i]
            first_face = mesh.faces[first]
            first_corners = [points[vertex] for vertex in first_face]
            for j in range(i + 1, count):
                second = bucket[j]
                pair = (first, second) if first < second else (second, first)
                if pair in found:
                    continue
                second_face = mesh.faces[second]
                if set(first_face) & set(second_face):
                    continue
                low, high = boxes[second]
                overlaps = all(
                    low[axis] <= boxes[first][1][axis] and boxes[first][0][axis] <= high[axis]
                    for axis in range(3)
                )
                if not overlaps:
                    continue
                second_corners = [points[vertex] for vertex in second_face]
                if triangles_intersect(first_corners, second_corners):
                    found.add(pair)
                    if len(found) >= limit:
                        return tuple(sorted(found))
    return tuple(sorted(found))


__all__ = ["RELATIVE_TOLERANCE", "Face", "find_self_intersections", "triangles_intersect"]
