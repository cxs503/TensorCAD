"""Classify the open boundary of an imported CAD surface and cap its planar cuts.

Why this module exists
----------------------
A publisher IGES hull is a set of trimmed patches, not a sewn solid: after
tessellation and welding the mesh still has a boundary, and that boundary is a
mixture of two very different things.

* a **planar cut** -- the file was cut by the CAD system at a plane (the symmetry
  plane of a half hull, the transom, a bow section).  The opening is *information*,
  not damage: capping it invents no surface, the cap lies exactly on the plane the
  file already implies, and the enclosed volume becomes well defined.
* a **defect gap** -- neighbouring patches overlap or do not abut, so the boundary
  wanders through the hull surface.  Nothing can be recovered here; fanning a
  triangle over it would invent metres of surface that the file never had.

A gate verdict that only counts free edges cannot tell these apart, and the two
call for opposite actions.  So this module measures every boundary loop, decides
which kind it is from a least-squares plane fit (a planar cut's vertices lie on
the plane to round-off; a defect gap's do not), and -- only for the first kind --
caps it.  Every number that went into a decision travels in the returned record,
so a caller can re-derive it instead of trusting a label.

The module is standard-library only and deterministic: loop order, ear order and
plane normals are all fixed by index order, not by a set iteration order.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any

__all__ = [
    "OpeningLoop",
    "OpeningReport",
    "boundary_loops",
    "classify_openings",
    "seal_planar_openings",
]

Vec3 = tuple[float, float, float]
Triangle = tuple[int, int, int]

#: A loop counts as a planar cut when every vertex lies within this fraction of the
#: loop's own diagonal from the least-squares plane.  1e-3 is two orders of magnitude
#: looser than round-off on a 250 m hull (mm on a metre-scale chord) and two orders
#: tighter than the metre-scale wander of a mis-registered patch seam.
DEFAULT_PLANAR_RATIO = 1e-3


# --------------------------------------------------------------------------- #
# small vector helpers (kept local: geomquality stays dependency-free)
# --------------------------------------------------------------------------- #
def _sub(a: Sequence[float], b: Sequence[float]) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _dot(a: Sequence[float], b: Sequence[float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a: Sequence[float], b: Sequence[float]) -> Vec3:
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def _norm(a: Sequence[float]) -> float:
    return math.sqrt(_dot(a, a))


def _normalise(a: Sequence[float]) -> Vec3:
    length = _norm(a)
    if length <= 0.0:
        return (0.0, 0.0, 0.0)
    return (a[0] / length, a[1] / length, a[2] / length)


def _centroid(points: Sequence[Vec3]) -> Vec3:
    count = float(len(points))
    return (sum(p[0] for p in points) / count,
            sum(p[1] for p in points) / count,
            sum(p[2] for p in points) / count)


def _extent(points: Sequence[Vec3]) -> float:
    """Longest edge of the axis-aligned box around ``points``."""
    lo = [min(p[axis] for p in points) for axis in range(3)]
    hi = [max(p[axis] for p in points) for axis in range(3)]
    return max(hi[axis] - lo[axis] for axis in range(3))


def _plane_normal(points: Sequence[Vec3]) -> Vec3:
    """Unit normal of the least-squares plane through ``points``.

    The plane's normal is the eigenvector of the point covariance with the smallest
    eigenvalue.  Inverse iteration converges to it: the power method applied to
    ``trace * I - C``, whose dominant eigenvector is ``C``'s smallest one.  Thirty
    steps are ample for a 3x3 symmetric matrix and cost nothing next to the loop walk.
    """
    if len(points) < 3:
        return (0.0, 0.0, 0.0)
    centre = _centroid(points)
    cov = [[0.0] * 3 for _ in range(3)]
    for point in points:
        delta = _sub(point, centre)
        for row in range(3):
            delta_row = delta[row]
            for column in range(3):
                cov[row][column] += delta_row * delta[column]
    trace = cov[0][0] + cov[1][1] + cov[2][2]

    def iterate(seed: Vec3) -> Vec3:
        vector = seed
        for _ in range(200):
            multiplied = (trace * vector[0] - (cov[0][0] * vector[0] + cov[0][1] * vector[1] + cov[0][2] * vector[2]),
                          trace * vector[1] - (cov[1][0] * vector[0] + cov[1][1] * vector[1] + cov[1][2] * vector[2]),
                          trace * vector[2] - (cov[2][0] * vector[0] + cov[2][1] * vector[1] + cov[2][2] * vector[2]))
            length = _norm(multiplied)
            if length <= 1e-300:
                return (0.0, 0.0, 0.0)
            vector = (multiplied[0] / length, multiplied[1] / length, multiplied[2] / length)
        return vector

    # Power iteration can stall on a start vector that is nearly orthogonal to the
    # dominant direction, so try the coordinate axes too and keep the best fit.
    best: Vec3 = (0.0, 0.0, 0.0)
    best_residual = float("inf")
    for seed in ((1.0, 1.0, 1.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)):
        candidate = iterate(seed)
        if candidate == (0.0, 0.0, 0.0):
            continue
        residual = max(abs(_dot(_sub(point, centre), candidate)) for point in points)
        if residual < best_residual:
            best, best_residual = candidate, residual
    return best


# --------------------------------------------------------------------------- #
# boundary loops
# --------------------------------------------------------------------------- #
def boundary_loops(vertices: Sequence[Sequence[float]],
                   faces: Sequence[Sequence[int]]) -> list[list[int]]:
    """Walk the mesh's free edges into closed vertex loops, deterministically.

    A free edge is one carried by exactly one face.  For a consistently oriented mesh
    those edges are already directed (the face's own winding), and a boundary loop is
    the cycle they form.  Loops shorter than three vertices, and chains that fail to
    close (which happens when an edge in the chain is non-manifold and shares three
    faces), are dropped: a cap needs a cycle.
    """
    counts: dict[tuple[int, int], int] = {}
    for face in faces:
        a, b, c = int(face[0]), int(face[1]), int(face[2])
        for start, end in ((a, b), (b, c), (c, a)):
            key = (start, end) if start < end else (end, start)
            counts[key] = counts.get(key, 0) + 1

    directed: list[tuple[int, int]] = []
    for face in faces:
        a, b, c = int(face[0]), int(face[1]), int(face[2])
        for start, end in ((a, b), (b, c), (c, a)):
            key = (start, end) if start < end else (end, start)
            if counts[key] == 1:
                directed.append((start, end))

    outgoing: dict[int, list[int]] = {}
    for start, end in directed:
        outgoing.setdefault(start, []).append(end)

    loops: list[list[int]] = []
    used: set[tuple[int, int]] = set()
    for start, end in directed:
        if (start, end) in used:
            continue
        loop = [start]
        current = end
        used.add((start, end))
        guard = len(directed) + 2
        while current != start and guard > 0:
            guard -= 1
            loop.append(current)
            candidates = [nxt for nxt in outgoing.get(current, ())
                          if (current, nxt) not in used]
            if not candidates:
                break
            nxt = candidates[0]
            used.add((current, nxt))
            current = nxt
        if current == start and len(loop) >= 3:
            loops.append(loop)
    return loops


# --------------------------------------------------------------------------- #
# classification
# --------------------------------------------------------------------------- #
class OpeningLoop:
    """One boundary loop, measured."""

    __slots__ = ("index", "vertices", "edges", "perimeter", "longest_edge",
                 "planarity_residual", "planarity_ratio", "normal", "kind", "capped",
                 "cap_area")

    def __init__(self, index: int, vertices: list[int]) -> None:
        self.index = index
        self.vertices = vertices
        self.edges = len(vertices)
        self.perimeter = 0.0
        self.longest_edge = 0.0
        self.planarity_residual = 0.0
        self.planarity_ratio = 0.0
        self.normal: Vec3 = (0.0, 0.0, 0.0)
        self.kind = "unknown"
        self.capped = False
        self.cap_area = 0.0

    def to_dict(self) -> dict[str, Any]:
        return {
            "index": self.index,
            "edges": self.edges,
            "perimeter": self.perimeter,
            "longest_edge": self.longest_edge,
            "planarity_residual": self.planarity_residual,
            "planarity_ratio": self.planarity_ratio,
            "normal": list(self.normal),
            "kind": self.kind,
            "capped": self.capped,
            "cap_area": self.cap_area,
        }


class OpeningReport:
    """The measured opening of one mesh, plus the verdict the measurements support."""

    __slots__ = ("loops", "free_edges", "planar_tolerance", "verdict")

    def __init__(self) -> None:
        self.loops: list[OpeningLoop] = []
        self.free_edges = 0
        self.planar_tolerance = 0.0
        self.verdict = "unknown"

    @property
    def planar_edges(self) -> int:
        return sum(loop.edges for loop in self.loops if loop.kind == "planar_cut")

    @property
    def defect_edges(self) -> int:
        return sum(loop.edges for loop in self.loops if loop.kind == "defect_gap")

    def settle(self) -> None:
        """Fix the verdict from the loop kinds; call once after measuring."""
        if self.free_edges == 0:
            self.verdict = "closed"
        elif self.defect_edges == 0:
            self.verdict = "closable_by_planar_caps"
        elif self.planar_edges == 0:
            self.verdict = "defect_gaps_present"
        else:
            self.verdict = "mixed_planar_cuts_and_defect_gaps"

    def to_dict(self) -> dict[str, Any]:
        return {
            "free_edges": self.free_edges,
            "planar_tolerance": self.planar_tolerance,
            "verdict": self.verdict,
            "planar_edges": self.planar_edges,
            "defect_edges": self.defect_edges,
            "planar_edge_fraction": (self.planar_edges / self.free_edges
                                    if self.free_edges else 1.0),
            "loops": [loop.to_dict() for loop in self.loops],
        }


def _measure_loop(loop: list[int], point_of: Sequence[Vec3], index: int,
                  planar_tolerance: float) -> OpeningLoop:
    measured = OpeningLoop(index, loop)
    points = [point_of[vertex] for vertex in loop]
    for position, start in enumerate(loop):
        end = loop[(position + 1) % len(loop)]
        length = _norm(_sub(point_of[end], point_of[start]))
        measured.perimeter += length
        measured.longest_edge = max(measured.longest_edge, length)
    normal = _plane_normal(points)
    measured.normal = normal
    centre = _centroid(points)
    residual = max(abs(_dot(_sub(point, centre), normal)) for point in points)
    measured.planarity_residual = residual
    span = _extent(points)
    measured.planarity_ratio = residual / span if span > 0.0 else 0.0
    measured.kind = "planar_cut" if residual <= planar_tolerance else "defect_gap"
    return measured


def classify_openings(vertices: Sequence[Sequence[float]],
                      faces: Sequence[Sequence[int]],
                      *,
                      planar_tolerance: float | None = None) -> OpeningReport:
    """Measure every boundary loop and label it ``planar_cut`` or ``defect_gap``.

    ``planar_tolerance`` is an absolute distance in model units.  When it is ``None``
    the tolerance is :data:`DEFAULT_PLANAR_RATIO` of the mesh's own largest extent, so
    a metre-scale hull and a millimetre-scale model are judged by the same relative
    standard.
    """
    point_of: list[Vec3] = [(float(vertex[0]), float(vertex[1]), float(vertex[2]))
                            for vertex in vertices]
    counts: dict[tuple[int, int], int] = {}
    for face in faces:
        a, b, c = int(face[0]), int(face[1]), int(face[2])
        for start, end in ((a, b), (b, c), (c, a)):
            key = (start, end) if start < end else (end, start)
            counts[key] = counts.get(key, 0) + 1
    free_edges = sum(1 for count in counts.values() if count == 1)

    report = OpeningReport()
    report.free_edges = free_edges
    if planar_tolerance is None:
        span = _extent(point_of) if point_of else 0.0
        report.planar_tolerance = DEFAULT_PLANAR_RATIO * span
    else:
        report.planar_tolerance = float(planar_tolerance)
    for index, loop in enumerate(boundary_loops(point_of, faces)):
        report.loops.append(_measure_loop(loop, point_of, index, report.planar_tolerance))
    # Report the biggest openings first: the sorted order is also the useful one.
    report.loops.sort(key=lambda loop: (-loop.edges, loop.index))
    for position, loop in enumerate(report.loops):
        loop.index = position
    report.settle()
    return report


# --------------------------------------------------------------------------- #
# capping
# --------------------------------------------------------------------------- #
def _triangulate_polygon(points2d: Sequence[tuple[float, float]]) -> list[tuple[int, int, int]] | None:
    """Ear-clip a simple polygon given in its own plane; ``None`` when it stalls.

    Convexity is tested with the polygon's own signed area sign, so the caller does not
    have to know the winding.  A stall (numerically collinear neighbourhood, a
    self-intersecting loop) returns ``None`` rather than a wrong triangulation.
    """
    count = len(points2d)
    if count < 3:
        return None
    area = 0.0
    for index in range(count):
        x0, y0 = points2d[index]
        x1, y1 = points2d[(index + 1) % count]
        area += x0 * y1 - x1 * y0
    if area == 0.0:
        return None
    sign = 1.0 if area > 0.0 else -1.0
    remaining = list(range(count))
    triangles: list[tuple[int, int, int]] = []
    guard = 2 * count * count
    while len(remaining) > 3 and guard > 0:
        guard -= 1
        clipped = False
        for position in range(len(remaining)):
            previous = remaining[position - 1]
            current = remaining[position]
            following = remaining[(position + 1) % len(remaining)]
            ax, ay = points2d[previous]
            bx, by = points2d[current]
            cx, cy = points2d[following]
            turn = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
            if sign * turn <= 0.0:
                continue
            contains = False
            for other in remaining:
                if other in (previous, current, following):
                    continue
                px, py = points2d[other]
                d1 = (bx - ax) * (py - ay) - (by - ay) * (px - ax)
                d2 = (cx - bx) * (py - by) - (cy - by) * (px - bx)
                d3 = (ax - cx) * (py - cy) - (ay - cy) * (px - cx)
                if sign * d1 >= 0.0 and sign * d2 >= 0.0 and sign * d3 >= 0.0:
                    contains = True
                    break
            if contains:
                continue
            triangles.append((previous, current, following))
            remaining.pop(position)
            clipped = True
            break
        if not clipped:
            return None
    if len(remaining) != 3:
        return None
    triangles.append((remaining[0], remaining[1], remaining[2]))
    return triangles


def seal_planar_openings(vertices: Sequence[Sequence[float]],
                         faces: Sequence[Sequence[int]],
                         *,
                         planar_tolerance: float | None = None) -> tuple[list[Vec3], list[Triangle], dict[str, Any]]:
    """Cap every ``planar_cut`` opening; leave ``defect_gap`` openings open and reported.

    A cap shares the loop's existing vertices -- no new surface is invented off the
    plane, and the plane is the one the file's own vertices define.  Triangles are
    appended with the loop's reverse winding, because a consistently oriented surface
    carries its boundary edge in the direction the cap must *not* use.
    """
    point_of: list[Vec3] = [(float(vertex[0]), float(vertex[1]), float(vertex[2]))
                            for vertex in vertices]
    working = [tuple(int(index) for index in face) for face in faces]
    report = classify_openings(point_of, working, planar_tolerance=planar_tolerance)
    added = 0
    capped_loops = 0
    for loop in report.loops:
        if loop.kind != "planar_cut":
            continue
        points = [point_of[vertex] for vertex in loop.vertices]
        normal = loop.normal
        centre = _centroid(points)
        # An in-plane basis: any vector not parallel to the normal, then a cross product.
        seed = (1.0, 0.0, 0.0) if abs(normal[0]) < 0.9 else (0.0, 1.0, 0.0)
        axis_u = _normalise(_cross(normal, seed))
        axis_v = _cross(normal, axis_u)
        points2d = [(_dot(_sub(point, centre), axis_u), _dot(_sub(point, centre), axis_v))
                    for point in points]
        triangles = _triangulate_polygon(points2d)
        if triangles is None:
            loop.capped = False
            continue
        area = 0.0
        for i, j, k in triangles:
            a, b, c = points[i], points[j], points[k]
            area += 0.5 * _norm(_cross(_sub(b, a), _sub(c, a)))
            loop_arities = (loop.vertices[i], loop.vertices[j], loop.vertices[k])
            # Reverse winding: the surface's boundary edge runs the other way round.
            working.append((loop_arities[0], loop_arities[2], loop_arities[1]))
            added += 1
        loop.capped = True
        loop.cap_area = area
        capped_loops += 1
    report.settle()
    record = {
        "applied": capped_loops > 0,
        "loops_capped": capped_loops,
        "triangles_added": added,
        "capped_area": sum(loop.cap_area for loop in report.loops if loop.capped),
        "classification": report.to_dict(),
    }
    return list(point_of), working, record