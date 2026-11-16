"""Open-boundary accounting: decide *what* each free edge of a tessellated CAD set is.

The import pipeline can say how many free edges a mesh has and how far each one is from
its nearest neighbour chord (:func:`tensorcad.geomquality.stitch.seam_report`), but a
free-edge *count* answers none of the questions a repair decision needs:

* is this opening a **seam** -- the neighbour's boundary passes within the coincidence
  budget, so sewing can close it (:mod:`.stitch`);
* is it a **patch overlap** -- the boundary lies on the *interior* of another face,
  which is the signature of untrimmed patches exported on top of each other; no amount
  of stitching can close that, the surface has to be trimmed (:mod:`.trim_mesh`);
* is it a **declared opening** -- on a symmetry plane, the deck line, a transom: geometry
  that is *supposed* to be open;
* or is it **unexplained** -- nothing within the budget, i.e. either a real gap of
  metres or a genuine hole.

The distinction is measurable and it is what turns "the file is not watertight" into a
closure budget in metres per cause.  The classes are exclusive and ordered: a chord on a
declared plane is reported as plane even if a neighbour chord is also nearby, because the
plane explains it.

Design rules:

* **No repair, no tolerance widening.**  This module only measures; it never moves a
  vertex.  The caller decides.
* **Length weighted.**  A budget in metres per class is what a modeller can act on; a
  count of chords is not, because chord lengths differ by three orders of magnitude on a
  real hull (a deck chord is metres, a stern knuckle chord is millimetres).
* **Pure Python, no runtime dependency.**  The core package keeps zero runtime
  dependencies, so the spatial hashes here are hand-rolled like the ones in
  :mod:`.stitch`.
"""
from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Sequence
from typing import Any

from .stitch import (
    ENDPOINT_FRACTION,
    Vec3,
    Triangle,
    _EdgeGrid,
    free_edges,
    project_on_segment,
)

#: The boundary lies on a neighbour's boundary chord: sewing can close it.
CLASS_SEAM = "seam"
#: The boundary lies on another face's interior: the patches overlap and need trimming.
CLASS_OVERLAP = "overlap"
#: Nothing within the budget: a genuine gap or a real hole.
CLASS_UNEXPLAINED = "unexplained"
#: Prefix used for a boundary that lies on a declared plane (``plane:y=0``).
CLASS_PLANE_PREFIX = "plane:"

CLASSES = (CLASS_SEAM, CLASS_OVERLAP, CLASS_UNEXPLAINED)


class _TriangleGrid:
    """Uniform spatial hash of triangles, for point-to-surface queries.

    Vertical faces of a ship hull are large and the chords that matter are short, so the
    cell is sized from the chord statistics (``cell`` argument) rather than from the
    triangles: one bucket may then hold many triangles, which is the cheap direction.
    """

    def __init__(self, vertices: Sequence[Vec3], faces: Sequence[Triangle],
                 cell: float) -> None:
        self.vertices = vertices
        self.cell = cell
        self.faces = list(faces)
        self.buckets: dict[tuple[int, int, int], list[int]] = defaultdict(list)
        self._large: list[int] = []
        for index, (a, b, c) in enumerate(self.faces):
            corners = (vertices[a], vertices[b], vertices[c])
            keys = {self._key(point) for point in corners}
            if len(keys) > 64:
                # A triangle spanning very many cells (a sliver along the hull) is kept
                # in an overflow list instead of being written into every cell.
                self._large.append(index)
                continue
            for key in keys:
                self.buckets[key].append(index)

    def _key(self, point: Sequence[float]) -> tuple[int, int, int]:
        cell = self.cell
        return (int(math.floor(point[0] / cell)),
                int(math.floor(point[1] / cell)),
                int(math.floor(point[2] / cell)))

    def near(self, point: Sequence[float], radius: float) -> list[int]:
        base = self._key(point)
        reach = max(1, int(math.ceil(radius / self.cell)))
        found: set[int] = set()
        for dx in range(-reach, reach + 1):
            for dy in range(-reach, reach + 1):
                for dz in range(-reach, reach + 1):
                    found.update(self.buckets.get((base[0] + dx, base[1] + dy,
                                                   base[2] + dz), ()))
        found.update(self._large)
        return sorted(found)


def point_triangle_distance(point: Sequence[float], a: Sequence[float],
                            b: Sequence[float], c: Sequence[float]) -> float:
    """Shortest distance from ``point`` to the closed triangle ``a-b-c``.

    Ericson, *Real-Time Collision Detection* §5.1.5 (closest point on a triangle),
    written out rather than delegated so the core package stays dependency free.
    Degenerate triangles fall back to their shortest side.
    """
    ab = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
    ac = (c[0] - a[0], c[1] - a[1], c[2] - a[2])
    ap = (point[0] - a[0], point[1] - a[1], point[2] - a[2])
    d1 = ab[0] * ap[0] + ab[1] * ap[1] + ab[2] * ap[2]
    d2 = ac[0] * ap[0] + ac[1] * ap[1] + ac[2] * ap[2]
    if d1 <= 0.0 and d2 <= 0.0:
        return _norm(ap)
    bp = (point[0] - b[0], point[1] - b[1], point[2] - b[2])
    d3 = ab[0] * bp[0] + ab[1] * bp[1] + ab[2] * bp[2]
    d4 = ac[0] * bp[0] + ac[1] * bp[1] + ac[2] * bp[2]
    if d3 >= 0.0 and d4 <= d3:
        return _norm(bp)
    vc = d1 * d4 - d3 * d2
    if vc <= 0.0 and d1 >= 0.0 and d3 <= 0.0:
        denominator = d1 - d3
        v = d1 / denominator if denominator else 0.0
        return _norm(_combine(point, a, ab, v))
    cp = (point[0] - c[0], point[1] - c[1], point[2] - c[2])
    d5 = ab[0] * cp[0] + ab[1] * cp[1] + ab[2] * cp[2]
    d6 = ac[0] * cp[0] + ac[1] * cp[1] + ac[2] * cp[2]
    if d6 >= 0.0 and d5 <= d6:
        return _norm(cp)
    vb = d5 * d2 - d1 * d6
    if vb <= 0.0 and d2 >= 0.0 and d6 <= 0.0:
        denominator = d2 - d6
        w = d2 / denominator if denominator else 0.0
        return _norm(_combine(point, a, ac, w))
    va = d3 * d6 - d5 * d4
    if va <= 0.0 and (d4 - d3) >= 0.0 and (d5 - d6) >= 0.0:
        denominator = (d4 - d3) + (d5 - d6)
        w = (d4 - d3) / denominator if denominator else 0.0
        edge = (c[0] - b[0], c[1] - b[1], c[2] - b[2])
        return _norm(_combine(point, b, edge, w))
    denominator = va + vb + vc
    if abs(denominator) <= 1e-300:
        # Coplanar-degenerate triangle: the closest point is on one of its sides.
        return min(_point_segment(point, a, b), _point_segment(point, b, c),
                   _point_segment(point, c, a))
    v = vb / denominator
    w = vc / denominator
    closest = (a[0] + ab[0] * v + ac[0] * w,
               a[1] + ab[1] * v + ac[1] * w,
               a[2] + ab[2] * v + ac[2] * w)
    return _norm((point[0] - closest[0], point[1] - closest[1], point[2] - closest[2]))


def _combine(point: Sequence[float], origin: Sequence[float], direction: Sequence[float],
             t: float) -> Vec3:
    return (origin[0] + direction[0] * t - point[0],
            origin[1] + direction[1] * t - point[1],
            origin[2] + direction[2] * t - point[2])


def _norm(vector: Sequence[float]) -> float:
    return math.sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2])


def _point_segment(point: Sequence[float], a: Sequence[float], b: Sequence[float]) -> float:
    distance, _ = project_on_segment(point, a, b)
    return distance


def _median_chord_length(vertices: Sequence[Vec3],
                         edges: Sequence[tuple[int, int]]) -> float:
    if not edges:
        return 0.0
    lengths = sorted(math.dist(vertices[a], vertices[b]) for a, b in edges)
    return lengths[len(lengths) // 2]


def _sample_points(vertices: Sequence[Vec3], a: int, b: int,
                   samples: int) -> list[Vec3]:
    """Interior sample points of a chord (endpoints excluded: they belong to the weld
    cases that :func:`tensorcad.geomquality.stitch.stitch_t_junctions` already covers)."""
    pa, pb = vertices[a], vertices[b]
    points = []
    for index in range(samples):
        t = (index + 1) / (samples + 1)
        points.append((pa[0] + (pb[0] - pa[0]) * t,
                       pa[1] + (pb[1] - pa[1]) * t,
                       pa[2] + (pb[2] - pa[2]) * t))
    return points


def _quantiles(values: Sequence[float]) -> dict[str, float]:
    if not values:
        return {key: math.inf for key in ("q00", "q10", "q50", "q90", "q100")}
    ordered = sorted(values)

    def pick(fraction: float) -> float:
        return ordered[min(len(ordered) - 1, int(round(fraction * (len(ordered) - 1))))]

    return {"q00": pick(0.0), "q10": pick(0.1), "q50": pick(0.5),
            "q90": pick(0.9), "q100": pick(1.0)}


def classify_open_boundary(vertices: Sequence[Vec3], faces: Sequence[Triangle], *,
                           tolerance: float,
                           planes: Sequence[tuple[str, tuple[float, float, float],
                                                  float]] = (),
                           samples: int = 3,
                           worst_limit: int = 8) -> dict[str, Any]:
    """Split the mesh's open boundary into seam / overlap / declared-plane / unexplained.

    Parameters
    ----------
    vertices, faces:
        The tessellated mesh.  It is not modified.
    tolerance:
        Coincidence budget in model units (metres).  A chord is a **seam** when one of
        its interior samples is within ``tolerance`` of a boundary chord it does not
        share a vertex with; it is an **overlap** when the nearest other surface is
        within ``tolerance`` but no *boundary* is (the chord lies inside another face,
        so the faces are stacked instead of butting).
    planes:
        Declared openings.  Each entry is ``(name, (origin, normal), plane_tolerance)``:
        a chord is reported as ``plane:<name>`` when **every** interior sample lies
        within ``plane_tolerance`` of the plane through ``origin`` normal to ``normal``.
        The plane is declarative: it says "this hull is a half breadth model cut at
        y = 0", and the audit then does not call that cut a defect.
    samples:
        Interior sample points per chord (3 is the default; 1 is a midpoint-only scan).
    worst_limit:
        How many of the longest unexplained chords to list, with coordinates, so a
        report can point at the geometry instead of only at a number.

    Returns
    -------
    A dict with the free-edge count, the open-boundary length, the same length split by
    class, distance quantiles for both probes, and the unexplained-chord witnesses.
    Both probes are sampled at interior points of the chord and both require a strictly
    interior foot point on the candidate (a neighbour touched only at a welded vertex is
    neither a seam nor an overlap).
    """
    if tolerance <= 0.0:
        raise ValueError(f"tolerance must be positive, got {tolerance!r}")
    if samples < 1:
        raise ValueError(f"samples must be >= 1, got {samples!r}")
    open_edges = free_edges(faces)
    if not open_edges:
        # A closed mesh: every class is present and zero, so a consumer can index the
        # result without a special case (the shapes must not depend on the mesh).
        empty = {name: {"edges": 0, "length": 0.0} for name in CLASSES}
        return {
            "free_edges": 0,
            "open_length": 0.0,
            "tolerance": tolerance,
            "samples": samples,
            "median_chord_length": 0.0,
            "classes": empty,
            "family": {"seam": dict(empty[CLASS_SEAM]),
                       "overlap": dict(empty[CLASS_OVERLAP]),
                       "declared_plane": {"edges": 0, "length": 0.0},
                       "unexplained": dict(empty[CLASS_UNEXPLAINED])},
            "neighbour_chord_distance": _quantiles([]),
            "other_surface_distance": _quantiles([]),
            "unexplained_worst": [],
        }

    median = _median_chord_length(vertices, open_edges)
    edge_grid = _EdgeGrid(vertices, open_edges, max(median, 4.0 * tolerance, 1e-12))
    triangle_grid = _TriangleGrid(vertices, faces, max(median, 4.0 * tolerance, 1e-12))
    carrier: dict[tuple[int, int], int] = {}
    for index, (a, b, c) in enumerate(faces):
        for u, v in ((a, b), (b, c), (c, a)):
            edge = (u, v) if u < v else (v, u)
            if edge not in carrier:
                carrier[edge] = index
    declared = []
    for entry in planes:
        name, (origin, normal), plane_tolerance = entry
        length = _norm(normal)
        if length <= 0.0:
            raise ValueError(f"plane {name!r} has a zero normal")
        unit = (normal[0] / length, normal[1] / length, normal[2] / length)
        declared.append((CLASS_PLANE_PREFIX + name, origin, unit, plane_tolerance))

    classes: dict[str, dict[str, Any]] = {name: {"edges": 0, "length": 0.0}
                                          for name in CLASSES}
    for plane_name, _origin, _normal, _tolerance in declared:
        classes.setdefault(plane_name, {"edges": 0, "length": 0.0})
    chord_distances: list[float] = []
    surface_distances: list[float] = []
    unexplained: list[tuple[float, Vec3, Vec3]] = []

    free_lookup = {edge: index for index, edge in enumerate(open_edges)}
    for edge in open_edges:
        a, b = edge
        length = math.dist(vertices[a], vertices[b])
        points = _sample_points(vertices, a, b, samples)
        # Probe 1: distance to a boundary chord this one does not touch.
        neighbour_chord = math.inf
        for point in points:
            for candidate in edge_grid.near(point, max(tolerance, 1e-12) * 4.0):
                other = open_edges[candidate]
                if free_lookup[other] == free_lookup[edge]:
                    continue
                # The sample must land *inside* the neighbour chord, not on one of its
                # endpoints.  Two chords that merely meet at a welded vertex are
                # neighbours, not a seam -- and on a coincident-but-differently-sampled
                # curve the long chord of one patch does share both endpoints with the
                # short chords of the other, so an endpoint test (rather than "shares a
                # vertex") is the discriminator that works on both.
                distance, parameter = project_on_segment(point, vertices[other[0]],
                                                         vertices[other[1]])
                if not ENDPOINT_FRACTION < parameter < 1.0 - ENDPOINT_FRACTION:
                    continue
                neighbour_chord = min(neighbour_chord, distance)
                if neighbour_chord <= 0.0:
                    break
        # Probe 2: distance to another face's surface (the carrier face excluded).
        other_surface = math.inf
        own = carrier.get(edge)
        for point in points:
            for candidate in triangle_grid.near(point, max(tolerance, 1e-12) * 4.0):
                if candidate == own:
                    continue
                face = faces[candidate]
                distance = point_triangle_distance(point, vertices[face[0]],
                                                   vertices[face[1]], vertices[face[2]])
                other_surface = min(other_surface, distance)
                if other_surface <= 0.0:
                    break
        chord_distances.append(neighbour_chord)
        surface_distances.append(other_surface)

        label = None
        for plane_name, origin, normal, plane_tolerance in declared:
            if all(abs((point[0] - origin[0]) * normal[0]
                       + (point[1] - origin[1]) * normal[1]
                       + (point[2] - origin[2]) * normal[2]) <= plane_tolerance
                   for point in points):
                label = plane_name
                break
        if label is None:
            if neighbour_chord <= tolerance:
                label = CLASS_SEAM
            elif other_surface <= tolerance:
                label = CLASS_OVERLAP
            else:
                label = CLASS_UNEXPLAINED
                unexplained.append((length, vertices[a], vertices[b]))
        entry = classes.setdefault(label, {"edges": 0, "length": 0.0})
        entry["edges"] += 1
        entry["length"] += length

    unexplained.sort(key=lambda item: -item[0])
    return {
        "free_edges": len(open_edges),
        "open_length": sum(classes[name]["length"] for name in classes),
        "tolerance": tolerance,
        "samples": samples,
        "median_chord_length": median,
        "classes": classes,
        "family": {
            "seam": dict(classes[CLASS_SEAM]),
            "overlap": dict(classes[CLASS_OVERLAP]),
            "declared_plane": {
                "edges": sum(entry["edges"] for name, entry in classes.items()
                             if name.startswith(CLASS_PLANE_PREFIX)),
                "length": sum(entry["length"] for name, entry in classes.items()
                              if name.startswith(CLASS_PLANE_PREFIX)),
            },
            "unexplained": dict(classes[CLASS_UNEXPLAINED]),
        },
        "neighbour_chord_distance": _quantiles(chord_distances),
        "other_surface_distance": _quantiles(surface_distances),
        "unexplained_worst": [
            {"length": length, "from": list(start), "to": list(end)}
            for length, start, end in unexplained[:worst_limit]
        ],
    }


def closure_verdict(audit: dict[str, Any], *, max_unexplained_length: float = 0.0,
                    max_overlap_length: float = 0.0) -> dict[str, Any]:
    """Judge an :func:`classify_open_boundary` result against a closure contract.

    A hull that is meant to be open (a half breadth model) declares its planes; the
    contract then says how many metres of boundary may remain *unexplained* (a real gap
    the modeller must close) and how many may be overlap (stacked patches the trimmer
    must cut).  Sewing a seam is the pipeline's job and never fails the contract -- the
    seam class is what :func:`tensorcad.geomquality.stitch.stitch_t_junctions` consumes.
    """
    unexplained = audit.get("family", {}).get("unexplained", {})
    overlap = audit.get("family", {}).get("overlap", {})
    unexplained_length = float(unexplained.get("length", 0.0))
    overlap_length = float(overlap.get("length", 0.0))
    passed = (unexplained_length <= max_unexplained_length
              and overlap_length <= max_overlap_length)
    return {
        "passed": passed,
        "unexplained_length": unexplained_length,
        "unexplained_budget": max_unexplained_length,
        "overlap_length": overlap_length,
        "overlap_budget": max_overlap_length,
        "verdict": "closed" if passed else "open",
    }