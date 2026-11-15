"""Seam stitching for tessellated CAD patch sets: T-junction resolution.

Why this exists
---------------
A CAD file hands over a hull as *N separate surfaces*.  Each surface is tessellated
in its own parameter space, so two surfaces that share a 3-D boundary curve
generally sample that curve at **different densities**: the coarse side gets one
long chord where the fine side has five short ones.  The result is not a gap in
the geometry -- the curve is the same curve -- but it is a gap in the *mesh*:
the coarse chord has no vertex where the fine chords do, so every edge along the
seam is used once instead of twice and the shell is reported open.

Measured on the publisher's :file:`KCS_hull_SVA.igs` (NMRI T2015, 55 patches):
753 free edges after a 1 mm weld, and **every free-edge endpoint lies exactly on
another free edge** (point-to-segment distance 0 at the 0th, 10th, 50th, 90th and
100th percentile, 2258 vertices tested).  A weld cannot fix that -- welding only
merges vertices that are *coincident*, and these are not; they lie in the middle
of a neighbour's chord.  What has to happen is the coarse chord must be **split at
the neighbour's vertices**, which is what this module does.

The distinction matters for an audit: a metre-wide gap between two different
surfaces is a modelling defect that no stitching may paper over, and this module
must leave it open.  Only a vertex that lies *on* the neighbour's chord (within
``tolerance``) is used, and the split never moves an existing vertex, so the
surface is unchanged and the closed shell encloses exactly the same volume.

Algorithm
---------
1. index the free (once-used) edges of the mesh in a uniform grid;
2. for every free-edge endpoint, query the grid for free edges it lies on, with
   the projection parameter strictly inside the chord (a vertex at an endpoint is
   a weld case, not a T-junction);
3. split the chord at all such vertices at once, re-triangulating the single face
   that owned the chord (a free edge belongs to exactly one face by definition);
4. repeat, because a split can expose a second-level T-junction on a shorter edge.

Splitting reuses the *existing* vertex index of the neighbour patch, so the seam
becomes topological immediately -- no tolerance-based vertex merging is involved
and the operation is idempotent.
"""
from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Sequence
from typing import Any

Vec3 = tuple[float, float, float]
Triangle = tuple[int, int, int]

#: A vertex closer than this fraction of its chord's length to a chord endpoint is
#: treated as an endpoint (weld case) instead of a split; splitting there would
#: create a sliver triangle with no area.
ENDPOINT_FRACTION = 1e-9


def edge_use_counts(faces: Sequence[Triangle]) -> dict[tuple[int, int], int]:
    """How many triangles use each undirected edge."""
    counts: dict[tuple[int, int], int] = defaultdict(int)
    for a, b, c in faces:
        for u, v in ((a, b), (b, c), (c, a)):
            counts[(u, v) if u < v else (v, u)] += 1
    return dict(counts)


def free_edges(faces: Sequence[Triangle]) -> list[tuple[int, int]]:
    """Edges used exactly once.  Their count is the mesh's free-edge count."""
    return [edge for edge, count in edge_use_counts(faces).items() if count == 1]


def _sub(a: Sequence[float], b: Sequence[float]) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _dot(a: Sequence[float], b: Sequence[float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _norm(a: Sequence[float]) -> float:
    return math.sqrt(_dot(a, a))


def project_on_segment(point: Sequence[float], a: Sequence[float],
                       b: Sequence[float]) -> tuple[float, float]:
    """Return ``(distance, t)``: distance from ``point`` to segment ``a-b``, and the
    parameter of the closest point (0 at ``a``, 1 at ``b``)."""
    ab = _sub(b, a)
    length2 = _dot(ab, ab)
    if length2 <= 0.0:
        return _norm(_sub(point, a)), 0.0
    t = _dot(_sub(point, a), ab) / length2
    clamped = 0.0 if t < 0.0 else (1.0 if t > 1.0 else t)
    closest = (a[0] + ab[0] * clamped, a[1] + ab[1] * clamped, a[2] + ab[2] * clamped)
    return _norm(_sub(point, closest)), t


class _EdgeGrid:
    """Uniform spatial hash of free edges, for point-to-chord queries."""

    def __init__(self, vertices: Sequence[Vec3], edges: Sequence[tuple[int, int]],
                 cell: float) -> None:
        self.vertices = vertices
        self.cell = cell
        self.buckets: dict[tuple[int, int, int], list[int]] = defaultdict(list)
        self.edges = list(edges)
        for index, (a, b) in enumerate(self.edges):
            for key in self._cells_between(vertices[a], vertices[b]):
                self.buckets[key].append(index)

    def _key(self, point: Sequence[float]) -> tuple[int, int, int]:
        cell = self.cell
        return (int(math.floor(point[0] / cell)),
                int(math.floor(point[1] / cell)),
                int(math.floor(point[2] / cell)))

    def _cells_between(self, a: Sequence[float], b: Sequence[float]) -> list[tuple[int, int, int]]:
        """Cells covered by the segment's bounding box, sampled along the segment.

        Sampling (rather than enumerating the whole box) keeps a long diagonal chord
        from touching O(n^3) cells, and is exact enough because a query only needs the
        cell that contains the chord *point* nearest the vertex -- and any candidate at
        distance <= ``tolerance`` is within one cell of the chord when ``cell`` is at
        least twice ``tolerance``.
        """
        keys = {self._key(a), self._key(b)}
        length = _norm(_sub(b, a))
        steps = min(int(length / (0.5 * self.cell)) + 1, 4096)
        for step in range(1, steps):
            t = step / steps
            keys.add(self._key((a[0] + (b[0] - a[0]) * t,
                                a[1] + (b[1] - a[1]) * t,
                                a[2] + (b[2] - a[2]) * t)))
        return list(keys)

    def near(self, point: Sequence[float], radius: float) -> list[int]:
        """Edge indices whose bucket is within ``radius`` of ``point``."""
        base = self._key(point)
        reach = max(1, int(math.ceil(radius / self.cell)))
        found: set[int] = set()
        for dx in range(-reach, reach + 1):
            for dy in range(-reach, reach + 1):
                for dz in range(-reach, reach + 1):
                    found.update(self.buckets.get((base[0] + dx, base[1] + dy,
                                                   base[2] + dz), ()))
        return sorted(found)


def _face_along_edge(faces: Sequence[Triangle], a: int, b: int) -> tuple[int, int] | None:
    """``(face index, 1 if the face traverses a->b else -1)`` for a once-used edge."""
    for index, face in enumerate(faces):
        for u, v in ((face[0], face[1]), (face[1], face[2]), (face[2], face[0])):
            if (u == a and v == b) or (u == b and v == a):
                return index, (1 if (u == a and v == b) else -1)
    return None


def _median_edge_length(vertices: Sequence[Vec3], edges: Sequence[tuple[int, int]]) -> float:
    lengths = sorted(_norm(_sub(vertices[b], vertices[a])) for a, b in edges)
    if not lengths:
        return 0.0
    return lengths[len(lengths) // 2]


def _third_vertex(face: Triangle, a: int, b: int) -> int | None:
    """Return the third vertex index of a triangle given two of its edge vertices.

    Returns ``None`` if the face is degenerate (all three vertices are ``a`` or ``b``),
    which can happen when a triangle has collapsed to a line or point.
    """
    for index in face:
        if index != a and index != b:
            return index
    return None


def stitch_t_junctions(vertices: Sequence[Vec3], faces: Sequence[Triangle], *,
                       tolerance: float, max_passes: int = 4,
                       max_split_fraction: float = 0.5) -> tuple[list[Vec3], list[Triangle], dict[str, Any]]:
    """Split free chords at neighbour vertices that lie on them.

    Parameters
    ----------
    vertices, faces:
        The tessellated mesh.  Vertices are returned unchanged (this stage never
        moves geometry, it only adds index-sharing).
    tolerance:
        A vertex is on a chord when its distance to the chord is ``<= tolerance``.
        Use a small physical length (the seam coincidence the CAD export guarantees),
        *not* a weld tolerance: a metre-wide gap must stay open.
    max_passes:
        A split on one seam can expose a T-junction on the shorter neighbour chords;
        each pass removes one level.
    max_split_fraction:
        Fail-closed ceiling on ``splits_total / len(faces)``.  A budget that is large
        against the actual seam offsets does not close the model -- it *cascades*: every
        snapped vertex becomes a new candidate for a neighbour chord, and the triangle
        count grows without closing the seam (measured on ``jbc-boss.igs``: a 50 mm budget
        produced 19 146 splits for 894 input triangles and 10 951 fragments).  Exceeding
        the ceiling aborts the stitch and records ``aborted`` instead of returning a
        shattered mesh.

    Returns
    -------
    ``(vertices, faces, record)`` where ``record`` carries the numbers that decided
    the operation, so a report can show what changed instead of claiming it worked.
    """
    out_vertices = list(vertices)
    out_faces = [tuple(face) for face in faces]
    record: dict[str, Any] = {
        "tolerance": tolerance,
        "free_edges_before": len(free_edges(out_faces)),
        "passes": [],
        "splits_total": 0,
    }
    if tolerance <= 0.0 or not out_faces:
        record["free_edges_after"] = record["free_edges_before"]
        record["skipped"] = "no tolerance or empty mesh"
        return out_vertices, out_faces, record

    for pass_index in range(max_passes):
        current_free = free_edges(out_faces)
        if not current_free:
            break
        median = _median_edge_length(out_vertices, current_free)
        cell = max(median, 4.0 * tolerance, 1e-12)
        grid = _EdgeGrid(out_vertices, current_free, cell)

        boundary_vertices = sorted({index for edge in current_free for index in edge})
        splits: dict[tuple[int, int], set[int]] = defaultdict(set)
        for vertex in boundary_vertices:
            point = out_vertices[vertex]
            for edge_index in grid.near(point, tolerance):
                a, b = current_free[edge_index]
                if vertex in (a, b):
                    continue
                distance, t = project_on_segment(point, out_vertices[a], out_vertices[b])
                if distance > tolerance:
                    continue
                chord = _norm(_sub(out_vertices[b], out_vertices[a]))
                if chord <= 0.0:
                    continue
                # The split must be strictly interior *in parameter*: a vertex that lands on
                # an endpoint of the chord is the weld case, and splitting there would build
                # a triangle of zero area.  The physical proximity that ``tolerance`` already
                # granted is not used to widen this margin, or short chords (a 5 cm chord is
                # ordinary at a stern knuckle) would lose every candidate to a
                # tolerance/length ratio above 0.5.
                if not ENDPOINT_FRACTION < t < 1.0 - ENDPOINT_FRACTION:
                    continue
                splits[(a, b)].add(vertex)

        if not splits:
            record["passes"].append({"pass": pass_index, "free_edges": len(current_free),
                                     "splits": 0, "note": "no vertex lies on a free chord"})
            break

        applied = 0
        edge_splits: list[dict[str, Any]] = []
        for (a, b), candidates in sorted(splits.items()):
            located = _face_along_edge(out_faces, a, b)
            if located is None:
                continue
            face_index, direction = located
            face = out_faces[face_index]
            third = _third_vertex(face, a, b)
            if third is None:
                # Degenerate face: all three vertices are a or b. Skip splitting
                # as there is no valid third vertex to form a triangle with.
                continue
            ordered = sorted(candidates,
                             key=lambda vertex: project_on_segment(out_vertices[vertex],
                                                                   out_vertices[a],
                                                                   out_vertices[b])[1])
            # Collapse vertices that project to (almost) the same parameter: they are the
            # same seam point exported twice, and splitting twice would make a sliver.
            chain: list[int] = []
            last_t = -1.0
            for vertex in ordered:
                _, t = project_on_segment(out_vertices[vertex], out_vertices[a], out_vertices[b])
                if chain and abs(t - last_t) <= ENDPOINT_FRACTION:
                    continue
                chain.append(vertex)
                last_t = t
            if not chain:
                continue
            path = [a, *chain, b] if direction == 1 else [b, *chain, a]
            replacement = [(path[i], path[i + 1], third) for i in range(len(path) - 1)]
            out_faces[face_index] = replacement[0]
            out_faces.extend(tuple(triangle) for triangle in replacement[1:])
            applied += len(chain)
            edge_splits.append({"edge": [a, b], "face": face_index, "vertices": chain})
        record["passes"].append({"pass": pass_index, "free_edges": len(current_free),
                                 "splits": applied, "detail": edge_splits[:12],
                                 "detail_truncated": max(0, len(edge_splits) - 12)})
        record["splits_total"] += applied
        if record["splits_total"] > max_split_fraction * len(faces):
            record["aborted"] = (f"splits_total {record['splits_total']} exceeds "
                                 f"budget {max_split_fraction * len(faces):.0f}")
            break

    record["free_edges_after"] = len(free_edges(out_faces))
    return out_vertices, out_faces, record


def seam_report(vertices: Sequence[Vec3], faces: Sequence[Triangle], *,
                tolerance: float) -> dict[str, Any]:
    """Measure how much of the open boundary is a T-junction seam.

    For every free-edge endpoint, find the nearest free chord it does not belong to.
    A population clustered at 0 is a stitching problem; a population at metres is a
    geometric gap that must stay open and be reported as such.
    """
    open_edges = free_edges(faces)
    if not open_edges:
        return {"free_edges": 0, "t_junctions": 0, "residual_free_endpoints": 0}
    median = _median_edge_length(vertices, open_edges)
    grid = _EdgeGrid(vertices, open_edges, max(median, 4.0 * tolerance, 1e-12))
    distances: list[float] = []
    for vertex in sorted({index for edge in open_edges for index in edge}):
        point = vertices[vertex]
        best = math.inf
        for edge_index in grid.near(point, max(median, 10.0 * tolerance)):
            a, b = open_edges[edge_index]
            if vertex in (a, b):
                continue
            distance, _ = project_on_segment(point, vertices[a], vertices[b])
            best = min(best, distance)
        if best < math.inf:
            distances.append(best)
    distances.sort()

    def percentile(fraction: float) -> float:
        if not distances:
            return math.inf
        position = min(len(distances) - 1, int(round(fraction * (len(distances) - 1))))
        return distances[position]

    t_junctions = sum(1 for value in distances if value <= tolerance)
    return {
        "free_edges": len(open_edges),
        "free_edge_endpoints": len(distances),
        "t_junctions": t_junctions,
        "residual_free_endpoints": len(distances) - t_junctions,
        "distance_to_nearest_free_chord_quantiles": {
            "q00": percentile(0.0), "q10": percentile(0.1), "q50": percentile(0.5),
            "q90": percentile(0.9), "q100": percentile(1.0),
        },
        "tolerance": tolerance,
    }
