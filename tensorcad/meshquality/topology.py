"""The half of mesh quality that is not per-element: closure, orientation, and faces.

``tensorcad.meshquality.metrics`` answers "how good is each cell". That is only half
of what a solver rejects a mesh for. The other half lives on the mesh as a whole and on
its faces:

* **closure** -- is the surface watertight, edge-manifold, vertex-manifold, and does it
  have duplicate faces? A closed surface is the precondition for a volume, a pressure
  load, or a boundary-element solve;
* **orientation** -- are the faces wound consistently, and does the consistent winding
  point *out* of the solid? An inverted but otherwise perfect hull meshes fine and then
  computes a negative displacement;
* **face coupling** -- for a volume mesh, the non-orthogonality of every internal face,
  the OpenFOAM ``checkMesh`` quantity. A previous repository audit found that mesh
  quality was measured in three vocabularies and none of them checked any of this.

Conventions (stated once, tested everywhere below):

* A surface is **watertight** when every undirected edge is used by exactly two faces.
  An edge used once is a *boundary* edge; used more than twice, a *non-manifold* edge.
* **Winding consistency** means the two faces on each interior edge traverse it in
  opposite directions. **Outward** means consistent winding *and* a positive signed
  volume from the divergence theorem, ``sum_i p0 . (p1 x p2) / 6``.
* **Euler characteristic** ``V - E + F`` uses only vertices an incident face touches;
  for a closed orientable shell ``chi = 2 - 2 * genus``.
* **Non-orthogonality** of an internal face is ``theta = acos(n_hat . d_hat)`` in
  degrees, where ``n_hat`` is the face unit normal flipped to point from the
  lower-indexed cell to the higher and ``d_hat`` is the unit centre-to-centre vector.
  ``0 deg`` is an orthogonal face; ``45 deg`` for a uniformly sheared hexahedral stack
  is the closed form used in the tests. Cell and face centres are node centroids, which
  is exact for the simplices and for affine hexahedra and is what makes the sheared
  stack an exact check.

Nothing here imports numpy or a mesh library, and nothing in the package core imports
this module's external oracle (:func:`tensorcad.meshquality.reference.trimesh_topology`).
"""
from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass
from typing import Sequence

from .metrics import _angle_deg, _cross, _sub
from .report import Check

__all__ = [
    "CELL_FACES",
    "FaceNonOrthogonality",
    "SurfaceTopology",
    "gate_surface",
    "non_orthogonality",
    "surface_topology",
]

Vec3 = tuple[float, float, float]

#: Face node tables, in the winding that points out of a positively oriented cell.
CELL_FACES: dict[str, tuple[tuple[int, ...], ...]] = {
    "tet": ((0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)),
    "hex": ((0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)),
}

#: Nodes per cell, which is not the nodes per face of the table above.
CELL_NODES = {"tet": 4, "hex": 8}


# --------------------------------------------------------------------------------------
# surface topology
# --------------------------------------------------------------------------------------
@dataclass(frozen=True)
class SurfaceTopology:
    """Whole-mesh topological facts about a triangle surface."""

    points: int
    faces: int
    used_points: int
    edges: int
    boundary_edges: int
    non_manifold_edges: int
    duplicate_faces: int
    non_manifold_vertices: int
    disconnected_components: int
    winding_consistent: bool
    signed_volume: float
    euler_characteristic: int

    @property
    def watertight(self) -> bool:
        """Every edge used exactly twice -- a closed surface."""
        return self.faces > 0 and self.boundary_edges == 0 and self.non_manifold_edges == 0

    @property
    def edge_manifold(self) -> bool:
        return self.faces > 0 and self.non_manifold_edges == 0

    @property
    def vertex_manifold(self) -> bool:
        return self.non_manifold_vertices == 0

    @property
    def outward(self) -> bool:
        """Closed, consistently wound, and enclosing a positive volume."""
        return self.watertight and self.winding_consistent and self.signed_volume > 0.0

    @property
    def genus(self) -> float:
        """Genus per shell; only meaningful for closed orientable surfaces."""
        if self.disconnected_components == 0:
            return math.nan
        return (2.0 * self.disconnected_components - self.euler_characteristic) / 2.0

    def problems(self) -> tuple[str, ...]:
        """Every defect, as a sentence a CI log can print."""
        found: list[str] = []
        if self.faces == 0:
            found.append("no faces")
            return tuple(found)
        if self.boundary_edges:
            found.append(f"{self.boundary_edges} boundary edge(s): surface is open")
        if self.non_manifold_edges:
            found.append(f"{self.non_manifold_edges} edge(s) used by more than two faces")
        if self.duplicate_faces:
            found.append(f"{self.duplicate_faces} duplicate face(s)")
        if self.non_manifold_vertices:
            found.append(f"{self.non_manifold_vertices} vertex/vertices with a pinched face fan")
        if not self.winding_consistent:
            found.append("inconsistent face winding")
        if self.watertight and self.winding_consistent and self.signed_volume <= 0.0:
            found.append(f"signed volume {self.signed_volume:.6g} <= 0: faces point inward")
        return tuple(found)

    def to_dict(self) -> dict[str, object]:
        return {
            "points": self.points,
            "faces": self.faces,
            "used_points": self.used_points,
            "edges": self.edges,
            "boundary_edges": self.boundary_edges,
            "non_manifold_edges": self.non_manifold_edges,
            "duplicate_faces": self.duplicate_faces,
            "non_manifold_vertices": self.non_manifold_vertices,
            "disconnected_components": self.disconnected_components,
            "winding_consistent": self.winding_consistent,
            "signed_volume": self.signed_volume,
            "euler_characteristic": self.euler_characteristic,
            "genus": self.genus,
            "watertight": self.watertight,
            "edge_manifold": self.edge_manifold,
            "vertex_manifold": self.vertex_manifold,
            "outward": self.outward,
            "problems": list(self.problems()),
        }


def _components(nodes: Sequence[int], adjacency: dict[int, set[int]]) -> int:
    """Number of connected components of a small undirected graph."""
    seen: set[int] = set()
    count = 0
    for start in nodes:
        if start in seen or start not in adjacency:
            continue
        count += 1
        stack = [start]
        seen.add(start)
        while stack:
            node = stack.pop()
            for other in adjacency[node]:
                if other not in seen:
                    seen.add(other)
                    stack.append(other)
    return count


def surface_topology(
    points: Sequence[Sequence[float]], triangles: Sequence[Sequence[int]]
) -> SurfaceTopology:
    """Closure, manifoldness, orientation, Euler characteristic and signed volume."""
    faces = [tuple(int(i) for i in face) for face in triangles]
    for face in faces:
        if len(face) != 3:
            raise ValueError(f"surface_topology needs triangles, got a face with {len(face)} nodes")
        for index in face:
            if index < 0 or index >= len(points):
                raise ValueError(f"face {face} indexes outside {len(points)} points")

    directed: dict[tuple[int, int], int] = defaultdict(int)
    face_edges: list[tuple[tuple[int, int], ...]] = []
    for face in faces:
        a, b, c = face
        edges = ((a, b), (b, c), (c, a))
        face_edges.append(edges)
        for edge in edges:
            directed[edge] += 1

    undirected_counts: dict[tuple[int, int], int] = defaultdict(int)
    for (a, b), count in directed.items():
        undirected_counts[(a, b) if a <= b else (b, a)] += count

    boundary = sum(1 for count in undirected_counts.values() if count == 1)
    non_manifold_edges = sum(1 for count in undirected_counts.values() if count > 2)

    # Winding: on an interior edge the two faces must traverse it in opposite senses.
    winding_consistent = all(
        directed[(a, b)] == directed[(b, a)]
        for (a, b) in undirected_counts
        if undirected_counts[(a, b)] == 2
    )

    keyed: dict[tuple[int, ...], int] = defaultdict(int)
    for face in faces:
        keyed[tuple(sorted(face))] += 1
    duplicate_faces = sum(count - 1 for count in keyed.values() if count > 1)

    # Vertex manifoldness: the faces around a vertex must form one fan when joined
    # through the edges that meet at that vertex.
    incident: dict[int, list[int]] = defaultdict(list)
    for index, face in enumerate(faces):
        for vertex in face:
            incident[vertex].append(index)
    non_manifold_vertices = 0
    for vertex, face_ids in incident.items():
        if len(face_ids) < 2:
            continue
        local = {fid: fid for fid in face_ids}

        def find(x: int) -> int:
            while local[x] != x:
                local[x] = local[local[x]]
                x = local[x]
            return x

        by_edge: dict[tuple[int, int], list[int]] = defaultdict(list)
        for fid in face_ids:
            face = faces[fid]
            for i in range(3):
                a, b = face[i], face[(i + 1) % 3]
                if a == vertex or b == vertex:
                    by_edge[(min(a, b), max(a, b))].append(fid)
        for sharing in by_edge.values():
            for other in sharing[1:]:
                local[find(other)] = find(sharing[0])
        roots = {find(fid) for fid in face_ids}
        if len(roots) > 1:
            non_manifold_vertices += 1

    # Shells: faces joined through shared edges (topological connectivity).
    edge_faces: dict[tuple[int, int], list[int]] = defaultdict(list)
    for index, edges in enumerate(face_edges):
        for a, b in edges:
            edge_faces[(min(a, b), max(a, b))].append(index)
    adjacency: dict[int, set[int]] = defaultdict(set)
    for sharing in edge_faces.values():
        for other in sharing[1:]:
            adjacency[sharing[0]].add(other)
            adjacency[other].add(sharing[0])
    components = _components(list(range(len(faces))), adjacency) if faces else 0

    used = sorted({index for face in faces for index in face})
    signed_volume = 0.0
    for face in faces:
        p0 = tuple(float(x) for x in points[face[0]])
        p1 = tuple(float(x) for x in points[face[1]])
        p2 = tuple(float(x) for x in points[face[2]])
        signed_volume += (
            p0[0] * (p1[1] * p2[2] - p1[2] * p2[1])
            + p0[1] * (p1[2] * p2[0] - p1[0] * p2[2])
            + p0[2] * (p1[0] * p2[1] - p1[1] * p2[0])
        ) / 6.0

    return SurfaceTopology(
        points=len(points),
        faces=len(faces),
        used_points=len(used),
        edges=len(undirected_counts),
        boundary_edges=boundary,
        non_manifold_edges=non_manifold_edges,
        duplicate_faces=duplicate_faces,
        non_manifold_vertices=non_manifold_vertices,
        disconnected_components=components,
        winding_consistent=bool(winding_consistent),
        signed_volume=signed_volume,
        euler_characteristic=len(used) - len(undirected_counts) + len(faces),
    )


def gate_surface(
    topology: SurfaceTopology,
    *,
    require_watertight: bool = True,
    require_manifold: bool = True,
    require_outward: bool = True,
    max_genus: float | None = None,
) -> tuple[Check, ...]:
    """Turn the topological facts into :class:`Check` records a gate report can carry."""

    def check(name: str, passed: bool, value: float, limit: str, detail: str) -> Check:
        return Check(
            name=name,
            kind="surface",
            metric=name.split(".")[-1],
            passed=bool(passed),
            value=float(value),
            limit=limit,
            detail=detail,
        )

    checks: list[Check] = []
    if require_watertight:
        checks.append(
            check(
                "surface.watertight",
                topology.watertight,
                1.0 if topology.watertight else 0.0,
                "== True",
                f"{topology.boundary_edges} boundary edge(s), "
                f"{topology.non_manifold_edges} non-manifold edge(s), "
                f"{topology.faces} face(s)",
            )
        )
    if require_manifold:
        checks.append(
            check(
                "surface.edge_manifold",
                topology.edge_manifold,
                1.0 if topology.edge_manifold else 0.0,
                "== True",
                f"{topology.non_manifold_edges} edge(s) with more than two faces",
            )
        )
        checks.append(
            check(
                "surface.vertex_manifold",
                topology.vertex_manifold,
                1.0 if topology.vertex_manifold else 0.0,
                "== True",
                f"{topology.non_manifold_vertices} pinched vertex/vertices",
            )
        )
    checks.append(
        check(
            "surface.winding",
            topology.winding_consistent,
            1.0 if topology.winding_consistent else 0.0,
            "== True",
            "opposite traversal on every interior edge",
        )
    )
    if require_outward:
        checks.append(
            check(
                "surface.outward",
                topology.outward,
                1.0 if topology.outward else 0.0,
                "== True",
                f"signed volume {topology.signed_volume:.9g} m^3 over "
                f"{topology.disconnected_components} shell(s)",
            )
        )
    if max_genus is not None:
        checks.append(
            check(
                "surface.genus",
                topology.watertight and topology.genus <= max_genus,
                topology.genus,
                f"<= {max_genus:g}",
                f"chi = {topology.euler_characteristic}, "
                f"{topology.disconnected_components} shell(s)",
            )
        )
    return tuple(checks)


# --------------------------------------------------------------------------------------
# internal-face non-orthogonality (OpenFOAM checkMesh)
# --------------------------------------------------------------------------------------
@dataclass(frozen=True)
class FaceNonOrthogonality:
    """Non-orthogonality of every internal face of a tetrahedral or hexahedral mesh."""

    kind: str
    cells: int
    internal_faces: int
    boundary_faces: int
    degenerate_faces: int
    angles: tuple[float, ...]

    @property
    def maximum(self) -> float:
        return max(self.angles) if self.angles else math.nan

    @property
    def mean(self) -> float:
        return sum(self.angles) / len(self.angles) if self.angles else math.nan

    def over(self, limit: float) -> tuple[float, ...]:
        """The faces worse than ``limit`` degrees, sorted worst first."""
        return tuple(sorted((a for a in self.angles if a > limit), reverse=True))

    def to_dict(self) -> dict[str, object]:
        return {
            "kind": self.kind,
            "cells": self.cells,
            "internal_faces": self.internal_faces,
            "boundary_faces": self.boundary_faces,
            "degenerate_faces": self.degenerate_faces,
            "max_deg": self.maximum,
            "mean_deg": self.mean,
            "angles_deg": list(self.angles),
        }


def _centroid(points: Sequence[Vec3], nodes: Sequence[int]) -> Vec3:
    n = float(len(nodes))
    return (
        sum(points[i][0] for i in nodes) / n,
        sum(points[i][1] for i in nodes) / n,
        sum(points[i][2] for i in nodes) / n,
    )


def _newell_normal(points: Sequence[Vec3]) -> Vec3:
    normal = (0.0, 0.0, 0.0)
    n = len(points)
    for i in range(n):
        a, b = points[i], points[(i + 1) % n]
        normal = (
            normal[0] + (a[1] - b[1]) * (a[2] + b[2]),
            normal[1] + (a[2] - b[2]) * (a[0] + b[0]),
            normal[2] + (a[0] - b[0]) * (a[1] + b[1]),
        )
    return normal


def non_orthogonality(
    points: Sequence[Sequence[float]], cells: Sequence[Sequence[int]], kind: str
) -> FaceNonOrthogonality:
    """Angle between each internal face normal and the line joining its two cell centres.

    ``kind`` is ``"tet"`` or ``"hex"``; the face tables are :data:`CELL_FACES`. Faces used
    by one cell are boundary faces and are counted, not measured (a boundary face has no
    second centre). A face used by three or more cells raises, because a non-conforming
    mesh has no unique centre-to-centre line and silently averaging it would hide exactly
    the defect this metric exists to find.
    """
    if kind not in CELL_FACES:
        raise ValueError(f"non_orthogonality supports {sorted(CELL_FACES)}, got {kind!r}")
    table = CELL_FACES[kind]
    pts: list[Vec3] = [tuple(float(x) for x in p) for p in points]  # type: ignore[misc]
    shapes = {len(cell) for cell in cells}
    expected = CELL_NODES[kind]
    if shapes and shapes != {expected}:
        raise ValueError(f"{kind} cells need {expected} nodes, got sizes {sorted(shapes)}")

    owners: dict[tuple[int, ...], list[int]] = defaultdict(list)
    face_order: dict[tuple[int, ...], tuple[int, ...]] = {}
    for index, cell in enumerate(cells):
        for node_set in table:
            face = tuple(int(cell[i]) for i in node_set)
            key = tuple(sorted(face))
            owners[key].append(index)
            face_order.setdefault(key, face)

    centres = [_centroid(pts, cell) for cell in cells]
    angles: list[float] = []
    degenerate = 0
    boundary = 0
    for key, sharing in owners.items():
        if len(sharing) == 1:
            boundary += 1
            continue
        if len(sharing) > 2:
            raise ValueError(f"face {key} is shared by {len(sharing)} cells: non-conforming mesh")
        low, high = sorted(sharing)
        # A face normal needs the polygon's cyclic order, not the sorted key: Newell on a
        # sorted node list can return a bow-tie normal (e.g. (2, 2, 0) for a unit square).
        face_normal = _newell_normal([pts[i] for i in face_order[key]])
        d = _sub(centres[high], centres[low])
        # Orient the normal from the lower-indexed cell towards the higher one.
        if sum(face_normal[i] * d[i] for i in range(3)) < 0.0:
            face_normal = (-face_normal[0], -face_normal[1], -face_normal[2])
        angle = _angle_deg(face_normal, d)
        if angle != angle:  # NaN: a collapsed face or two coincident cell centres
            degenerate += 1  # counted, not dropped: duplicate cells must fail a gate
        else:
            angles.append(angle)

    return FaceNonOrthogonality(
        kind=kind,
        cells=len(cells),
        internal_faces=len(angles),
        boundary_faces=boundary,
        degenerate_faces=degenerate,
        angles=tuple(angles),
    )


def gate_non_orthogonality(
    faces: FaceNonOrthogonality, maximum: float = 70.0
) -> tuple[Check, ...]:
    """One check: the worst internal face must not exceed ``maximum`` degrees."""
    worst = faces.maximum
    passed = bool(faces.internal_faces) and worst <= maximum and faces.degenerate_faces == 0
    if faces.degenerate_faces:
        detail = (
            f"{faces.degenerate_faces} degenerate internal face(s) (coincident cell centres "
            f"or collapsed faces); max {worst:.6g} deg over {faces.internal_faces} face(s)"
        )
    elif faces.internal_faces:
        detail = f"max {worst:.6g} deg over {faces.internal_faces} internal face(s)"
    else:
        detail = "no internal faces (single cell or disconnected cells)"
    return (
        Check(
            name="volume.max_non_orthogonality",
            kind="volume",
            metric="max_non_orthogonality",
            passed=passed,
            value=worst,
            limit=f"<= {maximum:g} deg",
            detail=detail,
        ),
    )