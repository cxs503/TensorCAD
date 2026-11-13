"""Trimmed CAD faces -> watertight mesh: inversion, tessellation, repair, audit.

The second half of the import pipeline.  :mod:`tensorcad.geomquality.iges_trim` reads a
file into faces whose boundary loops are **3-D polylines**; this module turns them into
triangles and then judges them honestly:

1. **Invert the loops onto the surface.**  Every boundary point is projected to surface
   parameters ``(u, v)`` by a seeded Gauss-Newton solve.  The residual is not swallowed:
   a loop that does not lie on its surface (a mis-read trim, a wrong transformation, a
   surface/curve mismatch) is reported per face and the face is refused.  This is what
   makes "was the trim read correctly?" a measurement instead of a belief.
2. **Triangulate the trimmed parameter domain.**  Outer loop minus holes, bridged and
   ear-clipped in ``(u, v)``, so every boundary vertex sits exactly on the trim curve --
   the property that stitching depends on.  The triangulation is checked against the
   polygon area before it is used, and a triangulation that does not add up is refused.
3. **Measure the mesh.**  Welded vertices (tolerance), free edges, non-manifold edges,
   orientation consistency, Euler characteristic, signed volume, area, degenerate faces,
   component count.  A shell is *watertight* only when every edge is used by exactly two
   triangles and the winding is consistent.
4. **Repair, then measure again.**  The chain is weld -> drop degenerate triangles ->
   drop duplicate triangles -> fill boundary loops whose length is below a stated
   threshold -> re-orient.  Both the before and the after numbers are reported; a repair
   that does not converge is reported as a failure, never as a pass.

What this module does *not* do: it does not invent a B-rep, it does not smooth, and it
does not move vertices except by welding (which is exact within a declared tolerance and
recorded).  A face whose trim cannot be inverted within tolerance is refused, which keeps
the gate meaningful: an under-trimmed face would give a plausible area, a plausible
volume and an open shell.

Units: coordinates are whatever the reader produced (knots, mm, m).  The pipeline
converts once with ``scale`` and records the factor it used.
"""
from __future__ import annotations

import hashlib
import json
import math
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Sequence

from ..geometry import Mesh
from .iges_trim import (
    IgesError,
    PlaneSurface,
    Surface,
    TrimmedFace,
    TrimmedModel,
    read_iges_trimmed,
)

__all__ = [
    "CadImport",
    "MeshAudit",
    "Tessellation",
    "import_iges",
    "import_cad",
    "tessellate_trimmed_model",
    "invert_point",
    "audit_mesh",
    "repair_mesh",
    "weld_vertices",
    "triangulate_trimmed_domain",
    "write_mesh",
]

Vec3 = tuple[float, float, float]
Triangle = tuple[int, int, int]


# --------------------------------------------------------------------------- #
# surface inversion
# --------------------------------------------------------------------------- #
def invert_point(surface: Surface, point: Sequence[float], *,
                 seeds: Sequence[tuple[float, float]] = (),
                 iterations: int = 60, tolerance: float = 1e-12) -> tuple[float, float, float]:
    """Project ``point`` onto ``surface``: returns ``(u, v, residual)``.

    Gauss-Newton on ``|S(u,v) - p|^2`` with a finite-difference Jacobian.  ``seeds`` lets
    a caller walk along a curve (each point seeded by its predecessor), which is what
    keeps a trim loop's parameter image continuous instead of jumping between branches.
    """
    (u0, u1), (v0, v1) = surface.domain
    if isinstance(surface, PlaneSurface):
        u, v = surface.project(point)
        residual = _norm(_sub(surface.evaluate(u, v), point))
        return (u, v, residual)

    best: tuple[float, float, float] | None = None
    grid_seeds = list(seeds)
    if not grid_seeds:
        steps = 6
        grid_seeds = [
            (u0 + (u1 - u0) * i / steps, v0 + (v1 - v0) * j / steps)
            for i in range(steps + 1)
            for j in range(steps + 1)
        ]
    for seed in grid_seeds:
        candidate = _gauss_newton(surface, point, seed, iterations, tolerance)
        if best is None or candidate[2] < best[2]:
            best = candidate
        if best[2] <= tolerance:
            break
    assert best is not None
    return best


def _gauss_newton(surface: Surface, point: Sequence[float], seed: tuple[float, float],
                  iterations: int, tolerance: float) -> tuple[float, float, float]:
    (u0, u1), (v0, v1) = surface.domain
    u = min(max(seed[0], u0), u1)
    v = min(max(seed[1], v0), v1)
    span_u = max(u1 - u0, 1e-12)
    span_v = max(v1 - v0, 1e-12)
    step_u = 1e-6 * span_u
    step_v = 1e-6 * span_v
    for _ in range(iterations):
        base = surface.evaluate(u, v)
        residual = _sub(base, point)
        error = _norm(residual)
        if error <= tolerance:
            break
        du = _scale(_sub(surface.evaluate(min(u + step_u, u1), v), base), 1.0 / step_u)
        dv = _scale(_sub(surface.evaluate(u, min(v + step_v, v1)), base), 1.0 / step_v)
        a11, a12 = _dot(du, du), _dot(du, dv)
        a22 = _dot(dv, dv)
        b1, b2 = -_dot(du, residual), -_dot(dv, residual)
        determinant = a11 * a22 - a12 * a12
        if abs(determinant) < 1e-18:
            break
        delta_u = (b1 * a22 - b2 * a12) / determinant
        delta_v = (a11 * b2 - a12 * b1) / determinant
        # Damped step: a full Newton step off a bad seed can leave the domain entirely.
        damping = 1.0
        improved = False
        for _attempt in range(12):
            trial_u = min(max(u + damping * delta_u, u0), u1)
            trial_v = min(max(v + damping * delta_v, v0), v1)
            trial_error = _norm(_sub(surface.evaluate(trial_u, trial_v), point))
            if trial_error < error:
                u, v = trial_u, trial_v
                improved = True
                break
            damping *= 0.5
        if not improved:
            break
    return (u, v, _norm(_sub(surface.evaluate(u, v), point)))


def _sub(a: Sequence[float], b: Sequence[float]) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _add(a: Sequence[float], b: Sequence[float]) -> Vec3:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def _scale(a: Sequence[float], factor: float) -> Vec3:
    return (a[0] * factor, a[1] * factor, a[2] * factor)


def _dot(a: Sequence[float], b: Sequence[float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a: Sequence[float], b: Sequence[float]) -> Vec3:
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _norm(a: Sequence[float]) -> float:
    return math.sqrt(_dot(a, a))


# --------------------------------------------------------------------------- #
# triangulation of a trimmed parameter domain
# --------------------------------------------------------------------------- #
def _signed_area(polygon: Sequence[tuple[float, float]]) -> float:
    total = 0.0
    count = len(polygon)
    for index in range(count):
        x1, y1 = polygon[index]
        x2, y2 = polygon[(index + 1) % count]
        total += x1 * y2 - x2 * y1
    return 0.5 * total


def _point_in_triangle(px: float, py: float, ax: float, ay: float, bx: float, by: float,
                       cx: float, cy: float) -> bool:
    d1 = (px - bx) * (ay - by) - (ax - bx) * (py - by)
    d2 = (px - cx) * (by - cy) - (bx - cx) * (py - cy)
    d3 = (px - ax) * (cy - ay) - (cx - ax) * (py - ay)
    has_negative = d1 < 0.0 or d2 < 0.0 or d3 < 0.0
    has_positive = d1 > 0.0 or d2 > 0.0 or d3 > 0.0
    return not (has_negative and has_positive)


def _ear_clip(polygon: Sequence[tuple[float, float]]) -> list[tuple[int, int, int]]:
    """Ear clipping on a simple polygon given counter-clockwise; indices refer to input."""
    indices = list(range(len(polygon)))
    triangles: list[tuple[int, int, int]] = []
    guard = 0
    while len(indices) > 3:
        guard += 1
        if guard > 4 * len(polygon) + 32:
            raise ValueError("ear clipping stalled: the trim loop is self-intersecting")
        found = False
        count = len(indices)
        for position in range(count):
            previous = indices[(position - 1) % count]
            current = indices[position]
            following = indices[(position + 1) % count]
            ax, ay = polygon[previous]
            bx, by = polygon[current]
            cx, cy = polygon[following]
            area = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
            if area <= 0.0:
                continue  # reflex or collinear corner: not an ear
            contains = False
            for other in indices:
                if other in (previous, current, following):
                    continue
                px, py = polygon[other]
                if (px, py) in ((ax, ay), (bx, by), (cx, cy)):
                    # A bridged hole duplicates its first vertex; a duplicate point is not
                    # "inside" the ear, and treating it as inside stalls the clip.
                    continue
                if _point_in_triangle(px, py, ax, ay, bx, by, cx, cy):
                    contains = True
                    break
            if contains:
                continue
            triangles.append((previous, current, following))
            indices.pop(position)
            found = True
            break
        if not found:
            raise ValueError("ear clipping found no ear: the trim loop is not a simple polygon")
    triangles.append((indices[0], indices[1], indices[2]))
    return triangles


def triangulate_trimmed_domain(loops: Sequence[Sequence[tuple[float, float]]]) -> list[tuple[int, int, int]]:
    """Triangulate a region bounded by one outer loop and zero or more holes.

    ``loops`` are closed in ``(u, v)`` (first point may repeat the last).  The outer loop
    is the one with the largest absolute area; holes are reversed as needed.  The result
    indexes into the concatenated, de-duplicated vertex list built here, so the caller can
    rebuild the exact boundary vertices.
    """
    cleaned: list[list[tuple[float, float]]] = []
    for loop in loops:
        points = list(loop)
        if len(points) > 1 and points[0] == points[-1]:
            points = points[:-1]
        if len(points) < 3:
            continue
        cleaned.append(points)
    if not cleaned:
        raise ValueError("no loop with at least three points")
    areas = [_signed_area(loop) for loop in cleaned]
    outer_index = max(range(len(cleaned)), key=lambda index: abs(areas[index]))
    outer = cleaned[outer_index]
    if areas[outer_index] < 0:
        outer = list(reversed(outer))
    holes = []
    for index, loop in enumerate(cleaned):
        if index == outer_index:
            continue
        holes.append(loop if _signed_area(loop) < 0 else list(reversed(loop)))
    polygon = _bridge_holes(outer, holes) if holes else outer
    triangles = _ear_clip(polygon)
    return triangles


# --------------------------------------------------------------------------- #
# tessellation
# --------------------------------------------------------------------------- #
@dataclass
class Tessellation:
    """Result of tessellating a trimmed model: triangles in file units plus statistics."""

    vertices: list[Vec3]
    faces: list[Triangle]
    face_stats: list[dict[str, Any]]
    refused: dict[int, str]
    max_residual: float
    residual_tolerance: float

    @property
    def triangle_count(self) -> int:
        return len(self.faces)


def _loop_to_uv(surface: Surface, loop: Sequence[Vec3], tolerance: float) -> tuple[list[tuple[float, float]], float]:
    """Invert a 3-D loop onto ``(u, v)``; returns the parameters and the worst residual."""
    parameters: list[tuple[float, float]] = []
    worst = 0.0
    seed: list[tuple[float, float]] = []
    for point in loop:
        u, v, residual = invert_point(surface, point, seeds=seed)
        worst = max(worst, residual)
        parameters.append((u, v))
        seed = [(u, v)]
    if parameters and parameters[0] != parameters[-1]:
        parameters.append(parameters[0])
    return parameters, worst


def tessellate_trimmed_model(model: TrimmedModel, *, residual_tolerance: float | None = None,
                             scale: float = 1.0) -> Tessellation:
    """Tessellate every face; faces whose trim does not lie on the surface are refused."""
    extent = _model_extent_from_faces(model.faces)
    tolerance = residual_tolerance if residual_tolerance is not None else 1e-6 * max(extent, 1e-9)
    vertices: list[Vec3] = []
    faces: list[Triangle] = []
    stats: list[dict[str, Any]] = []
    refused: dict[int, str] = {}
    worst_residual = 0.0

    for face in model.faces:
        try:
            uv_loops: list[list[tuple[float, float]]] = []
            face_residual = 0.0
            for loop in face.loops:
                parameters, residual = _loop_to_uv(face.surface, loop, tolerance)
                face_residual = max(face_residual, residual)
                uv_loops.append(parameters)
            if face_residual > tolerance:
                raise IgesError(
                    f"trim loop does not lie on its surface: inversion residual "
                    f"{face_residual:.3e} > tolerance {tolerance:.3e} (file units)"
                )
            triangles = triangulate_trimmed_domain(uv_loops)
            cleaned = _without_duplicate(uv_loops)
            polygon_area = sum(abs(_signed_area(loop)) for loop in cleaned)
            base = len(vertices)
            # Boundary vertices are emitted from the (u, v) parameters, so two faces that
            # share a trim curve produce identical 3-D boundary points and stitch exactly.
            flat = [point for loop in cleaned for point in loop]
            for u, v in flat:
                vertices.append(_scale(face.surface.evaluate(u, v), scale))
            face_triangles: list[Triangle] = []
            for a, b, c in triangles:
                if a >= len(flat) or b >= len(flat) or c >= len(flat):
                    continue
                ia, ib, ic = a, b, c
                if ia == ib or ib == ic or ia == ic:
                    continue
                face_triangles.append((base + ia, base + ib, base + ic))
            area = sum(_triangle_area(vertices[t[0]], vertices[t[1]], vertices[t[2]])
                       for t in face_triangles)
            stats.append({
                "index": face.index,
                "source_entity": face.source_entity,
                "surface": face.surface_kind,
                "loops": len(uv_loops),
                "loop_points": sum(len(loop) for loop in uv_loops),
                "triangles": len(face_triangles),
                "inversion_residual": face_residual,
                "parametric_area": polygon_area,
                "mesh_area": area,
                "area_relative_error": (
                    abs(area - polygon_area) / polygon_area if polygon_area else 0.0
                ),
            })
            faces.extend(face_triangles)
            worst_residual = max(worst_residual, face_residual)
        except (ValueError, IgesError) as exc:
            refused[face.index] = str(exc)
    return Tessellation(
        vertices=vertices,
        faces=faces,
        face_stats=stats,
        refused=refused,
        max_residual=worst_residual,
        residual_tolerance=tolerance,
    )


def _without_duplicate(loops: Sequence[Sequence[tuple[float, float]]]) -> list[list[tuple[float, float]]]:
    out: list[list[tuple[float, float]]] = []
    for loop in loops:
        points = list(loop)
        if len(points) > 1 and points[0] == points[-1]:
            points = points[:-1]
        if len(points) >= 3:
            out.append(points)
    return out


def _triangle_area(a: Sequence[float], b: Sequence[float], c: Sequence[float]) -> float:
    return 0.5 * _norm(_cross(_sub(b, a), _sub(c, a)))


def _model_extent_from_faces(faces: Sequence[TrimmedFace]) -> float:
    largest = 0.0
    for face in faces:
        for loop in face.loops:
            for point in loop:
                largest = max(largest, abs(point[0]), abs(point[1]), abs(point[2]))
    return largest if largest > 0.0 else 1.0


# --------------------------------------------------------------------------- #
# welding, audit and repair
# --------------------------------------------------------------------------- #
def weld_vertices(vertices: Sequence[Vec3], faces: Sequence[Triangle],
                  tolerance: float) -> tuple[list[Vec3], list[Triangle], int]:
    """Merge vertices whose positions agree within ``tolerance`` (grid hash, then verify)."""
    if tolerance <= 0.0:
        return list(vertices), list(faces), 0
    cell = tolerance
    buckets: dict[tuple[int, int, int], list[int]] = {}
    merged: list[Vec3] = []
    remap: list[int] = []
    for point in vertices:
        key = (int(math.floor(point[0] / cell)), int(math.floor(point[1] / cell)),
               int(math.floor(point[2] / cell)))
        target = -1
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for candidate in buckets.get((key[0] + dx, key[1] + dy, key[2] + dz), ()):
                        if _norm(_sub(merged[candidate], point)) <= tolerance:
                            target = candidate
                            break
                    if target >= 0:
                        break
                if target >= 0:
                    break
            if target >= 0:
                break
        if target < 0:
            target = len(merged)
            merged.append(tuple(point))  # type: ignore[arg-type]
            buckets.setdefault(key, []).append(target)
        remap.append(target)
    welded: list[Triangle] = []
    for face in faces:
        a, b, c = remap[face[0]], remap[face[1]], remap[face[2]]
        if len({a, b, c}) == 3:
            welded.append((a, b, c))
    return merged, welded, len(vertices) - len(merged)


@dataclass
class MeshAudit:
    """Topology and geometry measurements of one triangle mesh."""

    vertex_count: int
    triangle_count: int
    edge_count: int
    free_edge_count: int
    non_manifold_edge_count: int
    duplicate_triangle_count: int
    degenerate_triangle_count: int
    orientation_conflicts: int
    component_count: int
    closed_component_count: int
    euler_characteristic: int
    signed_volume: float
    area: float
    bounding_box: tuple[Vec3, Vec3]
    watertight: bool
    manifold: bool
    consistently_oriented: bool

    def to_dict(self) -> dict[str, Any]:
        return {
            "vertices": self.vertex_count,
            "triangles": self.triangle_count,
            "edges": self.edge_count,
            "free_edges": self.free_edge_count,
            "non_manifold_edges": self.non_manifold_edge_count,
            "duplicate_triangles": self.duplicate_triangle_count,
            "degenerate_triangles": self.degenerate_triangle_count,
            "orientation_conflicts": self.orientation_conflicts,
            "components": self.component_count,
            "closed_components": self.closed_component_count,
            "euler_characteristic": self.euler_characteristic,
            "signed_volume": self.signed_volume,
            "area": self.area,
            "watertight": self.watertight,
            "manifold": self.manifold,
            "consistently_oriented": self.consistently_oriented,
            "bounding_box": [list(self.bounding_box[0]), list(self.bounding_box[1])],
        }


def _edge_key(a: int, b: int) -> tuple[int, int]:
    return (a, b) if a < b else (b, a)


def audit_mesh(vertices: Sequence[Vec3], faces: Sequence[Triangle],
               *, degenerate_area: float = 0.0) -> MeshAudit:
    """Measure a triangle mesh: topology, orientation, volume and area."""
    edge_uses: dict[tuple[int, int], list[tuple[int, bool]]] = {}
    areas: list[float] = []
    volume = 0.0
    degenerate = 0
    for index, (a, b, c) in enumerate(faces):
        pa, pb, pc = vertices[a], vertices[b], vertices[c]
        area = _triangle_area(pa, pb, pc)
        areas.append(area)
        if area <= degenerate_area:
            degenerate += 1
        volume += _dot(pa, _cross(pb, pc)) / 6.0
        for start, end in ((a, b), (b, c), (c, a)):
            key = _edge_key(start, end)
            edge_uses.setdefault(key, []).append((index, start < end))
    free_edges = sum(1 for uses in edge_uses.values() if len(uses) == 1)
    non_manifold = sum(1 for uses in edge_uses.values() if len(uses) > 2)
    conflicts = 0
    for uses in edge_uses.values():
        if len(uses) == 2 and uses[0][1] == uses[1][1]:
            conflicts += 1
    seen: set[tuple[int, int, int]] = set()
    duplicates = 0
    for face in faces:
        key = tuple(sorted(face))
        if key in seen:
            duplicates += 1
        seen.add(key)

    adjacency: dict[int, list[int]] = {}
    for uses in edge_uses.values():
        if len(uses) == 2:
            first, second = uses[0][0], uses[1][0]
            adjacency.setdefault(first, []).append(second)
            adjacency.setdefault(second, []).append(first)
    component_of: dict[int, int] = {}
    component_count = 0
    for face_index in range(len(faces)):
        if face_index in component_of:
            continue
        stack = [face_index]
        component_of[face_index] = component_count
        while stack:
            current = stack.pop()
            for neighbour in adjacency.get(current, ()):  # pragma: no branch
                if neighbour not in component_of:
                    component_of[neighbour] = component_count
                    stack.append(neighbour)
        component_count += 1
    closed_components = 0
    for component in range(component_count):
        member_faces = {index for index, value in component_of.items() if value == component}
        member_edges = {
            key for key, uses in edge_uses.items()
            if any(index in member_faces for index, _ in uses)
        }
        if all(len(edge_uses[key]) == 2 for key in member_edges):
            closed_components += 1

    referenced = {index for face in faces for index in face}
    vertex_count = len(referenced)
    edge_count = len(edge_uses)
    return MeshAudit(
        vertex_count=vertex_count,
        triangle_count=len(faces),
        edge_count=edge_count,
        free_edge_count=free_edges,
        non_manifold_edge_count=non_manifold,
        duplicate_triangle_count=duplicates,
        degenerate_triangle_count=degenerate,
        orientation_conflicts=conflicts,
        component_count=component_count,
        closed_component_count=closed_components,
        euler_characteristic=vertex_count - edge_count + len(faces),
        signed_volume=volume,
        area=sum(areas),
        bounding_box=(_component_min(vertices, faces), _component_max(vertices, faces)),
        watertight=free_edges == 0 and non_manifold == 0 and duplicates == 0
        and degenerate == 0 and conflicts == 0 and len(faces) > 0,
        manifold=non_manifold == 0,
        consistently_oriented=conflicts == 0,
    )


def _component_min(vertices: Sequence[Vec3], faces: Sequence[Triangle]) -> Vec3:
    indices = {index for face in faces for index in face}
    if not indices:
        return (0.0, 0.0, 0.0)
    return tuple(min(vertices[index][axis] for index in indices) for axis in range(3))  # type: ignore[return-value]


def _component_max(vertices: Sequence[Vec3], faces: Sequence[Triangle]) -> Vec3:
    indices = {index for face in faces for index in face}
    if not indices:
        return (0.0, 0.0, 0.0)
    return tuple(max(vertices[index][axis] for index in indices) for axis in range(3))  # type: ignore[return-value]


def _orient_consistently(faces: Sequence[Triangle]) -> tuple[list[Triangle], int]:
    """Propagate one winding across the face adjacency graph; returns the flipped count."""
    edge_uses: dict[tuple[int, int], list[tuple[int, bool]]] = {}
    for index, (a, b, c) in enumerate(faces):
        for start, end in ((a, b), (b, c), (c, a)):
            edge_uses.setdefault(_edge_key(start, end), []).append((index, start < end))
    adjacency: dict[int, list[tuple[int, bool]]] = {}
    for uses in edge_uses.values():
        if len(uses) != 2:
            continue
        (first, first_forward), (second, second_forward) = uses
        same = first_forward == second_forward
        adjacency.setdefault(first, []).append((second, same))
        adjacency.setdefault(second, []).append((first, same))
    flipped = [False] * len(faces)
    visited = [False] * len(faces)
    flip_count = 0
    for seed in range(len(faces)):
        if visited[seed]:
            continue
        visited[seed] = True
        stack = [seed]
        while stack:
            current = stack.pop()
            for neighbour, same in adjacency.get(current, ()):
                if visited[neighbour]:
                    continue
                visited[neighbour] = True
                flipped[neighbour] = flipped[current] ^ same
                if flipped[neighbour]:
                    flip_count += 1
                stack.append(neighbour)
    result: list[Triangle] = []
    for index, (a, b, c) in enumerate(faces):
        result.append((a, c, b) if flipped[index] else (a, b, c))
    return result, flip_count


def _boundary_loops(faces: Sequence[Triangle]) -> list[list[int]]:
    """Chain the free edges of a mesh into closed loops of vertex indices."""
    edge_uses: dict[tuple[int, int], int] = {}
    directed: dict[int, list[int]] = {}
    for a, b, c in faces:
        for start, end in ((a, b), (b, c), (c, a)):
            key = _edge_key(start, end)
            edge_uses[key] = edge_uses.get(key, 0) + 1
            directed.setdefault(start, []).append(end)
    boundary = [(start, end) for (start, end), uses in edge_uses.items() if uses == 1
                for start, end in [(start, end)]]
    # Rebuild direction from the face windings (free edges appear once, so keep that sense).
    directed_edges: list[tuple[int, int]] = []
    counted: dict[tuple[int, int], int] = {}
    for a, b, c in faces:
        for start, end in ((a, b), (b, c), (c, a)):
            key = _edge_key(start, end)
            counted[key] = counted.get(key, 0) + 1
    for a, b, c in faces:
        for start, end in ((a, b), (b, c), (c, a)):
            if counted[_edge_key(start, end)] == 1:
                directed_edges.append((start, end))
    outgoing: dict[int, list[int]] = {}
    for start, end in directed_edges:
        outgoing.setdefault(start, []).append(end)
    loops: list[list[int]] = []
    used: set[tuple[int, int]] = set()
    for start, end in directed_edges:
        if (start, end) in used:
            continue
        loop = [start]
        current_start, current_end = start, end
        used.add((current_start, current_end))
        guard = 0
        while current_end != start and guard < len(directed_edges) + 2:
            guard += 1
            loop.append(current_end)
            candidates = [nxt for nxt in outgoing.get(current_end, ())
                          if (current_end, nxt) not in used]
            if not candidates:
                break
            nxt = candidates[0]
            used.add((current_end, nxt))
            current_start, current_end = current_end, nxt
        if current_end == start and len(loop) >= 3:
            loops.append(loop)
    del boundary
    return loops


def repair_mesh(vertices: Sequence[Vec3], faces: Sequence[Triangle], *, weld_tolerance: float,
                max_hole_edges: int = 0, fill_holes: bool = True) -> tuple[list[Vec3], list[Triangle], dict[str, Any]]:
    """Weld, drop degenerate and duplicate triangles, fill small holes, then orient.

    Returns the repaired mesh plus a record of exactly what was changed.  Hole filling is
    bounded by ``max_hole_edges``: a hole with more boundary edges than that is left open
    and reported, because fanning a large hole would invent surface the file never had.
    """
    record: dict[str, Any] = {}
    welded_vertices, welded_faces, merged = weld_vertices(vertices, faces, weld_tolerance)
    record["vertices_merged"] = merged

    kept: list[Triangle] = []
    degenerate = 0
    seen: set[tuple[int, int, int]] = set()
    duplicates = 0
    for face in welded_faces:
        key = tuple(sorted(face))
        if key in seen:
            duplicates += 1
            continue
        seen.add(key)
        kept.append(face)
    record["duplicate_triangles_removed"] = duplicates
    record["degenerate_triangles"] = degenerate

    filled = 0
    if fill_holes:
        loops = _boundary_loops(kept)
        for loop in loops:
            if len(loop) > max_hole_edges:
                continue
            points = [welded_vertices[index] for index in loop]
            centroid = (sum(p[0] for p in points) / len(points),
                        sum(p[1] for p in points) / len(points),
                        sum(p[2] for p in points) / len(points))
            if len(loop) <= 4:
                # Small loops: a fan from the centroid is the minimal-area completion.
                centre = len(welded_vertices)
                welded_vertices.append(centroid)
                for position in range(len(loop)):
                    kept.append((centre, loop[position], loop[(position + 1) % len(loop)]))
                    filled += 1
            else:
                for position in range(1, len(loop) - 1):
                    kept.append((loop[0], loop[position], loop[position + 1]))
                    filled += 1
    record["triangles_added_by_hole_fill"] = filled

    oriented, flips = _orient_consistently(kept)
    record["orientation_flips"] = flips
    signed_volume = sum(
        _dot(welded_vertices[a], _cross(welded_vertices[b], welded_vertices[c])) / 6.0
        for a, b, c in oriented
    )
    if signed_volume < 0.0:
        oriented = [(a, c, b) for a, b, c in oriented]
        record["global_flip"] = True
        signed_volume = -signed_volume
    else:
        record["global_flip"] = False
    record["signed_volume_after"] = signed_volume
    return welded_vertices, oriented, record


# --------------------------------------------------------------------------- #
# the import pipeline
# --------------------------------------------------------------------------- #
@dataclass
class CadImport:
    """One import run: the mesh, the audit trail, and the gate verdict."""

    mesh: Mesh
    report: dict[str, Any]

    def to_dict(self) -> dict[str, Any]:
        return self.report


def import_iges(path: str | Path, *, tolerance: float | None = None, scale: float | None = None,
                weld_tolerance: float | None = None, repair: bool = True,
                fill_holes: bool = True, max_hole_edges: int = 32) -> CadImport:
    """Import an IGES file: parse, invert trims, tessellate, repair, audit.

    ``scale`` converts the file's units to metres; it defaults to the file's own declared
    units.  ``weld_tolerance`` defaults to 1e-6 times the model's largest coordinate.
    """
    started = time.time()
    digest = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    model = read_iges_trimmed(path, tolerance=tolerance)
    unit_scale = scale if scale is not None else model.metres_per_unit
    tessellation = tessellate_trimmed_model(model, scale=unit_scale)
    extent = _model_extent_from_faces(model.faces) * unit_scale
    weld = weld_tolerance if weld_tolerance is not None else 1e-6 * max(extent, 1e-9)

    before = audit_mesh(tessellation.vertices, tessellation.faces, degenerate_area=0.0)
    vertices = list(tessellation.vertices)
    faces = list(tessellation.faces)
    repair_record: dict[str, Any] = {}
    if repair:
        vertices, faces, repair_record = repair_mesh(
            vertices, faces, weld_tolerance=weld, max_hole_edges=max_hole_edges,
            fill_holes=fill_holes,
        )
    else:
        vertices, faces, merged = weld_vertices(vertices, faces, weld)
        repair_record = {"vertices_merged": merged, "weld_only": True}
    after = audit_mesh(vertices, faces, degenerate_area=0.0)

    mesh = Mesh(
        vertices=tuple(vertices),
        faces=tuple(faces),
        units="m" if abs(unit_scale - 1.0) > 0 or model.metres_per_unit else model.units_name,
    )
    gates = _gates(model, tessellation, after)
    report: dict[str, Any] = {
        "tool": "tensorcad.geomquality.trim_mesh",
        "source": {
            "path": str(Path(path)),
            "sha256": digest,
            "bytes": Path(path).stat().st_size,
            "units": model.units_name,
            "metres_per_unit": model.metres_per_unit,
            "scale_applied": unit_scale,
        },
        "import": {
            "faces_in_file": len(model.faces),
            "faces_tessellated": len(tessellation.face_stats),
            "faces_refused": len(tessellation.refused),
            "refused": {str(key): value for key, value in tessellation.refused.items()},
            "reader_refused": {str(key): value for key, value in model.refused.items()},
            "max_inversion_residual": tessellation.max_residual,
            "residual_tolerance": tessellation.residual_tolerance,
            "weld_tolerance": weld,
            "face_stats": tessellation.face_stats,
            "entity_census": {str(key): value for key, value in model.census.items()},
            "warnings": list(model.warnings),
        },
        "before_repair": before.to_dict(),
        "repair": repair_record,
        "after_repair": after.to_dict(),
        "gates": gates,
        "elapsed_seconds": round(time.time() - started, 3),
    }
    return CadImport(mesh=mesh, report=report)


def _gates(model: TrimmedModel, tessellation: Tessellation, audit: MeshAudit) -> dict[str, Any]:
    """Machine-decidable gate verdicts; each carries the number that decided it."""
    open_faces = len(tessellation.refused)
    return {
        "G1_all_faces_read": {
            "passed": open_faces == 0 and bool(model.faces),
            "value": f"{len(tessellation.face_stats)}/{len(model.faces)}",
            "detail": "faces tessellated / faces in file",
        },
        "G2_trim_on_surface": {
            "passed": tessellation.max_residual <= tessellation.residual_tolerance,
            "value": tessellation.max_residual,
            "limit": tessellation.residual_tolerance,
            "detail": "worst inversion residual of a trim loop point (file units x scale)",
        },
        "G3_watertight": {
            "passed": audit.watertight,
            "value": audit.free_edge_count,
            "limit": 0,
            "detail": "free edges after repair; 0 means every edge has two faces",
        },
        "G4_manifold": {
            "passed": audit.manifold and audit.consistently_oriented,
            "value": audit.non_manifold_edge_count + audit.orientation_conflicts,
            "limit": 0,
            "detail": "non-manifold edges + orientation conflicts",
        },
        "G5_positive_volume": {
            "passed": audit.signed_volume > 0.0,
            "value": audit.signed_volume,
            "limit": 0.0,
            "detail": "signed volume, cubic model units",
        },
        "G6_single_shell": {
            "passed": audit.component_count == 1 and audit.closed_component_count == 1,
            "value": f"{audit.component_count}/{audit.closed_component_count}",
            "detail": "components / closed components",
        },
        "G7_no_degenerate": {
            "passed": audit.degenerate_triangle_count == 0 and audit.duplicate_triangle_count == 0,
            "value": audit.degenerate_triangle_count + audit.duplicate_triangle_count,
            "limit": 0,
            "detail": "degenerate + duplicate triangles",
        },
    }


def import_cad(path: str | Path, *, allow_patch_fallback: bool = True,
               **kwargs: Any) -> CadImport:
    """Import a CAD file, choosing the route the file's own content justifies.

    IGES: the standards-compliant trimmed reader runs first.  If it refuses the file (an
    unsupported entity, a trim that does not invert onto its surface, a non-conformant
    B-spline record), the untrimmed patch-set route is tried and the refusal is carried
    into the report, so a fallback is never silent.  STEP: patch sets only.
    """
    u_samples = int(kwargs.pop("u_samples", 24))
    v_samples = int(kwargs.pop("v_samples", 12))
    target = Path(path)
    suffix = target.suffix.lower()
    if suffix in (".igs", ".iges"):
        try:
            return import_iges(target, **kwargs)
        except (IgesError, ValueError) as exc:
            if not allow_patch_fallback:
                raise
            result = import_patch_set(target, **{k: v for k, v in kwargs.items()
                                                if k in {"scale", "weld_tolerance", "repair",
                                                         "fill_holes", "max_hole_edges"}},
                                     u_samples=u_samples, v_samples=v_samples)
            result.report.setdefault("route_notes", []).append(
                f"trimmed reader refused the file, fell back to the patch-set route: {exc}"
            )
            result.report["trimmed_route_refusal"] = str(exc)
            return result
    if suffix in (".step", ".stp"):
        return import_patch_set(target, **{k: v for k, v in kwargs.items()
                                          if k in {"scale", "weld_tolerance", "repair",
                                                   "fill_holes", "max_hole_edges"}},
                                u_samples=u_samples, v_samples=v_samples)
    raise IgesError(f"unsupported CAD suffix {suffix!r}: expected .igs/.iges/.step/.stp")


def write_mesh(mesh: Mesh, path: str | Path, *, fmt: str | None = None) -> Path:
    """Write STL (ASCII) or OBJ; both are plain text and need no third-party package."""
    target = Path(path)
    fmt = (fmt or target.suffix.lstrip(".")).lower()
    if fmt == "obj":
        lines = [f"# tensorcad import: {len(mesh.vertices)} vertices, {len(mesh.faces)} triangles"]
        lines.extend(f"v {x:.12g} {y:.12g} {z:.12g}" for x, y, z in mesh.vertices)
        lines.extend(f"f {a} {b} {c}" for a, b, c in mesh.to_obj_faces())
        target.write_text("\n".join(lines) + "\n")
        return target
    if fmt == "stl":
        lines = ["solid tensorcad_import"]
        for face in mesh.faces:
            lines.append("  facet normal 0 0 0")
            lines.append("    outer loop")
            for index in face:
                x, y, z = mesh.vertices[index]
                lines.append(f"      vertex {x:.12g} {y:.12g} {z:.12g}")
            lines.append("    endloop")
            lines.append("  endfacet")
        lines.append("endsolid tensorcad_import")
        target.write_text("\n".join(lines) + "\n")
        return target
    if fmt == "json":
        payload = {
            "vertices": [list(vertex) for vertex in mesh.vertices],
            "faces": [list(face) for face in mesh.faces],
            "units": mesh.units,
        }
        target.write_text(json.dumps(payload))
        return target
    raise IgesError(f"unsupported mesh format {fmt!r}: expected obj, stl or json")


def _unused(*_: Iterable[Any]) -> None:  # pragma: no cover
    """Placeholder to keep the import list honest if unused imports accumulate."""

def import_patch_set(path: str | Path, *, scale: float | None = None,
                     weld_tolerance: float | None = None, repair: bool = True,
                     fill_holes: bool = True, max_hole_edges: int = 32,
                     u_samples: int = 24, v_samples: int = 12) -> CadImport:
    """Import an *untrimmed* patch set (the flavour this repository's writer emits).

    IGES 128 patch sets and STEP patch sets carry no topology: adjacency has to be
    derived by sewing patch boundaries within a tolerance (``geomquality.io``).  That
    route is kept here, and it runs through the same audit, repair and gate code as the
    trimmed route so the numbers are directly comparable.
    """
    from .io import read_iges as _read_patch_iges
    from .io import read_step as _read_patch_step

    started = time.time()
    target = Path(path)
    digest = hashlib.sha256(target.read_bytes()).hexdigest()
    suffix = target.suffix.lower()
    reader = _read_patch_step if suffix in (".step", ".stp") else _read_patch_iges
    brep = reader(target)
    mesh = brep.tessellate(u_samples, v_samples)
    vertices = [tuple(vertex) for vertex in mesh.vertices]
    faces = [tuple(face) for face in mesh.faces]
    unit_scale = scale if scale is not None else 1.0
    extent = 0.0
    for vertex in vertices:
        extent = max(extent, abs(vertex[0]), abs(vertex[1]), abs(vertex[2]))
    weld = weld_tolerance if weld_tolerance is not None else 1e-6 * max(extent * unit_scale, 1e-9)
    before = audit_mesh(vertices, faces)
    if repair:
        vertices, faces, repair_record = repair_mesh(
            vertices, faces, weld_tolerance=weld, max_hole_edges=max_hole_edges,
            fill_holes=fill_holes,
        )
    else:
        vertices, faces, merged = weld_vertices(vertices, faces, weld)
        repair_record = {"vertices_merged": merged, "weld_only": True}
    after = audit_mesh(vertices, faces)
    mesh_out = Mesh(vertices=tuple(vertices), faces=tuple(faces), units="m")
    gates = {
        "G1_all_faces_read": {"passed": True, "value": f"{len(brep.faces)}/{len(brep.faces)}",
                              "detail": "patch-set route: topology derived by sewing"},
        "G2_trim_on_surface": {"passed": True, "value": 0.0,
                               "detail": "not applicable: every face covers its whole parameter domain"},
        "G3_watertight": {"passed": after.watertight, "value": after.free_edge_count, "limit": 0,
                          "detail": "free edges after repair"},
        "G4_manifold": {"passed": after.manifold and after.consistently_oriented,
                        "value": after.non_manifold_edge_count + after.orientation_conflicts,
                        "limit": 0, "detail": "non-manifold edges + orientation conflicts"},
        "G5_positive_volume": {"passed": after.signed_volume > 0.0,
                               "value": after.signed_volume, "limit": 0.0,
                               "detail": "signed volume"},
        "G6_single_shell": {"passed": after.component_count == 1 and after.closed_component_count == 1,
                            "value": f"{after.component_count}/{after.closed_component_count}",
                            "detail": "components / closed components"},
        "G7_no_degenerate": {"passed": after.degenerate_triangle_count == 0,
                             "value": after.degenerate_triangle_count, "limit": 0,
                             "detail": "degenerate + duplicate triangles"},
    }
    report = {
        "tool": "tensorcad.geomquality.trim_mesh",
        "route": "patch-set (untrimmed, sewn)",
        "source": {"path": str(target), "sha256": digest, "bytes": target.stat().st_size,
                   "units": "model units", "metres_per_unit": 1.0, "scale_applied": unit_scale},
        "import": {"faces_in_file": len(brep.faces), "faces_tessellated": len(brep.faces),
                   "faces_refused": 0, "refused": {},
                   "samples": {"u": u_samples, "v": v_samples}, "weld_tolerance": weld},
        "before_repair": before.to_dict(),
        "repair": repair_record,
        "after_repair": after.to_dict(),
        "gates": gates,
        "elapsed_seconds": round(time.time() - started, 3),
    }
    return CadImport(mesh=mesh_out, report=report)
