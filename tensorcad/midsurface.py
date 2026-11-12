"""Mid-surface extraction and defeaturing: a thin-walled solid becomes a shell sheet.

``GAP-CAD-013`` in ``docs/oss_inventory_cad/GAP-ANCHORS.md`` records that this package can
audit a surface (:mod:`tensorcad.geometry_audit`), hollow a solid (:mod:`tensorcad.shell`) and
repair a mesh (:mod:`tensorcad.mesh_repair`), but cannot turn a thin-walled *solid* into the
mid-surface sheet a shell finite-element model needs. Its acceptance criterion is "extract the
mid-surface of a plate of thickness ``t``: area error < 1 %, face count <= 1/10 of the source,
and feature suppression cuts the mesh size by >= 50 %".

**Why a ship needs this.** A double bottom, a bulkhead, a wing tank or a stiffened panel
arrives as a solid whose wall is a few millimetres thick inside a hull a few hundred metres
long. A solid-element model of that is unrunnable (element aspect ratios of 1e4); the
structural model is the mid-surface sheet plus one thickness per patch. Nothing in the package
produced that sheet.

**The method -- three steps, each one an independently checkable number.**

1. *Measure the wall, do not assume it.* For every triangle a ray is cast from its centroid
   along the inward normal (``-n``) through a uniform triangle grid; the first hit is the
   opposed wall if its normal is anti-parallel (``-dot(n_src, n_hit) >= opposition``) and the
   distance is within ``max_thickness``. The pair is kept only when it is **mutual** (the
   reverse ray returns the original face), which is what separates a real wall from a chance
   anti-parallel encounter between unrelated features. The result is a thickness map: paired
   thickness and partner per face, plus a named reason for every unpaired face.
2. *Group the faces into panels.* Faces are linked when they share a coincident vertex
   (positions are matched within ``weld_tolerance``) *and* their planes agree within
   ``angle_tolerance`` / ``offset_tolerance``; panels are the connected components. This is
   what removes spurious CAD seams: a face split into sixteen imported patches is one panel if
   the patches really do lie in one plane.
3. *Place the sheet.* Each pair of faces becomes **one** panel -- each unordered pair is used
   once, from its smaller-indexed face, so a hollow body does not produce two coincident
   sheets. The panel plane is offset inward by ``t/2``; its boundary starts as the convex hull
   of the projected panel vertices and is then clipped by the mid-planes of the *adjacent
   panels*. The clipping is what makes a cornered wall behave: a box shelled at ``t`` gives six
   patches meeting on the lines ``x = y = t/2``, i.e. the closed sheet of the eroded box, whose
   area is exactly ``6 (L - t)^2`` and whose enclosed volume is exactly ``(L - t)^3``. A
   neighbour that is not a panel -- the edge strip of a plate, whose opposed face is ``L`` away
   rather than ``t`` -- clips nothing, which is also the right answer: a plate's sheet runs out
   to the plate edge.

**Defeaturing.** :func:`defeature` re-runs steps 2-3 with relaxed tolerances and an optional
area floor: a speck of a few facets below ``min_panel_area`` disappears, and patches that a CAD
import left a hair apart are bridged within ``weld_tolerance``. Both operations are *measured*,
not asserted: :meth:`MidSurface.face_reduction` is the fraction of the source triangles the
sheet no longer needs, and :meth:`MidSurface.removed_area` is the area that dropping features
cost.

**What it refuses** (a wrongly extracted sheet still looks like a sheet, so the refusals are
part of the contract):

* an open or inward-oriented surface (run :func:`tensorcad.mesh.orient_outward` first);
* a body with no opposed thin walls at all -- a solid cube's faces meet faces whose normals
  point the same way, so nothing pairs, and the error names the qualified area fraction;
* a wall-thickness spread inside one panel larger than ``thickness_tolerance``: variable
  thickness is not a constant-offset sheet and is reported instead of averaged;
* a panel that is not convex once flattened: the hull would silently over-cover it, so the call
  fails with the measured overshoot.

**What it does not do.** No B-rep rim, no face removal and no stitching back into a solid
(``GAP-CAD-001``); the sheet of a closed thin-walled body is closed, the sheet of an open region
is a patch, and both are reported through :meth:`MidSurface.watertight`. Curved walls are
represented by their facets, so the mid-sheet of a faceted tube is exact for the *faceted*
solid and short of the analytic tube by ``O(h^2)``; the tests measure both. Thickness is
constant within a panel by construction, so a taper is refused rather than shelled as a wedge.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

from . import mesh as meshlib
from .geometry import Face, Mesh, Vec3

__all__ = [
    "MidSurface",
    "MidSurfaceError",
    "Panel",
    "ThicknessMap",
    "defeature",
    "extract_midsurface",
    "measure_thickness",
    "suggest_max_thickness",
]


class MidSurfaceError(ValueError):
    """Raised when no defensible mid-surface follows from the input."""


# --------------------------------------------------------------------------------------
# vector helpers (the core is dependency-free: no NumPy anywhere in this module)
# --------------------------------------------------------------------------------------


def _sub(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _add(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def _scale(a: Vec3, s: float) -> Vec3:
    return (a[0] * s, a[1] * s, a[2] * s)


def _dot(a: Vec3, b: Vec3) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a: Vec3, b: Vec3) -> Vec3:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _length(a: Vec3) -> float:
    return math.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])


def _unit(a: Vec3) -> Vec3:
    norm = _length(a)
    if norm <= 0.0:
        raise MidSurfaceError("a zero-length vector has no direction")
    return (a[0] / norm, a[1] / norm, a[2] / norm)


def _triangle_area_and_normal(mesh: Mesh, face: Face) -> tuple[float, Vec3]:
    a, b, c = (mesh.vertices[index] for index in face)
    cross = _cross(_sub(b, a), _sub(c, a))
    norm = _length(cross)
    if norm <= 0.0:
        return 0.0, (0.0, 0.0, 0.0)
    return 0.5 * norm, (cross[0] / norm, cross[1] / norm, cross[2] / norm)


def _centroid(mesh: Mesh, face: Face) -> Vec3:
    a, b, c = (mesh.vertices[index] for index in face)
    return (
        (a[0] + b[0] + c[0]) / 3.0,
        (a[1] + b[1] + c[1]) / 3.0,
        (a[2] + b[2] + c[2]) / 3.0,
    )


def _bounding_extents(mesh: Mesh) -> tuple[float, float, float]:
    xs = [point[0] for point in mesh.vertices]
    ys = [point[1] for point in mesh.vertices]
    zs = [point[2] for point in mesh.vertices]
    return (max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs))


def _scale_of(mesh: Mesh) -> float:
    scale = max(_bounding_extents(mesh))
    if scale <= 0.0:
        raise MidSurfaceError("the mesh has no extent, so there is no length scale in it")
    return scale


def _require_closed_outward(mesh: Mesh) -> None:
    if not mesh.faces:
        raise MidSurfaceError("the mesh has no faces")
    if not meshlib.is_watertight(mesh):
        raise MidSurfaceError(
            "the surface is not closed: a mid-surface is the bisector of two opposed walls, so "
            "it needs a solid. Repair the mesh (tensorcad.mesh.weld) first"
        )
    if not meshlib.is_edge_manifold(mesh):
        raise MidSurfaceError("the surface is closed but not edge-manifold (an edge has >2 faces)")
    volume = meshlib.total_signed_volume(mesh)
    if volume <= 0.0:
        raise MidSurfaceError(
            f"the surface is wound inward (signed volume {volume!r} <= 0): run "
            "tensorcad.mesh.orient_outward first, because the wall is measured along -n"
        )


def suggest_max_thickness(mesh: Mesh) -> float:
    """A defensible default for ``max_thickness``: ``1.5 x`` the smallest bounding extent.

    A sheet-like solid's wall thickness is its smallest extent, so a wall is accepted up to
    that scale and rejected beyond it -- which is exactly what separates a plate's two large
    faces (paired at ``t``) from its four edge strips (whose opposed partner is the far face,
    ``L`` away). It is a guess about the *input*, not a result, so it is a separate function
    and every report carries the value that was used.
    """
    return 1.5 * min(_bounding_extents(mesh))


# --------------------------------------------------------------------------------------
# ray casting: a uniform grid over the triangles, then Moller-Trumbore
# --------------------------------------------------------------------------------------


def _ray_triangle(origin: Vec3, direction: Vec3, a: Vec3, b: Vec3, c: Vec3) -> float | None:
    """Moller-Trumbore; returns the distance (``direction`` is a unit vector)."""
    edge1 = _sub(b, a)
    edge2 = _sub(c, a)
    pvec = _cross(direction, edge2)
    det = _dot(edge1, pvec)
    if abs(det) < 1e-18:
        return None
    inv = 1.0 / det
    tvec = _sub(origin, a)
    u = _dot(tvec, pvec) * inv
    if u < -1e-12 or u > 1.0 + 1e-12:
        return None
    qvec = _cross(tvec, edge1)
    v = _dot(direction, qvec) * inv
    if v < -1e-12 or u + v > 1.0 + 1e-12:
        return None
    t = _dot(edge2, qvec) * inv
    return t if t > 0.0 else None


class _RayGrid:
    """Uniform-cell bucket grid of the triangles; a ray collects a candidate superset."""

    def __init__(self, mesh: Mesh, cell: float) -> None:
        self.mesh = mesh
        self.cell = cell
        self.origin = (
            min(point[0] for point in mesh.vertices),
            min(point[1] for point in mesh.vertices),
            min(point[2] for point in mesh.vertices),
        )
        self.buckets: dict[tuple[int, int, int], list[int]] = {}
        for index, face in enumerate(mesh.faces):
            points = [mesh.vertices[corner] for corner in face]
            low = [min(point[axis] for point in points) for axis in range(3)]
            high = [max(point[axis] for point in points) for axis in range(3)]
            for key in self._cells(low, high):
                self.buckets.setdefault(key, []).append(index)

    def _cells(self, low: Sequence[float], high: Sequence[float]):
        low_key = [int(math.floor((low[axis] - self.origin[axis]) / self.cell)) for axis in range(3)]
        high_key = [int(math.floor((high[axis] - self.origin[axis]) / self.cell)) for axis in range(3)]
        for i in range(low_key[0], high_key[0] + 1):
            for j in range(low_key[1], high_key[1] + 1):
                for k in range(low_key[2], high_key[2] + 1):
                    yield (i, j, k)

    def candidates(self, origin: Vec3, direction: Vec3, reach: float) -> list[int]:
        end = _add(origin, _scale(direction, reach))
        pad = self.cell
        low = [min(origin[axis], end[axis]) - pad for axis in range(3)]
        high = [max(origin[axis], end[axis]) + pad for axis in range(3)]
        seen: set[int] = set()
        for key in self._cells(low, high):
            bucket = self.buckets.get(key)
            if bucket:
                seen.update(bucket)
        return list(seen)

    def cast(
        self, origin: Vec3, direction: Vec3, reach: float, skip: int
    ) -> tuple[float, int] | None:
        """Nearest intersection beyond the source triangle; ``(distance, face)``."""
        best: tuple[float, int] | None = None
        for index in self.candidates(origin, direction, reach):
            if index == skip:
                continue
            a, b, c = (self.mesh.vertices[corner] for corner in self.mesh.faces[index])
            hit = _ray_triangle(origin, direction, a, b, c)
            if hit is None or hit > reach:
                continue
            if best is None or hit < best[0]:
                best = (hit, index)
        return best


# --------------------------------------------------------------------------------------
# step 1: the thickness map
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class ThicknessMap:
    """Per-face wall thickness found by an inward ray, with a reason for every failure.

    ``thickness[i]`` is the distance from the centroid of face ``i`` to the opposed wall (or
    ``None``); ``partner[i]`` is the face that was hit; ``reason[i]`` is empty for a pair and a
    sentence otherwise. ``mutual[i]`` records whether the reverse ray came back, which is the
    difference between "this is a wall" and "this ray happened to meet an anti-parallel face".
    """

    thickness: tuple[float | None, ...]
    partner: tuple[int | None, ...]
    reason: tuple[str, ...]
    mutual: tuple[bool, ...]
    face_area: tuple[float, ...]
    max_thickness: float
    scale: float

    @property
    def paired_count(self) -> int:
        return sum(1 for value in self.thickness if value is not None)

    @property
    def qualified_faces(self) -> tuple[int, ...]:
        """Faces whose pair is mutual -- the faces a mid-surface can be taken from."""
        return tuple(index for index, ok in enumerate(self.mutual) if ok)

    def qualified_area(self) -> float:
        return sum(self.face_area[index] for index in self.qualified_faces)

    def total_area(self) -> float:
        return sum(self.face_area)

    def coverage(self) -> float:
        """Qualified area / total area: how much of the surface is a measurable wall."""
        total = self.total_area()
        return self.qualified_area() / total if total > 0.0 else 0.0

    def thickness_values(self) -> tuple[float, ...]:
        return tuple(value for value in self.thickness if value is not None)

    def as_dict(self) -> dict[str, object]:
        values = self.thickness_values()
        return {
            "faces": len(self.thickness),
            "paired_faces": self.paired_count,
            "qualified_faces": len(self.qualified_faces),
            "coverage": self.coverage(),
            "max_thickness": self.max_thickness,
            "min_thickness": min(values) if values else None,
            "mean_thickness": (sum(values) / len(values)) if values else None,
            "max_thickness_seen": max(values) if values else None,
            "unpaired_reasons": sorted({reason for reason in self.reason if reason}),
        }


def measure_thickness(
    mesh: Mesh,
    *,
    max_thickness: float | None = None,
    opposition: float = 0.5,
    require_mutual: bool = True,
) -> ThicknessMap:
    """Cast an inward ray from every face and pair it with the opposed wall it meets.

    ``opposition`` is the minimum ``-dot(n_source, n_hit)`` for a hit to count as the other side
    of a wall: ``0.5`` accepts walls up to 60 degrees off parallel -- a faceted cylinder
    (offsetting keeps the facet normals exactly anti-parallel) and a slightly conical wall -- and
    rejects a perpendicular neighbour (a box corner, ``dot = 0``). Pairs must be *mutual* by
    default, so a ray that reaches a genuinely different part of the body is reported rather
    than believed.
    """
    _require_closed_outward(mesh)
    scale = _scale_of(mesh)
    limit = suggest_max_thickness(mesh) if max_thickness is None else float(max_thickness)
    if not (limit > 0.0) or not math.isfinite(limit):
        raise MidSurfaceError(f"max_thickness must be positive and finite: {max_thickness!r}")
    areas: list[float] = []
    normals: list[Vec3] = []
    for face in mesh.faces:
        area, normal = _triangle_area_and_normal(mesh, face)
        areas.append(area)
        normals.append(normal)
    mean_edge = math.sqrt(sum(areas) / len(areas)) if areas else 0.0
    cell = max(limit, mean_edge, scale / 64.0)
    grid = _RayGrid(mesh, cell)
    epsilon = 1e-9 * scale

    forward_distance: list[float | None] = []
    forward_partner: list[int | None] = []
    forward_reason: list[str] = []
    for index, face in enumerate(mesh.faces):
        normal = normals[index]
        if _length(normal) == 0.0:
            forward_distance.append(None)
            forward_partner.append(None)
            forward_reason.append("degenerate face: no normal to cast along")
            continue
        origin = _sub(_centroid(mesh, face), _scale(normal, epsilon))
        hit = grid.cast(origin, _scale(normal, -1.0), limit + epsilon, index)
        if hit is None:
            forward_distance.append(None)
            forward_partner.append(None)
            forward_reason.append("no other face within max_thickness along -n")
            continue
        distance, other = hit
        alignment = -_dot(normal, normals[other])
        if alignment < opposition:
            forward_distance.append(None)
            forward_partner.append(None)
            forward_reason.append(
                f"nearest face {other} is not opposed (-dot = {alignment:+.3f} < {opposition})"
            )
            continue
        forward_distance.append(distance + epsilon)
        forward_partner.append(other)
        forward_reason.append("")

    mutual: list[bool] = []
    thickness: list[float | None] = []
    partner: list[int | None] = []
    reason: list[str] = []
    for index in range(len(mesh.faces)):
        other = forward_partner[index]
        if other is None:
            mutual.append(False)
            thickness.append(None)
            partner.append(None)
            reason.append(forward_reason[index])
            continue
        if require_mutual and forward_partner[other] != index:
            mutual.append(False)
            thickness.append(None)
            partner.append(other)
            reason.append(
                f"the pair is not mutual: face {other} pairs with {forward_partner[other]} "
                f"instead of {index}"
            )
            continue
        mutual.append(True)
        thickness.append(forward_distance[index])
        partner.append(other)
        reason.append("")

    return ThicknessMap(
        thickness=tuple(thickness),
        partner=tuple(partner),
        reason=tuple(reason),
        mutual=tuple(mutual),
        face_area=tuple(areas),
        max_thickness=limit,
        scale=scale,
    )


# --------------------------------------------------------------------------------------
# step 2: panels (connected, coplanar, mutually-paired face sets)
# --------------------------------------------------------------------------------------


class _UnionFind:
    def __init__(self, size: int) -> None:
        self.parent = list(range(size))

    def find(self, item: int) -> int:
        root = item
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[item] != root:
            self.parent[item], item = root, self.parent[item]
        return root

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[max(ra, rb)] = min(ra, rb)


def _vertex_representatives(mesh: Mesh, tolerance: float) -> list[int]:
    """Map every vertex to a canonical id shared by vertices within ``tolerance``."""
    cell = tolerance if tolerance > 0.0 else 1e-12
    buckets: dict[tuple[int, int, int], list[tuple[Vec3, int]]] = {}
    representatives: list[int] = []
    next_id = len(mesh.vertices)
    for point in mesh.vertices:
        key = (
            int(math.floor(point[0] / cell)),
            int(math.floor(point[1] / cell)),
            int(math.floor(point[2] / cell)),
        )
        found: int | None = None
        for dx in (-1, 0, 1):
            if found is not None:
                break
            for dy in (-1, 0, 1):
                if found is not None:
                    break
                for dz in (-1, 0, 1):
                    for other, representative in buckets.get(
                        (key[0] + dx, key[1] + dy, key[2] + dz), ()
                    ):
                        if _length(_sub(point, other)) <= tolerance:
                            found = representative
                            break
                    if found is not None:
                        break
        if found is None:
            found = next_id
            next_id += 1
            buckets.setdefault(key, []).append((point, found))
        representatives.append(found)
    return representatives


@dataclass(frozen=True)
class Panel:
    """One flat patch of the mid-surface: a mid-plane, a measured thickness, a polygon.

    ``faces`` are the source faces of the measured side (the smaller-indexed face of each mutual
    pair, so a hollow body yields one sheet and not two). ``polygon`` is the mid-surface patch
    in order, on the mid-plane ``normal . x = mid_offset``; ``area`` is the source patch area and
    ``sheet_area`` the polygon's, which differ only by the corner clipping.
    """

    normal: Vec3
    offset: float
    mid_offset: float
    thickness: float
    thickness_min: float
    thickness_max: float
    faces: tuple[int, ...]
    area: float
    sheet_area: float
    polygon: tuple[Vec3, ...]

    @property
    def thickness_spread(self) -> float:
        return self.thickness_max - self.thickness_min

    def as_dict(self) -> dict[str, object]:
        return {
            "normal": list(self.normal),
            "offset": self.offset,
            "mid_offset": self.mid_offset,
            "thickness": self.thickness,
            "thickness_spread": self.thickness_spread,
            "faces": len(self.faces),
            "area": self.area,
            "sheet_area": self.sheet_area,
            "polygon": [list(point) for point in self.polygon],
        }


def _panel_plane(
    mesh: Mesh, faces: Sequence[int], areas: Sequence[float], normals: Sequence[Vec3]
) -> tuple[Vec3, float, float]:
    """Area-weighted plane of a face set: ``(normal, offset, area)``."""
    accumulated = [0.0, 0.0, 0.0]
    weighted = [0.0, 0.0, 0.0]
    total = 0.0
    for face in faces:
        normal = normals[face]
        area = areas[face]
        total += area
        for axis in range(3):
            accumulated[axis] += 2.0 * area * normal[axis]
        centroid = _centroid(mesh, mesh.faces[face])
        for axis in range(3):
            weighted[axis] += area * centroid[axis]
    if total <= 0.0:
        raise MidSurfaceError("a panel of zero area has no plane")
    normal = _unit((accumulated[0], accumulated[1], accumulated[2]))
    centroid = (weighted[0] / total, weighted[1] / total, weighted[2] / total)
    return normal, _dot(normal, centroid), total


def _basis(normal: Vec3) -> tuple[Vec3, Vec3]:
    axis = (1.0, 0.0, 0.0) if abs(normal[0]) <= 0.9 else (0.0, 1.0, 0.0)
    u = _unit(_sub(axis, _scale(normal, _dot(normal, axis))))
    return u, _cross(normal, u)


def _convex_hull_2d(points: Sequence[tuple[float, float]]) -> list[tuple[float, float]]:
    ordered = sorted(set(points))
    if len(ordered) <= 2:
        return list(ordered)

    def cross(o, a, b) -> float:
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower: list[tuple[float, float]] = []
    for point in ordered:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0.0:
            lower.pop()
        lower.append(point)
    upper: list[tuple[float, float]] = []
    for point in reversed(ordered):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0.0:
            upper.pop()
        upper.append(point)
    return lower[:-1] + upper[:-1]


def _interpolate_2d(
    a: tuple[float, float], b: tuple[float, float], va: float, vb: float
) -> tuple[float, float]:
    denominator = va - vb
    if abs(denominator) < 1e-300:
        return a
    s = va / denominator
    return (a[0] + (b[0] - a[0]) * s, a[1] + (b[1] - a[1]) * s)


def _clip_half_plane(
    polygon: list[tuple[float, float]], a: float, b: float, c: float, sign: float
) -> list[tuple[float, float]]:
    """Sutherland-Hodgman against ``sign * (a x + b y + c) >= 0``."""
    if not polygon:
        return []
    result: list[tuple[float, float]] = []

    def value(point: tuple[float, float]) -> float:
        return sign * (a * point[0] + b * point[1] + c)

    previous = polygon[-1]
    previous_value = value(previous)
    for current in polygon:
        current_value = value(current)
        if current_value >= 0.0:
            if previous_value < 0.0:
                result.append(_interpolate_2d(previous, current, previous_value, current_value))
            result.append(current)
        elif previous_value >= 0.0:
            result.append(_interpolate_2d(previous, current, previous_value, current_value))
        previous, previous_value = current, current_value
    return result


def _shoelace(polygon: Sequence[tuple[float, float]]) -> float:
    total = 0.0
    for index, point in enumerate(polygon):
        nxt = polygon[(index + 1) % len(polygon)]
        total += point[0] * nxt[1] - nxt[0] * point[1]
    return 0.5 * total


def _centroid_of_faces(mesh: Mesh, faces: Sequence[int], areas: Sequence[float]) -> Vec3:
    total = 0.0
    accumulated = [0.0, 0.0, 0.0]
    for face in faces:
        weight = areas[face]
        total += weight
        centroid = _centroid(mesh, mesh.faces[face])
        for axis in range(3):
            accumulated[axis] += weight * centroid[axis]
    if total <= 0.0:
        return (0.0, 0.0, 0.0)
    return (accumulated[0] / total, accumulated[1] / total, accumulated[2] / total)


def _build_panel(
    mesh: Mesh,
    faces: Sequence[int],
    areas: Sequence[float],
    normals: Sequence[Vec3],
    thickness_map: ThicknessMap,
    clippers: Sequence[tuple[Vec3, float]],
    convexity_tolerance: float,
) -> Panel:
    normal, offset, area = _panel_plane(mesh, faces, areas, normals)
    side_thickness = [
        thickness_map.thickness[face] for face in faces if thickness_map.thickness[face] is not None
    ]
    if not side_thickness:
        raise MidSurfaceError("a panel without a measured thickness cannot be offset")
    thickness = sum(side_thickness) / len(side_thickness)
    mid_offset = offset - 0.5 * thickness
    u, v = _basis(normal)

    def project(point: Vec3) -> tuple[float, float]:
        shifted = _sub(point, _scale(normal, _dot(normal, point) - mid_offset))
        return (_dot(shifted, u), _dot(shifted, v))

    vertices = sorted({index for face in faces for index in mesh.faces[face]})
    flat = [project(mesh.vertices[index]) for index in vertices]
    hull = _convex_hull_2d(flat)
    if not hull:
        raise MidSurfaceError("a panel degenerates to no polygon once flattened")
    # The hull must be the panel: for a flat convex patch the projected triangle areas sum to
    # the same number. A bigger hull means the patch is not convex and the sheet would cover
    # material that is not there -- refused rather than approximated.
    hull_area = abs(_shoelace(hull))
    if area > 0.0:
        overshoot = hull_area / area - 1.0
        if overshoot > convexity_tolerance:
            raise MidSurfaceError(
                f"a panel is not convex once flattened: its convex hull is "
                f"{100.0 * overshoot:.3f}% larger than the {area!r} m^2 of its {len(faces)} "
                "faces. The mid-surface is only built from convex patches; split the panel"
            )

    centroid = project(_centroid_of_faces(mesh, faces, areas))
    for clip_normal, clip_offset in clippers:
        a = _dot(clip_normal, u)
        b = _dot(clip_normal, v)
        c = mid_offset * _dot(clip_normal, normal) - clip_offset
        hint = a * centroid[0] + b * centroid[1] + c
        if abs(hint) < 1e-13 * max(1.0, abs(clip_offset)):
            # The panel centre sits on the neighbour's mid-plane: fall back to the vertex that
            # is furthest from it, which lies on the panel's own side by construction.
            furthest = max(flat, key=lambda point: abs(a * point[0] + b * point[1] + c))
            hint = a * furthest[0] + b * furthest[1] + c
            if abs(hint) < 1e-300:
                continue
        sign = 1.0 if hint >= 0.0 else -1.0
        hull = _clip_half_plane(hull, a, b, c, sign)
        if len(hull) < 3:
            raise MidSurfaceError(
                "a panel is clipped away entirely by its neighbours' mid-planes: the panels "
                "overlap, so the thickness used for the offset is wrong"
            )
    if _shoelace(hull) < 0.0:
        hull.reverse()
    sheet_area = abs(_shoelace(hull))
    polygon = tuple(
        _add(_add(_scale(u, point[0]), _scale(v, point[1])), _scale(normal, mid_offset))
        for point in hull
    )
    return Panel(
        normal=normal,
        offset=offset,
        mid_offset=mid_offset,
        thickness=thickness,
        thickness_min=min(side_thickness),
        thickness_max=max(side_thickness),
        faces=tuple(sorted(faces)),
        area=area,
        sheet_area=sheet_area,
        polygon=polygon,
    )


# --------------------------------------------------------------------------------------
# the assembled result
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class MidSurface:
    """The sheet, its patches, and every number needed to judge whether to trust it."""

    mesh: Mesh
    panels: tuple[Panel, ...]
    thickness_map: ThicknessMap
    source: Mesh
    source_faces: int
    source_area: float
    thickness: float
    angle_tolerance: float
    offset_tolerance: float
    weld_tolerance: float
    thickness_tolerance: float
    convexity_tolerance: float
    min_coverage: float
    min_panel_area: float

    @property
    def sheet_faces(self) -> int:
        return len(self.mesh.faces)

    def face_reduction(self) -> float:
        """``1 - sheet_faces / source_faces``: the mesh saving the model gets."""
        if self.source_faces <= 0:
            return 0.0
        return 1.0 - len(self.mesh.faces) / self.source_faces

    def area(self) -> float:
        return meshlib.surface_area(self.mesh)

    def sheet_area(self) -> float:
        return sum(panel.sheet_area for panel in self.panels)

    def source_patch_area(self) -> float:
        return sum(panel.area for panel in self.panels)

    def removed_area(self) -> float:
        """Source area that entered no panel (unpaired walls, suppressed features)."""
        return self.thickness_map.total_area() - self.source_patch_area()

    def volume(self) -> float:
        """Signed volume of the sheet: zero for an open sheet, the eroded volume for a closed one."""
        return meshlib.total_signed_volume(self.mesh)

    def watertight(self) -> bool:
        return meshlib.is_watertight(self.mesh)

    def thickness_values(self) -> tuple[float, ...]:
        return tuple(panel.thickness for panel in self.panels)

    def as_dict(self) -> dict[str, object]:
        values = self.thickness_values()
        return {
            "source_faces": self.source_faces,
            "sheet_faces": self.sheet_faces,
            "face_reduction": self.face_reduction(),
            "panels": len(self.panels),
            "thickness": self.thickness,
            "thickness_min": min(values) if values else None,
            "thickness_max": max(values) if values else None,
            "source_area": self.source_area,
            "source_patch_area": self.source_patch_area(),
            "sheet_area": self.sheet_area(),
            "mesh_area": self.area(),
            "removed_area": self.removed_area(),
            "signed_volume": self.volume(),
            "watertight": self.watertight(),
            "panels_data": [panel.as_dict() for panel in self.panels],
            "thickness_map": self.thickness_map.as_dict(),
        }


def _sheet_mesh(panels: Sequence[Panel], units: str, weld_tolerance: float) -> Mesh:
    vertices: list[Vec3] = []
    faces: list[Face] = []
    for panel in panels:
        base = len(vertices)
        vertices.extend(panel.polygon)
        for index in range(1, len(panel.polygon) - 1):
            faces.append((base, base + index, base + index + 1))
    sheet = Mesh(vertices=tuple(vertices), faces=tuple(faces), units=units)
    if weld_tolerance > 0.0 and faces:
        sheet = meshlib.weld(sheet, tolerance=weld_tolerance)
    return sheet


def _group_faces(
    mesh: Mesh,
    reference: Sequence[int],
    planes: dict[int, tuple[Vec3, float]],
    *,
    angle_tolerance: float,
    offset_tolerance: float,
    weld_tolerance: float,
) -> list[list[int]]:
    """Connected, coplanar components of the ``reference`` faces."""
    if not reference:
        return []
    representatives = _vertex_representatives(mesh, weld_tolerance)
    faces_at_position: dict[int, list[int]] = {}
    for index in reference:
        for corner in mesh.faces[index]:
            faces_at_position.setdefault(representatives[corner], []).append(index)
    labels = {face: position for position, face in enumerate(reference)}
    union = _UnionFind(len(reference))
    # cos() is flat at 1.0 for small angles, so the comparison needs a floor at round-off:
    # two triangles of one facet produce normals that agree to ~1e-16, not to 1.0 exactly, and a
    # strict `dot < cos(tol)` test would split every quad into two panels.
    cosine = 1.0 - max(2e-16, 1.0 - math.cos(angle_tolerance))
    for faces in faces_at_position.values():
        if len(faces) < 2:
            continue
        # Every pair, not only pairs with the first face: a vertex is shared by faces of every
        # plane that meets there, so anchoring on any one of them would miss the others.
        for first in range(len(faces)):
            for second in range(first + 1, len(faces)):
                anchor = faces[first]
                other = faces[second]
                anchor_normal, anchor_offset = planes[anchor]
                other_normal, other_offset = planes[other]
                if _dot(anchor_normal, other_normal) < cosine:
                    continue
                if abs(anchor_offset - other_offset) > offset_tolerance:
                    continue
                union.union(labels[anchor], labels[other])
    groups: dict[int, list[int]] = {}
    for face in reference:
        groups.setdefault(union.find(labels[face]), []).append(face)
    return [sorted(group) for group in groups.values()]


def _panel_adjacency(
    mesh: Mesh, groups: Sequence[Sequence[int]], weld_tolerance: float
) -> list[set[int]]:
    representatives = _vertex_representatives(mesh, weld_tolerance)
    position_to_groups: dict[int, set[int]] = {}
    for group_index, group in enumerate(groups):
        for face in group:
            for corner in mesh.faces[face]:
                position_to_groups.setdefault(representatives[corner], set()).add(group_index)
    adjacency: list[set[int]] = [set() for _ in groups]
    for group_index, group in enumerate(groups):
        for face in group:
            for corner in mesh.faces[face]:
                adjacency[group_index].update(position_to_groups[representatives[corner]])
        adjacency[group_index].discard(group_index)
    return adjacency


def _build(
    mesh: Mesh,
    thickness_map: ThicknessMap,
    *,
    angle_tolerance: float,
    offset_tolerance: float,
    weld_tolerance: float,
    thickness_tolerance: float,
    convexity_tolerance: float,
    min_coverage: float,
    min_panel_area: float,
) -> MidSurface:
    areas = thickness_map.face_area
    normals = [_triangle_area_and_normal(mesh, face)[1] for face in mesh.faces]
    planes = {
        face: (normals[face], _dot(normals[face], _centroid(mesh, mesh.faces[face])))
        for face in range(len(mesh.faces))
    }
    # One reference face per mutual pair -- the smaller index, so a hollow body (its outer
    # faces first) produces the outer side of the sheet exactly once.
    reference = sorted(
        index
        for index in thickness_map.qualified_faces
        if thickness_map.partner[index] is not None and index < thickness_map.partner[index]  # type: ignore[operator]
    )
    total_area = thickness_map.total_area()
    if not reference:
        raise MidSurfaceError(
            "no opposed thin walls found: every face's inward ray either left the body or met "
            "a face whose normal points the same way (qualified area fraction "
            f"{thickness_map.coverage():.4f}). This is not a thin-walled solid, or "
            f"max_thickness = {thickness_map.max_thickness!r} is wrong for it"
        )
    # Coverage is the *qualified* (mutual-pair) area, not the reference subset: every pair is
    # used once for the sheet, so counting only the reference faces would report half of a fully
    # thin-walled body as "not a wall" and refuse it.
    coverage = thickness_map.coverage()
    if coverage < min_coverage:
        raise MidSurfaceError(
            f"only {100.0 * coverage:.2f}% of the surface is a measurable wall (min_coverage = "
            f"{100.0 * min_coverage:.2f}%): the body is mostly not sheet-like, so a mid-surface "
            "of it would describe something else"
        )
    groups = _group_faces(
        mesh,
        reference,
        planes,
        angle_tolerance=angle_tolerance,
        offset_tolerance=offset_tolerance,
        weld_tolerance=weld_tolerance,
    )
    adjacency = _panel_adjacency(mesh, groups, weld_tolerance)
    group_planes = [_panel_plane(mesh, group, areas, normals) for group in groups]
    group_thickness: list[float] = []
    for index, group in enumerate(groups):
        values = [
            thickness_map.thickness[face]
            for face in group
            if thickness_map.thickness[face] is not None
        ]
        mean = sum(values) / len(values) if values else 0.0
        if values and mean > 0.0:
            spread = max(values) - min(values)
            if spread > thickness_tolerance * mean:
                plane_normal = group_planes[index][0]
                raise MidSurfaceError(
                    f"the wall thickness inside one panel is not constant: {min(values)!r} to "
                    f"{max(values)!r} m (spread {100.0 * spread / mean:.2f}% of the mean, "
                    f"tolerance {100.0 * thickness_tolerance:.2f}%). A constant-offset sheet "
                    "cannot represent a taper; split the panel or loosen the panel grouping. "
                    f"Panel normal {tuple(round(value, 6) for value in plane_normal)}"
                )
        group_thickness.append(mean)

    kept = {
        index
        for index, group in enumerate(groups)
        if sum(areas[face] for face in group) >= min_panel_area
    }
    panels: list[Panel] = []
    for index, group in enumerate(groups):
        if index not in kept:
            continue
        clippers = [
            (group_planes[other][0], group_planes[other][1] - 0.5 * group_thickness[other])
            for other in sorted(adjacency[index])
            if other in kept
        ]
        panels.append(
            _build_panel(
                mesh,
                group,
                areas,
                normals,
                thickness_map,
                clippers,
                convexity_tolerance,
            )
        )
    if not panels:
        raise MidSurfaceError(
            f"min_panel_area = {min_panel_area!r} suppressed every panel; nothing is left"
        )
    sheet = _sheet_mesh(panels, mesh.units, weld_tolerance)
    values = [panel.thickness for panel in panels]
    return MidSurface(
        mesh=sheet,
        panels=tuple(panels),
        thickness_map=thickness_map,
        source=mesh,
        source_faces=len(mesh.faces),
        source_area=meshlib.surface_area(mesh),
        thickness=sum(values) / len(values),
        angle_tolerance=angle_tolerance,
        offset_tolerance=offset_tolerance,
        weld_tolerance=weld_tolerance,
        thickness_tolerance=thickness_tolerance,
        convexity_tolerance=convexity_tolerance,
        min_coverage=min_coverage,
        min_panel_area=min_panel_area,
    )


def extract_midsurface(
    mesh: Mesh,
    *,
    max_thickness: float | None = None,
    opposition: float = 0.5,
    angle_tolerance: float = 1e-9,
    offset_tolerance: float = 1e-9,
    weld_tolerance: float | None = None,
    thickness_tolerance: float = 0.05,
    convexity_tolerance: float = 1e-9,
    min_coverage: float = 0.5,
    min_panel_area: float = 0.0,
) -> MidSurface:
    """Extract the mid-surface sheet of a closed thin-walled solid.

    ``max_thickness`` bounds what counts as a wall (see :func:`suggest_max_thickness`);
    ``weld_tolerance`` defaults to ``1e-9 x`` the model size and is the distance at which two
    vertices are the same point, i.e. how wide a CAD seam is bridged; ``thickness_tolerance`` is
    the allowed relative spread of the measured wall inside one panel before a taper is
    refused; ``min_coverage`` is the fraction of the surface that must be a measurable wall.
    """
    thickness_map = measure_thickness(mesh, max_thickness=max_thickness, opposition=opposition)
    weld = 1e-9 * thickness_map.scale if weld_tolerance is None else float(weld_tolerance)
    return _build(
        mesh,
        thickness_map,
        angle_tolerance=angle_tolerance,
        offset_tolerance=offset_tolerance,
        weld_tolerance=weld,
        thickness_tolerance=thickness_tolerance,
        convexity_tolerance=convexity_tolerance,
        min_coverage=min_coverage,
        min_panel_area=min_panel_area,
    )


def defeature(
    midsurface: MidSurface,
    *,
    weld_tolerance: float | None = None,
    angle_tolerance: float | None = None,
    offset_tolerance: float | None = None,
    thickness_tolerance: float | None = None,
    convexity_tolerance: float | None = None,
    min_panel_area: float | None = None,
) -> MidSurface:
    """Re-extract with relaxed tolerances: bridge import seams, drop features below an area.

    Every argument defaults to the value the input used, so ``defeature(m, min_panel_area=a)``
    is exactly "suppress features smaller than ``a``" and nothing else. ``convexity_tolerance``
    has to be relaxed along with ``weld_tolerance``: bridging a gap makes the hull cover it, and
    that covered sliver is exactly what the convexity check measures. The result is a fresh
    :class:`MidSurface`; differencing its :meth:`MidSurface.face_reduction` and
    :meth:`MidSurface.removed_area` against the input reports what the simplification bought and
    cost. The numbers come from the caller, this function does not promise any.
    """
    return _build(
        midsurface.source,
        midsurface.thickness_map,
        angle_tolerance=midsurface.angle_tolerance if angle_tolerance is None else angle_tolerance,
        offset_tolerance=midsurface.offset_tolerance
        if offset_tolerance is None
        else offset_tolerance,
        weld_tolerance=midsurface.weld_tolerance if weld_tolerance is None else weld_tolerance,
        thickness_tolerance=midsurface.thickness_tolerance
        if thickness_tolerance is None
        else thickness_tolerance,
        convexity_tolerance=midsurface.convexity_tolerance
        if convexity_tolerance is None
        else convexity_tolerance,
        min_coverage=midsurface.min_coverage,
        min_panel_area=midsurface.min_panel_area if min_panel_area is None else min_panel_area,
    )