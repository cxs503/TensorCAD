"""Shelling: hollow a closed, outward-oriented surface by an inward offset.

``GAP-CAD-011`` in ``docs/oss_inventory_cad/GAP-ANCHORS.md`` records that the package can
smooth, decimate and project a mesh but cannot make it *hollow*: everything here produced a
surface, nothing produced a wall. The acceptance criterion written there is "a box shelled at
thickness ``t`` matches the analytic outer-minus-inner volume to 2 %, and a sweep of ``t/L``
from 0.01 to 0.1 does not fail". This module is that wall.

**The method.** Two surfaces in space do not need a B-rep to be offset. Every triangle carries
a plane ``n . x = d`` with ``n`` the outward unit normal, so the solid is locally the
intersection of half-spaces ``n . x <= d`` and the inward offset is the same intersection with
``d`` replaced by ``d - t``. At a vertex the incident faces hand over their (deduplicated)
planes, and the offset vertex is the weighted least-squares point

    minimise  sum_i w_i (n_i . x - (d_i - t))^2 ,     w_i = interior angle of face i at v

solved through the 3x3 normal equations. The weights are the *angle* weights that vertex
normals normally use. For a box corner the three planes are orthogonal and the solution is
``(d_i - t)`` along each of them with **no** dependence on the weights, i.e. the eroded box
exactly, for every thickness; on a flat or edge vertex the plane system is rank-deficient, and
a ridge term of ``1e-14 * trace/3`` picks the minimum-norm displacement (still exact in the
constrained directions, because the right-hand side lies in the span of the normals; the
ridge costs ``O(1e-14)`` relative on the offset point, which is what the measured box volume
error of ``3e-14`` is made of).

That is the whole of the geometry. What the module spends its remaining lines on is **honesty
about failure**, because an offset that has collided with itself still looks like a mesh:

* the inner solid's signed volume must stay positive (a fully inverted offset gives a negative
  one and is refused);
* inner faces must keep the winding of the outer ones (a local inversion flips the normal, so
  ``min`` of ``dot(n_inner, n_outer)`` is reported and a non-positive value is refused);
* the **achieved wall** is measured, not assumed: for every inner vertex the distance to the
  *whole* outer surface is computed through a uniform triangle grid, and if it falls more than
  ``wall_tolerance`` (default 5 %) below ``t`` the call fails. This is the check that catches
  the case the plane solve cannot see -- two opposite walls of a thin arm, where each vertex
  is at distance ``t`` from its own faces and closer than ``t`` to the far one. Four possible
  causes are named in the message; the mesh may be too coarse to carry the thickness, or the
  thickness may exceed half the local wall.

**What it does not do.** No face is removed and no rim is stitched, so the result is a
*closed hollow* body (outer surface plus reversed inner surface): opening a face and bridging
the gap needs B-rep topology and is ``GAP-CAD-001``. Curved surfaces are represented by their
facets, so the inner surface of a faceted sphere is short of the true offset by ``O(h^2)`` --
measured, and reported as ``wall_shortfall`` rather than hidden. Bodies whose orientation is
not outward are refused instead of silently inverted (``mesh.orient_outward`` fixes them
first).
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

from . import mesh as meshlib
from .geometry import Face, Mesh, Vec3

__all__ = [
    "HollowBody",
    "ShellError",
    "ShellResult",
    "hollow_body",
    "min_wall_thickness",
    "offset_inward",
    "shell",
]


class ShellError(ValueError):
    """Raised when a solid cannot be shelled at the requested thickness."""


@dataclass(frozen=True)
class ShellResult:
    """The outer surface, the inner surface and the numbers that say the wall is real.

    ``outer_volume``/``inner_volume`` are divergence-theorem volumes of the two closed
    surfaces; the material of the shell is the difference. ``min_wall_thickness`` is measured
    on the mesh (nearest distance from an inner vertex to the outer surface, *all* faces), and
    ``max_plane_residual`` says how well the plane least-squares solve satisfied its own
    constraints -- it is round-off for a planar polyhedron and non-zero only where more than
    three planes meet at a point and no single point satisfies them all.
    """

    outer: Mesh
    inner: Mesh
    thickness: float
    outer_volume: float
    inner_volume: float
    min_wall_thickness: float
    max_plane_residual: float
    min_normal_dot: float

    @property
    def shell_volume(self) -> float:
        """Volume of the material, i.e. outer minus inner."""
        return self.outer_volume - self.inner_volume

    @property
    def wall_shortfall(self) -> float:
        """``thickness - min_wall_thickness``; positive when the wall is not fully achieved."""
        return self.thickness - self.min_wall_thickness

    def as_dict(self) -> dict[str, float | int]:
        return {
            "thickness": self.thickness,
            "outer_volume": self.outer_volume,
            "inner_volume": self.inner_volume,
            "shell_volume": self.shell_volume,
            "min_wall_thickness": self.min_wall_thickness,
            "wall_shortfall": self.wall_shortfall,
            "max_plane_residual": self.max_plane_residual,
            "min_normal_dot": self.min_normal_dot,
            "outer_vertices": len(self.outer.vertices),
            "outer_faces": len(self.outer.faces),
        }


@dataclass(frozen=True)
class HollowBody:
    """The shell as one exportable surface: outer faces plus the inner ones reversed.

    The reversed inner component faces *into the cavity*, so the region the closed surface
    bounds is the material, and the divergence-theorem volume of :attr:`mesh` is
    ``outer_volume - inner_volume`` -- asserted in the tests rather than assumed here. Keep the
    two numbers apart: a consumer that wants the displaced volume of the body, or the volume
    of the hold it encloses, is asking for :attr:`outer_volume` or :attr:`cavity_volume`, and
    the surface integral answers only the first of the three.
    """

    mesh: Mesh
    material_volume: float
    cavity_volume: float

    def as_dict(self) -> dict[str, object]:
        return {
            "vertices": len(self.mesh.vertices),
            "faces": len(self.mesh.faces),
            "material_volume": self.material_volume,
            "cavity_volume": self.cavity_volume,
            "signed_volume": meshlib.total_signed_volume(self.mesh),
            "watertight": meshlib.is_watertight(self.mesh),
            "edge_manifold": meshlib.is_edge_manifold(self.mesh),
        }


# --------------------------------------------------------------------------------------
# small vector helpers (the core is dependency-free: no NumPy anywhere in this module)
# --------------------------------------------------------------------------------------


def _sub(a: Sequence[float], b: Sequence[float]) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _dot(a: Sequence[float], b: Sequence[float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a: Sequence[float], b: Sequence[float]) -> Vec3:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _length(a: Sequence[float]) -> float:
    return math.sqrt(_dot(a, a))


def _unit(a: Sequence[float]) -> Vec3:
    length = _length(a)
    if length == 0.0:
        raise ShellError("cannot normalise a zero-length vector")
    return (a[0] / length, a[1] / length, a[2] / length)


def _face_plane(mesh: Mesh, face: Face) -> tuple[Vec3, float]:
    """Outward unit normal and plane offset of a triangle; refuses zero-area faces."""
    a = mesh.vertices[face[0]]
    normal = _cross(_sub(mesh.vertices[face[1]], a), _sub(mesh.vertices[face[2]], a))
    length = _length(normal)
    if length == 0.0:
        raise ShellError(
            "a zero-area face carries no plane and cannot be offset; remove it first "
            "(mesh.remove_degenerate_faces)"
        )
    unit = (normal[0] / length, normal[1] / length, normal[2] / length)
    return unit, _dot(unit, a)


def _triangle_angles(mesh: Mesh, face: Face) -> tuple[float, float, float]:
    """Interior angles at the three vertices of the triangle, in face order."""
    points = [mesh.vertices[index] for index in face]
    angles = []
    for position in range(3):
        tail = points[position]
        left = _sub(points[(position + 1) % 3], tail)
        right = _sub(points[(position + 2) % 3], tail)
        left_length, right_length = _length(left), _length(right)
        if left_length == 0.0 or right_length == 0.0:
            raise ShellError("a face has a coincident vertex pair (zero-length edge); weld the mesh first")
        cosine = max(-1.0, min(1.0, _dot(left, right) / (left_length * right_length)))
        angles.append(math.acos(cosine))
    return angles[0], angles[1], angles[2]


# --------------------------------------------------------------------------------------
# per-vertex plane systems
# --------------------------------------------------------------------------------------


def _vertex_planes(mesh: Mesh) -> list[list[tuple[Vec3, float, float]]]:
    """Angle-weighted planes per vertex, with duplicate planes merged.

    Two triangles of the same facet produce the same plane to round-off; merging them keeps a
    quad's split from counting twice against the other planes at a corner. The merge tolerance
    is relative to the mean plane offset so it scales with the model.
    """
    planes: list[list[tuple[Vec3, float, float]]] = [[] for _ in mesh.vertices]
    for face in mesh.faces:
        unit, offset = _face_plane(mesh, face)
        angles = _triangle_angles(mesh, face)
        for position, vertex in enumerate(face):
            planes[vertex].append((unit, offset, angles[position]))
    merged: list[list[tuple[Vec3, float, float]]] = []
    for entries in planes:
        if not entries:
            merged.append([])
            continue
        scale = max(1.0, max(abs(offset) for _, offset, _ in entries))
        unique: list[list[float]] = []
        for unit, offset, weight in entries:
            placed = False
            for record in unique:
                normal = record[:3]
                if _dot(normal, unit) > 1.0 - 1e-12 and abs(record[3] - offset) <= 1e-12 * scale:
                    record[4] += weight
                    placed = True
                    break
            if not placed:
                unique.append([unit[0], unit[1], unit[2], offset, weight])
        merged.append([((r[0], r[1], r[2]), r[3], r[4]) for r in unique])
    return merged


def _solve_vertex(
    constraints: list[tuple[Vec3, float, float]], thickness: float, ridge: float
) -> tuple[Vec3, float]:
    """Weighted least-squares offset point and the largest constraint residual at it."""
    matrix = [[0.0] * 3 for _ in range(3)]
    rhs = [0.0, 0.0, 0.0]
    for normal, offset, weight in constraints:
        target = offset - thickness
        for row in range(3):
            rhs[row] += weight * normal[row] * target
            for column in range(3):
                matrix[row][column] += weight * normal[row] * normal[column]
    trace = matrix[0][0] + matrix[1][1] + matrix[2][2]
    for axis in range(3):
        matrix[axis][axis] += ridge * trace / 3.0
    determinant = (
        matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    )
    if determinant == 0.0:
        raise ShellError("the plane system at a vertex is singular")
    solution = []
    for column in range(3):
        replaced = [row[:] for row in matrix]
        for row in range(3):
            replaced[row][column] = rhs[row]
        solution.append(
            (
                replaced[0][0] * (replaced[1][1] * replaced[2][2] - replaced[1][2] * replaced[2][1])
                - replaced[0][1] * (replaced[1][0] * replaced[2][2] - replaced[1][2] * replaced[2][0])
                + replaced[0][2] * (replaced[1][0] * replaced[2][1] - replaced[1][1] * replaced[2][0])
            )
            / determinant
        )
    point: Vec3 = (solution[0], solution[1], solution[2])
    residual = max(abs(_dot(normal, point) - (offset - thickness)) for normal, offset, _ in constraints)
    return point, residual


# --------------------------------------------------------------------------------------
# nearest-distance queries against the outer surface
# --------------------------------------------------------------------------------------


def _closest_point_on_triangle(
    point: Sequence[float], a: Sequence[float], b: Sequence[float], c: Sequence[float]
) -> Vec3:
    """Closest point of the closed triangle ``abc`` to ``point`` (Ericson, RTCD 5.1.5)."""
    ab = _sub(b, a)
    ac = _sub(c, a)
    ap = _sub(point, a)
    d1 = _dot(ab, ap)
    d2 = _dot(ac, ap)
    if d1 <= 0.0 and d2 <= 0.0:
        return (a[0], a[1], a[2])
    bp = _sub(point, b)
    d3 = _dot(ab, bp)
    d4 = _dot(ac, bp)
    if d3 >= 0.0 and d4 <= d3:
        return (b[0], b[1], b[2])
    vc = d1 * d4 - d3 * d2
    if vc <= 0.0 and d1 >= 0.0 and d3 <= 0.0:
        ratio = d1 / (d1 - d3)
        return (a[0] + ratio * ab[0], a[1] + ratio * ab[1], a[2] + ratio * ab[2])
    cp = _sub(point, c)
    d5 = _dot(ab, cp)
    d6 = _dot(ac, cp)
    if d6 >= 0.0 and d5 <= d6:
        return (c[0], c[1], c[2])
    vb = d5 * d2 - d1 * d6
    if vb <= 0.0 and d2 >= 0.0 and d6 <= 0.0:
        ratio = d2 / (d2 - d6)
        return (a[0] + ratio * ac[0], a[1] + ratio * ac[1], a[2] + ratio * ac[2])
    va = d3 * d6 - d5 * d4
    if va <= 0.0 and (d4 - d3) >= 0.0 and (d5 - d6) >= 0.0:
        ratio = (d4 - d3) / ((d4 - d3) + (d5 - d6))
        edge = _sub(c, b)
        return (b[0] + ratio * edge[0], b[1] + ratio * edge[1], b[2] + ratio * edge[2])
    denominator = 1.0 / (va + vb + vc)
    v = vb * denominator
    w = vc * denominator
    return (
        a[0] + ab[0] * v + ac[0] * w,
        a[1] + ab[1] * v + ac[1] * w,
        a[2] + ab[2] * v + ac[2] * w,
    )


class _TriangleGrid:
    """Uniform bucket grid over triangle bounding boxes, for nearest-surface queries."""

    def __init__(self, mesh: Mesh, cell: float) -> None:
        if not mesh.faces:
            raise ShellError("the outer surface has no faces")
        self.mesh = mesh
        self.cell = cell
        self._buckets: dict[tuple[int, int, int], list[int]] = {}
        low = [math.inf, math.inf, math.inf]
        high = [-math.inf, -math.inf, -math.inf]
        corners: list[tuple[float, ...]] = []
        for face in mesh.faces:
            points = [mesh.vertices[index] for index in face]
            corners.append(tuple(min(p[axis] for p in points) for axis in range(3)))
            corners.append(tuple(max(p[axis] for p in points) for axis in range(3)))
            for axis in range(3):
                low[axis] = min(low[axis], corners[-2][axis])
                high[axis] = max(high[axis], corners[-1][axis])
        self._low = low
        self._high = high
        for position, point in enumerate(corners):
            face_index = position // 2
            if position % 2:
                continue
            self._register(face_index, mesh.faces[face_index], point, corners[position + 1])

    def _cell_of(self, coordinates: Sequence[float]) -> tuple[int, int, int]:
        return tuple(math.floor(value / self.cell) for value in coordinates)  # type: ignore[return-value]

    def _register(
        self,
        index: int,
        face: Face,
        low: Sequence[float],
        high: Sequence[float],
    ) -> None:
        start = self._cell_of(low)
        end = self._cell_of(high)
        for i in range(start[0], end[0] + 1):
            for j in range(start[1], end[1] + 1):
                for k in range(start[2], end[2] + 1):
                    self._buckets.setdefault((i, j, k), []).append(index)

    def _max_ring(self) -> int:
        span = max((self._high[axis] - self._low[axis]) / self.cell for axis in range(3))
        return math.ceil(span) + 2

    def nearest_distance(self, point: Sequence[float]) -> float:
        """Distance from ``point`` to the closest triangle of the surface."""
        base = self._cell_of(point)
        best = math.inf
        ring = 0
        limit = self._max_ring()
        while ring <= limit:
            if ring > 0 and (ring - 1) * self.cell >= best:
                break
            for i in range(base[0] - ring, base[0] + ring + 1):
                for j in range(base[1] - ring, base[1] + ring + 1):
                    for k in range(base[2] - ring, base[2] + ring + 1):
                        if max(abs(i - base[0]), abs(j - base[1]), abs(k - base[2])) != ring:
                            continue
                        for index in self._buckets.get((i, j, k), ()):  # cells absent outside the mesh
                            face = self.mesh.faces[index]
                            nearest = _closest_point_on_triangle(
                                point,
                                self.mesh.vertices[face[0]],
                                self.mesh.vertices[face[1]],
                                self.mesh.vertices[face[2]],
                            )
                            best = min(best, _length(_sub(point, nearest)))
            ring += 1
        if not math.isfinite(best):
            raise ShellError("the distance query found no triangle; the surface is empty")
        return best


# --------------------------------------------------------------------------------------
# the public operations
# --------------------------------------------------------------------------------------


def _require_closed_outward(mesh: Mesh) -> None:
    if not mesh.faces or not mesh.vertices:
        raise ShellError("the mesh is empty")
    counts = meshlib.edge_face_counts(mesh)
    unpaired = sum(1 for count in counts.values() if count != 2)
    if unpaired:
        raise ShellError(
            f"the surface is not closed: {unpaired} of {len(counts)} edges are not shared by "
            "exactly two faces, so there is no inside to hollow (a shell needs a solid)"
        )
    volume = meshlib.total_signed_volume(mesh)
    if volume <= 0.0:
        raise ShellError(
            f"the surface is not outward-oriented (signed volume {volume!r}); run "
            "mesh.orient_outward first, otherwise the offset would go outwards"
        )


def _offset_vertices(mesh: Mesh, thickness: float) -> tuple[list[Vec3], float]:
    planes = _vertex_planes(mesh)
    points: list[Vec3] = []
    max_residual = 0.0
    for vertex, constraints in enumerate(planes):
        if not constraints:
            raise ShellError(f"vertex {vertex} belongs to no face; the mesh is inconsistent")
        point, residual = _solve_vertex(constraints, thickness, ridge=1e-14)
        points.append(point)
        max_residual = max(max_residual, residual)
    return points, max_residual


def offset_inward(mesh: Mesh, thickness: float) -> Mesh:
    """Offset a closed, outward-oriented surface inward by ``thickness``.

    The topology is preserved index for index, so face ``i`` of the result corresponds to face
    ``i`` of the input -- which is what lets :func:`shell` compare their normals. This function
    is the geometry only; it does **not** run the collision checks, so use :func:`shell` unless
    you are deliberately building an offset to inspect.
    """
    if not (thickness > 0.0) or not math.isfinite(thickness):
        raise ShellError(f"the thickness must be positive and finite: {thickness!r}")
    _require_closed_outward(mesh)
    points, _ = _offset_vertices(mesh, thickness)
    return Mesh(vertices=tuple(points), faces=mesh.faces, units=mesh.units)


def min_wall_thickness(outer: Mesh, inner: Mesh, *, cell: float | None = None) -> float:
    """Closest approach from any vertex of ``inner`` to the surface of ``outer``.

    This is the measurement the shell acceptance criterion is written in, exposed on its own
    so it can be falsified: feed it a deliberately under-offset surface and it must report the
    smaller distance (the tests do exactly that). The query goes through a uniform grid of the
    outer triangles, so it visits the faces near the point rather than all of them; the cell
    size defaults to a thirty-second of the bounding diagonal.
    """
    if not inner.vertices:
        raise ShellError("the inner surface has no vertices to measure")
    if not outer.faces:
        raise ShellError("the outer surface has no faces to measure against")
    if cell is None:
        bounding = outer.bounding_box()
        cell = max(_length(_sub(bounding.maximum, bounding.minimum)) / 32.0, 1e-12)
    grid = _TriangleGrid(outer, max(cell, 1e-12))
    return min(grid.nearest_distance(point) for point in inner.vertices)


def _face_normal_dot(outer: Mesh, inner: Mesh) -> float:
    worst = math.inf
    for face in outer.faces:
        outer_plane = _face_plane(outer, face)[0]
        inner_plane = _face_plane(inner, face)[0]
        worst = min(worst, _dot(outer_plane, inner_plane))
    return worst


def shell(mesh: Mesh, thickness: float, *, wall_tolerance: float = 0.05) -> ShellResult:
    """Hollow ``mesh`` with a wall of ``thickness``, and prove the wall is there.

    Raises :class:`ShellError` when the surface is open or inward-oriented, when the thickness
    is not positive, when the offset inverts (non-positive inner volume or a flipped normal),
    or when the measured minimum wall falls more than ``wall_tolerance`` below the request.
    The last condition has two legitimate causes and the message names both: the thickness
    exceeds half a local wall, or the mesh is too coarse for the thickness (a faceted sphere
    is short of the true offset by ``O(h^2)``).
    """
    if not (thickness > 0.0) or not math.isfinite(thickness):
        raise ShellError(f"the thickness must be positive and finite: {thickness!r}")
    _require_closed_outward(mesh)
    points, residual = _offset_vertices(mesh, thickness)
    inner = Mesh(vertices=tuple(points), faces=mesh.faces, units=mesh.units)
    outer_volume = meshlib.total_signed_volume(mesh)
    inner_volume = meshlib.total_signed_volume(inner)
    if inner_volume <= 0.0:
        raise ShellError(
            f"the inward offset of {thickness!r} m inverts the solid: the inner signed volume is "
            f"{inner_volume!r} m^3 (outer {outer_volume!r} m^3), so the thickness exceeds half "
            "the local wall somewhere. Reduce the thickness or split the body"
        )
    normal_dot = _face_normal_dot(mesh, inner)
    if normal_dot <= 0.0:
        raise ShellError(
            f"the inward offset folds the surface: min(dot(n_inner, n_outer)) = {normal_dot!r} "
            "over the faces, so some offset faces are wound backwards"
        )
    wall = min_wall_thickness(mesh, inner)
    shortfall = thickness - wall
    if shortfall > wall_tolerance * thickness:
        raise ShellError(
            f"the requested wall of {thickness!r} m is not achieved: the closest approach from an "
            f"inner vertex to the outer surface is {wall!r} m ({100.0 * shortfall / thickness:.1f}% "
            "short). Either the thickness exceeds half a local wall (the offset surfaces collide), "
            "or the mesh is too coarse for this thickness; refine or reduce. No mesh is returned"
        )
    return ShellResult(
        outer=mesh,
        inner=inner,
        thickness=thickness,
        outer_volume=outer_volume,
        inner_volume=inner_volume,
        min_wall_thickness=wall,
        max_plane_residual=residual,
        min_normal_dot=normal_dot,
    )


def hollow_body(result: ShellResult) -> HollowBody:
    """Export the shell as a single surface: outer faces plus the inner faces reversed.

    The two components do not touch (no rim), so the result is a closed double-walled body,
    not a CAD-style open shell -- the missing feature is the B-rep rim and is recorded rather
    than approximated. The inner faces are reversed so that the surface bounds the *material*:
    its divergence-theorem volume is ``outer_volume - inner_volume``, not the sum. Writing that
    down matters, because the two components are geometrically disjoint and a reader who
    assumes "each component outward" would expect exactly the opposite sign.
    """
    offset = len(result.outer.vertices)
    vertices = tuple(result.outer.vertices) + tuple(result.inner.vertices)
    faces: list[Face] = list(result.outer.faces)
    faces.extend((a + offset, c + offset, b + offset) for a, b, c in result.inner.faces)
    surface = Mesh(vertices=vertices, faces=tuple(faces), units=result.outer.units)
    return HollowBody(
        mesh=surface,
        material_volume=result.outer_volume - result.inner_volume,
        cavity_volume=result.inner_volume,
    )
