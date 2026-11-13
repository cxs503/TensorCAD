"""Regularised boolean operations on closed triangle meshes (union/intersection/difference).

Why this module exists (GAP-CAD-003, ``docs/oss_inventory_cad/GAP-ANCHORS.md``)
-------------------------------------------------------------------------------
:mod:`tensorcad.csg` implements union/intersection/difference as **implicit-field
operations**: the operands are sampled on a grid and the result is a voxel field.  That is
the right tool for sculpting a hull envelope, but the wrong tool for cutting a hatch
opening in a deck: the answer is a *surface* with the topology of the cut, and a field
answer has to be re-meshed afterwards, which loses the part boundary wherever the grid
does not line up with it.

This module works on the other representation: closed edge-manifold triangle mesh in,
closed edge-manifold triangle mesh out.

Algorithm (surface classification, not BSP-leaf clipping)
--------------------------------------------------------
The classical BSP boolean of Naylor/Thibault (the ``csg.js`` formulation) classifies each
polygon against a node plane and resolves *coplanar* polygons by the direction they face.
That rule is exact when two coplanar faces are identical or disjoint, and wrong when they
partially overlap -- which is precisely the case of a hole cut flush with a face, i.e.
most deck/hatch problems.  Measured on a cube minus a coaxial through-bolt, the ``csg.js``
pipeline returns a surface with 300 boundary edges, so it is not used here.

Instead, each operand is decomposed into the cells of the *other* operand's
supporting-plane arrangement, and every cell is classified by two offset probes of a
winding-number point-in-solid test:

1. collect the distinct supporting planes of B (triangle planes, deduplicated at the
   model-relative tolerance);
2. split every triangle of A by those planes, keeping both sides.  Each resulting convex
   cell therefore lies entirely inside B, entirely outside B, or entirely *on* B's
   boundary -- the ambiguity that breaks BSP clipping cannot arise;
3. probe the cell centroid at ``+/- delta`` along the cell normal.  Four classes follow:
   outside, inside, and two boundary classes (material behind the cell vs. material in
   front of it);
4. keep the classes that belong to the requested set operation, reverse the second
   operand's interior surface for a difference, weld, and validate.

Contract
--------
The output is **validated, not assumed**: after welding it must be watertight,
edge-manifold, free of zero-area faces and outward-oriented with positive volume.
Otherwise :class:`MeshBooleanError` is raised with the offending counts -- an explicit
failure beats an invalid solid that a volume mesher will silently invert.

Tolerances are relative to the model scale (largest extent of the combined bounding box),
so a part modelled in millimetres and the same part in metres produce the same topology.
The probe distance :data:`PROBE_DISTANCE` is a *contractual* limit, stated rather than
hidden: features thinner than ``1e-6``x the model scale next to a boolean seam are outside
the guarantee, because no finite-precision predicate can separate a cell from the surface
it lies on.

Limitations stated rather than hidden: this is a *polyhedral* boolean.  The result is
exact up to the input tessellation, so the volume error against an analytic solid is the
tessellation error of the operands (measured in ``docs/mesh_boolean.md``), not a
boolean-specific error.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from .geometry import Face, Mesh, Vec3
from .mesh import (
    edge_face_counts,
    is_edge_manifold,
    is_watertight,
    total_signed_volume,
    zero_area_faces,
)

#: Distance tolerance, relative to the model scale, for point/plane classification.
RELATIVE_TOLERANCE = 1e-9

#: Offset used by the inside/outside probes, relative to the model scale.
PROBE_DISTANCE = 1e-6

#: Relative area below which an output triangle is considered a sliver and dropped.
AREA_TOLERANCE = 1e-14

OPERATIONS = ("union", "intersection", "difference")

# cell classes with respect to the other operand
_OUTSIDE = "outside"
_INSIDE = "inside"
_BEHIND = "behind"  # material behind the cell: same outward direction as the boundary
_FRONT = "front"  # material in front of the cell: opposite outward direction

#: which cells each operation keeps (see ``docs/mesh_boolean.md`` for the derivation)
_KEEP = {
    "union": {"first": (_OUTSIDE, _BEHIND), "second": (_OUTSIDE,), "flip_second": False},
    "intersection": {"first": (_INSIDE, _BEHIND), "second": (_INSIDE,), "flip_second": False},
    "difference": {"first": (_OUTSIDE, _FRONT), "second": (_INSIDE,), "flip_second": True},
}


class MeshBooleanError(ValueError):
    """Raised when a boolean cannot produce a valid closed manifold."""


# --------------------------------------------------------------------------- #
# vector helpers
# --------------------------------------------------------------------------- #
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


def _lerp(a: Vec3, b: Vec3, t: float) -> Vec3:
    return (a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t)


def _norm(a: Vec3) -> float:
    return math.sqrt(_dot(a, a))


# --------------------------------------------------------------------------- #
# planes and convex cells
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class _Plane:
    normal: Vec3
    offset: float  # {x : dot(normal, x) == offset}

    @staticmethod
    def through(a: Vec3, b: Vec3, c: Vec3) -> "_Plane | None":
        normal = _cross(_sub(b, a), _sub(c, a))
        length = _norm(normal)
        if length == 0.0:
            return None
        unit = (normal[0] / length, normal[1] / length, normal[2] / length)
        return _Plane(unit, _dot(unit, a))

    def distance(self, point: Vec3) -> float:
        return _dot(self.normal, point) - self.offset

    def key(self, eps: float) -> tuple[float, float, float, float]:
        # Planes agreeing within eps are the same supporting plane.  The normal sign is
        # kept: a face and the reverse side of a thin wall are different supporting planes
        # for the arrangement, and merging them would erase a cell boundary.
        return (
            round(self.normal[0] / eps),
            round(self.normal[1] / eps),
            round(self.normal[2] / eps),
            round(self.offset / eps),
        )


@dataclass(frozen=True)
class _Cell:
    vertices: tuple[Vec3, ...]
    normal: Vec3

    def centroid(self) -> Vec3:
        count = len(self.vertices)
        return (
            sum(vertex[0] for vertex in self.vertices) / count,
            sum(vertex[1] for vertex in self.vertices) / count,
            sum(vertex[2] for vertex in self.vertices) / count,
        )

    def flipped(self) -> "_Cell":
        return _Cell(tuple(reversed(self.vertices)), self.normal)


def _split_cell(cell: _Cell, plane: _Plane, eps: float) -> tuple[list[_Cell], list[_Cell]]:
    """Split ``cell`` by ``plane`` into (front, back); a coplanar cell returns ([], [])."""

    distances = [plane.distance(vertex) for vertex in cell.vertices]
    if all(distance >= -eps for distance in distances):
        # A cell lying *in* the cutting plane (all distances within eps) is kept, on the
        # front side.  Dropping it would delete the very surface the boolean must decide
        # about -- the flush-cut face -- and the offset probes classify it correctly
        # afterwards whoever's side of the plane it is filed under.
        return [cell], []
    if all(distance <= eps for distance in distances):
        return [], [cell]

    front: list[Vec3] = []
    back: list[Vec3] = []
    count = len(cell.vertices)
    for index in range(count):
        nxt = (index + 1) % count
        current, following = cell.vertices[index], cell.vertices[nxt]
        current_d, next_d = distances[index], distances[nxt]
        if current_d >= -eps:
            front.append(current)
        if current_d <= eps:
            back.append(current)
        if (current_d > eps and next_d < -eps) or (current_d < -eps and next_d > eps):
            t = current_d / (current_d - next_d)
            crossing = _lerp(current, following, t)
            front.append(crossing)
            back.append(crossing)

    out_front = [_Cell(tuple(front), cell.normal)] if len(front) >= 3 else []
    out_back = [_Cell(tuple(back), cell.normal)] if len(back) >= 3 else []
    return out_front, out_back


def _supporting_planes(mesh: Mesh, eps: float) -> list[_Plane]:
    planes: dict[tuple[float, float, float, float], _Plane] = {}
    for face in mesh.faces:
        a, b, c = (mesh.vertices[vertex] for vertex in face)
        plane = _Plane.through(a, b, c)
        if plane is None:
            raise MeshBooleanError(
                "operand has a degenerate (zero-area) triangle; clean it first with "
                "tensorcad.mesh.remove_degenerate_faces"
            )
        planes.setdefault(plane.key(eps), plane)
    return list(planes.values())


def _cells(mesh: Mesh, planes: list[_Plane], eps: float) -> list[_Cell]:
    cells: list[_Cell] = []
    for face in mesh.faces:
        a, b, c = (mesh.vertices[vertex] for vertex in face)
        plane = _Plane.through(a, b, c)
        if plane is None:
            continue
        queue = [_Cell((a, b, c), plane.normal)]
        for cutter in planes:
            if not queue:
                break
            front: list[_Cell] = []
            back: list[_Cell] = []
            for cell in queue:
                f, b_ = _split_cell(cell, cutter, eps)
                front.extend(f)
                back.extend(b_)
            queue = front + back
        cells.extend(queue)
    return cells


# --------------------------------------------------------------------------- #
# point-in-solid by winding number
# --------------------------------------------------------------------------- #
def _solid_angle(a: Vec3, b: Vec3, c: Vec3, p: Vec3) -> float:
    ax, ay, az = a[0] - p[0], a[1] - p[1], a[2] - p[2]
    bx, by, bz = b[0] - p[0], b[1] - p[1], b[2] - p[2]
    cx, cy, cz = c[0] - p[0], c[1] - p[1], c[2] - p[2]
    la = math.sqrt(ax * ax + ay * ay + az * az)
    lb = math.sqrt(bx * bx + by * by + bz * bz)
    lc = math.sqrt(cx * cx + cy * cy + cz * cz)
    if la == 0.0 or lb == 0.0 or lc == 0.0:
        return 0.0
    determinant = ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx)
    denominator = (
        la * lb * lc
        + (ax * bx + ay * by + az * bz) * lc
        + (bx * cx + by * cy + bz * cz) * la
        + (cx * ax + cy * ay + cz * az) * lb
    )
    return 2.0 * math.atan2(determinant, denominator)


def point_in_solid(mesh: Mesh, point: Vec3) -> bool:
    """True when ``point`` is inside the closed ``mesh``.

    The predicate is the winding number (sum of signed solid angles over ``4 pi``), so no
    ray has to be aimed -- there is no direction to be unlucky with -- and the decision is
    the sign of the sum, which is stable for points as close to the surface as the probe
    distance.
    """

    total = 0.0
    for face in mesh.faces:
        a, b, c = (mesh.vertices[vertex] for vertex in face)
        total += _solid_angle(a, b, c, point)
    return abs(total) > 2.0 * math.pi


def _classify(cell: _Cell, other: Mesh, delta: float) -> str:
    centroid = cell.centroid()
    normal = cell.normal
    plus = (
        centroid[0] + normal[0] * delta,
        centroid[1] + normal[1] * delta,
        centroid[2] + normal[2] * delta,
    )
    minus = (
        centroid[0] - normal[0] * delta,
        centroid[1] - normal[1] * delta,
        centroid[2] - normal[2] * delta,
    )
    in_plus = point_in_solid(other, plus)
    in_minus = point_in_solid(other, minus)
    if in_plus and in_minus:
        return _INSIDE
    if not in_plus and not in_minus:
        return _OUTSIDE
    return _BEHIND if in_minus else _FRONT


# --------------------------------------------------------------------------- #
# mesh helpers
# --------------------------------------------------------------------------- #
def _model_scale(*meshes: Mesh) -> float:
    lo = [math.inf, math.inf, math.inf]
    hi = [-math.inf, -math.inf, -math.inf]
    for mesh in meshes:
        for vertex in mesh.vertices:
            for axis in range(3):
                lo[axis] = min(lo[axis], vertex[axis])
                hi[axis] = max(hi[axis], vertex[axis])
    extent = max(hi[axis] - lo[axis] for axis in range(3))
    if not math.isfinite(extent) or extent <= 0.0:
        raise MeshBooleanError("operands have a degenerate (zero-extent) bounding box")
    return extent


def _check_operand(mesh: Mesh, name: str) -> None:
    if not mesh.faces:
        raise MeshBooleanError(f"{name} has no faces")
    if not is_watertight(mesh):
        raise MeshBooleanError(f"{name} is not closed (watertight); a boolean needs a solid")
    if not is_edge_manifold(mesh):
        raise MeshBooleanError(f"{name} is not edge-manifold; a boolean needs a solid")
    if zero_area_faces(mesh):
        raise MeshBooleanError(
            f"{name} has zero-area faces; run tensorcad.mesh.remove_degenerate_faces first"
        )
    if total_signed_volume(mesh) <= 0.0:
        raise MeshBooleanError(
            f"{name} has non-positive enclosed volume: its faces are not consistently "
            "outward-oriented (run tensorcad.mesh.orient_outward)"
        )


def _cells_to_mesh(cells: list[_Cell], eps: float, area_eps: float) -> Mesh:
    vertices: list[Vec3] = []
    index_of: dict[tuple[int, int, int], int] = {}
    faces: list[Face] = []
    seen: set[tuple[int, int, int]] = set()

    def vertex_id(point: Vec3) -> int:
        key = (round(point[0] / eps), round(point[1] / eps), round(point[2] / eps))
        found = index_of.get(key)
        if found is None:
            found = len(vertices)
            index_of[key] = found
            vertices.append(point)
        return found

    for cell in cells:
        ids: list[int] = []
        for vertex in cell.vertices:
            identifier = vertex_id(vertex)
            if not ids or ids[-1] != identifier:
                ids.append(identifier)
        if len(ids) > 2 and ids[0] == ids[-1]:
            ids.pop()
        if len(ids) < 3:
            continue
        for index in range(1, len(ids) - 1):
            a, b, c = ids[0], ids[index], ids[index + 1]
            if len({a, b, c}) < 3:
                continue
            a_pos, b_pos, c_pos = vertices[a], vertices[b], vertices[c]
            cross = _cross(_sub(b_pos, a_pos), _sub(c_pos, a_pos))
            if 0.5 * _norm(cross) <= area_eps:
                continue
            key = tuple(sorted((a, b, c)))
            if key in seen:
                continue
            seen.add(key)
            faces.append((a, b, c))
    return Mesh(tuple(vertices), tuple(faces))


def _stitch_t_junctions(
    mesh: Mesh, eps: float, area_eps: float = 0.0, max_passes: int = 8
) -> Mesh:
    """Split faces whose edge has another mesh vertex lying in its interior.

    Two operands are triangulated independently, so a seam can arrive with a T-junction:
    the cut through a side quad of a faceted cylinder introduces a vertex at the quad's
    *diagonal*, which the flat cut face on the other operand knows nothing about.  Welding
    cannot fix that -- the edges are collinear, not coincident -- so the face on the side
    that owns the long edge is split at the intruding vertex.  Left alone, those are
    exactly the boundary edges that make a boolean look successful and be unusable.
    """

    for _ in range(max_passes):
        counts = edge_face_counts(mesh)
        boundary = [edge for edge, count in counts.items() if count == 1]
        if not boundary:
            return mesh
        vertices = mesh.vertices
        splits: dict[tuple[int, int], list[int]] = {}
        for a, b in boundary:
            pa, pb = vertices[a], vertices[b]
            direction = _sub(pb, pa)
            length_sq = _dot(direction, direction)
            if length_sq <= 0.0:
                continue
            hits: list[tuple[float, int]] = []
            for index, point in enumerate(vertices):
                if index == a or index == b:
                    continue
                t = _dot(_sub(point, pa), direction) / length_sq
                if not 0.0 < t < 1.0:
                    continue
                projected = _lerp(pa, pb, t)
                if _norm(_sub(projected, point)) <= eps:
                    hits.append((t, index))
            if hits:
                splits[(min(a, b), max(a, b))] = [index for _, index in sorted(hits)]
        if not splits:
            return mesh

        faces: list[Face] = []
        for face in mesh.faces:
            polygon: list[int] = []
            count = len(face)
            for index in range(count):
                start, end = face[index], face[(index + 1) % count]
                polygon.append(start)
                key = (min(start, end), max(start, end))
                extra = splits.get(key)
                if extra:
                    forward = start < end
                    polygon.extend(extra if forward else list(reversed(extra)))
            if polygon and polygon[0] == polygon[-1]:
                polygon.pop()
            # a triangle whose edges were subdivided is still convex: fan it
            for index in range(1, len(polygon) - 1):
                face_out = (polygon[0], polygon[index], polygon[index + 1])
                if len(set(face_out)) < 3:
                    continue
                a, b, c = (vertices[i] for i in face_out)
                if 0.5 * _norm(_cross(_sub(b, a), _sub(c, a))) <= area_eps:
                    continue  # a stitch can collapse a sliver: drop it, never keep it
                faces.append(face_out)
        mesh = Mesh(vertices, tuple(faces))
    return mesh


def _validate(mesh: Mesh, operation: str) -> Mesh:
    if not mesh.faces:
        raise MeshBooleanError(f"{operation}: the result is empty")
    counts = edge_face_counts(mesh)
    boundary = sum(1 for count in counts.values() if count == 1)
    non_manifold = sum(1 for count in counts.values() if count > 2)
    slivers = len(zero_area_faces(mesh))
    if boundary or non_manifold or slivers or not is_watertight(mesh):
        raise MeshBooleanError(
            f"{operation}: result is not a closed edge-manifold surface "
            f"(faces={len(mesh.faces)}, boundary_edges={boundary}, "
            f"non_manifold_edges={non_manifold}, zero_area_faces={slivers}); "
            "refusing to return an invalid solid"
        )
    volume = total_signed_volume(mesh)
    if volume <= 0.0:
        raise MeshBooleanError(f"{operation}: result is inward-oriented (volume {volume!r})")
    return mesh


# --------------------------------------------------------------------------- #
# public API
# --------------------------------------------------------------------------- #
def boolean(
    first: Mesh,
    second: Mesh,
    operation: str,
    *,
    tolerance: float | None = None,
    probe: float | None = None,
) -> Mesh:
    """Return ``first <operation> second`` as a validated closed manifold mesh.

    ``operation`` is one of ``"union"``, ``"intersection"``, ``"difference"``.
    """

    if operation not in OPERATIONS:
        raise MeshBooleanError(f"unknown operation {operation!r}; expected one of {OPERATIONS}")
    _check_operand(first, "first operand")
    _check_operand(second, "second operand")

    scale = _model_scale(first, second)
    eps = (RELATIVE_TOLERANCE if tolerance is None else tolerance) * scale
    delta = (PROBE_DISTANCE if probe is None else probe) * scale
    area_eps = AREA_TOLERANCE * scale * scale
    rule = _KEEP[operation]

    cells_first = _cells(first, _supporting_planes(second, eps), eps)
    cells_second = _cells(second, _supporting_planes(first, eps), eps)

    kept: list[_Cell] = []
    for cell in cells_first:
        if _classify(cell, second, delta) in rule["first"]:
            kept.append(cell)
    for cell in cells_second:
        if _classify(cell, first, delta) in rule["second"]:
            kept.append(cell.flipped() if rule["flip_second"] else cell)

    assembled = _cells_to_mesh(kept, eps, area_eps)
    return _validate(_stitch_t_junctions(assembled, eps, area_eps), operation)


def mesh_union(first: Mesh, second: Mesh, *, tolerance: float | None = None) -> Mesh:
    """Boolean union ``first | second`` as a validated closed manifold."""

    return boolean(first, second, "union", tolerance=tolerance)


def mesh_intersection(first: Mesh, second: Mesh, *, tolerance: float | None = None) -> Mesh:
    """Boolean intersection ``first & second`` as a validated closed manifold."""

    return boolean(first, second, "intersection", tolerance=tolerance)


def mesh_difference(first: Mesh, second: Mesh, *, tolerance: float | None = None) -> Mesh:
    """Boolean difference ``first - second`` as a validated closed manifold."""

    return boolean(first, second, "difference", tolerance=tolerance)


__all__ = [
    "AREA_TOLERANCE",
    "MeshBooleanError",
    "OPERATIONS",
    "PROBE_DISTANCE",
    "RELATIVE_TOLERANCE",
    "boolean",
    "mesh_difference",
    "mesh_intersection",
    "mesh_union",
    "point_in_solid",
]