"""Repair of an open, non-manifold or inconsistently wound triangulation.

``docs/oss_inventory_cad/GAP-ANCHORS.md`` records GAP-CAD-009 as *partial*: the package
can weld, drop degenerate facets and *detect* self intersections, but after a defect is
detected there is no path back to a mesh a solver can use -- no hole filling, no
manifoldisation, no orientation repair, and no post-repair proof that the self-intersection
set is empty. This module is that path, and it is written the way the rest of the package
reports geometry: every step is measured, and a repair that would trade one defect for
another is **refused** with the reason attached instead of being returned as a success.

The steps, each with what it can and cannot do:

============================  ====================================================  ==========
step                          what it does                                          refuses
============================  ====================================================  ==========
``weld``                      merges vertices closer than the weld tolerance        --
``drop_degenerate_faces``     deletes zero-area facets                              --
``drop_duplicate_faces``      keeps one copy of coincident facets                   --
``manifoldize``               splits a vertex per face sheet, so a bow-tie edge     never (it
                              or a bow-tie vertex becomes two manifold sheets       only splits)
``orient_consistently``       floods a winding across manifold edges, then turns    a non-
                              each closed sheet outward                             orientable
                                                                                    sheet
``fill_holes``                closes each simple boundary loop with an ear-clipped  a loop that
                              planar fill over the loop's own vertices              is not a
                                                                                    simple loop
============================  ====================================================  ==========

The fill uses **the loop's existing vertices only**: no point is invented, so the patch
shares its boundary edges with the shell exactly and closure is a topological statement
rather than a tolerance. What that costs is shape: a hole in a curved surface is closed by a
planar patch, and the deviation from the plane is measured per hole
(``HoleFill.holes[i]["planar_deviation_m"]``), not hidden.

What the module deliberately does **not** do:

* It does not remove self intersections. That is a remeshing problem (fTetWild-class); the
  honest answer here is a refusal, and :func:`repair_mesh` returns the *original* mesh with
  ``accepted=False`` and the offending face pairs in the reason.
* It does not re-tessellate, so it cannot improve a facet that is simply too coarse.
* It does not equate edge-manifoldness with intersection-freeness: two sheets passing within
  the weld tolerance of each other are topologically separate after :func:`manifoldize`,
  which is what "manifold" means. The self-intersection count is the separate measurement
  that answers whether they still pierce each other, and it is part of the acceptance test.
"""
from __future__ import annotations

import math
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Any, Mapping, Sequence

from . import mesh as meshlib
from .geometry import Face, Mesh, Vec3
from .intersect import find_self_intersections

__all__ = [
    "DEFAULT_WELD_TOLERANCE",
    "DEFAULT_VOLUME_TOLERANCE",
    "DEFAULT_AREA_TOLERANCE",
    "DEFAULT_MAX_HOLE_EDGES",
    "MeshRepairError",
    "MeshMetrics",
    "HoleFill",
    "RepairStep",
    "RepairResult",
    "measure_mesh",
    "fill_holes",
    "boundary_loops",
    "manifoldize",
    "orient_consistently",
    "drop_duplicate_faces",
    "repair_mesh",
]

#: Vertices closer than this are the same point. Same default as ``geomquality``.
DEFAULT_WELD_TOLERANCE = 1e-6
#: A repair may not move the enclosed volume by more than this share of the reference
#: volume (``GAP-CAD-009`` asks for < 1 %).
DEFAULT_VOLUME_TOLERANCE = 0.01
#: ... nor the wetted area (a fill adds the patch's own area, so a reference is required).
DEFAULT_AREA_TOLERANCE = 0.01
#: A boundary loop longer than this many edges is not closed by a planar fill here: beyond
#: it the honest route is a remesher, not ear clipping.
DEFAULT_MAX_HOLE_EDGES = 512


class MeshRepairError(ValueError):
    """Raised when a repair step cannot be carried out honestly."""


# ---------------------------------------------------------------------------------------
# Measurement: the numbers before and after
# ---------------------------------------------------------------------------------------


@dataclass(frozen=True)
class MeshMetrics:
    """Everything the acceptance test reads, in one place."""

    vertices: int
    faces: int
    open_edges: int
    nonmanifold_edges: int
    components: int
    duplicate_faces: int
    winding_conflicts: int
    self_intersections: int
    signed_volume: float
    surface_area: float

    @property
    def closed(self) -> bool:
        return self.faces > 0 and self.open_edges == 0

    @property
    def manifold(self) -> bool:
        return self.nonmanifold_edges == 0

    @property
    def consistently_wound(self) -> bool:
        return self.winding_conflicts == 0

    @property
    def clean(self) -> bool:
        """Closed, edge-manifold, consistently wound, no duplicates, no intersections."""
        return (
            self.closed
            and self.manifold
            and self.consistently_wound
            and self.duplicate_faces == 0
            and self.self_intersections == 0
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "vertices": self.vertices,
            "faces": self.faces,
            "open_edges": self.open_edges,
            "nonmanifold_edges": self.nonmanifold_edges,
            "components": self.components,
            "duplicate_faces": self.duplicate_faces,
            "winding_conflicts": self.winding_conflicts,
            "self_intersections": self.self_intersections,
            "signed_volume_m3": self.signed_volume,
            "surface_area_m2": self.surface_area,
            "closed": self.closed,
            "manifold": self.manifold,
            "consistently_wound": self.consistently_wound,
            "clean": self.clean,
        }


def _directed_edges(mesh: Mesh) -> dict[tuple[int, int], list[tuple[int, int]]]:
    """Undirected edge -> ``(face index, direction)``.

    Direction is ``+1`` when the face traverses the edge from the smaller vertex id to the
    larger, ``-1`` otherwise.
    """
    out: dict[tuple[int, int], list[tuple[int, int]]] = defaultdict(list)
    for index, face in enumerate(mesh.faces):
        for k in range(3):
            a, b = face[k], face[(k + 1) % 3]
            out[meshlib.undirected_edge((a, b))].append((index, 1 if a < b else -1))
    return out


def _face_components(mesh: Mesh) -> tuple[tuple[int, ...], ...]:
    """Faces grouped into sheets joined through edges shared by exactly two faces."""
    parent = list(range(len(mesh.faces)))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for incident in _directed_edges(mesh).values():
        if len(incident) == 2:
            ra, rb = find(incident[0][0]), find(incident[1][0])
            if ra != rb:
                parent[rb] = ra
    groups: dict[int, list[int]] = defaultdict(list)
    for index in range(len(mesh.faces)):
        groups[find(index)].append(index)
    return tuple(tuple(sorted(group)) for group in groups.values())


def _winding_conflicts(mesh: Mesh) -> int:
    """Edges whose two incident faces traverse them in the *same* direction."""
    return sum(
        1
        for incident in _directed_edges(mesh).values()
        if len(incident) == 2 and incident[0][1] == incident[1][1]
    )


def measure_mesh(mesh: Mesh, *, self_intersection_limit: int = 16) -> MeshMetrics:
    """Read the metrics the acceptance test uses. Never mutates, never raises on defects."""
    counts = meshlib.edge_face_counts(mesh)
    seen = Counter(tuple(sorted(face)) for face in mesh.faces)
    return MeshMetrics(
        vertices=len(mesh.vertices),
        faces=len(mesh.faces),
        open_edges=sum(1 for count in counts.values() if count == 1),
        nonmanifold_edges=sum(1 for count in counts.values() if count > 2),
        components=len(_face_components(mesh)),
        duplicate_faces=sum(count - 1 for count in seen.values() if count > 1),
        winding_conflicts=_winding_conflicts(mesh),
        self_intersections=len(find_self_intersections(mesh, limit=self_intersection_limit)),
        signed_volume=meshlib.total_signed_volume(mesh),
        surface_area=meshlib.surface_area(mesh),
    )


# ---------------------------------------------------------------------------------------
# Steps
# ---------------------------------------------------------------------------------------


@dataclass(frozen=True)
class RepairStep:
    """One step of a repair log: what ran, how much it changed, and why it stopped."""

    name: str
    changed: int
    detail: str

    def to_dict(self) -> dict[str, Any]:
        return {"name": self.name, "changed": self.changed, "detail": self.detail}


def _rebuild_with_faces(mesh: Mesh, faces: Sequence[Face]) -> Mesh:
    """Reuse ``mesh``'s vertex list, drop the vertices no face uses, keep the metadata."""
    used: dict[int, int] = {}
    vertices: list[Vec3] = []
    out: list[Face] = []
    for face in faces:
        remapped = []
        for vertex in face:
            if vertex not in used:
                used[vertex] = len(vertices)
                vertices.append(mesh.vertices[vertex])
            remapped.append(used[vertex])
        out.append((remapped[0], remapped[1], remapped[2]))
    return Mesh(
        vertices=tuple(vertices),
        faces=tuple(out),
        units=mesh.units,
        coordinate_frame=mesh.coordinate_frame,
    )


def drop_duplicate_faces(mesh: Mesh) -> tuple[Mesh, RepairStep]:
    """Keep one copy of every coincident facet (the audit's ``DUPLICATE_FACE``)."""
    seen: set[tuple[int, int, int]] = set()
    keep: list[Face] = []
    dropped = 0
    for face in mesh.faces:
        key = tuple(sorted(face))
        if key in seen:
            dropped += 1
            continue
        seen.add(key)
        keep.append(face)
    if dropped == 0:
        return mesh, RepairStep("drop_duplicate_faces", 0, "no coincident facets")
    return (
        _rebuild_with_faces(mesh, keep),
        RepairStep(
            "drop_duplicate_faces",
            dropped,
            f"dropped {dropped} coincident facet(s); {len(keep)} kept",
        ),
    )


def manifoldize(mesh: Mesh) -> tuple[Mesh, RepairStep]:
    """Split vertices so that every edge is shared by at most two faces.

    Faces are grouped into sheets through edges already shared by exactly two faces (the
    manifold skeleton). A vertex used by more than one sheet -- a bow-tie vertex, or one end
    of an edge with three or more incident faces -- is duplicated once per sheet. The sheets
    then touch geometrically but not topologically, which is what "edge-manifold" means;
    :func:`measure_mesh`'s self-intersection count is the separate question of whether they
    still pass through each other.
    """
    if not mesh.faces:
        raise MeshRepairError("cannot manifoldise a mesh with no faces")
    components = _face_components(mesh)
    members: dict[int, int] = {}
    for sheet, group in enumerate(components):
        for face in group:
            members[face] = sheet

    owner: dict[tuple[int, int], int] = {}
    vertices: list[Vec3] = []
    faces: list[Face] = []
    for index, face in enumerate(mesh.faces):
        sheet = members[index]
        remapped = []
        for vertex in face:
            key = (vertex, sheet)
            existing = owner.get(key)
            if existing is None:
                existing = len(vertices)
                vertices.append(mesh.vertices[vertex])
                owner[key] = existing
            remapped.append(existing)
        faces.append((remapped[0], remapped[1], remapped[2]))

    result = Mesh(
        vertices=tuple(vertices),
        faces=tuple(faces),
        units=mesh.units,
        coordinate_frame=mesh.coordinate_frame,
    )
    leftover = sum(1 for count in meshlib.edge_face_counts(result).values() if count > 2)
    if leftover:
        raise MeshRepairError(
            f"manifoldisation left {leftover} non-manifold edge(s): the bow-tie is inside "
            "one sheet and needs a face-splitting repair, not a vertex split"
        )
    duplications = len(vertices) - len({vertex for face in mesh.faces for vertex in face})
    if duplications <= 0:
        return mesh, RepairStep("manifoldize", 0, "every vertex belongs to one sheet")
    return (
        result,
        RepairStep(
            "manifoldize",
            duplications,
            f"split {len(components)} sheet(s): {duplications} vertex cop(ies) created, "
            f"{len(meshlib.edge_face_counts(mesh))} edges -> {len(meshlib.edge_face_counts(result))}",
        ),
    )


def orient_consistently(mesh: Mesh) -> tuple[Mesh, RepairStep]:
    """Flood a winding across manifold edges, then turn every closed sheet outward.

    A sheet that is not orientable (a Möbius band) has no consistent winding and is
    refused: flipping it would only move the conflict around.
    """
    if not mesh.faces:
        raise MeshRepairError("cannot orient a mesh with no faces")
    flips = [1] * len(mesh.faces)
    visited = [False] * len(mesh.faces)
    conflicts = 0
    adjacency: dict[int, list[tuple[int, bool]]] = defaultdict(list)
    for incident in _directed_edges(mesh).values():
        if len(incident) == 2:
            (first, d_first), (second, d_second) = incident
            same = d_first == d_second
            adjacency[first].append((second, same))
            adjacency[second].append((first, same))

    for start in range(len(mesh.faces)):
        if visited[start]:
            continue
        visited[start] = True
        queue = [start]
        while queue:
            current = queue.pop()
            for neighbour, same in adjacency[current]:
                required = -flips[current] if same else flips[current]
                if not visited[neighbour]:
                    visited[neighbour] = True
                    flips[neighbour] = required
                    queue.append(neighbour)
                elif flips[neighbour] != required:
                    conflicts += 1
    if conflicts:
        raise MeshRepairError(
            f"windings cannot be made consistent: {conflicts} contradiction(s) -- the sheet "
            "is not orientable, or two sheets are glued along an edge"
        )

    # Flooding from an arbitrary start can invert a whole sheet; flip the sheet back when
    # that is the smaller change, so the reported count is the number of facets actually
    # wrong rather than the number the walk happened to visit.
    for group in _face_components(mesh):
        negative = sum(1 for index in group if flips[index] < 0)
        if negative * 2 > len(group):
            for index in group:
                flips[index] = -flips[index]

    flipped_faces = [
        face if flip > 0 else (face[0], face[2], face[1])
        for face, flip in zip(mesh.faces, flips)
    ]
    flipped = sum(1 for flip in flips if flip < 0)

    # Outward, per sheet: a positive signed volume. Sheets are recomputed after the flood
    # because that is the winding the outward test has to be applied to.
    result = _rebuild_with_faces(mesh, flipped_faces)
    outward: list[Face] = []
    turned = 0
    for group in _face_components(result):
        patch = [result.faces[index] for index in group]
        volume = meshlib.total_signed_volume(Mesh(vertices=result.vertices, faces=tuple(patch)))
        if volume < 0.0:
            patch = [(face[0], face[2], face[1]) for face in patch]
            turned += len(patch)
        outward.extend(patch)
    final = _rebuild_with_faces(result, outward)
    changed = flipped + turned
    if changed == 0:
        return mesh, RepairStep(
            "orient_consistently", 0, "winding already consistent and outward"
        )
    return (
        final,
        RepairStep(
            "orient_consistently",
            changed,
            f"flipped {flipped} facet(s) to agree across shared edges, {turned} to face "
            "outward",
        ),
    )


# ---------------------------------------------------------------------------------------
# Hole filling
# ---------------------------------------------------------------------------------------


@dataclass(frozen=True)
class HoleFill:
    """The refilled mesh plus one record per boundary loop the fill saw."""

    mesh: Mesh
    holes: tuple[Mapping[str, Any], ...] = ()
    refused: tuple[Mapping[str, Any], ...] = ()
    step: RepairStep | None = None

    @property
    def filled(self) -> int:
        return len(self.holes)

    def to_dict(self) -> dict[str, Any]:
        return {
            "filled": self.filled,
            "refused": len(self.refused),
            "holes": [dict(hole) for hole in self.holes],
            "refusals": [dict(item) for item in self.refused],
            "step": self.step.to_dict() if self.step else None,
        }


def _boundary_loops(mesh: Mesh) -> tuple[tuple[int, ...], ...]:
    """Directed boundary loops, each as the vertex cycle its owning faces walk.

    The loop is returned in the direction the *shell* walks it, so a fill patch has to
    traverse it backwards to agree with its neighbours.
    """
    counts = meshlib.edge_face_counts(mesh)
    boundary = {edge for edge, count in counts.items() if count == 1}
    if not boundary:
        return ()
    following: dict[int, list[int]] = defaultdict(list)
    remaining: set[tuple[int, int]] = set()
    for face in mesh.faces:
        for k in range(3):
            a, b = face[k], face[(k + 1) % 3]
            if meshlib.undirected_edge((a, b)) in boundary:
                following[a].append(b)
                remaining.add((a, b))
    for vertex, targets in following.items():
        if len(targets) > 1:
            raise MeshRepairError(
                f"boundary vertex {vertex} is walked by {len(targets)} loops; splitting it is "
                "a manifoldisation, not a fill"
            )
    loops: list[tuple[int, ...]] = []
    while remaining:
        start, second = min(remaining)
        loop = [start]
        current, nxt = start, second
        while True:
            remaining.discard((current, nxt))
            if nxt == loop[0]:
                break
            if nxt in loop:
                raise MeshRepairError(f"boundary walk closes a short cycle at vertex {nxt}")
            loop.append(nxt)
            targets = [value for value in following[nxt] if (nxt, value) in remaining]
            if not targets:
                raise MeshRepairError(f"boundary walk dead-ends at vertex {nxt}")
            current, nxt = nxt, targets[0]
        loops.append(tuple(loop))
    return tuple(loops)


def boundary_loops(mesh: Mesh) -> tuple[tuple[int, ...], ...]:
    """Public view of the directed boundary walk :func:`fill_holes` is built on.

    Each loop is returned in the direction the *shell* walks it, so a cap patch has to
    traverse it backwards to agree with its neighbours. Raises :class:`MeshRepairError` when
    a boundary vertex is walked by more than one loop or a walk dead-ends.
    """
    return _boundary_loops(mesh)


def _newell_normal(points: Sequence[Vec3]) -> Vec3:
    nx = ny = nz = 0.0
    count = len(points)
    for index in range(count):
        x0, y0, z0 = points[index]
        x1, y1, z1 = points[(index + 1) % count]
        nx += (y0 - y1) * (z0 + z1)
        ny += (z0 - z1) * (x0 + x1)
        nz += (x0 - x1) * (y0 + y1)
    return (nx, ny, nz)


def _dot(a: Sequence[float], b: Sequence[float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a: Sequence[float], b: Sequence[float]) -> Vec3:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _unit(vector: Sequence[float]) -> Vec3:
    length = math.sqrt(_dot(vector, vector))
    if length == 0.0:
        raise MeshRepairError("degenerate direction: zero-length vector")
    return (vector[0] / length, vector[1] / length, vector[2] / length)


def _signed_area_2d(points: Sequence[tuple[float, float]]) -> float:
    total = 0.0
    count = len(points)
    for index in range(count):
        x0, y0 = points[index]
        x1, y1 = points[(index + 1) % count]
        total += x0 * y1 - x1 * y0
    return 0.5 * total


def _cross_2d(o: tuple[float, float], a: tuple[float, float], b: tuple[float, float]) -> float:
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def _segments_meet_2d(p1, p2, p3, p4) -> bool:
    d1 = _cross_2d(p3, p4, p1)
    d2 = _cross_2d(p3, p4, p2)
    d3 = _cross_2d(p1, p2, p3)
    d4 = _cross_2d(p1, p2, p4)
    return ((d1 > 0.0) != (d2 > 0.0)) and ((d3 > 0.0) != (d4 > 0.0))


def _point_in_triangle_2d(point, a, b, c) -> bool:
    d1 = _cross_2d(a, b, point)
    d2 = _cross_2d(b, c, point)
    d3 = _cross_2d(c, a, point)
    has_negative = d1 < 0.0 or d2 < 0.0 or d3 < 0.0
    has_positive = d1 > 0.0 or d2 > 0.0 or d3 > 0.0
    return not (has_negative and has_positive)


def _ear_clip(points: Sequence[tuple[float, float]]) -> list[tuple[int, int, int]]:
    """Triangulate a simple polygon, in the polygon's own vertex order."""
    count = len(points)
    if count < 3:
        raise MeshRepairError(f"a polygon needs at least three vertices, got {count}")
    orientation = 1.0 if _signed_area_2d(points) > 0.0 else -1.0
    scale = max(1e-300, max(abs(coordinate) for point in points for coordinate in point))
    epsilon = 1e-12 * scale * scale
    remaining = list(range(count))
    triangles: list[tuple[int, int, int]] = []
    guard = 0
    while len(remaining) > 3:
        guard += 1
        if guard > 4 * count * count:
            raise MeshRepairError("ear clipping did not terminate: polygon is not simple")
        clipped = False
        for position in range(len(remaining)):
            previous = remaining[position - 1]
            current = remaining[position]
            nxt = remaining[(position + 1) % len(remaining)]
            a, b, c = points[previous], points[current], points[nxt]
            if orientation * _cross_2d(a, b, c) <= epsilon:
                continue  # reflex corner, not an ear
            blocked = False
            for other in remaining:
                if other in (previous, current, nxt):
                    continue
                if _point_in_triangle_2d(points[other], a, b, c):
                    blocked = True
                    break
            if blocked:
                continue
            triangles.append((previous, current, nxt))
            del remaining[position]
            clipped = True
            break
        if not clipped:
            raise MeshRepairError("no ear found: polygon is not simple")
    triangles.append((remaining[0], remaining[1], remaining[2]))
    return triangles


def _loop_is_simple(points2d: Sequence[tuple[float, float]]) -> bool:
    count = len(points2d)
    for i in range(count):
        a, b = points2d[i], points2d[(i + 1) % count]
        for j in range(i + 1, count):
            if j == i or (j + 1) % count == i or j == (i + 1) % count:
                continue
            c, d = points2d[j], points2d[(j + 1) % count]
            if _segments_meet_2d(a, b, c, d):
                return False
    return True


def _triangle_area(points: Sequence[Vec3], face: Face) -> float:
    ax, ay, az = points[face[0]]
    bx, by, bz = points[face[1]]
    cx, cy, cz = points[face[2]]
    ux, uy, uz = bx - ax, by - ay, bz - az
    vx, vy, vz = cx - ax, cy - ay, cz - az
    return 0.5 * math.sqrt(
        (uy * vz - uz * vy) ** 2 + (uz * vx - ux * vz) ** 2 + (ux * vy - uy * vx) ** 2
    )


def fill_holes(mesh: Mesh, *, max_hole_edges: int = DEFAULT_MAX_HOLE_EDGES) -> HoleFill:
    """Close every simple boundary loop with a planar fill over the loop's own vertices.

    Loops that cannot be closed honestly -- longer than ``max_hole_edges``, projecting onto
    a self-intersecting polygon, collinear, or producing a zero-area facet -- are left open
    and recorded in :attr:`HoleFill.refused` with the reason. Nothing is silently dropped.
    """
    if max_hole_edges < 3:
        raise MeshRepairError(f"max_hole_edges must be at least 3, got {max_hole_edges!r}")
    if not mesh.faces:
        raise MeshRepairError("cannot fill holes in a mesh with no faces")
    loops = _boundary_loops(mesh)
    if not loops:
        return HoleFill(mesh=mesh, step=RepairStep("fill_holes", 0, "no boundary edges"))

    vertices = list(mesh.vertices)
    faces = list(mesh.faces)
    holes: list[Mapping[str, Any]] = []
    refused: list[Mapping[str, Any]] = []
    for loop in loops:
        if len(loop) > max_hole_edges:
            refused.append(
                {
                    "loop_vertices": len(loop),
                    "reason": f"loop has {len(loop)} edges, above the {max_hole_edges} cap; a "
                    "planar fill is not a repair at that size",
                }
            )
            continue
        points = [vertices[vertex] for vertex in loop]
        normal = _newell_normal(points)
        length = math.sqrt(_dot(normal, normal))
        scale = max(max(abs(coordinate) for coordinate in point) for point in points) or 1.0
        if length <= 1e-12 * scale * scale:
            refused.append(
                {"loop_vertices": len(loop), "reason": "loop is collinear: no plane to fill in"}
            )
            continue
        unit_normal = _unit(normal)
        first_edge = tuple(points[1][axis] - points[0][axis] for axis in range(3))
        tangential = tuple(
            first_edge[axis] - _dot(first_edge, unit_normal) * unit_normal[axis]
            for axis in range(3)
        )
        try:
            tangent = _unit(tangential)
        except MeshRepairError:
            refused.append(
                {
                    "loop_vertices": len(loop),
                    "reason": "the loop's first edge is parallel to the loop plane normal",
                }
            )
            continue
        bitangent = _cross(unit_normal, tangent)
        projected = [(_dot(point, tangent), _dot(point, bitangent)) for point in points]
        if not _loop_is_simple(projected):
            refused.append(
                {
                    "loop_vertices": len(loop),
                    "reason": "the loop projects onto a self-intersecting polygon",
                }
            )
            continue
        try:
            # The shell walks the loop one way; the patch has to walk it the other way.
            triangles = _ear_clip(list(reversed(projected)))
        except MeshRepairError as error:
            refused.append({"loop_vertices": len(loop), "reason": str(error)})
            continue
        reversed_loop = list(reversed(loop))
        patch: list[Face] = [
            (reversed_loop[a], reversed_loop[b], reversed_loop[c]) for a, b, c in triangles
        ]
        degenerate = sum(1 for face in patch if _triangle_area(points, _remap(face, loop)) <= 0.0)
        if degenerate:
            refused.append(
                {
                    "loop_vertices": len(loop),
                    "reason": f"fill produced {degenerate} zero-area facet(s)",
                }
            )
            continue
        faces.extend(patch)
        plane_distance = max(
            abs(_dot(tuple(point[axis] - points[0][axis] for axis in range(3)), unit_normal))
            for point in points
        )
        loop_area = abs(_signed_area_2d(projected))
        holes.append(
            {
                "loop_vertices": len(loop),
                "triangles": len(patch),
                "loop_area_m2": loop_area,
                "planar_deviation_m": plane_distance,
                "planar_deviation_relative": plane_distance / math.sqrt(loop_area or 1.0),
            }
        )

    result = _rebuild_with_faces(
        Mesh(
            vertices=tuple(vertices),
            faces=tuple(faces),
            units=mesh.units,
            coordinate_frame=mesh.coordinate_frame,
        ),
        faces,
    )
    filled_triangles = sum(int(hole["triangles"]) for hole in holes)
    return HoleFill(
        mesh=result,
        holes=tuple(holes),
        refused=tuple(refused),
        step=RepairStep(
            "fill_holes",
            filled_triangles,
            f"{len(holes)} loop(s) closed with {filled_triangles} facet(s); "
            f"{len(refused)} refused",
        ),
    )


def _remap(face: Face, loop: Sequence[int]) -> Face:
    """Map a patch face's local (0..n-1) indices onto positions in the loop's vertex list."""
    return (loop.index(face[0]), loop.index(face[1]), loop.index(face[2]))


# ---------------------------------------------------------------------------------------
# The whole repair, with an acceptance test
# ---------------------------------------------------------------------------------------


@dataclass(frozen=True)
class RepairResult:
    """The repaired mesh (or the original, on refusal) plus the evidence for either."""

    mesh: Mesh
    accepted: bool
    reason: str
    steps: tuple[RepairStep, ...]
    before: MeshMetrics
    after: MeshMetrics
    volume_relative_change: float | None
    area_relative_change: float | None
    holes: HoleFill | None = None
    reference_volume: float | None = None
    reference_area: float | None = None
    tolerated_volume_change: float = DEFAULT_VOLUME_TOLERANCE
    tolerated_area_change: float = DEFAULT_AREA_TOLERANCE
    notes: tuple[str, ...] = field(default_factory=tuple)

    def to_dict(self) -> dict[str, Any]:
        return {
            "accepted": self.accepted,
            "reason": self.reason,
            "steps": [step.to_dict() for step in self.steps],
            "before": self.before.to_dict(),
            "after": self.after.to_dict(),
            "volume_relative_change": self.volume_relative_change,
            "area_relative_change": self.area_relative_change,
            "reference_volume_m3": self.reference_volume,
            "reference_area_m2": self.reference_area,
            "tolerated_volume_change": self.tolerated_volume_change,
            "tolerated_area_change": self.tolerated_area_change,
            "holes": self.holes.to_dict() if self.holes else None,
            "notes": list(self.notes),
        }


def _relative_change(before: float, after: float) -> float | None:
    if before == 0.0:
        return None
    return abs(after - before) / abs(before)


def repair_mesh(
    mesh: Mesh,
    *,
    reference_volume: float | None = None,
    reference_area: float | None = None,
    volume_tolerance: float = DEFAULT_VOLUME_TOLERANCE,
    area_tolerance: float = DEFAULT_AREA_TOLERANCE,
    weld_tolerance: float = DEFAULT_WELD_TOLERANCE,
    max_hole_edges: int = DEFAULT_MAX_HOLE_EDGES,
    self_intersection_limit: int = 16,
    heal_self_intersections: bool = False,
) -> RepairResult:
    """Run the repair chain and accept it only if the result is measurably better.

    Acceptance is a conjunction, and every clause is a number that can be printed:

    * the result is closed (no boundary edge), edge-manifold, consistently wound, and has
      no coincident facets;
    * the result has no self intersection. This module does not remove one by itself, so a
      mesh that arrives with one is refused rather than quietly passed on -- unless
      ``heal_self_intersections=True`` is passed, in which case :mod:`tensorcad.detangle`
      excavates the fold and caps the hole first and the refusal only stands if that fails;
    * the enclosed volume -- when a ``reference_volume`` is supplied, or when the input was
      already closed and its own volume can serve as the reference -- moved by no more than
      ``volume_tolerance``; likewise the area against ``reference_area``.

    A refusal returns the **original** mesh. It is not an exception: a caller that asked a
    repair to run wants to know whether the mesh is usable, and ``accepted=False`` with a
    reason is a more useful answer than a traceback.
    """
    if volume_tolerance < 0.0 or area_tolerance < 0.0:
        raise MeshRepairError("tolerances must be non-negative")
    before = measure_mesh(mesh, self_intersection_limit=self_intersection_limit)
    notes: list[str] = []
    steps: list[RepairStep] = []

    current = meshlib.weld(mesh, tolerance=weld_tolerance)
    if len(current.vertices) != len(mesh.vertices):
        steps.append(
            RepairStep(
                "weld",
                len(mesh.vertices) - len(current.vertices),
                f"merged {len(mesh.vertices) - len(current.vertices)} vertex(es) at "
                f"{weld_tolerance:g} m",
            )
        )
    trimmed = meshlib.remove_degenerate_faces(current, tolerance=0.0)
    if len(trimmed.faces) != len(current.faces):
        steps.append(
            RepairStep(
                "drop_degenerate_faces",
                len(current.faces) - len(trimmed.faces),
                f"dropped {len(current.faces) - len(trimmed.faces)} zero-area facet(s)",
            )
        )
    current = trimmed

    try:
        current, step = drop_duplicate_faces(current)
        steps.append(step)
        current, step = manifoldize(current)
        steps.append(step)
        current, step = orient_consistently(current)
        steps.append(step)
        holes = fill_holes(current, max_hole_edges=max_hole_edges)
        if holes.step:
            steps.append(holes.step)
        current = holes.mesh
        if holes.filled:
            current, step = orient_consistently(current)
            if step.changed:
                steps.append(step)
    except MeshRepairError as error:
        return RepairResult(
            mesh=mesh,
            accepted=False,
            reason=f"repair step refused: {error} -- original mesh returned",
            steps=tuple(steps),
            before=before,
            after=before,
            volume_relative_change=None,
            area_relative_change=None,
            reference_volume=reference_volume,
            reference_area=reference_area,
            tolerated_volume_change=volume_tolerance,
            tolerated_area_change=area_tolerance,
            notes=tuple(notes),
        )

    if heal_self_intersections and len(current.faces):
        probe = measure_mesh(current, self_intersection_limit=self_intersection_limit)
        if probe.self_intersections:
            from .detangle import heal_self_intersections as _heal_self_intersections

            healed = _heal_self_intersections(
                current,
                self_intersection_limit=max(self_intersection_limit, 64),
                reference_volume=reference_volume,
            )
            notes.extend(f"heal: {note}" for note in healed.notes)
            if healed.healed:
                steps.append(
                    RepairStep(
                        "heal_self_intersections",
                        healed.faces_removed,
                        f"{probe.self_intersections} penetrating pair(s) removed by cutting "
                        f"{healed.faces_removed} facet(s) in {healed.rings_used} ring(s) and "
                        f"capping {healed.holes.filled if healed.holes else 0} loop(s)",
                    )
                )
                current = healed.mesh
            else:
                notes.append(f"self-intersection heal refused: {healed.reason}")

    after = measure_mesh(current, self_intersection_limit=self_intersection_limit)

    if reference_volume is None and before.closed and before.signed_volume > 0.0:
        reference_volume = before.signed_volume
        notes.append("reference volume taken from the (closed) input mesh")
    elif reference_volume is None and not before.closed:
        notes.append(
            "volume check skipped: the input was open and no reference volume was supplied"
        )
    if reference_area is None:
        notes.append("area check skipped: a fill adds its own area, so a reference is required")

    volume_change = (
        _relative_change(reference_volume, after.signed_volume)
        if reference_volume is not None and after.closed
        else None
    )
    area_change = (
        _relative_change(reference_area, after.surface_area)
        if reference_area is not None
        else None
    )

    problems: list[str] = []
    if not after.closed:
        problems.append(f"still open: {after.open_edges} boundary edge(s)")
    if not after.manifold:
        problems.append(f"{after.nonmanifold_edges} non-manifold edge(s) left")
    if not after.consistently_wound:
        problems.append(f"{after.winding_conflicts} winding conflict(s) left")
    if after.duplicate_faces:
        problems.append(f"{after.duplicate_faces} coincident facet(s) left")
    if after.self_intersections:
        problems.append(
            f"{after.self_intersections} self-intersecting face pair(s); removing them is a "
            "remeshing problem, not a repair this module performs"
        )
    if volume_change is not None and volume_change > volume_tolerance:
        problems.append(
            f"volume moved {volume_change:.3%} against the reference (limit "
            f"{volume_tolerance:.3%})"
        )
    if area_change is not None and area_change > area_tolerance:
        problems.append(
            f"area moved {area_change:.3%} against the reference (limit {area_tolerance:.3%})"
        )
    refused_holes = len(holes.refused)
    if refused_holes:
        problems.append(
            f"{refused_holes} boundary loop(s) could not be filled: "
            + "; ".join(str(item["reason"]) for item in holes.refused)
        )

    if problems:
        return RepairResult(
            mesh=mesh,
            accepted=False,
            reason="; ".join(problems) + " -- original mesh returned",
            steps=tuple(steps),
            before=before,
            after=after,
            volume_relative_change=volume_change,
            area_relative_change=area_change,
            holes=holes,
            reference_volume=reference_volume,
            reference_area=reference_area,
            tolerated_volume_change=volume_tolerance,
            tolerated_area_change=area_tolerance,
            notes=tuple(notes),
        )

    changed = sum(step.changed for step in steps)
    if changed == 0:
        reason = "nothing to repair: the mesh was already closed, manifold and clean"
    else:
        reason = "; ".join(f"{step.name}: {step.detail}" for step in steps if step.changed)
    return RepairResult(
        mesh=current,
        accepted=True,
        reason=reason,
        steps=tuple(steps),
        before=before,
        after=after,
        volume_relative_change=volume_change,
        area_relative_change=area_change,
        holes=holes,
        reference_volume=reference_volume,
        reference_area=reference_area,
        tolerated_volume_change=volume_tolerance,
        tolerated_area_change=area_tolerance,
        notes=tuple(notes),
    )