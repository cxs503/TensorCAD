"""Generic, dependency-free mesh algorithms: topology checks and surface integrals.

These are the guarantees a geometry foundation actually has to provide before a
mesh is "solver-ready": closed, edge-manifold, consistently outward-oriented,
and free of zero-area faces. Everything here works on any :class:`Mesh`, so the
checks can also be run on meshes imported from a CAD kernel or a mesh file.

Conventions
-----------
* Faces are 0-based triangles.
* An edge is *manifold* when exactly two faces use it; an edge used once is a
  boundary edge (the mesh is open) and an edge used more than twice is
  non-manifold.
* Volume and centroid integrals use the divergence theorem over tetrahedra
  formed with the coordinate origin, so they are only meaningful for a closed
  surface; :func:`signed_volume` is nevertheless well defined for open ones and
  is what orientation checks rely on.
"""
from __future__ import annotations

import math
from collections import defaultdict, deque
from dataclasses import dataclass

from .geometry import Face, Mesh, Vec3


def face_edges(face: Face) -> tuple[tuple[int, int], ...]:
    """Directed edges of a triangle, in winding order."""
    a, b, c = face
    return ((a, b), (b, c), (c, a))


def undirected_edge(edge: tuple[int, int]) -> tuple[int, int]:
    a, b = edge
    return (a, b) if a <= b else (b, a)


def edge_face_counts(mesh: Mesh) -> dict[tuple[int, int], int]:
    counts: dict[tuple[int, int], int] = defaultdict(int)
    for face in mesh.faces:
        for edge in face_edges(face):
            counts[undirected_edge(edge)] += 1
    return dict(counts)


def boundary_edges(mesh: Mesh) -> tuple[tuple[int, int], ...]:
    """Edges belonging to exactly one face, i.e. holes in the surface."""
    edges = [edge for edge, count in edge_face_counts(mesh).items() if count == 1]
    return tuple(sorted(edges))


def is_watertight(mesh: Mesh) -> bool:
    """True when every edge is shared by exactly two faces (closed surface)."""
    if not mesh.faces:
        return False
    return all(count == 2 for count in edge_face_counts(mesh).values())


def is_edge_manifold(mesh: Mesh) -> bool:
    return all(count <= 2 for count in edge_face_counts(mesh).values())


def triangle_area(mesh: Mesh, face: Face) -> float:
    ax, ay, az = mesh.vertices[face[0]]
    bx, by, bz = mesh.vertices[face[1]]
    cx, cy, cz = mesh.vertices[face[2]]
    ux, uy, uz = bx - ax, by - ay, bz - az
    vx, vy, vz = cx - ax, cy - ay, cz - az
    return 0.5 * math.sqrt(
        (uy * vz - uz * vy) ** 2 + (uz * vx - ux * vz) ** 2 + (ux * vy - uy * vx) ** 2
    )


def triangle_areas(mesh: Mesh) -> tuple[float, ...]:
    return tuple(triangle_area(mesh, face) for face in mesh.faces)


def surface_area(mesh: Mesh) -> float:
    return math.fsum(triangle_areas(mesh))


def total_signed_volume(mesh: Mesh) -> float:
    """Divergence-theorem volume of the surface, positive when oriented outward."""
    return sum(_tet_volume(mesh, face) for face in mesh.faces) / 6.0


def _tet_volume(mesh: Mesh, face: Face) -> float:
    ax, ay, az = mesh.vertices[face[0]]
    bx, by, bz = mesh.vertices[face[1]]
    cx, cy, cz = mesh.vertices[face[2]]
    # a . (b x c)
    return ax * (by * cz - bz * cy) + ay * (bz * cx - bx * cz) + az * (bx * cy - by * cx)


def face_normal(mesh: Mesh, face: Face) -> Vec3:
    """Area-weighted normal (twice the area) pointing by the right-hand rule."""
    ax, ay, az = mesh.vertices[face[0]]
    bx, by, bz = mesh.vertices[face[1]]
    cx, cy, cz = mesh.vertices[face[2]]
    ux, uy, uz = bx - ax, by - ay, bz - az
    vx, vy, vz = cx - ax, cy - ay, cz - az
    return (uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx)


def face_centroid(mesh: Mesh, face: Face) -> Vec3:
    return tuple(sum(mesh.vertices[v][axis] for v in face) / 3.0 for axis in range(3))  # type: ignore[return-value]


def zero_area_faces(mesh: Mesh, tolerance: float = 0.0) -> tuple[int, ...]:
    """Indices of faces whose area is at or below ``tolerance``."""
    return tuple(index for index, area in enumerate(triangle_areas(mesh)) if area <= tolerance)


def remove_degenerate_faces(mesh: Mesh, tolerance: float = 0.0) -> Mesh:
    """Drop zero-area faces and compact the vertex list to what is still used."""
    keep = [
        index
        for index, area in enumerate(triangle_areas(mesh))
        if area > tolerance
    ]
    if len(keep) == len(mesh.faces):
        return mesh
    used: dict[int, int] = {}
    vertices: list[Vec3] = []
    for index in keep:
        for vertex in mesh.faces[index]:
            if vertex not in used:
                used[vertex] = len(vertices)
                vertices.append(mesh.vertices[vertex])
    faces = tuple(tuple(used[vertex] for vertex in mesh.faces[index]) for index in keep)
    return Mesh(
        vertices=tuple(vertices),
        faces=faces,  # type: ignore[arg-type]
        units=mesh.units,
        coordinate_frame=mesh.coordinate_frame,
    )


def weld(mesh: Mesh, tolerance: float = 1e-9) -> Mesh:
    """Merge vertices closer than ``tolerance``, rewriting face indices.

    Generation code should build shared corners exactly (see
    :class:`tensorcad.geometry.MeshBuilder`); this exists for *imported* meshes,
    where the file itself decides how close coincident points are.
    """
    if tolerance < 0:
        raise ValueError("tolerance must be non-negative")
    for point in mesh.vertices:
        if not all(math.isfinite(value) for value in point):
            raise ValueError("weld needs finite vertex coordinates")
    grid: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    # With tolerance == 0 the cell size is infinite, so every point lands in cell
    # (0, 0, 0) and merging degrades to exact equality.
    cell = tolerance if tolerance > 0 else float("inf")
    merged: list[int] = list(range(len(mesh.vertices)))
    for index, point in enumerate(mesh.vertices):
        bucket: tuple[int, int, int] = (
            math.floor(point[0] / cell),
            math.floor(point[1] / cell),
            math.floor(point[2] / cell),
        )
        for other in grid[bucket]:
            other_point = mesh.vertices[other]
            if all(abs(a - b) <= tolerance for a, b in zip(point, other_point, strict=False)):
                merged[index] = other
                break
        else:
            # Also look in the 26 neighbouring cells so clusters straddling a
            # cell boundary still merge.
            found = False
            if math.isfinite(cell):
                for dx in (-1, 0, 1):
                    for dy in (-1, 0, 1):
                        for dz in (-1, 0, 1):
                            neighbour: tuple[int, int, int] = (bucket[0] + dx, bucket[1] + dy, bucket[2] + dz)
                            for other in grid.get(neighbour, ()):
                                other_point = mesh.vertices[other]
                                if all(
                                    abs(p - q) <= tolerance for p, q in zip(point, other_point, strict=False)
                                ):
                                    merged[index] = other
                                    found = True
                                    break
                            if found:
                                break
                        if found:
                            break
                    if found:
                        break
            if not found:
                grid[bucket].append(index)
    remap = sorted({merged[vertex] for vertex in range(len(mesh.vertices))})
    compact = {old: new for new, old in enumerate(remap)}
    vertices = tuple(mesh.vertices[old] for old in remap)
    faces = tuple(tuple(compact[merged[vertex]] for vertex in face) for face in mesh.faces)
    return Mesh(
        vertices=vertices,
        faces=faces,  # type: ignore[arg-type]
        units=mesh.units,
        coordinate_frame=mesh.coordinate_frame,
    )


@dataclass(frozen=True)
class OrientationResult:
    """Outcome of :func:`orient_outward`."""

    mesh: Mesh
    components: int
    flipped_faces: int
    volume: float


def _edge_bits(face: Face) -> dict[tuple[int, int], int]:
    """Per undirected edge, +1 if the face traverses it low->high else -1."""
    bits: dict[tuple[int, int], int] = {}
    for u, v in face_edges(face):
        bits[(u, v) if u < v else (v, u)] = 1 if u < v else -1
    return bits


def orient_outward(mesh: Mesh) -> OrientationResult:
    """Make face winding consistent, then point it away from the enclosed volume.

    Pass one walks shared edges breadth-first and forces neighbours to traverse a
    common edge in opposite directions, which is the definition of consistent
    winding; mirrored or inverted patches are fixed by flipping whole faces, and
    the walk propagates across connected components only. Pass two flips each
    component whose signed volume is negative.

    Judging a component by its own volume is exact for a closed component. For an
    open one it is a documented heuristic, which is why callers must also check
    :func:`is_watertight` before trusting a volume.
    """
    raw: list[Face] = list(mesh.faces)
    bits = [_edge_bits(face) for face in raw]

    users: dict[tuple[int, int], list[int]] = defaultdict(list)
    for index, face in enumerate(raw):
        for edge in face_edges(face):
            users[undirected_edge(edge)].append(index)

    adjacency: dict[int, set[tuple[int, tuple[int, int]]]] = defaultdict(set)
    for edge, owners in users.items():
        for i in range(len(owners)):
            for j in range(i + 1, len(owners)):
                if owners[i] == owners[j]:
                    continue
                adjacency[owners[i]].add((owners[j], edge))
                adjacency[owners[j]].add((owners[i], edge))

    sign = [0] * len(raw)  # 0 = unvisited, +1 keep, -1 flip
    components: list[list[int]] = []
    for seed in range(len(raw)):
        if sign[seed]:
            continue
        sign[seed] = 1
        queue: deque[int] = deque([seed])
        component = [seed]
        while queue:
            current = queue.popleft()
            for other, edge in adjacency.get(current, ()):
                if sign[other]:
                    continue
                shared = bits[current][edge] * bits[other][edge]
                sign[other] = -sign[current] * shared
                queue.append(other)
                component.append(other)
        components.append(component)

    oriented: list[Face] = [
        raw[index][::-1] if sign[index] < 0 else raw[index] for index in range(len(raw))
    ]
    flipped = sum(1 for value in sign if value < 0)

    for component in components:
        volume = math.fsum(
            _tet_volume(mesh, oriented[index]) for index in component
        )
        if volume < 0:
            for index in component:
                oriented[index] = oriented[index][::-1]
                flipped += 1 if sign[index] > 0 else 0

    result = Mesh(
        vertices=mesh.vertices,
        faces=tuple(oriented),
        units=mesh.units,
        coordinate_frame=mesh.coordinate_frame,
    )
    return OrientationResult(
        mesh=result,
        components=len(components),
        flipped_faces=flipped,
        volume=total_signed_volume(result),
    )
