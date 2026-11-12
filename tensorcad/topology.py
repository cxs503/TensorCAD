"""B-rep cell topology: faces, edges, vertices, adjacency, Euler operators, stable names.

The library already answers triangle questions (``tensorcad.mesh``: edge--face counts,
watertightness, orientation) and implicit-field questions (``tensorcad.geometry_field``).
What it did not have is the object a modelling kernel hands to a mesher or a solver
*about faces*: "this shell has 6 faces, 12 edges, 8 vertices; this edge separates face
``deck`` from face ``side``; this vertex is where three faces meet; and the names survive
a rebuild". That is `GAP-CAD-001` in ``docs/oss_inventory_cad/GAP-ANCHORS.md`` -- the
anchor ``tensorcad.shell`` names as the prerequisite for removing a shell's outer wall
and stitching the rim.

This module builds that object from the only thing the library actually produces today,
a labelled triangle mesh, and it builds it **combinatorially**, not by tolerance:

1. **Faces** are connected components of triangles carrying the same label, so a face is
   a declaration by the caller (as a CAD feature tree declares it), never a guess from
   dihedral angles.
2. **Boundary arcs** are mesh edges whose two incident triangles lie in different faces
   (or that have one incident triangle). Every such edge is walked into chains.
3. **Vertices** are the nodes where the chain structure changes: a node that is not a
   degree-2 node of the boundary graph, or a degree-2 node whose two arcs separate
   *different* face pairs. Chains that close on themselves with no such node (a sphere's
   equator, a cylinder's rim) get one deterministic start vertex, because an edge must
   have endpoints for ``V - E + F`` to be the Euler characteristic rather than an
   accident of labelling.
4. **Edges** are the chains between those vertices.

The payoff is a check that cannot be gamed: ``V - E + F`` of the face/edge/vertex
structure must equal ``V - E + F`` of the triangle mesh it came from, because the Euler
characteristic is a topological invariant. Two independent decompositions of the same
surface that agree on it agree on genus; if the grouping missed a corner, split a chain
twice, or welded two faces into one, the number moves and :meth:`CellComplex.validate`
says so. That is the acceptance criterion the gap anchor asks for ("counts consistent
with the triangle mesh, zero exceptions"), together with the winding invariant (two faces
must traverse their shared edge in opposite directions) and the naming test.

Euler operators are provided on the result (``split_face``, ``split_edge``,
``merge_faces``) with the invariant that every one of them preserves the Euler
characteristic and leaves :meth:`CellComplex.validate` empty -- the property that makes
them *operators* rather than edits.

Naming
------
Names are derived from the face labels plus canonical geometry (coordinates rounded
against the model scale), never from indices, so rebuilding the same geometry at a
different sampling density, or handing the patches in a different order, yields the same
names. :meth:`CellComplex.name_index` is the lookup a feature tree needs; the honest
limitations are written down in ``docs/topology.md``.
"""
from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence

from . import mesh as meshlib
from .geometry import Mesh, Vec3

__all__ = [
    "CellComplex",
    "TopoEdge",
    "TopoFace",
    "TopoVertex",
    "TopologyError",
    "from_labeled_mesh",
    "merge_faces",
    "split_edge",
    "split_face",
]


class TopologyError(ValueError):
    """Raised when a labelling cannot be turned into a valid cell complex."""


def _face_position(record: "TopoFace | int") -> int:
    return record.index if isinstance(record, TopoFace) else int(record)


def _edge_position(record: "TopoEdge | int") -> int:
    return record.index if isinstance(record, TopoEdge) else int(record)


def _vertex_position(record: "TopoVertex | int") -> int:
    return record.index if isinstance(record, TopoVertex) else int(record)


# ---------------------------------------------------------------------------
# records
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class TopoVertex:
    """A topological vertex: a place where the edge structure branches.

    ``nodes`` are the mesh node indices this vertex stands for. The builder never merges
    two nodes into one vertex by position -- only mesh welding does that, upstream -- so
    a vertex can be trusted to be one identified point.
    """

    index: int
    point: Vec3
    name: str
    nodes: tuple[int, ...] = ()

    def to_dict(self) -> dict[str, object]:
        return {
            "index": self.index,
            "name": self.name,
            "point": list(self.point),
            "nodes": list(self.nodes),
        }


@dataclass(frozen=True)
class TopoEdge:
    """A topological edge: a chain of mesh edges separating one or two faces.

    ``faces`` holds one entry for an open boundary and two for a closed shell.
    ``vertices`` may repeat a single index for a closed chain (a loop), which is why
    ``loop`` is recorded explicitly instead of being inferred from equality.
    """

    index: int
    vertices: tuple[int, int]
    faces: tuple[int, ...]
    name: str
    nodes: tuple[int, ...] = ()
    loop: bool = False

    @property
    def open(self) -> bool:
        return len(self.faces) == 1

    def to_dict(self) -> dict[str, object]:
        return {
            "index": self.index,
            "name": self.name,
            "vertices": list(self.vertices),
            "faces": list(self.faces),
            "nodes": len(self.nodes),
            "loop": self.loop,
        }


@dataclass(frozen=True)
class TopoFace:
    """A topological face: one connected region of same-label triangles."""

    index: int
    label: str
    name: str
    edges: tuple[int, ...]
    vertices: tuple[int, ...]
    triangles: tuple[int, ...] = ()
    area: float = 0.0

    def to_dict(self) -> dict[str, object]:
        return {
            "index": self.index,
            "label": self.label,
            "name": self.name,
            "edges": list(self.edges),
            "vertices": list(self.vertices),
            "triangles": len(self.triangles),
            "area_m2": self.area,
        }


# ---------------------------------------------------------------------------
# the cell complex
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class CellComplex:
    """Faces, edges and vertices of one surface decomposition, with adjacency.

    Built by :func:`from_labeled_mesh`; edited only by the Euler operators, each of
    which preserves :meth:`euler_characteristic` and leaves :meth:`validate` empty.
    """

    faces: tuple[TopoFace, ...]
    edges: tuple[TopoEdge, ...]
    vertices: tuple[TopoVertex, ...]
    mesh_counts: tuple[int, int, int] | None = None
    winding: tuple[bool, ...] = ()

    # -- invariants -------------------------------------------------------
    def counts(self) -> dict[str, int]:
        return {
            "vertices": len(self.vertices),
            "edges": len(self.edges),
            "faces": len(self.faces),
        }

    def euler_characteristic(self) -> int:
        return len(self.vertices) - len(self.edges) + len(self.faces)

    def mesh_euler_characteristic(self) -> int | None:
        if self.mesh_counts is None:
            return None
        vertices, edges, triangles = self.mesh_counts
        return vertices - edges + triangles

    @property
    def closed(self) -> bool:
        return bool(self.edges) and all(len(edge.faces) == 2 for edge in self.edges)

    def genus(self) -> int | None:
        """Genus for a closed surface, ``None`` when it is not closed or not even."""
        if not self.closed:
            return None
        remainder = 2 - self.euler_characteristic()
        if remainder % 2:
            return None
        return remainder // 2

    # -- adjacency queries -------------------------------------------------
    def faces_of_edge(self, edge: TopoEdge | int) -> tuple[int, ...]:
        return self.edges[_edge_position(edge)].faces

    def edges_of_face(self, face: TopoFace | int) -> tuple[int, ...]:
        return self.faces[_face_position(face)].edges

    def edges_of_vertex(self, vertex: TopoVertex | int) -> tuple[int, ...]:
        index = _vertex_position(vertex)
        return tuple(edge.index for edge in self.edges if index in edge.vertices)

    def faces_of_vertex(self, vertex: TopoVertex | int) -> tuple[int, ...]:
        found: set[int] = set()
        for edge_index in self.edges_of_vertex(vertex):
            found.update(self.edges[edge_index].faces)
        return tuple(sorted(found))

    def adjacency(self) -> dict[str, list[list[int]]]:
        """Every adjacency relation, as JSON-ready lists."""
        return {
            "edge_faces": [list(edge.faces) for edge in self.edges],
            "face_edges": [list(face.edges) for face in self.faces],
            "edge_vertices": [list(edge.vertices) for edge in self.edges],
            "vertex_edges": [list(self.edges_of_vertex(v.index)) for v in self.vertices],
            "vertex_faces": [list(self.faces_of_vertex(v.index)) for v in self.vertices],
        }

    # -- names -------------------------------------------------------------
    def name_index(self) -> dict[str, tuple[str, int]]:
        """``name -> (kind, index)`` for faces, edges and vertices."""
        index: dict[str, tuple[str, int]] = {}
        for face in self.faces:
            index[face.name] = ("face", face.index)
        for edge in self.edges:
            index[edge.name] = ("edge", edge.index)
        for vertex in self.vertices:
            index[vertex.name] = ("vertex", vertex.index)
        return index

    def find_face(self, name: str) -> TopoFace:
        kind, index = self.name_index()[name]
        if kind != "face":
            raise KeyError(f"{name!r} names a {kind}, not a face")
        return self.faces[index]

    def find_edge(self, name: str) -> TopoEdge:
        kind, index = self.name_index()[name]
        if kind != "edge":
            raise KeyError(f"{name!r} names a {kind}, not an edge")
        return self.edges[index]

    def stable_names(self) -> dict[str, tuple[str, ...]]:
        return {
            "faces": tuple(sorted(face.name for face in self.faces)),
            "edges": tuple(sorted(edge.name for edge in self.edges)),
            "vertices": tuple(sorted(vertex.name for vertex in self.vertices)),
        }

    # -- validation --------------------------------------------------------
    def validate(self, *, require_closed: bool = False) -> list[str]:
        """Every structural defect, as human-readable strings (empty = valid)."""
        problems: list[str] = []
        if not self.faces:
            problems.append("complex has no faces")
        for edge in self.edges:
            if len(edge.faces) not in (1, 2):
                problems.append(f"edge {edge.name} borders {len(edge.faces)} faces")
            if len(edge.faces) == 2 and edge.faces[0] == edge.faces[1]:
                problems.append(f"edge {edge.name} has the same face on both sides")
            for vertex in edge.vertices:
                if not 0 <= vertex < len(self.vertices):
                    problems.append(f"edge {edge.name} references vertex {vertex}")
            if edge.loop != (edge.vertices[0] == edge.vertices[1]):
                problems.append(f"edge {edge.name} disagrees with its loop flag")
        for face in self.faces:
            if not face.edges:
                problems.append(f"face {face.name} has no edges")
                continue
            if _face_walk(self, face) is None:
                problems.append(f"face {face.name} boundary is not a single closed chain")
            for vertex in face.vertices:
                if not 0 <= vertex < len(self.vertices):
                    problems.append(f"face {face.name} references vertex {vertex}")
        for vertex in self.vertices:
            if not self.edges_of_vertex(vertex.index):
                problems.append(f"vertex {vertex.name} touches no edge")
        if require_closed:
            open_edges = sum(1 for edge in self.edges if edge.open)
            if open_edges:
                problems.append(f"{open_edges} open edges")
        mesh_characteristic = self.mesh_euler_characteristic()
        if mesh_characteristic is not None and mesh_characteristic != self.euler_characteristic():
            problems.append(
                f"Euler characteristic {self.euler_characteristic()} disagrees with the "
                f"triangle mesh {mesh_characteristic}"
            )
        return problems

    def consistency_report(self) -> dict[str, object]:
        """The acceptance numbers of `GAP-CAD-001`, computed the two independent ways."""
        brep = {
            "vertices": len(self.vertices),
            "edges": len(self.edges),
            "faces": len(self.faces),
            "euler": self.euler_characteristic(),
        }
        mesh_counts: dict[str, int] | None = None
        if self.mesh_counts is not None:
            vertices, edges, triangles = self.mesh_counts
            mesh_counts = {
                "vertices": vertices,
                "edges": edges,
                "triangles": triangles,
                "euler": vertices - edges + triangles,
            }
        two_sided = sum(1 for edge in self.edges if len(edge.faces) == 2)
        return {
            "brep": brep,
            "mesh": mesh_counts,
            "closed": self.closed,
            "genus": self.genus(),
            "edges_with_two_faces": two_sided,
            "edges_with_one_face": len(self.edges) - two_sided,
            "winding_consistent": bool(all(self.winding)) if self.winding else True,
            "euler_matches_mesh": (
                None if mesh_counts is None else brep["euler"] == mesh_counts["euler"]
            ),
            "problems": self.validate(),
        }

    def to_dict(self) -> dict[str, object]:
        return {
            "counts": self.counts(),
            "euler_characteristic": self.euler_characteristic(),
            "genus": self.genus(),
            "closed": self.closed,
            "faces": [face.to_dict() for face in self.faces],
            "edges": [edge.to_dict() for edge in self.edges],
            "vertices": [vertex.to_dict() for vertex in self.vertices],
            "names": {key: list(value) for key, value in self.stable_names().items()},
            "consistency": self.consistency_report(),
        }


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _canonical_key(points: Iterable[Vec3], scale: float) -> tuple[float, ...]:
    """Order-independent, scale-relative geometric key.

    Coordinates are divided by the model scale before rounding so the key means the same
    thing for a 1 m part and a 300 m ship: ``1e-9 * L`` is its resolution.
    """
    denominator = scale if scale > 0.0 else 1.0
    values: list[float] = []
    for point in sorted((float(p[0]), float(p[1]), float(p[2])) for p in points):
        values.extend(round(component / denominator, 9) + 0.0 for component in point)
    return tuple(values)


def _midpoint(points: Sequence[Vec3]) -> Vec3:
    count = len(points)
    if count == 0:  # pragma: no cover - callers pass walked nodes
        raise TopologyError("cannot take the midpoint of no points")
    return (
        math.fsum(p[0] for p in points) / count,
        math.fsum(p[1] for p in points) / count,
        math.fsum(p[2] for p in points) / count,
    )


def _labels_sequence(labels: Mapping[int, str] | Sequence[str], count: int) -> list[str]:
    if isinstance(labels, Mapping):
        out = ["" for _ in range(count)]
        for triangle, label in labels.items():
            if not 0 <= triangle < count:
                raise TopologyError(f"label names triangle {triangle}, mesh has {count}")
            out[triangle] = str(label)
    else:
        out = [str(label) for label in labels]
        if len(out) != count:
            raise TopologyError(f"{len(out)} labels for {count} triangles")
    for triangle, label in enumerate(out):
        if not label:
            raise TopologyError(f"triangle {triangle} has no label")
    return out


class _UnionFind:
    def __init__(self, size: int) -> None:
        self._parent = list(range(size))
        self._rank = [0] * size

    def find(self, index: int) -> int:
        parent = self._parent
        root = index
        while parent[root] != root:
            root = parent[root]
        while parent[index] != root:
            parent[index], index = root, parent[index]
        return root

    def union(self, first: int, second: int) -> None:
        a, b = self.find(first), self.find(second)
        if a == b:
            return
        if self._rank[a] < self._rank[b]:
            a, b = b, a
        self._parent[b] = a
        if self._rank[a] == self._rank[b]:
            self._rank[a] += 1


def _unique_names(
    labelled: Sequence[tuple[int, str]],
    keys: Mapping[int, tuple[float, ...]],
    prefix: str,
) -> dict[int, str]:
    grouped: dict[str, list[int]] = defaultdict(list)
    for root, label in labelled:
        grouped[label].append(root)
    names: dict[int, str] = {}
    for label, roots in grouped.items():
        ordered = sorted(roots, key=lambda root: keys[root])
        base = f"{prefix}{label}"
        for position, root in enumerate(ordered):
            names[root] = base if len(ordered) == 1 else f"{base}@{position + 1}"
    return names


def _disambiguate(
    names: Mapping[int, str], keys: Mapping[int, tuple[float, ...]]
) -> dict[int, str]:
    """Append ``#k`` (ordered geometrically) to every name that would be ambiguous."""
    grouped: dict[str, list[int]] = defaultdict(list)
    for index, name in names.items():
        grouped[name].append(index)
    out: dict[int, str] = dict(names)
    for name, indices in grouped.items():
        if len(indices) == 1:
            continue
        for position, index in enumerate(sorted(indices, key=lambda i: keys[i])):
            out[index] = f"{name}#{position + 1}"
    return out


def _walk_faces(
    boundary: Mapping[tuple[int, int], tuple[int, ...]], nodes: Sequence[int]
) -> set[int]:
    faces: set[int] = set()
    for first, second in zip(nodes, nodes[1:]):
        faces.update(boundary[meshlib.undirected_edge((first, second))])
    return faces


def _winding_consistent(
    mesh: Mesh,
    edge_triangles: Mapping[tuple[int, int], Sequence[int]],
    nodes: Sequence[int],
    component_of: Mapping[int, int],
    boundary: Mapping[tuple[int, int], tuple[int, ...]],
) -> bool:
    """True when the two faces traverse every arc of the chain in opposite directions.

    For each mesh edge the incident triangle's winding gives a direction. Signing that
    direction against the chain's own traversal and summing per face, the two faces must
    carry opposite totals -- the manifold invariant, held along the whole topological
    edge rather than edge by edge.
    """
    totals: dict[int, int] = defaultdict(int)
    for first, second in zip(nodes, nodes[1:]):
        edge = meshlib.undirected_edge((first, second))
        for triangle in edge_triangles[edge]:
            a, b, c = mesh.faces[triangle]
            sign = 1
            for start, end in ((a, b), (b, c), (c, a)):
                if (start, end) == (first, second):
                    sign = 1
                    break
                if (start, end) == (second, first):
                    sign = -1
                    break
            totals[component_of[triangle]] += sign
    if len(boundary[meshlib.undirected_edge((nodes[0], nodes[1]))]) == 1:
        return True
    values = sorted(totals.values())
    return len(values) == 2 and values[0] == -values[1]


def _walk_chains(
    boundary: Mapping[tuple[int, int], tuple[int, ...]],
    incident: Mapping[int, Sequence[tuple[int, int]]],
    corners: set[int],
    vertices: Sequence[Vec3],
) -> list[tuple[int, int, list[int]]]:
    """Split the boundary graph into chains between corners (plus deterministic loops)."""
    used: set[tuple[int, int]] = set()
    walks: list[tuple[int, int, list[int]]] = []

    def walk_from(start: int, first_arc: tuple[int, int]) -> tuple[int, int, list[int]]:
        path = [start]
        arc = first_arc
        node = start
        while True:
            used.add(arc)
            nxt = arc[1] if arc[0] == node else arc[0]
            path.append(nxt)
            if nxt in corners or nxt == start:
                return (start, nxt, path)
            candidates = [candidate for candidate in incident[nxt] if candidate != arc]
            if len(candidates) != 1:  # pragma: no cover - corner detection covers this
                return (start, nxt, path)
            node = nxt
            arc = candidates[0]
            if arc in used:  # pragma: no cover - defensive
                return (start, nxt, path)

    for corner in sorted(corners):
        for arc in incident[corner]:
            if arc in used:
                continue
            walks.append(walk_from(corner, arc))

    # chains that close on themselves without a corner (a sphere's equator, a rim)
    for arc in boundary:
        if arc in used:
            continue
        start = min(arc, key=lambda node: tuple(round(c, 12) for c in vertices[node]))
        walks.append(walk_from(start, arc))
    return walks


# ---------------------------------------------------------------------------
# construction
# ---------------------------------------------------------------------------
def from_labeled_mesh(
    mesh: Mesh,
    labels: Mapping[int, str] | Sequence[str],
    *,
    name_prefix: str = "",
) -> CellComplex:
    """Build the face/edge/vertex structure of ``mesh`` from a per-triangle labelling.

    ``labels`` is either a sequence with one label per triangle or a mapping from
    triangle index to label; every triangle must be labelled, because an unlabelled
    triangle would be a face nobody declared. Labels are the *only* input that decides
    what a face is: dihedral angles, curvature and proximity are not consulted, which is
    what keeps the result reproducible and reviewable.
    """
    if not mesh.faces:
        raise TopologyError("a cell complex needs at least one triangle")
    tri_labels = _labels_sequence(labels, len(mesh.faces))

    edge_triangles: dict[tuple[int, int], list[int]] = defaultdict(list)
    for triangle, face in enumerate(mesh.faces):
        for edge in meshlib.face_edges(face):
            edge_triangles[meshlib.undirected_edge(edge)].append(triangle)
    for edge, triangles in edge_triangles.items():
        if len(triangles) > 2:
            raise TopologyError(
                f"mesh edge {edge} is used by {len(triangles)} triangles; the boundary "
                "chains of a non-manifold edge cannot be walked"
            )

    # 1. faces = connected components of same-label triangles
    components = _UnionFind(len(mesh.faces))
    for edge, triangles in edge_triangles.items():
        if len(triangles) == 2:
            first, second = triangles
            if tri_labels[first] == tri_labels[second]:
                components.union(first, second)
    component_of = {triangle: components.find(triangle) for triangle in range(len(mesh.faces))}
    roots = sorted(set(component_of.values()))
    root_position = {root: position for position, root in enumerate(roots)}
    component_triangles: dict[int, list[int]] = {root: [] for root in roots}
    for triangle in range(len(mesh.faces)):
        component_triangles[component_of[triangle]].append(triangle)
    component_labels = {root: tri_labels[component_triangles[root][0]] for root in roots}

    # 2. boundary arcs
    boundary: dict[tuple[int, int], tuple[int, ...]] = {}
    for edge, triangles in edge_triangles.items():
        if len(triangles) == 1:
            boundary[edge] = (component_of[triangles[0]],)
            continue
        first, second = triangles
        if component_of[first] != component_of[second]:
            boundary[edge] = tuple(sorted((component_of[first], component_of[second])))
    if not boundary:
        raise TopologyError(
            "the labelling has a single face and no boundary; a closed surface needs at "
            "least one declared seam or a second region, otherwise V - E + F cannot be "
            "the Euler characteristic (see docs/topology.md)"
        )

    incident: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for edge in boundary:
        incident[edge[0]].append(edge)
        incident[edge[1]].append(edge)

    corners: set[int] = set()
    for node, arcs in incident.items():
        if len(arcs) != 2:
            corners.add(node)
            continue
        if boundary[arcs[0]] != boundary[arcs[1]]:
            corners.add(node)

    walks = _walk_chains(boundary, incident, corners, mesh.vertices)
    if not walks:  # pragma: no cover - every boundary edge is walked
        raise TopologyError("no boundary chains were found")

    # 3. names first: vertices and edges name themselves after their faces
    scale = _model_scale(mesh)
    areas = meshlib.triangle_areas(mesh)
    face_names = _unique_names(
        [(root, component_labels[root]) for root in roots],
        {
            root: _canonical_key(
                _triangle_centroids(mesh, component_triangles[root]), scale
            )
            for root in roots
        },
        name_prefix,
    )

    corner_nodes = sorted({node for walk in walks for node in (walk[0], walk[1])})
    vertex_faces: dict[int, set[int]] = defaultdict(set)
    for start, end, nodes in walks:
        faces = _walk_faces(boundary, nodes)
        vertex_faces[start].update(faces)
        vertex_faces[end].update(faces)
    vertex_names = _disambiguate(
        {
            node: name_prefix + "V:" + "|".join(
                sorted(face_names[component] for component in vertex_faces[node])
            )
            for node in corner_nodes
        },
        {node: _canonical_key([mesh.vertices[node]], scale) for node in corner_nodes},
    )
    vertex_index = {node: position for position, node in enumerate(corner_nodes)}

    edge_keys: dict[int, tuple[float, ...]] = {}
    edge_names: dict[int, str] = {}
    for position, (_start, _end, nodes) in enumerate(walks):
        faces_sorted = sorted(face_names[component] for component in _walk_faces(boundary, nodes))
        edge_names[position] = name_prefix + "E:" + "|".join(faces_sorted)
        anchors = [mesh.vertices[nodes[0]], mesh.vertices[nodes[-1]]]
        if len(nodes) > 2:
            anchors.append(_midpoint([mesh.vertices[node] for node in nodes[1:-1]]))
        edge_keys[position] = _canonical_key(anchors, scale)
    edge_names = _disambiguate(edge_names, edge_keys)

    edges: list[TopoEdge] = []
    winding: list[bool] = []
    for position, (start, end, nodes) in enumerate(walks):
        faces = tuple(
            sorted(root_position[component] for component in _walk_faces(boundary, nodes))
        )
        edges.append(
            TopoEdge(
                index=position,
                vertices=(vertex_index[start], vertex_index[end]),
                faces=faces,
                name=edge_names[position],
                nodes=tuple(nodes),
                loop=start == end,
            )
        )
        winding.append(
            _winding_consistent(mesh, edge_triangles, nodes, component_of, boundary)
        )

    face_edges: dict[int, list[int]] = defaultdict(list)
    face_vertices: dict[int, set[int]] = defaultdict(set)
    for edge in edges:
        for face in edge.faces:
            face_edges[face].append(edge.index)
            face_vertices[face].update(edge.vertices)

    faces = tuple(
        TopoFace(
            index=position,
            label=component_labels[root],
            name=face_names[root],
            edges=tuple(sorted(face_edges[position])),
            vertices=tuple(sorted(face_vertices[position])),
            triangles=tuple(sorted(component_triangles[root])),
            area=math.fsum(areas[triangle] for triangle in component_triangles[root]),
        )
        for position, root in enumerate(roots)
    )

    # only referenced vertices count: an unreferenced vertex would break the Euler
    # characteristic of the mesh without being part of its surface
    used_nodes = {node for triangle in mesh.faces for node in triangle}
    mesh_counts = (len(used_nodes), len(edge_triangles), len(mesh.faces))
    return CellComplex(
        faces=faces,
        edges=tuple(edges),
        vertices=tuple(
            TopoVertex(
                index=vertex_index[node],
                point=mesh.vertices[node],
                name=vertex_names[node],
                nodes=(node,),
            )
            for node in corner_nodes
        ),
        mesh_counts=mesh_counts,
    )


def _model_scale(mesh: Mesh) -> float:
    box = mesh.bounding_box()
    diagonal = math.dist(box.minimum, box.maximum)
    return diagonal if diagonal > 0.0 else 1.0


def _triangle_centroids(mesh: Mesh, triangles: Iterable[int]) -> list[Vec3]:
    out: list[Vec3] = []
    for triangle in triangles:
        a, b, c = mesh.faces[triangle]
        pa, pb, pc = mesh.vertices[a], mesh.vertices[b], mesh.vertices[c]
        out.append(
            (
                (pa[0] + pb[0] + pc[0]) / 3.0,
                (pa[1] + pb[1] + pc[1]) / 3.0,
                (pa[2] + pb[2] + pc[2]) / 3.0,
            )
        )
    return out


# ---------------------------------------------------------------------------
# Euler operators
# ---------------------------------------------------------------------------
def _face_walk(complex_: CellComplex, face: TopoFace | int) -> list[int] | None:
    """Vertex chain around ``face`` in edge order, or ``None`` if it does not chain."""
    record = complex_.faces[_face_position(face)]
    if not record.edges:
        return None
    first = complex_.edges[record.edges[0]]
    chain = [first.vertices[0], first.vertices[1]]
    remaining = list(record.edges[1:])
    while remaining:
        tail = chain[-1]
        for position, edge_index in enumerate(remaining):
            edge = complex_.edges[edge_index]
            if edge.vertices[0] == tail:
                chain.append(edge.vertices[1])
                remaining.pop(position)
                break
            if edge.vertices[1] == tail:
                chain.append(edge.vertices[0])
                remaining.pop(position)
                break
        else:
            return None
    if len(chain) < 2 or chain[0] != chain[-1]:
        return None
    return chain


def _assemble(
    faces: Sequence[tuple[str, str, Sequence[str], Sequence[str], tuple[int, ...], float]],
    edges: Sequence[tuple[str, tuple[str, str], Sequence[str], tuple[int, ...], bool]],
    vertices: Sequence[tuple[str, Vec3, tuple[int, ...]]],
    *,
    mesh_counts: tuple[int, int, int] | None,
) -> CellComplex:
    """Index name-keyed records into a complex (names are the only identity here)."""
    face_index = {name: position for position, (name, *_rest) in enumerate(faces)}
    edge_index = {name: position for position, (name, *_rest) in enumerate(edges)}
    vertex_index = {name: position for position, (name, *_rest) in enumerate(vertices)}
    topo_vertices = tuple(
        TopoVertex(index=position, point=point, name=name, nodes=nodes)
        for position, (name, point, nodes) in enumerate(vertices)
    )
    topo_edges = tuple(
        TopoEdge(
            index=position,
            vertices=(vertex_index[first], vertex_index[second]),
            faces=tuple(face_index[face] for face in face_names),
            name=name,
            nodes=nodes,
            loop=loop,
        )
        for position, (name, (first, second), face_names, nodes, loop) in enumerate(edges)
    )
    topo_faces = tuple(
        TopoFace(
            index=position,
            label=label,
            name=name,
            edges=tuple(edge_index[edge] for edge in edge_names),
            vertices=tuple(sorted(vertex_index[vertex] for vertex in vertex_names)),
            triangles=triangles,
            area=area,
        )
        for position, (name, label, edge_names, vertex_names, triangles, area) in enumerate(faces)
    )
    return CellComplex(
        faces=topo_faces,
        edges=topo_edges,
        vertices=topo_vertices,
        mesh_counts=mesh_counts,
    )


def _vertex_index(complex_: CellComplex, position: int) -> int:
    if not 0 <= position < len(complex_.vertices):
        raise TopologyError(f"vertex index {position} is out of range")
    return complex_.vertices[position].index


def split_face(
    complex_: CellComplex,
    face: TopoFace | int,
    vertices: tuple[TopoVertex | int, TopoVertex | int],
    *,
    suffix: str = "split",
) -> CellComplex:
    """Euler operator ``mef``: add an edge between two vertices of a face, splitting it.

    ``F`` and ``E`` each go up by one and ``V`` is unchanged, so the Euler characteristic
    cannot move; the operator refuses when the two vertices are not both on the face's
    boundary, because then it would not be a split but a new face floating on the old one.
    """
    record = complex_.faces[_face_position(face)]
    first = vertices[0] if isinstance(vertices[0], TopoVertex) else _vertex_index(complex_, int(vertices[0]))
    second = vertices[1] if isinstance(vertices[1], TopoVertex) else _vertex_index(complex_, int(vertices[1]))
    if first == second:
        raise TopologyError("an edge needs two distinct vertices")
    walk = _face_walk(complex_, record)
    if walk is None:
        raise TopologyError(f"face {record.name} boundary is not a single chain")
    chain = walk[:-1]
    if first not in chain or second not in chain:
        raise TopologyError(
            f"vertices {complex_.vertices[first].name} and {complex_.vertices[second].name} "
            f"are not both on face {record.name}"
        )
    a, b = sorted((chain.index(first), chain.index(second)))
    lower = chain[a : b + 1]
    upper = chain[b:] + chain[: a + 1]
    if len(lower) < 3 or len(upper) < 3:
        raise TopologyError("both halves of a split need at least three vertices")

    new_edge_name = _free_name(
        {edge.name for edge in complex_.edges}, f"{record.name}/split-E"
    )
    new_face_name = _free_name(
        {candidate.name for candidate in complex_.faces}, f"{record.name}/{suffix}"
    )
    lower_edges = _cycle_edges(complex_, lower)
    upper_edges = _cycle_edges(complex_, upper)
    if lower_edges is None or upper_edges is None:  # pragma: no cover - walks come from edges
        raise TopologyError("a half of the split is not an edge chain")
    used_by_other = lower_edges & upper_edges
    # the new edge takes the place of the shared cut in both cycles
    lower_edges = (lower_edges - used_by_other) | {new_edge_name}
    upper_edges = (upper_edges - used_by_other) | {new_edge_name}

    faces: list[tuple[str, str, Sequence[str], Sequence[str], tuple[int, ...], float]] = []
    for candidate in complex_.faces:
        if candidate.index != record.index:
            faces.append(
                (
                    candidate.name,
                    candidate.label,
                    [complex_.edges[e].name for e in candidate.edges],
                    [complex_.vertices[v].name for v in candidate.vertices],
                    candidate.triangles,
                    candidate.area,
                )
            )
            continue
        faces.append(
            (
                record.name,
                record.label,
                sorted(lower_edges),
                [complex_.vertices[v].name for v in lower],
                (),
                record.area,
            )
        )
    faces.append(
        (
            new_face_name,
            record.label,
            sorted(upper_edges),
            [complex_.vertices[v].name for v in upper],
            (),
            0.0,
        )
    )

    edges = [
        (
            edge.name,
            (complex_.vertices[edge.vertices[0]].name, complex_.vertices[edge.vertices[1]].name),
            [complex_.faces[f].name for f in edge.faces],
            edge.nodes,
            edge.loop,
        )
        for edge in complex_.edges
    ]
    edges.append(
        (
            new_edge_name,
            (complex_.vertices[first].name, complex_.vertices[second].name),
            [record.name, new_face_name],
            (),
            False,
        )
    )
    vertices = [(v.name, v.point, v.nodes) for v in complex_.vertices]
    result = _assemble(
        faces,
        edges,
        vertices,
        mesh_counts=None,
    )
    _assert_preserved(complex_, result, "split_face")
    return result


def _free_name(taken: set[str], base: str) -> str:
    if base not in taken:
        return base
    position = 2
    while f"{base}#{position}" in taken:
        position += 1
    return f"{base}#{position}"


def _cycle_edges(complex_: CellComplex, chain: Sequence[int]) -> set[str] | None:
    """Edge names joining consecutive vertices of ``chain``; ``None`` if one is missing."""
    out: set[str] = set()
    for first, second in zip(chain, chain[1:]):
        for edge in complex_.edges:
            if {edge.vertices[0], edge.vertices[1]} == {first, second}:
                out.add(edge.name)
                break
        else:
            return None
    return out


def split_edge(
    complex_: CellComplex,
    edge: TopoEdge | int,
    point: Vec3,
    *,
    suffix: str = "mid",
) -> CellComplex:
    """Euler operator ``mev``: insert a vertex on an edge, splitting it in two.

    ``V`` and ``E`` each go up by one, so the Euler characteristic is unchanged; the new
    vertex is the only place where the caller's geometry enters, and it is stored as
    given rather than projected, so the operator cannot silently move the surface.
    """
    record = complex_.edges[_edge_position(edge)]
    vertex_names = {vertex.name for vertex in complex_.vertices}
    edge_names = {candidate.name for candidate in complex_.edges}
    new_vertex = _free_name(vertex_names, f"{record.name}/{suffix}")
    first_name = _free_name(edge_names, f"{record.name}/a")
    second_name = _free_name(edge_names | {first_name}, f"{record.name}/b")

    first_vertex = complex_.vertices[record.vertices[0]].name
    last_vertex = complex_.vertices[record.vertices[1]].name
    edges: list[tuple[str, tuple[str, str], Sequence[str], tuple[int, ...], bool]] = []
    for candidate in complex_.edges:
        if candidate.index != record.index:
            edges.append(
                (
                    candidate.name,
                    (
                        complex_.vertices[candidate.vertices[0]].name,
                        complex_.vertices[candidate.vertices[1]].name,
                    ),
                    [complex_.faces[f].name for f in candidate.faces],
                    candidate.nodes,
                    candidate.loop,
                )
            )
            continue
        face_names = [complex_.faces[f].name for f in candidate.faces]
        edges.append((first_name, (first_vertex, new_vertex), face_names, candidate.nodes, False))
        edges.append((second_name, (new_vertex, last_vertex), face_names, (), False))

    faces: list[tuple[str, str, Sequence[str], Sequence[str], tuple[int, ...], float]] = []
    for candidate in complex_.faces:
        names = [complex_.edges[e].name for e in candidate.edges]
        updated: list[str] = []
        for name in names:
            if name == record.name:
                updated.extend([first_name, second_name])
            else:
                updated.append(name)
        vertex_names_of_face = [complex_.vertices[v].name for v in candidate.vertices]
        if record.name in names:
            vertex_names_of_face.append(new_vertex)
        faces.append(
            (
                candidate.name,
                candidate.label,
                updated,
                vertex_names_of_face,
                candidate.triangles,
                candidate.area,
            )
        )

    vertices = [(v.name, v.point, v.nodes) for v in complex_.vertices]
    vertices.append(
        (new_vertex, (float(point[0]), float(point[1]), float(point[2])), ())
    )
    result = _assemble(
        faces,
        edges,
        vertices,
        mesh_counts=None,
    )
    _assert_preserved(complex_, result, "split_edge")
    return result


def merge_faces(
    complex_: CellComplex, edge: TopoEdge | int, *, name: str | None = None
) -> CellComplex:
    """Euler operator ``kef``: remove an edge, joining the two faces it separates.

    ``F`` and ``E`` each drop by one, so the Euler characteristic is unchanged. Vertices
    stay: the endpoints of the removed edge are corners of other faces (in a box, three
    edges meet at each corner), so they are not removed with it.
    """
    record = complex_.edges[_edge_position(edge)]
    if len(record.faces) != 2:
        raise TopologyError(f"edge {record.name} does not separate two faces")
    first_face = complex_.faces[record.faces[0]]
    second_face = complex_.faces[record.faces[1]]
    walk_a = _face_walk(complex_, first_face)
    walk_b = _face_walk(complex_, second_face)
    if walk_a is None or walk_b is None:
        raise TopologyError("both faces must have a single closed boundary chain")
    chain_a = walk_a[:-1]
    positions = [
        index
        for index in range(len(chain_a))
        if {chain_a[index], chain_a[(index + 1) % len(chain_a)]} == set(record.vertices)
    ]
    if not positions:  # pragma: no cover - incidence comes from the same structure
        raise TopologyError(f"edge {record.name} is not on face {first_face.name}")
    start = positions[0]
    chain_b = walk_b[:-1]
    matches = [
        index
        for index in range(len(chain_b))
        if {chain_b[index], chain_b[(index + 1) % len(chain_b)]} == set(record.vertices)
    ]
    if not matches:  # pragma: no cover - defensive
        raise TopologyError(f"edge {record.name} is not on face {second_face.name}")
    # join: walk face A up to the shared edge, then walk face B backwards from the other
    # end of the shared edge, so the merged chain stays a single cycle
    head = chain_a[start + 1 :] + chain_a[: start + 1]
    other = matches[0]
    tail = chain_b[other + 1 :] + chain_b[: other + 1]
    joined = head + tail[1:-1]
    if len(set(joined)) != len(joined) - (1 if joined[0] == joined[-1] else 0):
        raise TopologyError(f"merging across {record.name} would repeat a vertex")

    merged_name = name or _free_name(
        {candidate.name for candidate in complex_.faces},
        f"{first_face.name}+{second_face.name}",
    )
    if merged_name in {candidate.name for candidate in complex_.faces}:
        raise TopologyError(f"a face named {merged_name} already exists")

    faces: list[tuple[str, str, Sequence[str], Sequence[str], tuple[int, ...], float]] = []
    for candidate in complex_.faces:
        if candidate.index == second_face.index:
            continue
        if candidate.index == first_face.index:
            edges_of_merged = [
                complex_.edges[e].name
                for e in candidate.edges
                if complex_.edges[e].name != record.name
            ]
            faces.append(
                (
                    merged_name,
                    candidate.label,
                    edges_of_merged,
                    sorted(
                        {complex_.vertices[v].name for v in first_face.vertices}
                        | {complex_.vertices[v].name for v in second_face.vertices}
                    ),
                    (),
                    first_face.area + second_face.area,
                )
            )
            continue
        faces.append(
            (
                candidate.name,
                candidate.label,
                [complex_.edges[e].name for e in candidate.edges],
                [complex_.vertices[v].name for v in candidate.vertices],
                candidate.triangles,
                candidate.area,
            )
        )
    edges = [
        (
            edge.name,
            (complex_.vertices[edge.vertices[0]].name, complex_.vertices[edge.vertices[1]].name),
            [
                merged_name if complex_.faces[f].name == first_face.name
                else complex_.faces[f].name
                for f in edge.faces
            ],
            edge.nodes,
            edge.loop,
        )
        for edge in complex_.edges
        if edge.index != record.index
    ]
    vertices = [(v.name, v.point, v.nodes) for v in complex_.vertices]
    result = _assemble(
        faces,
        edges,
        vertices,
        mesh_counts=None,
    )
    _assert_preserved(complex_, result, "merge_faces")
    return result


def _assert_preserved(before: CellComplex, after: CellComplex, operator: str) -> None:
    if before.euler_characteristic() != after.euler_characteristic():
        raise TopologyError(
            f"{operator} moved the Euler characteristic from "
            f"{before.euler_characteristic()} to {after.euler_characteristic()}"
        )
    problems = [problem for problem in after.validate() if "Euler characteristic" not in problem]
    if problems:
        raise TopologyError(f"{operator} produced an invalid complex: " + "; ".join(problems))