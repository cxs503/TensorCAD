"""Hull geometry quality audit: watertightness, fairness, hydrostatic consistency.

A hull that renders correctly is not a hull a solver or a hydrostatics table can
use.  This module is the machine gate between an exchanged model and the numbers
computed from it, ordered the way a shipyard checks a lines plan:

1. **Topology** (:func:`audit_watertightness`).  Free edges, seam gaps,
   non-manifold edges, duplicate facets, winding conflicts and a signed volume of
   the wrong sign.  Every finding is a :class:`GeometryDefect` with a magnitude in
   metres *and* as a fraction of ``L``, because "a 2 mm gap" is a defect on a 6 m
   model and is not one on a 300 m ship.
2. **Fairness** (:func:`audit_fairness`).  Discrete Gaussian curvature (angle
   deficit over mixed area) and the relative jump across every smooth interior
   edge.  A jump above :data:`CURVATURE_JUMP_THRESHOLD` (50 %) between neighbours
   is a curvature discontinuity: a dent, an unintended kink, a bad patch.  Designed
   creases (deck edge, transom knuckle, stem) are *skipped* by dihedral angle, not
   smoothed away, so a knuckle is not reported as a defect.
3. **Declared geometry** (:func:`audit_declared_geometry`).  Realised ``L/B/T``
   and ``Cb`` against the values the caller declared.  This runs *before* any
   hydrostatics: if the model is not the ship it claims to be, the displacement it
   integrates to is not evidence about that ship.
4. **Hydrostatic consistency** (:func:`audit_hydrostatic_consistency`).  The same
   hull is integrated twice by independent routes -- an offsets table (lines plan:
   half-breadths on a station x waterline grid, Simpson/trapezoid) and the triangle
   mesh (divergence theorem) -- and the differences in displacement volume, wetted
   area and ``Cb`` must stay inside :data:`HYDROSTATIC_TOLERANCE` (1 %).  A table
   built from *analytic* half-breadths is independent of the mesh; one built from
   the mesh itself (:func:`offsets_from_mesh`) shares the vertices but not the
   quadrature, which is the ordinary "lines plan versus 3-D model" check.

Negative cases are part of the contract, not an afterthought, so the module also
carries the injectors it is tested with: :func:`inject_seam_offset`,
:func:`inject_hole`, :func:`inject_flip_faces`, :func:`inject_duplicate_faces`,
:func:`inject_cb_offset` and :func:`inject_bump`.  They are small, they are the
only way to know the audit can see anything, and a pipeline that returns
``accepted`` for every one of them is worthless.

Gaussian curvature here is the *intrinsic* curvature, so it is zero on a
developable plate and large at a bilge: that is why the jump criterion is
relative, and why the absolute curvature is reported alongside it.
"""
from __future__ import annotations

import json
import math
import random
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Sequence

from . import mesh as meshlib
from .geometry import Mesh
from .hydrostatics import Hydrostatics, compute_hydrostatics

__all__ = [
    "CREASE_DEGREES",
    "CURVATURE_JUMP_THRESHOLD",
    "DECLARED_CB_TOLERANCE",
    "DIMENSION_TOLERANCE",
    "GAP_TOLERANCE_FRACTION",
    "HYDROSTATIC_TOLERANCE",
    "ConsistencyAudit",
    "CurvatureField",
    "CurvatureJumpRegion",
    "FairnessAudit",
    "GeometryAuditReport",
    "GeometryDefect",
    "HullSpec",
    "OffsetTable",
    "TableHydrostatics",
    "WatertightnessAudit",
    "audit_declared_geometry",
    "audit_fairness",
    "audit_hull_geometry",
    "audit_hydrostatic_consistency",
    "audit_watertightness",
    "gaussian_curvature",
    "inject_bump",
    "inject_cb_offset",
    "inject_duplicate_faces",
    "inject_flip_faces",
    "inject_hole",
    "inject_seam_offset",
    "integrate_offsets",
    "offsets_from_function",
    "offsets_from_mesh",
    "write_gates_json",
]

SCHEMA_AUDIT = "tensorcad.geometry_audit/1"
SCHEMA_GATES = "tensorcad.geometry_audit.gates/1"

#: A free edge further from the rest of the surface than this fraction of ``L`` is
#: a gap worth reporting: 0.1 % of the length (2 mm on a 2 m model, 30 mm on a
#: 30 m hull).  This is the acceptance figure the pipeline gates on.
GAP_TOLERANCE_FRACTION = 1e-3

#: Relative Gaussian-curvature jump ``|K_u - K_v| / mean(|K_u|,|K_v|)`` that marks
#: an edge as a curvature discontinuity.
CURVATURE_JUMP_THRESHOLD = 0.5

#: Edges whose dihedral angle exceeds this are *designed* creases (deck edge,
#: transom knuckle, stem).  Curvature jumps across them are expected and are not
#: defects; they are counted separately so the exclusion stays visible.
CREASE_DEGREES = 45.0

#: Number of edge rings around a designed crease that are absorbed by it.  A pinch
#: (the Series 60 ends taper to a line) produces a genuine, smooth but steep
#: curvature gradient next to the crease; without this dilation the field there
#: would be scored as an unfairness that no finer mesh can remove.
CREASE_DILATION_RINGS = 2

#: A curvature-jump region smaller than this fraction of the surface is treated as
#: discretisation noise rather than a defect.  1e-4 of the surface is roughly one
#: or two triangles of a 160x40 grid.
MIN_REGION_AREA_FRACTION = 1e-4

#: Allowed relative difference between the offsets-table and mesh integration.
HYDROSTATIC_TOLERANCE = 0.01

#: Allowed relative difference between realised and declared ``L/B/T``.
DIMENSION_TOLERANCE = 0.01

#: Allowed relative difference between realised and declared ``Cb``.  Wider than
#: the dimension band because ``Cb`` also carries the volume discretisation error.
DECLARED_CB_TOLERANCE = 0.02


# --------------------------------------------------------------------------- #
# defects and the watertightness audit
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class GeometryDefect:
    """One finding, with a size a reviewer can judge."""

    kind: str
    magnitude: float
    unit: str
    where: tuple[float, float, float] | None = None
    detail: str = ""
    #: ``magnitude`` divided by the reference length, when one was supplied.
    relative: float | None = None

    def to_dict(self) -> dict[str, Any]:
        data: dict[str, Any] = {
            "kind": self.kind,
            "magnitude": round(self.magnitude, 9),
            "unit": self.unit,
            "detail": self.detail,
        }
        if self.where is not None:
            data["where"] = [round(value, 6) for value in self.where]
        if self.relative is not None:
            data["relative"] = round(self.relative, 9)
        return data


@dataclass(frozen=True)
class WatertightnessAudit:
    """Topology verdict for one surface."""

    open_edges: int = 0
    boundary_components: int = 0
    max_gap: float = 0.0
    gap_tolerance: float = 0.0
    nonmanifold_edges: int = 0
    duplicate_faces: int = 0
    degenerate_faces: int = 0
    winding_conflicts: int = 0
    conflicting_area: float = 0.0
    signed_volume: float = 0.0
    closed: bool = False
    defects: tuple[GeometryDefect, ...] = ()
    length: float | None = None
    patch_seams: int = 0
    face_components: int = 1

    @property
    def passed(self) -> bool:
        return not self.defects

    @property
    def max_gap_relative(self) -> float | None:
        if not self.length:
            return None
        return self.max_gap / float(self.length)

    def to_dict(self) -> dict[str, Any]:
        return {
            "passed": self.passed,
            "closed": self.closed,
            "open_edges": self.open_edges,
            "boundary_components": self.boundary_components,
            "max_gap_m": round(self.max_gap, 9),
            "gap_tolerance_m": round(self.gap_tolerance, 9),
            "max_gap_relative_to_length": None
            if self.max_gap_relative is None
            else round(self.max_gap_relative, 9),
            "nonmanifold_edges": self.nonmanifold_edges,
            "duplicate_faces": self.duplicate_faces,
            "degenerate_faces": self.degenerate_faces,
            "winding_conflicts": self.winding_conflicts,
            "conflicting_area_m2": round(self.conflicting_area, 9),
            "signed_volume_m3": round(self.signed_volume, 9),
            "patch_seams": self.patch_seams,
            "face_components": self.face_components,
            "defects": [defect.to_dict() for defect in self.defects],
        }


def _distance(a: Sequence[float], b: Sequence[float]) -> float:
    return math.dist(tuple(a), tuple(b))


def _boundary_components(mesh: Mesh) -> list[tuple[tuple[int, int], ...]]:
    """Group free (one-face) edges into chains/loops of consecutive free edges."""
    free = [edge for edge, count in meshlib.edge_face_counts(mesh).items() if count == 1]
    if not free:
        return []
    adjacency: dict[int, list[int]] = defaultdict(list)
    for index, (u, v) in enumerate(free):
        adjacency[u].append(index)
        adjacency[v].append(index)
    seen: set[int] = set()
    components: list[tuple[tuple[int, int], ...]] = []
    for start in range(len(free)):
        if start in seen:
            continue
        stack = [start]
        seen.add(start)
        group: list[tuple[int, int]] = []
        while stack:
            index = stack.pop()
            group.append(free[index])
            for vertex in free[index]:
                for neighbour in adjacency[vertex]:
                    if neighbour not in seen:
                        seen.add(neighbour)
                        stack.append(neighbour)
        components.append(tuple(group))
    return components


def _component_vertices(component: Sequence[tuple[int, int]]) -> list[int]:
    return list(dict.fromkeys(vertex for edge in component for vertex in edge))


def _component_opening(
    mesh: Mesh, component: Sequence[tuple[int, int]]
) -> tuple[float, tuple[float, float, float]]:
    """Opening size of one boundary component, and its centroid.

    A closed ring's opening is the largest chord between vertices at least a third
    of the ring apart; an open chain's opening is the distance between its two
    ends.  Both answer the same question -- *how big is the hole* -- without
    needing correspondence between the two sides of a gap.
    """
    vertices = _component_vertices(component)
    points = [mesh.vertices[index] for index in vertices]
    n = len(vertices)
    centroid = (
        sum(point[0] for point in points) / n,
        sum(point[1] for point in points) / n,
        sum(point[2] for point in points) / n,
    )
    if n == 1:
        return 0.0, centroid
    if n == 2:
        return _distance(points[0], points[1]), centroid
    # Ring or chain: decided by whether the two ends of the free-edge path coincide.
    ends_coincide = _distance(points[0], points[-1]) <= 1e-9 + 1e-6 * max(
        _distance(points[0], points[1]), _distance(points[-1], points[-2])
    )
    if not ends_coincide:
        return _distance(points[0], points[-1]), centroid
    stride = max(1, n // 3)
    best = 0.0
    for i in range(n):
        for j in range(i + stride, min(i + n - stride + 1, n)):
            best = max(best, _distance(points[i], points[j % n]))
    return best, centroid


def _min_component_distance(
    mesh: Mesh, left: Sequence[tuple[int, int]], right: Sequence[tuple[int, int]]
) -> float:
    left_vertices = _component_vertices(left)
    right_vertices = _component_vertices(right)
    if len(left_vertices) > 200:
        left_vertices = left_vertices[:: max(1, len(left_vertices) // 200)]
    if len(right_vertices) > 200:
        right_vertices = right_vertices[:: max(1, len(right_vertices) // 200)]
    best = math.inf
    for i in left_vertices:
        a = mesh.vertices[i]
        for j in right_vertices:
            best = min(best, _distance(a, mesh.vertices[j]))
    return best


def _winding_conflicts(mesh: Mesh) -> tuple[int, float]:
    """Edges whose two incident faces traverse them in the *same* direction."""
    directions: dict[tuple[int, int], list[tuple[int, int]]] = defaultdict(list)
    for face_index, face in enumerate(mesh.faces):
        for u, v in meshlib.face_edges(face):
            directions[meshlib.undirected_edge((u, v))].append((face_index, (u, v)))
    conflicts = 0
    flagged_faces: dict[int, float] = {}
    for entries in directions.values():
        if len(entries) != 2:
            continue
        (face_a, (a0, a1)), (face_b, (b0, b1)) = entries
        if a0 == b0 and a1 == b1:
            conflicts += 1
            for index in (face_a, face_b):
                flagged_faces[index] = meshlib.triangle_area(mesh, mesh.faces[index])
    return conflicts, math.fsum(flagged_faces.values())


def _face_components(mesh: Mesh) -> int:
    """Number of connected components over the faces (vertices as the link)."""
    parent: dict[int, int] = {}

    def find(vertex: int) -> int:
        root = vertex
        while parent[root] != root:
            root = parent[root]
        while parent[vertex] != root:
            parent[vertex], vertex = root, parent[vertex]
        return root

    for face in mesh.faces:
        for vertex in face:
            parent.setdefault(vertex, vertex)
        for u, v in ((face[0], face[1]), (face[1], face[2])):
            root_u, root_v = find(u), find(v)
            if root_u != root_v:
                parent[root_u] = root_v
    return len({find(vertex) for vertex in parent})


def audit_watertightness(
    mesh: Mesh,
    *,
    length: float | None = None,
    gap_tolerance_fraction: float = GAP_TOLERANCE_FRACTION,
    expected_components: int = 1,
    patch_ids: Sequence[str] | None = None,
) -> WatertightnessAudit:
    """Check that the surface can bound a volume, and measure how badly it fails."""
    if not mesh.faces:
        raise ValueError("cannot audit a mesh with no faces")
    box = mesh.bounding_box()
    reference_length = float(length) if length is not None else float(box.sizes[0])
    if reference_length <= 0:
        raise ValueError(f"reference length must be positive, got {reference_length!r}")
    tolerance = gap_tolerance_fraction * reference_length

    counts = meshlib.edge_face_counts(mesh)
    open_edges = sum(1 for count in counts.values() if count == 1)
    nonmanifold = sum(1 for count in counts.values() if count > 2)
    degenerate = len(meshlib.zero_area_faces(mesh, tolerance=0.0))
    seen: dict[tuple[int, int, int], int] = defaultdict(int)
    for face in mesh.faces:
        seen[tuple(sorted(face))] += 1
    duplicates = sum(count - 1 for count in seen.values() if count > 1)
    conflicts, conflicting_area = _winding_conflicts(mesh)
    signed_volume = meshlib.total_signed_volume(mesh)

    defects: list[GeometryDefect] = []
    components = _boundary_components(mesh)
    openings = [_component_opening(mesh, component) for component in components]
    max_gap = 0.0
    seams = 0
    for index, (opening, where) in enumerate(openings):
        nearest_other = math.inf
        for other_index, other in enumerate(components):
            if other_index == index:
                continue
            nearest_other = min(
                nearest_other, _min_component_distance(mesh, components[index], other)
            )
        if math.isfinite(nearest_other):
            seams += 1
            # The separation between pieces *is* the gap: it is the distance the
            # surface has to move to close.  The component's own opening (ring span)
            # only matters when the piece is alone.
            gap = nearest_other
            kind = "seam_gap"
            detail = (
                f"surface splits into {len(components)} pieces; the nearest piece is "
                f"{nearest_other:.6g} m away while this component's own opening is "
                f"{opening:.6g} m"
            )
        else:
            gap = opening
            kind = "open_boundary"
            detail = (
                f"single boundary component with "
                f"{len(_component_vertices(components[index]))} free-edge vertices"
            )
        max_gap = max(max_gap, gap)
        if gap > tolerance:
            defects.append(
                GeometryDefect(
                    kind=kind,
                    magnitude=gap,
                    unit="m",
                    where=where,
                    detail=detail,
                    relative=gap / reference_length,
                )
            )
    if nonmanifold:
        defects.append(
            GeometryDefect(
                kind="nonmanifold_edge",
                magnitude=float(nonmanifold),
                unit="count",
                detail="an edge is shared by more than two faces: overlap or pinch",
            )
        )
    if duplicates:
        defects.append(
            GeometryDefect(
                kind="duplicate_facets",
                magnitude=float(duplicates),
                unit="count",
                detail="identical facets repeated: overlapping shells inside the solid",
            )
        )
    if degenerate:
        defects.append(
            GeometryDefect(
                kind="degenerate_facets",
                magnitude=float(degenerate),
                unit="count",
                detail="zero-area facets carry a normal but no area",
            )
        )
    if conflicts:
        defects.append(
            GeometryDefect(
                kind="normal_inconsistency",
                magnitude=float(conflicts),
                unit="count",
                detail=(
                    f"{conflicts} edges are traversed in the same direction by both incident "
                    f"faces ({conflicting_area:.6g} m^2 of flipped surface): the normals are "
                    "not consistent"
                ),
            )
        )
    face_components = _face_components(mesh)
    if face_components != expected_components:
        defects.append(
            GeometryDefect(
                kind="component_count",
                magnitude=float(face_components),
                unit="count",
                detail=(
                    f"{face_components} disconnected shells for an expected "
                    f"{expected_components}: a detached or duplicated solid is not one hull"
                ),
                relative=float(face_components - expected_components),
            )
        )
    if components and signed_volume <= 0:
        defects.append(
            GeometryDefect(
                kind="inverted_orientation",
                magnitude=signed_volume,
                unit="m^3",
                detail="signed volume is not positive: the surface is oriented inward",
            )
        )
    return WatertightnessAudit(
        open_edges=open_edges,
        boundary_components=len(components),
        max_gap=max_gap,
        gap_tolerance=tolerance,
        nonmanifold_edges=nonmanifold,
        duplicate_faces=duplicates,
        degenerate_faces=degenerate,
        winding_conflicts=conflicts,
        conflicting_area=conflicting_area,
        signed_volume=signed_volume,
        closed=not open_edges and not nonmanifold,
        defects=tuple(defects),
        length=reference_length,
        patch_seams=seams,
        face_components=face_components,
    )


# --------------------------------------------------------------------------- #
# fairness: discrete Gaussian curvature and its jumps
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class CurvatureField:
    """Per-vertex discrete Gaussian curvature (angle deficit over mixed area)."""

    curvature: tuple[float, ...]
    mixed_area: tuple[float, ...]
    interior: tuple[bool, ...]
    normals: tuple[tuple[float, float, float], ...]
    mean_abs_curvature: float

    def curvature_at(self, vertex: int) -> float:
        return self.curvature[vertex]

    def to_dict(self) -> dict[str, Any]:
        return {
            "vertex_count": len(self.curvature),
            "interior_vertices": sum(1 for flag in self.interior if flag),
            "mean_abs_curvature_1_per_m2": round(self.mean_abs_curvature, 9),
        }


def _angles_and_mixed_areas(
    a: Sequence[float], b: Sequence[float], c: Sequence[float]
) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    """Corner angles and Meyer-style mixed-area contributions for one triangle."""

    def squared(p: Sequence[float], q: Sequence[float]) -> float:
        return (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 + (p[2] - q[2]) ** 2

    def angle(p: Sequence[float], q: Sequence[float], r: Sequence[float]) -> float:
        ux, uy, uz = (q[0] - p[0], q[1] - p[1], q[2] - p[2])
        vx, vy, vz = (r[0] - p[0], r[1] - p[1], r[2] - p[2])
        cross = (uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx)
        return math.atan2(
            math.sqrt(cross[0] ** 2 + cross[1] ** 2 + cross[2] ** 2),
            ux * vx + uy * vy + uz * vz,
        )

    alpha = angle(a, b, c)
    beta = angle(b, c, a)
    gamma = angle(c, a, b)
    lab = squared(a, b)
    lbc = squared(b, c)
    lca = squared(c, a)
    area = 0.5 * math.sqrt(max(lab * lca - ((lab + lca - lbc) / 2.0) ** 2, 0.0))
    right = 0.5 * math.pi
    if alpha < right and beta < right and gamma < right:
        mixed_a = (1.0 / 8.0) * (lab / math.tan(gamma) + lca / math.tan(beta))
        mixed_b = (1.0 / 8.0) * (lab / math.tan(gamma) + lbc / math.tan(alpha))
        mixed_c = (1.0 / 8.0) * (lca / math.tan(beta) + lbc / math.tan(alpha))
    elif alpha >= right:
        mixed_a, mixed_b, mixed_c = 0.5 * area, 0.25 * area, 0.25 * area
    elif beta >= right:
        mixed_a, mixed_b, mixed_c = 0.25 * area, 0.5 * area, 0.25 * area
    else:
        mixed_a, mixed_b, mixed_c = 0.25 * area, 0.25 * area, 0.5 * area
    return (alpha, beta, gamma), (mixed_a, mixed_b, mixed_c)


def gaussian_curvature(mesh: Mesh) -> CurvatureField:
    """Discrete Gaussian curvature ``K = (2*pi - sum(theta)) / A_mixed`` per vertex."""
    count = len(mesh.vertices)
    angle_sums = [0.0] * count
    mixed = [0.0] * count
    normals = [[0.0, 0.0, 0.0] for _ in range(count)]
    for face in mesh.faces:
        a, b, c = (mesh.vertices[index] for index in face)
        (alpha, beta, gamma), (mixed_a, mixed_b, mixed_c) = _angles_and_mixed_areas(a, b, c)
        angle_sums[face[0]] += alpha
        angle_sums[face[1]] += beta
        angle_sums[face[2]] += gamma
        mixed[face[0]] += mixed_a
        mixed[face[1]] += mixed_b
        mixed[face[2]] += mixed_c
        nx, ny, nz = meshlib.face_normal(mesh, face)
        for index in face:
            normals[index][0] += nx
            normals[index][1] += ny
            normals[index][2] += nz
    boundary = {vertex for edge in meshlib.boundary_edges(mesh) for vertex in edge}
    curvature = [float("nan")] * count
    interior = [False] * count
    for index in range(count):
        if index in boundary or mixed[index] <= 0.0:
            continue
        curvature[index] = (2.0 * math.pi - angle_sums[index]) / mixed[index]
        interior[index] = True
    unit_normals: list[tuple[float, float, float]] = []
    for vector in normals:
        magnitude = math.sqrt(vector[0] ** 2 + vector[1] ** 2 + vector[2] ** 2)
        if magnitude <= 0.0:
            unit_normals.append((0.0, 0.0, 1.0))
        else:
            unit_normals.append((vector[0] / magnitude, vector[1] / magnitude, vector[2] / magnitude))
    values = [value for value, flag in zip(curvature, interior) if flag and math.isfinite(value)]
    mean_abs = math.fsum(abs(value) for value in values) / len(values) if values else 0.0
    return CurvatureField(
        curvature=tuple(curvature),
        mixed_area=tuple(mixed),
        interior=tuple(interior),
        normals=tuple(unit_normals),
        mean_abs_curvature=mean_abs,
    )


def _unit(vector: Sequence[float]) -> tuple[float, float, float]:
    magnitude = math.sqrt(sum(component * component for component in vector))
    if magnitude <= 0.0:
        return (0.0, 0.0, 1.0)
    return (vector[0] / magnitude, vector[1] / magnitude, vector[2] / magnitude)


def _edge_face_table(mesh: Mesh) -> dict[tuple[int, int], list[int]]:
    """Edge -> incident face indices (the two-face case is what the crease check uses)."""
    table: dict[tuple[int, int], list[int]] = defaultdict(list)
    for face_index, face in enumerate(mesh.faces):
        for u, v in meshlib.face_edges(face):
            table[meshlib.undirected_edge((u, v))].append(face_index)
    return table


def dihedral_angle_degrees(
    mesh: Mesh,
    edge: tuple[int, int],
    *,
    normals: Sequence[Sequence[float]] | None = None,
    incident: Sequence[int] | None = None,
) -> float:
    """Angle between the two faces sharing ``edge`` (0 = flat, 180 = folded back).

    ``normals``/``incident`` let a caller that already built the edge table pass it
    in; without them this is an O(F) scan, which is fine for one edge and hopeless
    for a whole mesh.
    """
    if incident is None:
        faces = [face for face in mesh.faces if edge[0] in face and edge[1] in face]
        incident = list(range(len(faces)))
        face_list = faces
    else:
        face_list = None  # type: ignore[assignment]
    if len(incident) != 2:
        return 0.0
    if normals is None:
        normals = [
            _unit(meshlib.face_normal(mesh, mesh.faces[index])) for index in incident
        ]
        unit_a, unit_b = normals[0], normals[1]
    else:
        unit_a, unit_b = normals[incident[0]], normals[incident[1]]
    cosine = sum(a * b for a, b in zip(unit_a, unit_b))
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


@dataclass(frozen=True)
class CurvatureJumpRegion:
    """A connected group of vertices where the curvature field breaks."""

    vertex_count: int
    area: float
    area_fraction: float
    centroid: tuple[float, float, float]
    max_ratio: float
    curvature_min: float
    curvature_max: float
    count_as_defect: bool

    def to_dict(self) -> dict[str, Any]:
        return {
            "vertex_count": self.vertex_count,
            "area_m2": round(self.area, 9),
            "area_fraction": round(self.area_fraction, 9),
            "centroid": [round(value, 6) for value in self.centroid],
            "max_jump_ratio": round(self.max_ratio, 6),
            "gaussian_curvature_min_1_per_m2": round(self.curvature_min, 9),
            "gaussian_curvature_max_1_per_m2": round(self.curvature_max, 9),
            "counts_as_defect": self.count_as_defect,
        }


@dataclass(frozen=True)
class FairnessAudit:
    """Curvature-continuity verdict for one surface."""

    threshold: float = CURVATURE_JUMP_THRESHOLD
    vertex_count: int = 0
    interior_vertices: int = 0
    edges_considered: int = 0
    crease_edges_skipped: int = 0
    crease_edges: int = 0
    crease_vertices: int = 0
    mean_abs_curvature: float = 0.0
    max_jump_ratio: float = 0.0
    flagged_edges: int = 0
    flagged_vertices: int = 0
    flagged_area_fraction: float = 0.0
    regions: tuple[CurvatureJumpRegion, ...] = ()

    @property
    def passed(self) -> bool:
        return not any(region.count_as_defect for region in self.regions)

    @property
    def defects(self) -> tuple[GeometryDefect, ...]:
        findings: list[GeometryDefect] = []
        for region in self.regions:
            if not region.count_as_defect:
                continue
            findings.append(
                GeometryDefect(
                    kind="curvature_discontinuity",
                    magnitude=region.max_ratio,
                    unit="relative_jump",
                    where=region.centroid,
                    detail=(
                        f"{region.vertex_count} vertices over {region.area:.6g} m^2 "
                        f"({100.0 * region.area_fraction:.3f}% of the surface) carry a Gaussian "
                        f"curvature jump above {100.0 * self.threshold:.0f}%"
                    ),
                    relative=region.area_fraction,
                )
            )
        return tuple(findings)

    def to_dict(self) -> dict[str, Any]:
        return {
            "passed": self.passed,
            "threshold": self.threshold,
            "vertex_count": self.vertex_count,
            "interior_vertices": self.interior_vertices,
            "edges_considered": self.edges_considered,
            "crease_edges_skipped": self.crease_edges_skipped,
            "crease_edges": self.crease_edges,
            "crease_vertices": self.crease_vertices,
            "mean_abs_gaussian_curvature_1_per_m2": round(self.mean_abs_curvature, 9),
            "max_jump_ratio": round(self.max_jump_ratio, 6),
            "flagged_edges": self.flagged_edges,
            "flagged_vertices": self.flagged_vertices,
            "flagged_area_fraction": round(self.flagged_area_fraction, 9),
            "regions": [region.to_dict() for region in self.regions],
        }


def _dilate(mesh: Mesh, seed: set[int], rings: int) -> set[int]:
    """Graph distance (over mesh edges) at most ``rings`` from ``seed``."""
    if not seed or rings <= 0:
        return set(seed)
    adjacency: dict[int, list[int]] = defaultdict(list)
    for edge in meshlib.edge_face_counts(mesh):
        u, v = edge
        adjacency[u].append(v)
        adjacency[v].append(u)
    reached = set(seed)
    frontier = set(seed)
    for _ in range(int(rings)):
        nxt: set[int] = set()
        for vertex in frontier:
            for neighbour in adjacency[vertex]:
                if neighbour not in reached:
                    nxt.add(neighbour)
        reached |= nxt
        frontier = nxt
        if not frontier:
            break
    return reached


def audit_fairness(
    mesh: Mesh,
    *,
    threshold: float = CURVATURE_JUMP_THRESHOLD,
    crease_degrees: float = CREASE_DEGREES,
    min_region_area_fraction: float = MIN_REGION_AREA_FRACTION,
    crease_dilation: int = CREASE_DILATION_RINGS,
    field: CurvatureField | None = None,
) -> FairnessAudit:
    """Flag edges whose relative Gaussian-curvature jump exceeds ``threshold``.

    Designed creases (dihedral above ``crease_degrees``) are excluded from the
    comparison and counted in ``crease_edges_skipped``: the deck edge of a closed
    solid is a genuine curvature break and is not a defect.  The denominator is
    floored at a tenth of the mean absolute curvature so that a nearly developable
    patch does not report an unbounded jump from one triangle of noise.
    """
    curvature = field or gaussian_curvature(mesh)
    floor = 0.1 * curvature.mean_abs_curvature
    if floor <= 0.0:
        floor = 1e-12
    normals = [_unit(meshlib.face_normal(mesh, face)) for face in mesh.faces]
    edge_faces = _edge_face_table(mesh)
    # Gaussian curvature is not defined *at* a designed crease: the angle deficit
    # there mixes two surfaces (hull and deck cap, or the two halves of a
    # zero-breadth stem).  Those vertices are excluded from the comparison and
    # reported, rather than being allowed to dominate the flag list.
    crease_vertices: set[int] = set()
    crease_edges = 0
    for edge, incident in edge_faces.items():
        if len(incident) != 2:
            continue
        if (
            dihedral_angle_degrees(mesh, edge, normals=normals, incident=incident)
            > crease_degrees
        ):
            crease_edges += 1
            crease_vertices.update(edge)
    excluded_vertices = _dilate(mesh, crease_vertices, crease_dilation)
    flagged_pairs: list[tuple[int, int, float]] = []
    considered = 0
    skipped = 0
    max_ratio = 0.0
    for edge, incident in edge_faces.items():
        if len(incident) != 2:
            continue
        u, v = edge
        if not (curvature.interior[u] and curvature.interior[v]):
            continue
        ku, kv = curvature.curvature[u], curvature.curvature[v]
        if not (math.isfinite(ku) and math.isfinite(kv)):
            continue
        if u in excluded_vertices or v in excluded_vertices:
            skipped += 1
            continue
        considered += 1
        denominator = max(0.5 * (abs(ku) + abs(kv)), floor)
        ratio = abs(ku - kv) / denominator
        max_ratio = max(max_ratio, ratio)
        if ratio > threshold:
            flagged_pairs.append((u, v, ratio))

    flagged_vertices = sorted({vertex for u, v, _ in flagged_pairs for vertex in (u, v)})
    total_area = math.fsum(curvature.mixed_area)

    parent: dict[int, int] = {vertex: vertex for vertex in flagged_vertices}

    def find(vertex: int) -> int:
        root = vertex
        while parent[root] != root:
            root = parent[root]
        while parent[vertex] != root:
            parent[vertex], vertex = root, parent[vertex]
        return root

    for u, v, _ in flagged_pairs:
        root_u, root_v = find(u), find(v)
        if root_u != root_v:
            parent[root_u] = root_v
    groups: dict[int, list[int]] = defaultdict(list)
    for vertex in flagged_vertices:
        groups[find(vertex)].append(vertex)
    ratio_of_vertex: dict[int, float] = defaultdict(float)
    for u, v, ratio in flagged_pairs:
        ratio_of_vertex[u] = max(ratio_of_vertex[u], ratio)
        ratio_of_vertex[v] = max(ratio_of_vertex[v], ratio)

    regions: list[CurvatureJumpRegion] = []
    for members in groups.values():
        area = math.fsum(curvature.mixed_area[vertex] for vertex in members)
        centroid = (
            sum(mesh.vertices[vertex][0] for vertex in members) / len(members),
            sum(mesh.vertices[vertex][1] for vertex in members) / len(members),
            sum(mesh.vertices[vertex][2] for vertex in members) / len(members),
        )
        fraction = area / total_area if total_area > 0 else 0.0
        values = [curvature.curvature[vertex] for vertex in members]
        regions.append(
            CurvatureJumpRegion(
                vertex_count=len(members),
                area=area,
                area_fraction=fraction,
                centroid=centroid,
                max_ratio=max(ratio_of_vertex[vertex] for vertex in members),
                curvature_min=min(values),
                curvature_max=max(values),
                count_as_defect=fraction >= min_region_area_fraction,
            )
        )
    regions.sort(key=lambda region: region.area, reverse=True)
    flagged_area = math.fsum(curvature.mixed_area[vertex] for vertex in flagged_vertices)
    return FairnessAudit(
        threshold=threshold,
        vertex_count=len(mesh.vertices),
        interior_vertices=sum(1 for flag in curvature.interior if flag),
        edges_considered=considered,
        crease_edges_skipped=skipped,
        crease_edges=crease_edges,
        crease_vertices=len(excluded_vertices),
        mean_abs_curvature=curvature.mean_abs_curvature,
        max_jump_ratio=max_ratio,
        flagged_edges=len(flagged_pairs),
        flagged_vertices=len(flagged_vertices),
        flagged_area_fraction=flagged_area / total_area if total_area > 0 else 0.0,
        regions=tuple(regions),
    )# --------------------------------------------------------------------------- #
# offsets table (lines plan) integration
# --------------------------------------------------------------------------- #
def _simpson_or_trapezoid(values: Sequence[float], spacing: float) -> float:
    """Composite Simpson when the sample count is odd, trapezoid otherwise.

    Both are second order on a smooth integrand; Simpson is used where the grid
    allows it because the volume gate is tight (1 %) and the hull forms here have
    zero slope at the ends, where the trapezoid error is largest.
    """
    n = len(values)
    if n < 2:
        return 0.0
    if spacing <= 0.0:
        raise ValueError(f"spacing must be positive, got {spacing!r}")
    if n % 2 == 1:
        total = values[0] + values[-1]
        total += 4.0 * math.fsum(values[1:-1:2])
        total += 2.0 * math.fsum(values[2:-1:2])
        return total * spacing / 3.0
    total = 0.5 * (values[0] + values[-1]) + math.fsum(values[1:-1])
    return total * spacing


@dataclass(frozen=True)
class OffsetTable:
    """Half-breadths on a station x waterline grid: a lines plan as numbers.

    ``stations`` are longitudinal positions, ``levels`` heights above the keel and
    ``half_breadths[i][j]`` the half-breadth at station ``i``, level ``j``.  Units
    are metres.  ``source`` records whether the table came from analytic
    half-breadths or was digitised off the mesh, because that decides how
    independent the cross-check against the mesh actually is.
    """

    stations: tuple[float, ...]
    levels: tuple[float, ...]
    half_breadths: tuple[tuple[float, ...], ...]
    length: float
    beam: float
    draft: float
    source: str = "unspecified"

    def __post_init__(self) -> None:
        if len(self.stations) < 3 or len(self.levels) < 3:
            raise ValueError("an offsets table needs at least 3 stations and 3 levels")
        if len(self.half_breadths) != len(self.stations):
            raise ValueError("half_breadths must have one row per station")
        for row in self.half_breadths:
            if len(row) != len(self.levels):
                raise ValueError("every station row must have one entry per level")
        for name, values in (("stations", self.stations), ("levels", self.levels)):
            steps = [b - a for a, b in zip(values, values[1:])]
            if min(steps) <= 0.0:
                raise ValueError(f"{name} must be strictly increasing")

    @property
    def n_stations(self) -> int:
        return len(self.stations)

    @property
    def n_levels(self) -> int:
        return len(self.levels)

    @property
    def station_spacing(self) -> float:
        return self.stations[1] - self.stations[0]

    @property
    def level_spacing(self) -> float:
        return self.levels[1] - self.levels[0]

    def offset(self, i: int, j: int) -> float:
        return self.half_breadths[i][j]

    def to_dict(self) -> dict[str, Any]:
        return {
            "source": self.source,
            "stations": len(self.stations),
            "levels": len(self.levels),
            "length_m": round(self.length, 9),
            "beam_m": round(self.beam, 9),
            "draft_m": round(self.draft, 9),
            "max_half_breadth_m": round(max(max(row) for row in self.half_breadths), 9),
        }


def offsets_from_function(
    half_breadth: Callable[[float, float], float],
    *,
    length: float,
    beam: float,
    draft: float,
    n_stations: int = 201,
    n_levels: int = 51,
    source: str = "analytic half-breadths",
) -> OffsetTable:
    """Sample analytic half-breadths onto a station x waterline grid.

    ``half_breadth(x, z)`` takes ``x`` along the hull in ``[0, length]`` and ``z``
    above the keel in ``[0, draft]``.  The counts are forced odd so the composite
    Simpson rule applies in both directions.
    """
    if n_stations % 2 == 0:
        n_stations += 1
    if n_levels % 2 == 0:
        n_levels += 1
    stations = tuple(length * i / (n_stations - 1) for i in range(n_stations))
    levels = tuple(draft * j / (n_levels - 1) for j in range(n_levels))
    rows: list[tuple[float, ...]] = []
    for x in stations:
        row = []
        for z in levels:
            value = float(half_breadth(x, z))
            if not math.isfinite(value) or value < 0.0:
                value = 0.0
            row.append(value)
        rows.append(tuple(row))
    return OffsetTable(
        stations=stations,
        levels=levels,
        half_breadths=tuple(rows),
        length=length,
        beam=beam,
        draft=draft,
        source=source,
    )


def offsets_from_mesh(
    mesh: Mesh,
    *,
    waterline_z: float | None = None,
    n_stations: int = 81,
    n_levels: int = 41,
    source: str = "digitised off the mesh",
) -> OffsetTable:
    """Digitise an offsets table off a mesh: the "3-D model -> lines plan" step.

    Each station is cut with :func:`tensorcad.bonjean.section_loops`, and at every
    waterline the largest starboard half-breadth of that section is taken.  This
    only describes a hull whose half-breadth is single valued in ``z`` below the
    waterline; a bulbous bow or an overhanging flare needs a real loft, and the
    table is a *reconstruction* of the mesh rather than an independent source.
    """
    from .bonjean import section_loops

    box = mesh.bounding_box()
    top = float(box.maximum[2] if waterline_z is None else waterline_z)
    keel = float(box.minimum[2])
    draft = top - keel
    if draft <= 0.0:
        raise ValueError(f"waterline {top!r} is not above the keel {keel!r}")
    length = float(box.sizes[0])
    beam = float(box.sizes[1])
    if n_stations % 2 == 0:
        n_stations += 1
    if n_levels % 2 == 0:
        n_levels += 1
    # Sample the ends exactly: Simpson weights them, and a half-cell margin at the
    # stem is a first-order loss of wetted area (it showed up as a -1 % bias).
    stations = tuple(
        float(box.minimum[0]) + length * i / (n_stations - 1) for i in range(n_stations)
    )
    levels = tuple(draft * j / (n_levels - 1) for j in range(n_levels))
    rows: list[tuple[float, ...]] = []
    for station in stations:
        try:
            loops = section_loops(mesh, station)
        except Exception:  # pragma: no cover - a degenerate cut is reported, not fatal
            loops = ()
        polygon: list[tuple[float, float]] = [point for loop in loops for point in loop]
        row: list[float] = []
        for level in levels:
            height = level + keel
            best = 0.0
            for y, z in polygon:
                if y > best and abs(z - height) <= 1e-12:
                    best = y
            for (y1, z1), (y2, z2) in zip(polygon, polygon[1:] + polygon[:1]):
                low, high = (z1, z2) if z1 <= z2 else (z2, z1)
                if low <= height <= high and abs(z2 - z1) > 1e-15:
                    fraction = (height - z1) / (z2 - z1)
                    y = y1 + fraction * (y2 - y1)
                    best = max(best, y)
            row.append(max(0.0, best))
        rows.append(tuple(row))
    return OffsetTable(
        stations=stations,
        levels=levels,
        half_breadths=tuple(rows),
        length=length,
        beam=beam,
        draft=draft,
        source=source,
    )


@dataclass(frozen=True)
class TableHydrostatics:
    """Hydrostatics integrated from an offsets table, with the quadrature recorded."""

    displacement_volume: float
    wetted_surface_area: float
    waterplane_area: float
    block_coefficient: float
    prismatic_coefficient: float
    maximum_section_area: float
    section_areas: tuple[float, ...]
    source: str = "offsets table"

    def to_dict(self) -> dict[str, Any]:
        return {
            "displacement_volume_m3": round(self.displacement_volume, 9),
            "wetted_surface_area_m2": round(self.wetted_surface_area, 9),
            "waterplane_area_m2": round(self.waterplane_area, 9),
            "block_coefficient": round(self.block_coefficient, 6),
            "prismatic_coefficient": round(self.prismatic_coefficient, 6),
            "maximum_section_area_m2": round(self.maximum_section_area, 9),
            "source": self.source,
        }


def _slopes(table: OffsetTable) -> tuple[list[list[float]], list[list[float]]]:
    """One-sided (boundary) or central (interior) finite differences of the table."""
    dx = table.station_spacing
    dz = table.level_spacing
    ns, nl = table.n_stations, table.n_levels
    d_dx = [[0.0] * nl for _ in range(ns)]
    d_dz = [[0.0] * nl for _ in range(ns)]
    for i in range(ns):
        for j in range(nl):
            if i == 0:
                d_dx[i][j] = (table.offset(1, j) - table.offset(0, j)) / dx
            elif i == ns - 1:
                d_dx[i][j] = (table.offset(ns - 1, j) - table.offset(ns - 2, j)) / dx
            else:
                d_dx[i][j] = (table.offset(i + 1, j) - table.offset(i - 1, j)) / (2.0 * dx)
            if j == 0:
                d_dz[i][j] = (table.offset(i, 1) - table.offset(i, 0)) / dz
            elif j == nl - 1:
                d_dz[i][j] = (table.offset(i, nl - 1) - table.offset(i, nl - 2)) / dz
            else:
                d_dz[i][j] = (table.offset(i, j + 1) - table.offset(i, j - 1)) / (2.0 * dz)
    return d_dx, d_dz


def integrate_offsets(table: OffsetTable) -> TableHydrostatics:
    """Integrate volume, waterplane and wetted area from an offsets table.

    Volume: ``2 * int A_s(x) dx`` with the section areas ``A_s = 2 * int y dz``,
    Simpson in both directions.  Wetted area: ``2 * int int sqrt(1 + y_x^2 +
    y_z^2) dx dz``, evaluated per table cell from the corner slopes.  The factor 2
    is the mirror plane, so a *half* model in the table is the right input.
    """
    dz = table.level_spacing
    dx = table.station_spacing
    section_areas = [
        2.0 * _simpson_or_trapezoid(row, dz) for row in table.half_breadths
    ]
    volume = _simpson_or_trapezoid(section_areas, dx)
    waterplane = 2.0 * _simpson_or_trapezoid([row[-1] for row in table.half_breadths], dx)
    d_dx, d_dz = _slopes(table)
    wetted = 0.0
    ns, nl = table.n_stations, table.n_levels

    def integrand(i: int, j: int) -> float:
        sx = d_dx[i][j]
        sz = d_dz[i][j]
        return math.sqrt(1.0 + sx * sx + sz * sz)

    for i in range(ns - 1):
        for j in range(nl - 1):
            # Trapezoid over the cell: the four corner slopes, then the cell area.
            value = 0.25 * (
                integrand(i, j) + integrand(i + 1, j) + integrand(i, j + 1) + integrand(i + 1, j + 1)
            )
            wetted += value * dx * dz
    wetted *= 2.0
    maximum_section = max(section_areas) if section_areas else 0.0
    block = volume / (table.length * table.beam * table.draft)
    prismatic = volume / (table.length * maximum_section) if maximum_section > 0 else float("nan")
    return TableHydrostatics(
        displacement_volume=volume,
        wetted_surface_area=wetted,
        waterplane_area=waterplane,
        block_coefficient=block,
        prismatic_coefficient=prismatic,
        maximum_section_area=maximum_section,
        section_areas=tuple(section_areas),
        source=table.source,
    )


# --------------------------------------------------------------------------- #
# declared geometry gate (Cb / L / B / T first)
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class HullSpec:
    """What the caller says the model is."""

    length: float
    beam: float
    draft: float
    block_coefficient: float | None = None
    name: str = ""

    def __post_init__(self) -> None:
        for label, value in (("length", self.length), ("beam", self.beam), ("draft", self.draft)):
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{label} must be positive and finite, got {value!r}")
        if self.block_coefficient is not None and not 0.0 < self.block_coefficient <= 1.0:
            raise ValueError(
                f"block_coefficient must lie in (0, 1], got {self.block_coefficient!r}"
            )

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "length_m": self.length,
            "beam_m": self.beam,
            "draft_m": self.draft,
            "block_coefficient": self.block_coefficient,
        }


@dataclass(frozen=True)
class DeclaredGeometryAudit:
    """Realised dimensions and ``Cb`` against the declaration."""

    spec: HullSpec
    realised_length: float
    realised_beam: float
    realised_draft: float
    realised_block_coefficient: float
    dimension_tolerance: float
    cb_tolerance: float
    defects: tuple[GeometryDefect, ...] = ()

    @property
    def passed(self) -> bool:
        return not self.defects

    def relative(self, name: str) -> float:
        table = {
            "length": (self.realised_length, self.spec.length),
            "beam": (self.realised_beam, self.spec.beam),
            "draft": (self.realised_draft, self.spec.draft),
            "block_coefficient": (self.realised_block_coefficient, self.spec.block_coefficient),
        }
        actual, declared = table[name]
        if declared in (None, 0.0):
            return float("nan")
        return (actual - float(declared)) / float(declared)

    def to_dict(self) -> dict[str, Any]:
        return {
            "passed": self.passed,
            "spec": self.spec.to_dict(),
            "realised": {
                "length_m": round(self.realised_length, 9),
                "beam_m": round(self.realised_beam, 9),
                "draft_m": round(self.realised_draft, 9),
                "block_coefficient": round(self.realised_block_coefficient, 9),
            },
            "relative_deviation": {
                name: round(self.relative(name), 9)
                for name in ("length", "beam", "draft", "block_coefficient")
            },
            "dimension_tolerance": self.dimension_tolerance,
            "cb_tolerance": self.cb_tolerance,
            "defects": [defect.to_dict() for defect in self.defects],
        }


def audit_declared_geometry(
    mesh: Mesh,
    spec: HullSpec,
    *,
    waterline_z: float | None = None,
    volume: float | None = None,
    dimension_tolerance: float = DIMENSION_TOLERANCE,
    cb_tolerance: float = DECLARED_CB_TOLERANCE,
) -> DeclaredGeometryAudit:
    """Compare realised ``L/B/T/Cb`` with the declaration, before any hydrostatics."""
    box = mesh.bounding_box()
    top = float(box.maximum[2] if waterline_z is None else waterline_z)
    keel = float(box.minimum[2])
    realised = {
        "length": float(box.sizes[0]),
        "beam": float(box.sizes[1]),
        "draft": top - keel,
    }
    if volume is None:
        volume = abs(meshlib.total_signed_volume(mesh))
    displacement = abs(float(volume))
    cb = displacement / (realised["length"] * realised["beam"] * realised["draft"])
    declared = {
        "length": spec.length,
        "beam": spec.beam,
        "draft": spec.draft,
    }
    defects: list[GeometryDefect] = []
    for name, value in realised.items():
        deviation = (value - declared[name]) / declared[name]
        if abs(deviation) > dimension_tolerance:
            defects.append(
                GeometryDefect(
                    kind=f"{name}_mismatch",
                    magnitude=value,
                    unit="m",
                    detail=(
                        f"realised {name} {value:.6g} m against declared {declared[name]:.6g} m: "
                        f"{100.0 * deviation:+.2f}% (tolerance {100.0 * dimension_tolerance:.2f}%)"
                    ),
                    relative=deviation,
                )
            )
    if spec.block_coefficient is not None:
        deviation = (cb - spec.block_coefficient) / spec.block_coefficient
        if abs(deviation) > cb_tolerance:
            defects.append(
                GeometryDefect(
                    kind="block_coefficient_mismatch",
                    magnitude=cb,
                    unit="dimensionless",
                    detail=(
                        f"realised Cb {cb:.6f} against declared {spec.block_coefficient:.6f}: "
                        f"{100.0 * deviation:+.2f}% (tolerance {100.0 * cb_tolerance:.2f}%)"
                    ),
                    relative=deviation,
                )
            )
    return DeclaredGeometryAudit(
        spec=spec,
        realised_length=realised["length"],
        realised_beam=realised["beam"],
        realised_draft=realised["draft"],
        realised_block_coefficient=cb,
        dimension_tolerance=dimension_tolerance,
        cb_tolerance=cb_tolerance,
        defects=tuple(defects),
    )


# --------------------------------------------------------------------------- #
# hydrostatic consistency: table vs mesh
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class ConsistencyAudit:
    """Agreement between the offsets-table and mesh-integrated hydrostatics."""

    table: TableHydrostatics
    mesh: Hydrostatics
    tolerance: float
    differences: dict[str, float]

    @property
    def passed(self) -> bool:
        return all(abs(value) <= self.tolerance for value in self.differences.values())

    @property
    def defects(self) -> tuple[GeometryDefect, ...]:
        findings: list[GeometryDefect] = []
        for name, value in self.differences.items():
            if abs(value) <= self.tolerance:
                continue
            findings.append(
                GeometryDefect(
                    kind="hydrostatic_inconsistency",
                    magnitude=value,
                    unit="relative",
                    detail=(
                        f"{name}: offsets-table and mesh integrals differ by {100.0 * value:+.3f}% "
                        f"(tolerance {100.0 * self.tolerance:.2f}%); table uses {self.table.source}"
                    ),
                    relative=value,
                )
            )
        return tuple(findings)

    def to_dict(self) -> dict[str, Any]:
        return {
            "passed": self.passed,
            "tolerance": self.tolerance,
            "table_source": self.table.source,
            "compare": {
                name: {
                    "offsets_table": round(getattr(self.table, name), 9),
                    "mesh": round(getattr(self.mesh, name), 9),
                    "relative_difference": round(value, 9),
                }
                for name, value in self.differences.items()
            },
            "table": self.table.to_dict(),
            "mesh": self.mesh.to_dict(),
        }


_COMPARED_KEYS = (
    "displacement_volume",
    "wetted_surface_area",
    "block_coefficient",
    "waterplane_area",
    "prismatic_coefficient",
)


def audit_hydrostatic_consistency(
    table_hydrostatics: TableHydrostatics,
    mesh_hydrostatics: Hydrostatics,
    *,
    tolerance: float = HYDROSTATIC_TOLERANCE,
) -> ConsistencyAudit:
    """Compare two independent integrations of the same hull."""
    differences: dict[str, float] = {}
    for key in _COMPARED_KEYS:
        table_value = getattr(table_hydrostatics, key)
        mesh_value = getattr(mesh_hydrostatics, key)
        if not (math.isfinite(table_value) and math.isfinite(mesh_value)):
            continue
        reference = max(abs(mesh_value), abs(table_value))
        if reference <= 1e-12:
            continue
        differences[key] = (table_value - mesh_value) / mesh_value
    return ConsistencyAudit(
        table=table_hydrostatics,
        mesh=mesh_hydrostatics,
        tolerance=tolerance,
        differences=differences,
    )


# --------------------------------------------------------------------------- #
# the pipeline
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class GeometryAuditReport:
    """One audit pass over one hull model."""

    watertight: WatertightnessAudit
    fairness: FairnessAudit | None
    declared: DeclaredGeometryAudit | None
    consistency: ConsistencyAudit | None
    mesh_hydrostatics: Hydrostatics | None = None
    tables: tuple[OffsetTable, ...] = ()
    stages_run: tuple[str, ...] = ()
    schema: str = SCHEMA_AUDIT

    @property
    def accepted(self) -> bool:
        for stage in (self.watertight, self.fairness, self.declared, self.consistency):
            if stage is not None and not stage.passed:
                return False
        return True

    @property
    def defects(self) -> tuple[GeometryDefect, ...]:
        found: list[GeometryDefect] = list(self.watertight.defects)
        for stage in (self.fairness, self.declared, self.consistency):
            if stage is not None:
                found.extend(stage.defects)
        return tuple(found)

    @property
    def verdict(self) -> str:
        return "accepted" if self.accepted else "rejected"

    def to_dict(self) -> dict[str, Any]:
        data: dict[str, Any] = {
            "schema": self.schema,
            "verdict": self.verdict,
            "accepted": self.accepted,
            "stages_run": list(self.stages_run),
            "defect_count": len(self.defects),
            "defects": [defect.to_dict() for defect in self.defects],
            "watertightness": self.watertight.to_dict(),
        }
        if self.fairness is not None:
            data["fairness"] = self.fairness.to_dict()
        if self.declared is not None:
            data["declared_geometry"] = self.declared.to_dict()
        if self.consistency is not None:
            data["hydrostatic_consistency"] = self.consistency.to_dict()
        if self.mesh_hydrostatics is not None:
            data["mesh_hydrostatics"] = self.mesh_hydrostatics.to_dict()
        if self.tables:
            data["offsets_tables"] = [table.to_dict() for table in self.tables]
        return data

    def to_json(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(self.to_dict(), indent=2) + "\n", encoding="utf-8")
        return target


def audit_hull_geometry(
    mesh: Mesh,
    *,
    spec: HullSpec | None = None,
    offsets: OffsetTable | None = None,
    half_breadth: Callable[[float, float], float] | None = None,
    waterline_z: float | None = None,
    water_density: float = 1025.0,
    gap_tolerance_fraction: float = GAP_TOLERANCE_FRACTION,
    curvature_threshold: float = CURVATURE_JUMP_THRESHOLD,
    hydrostatic_tolerance: float = HYDROSTATIC_TOLERANCE,
    dimension_tolerance: float = DIMENSION_TOLERANCE,
    cb_tolerance: float = DECLARED_CB_TOLERANCE,
    check_fairness: bool = True,
) -> GeometryAuditReport:
    """Run the full audit: topology, fairness, declared geometry, consistency.

    The order matters and is enforced: a surface that does not bound a volume is
    not integrated at all, because the displacement of an open shell is not a
    number about the ship.  A rejected report therefore carries the topology
    defects that stopped it and ``stages_run`` says which stages were reached.

    ``offsets`` (or ``half_breadth``, sampled onto a table) is the lines-plan
    route.  Without it the table is digitised off the mesh itself
    (:func:`offsets_from_mesh`), which still cross-checks the quadrature but not
    the geometry; the report records which of the two it used.
    """
    if not mesh.faces:
        raise ValueError("cannot audit a mesh with no faces")
    box = mesh.bounding_box()
    reference_length = spec.length if spec is not None else float(box.sizes[0])
    stages: list[str] = []

    watertight = audit_watertightness(
        mesh, length=reference_length, gap_tolerance_fraction=gap_tolerance_fraction
    )
    stages.append("watertightness")
    tables: list[OffsetTable] = []
    if not watertight.passed:
        return GeometryAuditReport(
            watertight=watertight,
            fairness=None,
            declared=None,
            consistency=None,
            stages_run=tuple(stages),
        )

    fairness = None
    if check_fairness:
        # A fairness failure is recorded but does not stop the pipeline: a sharp
        # bilge does not invalidate the displacement it encloses, and suppressing
        # the hydrostatics would only hide the geometry finding.  Topology and the
        # declared-geometry gate *do* stop it, for the reasons given above.
        fairness = audit_fairness(mesh, threshold=curvature_threshold)
        stages.append("fairness")

    top = float(box.maximum[2] if waterline_z is None else waterline_z)
    mesh_hydro = compute_hydrostatics(mesh, top, water_density=water_density, closed=True)
    stages.append("mesh_hydrostatics")

    declared = None
    if spec is not None:
        declared = audit_declared_geometry(
            mesh,
            spec,
            waterline_z=top,
            volume=mesh_hydro.displacement_volume,
            dimension_tolerance=dimension_tolerance,
            cb_tolerance=cb_tolerance,
        )
        stages.append("declared_geometry")
        if not declared.passed:
            return GeometryAuditReport(
                watertight=watertight,
                fairness=fairness,
                declared=declared,
                consistency=None,
                mesh_hydrostatics=mesh_hydro,
                stages_run=tuple(stages),
            )

    if offsets is None and half_breadth is not None:
        offsets = offsets_from_function(
            half_breadth,
            length=float(box.sizes[0]),
            beam=float(box.sizes[1]),
            draft=top - float(box.minimum[2]),
        )
    if offsets is None:
        offsets = offsets_from_mesh(mesh, waterline_z=top)
    tables.append(offsets)
    table_hydro = integrate_offsets(offsets)
    consistency = audit_hydrostatic_consistency(
        table_hydro, mesh_hydro, tolerance=hydrostatic_tolerance
    )
    stages.append("hydrostatic_consistency")
    return GeometryAuditReport(
        watertight=watertight,
        fairness=fairness,
        declared=declared,
        consistency=consistency,
        mesh_hydrostatics=mesh_hydro,
        tables=tuple(tables),
        stages_run=tuple(stages),
    )


# --------------------------------------------------------------------------- #
# injectors: the only way to know the audit can see anything
# --------------------------------------------------------------------------- #
def inject_seam_offset(
    mesh: Mesh,
    *,
    delta: float,
    x_split: float,
    axis: int = 0,
    snap: bool = True,
) -> tuple[Mesh, tuple[str, ...]]:
    """Separate the two patches at ``x_split`` by ``delta`` along ``axis``.

    The seam vertices are *duplicated* for the aft patch before it is translated,
    so the two pieces become open shells whose free boundaries sit exactly ``delta``
    apart.  Simply translating shared vertices would keep the solid closed and the
    defect invisible -- which is precisely the case that makes this injector worth
    having: it is the only way to be sure the gap detector is not measuring its own
    bookkeeping.

    ``snap`` moves the cut onto the nearest vertex plane, which is what makes the
    offset an exact ground truth; ``snap=False`` reproduces a staircase tear.

    The default offset is longitudinal (``axis=0``) because it is the only
    direction in which the two free boundaries are *exactly* ``delta`` apart:
    a lateral offset near a narrow section can bring the port flank of the aft
    piece next to the starboard flank of the fore piece, and the pieces then
    genuinely touch, so the smallest separation is legitimately smaller than
    ``delta``.  Both cases are defects; only the longitudinal one is an exact
    ground truth for scoring.

    Returns the defected mesh and the per-face patch ids (``aft``/``fore``).
    """
    cut = float(x_split)
    if snap:
        # Move the cut onto the nearest plane of vertices.  A cut that falls between
        # two mesh rings tears the surface along a staircase, and then the smallest
        # separation is |h - delta| (one station spacing minus the offset) rather
        # than the offset: real, but not a usable ground truth.  With the cut on a
        # ring the two pieces are exactly ``delta`` apart.
        planes = sorted({vertex[axis] for vertex in mesh.vertices})
        cut = min(planes, key=lambda value: abs(value - x_split))
    aft = [
        index for index, face in enumerate(mesh.faces)
        if meshlib.face_centroid(mesh, face)[axis] < cut
    ]
    if not aft or len(aft) == len(mesh.faces):
        raise ValueError(f"x_split {x_split!r} does not divide the mesh into two patches")
    aft_set = set(aft)
    fore_vertices = {vertex for index, face in enumerate(mesh.faces)
                     if index not in aft_set for vertex in face}
    vertices = list(mesh.vertices)
    offset = [0.0, 0.0, 0.0]
    # The aft patch is moved *away* from the fore patch (negative along the axis), so
    # the defect is a clean gap rather than an interpenetration.
    offset[axis] = -float(delta)
    copy_of: dict[int, int] = {}
    faces: list[tuple[int, int, int]] = []
    for index, face in enumerate(mesh.faces):
        if index not in aft_set:
            faces.append(face)
            continue
        moved: list[int] = []
        for vertex in face:
            if vertex not in copy_of:
                point = vertices[vertex]
                copy_of[vertex] = len(vertices)
                vertices.append(
                    (point[0] + offset[0], point[1] + offset[1], point[2] + offset[2])
                )
                if vertex not in fore_vertices:
                    # An interior vertex of the aft patch has no seam to open; reusing
                    # its copy keeps the patch connected.
                    vertices[copy_of[vertex]] = point if delta == 0.0 else vertices[copy_of[vertex]]
            moved.append(copy_of[vertex])
        faces.append((moved[0], moved[1], moved[2]))
    shifted = Mesh(
        vertices=tuple(vertices), faces=tuple(faces), units=mesh.units,
        coordinate_frame=mesh.coordinate_frame,
    )
    patch_ids = tuple("aft" if index in aft_set else "fore" for index in range(len(mesh.faces)))
    return shifted, patch_ids


def inject_hole(
    mesh: Mesh,
    *,
    x_centre: float,
    z_centre: float,
    x_size: float,
    z_size: float,
) -> Mesh:
    """Delete the facets inside an ``x/z`` window, leaving one open boundary."""
    kept: list[tuple[int, int, int]] = []
    for face in mesh.faces:
        cx, _, cz = meshlib.face_centroid(mesh, face)
        inside = abs(cx - x_centre) <= 0.5 * x_size and abs(cz - z_centre) <= 0.5 * z_size
        if not inside:
            kept.append(face)
    if len(kept) == len(mesh.faces):
        raise ValueError("the requested window removed no facets")
    return Mesh(vertices=mesh.vertices, faces=tuple(kept), units=mesh.units,
                coordinate_frame=mesh.coordinate_frame)


def inject_flip_faces(mesh: Mesh, *, count: int, seed: int = 0) -> Mesh:
    """Reverse the winding of ``count`` random facets: inconsistent normals."""
    if count <= 0:
        raise ValueError("count must be positive")
    rng = random.Random(seed)
    chosen = set(rng.sample(range(len(mesh.faces)), min(count, len(mesh.faces))))
    faces = [
        (face[0], face[2], face[1]) if index in chosen else face
        for index, face in enumerate(mesh.faces)
    ]
    return Mesh(vertices=mesh.vertices, faces=tuple(faces), units=mesh.units,
                coordinate_frame=mesh.coordinate_frame)


def inject_duplicate_faces(mesh: Mesh, *, count: int, seed: int = 0) -> Mesh:
    """Append ``count`` copies of random facets: overlapping shells."""
    if count <= 0:
        raise ValueError("count must be positive")
    rng = random.Random(seed)
    extra = [mesh.faces[index] for index in rng.sample(range(len(mesh.faces)),
                                                       min(count, len(mesh.faces)))]
    return Mesh(vertices=mesh.vertices, faces=mesh.faces + tuple(extra), units=mesh.units,
                coordinate_frame=mesh.coordinate_frame)


def inject_cb_offset(
    mesh: Mesh,
    *,
    factor: float,
    waterline_z: float | None = None,
) -> Mesh:
    """Displace ``Cb`` by ``factor`` while leaving ``L``, ``B`` and ``T`` untouched.

    Scaling the whole beam would trip the beam check first and prove nothing about
    the ``Cb`` gate.  Instead the half-breadths are reweighted with a profile that
    vanishes at the waterline -- ``y' = y * (1 + k * w(z))`` with ``w`` quadratic in
    depth -- so the deck edge keeps the declared beam and the added volume goes
    below it.  Volume is linear in ``k``, so one probe evaluation solves for ``k``
    exactly; the realised block coefficient then moves by ``factor`` and by nothing
    else.
    """
    if factor <= 0.0:
        raise ValueError("factor must be positive")
    box = mesh.bounding_box()
    top = float(box.maximum[2] if waterline_z is None else waterline_z)
    draft = top - float(box.minimum[2])
    if draft <= 0.0:
        raise ValueError("waterline is not above the keel")
    base = abs(meshlib.total_signed_volume(mesh))

    def scaled(k: float) -> Mesh:
        vertices = []
        for x, y, z in mesh.vertices:
            depth_fraction = min(max((top - z) / draft, 0.0), 1.0)
            weight = depth_fraction * depth_fraction
            vertices.append((x, y * (1.0 + k * weight), z))
        return Mesh(vertices=tuple(vertices), faces=mesh.faces, units=mesh.units,
                    coordinate_frame=mesh.coordinate_frame)

    probe = abs(meshlib.total_signed_volume(scaled(1.0)))
    delta = probe - base
    if abs(delta) <= 1e-15:
        raise ValueError("the mesh has no volume to redistribute")
    k = (factor - 1.0) * base / delta
    return scaled(k)


def inject_bump(
    mesh: Mesh,
    *,
    centre: Sequence[float],
    amplitude: float,
    sigma: float,
    field: CurvatureField | None = None,
) -> Mesh:
    """Push vertices along their normals with a Gaussian profile: a dent or blister.

    ``amplitude``/``sigma`` control how sharp the curvature break is; the audit is
    scored on whether the resulting ``|dK|/K`` region is marked, and the test that
    drives it measures ``|dK|/K`` before and after so the ground truth is not
    assumed.
    """
    if sigma <= 0.0:
        raise ValueError("sigma must be positive")
    normals = (field or gaussian_curvature(mesh)).normals
    centre_point = tuple(float(value) for value in centre)
    vertices = []
    for index, point in enumerate(mesh.vertices):
        distance = math.dist(point, centre_point)
        weight = math.exp(-0.5 * (distance / sigma) ** 2)
        normal = normals[index]
        vertices.append(
            (
                point[0] + amplitude * weight * normal[0],
                point[1] + amplitude * weight * normal[1],
                point[2] + amplitude * weight * normal[2],
            )
        )
    return Mesh(vertices=tuple(vertices), faces=mesh.faces, units=mesh.units,
                coordinate_frame=mesh.coordinate_frame)


# --------------------------------------------------------------------------- #
# gates artifact
# --------------------------------------------------------------------------- #
def write_gates_json(
    path: str | Path,
    gates: Sequence[dict[str, Any]],
    *,
    pipeline: str = "tensorcad.geometry_audit",
    mission: str = "",
    extra: dict[str, Any] | None = None,
) -> Path:
    """Write the machine-readable gate table (one row per acceptance criterion)."""
    rows = []
    for gate in gates:
        row = dict(gate)
        row.setdefault("passed", False)
        rows.append(row)
    payload: dict[str, Any] = {
        "schema": SCHEMA_GATES,
        "pipeline": pipeline,
        "mission": mission,
        "generated_by": "tensorcad.geometry_audit.write_gates_json",
        "all_passed": all(bool(row["passed"]) for row in rows),
        "gate_count": len(rows),
        "gates": rows,
    }
    if extra:
        payload["extra"] = extra
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return target