"""B-rep topology for CAD surface patch sets: edge-face graph and watertightness.

The object a shipyard hands over is a set of trimmed surfaces with *sewn* topology;
the object a solver can use is a closed, orientable, consistently wound shell.  The
difference between the two is invisible in a viewer, so the pipeline in
:mod:`tensorcad.geomquality` decides it here, on the adjacency graph, before any
hydrostatic number is computed.

Model
-----
A :class:`BRep` is ``faces`` plus ``edges``.  Each face carries its supporting
surface and an ordered boundary loop of :class:`EdgeRef` (edge index + traversal
direction).  This is the same graph a real B-rep kernel keeps, only without the
trimmed-pcurve detail: an edge is used once by one face when the model is open,
twice with opposite directions when it is a manifold interior edge, and three or
more times when it is non-manifold.

Two ways to obtain the graph:

``BRep(faces=..., edges=...)``
    explicit topology, as read from STEP (``edge_curve`` / ``edge_loop``), where the
    file itself declares which faces share which edge.
:func:`sew`
    the graph is *derived from geometry*: each face contributes its four boundary
    isocurves and curves that coincide within a tolerance are welded into one edge.
    This is what an IGES export of untrimmed patches gives you, and it is where a
    seam crack (two boundaries that nearly but not exactly coincide) shows up: the
    two curves stay distinct edges, so the shell has free edges *and* a measurable
    gap.  Sewing is therefore an audit, not a repair.

Every finding is a :class:`Defect` carrying a code, a magnitude in metres, the same
magnitude as a fraction of ``L`` and a repair hint, because "a 2 mm gap" is a defect
on a 6 m model and is not one on a 300 m ship.
"""
from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from typing import Any

from ..geometry import Mesh, MeshBuilder
from ..nurbs import NURBSSurface
from ._eval import evaluate_fast

Vec3 = tuple[float, float, float]

#: Absolute weld tolerance in model units (metres): two boundary samples closer than
#: this are the same point as far as CAD sewing is concerned.
DEFAULT_TOLERANCE = 1e-6

#: A crack whose widest separation stays under this fraction of ``L`` is a modelling
#: artefact of the export, not a defect; above it the shell leaks.
GAP_TOLERANCE_FRACTION = 1e-5

#: How many samples per boundary curve are compared when sewing / measuring gaps.
BOUNDARY_SAMPLES = 33

#: An edge shorter than ``DEGENERATE_FRACTION * L`` is a collapsed boundary (keel line
#: end, bow apex).  Real B-rep kernels keep those as seam edges: they bound no area and
#: cannot leak.  Above this and below the weld tolerance, the edge is a defect.
DEGENERATE_FRACTION = 1e-12

#: Repair hints keyed by defect code.  Kept next to the detector that emits the code
#: so a report can always say what to do about it, not only that something is wrong.
REPAIR_HINTS: dict[str, str] = {
    "FREE_EDGE": "sew the two coincident boundaries, or cap the opening with a face",
    "OPEN_TOP": "close the waterline: add the deck / waterplane face at z = T",
    "SEAM_GAP": "snap the seam vertices to a common tolerance and re-sew the patch pair",
    "NON_MANIFOLD": "remove the extra face at the edge, or split the shell at the T-junction",
    "ORIENTATION": "flip the offending face so its loop is anticlockwise about the outward normal",
    "DUPLICATE_FACE": "delete the redundant face; it contributes displacement twice",
    "SLIVER": "collapse the degenerate face to an edge or rebuild the patch with a coarser net",
    "MICRO_EDGE": "merge the two vertices of the micro edge (tolerance weld)",
    "SELF_INTERSECTION": "reduce the offending control point; the patch normal reverses inside the domain",
    "TANGENT_JUMP": "make the two patch tangents collinear at the seam (G1) before fairing",
    "CURVATURE_JUMP": "re-fair the seam: match the second derivative (G2) or accept the knuckle explicitly",
    "ROUGHNESS": "low-pass the control net; a high-frequency ripple is not a design feature",
    "CB_MISMATCH": "re-scale the sections to the declared Cb, or correct the declared value",
    "UNIT_SCALE": "apply the declared unit factor to the whole model (mm/m confusion)",
    "NEGATIVE_VOLUME": "reverse the shell orientation: the signed volume is negative",
}

#: Defect codes considered fatal for the watertightness score.
TOPOLOGY_CODES = frozenset(
    {
        "FREE_EDGE",
        "OPEN_TOP",
        "SEAM_GAP",
        "NON_MANIFOLD",
        "ORIENTATION",
        "DUPLICATE_FACE",
        "SLIVER",
        "MICRO_EDGE",
        "NEGATIVE_VOLUME",
    }
)


@dataclass(frozen=True)
class Defect:
    """One machine-readable finding.

    ``magnitude`` is in model units (metres for a hull), ``relative`` is the same
    quantity divided by the reference length ``L`` when the audit knows it, and
    ``repair`` is copied from :data:`REPAIR_HINTS`.
    """

    code: str
    message: str
    magnitude: float = 0.0
    relative: float | None = None
    where: dict[str, Any] = field(default_factory=dict)
    repair: str = ""

    def __post_init__(self) -> None:
        if not self.repair:
            object.__setattr__(self, "repair", REPAIR_HINTS.get(self.code, ""))

    def to_dict(self) -> dict[str, Any]:
        return {
            "code": self.code,
            "message": self.message,
            "magnitude": self.magnitude,
            "relative": self.relative,
            "where": dict(self.where),
            "repair": self.repair,
        }


@dataclass(frozen=True)
class Edge:
    """A boundary curve: a polyline plus its arc length, in model units."""

    index: int
    points: tuple[Vec3, ...]

    @property
    def length(self) -> float:
        return sum(
            math.dist(self.points[i], self.points[i + 1])
            for i in range(len(self.points) - 1)
        )

    @property
    def start(self) -> Vec3:
        return self.points[0]

    @property
    def end(self) -> Vec3:
        return self.points[-1]

    def reversed(self) -> Edge:
        return Edge(index=self.index, points=tuple(reversed(self.points)))

    def to_dict(self) -> dict[str, Any]:
        return {
            "index": self.index,
            "points": len(self.points),
            "length": self.length,
            "start": self.start,
            "end": self.end,
        }


@dataclass(frozen=True)
class EdgeRef:
    """One use of an edge in a face boundary loop."""

    edge: int
    forward: bool = True


@dataclass(frozen=True)
class Face:
    """A supporting surface plus its boundary loop.

    ``same_sense`` is the STEP flag: ``False`` means the face normal is the *reverse*
    of ``S_u x S_v``.  ``full_domain`` says whether the loop covers the whole
    rectangular parameter domain -- only such faces can be integrated by patch
    quadrature; anything else falls back to tessellation.
    """

    index: int
    surface: NURBSSurface
    boundary: tuple[EdgeRef, ...]
    same_sense: bool = True
    name: str = ""
    #: ``False`` for a face that is not part of the wetted surface (e.g. the deck
    #: patch used to close the waterplane).
    wetted: bool = True

    @property
    def domain(self) -> tuple[tuple[float, float], tuple[float, float]]:
        return self.surface.domain

    def corners(self) -> tuple[Vec3, Vec3, Vec3, Vec3]:
        """``(u0,v0)``, ``(u1,v0)``, ``(u1,v1)``, ``(u0,v1)`` -- the loop order."""
        (u0, u1), (v0, v1) = self.domain
        evaluate = self.surface.evaluate
        return (evaluate(u0, v0), evaluate(u1, v0), evaluate(u1, v1), evaluate(u0, v1))

    def loop_points(self) -> tuple[Vec3, ...]:
        """Boundary loop as four polyline segments, in loop order."""
        (u0, u1), (v0, v1) = self.domain
        n = 8

        def evaluate(u: float, v: float) -> Vec3:
            return evaluate_fast(self.surface, u, v)

        segments = []
        for i in range(4):
            segment = []
            for k in range(n + 1):
                t = k / n
                if i == 0:
                    point = evaluate(u0 + (u1 - u0) * t, v0)
                elif i == 1:
                    point = evaluate(u1, v0 + (v1 - v0) * t)
                elif i == 2:
                    point = evaluate(u1 - (u1 - u0) * t, v1)
                else:
                    point = evaluate(u0, v1 - (v1 - v0) * t)
                segment.append(point)
            segments.append(tuple(segment))
        return tuple(segments)

    def to_dict(self) -> dict[str, Any]:
        return {
            "index": self.index,
            "name": self.name,
            "degree": self.surface.degree,
            "control_net": (len(self.surface.control_net), len(self.surface.control_net[0])),
            "rational": self.surface.rational,
            "same_sense": self.same_sense,
            "wetted": self.wetted,
            "edges": [ref.edge for ref in self.boundary],
        }


@dataclass
class BRep:
    """Faces plus the shared edge pool: a surface model with explicit topology."""

    faces: tuple[Face, ...]
    edges: tuple[Edge, ...]
    source: str = ""
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def edge_count(self) -> int:
        return len(self.edges)

    def edge_uses(self) -> dict[int, list[tuple[int, bool]]]:
        """``edge index -> [(face index, forward)]``: the adjacency graph itself."""
        uses: dict[int, list[tuple[int, bool]]] = {edge.index: [] for edge in self.edges}
        for face in self.faces:
            for ref in face.boundary:
                uses.setdefault(ref.edge, []).append((face.index, ref.forward))
        return uses

    def bounding_box(self) -> tuple[Vec3, Vec3]:
        points = [
            point
            for face in self.faces
            for segment in face.loop_points()
            for point in segment
        ]
        if not points:
            raise ValueError("BRep has no faces")
        lo = tuple(min(p[i] for p in points) for i in range(3))
        hi = tuple(max(p[i] for p in points) for i in range(3))
        return lo, hi  # type: ignore[return-value]

    @property
    def length(self) -> float:
        """Model length along x: the reference ``L`` every relative figure uses."""
        lo, hi = self.bounding_box()
        return hi[0] - lo[0]

    def tessellate(self, u_samples: int = 16, v_samples: int = 8) -> Mesh:
        """Sample every face into triangles, honouring ``same_sense``.

        Vertices are welded by exact position, so a seam that closes bit-for-bit
        closes in the mesh too.  Faces that only meet within a tolerance stay
        separate on purpose: the mesh is evidence, not a repair.
        """
        builder = MeshBuilder()
        for face in self.faces:
            (u0, u1), (v0, v1) = face.domain
            us = [u0 + (u1 - u0) * i / u_samples for i in range(u_samples + 1)]
            vs = [v0 + (v1 - v0) * j / v_samples for j in range(v_samples + 1)]
            us[0], us[-1] = u0, u1
            vs[0], vs[-1] = v0, v1
            grid = [[builder.add_vertex(evaluate_fast(face.surface, u, v)) for v in vs] for u in us]
            for i in range(u_samples):
                for j in range(v_samples):
                    quad = (grid[i][j], grid[i + 1][j], grid[i + 1][j + 1], grid[i][j + 1])
                    if face.same_sense:
                        builder.add_quad(*quad)
                    else:
                        builder.add_quad(quad[0], quad[3], quad[2], quad[1])
        return builder.build()

    def to_dict(self) -> dict[str, Any]:
        return {
            "source": self.source,
            "faces": [face.to_dict() for face in self.faces],
            "edges": len(self.edges),
            "metadata": dict(self.metadata),
        }


# --------------------------------------------------------------------------- #
# construction
# --------------------------------------------------------------------------- #
def _boundary_polyline(face_parts: Sequence[tuple[NURBSSurface, bool]], segment: int,
                       samples: int) -> tuple[Vec3, ...]:
    surface = face_parts[0][0]
    (u0, u1), (v0, v1) = surface.domain
    points = []
    for k in range(samples + 1):
        t = k / samples
        if segment == 0:
            point = surface.evaluate(u0 + (u1 - u0) * t, v0)
        elif segment == 1:
            point = surface.evaluate(u1, v0 + (v1 - v0) * t)
        elif segment == 2:
            point = surface.evaluate(u1 - (u1 - u0) * t, v1)
        else:
            point = surface.evaluate(u0, v1 - (v1 - v0) * t)
        points.append(point)
    return tuple(points)


def _curve_deviation(a: Sequence[Vec3], b: Sequence[Vec3]) -> float:
    """Largest pointwise separation between two equally sampled curves."""
    if len(a) != len(b):
        raise ValueError("curve deviation needs equal sample counts")
    return max(math.dist(p, q) for p, q in zip(a, b))


def _reverse_deviation(a: Sequence[Vec3], b: Sequence[Vec3]) -> float:
    return max(math.dist(p, q) for p, q in zip(a, reversed(b)))


def sew(
    surfaces: Iterable[NURBSSurface],
    *,
    tolerance: float = DEFAULT_TOLERANCE,
    samples: int = BOUNDARY_SAMPLES,
    names: Sequence[str] | None = None,
    same_sense: Sequence[bool] | None = None,
    wetted: Sequence[bool] | None = None,
    source: str = "sew",
    metadata: dict[str, Any] | None = None,
) -> BRep:
    """Build a B-rep from untrimmed patch surfaces by welding coincident boundaries.

    Each surface contributes four boundary isocurves.  A curve that matches an
    already collected edge -- forward or reversed, within ``tolerance`` -- becomes a
    *use* of that edge; anything else starts a new edge and therefore, if it is never
    matched, ends up as a free edge in :func:`audit_topology`.
    """
    surface_list = list(surfaces)
    sample_count = samples
    edges: list[Edge] = []
    edge_points: list[tuple[Vec3, ...]] = []
    faces: list[Face] = []
    flags = list(same_sense) if same_sense is not None else [True] * len(surface_list)
    if len(flags) != len(surface_list):
        raise ValueError("same_sense must have one entry per surface")
    wet = list(wetted) if wetted is not None else [True] * len(surface_list)
    if len(wet) != len(surface_list):
        raise ValueError("wetted must have one entry per surface")

    for face_index, surface in enumerate(surface_list):
        refs: list[EdgeRef] = []
        for segment in range(4):
            curve = _boundary_polyline([(surface, True)], segment, sample_count)
            match: tuple[int, bool] | None = None
            for edge_index, existing in enumerate(edge_points):
                if abs(math.dist(existing[0], curve[0]) ) <= tolerance and _curve_deviation(
                    existing, curve
                ) <= tolerance:
                    match = (edge_index, True)
                    break
                if abs(math.dist(existing[0], curve[-1])) <= tolerance and _reverse_deviation(
                    existing, curve
                ) <= tolerance:
                    match = (edge_index, False)
                    break
            if match is None:
                edge_index = len(edge_points)
                edge_points.append(curve)
                edges.append(Edge(index=edge_index, points=curve))
            else:
                edge_index = match[0]
            refs.append(EdgeRef(edge=edge_index, forward=match[1] if match else True))
        if not flags[face_index]:
            # The face normal is the reverse of S_u x S_v, so its outward loop is the
            # clockwise one: reverse the traversal and every edge direction.
            refs = [EdgeRef(edge=ref.edge, forward=not ref.forward) for ref in reversed(refs)]
        faces.append(
            Face(
                index=face_index,
                surface=surface,
                boundary=tuple(refs),
                same_sense=flags[face_index],
                name=(names[face_index] if names is not None else f"face{face_index}"),
                wetted=wet[face_index],
            )
        )
    return BRep(
        faces=tuple(faces),
        edges=tuple(edges),
        source=source,
        metadata=dict(metadata or {}),
    )


# --------------------------------------------------------------------------- #
# topology audit
# --------------------------------------------------------------------------- #
def _dot(a: Sequence[float], b: Sequence[float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _loop_normal(brep: BRep, face: Face) -> Vec3:
    """Newell normal of the *stored* boundary loop, walked through the edge pool."""
    edge_by_index = {edge.index: edge for edge in brep.edges}
    points: list[Vec3] = []
    for reference in face.boundary:
        edge = edge_by_index.get(reference.edge)
        if edge is None:
            continue
        polyline = edge.points if reference.forward else tuple(reversed(edge.points))
        points.extend(polyline if not points else polyline[1:])
    nx = ny = nz = 0.0
    for i, current in enumerate(points):
        nxt = points[(i + 1) % len(points)]
        nx += (current[1] - nxt[1]) * (current[2] + nxt[2])
        ny += (current[2] - nxt[2]) * (current[0] + nxt[0])
        nz += (current[0] - nxt[0]) * (current[1] + nxt[1])
    length = math.sqrt(nx * nx + ny * ny + nz * nz)
    if length == 0.0:
        return (0.0, 0.0, 0.0)
    return (nx / length, ny / length, nz / length)


def orient_shell(brep: BRep, *, flip_negative_volume: bool = True) -> BRep:
    """Propagate a consistent orientation across the adjacency graph.

    IGES carries no orientation at all, so a shell sewn from IGES surfaces arrives with
    arbitrary ``same_sense`` flags.  The propagation is the standard one: walk the face
    graph and require every interior edge to be traversed in *opposite* directions by its
    two faces; then, if the closed shell integrates to a negative volume, flip all of it.
    """
    uses = brep.edge_uses()
    adjacency: dict[int, list[tuple[int, int, int]]] = {face.index: [] for face in brep.faces}
    for edge_index, entries in uses.items():
        if len(entries) != 2:
            continue
        (face_a, forward_a), (face_b, forward_b) = entries
        adjacency[face_a].append((face_b, int(forward_a), int(forward_b)))
        adjacency[face_b].append((face_a, int(forward_b), int(forward_a)))
    orientation: dict[int, int] = {}
    for face in brep.faces:
        if face.index in orientation:
            continue
        orientation[face.index] = 1
        queue = [face.index]
        while queue:
            current = queue.pop()
            for neighbour, forward_current, forward_neighbour in adjacency[current]:
                if neighbour in orientation:
                    continue
                # forward_current * orient_current = -forward_neighbour * orient_neighbour
                direction_current = (1 if forward_current == 1 else -1) * orientation[current]
                orientation[neighbour] = -direction_current * (1 if forward_neighbour == 1 else -1)
                queue.append(neighbour)
    faces = []
    for face in brep.faces:
        flip = orientation.get(face.index, 1) < 0
        if flip:
            faces.append(
                Face(
                    index=face.index,
                    surface=face.surface,
                    boundary=tuple(
                        EdgeRef(edge=reference.edge, forward=not reference.forward)
                        for reference in reversed(face.boundary)
                    ),
                    same_sense=not face.same_sense,
                    name=face.name,
                    wetted=face.wetted,
                )
            )
        else:
            faces.append(face)
    oriented = BRep(faces=tuple(faces), edges=brep.edges, source=brep.source + "+oriented",
                    metadata=dict(brep.metadata))
    if flip_negative_volume:
        from .. import mesh as meshlib

        if meshlib.total_signed_volume(oriented.tessellate(8, 4)) < 0:
            flipped = [
                Face(
                    index=face.index,
                    surface=face.surface,
                    boundary=tuple(
                        EdgeRef(edge=reference.edge, forward=not reference.forward)
                        for reference in reversed(face.boundary)
                    ),
                    same_sense=not face.same_sense,
                    name=face.name,
                    wetted=face.wetted,
                )
                for face in oriented.faces
            ]
            oriented = BRep(
                faces=tuple(flipped),
                edges=oriented.edges,
                source=oriented.source + "+outward",
                metadata=dict(oriented.metadata),
            )
    return oriented


def _defect(
    code: str, message: str, magnitude: float, length: float | None, **where: Any
) -> Defect:
    relative = (magnitude / length) if (length and length > 0) else None
    return Defect(
        code=code,
        message=message,
        magnitude=magnitude,
        relative=relative,
        where=where,
    )


@dataclass
class TopologyAudit:
    """Result of the edge-face graph check: the watertightness score is 0 or 1."""

    defects: tuple[Defect, ...]
    free_edges: tuple[int, ...]
    non_manifold_edges: tuple[int, ...]
    coincident_seams: tuple[tuple[int, int, float], ...]
    signed_volume: float
    length: float
    edge_uses: dict[int, int]
    degenerate_edges: tuple[int, ...] = ()

    @property
    def watertight(self) -> bool:
        """1 iff the shell is closed, orientable and positively oriented."""
        return not any(defect.code in TOPOLOGY_CODES for defect in self.defects)

    @property
    def score(self) -> float:
        return 1.0 if self.watertight else 0.0

    @property
    def max_gap(self) -> float:
        gaps = [defect.magnitude for defect in self.defects if defect.code == "SEAM_GAP"]
        return max(gaps) if gaps else 0.0

    def to_dict(self) -> dict[str, Any]:
        return {
            "watertight": self.watertight,
            "score": self.score,
            "signed_volume": self.signed_volume,
            "length": self.length,
            "free_edges": list(self.free_edges),
            "degenerate_edges": list(self.degenerate_edges),
            "non_manifold_edges": list(self.non_manifold_edges),
            "coincident_seams": [
                {"edge_a": a, "edge_b": b, "deviation": d}
                for a, b, d in self.coincident_seams
            ],
            "edge_use_histogram": _histogram(self.edge_uses.values()),
            "defects": [defect.to_dict() for defect in self.defects],
        }


def _histogram(values: Iterable[int]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for value in values:
        key = str(value)
        counts[key] = counts.get(key, 0) + 1
    return dict(sorted(counts.items()))


def _free_edge_chains(
    free: set[int], edges: dict[int, Edge], tolerance: float
) -> list[list[int]]:
    """Group free edges into connected chains by shared endpoints."""
    remaining = set(free)
    chains: list[list[int]] = []
    while remaining:
        seed = remaining.pop()
        chain = [seed]
        grew = True
        while grew:
            grew = False
            for index in sorted(remaining):
                for anchor in (chain[0], chain[-1]):
                    edge = edges[anchor]
                    other = edges[index]
                    if (
                        min(
                            math.dist(edge.start, other.start),
                            math.dist(edge.start, other.end),
                            math.dist(edge.end, other.start),
                            math.dist(edge.end, other.end),
                        )
                        <= max(tolerance, 1e-9)
                    ):
                        chain.append(index)
                        remaining.discard(index)
                        grew = True
                        break
                if grew:
                    break
        chains.append(chain)
    return chains


def audit_topology(
    brep: BRep,
    *,
    tolerance: float = DEFAULT_TOLERANCE,
    gap_tolerance: float | None = None,
    length: float | None = None,
    expected_length: float | None = None,
    declared_cb: float | None = None,
    reference_cb: float | None = None,
    tessellation: tuple[int, int] = (8, 4),
) -> TopologyAudit:
    """Check the adjacency graph and return every finding plus a 0/1 score.

    Parameters that describe the *declared* model (``expected_length``,
    ``declared_cb``) turn the audit into a consistency check too: a model that
    integrates to a different block coefficient than the one on its drawing is a
    geometry defect even when the shell is closed.
    """
    uses = brep.edge_uses()
    edge_by_index = {edge.index: edge for edge in brep.edges}
    reference_length = length if length is not None else brep.length
    gap_limit = (
        gap_tolerance
        if gap_tolerance is not None
        else GAP_TOLERANCE_FRACTION * reference_length
    )
    defects: list[Defect] = []

    degenerate_limit = max(DEGENERATE_FRACTION * reference_length, 1e-15)
    degenerate = tuple(
        sorted(
            index for index, edge in edge_by_index.items() if edge.length <= degenerate_limit
        )
    )
    free_edges = tuple(
        sorted(
            index
            for index, entries in uses.items()
            if len(entries) == 1 and edge_by_index[index].length > degenerate_limit
        )
    )
    non_manifold = tuple(sorted(index for index, entries in uses.items() if len(entries) > 2))

    # 1. free edges, classified into holes (closed loops) and open boundaries.
    if free_edges:
        chains = _free_edge_chains(set(free_edges), edge_by_index, tolerance)
        total_length = 0.0
        for chain in chains:
            points: list[Vec3] = []
            for index in chain:
                edge = edge_by_index[index]
                if points and math.dist(points[-1], edge.start) > math.dist(
                    points[-1], edge.end
                ):
                    edge = edge.reversed()
                points.extend(edge.points if not points else edge.points[1:])
            chain_length = sum(
                math.dist(points[i], points[i + 1]) for i in range(len(points) - 1)
            )
            total_length += chain_length
            closure = math.dist(points[0], points[-1])
            top = max(point[2] for point in points)
            bottom = min(point[2] for point in points)
            if top - bottom <= degenerate_limit and len(points) > 2:
                # An open boundary lying in a plane of constant z is the waterline: the
                # deck is missing, which is a different repair from a hole in the shell.
                defects.append(
                    _defect(
                        "OPEN_TOP",
                        "open boundary at constant z: the waterline/deck is not closed",
                        chain_length,
                        reference_length,
                        chain=chain,
                        z=top,
                    )
                )
            elif closure <= gap_limit:
                defects.append(
                    _defect(
                        "FREE_EDGE",
                        f"closed boundary loop of {len(chain)} free edges: a hole",
                        chain_length,
                        reference_length,
                        chain=chain,
                        closure=closure,
                    )
                )
            elif top - bottom <= gap_limit:
                defects.append(
                    _defect(
                        "OPEN_TOP",
                        "open boundary at constant z: the waterline/deck is not closed",
                        chain_length,
                        reference_length,
                        chain=chain,
                        z=top,
                    )
                )
            else:
                defects.append(
                    _defect(
                        "FREE_EDGE",
                        f"open boundary chain of {len(chain)} edges, closure {closure:.3g}",
                        chain_length,
                        reference_length,
                        chain=chain,
                        closure=closure,
                    )
                )

    # 2. non-manifold edges.
    for index in non_manifold:
        edge = edge_by_index[index]
        defects.append(
            _defect(
                "NON_MANIFOLD",
                f"edge used by {len(uses[index])} faces: the shell is not a 2-manifold there",
                edge.length,
                reference_length,
                edge=index,
                faces=[face for face, _ in uses[index]],
            )
        )

    # 3. orientation: an interior edge must be traversed once in each direction.
    conflicts = [
        index
        for index, entries in uses.items()
        if len(entries) == 2 and entries[0][1] == entries[1][1]
    ]
    for index in sorted(conflicts):
        defects.append(
            _defect(
                "ORIENTATION",
                "both faces traverse the shared edge in the same direction",
                edge_by_index[index].length,
                reference_length,
                edge=index,
                faces=[face for face, _ in uses[index]],
            )
        )

    # 3b. each face's declared orientation must agree with the winding of its own loop.
    for face in brep.faces:
        try:
            normal = face.surface.normal(
                0.5 * (face.domain[0][0] + face.domain[0][1]),
                0.5 * (face.domain[1][0] + face.domain[1][1]),
            )
        except Exception:
            continue
        effective = normal if face.same_sense else tuple(-value for value in normal)
        loop_normal = _loop_normal(brep, face)
        if _dot(effective, loop_normal) < 0.0:
            defects.append(
                _defect(
                    "ORIENTATION",
                    f"face '{face.name}' is wound against its declared orientation",
                    abs(_dot(effective, loop_normal)),
                    reference_length,
                    face=face.index,
                )
            )

    # 4. coincident but unsewn seams: two distinct boundary curves that nearly meet.
    coincident: list[tuple[int, int, float]] = []
    free_list = list(free_edges)
    for i, a_index in enumerate(free_list):
        a = edge_by_index[a_index]
        for b_index in free_list[i + 1 :]:
            b = edge_by_index[b_index]
            if min(
                math.dist(a.start, b.start), math.dist(a.start, b.end)
            ) > gap_limit * 1e3 + 1e-9:
                continue
            forward = _reverse_deviation(a.points, b.points) if len(a.points) == len(b.points) else math.inf
            backward = _curve_deviation(a.points, b.points) if len(a.points) == len(b.points) else math.inf
            deviation = min(forward, backward)
            if deviation <= gap_limit:
                coincident.append((a_index, b_index, deviation))
    for a_index, b_index, deviation in coincident:
        defects.append(
            _defect(
                "SEAM_GAP",
                f"unsewn seam between edges {a_index} and {b_index}: "
                f"max separation {deviation:.3g} m",
                deviation,
                reference_length,
                edge_a=a_index,
                edge_b=b_index,
            )
        )

    # 5. micro edges and slivers.
    for edge in brep.edges:
        if degenerate_limit < edge.length <= tolerance:
            defects.append(
                _defect(
                    "MICRO_EDGE",
                    "boundary curve shorter than the weld tolerance",
                    edge.length,
                    reference_length,
                    edge=edge.index,
                )
            )
    areas = []
    for face in brep.faces:
        area = _face_area_estimate(face)
        areas.append(area)
        sliver_limit = (1e-4 * reference_length) ** 2
        if area <= sliver_limit:
            defects.append(
                _defect(
                    "SLIVER",
                    f"face area {area:.3g} m^2 is below the sliver limit",
                    area,
                    reference_length,
                    face=face.index,
                )
            )

    # 6. duplicate faces: identical supporting surface *and* identical domain.
    seen: dict[tuple[Any, ...], int] = {}
    for face in brep.faces:
        key = (
            face.surface.degree,
            tuple(tuple(point) for row in face.surface.control_net for point in row),
            face.surface.knots,
            face.domain,
        )
        digest = hash(key)
        if digest in seen:
            defects.append(
                _defect(
                    "DUPLICATE_FACE",
                    f"face {face.index} repeats face {seen[digest]}: displacement counted twice",
                    face.surface.area(order=6),
                    reference_length,
                    face=face.index,
                    duplicate_of=seen[digest],
                )
            )
        else:
            seen[digest] = face.index

    # 7. signed volume from the tessellation of every face (origin-tet divergence).
    from .. import mesh as meshlib

    volume = meshlib.total_signed_volume(brep.tessellate(*tessellation))
    if volume < 0:
        defects.append(
            _defect(
                "NEGATIVE_VOLUME",
                f"signed volume {volume:.6g} m^3 is negative: the shell is inside out",
                abs(volume),
                reference_length,
            )
        )

    # 8. declared-vs-realised checks (only when the caller declares something).
    if expected_length is not None and expected_length > 0:
        ratio = reference_length / expected_length
        if not 0.98 <= ratio <= 1.02:
            defects.append(
                _defect(
                    "UNIT_SCALE",
                    f"model length {reference_length:.6g} vs declared {expected_length:.6g} "
                    f"(factor {ratio:.6g})",
                    abs(reference_length - expected_length),
                    expected_length,
                )
            )
    if declared_cb is not None and reference_cb is not None and declared_cb > 0:
        deviation = abs(reference_cb - declared_cb) / declared_cb
        if deviation > 0.005:
            defects.append(
                Defect(
                    code="CB_MISMATCH",
                    message=(
                        f"realised Cb {reference_cb:.5f} vs declared {declared_cb:.5f} "
                        f"({100 * deviation:.2f} %)"
                    ),
                    magnitude=abs(reference_cb - declared_cb),
                    relative=deviation,
                )
            )

    return TopologyAudit(
        defects=tuple(defects),
        free_edges=free_edges,
        non_manifold_edges=non_manifold,
        degenerate_edges=degenerate,
        coincident_seams=tuple(coincident),
        signed_volume=volume,
        length=reference_length,
        edge_uses={index: len(entries) for index, entries in uses.items()},
    )


def _face_area_estimate(face: Face, u_samples: int = 8, v_samples: int = 4) -> float:
    """Triangle-sum area from a coarse sample grid: a sliver test, not a metrology."""
    (u0, u1), (v0, v1) = face.domain
    grid = [
        [
            evaluate_fast(face.surface, u0 + (u1 - u0) * i / u_samples, v0 + (v1 - v0) * j / v_samples)
            for j in range(v_samples + 1)
        ]
        for i in range(u_samples + 1)
    ]
    total = 0.0
    for i in range(u_samples):
        for j in range(v_samples):
            a, b, c, d = grid[i][j], grid[i + 1][j], grid[i + 1][j + 1], grid[i][j + 1]
            for p, q, r in ((a, b, c), (a, c, d)):
                u = (q[0] - p[0], q[1] - p[1], q[2] - p[2])
                v = (r[0] - p[0], r[1] - p[1], r[2] - p[2])
                cross = (
                    u[1] * v[2] - u[2] * v[1],
                    u[2] * v[0] - u[0] * v[2],
                    u[0] * v[1] - u[1] * v[0],
                )
                total += 0.5 * math.sqrt(cross[0] ** 2 + cross[1] ** 2 + cross[2] ** 2)
    return total


def merge_breps(*breps: BRep, source: str = "merged") -> BRep:
    """Concatenate B-reps, keeping edges distinct unless their polylines coincide."""
    faces: list[Face] = []
    edges: list[Edge] = []
    points: list[tuple[Vec3, ...]] = []
    tolerance = DEFAULT_TOLERANCE
    for brep in breps:
        remap: dict[int, int] = {}
        for edge in brep.edges:
            target = None
            for index, existing in enumerate(points):
                if _curve_deviation(existing, edge.points) <= tolerance or _reverse_deviation(
                    existing, edge.points
                ) <= tolerance:
                    target = index
                    break
            if target is None:
                target = len(points)
                points.append(edge.points)
                edges.append(Edge(index=target, points=edge.points))
            remap[edge.index] = target
        for face in brep.faces:
            faces.append(
                Face(
                    index=len(faces),
                    surface=face.surface,
                    boundary=tuple(
                        EdgeRef(edge=remap[ref.edge], forward=ref.forward)
                        for ref in face.boundary
                    ),
                    same_sense=face.same_sense,
                    name=face.name,
                    wetted=face.wetted,
                )
            )
    return BRep(faces=tuple(faces), edges=tuple(edges), source=source)


def mirror(brep: BRep, *, axis: int = 1, plane: float = 0.0) -> BRep:
    """Mirror a half model about ``x_axis = plane`` and weld the symmetry seam.

    A reflection reverses orientation, so the copied faces carry ``same_sense=False``
    *and* a reversed boundary loop: the two faces then traverse their shared seam edge
    in opposite directions, which is what the manifold test requires.
    """
    from ..nurbs import NURBSSurface

    def mirrored_surface(surface: NURBSSurface) -> NURBSSurface:
        def flip(point: Sequence[float]) -> tuple[float, float, float]:
            coordinates = list(point)
            coordinates[axis] = 2.0 * plane - coordinates[axis]
            return (coordinates[0], coordinates[1], coordinates[2])

        net = tuple(tuple(flip(point) for point in row) for row in surface.control_net)
        return NURBSSurface(
            degree=surface.degree,
            control_net=net,
            weights=surface.weights,
            knots=surface.knots,
        )

    faces = [
        Face(
            index=len(brep.faces) + index,
            surface=mirrored_surface(face.surface),
            boundary=tuple(
                EdgeRef(edge=ref.edge, forward=not ref.forward)
                for ref in reversed(face.boundary)
            ),
            same_sense=not face.same_sense,
            name=face.name + "_mirror",
            wetted=face.wetted,
        )
        for index, face in enumerate(brep.faces)
    ]
    def flip_edge(edge: Edge) -> Edge:
        points = []
        for point in edge.points:
            coordinates = list(point)
            coordinates[axis] = 2.0 * plane - coordinates[axis]
            points.append((coordinates[0], coordinates[1], coordinates[2]))
        return Edge(index=edge.index, points=tuple(points))

    extra = BRep(
        faces=tuple(faces),
        edges=tuple(flip_edge(edge) for edge in brep.edges),
        source=brep.source + "+mirror",
    )
    return merge_breps(brep, extra, source=brep.source + "+mirror")
