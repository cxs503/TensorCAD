"""Removal of local self-intersections: excavate the fold, then cap the hole.

Why this module exists
----------------------
:mod:`tensorcad.mesh_repair` closes, manifoldises, re-winds and fills holes, and it measures
self intersections so that a fold cannot slip through the acceptance test. What it could not
do is *remove* one: ``repair_mesh`` documents the defect as the thing it refuses, and a mesh
that arrives with an intersection is handed back unchanged with ``accepted=False``. The
parameterised hull, a filleted bracket, an exported CSG result and every imported STL can all
carry one, so "detected and refused" was the end of the route (``GAP-CAD-009``).

What it does
------------
The repair is deliberately *local*. An intersection is only ever reported between two facets
that do not share a vertex, so a fold such as a vertex poked through the far side of a shell
appears as a small set of penetrating facet pairs. The route is:

1. ``weld`` the input (a fold is often only a duplicate vertex apart) and measure it;
2. take the facets named by :func:`tensorcad.intersect.find_self_intersections` as seeds and
   dilate them by ``rings`` layers of edge adjacency -- the *excavation set*;
3. delete the excavation set and re-measure: any residual penetrating pair means the band was
   wider than the excavation, so grow one more ring and try again, up to ``max_rings``;
4. separately per original shell, keep the largest surviving piece and drop the shards that
   the cut detached (recorded, never silent);
5. close the resulting boundary loops with the planar fill of
   :func:`tensorcad.mesh_repair.fill_holes`, then ``orient_outward``;
6. accept only if the result is closed, edge-manifold, consistently wound, carries **no**
   self intersection, has positive volume, and -- when a ``reference_volume`` is supplied --
   moves the enclosed volume by no more than ``volume_tolerance``.

The answer is a :class:`DetangleResult`. On refusal it carries the **original** mesh, exactly
like ``repair_mesh``: a caller wants to know whether the mesh is usable, and
``healed=False`` with a printed reason is more useful than an exception.

Limits (say them out loud)
--------------------------
* This is not a boolean re-triangulation. Two sheets that genuinely cross over a wide area --
  not a local fold -- will either be capped into a membrane or refused with the residual pair
  count; the mesh is closed either way, but the shape at the crossing is not reconstructed.
* The cap is planar, so a hole spanning a strongly curved region is closed by a flat patch.
  That is the same trade ``fill_holes`` already makes, and it is why the volume tolerance is
  part of the acceptance test rather than an afterthought.
* A hole whose boundary loop is longer than ``max_hole_edges`` is refused, not patched.
* The input must be closed. Excavating a mesh that already has boundary edges would make the
  repair's own holes indistinguishable from the mesh's, so that case is refused up front.
"""
from __future__ import annotations

import math
from collections import defaultdict, deque
from dataclasses import dataclass, field
from typing import Any, Mapping, Sequence

from . import mesh as meshlib
from .geometry import Face, Mesh
from .intersect import find_self_intersections
from .mesh_repair import (
    DEFAULT_MAX_HOLE_EDGES,
    DEFAULT_VOLUME_TOLERANCE,
    DEFAULT_WELD_TOLERANCE,
    HoleFill,
    MeshMetrics,
    MeshRepairError,
    RepairStep,
    boundary_loops,
    fill_holes,
    measure_mesh,
)

__all__ = [
    "DEFAULT_MAX_RINGS",
    "DEFAULT_SELF_INTERSECTION_LIMIT",
    "DEFAULT_MIN_VOLUME_TOLERANCE",
    "DetangleResult",
    "ExcavationBand",
    "heal_self_intersections",
    "excavation_faces",
]

#: How many adjacency rings may be cut away before the repair gives up on a fold.
DEFAULT_MAX_RINGS = 4
#: Detection cap. A fold wider than this many pairs is not a *local* fold any more.
DEFAULT_SELF_INTERSECTION_LIMIT = 64
#: Floor on the volume check when the caller does not supply a reference volume: below this
#: the excavation is not local any more, whatever the reference.
DEFAULT_MIN_VOLUME_TOLERANCE = 0.05


# ---------------------------------------------------------------------------------------
# topology helpers (face-level, on the original vertex numbering)
# ---------------------------------------------------------------------------------------


def _face_neighbours(mesh: Mesh) -> dict[int, set[int]]:
    """Face adjacency across shared (undirected) edges."""
    owners: dict[tuple[int, int], list[int]] = defaultdict(list)
    for index, face in enumerate(mesh.faces):
        for edge in meshlib.face_edges(face):
            owners[meshlib.undirected_edge(edge)].append(index)
    neighbours: dict[int, set[int]] = {index: set() for index in range(len(mesh.faces))}
    for faces in owners.values():
        for i in range(len(faces)):
            for j in range(i + 1, len(faces)):
                if faces[i] != faces[j]:
                    neighbours[faces[i]].add(faces[j])
                    neighbours[faces[j]].add(faces[i])
    return neighbours


def _dilate(seeds: Sequence[int], neighbours: Mapping[int, set[int]], rings: int) -> set[int]:
    """Seeds plus ``rings`` layers of edge adjacency."""
    selected = set(seeds)
    frontier = set(seeds)
    for _ in range(rings):
        grown: set[int] = set()
        for index in frontier:
            grown |= neighbours[index] - selected
        if not grown:
            break
        selected |= grown
        frontier = grown
    return selected


def _components(indices: Sequence[int], neighbours: Mapping[int, set[int]]) -> list[list[int]]:
    """Connected components of the face subset ``indices``."""
    remaining = set(indices)
    out: list[list[int]] = []
    while remaining:
        start = min(remaining)
        remaining.discard(start)
        group = [start]
        queue = deque([start])
        while queue:
            current = queue.popleft()
            for other in neighbours[current]:
                if other in remaining:
                    remaining.discard(other)
                    group.append(other)
                    queue.append(other)
        out.append(sorted(group))
    return out


def _submesh(mesh: Mesh, face_indices: Sequence[int]) -> Mesh:
    """The mesh restricted to ``face_indices`` (kept in order). Vertices are not renumbered."""
    return Mesh(
        vertices=mesh.vertices,
        faces=tuple(mesh.faces[index] for index in face_indices),
        units=mesh.units,
        coordinate_frame=mesh.coordinate_frame,
    )


def _fan_cap(mesh: Mesh, loops: Sequence[Sequence[int]]) -> tuple[Mesh, int]:
    """Cap each loop with a triangle fan from the loop's own centroid.

    The fallback for a loop whose projection onto its own plane is not a simple polygon --
    an excavated loop around a fold usually is not, because the facet normals it walks over
    turn through more than 90 degrees. A centroid fan has no such requirement: it needs only
    the loop itself, and every triangle it adds is non-degenerate as long as the loop is not
    collinear with its centroid.
    """
    vertices = list(mesh.vertices)
    faces = list(mesh.faces)
    capped = 0
    for loop in loops:
        points = [mesh.vertices[index] for index in loop]
        centroid = tuple(
            sum(point[axis] for point in points) / len(points) for axis in range(3)
        )
        centre = len(vertices)
        vertices.append(centroid)
        walk = list(reversed(list(loop)))
        for k in range(len(walk)):
            face = (walk[k], walk[(k + 1) % len(walk)], centre)
            a, b, c = (vertices[index] for index in face)
            ux, uy, uz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
            vx, vy, vz = c[0] - a[0], c[1] - a[1], c[2] - a[2]
            nx, ny, nz = uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx
            if math.sqrt(nx * nx + ny * ny + nz * nz) <= 0.0:
                raise MeshRepairError(
                    f"the centroid fan over a {len(loop)}-edge loop produced a zero-area facet"
                )
            faces.append(face)
        capped += 1
    return (
        Mesh(
            vertices=tuple(vertices),
            faces=tuple(faces),
            units=mesh.units,
            coordinate_frame=mesh.coordinate_frame,
        ),
        capped,
    )


# ---------------------------------------------------------------------------------------
# the result
# ---------------------------------------------------------------------------------------


@dataclass(frozen=True)
class ExcavationBand:
    """One attempt at cutting the fold away, with the residual it left behind."""

    rings: int
    faces_removed: int
    residual_pairs: int
    capped: bool = False

    def to_dict(self) -> dict[str, Any]:
        return {
            "rings": self.rings,
            "faces_removed": self.faces_removed,
            "residual_pairs": self.residual_pairs,
            "capped": self.capped,
        }


@dataclass(frozen=True)
class DetangleResult:
    """The healed mesh (or the original, on refusal) plus the evidence for either."""

    mesh: Mesh
    healed: bool
    reason: str
    before: MeshMetrics
    after: MeshMetrics
    pairs_found: int
    bands: tuple[ExcavationBand, ...] = ()
    rings_used: int = 0
    faces_removed: int = 0
    fragments_dropped: int = 0
    dropped_faces: int = 0
    holes: HoleFill | None = None
    volume_before: float = 0.0
    volume_after: float = 0.0
    volume_relative_change: float | None = None
    reference_volume: float | None = None
    tolerated_volume_change: float = DEFAULT_VOLUME_TOLERANCE
    notes: tuple[str, ...] = field(default_factory=tuple)

    @property
    def changed(self) -> bool:
        return self.healed and self.faces_removed > 0

    def to_dict(self) -> dict[str, Any]:
        return {
            "healed": self.healed,
            "changed": self.changed,
            "reason": self.reason,
            "pairs_found": self.pairs_found,
            "rings_used": self.rings_used,
            "faces_removed": self.faces_removed,
            "fragments_dropped": self.fragments_dropped,
            "dropped_faces": self.dropped_faces,
            "bands": [band.to_dict() for band in self.bands],
            "before": self.before.to_dict(),
            "after": self.after.to_dict(),
            "volume_before_m3": self.volume_before,
            "volume_after_m3": self.volume_after,
            "volume_relative_change": self.volume_relative_change,
            "reference_volume_m3": self.reference_volume,
            "tolerated_volume_change": self.tolerated_volume_change,
            "holes": self.holes.to_dict() if self.holes else None,
            "notes": list(self.notes),
        }


def _relative_change(before: float, after: float) -> float | None:
    if before == 0.0:
        return None
    return abs(after - before) / abs(before)


@dataclass(frozen=True)
class _Attempt:
    """One cut-and-cap attempt, accepted or not."""

    ok: bool
    reason: str
    mesh: Mesh | None = None
    after: MeshMetrics | None = None
    holes: HoleFill | None = None
    fragments_dropped: int = 0
    dropped_faces: int = 0
    volume_after: float = 0.0
    volume_relative_change: float | None = None
    reference_volume: float | None = None
    notes: tuple[str, ...] = ()


def _cap_attempt(
    welded: Mesh,
    keep_index: Sequence[int],
    neighbours: Mapping[int, set[int]],
    shell_of: Mapping[int, int],
    *,
    max_hole_edges: int,
    self_intersection_limit: int,
    reference_volume: float | None,
    volume_before: float,
) -> _Attempt:
    """Cut the kept facets out, drop shards, cap every hole, and measure the result."""
    notes: list[str] = []

    pieces: dict[int, list[list[int]]] = defaultdict(list)
    for piece in _components(keep_index, neighbours):
        pieces[shell_of[piece[0]]].append(piece)
    keep: list[int] = []
    fragments_dropped = 0
    dropped_faces = 0
    for shell in sorted(pieces):
        survivors = pieces[shell]
        biggest = max(len(piece) for piece in survivors)
        for piece in survivors:
            if len(piece) < biggest:
                fragments_dropped += 1
                dropped_faces += len(piece)
                notes.append(
                    f"dropped a {len(piece)}-facet fragment detached from shell {shell} "
                    f"(the largest surviving piece of that shell keeps {biggest} facets)"
                )
            else:
                keep.extend(piece)
    keep.sort()
    if not keep:
        return _Attempt(
            ok=False,
            reason="refusing to excavate: no shell survives the cut",
            fragments_dropped=fragments_dropped,
            dropped_faces=dropped_faces,
            notes=tuple(notes),
        )

    holes: HoleFill = fill_holes(_submesh(welded, keep), max_hole_edges=max_hole_edges)
    if holes.refused:
        try:
            loops = boundary_loops(holes.mesh)
            fanned_mesh, fanned_loops = _fan_cap(holes.mesh, loops)
        except MeshRepairError as error:
            return _Attempt(
                ok=False,
                reason=(
                    f"the {len(holes.refused)} excavated boundary loop(s) could not be capped: "
                    f"the planar fill refused them ({holes.refused[0].get('reason')}) and the "
                    f"centroid fan could not close them either ({error})"
                ),
                holes=holes,
                fragments_dropped=fragments_dropped,
                dropped_faces=dropped_faces,
                notes=tuple(notes),
            )
        notes.append(
            f"the planar fill refused {len(holes.refused)} loop(s) "
            f"({holes.refused[0].get('reason')}); {fanned_loops} of them were closed with a "
            "centroid fan instead"
        )
        holes = HoleFill(
            mesh=fanned_mesh,
            holes=holes.holes
            + tuple(
                {"loop_vertices": len(loop), "triangles": len(loop), "fill": "centroid_fan"}
                for loop in loops
            ),
            step=RepairStep(
                "fill_holes",
                fanned_loops,
                f"{fanned_loops} loop(s) closed with a centroid fan after the planar fill "
                "was refused",
            ),
        )
    try:
        leftovers = boundary_loops(holes.mesh)
    except MeshRepairError as error:  # pragma: no cover - a manifoldisation failure
        return _Attempt(
            ok=False,
            reason=f"the capped mesh still has an unwalkable boundary: {error}",
            holes=holes,
            fragments_dropped=fragments_dropped,
            dropped_faces=dropped_faces,
            notes=tuple(notes),
        )
    if leftovers:
        return _Attempt(
            ok=False,
            reason=(
                f"{len(leftovers)} excavated boundary loop(s) stay open after the fill "
                f"(the {max_hole_edges}-edge cap limit is the usual cause)"
            ),
            holes=holes,
            fragments_dropped=fragments_dropped,
            dropped_faces=dropped_faces,
            notes=tuple(notes),
        )

    oriented = meshlib.orient_outward(holes.mesh)
    after = measure_mesh(oriented.mesh, self_intersection_limit=self_intersection_limit)
    volume_after = meshlib.total_signed_volume(oriented.mesh)
    reference = reference_volume if reference_volume is not None else volume_before
    change = _relative_change(reference, volume_after)
    base = dict(
        after=after,
        holes=holes,
        fragments_dropped=fragments_dropped,
        dropped_faces=dropped_faces,
        volume_after=volume_after,
        volume_relative_change=change,
        reference_volume=reference,
        notes=tuple(notes),
    )
    if not (after.closed and after.manifold and after.consistently_wound):
        return _Attempt(
            ok=False,
            reason=(
                "the capped mesh is not clean "
                f"(open edges {after.open_edges}, non-manifold edges "
                f"{after.nonmanifold_edges}, winding conflicts {after.winding_conflicts})"
            ),
            **base,
        )
    if after.self_intersections:
        return _Attempt(
            ok=False,
            reason=(
                f"the cap introduces {after.self_intersections} new penetrating pair(s); the "
                "excavated region is too curved for the fill available here"
            ),
            **base,
        )
    return _Attempt(ok=True, reason="clean after cut and cap", mesh=oriented.mesh, **base)


def excavation_faces(
    mesh: Mesh, *, rings: int = 1, self_intersection_limit: int = DEFAULT_SELF_INTERSECTION_LIMIT
) -> frozenset[int]:
    """The facets that would be cut for ``rings`` layers around the detected folds.

    Exposed so a caller can inspect (or plot) the band *before* committing to a repair.
    """
    if rings < 1:
        raise ValueError("rings must be at least 1")
    pairs = find_self_intersections(mesh, limit=self_intersection_limit)
    if not pairs:
        return frozenset()
    seeds = {index for pair in pairs for index in pair}
    return frozenset(_dilate(sorted(seeds), _face_neighbours(mesh), rings))


# ---------------------------------------------------------------------------------------
# the repair
# ---------------------------------------------------------------------------------------


def heal_self_intersections(
    mesh: Mesh,
    *,
    max_rings: int = DEFAULT_MAX_RINGS,
    self_intersection_limit: int = DEFAULT_SELF_INTERSECTION_LIMIT,
    max_hole_edges: int = DEFAULT_MAX_HOLE_EDGES,
    reference_volume: float | None = None,
    volume_tolerance: float | None = None,
    weld_tolerance: float = DEFAULT_WELD_TOLERANCE,
) -> DetangleResult:
    """Remove local self intersections by cutting the fold away and capping the hole.

    ``volume_tolerance`` defaults to :data:`DEFAULT_VOLUME_TOLERANCE` (1 %, the ``GAP-CAD-009``
    budget) when a ``reference_volume`` is supplied, and to
    :data:`DEFAULT_MIN_VOLUME_TOLERANCE` (5 %) otherwise -- without a reference the input's own
    signed volume is the only yardstick, and on a self-intersecting mesh that value is itself
    the thing under suspicion, so it is checked more loosely and the note says so.
    """
    if max_rings < 1:
        raise ValueError("max_rings must be at least 1")
    if max_hole_edges < 3:
        raise ValueError("max_hole_edges must be at least 3")
    if volume_tolerance is not None and volume_tolerance < 0.0:
        raise ValueError("volume_tolerance must be non-negative")
    notes: list[str] = []

    welded = meshlib.weld(mesh, tolerance=weld_tolerance)
    if len(welded.vertices) != len(mesh.vertices):
        notes.append(
            f"welded {len(mesh.vertices) - len(welded.vertices)} coincident vertex (or "
            f"vertices) at {weld_tolerance:g} m before measuring"
        )
    before = measure_mesh(welded, self_intersection_limit=self_intersection_limit)
    volume_before = meshlib.total_signed_volume(welded)
    effective_tolerance = (
        volume_tolerance
        if volume_tolerance is not None
        else (
            DEFAULT_VOLUME_TOLERANCE
            if reference_volume is not None
            else DEFAULT_MIN_VOLUME_TOLERANCE
        )
    )

    def refuse(reason: str, holes: HoleFill | None = None) -> DetangleResult:
        return DetangleResult(
            mesh=mesh,
            healed=False,
            reason=reason,
            before=before,
            after=before,
            pairs_found=before.self_intersections,
            holes=holes,
            volume_before=volume_before,
            volume_after=volume_before,
            volume_relative_change=None,
            reference_volume=reference_volume,
            tolerated_volume_change=effective_tolerance,
            notes=tuple(notes),
        )

    pairs = find_self_intersections(welded, limit=self_intersection_limit)
    if not pairs:
        return refuse(
            "nothing to heal: no penetrating facet pair within the detection limit of "
            f"{self_intersection_limit} (the input's own upper bound is "
            f"{before.self_intersections})"
        )
    if not before.closed:
        return refuse(
            f"refusing to excavate a mesh with {before.open_edges} boundary edge(s): the "
            "repair's own holes would be indistinguishable from the mesh's, so a closed input "
            "is required"
        )

    neighbours = _face_neighbours(welded)
    shells = _components(list(range(len(welded.faces))), neighbours)
    shell_of: dict[int, int] = {}
    for shell_index, shell in enumerate(shells):
        for face_index in shell:
            shell_of[face_index] = shell_index

    if reference_volume is None:
        notes.append(
            "no reference volume supplied: the input's own signed volume is the yardstick, at "
            f"the looser {effective_tolerance:.1%} budget, because on a self-intersecting mesh "
            "that volume is itself the suspect quantity"
        )
    seeds = sorted({index for pair in pairs for index in pair})
    bands: list[ExcavationBand] = []
    last_reason = "the excavation never produced a cappable surface"
    stopped: _Attempt | None = None
    for rings in range(1, max_rings + 1):
        removed = _dilate(seeds, neighbours, rings)
        kept = [index for index in range(len(welded.faces)) if index not in removed]
        if not kept:
            return refuse("refusing to excavate: every facet is inside the excavation band")
        residual = find_self_intersections(_submesh(welded, kept), limit=self_intersection_limit)
        if residual:
            bands.append(
                ExcavationBand(
                    rings=rings, faces_removed=len(removed), residual_pairs=len(residual)
                )
            )
            last_reason = (
                f"{len(residual)} penetrating pair(s) survive after cutting {rings} ring(s) "
                f"({len(removed)} facet(s)); the fold is wider than that excavation"
            )
            continue
        attempt = _cap_attempt(
            welded,
            kept,
            neighbours,
            shell_of,
            max_hole_edges=max_hole_edges,
            self_intersection_limit=self_intersection_limit,
            reference_volume=reference_volume,
            volume_before=volume_before,
        )
        bands.append(
            ExcavationBand(
                rings=rings,
                faces_removed=len(removed),
                residual_pairs=0,
                capped=attempt.ok,
            )
        )
        last_reason = attempt.reason
        notes.extend(attempt.notes)
        if not attempt.ok:
            continue
        if (
            attempt.volume_relative_change is not None
            and attempt.volume_relative_change > effective_tolerance
        ):
            stopped = attempt
            last_reason = (
                f"the repair moves the enclosed volume by {attempt.volume_relative_change:.3%}, "
                f"above the {effective_tolerance:.1%} budget ({volume_before:.6f} -> "
                f"{attempt.volume_after:.6f} m^3): the cut was not local enough to be called a "
                "repair, and a wider cut only removes more material"
            )
            break
        return DetangleResult(
            mesh=attempt.mesh,
            healed=True,
            reason=(
                f"cut {len(removed)} facet(s) in {rings} ring(s), dropped "
                f"{attempt.fragments_dropped} fragment(s), capped {attempt.holes.filled} loop(s); "
                f"{len(pairs)} penetrating pair(s) -> {attempt.after.self_intersections}"
            ),
            before=before,
            after=attempt.after,
            pairs_found=len(pairs),
            bands=tuple(bands),
            rings_used=rings,
            faces_removed=len(removed),
            fragments_dropped=attempt.fragments_dropped,
            dropped_faces=attempt.dropped_faces,
            holes=attempt.holes,
            volume_before=volume_before,
            volume_after=attempt.volume_after,
            volume_relative_change=attempt.volume_relative_change,
            reference_volume=attempt.reference_volume,
            tolerated_volume_change=effective_tolerance,
            notes=tuple(notes),
        )

    # A refusal still carries the best attempt's measurements -- the reason quotes those
    # numbers, so they belong in the structured record too, not only in the sentence.
    return DetangleResult(
        mesh=mesh,
        healed=False,
        reason=last_reason,
        before=before,
        after=stopped.after if stopped and stopped.after is not None else before,
        pairs_found=len(pairs),
        bands=tuple(bands),
        rings_used=bands[-1].rings if bands else 0,
        faces_removed=bands[-1].faces_removed if bands else 0,
        fragments_dropped=stopped.fragments_dropped if stopped else 0,
        dropped_faces=stopped.dropped_faces if stopped else 0,
        holes=stopped.holes if stopped else None,
        volume_before=volume_before,
        volume_after=stopped.volume_after if stopped else volume_before,
        volume_relative_change=stopped.volume_relative_change if stopped else None,
        reference_volume=(stopped.reference_volume if stopped else None) or reference_volume,
        tolerated_volume_change=effective_tolerance,
        notes=tuple(notes),
    )
