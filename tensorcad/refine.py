"""Target-size-field mesh refinement: the half of remeshing conditioning was missing.

:mod:`tensorcad.mesh_conditioning` moves vertices along the surface (Taubin), puts them
back on it (projection) and removes them (quadric-error decimation). Every one of those
steps *coarsens*: after the pipeline the mesh has the same vertices it started with or
fewer. A solver that needs a finer discretisation than the sampler produced -- a BEM
panel model that wants a tenth of a wavelength near the waterline, a FEM surface that
wants elements below the local radius of curvature -- had no route back up, and the
README says so ("remeshing to a target size field: mesh conditioning only coarsens").

This module is that route. It splits triangles, it never moves an existing vertex, and
it reports what splitting did to the geometry:

* :func:`refine` -- split every edge longer than the local target size until none is
  left, and return the mesh together with the before/after measurements;
* :class:`RefinementReport` -- vertices, faces, edge-length extremes, worst angle, worst
  aspect, volume, area, centroid, watertightness, edge-manifoldness, Euler
  characteristic, and the projection residuals when a field is supplied;
* :class:`StageMetrics` -- the same block for one stage, so the report can be compared
  against ``measure_mesh`` from the conditioning pipeline.

**Why the refinement is exact.** A new vertex is the midpoint of an existing edge, so
every sub-triangle lies in the plane of its parent, and the signed volume, the surface
area and the volume centroid are invariant under the split -- not approximately, to
round-off. That is asserted in the tests at 1e-13 relative, and it is the honest reading
of the operation: **refinement adds resolution, it does not add accuracy.** A coarse
faceted hull stays equally wrong at four times the face count. Accuracy comes from the
optional ``project`` argument, which moves the *new* vertices onto the analytic surface
with a bounded Newton step; the report then carries the residual ``|f|`` before and
after and the displacement that was paid for it.

**Conformity.** Edges are marked by length, and every face is split according to how
many of its three edges are marked (the 0/1/2/3-edge cases are the classical red-green
patterns). A marked edge has exactly one midpoint and both of its faces use it, so no
hanging node can appear: watertightness, edge-manifoldness and the Euler characteristic
of the input survive the operation. A mesh that is *not* watertight stays not watertight
in the same way -- the report shows both states rather than silently welding cracks
(coincident input vertices are welded by the builder and the count is reported as
``welded_positions``, because that is a topology change worth knowing about).

**What it is not.** No curvature sensor (the target is the caller's field, uniform or
``(x, y, z) -> h``, so a wavelength or a curvature estimate is the caller's business),
no anisotropic splitting, no surface re-parametrisation, and no projection of the
*existing* vertices. Splitting alone never improves geometric accuracy, and the
docstring of :func:`refine` says so where a caller will read it.

Only the Python standard library is required. A field to project onto is duck-typed on
:class:`tensorcad.geometry_field.GeometryField` (``signed_distance`` and ``gradient``),
so no import of the field module is needed here.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Any, Protocol

from . import mesh as meshlib
from .geometry import Face, Mesh, MeshBuilder, Vec3

__all__ = [
    "RefinementError",
    "RefinementReport",
    "StageMetrics",
    "edge_lengths",
    "refine",
]

EdgeKey = tuple[int, int]


class RefinementError(ValueError):
    """Raised when the requested refinement cannot be performed."""


class _Projector(Protocol):
    """The slice of :class:`GeometryField` that projection needs."""

    def signed_distance(self, points: Sequence[Sequence[float]]) -> list[float]: ...

    def gradient(self, points: Sequence[Sequence[float]]) -> list[tuple[float, float, float]]: ...


# --------------------------------------------------------------------------- #
# Small vector helpers (standard library only)
# --------------------------------------------------------------------------- #


def _sub(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _cross(a: Vec3, b: Vec3) -> Vec3:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _norm(a: Vec3) -> float:
    return math.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])


def _midpoint(a: Vec3, b: Vec3) -> Vec3:
    """Exact midpoint; the same expression is used for both faces of an edge."""
    return ((a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5, (a[2] + b[2]) * 0.5)


def _edge_key(a: int, b: int) -> EdgeKey:
    return (a, b) if a < b else (b, a)


def _finite_positive(name: str, value: Any) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise RefinementError(f"{name} must be a finite positive number, got {value!r}")
    return number


# --------------------------------------------------------------------------- #
# Measurements
# --------------------------------------------------------------------------- #


def _volume_centroid(mesh: Mesh) -> Vec3:
    """Volume centroid of a closed surface, by signed tetrahedra about the origin."""
    total = 0.0
    sx = sy = sz = 0.0
    for face in mesh.faces:
        a, b, c = (mesh.vertices[face[0]], mesh.vertices[face[1]], mesh.vertices[face[2]])
        det = (
            a[0] * (b[1] * c[2] - b[2] * c[1])
            + a[1] * (b[2] * c[0] - b[0] * c[2])
            + a[2] * (b[0] * c[1] - b[1] * c[0])
        )
        total += det
        sx += det * (a[0] + b[0] + c[0]) * 0.25
        sy += det * (a[1] + b[1] + c[1]) * 0.25
        sz += det * (a[2] + b[2] + c[2]) * 0.25
    if abs(total) <= 1e-300:
        return (math.nan, math.nan, math.nan)
    return (sx / total, sy / total, sz / total)


def _face_angles_and_aspect(mesh: Mesh) -> tuple[float, float]:
    """Smallest angle in degrees and largest longest-edge/shortest-altitude ratio."""
    min_angle = math.inf
    max_aspect = 0.0
    for face in mesh.faces:
        a, b, c = (mesh.vertices[face[0]], mesh.vertices[face[1]], mesh.vertices[face[2]])
        for first, second, opposite in ((a, b, c), (b, c, a), (c, a, b)):
            u = _norm(_sub(second, first))
            v = _norm(_sub(opposite, first))
            w = _norm(_sub(opposite, second))
            if u <= 0.0 or v <= 0.0 or w <= 0.0:
                continue
            cosine = (u * u + v * v - w * w) / (2.0 * u * v)
            angle = math.degrees(math.acos(min(1.0, max(-1.0, cosine))))
            min_angle = min(min_angle, angle)
            area2 = _norm(_cross(_sub(b, a), _sub(c, a)))
            if area2 > 0.0:
                longest = max(u, v, w)
                # area = 0.5 * longest * altitude -> longest altitude ratio
                max_aspect = max(max_aspect, longest * longest / area2)
    if not math.isfinite(min_angle):
        min_angle = math.nan
    return min_angle, max_aspect


def edge_lengths(mesh: Mesh) -> dict[EdgeKey, float]:
    """Length of every undirected edge, keyed by its sorted vertex pair."""
    lengths: dict[EdgeKey, float] = {}
    for face in mesh.faces:
        for first, second in (
            (face[0], face[1]),
            (face[1], face[2]),
            (face[2], face[0]),
        ):
            key = _edge_key(first, second)
            if key not in lengths:
                lengths[key] = _norm(_sub(mesh.vertices[key[1]], mesh.vertices[key[0]]))
    return lengths


@dataclass(frozen=True)
class StageMetrics:
    """One mesh state, measured the way the conditioning pipeline measures it."""

    vertices: int
    faces: int
    max_edge: float
    mean_edge: float
    min_angle_deg: float
    max_aspect: float
    volume: float
    area: float
    centroid: Vec3
    watertight: bool
    edge_manifold: bool
    euler_characteristic: int

    def to_dict(self) -> dict[str, Any]:
        return {
            "vertices": self.vertices,
            "faces": self.faces,
            "max_edge": self.max_edge,
            "mean_edge": self.mean_edge,
            "min_angle_deg": self.min_angle_deg,
            "max_aspect": self.max_aspect,
            "volume": self.volume,
            "area": self.area,
            "centroid": self.centroid,
            "watertight": self.watertight,
            "edge_manifold": self.edge_manifold,
            "euler_characteristic": self.euler_characteristic,
        }


def _stage_metrics(mesh: Mesh) -> StageMetrics:
    lengths = edge_lengths(mesh)
    min_angle, max_aspect = _face_angles_and_aspect(mesh)
    counts = meshlib.edge_face_counts(mesh)
    values = tuple(lengths.values())
    return StageMetrics(
        vertices=len(mesh.vertices),
        faces=len(mesh.faces),
        max_edge=max(values) if values else math.nan,
        mean_edge=math.fsum(values) / len(values) if values else math.nan,
        min_angle_deg=min_angle,
        max_aspect=max_aspect,
        volume=meshlib.total_signed_volume(mesh),
        area=meshlib.surface_area(mesh),
        centroid=_volume_centroid(mesh),
        watertight=meshlib.is_watertight(mesh),
        edge_manifold=all(count <= 2 for count in counts.values()),
        euler_characteristic=len(mesh.vertices) - len(counts) + len(mesh.faces),
    )


# --------------------------------------------------------------------------- #
# Target size field
# --------------------------------------------------------------------------- #


def _resolve_target(target: float | Callable[[Vec3], float]) -> Callable[[Vec3], float]:
    if callable(target):
        probe = target((0.0, 0.0, 0.0))
        _finite_positive("target size field", probe)

        def evaluate(point: Vec3) -> float:
            return _finite_positive("target size field", target(point))

        return evaluate
    constant = _finite_positive("target_edge_length", target)
    return lambda point: constant


# --------------------------------------------------------------------------- #
# Splitting
# --------------------------------------------------------------------------- #


def _split_face(
    builder: MeshBuilder,
    face: Face,
    midpoints: list[int | None],
    flags: list[bool],
) -> None:
    """Emit the red-green pattern for one face.

    ``midpoints[i]`` is the new-vertex index of edge ``i`` (0: a-b, 1: b-c, 2: c-a)
    and ``flags[i]`` says whether that edge was marked. All patterns keep the parent
    winding, so the output stays consistently oriented.
    """
    a, b, c = face
    mask = (1 if flags[0] else 0) | (2 if flags[1] else 0) | (4 if flags[2] else 0)
    if mask == 0:
        builder.add_triangle(a, b, c)
        return
    p = midpoints[0] if flags[0] else None
    q = midpoints[1] if flags[1] else None
    r = midpoints[2] if flags[2] else None
    if mask == 7:
        assert p is not None and q is not None and r is not None
        builder.add_triangle(a, p, r)
        builder.add_triangle(p, b, q)
        builder.add_triangle(r, q, c)
        builder.add_triangle(p, q, r)
    elif mask == 1:
        assert p is not None
        builder.add_triangle(a, p, c)
        builder.add_triangle(p, b, c)
    elif mask == 2:
        assert q is not None
        builder.add_triangle(a, b, q)
        builder.add_triangle(a, q, c)
    elif mask == 4:
        assert r is not None
        builder.add_triangle(a, b, r)
        builder.add_triangle(r, b, c)
    elif mask == 3:
        assert p is not None and q is not None
        builder.add_triangle(p, b, q)
        builder.add_triangle(a, p, q)
        builder.add_triangle(a, q, c)
    elif mask == 6:
        assert q is not None and r is not None
        builder.add_triangle(a, b, q)
        builder.add_triangle(a, q, r)
        builder.add_triangle(c, r, q)
    elif mask == 5:
        assert p is not None and r is not None
        builder.add_triangle(a, p, r)
        builder.add_triangle(p, b, c)
        builder.add_triangle(p, c, r)
    else:  # pragma: no cover - masks are exhaustive above
        raise RefinementError(f"unhandled edge mask {mask}")


def _split_pass(
    mesh: Mesh,
    target: Callable[[Vec3], float],
) -> tuple[Mesh, list[float], int]:
    """One conforming split of every edge longer than its local target.

    Returns the new mesh, the parent edge length of every midpoint in the order the
    builder created them, and the number of coincident input positions it welded.
    """
    builder = MeshBuilder(mesh.units, mesh.coordinate_frame)
    original = [builder.add_vertex(point) for point in mesh.vertices]
    welded = len(mesh.vertices) - builder.vertex_count
    midpoint_index: dict[EdgeKey, int] = {}
    midpoint_lengths: list[float] = []

    def midpoint(key: EdgeKey, length: float) -> int:
        index = midpoint_index.get(key)
        if index is None:
            point = _midpoint(mesh.vertices[key[0]], mesh.vertices[key[1]])
            index = builder.add_vertex(point)
            midpoint_index[key] = index
            midpoint_lengths.append(length)
        return index

    for face in mesh.faces:
        keys = (
            _edge_key(face[0], face[1]),
            _edge_key(face[1], face[2]),
            _edge_key(face[2], face[0]),
        )
        flags: list[bool] = []
        indices: list[int | None] = []
        for key in keys:
            first, second = mesh.vertices[key[0]], mesh.vertices[key[1]]
            length = _norm(_sub(second, first))
            limit = target(_midpoint(first, second))
            if length > limit:
                flags.append(True)
                indices.append(midpoint(key, length))
            else:
                flags.append(False)
                indices.append(None)
        _split_face(
            builder,
            (original[face[0]], original[face[1]], original[face[2]]),
            indices,
            flags,
        )
    return builder.build(), midpoint_lengths, welded


# --------------------------------------------------------------------------- #
# Projection onto an analytic surface
# --------------------------------------------------------------------------- #


def _project_new_vertices(
    field: _Projector,
    vertices: list[Vec3],
    start: int,
    iterations: int,
    bounds: Sequence[float],
) -> tuple[list[Vec3], dict[str, float | int]]:
    """Newton-project ``vertices[start:]`` onto ``field``'s zero level set.

    Each step is ``p - f(p)/|grad f|^2 * grad f``, clipped to ``bounds[i]`` (a
    multiple of that vertex's own parent edge) and backtracked when it fails to reduce
    ``|f|``; a vertex whose step cannot be accepted is left where it was and the report
    counts it.
    """
    points = list(vertices)
    targets = points[start:]
    if not targets:
        return points, {
            "vertices_projected": 0,
            "max_displacement": 0.0,
            "mean_displacement": 0.0,
            "max_residual_before": 0.0,
            "max_residual_after": 0.0,
            "stalled": 0,
        }
    before = [abs(value) for value in field.signed_distance(targets)]
    displacement = [0.0] * len(targets)
    stalled = 0
    for _ in range(max(0, iterations)):
        values = field.signed_distance(targets)
        gradients = field.gradient(targets)
        candidates: list[Vec3] = []
        for index, (point, value, gradient) in enumerate(zip(targets, values, gradients, strict=True)):
            norm2 = gradient[0] ** 2 + gradient[1] ** 2 + gradient[2] ** 2
            if norm2 <= 1e-300 or not math.isfinite(norm2):
                candidates.append(point)
                continue
            step = -value / norm2
            length = abs(step) * math.sqrt(norm2)
            max_displacement = bounds[index] if index < len(bounds) else 0.0
            if length > max_displacement > 0.0:
                step *= max_displacement / length
            candidates.append(
                (
                    point[0] + step * gradient[0],
                    point[1] + step * gradient[1],
                    point[2] + step * gradient[2],
                )
            )
        candidate_values = [abs(value) for value in field.signed_distance(candidates)]
        for index, (point, candidate, candidate_value) in enumerate(
            zip(targets, candidates, candidate_values, strict=True)
        ):
            if candidate_value <= abs(values[index]):
                displacement[index] += _norm(_sub(candidate, point))
                targets[index] = candidate
            else:
                stalled += 1
    after = [abs(value) for value in field.signed_distance(targets)]
    points[start:] = targets
    return points, {
        "vertices_projected": len(targets),
        "max_displacement": max(displacement) if displacement else 0.0,
        "mean_displacement": math.fsum(displacement) / len(displacement) if displacement else 0.0,
        "max_residual_before": max(before) if before else 0.0,
        "max_residual_after": max(after) if after else 0.0,
        "stalled": stalled,
    }


# --------------------------------------------------------------------------- #
# The report and the entry point
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class RefinementReport:
    """What the refinement did, measured on both ends.

    ``converged`` means *every* edge is at or below its local target (the tolerance
    below is a relative slack on that comparison, not a licence to stop early); when
    ``max_passes`` runs out first the report says ``False`` and still prints the
    largest remaining edge, so a caller cannot mistake a truncated run for a finished
    one.
    """

    passes: int
    converged: bool
    target_edge_length: float | None
    target_callable: bool
    target_range: tuple[float, float] | None
    tolerance: float
    before: StageMetrics
    after: StageMetrics
    projection: dict[str, float | int] | None
    welded_positions: int

    @property
    def volume_relative_change(self) -> float:
        if self.before.volume == 0.0:
            return math.nan
        return (self.after.volume - self.before.volume) / self.before.volume

    @property
    def area_relative_change(self) -> float:
        if self.before.area == 0.0:
            return math.nan
        return (self.after.area - self.before.area) / self.before.area

    @property
    def centroid_shift(self) -> float:
        before, after = self.before.centroid, self.after.centroid
        if not all(math.isfinite(value) for value in before + after):
            return math.nan
        return _norm(_sub(after, before))

    def to_dict(self) -> dict[str, Any]:
        return {
            "passes": self.passes,
            "converged": self.converged,
            "target_edge_length": self.target_edge_length,
            "target_callable": self.target_callable,
            "target_range": self.target_range,
            "tolerance": self.tolerance,
            "before": self.before.to_dict(),
            "after": self.after.to_dict(),
            "projection": self.projection,
            "welded_positions": self.welded_positions,
            "volume_relative_change": self.volume_relative_change,
            "area_relative_change": self.area_relative_change,
            "centroid_shift": self.centroid_shift,
        }


def refine(
    mesh: Mesh,
    target: float | Callable[[Vec3], float],
    *,
    max_passes: int | None = None,
    tolerance: float = 0.0,
    project: _Projector | None = None,
    project_iterations: int = 2,
    max_faces: int = 4_000_000,
) -> tuple[Mesh, RefinementReport]:
    """Split every edge longer than its local target, conformingly, and report it.

    ``target`` is either a length in the mesh's units or a callable
    ``(x, y, z) -> h`` -- the size field. An edge is split when
    ``length > h(midpoint) * (1 + tolerance)``; the tolerance exists so a caller can
    stop one hair above the target instead of chasing float noise.

    Splitting is **geometry-neutral**: new vertices are edge midpoints, so volume,
    area and centroid are invariant to round-off and the only thing that changes is
    the resolution. Pass ``project`` (any object with ``signed_distance`` and
    ``gradient``, e.g. a :class:`tensorcad.geometry_field.GeometryField`) to move the
    new vertices onto the analytic surface; the report then carries the residual
    before/after and the displacement, and *that* is where accuracy comes from.

    ``max_passes`` defaults to what the worst edge needs (``log2(max_edge/min h)``)
    plus two; ``max_faces`` refuses a refinement that would blow up rather than
    exhausting memory.
    """
    if not isinstance(mesh, Mesh):
        raise RefinementError(f"expected a Mesh, got {type(mesh).__name__}")
    if not mesh.faces:
        raise RefinementError("cannot refine a mesh with no faces")
    if tolerance < 0.0 or not math.isfinite(tolerance):
        raise RefinementError(f"tolerance must be finite and non-negative, got {tolerance!r}")
    for index, point in enumerate(mesh.vertices):
        if not all(math.isfinite(axis) for axis in point):
            raise RefinementError(f"mesh vertex {index} is not finite: {point!r}")
    size_field = _resolve_target(target)
    slack = 1.0 + tolerance

    lengths = edge_lengths(mesh)
    worst = max(lengths.values())
    sampled = [size_field(point) for point in mesh.vertices]
    smallest = min(sampled)
    largest = max(sampled)
    if max_passes is None:
        # One pass halves a marked edge at best, so this many passes is enough for
        # the uniform case; adaptive fields can need a couple more and the loop
        # still stops as soon as nothing is marked.
        needed = math.ceil(math.log2(worst / smallest)) + 2 if worst > smallest else 1
        max_passes = max(1, min(64, needed))
    if max_passes < 1:
        raise RefinementError(f"max_passes must be at least 1, got {max_passes!r}")

    before = _stage_metrics(mesh)
    current = mesh
    passes = 0
    converged = False
    welded = 0
    projection: dict[str, float | int] | None = None
    for _ in range(max_passes):
        passes += 1
        start = len(current.vertices)
        split_mesh, parent_lengths, new_welded = _split_pass(current, size_field)
        welded += new_welded
        if not parent_lengths:
            current = split_mesh
            converged = True
            break
        if len(split_mesh.faces) > max_faces:
            raise RefinementError(
                f"refinement would produce {len(split_mesh.faces)} faces, above max_faces={max_faces}"
            )
        if project is not None:
            points = list(split_mesh.vertices)
            if len(points) > start:
                # Each new vertex may travel at most a few of its own parent edges;
                # that is the bound that stops a Newton step on a bad gradient from
                # flinging a vertex across the model.
                if len(points) == start + len(parent_lengths):
                    bounds = [4.0 * length for length in parent_lengths]
                else:  # pragma: no cover - only if the builder welded midpoints away
                    bounds = [4.0 * largest] * (len(points) - start)
                points, step_report = _project_new_vertices(
                    project,
                    points,
                    start,
                    project_iterations,
                    bounds,
                )
                split_mesh = Mesh(
                    vertices=tuple(points),
                    faces=split_mesh.faces,
                    units=split_mesh.units,
                    coordinate_frame=split_mesh.coordinate_frame,
                )
                if projection is None:
                    projection = step_report
                else:
                    for key, value in step_report.items():
                        if isinstance(value, float) and key.startswith("max_"):
                            projection[key] = max(float(projection.get(key, 0.0)), value)
                        elif key == "vertices_projected" or key == "stalled":
                            projection[key] = int(projection.get(key, 0)) + int(value)
        current = split_mesh
    if not converged:
        # The pass budget ran out. Re-check honestly: a caller must not read a
        # truncated run as a finished one, so the largest remaining edge is in
        # ``after.max_edge`` and ``converged`` stays False unless the target was met.
        converged = all(
            length <= size_field(_midpoint(current.vertices[key[0]], current.vertices[key[1]])) * slack
            for key, length in edge_lengths(current).items()
        )

    after = _stage_metrics(current)
    report = RefinementReport(
        passes=passes,
        converged=converged,
        target_edge_length=None if callable(target) else float(target),
        target_callable=bool(callable(target)),
        target_range=(smallest, largest),
        tolerance=tolerance,
        before=before,
        after=after,
        projection=projection,
        welded_positions=welded,
    )
    return current, report
