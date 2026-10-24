"""Mesh conditioning: shrink-free smoothing and topology-preserving decimation.

Marching tetrahedra is a *volumetric* sampler: it guarantees a closed surface with the
right sign, and says nothing about element shape. On a slender body with knife edges
(a NACA leading edge, a sail's elliptical cap) it produces what it must -- a few
percent of triangles thinner than 2 degrees and aspect ratios in the thousands. Those
are warnings, not errors, and this module is what turns them into a mesh a
finite-element code will accept.

Two operations, in this order, and one rule that governs both:

* :func:`taubin_smooth` -- alternating lambda/mu Laplacian sweeps. The lambda pass
  shrinks, the mu pass (|mu| > lambda) expands, so the pair removes high-frequency
  noise without the volume loss a plain Laplacian would cause.
* :func:`decimate` -- quadric-error-metric edge collapse. Each vertex carries the
  summed plane quadric of its incident faces; collapsing an edge is allowed only when
  it cannot change the topology (the *link condition*: the two one-rings meet in
  exactly the two opposite vertices of the shared edge) and cannot flip or degenerate
  any surviving triangle. Volume is kept because the collapse placement is the
  quadric's minimiser, not the midpoint.

The rule: neither operation is allowed to move the surface off the geometry it came
from. When the caller passes the originating ``GeometryField``, every replacement
position -- each collapse placement and every smoothing sweep -- is projected back onto
the field's zero set with a Newton step. That is what separates "make the triangles
prettier" from "make the shape different": the shape stays the field's, and the
conditioning only re-arranges how it is cut into triangles.

Everything is verified after the fact (:func:`measure_mesh`), because the interesting
failure mode of a decimator is a mesh that looks fine and has a hole in it. The
verification is not a formality: :func:`condition_mesh` returns the measurements before
and after, and refuses to return a mesh that lost watertightness, manifoldness, the
Euler characteristic, or more than a declared share of its volume and area.

Only NumPy is needed (imported lazily, like every other sampling path), and the
topology code is plain Python so a caller can reason about it.
"""
from __future__ import annotations

import heapq
import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

from .geometry import Face, Mesh
from .geometry_field import GeometryField
from .mesh import (
    edge_face_counts,
    is_edge_manifold,
    is_watertight,
    surface_area,
    total_signed_volume,
    triangle_areas,
)
from .validation import Check, measure_quality

__all__ = [
    "ConditionedMesh",
    "MeshConditioningError",
    "condition_mesh",
    "decimate",
    "measure_mesh",
    "project_to_field",
    "taubin_smooth",
]


class MeshConditioningError(ValueError):
    """Raised when conditioning cannot preserve what it promises to preserve."""


# --------------------------------------------------------------------------- #
# Measuring: one record for before and after
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class MeshMeasurements:
    """What a conditioning step is judged on, measured from the mesh itself."""

    vertices: int
    triangles: int
    watertight: bool
    edge_manifold: bool
    signed_volume: float
    surface_area: float
    min_angle_degrees: float
    max_aspect_ratio: float
    sliver_share: float
    zero_area_faces: int

    def to_dict(self) -> dict[str, Any]:
        return {
            "vertices": self.vertices,
            "triangles": self.triangles,
            "watertight": self.watertight,
            "edge_manifold": self.edge_manifold,
            "signed_volume": self.signed_volume,
            "surface_area": self.surface_area,
            "min_angle_degrees": self.min_angle_degrees,
            "max_aspect_ratio": self.max_aspect_ratio,
            "sliver_share": self.sliver_share,
            "zero_area_faces": self.zero_area_faces,
        }

    @property
    def passed(self) -> bool:
        return self.watertight and self.edge_manifold and self.zero_area_faces == 0


def measure_mesh(
    mesh: Mesh,
    *,
    min_angle_degrees: float = 2.0,
    max_aspect_ratio: float = 25.0,
) -> MeshMeasurements:
    """Measure the topology and the element quality of ``mesh`` in one pass.

    The three quality numbers come from :func:`tensorcad.validation.measure_quality`,
    the same function behind :func:`tensorcad.validation.inspect_mesh_quality`, so the
    mesh a caller conditions is judged by exactly the same thresholds as the mesh the
    sampler produced -- there is no second, laxer definition of "good triangle" here.
    """
    quality = measure_quality(
        mesh, min_angle_degrees=min_angle_degrees, max_aspect_ratio=max_aspect_ratio
    )
    zero_area = sum(1 for area in triangle_areas(mesh) if area <= 0.0)
    return MeshMeasurements(
        vertices=len(mesh.vertices),
        triangles=len(mesh.faces),
        watertight=is_watertight(mesh),
        edge_manifold=is_edge_manifold(mesh),
        signed_volume=total_signed_volume(mesh),
        surface_area=surface_area(mesh),
        min_angle_degrees=quality.min_angle_degrees,
        max_aspect_ratio=quality.max_aspect_ratio,
        sliver_share=quality.sliver_share,
        zero_area_faces=zero_area,
    )


# --------------------------------------------------------------------------- #
# Projection: keep the surface on the field it came from
# --------------------------------------------------------------------------- #


def vertex_normals(mesh: Mesh) -> Any:
    """Area-weighted unit vertex normals, as an ``(n, 3)`` array.

    The mesh's own estimate of the outward direction. It is not the field's gradient
    -- but it is free (one pass over the faces, vectorised) and it is within a few
    degrees of the gradient on everything except a knife edge, which is why
    :func:`project_to_field` uses it as the Newton direction.
    """
    numpy = _require_numpy()
    if not mesh.vertices:
        return numpy.zeros((0, 3), dtype=numpy.float64)
    vertices = numpy.asarray(mesh.vertices, dtype=numpy.float64)
    faces = numpy.asarray(mesh.faces, dtype=numpy.int64)
    if not len(faces):
        return numpy.zeros_like(vertices)
    p0, p1, p2 = vertices[faces[:, 0]], vertices[faces[:, 1]], vertices[faces[:, 2]]
    weighted = numpy.cross(p1 - p0, p2 - p0)
    total = numpy.zeros_like(vertices)
    for corner in range(3):
        numpy.add.at(total, faces[:, corner], weighted)
    lengths = numpy.linalg.norm(total, axis=1)
    safe = numpy.where(lengths > 0.0, lengths, 1.0)
    return total / safe[:, None]


def project_to_field(
    mesh: Mesh,
    field: GeometryField,
    *,
    iterations: int = 2,
    max_step: float | None = None,
    method: str = "normal",
) -> tuple[Mesh, float]:
    """Newton-project every vertex onto ``field``'s zero set.

    ``v <- v - f(v) grad f(v) / |grad f(v)|^2`` is the exact step for a signed distance
    (it lands in one iteration) and a contraction for the normalised implicit fields in
    this repository. The step is clamped to ``max_step`` (default: a tenth of the mesh's
    bounding diagonal) so a vertex that started far away cannot jump across the solid.

    Two directions are available, and the difference is a real trade-off rather than a
    detail:

    * ``"normal"`` (default) steps along each vertex's own area-weighted normal,
      ``v <- v - f(v) n_v``. One field evaluation per sweep, no gradient: on a signed
      distance this is exact, on a normalised implicit function it is accurate to the
      angle between ``n_v`` and the true normal, which is second order. Measured on the
      SUBOFF union (70,378 vertices): 0.8 s per sweep, leaving a mean residual of
      7.5 µm and a worst residual of 5.4 mm at the few crease vertices where the mesh
      normal is genuinely ambiguous.
    * ``"gradient"`` steps along ``grad f``, which is exact wherever the field has an
      analytic gradient and *far* more expensive where it does not: this repository's
      appendage fields estimate their gradient by central differences, so this costs
      six field evaluations per vertex per sweep (30 s per sweep on the same mesh,
      against 0.8 s).

    Returns the projected mesh and the largest distance any vertex was moved, which is
    the number to quote when asked how much the projection changed the geometry.
    """
    numpy = _require_numpy()
    batch = int(iterations)
    if batch < 0:
        raise MeshConditioningError(f"iterations must be non-negative, got {iterations!r}")
    if method not in ("normal", "gradient"):
        raise MeshConditioningError(f"method must be 'normal' or 'gradient', got {method!r}")
    if not mesh.vertices:
        return mesh, 0.0
    points = numpy.asarray(mesh.vertices, dtype=numpy.float64)
    box = mesh.bounding_box()
    diagonal = math.sqrt(sum(size * size for size in box.sizes))
    if max_step is None:
        # The step is clamped per vertex to that vertex's own mean edge length. A
        # global clamp (a fraction of the bounding diagonal) is worthless here: on a
        # 4.4 m model it allows a 0.5 m jump, and a handful of crease vertices where
        # the mesh normal is nearly perpendicular to the true one then take it. That
        # is not hypothetical -- it inflated the SUBOFF union's surface area by 3.6 %
        # while leaving every triangle legal, which is exactly the failure the volume
        # and area gates in `condition_mesh` exist to catch.
        local = _local_edge_lengths(mesh)
        limit = numpy.minimum(local, 0.1 * diagonal)
    else:
        limit = float(max_step)
        if limit <= 0.0:
            raise MeshConditioningError(f"max_step must be positive, got {max_step!r}")
        limit = numpy.full(len(points), limit)
    original = points.copy()
    for _ in range(batch):
        values = numpy.asarray(field.signed_distance(points), dtype=numpy.float64)
        if method == "normal":
            direction = vertex_normals(_rebuild(mesh, points))
        else:
            direction = numpy.asarray(field.gradient(points), dtype=numpy.float64)
            denominator = numpy.einsum("ij,ij->i", direction, direction)
            direction = direction / numpy.where(denominator > 1e-30, denominator, 1.0)[:, None]
        step = -values[:, None] * direction
        lengths = numpy.linalg.norm(step, axis=1)
        scale = numpy.where(lengths > limit, limit / numpy.maximum(lengths, 1e-30), 1.0)
        step = step * scale[:, None]
        # Backtracking line search. A Newton step along a *mesh* normal is only
        # approximately the field's steepest-descent direction, and on a thin blade
        # (the SUBOFF fins are 0.1 ft thick at the root) a full step can land on the
        # other side of the blade, where the sign flips and the next sweep pushes it
        # back. Left alone that oscillation random-walks vertices off the surface --
        # measured: 20 sweeps moved individual vertices 120 mm and inflated the area by
        # 3 %. Accepting a step only when it reduces |f| makes the projection monotone
        # and costs at most four extra batched evaluations.
        moved = numpy.abs(step).max(axis=1) > 0.0
        if numpy.any(moved):
            best_points = points.copy()
            best_values = values.copy()
            active = moved
            trial_step = step
            for _ in range(4):
                if not numpy.any(active):
                    break
                trial = points + trial_step * active[:, None]
                trial_values = numpy.asarray(
                    field.signed_distance(trial), dtype=numpy.float64
                )
                better = active & (numpy.abs(trial_values) < numpy.abs(best_values))
                best_points[better] = trial[better]
                best_values[better] = trial_values[better]
                active = active & ~better
                trial_step = trial_step * 0.5
            points = best_points
    moved = numpy.linalg.norm(points - original, axis=1)
    return _rebuild(mesh, points), float(moved.max()) if len(moved) else 0.0


# --------------------------------------------------------------------------- #
# Smoothing
# --------------------------------------------------------------------------- #


def _neighbour_csr(vertex_count: int, faces: Sequence[Face]) -> tuple[Any, Any]:
    """Unique neighbour lists as (indptr, indices), the form NumPy likes."""
    numpy = _require_numpy()
    loops = numpy.asarray(faces, dtype=numpy.int64)
    sources = numpy.concatenate([loops[:, 0], loops[:, 1], loops[:, 2]])
    targets = numpy.concatenate([loops[:, 1], loops[:, 2], loops[:, 0]])
    # Both directions: a triangle contributes u->v and v->u for each of its edges.
    sources, targets = (
        numpy.concatenate([sources, targets]),
        numpy.concatenate([targets, sources]),
    )
    keys = numpy.unique(sources * vertex_count + targets)
    sources = keys // vertex_count
    targets = keys % vertex_count
    counts = numpy.bincount(sources, minlength=vertex_count)
    indptr = numpy.concatenate([[0], numpy.cumsum(counts)]).astype(numpy.int64)
    return indptr, targets.astype(numpy.int64)


def taubin_smooth(
    mesh: Mesh,
    *,
    iterations: int = 10,
    lam: float = 0.5,
    mu: float = -0.53,
    field: GeometryField | None = None,
    projection_iterations: int = 1,
    projection_method: str = "normal",
) -> tuple[Mesh, float]:
    """Taubin lambda/mu smoothing, optionally re-projected onto ``field`` each sweep.

    A single Laplacian pass (lambda) is a low-pass filter that also shrinks the shape by
    an amount proportional to the curvature; the alternating mu pass with ``|mu| >
    lam`` is its inverse at high frequency, so the pair is a band-pass and the volume
    survives. With ``field`` the sweep is followed by a projection, which makes the
    smoothing *tangential*: vertices slide along the surface but cannot leave it.

    The filter is defined for a sampled surface, and it assumes the mesh is dense
    relative to the features it carries. On a very coarse or highly irregular mesh --
    the extreme is a twelve-triangle cube, whose corner modes the mu pass amplifies --
    a sweep can grow the solid instead of relaxing it. There is no internal gate for
    that here; callers that care pass a ``field`` or run the step through
    :func:`condition_mesh`, which refuses to return a mesh whose volume or area moved
    beyond the declared tolerance.

    Returns the smoothed mesh and the largest vertex displacement, so the caller can
    report how much it moved rather than trusting that it did not.
    """
    numpy = _require_numpy()
    if iterations < 0:
        raise MeshConditioningError(f"iterations must be non-negative, got {iterations!r}")
    if not (0.0 < lam < 1.0):
        raise MeshConditioningError(f"lam must be in (0, 1), got {lam!r}")
    if not (-1.0 < mu < 0.0):
        raise MeshConditioningError(f"mu must be in (-1, 0), got {mu!r}")
    if iterations == 0 or not mesh.faces:
        return mesh, 0.0
    points = numpy.asarray(mesh.vertices, dtype=numpy.float64)
    original = points.copy()
    indptr, indices = _neighbour_csr(len(mesh.vertices), mesh.faces)
    counts = numpy.diff(indptr)
    safe_counts = numpy.where(counts > 0, counts, 1)

    def average(values: Any) -> Any:
        sums = numpy.add.reduceat(values[indices], indptr[:-1], axis=0)
        return sums / safe_counts[:, None]

    for _ in range(iterations):
        for factor in (lam, mu):
            points = points + factor * (average(points) - points)
        if field is not None:
            projected, _ = project_to_field(
                _rebuild(mesh, points),
                field,
                iterations=projection_iterations,
                method=projection_method,
            )
            points = numpy.asarray(projected.vertices, dtype=numpy.float64)
    moved = numpy.linalg.norm(points - original, axis=1)
    return _rebuild(mesh, points), float(moved.max()) if len(moved) else 0.0


# --------------------------------------------------------------------------- #
# Decimation: quadric-error edge collapse
# --------------------------------------------------------------------------- #


class _Quadrics:
    """Per-vertex plane quadrics, as the ten unique entries of a symmetric 4x4.

    ``Q(p) = a11 x^2 + a22 y^2 + a33 z^2 + 2 a12 xy + 2 a13 xz + 2 a23 yz
    + 2 a14 x + 2 a24 y + 2 a34 z + a44`` at the homogeneous point ``(x, y, z, 1)``,
    which is the squared distance to the face planes summed over the vertex's
    incident faces and weighted by their areas.

    The entries live in a list of tuples rather than a NumPy array because the
    decimator touches them once per candidate edge -- hundreds of thousands of times
    in pure Python -- and a 3x3 Cramer solve on ten floats is an order of magnitude
    cheaper there than ``numpy.linalg.solve`` on a 4x4.
    """

    def __init__(self, vertices: Any, faces: Any) -> None:
        numpy = _require_numpy()
        p0 = vertices[faces[:, 0]]
        p1 = vertices[faces[:, 1]]
        p2 = vertices[faces[:, 2]]
        normals = numpy.cross(p1 - p0, p2 - p0)
        area2 = numpy.linalg.norm(normals, axis=1)
        good = area2 > 1e-300
        unit = normals[good] / area2[good, None]
        centroids = (p0[good] + p1[good] + p2[good]) / 3.0
        offset = -numpy.einsum("ij,ij->i", unit, centroids)
        weight = 0.5 * area2[good]
        contributions = numpy.stack(
            [
                unit[:, 0] * unit[:, 0],
                unit[:, 0] * unit[:, 1],
                unit[:, 0] * unit[:, 2],
                unit[:, 0] * offset,
                unit[:, 1] * unit[:, 1],
                unit[:, 1] * unit[:, 2],
                unit[:, 1] * offset,
                unit[:, 2] * unit[:, 2],
                unit[:, 2] * offset,
                offset * offset,
            ],
            axis=1,
        ) * weight[:, None]
        total = numpy.zeros((len(vertices), 10), dtype=numpy.float64)
        for corner in range(3):
            numpy.add.at(total, faces[good][:, corner], contributions)
        self.data: list[tuple[float, ...]] = [tuple(row) for row in total]

    def add(self, target: int, source: int) -> None:
        left = self.data[target]
        right = self.data[source]
        self.data[target] = tuple(a + b for a, b in zip(left, right, strict=True))

    @staticmethod
    def evaluate(quadric: Sequence[float], point: Sequence[float]) -> float:
        (a11, a12, a13, a14, a22, a23, a24, a33, a34, a44) = quadric
        x, y, z = point
        return (
            a11 * x * x
            + a22 * y * y
            + a33 * z * z
            + 2.0 * (a12 * x * y + a13 * x * z + a23 * y * z + a14 * x + a24 * y + a34 * z)
            + a44
        )

    def optimal(
        self,
        u: int,
        v: int,
        positions: Sequence[Any],
    ) -> tuple[tuple[float, float, float], float]:
        """The merged quadric's minimiser -- the collapse placement -- and its error."""
        left = self.data[u]
        right = self.data[v]
        merged = tuple(a + b for a, b in zip(left, right, strict=True))
        a11, a12, a13, a14, a22, a23, a24, a33, a34, _a44 = merged
        determinant = (
            a11 * (a22 * a33 - a23 * a23)
            - a12 * (a12 * a33 - a23 * a13)
            + a13 * (a12 * a23 - a22 * a13)
        )
        scale = max(abs(a11), abs(a22), abs(a33), 1e-300)
        if abs(determinant) > 1e-10 * scale**3:
            inv = 1.0 / determinant
            bx, by, bz = -a14, -a24, -a34
            x = inv * (
                bx * (a22 * a33 - a23 * a23)
                - a12 * (by * a33 - a23 * bz)
                + a13 * (by * a23 - a22 * bz)
            )
            y = inv * (
                a11 * (by * a33 - bz * a23)
                - bx * (a12 * a33 - a13 * a23)
                + a13 * (a12 * bz - by * a13)
            )
            z = inv * (
                a11 * (a22 * bz - by * a23)
                - a12 * (a12 * bz - a13 * by)
                + bx * (a12 * a23 - a22 * a13)
            )
        else:
            # A flat or symmetric neighbourhood leaves the normal direction
            # unconstrained; the midpoint is the standard fallback and is what keeps
            # the collapse from extrapolating off a plane.
            x = 0.5 * (positions[u][0] + positions[v][0])
            y = 0.5 * (positions[u][1] + positions[v][1])
            z = 0.5 * (positions[u][2] + positions[v][2])
        point = (float(x), float(y), float(z))
        # Keep the placement inside the edge's own neighbourhood. An ill-conditioned
        # quadric (two nearly parallel planes meeting a third at a shallow angle) can
        # put its minimiser far outside the edge, and a vertex that jumps 80 mm while
        # its neighbours stay put is how a decimator tears a surface even when every
        # triangle it checks is legal. Outside that ball the midpoint is used instead.
        mid = (
            0.5 * (positions[u][0] + positions[v][0]),
            0.5 * (positions[u][1] + positions[v][1]),
            0.5 * (positions[u][2] + positions[v][2]),
        )
        reach = math.dist(positions[u], positions[v])
        if math.dist(point, mid) > reach:
            point = mid
        return point, max(self.evaluate(merged, point), 0.0)


def decimate(
    mesh: Mesh,
    *,
    target_faces: int | None = None,
    max_edge_length: float | None = None,
    field: GeometryField | None = None,
    length_weight: float = 0.0,
    max_aspect_ratio: float = 25.0,
    projection_iterations: int = 1,
) -> tuple[Mesh, dict[str, int]]:
    """Collapse edges under a quadric-error metric, preserving the topology.

    Two stopping rules, both optional and combinable:

    * ``max_edge_length`` -- only edges at most this long may collapse. This is the
      sliver killer: a sliver *is* a short edge between two long ones, so removing short
      edges removes slivers by construction.
    * ``target_faces`` -- stop once the face count reaches it (uniform coarsening).

    An edge is only collapsed when

    * the link condition holds -- ``N(u) & N(v)`` is exactly the opposite vertices of
      the two shared faces -- which is what keeps the surface a manifold and keeps the
      Euler characteristic;
    * every surviving triangle keeps its orientation (strictly positive area, normal
      not reversed), which is what keeps the collapse from folding the surface;
    * the collapse creates no duplicate face;
    * no surviving triangle becomes worse than ``max_aspect_ratio`` unless it was
      already worse (a collapse may repair a bad triangle but may not create one);
    * when ``field`` is given, every vertex of the result is projected onto the field
      once at the end. Projecting each placement *during* the collapse would be more
      precise and much slower -- a single-point field query costs ~4 ms against
      13 microseconds amortised in a batch -- and the difference is below the quadric
      placement's own error.

    Returns the decimated mesh and a tally, including the collapses that were
    *rejected* and why: a decimator that reports only its successes cannot be audited.
    """
    numpy = _require_numpy()
    if target_faces is not None and target_faces < 4:
        raise MeshConditioningError(f"target_faces must be at least 4, got {target_faces!r}")
    if max_edge_length is not None and max_edge_length <= 0.0:
        raise MeshConditioningError(f"max_edge_length must be positive, got {max_edge_length!r}")
    if length_weight < 0.0:
        raise MeshConditioningError(f"length_weight must be non-negative, got {length_weight!r}")
    if target_faces is None and max_edge_length is None:
        raise MeshConditioningError("decimate needs target_faces or max_edge_length")

    if not mesh.faces:
        return mesh, {
            "collapsed": 0,
            "rejected_link": 0,
            "rejected_flip": 0,
            "rejected_quality": 0,
            "projected": 0,
        }
    faces: list[list[int]] = [list(face) for face in mesh.faces]
    alive = [True] * len(faces)
    vertex_count = len(mesh.vertices)
    positions = numpy.asarray(mesh.vertices, dtype=numpy.float64).copy()
    position_list: list[tuple[float, float, float]] = [
        (float(row[0]), float(row[1]), float(row[2])) for row in positions
    ]
    vertex_faces: list[set[int]] = [set() for _ in range(vertex_count)]
    for index, face in enumerate(faces):
        for vertex in face:
            vertex_faces[vertex].add(index)
    versions = [0] * vertex_count
    quadrics = _Quadrics(positions, numpy.asarray(mesh.faces, dtype=numpy.int64))
    stats = {
        "collapsed": 0,
        "rejected_link": 0,
        "rejected_flip": 0,
        "rejected_quality": 0,
        "projected": 0,
    }
    face_count = len(faces)
    if face_count == 0:
        return mesh, stats

    box = mesh.bounding_box()
    diagonal = math.sqrt(sum(size * size for size in box.sizes))
    area_floor = (1e-7 * diagonal) ** 2

    def neighbours(vertex: int) -> set[int]:
        found: set[int] = set()
        for face_index in vertex_faces[vertex]:
            found.update(faces[face_index])
        found.discard(vertex)
        return found

    def candidate(u: int, v: int) -> tuple[float, tuple[float, float, float]] | None:
        """Cost and placement of collapsing ``(u, v)`` into ``u``, or ``None``."""
        shared = vertex_faces[u] & vertex_faces[v]
        if len(shared) != 2:
            return None
        opposite: set[int] = set()
        for face_index in shared:
            opposite.update(vertex for vertex in faces[face_index] if vertex != u and vertex != v)
        if len(opposite) != 2 or (neighbours(u) & neighbours(v)) != opposite:
            return None
        length = math.dist(position_list[u], position_list[v])
        if max_edge_length is not None and length > max_edge_length:
            return None
        point, error = quadrics.optimal(u, v, position_list)
        return error + length_weight * length * length, point

    def legal(u: int, v: int, point: tuple[float, float, float]) -> bool:
        """No flip, no degenerate face, no duplicate face, no worse triangle."""
        shared = vertex_faces[u] & vertex_faces[v]
        survivors = (vertex_faces[u] | vertex_faces[v]) - shared
        created: set[tuple[int, int, int]] = set()
        for face_index in survivors:
            face = faces[face_index]
            old = (position_list[face[0]], position_list[face[1]], position_list[face[2]])
            new: tuple[tuple[float, float, float], ...] = tuple(
                point if index in (u, v) else position_list[index] for index in face
            )
            old_normal = _cross(_sub(old[1], old[0]), _sub(old[2], old[0]))
            new_normal = _cross(_sub(new[1], new[0]), _sub(new[2], new[0]))
            new_area2 = math.sqrt(sum(component * component for component in new_normal))
            if new_area2 <= area_floor:
                stats["rejected_flip"] += 1
                return False
            old_area2 = math.sqrt(sum(component * component for component in old_normal))
            if (
                old_area2 > 0.0
                and sum(a * b for a, b in zip(old_normal, new_normal, strict=True))
                <= 0.3 * old_area2 * new_area2
            ):
                stats["rejected_flip"] += 1
                return False
            new_aspect = _aspect_ratio(new, new_normal, new_area2)
            if new_aspect > max_aspect_ratio and new_aspect > _aspect_ratio(old, old_normal, old_area2):
                # A collapse may repair a bad triangle but may not create one: without
                # this gate the plain quadric minimiser happily trades a sliver for a
                # needle, which is how a decimator makes a mesh *look* better while
                # making it worse.
                stats["rejected_quality"] += 1
                return False
            tag = tuple(sorted(u if index in (u, v) else index for index in face))
            created.add((tag[0], tag[1], tag[2]))
        if len(created) != len(survivors):
            stats["rejected_flip"] += 1
            return False
        return True

    def apply(u: int, v: int, point: tuple[float, float, float]) -> None:
        nonlocal face_count
        shared = vertex_faces[u] & vertex_faces[v]
        for face_index in shared:
            alive[face_index] = False
            for vertex in faces[face_index]:
                vertex_faces[vertex].discard(face_index)
            face_count -= 1
        position_list[u] = point
        positions[u] = point
        for face_index in list(vertex_faces[v]):
            faces[face_index] = [u if index == v else index for index in faces[face_index]]
            vertex_faces[u].add(face_index)
        vertex_faces[v] = set()
        quadrics.add(u, v)
        stats["collapsed"] += 1

    def push(u: int, v: int) -> None:
        result = candidate(u, v)
        if result is not None:
            heapq.heappush(heap, (result[0], u, v, versions[u], versions[v]))

    heap: list[tuple[float, int, int, int, int]] = []
    for u, v in edge_face_counts(mesh):
        if u != v:
            push(min(u, v), max(u, v))
    heapq.heapify(heap)

    while heap:
        if target_faces is not None and face_count <= target_faces:
            break
        _, u, v, version_u, version_v = heapq.heappop(heap)
        if not vertex_faces[u] or not vertex_faces[v]:
            continue
        if version_u != versions[u] or version_v != versions[v]:
            continue
        result = candidate(u, v)
        if result is None:
            stats["rejected_link"] += 1
            continue
        _, point = result
        if not legal(u, v, point):
            continue
        apply(u, v, point)
        around: set[int] = set()
        for face_index in vertex_faces[u]:
            around.update(faces[face_index])
        around.discard(u)
        around.discard(v)
        # Bump every touched vertex *once*, then re-queue: bumping u inside the loop
        # below would invalidate the entries queued moments earlier in the same loop,
        # and the queue would drain long before the target face count.
        versions[u] += 1
        for vertex in around:
            versions[vertex] += 1
        for vertex in around:
            push(min(u, vertex), max(u, vertex))

    keep = [index for index in range(len(faces)) if alive[index]]
    compact: dict[int, int] = {}
    vertices: list[tuple[float, float, float]] = []
    rebuilt: list[tuple[int, int, int]] = []
    for face_index in keep:
        corners = []
        for vertex in faces[face_index]:
            if vertex not in compact:
                compact[vertex] = len(vertices)
                vertices.append(position_list[vertex])
            corners.append(compact[vertex])
        rebuilt.append((corners[0], corners[1], corners[2]))
    reduced = Mesh(
        vertices=tuple(vertices),
        faces=tuple(rebuilt),
        units=mesh.units,
        coordinate_frame=mesh.coordinate_frame,
    )
    if field is not None and stats["collapsed"] and reduced.faces:
        reduced, _ = project_to_field(reduced, field, iterations=projection_iterations)
        stats["projected"] = len(reduced.vertices)
    return reduced, stats


def _local_edge_lengths(mesh: Mesh) -> Any:
    """Mean incident edge length per vertex -- the natural step scale of a mesh."""
    numpy = _require_numpy()
    if not mesh.faces:
        return numpy.zeros(len(mesh.vertices), dtype=numpy.float64)
    vertices = numpy.asarray(mesh.vertices, dtype=numpy.float64)
    faces = numpy.asarray(mesh.faces, dtype=numpy.int64)
    total = numpy.zeros(len(vertices), dtype=numpy.float64)
    counts = numpy.zeros(len(vertices), dtype=numpy.float64)
    for a, b in ((0, 1), (1, 2), (2, 0)):
        first = faces[:, a]
        second = faces[:, b]
        lengths = numpy.linalg.norm(vertices[first] - vertices[second], axis=1)
        numpy.add.at(total, first, lengths)
        numpy.add.at(total, second, lengths)
        numpy.add.at(counts, first, 1.0)
        numpy.add.at(counts, second, 1.0)
    return total / numpy.where(counts > 0.0, counts, 1.0)


def _aspect_ratio(
    triangle: Sequence[Sequence[float]],
    normal: Sequence[float],
    area2: float,
) -> float:
    """Longest edge over twice the height, i.e. the same measure ``validation`` uses."""
    longest = max(
        math.dist(triangle[0], triangle[1]),
        math.dist(triangle[1], triangle[2]),
        math.dist(triangle[2], triangle[0]),
    )
    if area2 <= 0.0:  # pragma: no cover - callers check the area first
        return math.inf
    return longest * longest / area2


def _rebuild(mesh: Mesh, points: Any) -> Mesh:
    numpy = _require_numpy()
    vertices = tuple(
        (float(row[0]), float(row[1]), float(row[2])) for row in numpy.asarray(points)
    )
    return Mesh(
        vertices=vertices,
        faces=mesh.faces,
        units=mesh.units,
        coordinate_frame=mesh.coordinate_frame,
    )


def _require_numpy() -> Any:
    try:
        import numpy
    except ImportError as exc:  # pragma: no cover - numpy is a test/optional dependency
        raise MeshConditioningError(
            "mesh conditioning needs NumPy; install the `test` extra"
        ) from exc
    return numpy


def _sub(a: Sequence[float], b: Sequence[float]) -> tuple[float, float, float]:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _cross(a: Sequence[float], b: Sequence[float]) -> tuple[float, float, float]:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


@dataclass(frozen=True)
class ConditionedMesh:
    """A conditioned mesh plus the measurements that justify handing it on."""

    mesh: Mesh
    before: MeshMeasurements
    after: MeshMeasurements
    collapsed_edges: int
    rejected_link: int
    rejected_flip: int
    max_projection_move: float
    max_smoothing_move: float
    checks: tuple[Check, ...] = ()

    @property
    def passed(self) -> bool:
        return self.before.watertight == self.after.watertight and self.after.passed

    def summary(self) -> dict[str, Any]:
        return {
            "schema": "tensorcad.mesh.conditioning/1",
            "before": self.before.to_dict(),
            "after": self.after.to_dict(),
            "collapsed_edges": self.collapsed_edges,
            "rejected_link": self.rejected_link,
            "rejected_flip": self.rejected_flip,
            "max_projection_move": self.max_projection_move,
            "max_smoothing_move": self.max_smoothing_move,
            "passed": self.passed,
        }


def condition_mesh(
    mesh: Mesh,
    *,
    field: GeometryField | None = None,
    target_faces: int | None = None,
    max_edge_length: float | None = None,
    smooth_iterations: int = 10,
    lam: float = 0.5,
    mu: float = -0.53,
    projection_iterations: int = 2,
    projection_method: str = "normal",
    volume_tolerance: float = 0.005,
    area_tolerance: float = 0.01,
) -> ConditionedMesh:
    """Decimate, then smooth, then check that the shape survived.

    The order matters: collapsing first removes the slivers that make smoothing
    ill-conditioned (a needle's Laplacian is dominated by its short edge), and
    smoothing afterwards repairs the triangles the collapses left behind.

    Refuses to return a mesh that lost watertightness, edge-manifoldness, the Euler
    characteristic, more than ``volume_tolerance`` of its signed volume, or more than
    ``area_tolerance`` of its surface area. Those are not warnings here: an exported
    handoff record that silently changed the shape is worse than no handoff record.
    """
    before = measure_mesh(mesh)
    result = mesh
    collapsed = rejected_link = rejected_flip = 0
    projection_move = 0.0
    smoothing_move = 0.0

    if target_faces is not None or max_edge_length is not None:
        result, stats = decimate(
            result,
            target_faces=target_faces,
            max_edge_length=max_edge_length,
            field=field,
        )
        collapsed = stats["collapsed"]
        rejected_link = stats["rejected_link"]
        rejected_flip = stats["rejected_flip"]

    if smooth_iterations > 0:
        result, smoothing_move = taubin_smooth(
            result,
            iterations=smooth_iterations,
            lam=lam,
            mu=mu,
            field=field,
            projection_method=projection_method,
        )
    if field is not None:
        result, projection_move = project_to_field(
            result, field, iterations=projection_iterations, method=projection_method
        )

    after = measure_mesh(result)
    checks: list[Check] = []
    checks.append(
        Check("conditioned_watertight", after.watertight == before.watertight and after.watertight,
              "" if after.watertight else "the conditioned mesh is not closed", severity="error")
    )
    checks.append(
        Check("conditioned_manifold", after.edge_manifold == before.edge_manifold and after.edge_manifold,
              "" if after.edge_manifold else "the conditioned mesh is not edge-manifold",
              severity="error")
    )
    volume_change = abs(after.signed_volume / before.signed_volume - 1.0) if before.signed_volume else 0.0
    area_change = abs(after.surface_area / before.surface_area - 1.0) if before.surface_area else 0.0
    checks.append(
        Check(
            "conditioned_volume",
            volume_change <= volume_tolerance,
            f"volume changed {100 * volume_change:+.4f} % (limit {100 * volume_tolerance:g} %)",
            severity="error",
        )
    )
    checks.append(
        Check(
            "conditioned_area",
            area_change <= area_tolerance,
            f"area changed {100 * area_change:+.4f} % (limit {100 * area_tolerance:g} %)",
            severity="error",
        )
    )
    conditioned = ConditionedMesh(
        mesh=result,
        before=before,
        after=after,
        collapsed_edges=collapsed,
        rejected_link=rejected_link,
        rejected_flip=rejected_flip,
        max_projection_move=projection_move,
        max_smoothing_move=smoothing_move,
        checks=tuple(checks),
    )
    if not all(check.passed for check in checks):
        failures = "; ".join(check.detail for check in checks if not check.passed)
        raise MeshConditioningError(f"conditioning changed the mesh: {failures}")
    return conditioned
