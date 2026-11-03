"""Curvature and surface fairness.

Every other check in this package asks a *topological* or *metric* question -- is
the surface closed, is it outward oriented, is it the right volume, are its
triangles thin. None of them asks whether the surface is **fair**: whether its
curvature is finite, continuous and free of the accidental kinks that make a hull
expensive to build and noisy in a resistance or added-mass solve. ``hull.py`` says
so in its own docstring ("nothing here checks fairness"), and the mesh conditioning
document ends up reporting slivers along sharp edges without ever measuring the
curvature that produces them.

This module answers the question in two independent ways, because the answer has
to survive the fact that a mesh is a discretisation of something else:

* :func:`principal_curvatures` reads the curvature off the *analytic* field. The
  zero set of a differentiable field ``f`` has unit normal ``n = grad f / |grad f|``
  and its shape operator is the tangential part of ``Hess f / |grad f|``::

      M = P (Hess f) P / |grad f|,        P = I - n n^T

  whose two non-zero eigenvalues are the principal curvatures, positive where the
  surface is convex in the direction of ``grad f`` (a sphere of radius ``R`` is
  ``1/R``, ``1/R``). Because ``P`` kills the normal direction, the answer is
  independent of whether the field is a true signed distance or only an implicit
  function, which is what makes it usable on the ellipsoid and on the learned
  fields as well as on the exact primitives.
* :func:`vertex_mean_curvature` reads it off the *mesh* with the cotangent
  Laplacian, which is the only estimate a caller has once a surface has been
  extracted, exported or conditioned. It converges to the analytic value, and
  :func:`assess_fairness` wraps it together with the per-edge dihedral angles into
  a single report: creases (G0/G1 breaks) are counted instead of being smoothed
  away silently.

Why both: a field is exact where it exists but a solver consumes a mesh, and a mesh
estimate is what a downstream consumer can recompute. If the two disagree the
disagreement is information -- it is the discretisation error of the mesh, and
``examples/fairness`` measures it rather than asserting it away.

Only the standard library is used. Evaluation is a pure Python loop, which is the
right trade at these sizes: curvature is a diagnostic over a few thousand probe
points, not a per-vertex pass over a million-triangle export.
"""
from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import Any

from .geometry import Mesh
from .geometry_field import GeometryField

Point = tuple[float, float, float]

#: Relative probe step used for the finite-difference Hessian, in units of the
#: field's bounding-box diagonal. Central differences of the gradient are ``O(h^2)``
#: accurate and lose ``~eps/|f| / h`` to cancellation, so ``1e-4`` sits in the
#: middle of the usable window by five orders of magnitude on either side.
DEFAULT_STEP_FRACTION = 1e-4

#: Dihedral angle above which an interior edge is classified as a *designed* crease
#: rather than a discretisation artefact. 15 degrees is the usual ship-form
#: convention: a knuckle or a transom edge is far above it, and an under-resolved
#: smooth surface is far below it on any grid this package produces.
DEFAULT_CREASE_DEGREES = 15.0

SCHEMA_CURVATURE = "tensorcad.curvature.sample/1"
SCHEMA_FAIRNESS = "tensorcad.curvature.fairness/1"


class CurvatureError(ValueError):
    """Raised when curvature cannot be defined at the requested point."""


def _checked_points(points: Iterable[Sequence[float]]) -> list[Point]:
    result: list[Point] = []
    for index, point in enumerate(points):
        values = tuple(float(value) for value in point)
        if len(values) != 3:
            raise CurvatureError(f"point {index} has {len(values)} coordinates, expected 3")
        for axis, value in enumerate(values):
            if not math.isfinite(value):
                raise CurvatureError(f"point {index} coordinate {axis} is not finite: {value!r}")
        result.append((values[0], values[1], values[2]))
    return result


def _dot(a: Sequence[float], b: Sequence[float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a: Sequence[float], b: Sequence[float]) -> Point:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _unit(vector: Sequence[float], where: str) -> Point:
    length = math.sqrt(_dot(vector, vector))
    if not math.isfinite(length) or length <= 0.0:
        raise CurvatureError(f"gradient vanished at {where}; the shape operator is undefined there")
    return (vector[0] / length, vector[1] / length, vector[2] / length)


def _tangent_basis(normal: Sequence[float]) -> tuple[Point, Point]:
    """Two orthonormal vectors spanning the plane perpendicular to ``normal``."""
    axis = min(range(3), key=lambda index: abs(normal[index]))
    seed = [0.0, 0.0, 0.0]
    seed[axis] = 1.0
    t1 = _unit(_cross(normal, seed), "tangent basis")
    t2 = _cross(normal, t1)
    return t1, t2


@dataclass(frozen=True)
class CurvatureSample:
    """Principal curvatures of a level set at one point, in ``1/m``.

    ``kappa_max``/``kappa_min`` are ordered by value and signed with respect to the
    field's own normal (``grad f / |grad f|``), so a surface that is convex when
    seen from the outside of the field reports positive curvatures. ``step_error``
    is the change between the two probe steps used, which is the honest error
    indicator for a finite-difference estimate: on the exact primitives it comes
    out at the round-off floor, and on a field with a coarse gradient it is the
    only number that says so.
    """

    point: Point
    normal: Point
    kappa_max: float
    kappa_min: float
    gradient_norm: float
    step: float
    step_error: float

    @property
    def mean_curvature(self) -> float:
        """``H = (kappa_max + kappa_min) / 2``."""
        return 0.5 * (self.kappa_max + self.kappa_min)

    @property
    def gaussian_curvature(self) -> float:
        """``K = kappa_max * kappa_min``."""
        return self.kappa_max * self.kappa_min

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA_CURVATURE,
            "point": list(self.point),
            "normal": list(self.normal),
            "kappa_max_per_m": self.kappa_max,
            "kappa_min_per_m": self.kappa_min,
            "mean_curvature_per_m": self.mean_curvature,
            "gaussian_curvature_per_m2": self.gaussian_curvature,
            "gradient_norm": self.gradient_norm,
            "step": self.step,
            "step_error_per_m": self.step_error,
        }


def _batch_gradients(
    field: GeometryField, probes: Sequence[Point], step: float
) -> list[tuple[Point, list[list[float]]]]:
    """Gradient at each probe and the central-difference Hessian around it.

    The 7 gradients per probe go out in **one** call: fields evaluate a batch
    vectorised (the analytic ones in a tight loop, the learned ones through torch),
    so a per-probe call is dominated by interpreter overhead -- measured at roughly
    an order of magnitude on the AFF-8 field, which is why the batching is here and
    not in the caller.
    """
    batch: list[Point] = []
    for point in probes:
        batch.append(point)
        for axis in range(3):
            plus = list(point)
            minus = list(point)
            plus[axis] += step
            minus[axis] -= step
            batch.append((plus[0], plus[1], plus[2]))
            batch.append((minus[0], minus[1], minus[2]))
    gradients = field.gradient(batch)

    result: list[tuple[Point, list[list[float]]]] = []
    for index, point in enumerate(probes):
        base = 7 * index
        centre = gradients[base]
        hessian = [[0.0, 0.0, 0.0] for _ in range(3)]
        for axis in range(3):
            g_plus = gradients[base + 1 + 2 * axis]
            g_minus = gradients[base + 2 + 2 * axis]
            for row in range(3):
                hessian[row][axis] = (g_plus[row] - g_minus[row]) / (2.0 * step)
        # Symmetrise: the analytic Hessian is symmetric, the difference one is only
        # up to O(h^2), and the 2x2 eigenproblem below assumes symmetry.
        for i in range(3):
            for j in range(i + 1, 3):
                mean = 0.5 * (hessian[i][j] + hessian[j][i])
                hessian[i][j] = mean
                hessian[j][i] = mean
        result.append((tuple(centre), hessian))
    return result


def _curvatures_at(
    gradient: Sequence[float], hessian: Sequence[Sequence[float]], point: Point
) -> tuple[Point, float, float, float]:
    """``(normal, gradient_norm, kappa_max, kappa_min)`` from one gradient + Hessian."""
    normal = _unit(gradient, f"{point!r}")
    length = math.sqrt(_dot(gradient, gradient))
    t1, t2 = _tangent_basis(normal)

    def quadratic(a: Point, b: Point) -> float:
        return sum(a[i] * hessian[i][j] * b[j] for i in range(3) for j in range(3))

    m11 = quadratic(t1, t1) / length
    m12 = quadratic(t1, t2) / length
    m22 = quadratic(t2, t2) / length
    trace = m11 + m22
    determinant = m11 * m22 - m12 * m12
    discriminant = max(0.0, 0.25 * trace * trace - determinant)
    root = math.sqrt(discriminant)
    return normal, length, 0.5 * trace + root, 0.5 * trace - root


def principal_curvatures(
    field: GeometryField,
    points: Iterable[Sequence[float]],
    *,
    step: float | None = None,
    refine: bool = True,
) -> list[CurvatureSample]:
    """Principal curvatures of ``field``'s zero set at each probe point.

    The probe points are expected to lie *on* the surface. A point off the surface
    still returns the curvature of the level set that passes through it, which is
    the right answer for a signed distance and a meaningless one for an implicit
    function whose level sets do not resemble the shape -- so callers that need the
    surface should project first (``to_mesh`` or ``project_to_field`` does that).

    ``step`` defaults to :data:`DEFAULT_STEP_FRACTION` of the field's bounding-box
    diagonal. With ``refine`` the probe is repeated at half the step and the two
    results are Richardson-extrapolated (the difference estimator is second order),
    with ``step_error`` reporting how far the coarse probe had to move.
    """
    probes = _checked_points(points)
    if step is None:
        bounds = field.bounds()
        diagonal = 1.0 if bounds is None else _dot(
            [bounds.maximum[i] - bounds.minimum[i] for i in range(3)],
            [bounds.maximum[i] - bounds.minimum[i] for i in range(3)],
        ) ** 0.5
        step = DEFAULT_STEP_FRACTION * diagonal if diagonal > 0.0 else DEFAULT_STEP_FRACTION
    if not (math.isfinite(step) and step > 0.0):
        raise CurvatureError(f"step must be positive and finite, got {step!r}")

    coarse = _batch_gradients(field, probes, step)
    fine = _batch_gradients(field, probes, 0.5 * step) if refine else []

    samples: list[CurvatureSample] = []
    for index, point in enumerate(probes):
        gradient, hessian = coarse[index]
        normal, length, kappa_max, kappa_min = _curvatures_at(gradient, hessian, point)
        if not refine:
            samples.append(
                CurvatureSample(point, normal, kappa_max, kappa_min, length, step, float("nan"))
            )
            continue
        _, fine_hessian = fine[index]
        _, _, fine_max, fine_min = _curvatures_at(gradient, fine_hessian, point)
        # Richardson for a second-order estimator: 4/3 of the fine value minus 1/3
        # of the coarse one removes the leading O(h^2) term.
        extrapolated_max = (4.0 * fine_max - kappa_max) / 3.0
        extrapolated_min = (4.0 * fine_min - kappa_min) / 3.0
        error = max(
            abs(kappa_max - fine_max),
            abs(kappa_min - fine_min),
            abs(extrapolated_max - fine_max),
            abs(extrapolated_min - fine_min),
        )
        samples.append(
            CurvatureSample(
                point, normal, extrapolated_max, extrapolated_min, length, 0.5 * step, error
            )
        )
    return samples


def vertex_mean_curvature(
    mesh: Mesh,
    *,
    clamp_cotangents: bool = True,
) -> list[float]:
    """Mean curvature per vertex from the cotangent Laplacian, in ``1/m``.

    ``K_i = (1 / (2 A_i)) * sum_j (cot a_ij + cot b_ij) (v_i - v_j)`` over the
    one-ring, with ``A_i`` the *mixed* (Voronoi) area and the signed mean curvature
    read off as ``H_i = 0.5 * (K_i . n_i)`` for outward vertex normals -- the sign
    convention is the field one, so a convex surface is positive.

    Obtuse triangles make a cotangent negative and, worse, can make their sum zero,
    which blows the estimate up on exactly the sliver field this package complains
    about elsewhere. ``clamp_cotangents`` clamps each cotangent from below, the
    standard guard: it costs accuracy on a well-shaped mesh (nothing to clamp) and
    turns a singular value into a finite one on a bad one.

    The estimate is second-order accurate on a regular triangulation, and
    ``examples/fairness/run_example.py`` measures its convergence on a sphere
    instead of assuming it.
    """
    vertices = mesh.vertices
    if not mesh.faces:
        raise CurvatureError("mesh has no faces")
    cotangent_floor = -1e3 if clamp_cotangents else -float("inf")

    # edge -> (angles opposite the edge) for each incident face
    edge_angles: dict[tuple[int, int], list[float]] = defaultdict(list)
    vertex_angle_sum = [0.0] * len(vertices)
    vertex_area_sum = [0.0] * len(vertices)
    for face in mesh.faces:
        a, b, c = face
        for i, j, k in ((a, b, c), (b, c, a), (c, a, b)):
            # angle at k, opposite the edge (i, j)
            u = [vertices[i][axis] - vertices[k][axis] for axis in range(3)]
            v = [vertices[j][axis] - vertices[k][axis] for axis in range(3)]
            cross = _cross(u, v)
            sin = math.sqrt(_dot(cross, cross))
            cos = _dot(u, v)
            angle = math.atan2(sin, cos)
            vertex_angle_sum[k] += angle
            edge_angles[(i, j) if i <= j else (j, i)].append(angle)
    for face in mesh.faces:
        a, b, c = face
        u = [vertices[b][axis] - vertices[a][axis] for axis in range(3)]
        v = [vertices[c][axis] - vertices[a][axis] for axis in range(3)]
        cross = _cross(u, v)
        area = 0.5 * math.sqrt(_dot(cross, cross))
        for index in (a, b, c):
            vertex_area_sum[index] += area / 3.0

    laplacian = [[0.0, 0.0, 0.0] for _ in vertices]
    for (i, j), angles in edge_angles.items():
        weight = sum(1.0 / math.tan(angle) if abs(math.tan(angle)) > 1e-12 else 0.0 for angle in angles)
        weight = max(weight, cotangent_floor)
        delta = [vertices[i][axis] - vertices[j][axis] for axis in range(3)]
        for axis in range(3):
            laplacian[i][axis] += weight * delta[axis]
            laplacian[j][axis] -= weight * delta[axis]

    normals = _vertex_normals(mesh)
    result: list[float] = []
    for index, area in enumerate(vertex_area_sum):
        if area <= 0.0:
            result.append(float("nan"))
            continue
        vector = [laplacian[index][axis] / (2.0 * area) for axis in range(3)]
        result.append(0.5 * _dot(vector, normals[index]))
    return result


def _vertex_normals(mesh: Mesh) -> list[Point]:
    """Area-weighted vertex normals, normalised (outward for an oriented mesh)."""
    accum = [[0.0, 0.0, 0.0] for _ in mesh.vertices]
    for face in mesh.faces:
        a, b, c = face
        u = [mesh.vertices[b][axis] - mesh.vertices[a][axis] for axis in range(3)]
        v = [mesh.vertices[c][axis] - mesh.vertices[a][axis] for axis in range(3)]
        normal = _cross(u, v)
        for index in face:
            for axis in range(3):
                accum[index][axis] += normal[axis]
    normals: list[Point] = []
    for vector in accum:
        length = math.sqrt(_dot(vector, vector))
        if length <= 0.0:
            normals.append((0.0, 0.0, 0.0))
        else:
            normals.append((vector[0] / length, vector[1] / length, vector[2] / length))
    return normals


def _percentile(sorted_values: Sequence[float], fraction: float) -> float:
    if not sorted_values:
        return float("nan")
    if len(sorted_values) == 1:
        return float(sorted_values[0])
    position = fraction * (len(sorted_values) - 1)
    lower = int(math.floor(position))
    upper = min(lower + 1, len(sorted_values) - 1)
    weight = position - lower
    return float(sorted_values[lower] * (1.0 - weight) + sorted_values[upper] * weight)


def _distribution(values: Sequence[float]) -> dict[str, float]:
    finite = sorted(value for value in values if math.isfinite(value))
    if not finite:
        return {"count": 0, "mean": float("nan"), "p50": float("nan"),
                "p95": float("nan"), "max": float("nan")}
    return {
        "count": len(finite),
        "mean": math.fsum(finite) / len(finite),
        "p50": _percentile(finite, 0.5),
        "p95": _percentile(finite, 0.95),
        "max": float(finite[-1]),
    }


@dataclass(frozen=True)
class FairnessReport:
    """How smooth a surface is, in the two ways that matter to a consumer.

    *Creases* are the G1 story: the distribution of dihedral angles over interior
    edges, and how many of them exceed the crease threshold. A designed knuckle and
    an under-resolved smooth patch look identical edge by edge, so the report keeps
    the distribution rather than only the count.
    *Curvature jumps* are the G2 story, and they need a field: with one, the mean
    curvature is evaluated at the two face centroids either side of an interior edge
    and their difference is the jump the geometry actually has, independent of how
    the mesh was triangulated.

    ``creases`` is the check a modeller wants: no interior edge above the threshold
    that is not meant to be there. It is reported, not enforced -- the AFF-8
    appendage set has genuine 0-degree knife edges, so a gate on it would reject a
    correct model.

    The curvature summary is reported twice on purpose. ``mesh_curvature_per_m``
    comes from :func:`vertex_mean_curvature` and is what a consumer of the mesh can
    recompute; its ``area_weighted_mean`` is the integral ``(1/A) * int H dA`` and
    its median/p95 describe the pointwise field, which on a raw marching-tet mesh is
    dominated by slivers. ``field_curvature_per_m`` is the analytic answer at the
    face centroids and is the reference the mesh number has to be judged against.
    """

    triangles: int
    vertices: int
    interior_edges: int
    boundary_edges: int
    watertight: bool
    edge_manifold: bool
    crease_angle_degrees: float
    dihedral_degrees: dict[str, float]
    crease_edges: int
    crease_length_fraction: float
    mesh_curvature_per_m: dict[str, float]
    field_curvature_per_m: dict[str, float] | None
    curvature_jump_per_m: dict[str, float] | None

    @property
    def minimum_angle_degrees(self) -> float:
        """Smallest dihedral angle seen; ``0`` means the surface has a knife edge."""
        return float(self.dihedral_degrees["min"]) if "min" in self.dihedral_degrees else float("nan")

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA_FAIRNESS,
            "triangles": self.triangles,
            "vertices": self.vertices,
            "interior_edges": self.interior_edges,
            "boundary_edges": self.boundary_edges,
            "watertight": self.watertight,
            "edge_manifold": self.edge_manifold,
            "crease_angle_degrees": self.crease_angle_degrees,
            "dihedral_degrees": dict(self.dihedral_degrees),
            "crease_edges": self.crease_edges,
            "crease_length_fraction": self.crease_length_fraction,
            "mesh_curvature_per_m": dict(self.mesh_curvature_per_m),
            "field_curvature_per_m": (
                None if self.field_curvature_per_m is None else dict(self.field_curvature_per_m)
            ),
            "curvature_jump_per_m": (
                None if self.curvature_jump_per_m is None else dict(self.curvature_jump_per_m)
            ),
        }


def assess_fairness(
    mesh: Mesh,
    *,
    field: GeometryField | None = None,
    crease_angle_degrees: float = DEFAULT_CREASE_DEGREES,
) -> FairnessReport:
    """Dihedral-angle and curvature-jump statistics for ``mesh``.

    Supply ``field`` when the mesh was extracted from one: the curvature jump across
    an edge is then a property of the geometry rather than of the triangulation.
    Without it the report is dihedral-only, which is still enough to tell a sharp
    edge from a smooth patch.
    """
    if not (math.isfinite(crease_angle_degrees) and 0.0 < crease_angle_degrees < 180.0):
        raise CurvatureError(
            f"crease_angle_degrees must lie in (0, 180), got {crease_angle_degrees!r}"
        )
    if not mesh.faces:
        raise CurvatureError("mesh has no faces")

    edge_faces: dict[tuple[int, int], list[int]] = defaultdict(list)
    for index, face in enumerate(mesh.faces):
        for i, j in ((face[0], face[1]), (face[1], face[2]), (face[2], face[0])):
            edge_faces[(i, j) if i <= j else (j, i)].append(index)

    normals: list[Point] = []
    face_areas: list[float] = []
    for index, face in enumerate(mesh.faces):
        u = [mesh.vertices[face[1]][axis] - mesh.vertices[face[0]][axis] for axis in range(3)]
        v = [mesh.vertices[face[2]][axis] - mesh.vertices[face[0]][axis] for axis in range(3)]
        cross = _cross(u, v)
        normals.append(_unit(cross, f"face {index}"))  # raises on a degenerate face
        face_areas.append(0.5 * math.sqrt(_dot(cross, cross)))

    dihedrals: list[float] = []
    crease_edges = 0
    crease_length = 0.0
    interior_length = 0.0
    interior_edges = 0
    interior_pairs: list[tuple[int, int]] = []
    for edge, faces in edge_faces.items():
        if len(faces) != 2:
            continue
        interior_edges += 1
        edge_vector = [mesh.vertices[edge[1]][axis] - mesh.vertices[edge[0]][axis] for axis in range(3)]
        edge_length = math.sqrt(_dot(edge_vector, edge_vector))
        interior_length += edge_length
        cosine = max(-1.0, min(1.0, _dot(normals[faces[0]], normals[faces[1]])))
        angle = math.degrees(math.acos(cosine))
        dihedrals.append(angle)
        if angle > crease_angle_degrees:
            crease_edges += 1
            crease_length += edge_length
        interior_pairs.append((faces[0], faces[1]))

    distribution = _distribution(dihedrals)
    distribution["min"] = float(min(dihedrals)) if dihedrals else float("nan")

    mesh_curvature = vertex_mean_curvature(mesh)
    mesh_weights = [0.0] * len(mesh.vertices)
    for face, area in zip(mesh.faces, face_areas):
        for index in face:
            mesh_weights[index] += area / 3.0
    mesh_curvature_summary = _distribution([abs(value) for value in mesh_curvature])
    total_weight = math.fsum(mesh_weights)
    weighted_pairs = [
        (value, weight)
        for value, weight in zip(mesh_curvature, mesh_weights)
        if math.isfinite(value) and weight > 0.0
    ]
    mesh_curvature_summary["area_weighted_mean"] = (
        math.fsum(value * weight for value, weight in weighted_pairs)
        / math.fsum(weight for _, weight in weighted_pairs)
        if weighted_pairs and total_weight > 0.0
        else float("nan")
    )

    curvature_summary: dict[str, float] | None = None
    jump_summary: dict[str, float] | None = None
    if field is not None:
        centroids = [
            tuple(
                (mesh.vertices[face[0]][axis] + mesh.vertices[face[1]][axis] + mesh.vertices[face[2]][axis]) / 3.0
                for axis in range(3)
            )
            for face in mesh.faces
        ]
        samples = principal_curvatures(field, centroids, refine=False)
        mean_curvatures = [sample.mean_curvature for sample in samples]
        curvature_summary = _distribution([abs(value) for value in mean_curvatures])
        jumps = [
            abs(mean_curvatures[first] - mean_curvatures[second])
            for first, second in interior_pairs
        ]
        jump_summary = _distribution(jumps)

    boundary = sum(1 for faces in edge_faces.values() if len(faces) == 1)
    return FairnessReport(
        triangles=len(mesh.faces),
        vertices=len(mesh.vertices),
        interior_edges=interior_edges,
        boundary_edges=boundary,
        watertight=boundary == 0 and interior_edges > 0,
        edge_manifold=all(len(faces) <= 2 for faces in edge_faces.values()),
        crease_angle_degrees=crease_angle_degrees,
        dihedral_degrees=distribution,
        crease_edges=crease_edges,
        crease_length_fraction=(crease_length / interior_length) if interior_length > 0.0 else float("nan"),
        mesh_curvature_per_m=mesh_curvature_summary,
        field_curvature_per_m=curvature_summary,
        curvature_jump_per_m=jump_summary,
    )


__all__ = [
    "CurvatureError",
    "CurvatureSample",
    "DEFAULT_CREASE_DEGREES",
    "DEFAULT_STEP_FRACTION",
    "FairnessReport",
    "assess_fairness",
    "principal_curvatures",
    "vertex_mean_curvature",
]