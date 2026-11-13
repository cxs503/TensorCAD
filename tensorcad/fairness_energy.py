"""Curvature *fields* and smoothing energy: what a surface's curvature costs.

``curvature.py`` answers "what is the curvature here?" for one probe point and
summarises it as a distribution; ``geometry_audit.py`` turns the discrete Gaussian
curvature into an audit. Neither answers the question a hull designer actually
asks: **how much curvature does this shape carry in total, and is that number
stable?** That is the smoothing (fairness) energy, and it is what this module
adds:

    E_bending = int_S (kappa1**2 + kappa2**2) dA      thin-plate / fairness energy
    E_willmore = int_S H**2 dA                        Willmore integrand
    E_gauss = int_S K dA                              Gauss-Bonnet charge

All three are dimensionless -- ``dA`` carries ``1/m**2`` -- which is what makes
them comparable between a 10 m model and a 300 m ship, and it is the first claim
this module is tested against (uniform scaling must leave the number alone).

They are not independent. With ``H = (k1 + k2) / 2`` and ``K = k1 k2``::

    k1**2 + k2**2 = 4 H**2 - 2 K

so ``E_bending = 4 E_willmore - 2 E_gauss`` cell by cell. On the field path the
identity is a round-off check of the eigensolver; on the mesh path it is the
*definition* (the mesh has H and K, not k1 and k2), and the residual is reported
rather than hidden.

Two reference anchors make the numbers falsifiable without a third-party tool:

* a closed genus-0 surface has ``E_gauss = 4*pi`` exactly (Gauss-Bonnet), and the
  discrete angle-defect estimator reproduces that to round-off by construction, so
  it validates the quadrature weights rather than the geometry;
* a sphere of radius ``R`` has ``E_bending = 8*pi`` and ``E_willmore = 4*pi`` for
  every ``R`` (Willmore's inequality is saturated only by the round sphere).

For a hull the interesting part is *where the energy sits*. A designed knuckle
(keel line, transom edge, deck edge) has an infinite curvature in the continuum
limit, so a pointwise mesh estimate of it grows like ``1/h`` and the total energy
diverges under refinement -- which is a property of the shape, not a bug in the
estimator. ``crease_angle_degrees`` therefore splits the surface: cells touching a
crease or an open boundary are summed separately (:class:`SmoothingEnergy.crease_energy`)
from the smooth interior, and only the latter is expected to converge. Reporting
the split ratio is the diagnosis: it says how much of the surface's curvature cost
is a design decision and how much is a fair surface.

Numerical stability is a design constraint, not an afterthought:

* every accumulation goes through :func:`math.fsum`, which is exactly rounded, so
  the result is bit-for-bit independent of face order (tested by shuffling);
* the ``condition`` field is ``sum|integrand*dA| / |sum integrand*dA|``: 1 for the
  positive denstities, and as large as ``sum|K|dA / 4*pi`` for the signed charge,
  which is the honest measure of how much cancellation the sum survived;
* non-positive cells are counted (``nonpositive_fraction``) instead of clipped
  silently, because ``4 H**2 - 2 K < 0`` is a discretisation artefact and hiding it
  would make the energy look smoother than the mesh is.

Standard library only, like the rest of the geometry stack.
"""
from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Any

from .curvature import CurvatureError, principal_curvatures, vertex_mean_curvature
from .geometry import Mesh
from .geometry_audit import gaussian_curvature
from .geometry_field import GeometryField

SCHEMA_FIELD = "tensorcad.fairness_energy.field/1"
SCHEMA_ENERGY = "tensorcad.fairness_energy.energy/1"
SCHEMA_REPORT = "tensorcad.fairness_energy.report/1"

#: The three integrands, keyed by the name the API takes. ``bending`` is the
#: thin-plate (``k1**2 + k2**2``) energy, ``willmore`` is ``H**2``, and ``gaussian``
#: is the signed Gauss-Bonnet charge ``K``.
ENERGY_KINDS: tuple[str, ...] = ("bending", "willmore", "gaussian")


def _mean(values: Sequence[float]) -> float:
    finite = [value for value in values if math.isfinite(value)]
    return math.fsum(finite) / len(finite) if finite else float("nan")


def _integrand(kind: str, kappa_max: float, kappa_min: float, h: float, k: float) -> float:
    if kind == "bending":
        return kappa_max * kappa_max + kappa_min * kappa_min
    if kind == "willmore":
        return h * h
    if kind == "gaussian":
        return k
    raise CurvatureError(f"unknown energy kind {kind!r}; expected one of {ENERGY_KINDS}")


# --------------------------------------------------------------------------- #
# curvature field
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class MeshCurvatureField:
    """Mean, Gaussian and principal curvature at every vertex of a mesh.

    ``kappa_max``/``kappa_min`` come from inverting ``H`` and ``K``::

        k1,2 = H +- sqrt(H**2 - K)

    The discriminant is non-negative in the exact geometry and can go negative on a
    mesh, because the cotangent Laplacian and the angle-defect estimator are two
    different discretisations. ``negative_discriminant_vertices`` counts those
    vertices and their curvatures are reported as ``nan`` (never clipped to zero:
    a clipped vertex would look perfectly umbilic). Boundary vertices are marked
    ``interior=False`` for the same reason -- a one-ring that is not closed cannot
    define a curvature.
    """

    mean: tuple[float, ...]
    gaussian: tuple[float, ...]
    kappa_max: tuple[float, ...]
    kappa_min: tuple[float, ...]
    area: tuple[float, ...]
    interior: tuple[bool, ...]
    negative_discriminant_vertices: int

    @property
    def curvedness(self) -> tuple[float, ...]:
        """``sqrt((k1**2 + k2**2) / 2)`` per vertex; ``nan`` where undefined."""
        return tuple(
            math.sqrt(0.5 * (first * first + second * second))
            if math.isfinite(first) and math.isfinite(second)
            else float("nan")
            for first, second in zip(self.kappa_max, self.kappa_min, strict=True)
        )

    def to_dict(self) -> dict[str, Any]:
        finite_h = [value for value in self.mean if math.isfinite(value)]
        finite_k = [value for value in self.gaussian if math.isfinite(value)]
        positive = [value for value in self.gaussian if math.isfinite(value) and value > 0.0]
        return {
            "schema": SCHEMA_FIELD,
            "vertex_count": len(self.mean),
            "interior_vertices": sum(1 for flag in self.interior if flag),
            "negative_discriminant_vertices": self.negative_discriminant_vertices,
            "mean_curvature_per_m": {
                "mean": _mean(finite_h),
                "abs_mean": _mean([abs(value) for value in finite_h]),
                "max_abs": max((abs(value) for value in finite_h), default=float("nan")),
            },
            "gaussian_curvature_per_m2": {
                "mean": _mean(finite_k),
                "abs_mean": _mean([abs(value) for value in finite_k]),
                "elliptic_fraction": (
                    len(positive) / len(finite_k) if finite_k else float("nan")
                ),
            },
        }


def principal_curvature_field(mesh: Mesh) -> MeshCurvatureField:
    """Mean/Gaussian/principal curvature at every vertex of ``mesh``.

    ``H`` from the cotangent Laplacian (:func:`tensorcad.curvature.vertex_mean_curvature`)
    and ``K`` from the angle deficit over the mixed area
    (:func:`tensorcad.geometry_audit.gaussian_curvature`) are two *independent*
    estimators -- a divergence of a gradient and a holonomy -- and the principal
    curvatures are the pair that reproduces both. Where they disagree enough to make
    ``H**2 - K`` negative the vertex is reported as ``nan`` and counted.
    """
    if not mesh.faces:
        raise CurvatureError("mesh has no faces")
    mean = vertex_mean_curvature(mesh)
    discrete = gaussian_curvature(mesh)
    kappa_max: list[float] = []
    kappa_min: list[float] = []
    negative = 0
    for h, k, flag in zip(mean, discrete.curvature, discrete.interior, strict=True):
        if not flag or not (math.isfinite(h) and math.isfinite(k)):
            kappa_max.append(float("nan"))
            kappa_min.append(float("nan"))
            continue
        discriminant = h * h - k
        if discriminant < 0.0:
            negative += 1
            kappa_max.append(float("nan"))
            kappa_min.append(float("nan"))
            continue
        root = math.sqrt(discriminant)
        kappa_max.append(h + root)
        kappa_min.append(h - root)
    return MeshCurvatureField(
        mean=tuple(mean),
        gaussian=tuple(discrete.curvature),
        kappa_max=tuple(kappa_max),
        kappa_min=tuple(kappa_min),
        area=tuple(discrete.mixed_area),
        interior=tuple(discrete.interior),
        negative_discriminant_vertices=negative,
    )


# --------------------------------------------------------------------------- #
# crease / boundary classification
# --------------------------------------------------------------------------- #
def _face_edges(face: Sequence[int]) -> tuple[tuple[int, int], ...]:
    return (
        (face[0], face[1]) if face[0] <= face[1] else (face[1], face[0]),
        (face[1], face[2]) if face[1] <= face[2] else (face[2], face[1]),
        (face[2], face[0]) if face[2] <= face[0] else (face[0], face[2]),
    )


def _face_area(mesh: Mesh, face: Sequence[int]) -> float:
    a, b, c = (mesh.vertices[index] for index in face)
    u = [b[axis] - a[axis] for axis in range(3)]
    v = [c[axis] - a[axis] for axis in range(3)]
    cross = (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    )
    return 0.5 * math.sqrt(cross[0] ** 2 + cross[1] ** 2 + cross[2] ** 2)


def _face_centroid(mesh: Mesh, face: Sequence[int]) -> tuple[float, float, float]:
    a, b, c = (mesh.vertices[index] for index in face)
    return (
        (a[0] + b[0] + c[0]) / 3.0,
        (a[1] + b[1] + c[1]) / 3.0,
        (a[2] + b[2] + c[2]) / 3.0,
    )


def _crease_and_boundary_faces(
    mesh: Mesh, crease_angle_degrees: float | None
) -> tuple[set[int], set[int]]:
    """Faces touching an open boundary, and faces touching a crease above threshold.

    "Touching" is taken **at vertex level, one ring wide**. Excluding only the faces
    that carry the crease edge is not enough: the cotangent estimate at the vertices
    on either side of a crease is ``O(1/h)`` in its own right, so a face that merely
    *shares a vertex* with a crease still inherits the blow-up. The filter therefore
    flags the vertices of every crease/boundary edge and excludes every face that
    uses one of them, which turns a ``1/h`` divergence into a converged sum (the
    excluded band has ``O(h)`` width, so its area share vanishes under refinement).
    """
    boundary_faces: set[int] = set()
    crease_faces: set[int] = set()
    edge_faces: dict[tuple[int, int], list[int]] = defaultdict(list)
    for index, face in enumerate(mesh.faces):
        for edge in _face_edges(face):
            edge_faces[edge].append(index)
    normals: dict[int, tuple[float, float, float]] = {}
    if crease_angle_degrees is not None:
        for index, face in enumerate(mesh.faces):
            a, b, c = (mesh.vertices[vertex] for vertex in face)
            u = [b[axis] - a[axis] for axis in range(3)]
            v = [c[axis] - a[axis] for axis in range(3)]
            cross = (
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            )
            length = math.sqrt(cross[0] ** 2 + cross[1] ** 2 + cross[2] ** 2)
            normals[index] = (
                (cross[0] / length, cross[1] / length, cross[2] / length)
                if length > 0.0
                else (0.0, 0.0, 0.0)
            )
    sharp_vertices: set[int] = set()
    for edge, faces in edge_faces.items():
        if len(faces) == 1:
            sharp_vertices.update(edge)
            continue
        if len(faces) != 2 or crease_angle_degrees is None:
            continue
        first, second = normals[faces[0]], normals[faces[1]]
        cosine = max(-1.0, min(1.0, sum(a * b for a, b in zip(first, second, strict=True))))
        if math.degrees(math.acos(cosine)) > crease_angle_degrees:
            sharp_vertices.update(edge)
    for index, face in enumerate(mesh.faces):
        if any(vertex in sharp_vertices for vertex in face):
            crease_faces.add(index)
    boundary_faces = {index for index, face in enumerate(mesh.faces)
                      if any(vertex in sharp_vertices for vertex in face)}
    boundary_faces = set(boundary_faces & crease_faces)
    return crease_faces, boundary_faces


# --------------------------------------------------------------------------- #
# smoothing energy
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class SmoothingEnergy:
    """``int_S integrand dA`` for one curvature density, with its conditioning.

    ``smooth_energy`` sums cells that touch neither a crease above the requested
    threshold nor an open boundary; ``crease_energy`` sums the cells that do, and
    ``total_energy`` is the two together -- reported as three separate numbers so a
    divergent crease contribution cannot hide inside a converged one.
    """

    kind: str
    source: str
    smooth_energy: float
    crease_energy: float
    total_energy: float
    smooth_area: float
    crease_area: float
    cells: int
    smooth_cells: int
    integrand: dict[str, float]
    nonpositive_fraction: float
    condition: float
    crease_angle_degrees: float | None
    identity_residual: float | None = None
    summation: str = "math.fsum"
    notes: tuple[str, ...] = field(default_factory=tuple)

    @property
    def area(self) -> float:
        return self.smooth_area + self.crease_area

    @property
    def smooth_per_unit_area(self) -> float:
        return self.smooth_energy / self.smooth_area if self.smooth_area > 0.0 else float("nan")

    @property
    def crease_area_fraction(self) -> float:
        return self.crease_area / self.area if self.area > 0.0 else float("nan")

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA_ENERGY,
            "kind": self.kind,
            "source": self.source,
            "smooth_energy": self.smooth_energy,
            "crease_energy": self.crease_energy,
            "total_energy": self.total_energy,
            "smooth_area": self.smooth_area,
            "crease_area": self.crease_area,
            "smooth_per_unit_area": self.smooth_per_unit_area,
            "crease_area_fraction": self.crease_area_fraction,
            "cells": self.cells,
            "smooth_cells": self.smooth_cells,
            "integrand": dict(self.integrand),
            "nonpositive_fraction": self.nonpositive_fraction,
            "condition": self.condition,
            "crease_angle_degrees": self.crease_angle_degrees,
            "identity_residual": self.identity_residual,
            "summation": self.summation,
            "notes": list(self.notes),
        }


def smoother(
    *,
    kind: str,
    source: str,
    terms: Sequence[tuple[float, float, bool]],
    crease_angle_degrees: float | None,
    integrand_values: Sequence[float],
    identity_residual: float | None,
    notes: Sequence[str] = (),
) -> SmoothingEnergy:
    """Assemble a :class:`SmoothingEnergy` from ``(integrand, area, is_smooth)`` terms.

    Kept separate from the two evaluators so the accumulation policy lives in one
    place: exactly rounded sums, three separate buckets (smooth, crease, total) and
    no clipping.
    """
    smooth = [(value, area) for value, area, flag in terms if flag]
    crease = [(value, area) for value, area, flag in terms if not flag]
    smooth_energy = math.fsum(value * area for value, area in smooth)
    crease_energy = math.fsum(value * area for value, area in crease)
    total_energy = math.fsum((smooth_energy, crease_energy))
    absolute = math.fsum(abs(value * area) for value, area, _ in terms)
    finite_values = [value for value in integrand_values if math.isfinite(value)]
    negatives = sum(1 for value in finite_values if value <= 0.0)
    return SmoothingEnergy(
        kind=kind,
        source=source,
        smooth_energy=smooth_energy,
        crease_energy=crease_energy,
        total_energy=total_energy,
        smooth_area=math.fsum(area for _, area in smooth),
        crease_area=math.fsum(area for _, area in crease),
        cells=len(terms),
        smooth_cells=len(smooth),
        integrand={
            "count": len(finite_values),
            "mean": _mean(finite_values),
            "min": min(finite_values, default=float("nan")),
            "max": max(finite_values, default=float("nan")),
        },
        nonpositive_fraction=(negatives / len(finite_values)) if finite_values else float("nan"),
        condition=(absolute / abs(total_energy)) if total_energy != 0.0 else float("inf"),
        crease_angle_degrees=crease_angle_degrees,
        identity_residual=identity_residual,
        notes=tuple(notes),
    )


def smoothing_energy(
    mesh: Mesh,
    *,
    kind: str = "bending",
    field: GeometryField | None = None,
    crease_angle_degrees: float | None = None,
) -> SmoothingEnergy:
    """Smoothing energy of ``mesh``, read off ``field`` when it is supplied.

    * ``field`` given: the principal curvatures are exact at each face centroid
      (``refine=False``, because Richardson at a centroid buys round-off and costs
      7 field evaluations), and the energy is a pure quadrature of the analytic
      surface.
    * no field: ``H`` from the cotangent Laplacian and ``K`` from the angle defect,
      averaged over the three vertices of the cell. This is what a consumer of an
      exported mesh can recompute.
    * ``crease_angle_degrees``: cells touching an edge sharper than this, or an open
      boundary edge, are reported separately. ``None`` treats every cell as smooth
      (only meaningful on a closed smooth mesh -- the sphere test does that and the
      Wigley test shows why a crease threshold is needed).
    """
    if kind not in ENERGY_KINDS:
        raise CurvatureError(f"unknown energy kind {kind!r}; expected one of {ENERGY_KINDS}")
    if not mesh.faces:
        raise CurvatureError("mesh has no faces")
    if crease_angle_degrees is not None and not (
        math.isfinite(crease_angle_degrees) and 0.0 < crease_angle_degrees < 180.0
    ):
        raise CurvatureError(
            f"crease_angle_degrees must lie in (0, 180), got {crease_angle_degrees!r}"
        )

    crease_faces, boundary_faces = _crease_and_boundary_faces(mesh, crease_angle_degrees)
    sharp = crease_faces | boundary_faces

    terms: list[tuple[float, float, bool]] = []
    integrand_values: list[float] = []
    identity_residual = None

    if field is not None:
        centroids = [_face_centroid(mesh, face) for face in mesh.faces]
        samples = principal_curvatures(field, centroids, refine=False)
        bending_terms: list[tuple[float, float, bool]] = []
        willmore_terms: list[tuple[float, float, bool]] = []
        gaussian_terms: list[tuple[float, float, bool]] = []
        residuals: list[float] = []
        for index, (face, sample) in enumerate(zip(mesh.faces, samples, strict=True)):
            area = _face_area(mesh, face)
            h = sample.mean_curvature
            bending = sample.kappa_max**2 + sample.kappa_min**2
            willmore = h * h
            gauss = sample.gaussian_curvature
            residuals.append(abs(bending - (4.0 * willmore - 2.0 * gauss)))
            flag = index not in sharp
            bending_terms.append((bending, area, flag))
            willmore_terms.append((willmore, area, flag))
            gaussian_terms.append((gauss, area, flag))
            if kind == "bending":
                integrand_values.append(bending)
            elif kind == "willmore":
                integrand_values.append(willmore)
            else:
                integrand_values.append(gauss)
        identity_residual = max(residuals)
        terms = {
            "bending": bending_terms,
            "willmore": willmore_terms,
            "gaussian": gaussian_terms,
        }[kind]
        source = "field"
        notes = (
            "principal curvatures evaluated exactly at face centroids",
            "identity |(k1^2+k2^2) - (4H^2-2K)| reported as identity_residual",
        )
    else:
        # The mesh path deliberately does **not** go through the (H, K) -> (k1, k2)
        # inversion. On a real mesh the cotangent mean curvature and the angle-defect
        # Gaussian curvature are two different discretisations, and their pointwise
        # inconsistency is routinely large enough to make H**2 - K negative (on a
        # conditioned marching-tet sphere it is negative at half the vertices), which
        # would silently drop most of the surface from the integral. The densities
        # below only need H and K, so they are used directly:
        #     willmore = H**2,  gaussian = K,  bending = 4 H**2 - 2 K.
        discrete = gaussian_curvature(mesh)
        mean = vertex_mean_curvature(mesh)
        for index, face in enumerate(mesh.faces):
            corners_h = [mean[vertex] for vertex in face]
            corners_k = [discrete.curvature[vertex] for vertex in face]
            if any(not math.isfinite(value) for value in corners_h + corners_k):
                continue
            h = math.fsum(corners_h) / 3.0
            k = math.fsum(corners_k) / 3.0
            value = {
                "bending": 4.0 * h * h - 2.0 * k,
                "willmore": h * h,
                "gaussian": k,
            }[kind]
            area = _face_area(mesh, face)
            terms.append((value, area, index not in sharp))
            integrand_values.append(value)
        source = "mesh"
        notes = (
            "H from cotangent Laplacian, K from angle defect; mean of the three vertices",
            f"cells dropped for undefined vertex curvature: {len(mesh.faces) - len(terms)}",
            "bending density taken as 4H^2-2K, not via the H/K -> k1,k2 inversion",
        )

    return smoother(
        kind=kind,
        source=source,
        terms=terms,
        crease_angle_degrees=crease_angle_degrees,
        integrand_values=integrand_values,
        identity_residual=identity_residual,
        notes=notes,
    )


@dataclass(frozen=True)
class FairnessEnergyReport:
    """All three energies of one surface, plus the identity that ties them together."""

    bending: SmoothingEnergy
    willmore: SmoothingEnergy
    gaussian: SmoothingEnergy
    source: str

    @property
    def identity_gap(self) -> float:
        """``E_bending - (4 E_willmore - 2 E_gauss)`` on the smooth part.

        Zero to round-off on the field path (same eigenvalues), and the
        interpolation residual on the mesh path -- which is a genuine measure of how
        far the two discrete operators (cotangent, angle defect) are from being a
        consistent pair.
        """
        return self.bending.smooth_energy - (
            4.0 * self.willmore.smooth_energy - 2.0 * self.gaussian.smooth_energy
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA_REPORT,
            "source": self.source,
            "bending": self.bending.to_dict(),
            "willmore": self.willmore.to_dict(),
            "gaussian": self.gaussian.to_dict(),
            "identity_gap": self.identity_gap,
            "willmore_inequality_slack": self.willmore.smooth_energy - 4.0 * math.pi,
        }


def fairness_energy_report(
    mesh: Mesh,
    *,
    field: GeometryField | None = None,
    crease_angle_degrees: float | None = None,
) -> FairnessEnergyReport:
    """``E_bending``, ``E_willmore`` and ``E_gauss`` for ``mesh`` in one pass each."""
    energies = {
        kind: smoothing_energy(
            mesh, kind=kind, field=field, crease_angle_degrees=crease_angle_degrees
        )
        for kind in ENERGY_KINDS
    }
    return FairnessEnergyReport(
        bending=energies["bending"],
        willmore=energies["willmore"],
        gaussian=energies["gaussian"],
        source="field" if field is not None else "mesh",
    )


__all__ = [
    "ENERGY_KINDS",
    "FairnessEnergyReport",
    "MeshCurvatureField",
    "SCHEMA_ENERGY",
    "SCHEMA_FIELD",
    "SCHEMA_REPORT",
    "SmoothingEnergy",
    "fairness_energy_report",
    "principal_curvature_field",
    "smoothing_energy",
]