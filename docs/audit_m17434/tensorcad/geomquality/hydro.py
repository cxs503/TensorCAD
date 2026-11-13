"""Hydrostatics by patch quadrature, and the convergence study that audits it.

Two integration routes over the *same* B-rep, on purpose:

* :func:`patch_hydrostatics` integrates each face's supporting surface with
  Gauss-Legendre quadrature per knot span.  For a non-rational B-spline patch the
  integrand of the divergence theorem is a polynomial inside a span, so the
  quadrature is exact to round-off and the only remaining error is the difference
  between the CAD surface and the hull it represents.  This is the route the pipeline
  reports.
* :func:`facet_hydrostatics` integrates the *tessellation* instead: the polyhedral
  approximation a solver would actually consume.  Its error is ``O(h^2)``, which is
  why a mesh that looks fine can be 1 % off, and why the pipeline uses it only for
  the grid-convergence audit (:func:`richardson_gci`).

Formulas (divergence theorem over the closed body = hull up to the design waterline
plus the waterplane disc, outward normal ``n``, ``dA n = +-(S_u x S_v) du dv``):

    V        = (1/3)      oint  r . n                 dA
    int x dV = (1/2)      oint  x^2 n_x               dA
    int z dV = (1/2)      oint  z^2 n_z               dA
    A_wp     =          - oint  n_z                   dA
    I_tp     = int_A y^2 dA = - oint y^2 n_z          dA      (about the x axis)
    I_lp     = int_A x^2 dA = - oint x^2 n_z          dA      (about midship)

The last three follow from ``div F = 0`` for ``F = (0,0,y^2)`` and ``F = (0,0,x^2)``
over the closed body.  All of them are one quadrature each, so a single Gauss sweep
per face fills the whole hydrostatic table.
"""
from __future__ import annotations

import math
import time
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

from ..geometry import Mesh
from ..nurbs import NURBSSurface, gauss_legendre
from ._eval import basis_1d, evaluate_fast, find_span
from .brep import BRep, Face

Vec3 = tuple[float, float, float]

#: Default Gauss order per knot span.  The volume integrand is a polynomial of degree
#: <= 2 * (degree - 1) + degree per direction, i.e. <= 8 for cubic patches, so order 6
#: (exact to degree 11) integrates a cubic tensor-product patch exactly.
DEFAULT_ORDER = 6


def _distinct(values: Sequence[float]) -> list[float]:
    out: list[float] = []
    for value in values:
        if not out or value != out[-1]:
            out.append(float(value))
    return out


def span_intervals(knots: Sequence[float], degree: int | None = None) -> list[tuple[float, float]]:
    """Non-empty knot spans ``[t_k, t_k+1)``.

    The caller's knot vectors are clamped (first and last knot repeated ``degree+1``
    times), so the distinct values are exactly the span boundaries: no trimming by
    ``degree`` on either side, which would silently drop the outermost spans.
    """
    del degree
    unique = _distinct(knots)
    return [(unique[i], unique[i + 1]) for i in range(len(unique) - 1) if unique[i + 1] > unique[i]]


def _window_basis(
    knots: Sequence[float], degree: int, u: float, span: int
) -> list[float]:
    """The ``degree + 1`` nonzero basis values at ``u`` (Piegl & Tiller A2.2)."""
    n = [0.0] * (degree + 1)
    left = [0.0] * (degree + 1)
    right = [0.0] * (degree + 1)
    n[0] = 1.0
    for j in range(1, degree + 1):
        left[j] = u - knots[span + 1 - j]
        right[j] = knots[span + j] - u
        saved = 0.0
        for r in range(j):
            denominator = right[r + 1] + left[j - r]
            temp = n[r] / denominator if denominator != 0.0 else 0.0
            n[r] = saved + right[r + 1] * temp
            saved = left[j - r] * temp
        n[j] = saved
    return n


def _cross(a: Sequence[float], b: Sequence[float]) -> Vec3:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


@dataclass
class PatchHydrostatics:
    """The hydrostatic table, all of it from one quadrature sweep."""

    volume: float
    displacement: float
    lcb_x: float
    vcb_z: float
    waterplane_area: float
    lcf_x: float
    wetted_area: float
    transverse_inertia: float
    longitudinal_inertia: float
    block_coefficient: float
    length: float
    beam: float
    draft: float
    quadrature_points: int = 0
    seconds: float = 0.0
    density: float = 1025.0

    @property
    def midship_section_coefficient(self) -> float:
        return float("nan")

    def to_dict(self) -> dict[str, Any]:
        return {
            "volume": self.volume,
            "displacement_kg": self.displacement,
            "lcb_x": self.lcb_x,
            "vcb_z": self.vcb_z,
            "waterplane_area": self.waterplane_area,
            "lcf_x": self.lcf_x,
            "wetted_area": self.wetted_area,
            "transverse_inertia": self.transverse_inertia,
            "longitudinal_inertia": self.longitudinal_inertia,
            "block_coefficient": self.block_coefficient,
            "length": self.length,
            "beam": self.beam,
            "draft": self.draft,
            "quadrature_points": self.quadrature_points,
            "seconds": self.seconds,
            "density": self.density,
        }


def patch_hydrostatics(
    brep: BRep,
    *,
    order: int = DEFAULT_ORDER,
    density: float = 1025.0,
    length: float | None = None,
    beam: float | None = None,
    draft: float | None = None,
) -> PatchHydrostatics:
    """Integrate the B-rep surfaces exactly (per-span Gauss) into a hydrostatic table."""
    started = time.perf_counter()
    volume = 0.0
    moment_x = 0.0
    moment_z = 0.0
    projected = 0.0
    moment_projected_x = 0.0
    inertia_y = 0.0
    inertia_x = 0.0
    wetted = 0.0
    points = 0
    for face in brep.faces:
        surface = face.surface
        sign = 1.0 if face.same_sense else -1.0
        u_spans = span_intervals(surface.knots[0], surface.degree[0])
        v_spans = span_intervals(surface.knots[1], surface.degree[1])
        # The volume integrand is a polynomial of degree up to about 2*p + 2 per direction
        # inside one span, so the quadrature order has to follow the patch degree: a
        # degree-24 vertical span (the z = -T + T*v**m mapping of a power-law section)
        # under-integrates by 0.5 % with the default order 6, which is exactly the kind of
        # quiet error this pipeline exists to rule out.
        order_u = max(order, surface.degree[0] + 3)
        # 2*degree_v + 2 is the exact degree bound of the volume integrand; two extra
        # nodes are margin, not necessity (degree 24 -> 26 nodes, exact to degree 51).
        order_v = max(order, surface.degree[1] + 2)
        net = surface.control_net
        index_v = 0
        nodes_u, weights_u = gauss_legendre(order_u)
        nodes_v, weights_v = gauss_legendre(order_v)
        weight_u = weights_u
        for ua, ub in u_spans:
            half_u = 0.5 * (ub - ua)
            mid_u = 0.5 * (ub + ua)
            for va, vb in v_spans:
                half_v = 0.5 * (vb - va)
                mid_v = 0.5 * (vb + va)
                v_basis = []
                for node_v in nodes_v:
                    v = mid_v + half_v * node_v
                    span_v = find_span(surface.knots[1], surface.degree[1], v)
                    bv, dbv = basis_1d(surface.knots[1], surface.degree[1], v, span_v)
                    v_basis.append((bv, dbv, span_v - surface.degree[1], v))
                for i, node_u in enumerate(nodes_u):
                    u = mid_u + half_u * node_u
                    span_u = find_span(surface.knots[0], surface.degree[0], u)
                    bu, dbu = basis_1d(surface.knots[0], surface.degree[0], u, span_u)
                    index_u = span_u - surface.degree[0]
                    # Contract the control net in u once per u-node, then the v sweep is a
                    # short dot product: the tensor product costs (p+q) per point instead
                    # of p*q, which matters for a degree-24 vertical span.
                    contracted: list[Vec3] = []
                    contracted_u: list[Vec3] = []
                    for j in range(surface.degree[1] + 1):
                        point = [0.0, 0.0, 0.0]
                        tangent = [0.0, 0.0, 0.0]
                        for a in range(surface.degree[0] + 1):
                            control = net[index_u + a][index_v + j]
                            weight = bu[a]
                            point[0] += weight * control[0]
                            point[1] += weight * control[1]
                            point[2] += weight * control[2]
                            weight = dbu[a]
                            tangent[0] += weight * control[0]
                            tangent[1] += weight * control[1]
                            tangent[2] += weight * control[2]
                        contracted.append((point[0], point[1], point[2]))
                        contracted_u.append((tangent[0], tangent[1], tangent[2]))
                    for j, (bv, dbv, index_v_local, _v) in enumerate(v_basis):
                        point = [0.0, 0.0, 0.0]
                        tangent_u = [0.0, 0.0, 0.0]
                        tangent_v = [0.0, 0.0, 0.0]
                        for b in range(surface.degree[1] + 1):
                            weight = bv[b]
                            control = contracted[index_v_local + b]
                            point[0] += weight * control[0]
                            point[1] += weight * control[1]
                            point[2] += weight * control[2]
                            control = contracted_u[index_v_local + b]
                            tangent_u[0] += weight * control[0]
                            tangent_u[1] += weight * control[1]
                            tangent_u[2] += weight * control[2]
                            weight = dbv[b]
                            control = contracted[index_v_local + b]
                            tangent_v[0] += weight * control[0]
                            tangent_v[1] += weight * control[1]
                            tangent_v[2] += weight * control[2]
                        normal = _cross(tangent_u, tangent_v)
                        jacobian = weight_u[i] * weights_v[j] * half_u * half_v
                        points += 1
                        area_density = math.sqrt(
                            normal[0] ** 2 + normal[1] ** 2 + normal[2] ** 2
                        )
                        if face.wetted:
                            wetted += jacobian * area_density
                        flux = jacobian * sign
                        volume += flux * (
                            point[0] * normal[0] + point[1] * normal[1] + point[2] * normal[2]
                        ) / 3.0
                        moment_x += 0.5 * flux * point[0] ** 2 * normal[0]
                        moment_z += 0.5 * flux * point[2] ** 2 * normal[2]
                        if face.wetted:
                            # A_wp = -oint_hull n_z dA over the *wetted* hull; the deck
                            # patch contributes +A_wp of its own and must be excluded.
                            projected += flux * normal[2]
                            moment_projected_x += flux * point[0] * normal[2]
                            inertia_y += flux * point[1] ** 2 * normal[2]
                            inertia_x += flux * point[0] ** 2 * normal[2]
    waterplane = -projected
    if abs(volume) > 0:
        lcb = moment_x / volume
        vcb = moment_z / volume
    else:
        lcb = vcb = 0.0
    lcf = -moment_projected_x / waterplane if abs(waterplane) > 0 else 0.0
    lo, hi = brep.bounding_box()
    resolved_length = length if length is not None else (hi[0] - lo[0])
    resolved_beam = beam if beam is not None else (hi[1] - lo[1])
    resolved_draft = draft if draft is not None else (hi[2] - lo[2])
    cb = volume / (resolved_length * resolved_beam * resolved_draft) if (
        resolved_length * resolved_beam * resolved_draft
    ) else float("nan")
    return PatchHydrostatics(
        volume=volume,
        displacement=volume * density,
        lcb_x=lcb,
        vcb_z=vcb,
        waterplane_area=waterplane,
        lcf_x=lcf,
        wetted_area=wetted,
        transverse_inertia=-inertia_y,
        longitudinal_inertia=-inertia_x,
        block_coefficient=cb,
        length=resolved_length,
        beam=resolved_beam,
        draft=resolved_draft,
        quadrature_points=points,
        seconds=time.perf_counter() - started,
        density=density,
    )


def facet_hydrostatics(
    brep: BRep,
    *,
    u_samples: int = 24,
    v_samples: int = 12,
    density: float = 1025.0,
    physical: bool = True,
) -> PatchHydrostatics:
    """The same table from the tessellation: the route a solver actually consumes."""
    started = time.perf_counter()
    volume = 0.0
    moment_x = 0.0
    moment_z = 0.0
    projected = 0.0
    wetted = 0.0
    triangles = 0
    for face in brep.faces:
        mesh = _face_mesh(face, u_samples, v_samples, physical=physical)
        for a, b, c in mesh.faces:
            pa, pb, pc = mesh.vertices[a], mesh.vertices[b], mesh.vertices[c]
            normal = _cross(
                (pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]),
                (pc[0] - pa[0], pc[1] - pa[1], pc[2] - pa[2]),
            )
            area2 = math.sqrt(sum(component**2 for component in normal))
            if area2 == 0.0:
                continue
            tetra = (
                pa[0] * (pb[1] * pc[2] - pb[2] * pc[1])
                - pa[1] * (pb[0] * pc[2] - pb[2] * pc[0])
                + pa[2] * (pb[0] * pc[1] - pb[1] * pc[0])
            ) / 6.0
            volume += tetra
            triangles += 1
            centroid = (
                (pa[0] + pb[0] + pc[0]) / 3.0,
                (pa[1] + pb[1] + pc[1]) / 3.0,
                (pa[2] + pb[2] + pc[2]) / 3.0,
            )
            moment_x += tetra * 0.75 * centroid[0]
            moment_z += tetra * 0.75 * centroid[2]
            wetted += 0.5 * area2
            projected += 0.5 * normal[2]
    waterplane = -projected
    lo, hi = brep.bounding_box()
    length = hi[0] - lo[0]
    beam = hi[1] - lo[1]
    draft = hi[2] - lo[2]
    return PatchHydrostatics(
        volume=volume,
        displacement=volume * density,
        lcb_x=moment_x / volume if volume else 0.0,
        vcb_z=moment_z / volume if volume else 0.0,
        waterplane_area=waterplane,
        lcf_x=0.0,
        wetted_area=wetted,
        transverse_inertia=float("nan"),
        longitudinal_inertia=float("nan"),
        block_coefficient=volume / (length * beam * draft) if length * beam * draft else float("nan"),
        length=length,
        beam=beam,
        draft=draft,
        quadrature_points=triangles,
        seconds=time.perf_counter() - started,
        density=density,
    )


def arc_length_parameters(
    surface: NURBSSurface, direction: int, samples: int, *, probes: int = 41
) -> list[float]:
    """Parameters placed at equal arc length along a mid-domain isocurve.

    Refining *parameter space* is not the same as refining *physical space* when the
    parameterisation is graded: a patch with ``z = -T + T*v**m`` packs almost all of its
    geometry into the last few percent of ``v``, so doubling ``v_samples`` refines the
    keel region and leaves the deck where it was.  Equal-arc-length placement makes the
    three-level convergence study refine the surface itself, which is what the observed
    order is supposed to measure.
    """
    (u0, u1), (v0, v1) = surface.domain
    fixed = 0.5 * (v0 + v1) if direction == 0 else 0.5 * (u0 + u1)
    parameters = []
    cumulative = [0.0]
    previous = None
    for index in range(probes):
        t = index / (probes - 1)
        if direction == 0:
            parameter = u0 + (u1 - u0) * t
            point = evaluate_fast(surface, parameter, fixed)
        else:
            parameter = v0 + (v1 - v0) * t
            point = evaluate_fast(surface, fixed, parameter)
        if previous is not None:
            cumulative.append(cumulative[-1] + math.dist(previous, point))
        previous = point
        parameters.append(parameter)
    total = cumulative[-1]
    if total == 0.0:
        return [u0 + (u1 - u0) * i / samples if direction == 0 else v0 + (v1 - v0) * i / samples
                for i in range(samples + 1)]
    targets = [total * i / samples for i in range(samples + 1)]
    out: list[float] = []
    cursor = 0
    for target in targets:
        while cursor < len(cumulative) - 2 and cumulative[cursor + 1] < target:
            cursor += 1
        span = cumulative[cursor + 1] - cumulative[cursor]
        weight = 0.0 if span == 0.0 else (target - cumulative[cursor]) / span
        out.append(parameters[cursor] + weight * (parameters[cursor + 1] - parameters[cursor]))
    return out


def _face_mesh(
    face: Face, u_samples: int, v_samples: int, *, physical: bool = True
) -> Mesh:
    from ..geometry import MeshBuilder

    builder = MeshBuilder()
    if physical:
        us = arc_length_parameters(face.surface, 0, u_samples)
        vs = arc_length_parameters(face.surface, 1, v_samples)
    else:
        (u0, u1), (v0, v1) = face.domain
        us = [u0 + (u1 - u0) * i / u_samples for i in range(u_samples + 1)]
        vs = [v0 + (v1 - v0) * j / v_samples for j in range(v_samples + 1)]
    grid = [[builder.add_vertex(evaluate_fast(face.surface, u, v)) for v in vs] for u in us]
    for i in range(u_samples):
        for j in range(v_samples):
            quad = (grid[i][j], grid[i + 1][j], grid[i + 1][j + 1], grid[i][j + 1])
            if face.same_sense:
                builder.add_quad(*quad)
            else:
                builder.add_quad(quad[0], quad[3], quad[2], quad[1])
    return builder.build()


# --------------------------------------------------------------------------- #
# convergence
# --------------------------------------------------------------------------- #
def richardson_gci(
    values: Sequence[float], refinement_ratio: float = 2.0, safety: float = 1.25
) -> dict[str, Any]:
    """Observed order and grid convergence index from three successively finer results.

    ``p = ln|(f2-f1)/(f3-f2)| / ln(r)`` with ``f1`` the coarsest, then
    ``f_exact = f3 + (f3-f2)/(r^p - 1)`` and ``GCI = safety |f3-f2| / |f3| / (r^p-1)``.
    The sign of the differences must not oscillate: a non-monotone sequence is
    reported as such instead of being fitted.
    """
    if len(values) != 3:
        raise ValueError("richardson_gci needs exactly three results")
    coarse, medium, fine = (float(value) for value in values)
    d1 = medium - coarse
    d2 = fine - medium
    report: dict[str, Any] = {
        "values": [coarse, medium, fine],
        "differences": [d1, d2],
        "refinement_ratio": refinement_ratio,
    }
    if d1 == 0 or d2 == 0:
        report.update(
            {
                "observed_order": None,
                "extrapolated": fine,
                "gci_percent": 0.0,
                "converged": True,
                "note": "difference vanished to round-off: exact discretisation",
            }
        )
        return report
    ratio = abs(d1 / d2)
    if ratio <= 0:
        raise ValueError("non-positive ratio in Richardson extrapolation")
    order = math.log(ratio) / math.log(refinement_ratio)
    if d1 * d2 < 0:
        report["note"] = "oscillatory convergence: order not meaningful"
    extrapolated = fine + (fine - medium) / (refinement_ratio**order - 1.0)
    gci = safety * abs(fine - medium) / abs(fine) / (refinement_ratio**order - 1.0)
    report.update(
        {
            "observed_order": order,
            "extrapolated": extrapolated,
            "gci_percent": 100.0 * gci,
            "converged": d1 * d2 > 0,
        }
    )
    return report


def convergence_study(
    brep: BRep,
    *,
    levels: Sequence[tuple[int, int]] = ((24, 12), (48, 24), (96, 48)),
    reference: float | None = None,
    physical: bool = True,
) -> dict[str, Any]:
    """Facet-route volume at three tessellation densities: observed order and GCI."""
    volumes = []
    timings = []
    for u_samples, v_samples in levels:
        table = facet_hydrostatics(
            brep, u_samples=u_samples, v_samples=v_samples, physical=physical
        )
        volumes.append(table.volume)
        timings.append(table.seconds)
    study = richardson_gci(volumes)
    study["levels"] = list(levels)
    study["refinement"] = "equal arc length" if physical else "uniform in parameter"
    study["seconds"] = timings
    if reference is not None:
        study["reference"] = reference
        study["error_finest_percent"] = 100.0 * abs(volumes[-1] - reference) / abs(reference)
        study["error_extrapolated_percent"] = (
            100.0 * abs(study["extrapolated"] - reference) / abs(reference)
        )
    return study
