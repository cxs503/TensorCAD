"""p-order quadrature on parametric surfaces, and the curved element it exports.

Why this module exists
----------------------
A solver that consumes a NURBS patch needs four numbers per quadrature point -- the
position ``x``, the weight ``w``, the unit normal ``n`` and the surface Jacobian
``|S_u x S_v|`` -- *and* it needs to know the polynomial order those numbers are
accurate to. This module produces both, and it measures the second one instead of
asserting it:

* A ``p``-point Gauss-Legendre rule integrates polynomials of degree ``2p - 1`` exactly
  and degree ``2p`` not at all (:func:`gauss_exactness_probe` reports both sides, so the
  rule's exactness order is pinned from above and below rather than claimed).
* On a curved patch the quadrature error of the area, the volume and the curvature
  integrals must decay at least as fast as that exactness order under refinement
  (:func:`convergence_study` plus :func:`observed_orders`).
* The negative control is the same machinery applied to the *faceted* patch
  (:func:`tessellation_convergence`): a patch exported as flat triangles converges at
  second order no matter what ``p`` is, which is exactly the failure the p-order claims
  are supposed to rule out.

Named limitations
-----------------
* The rule is tensor-product Gauss-Legendre on the reference square. A simplex rule for
  triangular elements is deliberately absent: a tensor-product NURBS patch has a
  quadrilateral parameter domain, and exposing a "triangle rule" here would invite a
  quadrature that no patch in this repository can back.
* The integration is per knot span (each span of a degree-``p`` patch is one Bezier
  segment); crossing a knot is never done with a single rule, because the integrand is
  only piecewise smooth there. That is also why the reported "order" is a property of
  the *rule on one span* plus the subdivision count, not of the whole patch.
* "Observed order" is the pairwise ratio ``log2(e_{n} / e_{2n})``; once the error is at
  the round-off floor the ratio is noise, so :func:`observed_orders` returns ``None`` for
  such a pair instead of manufacturing a number out of rounding error.
* The faceted integrator is exact for a flat triangle (that is the point of a control):
  it is not meant to be a good area estimator for curved geometry.
"""

from __future__ import annotations

import itertools
import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Any

from ..geometry import Mesh, Vec3
from ..nurbs import NURBSSurface, gauss_legendre

Vec2 = tuple[float, float]

#: A scalar field on the parameter domain, ``f(u, v) -> float``.
Field = Callable[[float, float], float]

#: Errors below this fraction of the integral are treated as round-off, not as signal.
NOISE_FLOOR = 1e-13

__all__ = [
    "NOISE_FLOOR",
    "CurvedElement",
    "Field",
    "convergence_study",
    "curved_elements",
    "faceted_area_and_volume",
    "gauss_exactness_probe",
    "gauss_rule",
    "integrate_surface",
    "knot_spans",
    "mesh_area_and_volume",
    "observed_orders",
    "surface_area",
    "surface_volume",
    "tessellation_convergence",
]


# ----------------------------------------------------------------------------------
# rules and spans
# ----------------------------------------------------------------------------------
def gauss_rule(order: int) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """``order``-point Gauss-Legendre nodes and weights on ``[-1, 1]``."""
    if not isinstance(order, int) or order < 1:
        raise ValueError(f"quadrature order must be a positive integer, got {order!r}")
    return gauss_legendre(order)


def knot_spans(knots: Sequence[float], degree: int) -> tuple[Vec2, ...]:
    """The distinct non-degenerate spans of a knot vector, in order."""
    if degree < 0:
        raise ValueError(f"degree must be non-negative, got {degree!r}")
    if len(knots) < 2 * (degree + 1):
        raise ValueError(f"a degree-{degree} clamped knot vector needs at least "
                         f"{2 * (degree + 1)} entries, got {len(knots)}")
    low = knots[degree]
    high = knots[len(knots) - degree - 1]
    values = sorted({float(value) for value in knots if low <= value <= high})
    spans: list[Vec2] = []
    for left, right in itertools.pairwise(values):
        if right - left > 0.0:
            spans.append((left, right))
    if not spans:
        raise ValueError("knot vector has no positive-length span")
    return tuple(spans)


# ----------------------------------------------------------------------------------
# integration
# ----------------------------------------------------------------------------------
def _check_subdivisions(subdivisions: tuple[int, int]) -> tuple[int, int]:
    nu, nv = int(subdivisions[0]), int(subdivisions[1])
    if nu < 1 or nv < 1:
        raise ValueError(f"subdivisions must be positive, got {subdivisions!r}")
    return nu, nv


def _cells(
    surface: NURBSSurface, subdivisions: tuple[int, int]
) -> list[tuple[Vec2, Vec2]]:
    nu, nv = _check_subdivisions(subdivisions)
    cells: list[tuple[Vec2, Vec2]] = []
    for ua, ub in knot_spans(surface.knots[0], surface.degree[0]):
        for va, vb in knot_spans(surface.knots[1], surface.degree[1]):
            hu = (ub - ua) / nu
            hv = (vb - va) / nv
            for iu in range(nu):
                for iv in range(nv):
                    cells.append(
                        ((ua + iu * hu, ua + (iu + 1) * hu), (va + iv * hv, va + (iv + 1) * hv))
                    )
    return cells


def integrate_surface(
    surface: NURBSSurface,
    field: Field,
    *,
    order: int = 8,
    subdivisions: tuple[int, int] = (1, 1),
) -> float:
    """``int int f(u, v) du dv`` over the patch, Gauss-Legendre per span and cell."""
    return float(integrate_surface_report(
        surface, field, order=order, subdivisions=subdivisions
    )["value"])


def integrate_surface_report(
    surface: NURBSSurface,
    field: Field,
    *,
    order: int = 8,
    subdivisions: tuple[int, int] = (1, 1),
) -> dict[str, Any]:
    """The integral plus how it was computed (order, cells, evaluations)."""
    nodes, weights = gauss_rule(order)
    cells = _cells(surface, subdivisions)
    total = 0.0
    evaluations = 0
    for (ua, ub), (va, vb) in cells:
        u_half = 0.5 * (ub - ua)
        v_half = 0.5 * (vb - va)
        scale = u_half * v_half
        for node_u, weight_u in zip(nodes, weights, strict=True):
            u = u_half * (node_u + 1.0) + ua
            for node_v, weight_v in zip(nodes, weights, strict=True):
                v = v_half * (node_v + 1.0) + va
                total += weight_u * weight_v * scale * field(u, v)
                evaluations += 1
    return {
        "value": total,
        "order": order,
        "subdivisions": list(subdivisions),
        "cells": len(cells),
        "evaluations": evaluations,
    }


def area_element(surface: NURBSSurface, u: float, v: float) -> float:
    """``|S_u x S_v|``: the surface Jacobian at ``(u, v)``."""
    _, _, cross = surface.jacobian(u, v)
    return math.sqrt(cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2])


def volume_element(surface: NURBSSurface, u: float, v: float) -> float:
    """Integrand of the divergence-theorem volume ``(1/3) x . (S_u x S_v)``.

    Summing this over a *closed* patch gives the enclosed volume; on an open patch it is
    the volume of the cone over the patch with apex at the origin and carries no physical
    meaning on its own.
    """
    point = surface.evaluate(u, v)
    _, _, cross = surface.jacobian(u, v)
    return (point[0] * cross[0] + point[1] * cross[1] + point[2] * cross[2]) / 3.0


def surface_area(
    surface: NURBSSurface, *, order: int = 8, subdivisions: tuple[int, int] = (1, 1)
) -> float:
    """Patch area ``int int |S_u x S_v| du dv``."""
    return integrate_surface(
        surface, lambda u, v: area_element(surface, u, v), order=order, subdivisions=subdivisions
    )


def surface_volume(
    surface: NURBSSurface, *, order: int = 8, subdivisions: tuple[int, int] = (1, 1)
) -> float:
    """Enclosed volume of a closed patch, by the divergence theorem."""
    return integrate_surface(
        surface, lambda u, v: volume_element(surface, u, v), order=order, subdivisions=subdivisions
    )


# ----------------------------------------------------------------------------------
# convergence bookkeeping
# ----------------------------------------------------------------------------------
def observed_orders(
    levels: Sequence[float], errors: Sequence[float], *, floor: float = NOISE_FLOOR
) -> list[float | None]:
    """Pairwise ``log(e_i / e_{i+1}) / log(h_i / h_{i+1})``, ``None`` at the noise floor.

    ``levels`` are the refinement parameters (Gauss order or tessellation count) and are
    *not* assumed to double, so the denominator is computed rather than hard-coded to
    ``log 2``. A pair whose finer error is at the noise floor cannot support a statement
    about convergence and yields ``None``.
    """
    if len(levels) != len(errors):
        raise ValueError(f"{len(levels)} levels against {len(errors)} errors")
    orders: list[float | None] = []
    for index in range(len(levels) - 1):
        fine, coarse = errors[index + 1], errors[index]
        if fine <= floor or coarse <= 0.0 or fine <= 0.0:
            orders.append(None)
            continue
        ratio = math.log(coarse / fine)
        span = math.log(levels[index + 1] / levels[index])
        orders.append(ratio / span)
    return orders


def fitted_order(
    levels: Sequence[float], errors: Sequence[float], *, floor: float = NOISE_FLOOR
) -> float | None:
    """Least-squares slope of ``log e`` against ``log level`` above the noise floor."""
    points = [
        (math.log(level), math.log(error))
        for level, error in zip(levels, errors, strict=True)
        if error > floor
    ]
    if len(points) < 2:
        return None
    mean_x = sum(x for x, _ in points) / len(points)
    mean_y = sum(y for _, y in points) / len(points)
    numerator = sum((x - mean_x) * (y - mean_y) for x, y in points)
    denominator = sum((x - mean_x) ** 2 for x, _ in points)
    if denominator == 0.0:
        return None
    return -numerator / denominator


def convergence_study(
    surface: NURBSSurface,
    field: Field,
    orders: Sequence[int] = (2, 3, 4, 6, 8, 12, 16),
    *,
    reference_order: int = 40,
    subdivisions: tuple[int, int] = (1, 1),
    exact: float | None = None,
) -> dict[str, Any]:
    """Quadrature error against Gauss order, with the observed convergence order.

    The reference is the same rule at ``reference_order`` unless ``exact`` (an analytic
    value) is supplied. Relative errors are taken against that reference, and the
    observed order is reported both pairwise and as a fitted slope over the orders whose
    error is above :data:`NOISE_FLOOR`.
    """
    reference = (
        float(exact)
        if exact is not None
        else integrate_surface(surface, field, order=reference_order, subdivisions=subdivisions)
    )
    if reference == 0.0:
        raise ValueError("reference integral is zero; relative errors are undefined")
    rows: list[dict[str, Any]] = []
    for order in orders:
        value = integrate_surface(surface, field, order=order, subdivisions=subdivisions)
        rows.append(
            {
                "order": order,
                "value": value,
                "relative_error": abs(value - reference) / abs(reference),
            }
        )
    levels = [row["order"] for row in rows]
    errors = [row["relative_error"] for row in rows]
    return {
        "reference": reference,
        "reference_order": None if exact is not None else reference_order,
        "subdivisions": list(subdivisions),
        "rows": rows,
        "pairwise_orders": observed_orders(levels, errors),
        "fitted_order": fitted_order(levels, errors),
    }


def gauss_exactness_probe(order: int) -> dict[str, Any]:
    """Degree-of-exactness of the ``order``-point rule, measured on ``[-1, 1]``.

    ``int_{-1}^{1} x^k dx`` is ``2 / (k + 1)`` for even ``k`` and zero for odd ``k``. The
    rule is exact at ``k = 2p - 1`` and wrong at ``k = 2p`` for every ``p >= 1``; both
    numbers are returned, so the exactness order is bounded from both sides.
    """
    nodes, weights = gauss_rule(order)

    def rule(power: int) -> float:
        return sum(weight * node**power for node, weight in zip(nodes, weights, strict=True))

    def exact(power: int) -> float:
        return 2.0 / (power + 1) if power % 2 == 0 else 0.0

    at_limit = abs(rule(2 * order - 1) - exact(2 * order - 1))
    beyond = abs(rule(2 * order) - exact(2 * order))
    return {
        "order": order,
        "degree_exact": 2 * order - 1,
        "error_at_degree_exact": at_limit,
        "degree_beyond": 2 * order,
        "error_one_degree_beyond": beyond,
        "exact_at_limit": at_limit <= 1e-14,
        "fails_beyond_limit": beyond > 1e-12,
    }


# ----------------------------------------------------------------------------------
# isoparametric curved elements
# ----------------------------------------------------------------------------------
@dataclass(frozen=True)
class CurvedElement:
    """One quadrature cell of an isoparametric ``p``-order surface element.

    The weights are already in surface measure: ``weights[i] = w_ref[i] * |J|[i] * scale``
    where ``scale`` converts the reference square ``[-1, 1]^2`` to the cell. So
    ``sum(weights)`` is the element area, ``sum`` over all elements is the patch area, and
    ``sum(weights * x . n / 3)`` is the element's divergence-theorem volume contribution.
    A solver consuming this needs no Jacobian of its own for surface integrals; it does
    get ``jacobian_det`` and ``reference_weights`` separately in case it wants to place a
    volume element on the surface.
    """

    element_id: int
    degree: tuple[int, int]
    span: tuple[Vec2, Vec2]
    quadrature_order: int
    metric_scale: float
    uv: tuple[Vec2, ...]
    points: tuple[Vec3, ...]
    reference_weights: tuple[float, ...]
    jacobian_det: tuple[float, ...]
    normals: tuple[Vec3, ...]

    @property
    def weights(self) -> tuple[float, ...]:
        """Quadrature weights in surface measure (``sum`` is the element area)."""
        return tuple(
            reference * scale * jac for reference, scale, jac in zip(
                self.reference_weights, [self.metric_scale] * len(self.reference_weights),
                self.jacobian_det, strict=True
            )
        )

    @property
    def area(self) -> float:
        return math.fsum(self.weights)

    @property
    def volume_contribution(self) -> float:
        return math.fsum(
            weight * (point[0] * normal[0] + point[1] * normal[1] + point[2] * normal[2]) / 3.0
            for weight, point, normal in zip(self.weights, self.points, self.normals, strict=True)
        )

    @property
    def centroid(self) -> Vec3:
        weights = self.weights
        total = math.fsum(weights)
        if total == 0.0:
            raise ValueError(f"element {self.element_id} has zero area")
        return (
            math.fsum(w * point[0] for w, point in zip(weights, self.points, strict=True)) / total,
            math.fsum(w * point[1] for w, point in zip(weights, self.points, strict=True)) / total,
            math.fsum(w * point[2] for w, point in zip(weights, self.points, strict=True)) / total,
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "element_id": self.element_id,
            "degree": list(self.degree),
            "span": [[self.span[0][0], self.span[0][1]], [self.span[1][0], self.span[1][1]]],
            "quadrature_order": self.quadrature_order,
            "metric_scale": self.metric_scale,
            "area": self.area,
            "volume_contribution": self.volume_contribution,
            "centroid": list(self.centroid),
            "points": len(self.points),
        }


def curved_elements(
    surface: NURBSSurface,
    *,
    order: int = 4,
    subdivisions: tuple[int, int] = (1, 1),
) -> tuple[CurvedElement, ...]:
    """The full isoparametric element set of a patch at quadrature ``order``.

    One element per knot span per subdivision; ``order`` must be at least ``p`` so that
    the element's own geometry (degree ``p``) is integrated consistently, and the caller
    is expected to report ``degree`` alongside the element so a consumer can check that.
    """
    p, q = surface.degree
    if order < 1:
        raise ValueError(f"order must be positive, got {order!r}")
    nodes, weights = gauss_rule(order)
    elements: list[CurvedElement] = []
    for (ua, ub), (va, vb) in _cells(surface, subdivisions):
        u_half = 0.5 * (ub - ua)
        v_half = 0.5 * (vb - va)
        scale = u_half * v_half
        uv: list[Vec2] = []
        points: list[Vec3] = []
        jacobians: list[float] = []
        normals: list[Vec3] = []
        reference: list[float] = []
        for node_u, weight_u in zip(nodes, weights, strict=True):
            u = u_half * (node_u + 1.0) + ua
            for node_v, weight_v in zip(nodes, weights, strict=True):
                v = v_half * (node_v + 1.0) + va
                uv.append((u, v))
                points.append(surface.evaluate(u, v))
                jacobians.append(area_element(surface, u, v))
                normals.append(surface.normal(u, v))
                reference.append(weight_u * weight_v)
        elements.append(
            CurvedElement(
                element_id=len(elements),
                degree=(p, q),
                span=((ua, ub), (va, vb)),
                quadrature_order=order,
                metric_scale=scale,
                uv=tuple(uv),
                points=tuple(points),
                reference_weights=tuple(reference),
                jacobian_det=tuple(jacobians),
                normals=tuple(normals),
            )
        )
    return tuple(elements)


# ----------------------------------------------------------------------------------
# meshes: the faceted negative control
# ----------------------------------------------------------------------------------
def mesh_area_and_volume(mesh: Mesh) -> tuple[float, float]:
    """Total area and signed enclosed volume of a triangle mesh.

    The volume is the sum of signed tetrahedra from the origin; its sign is the mesh's
    winding, so an inward-oriented shell reports a negative number rather than a
    plausible-looking positive one.
    """
    area = 0.0
    volume = 0.0
    for face in mesh.faces:
        a = mesh.vertices[face[0]]
        b = mesh.vertices[face[1]]
        c = mesh.vertices[face[2]]
        ab = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
        ac = (c[0] - a[0], c[1] - a[1], c[2] - a[2])
        cross = (
            ab[1] * ac[2] - ab[2] * ac[1],
            ab[2] * ac[0] - ab[0] * ac[2],
            ab[0] * ac[1] - ab[1] * ac[0],
        )
        area += 0.5 * math.sqrt(cross[0] ** 2 + cross[1] ** 2 + cross[2] ** 2)
        volume += (
            a[0] * (b[1] * c[2] - b[2] * c[1])
            - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0])
        ) / 6.0
    return area, volume


def faceted_area_and_volume(
    surface: NURBSSurface, samples: int
) -> tuple[float, float]:
    """Area and volume of a linearized (``samples x samples`` facet) copy of the patch."""
    if samples < 1:
        raise ValueError(f"samples must be positive, got {samples!r}")
    return mesh_area_and_volume(surface.tessellate(samples, samples).mesh)


def tessellation_convergence(
    surface: NURBSSurface,
    samples: Sequence[int] = (8, 16, 32, 64),
    *,
    exact_area: float,
    exact_volume: float,
) -> dict[str, Any]:
    """Observed order of a *faceted* patch: the negative control for the p-order claims.

    A chord replaces an arc to first order in the facet size while the area deficit is
    quadratic, so the faceted area and volume are expected to converge at order 2 -- and
    to be nowhere near a ``1e-6`` relative volume at any affordable facet count. Both
    statements are computed here so a gate can require them instead of trusting them.
    """
    rows: list[dict[str, Any]] = []
    for count in samples:
        area, volume = faceted_area_and_volume(surface, count)
        rows.append(
            {
                "samples": count,
                "facets": len(surface.tessellate(count, count).mesh.faces),
                "area": area,
                "volume": volume,
                "area_relative_error": abs(area - exact_area) / abs(exact_area),
                "volume_relative_error": abs(volume - exact_volume) / abs(exact_volume),
            }
        )
    levels = [row["samples"] for row in rows]
    areas = [row["area_relative_error"] for row in rows]
    volumes = [row["volume_relative_error"] for row in rows]
    return {
        "rows": rows,
        "area_pairwise_orders": observed_orders(levels, areas),
        "volume_pairwise_orders": observed_orders(levels, volumes),
        "area_fitted_order": fitted_order(levels, areas),
        "volume_fitted_order": fitted_order(levels, volumes),
    }
