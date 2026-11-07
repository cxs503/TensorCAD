"""Fast tensor-product evaluation for the geometry-quality pipeline.

``NURBSSurface.evaluate`` rebuilds and differentiates the derivative control nets on
every call, which is the right thing for a one-off evaluation and far too slow for the
quadrature and tessellation sweeps this package runs (a 96-span, degree-10 patch costs
milliseconds per point).  Inside one knot span the nonzero basis functions are a
``degree + 1`` window, so the point, ``S_u`` and ``S_v`` are a tensor contraction over
that window -- a few hundred flops, no allocation of derivative nets.

Non-rational surfaces take this path; a rational one falls back to the reference
implementation, because the quotient rule needs the homogeneous derivatives and that
case is not on the hot path here.
"""
from __future__ import annotations

import math
from collections.abc import Sequence

from ..nurbs import NURBSSurface

Vec3 = tuple[float, float, float]


def find_span(knots: Sequence[float], degree: int, u: float) -> int:
    """Index of the knot span containing ``u`` (clamped at the last span)."""
    count = len(knots) - degree - 1
    if u >= knots[-1]:
        return count - 1
    low, high = degree, count
    while high - low > 1:
        mid = (low + high) // 2
        if u < knots[mid]:
            high = mid
        else:
            low = mid
    return low


def window_basis(knots: Sequence[float], degree: int, u: float, span: int) -> list[float]:
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


def basis_1d(
    knots: Sequence[float], degree: int, u: float, span: int
) -> tuple[list[float], list[float]]:
    """Values and first derivatives of the nonzero basis functions at ``u``.

    ``N'_{i,p} = p/(t_{i+p}-t_i) N_{i,p-1} - p/(t_{i+p+1}-t_{i+1}) N_{i+1,p-1}`` with the
    degree-``p-1`` window from the same span; terms outside the window are zero, which is
    what makes the boundary cases exact rather than clipped.
    """
    values = window_basis(knots, degree, u, span)
    if degree == 0:
        return values, [0.0]
    lower = window_basis(knots, degree - 1, u, span)
    derivatives: list[float] = []
    for a in range(degree + 1):
        i = span - degree + a
        term = 0.0
        left = knots[i + degree] - knots[i]
        if left > 0 and a >= 1:
            term += degree * lower[a - 1] / left
        right = knots[i + degree + 1] - knots[i + 1]
        if right > 0 and a <= degree - 1:
            term -= degree * lower[a] / right
        derivatives.append(term)
    return values, derivatives


def eval_tensor(
    surface: NURBSSurface, u: float, v: float, span_u: int, span_v: int
) -> tuple[Vec3, Vec3, Vec3]:
    """``(S, S_u, S_v)`` for a non-rational tensor-product patch, by contraction."""
    net = surface.control_net
    degree_u, degree_v = surface.degree
    bu, dbu = basis_1d(surface.knots[0], degree_u, u, span_u)
    bv, dbv = basis_1d(surface.knots[1], degree_v, v, span_v)
    index_u = span_u - degree_u
    index_v = span_v - degree_v
    point = [0.0, 0.0, 0.0]
    tangent_u = [0.0, 0.0, 0.0]
    tangent_v = [0.0, 0.0, 0.0]
    for a in range(degree_u + 1):
        row = net[index_u + a]
        for b in range(degree_v + 1):
            control = row[index_v + b]
            weight = bu[a] * bv[b]
            point[0] += weight * control[0]
            point[1] += weight * control[1]
            point[2] += weight * control[2]
            weight = dbu[a] * bv[b]
            tangent_u[0] += weight * control[0]
            tangent_u[1] += weight * control[1]
            tangent_u[2] += weight * control[2]
            weight = bu[a] * dbv[b]
            tangent_v[0] += weight * control[0]
            tangent_v[1] += weight * control[1]
            tangent_v[2] += weight * control[2]
    return (
        (point[0], point[1], point[2]),
        (tangent_u[0], tangent_u[1], tangent_u[2]),
        (tangent_v[0], tangent_v[1], tangent_v[2]),
    )


def evaluate_fast(surface: NURBSSurface, u: float, v: float) -> Vec3:
    """``S(u, v)``: fast path for non-rational patches, reference path otherwise."""
    if surface.rational:
        return surface.evaluate(u, v)
    return eval_tensor(
        surface,
        u,
        v,
        find_span(surface.knots[0], surface.degree[0], u),
        find_span(surface.knots[1], surface.degree[1], v),
    )[0]


def jacobian_fast(surface: NURBSSurface, u: float, v: float) -> tuple[Vec3, Vec3, Vec3]:
    """``(S_u, S_v, S_u x S_v)`` on the fast path when the patch is non-rational."""
    if surface.rational:
        return surface.jacobian(u, v)
    _, tangent_u, tangent_v = eval_tensor(
        surface,
        u,
        v,
        find_span(surface.knots[0], surface.degree[0], u),
        find_span(surface.knots[1], surface.degree[1], v),
    )
    normal = (
        tangent_u[1] * tangent_v[2] - tangent_u[2] * tangent_v[1],
        tangent_u[2] * tangent_v[0] - tangent_u[0] * tangent_v[2],
        tangent_u[0] * tangent_v[1] - tangent_u[1] * tangent_v[0],
    )
    return tangent_u, tangent_v, normal


def cross_norm(normal: Sequence[float]) -> float:
    return math.sqrt(normal[0] ** 2 + normal[1] ** 2 + normal[2] ** 2)
