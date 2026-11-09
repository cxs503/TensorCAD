"""Differential geometry of a NURBS patch: normals, fundamental forms, curvature.

Second derivatives are the crux. :meth:`NURBSSurface.derivatives` differentiates the
*control net*, so a direction whose interior knot multiplicity equals the degree (the
three-arc rational circle of :func:`tensorcad.nurbs.sphere_surface`, for example) has no
differentiable second derivative in that direction -- and rightly refuses to invent one.
The standard remedy is Bezier extraction: insert knots until the span of interest is a
single Bezier segment and differentiate that segment. Knot insertion changes the control
net without changing the geometry or the parametrisation, so the extracted segment
evaluates to the same point at the same parameter and its derivatives are exact.
:func:`restrict_patch` implements the extraction and :func:`differential_geometry` uses
it, so ``S_uu``, ``S_uv`` and ``S_vv`` here are analytic rather than finite differences.

That matters for the accuracy claim this module exists to support. A finite-difference
"analytic" curvature would make the finite-difference cross-check in
:func:`finite_difference_geometry` circular; with exact extraction the two routes are
genuinely independent (one differentiates a rebuilt control net, the other differentiates
the tangent field numerically), and their agreement is evidence that both are right.

Named limitations
-----------------
* Only axis-aligned, tensor-product patches: no trimmed or multi-patch geometry, so no
  curvature is averaged across a patch seam and a multi-patch hull must be handled patch
  by patch.
* A direction of degree 1 has no analytic second derivative through this route; rather
  than return a zero that would silently be right in one parameterisation and wrong in a
  rational one, :func:`differential_geometry` raises and points at
  :func:`finite_difference_geometry`.
* Degenerate (pole) rows are excluded from :func:`curvature_scan` for the obvious reason:
  the parametrisation is singular there and curvature is undefined, not infinite.
* :func:`ellipsoid_principal_curvatures` is the closed-form truth for an axis-aligned
  ellipsoid only. It is deliberately independent of the NURBS route (it differentiates
  the implicit function), which is what makes the comparison meaningful.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

from ..nurbs import NURBSError, NURBSSurface, find_span, sphere_surface
from .elliptic import triaxial_ellipsoid_area
from .quadrature import knot_spans

Vec3 = tuple[float, float, float]
Vec2 = tuple[float, float]

__all__ = [
    "SurfacePoint",
    "curvature_scan",
    "differential_geometry",
    "ellipsoid_normal",
    "ellipsoid_principal_curvatures",
    "ellipsoid_surface",
    "finite_difference_geometry",
    "half_ellipsoid_patch",
    "hydrostatics_reference",
    "insert_knot",
    "knot_multiplicity",
    "restrict_patch",
]


# ----------------------------------------------------------------------------------
# small vector algebra (kept local so the module has no runtime dependency)
# ----------------------------------------------------------------------------------
def _dot(a: Sequence[float], b: Sequence[float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _sub(a: Sequence[float], b: Sequence[float]) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _cross(a: Sequence[float], b: Sequence[float]) -> Vec3:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _norm(a: Sequence[float]) -> float:
    return math.sqrt(_dot(a, a))


def _scale(a: Sequence[float], factor: float) -> Vec3:
    return (a[0] * factor, a[1] * factor, a[2] * factor)


# ----------------------------------------------------------------------------------
# exact Bezier extraction
# ----------------------------------------------------------------------------------
def knot_multiplicity(knots: Sequence[float], value: float) -> int:
    """How many times ``value`` occurs in the knot vector (exact float comparison)."""
    return sum(1 for knot in knots if knot == value)


def _homogeneous_curves(surface: NURBSSurface, axis: int) -> list[list[tuple[float, ...]]]:
    """Homogeneous control polygons along ``axis``, one per line of the other index."""
    rows = len(surface.control_net)
    columns = len(surface.control_net[0])
    curves: list[list[tuple[float, ...]]] = []
    if axis == 0:
        for j in range(columns):
            curves.append([
                (
                    surface.weights[i][j] * surface.control_net[i][j][0],
                    surface.weights[i][j] * surface.control_net[i][j][1],
                    surface.weights[i][j] * surface.control_net[i][j][2],
                    surface.weights[i][j],
                )
                for i in range(rows)
            ])
    else:
        for i in range(rows):
            curves.append([
                (
                    surface.weights[i][j] * surface.control_net[i][j][0],
                    surface.weights[i][j] * surface.control_net[i][j][1],
                    surface.weights[i][j] * surface.control_net[i][j][2],
                    surface.weights[i][j],
                )
                for j in range(columns)
            ])
    return curves


def insert_knot(
    surface: NURBSSurface, axis: int, value: float, times: int = 1
) -> NURBSSurface:
    """Insert ``value`` ``times`` times into the ``axis`` knot vector (Boehm A5.1).

    Knot insertion preserves both geometry and parametrisation: the returned surface
    evaluates to the same point at the same parameter. It is the enabling operation for
    exact Bezier extraction, and :func:`restrict_patch` uses it that way.
    """
    if axis not in (0, 1):
        raise ValueError(f"axis must be 0 or 1, got {axis!r}")
    if times < 0:
        raise ValueError(f"times must be non-negative, got {times!r}")
    if times == 0:
        return surface
    degree = surface.degree[axis]
    knots = list(surface.knots[axis])
    low = knots[degree]
    high = knots[len(knots) - degree - 1]
    if not (low <= value <= high):
        raise NURBSError(f"knot {value!r} is outside the domain [{low!r}, {high!r}]")

    curves = _homogeneous_curves(surface, axis)
    for _ in range(times):
        span = find_span(tuple(knots), degree, value)
        updated: list[list[tuple[float, ...]]] = []
        for points in curves:
            head = points[: span - degree + 1]
            blended: list[tuple[float, ...]] = []
            for index in range(span - degree + 1, span + 1):
                denominator = knots[index + degree] - knots[index]
                alpha = 0.0 if denominator == 0.0 else (value - knots[index]) / denominator
                left = points[index]
                right = points[index - 1]
                blended.append(
                    tuple(alpha * left[k] + (1.0 - alpha) * right[k] for k in range(4))
                )
            updated.append(head + blended + points[span:])
        knots.insert(span + 1, value)
        curves = updated

    length = len(curves[0])
    count = len(curves)
    if axis == 0:
        net = tuple(
            tuple(
                (
                    curves[j][i][0] / curves[j][i][3],
                    curves[j][i][1] / curves[j][i][3],
                    curves[j][i][2] / curves[j][i][3],
                )
                for j in range(count)
            )
            for i in range(length)
        )
        weights = tuple(tuple(curves[j][i][3] for j in range(count)) for i in range(length))
        knot_pair = (tuple(knots), surface.knots[1])
    else:
        net = tuple(
            tuple(
                (
                    curves[i][j][0] / curves[i][j][3],
                    curves[i][j][1] / curves[i][j][3],
                    curves[i][j][2] / curves[i][j][3],
                )
                for j in range(length)
            )
            for i in range(count)
        )
        weights = tuple(tuple(curves[i][j][3] for j in range(length)) for i in range(count))
        knot_pair = (surface.knots[0], tuple(knots))
    return NURBSSurface(
        degree=surface.degree, knots=knot_pair, control_net=net, weights=weights
    )


def split_patch(
    surface: NURBSSurface, axis: int, value: float
) -> tuple[NURBSSurface, NURBSSurface]:
    """Split the patch at an isoparametric ``value`` into a left and a right part.

    The two parts share the boundary control point row, their parameters are unchanged, and
    each is a clamped B-spline surface whose domain is one side of ``value``. This is the
    primitive behind :func:`restrict_patch`: splitting at both bounds and keeping the
    middle piece gives an exact sub-patch without any index arithmetic that would have to
    be trusted.
    """
    if axis not in (0, 1):
        raise ValueError(f"axis must be 0 or 1, got {axis!r}")
    degree = surface.degree[axis]
    knots_of = surface.knots[axis]
    low = knots_of[degree]
    high = knots_of[len(knots_of) - degree - 1]
    if not (low < value < high):
        raise NURBSError(
            f"split value {value!r} must be strictly inside the domain [{low!r}, {high!r}]"
        )
    work = surface
    multiplicity = knot_multiplicity(work.knots[axis], value)
    if multiplicity < degree:
        work = insert_knot(work, axis, value, degree - multiplicity)
    knots = list(work.knots[axis])
    span = find_span(tuple(knots), degree, value)
    total = len(work.control_net) if axis == 0 else len(work.control_net[0])
    left_count = span - degree + 1
    right_start = span - degree
    if left_count < 1 or right_start >= total:
        raise NURBSError(f"splitting at {value!r} produced an empty part")
    left_knots = [knot for knot in knots[:span] if knot != value] + [value] * (degree + 1)
    right_knots = [value] * (degree + 1) + [knot for knot in knots[span + 1:] if knot != value]
    if axis == 0:
        left = NURBSSurface(
            degree=work.degree,
            knots=(tuple(left_knots), work.knots[1]),
            control_net=tuple(work.control_net[i] for i in range(left_count)),
            weights=tuple(work.weights[i] for i in range(left_count)),
        )
        right = NURBSSurface(
            degree=work.degree,
            knots=(tuple(right_knots), work.knots[1]),
            control_net=tuple(work.control_net[i] for i in range(right_start, total)),
            weights=tuple(work.weights[i] for i in range(right_start, total)),
        )
    else:
        left = NURBSSurface(
            degree=work.degree,
            knots=(work.knots[0], tuple(left_knots)),
            control_net=tuple(row[:left_count] for row in work.control_net),
            weights=tuple(row[:left_count] for row in work.weights),
        )
        right = NURBSSurface(
            degree=work.degree,
            knots=(work.knots[0], tuple(right_knots)),
            control_net=tuple(row[right_start:] for row in work.control_net),
            weights=tuple(row[right_start:] for row in work.weights),
        )
    return left, right


def restrict_patch(surface: NURBSSurface, u_bounds: Vec2, v_bounds: Vec2) -> NURBSSurface:
    """The sub-patch over ``u_bounds x v_bounds``, exactly and with the same parameters.

    The bounds need not be knots. The restriction is obtained by splitting at each bound
    and keeping the middle piece, so ``S_restricted(u, v)`` equals ``S(u, v)`` to round-off
    for every parameter inside the restriction -- which the tests check directly rather
    than infer from the knot bookkeeping.
    """
    work = surface
    for axis, (lower, upper) in ((0, u_bounds), (1, v_bounds)):
        knots = work.knots[axis]
        degree = work.degree[axis]
        domain_low = knots[degree]
        domain_high = knots[len(knots) - degree - 1]
        if not (domain_low <= lower < upper <= domain_high):
            raise NURBSError(
                f"bounds {(lower, upper)!r} are not a positive interval inside "
                f"[{domain_low!r}, {domain_high!r}]"
            )
        if upper < domain_high:
            work, _ = split_patch(work, axis, upper)
        if lower > domain_low:
            _, work = split_patch(work, axis, lower)
    return work


def _span_of(knots: Sequence[float], degree: int, value: float) -> Vec2:
    for left, right in knot_spans(knots, degree):
        if left <= value <= right:
            return (left, right)
    raise NURBSError(f"parameter {value!r} is not in any span")


def span_patch(surface: NURBSSurface, u: float, v: float) -> NURBSSurface:
    """The single Bezier segment containing ``(u, v)``, as its own surface."""
    return restrict_patch(
        surface,
        _span_of(surface.knots[0], surface.degree[0], u),
        _span_of(surface.knots[1], surface.degree[1], v),
    )


# ----------------------------------------------------------------------------------
# fundamental forms
# ----------------------------------------------------------------------------------
@dataclass(frozen=True)
class SurfacePoint:
    """Everything a solver needs to know about the surface at one parameter pair."""

    u: float
    v: float
    point: Vec3
    du: Vec3
    dv: Vec3
    duu: Vec3
    duv: Vec3
    dvv: Vec3
    normal: Vec3
    area_element: float
    first_form: Vec3  # (E, F, G)
    second_form: Vec3  # (L, M, N)
    mean_curvature: float
    gauss_curvature: float
    principal_curvatures: tuple[float, float]
    second_derivatives: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "u": self.u,
            "v": self.v,
            "point": list(self.point),
            "normal": list(self.normal),
            "area_element": self.area_element,
            "first_form": list(self.first_form),
            "second_form": list(self.second_form),
            "mean_curvature": self.mean_curvature,
            "gauss_curvature": self.gauss_curvature,
            "principal_curvatures": list(self.principal_curvatures),
            "second_derivatives": self.second_derivatives,
        }


def _assemble(
    u: float,
    v: float,
    point: Vec3,
    du: Vec3,
    dv: Vec3,
    duu: Vec3,
    duv: Vec3,
    dvv: Vec3,
    method: str,
) -> SurfacePoint:
    cross = _cross(du, dv)
    jacobian = _norm(cross)
    if jacobian == 0.0:
        raise NURBSError(f"degenerate parametrisation at (u, v)=({u!r}, {v!r})")
    normal = _scale(cross, 1.0 / jacobian)
    e_form = _dot(du, du)
    f_form = _dot(du, dv)
    g_form = _dot(dv, dv)
    metric = e_form * g_form - f_form * f_form
    if metric <= 0.0:
        raise NURBSError(f"degenerate metric at (u, v)=({u!r}, {v!r})")
    l_form = _dot(duu, normal)
    m_form = _dot(duv, normal)
    n_form = _dot(dvv, normal)
    mean = (e_form * n_form - 2.0 * f_form * m_form + g_form * l_form) / (2.0 * metric)
    gauss = (l_form * n_form - m_form * m_form) / metric
    root = math.sqrt(max(mean * mean - gauss, 0.0))
    return SurfacePoint(
        u=u,
        v=v,
        point=point,
        du=du,
        dv=dv,
        duu=duu,
        duv=duv,
        dvv=dvv,
        normal=normal,
        area_element=jacobian,
        first_form=(e_form, f_form, g_form),
        second_form=(l_form, m_form, n_form),
        mean_curvature=mean,
        gauss_curvature=gauss,
        principal_curvatures=(mean + root, mean - root),
        second_derivatives=method,
    )


def differential_geometry(surface: NURBSSurface, u: float, v: float) -> SurfacePoint:
    """Analytic first and second fundamental forms at ``(u, v)``.

    The second derivatives come from the Bezier segment containing the parameter. Both
    degrees must be at least 2; a linear direction has no second derivative through the
    knot-insertion route and the caller is sent to :func:`finite_difference_geometry`
    instead of being handed a zero that may be an artefact of a non-rational
    parameterisation.
    """
    local = span_patch(surface, u, v)
    p, q = local.degree
    if p < 2 or q < 2:
        raise NURBSError(
            f"analytic second derivatives need degree >= 2 in both directions, got {p}x{q}"
            "; use finite_difference_geometry for a linear patch"
        )
    derivatives = local.derivatives(u, v, 2)
    return _assemble(
        u,
        v,
        derivatives[(0, 0)],
        derivatives[(1, 0)],
        derivatives[(0, 1)],
        derivatives[(2, 0)],
        derivatives[(1, 1)],
        derivatives[(0, 2)],
        "analytic",
    )


def finite_difference_geometry(
    surface: NURBSSurface, u: float, v: float, *, step: float = 1e-5
) -> SurfacePoint:
    """The same measures with the second derivatives taken by central differences.

    Position, first derivatives and therefore the normal and ``|J|`` stay exact -- only
    the curvature route is numerical -- so agreement with
    :func:`differential_geometry` tests the second-derivative supply chain rather than the
    whole module against itself.
    """
    (u_low, u_high), (v_low, v_high) = surface.domain
    h = min(step, 0.25 * (u - u_low), 0.25 * (u_high - u), 0.25 * (v - v_low), 0.25 * (v_high - v))
    if h <= 0.0:
        raise NURBSError(
            f"finite differences need an interior parameter, got (u, v)=({u!r}, {v!r})"
        )
    base = surface.derivatives(u, v, 1)
    du_plus = surface.derivatives(u + h, v, 1)[(1, 0)]
    du_minus = surface.derivatives(u - h, v, 1)[(1, 0)]
    du_v_plus = surface.derivatives(u, v + h, 1)[(1, 0)]
    du_v_minus = surface.derivatives(u, v - h, 1)[(1, 0)]
    dv_plus = surface.derivatives(u, v + h, 1)[(0, 1)]
    dv_minus = surface.derivatives(u, v - h, 1)[(0, 1)]
    return _assemble(
        u,
        v,
        base[(0, 0)],
        base[(1, 0)],
        base[(0, 1)],
        _scale(_sub(du_plus, du_minus), 1.0 / (2.0 * h)),
        _scale(_sub(du_v_plus, du_v_minus), 1.0 / (2.0 * h)),
        _scale(_sub(dv_plus, dv_minus), 1.0 / (2.0 * h)),
        "finite_difference",
    )


def curvature_scan(
    surface: NURBSSurface,
    *,
    samples: Sequence[Vec2] = (),
    grid: int = 1,
    u_margin: float = 0.02,
    v_margin: float = 0.02,
) -> list[SurfacePoint]:
    """Analytic curvature on a parameter grid (or an explicit sample list).

    Margins keep the scan strictly inside the domain so that a pole or a degenerate
    boundary row -- where curvature is undefined rather than large -- is never sampled.
    """
    (u_low, u_high), (v_low, v_high) = surface.domain
    u_span = u_high - u_low
    v_span = v_high - v_low
    parameters: list[Vec2] = list(samples)
    if grid > 1:
        for i in range(grid):
            u = u_low + u_span * (u_margin + (1.0 - 2.0 * u_margin) * i / (grid - 1))
            for j in range(grid):
                v = v_low + v_span * (v_margin + (1.0 - 2.0 * v_margin) * j / (grid - 1))
                parameters.append((u, v))
    if not parameters:
        raise ValueError("curvature_scan needs explicit samples or grid > 1")
    return [differential_geometry(surface, u, v) for u, v in parameters]


# ----------------------------------------------------------------------------------
# the exact ellipsoid: NURBS representation and closed-form curvature
# ----------------------------------------------------------------------------------
def ellipsoid_surface(a: float, b: float, c: float, *, arcs: int = 3) -> NURBSSurface:
    """An exact axis-aligned ellipsoid with semi-axes ``a, b, c``, as one rational patch.

    Scaling the control points of a rational sphere scales the surface: the quotient
    ``sum(N w P) / sum(N w)`` is linear in ``P``, so applying ``diag(a, b, c)`` to every
    control point applies it to the surface and leaves the weights alone. The result is
    the true ellipsoid (``(x/a)^2 + (y/b)^2 + (z/c)^2 = 1`` to round-off), not a fit.
    """
    for name, value in (("a", a), ("b", b), ("c", c)):
        if value <= 0.0:
            raise NURBSError(f"semi-axis {name} must be positive, got {value!r}")
    base = sphere_surface(1.0, arcs=arcs)
    net = tuple(
        tuple((a * point[0], b * point[1], c * point[2]) for point in row)
        for row in base.control_net
    )
    return NURBSSurface(
        degree=base.degree, knots=base.knots, control_net=net, weights=base.weights
    )


def half_ellipsoid_patch(
    a: float, b: float, c: float, *, arcs: int = 3
) -> NURBSSurface:
    """The ``z >= 0`` half of the ellipsoid: a hull-shaped patch with the deck plane open.

    ``v`` runs from the south pole (0) through the equator (1) to the north pole (2), so
    the upper half is ``v in [1, 2]``. Its area is half the closed area by symmetry and
    its divergence-theorem volume is ``2/3 pi a b c``.
    """
    return restrict_patch(ellipsoid_surface(a, b, c, arcs=arcs), (0.0, 1.0), (1.0, 2.0))


def ellipsoid_normal(a: float, b: float, c: float, point: Sequence[float]) -> Vec3:
    """Outward unit normal of the ellipsoid at ``point`` (gradient of the implicit form)."""
    gradient = (
        point[0] / (a * a),
        point[1] / (b * b),
        point[2] / (c * c),
    )
    length = _norm(gradient)
    if length == 0.0:
        raise NURBSError("ellipsoid normal is undefined at the centre")
    return _scale(gradient, 1.0 / length)


def ellipsoid_principal_curvatures(
    a: float, b: float, c: float, point: Sequence[float]
) -> tuple[float, float]:
    """Closed-form principal curvatures of the ellipsoid at ``point``.

    For ``F = sum (x_i / a_i)^2 - 1`` the shape operator on the tangent plane is the
    Hessian of ``F`` restricted to that plane and divided by ``|grad F|``. With an
    orthonormal tangent basis the 2x2 second fundamental form is symmetric, so its
    eigenvalues are available in closed form -- no iteration, no finite differences, and
    no reuse of the NURBS machinery whose answer it is checking.
    """
    normal = ellipsoid_normal(a, b, c, point)
    gradient_length = 2.0 * math.sqrt(
        (point[0] / (a * a)) ** 2 + (point[1] / (b * b)) ** 2 + (point[2] / (c * c)) ** 2
    )
    # Orthonormal tangent basis: take the axis least aligned with the normal.
    axes = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
    reference = min(axes, key=lambda axis: abs(_dot(axis, normal)))
    t1 = _cross(normal, reference)
    t1 = _scale(t1, 1.0 / _norm(t1))
    t2 = _cross(normal, t1)
    # II(t_i, t_j) = t_i^T Hess F t_j / |grad F|, Hess F = 2 diag(1/a^2, 1/b^2, 1/c^2).
    curvatures = (2.0 / (a * a), 2.0 / (b * b), 2.0 / (c * c))
    b11 = sum(t1[k] * t1[k] * curvatures[k] for k in range(3)) / gradient_length
    b22 = sum(t2[k] * t2[k] * curvatures[k] for k in range(3)) / gradient_length
    b12 = sum(t1[k] * t2[k] * curvatures[k] for k in range(3)) / gradient_length
    half_trace = 0.5 * (b11 + b22)
    root = math.sqrt(max(0.25 * (b11 - b22) ** 2 + b12 * b12, 0.0))
    return (half_trace + root, half_trace - root)


def hydrostatics_reference(a: float, b: float, c: float) -> dict[str, Any]:
    """Closed-form hydrostatics of the half-ellipsoid hull (analytic reference for Cb/S).

    Hull: ``L = 2a``, ``B = 2b``, draft ``T = c``, waterplane ``z = 0``. Volume
    ``2/3 pi a b c``, waterplane area ``pi a b``, wetted surface half the ellipsoid area
    (the deck is not wetted). ``Cb = pi/6`` and ``Cp = 1/3`` follow from the volume
    integrals, and the wetted surface is the one quantity that needs a real elliptic
    integral -- which is why the cross-check in
    ``benchmarks/verified/cad_highorder_geometry`` is worth running.
    """
    volume = 2.0 / 3.0 * math.pi * a * b * c
    length = 2.0 * a
    breadth = 2.0 * b
    draft = c
    waterplane = math.pi * a * b
    midship = math.pi * b * c
    total_surface = triaxial_ellipsoid_area(a, b, c)
    return {
        "L": length,
        "B": breadth,
        "T": draft,
        "volume": volume,
        "Cb": volume / (length * breadth * draft),
        "Cp": volume / (midship * length),
        "Cwp": waterplane / (length * breadth),
        "waterplane_area": waterplane,
        "wetted_surface": 0.5 * total_surface,
        "total_surface": total_surface,
    }
