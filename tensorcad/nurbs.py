"""Rational B-spline (NURBS) curves and surfaces: exact arcs, area, tessellation.

Why this module exists
----------------------
The geometry service had sampled hulls, implicit fields and triangle meshes, but no
*parametric* surface of record. A parametric surface is the only representation that
carries its own analytic structure: a circular arc is exact (not a 96-gon), the
lateral area is ``int int |S_u x S_v| du dv`` (not a sum over triangles), and the
continuity of the parametrisation across a knot is a property of the *data*, not an
artefact of the sampling density. This module is the minimal such surface: rational
quadratic arcs, tensor-product cylinders, exact derivatives, parameter-derivative
continuity checks, and a tessellator whose faces all carry a consistent, outward
normal.

Coordinate frame (identical to the rest of the package): ``x`` longitudinal, ``y``
transverse starboard-positive, ``z`` vertical up.

What is checked, and how
------------------------
* **Analytic anchors.** An arc of radius ``r`` and sweep ``theta`` has length
  ``r * theta``; a cylinder of radius ``r`` and height ``h`` has lateral area
  ``2 pi r h`` and volume ``pi r^2 h``. :meth:`NURBSSurface.area` integrates with
  Gauss-Legendre *per knot span*, so the parametric area agrees with the closed form
  to round-off; the tessellated area converges to it as the grid is refined, and that
  convergence is reported rather than assumed.
* **Parameter-derivative continuity.** At every interior knot the derivatives of
  order ``k`` are evaluated from the span on the left and from the span on the right;
  their difference is the jump. The declared continuity of a B-spline is ``p - m``
  (degree minus knot multiplicity). The measured jump says whether the implementation
  delivers it, and for rational arcs it exposes the classical result that the
  three-arc circle is ``C1`` in parameter although its multiplicities only promise
  ``C0``.
* **Normal consistency.** Every tessellated face is compared with the analytic
  outward normal at the face's parametric centre: no face may point inward and the
  angle deviation must decay with the grid. A closed tessellation must additionally
  be watertight (every edge shared by exactly two faces), which is what makes an STL
  export usable by a solver.

Only the Python standard library is required: the exactness claims are arithmetic
claims, not library claims, and no CAD kernel should be a prerequisite for auditing
them.

Named limitations
-----------------
* A single quadratic arc cannot sweep more than ``pi``; :func:`circle_curve` stitches
  ``arcs`` of them (three by default) and the joins are only ``C0`` in the knot
  bookkeeping (measured ``C1``, see above).
* A derivative of order ``p - m + 1`` across a knot of multiplicity ``m`` does not
  exist as a function (it is a distribution there), so :meth:`NURBSCurve.derivatives`
  and :meth:`NURBSSurface.derivatives` **raise** for it instead of returning a
  plausible-looking number; :func:`curve_derivative_continuity` clamps the order it
  measures to what is representable and reports ``max_measurable_order``.
* The surface type is a tensor-product patch with a rectangular control net: no
  trimmed patches, no boolean operations, no non-tensor-product topology. A hull with
  a bulbous bow would need a trimmed multi-patch model, which is out of scope here.
* A *closed* patch with singular pole rows is supported: :func:`sphere_surface` is one
  patch whose two pole rings pinch to a vertex, and its tessellation alone is a
  watertight, outward, positive-volume shell (anchored in
  ``docs/nurbs_pole_surface.md``). What is still absent is topology beyond that patch:
  no multi-patch stitching and no trimmed seam.
* The PC6 compatibility block carries the *representation* and the L/B/Cb gate; the
  PC6 hull surface itself is still lofted by :mod:`tensorcad.pc6` from analytic
  stations, not re-fitted as one NURBS patch.
"""
from __future__ import annotations

import json
import math
from collections.abc import Sequence
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path
from typing import Any

from .geometry import Mesh, MeshBuilder, Vec3

#: Schema tag of the additive parametric-surface block (see :func:`pc6_compat_block`).
NURBS_SCHEMA = "tensorcad.nurbs.surface/1.0"

TWO_PI = 2.0 * math.pi

#: ``cos(45 deg)``: the middle weight of an exact 90-degree rational quadratic arc.
SQRT_HALF = math.sqrt(0.5)

__all__ = [
    "ContinuityRecord",
    "NURBSCurve",
    "NURBSError",
    "NURBSSurface",
    "Tessellation",
    "circle_curve",
    "circular_arc",
    "curve_derivative_continuity",
    "cylinder_solid_mesh",
    "cylinder_surface",
    "derivative_net",
    "find_span",
    "gauss_legendre",
    "mesh_solid_report",
    "pc6_compat_block",
    "spline_derivatives",
    "spline_point",
    "surface_derivative_continuity",
    "tensor_point",
    "tessellation_report",
    "validate_knots",
    "validate_weights",
]


class NURBSError(ValueError):
    """Raised when data is not a valid rational B-spline, or an anchor fails."""


# ---------------------------------------------------------------------------
# validation (the campaign's negative controls live here)
# ---------------------------------------------------------------------------
def _as_floats(values: Sequence[float], label: str) -> tuple[float, ...]:
    out = tuple(float(value) for value in values)
    for index, value in enumerate(out):
        if not math.isfinite(value):
            raise NURBSError(f"{label}[{index}] is not finite: {value!r}")
    return out


def validate_knots(
    knots: Sequence[float],
    degree: int,
    control_count: int | None = None,
    *,
    clamped: bool = True,
) -> tuple[float, ...]:
    """Validate a knot vector and return it as a tuple of floats.

    Rejected, each naming the offending index:

    * a non-integer or non-positive degree;
    * a wrong number of knots (a degree ``p`` spline with ``n + 1`` control points
      needs exactly ``n + p + 2`` knots);
    * a **non-monotone** knot vector -- knot values must be non-decreasing, so one
      decrease makes the basis functions undefined and every downstream integral
      meaningless (negative control 1);
    * a knot whose multiplicity exceeds ``p + 1``;
    * an empty parameter domain (``U[p] >= U[n+1]``);
    * (when ``clamped``) ends that are not ``p + 1``-fold, because an unclamped end
      produces a curve that misses its end control points and silently distorts the
      length or the area.
    """
    if not isinstance(degree, int) or isinstance(degree, bool) or degree < 1:
        raise NURBSError(f"degree must be an integer >= 1, got {degree!r}")
    values = _as_floats(knots, "knots")
    if control_count is not None:
        expected = control_count + degree + 1
        if len(values) != expected:
            raise NURBSError(
                f"a degree-{degree} spline with {control_count} control points needs "
                f"{expected} knots, got {len(values)}"
            )
    if len(values) < 2 * (degree + 1):
        raise NURBSError(f"knot vector of length {len(values)} is too short for degree {degree}")
    multiplicity = 1
    for index in range(1, len(values)):
        if values[index] < values[index - 1]:
            raise NURBSError(
                f"knot vector must be non-decreasing: knots[{index}]={values[index]!r} < "
                f"knots[{index - 1}]={values[index - 1]!r}"
            )
        multiplicity = multiplicity + 1 if values[index] == values[index - 1] else 1
        if multiplicity > degree + 1:
            raise NURBSError(
                f"knot value {values[index]!r} appears {multiplicity} times, more than "
                f"degree + 1 = {degree + 1}"
            )
    n = len(values) - degree - 2
    if not values[degree] < values[n + 1]:
        raise NURBSError(
            f"empty parameter domain: U[p]={values[degree]!r} is not below "
            f"U[n+1]={values[n + 1]!r}"
        )
    if clamped:
        for offset in range(degree):
            if values[offset] != values[degree]:
                raise NURBSError(
                    f"knot vector is not clamped at the start: knots[{offset}]="
                    f"{values[offset]!r} != knots[p]={values[degree]!r}"
                )
            index = n + 2 + offset
            if values[index] != values[n + 1]:
                raise NURBSError(
                    f"knot vector is not clamped at the end: knots[{index}]="
                    f"{values[index]!r} != knots[n+1]={values[n + 1]!r}"
                )
    return values


def validate_weights(
    weights: Sequence[float],
    control_count: int | None = None,
) -> tuple[float, ...]:
    """Validate homogeneous weights and return them as a tuple of floats.

    A rational B-spline needs **strictly positive** weights: a zero weight divides
    the projective point by zero, and a negative weight folds the patch inside out at
    some parameter value (the point at infinity enters the domain). Both are rejected
    here instead of producing a silently wrong area (negative control 2).
    """
    values = _as_floats(weights, "weights")
    if control_count is not None and len(values) != control_count:
        raise NURBSError(f"{len(values)} weights for {control_count} control points")
    for index, value in enumerate(values):
        if value <= 0.0:
            raise NURBSError(
                f"weights[{index}] = {value!r}: rational B-spline weights must be strictly "
                f"positive (zero divides by zero, a negative weight folds the patch)"
            )
    return values


# ---------------------------------------------------------------------------
# basis functions
# ---------------------------------------------------------------------------
def find_span(knots: Sequence[float], degree: int, u: float) -> int:
    """Index of the knot span containing ``u`` (KnotSpan, Piegl & Tiller A2.1).

    ``u`` must lie in the closed domain ``[U[p], U[n+1]]``. Outside it every basis
    function is zero and the surface would evaluate to ``0/0``, so an out-of-domain
    parameter raises instead of clamping silently.
    """
    n = len(knots) - degree - 2
    low, high = knots[degree], knots[n + 1]
    if u < low or u > high:
        raise NURBSError(f"parameter {u!r} is outside the domain [{low!r}, {high!r}]")
    if u == high:
        return n
    if u == low:
        return degree
    lo, hi = degree, n + 1
    mid = (lo + hi) // 2
    while u < knots[mid] or u >= knots[mid + 1]:
        if u < knots[mid]:
            hi = mid
        else:
            lo = mid
        mid = (lo + hi) // 2
    return mid


def _basis_values(span: int, u: float, degree: int, knots: Sequence[float]) -> list[float]:
    """The ``degree + 1`` nonzero basis functions on ``span`` (Piegl & Tiller A2.2)."""
    left = [0.0] * (degree + 1)
    right = [0.0] * (degree + 1)
    values = [1.0] + [0.0] * degree
    for j in range(1, degree + 1):
        left[j] = u - knots[span + 1 - j]
        right[j] = knots[span + j] - u
        saved = 0.0
        for r in range(j):
            denominator = right[r + 1] + left[j - r]
            temp = values[r] / denominator if denominator != 0.0 else 0.0
            values[r] = saved + right[r + 1] * temp
            saved = left[j - r] * temp
        values[j] = saved
    return values


def derivative_net(
    points: Sequence[Sequence[float]],
    degree: int,
    knots: Sequence[float],
) -> tuple[tuple[tuple[float, ...], ...], int, tuple[float, ...]]:
    """Control net of the derivative curve/spline (Piegl & Tiller A3.4).

    Differentiating a B-spline is exact and stays inside the representation::

        D_i = p (P_{i+1} - P_i) / (U_{i+p+1} - U_{i+1}),
        knots(D) = U[1:-1], degree(D) = p - 1

    Applying this ``k`` times gives the ``k``-th derivative spline. Building the
    derivative *net* (instead of differentiating the basis functions in place) keeps
    every step a convex combination of control points with a positive denominator,
    which is why it cannot produce the silent zeros a hand-rolled derivative table
    can. The test suite checks it against central differences.
    """
    if degree < 1:
        raise NURBSError(f"cannot differentiate a degree-{degree} spline")
    width = len(points[0])
    out: list[tuple[float, ...]] = []
    for i in range(len(points) - 1):
        denominator = knots[i + degree + 1] - knots[i + 1]
        if denominator == 0.0:
            raise NURBSError(
                f"zero-length knot interval U[{i + degree + 1}] - U[{i + 1}] while "
                f"differentiating"
            )
        factor = degree / denominator
        out.append(
            tuple(factor * (points[i + 1][c] - points[i][c]) for c in range(width))
        )
    return tuple(out), degree - 1, tuple(knots[1:-1])


def spline_point(
    points: Sequence[Sequence[float]],
    degree: int,
    knots: Sequence[float],
    u: float,
    *,
    span: int | None = None,
) -> tuple[float, ...]:
    """Evaluate a (vector-valued) B-spline at ``u``; ``span`` selects the branch."""
    width = len(points[0])
    if degree < 0:
        return tuple(0.0 for _ in range(width))
    used = find_span(knots, max(degree, 0), u) if span is None else span
    values = _basis_values(used, u, degree, knots)
    base = used - degree
    return tuple(
        sum(values[i] * points[base + i][c] for i in range(degree + 1)) for c in range(width)
    )


def _dehomogenise(point: Sequence[float]) -> tuple[float, float, float]:
    """3-D point of a homogeneous control point; a zero weight is an error, not a NaN."""
    weight = point[3]
    if weight == 0.0:
        raise NURBSError("degree elevation produced a zero weight")
    return (point[0] / weight, point[1] / weight, point[2] / weight)


def _uniform_weight(points: Sequence[Sequence[float]]) -> float | None:
    """The common weight if every homogeneous point carries the same one, else ``None``.

    A net whose weights are all equal represents the same point set as the
    non-rational net with unit weights (scaling every homogeneous coordinate by a
    constant leaves the rational quotient unchanged). The elevation solve returns
    that constant as ``w +/- 4e-16``, which would silently switch the result onto
    the rational evaluation path -- so the caller restores the exact constant.
    """
    first = points[0][3]
    return first if all(point[3] == first for point in points) else None


def _rescale_to_uniform(points: Sequence[Sequence[float]], weight: float) -> list[list[float]]:
    """Rescale homogeneous points so every one carries exactly ``weight``.

    Homogeneous scaling (all four components by a common factor) does not move the
    point, so this restores the exact constant that the solve returned as
    ``weight +/- 4e-16`` without touching the geometry.
    """
    out: list[list[float]] = []
    for point in points:
        factor = weight / point[3]
        out.append([value * factor for value in point])
    return out


def _elevate_homogeneous(
    points: Sequence[Sequence[float]],
    degree: int,
    knots: Sequence[float],
    times: int,
    *,
    label: str = "curve",
) -> tuple[int, tuple[float, ...], tuple[tuple[float, ...], ...]]:
    """Degree-elevate a **homogeneous** control polygon; shared by curves and surfaces.

    ``points`` are 4-component homogeneous control points ``(w x, w y, w z, w)``.
    The returned net spans the same geometry in the degree ``degree + times``
    spline space, with every knot multiplicity grown by ``times`` (ISO 10303-42
    / Piegl & Tiller A5.9), so the original continuity is preserved and the
    representation grows by the minimum ``times * (#spans)`` points.

    Working in homogeneous coordinates is what makes this exact for rational
    geometry: elevating the 3-D net and carrying the old weights across is a
    different, generally wrong curve. Collocating at the Greville abscissae of
    the new knot vector gives a banded, totally positive system, so the
    elimination needs no pivoting.
    """
    low = knots[degree]
    high = knots[len(points)]
    target = degree + times
    unique = sorted({value for value in knots if low < value < high})
    multiplicity = {value: sum(1 for knot in knots if knot == value) for value in unique}
    disconnected = [value for value in unique if multiplicity[value] > degree]
    if disconnected:
        # A multiplicity of degree + 1 is not a C0 join but a break: the two sides
        # share no control point, so the single-patch collocation matrix is
        # singular by construction. Refuse with the fix instead of a generic
        # singular-matrix message.
        raise NURBSError(
            f"degree elevation needs a C0-or-better {label}; knot "
            f"{disconnected[0]!r} has multiplicity {multiplicity[disconnected[0]]} "
            f"> degree {degree}, so the {label} is discontinuous there -- for a curve "
            "use NURBSCurve.split_at_discontinuities(), or "
            "NURBSCurve.elevate_degree_piecewise(times) which splits and elevates in "
            "one step"
        )

    elevated_knots: list[float] = [low] * (target + 1)
    for value in unique:
        elevated_knots.extend([value] * (multiplicity[value] + times))
    elevated_knots.extend([high] * (target + 1))
    count = len(elevated_knots) - target - 1
    greville = [sum(elevated_knots[i + 1 : i + target + 1]) / target for i in range(count)]

    values = [spline_point(points, degree, knots, u) for u in greville]

    matrix = [[0.0] * count for _ in range(count)]
    for row, u in enumerate(greville):
        span = find_span(elevated_knots, target, u)
        basis = _basis_values(span, u, target, elevated_knots)
        for offset, value in enumerate(basis):
            matrix[row][span - target + offset] = value

    # Solve the 4 homogeneous components with one elimination: same matrix, four
    # right-hand sides. The collocation matrix at Greville points is totally
    # positive, hence every leading principal minor is positive -- no pivoting.
    solutions = [[values[row][axis] for row in range(count)] for axis in range(4)]
    for pivot in range(count):
        diagonal = matrix[pivot][pivot]
        if diagonal == 0.0:
            raise NURBSError(
                "degree elevation failed: the Greville collocation matrix is singular"
            )
        for row in range(pivot + 1, count):
            factor = matrix[row][pivot]
            if factor == 0.0:
                continue
            factor /= diagonal
            for column in range(pivot, count):
                matrix[row][column] -= factor * matrix[pivot][column]
            for axis in range(4):
                solutions[axis][row] -= factor * solutions[axis][pivot]
    coefficients = [[0.0] * count for _ in range(4)]
    for pivot in range(count - 1, -1, -1):
        for axis in range(4):
            total = solutions[axis][pivot]
            for column in range(pivot + 1, count):
                total -= matrix[pivot][column] * coefficients[axis][column]
            coefficients[axis][pivot] = total / matrix[pivot][pivot]

    net = tuple(
        (
            coefficients[0][index],
            coefficients[1][index],
            coefficients[2][index],
            coefficients[3][index],
        )
        for index in range(count)
    )
    return target, tuple(elevated_knots), net


def spline_derivatives(
    points: Sequence[Sequence[float]],
    degree: int,
    knots: Sequence[float],
    u: float,
    order: int,
    *,
    span: int | None = None,
) -> tuple[tuple[float, ...], ...]:
    """Derivatives ``0..order`` of a B-spline at ``u`` (exact, via derivative nets)."""
    if order < 0:
        raise NURBSError(f"derivative order must be >= 0, got {order}")
    out: list[tuple[float, ...]] = []
    net = tuple(tuple(point) for point in points)
    deg = degree
    kts = tuple(knots)
    used = span
    for k in range(order + 1):
        if k > 0:
            net, deg, kts = derivative_net(net, deg, kts)
            if used is not None:
                used -= 1
        out.append(spline_point(net, deg, kts, u, span=used))
    return tuple(out)


def _differentiate_axis(
    net: Sequence[Sequence[Sequence[float]]],
    axis: int,
    degree: int,
    knots: Sequence[float],
) -> tuple[tuple[tuple[tuple[float, ...], ...], ...], int, tuple[float, ...]]:
    """Differentiate a control net along one axis (``0`` = rows/u, ``1`` = columns/v)."""
    if axis == 0:
        columns = list(zip(*net, strict=True))
        reduced = [derivative_net(column, degree, knots)[0] for column in columns]
        return tuple(zip(*reduced, strict=True)), degree - 1, tuple(knots[1:-1])
    rows = [derivative_net(row, degree, knots)[0] for row in net]
    return tuple(rows), degree - 1, tuple(knots[1:-1])


def tensor_point(
    net: Sequence[Sequence[Sequence[float]]],
    degrees: tuple[int, int],
    knots: tuple[Sequence[float], Sequence[float]],
    u: float,
    v: float,
    spans: tuple[int | None, int | None] = (None, None),
) -> tuple[float, ...]:
    """Evaluate a tensor-product B-spline patch of vector-valued control points."""
    p, q = degrees
    knot_u, knot_v = knots
    width = len(net[0][0])
    if p < 0 or q < 0:
        return tuple(0.0 for _ in range(width))
    span_u = find_span(knot_u, max(p, 0), u) if spans[0] is None else spans[0]
    span_v = find_span(knot_v, max(q, 0), v) if spans[1] is None else spans[1]
    nu = _basis_values(span_u, u, p, knot_u)
    nv = _basis_values(span_v, v, q, knot_v)
    base_u, base_v = span_u - p, span_v - q
    accumulator = [0.0] * width
    for i in range(p + 1):
        if nu[i] == 0.0:
            continue
        row = net[base_u + i]
        for j in range(q + 1):
            weight = nu[i] * nv[j]
            if weight == 0.0:
                continue
            point = row[base_v + j]
            for c in range(width):
                accumulator[c] += weight * point[c]
    return tuple(accumulator)


# ---------------------------------------------------------------------------
# Gauss-Legendre quadrature
# ---------------------------------------------------------------------------
def gauss_legendre(order: int) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Nodes and weights of ``order``-point Gauss-Legendre quadrature on [-1, 1].

    Newton iteration on the Legendre polynomial (Abramowitz & Stegun 25.4.29): a
    handful of steps reach round-off for the orders used here, and it keeps numpy out
    of the core package.
    """
    if order < 1:
        raise NURBSError(f"quadrature order must be >= 1, got {order}")
    nodes: list[float] = []
    weights: list[float] = []
    half = (order + 1) // 2
    for i in range(1, half + 1):
        x = math.cos(math.pi * (i - 0.25) / (order + 0.5))
        derivative = 1.0
        for _ in range(100):
            p0, p1 = 1.0, x
            for j in range(2, order + 1):
                p0, p1 = p1, ((2 * j - 1) * x * p1 - (j - 1) * p0) / j
            derivative = order * (x * p1 - p0) / (x * x - 1.0)
            step = p1 / derivative
            x -= step
            if abs(step) < 1e-16:
                break
        weights.append(2.0 / ((1.0 - x * x) * derivative * derivative))
        nodes.append(x)
    pairs = list(zip(nodes, weights, strict=True))
    mirrored = [(-node, weight) for node, weight in pairs if order % 2 == 0 or abs(node) > 1e-15]
    pairs.extend(mirrored)
    pairs.sort()
    return tuple(pair[0] for pair in pairs), tuple(pair[1] for pair in pairs)


# ---------------------------------------------------------------------------
# small vector helpers (plain tuples, matching geometry.Mesh)
# ---------------------------------------------------------------------------
def _norm(vector: Sequence[float]) -> float:
    return math.sqrt(sum(value * value for value in vector))


def _unit(vector: Sequence[float]) -> tuple[float, float, float]:
    length = _norm(vector)
    if length == 0.0:
        raise NURBSError("cannot normalise a zero vector")
    return (vector[0] / length, vector[1] / length, vector[2] / length)


def _sub(a: Sequence[float], b: Sequence[float]) -> tuple[float, float, float]:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _cross(a: Sequence[float], b: Sequence[float]) -> tuple[float, float, float]:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _dot(a: Sequence[float], b: Sequence[float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _span_intervals(
    knots: Sequence[float], degree: int, count: int
) -> tuple[tuple[float, float], ...]:
    """The non-empty knot spans inside the parameter domain, in increasing order."""
    intervals: list[tuple[float, float]] = []
    for index in range(degree, count):
        start, end = knots[index], knots[index + 1]
        if end > start:
            intervals.append((start, end))
    return tuple(intervals)


def _interior_knots(knots: Sequence[float], degree: int, count: int) -> tuple[float, ...]:
    low, high = knots[degree], knots[count]
    return tuple(sorted({value for value in knots if low < value < high}))


def _span_before(knots: Sequence[float], degree: int, value: float) -> int:
    """Index of the span lying left of an interior knot value (one-sided limit)."""
    first = min(index for index, knot in enumerate(knots) if knot == value)
    span = first - 1
    if span < degree:
        raise NURBSError(f"no span left of interior knot {value!r}")
    return span


# ---------------------------------------------------------------------------
# curve
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class NURBSCurve:
    """A rational B-spline curve in 3-D.

    ``control_points`` and ``weights`` describe the homogeneous control points
    ``w_i (x_i, y_i, z_i)``. ``weights`` defaults to all ones, which degenerates to a
    polynomial B-spline (``rational`` is then ``False``).
    """

    degree: int
    knots: tuple[float, ...]
    control_points: tuple[Vec3, ...]
    weights: tuple[float, ...] = ()

    def __post_init__(self) -> None:
        degree = int(self.degree)
        points = tuple(
            (float(point[0]) + 0.0, float(point[1]) + 0.0, float(point[2]) + 0.0)
            for point in self.control_points
        )
        if not points:
            raise NURBSError("a curve needs at least one control point")
        for index, point in enumerate(points):
            if not all(math.isfinite(value) for value in point):
                raise NURBSError(f"control_points[{index}] is not finite: {point!r}")
        weights = tuple(float(weight) for weight in (self.weights or tuple(1.0 for _ in points)))
        knots = validate_knots(self.knots, degree, len(points))
        weights = validate_weights(weights, len(points))
        object.__setattr__(self, "degree", degree)
        object.__setattr__(self, "control_points", points)
        object.__setattr__(self, "weights", weights)
        object.__setattr__(self, "knots", knots)

    # -- introspection ------------------------------------------------------
    @property
    def rational(self) -> bool:
        return any(weight != 1.0 for weight in self.weights)

    @property
    def domain(self) -> tuple[float, float]:
        return (self.knots[self.degree], self.knots[len(self.control_points)])

    @property
    def interior_knots(self) -> tuple[float, ...]:
        return _interior_knots(self.knots, self.degree, len(self.control_points))

    # -- evaluation ---------------------------------------------------------
    def _homogeneous(
        self, u: float, order: int, span: int | None
    ) -> tuple[tuple[tuple[float, float, float, float], ...], int]:
        """Homogeneous derivatives ``(w x, w y, w z, w)`` up to ``order``."""
        points = tuple(
            (weight * point[0], weight * point[1], weight * point[2], weight)
            for point, weight in zip(self.control_points, self.weights, strict=True)
        )
        raw = spline_derivatives(points, self.degree, self.knots, u, order, span=span)
        used = find_span(self.knots, self.degree, u) if span is None else span
        return tuple(raw), used

    def derivatives(
        self, u: float, order: int = 1, *, span: int | None = None
    ) -> tuple[Vec3, ...]:
        """Euclidean derivatives up to ``order`` (rational quotient recursion A3.2)."""
        homogeneous, _ = self._homogeneous(u, order, span)
        w0 = homogeneous[0][3]
        if w0 == 0.0:
            raise NURBSError(f"degenerate evaluation at u={u!r}: weight sums to zero")
        out: list[Vec3] = [
            (homogeneous[0][0] / w0, homogeneous[0][1] / w0, homogeneous[0][2] / w0)
        ]
        for k in range(1, order + 1):
            vector = list(homogeneous[k][:3])
            for i in range(1, k + 1):
                factor = math.comb(k, i) * homogeneous[i][3]
                previous = out[k - i]
                for axis in range(3):
                    vector[axis] -= factor * previous[axis]
            out.append((vector[0] / w0, vector[1] / w0, vector[2] / w0))
        return tuple(out)

    def evaluate(self, u: float) -> Vec3:
        return self.derivatives(u, 0)[0]

    def tangent(self, u: float) -> Vec3:
        return _unit(self.derivatives(u, 1)[1])

    # -- integrals ----------------------------------------------------------
    def arc_length(self, order: int = 24) -> float:
        """Arc length ``int |C'(u)| du``, Gauss-Legendre per knot span."""
        nodes, weights = gauss_legendre(order)
        total = 0.0
        for low, high in _span_intervals(self.knots, self.degree, len(self.control_points)):
            half = 0.5 * (high - low)
            for node, weight in zip(nodes, weights, strict=True):
                total += weight * half * _norm(self.derivatives(half * (node + 1.0) + low, 1)[1])
        return total

    def polyline(self, samples: int = 256) -> tuple[Vec3, ...]:
        if samples < 2:
            raise NURBSError(f"samples must be >= 2, got {samples}")
        low, high = self.domain
        step = (high - low) / (samples - 1)
        return tuple(self.evaluate(low + index * step) for index in range(samples))


    def elevate_degree(self, times: int = 1) -> NURBSCurve:
        """A curve of degree ``degree + times`` tracing exactly the same points.

        Standard CAD operation (ISO 10303-42 degree elevation; Piegl & Tiller A5.9,
        exposed by OCCT as ``Geom_BSplineCurve::IncreaseDegree``). A curve cannot be
        stitched to a patch of higher degree without it, so this is the curve-level
        half of the assembly operation the surface service lacks.

        The representation grows minimally and keeps the original continuity: a knot
        of multiplicity ``m`` becomes ``m + times``, so an interior ``C2`` join stays
        ``C2`` in the larger degree. The control points are the unique coefficients of
        the elevated curve in that space, obtained by collocation at the Greville
        abscissae of the new knot vector (the collocation matrix is banded and totally
        positive, so elimination needs no pivoting). Because the solve is done in
        homogeneous coordinates the operation is exact for rational curves too:
        ``elevated.evaluate(u) == self.evaluate(u)`` to round-off at the same ``u``.
        """
        if not isinstance(times, int) or isinstance(times, bool):
            raise NURBSError(f"times must be a non-negative integer, got {times!r}")
        if times < 0:
            raise NURBSError(f"times must be a non-negative integer, got {times!r}")
        if times == 0:
            return self

        homogeneous = tuple(
            (weight * point[0], weight * point[1], weight * point[2], weight)
            for point, weight in zip(self.control_points, self.weights, strict=True)
        )
        target, knots, net = _elevate_homogeneous(homogeneous, self.degree, self.knots, times)
        uniform = _uniform_weight(homogeneous)
        weights = (
            tuple(point[3] for point in net)
            if uniform is None
            else tuple(uniform for _ in net)
        )
        return NURBSCurve(
            degree=target,
            knots=knots,
            control_points=tuple(_dehomogenise(point) for point in net),
            weights=weights,
        )

    # -- piecewise (break-aware) operations ---------------------------------
    def split_at_discontinuities(self) -> tuple[NURBSCurve, ...]:
        """The maximal ``C0``-or-better pieces, split where the curve is broken.

        An interior knot of multiplicity ``degree + 1`` is not a ``C0`` join: the two
        sides share no control point, so a single-curve collocation matrix is singular
        there and :meth:`elevate_degree` (rightly) refuses. This method is the
        splitting step that refusal asks for; :meth:`elevate_degree_piecewise`
        composes the two.

        The pieces are exact, not refitted: because a break already carries
        multiplicity ``degree + 1``, the curve is the direct sum of its pieces -- the
        i-th piece's control points are a contiguous window of the original polygon and
        its knots the corresponding sub-range of the original knot vector. Nothing is
        inserted, removed or re-solved, so every interior knot (and therefore every
        ``C0``/``C1`` join *inside* a piece) is preserved verbatim, parameters are
        unchanged (the i-th piece is parameterised on the knot interval it came from),
        and rational weights travel untouched. A curve with no break returns
        ``(self,)``: the identity, not a copy.

        Only true breaks (multiplicity ``degree + 1``) split the curve. Note that the
        *whole* broken curve is point-evaluable but not differentiable across a break
        (the derivative net divides by a zero-length knot interval), which is exactly
        why ``arc_length`` and ``curvature`` need the pieces this method returns.
        """
        degree = self.degree
        low, high = self.domain
        interior = sorted({value for value in self.knots if low < value < high})
        breaks = [value for value in interior if _multiplicity(self.knots, value) > degree]
        if not breaks:
            return (self,)

        knots = self.knots
        bounds = [low, *breaks, high]
        pieces: list[NURBSCurve] = []
        for left, right in pairwise(bounds):
            last_left = max(index for index, value in enumerate(knots) if value == left)
            piece_knots = [left] * (degree + 1)
            piece_knots.extend(value for value in knots if left < value < right)
            piece_knots.extend([right] * (degree + 1))
            count = len(piece_knots) - degree - 1
            start = last_left - degree
            pieces.append(
                NURBSCurve(
                    degree=degree,
                    knots=tuple(piece_knots),
                    control_points=self.control_points[start : start + count],
                    weights=self.weights[start : start + count],
                )
            )
        return tuple(pieces)

    def elevate_degree_piecewise(self, times: int = 1) -> tuple[NURBSCurve, ...]:
        """Elevate every continuous piece of a curve that is allowed to have breaks.

        The piecewise counterpart of :meth:`elevate_degree`: split at every
        multiplicity ``degree + 1`` knot (:meth:`split_at_discontinuities`) and elevate
        each piece. The result is a tuple of segments, because a break stays a break at
        every degree -- no single spline of degree ``degree + times`` represents a
        broken curve. Every segment is exact, so the union traces the original curve
        point-for-point while each segment is now usable by operations that need one
        connected control polygon (patch-degree alignment against a neighbour,
        curvature, collocation).
        """
        if not isinstance(times, int) or isinstance(times, bool) or times < 0:
            raise NURBSError(f"times must be a non-negative integer, got {times!r}")
        return tuple(piece.elevate_degree(times) for piece in self.split_at_discontinuities())

    # -- continuity ---------------------------------------------------------
    def continuity(self, order: int = 2) -> dict[str, Any]:
        records = curve_derivative_continuity(self, order)
        return {
            "kind": "curve",
            "degree": self.degree,
            "rational": self.rational,
            "knots": [record.to_dict() for record in records],
            "declared_minimum_order": min(
                (record.declared_order for record in records), default=self.degree
            ),
            "passed": all(record.passed for record in records),
        }


# ---------------------------------------------------------------------------
# continuity records
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class ContinuityRecord:
    """One interior knot: what continuity is declared, and what is measured."""

    direction: str
    knot: float
    multiplicity: int
    degree: int
    declared_order: int
    jumps: tuple[float, ...]
    scales: tuple[float, ...]
    max_measurable_order: int = 0

    @property
    def passed(self) -> bool:
        """True when every derivative up to the declared order is continuous.

        The tolerance is relative to the magnitude of the derivative itself: an
        absolute bound on a curvature-scale derivative (1e-3 m^-2) would be a
        meaningless test, so each order is scaled by its own magnitude.
        """
        for order in range(min(self.declared_order, len(self.jumps) - 1) + 1):
            if self.jumps[order] > 1e-9 * max(self.scales[order], 1.0):
                return False
        return True

    @property
    def measured_continuous_order(self) -> int:
        """Highest derivative order the one-sided measurement finds continuous.

        For the three-arc circle this is 1 while the B-spline bookkeeping only
        promises 0: the measurement is what turns "the joins look fine" into a number.
        """
        achieved = -1
        for order, jump in enumerate(self.jumps):
            if jump <= 1e-9 * max(self.scales[order], 1.0):
                achieved = order
            else:
                break
        return achieved

    def to_dict(self) -> dict[str, Any]:
        return {
            "direction": self.direction,
            "knot": self.knot,
            "multiplicity": self.multiplicity,
            "degree": self.degree,
            "declared_order": self.declared_order,
            "max_measurable_order": self.max_measurable_order,
            "measured_continuous_order": self.measured_continuous_order,
            "jumps": list(self.jumps),
            "scales": list(self.scales),
            "passed": self.passed,
        }


def _multiplicity(knots: Sequence[float], value: float) -> int:
    return sum(1 for knot in knots if knot == value)


def curve_derivative_continuity(
    curve: NURBSCurve, order: int = 2
) -> tuple[ContinuityRecord, ...]:
    """One-sided derivative comparison at every interior knot of ``curve``.

    The jump of order ``k`` is ``|C^(k)(knot+0) - C^(k)(knot-0)|`` computed from the
    analytic derivatives of the two neighbouring spans, so it is a property of the
    surface data and not of a finite-difference step size.
    """
    interior = curve.interior_knots
    if not interior:
        return ()
    # A derivative of order p - m + 1 does not exist as a *function* across a knot of
    # multiplicity m (it is a distribution there), so the measurement is clamped to the
    # highest order that is representable.
    highest = min(
        curve.degree - _multiplicity(curve.knots, value) + 1 for value in interior
    )
    order = max(min(order, highest), 0)
    records: list[ContinuityRecord] = []
    for value in interior:
        multiplicity = _multiplicity(curve.knots, value)
        right_span = find_span(curve.knots, curve.degree, value)
        left_span = _span_before(curve.knots, curve.degree, value)
        left = curve.derivatives(value, order, span=left_span)
        right = curve.derivatives(value, order, span=right_span)
        records.append(
            ContinuityRecord(
                direction="u",
                knot=value,
                multiplicity=multiplicity,
                degree=curve.degree,
                declared_order=max(curve.degree - multiplicity, 0),
                jumps=tuple(_norm(_sub(right[k], left[k])) for k in range(order + 1)),
                scales=tuple(max(_norm(left[k]), _norm(right[k])) for k in range(order + 1)),
                max_measurable_order=highest,
            )
        )
    return tuple(records)

# ---------------------------------------------------------------------------
# surface
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class NURBSSurface:
    """A tensor-product rational B-spline surface.

    ``control_net[i][j]`` is the control point at ``u``-index ``i`` and ``v``-index
    ``j``; ``weights`` shares that layout. The surface is the projective image of the
    4-D tensor-product B-spline through ``(w x, w y, w z, w)``.
    """

    degree: tuple[int, int]
    knots: tuple[tuple[float, ...], tuple[float, ...]]
    control_net: tuple[tuple[Vec3, ...], ...]
    weights: tuple[tuple[float, ...], ...] = ()

    def __post_init__(self) -> None:
        degree = (int(self.degree[0]), int(self.degree[1]))
        net = tuple(
            tuple((float(point[0]) + 0.0, float(point[1]) + 0.0, float(point[2]) + 0.0)
                  for point in row)
            for row in self.control_net
        )
        if not net or not net[0]:
            raise NURBSError("control net is empty")
        width = len(net[0])
        for row_index, row in enumerate(net):
            if len(row) != width:
                raise NURBSError(
                    f"control net row {row_index} has {len(row)} points, expected {width}"
                )
            for column_index, point in enumerate(row):
                if not all(math.isfinite(value) for value in point):
                    raise NURBSError(
                        f"control_net[{row_index}][{column_index}] is not finite: {point!r}"
                    )
        weights = self.weights or tuple(tuple(1.0 for _ in row) for row in net)
        if len(weights) != len(net) or any(len(row) != width for row in weights):
            raise NURBSError(
                f"weights layout does not match the control net {len(net)}x{width}"
            )
        checked = tuple(
            validate_weights(tuple(row), width) for row in weights
        )
        knots_u = validate_knots(self.knots[0], degree[0], len(net))
        knots_v = validate_knots(self.knots[1], degree[1], width)
        object.__setattr__(self, "degree", degree)
        object.__setattr__(self, "control_net", net)
        object.__setattr__(self, "weights", checked)
        object.__setattr__(self, "knots", (knots_u, knots_v))

    # -- introspection ------------------------------------------------------
    @property
    def rational(self) -> bool:
        return any(weight != 1.0 for row in self.weights for weight in row)

    @property
    def domain(self) -> tuple[tuple[float, float], tuple[float, float]]:
        return (
            (self.knots[0][self.degree[0]], self.knots[0][len(self.control_net)]),
            (self.knots[1][self.degree[1]], self.knots[1][len(self.control_net[0])]),
        )

    @property
    def interior_knots(self) -> tuple[tuple[float, ...], tuple[float, ...]]:
        return (
            _interior_knots(self.knots[0], self.degree[0], len(self.control_net)),
            _interior_knots(self.knots[1], self.degree[1], len(self.control_net[0])),
        )

    # -- evaluation ---------------------------------------------------------
    def derivatives(
        self,
        u: float,
        v: float,
        order: int = 2,
        *,
        spans: tuple[int | None, int | None] | None = None,
    ) -> dict[tuple[int, int], Vec3]:
        """Euclidean partial derivatives ``S^(k,l)`` for ``k + l <= order`` (A3.6).

        Both directions are differentiated by building the derivative control net
        (exact and structure-preserving), the tensor product is evaluated with the
        verified basis functions, and the rational quotient is then unwound with the
        standard recursion. So the ``(1, 0)`` entry is the true ``dS/du``, not a
        finite-difference approximation of it.
        """
        p, q = self.degree
        span_u, span_v = spans if spans is not None else (None, None)
        net: tuple[tuple[tuple[float, ...], ...], ...] = tuple(
            tuple(
                (
                    self.weights[i][j] * self.control_net[i][j][0],
                    self.weights[i][j] * self.control_net[i][j][1],
                    self.weights[i][j] * self.control_net[i][j][2],
                    self.weights[i][j],
                )
                for j in range(len(self.control_net[i]))
            )
            for i in range(len(self.control_net))
        )
        levels_u: list[tuple[object, int, tuple[float, ...]]] = [(net, p, tuple(self.knots[0]))]
        for _ in range(1, min(order, p) + 1):
            reduced, degree, reduced_knots = _differentiate_axis(
                levels_u[-1][0], 0, levels_u[-1][1], levels_u[-1][2]
            )
            levels_u.append((reduced, degree, reduced_knots))
        homogeneous: dict[tuple[int, int], tuple[float, float, float, float]] = {}
        for k in range(order + 1):
            if k > p:
                for l in range(order - k + 1):
                    homogeneous[(k, l)] = (0.0, 0.0, 0.0, 0.0)
                continue
            base_net, base_degree, base_knots = levels_u[k]
            levels_v: list[tuple[object, int, tuple[float, ...]]] = [
                (base_net, q, tuple(self.knots[1]))
            ]
            for _ in range(1, min(order - k, q) + 1):
                reduced, degree, reduced_knots = _differentiate_axis(
                    levels_v[-1][0], 1, levels_v[-1][1], levels_v[-1][2]
                )
                levels_v.append((reduced, degree, reduced_knots))
            for l in range(order - k + 1):
                if l > q:
                    homogeneous[(k, l)] = (0.0, 0.0, 0.0, 0.0)
                    continue
                level_net, level_degree, level_knots = levels_v[l]
                reduced_span_u = None if span_u is None else span_u - k
                reduced_span_v = None if span_v is None else span_v - l
                homogeneous[(k, l)] = tensor_point(
                    level_net,
                    (base_degree, level_degree),
                    (base_knots, level_knots),
                    u,
                    v,
                    (reduced_span_u, reduced_span_v),
                )
        w0 = homogeneous[(0, 0)][3]
        if w0 == 0.0:
            raise NURBSError(f"degenerate evaluation at (u, v)=({u!r}, {v!r})")
        out: dict[tuple[int, int], Vec3] = {}
        for k in range(order + 1):
            for l in range(order + 1 - k):
                vector = list(homogeneous[(k, l)][:3])
                for j in range(1, l + 1):
                    factor = math.comb(l, j) * homogeneous[(0, j)][3]
                    for axis in range(3):
                        vector[axis] -= factor * out[(k, l - j)][axis]
                for i in range(1, k + 1):
                    factor = math.comb(k, i) * homogeneous[(i, 0)][3]
                    for axis in range(3):
                        vector[axis] -= factor * out[(k - i, l)][axis]
                    for j in range(1, l + 1):
                        factor = math.comb(k, i) * math.comb(l, j) * homogeneous[(i, j)][3]
                        for axis in range(3):
                            vector[axis] -= factor * out[(k - i, l - j)][axis]
                out[(k, l)] = (vector[0] / w0, vector[1] / w0, vector[2] / w0)
        return out

    def evaluate(self, u: float, v: float) -> Vec3:
        return self.derivatives(u, v, 0)[(0, 0)]

    def elevate_degree(self, times_u: int = 0, times_v: int = 0) -> NURBSSurface:
        """A surface of degree ``degree + (times_u, times_v)`` tracing the same patch.

        The tensor-product counterpart of :meth:`NURBSCurve.elevate_degree`: each
        fixed-``v`` column is a curve in ``u`` and each fixed-``u`` row is a curve
        in ``v``, so the patch is elevated one direction at a time. Because the
        operation is linear in homogeneous control points, elevating ``u`` then
        ``v`` equals elevating ``v`` then ``u`` (asserted in
        ``tests/test_nurbs_surface_degree_elevation.py``).

        This is the surface-level half of patch-degree alignment: patches of an
        assembly that disagree in degree cannot share one edge representation
        without it. Trimming and stitching remain outside this module.

        Everything the curve version promises carries over per direction: the
        surfaces agree to round-off at the same ``(u, v)``, an interior knot of
        multiplicity ``m`` becomes ``m + times`` (so a ``C1`` join stays ``C1``),
        the control net grows by the minimum ``times * (#spans)`` points, and the
        solve happens in homogeneous coordinates, so rational patches (the exact
        cylinder and sphere anchors) stay exact.
        """
        for name, value in (("times_u", times_u), ("times_v", times_v)):
            if not isinstance(value, int) or isinstance(value, bool):
                raise NURBSError(f"{name} must be a non-negative integer, got {value!r}")
            if value < 0:
                raise NURBSError(f"{name} must be a non-negative integer, got {value!r}")
        if times_u == 0 and times_v == 0:
            return self

        width = len(self.control_net[0])
        net = [
            [
                (
                    self.weights[i][j] * self.control_net[i][j][0],
                    self.weights[i][j] * self.control_net[i][j][1],
                    self.weights[i][j] * self.control_net[i][j][2],
                    self.weights[i][j],
                )
                for j in range(width)
            ]
            for i in range(len(self.control_net))
        ]
        degree = [self.degree[0], self.degree[1]]
        knots = [self.knots[0], self.knots[1]]
        uniform = _uniform_weight([point for row in net for point in row])

        if times_u:
            columns = [[net[i][j] for i in range(len(net))] for j in range(width)]
            elevated = [
                _elevate_homogeneous(column, degree[0], knots[0], times_u, label="u-column")
                for column in columns
            ]
            degree[0] = elevated[0][0]
            knots[0] = elevated[0][1]
            net = [
                [elevated[j][2][i] for j in range(width)] for i in range(len(elevated[0][2]))
            ]
            if uniform is not None:
                net = [_rescale_to_uniform(row, uniform) for row in net]

        if times_v:
            rows = [
                _elevate_homogeneous(row, degree[1], knots[1], times_v, label="v-row")
                for row in net
            ]
            degree[1] = rows[0][0]
            knots[1] = rows[0][1]
            net = [list(row[2]) for row in rows]
            if uniform is not None:
                net = [_rescale_to_uniform(row, uniform) for row in net]

        control_net = tuple(tuple(_dehomogenise(point) for point in row) for row in net)
        weights = (
            tuple(tuple(point[3] for point in row) for row in net)
            if uniform is None
            else tuple(tuple(uniform for _ in row) for row in net)
        )
        return NURBSSurface(
            degree=(degree[0], degree[1]),
            knots=(knots[0], knots[1]),
            control_net=control_net,
            weights=weights,
        )

    def jacobian(self, u: float, v: float) -> tuple[Vec3, Vec3, Vec3]:
        """``(S_u, S_v, S_u x S_v)``; the third entry is the normal density vector."""
        ders = self.derivatives(u, v, 1)
        return ders[(1, 0)], ders[(0, 1)], _cross(ders[(1, 0)], ders[(0, 1)])

    def normal(self, u: float, v: float) -> Vec3:
        """Unit normal ``S_u x S_v / |S_u x S_v|``; a zero-area patch is rejected."""
        _, _, cross = self.jacobian(u, v)
        length = _norm(cross)
        if length == 0.0:
            raise NURBSError(f"degenerate patch at (u, v)=({u!r}, {v!r}): zero Jacobian")
        return (cross[0] / length, cross[1] / length, cross[2] / length)

    # -- integrals ----------------------------------------------------------
    def area(self, order: int = 24) -> float:
        """Surface area ``int int |S_u x S_v| du dv``, Gauss-Legendre per knot span."""
        nodes, weights = gauss_legendre(order)
        total = 0.0
        for ua, ub in _span_intervals(self.knots[0], self.degree[0], len(self.control_net)):
            u_half = 0.5 * (ub - ua)
            for va, vb in _span_intervals(
                self.knots[1], self.degree[1], len(self.control_net[0])
            ):
                v_half = 0.5 * (vb - va)
                for node_u, weight_u in zip(nodes, weights, strict=True):
                    u = u_half * (node_u + 1.0) + ua
                    for node_v, weight_v in zip(nodes, weights, strict=True):
                        v = v_half * (node_v + 1.0) + va
                        _, _, cross = self.jacobian(u, v)
                        total += weight_u * weight_v * u_half * v_half * _norm(cross)
        return total

    def area_convergence(self, orders: Sequence[int] = (2, 3, 4, 6, 8, 12, 16, 24)) -> dict[str, Any]:
        """Quadrature order vs integrated area: the whole error curve, not one number."""
        reference = self.area(max(max(orders), 48))
        rows = []
        for order in orders:
            value = self.area(order)
            rows.append(
                {
                    "order": order,
                    "area": value,
                    "relative_error": abs(value - reference) / reference if reference else math.inf,
                }
            )
        return {"reference_area": reference, "rows": rows}

    # -- continuity ---------------------------------------------------------
    def continuity(self, order: int = 1) -> dict[str, Any]:
        records = surface_derivative_continuity(self, order)
        return {
            "kind": "surface",
            "degree": list(self.degree),
            "rational": self.rational,
            "u": [record.to_dict() for record in records["u"]],
            "v": [record.to_dict() for record in records["v"]],
            "passed": all(record.passed for record in records["u"] + records["v"]),
        }

    # -- tessellation -------------------------------------------------------
    def singular_boundary_point(self, axis: str, *, upper: bool = False) -> Vec3 | None:
        """The point a coincident boundary control row collapses to, ``None`` if none.

        A row of the control net whose control points are *all the same point*
        forces the corresponding boundary curve of the patch to that point: the
        north/south pole of a sphere, the apex of a cone. The parametrisation is
        singular there -- every value of the other parameter maps to that one
        point -- and evaluating it sample by sample gives values that agree only
        to the last bit, so exact-position vertex welding would leave a ragged
        pole ring. Returning the control point lets :meth:`tessellate` key the
        whole row on a single vertex, which is a fact about the control net
        rather than a tolerance.
        """
        if axis not in ("u", "v"):
            raise NURBSError(f"axis must be 'u' or 'v', got {axis!r}")
        if axis == "v":
            index = len(self.control_net[0]) - 1 if upper else 0
            row = [self.control_net[i][index] for i in range(len(self.control_net))]
        else:
            row = list(self.control_net[-1 if upper else 0])
        first = row[0]
        return first if all(point == first for point in row) else None

    def tessellate(self, u_samples: int = 64, v_samples: int = 8) -> Tessellation:
        """Sample the surface into a triangle mesh with consistent winding.

        Each quad ``(u_i, v_j) - (u_i+1, v_j) - (u_i+1, v_j+1) - (u_i, v_j+1)`` is
        traversed so that its normal is ``S_u x S_v``. Vertices are shared by exact
        position, so a seam that closes in the geometry (a full revolution, where
        ``S(u_min, v) == S(u_max, v)`` to the last bit) also closes in the topology --
        the watertightness check then has something real to say.

        A *singular* boundary row (every control point of the boundary row is the
        same point: a sphere's poles, a cone's apex) is keyed on that one vertex
        instead of on the per-sample evaluation, because the evaluation only agrees
        with itself to the last bit, not bit-for-bit (see
        :meth:`singular_boundary_point`). Without that, the pole quads do not
        collapse and the shell is not watertight at any density.
        """
        if u_samples < 1 or v_samples < 1:
            raise NURBSError(
                f"tessellation needs at least 1 x 1 samples, got {u_samples} x {v_samples}"
            )
        (u_low, u_high), (v_low, v_high) = self.domain
        u_values = [u_low + (u_high - u_low) * index / u_samples for index in range(u_samples + 1)]
        v_values = [v_low + (v_high - v_low) * index / v_samples for index in range(v_samples + 1)]
        u_values[0], u_values[-1] = u_low, u_high
        v_values[0], v_values[-1] = v_low, v_high
        builder = MeshBuilder()
        pole_v_low = self.singular_boundary_point("v")
        pole_v_high = self.singular_boundary_point("v", upper=True)
        pole_u_low = self.singular_boundary_point("u")
        pole_u_high = self.singular_boundary_point("u", upper=True)
        last_u = len(u_values) - 1
        last_v = len(v_values) - 1
        grid: list[list[int]] = []
        for i, u in enumerate(u_values):
            row: list[int] = []
            pole_u = pole_u_low if i == 0 else (pole_u_high if i == last_u else None)
            for j, v in enumerate(v_values):
                if pole_u is not None:
                    row.append(builder.add_vertex(pole_u))
                elif j == 0 and pole_v_low is not None:
                    row.append(builder.add_vertex(pole_v_low))
                elif j == last_v and pole_v_high is not None:
                    row.append(builder.add_vertex(pole_v_high))
                else:
                    row.append(builder.add_vertex(self.evaluate(u, v)))
            grid.append(row)
        faces_per_quad: list[tuple[int, ...]] = []
        emitted = 0
        for i in range(u_samples):
            for j in range(v_samples):
                produced = builder.add_quad(
                    grid[i][j], grid[i + 1][j], grid[i + 1][j + 1], grid[i][j + 1]
                )
                faces_per_quad.append(tuple(range(emitted, emitted + produced)))
                emitted += produced
        return Tessellation(
            mesh=builder.build(),
            u=tuple(u_values),
            v=tuple(v_values),
            degree=self.degree,
            degenerate_skipped=builder.degenerate_skipped,
            faces_per_quad=tuple(faces_per_quad),
        )


def surface_derivative_continuity(
    surface: NURBSSurface, order: int = 1
) -> dict[str, tuple[ContinuityRecord, ...]]:
    """Continuity of the surface derivatives across the interior knots of U and V.

    For each interior knot the derivative is evaluated twice at the same parameter
    value, once from the span on its left and once from the span on its right. In the
    other direction the mid-parameter is used, so the comparison is uncontaminated.
    """
    out: dict[str, tuple[ContinuityRecord, ...]] = {}
    (u_low, u_high), (v_low, v_high) = surface.domain
    for axis, label in ((0, "u"), (1, "v")):
        knots = surface.knots[axis]
        degree = surface.degree[axis]
        low = u_low if axis == 0 else v_low
        high = u_high if axis == 0 else v_high
        interior = sorted({knot for knot in knots if low < knot < high})
        if not interior:
            out[label] = ()
            continue
        highest = min(degree - _multiplicity(knots, knot) + 1 for knot in interior)
        axis_order = max(min(order, highest), 0)
        records: list[ContinuityRecord] = []
        for value in interior:
            multiplicity = _multiplicity(knots, value)
            u = value if axis == 0 else 0.5 * (u_low + u_high)
            v = 0.5 * (v_low + v_high) if axis == 0 else value
            left_span = _span_before(knots, degree, value)
            right_span = find_span(knots, degree, value)
            spans_left = (left_span, None) if axis == 0 else (None, left_span)
            spans_right = (right_span, None) if axis == 0 else (None, right_span)
            left = surface.derivatives(u, v, axis_order, spans=spans_left)
            right = surface.derivatives(u, v, axis_order, spans=spans_right)
            key = (1, 0) if axis == 0 else (0, 1)
            jumps = [_norm(_sub(right[(0, 0)], left[(0, 0)]))]
            scales = [max(_norm(left[(0, 0)]), _norm(right[(0, 0)]))]
            for k in range(1, axis_order + 1):
                component = key if k == 1 else (k, 0) if axis == 0 else (0, k)
                jumps.append(_norm(_sub(right[component], left[component])))
                scales.append(max(_norm(left[component]), _norm(right[component])))
            records.append(
                ContinuityRecord(
                    direction=label,
                    knot=value,
                    multiplicity=multiplicity,
                    degree=degree,
                    declared_order=max(degree - multiplicity, 0),
                    jumps=tuple(jumps),
                    scales=tuple(scales),
                    max_measurable_order=highest,
                )
            )
        out[label] = tuple(records)
    return out


# ---------------------------------------------------------------------------
# exact arcs, the cylinder, and the closed solid used by the watertightness check
# ---------------------------------------------------------------------------
def circular_arc(
    radius: float,
    sweep: float = 0.5 * math.pi,
    *,
    start_angle: float = 0.0,
    z: float = 0.0,
) -> NURBSCurve:
    """The *exact* rational quadratic representation of a circular arc.

    Control points are the start, the intersection of the end tangents and the end,
    with weights ``(1, cos(sweep / 2), 1)``. A conic is thus represented, not
    approximated: the point at any parameter lies on the circle to round-off.
    ``sweep`` must not exceed ``pi`` because a single quadratic arc folds beyond a half
    circle, which is why :func:`circle_curve` stitches three 120 degree arcs.
    """
    if radius <= 0.0:
        raise NURBSError(f"radius must be positive, got {radius!r}")
    if not 0.0 < abs(sweep) <= math.pi * (1.0 + 1e-12):
        raise NURBSError(f"sweep must lie in (0, pi], got {sweep!r}")
    half = 0.5 * sweep
    cos_half = math.cos(half)
    if cos_half <= 0.0:
        raise NURBSError("a 180 degree arc needs the limiting tangent intersection")
    start = start_angle
    mid = start_angle + half
    end = start_angle + sweep
    p0 = (radius * math.cos(start), radius * math.sin(start), z)
    p1 = (radius / cos_half * math.cos(mid), radius / cos_half * math.sin(mid), z)
    p2 = (radius * math.cos(end), radius * math.sin(end), z)
    return NURBSCurve(
        degree=2,
        knots=(0.0, 0.0, 0.0, 1.0, 1.0, 1.0),
        control_points=(p0, p1, p2),
        weights=(1.0, cos_half, 1.0),
    )


def circle_curve(radius: float = 1.0, *, z: float = 0.0, arcs: int = 3) -> NURBSCurve:
    """A full circle as ``arcs`` exact rational quadratic arcs, closed at the seam.

    Interior joints carry multiplicity 2 (``C0`` by the B-spline bookkeeping), yet the
    rational arcs are ``C1`` there: a symmetric arc has the same parametric speed at
    both of its ends, so the join matches. :func:`curve_derivative_continuity`
    measures that instead of taking it on faith.
    """
    if radius <= 0.0:
        raise NURBSError(f"radius must be positive, got {radius!r}")
    if arcs < 2:
        raise NURBSError(f"a closed circle needs at least 2 arcs, got {arcs}")
    sweep = TWO_PI / arcs
    half = 0.5 * sweep
    cos_half = math.cos(half)
    angles = [index * sweep for index in range(arcs + 1)]
    ring = [(radius * math.cos(angle), radius * math.sin(angle), z) for angle in angles[:-1]]
    ring.append(ring[0])  # exact seam closure: the last point *is* the first tuple
    points: list[Vec3] = []
    weights: list[float] = []
    for index in range(arcs):
        mid = angles[index] + half
        points.append(ring[index])
        points.append((radius / cos_half * math.cos(mid), radius / cos_half * math.sin(mid), z))
        weights.extend((1.0, cos_half))
    points.append(ring[-1])
    weights.append(1.0)
    knots = [0.0, 0.0, 0.0]
    for index in range(1, arcs):
        knots.extend((index / arcs, index / arcs))
    knots.extend((1.0, 1.0, 1.0))
    return NURBSCurve(
        degree=2,
        knots=tuple(knots),
        control_points=tuple(points),
        weights=tuple(weights),
    )


def cylinder_surface(
    radius: float,
    height: float,
    *,
    base_z: float = 0.0,
    arcs: int = 3,
) -> NURBSSurface:
    """An exact tensor-product cylinder: ``u`` is the exact circle, ``v`` is linear.

    Because both ``v``-rows share the same rational weights, the quotient leaves
    ``z = base_z + v * height`` untouched and the surface is the true cylinder. Its
    lateral area is therefore ``2 pi r h`` to round-off, which is the analytic anchor
    the acceptance criterion uses.
    """
    if radius <= 0.0 or height <= 0.0:
        raise NURBSError(f"radius and height must be positive, got {radius!r}, {height!r}")
    curve = circle_curve(radius, z=base_z, arcs=arcs)
    base_net = curve.control_points
    top_net = tuple((x, y, z + height) for x, y, z in base_net)
    net = tuple((base_net[i], top_net[i]) for i in range(len(base_net)))
    weights = tuple((curve.weights[i], curve.weights[i]) for i in range(len(base_net)))
    return NURBSSurface(
        degree=(2, 1),
        knots=(curve.knots, (0.0, 0.0, 1.0, 1.0)),
        control_net=net,
        weights=weights,
    )


def sphere_surface(radius: float, *, arcs: int = 3) -> NURBSSurface:
    """An exact rational biquadratic sphere with two *singular* pole rings.

    ``u`` is the exact circle (:func:`circle_curve`) and ``v`` is a semicircle in
    the meridian half-plane built from two exact 90-degree rational quadratic
    arcs (weights ``1, sqrt(1/2), 1``), so the patch is the true sphere rather
    than a fit: every evaluation satisfies ``|S(u, v)| == radius`` to round-off,
    which is the anchor the tests use. The whole sphere is one patch -- it needs
    no caps, unlike :func:`cylinder_surface`.

    The meridian runs from the **south** pole to the **north** pole, so
    ``S_u x S_v`` points *outward* and a tessellation of the surface alone is a
    closed shell with positive signed volume. The pole ordering decides that
    sign: a north-to-south meridian produces the same closed shell with every
    face consistently *inward* and a negative volume (both are documented in
    ``docs/nurbs_pole_surface.md``).

    Both ``v`` ends are singular rows of the control net -- every ``u`` maps to
    the same point -- which :meth:`NURBSSurface.singular_boundary_point` reports
    and :meth:`NURBSSurface.tessellate` welds to a single vertex. The area is
    ``4 pi r^2`` and the volume ``4/3 pi r^3``; ``arcs`` only controls how the
    full circle is stitched (three 120-degree arcs by default).
    """
    if radius <= 0.0:
        raise NURBSError(f"radius must be positive, got {radius!r}")
    meridian = ((0.0, -1.0), (1.0, -1.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))
    meridian_weights = (1.0, SQRT_HALF, 1.0, SQRT_HALF, 1.0)
    meridian_knots = (0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 2.0)
    circle = circle_curve(1.0, arcs=arcs)
    # net[i][j]: i runs over u, j over v.  Scaling the unit circle's control
    # polygon by the meridian control radius r_j makes the rational quotient
    # separate, so S = radius * (r(v) X(u), r(v) Y(u), z(v)) exactly.
    net = tuple(
        tuple((radius * r * point[0], radius * r * point[1], radius * z) for r, z in meridian)
        for point in circle.control_points
    )
    weights = tuple(
        tuple(circle.weights[i] * meridian_weights[j] for j in range(len(meridian)))
        for i in range(len(circle.control_points))
    )
    return NURBSSurface(
        degree=(2, 2),
        knots=(circle.knots, meridian_knots),
        control_net=net,
        weights=weights,
    )


@dataclass(frozen=True)
class Tessellation:
    """A parametric sample grid and the mesh it produced."""

    mesh: Mesh
    u: tuple[float, ...]
    v: tuple[float, ...]
    degree: tuple[int, int]
    degenerate_skipped: int = 0
    faces_per_quad: tuple[tuple[int, ...], ...] = ()

    def faces_of_quad(self, quad: tuple[int, int]) -> tuple[int, ...]:
        """Triangle indices the tessellator emitted for the quad ``(i, j)``."""
        i, j = quad
        width = max(len(self.v) - 1, 1)  # one entry per quad column
        index = i * width + j
        if not self.faces_per_quad:
            return ()
        return self.faces_per_quad[index]

    @property
    def quad_count(self) -> int:
        return (len(self.u) - 1) * (len(self.v) - 1)

    def to_dict(self) -> dict[str, Any]:
        return {
            "u_samples": len(self.u),
            "v_samples": len(self.v),
            "quads": self.quad_count,
            "triangles": len(self.mesh.faces),
            "vertices": len(self.mesh.vertices),
            "degenerate_quads_skipped": self.degenerate_skipped,
            "degree": list(self.degree),
        }


def cylinder_solid_mesh(
    radius: float,
    height: float,
    *,
    base_z: float = 0.0,
    u_samples: int = 192,
    v_samples: int = 8,
) -> Mesh:
    """Tessellate the *closed* cylinder (lateral surface plus two caps).

    The caps are triangle fans whose rim vertices are the very coordinates the lateral
    surface evaluates at ``v = 0`` and ``v = 1``, so the shared positions collapse to
    shared indices and the shell is watertight by construction. The volume of the
    resulting solid converges to ``pi r^2 h`` from below as the rim is refined, which
    is reported as a convergence curve rather than as a single tolerance.
    """
    surface = cylinder_surface(radius, height, base_z=base_z)
    tess = surface.tessellate(u_samples, v_samples)
    builder = MeshBuilder()
    u_values = list(tess.u)
    v_values = list(tess.v)
    nodes = [
        [builder.add_vertex(surface.evaluate(u, v)) for v in v_values] for u in u_values
    ]
    for i in range(len(u_values) - 1):
        for j in range(len(v_values) - 1):
            builder.add_quad(nodes[i][j], nodes[i + 1][j], nodes[i + 1][j + 1], nodes[i][j + 1])
    bottom = builder.add_vertex((0.0, 0.0, base_z))
    top = builder.add_vertex((0.0, 0.0, base_z + height))
    for i in range(len(u_values) - 1):
        low_a, low_b = nodes[i][0], nodes[i + 1][0]
        builder.add_triangle(bottom, low_b, low_a)
        up_a, up_b = nodes[i][-1], nodes[i + 1][-1]
        builder.add_triangle(top, up_a, up_b)
    return builder.build()

# ---------------------------------------------------------------------------
# tessellation audit: watertightness + normal consistency
# ---------------------------------------------------------------------------
def tessellation_report(
    surface: NURBSSurface,
    tessellation: Tessellation,
    *,
    area_reference: float | None = None,
) -> dict[str, Any]:
    """Audit a tessellation: topology, area error, and normal consistency.

    *Topology* is checked on the mesh, not on the intent: every edge must be shared by
    exactly two faces, no edge may be shared by three, and the shell must have exactly
    one connected component (a duplicated seam ring would show up as two).

    *Normals* are compared face by face against the analytic unit normal of the
    surface at the quad's parametric centre. ``inward_faces`` counts faces whose normal
    has a negative dot product with that reference, i.e. faces the solver would see
    from the wrong side; ``max_angle_deg`` is the worst deviation, which must shrink as
    ``O(h^2)`` for a smooth patch.

    *Area* compares the triangle area with the parametric integral: the difference is
    the discretisation error of the tessellation, and it converges from below.
    """
    from . import mesh as meshlib

    mesh = tessellation.mesh
    counts = meshlib.edge_face_counts(mesh)
    components = _component_count(mesh)
    reference_area = surface.area() if area_reference is None else area_reference
    tessellated_area = meshlib.surface_area(mesh)

    u_values = tessellation.u
    v_values = tessellation.v
    worst_angle = 0.0
    angle_sum = 0.0
    face_count = 0
    inward = 0
    flipped = 0
    previous_normal: Vec3 | None = None
    for i in range(len(u_values) - 1):
        for j in range(len(v_values) - 1):
            centre_u = 0.5 * (u_values[i] + u_values[i + 1])
            centre_v = 0.5 * (v_values[j] + v_values[j + 1])
            reference = surface.normal(centre_u, centre_v)
            quad = (i, j)
            for face_index in tessellation.faces_of_quad(quad):
                # face_normal returns the raw cross product; the comparison is
                # against a unit reference, so normalise it first.
                geometric = _unit(meshlib.face_normal(mesh, mesh.faces[face_index]))
                dot = _dot(geometric, reference)
                if dot < 0.0:
                    inward += 1
                if previous_normal is not None and _dot(geometric, previous_normal) < 0.0:
                    flipped += 1
                previous_normal = geometric
                angle = math.degrees(math.acos(max(-1.0, min(1.0, dot))))
                worst_angle = max(worst_angle, angle)
                angle_sum += angle
                face_count += 1

    return {
        "triangles": len(mesh.faces),
        "vertices": len(mesh.vertices),
        "watertight": meshlib.is_watertight(mesh),
        "edge_manifold": meshlib.is_edge_manifold(mesh),
        "boundary_edges": len(meshlib.boundary_edges(mesh)),
        "non_manifold_edges": sum(1 for count in counts.values() if count > 2),
        "components": components,
        "parametric_area_m2": reference_area,
        "tessellated_area_m2": tessellated_area,
        "area_relative_error": (
            (tessellated_area - reference_area) / reference_area if reference_area else math.inf
        ),
        "faces_checked": face_count,
        "inward_faces": inward,
        "orientation_flips": flipped,
        "max_normal_angle_deg": worst_angle,
        "mean_normal_angle_deg": angle_sum / face_count if face_count else math.inf,
        "normal_reference": "S_u x S_v at the quad parametric centre",
    }


def _component_count(mesh: Mesh) -> int:
    """Number of face-connected shells (union-find over shared vertices)."""
    parent = list(range(len(mesh.vertices)))

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    for face in mesh.faces:
        first = find(face[0])
        for vertex in face[1:]:
            other = find(vertex)
            if other != first:
                parent[other] = first
    shells = {find(vertex) for face in mesh.faces for vertex in face}
    return len(shells)


def mesh_solid_report(
    mesh: Mesh,
    *,
    volume_reference: float,
    area_reference: float | None = None,
) -> dict[str, Any]:
    """Audit a closed tessellation against its analytic volume (and optionally area).

    Parallel to :func:`tessellation_report`, but for a *solid*: the volume is what a
    hydrostatics solver reads, so the sign, the watertightness and the convergence of
    ``V_mesh / V_exact`` are reported together. A positive signed volume with every
    edge shared by two faces is the pair of facts an STL export must have before it can
    be called a body.
    """
    from . import mesh as meshlib

    counts = meshlib.edge_face_counts(mesh)
    volume = meshlib.total_signed_volume(mesh)
    result: dict[str, Any] = {
        "triangles": len(mesh.faces),
        "vertices": len(mesh.vertices),
        "watertight": meshlib.is_watertight(mesh),
        "edge_manifold": meshlib.is_edge_manifold(mesh),
        "boundary_edges": len(meshlib.boundary_edges(mesh)),
        "non_manifold_edges": sum(1 for count in counts.values() if count > 2),
        "components": _component_count(mesh),
        "signed_volume_m3": volume,
        "volume_reference_m3": volume_reference,
        "volume_relative_error": (volume - volume_reference) / volume_reference,
    }
    if area_reference is not None:
        area = meshlib.surface_area(mesh)
        result["area_m2"] = area
        result["area_reference_m2"] = area_reference
        result["area_relative_error"] = (area - area_reference) / area_reference
    orientation = meshlib.orient_outward(mesh)
    result["orientation_flipped_faces"] = orientation.flipped_faces
    result["orientation_consistent"] = orientation.flipped_faces == 0
    return result


# ---------------------------------------------------------------------------
# PC6 contract compatibility (additive block, P1's schema untouched)
# ---------------------------------------------------------------------------
def _geometry_gate(
    parameters: Any, declared: dict[str, float] | None
) -> dict[str, Any]:
    """Geometry-first gate: L, B, T and Cb must agree before a surface is published.

    The campaign rule is that geometry is verified before algorithms are trusted, so
    this gate recomputes the closed-form form coefficients of the P1 hull and, when the
    caller *declares* values (for example from an incoming contract), compares them
    against the recomputed ones. A mismatch is an error, not a warning: publishing a
    surface block whose stated Cb disagrees with its own parameters is how two solvers
    end up with two hulls.
    """
    from .pc6 import analytic_form

    form = analytic_form(parameters)
    values = {
        "length_m": float(parameters.length),
        "beam_m": float(parameters.beam),
        "draft_m": float(parameters.draft),
        "displacement_volume_m3": form.volume,
        "block_coefficient": form.block_coefficient,
        "prismatic_coefficient": form.prismatic_coefficient,
        "midship_coefficient": form.midship_coefficient,
    }
    checks: dict[str, bool] = {
        "length_positive": values["length_m"] > 0.0,
        "beam_positive": values["beam_m"] > 0.0,
        "draft_positive": values["draft_m"] > 0.0,
        "block_coefficient_in_unit_interval": 0.0 < form.block_coefficient < 1.0,
        "volume_consistent_with_LBT": abs(
            form.volume / (values["length_m"] * values["beam_m"] * values["draft_m"])
            - form.block_coefficient
        )
        < 1e-12,
    }
    if declared:
        for name, value in declared.items():
            if name not in values:
                raise NURBSError(
                    f"declared geometry field {name!r} is not part of the gate "
                    f"({sorted(values)})"
                )
            reference = values[name]
            if reference == 0.0:
                raise NURBSError(f"cannot compare declared {name!r}: reference is zero")
            deviation = abs(value - reference) / abs(reference)
            checks[f"declared_{name}_within_1e-3"] = deviation <= 1e-3
            if deviation > 1e-3:
                raise NURBSError(
                    f"declared {name} = {value!r} deviates {deviation:.3e} from the "
                    f"recomputed {reference!r}; refusing to publish the surface block"
                )
    return {"values": values, "checks": checks, "passed": all(checks.values())}


def pc6_compat_block(
    surface: NURBSSurface,
    *,
    parameters: Any = None,
    declared: dict[str, float] | None = None,
    tessellation: Tessellation | None = None,
) -> dict[str, Any]:
    """An **additive** parametric-surface block for P1's ``geom_pc6.json`` contract.

    Compatibility is by construction and by test, not by promise:

    * the block re-uses P1's schema tag (``tensorcad.pc6.geom/1.0``) verbatim, so the
      file keeps one schema and the loader keeps one branch;
    * the block is *additive*: :func:`tensorcad.pc6_contract.load_contract` reads
      ``schema``, ``parameters``, ``artifacts``, ``geometry`` and ``acceptance`` only,
      so an extra key neither changes nor loosens the existing checks -- the test
      suite re-parses the parameters through P1's own ``_parameters_from_contract``
      with the block appended;
    * the geometry gate below recomputes L, B, T and Cb from P1's closed form and
      refuses to publish when a declared value deviates by more than 1e-3.
    """
    from .pc6 import PC6Parameters
    from .pc6_contract import CONTRACT_SCHEMA

    parameters = PC6Parameters() if parameters is None else parameters
    tessellation = tessellation or surface.tessellate(48, 6)
    report = tessellation_report(surface, tessellation)
    gate = _geometry_gate(parameters, declared)
    return {
        "schema": CONTRACT_SCHEMA,
        "representation": NURBS_SCHEMA,
        "additive": True,
        "note": (
            "Parametric (NURBS) representation of one surface of the same body. "
            "Additive to the P1 contract: no existing key is modified or relaxed."
        ),
        "surface": {
            "degree": list(surface.degree),
            "knots": [list(knots) for knots in surface.knots],
            "weights": [list(row) for row in surface.weights],
            "control_net": [[list(point) for point in row] for row in surface.control_net],
            "units": "m",
        },
        "metrics": {
            "parametric_area_m2": report["parametric_area_m2"],
            "tessellated_area_m2": report["tessellated_area_m2"],
            "area_relative_error": report["area_relative_error"],
            "watertight": report["watertight"],
            "max_normal_angle_deg": report["max_normal_angle_deg"],
            "inward_faces": report["inward_faces"],
        },
        "geometry_gate": gate,
        "continuity": surface.continuity(1),
    }


def write_block(block: dict[str, Any], path: str | Path) -> Path:
    """Write a compatibility block as JSON (stable key order, no NaN)."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(block, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return target
