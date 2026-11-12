"""Skinning (lofting) a NURBS surface through an arbitrary sequence of section curves.

Why this exists
---------------
:func:`tensorcad.parametric_hull.loft_nurbs_surface` skins a hull from a sequence of
:class:`~tensorcad.parametric_hull.CrossSection` objects: every section is a ``(y, z)``
polyline living in its own plane ``x = const`` and re-expressed in one half-breadth
frame.  That is a loft *of hull stations*, not a loft: it cannot take a general 3-D
NURBS curve, it cannot take sections that are tilted, non-planar or not mutually
parallel, and it forces every section through a single frame.
``GAP-CAD-016`` in ``docs/oss_inventory_cad/GAP-ANCHORS.md`` asks for the general
operation -- an arbitrary section sequence -- with a machine-decidable acceptance:

* the surface passes through every section (deviation ``< 1e-9 L``), and
* the joins in the loft direction are ``G1``/``G2`` and that is *verifiable*.

The method here is the classical **exact** skin (Piegl & Tiller, ch. 10), not a fit:

1. **Unify the sections.**  Every section is knot-normalised to ``[0, 1]``, degree-
   elevated to the highest degree in the sequence, and every interior knot of the union
   knot vector is inserted into it (Boehm).  All sections then share one degree ``q`` and
   one knot vector, hence one control point count -- and they are still the *same
   geometry*, because knot normalisation, elevation and insertion are exact
   re-representations, not approximations.
2. **Interpolate the control nets.**  For every transverse index ``j`` the ``n_u``
   homogeneous control points ``(w x, w y, w z, w)`` are interpolated along the loft
   direction on the mean chord-length parametrisation (Piegl & Tiller eqs. 9.1--9.8).
   At ``u = u_i`` that interpolation basis is the Kronecker delta, so the surface
   restricted to the row is *the section itself*: same degree, same knots, same weights.
3. **One patch, no stitching.**  With ``degree_u = 3`` and simple interior knots the
   skin is ``C2`` in the loft direction.  :meth:`Loft.measure` checks that two-sided at
   every interior knot, and the same measurement applied to a deliberately reduced
   ``degree_u = 2`` skin reports the finite second-derivative jump -- so the check can
   fail, it is not a tautology.

Why exactness matters: a sampled/fitted skin through exact rational circles carries a
``O(h^4)`` interpolation error between the sample rows; the exact skin *is* the cylinder,
so the measured radial deviation is round-off (``docs/loft.md``).

What is *not* done
------------------
No section reparametrisation (each keeps its own ``v`` parametrisation, so the knot union
is a representation union, not a geometric correspondence), no self-intersection or fold
detection, no trimming.  A section with a genuine break (interior knot of multiplicity
``degree + 1``) is rejected with a named reason instead of being silently elevated --
split it with :meth:`tensorcad.nurbs.NURBSCurve.split_at_discontinuities` first.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

from .nurbs import NURBSCurve, NURBSSurface, derivative_net, find_span
from .parametric_hull import interpolate_control_points

Vec3 = tuple[float, float, float]
Homogeneous = tuple[float, float, float, float]

__all__ = [
    "Loft",
    "LoftError",
    "LoftMeasurement",
    "knot_union",
    "loft_sections",
    "loft_surface",
    "section_deviation",
]


class LoftError(ValueError):
    """Raised when a section sequence cannot be skinned (named reason in the message)."""


# ---------------------------------------------------------------------------
# Vector helpers
# ---------------------------------------------------------------------------
def _sub(a: Sequence[float], b: Sequence[float]) -> tuple[float, ...]:
    return tuple(x - y for x, y in zip(a, b, strict=True))


def _norm(a: Sequence[float]) -> float:
    return math.sqrt(sum(value * value for value in a))


def _distance(a: Sequence[float], b: Sequence[float]) -> float:
    return _norm(_sub(a, b))


def _homogeneous(curve: NURBSCurve) -> tuple[Homogeneous, ...]:
    return tuple(
        (weight * point[0], weight * point[1], weight * point[2], weight)
        for point, weight in zip(curve.control_points, curve.weights, strict=True)
    )


def _dehomogenise(point: Homogeneous) -> Vec3:
    weight = point[3]
    if weight <= 0.0:
        raise LoftError(f"non-positive weight {weight!r} in the control net")
    return (point[0] / weight, point[1] / weight, point[2] / weight)


# ---------------------------------------------------------------------------
# Exact re-representation: normalisation, knot insertion, degree elevation
# ---------------------------------------------------------------------------
def _normalised_knots(knots: Sequence[float]) -> tuple[float, ...]:
    """Clamp-normalise a knot vector to ``[0, 1]`` (a parameter change, not a shape change)."""
    if knots[0] == 0.0 and knots[-1] == 1.0:
        return tuple(knots)
    low, high = knots[0], knots[-1]
    span = high - low
    if span <= 0.0:
        raise LoftError(f"degenerate knot span {low!r}..{high!r}")
    return tuple((value - low) / span for value in knots)


def _multiplicity(knots: Sequence[float], value: float) -> int:
    return sum(1 for knot in knots if knot == value)


def _insert_knot_once(
    degree: int,
    knots: Sequence[float],
    points: Sequence[Homogeneous],
    u: float,
) -> tuple[tuple[float, ...], tuple[Homogeneous, ...]]:
    """Boehm's single knot insertion in homogeneous coordinates (Piegl & Tiller A5.1).

    The refinement is exact: the inserted knot adds a control point without moving the
    curve, so a rational section stays the same geometry.
    """
    p = degree
    n = len(points) - 1
    if not (knots[p] < u < knots[n + 1]):
        raise LoftError(f"cannot insert knot {u!r}: not strictly inside the domain")
    span = find_span(knots, p, u)
    s = _multiplicity(knots, u)
    if s >= p:
        raise LoftError(f"knot {u!r} already has multiplicity {s} > degree {p}")

    out: list[Homogeneous | None] = [None] * (n + 2)
    for i in range(0, span - p + 1):
        out[i] = points[i]
    for i in range(span - s, n + 1):
        out[i + 1] = points[i]
    reduced = [list(points[i]) for i in range(span - p, span - s + 1)]
    left = span - p + 1
    for i in range(0, p - s):
        denominator = knots[i + span + 1] - knots[left + i]
        if denominator == 0.0:
            raise LoftError(f"degenerate denominator inserting knot {u!r}")
        alpha = (u - knots[left + i]) / denominator
        previous, following = reduced[i], reduced[i + 1]
        reduced[i] = [
            (1.0 - alpha) * previous[axis] + alpha * following[axis] for axis in range(4)
        ]
    out[left] = tuple(reduced[0])
    # A5.1 with r - s = 1 (a single insertion): Qw[k + r - j - s] = Rw[p - j - s]
    out[span - s] = tuple(reduced[p - 1 - s])
    for i in range(left + 1, span - s):
        out[i] = tuple(reduced[i - left])
    for i in range(0, left):
        out[i] = points[i]
    for i in range(span - s + 1, n + 2):
        out[i] = points[i - 1]

    filled: list[Homogeneous] = []
    for index, value in enumerate(out):
        if value is None:  # pragma: no cover - defensive
            raise LoftError(f"knot insertion left control point {index} undefined")
        filled.append(value)
    new_knots = list(knots)
    new_knots.insert(span + 1, u)
    return tuple(new_knots), tuple(filled)


def _insert_knot(
    degree: int,
    knots: Sequence[float],
    points: Sequence[Homogeneous],
    u: float,
    times: int,
) -> tuple[tuple[float, ...], tuple[Homogeneous, ...]]:
    current_knots, current_points = tuple(knots), tuple(points)
    for _ in range(times):
        current_knots, current_points = _insert_knot_once(
            degree, current_knots, current_points, u
        )
    return current_knots, current_points


def _interior_breaks(knots: Sequence[float], degree: int) -> list[float]:
    """Interior knot values whose multiplicity is ``degree + 1`` (a genuine break)."""
    interior = sorted({value for value in knots[degree + 1 : -(degree + 1)]})
    return [value for value in interior if _multiplicity(knots, value) >= degree + 1]


def _elevate_section(curve: NURBSCurve, times: int) -> NURBSCurve:
    if times == 0:
        return curve
    try:
        return curve.elevate_degree(times)
    except Exception as error:  # pragma: no cover - surfaced as LoftError
        raise LoftError(f"degree elevation failed: {error}") from error


def knot_union(
    sections: Sequence[NURBSCurve], *, degree: int | None = None
) -> tuple[int, tuple[float, ...]]:
    """Common ``(degree, knot vector)`` of the exact skin of ``sections``.

    The degree is the highest section degree (or ``degree`` when larger); the knot
    vector is the union of the normalised section knot vectors with the maximum
    multiplicity of each interior value.  Inserting the difference into every section
    makes them share this representation exactly.
    """
    if not sections:
        raise LoftError("no sections")
    target = max(curve.degree for curve in sections)
    if degree is not None:
        if degree < 1:
            raise LoftError(f"degree must be >= 1, got {degree}")
        target = max(target, degree)
    wanted: dict[float, int] = {}
    for index, curve in enumerate(sections):
        elevated = _elevate_section(curve, target - curve.degree)
        knots = _normalised_knots(elevated.knots)
        breaks = _interior_breaks(knots, elevated.degree)
        if breaks:
            raise LoftError(
                f"section {index} has a break at v={breaks[0]!r} (multiplicity "
                f"{_multiplicity(knots, breaks[0])} = degree + 1): split it with "
                "NURBSCurve.split_at_discontinuities() before lofting"
            )
        for value in knots[elevated.degree + 1 : -(elevated.degree + 1)]:
            wanted[value] = max(wanted.get(value, 0), _multiplicity(knots, value))
    combined: list[float] = [0.0] * (target + 1)
    for value in sorted(wanted):
        combined.extend([value] * wanted[value])
    combined.extend([1.0] * (target + 1))
    return target, tuple(combined)


def _unified_sections(
    sections: Sequence[NURBSCurve], degree: int, knots: tuple[float, ...]
) -> tuple[tuple[Homogeneous, ...], ...]:
    """Every section re-expressed exactly on ``(degree, knots)``."""
    unified: list[tuple[Homogeneous, ...]] = []
    for index, curve in enumerate(sections):
        section = _elevate_section(curve, degree - curve.degree)
        section_knots = _normalised_knots(section.knots)
        points = _homogeneous(section)
        new_knots = section_knots
        for value in knots:
            if value in (0.0, 1.0):
                continue
            missing = _multiplicity(knots, value) - _multiplicity(new_knots, value)
            if missing > 0:
                new_knots, points = _insert_knot(degree, new_knots, points, value, missing)
        if tuple(new_knots) != tuple(knots):
            raise LoftError(
                f"section {index} could not be brought onto the union knot vector"
            )
        expected = len(knots) - degree - 1
        if len(points) != expected:
            raise LoftError(
                f"section {index} has {len(points)} control points, expected {expected}"
            )
        unified.append(points)
    return tuple(unified)


# ---------------------------------------------------------------------------
# Exact one-sided u-derivatives (independent of the v-direction knot multiplicity)
# ---------------------------------------------------------------------------
def _basis_row(
    knots: Sequence[float], degree: int, u: float, n_ctrl: int, span: int | None = None
) -> list[float]:
    """Non-zero B-spline basis values at ``u``, written into a full-length row (P&T A2.1).

    ``span`` forces the polynomial piece, which is how the one-sided limit at an interior
    knot is taken: ``find_span`` returns the piece to the *right* of the knot, and the
    piece to the left is one less.
    """
    if span is None:
        span = find_span(knots, degree, u)
    span = max(degree, min(span, n_ctrl - 1))
    values = [0.0] * (degree + 1)
    values[0] = 1.0
    left = [0.0] * (degree + 1)
    right = [0.0] * (degree + 1)
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
    row = [0.0] * n_ctrl
    for j in range(degree + 1):
        row[span - degree + j] = values[j]
    return row


def _eval_homogeneous_net(
    net: Sequence[Sequence[Homogeneous]],
    degree: tuple[int, int],
    knots: tuple[Sequence[float], Sequence[float]],
    u: float,
    v: float,
    *,
    span_u: int | None = None,
) -> Homogeneous:
    row_u = _basis_row(knots[0], degree[0], u, len(net), span=span_u)
    row_v = _basis_row(knots[1], degree[1], v, len(net[0]))
    total = [0.0, 0.0, 0.0, 0.0]
    for i, weight_u in enumerate(row_u):
        if weight_u == 0.0:
            continue
        for j, weight_v in enumerate(row_v):
            if weight_v == 0.0:
                continue
            factor = weight_u * weight_v
            point = net[i][j]
            for axis in range(4):
                total[axis] += factor * point[axis]
    return (total[0], total[1], total[2], total[3])


def _u_derivatives(
    surface: NURBSSurface, u: float, v: float, *, span_u: int | None = None
) -> tuple[Vec3, Vec3, Vec3 | None]:
    """``(S, dS/du, d2S/du2)`` exactly, by differentiating the homogeneous net in ``u``.

    ``NURBSSurface.derivatives`` differentiates both directions up to the requested
    order, which a section-representation union can forbid in ``v`` (an exact circle has
    interior knots of multiplicity 2, so its *second* ``v``-derivative net is not a valid
    B-spline net).  The loft direction is where the continuity claim lives, so it is
    differentiated on its own here -- with the same exact rational quotient recursion
    (Piegl & Tiller A3.2).
    """
    degree_u, degree_v = surface.degree
    knots_u, knots_v = surface.knots
    net = tuple(
        tuple(
            (
                surface.weights[i][j] * surface.control_net[i][j][0],
                surface.weights[i][j] * surface.control_net[i][j][1],
                surface.weights[i][j] * surface.control_net[i][j][2],
                surface.weights[i][j],
            )
            for j in range(len(surface.control_net[i]))
        )
        for i in range(len(surface.control_net))
    )

    def differentiate(current: tuple[tuple[Homogeneous, ...], ...], degree: int,
                      knots: Sequence[float]) -> tuple[tuple[tuple[Homogeneous, ...], ...], int,
                                                        tuple[float, ...]]:
        columns = [
            derivative_net([current[i][j] for i in range(len(current))], degree, knots)[0]
            for j in range(len(current[0]))
        ]
        new_knots = tuple(knots[1:-1])
        new_net = tuple(
            tuple(columns[j][i] for j in range(len(current[0])))
            for i in range(len(columns[0]))
        )
        return new_net, degree - 1, new_knots

    first, degree_1, knots_1 = differentiate(net, degree_u, knots_u)
    if degree_1 >= 1:
        second, degree_2, knots_2 = differentiate(first, degree_1, knots_1)
    else:
        second, degree_2, knots_2 = None, -1, None

    a0 = _eval_homogeneous_net(
        net, (degree_u, degree_v), (knots_u, knots_v), u, v, span_u=span_u
    )
    a1 = _eval_homogeneous_net(
        first, (degree_1, degree_v), (knots_1, knots_v), u, v,
        span_u=None if span_u is None else span_u - 1,
    )
    if second is not None:
        a2 = _eval_homogeneous_net(
            second, (degree_2, degree_v), (knots_2, knots_v), u, v,
            span_u=None if span_u is None else span_u - 2,
        )
    else:
        a2 = None

    weight = a0[3]
    if weight <= 0.0:
        raise LoftError(f"degenerate evaluation at (u={u!r}, v={v!r}): weight {weight!r}")
    position = (a0[0] / weight, a0[1] / weight, a0[2] / weight)
    weight_u = a1[3]
    first_vec = tuple((a1[axis] - weight_u * position[axis]) / weight for axis in range(3))
    if a2 is None:
        return position, first_vec, None
    weight_uu = a2[3]
    second_vec = tuple(
        (a2[axis] - 2.0 * weight_u * first_vec[axis] - weight_uu * position[axis]) / weight
        for axis in range(3)
    )
    return position, first_vec, second_vec


# ---------------------------------------------------------------------------
# Measurement
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class LoftMeasurement:
    """What the returned surface actually does, measured on the surface itself."""

    n_sections: int
    n_v: int
    degree: tuple[int, int]
    scale: float
    u_parameters: tuple[float, ...]
    interior_knots_u: tuple[float, ...]
    interpolation_defect: tuple[float, ...]
    deviation: tuple[float, ...]
    relative_interpolation_defect: float
    relative_deviation: float
    second_derivative_jump: tuple[float, ...]
    first_derivative_jump: tuple[float, ...]
    relative_second_derivative_jump: float
    relative_first_derivative_jump: float
    measured_continuous_order: int

    def to_dict(self) -> dict[str, object]:
        return {
            "n_sections": self.n_sections,
            "n_v": self.n_v,
            "degree": list(self.degree),
            "scale": self.scale,
            "u_parameters": list(self.u_parameters),
            "interior_knots_u": list(self.interior_knots_u),
            "interpolation_defect": list(self.interpolation_defect),
            "deviation": list(self.deviation),
            "relative_interpolation_defect": self.relative_interpolation_defect,
            "relative_deviation": self.relative_deviation,
            "second_derivative_jump": list(self.second_derivative_jump),
            "first_derivative_jump": list(self.first_derivative_jump),
            "relative_second_derivative_jump": self.relative_second_derivative_jump,
            "relative_first_derivative_jump": self.relative_first_derivative_jump,
            "measured_continuous_order": self.measured_continuous_order,
        }


def _nearest_on_isocurve(
    surface: NURBSSurface, u: float, point: Sequence[float], *, grid: int = 96
) -> float:
    """Distance from ``point`` to ``S(u, .)`` by a grid scan plus golden-section refine."""
    v0, v1 = surface.domain[1]
    best, best_v = math.inf, v0
    for index in range(grid + 1):
        v = v0 + (v1 - v0) * index / grid
        distance = _distance(point, surface.evaluate(u, v))
        if distance < best:
            best, best_v = distance, v
    step = (v1 - v0) / grid
    low = max(v0, best_v - step)
    high = min(v1, best_v + step)
    ratio = 0.5 * (math.sqrt(5.0) - 1.0)
    for _ in range(80):
        a = high - ratio * (high - low)
        b = low + ratio * (high - low)
        if _distance(point, surface.evaluate(u, a)) < _distance(point, surface.evaluate(u, b)):
            high = b
        else:
            low = a
    for candidate in (low, high, 0.5 * (low + high)):
        best = min(best, _distance(point, surface.evaluate(u, candidate)))
    return best


def section_deviation(
    surface: NURBSSurface, u: float, curve: NURBSCurve, *, samples: int = 65
) -> tuple[float, float]:
    """``(max, mean)`` distance from ``curve`` to the isoparametric line ``S(u, .)``.

    Measured against a dense sampling of the *curve* and a refined point-to-curve search
    on the surface, so this is a geometric statement, not a re-evaluation of the data the
    solve was handed.
    """
    v0, v1 = curve.domain
    worst, total = 0.0, 0.0
    for index in range(samples):
        v = v0 + (v1 - v0) * index / (samples - 1)
        distance = _nearest_on_isocurve(surface, u, curve.evaluate(v))
        worst = max(worst, distance)
        total += distance
    return worst, total / samples


# ---------------------------------------------------------------------------
# The loft
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Loft:
    """A skinned patch plus the data it was built from (so it can be measured)."""

    surface: NURBSSurface
    sections: tuple[NURBSCurve, ...]
    control_rows: tuple[tuple[Homogeneous, ...], ...]
    u_parameters: tuple[float, ...]

    @property
    def n_sections(self) -> int:
        return len(self.sections)

    @property
    def n_v(self) -> int:
        return len(self.control_rows[0])

    def _scale(self) -> float:
        low = [math.inf] * 3
        high = [-math.inf] * 3
        for row in self.control_rows:
            for point in row:
                position = _dehomogenise(point)
                for axis in range(3):
                    low[axis] = min(low[axis], position[axis])
                    high[axis] = max(high[axis], position[axis])
        return math.sqrt(sum((high[axis] - low[axis]) ** 2 for axis in range(3)))

    def measure(self, *, samples: int = 65, deviation_samples: int = 33) -> LoftMeasurement:
        """Measure the exact-interpolation defect, geometric deviation and continuity."""
        surface = self.surface
        degree_u = surface.degree[0]
        scale = self._scale()
        if scale <= 0.0:
            raise LoftError("degenerate section set: zero bounding-box diagonal")

        defects: list[float] = []
        deviations: list[float] = []
        for index, curve in enumerate(self.sections):
            u = self.u_parameters[index]
            v0, v1 = curve.domain
            defects.append(
                max(
                    _distance(
                        curve.evaluate(v0 + (v1 - v0) * j / (samples - 1)),
                        surface.evaluate(u, j / (samples - 1)),
                    )
                    for j in range(samples)
                )
            )
            worst_deviation, _ = section_deviation(
                surface, u, curve, samples=deviation_samples
            )
            deviations.append(worst_deviation)

        knots_u = surface.knots[0]
        interior = sorted({value for value in knots_u[degree_u + 1 : -(degree_u + 1)]})
        v_mid = 0.5 * (surface.domain[1][0] + surface.domain[1][1])
        second_jumps: list[float] = []
        first_jumps: list[float] = []
        for knot in interior:
            # one-sided limits: the piece right of the knot is what find_span returns,
            # the piece left of it is one span earlier.  Taking them exactly (no delta)
            # is what makes a C2 patch read as round-off instead of O(delta).
            right_span = find_span(knots_u, degree_u, knot)
            _, left_first, left_second = _u_derivatives(
                surface, knot, v_mid, span_u=right_span - 1
            )
            _, right_first, right_second = _u_derivatives(
                surface, knot, v_mid, span_u=right_span
            )
            first_jumps.append(_distance(left_first, right_first))
            if left_second is not None and right_second is not None:
                second_jumps.append(_distance(left_second, right_second))

        tolerance = 1e-9 * scale
        if not interior:
            measured = degree_u  # no join to fail
        elif second_jumps and all(jump <= tolerance for jump in second_jumps):
            measured = 2
        elif all(jump <= tolerance for jump in first_jumps):
            measured = 1
        else:
            measured = 0
        measured = min(measured, degree_u)

        return LoftMeasurement(
            n_sections=self.n_sections,
            n_v=self.n_v,
            degree=surface.degree,
            scale=scale,
            u_parameters=self.u_parameters,
            interior_knots_u=tuple(interior),
            interpolation_defect=tuple(defects),
            deviation=tuple(deviations),
            relative_interpolation_defect=max(defects) / scale,
            relative_deviation=max(deviations) / scale,
            second_derivative_jump=tuple(second_jumps),
            first_derivative_jump=tuple(first_jumps),
            relative_second_derivative_jump=(
                max(second_jumps) / scale if second_jumps else 0.0
            ),
            relative_first_derivative_jump=(
                max(first_jumps) / scale if first_jumps else 0.0
            ),
            measured_continuous_order=measured,
        )


def _mean_chord_parameters(columns: Sequence[Sequence[Sequence[float]]]) -> tuple[float, ...]:
    """Mean normalised chord-length parametrisation shared by the control columns."""
    n = len(columns[0])
    if n < 2:
        raise LoftError("need at least 2 sections to parametrise the loft direction")
    accumulator = [0.0] * n
    counted = 0
    for column in columns:
        if len(column) != n:
            raise LoftError("control columns must have equal length")
        cumulative = [0.0]
        for index in range(1, n):
            cumulative.append(
                cumulative[-1] + _distance(column[index][:3], column[index - 1][:3])
            )
        total = cumulative[-1]
        if total <= 1e-30:
            continue
        for index in range(n):
            accumulator[index] += cumulative[index] / total
        counted += 1
    if counted == 0:
        return tuple(index / (n - 1) for index in range(n))
    return tuple(value / counted for value in accumulator)


def _averaged_knots(degree: int, params: Sequence[float]) -> tuple[float, ...]:
    """Averaging knot vector ``U_{j+p} = (u_{j-1} + ... + u_{j+p-2}) / p`` (Piegl 9.8)."""
    n_ctrl = len(params)
    if n_ctrl < degree + 1:
        raise LoftError(f"{n_ctrl} data points insufficient for degree {degree}")
    last = n_ctrl - 1
    knots = [0.0] * (degree + 1)
    for j in range(1, last - degree + 1):
        knots.append(sum(params[j + i - 1] for i in range(1, degree + 1)) / degree)
    knots.extend([1.0] * (degree + 1))
    return tuple(knots)


def loft_sections(
    sections: Sequence[NURBSCurve], *, degree_u: int = 3, degree_v: int | None = None
) -> Loft:
    """Skin a NURBS surface **through** ``sections`` (ordered along the loft direction).

    Parameters
    ----------
    sections : sequence of NURBSCurve
        Arbitrary 3-D rational curves.  They need not be planar, parallel, equally
        spaced, equally parametrised or of equal degree.
    degree_u : int
        Degree along the loft direction (clamped by ``n_sections - 1``; a loft of two
        sections is ruled whatever is asked for, and ``Loft.surface.degree`` says so).
    degree_v : int, optional
        Minimum transverse degree; the actual transverse degree is the highest section
        degree elevated to at least this value.

    Raises
    ------
    LoftError
        Fewer than two sections, a non-curve section, a non-positive degree, two
        coincident consecutive sections (degenerate strip), or a section with a genuine
        break.  The message names the failing index.
    """
    if len(sections) < 2:
        raise LoftError(f"need at least 2 sections to loft, got {len(sections)}")
    if degree_u < 1:
        raise LoftError(f"degree_u must be >= 1, got {degree_u}")
    for index, section in enumerate(sections):
        if not isinstance(section, NURBSCurve):
            raise LoftError(
                f"section {index} is {type(section).__name__}, expected "
                "tensorcad.nurbs.NURBSCurve"
            )

    transverse_degree, knots_v = knot_union(sections, degree=degree_v)
    rows = _unified_sections(sections, transverse_degree, knots_v)
    n_u = len(rows)
    n_v = len(rows[0])
    degree_u_actual = min(degree_u, n_u - 1)

    scale = 0.0
    for row in rows:
        for point in row:
            scale = max(scale, _norm(_dehomogenise(point)))
    if scale <= 0.0:
        raise LoftError("degenerate section set: zero extent")
    for index in range(n_u - 1):
        gap = max(_distance(rows[index][j], rows[index + 1][j]) for j in range(n_v))
        if gap <= 1e-14 * scale:
            raise LoftError(
                f"sections {index} and {index + 1} coincide (max control-net gap "
                f"{gap!r}); the loft strip is degenerate -- drop the duplicate"
            )

    # A skin of polynomial sections is polynomial.  Solving for the homogeneous weights
    # as well would return 1 +/- 1e-16 and turn an exactly polynomial patch into a
    # rational quotient of two nearly-equal sums -- so the plain 3-D solve is used when
    # every section is polynomial, and the weights are exactly one.
    polynomial = not any(section.rational for section in sections)
    if polynomial:
        columns = [
            [tuple(point[:3]) for point in (rows[i][j] for i in range(n_u))]
            for j in range(n_v)
        ]
    else:
        columns = [[rows[i][j] for i in range(n_u)] for j in range(n_v)]
    u_parameters = _mean_chord_parameters(columns)
    if any(
        u_parameters[index + 1] - u_parameters[index] <= 0.0
        for index in range(len(u_parameters) - 1)
    ):
        raise LoftError(f"loft parameters are not strictly increasing: {u_parameters!r}")
    knots_u = _averaged_knots(degree_u_actual, u_parameters)
    column_control = [
        interpolate_control_points(
            [tuple(value for value in point) for point in column],
            degree_u_actual,
            params=u_parameters,
            knots=knots_u,
        )[0]
        for column in columns
    ]

    if polynomial:
        control_net = tuple(
            tuple(tuple(column_control[j][i]) for j in range(n_v)) for i in range(n_u)
        )
        weights = tuple(tuple(1.0 for _ in range(n_v)) for _ in range(n_u))
    else:
        control_net = tuple(
            tuple(_dehomogenise(column_control[j][i]) for j in range(n_v))
            for i in range(n_u)
        )
        weights = tuple(
            tuple(float(column_control[j][i][3]) for j in range(n_v)) for i in range(n_u)
        )
    surface = NURBSSurface(
        degree=(degree_u_actual, transverse_degree),
        knots=(knots_u, knots_v),
        control_net=control_net,
        weights=weights,
    )
    return Loft(
        surface=surface,
        sections=tuple(sections),
        control_rows=rows,
        u_parameters=u_parameters,
    )


def loft_surface(
    sections: Sequence[NURBSCurve], *, degree_u: int = 3, degree_v: int | None = None
) -> NURBSSurface:
    """Convenience wrapper around :func:`loft_sections` returning only the surface."""
    return loft_sections(sections, degree_u=degree_u, degree_v=degree_v).surface