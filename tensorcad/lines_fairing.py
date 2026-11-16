"""Three-direction fairing of an offset table: curvature energy, not eyeballing.

A ship's lines are handed over as a table of half-breadths ``y(x, z)`` on a
``station x waterline`` grid (``tensorcad.offset_loft.ShipOffsetTable``).  Before that
table can be lofted, ballasted or meshed it has to be **fair**: the designer's three
families of curves -- waterlines (``z`` const), stations (``x`` const) and buttocks
(``y`` const, the vertical-plane lines) -- must be smooth, which in practice means
their curvature may change sign only where the shape really has an inflection.

``GAP-CAD-C3-LINES-FAIRING`` asks for that operation with a machine-decidable
acceptance::

    fairing deviation < the given tolerance
    AND the curvature sign-flip count goes monotonically down, on >= 2 ship tables
    AND agreement with a manual-fairing reference

This module implements it as a **variational problem on the three families at once**.
Write the table in non-dimensional coordinates ``xi = x / L``, ``zeta = z / T``,
``eta = 2 y / B`` and let ``eta0`` be the tabulated ordinates.  The faired table
``eta = eta0 + delta`` minimises

.. code-block:: text

    J(delta) = sum_ij delta_ij^2
             + lam_wl * || D_xx eta ||^2_Q      waterlines   (curves z = const)
             + lam_st * || D_zz eta ||^2_Q      stations     (curves x = const)
             + lam_bt * || D_tt eta ||^2_Q      buttocks     (curves y = const)

The three penalty terms are the curvature energy of the three families, and the reason
they can be written on one grid is the level-set identity

.. code-block:: text

    kappa_buttock = -eta_tt / |grad eta|          (t = tangent of the level curve)

so the curvature of a buttock is ``eta``'s second derivative *along the buttock
tangent*, which :class:`_TangentOperator` discretises with bilinear interpolation at
``p +- h t``.  ``|grad eta| > 0`` away from the extrema, so **sign flips of ``eta_tt``
and sign flips of the buttock's geometric curvature are the same count** -- which makes
the objective and the acceptance metric the same quantity instead of two things that
merely correlate.  ``tests/test_lines_fairing.py`` pins that identity on an analytic
surface, and the study reports the flips measured on the *contoured* buttock polylines,
not on the operator.

Three design decisions worth stating up front:

* **Every reported number describes the delivered table.**  An outer trial that raises the
  flip count is rejected, so the last trial can describe ordinates that were discarded;
  :attr:`FairingResult.flips_after` is therefore re-measured on ``table_after`` and the
  rejected trial stays in the audit trail (``iterations``, each with ``accepted``).
* **The tolerance is a budget, not a check.**  ``lam`` is solved for by bisection so
  that the correction uses ``98 %`` of the deviation budget and never exceeds it; the
  final correction is additionally scaled down if the bisection landed high
  (:attr:`FairingIteration.contraction`).  A fairing that quietly walks the table away
  from the ship is worthless, so the deviation is bounded by construction.
* **Fixed edges.**  The perpendiculars (``xi = 0, 1``), the keel (``zeta = 0``) and
  every ordinate that is already ``0`` are frozen: those values belong to the principal
  dimensions, not to the fairing.  The watertight closure of the sections downstream
  depends on it.
* **A buttock is a level set, and a level set can have several branches.** With
  tumblehome -- any section whose half-breadth falls again above its maximum breadth --
  a level ``eta = const`` is crossed twice, and the two crossings are two different
  buttocks.  :func:`buttock_curves` traces the crossings station to station and returns
  one polyline per branch, splitting where a level set disappears for a stretch; the
  levels themselves come from the table (:func:`buttock_levels`) rather than being fixed
  fractions of the global peak, because a near-deck level spans only a few stations and
  used to drop out of the flip count altogether.  Both are reported per level in
  :func:`buttock_coverage`, so the third family cannot shrink to nothing unnoticed.
* **Monotonicity is enforced, and rejections are reported.**  An outer iteration is
  kept only if the total sign-flip count does not rise; a rejected step leaves the
  table untouched and increments :attr:`FairingResult.rejected_steps`.  The
  per-iteration table in :meth:`FairingResult.to_dict` is the evidence, so a run that
  never improves cannot masquerade as a converging one.

Standard library only, like the rest of the geometry stack (the solver is a
conjugate-gradient iteration on the normal equations; the grid is a few hundred nodes,
so the matrix is never formed).
"""
from __future__ import annotations

import argparse
import json
import math
import random
import sys
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .offset_loft import ShipOffsetTable, table_integrals

SCHEMA_REPORT = "tensorcad.lines_fairing.report/1"

#: Direction keys, in the order the objective sums them.
DIRECTIONS: tuple[str, str, str] = ("waterlines", "stations", "buttocks")

#: Relative threshold below which a discrete curvature is treated as zero when the
#: sign-flip count is taken.  Without it round-off at an inflection counts as a flip.
_SIGN_EPS = 1e-3

#: Fraction of the deviation budget the bisection aims for (leaves head-room for the
#: exact solution to be scaled back instead of overrunning the tolerance).
_BUDGET_FRACTION = 0.98

#: Share of the zeta-direction bending energy at which one waterline row is read as a
#: *design feature* rather than as unfairness.  On a full ``zeta**(1/m)`` section the
#: near-keel row owns 99-100 % of that energy (the section is nearly vertical there, so
#: its second difference in ``zeta`` is singular), while on a parabolic hull no row owns
#: more than 10 %.  The gap is two orders of magnitude wide, so the threshold does not
#: need tuning -- see :func:`feature_waterlines`.
_FEATURE_SHARE = 0.5


class FairingError(ValueError):
    """The table or the request cannot be faired as asked."""


# ---------------------------------------------------------------------------
# small dense helpers -- the whole problem is a few hundred unknowns
# ---------------------------------------------------------------------------


def _spacings(t: Sequence[float]) -> list[float]:
    return [t[i + 1] - t[i] for i in range(len(t) - 1)]


def second_differences(values: Sequence[float], t: Sequence[float]) -> list[float]:
    """Non-uniform three-point second derivative ``d2v/dt2`` at the interior nodes.

    ``2 ((v[i+1]-v[i])/h[i+1] - (v[i]-v[i-1])/h[i]) / (h[i] + h[i+1])`` -- exact for a
    quadratic on any grid.  The end values are ``0``: their curvature is not defined by
    three points, and the callers freeze those nodes anyway.
    """
    n = len(values)
    if n != len(t):
        raise FairingError(f"second_differences: {n} values against {len(t)} abscissae")
    if n < 3:
        return [0.0] * n
    h = _spacings(t)
    out = [0.0] * n
    for i in range(1, n - 1):
        hi, hi1 = h[i - 1], h[i]
        out[i] = 2.0 * (
            (values[i + 1] - values[i]) / hi1 - (values[i] - values[i - 1]) / hi
        ) / (hi + hi1)
    return out


def curvature_sign_flips(
    values: Sequence[float],
    t: Sequence[float],
    *,
    relative_floor: float = _SIGN_EPS,
) -> int:
    """Number of curvature sign changes along one curve.

    The curvature is taken as :func:`second_differences`; for a planar graph that is the
    curvature up to the positive factor ``(1 + v'^2)**-1.5``, so the *signs* -- and
    therefore the flip count -- are exact, and so is the count on a parametric curve
    whose parameter is monotone.  ``relative_floor`` expresses the dead band as a
    fraction of the curve's own curvature scale, so the count is invariant under a
    change of units.
    """
    n = len(values)
    if n < 3:
        return 0
    d2 = second_differences(values, t)
    span = t[-1] - t[0]
    scale = max(abs(v) for v in values) / (span * span) if span > 0.0 else 0.0
    floor = relative_floor * scale
    flips = 0
    previous = 0
    for value in d2:
        sign = 0 if abs(value) <= floor else (1 if value > 0.0 else -1)
        if sign == 0:
            continue
        if previous != 0 and sign != previous:
            flips += 1
        previous = sign
    return flips


def _cross(a: Sequence[float], b: Sequence[float]) -> float:
    return a[0] * b[1] - a[1] * b[0]


def polyline_sign_flips(
    points: Sequence[Sequence[float]], *, relative_floor: float = _SIGN_EPS
) -> int:
    """Curvature sign flips of a 2-D polyline, via the turning cross product.

    Used for the buttocks, where the curve is a *level set* extracted by contouring:
    the discrete turning sense ``cross(p[i]-p[i-1], p[i+1]-p[i])`` has the sign of the
    curvature, so counting its changes measures the same thing as
    :func:`curvature_sign_flips` without needing a parameterisation.
    """
    turns: list[float] = []
    for i in range(1, len(points) - 1):
        first = (points[i][0] - points[i - 1][0], points[i][1] - points[i - 1][1])
        second = (points[i + 1][0] - points[i][0], points[i + 1][1] - points[i][1])
        length_first = math.hypot(*first)
        length_second = math.hypot(*second)
        if length_first == 0.0 or length_second == 0.0:
            continue
        turns.append(_cross(first, second) / (length_first * length_second))
    if not turns:
        return 0
    floor = relative_floor * max(abs(v) for v in turns)
    flips = 0
    previous = 0
    for turn in turns:
        sign = 0 if abs(turn) <= floor else (1 if turn > 0.0 else -1)
        if sign == 0:
            continue
        if previous != 0 and sign != previous:
            flips += 1
        previous = sign
    return flips


def _conjugate_gradient(
    matvec: Any,
    rhs: Sequence[float],
    start: Sequence[float] | None = None,
    *,
    tolerance: float = 1e-10,
    max_iterations: int = 400,
) -> tuple[list[float], int, float]:
    """CG on a symmetric positive-definite system; returns ``(x, iterations, residual)``."""
    n = len(rhs)
    x = list(start) if start is not None else [0.0] * n
    ax = matvec(x)
    r = [rhs[i] - ax[i] for i in range(n)]
    p = list(r)
    rr = math.fsum(v * v for v in r)
    reference = math.sqrt(math.fsum(v * v for v in rhs)) or 1.0
    if math.sqrt(rr) <= tolerance * reference:
        return x, 0, math.sqrt(rr) / reference
    iterations = 0
    for iterations in range(1, max_iterations + 1):
        ap = matvec(p)
        pap = math.fsum(pi * ai for pi, ai in zip(p, ap))
        if not (pap > 0.0):
            break
        alpha = rr / pap
        for i in range(n):
            x[i] += alpha * p[i]
            r[i] -= alpha * ap[i]
        rr_new = math.fsum(v * v for v in r)
        if math.sqrt(rr_new) <= tolerance * reference:
            rr = rr_new
            break
        beta = rr_new / rr
        rr = rr_new
        for i in range(n):
            p[i] = r[i] + beta * p[i]
    return x, iterations, math.sqrt(rr) / reference


# ---------------------------------------------------------------------------
# the three curvature operators, all on one flat grid
# ---------------------------------------------------------------------------


def _divided_difference_coefficients(
    abscissa: Sequence[float], index: int, order: int
) -> list[tuple[int, float]] | None:
    """Coefficients of ``order! * f[x_{i-k}, ..., x_{i+1}]`` -- the ``order``-th derivative.

    Newton's divided difference is built as a linear combination of the stencil values,
    so the same code gives the non-uniform weights the clustered waterline grids need.
    The stencil is ``index - order + 1 .. index + 1`` (backward-biased, so the operator
    stays inside the grid whenever the index does).  ``None`` when the stencil leaves the
    curve -- the caller then contributes no row rather than a one-sided guess.
    """
    lo = index - order + 1
    hi = index + 1
    if lo < 0 or hi > len(abscissa) - 1:
        return None
    nodes = list(range(lo, hi + 1))

    def weights(columns: Sequence[int]) -> dict[int, float]:
        if len(columns) == 2:
            a, b = columns
            gap = abscissa[b] - abscissa[a]
            return {a: -1.0 / gap, b: 1.0 / gap}
        head = weights(columns[1:])
        tail = weights(columns[:-1])
        scale = abscissa[columns[-1]] - abscissa[columns[0]]
        combined: dict[int, float] = {}
        for key, value in head.items():
            combined[key] = combined.get(key, 0.0) + value / scale
        for key, value in tail.items():
            combined[key] = combined.get(key, 0.0) - value / scale
        return combined

    factor = float(math.factorial(order))
    return sorted((key, factor * value) for key, value in weights(nodes).items())


class _AxisOperator:
    """``d^order/dxi^order`` or ``d^order/dzeta^order`` on the grid, with quadrature weights.

    ``order = 2`` is the classical thin-plate **curvature energy** ``integral kappa^2 ds``
    (what the gap names); ``order = 3`` is the **curvature-variation energy**
    ``integral (dkappa/ds)^2 ds``, which a parabola satisfies exactly and which therefore
    does not fight a legitimately highly curved region such as the ``zeta**(1/10)``
    section of a VLCC.  Both are offered because the difference is not academic: see
    ``docs/lines_fairing.md`` for the measurement that made ``order = 3`` the default.
    """

    def __init__(
        self,
        label: str,
        n: int,
        m: int,
        abscissa: Sequence[float],
        axis: str,
        *,
        order: int = 2,
        boundary_margin: int = 1,
    ) -> None:
        if axis not in ("station", "waterline"):
            raise FairingError(f"_AxisOperator: unknown axis {axis!r}")
        if order not in (2, 3):
            raise FairingError(f"_AxisOperator: order must be 2 or 3, got {order!r}")
        if boundary_margin < 0:
            raise FairingError(f"_AxisOperator: boundary_margin must be >= 0, got {boundary_margin}")
        self.label = label
        self.n = n
        self.m = m
        self.axis = axis
        self.order = order
        self.size = n * m
        self.rows: list[tuple[int, list[tuple[int, float]], float]] = []
        for k in range(boundary_margin, len(abscissa) - boundary_margin):
            coefficients = _divided_difference_coefficients(abscissa, k, order)
            if coefficients is None:
                continue
            # Trapezoid measure of the stencil: the interval spanned divided by the
            # number of sub-intervals, so ``sum q * (D**order)**2`` approximates the
            # integral of the squared derivative along the curve.
            weight = (abscissa[coefficients[-1][0]] - abscissa[coefficients[0][0]]) / order
            for c in range(m if axis == "station" else n):
                node = k * m + c if axis == "station" else c * m + k
                entries = [
                    (i * m + c if axis == "station" else c * m + i, value)
                    for i, value in coefficients
                ]
                self.rows.append((node, entries, weight))

    def evaluate(self, field: Sequence[float]) -> list[float]:
        out = [0.0] * self.size
        for node, coefficients, _ in self.rows:
            out[node] = math.fsum(value * field[index] for index, value in coefficients)
        return out

    def transpose_of_weighted(self, residual: Sequence[float]) -> list[float]:
        """``A^T Q r`` -- the gradient contribution of this direction."""
        acc = [0.0] * self.size
        for node, coefficients, weight in self.rows:
            scaled = weight * residual[node]
            if scaled == 0.0:
                continue
            for index, value in coefficients:
                acc[index] += value * scaled
        return acc

    def bend_energy(self, field: Sequence[float]) -> float:
        out = self.evaluate(field)
        return math.fsum(weight * out[node] ** 2 for node, _, weight in self.rows)


class _TangentOperator:
    """``d^order/dtau^order`` along the buttock tangent: the curvature of the level curves.

    Order 2 is built from the **tangential quadratic form** ``eta_tt = t^T H t`` with
    ``t = (-eta_zeta, eta_xi)/|grad eta|`` and the Hessian read off the grid by central
    divided differences -- exact for a quadratic surface, and therefore exact on the
    grid polynomials the axis operators also reproduce.  (An earlier version of this
    class sampled ``eta`` at ``p +- h t`` through a bilinear interpolant; bilinear
    interpolation is not quadratic-exact, and the second difference of its ``O(h**2)``
    error divided by ``h**2`` is ``O(1)``, which measured 100 % error on a quadratic
    surface.  That is why the form used here is the analytic one.)

    Order 3 -- the variation of the buttock curvature -- keeps the sampled route with
    the uniform-spacing four-point stencil, at the cost of that lower order accuracy;
    it is offered for the study, not as the default.

    Nodes whose stencil would leave the grid, or sit where ``grad eta`` vanishes (an
    extremum has no level-curve tangent), contribute no row -- counted in
    :attr:`nodes_without_row` rather than silently penalised.
    """

    #: Uniform-spacing stencil for ``d3/ds3``: divide by ``h**3`` and dot with the samples.
    _STENCIL_3: tuple[tuple[int, ...], tuple[float, ...]] = ((-1, 0, 1, 2), (-1.0, 3.0, -3.0, 1.0))

    def __init__(
        self,
        label: str,
        xi: Sequence[float],
        zeta: Sequence[float],
        eta: Sequence[Sequence[float]],
        *,
        order: int = 2,
        step_fraction: float = 0.5,
        boundary_margin: int = 1,
    ) -> None:
        if order not in (2, 3):
            raise FairingError(f"_TangentOperator: order must be 2 or 3, got {order!r}")
        self.label = label
        self.order = order
        self.n = len(xi)
        self.m = len(zeta)
        self.size = self.n * self.m
        self.rows: list[tuple[int, list[tuple[int, float]], float]] = []
        self.nodes_without_row = 0
        self.measure = 0.5 * (
            (xi[-1] - xi[0]) / max(self.n - 1, 1) + (zeta[-1] - zeta[0]) / max(self.m - 1, 1)
        )
        base_step = step_fraction * min(min(_spacings(xi)), min(_spacings(zeta)))
        for i in range(self.n):
            for j in range(self.m):
                margin = boundary_margin
                if (
                    i < margin
                    or i > self.n - 1 - margin
                    or j < margin
                    or j > self.m - 1 - margin
                ):
                    self.nodes_without_row += 1
                    continue
                d_xi = (eta[i + 1][j] - eta[i - 1][j]) / (xi[i + 1] - xi[i - 1])
                d_zeta = (eta[i][j + 1] - eta[i][j - 1]) / (zeta[j + 1] - zeta[j - 1])
                gradient = math.hypot(d_xi, d_zeta)
                if not (gradient > 1e-12):
                    self.nodes_without_row += 1
                    continue
                tau = (-d_zeta / gradient, d_xi / gradient)
                if order == 2:
                    row = self._quadratic_form_row(xi, zeta, i, j, tau)
                else:
                    row = self._sampled_row(xi, zeta, i, j, tau, base_step)
                if row is None:
                    self.nodes_without_row += 1
                else:
                    self.rows.append(row)

    def _quadratic_form_row(
        self, xi: Sequence[float], zeta: Sequence[float], i: int, j: int, tau: Sequence[float]
    ) -> tuple[int, list[tuple[int, float]], float] | None:
        along_xi = _divided_difference_coefficients(xi, i, 2)
        along_zeta = _divided_difference_coefficients(zeta, j, 2)
        if along_xi is None or along_zeta is None:
            return None
        node = i * self.m + j
        collected: dict[int, float] = {}
        tx, tz = tau
        for index, value in along_xi:
            collected[index * self.m + j] = collected.get(index * self.m + j, 0.0) + tx * tx * value
        for index, value in along_zeta:
            collected[i * self.m + index] = collected.get(i * self.m + index, 0.0) + tz * tz * value
        if abs(2.0 * tx * tz) > 0.0:
            span_xi = xi[i + 1] - xi[i - 1]
            span_zeta = zeta[j + 1] - zeta[j - 1]
            factor = 2.0 * tx * tz / (span_xi * span_zeta)
            corners = (
                ((i + 1) * self.m + j + 1, 1.0),
                ((i + 1) * self.m + j - 1, -1.0),
                ((i - 1) * self.m + j + 1, -1.0),
                ((i - 1) * self.m + j - 1, 1.0),
            )
            for index, sign in corners:
                collected[index] = collected.get(index, 0.0) + sign * factor
        return node, sorted(collected.items()), self.measure

    def _sampled_row(
        self,
        xi: Sequence[float],
        zeta: Sequence[float],
        i: int,
        j: int,
        tau: Sequence[float],
        base_step: float,
    ) -> tuple[int, list[tuple[int, float]], float] | None:
        node = i * self.m + j
        offsets, coefficients = self._STENCIL_3
        for factor in (1.0, 0.5, 0.25):
            h = base_step * factor
            collected: dict[int, float] = {}
            ok = True
            for offset, coefficient in zip(offsets, coefficients):
                sample = self._bilinear(
                    xi, zeta, xi[i] + offset * h * tau[0], zeta[j] + offset * h * tau[1]
                )
                if sample is None:
                    ok = False
                    break
                for index, weight in sample:
                    collected[index] = collected.get(index, 0.0) + coefficient * weight / h**3
            if not ok:
                continue
            return node, sorted(collected.items()), h
        return None

    def _bilinear(
        self, xi: Sequence[float], zeta: Sequence[float], x: float, z: float
    ) -> list[tuple[int, float]] | None:
        if not (xi[0] <= x <= xi[-1]) or not (zeta[0] <= z <= zeta[-1]):
            return None
        i = _locate(xi, x)
        j = _locate(zeta, z)
        if i < 0 or j < 0:
            return None
        sx = (x - xi[i]) / (xi[i + 1] - xi[i])
        sz = (z - zeta[j]) / (zeta[j + 1] - zeta[j])
        return [
            (i * self.m + j, (1.0 - sx) * (1.0 - sz)),
            (i * self.m + j + 1, (1.0 - sx) * sz),
            ((i + 1) * self.m + j, sx * (1.0 - sz)),
            ((i + 1) * self.m + j + 1, sx * sz),
        ]

    def evaluate(self, field: Sequence[float]) -> list[float]:
        out = [0.0] * self.size
        for node, coefficients, _ in self.rows:
            out[node] = math.fsum(c * field[index] for index, c in coefficients)
        return out

    def transpose_of_weighted(self, residual: Sequence[float]) -> list[float]:
        acc = [0.0] * self.size
        for node, coefficients, weight in self.rows:
            scaled = weight * residual[node]
            if scaled == 0.0:
                continue
            for index, coefficient in coefficients:
                acc[index] += coefficient * scaled
        return acc

    def bend_energy(self, field: Sequence[float]) -> float:
        out = self.evaluate(field)
        return math.fsum(weight * out[node] ** 2 for node, _, weight in self.rows)


def _locate(abscissa: Sequence[float], value: float) -> int:
    """Index ``i`` with ``abscissa[i] <= value <= abscissa[i+1]`` (``-1`` if none)."""
    if value < abscissa[0] or value > abscissa[-1]:
        return -1
    low, high = 0, len(abscissa) - 1
    while high - low > 1:
        middle = (low + high) // 2
        if abscissa[middle] <= value:
            low = middle
        else:
            high = middle
    return low

# ---------------------------------------------------------------------------
# the grid, the freeze mask and the problem
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class _Grid:
    """The table in non-dimensional coordinates, plus what the fairing may touch."""

    xi: tuple[float, ...]
    zeta: tuple[float, ...]
    eta: tuple[tuple[float, ...], ...]
    half_beam: float
    draft: float

    @property
    def n(self) -> int:
        return len(self.xi)

    @property
    def m(self) -> int:
        return len(self.zeta)

    def flat(self) -> list[float]:
        return [value for row in self.eta for value in row]

    def frozen(self) -> list[bool]:
        """Perpendiculars, keel and every already-zero ordinate stay where they are."""
        mask: list[bool] = []
        for i in range(self.n):
            for j in range(self.m):
                value = self.eta[i][j]
                mask.append(
                    i == 0 or i == self.n - 1 or j == 0 or j == self.m - 1 or value <= 0.0
                )
        return mask


def _grid_of(table: ShipOffsetTable) -> _Grid:
    half_beam = 0.5 * table.beam
    if not (half_beam > 0.0):
        raise FairingError(f"{table.name}: degenerate table, half-beam is {half_beam!r}")
    length = table.length
    draft = table.design_draft
    if not (length > 0.0 and draft > 0.0):
        raise FairingError(
            f"{table.name}: length/draft must be positive, got {length!r}/{draft!r}"
        )
    xi = tuple(x / length for x in table.stations)
    zeta = tuple(z / draft for z in table.waterlines)
    eta = tuple(tuple(y / half_beam for y in row) for row in table.half_breadths)
    return _Grid(xi=xi, zeta=zeta, eta=eta, half_beam=half_beam, draft=draft)


def _table_from(grid: _Grid, table: ShipOffsetTable, flat: Sequence[float]) -> ShipOffsetTable:
    m = grid.m
    rows = tuple(
        tuple(flat[i * m + j] * grid.half_beam for j in range(m)) for i in range(grid.n)
    )
    return ShipOffsetTable(
        name=table.name,
        length=table.length,
        stations=table.stations,
        waterlines=table.waterlines,
        half_breadths=rows,
        design_draft=table.design_draft,
        provenance=table.provenance,
        units=table.units,
        waterline_exponent=table.waterline_exponent,
    )


def _operators(
    grid: _Grid,
    eta: Sequence[Sequence[float]],
    *,
    order: int = 3,
    boundary_margin: int = 1,
) -> tuple[Any, Any, Any]:
    waterlines = _AxisOperator(
        "waterlines", grid.n, grid.m, grid.xi, axis="station",
        order=order, boundary_margin=boundary_margin,
    )  # d^order/dxi^order at fixed zeta
    stations = _AxisOperator(
        "stations", grid.n, grid.m, grid.zeta, axis="waterline",
        order=order, boundary_margin=boundary_margin,
    )  # d^order/dzeta^order at fixed xi
    buttocks = _TangentOperator(
        "buttocks", grid.xi, grid.zeta, eta, order=order, boundary_margin=boundary_margin
    )
    return waterlines, stations, buttocks


def _free_indices(frozen: Sequence[bool]) -> list[int]:
    return [index for index, flag in enumerate(frozen) if not flag]


def _energy_share_per_waterline(grid: _Grid, *, order: int = 2) -> dict[int, float]:
    """Share of the ``zeta``-curvature bending energy anchored on each waterline row.

    The rows of the ``stations`` operator are the curvature of the section profiles, so
    the ones anchored at waterline ``j`` are exactly the ones that see the near-keel band.
    Normalising by the total makes the number grid- and scale-free.
    """
    operator = _AxisOperator(
        "stations", grid.n, grid.m, grid.zeta, axis="waterline", order=order
    )
    residual = operator.evaluate(grid.flat())
    energy: dict[int, float] = {}
    for node, _coefficients, weight in operator.rows:
        index = node % grid.m
        energy[index] = energy.get(index, 0.0) + weight * residual[node] ** 2
    total = math.fsum(energy.values())
    if not total > 0.0:
        return {index: 0.0 for index in range(grid.m)}
    return {index: value / total for index, value in energy.items()}


def feature_waterlines(
    table: ShipOffsetTable, *, threshold: float = _FEATURE_SHARE, order: int = 2
) -> tuple[int, ...]:
    """Waterline rows whose section curvature is a **design feature**, not unfairness.

    A power-law section ``zeta**(1/m)`` (Series 60 ``m = 24``, the KVLCC2-class hull
    ``m = 10``) turns the curvature energy of the near-keel band into an unbounded number:
    the fairing then spends the whole max-norm deviation budget on that band and the rest
    of the table barely moves (the limitation measured in ``docs/lines_fairing.md``
    before this helper existed).  Declaring the band a *feature* -- the language a
    B-spline fairing uses for a crease -- removes its rows from the energy, so the budget
    is spent on the smooth part of the surface.

    The detection is one number per row: the share of the ``zeta``-curvature energy the
    row carries.  A row above ``threshold`` is reported.  The keel row (index ``0``) and
    the deck row (index ``m-1``) are never reported: they are frozen by
    :meth:`_Grid.frozen`, so they are not unknowns of the problem in the first place.
    """
    if not (0.0 < threshold <= 1.0):
        raise FairingError(
            f"feature_waterlines: threshold must be in (0, 1], got {threshold!r}"
        )
    grid = _grid_of(table)
    share = _energy_share_per_waterline(grid, order=order)
    return tuple(
        index for index in range(1, grid.m - 1) if share.get(index, 0.0) > threshold
    )


def _resolve_features(
    features: Sequence[int] | str, grid: _Grid, *, order: int
) -> tuple[int, ...]:
    """Turn the ``features`` request of :func:`fair_offset_table` into row indices."""
    if isinstance(features, str):
        if features in ("", "off", "none"):
            return ()
        if features != "auto":
            raise FairingError(
                "fair_offset_table: features must be a sequence of waterline indices, "
                f"'off' or 'auto', got {features!r}"
            )
        share = _energy_share_per_waterline(grid, order=order)
        return tuple(
            index for index in range(1, grid.m - 1) if share.get(index, 0.0) > _FEATURE_SHARE
        )
    resolved: list[int] = []
    for index in features:
        if not isinstance(index, int) or isinstance(index, bool):
            raise FairingError(
                f"fair_offset_table: feature waterline indices must be integers, got {index!r}"
            )
        if not 0 <= index < grid.m:
            raise FairingError(
                f"fair_offset_table: feature waterline index {index} outside the grid "
                f"(0 <= index < {grid.m})"
            )
        if index not in resolved:
            resolved.append(index)
    return tuple(sorted(resolved))


def _drop_feature_rows(
    operators: Sequence[Any], features: Sequence[int], m: int
) -> int:
    """Remove the ``zeta``-direction rows anchored on ``features``; return how many.

    Only the ``stations`` (``d^n/dzeta^n``) and ``buttocks`` (level-curve curvature)
    families are touched.  The ``waterlines`` family differentiates along ``x`` at fixed
    ``zeta``, where a power-law section is smooth: excluding it there would delete real
    evidence rather than a design feature.
    """
    if not features:
        return 0
    excluded = set(features)
    dropped = 0
    for operator in operators:
        if operator.label == "waterlines":
            continue
        kept = [row for row in operator.rows if (row[0] % m) not in excluded]
        dropped += len(operator.rows) - len(kept)
        operator.rows = kept
    return dropped


def _solve_with_penalty(
    eta0: Sequence[float],
    operators: Sequence[Any],
    lam: Sequence[float],
    free: Sequence[int],
    size: int,
    *,
    start: Sequence[float] | None = None,
    tolerance: float = 1e-9,
    max_iterations: int = 400,
) -> tuple[list[float], int, float]:
    """Solve ``min ||delta||^2 + sum_d lam_d ||D_d (eta0 + delta)||^2_Q`` on the free nodes.

    The data weight is ``1`` per node, so ``lam`` carries the units of the penalty: a
    larger ``lam`` buys fairness with deviation, which is exactly the trade the
    tolerance budget is about.
    """
    position = [-1] * size
    for local, index in enumerate(free):
        position[index] = local

    def scatter(reduced: Sequence[float]) -> list[float]:
        full = [0.0] * size
        for local, index in enumerate(free):
            full[index] = reduced[local]
        return full

    rhs = [0.0] * len(free)
    for weight, operator in zip(lam, operators):
        if weight == 0.0:
            continue
        gradient = operator.transpose_of_weighted(operator.evaluate(eta0))
        for local, index in enumerate(free):
            rhs[local] -= weight * gradient[index]

    def matvec(reduced: Sequence[float]) -> list[float]:
        full = scatter(reduced)
        out = list(reduced)
        for weight, operator in zip(lam, operators):
            if weight == 0.0:
                continue
            gradient = operator.transpose_of_weighted(operator.evaluate(full))
            for local, index in enumerate(free):
                out[local] += weight * gradient[index]
        return out

    reduced, iterations, residual = _conjugate_gradient(
        matvec, rhs, start, tolerance=tolerance, max_iterations=max_iterations
    )
    return scatter(reduced), iterations, residual


# ---------------------------------------------------------------------------
# the measurable evidence
# ---------------------------------------------------------------------------


def _flips_per_direction(grid: _Grid, flat: Sequence[float]) -> dict[str, int]:
    eta = [[flat[i * grid.m + j] for j in range(grid.m)] for i in range(grid.n)]
    waterlines = sum(
        curvature_sign_flips([eta[i][j] for i in range(grid.n)], grid.xi)
        for j in range(grid.m)
    )
    stations = sum(
        curvature_sign_flips([eta[i][j] for j in range(grid.m)], grid.zeta)
        for i in range(grid.n)
    )
    buttocks = sum(
        polyline_sign_flips(curve)
        for curve in buttock_curves(grid.eta, grid.xi, grid.zeta, levels=_buttock_levels(grid))
    )
    return {"waterlines": waterlines, "stations": stations, "buttocks": buttocks}


#: Fallback level fractions (of the peak ordinate) used when no candidate level of
#: :func:`buttock_levels` is straddled by enough stations.
_FALLBACK_LEVEL_FRACTIONS: tuple[float, ...] = (0.2, 0.35, 0.5, 0.65, 0.8)

#: Number of candidate levels the adaptive search scans, and their range as fractions of
#: the peak ordinate.  The bottom of the range is above ``0`` because the keel row is
#: frozen and a level at ``eta = 0`` degenerates onto it.
_LEVEL_CANDIDATES = 61
_LEVEL_RANGE = (0.05, 0.95)

#: Default coverage floor of a chosen level, as a fraction of the station count: a buttock
#: has to carry at least this many samples for its flip count to mean anything.
_MIN_COVERAGE = 0.3

#: Largest ``zeta`` step (as a fraction of the draft) at which two crossings of the same
#: level are still the same branch.  Above it the branch ends and a new one starts: a
#: level set that jumps is a gap, not a chord, and joining the two would invent curvature.
_TRACE_JUMP = 0.35


def _straddling_stations(eta: Sequence[Sequence[float]], level: float) -> int:
    """How many stations have ``level`` strictly inside their ``zeta`` profile."""
    return sum(1 for row in eta if min(row) < level < max(row))


def _longest_branch(
    eta: Sequence[Sequence[float]],
    xi: Sequence[float],
    zeta: Sequence[float],
    level: float,
) -> int:
    """Samples on the longest traced branch of the level set -- the coverage of ``level``.

    A level that only grazes a section (``level`` just under the local maximum breadth) is
    crossed at every station there, but at a ``zeta`` that jumps between neighbours: the
    traced branches are then 1, 5 and 1 points long, not one 7-point buttock.  Counting
    the traced branches instead of the straddled stations is what keeps such a level from
    passing the coverage floor on the strength of stations it does not actually join.
    """
    return max((len(branch) for branch in _trace_branches(eta, xi, zeta, level)), default=0)


def buttock_levels(
    eta: Sequence[Sequence[float]],
    xi: Sequence[float],
    zeta: Sequence[float],
    *,
    count: int = 5,
    min_samples: int | None = None,
) -> tuple[float, ...]:
    """Levels ``eta = const`` at which the buttock family is contoured, chosen from the table.

    A buttock exists only over the length where that half-breadth is reached, so a level
    near the peak spans a handful of stations and contributes nothing to the flip count
    once the ``>= 3`` point floor is applied.  Fixed fractions of the *global* peak
    therefore under-cover the third family on any table whose maximum breadth lives over
    a short stretch: on the built-in forms' 21 x 17 grids most of the five nominal levels
    came back with one or two points.  The levels are instead chosen from the data --
    scan :data:`_LEVEL_CANDIDATES` fractions of the peak in :data:`_LEVEL_RANGE`, keep
    those whose longest traced branch carries at least ``min_samples`` samples (default
    ``max(5, round(_MIN_COVERAGE * n_stations))``, i.e. a third of the stations), and spread
    ``count`` of the survivors evenly over the surviving range -- so the number of stations behind the
    metric is a property of the report, not of the caller's luck.

    ``xi`` and ``zeta`` are accepted (and length-checked) because the caller's grid is
    the definition of the table; the selection itself reads only ``eta``.
    """
    rows = [list(row) for row in eta]
    if not rows:
        return ()
    width = len(zeta)
    if any(len(row) != width for row in rows):
        raise FairingError("buttock_levels: eta is not a (station x waterline) grid")
    if len(xi) != len(rows):
        raise FairingError(
            f"buttock_levels: {len(rows)} stations against {len(xi)} abscissae"
        )
    if count < 1:
        raise FairingError(f"buttock_levels: count must be >= 1, got {count!r}")
    floor = max(5, round(_MIN_COVERAGE * len(rows))) if min_samples is None else min_samples
    if floor < 3:
        raise FairingError(
            f"buttock_levels: min_samples must be >= 3 (a flip count needs 3 points), "
            f"got {floor!r}"
        )
    peak = max((max(row) for row in rows if row), default=0.0)
    if not (peak > 0.0):
        return ()
    low, high = _LEVEL_RANGE
    candidates = [
        peak * (low + (high - low) * k / (_LEVEL_CANDIDATES - 1))
        for k in range(_LEVEL_CANDIDATES)
    ]
    eligible = [
        level for level in candidates
        if _longest_branch(rows, xi, zeta, level) >= floor
    ]
    if not eligible:
        eligible = [
            peak * fraction
            for fraction in _FALLBACK_LEVEL_FRACTIONS
            if _straddling_stations(rows, peak * fraction) >= 3
        ]
    if not eligible:
        return ()
    if len(eligible) <= count:
        return tuple(eligible)
    picked: list[float] = []
    for k in range(count):
        index = round(k * (len(eligible) - 1) / (count - 1))
        level = eligible[index]
        if not picked or level > picked[-1]:
            picked.append(level)
    return tuple(picked)


def _buttock_levels(grid: _Grid, count: int = 5) -> tuple[float, ...]:
    """Adaptive levels of :func:`buttock_levels` for a normalised grid."""
    return buttock_levels(grid.eta, grid.xi, grid.zeta, count=count)


def _level_crossings(
    profile: Sequence[float], zeta: Sequence[float], level: float
) -> dict[int, float]:
    """Every ``zeta`` at which the piecewise-linear profile crosses ``level``.

    Keyed by the *segment* ``[j, j + 1]`` the crossing falls in.  That key is what makes
    the branches of a level set separable: along one branch a level is crossed between
    the same pair of waterlines from station to station, so grouping the crossings of a
    level by their key recovers the branches without fitting or sorting anything.  The
    old contour bailed out whenever ``level`` was not bracketed by the endpoints of the
    profile, which silently dropped every level that a tumblehome section crosses twice.
    Exact node hits are assigned to the segment that starts at the node, the last node to
    the segment that ends there, so a crossing is produced once.
    """
    marks: dict[int, float] = {}
    last = len(profile) - 1
    for j in range(last):
        below, above = profile[j], profile[j + 1]
        if below == above:
            continue
        if below == level:
            key, value = (j if j > 0 else 0), zeta[j]
        elif above == level:
            key, value = (j + 1 if j + 1 < last else j), zeta[j + 1]
        elif (below - level) * (above - level) < 0.0:
            key = j
            value = zeta[j] + (level - below) / (above - below) * (zeta[j + 1] - zeta[j])
        else:
            continue
        marks.setdefault(key, value)
    return marks


def _trace_branches(
    eta: Sequence[Sequence[float]],
    xi: Sequence[float],
    zeta: Sequence[float],
    level: float,
    *,
    samples: int = 0,
) -> list[list[tuple[float, float]]]:
    """Branches of the level set ``eta = level``, traced station by station.

    At every station the level is crossed once per branch; a crossing is attached to the
    open branch whose last sample is nearest in ``zeta``, and a branch that no crossing
    matches simply ends (and a later crossing starts a new branch).  Continuity tracking
    -- rather than grouping the crossings by the waterline segment they fall in -- is what
    keeps a buttock that runs across several segments one curve: on the Wigley table the
    topmost level crosses segment 12 amidships and segment 14 near the ends, and grouping
    by segment cuts one 5-point buttock into a branch of 3 and an orphan of 2.  A crossing
    further than :data:`_TRACE_JUMP` (fraction of the ``zeta`` range) from every open
    branch is treated as a new branch, so a jump in ``zeta`` is never a chord.
    """
    span = (zeta[-1] - zeta[0]) if len(zeta) > 1 else 0.0
    max_jump = _TRACE_JUMP * span
    open_branches: list[list[tuple[float, float]]] = []
    done: list[list[tuple[float, float]]] = []
    for i, x in enumerate(xi):
        crossings = sorted(_level_crossings(eta[i], zeta, level).values())
        taken = [False] * len(open_branches)
        for value in crossings:
            best: int | None = None
            best_distance = math.inf
            for index, branch in enumerate(open_branches):
                if taken[index]:
                    continue
                distance = abs(branch[-1][1] - value)
                if distance <= max_jump and distance < best_distance:
                    best, best_distance = index, distance
            if best is None:
                open_branches.append([(x, value)])
                taken.append(True)
            else:
                open_branches[best].append((x, value))
                taken[best] = True
        keep = [branch for index, branch in enumerate(open_branches) if taken[index]]
        done.extend(branch for index, branch in enumerate(open_branches) if not taken[index])
        open_branches = keep
    done.extend(open_branches)
    if samples > 0:
        done = [_resample(branch, samples) if len(branch) >= 3 else branch for branch in done]
    return done


def buttock_curves(
    eta: Sequence[Sequence[float]],
    xi: Sequence[float],
    zeta: Sequence[float],
    *,
    levels: Iterable[float],
    samples: int = 0,
    min_points: int = 3,
) -> list[list[tuple[float, float]]]:
    """Contour the ``eta`` grid at each level: the buttock ``zeta(xi)`` polylines.

    Inversion is by linear interpolation on the piecewise-linear profile at each station,
    so a buttock appears only over the length where that half-breadth exists -- which is
    the shape of a real buttock plan.  A level set may have **more than one branch**: a
    section that is widest below the deck (tumblehome, or any hull whose half-breadth
    falls again above the maximum breadth) crosses a level twice, and the two crossings
    are different buttocks.  Every branch is returned as its own polyline (see
    :func:`_trace_branches`), and a polyline needs ``min_points`` samples to be returned
    at all.  ``samples > 0`` resamples each branch on a uniform ``xi`` grid (the default
    keeps the tabulated stations).
    """
    curves: list[list[tuple[float, float]]] = []
    for level in levels:
        for branch in _trace_branches(eta, xi, zeta, level, samples=samples):
            if len(branch) >= min_points:
                curves.append(branch)
    return curves


def buttock_coverage(
    eta: Sequence[Sequence[float]],
    xi: Sequence[float],
    zeta: Sequence[float],
    *,
    levels: Sequence[float] | None = None,
    min_points: int = 3,
) -> dict[str, Any]:
    """How much of the third family the metric actually sees, per level.

    A flip count is only as meaningful as the curves behind it, so the report states, for
    every level: how many stations straddle it, how many branches it has, and how many
    samples each returned branch carries.  A level whose branches are all shorter than
    ``min_points`` contributes nothing and shows up here as ``branches: 0`` instead of
    passing as a silent zero -- which is exactly how the third family lost most of its
    weight in the earlier report of this gap.
    """
    chosen = tuple(levels) if levels is not None else buttock_levels(eta, xi, zeta)
    peak = max((max(row) for row in eta if row), default=0.0)
    detail: list[dict[str, Any]] = []
    curves_total = 0
    for level in chosen:
        kept = [
            branch
            for branch in _trace_branches(eta, xi, zeta, level)
            if len(branch) >= min_points
        ]
        curves_total += len(kept)
        detail.append(
            {
                "level": level,
                "fraction_of_peak": (level / peak) if peak > 0.0 else None,
                "stations_straddling": _straddling_stations(eta, level),
                "branches": len(kept),
                "samples": [len(branch) for branch in kept],
            }
        )
    return {
        "levels": list(chosen),
        "curves": curves_total,
        "stations_total": len(xi),
        "detail": detail,
    }

def _resample(points: Sequence[tuple[float, float]], count: int) -> list[tuple[float, float]]:
    first, last = points[0][0], points[-1][0]
    out: list[tuple[float, float]] = []
    for k in range(count):
        target = first + (last - first) * k / (count - 1)
        for index in range(len(points) - 1):
            x0, z0 = points[index]
            x1, z1 = points[index + 1]
            if x0 <= target <= x1 and x1 > x0:
                out.append((target, z0 + (z1 - z0) * (target - x0) / (x1 - x0)))
                break
    return out


def _direction_metrics(
    grid_before: _Grid,
    grid_after: _Grid,
    operators: Sequence[Any],
    flat_before: Sequence[float],
    flat_after: Sequence[float],
) -> tuple[DirectionMetrics, ...]:
    flips_before = _flips_per_direction(grid_before, flat_before)
    flips_after = _flips_per_direction(grid_after, flat_after)
    metrics: list[DirectionMetrics] = []
    for key, operator in zip(DIRECTIONS, operators):
        if key == "buttocks":
            curves = len(buttock_curves(grid_before.eta, grid_before.xi, grid_before.zeta,
                                       levels=_buttock_levels(grid_before)))
        elif key == "waterlines":
            curves = grid_before.m
        else:
            curves = grid_before.n
        metrics.append(
            DirectionMetrics(
                direction=key,
                curves=curves,
                flips_before=flips_before[key],
                flips_after=flips_after[key],
                energy_before=operator.bend_energy(flat_before),
                energy_after=operator.bend_energy(flat_after),
            )
        )
    return tuple(metrics)


@dataclass(frozen=True)
class DirectionMetrics:
    """One family's evidence: curvature sign flips and bending energy, before/after."""

    direction: str
    curves: int
    flips_before: int
    flips_after: int
    energy_before: float
    energy_after: float

    @property
    def flips_reduced(self) -> bool:
        return self.flips_after < self.flips_before

    def to_dict(self) -> dict[str, Any]:
        return {
            "direction": self.direction,
            "curves": self.curves,
            "flips_before": self.flips_before,
            "flips_after": self.flips_after,
            "flips_delta": self.flips_after - self.flips_before,
            "bending_energy_before": self.energy_before,
            "bending_energy_after": self.energy_after,
            "energy_ratio": (self.energy_after / self.energy_before)
            if self.energy_before > 0.0
            else None,
        }


@dataclass(frozen=True)
class FairingIteration:
    """What one outer iteration achieved, and whether it was kept."""

    iteration: int
    penalty_waterlines: float
    penalty_stations: float
    penalty_buttocks: float
    deviation_max: float
    contraction: float
    cg_iterations: int
    total_flips_before: int
    total_flips_after: int
    metrics: tuple[DirectionMetrics, ...]
    accepted: bool
    feature_waterlines: tuple[int, ...] = ()
    penalty_rows_dropped: int = 0

    def to_dict(self) -> dict[str, Any]:
        return {
            "iteration": self.iteration,
            "penalty": {
                "waterlines": self.penalty_waterlines,
                "stations": self.penalty_stations,
                "buttocks": self.penalty_buttocks,
            },
            "deviation_max_m": self.deviation_max,
            "contraction": self.contraction,
            "cg_iterations": self.cg_iterations,
            "total_flips_before": self.total_flips_before,
            "total_flips_after": self.total_flips_after,
            "accepted": self.accepted,
            "feature_waterlines": list(self.feature_waterlines),
            "penalty_rows_dropped": self.penalty_rows_dropped,
            "directions": [metric.to_dict() for metric in self.metrics],
        }


@dataclass(frozen=True)
class ReferenceComparison:
    """Distance to a smoother table taken as the manual-fairing stand-in."""

    name: str
    max_before: float
    max_after: float
    rms_before: float
    rms_after: float

    @property
    def max_improvement(self) -> float | None:
        return (self.max_before / self.max_after) if self.max_after > 0.0 else None

    @property
    def rms_improvement(self) -> float | None:
        return (self.rms_before / self.rms_after) if self.rms_after > 0.0 else None

    def to_dict(self) -> dict[str, Any]:
        return {
            "reference": self.name,
            "max_before_m": self.max_before,
            "max_after_m": self.max_after,
            "rms_before_m": self.rms_before,
            "rms_after_m": self.rms_after,
            "max_improvement": (self.max_before / self.max_after)
            if self.max_after > 0.0
            else None,
            "rms_improvement": (self.rms_before / self.rms_after)
            if self.rms_after > 0.0
            else None,
        }


@dataclass(frozen=True)
class FairingResult:
    """The faired table plus the whole audit trail."""

    table_before: ShipOffsetTable
    table_after: ShipOffsetTable
    tolerance: float
    iterations: tuple[FairingIteration, ...]
    rejected_steps: int
    reference: ReferenceComparison | None = None
    feature_waterlines: tuple[int, ...] = ()
    delivered_flips: dict[str, int] | None = None

    # -- summary numbers -----------------------------------------------------
    @property
    def deviation_max(self) -> float:
        """Largest ordinate change, in metres (the number the tolerance bounds)."""
        return max(
            (abs(after - before) for row_before, row_after in zip(
                self.table_before.half_breadths, self.table_after.half_breadths
            ) for before, after in zip(row_before, row_after)),
            default=0.0,
        )

    @property
    def deviation_argmax(self) -> tuple[int, int]:
        """Grid cell that absorbed the largest ordinate change.

        Reported because it is the diagnosis of a weak run: when this lands in the
        near-keel band of a power-law section, the tolerance budget is being spent on a
        *design feature* and the rest of the table barely moves -- see
        ``docs/lines_fairing.md`` and :func:`fairing_study`.
        """
        best = (-1.0, 0, 0)
        for i, (row_before, row_after) in enumerate(
            zip(self.table_before.half_breadths, self.table_after.half_breadths)
        ):
            for j, (before, after) in enumerate(zip(row_before, row_after)):
                if abs(after - before) > best[0]:
                    best = (abs(after - before), i, j)
        return best[1], best[2]

    @property
    def deviation_fraction_of_half_beam(self) -> float:
        half_beam = 0.5 * self.table_before.beam
        return self.deviation_max / half_beam if half_beam else float("nan")

    @property
    def within_tolerance(self) -> bool:
        return self.deviation_max <= self.tolerance * (1.0 + 1e-12)

    @property
    def flips_before_by_direction(self) -> dict[str, int]:
        """Per-family flip count of the table as it was handed in."""
        return {metric.direction: metric.flips_before for metric in self.iterations[0].metrics}

    @property
    def flips_after_by_direction(self) -> dict[str, int]:
        """Per-family flip count of :attr:`table_after` -- measured on the delivered table.

        Deliberately **not** ``self.iterations[-1].metrics``: the last outer trial may have
        been rejected (its flip count rose), and then that trial describes a table which was
        thrown away.  Every reported number has to describe the ordinates the caller gets
        back, otherwise ``flips_after`` can disagree with ``table_after.digest()`` -- which
        is exactly what a rejected step used to do (published Series 60 / KVLCC2 rows
        reported flip counts of discarded tables).
        """
        if self.delivered_flips is not None:
            return dict(self.delivered_flips)
        return _flips_per_direction(*_grid_and_flat(self.table_after))

    @property
    def flips_before(self) -> int:
        return sum(self.flips_before_by_direction.values())

    @property
    def flips_after(self) -> int:
        return sum(self.flips_after_by_direction.values())

    @property
    def last_trial_accepted(self) -> bool:
        """Whether the final outer trial survived the accept rule."""
        return bool(self.iterations) and self.iterations[-1].accepted

    @property
    def monotone(self) -> bool:
        """The delivered total flip count never rose along the accepted path."""
        sequence = [self.iterations[0].total_flips_before]
        for iteration in self.iterations:
            if iteration.accepted:
                sequence.append(iteration.total_flips_after)
        sequence.append(self.flips_after)
        return all(later <= earlier for earlier, later in zip(sequence, sequence[1:]))

    @property
    def flips_reduced(self) -> bool:
        return self.flips_after < self.flips_before

    def accepted(self) -> bool:
        """The acceptance of ``GAP-CAD-C3-LINES-FAIRING``, as a boolean."""
        return self.within_tolerance and self.monotone and self.flips_reduced

    def hydrostatics(self) -> dict[str, Any]:
        before = table_integrals(self.table_before)
        after = table_integrals(self.table_after)
        payload: dict[str, Any] = {"before": before.to_dict(), "after": after.to_dict()}
        payload["relative_change"] = {
            key: ((getattr(after, key) - getattr(before, key)) / getattr(before, key))
            if getattr(before, key)
            else None
            for key in ("volume", "waterplane_area", "lcb", "vcb")
        }
        return payload

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA_REPORT,
            "table": self.table_before.name,
            "grid": {
                "n_stations": self.table_before.n_stations,
                "n_waterlines": self.table_before.n_waterlines,
                "length_m": self.table_before.length,
                "half_beam_m": 0.5 * self.table_before.beam,
                "draft_m": self.table_before.design_draft,
                "digest_before": self.table_before.digest(),
                "digest_after": self.table_after.digest(),
            },
            "tolerance_m": self.tolerance,
            "deviation_max_m": self.deviation_max,
            "deviation_fraction_of_half_beam": self.deviation_fraction_of_half_beam,
            "deviation_argmax": {
                "station_index": self.deviation_argmax[0],
                "waterline_index": self.deviation_argmax[1],
                "station_m": self.table_before.stations[self.deviation_argmax[0]],
                "waterline_m": self.table_before.waterlines[self.deviation_argmax[1]],
            },
            "within_tolerance": self.within_tolerance,
            "feature_waterlines": list(self.feature_waterlines),
            "rejected_steps": self.rejected_steps,
            "flips_before": self.flips_before,
            "flips_after": self.flips_after,
            "monotone": self.monotone,
            "flips_reduced": self.flips_reduced,
            "flips_before_by_direction": self.flips_before_by_direction,
            "flips_after_by_direction": self.flips_after_by_direction,
            "last_trial_accepted": self.last_trial_accepted,
            "accepted": self.accepted(),
            "iterations": [iteration.to_dict() for iteration in self.iterations],
            "reference": self.reference.to_dict() if self.reference else None,
            "hydrostatics": self.hydrostatics(),
        }


def _deviation(flat: Sequence[float], reference: Sequence[float], free: Sequence[int]) -> float:
    return max((abs(flat[index] - reference[index]) for index in free), default=0.0)


def _search_penalty(
    origin: Sequence[float],
    operators: Sequence[Any],
    weights: Sequence[float],
    free: Sequence[int],
    size: int,
    budget: float,
    penalty_range: tuple[float, float],
    bisection_steps: int,
    cg_tolerance: float,
    cg_max_iterations: int,
    warm: Sequence[float] | None,
) -> tuple[list[float], float, int, float]:
    """Largest common penalty whose correction fits ``budget``; returns the fit.

    The deviation grows monotonically with the penalty (the data term pulls harder as
    ``lam`` falls), so the search is a bracketing sweep -- decades at first, then a
    geometric bisection -- and every trial is warm-started from the previous one.
    """
    cache: dict[float, tuple[list[float], float, int]] = {}

    def trial(log_value: float) -> tuple[list[float], float, int]:
        key = round(log_value, 12)
        if key in cache:
            return cache[key]
        penalty = 10.0 ** log_value
        lam = tuple(weight * penalty for weight in weights)
        correction, cg_iterations, _ = _solve_with_penalty(
            origin, operators, lam, free, size,
            start=warm, tolerance=cg_tolerance, max_iterations=cg_max_iterations,
        )
        candidate = [origin[index] + correction[index] for index in range(size)]
        result = (candidate, _deviation(candidate, origin, free), cg_iterations)
        cache[key] = result
        return result

    low_limit, high_limit = (math.log10(value) for value in penalty_range)
    log_value = 0.0
    candidate, deviation, cg_iterations = trial(log_value)
    if deviation <= budget:
        log_low, log_high = log_value, None
        step = 2.0
        while log_high is None and log_value < high_limit:
            advanced = min(high_limit, log_value + step)
            if advanced == log_value:
                break
            log_value = advanced
            step *= 1.5
            candidate_high, deviation_high, iterations_high = trial(log_value)
            if deviation_high > budget:
                log_high = log_value
            else:
                log_low, candidate, deviation, cg_iterations = (
                    log_value, candidate_high, deviation_high, iterations_high,
                )
        if log_high is None:
            log_high = min(high_limit, log_value + step)
    else:
        log_high = log_value
        log_low = None
        step = -2.0
        while log_low is None and log_value > low_limit:
            advanced = max(low_limit, log_value + step)
            if advanced == log_value:
                break
            log_value = advanced
            step *= 1.5
            candidate_low, deviation_low, iterations_low = trial(log_value)
            if deviation_low <= budget:
                log_low = log_value
                candidate, deviation, cg_iterations = candidate_low, deviation_low, iterations_low
            else:
                log_high = log_value
        if log_low is None:
            raise FairingError(
                "fair_offset_table: no penalty in "
                f"[{low_limit:g}, {high_limit:g}] keeps the deviation inside the tolerance; "
                "the requested tolerance is below the round-off of the algebra"
            )

    for _ in range(bisection_steps):
        middle = 0.5 * (log_low + log_high)
        candidate_mid, deviation_mid, iterations_mid = trial(middle)
        if deviation_mid <= budget:
            log_low = middle
            candidate, deviation, cg_iterations = candidate_mid, deviation_mid, iterations_mid
        else:
            log_high = middle
    return candidate, deviation, cg_iterations, log_low


def fair_offset_table(
    table: ShipOffsetTable,
    *,
    tolerance: float,
    iterations: int = 3,
    order: int = 2,
    weights: tuple[float, float, float] = (1.0, 1.0, 1.0),
    bisection_steps: int = 8,
    cg_tolerance: float = 1e-8,
    cg_max_iterations: int = 240,
    penalty_range: tuple[float, float] = (1e-16, 1e2),
    reference: ShipOffsetTable | None = None,
    reference_name: str = "reference",
    features: Sequence[int] | str = (),
) -> FairingResult:
    """Fair ``table`` on all three families inside ``tolerance`` (metres).

    ``tolerance`` is the largest ordinate change allowed, in metres; the penalty is
    bisected in log space so the correction uses :data:`_BUDGET_FRACTION` of it.
    ``order`` selects the energy: ``2`` is the curvature energy ``integral kappa^2`` the
    gap names, ``3`` the curvature-variation energy ``integral (dkappa/ds)^2``.  An
    outer iteration is kept only when the total curvature sign-flip count does not rise
    (see the module docstring); ``weights`` fixes the relative cost of the three
    families and ``iterations`` bounds the outer sweep, because the buttock operator
    depends on the table it is evaluated on.

    ``features`` declares waterline rows whose section curvature is a *design feature*
    (the near-keel band of a ``zeta**(1/m)`` section): their ``zeta``-direction rows are
    removed from the energy, so the deviation budget is spent on the smooth part of the
    table instead of on the singular band.  ``"auto"`` detects the rows with
    :func:`feature_waterlines`; a sequence fixes them by hand; the default ``()`` leaves
    the energy exactly as it was before the option existed (the same table comes out).
    """
    if not (tolerance > 0.0):
        raise FairingError(f"fair_offset_table: tolerance must be positive, got {tolerance!r}")
    if iterations < 1:
        raise FairingError(f"fair_offset_table: iterations must be >= 1, got {iterations!r}")
    if order not in (2, 3):
        raise FairingError(f"fair_offset_table: order must be 2 or 3, got {order!r}")
    if len(weights) != 3 or any(w < 0.0 for w in weights) or not any(w > 0.0 for w in weights):
        raise FairingError(f"fair_offset_table: weights must be 3 non-negative numbers with a "
                           f"positive one, got {weights!r}")
    if isinstance(features, str) and features not in ("", "off", "none", "auto"):
        raise FairingError(
            "fair_offset_table: features must be a sequence of waterline indices, "
            f"'off' or 'auto', got {features!r}"
        )

    grid = _grid_of(table)
    origin = grid.flat()
    size = len(origin)
    frozen = grid.frozen()
    free = _free_indices(frozen)
    if not free:
        raise FairingError(f"{table.name}: every ordinate is frozen, nothing to fair")
    budget = _BUDGET_FRACTION * tolerance / grid.half_beam

    current = list(origin)
    records: list[FairingIteration] = []
    rejected = 0
    warm: list[float] | None = None
    for outer in range(1, iterations + 1):
        current_grid = _Grid(
            xi=grid.xi, zeta=grid.zeta,
            eta=tuple(tuple(current[i * grid.m + j] for j in range(grid.m))
                      for i in range(grid.n)),
            half_beam=grid.half_beam, draft=grid.draft,
        )
        operators = _operators(current_grid, current_grid.eta, order=order)
        active_features = _resolve_features(features, current_grid, order=order)
        rows_dropped = _drop_feature_rows(operators, active_features, grid.m)
        flips_now = _flips_per_direction(current_grid, current)
        total_before = sum(flips_now.values())

        # Bracket the penalty, then bisect: the largest lam whose deviation still fits
        # the budget.  Bracketing (rather than a fixed decade range) is what makes the
        # search independent of ``order`` and of the grid: a third-difference energy of
        # a jittered table is ~h^-2 times larger than a second-difference one, so a fixed
        # range would silently pick the wrong end for one of the two.
        candidate, deviation, cg_iterations, log_penalty = _search_penalty(
            origin, operators, weights, free, size, budget,
            penalty_range, bisection_steps, cg_tolerance, cg_max_iterations, warm,
        )
        contraction = 1.0
        if deviation > budget:
            contraction = budget / deviation
            candidate = [
                origin[index] + contraction * (candidate[index] - origin[index])
                for index in range(size)
            ]
            deviation = _deviation(candidate, origin, free)

        candidate_grid = _Grid(
            xi=grid.xi, zeta=grid.zeta,
            eta=tuple(tuple(candidate[i * grid.m + j] for j in range(grid.m))
                      for i in range(grid.n)),
            half_beam=grid.half_beam, draft=grid.draft,
        )
        metrics = _direction_metrics(current_grid, candidate_grid, operators, current, candidate)
        total_after = sum(metric.flips_after for metric in metrics)
        accepted = total_after <= total_before
        records.append(
            FairingIteration(
                iteration=outer,
                penalty_waterlines=weights[0] * 10.0 ** log_penalty,
                penalty_stations=weights[1] * 10.0 ** log_penalty,
                penalty_buttocks=weights[2] * 10.0 ** log_penalty,
                deviation_max=deviation * grid.half_beam,
                contraction=contraction,
                cg_iterations=cg_iterations,
                total_flips_before=total_before,
                total_flips_after=total_after,
                metrics=metrics,
                accepted=accepted,
                feature_waterlines=tuple(active_features),
                penalty_rows_dropped=rows_dropped,
            )
        )
        if accepted:
            current = candidate
            if total_after == total_before:
                # The flip count has stopped moving: further sweeps would reproduce this
                # step, so the run stops here and the record shows the plateau.
                break
        else:
            rejected += 1
            break

    after = _table_from(grid, table, current)
    # The verdict has to describe the table that is returned, not the last trial: a
    # rejected trial (flip count rose) leaves ``table_after`` at the previous iterate.
    delivered_flips = _flips_per_direction(*_grid_and_flat(after))
    comparison = None
    if reference is not None:
        comparison = compare_to_reference(table, after, reference, name=reference_name)
    used = next(
        (record.feature_waterlines for record in reversed(records) if record.accepted),
        records[0].feature_waterlines if records else (),
    )
    return FairingResult(
        table_before=table,
        table_after=after,
        tolerance=tolerance,
        iterations=tuple(records),
        rejected_steps=rejected,
        reference=comparison,
        feature_waterlines=tuple(used),
        delivered_flips=delivered_flips,
    )


def _rms(values: Iterable[float]) -> float:
    collected = list(values)
    return math.sqrt(math.fsum(v * v for v in collected) / len(collected)) if collected else 0.0


def compare_to_reference(
    table: ShipOffsetTable,
    faired: ShipOffsetTable,
    reference: ShipOffsetTable,
    *,
    name: str = "reference",
) -> ReferenceComparison:
    """``faired`` and ``table`` against a smoother ``reference``, ordinate by ordinate.

    This is the "manual fairing" comparison of the acceptance: the reference is the
    table a designer would accept, so what matters is whether fairing moved the table
    *towards* it.  All three tables must share a grid; the reference is not required to
    satisfy ``tolerance``.
    """
    for other, label in ((faired, "faired"), (reference, "reference")):
        if (len(other.stations) != len(table.stations)
                or len(other.waterlines) != len(table.waterlines)):
            raise FairingError(
                f"compare_to_reference: {label} grid "
                f"{len(other.stations)}x{len(other.waterlines)} against "
                f"{len(table.stations)}x{len(table.waterlines)}"
            )
    before = [
        abs(b - r)
        for row_b, row_r in zip(table.half_breadths, reference.half_breadths)
        for b, r in zip(row_b, row_r)
    ]
    after = [
        abs(a - r)
        for row_a, row_r in zip(faired.half_breadths, reference.half_breadths)
        for a, r in zip(row_a, row_r)
    ]
    return ReferenceComparison(
        name=name,
        max_before=max(before, default=0.0),
        max_after=max(after, default=0.0),
        rms_before=_rms(before),
        rms_after=_rms(after),
    )


def perturb_table(
    table: ShipOffsetTable,
    *,
    amplitude: float,
    seed: int = 0,
    only_free: bool = True,
) -> ShipOffsetTable:
    """A deliberately *unfair* copy of ``table``: ordinates jittered by ``+-amplitude`` m.

    A measured or hand-drafted table is not smooth at the grid scale, and that is what
    the fairing is supposed to remove; this helper manufactures the defect so the
    removal can be measured.  ``random.Random`` is seeded, so the study is reproducible
    bit for bit.  With ``only_free`` the perpendiculars, the keel and the zero ordinates
    are left alone -- those are principal dimensions, not noise.
    """
    if amplitude < 0.0:
        raise FairingError(f"perturb_table: amplitude must be non-negative, got {amplitude!r}")
    grid = _grid_of(table)
    frozen = grid.frozen()
    generator = random.Random(seed)
    m = grid.m
    rows: list[tuple[float, ...]] = []
    for i in range(grid.n):
        row: list[float] = []
        for j in range(m):
            index = i * m + j
            base = table.half_breadths[i][j]
            if only_free and frozen[index]:
                row.append(base)
            else:
                row.append(max(0.0, base + generator.uniform(-amplitude, amplitude)))
        rows.append(tuple(row))
    return ShipOffsetTable(
        name=f"{table.name} (+-{amplitude:g} m ordinate jitter, seed {seed})",
        length=table.length,
        stations=table.stations,
        waterlines=table.waterlines,
        half_breadths=tuple(rows),
        design_draft=table.design_draft,
        provenance=table.provenance,
        units=table.units,
        waterline_exponent=table.waterline_exponent,
    )


# ---------------------------------------------------------------------------
# the study: three built-in tables, each deliberately made unfair first
# ---------------------------------------------------------------------------

DEFAULT_FORMS: tuple[str, ...] = ("wigley_iii", "series60_cb060", "kvlcc2_class")


def fairing_study(
    *,
    forms: Sequence[str] = DEFAULT_FORMS,
    n_stations: int = 21,
    n_waterlines: int = 17,
    noise_fraction: float = 0.004,
    tolerance_fractions: Sequence[float] = (0.004, 0.008, 0.016, 0.032),
    iterations: int = 3,
    order: int = 2,
    seed: int = 20261114,
    waterline_exponent: float = 1.0,
    feature_modes: Sequence[str] = ("off", "auto"),
) -> dict[str, Any]:
    """Fair every built-in form over a sweep of tolerances and report the evidence.

    ``noise_fraction`` and every entry of ``tolerance_fractions`` are fractions of the
    half-beam, so the experiment is scale-free: the same settings mean the same thing on
    a 4 m Wigley and on a 320 m tanker.

    The sweep is the point.  A single tolerance cannot distinguish "the fairing does not
    work" from "the tolerance is too tight", and on a full ``zeta**(1/m)`` section the
    two look identical -- the near-keel rows cost so much curvature energy that a
    max-norm deviation budget is spent there and the rest of the table barely moves.
    Sweeping the budget shows the flip count following it down, and
    :attr:`FairingResult.deviation_argmax` names the cell that is doing the absorbing.

    The analytic ordinates of each form play the manual-fairing reference.  They are
    smooth to machine precision by construction; no human-faired table exists in this
    repository.  That substitution is the study's main assumption and is restated in
    ``docs/lines_fairing.md``.

    ``feature_modes`` runs the whole sweep once per mode and keeps both: ``"off"`` is the
    plain energy and ``"auto"`` declares the detected design-feature rows
    (:func:`feature_waterlines`), so the effect of the declaration is measured on the same
    jittered tables rather than asserted.  The sweeps live in ``case["sweeps"][mode]``
    (the payload written before the option existed is ``case["sweeps"]["off"]``, and
    ``case["accepted_at"]`` still reports that mode).
    """
    from .offset_forms import offset_table

    unknown = [mode for mode in feature_modes if mode not in ("off", "auto")]
    if unknown:
        raise FairingError(
            f"fairing_study: unknown feature mode(s) {unknown!r}; use 'off' and/or 'auto'"
        )
    if not feature_modes:
        raise FairingError("fairing_study: feature_modes must name at least one mode")

    cases: list[dict[str, Any]] = []
    for key in forms:
        clean = offset_table(
            key, n_stations=n_stations, n_waterlines=n_waterlines,
            waterline_exponent=waterline_exponent,
        )
        half_beam = 0.5 * clean.beam
        unfair = perturb_table(clean, amplitude=noise_fraction * half_beam, seed=seed)
        clean_flips = _flips_per_direction(*_grid_and_flat(clean))
        unfair_flips = _flips_per_direction(*_grid_and_flat(unfair))
        # The third family is contoured, not tabulated, so its weight is only as real as
        # the number of samples behind its curves.  Record that, per level.
        grid_unfair = _grid_of(unfair)
        coverage = buttock_coverage(grid_unfair.eta, grid_unfair.xi, grid_unfair.zeta)
        sweeps: dict[str, list[dict[str, Any]]] = {}
        for mode in feature_modes:
            request: Sequence[int] | str = "auto" if mode == "auto" else ()
            sweep: list[dict[str, Any]] = []
            for fraction in tolerance_fractions:
                result = fair_offset_table(
                    unfair,
                    tolerance=fraction * half_beam,
                    iterations=iterations,
                    order=order,
                    reference=clean,
                    reference_name=f"{key} analytic ordinates (manual-fairing stand-in)",
                    features=request,
                )
                payload = result.to_dict()
                payload["tolerance_fraction_of_half_beam"] = fraction
                payload["feature_mode"] = mode
                payload["flips_by_direction"] = {
                    direction: {
                        "before": result.flips_before_by_direction[direction],
                        "after": result.flips_after_by_direction.get(direction),
                    }
                    for direction in result.flips_before_by_direction
                }
                sweep.append(payload)
            sweeps[mode] = sweep
        cases.append(
            {
                "form": key,
                "name": clean.name,
                "grid": {
                    "n_stations": clean.n_stations,
                    "n_waterlines": clean.n_waterlines,
                    "length_m": clean.length,
                    "half_beam_m": half_beam,
                    "draft_m": clean.design_draft,
                    "waterline_exponent": waterline_exponent,
                },
                "settings": {
                    "noise_amplitude_m": noise_fraction * half_beam,
                    "noise_fraction_of_half_beam": noise_fraction,
                    "iterations": iterations,
                    "order": order,
                    "seed": seed,
                    "feature_modes": list(feature_modes),
                },
                "analytic_flips": clean_flips,
                "unfair_flips": unfair_flips,
                "buttock_coverage": coverage,
                "sweeps": sweeps,
                "accepted_at": [
                    entry["tolerance_fraction_of_half_beam"]
                    for entry in sweeps.get("off", next(iter(sweeps.values())))
                    if entry["accepted"]
                ],
                "accepted_at_by_mode": {
                    mode: [
                        entry["tolerance_fraction_of_half_beam"]
                        for entry in sweep
                        if entry["accepted"]
                    ]
                    for mode, sweep in sweeps.items()
                },
            }
        )
    return {
        "schema": SCHEMA_REPORT,
        "study": "three-direction offset-table fairing (curvature-energy minimisation)",
        "criterion": (
            "deviation <= tolerance, curvature sign-flip count non-increasing over the "
            "accepted iterations and strictly lower overall, on >= 2 hull tables"
        ),
        "order": order,
        "feature_modes": list(feature_modes),
        "cases": cases,
        "cases_total": len(cases),
        "cases_meeting_criterion_at_a_tolerance": sum(1 for case in cases if case["accepted_at"]),
        "cases_meeting_criterion_at_a_tolerance_by_mode": {
            mode: sum(1 for case in cases if case["accepted_at_by_mode"][mode])
            for mode in feature_modes
        },
    }


def _grid_and_flat(table: ShipOffsetTable) -> tuple[_Grid, list[float]]:
    grid = _grid_of(table)
    return grid, grid.flat()


# ---------------------------------------------------------------------------
# command line
# ---------------------------------------------------------------------------


def _plot(result: FairingResult, path: Path, *, form: str) -> Path:
    """Waterlines, stations and one buttock family before/after, plus the flip history."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    before = result.table_before
    after = result.table_after
    half_beam = 0.5 * before.beam
    figure, axes = plt.subplots(2, 2, figsize=(11.0, 8.0))
    (ax_wl, ax_st), (ax_bt, ax_hist) = axes

    picks = sorted({0, before.n_waterlines // 4, before.n_waterlines // 2,
                    3 * before.n_waterlines // 4, before.n_waterlines - 1})
    for j in picks:
        ax_wl.plot(before.stations, [row[j] for row in before.half_breadths], "o",
                   ms=2.5, color="tab:gray", alpha=0.8)
        ax_wl.plot(after.stations, [row[j] for row in after.half_breadths], "-",
                   lw=1.4, label=f"z = {before.waterlines[j]:.2f} m")
    ax_wl.set_xlabel("x [m]")
    ax_wl.set_ylabel("y [m]")
    ax_wl.set_title("waterlines: ordinates (dots, before) vs faired (lines)")
    ax_wl.legend(fontsize=7)

    station_picks = sorted({before.n_stations // 2, before.n_stations // 3,
                            2 * before.n_stations // 3})
    for i in station_picks:
        ax_st.plot(before.waterlines, before.half_breadths[i], "o", ms=2.5,
                   color="tab:gray", alpha=0.8)
        ax_st.plot(after.waterlines, after.half_breadths[i], "-", lw=1.4,
                   label=f"x = {before.stations[i]:.1f} m")
    ax_st.set_xlabel("z [m]")
    ax_st.set_ylabel("y [m]")
    ax_st.set_title("stations")
    ax_st.legend(fontsize=7)

    grid_before = _grid_of(before)
    grid_after = _grid_of(after)
    levels = _buttock_levels(grid_before)
    for level in levels:
        # Plot every branch of the level set on both tables.  Zipping one "before" list
        # against one "after" list would mislabel the figure as soon as the fairing
        # changes how many branches a level has -- which it does as soon as a section
        # loses its tumblehome.
        for curve in buttock_curves(grid_before.eta, grid_before.xi, grid_before.zeta,
                                    levels=(level,)):
            ax_bt.plot([p[0] for p in curve], [p[1] for p in curve], "o--",
                       ms=3, lw=0.8, color="tab:gray", alpha=0.9)
        for curve in buttock_curves(grid_after.eta, grid_after.xi, grid_after.zeta,
                                    levels=(level,)):
            ax_bt.plot([p[0] for p in curve], [p[1] for p in curve], "-", lw=1.4,
                       label=f"y = {level * half_beam:.3f} m")
    ax_bt.set_xlabel("x / L")
    ax_bt.set_ylabel("z / T")
    ax_bt.set_title("buttocks (level sets), dashed = before")
    ax_bt.legend(fontsize=7)

    # The histogram follows the accepted path and ends on the delivered table; a rejected
    # trial is drawn off the path, because it is not a result the caller ever sees.
    totals = [result.iterations[0].total_flips_before]
    labels = ["before"]
    rejected: list[tuple[int, int, int]] = []
    for iteration in result.iterations:
        if iteration.accepted:
            totals.append(iteration.total_flips_after)
            labels.append(f"iter {iteration.iteration}")
        else:
            rejected.append((len(totals) - 1, iteration.total_flips_after, iteration.iteration))
    totals.append(result.flips_after)
    labels.append("delivered")
    ax_hist.plot(range(len(totals)), totals, "o-")
    for position, value, index in rejected:
        ax_hist.plot([position], [value], "x", color="tab:red", ms=9,
                     label=f"iter {index} rejected")
    if rejected:
        ax_hist.legend(fontsize=7)
    ax_hist.set_xticks(range(len(totals)))
    ax_hist.set_xticklabels(labels, fontsize=8)
    ax_hist.set_ylabel("curvature sign flips (three families)")
    ax_hist.set_title(f"{form}: deviation {result.deviation_max:.4f} m "
                      f"<= {result.tolerance:.4f} m")
    ax_hist.grid(alpha=0.3)

    figure.suptitle(f"three-direction fairing of {before.name}")
    figure.tight_layout()
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(target, dpi=150)
    plt.close(figure)
    return target


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="tensorcad-lines-fairing",
        description="Three-direction fairing of an offset table (curvature-energy minimisation).",
    )
    parser.add_argument("--form", default="wigley_iii",
                        help="built-in hull form to tabulate (offset_forms.builtin_forms)")
    parser.add_argument("--stations", type=int, default=21)
    parser.add_argument("--waterlines", type=int, default=17)
    parser.add_argument("--waterline-exponent", type=float, default=1.0,
                        help="clustering exponent of the waterline grid")
    parser.add_argument("--noise-fraction", type=float, default=0.004,
                        help="ordinate jitter amplitude as a fraction of the half-beam")
    parser.add_argument("--tolerance-fraction", type=float, default=0.008,
                        help="deviation budget as a fraction of the half-beam (single case)")
    parser.add_argument("--tolerance-fractions", default="0.004,0.008,0.016,0.032",
                        help="comma-separated deviation budgets for --study")
    parser.add_argument("--order", type=int, default=2, choices=(2, 3),
                        help="2 = curvature energy, 3 = curvature-variation energy")
    parser.add_argument("--iterations", type=int, default=3)
    parser.add_argument("--seed", type=int, default=20261114)
    parser.add_argument("--features", default="off", choices=("off", "auto"),
                        help="declare design-feature waterline rows (single case); "
                             "'auto' uses the energy-share detection")
    parser.add_argument("--feature-modes", default="off,auto",
                        help="comma-separated feature modes for --study (off, auto)")
    parser.add_argument("--study", action="store_true",
                        help="run all built-in forms instead of one")
    parser.add_argument("--json", type=Path, help="write the report JSON here")
    parser.add_argument("--csv", type=Path, help="write the faired table here (CSV)")
    parser.add_argument("--plot", type=Path, help="write the figure here (needs matplotlib)")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    from .offset_forms import offset_table

    args = build_parser().parse_args(argv)
    if args.study:
        report = fairing_study(
            n_stations=args.stations,
            n_waterlines=args.waterlines,
            noise_fraction=args.noise_fraction,
            tolerance_fractions=tuple(
                float(value) for value in args.tolerance_fractions.split(",") if value.strip()
            ),
            order=args.order,
            iterations=args.iterations,
            seed=args.seed,
            waterline_exponent=args.waterline_exponent,
            feature_modes=tuple(
                mode.strip() for mode in args.feature_modes.split(",") if mode.strip()
            ),
        )
        text = json.dumps(report, indent=2, sort_keys=True)
        if args.json:
            args.json.parent.mkdir(parents=True, exist_ok=True)
            args.json.write_text(text + "\n", encoding="utf-8")
        print(text)
        return 0 if report["cases_meeting_criterion_at_a_tolerance"] == report["cases_total"] else 1

    clean = offset_table(
        args.form, n_stations=args.stations, n_waterlines=args.waterlines,
        waterline_exponent=args.waterline_exponent,
    )
    half_beam = 0.5 * clean.beam
    unfair = perturb_table(clean, amplitude=args.noise_fraction * half_beam, seed=args.seed)
    result = fair_offset_table(
        unfair,
        tolerance=args.tolerance_fraction * half_beam,
        iterations=args.iterations,
        order=args.order,
        reference=clean,
        reference_name=f"{args.form} analytic ordinates (manual-fairing stand-in)",
        features=args.features,
    )
    payload = result.to_dict()
    payload["form"] = args.form
    text = json.dumps(payload, indent=2, sort_keys=True)
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(text + "\n", encoding="utf-8")
    if args.csv:
        from .offset_loft import write_offset_csv

        args.csv.parent.mkdir(parents=True, exist_ok=True)
        write_offset_csv(result.table_after, args.csv)
    if args.plot:
        _plot(result, args.plot, form=args.form)
    print(f"{args.form}: "
          f"flips {result.flips_before} -> {result.flips_after}, "
          f"deviation {result.deviation_max:.4f} m (budget {result.tolerance:.4f} m), "
          f"accepted={result.accepted()}")
    return 0 if result.accepted() else 1


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
