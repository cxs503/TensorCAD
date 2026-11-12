"""Sketch constraint solving: a 2-D parametric sketch plus constraints, solved by Newton.

``GAP-CAD-006`` in ``docs/oss_inventory_cad/GAP-ANCHORS.md`` records this as the CAD
feature that was *absent*: ``parameters.py`` validates parameter domains, but nothing
turned "these two lines are horizontal and 4 m long" into coordinates, and nothing
reported how many degrees of freedom a sketch has left. The machine-decidable acceptance
criterion attached to the gap is: ten fixtures including under- and over-constrained
ones, 100 % convergence, the correct degree-of-freedom count when a sketch is
under-constrained, and positions within ``1e-10`` of the analytic value.

**The formulation.** A sketch is a parameter vector ``x`` (two per point, one per circle
radius) and a list of scalar constraints ``r(x) = 0``. Every constraint is a residual
written once in ordinary arithmetic and differentiated by forward-mode AD
(:mod:`tensorcad.autodiff`), so the Jacobian ``J = dr/dx`` is exact rather than a finite
difference -- which matters, because a differenced Jacobian stalls Newton at the
difference's truncation error and the ``1e-10`` acceptance bound is below it.

**Why Levenberg-Marquardt rather than plain Newton.** A sketch is generally *not* a
square, non-singular system: under-constrained sketches (fewer independent rows than
parameters) and over-constrained ones (redundant rows) are the normal cases, not errors.
The damped normal equations

    (J^T J + lambda diag(J^T J + 1)) dx = -J^T r

exist for any rank and hand back the minimum-norm step near singular directions, so a
sketch with free parameters converges instead of dividing by a zero pivot. ``lambda`` is
adapted by trial: a step that lowers ``sum r^2`` is accepted and the damping is cut by
three, otherwise the damping is raised by ten and the step retried. That is the whole
solver; there is no line search, because the Levenberg-Marquardt trial *is* the search.

**Diagnostics are the second half of the deliverable.** After convergence the Jacobian is
row-normalised (so the rank tolerance does not depend on whether the sketch is drawn in
metres or millimetres) and reduced to row echelon form. Then

* ``dof = n_parameters - rank(J)`` -- the number of independent motions left;
* a consistent but non-independent row set (``rows > rank``, residuals driven to zero) is
  reported as *redundant*, not as a failure, because drawing the same dimension twice is
  a modelling convenience, not an error;
* residual rows that cannot be driven to zero are *inconsistent* -- the over-constrained
  case where two dimensions disagree -- and the offending constraints are named.

The null-space basis of the Jacobian is returned too, so "3 degrees of freedom" comes
with the three motions: for a line with only a horizontal constraint those are the two
translations plus the split of the two endpoints, which is exactly what a designer needs
to see to know what to dimension next.

**What this does not do.** No 3-D constraints, no inequality/dimension-scheme ordering, no
constraint priority when a system is inconsistent beyond naming the blocking rows, and no
global solve: Newton is local, so a mirror-image branch is selected by the starting guess
(this is why every fixture states its initial coordinates). Tangency is written with a
side convention (:class:`LineCircleTangency`) because an unsigned distance is not
differentiable where it vanishes.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import ClassVar, Iterable, Sequence

from .autodiff import Dual, atan2, sqrt

__all__ = [
    "Angle",
    "Coincident",
    "Constraint",
    "Distance",
    "EqualLength",
    "Fix",
    "FixX",
    "FixY",
    "Horizontal",
    "InconsistentConstraintsError",
    "LineCircleTangency",
    "Parallel",
    "Perpendicular",
    "PointLineDistance",
    "PointOnCircle",
    "PointOnLine",
    "Radius",
    "Sketch",
    "SketchCircle",
    "SketchError",
    "SketchLine",
    "SketchPoint",
    "SolveResult",
    "Vertical",
]


class SketchError(Exception):
    """A sketch that cannot be solved as asked."""


class InconsistentConstraintsError(SketchError):
    """The constraint set has no solution from this starting point.

    Carries the residual that could not be removed and the constraints that own it, so a
    caller can say *which* two dimensions disagree instead of reporting "solve failed".
    """

    def __init__(self, message: str, *, blocking: Sequence[str], max_residual: float,
                 result: "SolveResult | None" = None) -> None:
        super().__init__(message)
        self.blocking = tuple(blocking)
        self.max_residual = max_residual
        self.result = result


# ---------------------------------------------------------------------------------------
# entities
# ---------------------------------------------------------------------------------------
@dataclass(frozen=True)
class SketchPoint:
    """A free point: parameters ``(2 * index, 2 * index + 1)`` in the vector."""

    index: int
    name: str


@dataclass(frozen=True)
class SketchLine:
    """A line *through* two points, not a separate entity: constraining one endpoint
    moves the shared point, which is what keeps a closed profile closed."""

    start: SketchPoint
    end: SketchPoint
    name: str = ""


@dataclass(frozen=True)
class SketchCircle:
    """A circle is a centre point plus one radius parameter."""

    center: SketchPoint
    radius_index: int
    name: str = ""

    @property
    def label(self) -> str:
        return self.name or f"circle@{self.radius_index}"


class _Eval:
    """The parameter vector as dual numbers, with entity accessors."""

    __slots__ = ("params",)

    def __init__(self, values: Sequence[float]) -> None:
        self.params = [Dual.variable(i, v) for i, v in enumerate(values)]

    def xy(self, point: SketchPoint) -> tuple[Dual, Dual]:
        base = 2 * point.index
        return self.params[base], self.params[base + 1]

    def radius(self, circle: SketchCircle) -> Dual:
        return self.params[circle.radius_index]


def _direction(line: SketchLine, ev: _Eval) -> tuple[Dual, Dual]:
    ax, ay = ev.xy(line.start)
    bx, by = ev.xy(line.end)
    return bx - ax, by - ay


def _length(dx: Dual, dy: Dual) -> Dual:
    return sqrt(dx * dx + dy * dy)


def _cross(a: tuple[Dual, Dual], b: tuple[Dual, Dual]) -> Dual:
    return a[0] * b[1] - a[1] * b[0]


def _dot(a: tuple[Dual, Dual], b: tuple[Dual, Dual]) -> Dual:
    return a[0] * b[0] + a[1] * b[1]


# ---------------------------------------------------------------------------------------
# constraints
# ---------------------------------------------------------------------------------------
class Constraint:
    """One or more residual rows. ``label`` is what a failure message will print."""

    kind: ClassVar[str] = "constraint"

    def rows(self, ev: _Eval) -> list[Dual]:  # pragma: no cover - interface
        raise NotImplementedError

    def label(self) -> str:
        name = getattr(self, "name", "")
        return name or self.kind

    def __str__(self) -> str:  # pragma: no cover - diagnostics
        return self.label()


@dataclass(frozen=True)
class Coincident(Constraint):
    """Two points occupy the same location (a shared corner drawn as two points)."""

    first: SketchPoint
    second: SketchPoint
    name: str = ""
    kind: ClassVar[str] = "coincident"

    def rows(self, ev: _Eval) -> list[Dual]:
        ax, ay = ev.xy(self.first)
        bx, by = ev.xy(self.second)
        return [ax - bx, ay - by]


@dataclass(frozen=True)
class Fix(Constraint):
    """Pin a point: the gauge choice that removes rigid-body freedom."""

    point: SketchPoint
    x: float
    y: float
    name: str = ""
    kind: ClassVar[str] = "fix"

    def rows(self, ev: _Eval) -> list[Dual]:
        px, py = ev.xy(self.point)
        return [px - self.x, py - self.y]


@dataclass(frozen=True)
class FixX(Constraint):
    point: SketchPoint
    x: float
    name: str = ""
    kind: ClassVar[str] = "fix-x"

    def rows(self, ev: _Eval) -> list[Dual]:
        px, _ = ev.xy(self.point)
        return [px - self.x]


@dataclass(frozen=True)
class FixY(Constraint):
    point: SketchPoint
    y: float
    name: str = ""
    kind: ClassVar[str] = "fix-y"

    def rows(self, ev: _Eval) -> list[Dual]:
        _, py = ev.xy(self.point)
        return [py - self.y]


@dataclass(frozen=True)
class Distance(Constraint):
    """Distance between two points equals ``value`` (one equation, two branches)."""

    first: SketchPoint
    second: SketchPoint
    value: float
    name: str = ""
    kind: ClassVar[str] = "distance"

    def rows(self, ev: _Eval) -> list[Dual]:
        ax, ay = ev.xy(self.first)
        bx, by = ev.xy(self.second)
        return [_length(bx - ax, by - ay) - self.value]


@dataclass(frozen=True)
class Horizontal(Constraint):
    line: SketchLine
    name: str = ""
    kind: ClassVar[str] = "horizontal"

    def rows(self, ev: _Eval) -> list[Dual]:
        _, dy = _direction(self.line, ev)
        return [dy]


@dataclass(frozen=True)
class Vertical(Constraint):
    line: SketchLine
    name: str = ""
    kind: ClassVar[str] = "vertical"

    def rows(self, ev: _Eval) -> list[Dual]:
        dx, _ = _direction(self.line, ev)
        return [dx]


@dataclass(frozen=True)
class PointOnLine(Constraint):
    """Point lies on the infinite line through the two endpoints.

    The residual is the *signed* perpendicular distance ``cross / |d|``: it is a length,
    so its tolerance and the rank of its row are both expressed in metres, and it stays
    bounded far from the line (unlike a raw cross product, which grows with the line
    length and makes the rank tolerance scale dependent).
    """

    point: SketchPoint
    line: SketchLine
    name: str = ""
    kind: ClassVar[str] = "point-on-line"

    def rows(self, ev: _Eval) -> list[Dual]:
        return [_signed_distance(self.point, self.line, ev)]


def _signed_distance(point: SketchPoint, line: SketchLine, ev: _Eval) -> Dual:
    px, py = ev.xy(point)
    ax, ay = ev.xy(line.start)
    dx, dy = _direction(line, ev)
    length = _length(dx, dy)
    return (dx * (py - ay) - dy * (px - ax)) / length


@dataclass(frozen=True)
class PointLineDistance(Constraint):
    """Signed perpendicular distance from a point to a line equals ``value``."""

    point: SketchPoint
    line: SketchLine
    value: float
    name: str = ""
    kind: ClassVar[str] = "point-line-distance"

    def rows(self, ev: _Eval) -> list[Dual]:
        return [_signed_distance(self.point, self.line, ev) - self.value]


@dataclass(frozen=True)
class Parallel(Constraint):
    """Directions are parallel; with unit directions the residual is ``sin(angle)``."""

    first: SketchLine
    second: SketchLine
    name: str = ""
    kind: ClassVar[str] = "parallel"

    def rows(self, ev: _Eval) -> list[Dual]:
        d1 = _direction(self.first, ev)
        d2 = _direction(self.second, ev)
        return [_cross(d1, d2) / (_length(*d1) * _length(*d2))]


@dataclass(frozen=True)
class Perpendicular(Constraint):
    """Directions are perpendicular; with unit directions the residual is ``cos(angle)``."""

    first: SketchLine
    second: SketchLine
    name: str = ""
    kind: ClassVar[str] = "perpendicular"

    def rows(self, ev: _Eval) -> list[Dual]:
        d1 = _direction(self.first, ev)
        d2 = _direction(self.second, ev)
        return [_dot(d1, d2) / (_length(*d1) * _length(*d2))]


@dataclass(frozen=True)
class Angle(Constraint):
    """Signed angle from ``first`` to ``second`` in degrees.

    Signed and measured with ``atan2``, so the residual is continuous through the right
    angle but wraps by ``2 pi``; the branch is the one the starting geometry is on.
    """

    first: SketchLine
    second: SketchLine
    degrees: float
    name: str = ""
    kind: ClassVar[str] = "angle"

    def rows(self, ev: _Eval) -> list[Dual]:
        d1 = _direction(self.first, ev)
        d2 = _direction(self.second, ev)
        return [atan2(_cross(d1, d2), _dot(d1, d2)) - math.radians(self.degrees)]


@dataclass(frozen=True)
class EqualLength(Constraint):
    first: SketchLine
    second: SketchLine
    name: str = ""
    kind: ClassVar[str] = "equal-length"

    def rows(self, ev: _Eval) -> list[Dual]:
        return [_length(*_direction(self.first, ev)) - _length(*_direction(self.second, ev))]


@dataclass(frozen=True)
class Radius(Constraint):
    circle: SketchCircle
    value: float
    name: str = ""
    kind: ClassVar[str] = "radius"

    def rows(self, ev: _Eval) -> list[Dual]:
        return [ev.radius(self.circle) - self.value]


@dataclass(frozen=True)
class PointOnCircle(Constraint):
    point: SketchPoint
    circle: SketchCircle
    name: str = ""
    kind: ClassVar[str] = "point-on-circle"

    def rows(self, ev: _Eval) -> list[Dual]:
        px, py = ev.xy(self.point)
        cx, cy = ev.xy(self.circle.center)
        return [_length(px - cx, py - cy) - ev.radius(self.circle)]


@dataclass(frozen=True)
class LineCircleTangency(Constraint):
    """Line tangent to a circle, with the side chosen explicitly.

    ``side="left"`` means the centre is on the left of the ``start -> end`` direction, so
    the signed distance equals ``+r``; ``side="right"`` is the mirror. An unsigned
    ``|distance| - r`` would be differentiable everywhere except exactly at the solution,
    where its kink is the one place the solver must pass through.
    """

    line: SketchLine
    circle: SketchCircle
    side: str = "left"
    name: str = ""
    kind: ClassVar[str] = "line-circle-tangency"

    def __post_init__(self) -> None:
        if self.side not in ("left", "right"):
            raise ValueError(f"side must be 'left' or 'right', got {self.side!r}")

    def rows(self, ev: _Eval) -> list[Dual]:
        sign = 1.0 if self.side == "left" else -1.0
        return [_signed_distance(self.circle.center, self.line, ev) - sign * ev.radius(self.circle)]


# ---------------------------------------------------------------------------------------
# the sketch and the solver
# ---------------------------------------------------------------------------------------
@dataclass(frozen=True)
class SolveResult:
    """Outcome of :meth:`Sketch.solve`, including the degrees of freedom that remain."""

    #: ``"fully_constrained"`` (residual inside ``tol`` and no freedom left),
    #: ``"under_constrained"`` (residual inside ``tol``, freedom left -- ``null_space``
    #: says which), ``"not_converged"`` (residual above ``tol`` but below
    #: ``residual_tol``: the structural verdict is withheld because these coordinates do
    #: not satisfy the sketch), or ``"inconsistent"`` (residual above ``residual_tol``).
    status: str
    dof: int
    rank: int
    parameters: int
    residual_rows: int
    redundant_rows: int
    max_residual: float
    residuals: tuple[float, ...]
    constraint_labels: tuple[str, ...]
    iterations: int
    converged: bool
    null_space: tuple[tuple[float, ...], ...]
    history: tuple[float, ...]
    values: tuple[float, ...]
    positions: dict[str, tuple[float, float]]
    radii: dict[str, float]
    blocking: tuple[str, ...] = field(default_factory=tuple)
    #: True when a row had no derivative at the caller's start and the solve was restarted
    #: from :data:`SINGULAR_START_NUDGE` away from it (separate coincident points).
    nudged_from_singular_start: bool = False

    @property
    def fully_constrained(self) -> bool:
        return self.status == "fully_constrained"

    @property
    def under_constrained(self) -> bool:
        return self.status == "under_constrained"

    @property
    def inconsistent(self) -> bool:
        return self.status == "inconsistent"

    @property
    def not_converged(self) -> bool:
        """The solve stopped above ``tol`` without being provably inconsistent.

        The coordinates are the best the iteration reached, *not* a solution: ``dof`` is
        the local rank there. A caller that needs a solution must treat this as a
        failure -- ``fully_constrained`` can no longer be true here.
        """
        return self.status == "not_converged"

    @property
    def solved(self) -> bool:
        """Residual inside ``tol``: the coordinates satisfy every constraint row."""
        return self.converged


def _evaluate(sketch: "Sketch", values: Sequence[float]) -> tuple[list[float], list[list[float]]]:
    """Residual vector and dense Jacobian at ``values`` (one pass over the constraints)."""
    ev = _Eval(values)
    residuals: list[float] = []
    jacobian: list[list[float]] = []
    n = len(values)
    for constraint in sketch.constraints:
        for row in constraint.rows(ev):
            residuals.append(row.value)
            dense = [0.0] * n
            for index, value in row.grad.items():
                dense[index] = value
            jacobian.append(dense)
    return residuals, jacobian


#: Deterministic displacement used when a row has no derivative at the starting point.
#: A ``Distance`` between two points that both start at the origin is a length row whose
#: gradient is ``d / |d|`` with ``|d| = 0``: the geometry is well posed, only the *start*
#: is not, so the solver steps off the singular point instead of refusing the sketch.
#: ``1e-6`` is nine orders above the default convergence tolerance and four orders below
#: anything a drawing dimension would notice, so the solution it reaches is the same one
#: the caller would have got from a non-degenerate start.
SINGULAR_START_NUDGE = 1e-6


def _nudge_from_singular(values: Sequence[float]) -> list[float]:
    """Step off a singular start, deterministically and symmetrically.

    The displacements are ``+1, -2, +3, -4, ...`` times :data:`SINGULAR_START_NUDGE`: every
    parameter gets a *distinct* magnitude, so any two parameters that hold numerically
    equal values -- the coincidence that made the row singular -- are separated whatever
    pair they are. (An alternating *sign* alone does not do this: the two coordinates of
    one point and of the next point would get the same pattern and move together.) The
    magnitudes grow linearly with the parameter index, which the module's own cost note
    keeps small: this solver is for sketches with tens of parameters, so the largest
    displacement is a few tens of micrometres.
    """
    return [
        value + SINGULAR_START_NUDGE * (index + 1) * (1.0 if index % 2 == 0 else -1.0)
        for index, value in enumerate(values)
    ]


def _singular_label(sketch: "Sketch", values: Sequence[float]) -> str:
    """Label of the first constraint whose rows cannot be evaluated at ``values``."""
    ev = _Eval(list(values))
    for constraint in sketch.constraints:
        try:
            constraint.rows(ev)
        except ZeroDivisionError:
            return constraint.label()
    return "unattributed"


def _singular_message(sketch: "Sketch", values: Sequence[float], cause: Exception) -> str:
    return (
        f"constraint {_singular_label(sketch, values)!r} has no derivative at these "
        f"coordinates: a zero-length line has no direction for a length or direction row, "
        f"and a zero-radius circle has no tangent ({cause}). Give the entities distinct "
        "coordinates, or state the constraint without a direction."
    )


def _normal_equations(jacobian: Sequence[Sequence[float]], residuals: Sequence[float],
                      damping: float) -> tuple[list[list[float]], list[float]]:
    """``A = J^T J + lambda diag(J^T J + 1)`` and ``b = -J^T r``.

    The ``+1`` inside the diagonal keeps a parameter that no residual touches from
    producing a zero pivot: its damped entry is ``lambda``, so its step is zero, which is
    the right answer (nothing constrains it) rather than a division by zero.
    """
    n = len(jacobian[0]) if jacobian else 0
    a = [[0.0] * n for _ in range(n)]
    b = [0.0] * n
    for row, residual in zip(jacobian, residuals):
        for i in range(n):
            ri = row[i]
            if ri == 0.0:
                continue
            b[i] -= ri * residual
            ai = a[i]
            for j in range(i, n):
                ai[j] += ri * row[j]
    for i in range(n):
        for j in range(i):
            a[i][j] = a[j][i]
        a[i][i] += damping * (a[i][i] + 1.0)
    return a, b


class _Singular(Exception):
    """The damped normal equations were not positive definite; raise the damping."""


def _cholesky_solve(a: Sequence[Sequence[float]], b: Sequence[float]) -> list[float]:
    """Solve ``A x = b`` for symmetric positive definite ``A`` by Cholesky."""
    n = len(b)
    if n == 0:
        return []
    lower = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            s = a[i][j] - sum(lower[i][k] * lower[j][k] for k in range(j))
            if i == j:
                if s <= 0.0:
                    raise _Singular
                lower[i][j] = math.sqrt(s)
            else:
                lower[i][j] = s / lower[j][j]
    y = [0.0] * n
    for i in range(n):
        y[i] = (b[i] - sum(lower[i][k] * y[k] for k in range(i))) / lower[i][i]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - sum(lower[k][i] * x[k] for k in range(i + 1, n))) / lower[i][i]
    return x


def _row_echelon(jacobian: Sequence[Sequence[float]], tolerance: float):
    """Reduced row echelon form of the row-normalised Jacobian.

    Row normalisation first: a row's scale is an artefact of the units in the residual
    (a ``distance`` row is ``O(1)`` per metre, a ``fix`` row is ``1``), and comparing a
    pivot against a fixed absolute tolerance without it makes the reported rank depend on
    the drawing's units. Row operations preserve the null space, so the basis returned
    from this matrix is the null space of the original Jacobian.
    """
    matrix: list[list[float]] = []
    for row in jacobian:
        norm = math.sqrt(sum(v * v for v in row))
        matrix.append([v / norm for v in row] if norm > 0.0 else list(row))
    rows = len(matrix)
    cols = len(matrix[0]) if rows else 0
    pivots: list[int] = []
    current = 0
    for column in range(cols):
        best = max(range(current, rows), key=lambda i: abs(matrix[i][column]), default=None)
        if best is None or abs(matrix[best][column]) <= tolerance:
            continue
        matrix[current], matrix[best] = matrix[best], matrix[current]
        pivot = matrix[current][column]
        matrix[current] = [v / pivot for v in matrix[current]]
        for i in range(rows):
            if i != current and matrix[i][column] != 0.0:
                factor = matrix[i][column]
                matrix[i] = [vi - factor * vp for vi, vp in zip(matrix[i], matrix[current])]
        pivots.append(column)
        current += 1
        if current == rows:
            break
    return current, matrix, pivots


def _null_space(rref: Sequence[Sequence[float]], pivots: Sequence[int], cols: int) -> tuple[tuple[float, ...], ...]:
    """Basis of the null space, one vector per free column of the reduced Jacobian."""
    free = [c for c in range(cols) if c not in set(pivots)]
    basis: list[tuple[float, ...]] = []
    for column in free:
        vector = [0.0] * cols
        vector[column] = 1.0
        for row_index, pivot_column in enumerate(pivots):
            vector[pivot_column] = -rref[row_index][column]
        basis.append(tuple(vector))
    return tuple(basis)


class Sketch:
    """A 2-D parametric sketch: points, lines through them, circles, and constraints."""

    def __init__(self) -> None:
        self._values: list[float] = []
        self.points: list[SketchPoint] = []
        self.lines: list[SketchLine] = []
        self.circles: list[SketchCircle] = []
        self.constraints: list[Constraint] = []

    # -- construction -----------------------------------------------------------------
    def point(self, x: float, y: float, name: str | None = None) -> SketchPoint:
        index = len(self.points)
        self._values.extend((float(x), float(y)))
        created = SketchPoint(index, name or f"P{index}")
        self.points.append(created)
        return created

    def line(self, start: SketchPoint, end: SketchPoint, name: str | None = None) -> SketchLine:
        created = SketchLine(start, end, name or f"L{len(self.lines)}")
        self.lines.append(created)
        return created

    def circle(self, center: SketchPoint, radius: float, name: str | None = None) -> SketchCircle:
        index = len(self._values)
        self._values.append(float(radius))
        created = SketchCircle(center, index, name or f"C{len(self.circles)}")
        self.circles.append(created)
        return created

    def add(self, *constraints: Constraint) -> "Sketch":
        for constraint in constraints:
            if isinstance(constraint, (list, tuple)):
                self.add(*constraint)
            else:
                self.constraints.append(constraint)
        return self

    constrain = add

    # -- state ------------------------------------------------------------------------
    @property
    def parameter_count(self) -> int:
        return len(self._values)

    @property
    def constraint_labels(self) -> tuple[str, ...]:
        labels: list[str] = []
        for constraint in self.constraints:
            labels.extend([constraint.label()] * len(constraint.rows(_Eval(self._values))))
        return tuple(labels)

    def values(self) -> tuple[float, ...]:
        return tuple(self._values)

    def set_values(self, values: Sequence[float]) -> None:
        if len(values) != len(self._values):
            raise SketchError(
                f"expected {len(self._values)} parameters, got {len(values)}; a vector "
                "from another sketch cannot be applied here"
            )
        self._values = [float(v) for v in values]

    def position(self, point: SketchPoint) -> tuple[float, float]:
        base = 2 * point.index
        return self._values[base], self._values[base + 1]

    def set_position(self, point: SketchPoint, x: float, y: float) -> None:
        base = 2 * point.index
        self._values[base] = float(x)
        self._values[base + 1] = float(y)

    def radius(self, circle: SketchCircle) -> float:
        return self._values[circle.radius_index]

    def set_radius(self, circle: SketchCircle, radius: float) -> None:
        self._values[circle.radius_index] = float(radius)

    def residuals(self, values: Sequence[float] | None = None) -> tuple[float, ...]:
        """Constraint residuals at ``values`` (default: the current coordinates).

        Public because a consumer downstream of a solve needs to re-check the answer for
        itself: pointing this at the returned coordinates and finding zero is an
        independent verification, not a repetition of the solver's own bookkeeping.
        """
        evaluated = self._values if values is None else list(values)
        try:
            residuals, _ = _evaluate(self, evaluated)
        except ZeroDivisionError as exc:
            raise SketchError(_singular_message(self, evaluated, exc)) from exc
        return tuple(residuals)

    def jacobian(self, values: Sequence[float] | None = None) -> tuple[tuple[float, ...], ...]:
        """Dense constraint Jacobian ``dr/dx`` at ``values`` (exact, by forward-mode AD)."""
        evaluated = self._values if values is None else list(values)
        try:
            _, jacobian = _evaluate(self, evaluated)
        except ZeroDivisionError as exc:
            raise SketchError(_singular_message(self, evaluated, exc)) from exc
        return tuple(tuple(row) for row in jacobian)

    def positions(self) -> dict[str, tuple[float, float]]:
        return {point.name: self.position(point) for point in self.points}

    def radii(self) -> dict[str, float]:
        return {circle.label: self.radius(circle) for circle in self.circles}

    def perturb(self, scale: float, rng) -> None:
        """Jitter every point (and radius) by ``uniform(-scale, scale)``.

        Used for the convergence study: a solver that only works from the coordinates its
        own fixture author typed is not a solver. Deterministic for a seeded ``rng``.
        """
        for i in range(len(self._values)):
            self._values[i] += rng.uniform(-scale, scale)

    # -- solve ------------------------------------------------------------------------
    def solve(
        self,
        *,
        initial: Sequence[float] | None = None,
        tol: float = 1e-12,
        residual_tol: float = 1e-9,
        max_iterations: int = 100,
        strict: bool = True,
        rank_tolerance: float = 1e-9,
    ) -> SolveResult:
        """Drive the constraints to zero and report rank, DOF and any blocking rows.

        ``tol`` is the convergence threshold on the largest residual; ``residual_tol`` is
        the threshold above which an unconverged residual is called inconsistent rather
        than merely unfinished. ``strict`` raises
        :class:`InconsistentConstraintsError` on an inconsistent system; with
        ``strict=False`` the same information comes back in the result, which is what the
        fixtures use to assert on the failure mode itself.
        """
        x = [float(v) for v in (self._values if initial is None else initial)]
        if len(x) != len(self._values):
            raise SketchError(
                f"initial vector has {len(x)} parameters, the sketch has {len(self._values)}"
            )

        damping = 1e-8
        iterations = 0
        nudged_from_singular_start = False
        try:
            residuals, jacobian = _evaluate(self, x)
        except ZeroDivisionError as first:
            # The start is on the singular set of a length row (coincident points, a
            # zero-length line). The *geometry* may still be well posed -- a distance
            # between two points that happen to start together is an ordinary request --
            # so move deterministically off the singular point once and retry.
            x = _nudge_from_singular(x)
            nudged_from_singular_start = True
            try:
                residuals, jacobian = _evaluate(self, x)
            except ZeroDivisionError as second:
                raise SketchError(_singular_message(self, x, second)) from first
        max_residual = max((abs(v) for v in residuals), default=0.0)
        history: list[float] = [max_residual]
        while max_residual > tol and iterations < max_iterations:
            iterations += 1
            cost = sum(v * v for v in residuals)
            step: list[float] | None = None
            for _ in range(40):
                a, b = _normal_equations(jacobian, residuals, damping)
                try:
                    candidate = _cholesky_solve(a, b)
                except _Singular:
                    damping *= 10.0
                    continue
                trial = [xi + di for xi, di in zip(x, candidate)]
                try:
                    trial_residuals, _ = _evaluate(self, trial)
                except ZeroDivisionError:
                    # The step collapsed a line (or a radius) onto the singular set. That
                    # step cannot be scored, so treat it as a rejected trial and damp.
                    damping *= 10.0
                    continue
                if sum(v * v for v in trial_residuals) < cost:
                    x = trial
                    step = candidate
                    damping = max(damping / 3.0, 1e-12)
                    break
                damping *= 10.0
            if step is None:  # no damping level improved the cost: a local minimum
                break
            residuals, jacobian = _evaluate(self, x)
            max_residual = max((abs(v) for v in residuals), default=0.0)
            history.append(max_residual)

        residuals, jacobian = _evaluate(self, x)
        max_residual = max((abs(v) for v in residuals), default=0.0)
        converged = max_residual <= tol
        rank, rref, pivots = _row_echelon(jacobian, rank_tolerance)
        parameters = len(x)
        rows = len(residuals)
        dof = parameters - rank

        labels: list[str] = []
        for constraint in self.constraints:
            labels.extend([constraint.label()] * len(constraint.rows(_Eval(x))))
        # An inconsistent set names the rows that cannot be driven home; a merely
        # unfinished solve names the rows that are still above the convergence tolerance.
        blocking_threshold = residual_tol if max_residual > residual_tol else tol
        blocking = tuple(sorted({label for label, value in zip(labels, residuals)
                                  if abs(value) > blocking_threshold}))

        if max_residual > residual_tol:
            status = "inconsistent"
        elif not converged:
            # dof is the rank of the Jacobian *at these coordinates*, and those
            # coordinates do not satisfy the sketch. Reporting "fully_constrained" here
            # would hand a consumer geometry that misses its dimensions, so the
            # structural verdict is withheld until the residual is inside ``tol``.
            status = "not_converged"
        elif dof > 0:
            status = "under_constrained"
        else:
            status = "fully_constrained"

        self._values = list(x)
        result = SolveResult(
            status=status,
            dof=dof,
            rank=rank,
            parameters=parameters,
            residual_rows=rows,
            redundant_rows=rows - rank,
            max_residual=max_residual,
            residuals=tuple(residuals),
            constraint_labels=tuple(labels),
            iterations=iterations,
            converged=converged,
            null_space=_null_space(rref, pivots, parameters),
            history=tuple(history),
            values=tuple(x),
            positions=self.positions(),
            radii=self.radii(),
            blocking=blocking,
            nudged_from_singular_start=nudged_from_singular_start,
        )
        if status == "inconsistent" and strict:
            raise InconsistentConstraintsError(
                f"constraints are inconsistent from this start: max residual "
                f"{max_residual:.3e} after {iterations} iterations; blocking rows: "
                f"{', '.join(blocking) or 'unattributed'}. Either the dimensions disagree "
                "or the starting geometry is on the wrong branch.",
                blocking=blocking,
                max_residual=max_residual,
                result=result,
            )
        return result
