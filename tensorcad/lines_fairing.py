"""Three-directional fairing of an offset table (``GAP-CAD-C3-LINES-FAIRING``).

A ship's lines are handed over as an **offset table** -- half-breadths ``y`` on a grid of
stations ``x`` and waterlines ``z`` -- and a table straight from the model basin is never
fair: digitiser scatter, half-millimetre rounding and the inevitable small kinks from
interpolating between waterlines all show up as curvature that reverses sign from one
ordinate to the next.  A designer fixes that by hand, re-drawing every **waterline**
(``y(x)`` at constant ``z``), every **station** (``y(z)`` at constant ``x``) and every
**buttock** (``z(x)`` at constant ``y``) until the curvature runs smoothly -- the
"three-directional fairing" of the classical lines plan (Watson 1998, *Practical Ship
Design*, ch. 3; Letcher 2009, ch. 6). This module does that pass numerically and, more
importantly, **measures** whether it worked.

What is being minimised
-----------------------
Per curve, the fairing energy is the **variation of curvature**, i.e. the third derivative
in the curve's own normalised parameter ``t in [0, 1]``::

    E[f] = sum_i (f_i - y_i)**2  +  lambda * sum_k (f'''_k)**2 * dt

``f''' = 0`` is a parabola, so quadratics -- the shapes a waterline and a section spend
most of their length approximating -- are in the *null space* of the penalty.  That is
deliberate: a clean parabolic table is a fixed point of this operator (measured deviation
below ``1e-12`` m on the Wigley table), whereas an operator built on ``f''`` would flatten
a perfectly good waterline towards a straight line.  ``lambda`` is not a hand-picked
number: it is found by bisection as the **largest** smoothing weight whose ordinates stay
inside the tolerance band, i.e. the smoothest curve the designer's tolerance allows.

The constraint, not the penalty, is what the acceptance asks for
----------------------------------------------------------------
``GAP-CAD-C3-LINES-FAIRING`` accepts on two machine-decidable clauses:

* every faired ordinate stays within a stated tolerance of the table it came from
  (``max |faired - input| <= offset_tolerance``, enforced by projection after every pass
  and re-asserted at the end), and
* the number of curvature sign reversals per family falls monotonically, on at least two
  hulls.

The second clause is a property of the *driver*, not of the solver: a smoothing pass can
locally overshoot and add a flip.  :func:`fair_offset_table` therefore accepts a pass only
if the total flip count does not increase (lexicographic ``(flips, energy)``), halves the
relaxation and retries otherwise, and stops when no relaxed step improves -- so the
recorded flip history is non-increasing **by construction and testable as such**.

Sign reversals are counted on the discrete circumscribed-circle curvature
(``kappa = 2 * cross(a, b) / (|a| |b| |c|)``), ignoring values below
``curvature_floor`` times the curve's maximum ``|kappa|`` so that round-off on a straight
run is not mistaken for an inflection.

The three directions
--------------------
``waterline``  rows: ``y(x)`` at fixed ``z``, ordinate faired directly.
``station``    columns: ``y(z)`` at fixed ``x``, ordinate faired directly.
``buttock``    iso-``y`` contours: the depth ``z(x)`` at which the half-breadth equals a
               level ``y_k``.  Each section is inverted for the crossing depth, ``z(x)``
               is faired, and the depth correction is mapped back to an ordinate through
               the local section slope (``dy = -slope * dz / weight``), scaled so the
               ordinate change cannot exceed the tolerance.  A buttock correction that
               would make a section non-monotone in ``z`` is refused and counted
               (``monotone_refusals``) rather than silently applied.

What is *not* done
------------------
No automatic topology change (no station/waterline insertion or deletion), no fairness
across the centre plane or the deck edge beyond what the table holds, no bulb or
transom-specific treatment (those need a patch model, not a single table), and no claim
that a formula hull's buttock lines are the fairing target -- the buttock pass has its own
band and its own audit so a reader can see exactly what it moved.  The "manual fairing
reference" here is (a) the analytic surface for the formula hulls and (b) an independent
local-polynomial smoother (:func:`hand_fairing_reference`); no human designer's marked-up
table is in this repository, and the docs say so.

Standard library only, like the rest of the geometry stack.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Sequence

from .offset_loft import OffsetProvenance, ShipOffsetTable, table_integrals

__all__ = [
    "FAIRING_FAMILIES",
    "CurveFairness",
    "FairingError",
    "FairingResult",
    "FairingSettings",
    "FamilyFairness",
    "FairingPass",
    "curvature_sign_flips",
    "fair_offset_table",
    "hand_fairing_reference",
    "line_bending_energy",
    "signed_curvature",
]

SCHEMA_FAIRING = "tensorcad.lines_fairing.result/1"

#: The three families, in the order the driver sweeps them.
FAIRING_FAMILIES: tuple[str, ...] = ("waterline", "station", "buttock")

_BIG_WEIGHT = 1.0e8


class FairingError(ValueError):
    """Raised when a table or a curve cannot be faired, with a named reason."""


# --------------------------------------------------------------------------- #
# settings
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class FairingSettings:
    """Tolerances and effort knobs; every default is pinned by a test.

    ``offset_tolerance`` is in metres and is the one number the acceptance names: no
    ordinate may move further than this from the table it came from.  ``band_fraction``
    keeps the *solver* strictly inside that band so that the final projection is a
    formality instead of a silent clamp.  ``smoothing_search`` is the ``log10 lambda``
    bracket the bisection searches; ``search_steps`` bisection halvings per curve.
    """

    offset_tolerance: float = 2.0e-3
    band_fraction: float = 0.95
    smoothing_search: tuple[float, float] = (-12.0, -1.0)
    search_steps: int = 8
    clamping_iterations: int = 8
    max_passes: int = 6
    relax: float = 1.0
    #: smoothing kernel: ``"local"`` (moving least-squares quadratic, shape-preserving) or
    #: ``"energy"`` (global third-derivative minimisation inside the band).  The measured
    #: comparison is in ``docs/lines_fairing.md``: the energy kernel over-fairs sections
    #: that are steeper than the grid can resolve.
    kernel: str = "local"
    local_window: int = 5
    #: a point whose local-quadratic residual exceeds this many robust standard deviations
    #: is treated as *shape* (a knuckle, a turn of the bilge, the first waterline of a
    #: near-rectangular section) and left where the designer put it.
    shape_threshold: float = 3.0
    families: tuple[str, ...] = FAIRING_FAMILIES
    curvature_floor: float = 0.02
    buttock_levels: int = 6
    #: share of ``offset_tolerance`` the buttock pass may spend.  A buttock is a *derived*
    #: contour of the surface, not a measured ordinate, and its noise is amplified by
    #: ``1 / (dy/dz)`` -- so it gets a minority of the budget by construction (measured:
    #: spending the full band on the buttocks buys no flips and costs accuracy; see
    #: ``docs/lines_fairing.md``).
    buttock_band_fraction: float = 0.25
    min_curve_points: int = 5

    def __post_init__(self) -> None:
        if not (self.offset_tolerance > 0.0):
            raise FairingError(
                f"offset_tolerance must be positive, got {self.offset_tolerance!r}"
            )
        if not (0.0 < self.band_fraction <= 1.0):
            raise FairingError(f"band_fraction must be in (0, 1], got {self.band_fraction!r}")
        for family in self.families:
            if family not in FAIRING_FAMILIES:
                raise FairingError(
                    f"unknown family {family!r}; expected a subset of {FAIRING_FAMILIES}"
                )
        if not (0.0 < self.relax <= 1.0):
            raise FairingError(f"relax must be in (0, 1], got {self.relax!r}")
        if self.kernel not in ("local", "energy"):
            raise FairingError(f"kernel must be 'local' or 'energy', got {self.kernel!r}")
        if self.local_window < 3 or self.local_window % 2 == 0:
            raise FairingError(
                f"local_window must be odd and >= 3, got {self.local_window!r}"
            )
        low, high = self.smoothing_search
        if not low < high:
            raise FairingError(
                f"smoothing_search must increase, got ({low!r}, {high!r})"
            )


# --------------------------------------------------------------------------- #
# curvature of a polyline
# --------------------------------------------------------------------------- #
def signed_curvature(coords: Sequence[float], values: Sequence[float]) -> tuple[float, ...]:
    """Discrete circumscribed-circle curvature at every interior point.

    ``kappa_i = 2 * cross(p_i - p_{i-1}, p_{i+1} - p_i) / (|a| |b| |c|)`` with the three
    points taken in the ``(coord, value)`` plane; the sign is positive where the curve
    turns towards ``+value`` (``value'' > 0``).  A degenerate triple (repeated coordinate
    or a straight run of zero length) yields ``nan`` instead of an infinity, and callers
    skip those.
    """
    if len(coords) != len(values):
        raise FairingError(
            f"curve: {len(coords)} coordinates against {len(values)} values"
        )
    out: list[float] = []
    for i in range(1, len(coords) - 1):
        ax = coords[i] - coords[i - 1]
        ay = values[i] - values[i - 1]
        bx = coords[i + 1] - coords[i]
        by = values[i + 1] - values[i]
        cx = coords[i + 1] - coords[i - 1]
        cy = values[i + 1] - values[i - 1]
        la = math.hypot(ax, ay)
        lb = math.hypot(bx, by)
        lc = math.hypot(cx, cy)
        if la <= 0.0 or lb <= 0.0 or lc <= 0.0:
            out.append(float("nan"))
            continue
        cross = ax * by - ay * bx
        out.append(2.0 * cross / (la * lb * lc))
    return tuple(out)


def curvature_sign_flips(
    coords: Sequence[float],
    values: Sequence[float],
    *,
    floor_ratio: float = 0.02,
) -> int:
    """Number of curvature sign reversals along a polyline.

    Curvatures below ``floor_ratio * max|kappa|`` are treated as zero and neither start
    nor end a run, so a straight stretch between two arcs is one reversal and round-off on
    a straight run is none.  A curve too short to carry an interior curvature returns 0.
    """
    kappas = signed_curvature(coords, values)
    finite = [k for k in kappas if math.isfinite(k)]
    if len(finite) < 2:
        return 0
    scale = max(abs(k) for k in finite)
    if scale <= 0.0:
        return 0
    threshold = floor_ratio * scale
    flips = 0
    previous = 0
    for kappa in kappas:
        if not math.isfinite(kappa) or abs(kappa) <= threshold:
            continue
        sign = 1 if kappa > 0.0 else -1
        if previous and sign != previous:
            flips += 1
        previous = sign
    return flips


def line_bending_energy(coords: Sequence[float], values: Sequence[float]) -> float:
    """``sum kappa**2 * ds`` along the polyline (1/m: converges under refinement).

    This is the 1-D strain energy a fair line minimises; it is reported before and after
    fairing so the improvement is a number rather than an impression.
    """
    kappas = signed_curvature(coords, values)
    total = 0.0
    pieces: list[float] = []
    for index, kappa in enumerate(kappas, start=1):
        if not math.isfinite(kappa):
            continue
        ds = 0.5 * (
            math.hypot(coords[index] - coords[index - 1], values[index] - values[index - 1])
            + math.hypot(
                coords[index + 1] - coords[index], values[index + 1] - values[index]
            )
        )
        pieces.append(kappa * kappa * ds)
    return math.fsum(pieces)


# --------------------------------------------------------------------------- #
# banded Cholesky (the smoothing system is 7-diagonal: bandwidth 3)
# --------------------------------------------------------------------------- #
def _cholesky_banded_solve(
    matrix: Sequence[Sequence[float]], rhs: Sequence[float], bandwidth: int
) -> list[float]:
    """Solve ``A x = b`` for symmetric positive-definite ``A`` of half-bandwidth ``p``.

    ``A`` is given as dense rows (only the band is read).  Cholesky keeps the same
    bandwidth, so the factorisation costs ``O(n p**2)`` instead of ``O(n**3)`` -- which
    matters because the curvature bisection solves this system dozens of times per curve.
    """
    n = len(rhs)
    p = bandwidth
    factor = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(max(0, i - p), i + 1):
            total = matrix[i][j]
            for k in range(max(0, j - p), j):
                total -= factor[i][k] * factor[j][k]
            if i == j:
                if total <= 0.0:
                    raise FairingError(
                        f"smoothing system is not positive definite at row {i} "
                        f"(pivot {total!r}); the tolerance band or the curve is degenerate"
                    )
                factor[i][i] = math.sqrt(total)
            else:
                factor[i][j] = total / factor[j][j]
    forward = [0.0] * n
    for i in range(n):
        total = rhs[i]
        for k in range(max(0, i - p), i):
            total -= factor[i][k] * forward[k]
        forward[i] = total / factor[i][i]
    solution = [0.0] * n
    for i in range(n - 1, -1, -1):
        total = forward[i]
        for k in range(i + 1, min(n, i + p + 1)):
            total -= factor[k][i] * solution[k]
        solution[i] = total / factor[i][i]
    return solution


# --------------------------------------------------------------------------- #
# one curve
# --------------------------------------------------------------------------- #
def _third_difference_operator(n: int) -> tuple[float, list[tuple[int, tuple[float, ...]]]]:
    """Third-difference rows on the normalised parameter ``t`` and the panel weight.

    With ``h = 1 / (n - 1)`` the third derivative is approximated by
    ``(f_{k+3} - 3 f_{k+2} + 3 f_{k+1} - f_k) / h**3`` and the penalty is its square
    integrated over the panel (weight ``h``).  In ``t`` units ``f'''`` has the units of
    the ordinate, so ``lambda`` is dimensionless and comparable between a 5 m model and a
    320 m tanker.
    """
    h = 1.0 / (n - 1)
    scale = 1.0 / (h * h * h)
    rows: list[tuple[int, tuple[float, ...]]] = []
    for k in range(n - 3):
        rows.append((k, (-scale, 3.0 * scale, -3.0 * scale, scale)))
    return h, rows


def _solve_with_band(
    values: Sequence[float],
    rows: Sequence[tuple[int, tuple[float, ...]]],
    panel: float,
    lam: float,
    band: float,
    iterations: int,
) -> list[float]:
    """``min sum (f - y)**2 + lam * sum (f''')**2 dt`` subject to ``|f - y| <= band``.

    The band is handled as an active set: pinned points enter the normal equations with a
    large weight at the band edge, the system is re-solved, and any new violation is
    added.  The returned vector is finally projected onto the band, so the constraint
    holds exactly (not just to the weight).
    """
    n = len(values)
    pinned: dict[int, float] = {}
    solution = list(values)
    for _ in range(iterations):
        weights = [1.0] * n
        rhs = [float(value) for value in values]
        for index, target in pinned.items():
            weights[index] = _BIG_WEIGHT
            rhs[index] = _BIG_WEIGHT * target
        matrix = [[0.0] * n for _ in range(n)]
        for i in range(n):
            matrix[i][i] += weights[i]
        for start, coefficients in rows:
            coefficient = lam * panel
            for a in range(4):
                ca = coefficient * coefficients[a]
                if ca == 0.0:
                    continue
                for b in range(4):
                    matrix[start + a][start + b] += ca * coefficients[b]
        solution = _cholesky_banded_solve(matrix, rhs, 3)
        violations = [
            index
            for index in range(n)
            if abs(solution[index] - values[index]) > band and index not in pinned
        ]
        if not violations:
            break
        for index in violations:
            edge = values[index] + math.copysign(band, solution[index] - values[index])
            pinned[index] = edge
    return [
        min(max(value, values[i] - band), values[i] + band)
        for i, value in enumerate(solution)
    ]


def _local_quadratic(
    coords: Sequence[float], values: Sequence[float], window: int
) -> list[float]:
    """Moving least-squares quadratic: the value of the local fit at every interior point.

    Over a ``window``-point stencil the fitted value is the minimiser of the local data
    misfit over ``span{1, t, t**2}``, so it reproduces any quadratic exactly (a Wigley
    waterline comes back unchanged to round-off) and averages away the grid-scale part of
    a digitised table.  Edge points are left as they are: a one-sided stencil is an
    extrapolation, not a fit.
    """
    n = len(values)
    half = window // 2
    out = [float(value) for value in values]
    for i in range(half, n - half):
        centre = coords[i]
        xs = [coords[k] - centre for k in range(i - half, i + half + 1)]
        ys = [values[k] for k in range(i - half, i + half + 1)]
        matrix = [
            [math.fsum(x ** (p + q) for x in xs) for q in range(3)] for p in range(3)
        ]
        rhs = [math.fsum((x ** p) * y for x, y in zip(xs, ys, strict=True)) for p in range(3)]
        try:
            coefficients = _solve_dense(matrix, rhs)
        except FairingError:
            continue
        out[i] = coefficients[0]
    return out


def _fair_curve(
    coords: Sequence[float],
    values: Sequence[float],
    settings: FairingSettings,
    *,
    tolerance: float | None = None,
) -> list[float]:
    """Fair one curve inside an ordinate band; returns the faired values.

    ``tolerance`` overrides the settings' band for a family whose ordinate is not the
    faired quantity (buttock: the faired quantity is a depth).

    ``kernel="local"``: one moving least-squares quadratic pass, projected onto the band.
    ``kernel="energy"``: the largest third-derivative weight whose solution stays in the
    band, found by geometric bisection (the classical smoothing spline).
    """
    n = len(values)
    if n < settings.min_curve_points:
        return [float(value) for value in values]
    band = settings.band_fraction * (
        settings.offset_tolerance if tolerance is None else tolerance
    )
    if max(values) - min(values) <= 1.0e-12:
        return [float(value) for value in values]  # a flat curve has nothing to fair
    fit_curve = _local_quadratic(coords, values, settings.local_window)
    sigma = _residual_scale(values, fit_curve, settings.local_window)
    if not (sigma > 0.0):
        return [float(value) for value in values]
    shape_cut = settings.shape_threshold * sigma
    if settings.kernel == "local":
        candidate = fit_curve
    else:
        candidate, band = _energy_kernel(values, settings, band, n)
    out: list[float] = []
    for i, value in enumerate(values):
        moved = candidate[i] - value
        if abs(fit_curve[i] - value) > shape_cut or moved == 0.0:
            out.append(float(value))
            continue
        out.append(min(max(candidate[i], value - band), value + band))
    return out


def _energy_kernel(
    values: Sequence[float],
    settings: FairingSettings,
    band: float,
    n: int,
) -> tuple[list[float], float]:
    """Largest third-derivative weight whose solution stays inside the band."""
    panel, rows = _third_difference_operator(n)
    if not rows:
        return [float(value) for value in values], band
    low, high = settings.smoothing_search
    best = _solve_with_band(
        values, rows, panel, 10.0 ** low, band, settings.clamping_iterations
    )
    if _max_deviation(best, values) > band:
        return [float(value) for value in values], band
    for _ in range(settings.search_steps):
        middle = 0.5 * (low + high)
        candidate = _solve_with_band(
            values, rows, panel, 10.0 ** middle, band, settings.clamping_iterations
        )
        if _max_deviation(candidate, values) <= band:
            low = middle
            best = candidate
        else:
            high = middle
    return best, band


def _residual_scale(
    values: Sequence[float], fitted: Sequence[float], window: int
) -> float:
    """Robust noise scale: ``1.4826 * median |y - local quadratic|`` over interior points.

    The median is what makes this usable on ship data: a section with a sharp turn of the
    bilge puts a large residual on a handful of points, and a mean would call the whole
    curve noisy and smooth the knuckle away.  With the median, a few shape points do not
    move the scale, and the driver can then recognise them (:data:`FairingSettings.shape_threshold`)
    and leave them alone.
    """
    half = window // 2
    residuals = [
        abs(values[i] - fitted[i]) for i in range(half, len(values) - half)
    ]
    if not residuals:
        return 0.0
    residuals.sort()
    median = residuals[len(residuals) // 2]
    scale = 1.4826 * median
    peak = max((abs(value) for value in values), default=0.0)
    return scale if scale > 1.0e-12 * max(peak, 1.0) else 0.0


def _max_deviation(faired: Sequence[float], original: Sequence[float]) -> float:
    return max((abs(a - b) for a, b in zip(faired, original, strict=True)), default=0.0)


# --------------------------------------------------------------------------- #
# audits
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class CurveFairness:
    """What one curve looked like before and after the pass that touched it."""

    family: str
    index: int
    label: str
    points: int
    max_deviation: float
    flips_before: int
    flips_after: int
    energy_before: float
    energy_after: float
    skipped: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "family": self.family,
            "index": self.index,
            "label": self.label,
            "points": self.points,
            "max_deviation": self.max_deviation,
            "flips_before": self.flips_before,
            "flips_after": self.flips_after,
            "bending_before": self.energy_before,
            "bending_after": self.energy_after,
            "skipped": self.skipped,
        }


@dataclass(frozen=True)
class FamilyFairness:
    """Aggregate of :class:`CurveFairness` over one direction."""

    family: str
    curves: tuple[CurveFairness, ...]
    monotone_refusals: int = 0

    @property
    def evaluated(self) -> tuple[CurveFairness, ...]:
        return tuple(curve for curve in self.curves if not curve.skipped)

    @property
    def flips_before(self) -> int:
        return sum(curve.flips_before for curve in self.evaluated)

    @property
    def flips_after(self) -> int:
        return sum(curve.flips_after for curve in self.evaluated)

    @property
    def max_deviation(self) -> float:
        return max((curve.max_deviation for curve in self.evaluated), default=0.0)

    @property
    def energy_before(self) -> float:
        return math.fsum(curve.energy_before for curve in self.evaluated)

    @property
    def energy_after(self) -> float:
        return math.fsum(curve.energy_after for curve in self.evaluated)

    def to_dict(self) -> dict[str, Any]:
        return {
            "family": self.family,
            "curves": len(self.curves),
            "evaluated": len(self.evaluated),
            "skipped": [curve.label for curve in self.curves if curve.skipped],
            "flips_before": self.flips_before,
            "flips_after": self.flips_after,
            "max_deviation": self.max_deviation,
            "bending_before": self.energy_before,
            "bending_after": self.energy_after,
            "monotone_refusals": self.monotone_refusals,
        }


@dataclass(frozen=True)
class FairingPass:
    """One sweep over the three families, with the accept/reject decision recorded."""

    index: int
    relaxation: float
    accepted: bool
    max_update: float
    total_flips: int
    total_energy: float
    families: tuple[FamilyFairness, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "index": self.index,
            "relaxation": self.relaxation,
            "accepted": self.accepted,
            "max_update": self.max_update,
            "total_flips": self.total_flips,
            "total_bending": self.total_energy,
            "families": [family.to_dict() for family in self.families],
        }


@dataclass(frozen=True)
class FairingResult:
    """The faired table plus everything needed to contest the result."""

    input_table: ShipOffsetTable
    table: ShipOffsetTable
    settings: FairingSettings
    passes: tuple[FairingPass, ...]
    max_deviation: float
    converged: bool
    volume_before: float
    volume_after: float
    notes: tuple[str, ...] = field(default_factory=tuple)

    @property
    def flips(self) -> dict[str, tuple[int, ...]]:
        """Flip total per enabled family: pass 0 (input) then every accepted pass."""
        initial = self._initial_flips
        out: dict[str, tuple[int, ...]] = {}
        for family in self.settings.families:
            history = [initial[family]]
            for entry in self.passes:
                if not entry.accepted:
                    continue
                match = [audit for audit in entry.families if audit.family == family]
                if match:
                    history.append(match[0].flips_after)
            out[family] = tuple(history)
        return out

    @property
    def _initial_flips(self) -> dict[str, int]:
        return {
            family: _family_flips(self.input_table, family, self.settings)
            for family in self.settings.families
        }

    @property
    def total_flips(self) -> tuple[int, ...]:
        """Table-level flip total: pass 0 (input) then every accepted pass."""
        history = [sum(self._initial_flips.values())]
        history.extend(fairing.total_flips for fairing in self.passes if fairing.accepted)
        return tuple(history)

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA_FAIRING,
            "table": self.table.name,
            "offset_tolerance": self.settings.offset_tolerance,
            "band_fraction": self.settings.band_fraction,
            "families": list(self.settings.families),
            "max_deviation": self.max_deviation,
            "tolerance_respected": self.max_deviation
            <= self.settings.offset_tolerance * (1.0 + 1e-12),
            "converged": self.converged,
            "flips": {key: list(value) for key, value in self.flips.items()},
            "total_flips": list(self.total_flips),
            "volume_before": self.volume_before,
            "volume_after": self.volume_after,
            "volume_change": self.volume_after - self.volume_before,
            "passes": [fairing.to_dict() for fairing in self.passes],
            "notes": list(self.notes),
        }


# --------------------------------------------------------------------------- #
# curve enumeration for the three families
# --------------------------------------------------------------------------- #
def _sections(table: ShipOffsetTable) -> list[tuple[float, ...]]:
    return [tuple(row) for row in table.half_breadths]


def _family_flips(table: ShipOffsetTable, family: str, settings: FairingSettings) -> int:
    if family == "waterline":
        grid = _sections(table)
        total = 0
        for j in range(len(table.waterlines)):
            values = [grid[i][j] for i in range(len(table.stations))]
            total += curvature_sign_flips(
                table.stations, values, floor_ratio=settings.curvature_floor
            )
        return total
    if family == "station":
        total = 0
        for row in table.half_breadths:
            total += curvature_sign_flips(
                table.waterlines, row, floor_ratio=settings.curvature_floor
            )
        return total
    total = 0
    for level, coords, depths in _buttock_curves(table, settings):
        del level
        total += curvature_sign_flips(coords, depths, floor_ratio=settings.curvature_floor)
    return total


def _section_crossing(
    z_levels: Sequence[float], row: Sequence[float], target: float
) -> tuple[float, int, float] | None:
    """Depth where a section reaches ``target`` half-breadth, the bracketing lower index,
    and the local ``dy/dz``.  ``None`` when the section never reaches it."""
    if target <= 0.0 or target > max(row):
        return None
    for j in range(len(row) - 1):
        low, high = row[j], row[j + 1]
        if low <= target <= high and high > low:
            fraction = (target - low) / (high - low)
            depth = z_levels[j] + fraction * (z_levels[j + 1] - z_levels[j])
            slope = (high - low) / (z_levels[j + 1] - z_levels[j])
            return depth, j, slope
    return None


def _buttock_curves(
    table: ShipOffsetTable, settings: FairingSettings
) -> list[tuple[float, tuple[float, ...], tuple[float, ...]]]:
    """``(level, stations, depths)`` of the iso-half-breadth contours that exist.

    The levels are ``k / (levels + 1)`` of the widest ordinate in the table, so they are
    hull-independent and always inside the surface.  Only contours reaching
    ``min_curve_points`` stations are returned: a buttock that exists over three stations
    is not a curve and fairing it would be fitting noise.
    """
    peak = max((max(row) for row in table.half_breadths), default=0.0)
    if peak <= 0.0:
        return []
    out: list[tuple[float, tuple[float, ...], tuple[float, ...]]] = []
    for k in range(1, settings.buttock_levels + 1):
        target = peak * k / (settings.buttock_levels + 1)
        stations: list[float] = []
        depths: list[float] = []
        for x, row in zip(table.stations, table.half_breadths, strict=True):
            crossing = _section_crossing(table.waterlines, row, target)
            if crossing is None:
                continue
            stations.append(x)
            depths.append(crossing[0])
        if len(stations) >= settings.min_curve_points:
            out.append((target, tuple(stations), tuple(depths)))
    return out


# --------------------------------------------------------------------------- #
# the three passes
# --------------------------------------------------------------------------- #
def _audit_curve(
    family: str,
    index: int,
    label: str,
    coords: Sequence[float],
    before: Sequence[float],
    after: Sequence[float],
    settings: FairingSettings,
    skipped: str = "",
    deviation: float | None = None,
) -> CurveFairness:
    """Summarise one curve.

    ``deviation`` overrides the measured ``|after - before|`` when the faired quantity is
    not the ordinate (buttock: the curve is a depth, so the audit reports the ordinate
    movement the correction implies, which is the quantity the tolerance bounds).
    """
    floor = settings.curvature_floor
    return CurveFairness(
        family=family,
        index=index,
        label=label,
        points=len(coords),
        max_deviation=(
            _max_deviation(after, before) if deviation is None else float(deviation)
        ),
        flips_before=curvature_sign_flips(coords, before, floor_ratio=floor),
        flips_after=curvature_sign_flips(coords, after, floor_ratio=floor),
        energy_before=line_bending_energy(coords, before),
        energy_after=line_bending_energy(coords, after),
        skipped=skipped,
    )


def _pass_waterlines(
    grid: list[list[float]], table: ShipOffsetTable, settings: FairingSettings, relax: float
) -> tuple[list[list[float]], FamilyFairness]:
    new = [list(row) for row in grid]
    curves: list[CurveFairness] = []
    for j in range(len(table.waterlines)):
        before = [grid[i][j] for i in range(len(grid))]
        faired = _fair_curve(table.stations, before, settings)
        for i, value in enumerate(faired):
            new[i][j] = grid[i][j] + relax * (value - grid[i][j])
        curves.append(
            _audit_curve(
                "waterline", j, f"z={table.waterlines[j]:.3f}", table.stations, before,
                [new[i][j] for i in range(len(new))], settings,
            )
        )
    return new, FamilyFairness("waterline", tuple(curves))


def _pass_stations(
    grid: list[list[float]], table: ShipOffsetTable, settings: FairingSettings, relax: float
) -> tuple[list[list[float]], FamilyFairness]:
    new = [list(row) for row in grid]
    curves: list[CurveFairness] = []
    for i in range(len(grid)):
        before = list(grid[i])
        faired = _fair_curve(table.waterlines, before, settings)
        new[i] = [before[k] + relax * (faired[k] - before[k]) for k in range(len(before))]
        curves.append(
            _audit_curve(
                "station", i, f"x={table.stations[i]:.3f}", table.waterlines, before,
                new[i], settings,
            )
        )
    return new, FamilyFairness("station", tuple(curves))


def _pass_buttocks(
    grid: list[list[float]], table: ShipOffsetTable, settings: FairingSettings, relax: float
) -> tuple[list[list[float]], FamilyFairness]:
    """Fair the iso-``y`` contours and push the depth correction back onto ordinates.

    For a contour at level ``y_k`` the crossing depth at station ``i`` moves by ``dz``;
    along the section the ordinate at the bracketing waterline must move by
    ``-slope * dz / w`` (``w`` the linear interpolation weight of that node at the
    crossing) for the crossing to land there, and the per-curve band is the tolerance
    divided by the largest section slope, so ``|dy| <= tolerance`` holds by construction.
    A correction that would drive an ordinate negative, and any section that would stop
    being monotone in ``z``, is refused and reported (``monotone_refusals``) instead of
    being applied silently.
    """
    new = [list(row) for row in grid]
    curves: list[CurveFairness] = []
    refusals = 0
    for index, (level, stations, depths) in enumerate(_buttock_curves(table, settings)):
        moved = 0.0
        brackets: list[tuple[int, int, float, float]] = []  # station, node, weight, slope
        for x in stations:
            i = min(range(len(table.stations)), key=lambda k: abs(table.stations[k] - x))
            crossing = _section_crossing(table.waterlines, grid[i], level)
            if crossing is None:
                brackets.append((-1, -1, 0.0, float("inf")))
                continue
            _, j, slope = crossing
            span = max(table.waterlines[j + 1] - table.waterlines[j], 1.0e-12)
            weight = 1.0 - abs(crossing[0] - table.waterlines[j]) / span
            if weight >= 0.5:
                brackets.append((i, j, weight, slope))
            else:
                brackets.append((i, j + 1, 1.0 - weight, slope))
        slopes = [entry[3] for entry in brackets if math.isfinite(entry[3]) and entry[3] > 0.0]
        if len(slopes) < settings.min_curve_points:
            continue
        budget = settings.buttock_band_fraction * settings.offset_tolerance
        depth_tolerance = budget / max(slopes)
        faired = _fair_curve(stations, depths, settings, tolerance=depth_tolerance)
        for position, (i, node, weight, slope) in enumerate(brackets):
            if i < 0 or node < 0 or not math.isfinite(slope) or slope <= 0.0:
                continue
            dz = relax * (faired[position] - depths[position])
            if dz == 0.0:
                continue
            delta = -slope * dz / max(weight, 1.0e-3)
            delta = min(max(delta, -budget), budget)
            if new[i][node] + delta < 0.0:
                refusals += 1
                continue
            new[i][node] += delta
            moved = max(moved, abs(delta))
        curves.append(
            _audit_curve(
                "buttock",
                index,
                f"y={level:.3f}",
                stations,
                depths,
                [depths[p] + relax * (faired[p] - depths[p]) for p in range(len(depths))],
                settings,
                deviation=moved,
            )
        )
    for row in new:
        for k in range(1, len(row)):
            if row[k] < row[k - 1]:
                row[k] = row[k - 1]
                refusals += 1
    return new, FamilyFairness("buttock", tuple(curves), monotone_refusals=refusals)


# --------------------------------------------------------------------------- #
# driver
# --------------------------------------------------------------------------- #
def _snapshot(table: ShipOffsetTable, grid: Sequence[Sequence[float]]) -> ShipOffsetTable:
    return ShipOffsetTable(
        name=f"{table.name}-faired",
        length=table.length,
        stations=table.stations,
        waterlines=table.waterlines,
        half_breadths=tuple(tuple(float(value) for value in row) for row in grid),
        design_draft=table.design_draft,
        provenance=OffsetProvenance(
            source=f"{table.name} (faired in place by tensorcad.lines_fairing)",
            method=(
                "three-directional curvature-variation fairing inside an ordinate band; "
                "see tensorcad/lines_fairing.py and docs/lines_fairing.md"
            ),
            url="https://github.com/cxs503/TensorCAD/blob/main/tensorcad/lines_fairing.py",
            licence=table.provenance.licence,
            official=False,
            notes=f"derived from: {table.provenance.source}",
        ),
        units=table.units,
        waterline_exponent=table.waterline_exponent,
    )


def _project(grid: Sequence[Sequence[float]], original: Sequence[Sequence[float]], band: float) -> list[list[float]]:
    # a half-breadth is a length: the band is intersected with y >= 0
    return [
        [
            min(max(value, max(0.0, source - band)), source + band)
            for value, source in zip(row, source_row, strict=True)
        ]
        for row, source_row in zip(grid, original, strict=True)
    ]


def _total(grid: list[list[float]], table: ShipOffsetTable, settings: FairingSettings) -> tuple[int, float]:
    snapshot = _snapshot(table, grid)
    flips = sum(_family_flips(snapshot, family, settings) for family in settings.families)
    energy = 0.0
    for j in range(len(table.waterlines)):
        energy += line_bending_energy(
            table.stations, [grid[i][j] for i in range(len(grid))]
        )
    for row in grid:
        energy += line_bending_energy(table.waterlines, row)
    return flips, energy


def fair_offset_table(
    table: ShipOffsetTable, settings: FairingSettings | None = None
) -> FairingResult:
    """Fair ``table`` in all three directions and report what moved.

    The driver sweeps the enabled families once per pass.  A pass is committed only if the
    lexicographic objective ``(total flips, total bending energy)`` does not get worse;
    otherwise the sweep is retried with half the relaxation (up to ``max_passes``), and if
    no relaxed sweep helps the run stops.  The recorded flip history is therefore
    non-increasing by construction, which is exactly the acceptance clause, and a stalled
    run is reported as ``converged=False`` rather than papered over.  The guarantee is on
    the *total*: one family may gain a reversal while another loses two (the per-family
    histories are reported so that trade is visible, never hidden).
    """
    settings = settings or FairingSettings()
    original = [list(row) for row in table.half_breadths]
    grid = [list(row) for row in original]
    band = settings.offset_tolerance
    passes: list[FairingPass] = []
    notes: list[str] = []
    converged = False
    flips, energy = _total(grid, table, settings)

    for index in range(1, settings.max_passes + 1):
        relaxation = settings.relax
        accepted = False
        attempt = 0
        while not accepted and attempt < 3:
            attempt += 1
            candidate = [list(row) for row in grid]
            families: list[FamilyFairness] = []
            if "waterline" in settings.families:
                candidate, audit = _pass_waterlines(candidate, table, settings, relaxation)
                families.append(audit)
            if "station" in settings.families:
                candidate, audit = _pass_stations(candidate, table, settings, relaxation)
                families.append(audit)
            if "buttock" in settings.families:
                candidate, audit = _pass_buttocks(candidate, table, settings, relaxation)
                families.append(audit)
            candidate = _project(candidate, original, band)
            candidate_flips, candidate_energy = _total(candidate, table, settings)
            max_update = max(
                (
                    abs(candidate[i][j] - grid[i][j])
                    for i in range(len(grid))
                    for j in range(len(grid[i]))
                ),
                default=0.0,
            )
            improve = (candidate_flips, candidate_energy) < (flips, energy)
            if improve:
                passes.append(
                    FairingPass(
                        index=index,
                        relaxation=relaxation,
                        accepted=True,
                        max_update=max_update,
                        total_flips=candidate_flips,
                        total_energy=candidate_energy,
                        families=tuple(families),
                    )
                )
                grid = candidate
                flips, energy = candidate_flips, candidate_energy
                accepted = True
            else:
                passes.append(
                    FairingPass(
                        index=index,
                        relaxation=relaxation,
                        accepted=False,
                        max_update=max_update,
                        total_flips=candidate_flips,
                        total_energy=candidate_energy,
                        families=tuple(families),
                    )
                )
                relaxation *= 0.5
                if max_update <= 0.0:
                    break
        if not accepted:
            notes.append(
                f"pass {index}: no relaxed sweep improved (flips, energy); stopped at "
                f"{flips} flips"
            )
            break
        if max_update == 0.0:
            converged = True
            notes.append(f"pass {index}: table already fair, nothing moved")
            break
    else:
        notes.append(
            f"max_passes={settings.max_passes} reached with flips={flips}; "
            "the tolerance band is the likely limiter"
        )

    faired = _snapshot(table, grid)
    before = table_integrals(table)
    after = table_integrals(faired)
    max_deviation = max(
        (
            abs(grid[i][j] - original[i][j])
            for i in range(len(grid))
            for j in range(len(grid[i]))
        ),
        default=0.0,
    )
    if max_deviation > band * (1.0 + 1e-12):
        raise FairingError(
            f"fairing left the tolerance band: {max_deviation!r} > {band!r}"
        )
    return FairingResult(
        input_table=table,
        table=faired,
        settings=settings,
        passes=tuple(passes),
        max_deviation=max_deviation,
        converged=converged or flips == 0,
        volume_before=before.volume,
        volume_after=after.volume,
        notes=tuple(notes),
    )


# --------------------------------------------------------------------------- #
# an independent reference smoother
# --------------------------------------------------------------------------- #
def hand_fairing_reference(
    table: ShipOffsetTable, *, window: int = 5, passes: int = 1
) -> tuple[tuple[tuple[float, ...], ...], tuple[tuple[float, ...], ...]]:
    """Local-quadratic (Savitzky-Golay) smoothing of the same table, in the same directions.

    This exists so the energy-minimising fairing can be compared against an algorithm that
    shares none of its machinery: a moving least-squares quadratic fit to a
    ``window``-point stencil, applied to the waterlines and the stations, with the edges
    left alone.  It is *not* a human designer's hand-faired table -- no marked-up lines
    plan is available here -- but it is an independent estimate of "where the fair
    ordinates are", and the agreement between the two is reported.

    Returns ``(faired_rows, weights)``: the smoothed ordinates and ``1.0`` for a point that
    the smoother actually touched, ``0.0`` for an edge point it left as it found it.
    """
    if window < 3 or window % 2 == 0:
        raise FairingError(f"window must be odd and >= 3, got {window!r}")
    half = window // 2
    grid = [list(row) for row in table.half_breadths]
    touched = [[0.0] * len(grid[0]) for _ in grid]

    def smooth(values: Sequence[float], coords: Sequence[float]) -> list[float]:
        n = len(values)
        out = list(values)
        for i in range(half, n - half):
            xs = coords[i - half : i + half + 1]
            ys = values[i - half : i + half + 1]
            # least squares quadratic through the stencil, evaluated at x_i
            moment = [[sum((x - xs[half]) ** (p + q) for x in xs) for q in range(3)] for p in range(3)]
            rhs = [sum((x - xs[half]) ** p * y for x, y in zip(xs, ys, strict=True)) for p in range(3)]
            coefficients = _solve_dense(moment, rhs)
            out[i] = coefficients[0]
        return out

    for _ in range(passes):
        for j in range(len(table.waterlines)):
            column = [grid[i][j] for i in range(len(grid))]
            smoothed = smooth(column, table.stations)
            for i, value in enumerate(smoothed):
                grid[i][j] = value
                if half <= i < len(grid) - half:
                    touched[i][j] = 1.0
        for i, row in enumerate(grid):
            smoothed = smooth(row, table.waterlines)
            grid[i] = smoothed
            for j in range(len(row)):
                if half <= j < len(row) - half:
                    touched[i][j] = 1.0
    return (
        tuple(tuple(value for value in row) for row in grid),
        tuple(tuple(flag for flag in row) for row in touched),
    )


def _solve_dense(matrix: Sequence[Sequence[float]], rhs: Sequence[float]) -> list[float]:
    """Gaussian elimination with partial pivoting (3x3 systems here)."""
    n = len(rhs)
    a = [list(row) + [float(value)] for row, value in zip(matrix, rhs, strict=True)]
    for column in range(n):
        pivot = max(range(column, n), key=lambda r: abs(a[r][column]))
        if abs(a[pivot][column]) <= 1.0e-18:
            raise FairingError("singular stencil in hand_fairing_reference")
        a[column], a[pivot] = a[pivot], a[column]
        for r in range(column + 1, n):
            factor = a[r][column] / a[column][column]
            for c in range(column, n + 1):
                a[r][c] -= factor * a[column][c]
    solution = [0.0] * n
    for r in range(n - 1, -1, -1):
        total = a[r][n] - sum(a[r][c] * solution[c] for c in range(r + 1, n))
        solution[r] = total / a[r][r]
    return solution