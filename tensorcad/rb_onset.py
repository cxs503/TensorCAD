"""Rayleigh-Benard onset: linear-stability reference value and Ra-bisection harness.

Why this module exists
----------------------
Case-pool entry ``LBM-RB-ONSET-RA1708`` asks for the *critical Rayleigh number*
at which a Newtonian layer heated from below, confined between two infinite
parallel plates, loses its conduction state.  The headline reference value is
Chandrasekhar's rigid-rigid result

    Ra_c = 1707.762   at   a_c = 3.11632                          (exp 1)

and the stress-free companion, which is *exact* for a sine mode,

    Ra_c = 27*pi^4/4 = 657.5114   at   a_c = pi/sqrt(2) = 2.2214  (exp 2)

TensorCAD's role in the ecosystem is the geometry / reference-data side of a
benchmark -- exactly what it does for the lid-driven cavity
(:mod:`tensorcad.cavity_mesh`) and for Kovasznay flow (:mod:`tensorcad.kovasznay`);
the flow solver (LBM / FVM / DG) lives in the solver repositories.  So this
module supplies

1. an **independent, dependency-free** reference computation of ``Ra_c`` and the
   critical wavenumber, by the criterion the case-pool acceptance names --
   *Ra bisection scan driven by the stationary linear-stability criterion*
   (:func:`critical_rayleigh`, :func:`onset`, :func:`neutral_curve`);
2. a **solver-agnostic Ra-bisection harness** (:func:`bisect_onset`) that takes a
   growth-rate probe supplied by a solver and returns the measured onset, plus a
   log-linear growth-rate estimator (:func:`growth_rate`) for the probe side;
3. the **cell specification** a 2-D solver consumes (:func:`rb_cell`): aspect
   ratio, boundary-condition tags, mesh family, probe locations and the monitor
   contract (perturbation energy, Nusselt number);
4. the **critical eigenfunction** (temperature / stream-function depth profiles)
   so a solver's mode shape can be compared, not just its eigenvalue
   (:func:`onset_mode`, :func:`chebyshev_interpolate`).

Theory in one paragraph
-----------------------
Linearise the Boussinesq equations about the conduction state ``T = 1 - z``
(length ``d``, temperature drop ``dT``, bottom plate at ``z = 0``), take normal
modes ``~ exp(i a x)`` and scale temperature so the buoyancy term reads ``Ra``;
the stationary (marginal, exchange-of-stability) problem is

    (D^2 - a^2)^2 Psi - Ra Theta = 0                                (1a)
    (D^2 - a^2) Theta + a^2 Psi  = 0                                (1b)

with ``D = d/dz``.  Eliminating ``Theta`` recovers the classical sixth-order form
``(D^2 - a^2)^3 Psi = -Ra a^2 Psi`` (Chandrasekhar 1961, ch. II).  Boundary
conditions on impermeable isothermal plates: ``Theta = 0`` everywhere on the
plates together with ``Psi = D Psi = 0`` for **rigid** (no-slip) plates or
``Psi = D^2 Psi = 0`` for **free** (stress-free) plates.  The conduction state is
linearly stable iff (1) has only the trivial solution, i.e. iff
``det M(a, Ra) != 0``; the neutral curve ``Ra(a)`` is the smallest positive
``Ra`` with ``det M = 0`` and ``Ra_c = min_a Ra(a)``.  ``Ra`` enters (1a)
linearly, so ``det M`` is a polynomial in ``Ra`` and the first sign change of
``det M`` while scanning ``Ra`` upward is the marginal state of the lowest mode.

Discretisation
--------------
Chebyshev-Lobatto collocation in ``z in [0, 1]`` on the coupled system (1a)-(1b);
the four ``Psi`` conditions and the two ``Theta`` conditions replace the operator
rows adjacent to each wall ("boundary bordering").  Verified against (exp 2):
``Ra_c(a_c) = 657.511364479`` (relative error ``3e-13`` at ``order = 16``);
against (exp 1) the rigid-rigid value converges to ``1707.76178`` (relative
error ``1.3e-7``; see the table in ``docs/rb_onset.md``).  The determinant is
accumulated in log space because ``det`` of the bordered matrix (entries up to
``order^8``) overflows double precision -- see :func:`determinant_sign_log`.

Reference
---------
Chandrasekhar, S. (1961). *Hydrodynamic and Hydromagnetic Stability*, Oxford
University Press, ch. II (Benard problem), table III (rigid-rigid
``Ra_c = 1707.762``).  Drazin & Reid (1981), *Hydrodynamic Stability*, ch. II.
de Vahl Davis (1983), *Int. J. Numer. Methods Fluids* 3, 249 (``Nu`` benchmarks
used by the supercritical check ``Ra = 1e4, 1e5``).
"""
from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

__all__ = [
    "CASE_SPEC",
    "CHANDRASEKHAR_1961",
    "FREE_FREE_ANALYTIC",
    "NONDIMENSIONALISATION",
    "RB_ONSET_ACCEPTANCE",
    "OnsetMode",
    "OnsetScan",
    "OnsetVerdict",
    "RbCell",
    "bisect_onset",
    "chebyshev_differentiation_matrix",
    "chebyshev_interpolate",
    "chebyshev_lobatto",
    "critical_rayleigh",
    "determinant_sign_log",
    "growth_rate",
    "linear_stability_matrix",
    "neutral_curve",
    "nusselt_number",
    "onset",
    "onset_mode",
    "rb_cell",
    "write_bundle",
]

# --------------------------------------------------------------------------- #
# Reference values
# --------------------------------------------------------------------------- #

#: Chandrasekhar (1961) table III: rigid (no-slip), isothermal plates.
CHANDRASEKHAR_1961: dict[str, Any] = {
    "walls": "rigid-rigid",
    "rayleigh_critical": 1707.762,
    "wavenumber_critical": 3.11632,
    "source": "Chandrasekhar 1961, Hydrodynamic and Hydromagnetic Stability, ch. II, table III",
    "doi": "10.1016/0368-3281(63)90025-1",
}

#: Exact closed form for stress-free, isothermal plates (sine mode, n = 1).
FREE_FREE_ANALYTIC: dict[str, Any] = {
    "walls": "free-free",
    "rayleigh_critical": 27.0 * math.pi**4 / 4.0,  # 657.5113644795163
    "wavenumber_critical": math.pi / math.sqrt(2.0),  # 2.221441469079183
    "source": "Ra = (n^2 pi^2 + a^2)^3 / a^2 minimised at a^2 = n^2 pi^2 / 2 -> 27 n^4 pi^4 / 4",
}

#: Acceptance of case-pool entry ``LBM-RB-ONSET-RA1708`` (tolerance on Ra_c).
RB_ONSET_ACCEPTANCE: dict[str, Any] = {
    "case": "LBM-RB-ONSET-RA1708",
    "metric": "critical Rayleigh number judged by an Ra bisection scan",
    "reference": 1707.762,
    "reference_source": "Chandrasekhar 1961 ch. II (rigid-rigid)",
    "tolerance_relative": 0.02,
    "wavenumber_reference": 3.11632,
    "wavenumber_tolerance_relative": 0.02,
}

#: Nondimensionalisation of the cell (documented once, referenced everywhere).
NONDIMENSIONALISATION: dict[str, str] = {
    "length": "d (plate separation)",
    "temperature": "dT = T_hot - T_cold = 1 (bottom plate hot, z = 0)",
    "velocity": "kappa / d",
    "time": "d^2 / kappa",
    "rayleigh": "Ra = g beta dT d^3 / (kappa nu)",
    "note": "the stationary onset is independent of the Prandtl number",
}

_FREE = "free"
_RIGID = "rigid"
_WALLS = {"rigid": _RIGID, "free": _FREE, "rigid-rigid": _RIGID, "free-free": _FREE}

#: Machine-readable benchmark specification (mirrors the case-pool ledger entry).
CASE_SPEC: dict[str, Any] = {
    "case": RB_ONSET_ACCEPTANCE["case"],
    "title": "Rayleigh-Benard onset: critical Rayleigh number (rigid-rigid plates)",
    "reference": CHANDRASEKHAR_1961["rayleigh_critical"],
    "reference_wavenumber": CHANDRASEKHAR_1961["wavenumber_critical"],
    "reference_source": CHANDRASEKHAR_1961["source"],
    "doi_secondary": CHANDRASEKHAR_1961["doi"],
    "tolerance_relative": RB_ONSET_ACCEPTANCE["tolerance_relative"],
    "method": "Ra bisection scan with the stationary linear-stability criterion",
    "nondimensionalisation": NONDIMENSIONALISATION,
    "analytic_sanity_check": {
        "walls": "free-free",
        "rayleigh_critical": FREE_FREE_ANALYTIC["rayleigh_critical"],
        "wavenumber_critical": FREE_FREE_ANALYTIC["wavenumber_critical"],
        "source": FREE_FREE_ANALYTIC["source"],
    },
}


def _normalise_walls(walls: str) -> str:
    try:
        return _WALLS[str(walls).lower()]
    except KeyError:
        raise ValueError(f"unknown walls {walls!r}; expected 'rigid' or 'free'") from None


# --------------------------------------------------------------------------- #
# Chebyshev machinery (standard library only)
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class _Chebyshev:
    """Chebyshev-Lobatto nodes on ``[0, 1]`` and derivative operators."""

    order: int
    nodes: tuple[float, ...]
    d1: tuple[tuple[float, ...], ...]
    d2: tuple[tuple[float, ...], ...]


def chebyshev_lobatto(order: int) -> list[float]:
    """Chebyshev-Lobatto nodes ``0.5 * (1 - cos(pi j / order))`` in ``[0, 1]``.

    Ordered from the hot plate (``z = 0``) upward to the cold plate (``z = 1``).
    """
    if order < 1:
        raise ValueError("order must be >= 1")
    return [0.5 * (1.0 - math.cos(math.pi * j / order)) for j in range(order + 1)]


def chebyshev_differentiation_matrix(order: int) -> list[list[float]]:
    """First-derivative matrix ``D_z`` for :func:`chebyshev_lobatto` nodes on [0, 1].

    Trefethen (2000) *Spectral Methods in MATLAB*, ``cheb.m``, chained to ``[0, 1]``
    through ``z = (1 - xi) / 2``, i.e. ``d/dz = -2 d/dxi``.
    """
    if order < 1:
        raise ValueError("order must be >= 1")
    n = order
    xi = [math.cos(math.pi * j / n) for j in range(n + 1)]
    c = [2.0 if j in (0, n) else 1.0 for j in range(n + 1)]
    d = [[0.0] * (n + 1) for _ in range(n + 1)]
    for i in range(n + 1):
        for j in range(n + 1):
            if i != j:
                d[i][j] = (c[i] / c[j]) * ((-1.0) ** (i + j)) / (xi[i] - xi[j])
    for i in range(1, n):
        d[i][i] = -xi[i] / (2.0 * (1.0 - xi[i] ** 2))
    d[0][0] = (2.0 * n * n + 1.0) / 6.0
    d[n][n] = -(2.0 * n * n + 1.0) / 6.0
    # z = (1 - xi) / 2, so d/dz = -2 d/dxi: the node ordering is reversed by the map
    return [[-2.0 * v for v in row] for row in d]


def _matmul(a: Sequence[Sequence[float]], b: Sequence[Sequence[float]]) -> list[list[float]]:
    rows, inner, cols = len(a), len(b), len(b[0])
    return [[sum(a[i][t] * b[t][j] for t in range(inner)) for j in range(cols)] for i in range(rows)]


_CHEB_CACHE: dict[int, _Chebyshev] = {}


def _chebyshev(order: int) -> _Chebyshev:
    cached = _CHEB_CACHE.get(order)
    if cached is not None:
        return cached
    d1 = chebyshev_differentiation_matrix(order)
    d2 = _matmul(d1, d1)
    table = _Chebyshev(
        order=order,
        nodes=tuple(chebyshev_lobatto(order)),
        d1=tuple(tuple(r) for r in d1),
        d2=tuple(tuple(r) for r in d2),
    )
    _CHEB_CACHE[order] = table
    return table


def chebyshev_interpolate(nodes: Sequence[float], values: Sequence[float], z: float) -> float:
    """Barycentric interpolation of ``values`` sampled on Chebyshev-Lobatto nodes.

    Second barycentric formula with weights ``(-1)^j * (1/2 at the ends)``
    (Berrut & Trefethen 2004).  Exact for polynomials of degree ``<= len(nodes)-1``.
    """
    n = len(nodes) - 1
    if len(values) != len(nodes):
        raise ValueError("nodes and values must have the same length")
    w = [(-1.0) ** j * (0.5 if j in (0, n) else 1.0) for j in range(n + 1)]
    num = 0.0
    den = 0.0
    for j in range(n + 1):
        diff = z - nodes[j]
        if abs(diff) < 1e-15:
            return float(values[j])
        term = w[j] / diff
        num += term * values[j]
        den += term
    return num / den


# --------------------------------------------------------------------------- #
# Linear stability system
# --------------------------------------------------------------------------- #


def linear_stability_matrix(
    wavenumber: float, rayleigh: float, order: int = 32, walls: str = "rigid"
) -> list[list[float]]:
    """Bordered collocation matrix ``M(a, Ra)`` of the stationary problem (1).

    Unknowns are ordered ``[Psi_0 .. Psi_order, Theta_0 .. Theta_order]``; a
    non-trivial stationary perturbation exists iff ``det M = 0``.
    """
    wall = _normalise_walls(walls)
    if order < 6:
        raise ValueError("order must be >= 6")
    if wavenumber <= 0.0:
        raise ValueError("wavenumber must be positive")
    table = _chebyshev(order)
    n = order
    size = 2 * (n + 1)
    a2 = wavenumber * wavenumber
    m = [[0.0] * size for _ in range(size)]
    # L = D2 - a^2 I: applied to Psi in rows 0..n, to Theta in rows n+1..2n+1
    l_op = [[table.d2[i][j] - (a2 if i == j else 0.0) for j in range(n + 1)] for i in range(n + 1)]
    l2 = _matmul(l_op, l_op)
    for i in range(n + 1):
        for j in range(n + 1):
            m[i][j] += l2[i][j]  # (1a): (D^2 - a^2)^2 Psi
            m[n + 1 + i][j] += a2 if i == j else 0.0  # (1b): + a^2 Psi
            m[n + 1 + i][n + 1 + j] += l_op[i][j]  # (1b): (D^2 - a^2) Theta
        m[i][n + 1 + i] += -rayleigh  # (1a): - Ra Theta

    def clear(row: int) -> None:
        for j in range(size):
            m[row][j] = 0.0

    # Psi = 0 at both plates (replaces the (1a) rows sitting on the walls)
    clear(0)
    m[0][0] = 1.0
    clear(n)
    m[n][n] = 1.0
    # kinematic condition on the plates: D Psi (rigid) or D^2 Psi (free)
    wall_op = table.d1 if wall == _RIGID else table.d2
    for row, src in ((1, 0), (n - 1, n)):
        clear(row)
        for j in range(n + 1):
            m[row][j] = wall_op[src][j]
    # Theta = 0 at both plates (replaces the (1b) rows sitting on the walls)
    clear(n + 1)
    m[n + 1][n + 1] = 1.0
    clear(2 * n + 1)
    m[2 * n + 1][2 * n + 1] = 1.0
    return m


def determinant_sign_log(matrix: Sequence[Sequence[float]]) -> tuple[int, float]:
    """``(sign, log|det|)`` by LU with partial pivoting.

    The bordered matrix carries fourth-order Chebyshev entries (``~ order^8``),
    so ``det`` overflows double precision for ``order >= 20``; the determinant is
    accumulated as ``sum log|u_ii|`` with the sign tracked separately.  This is
    the only quantity the stability criterion needs.
    """
    size = len(matrix)
    a = [list(row) for row in matrix]
    sign = 1
    log_abs = 0.0
    for k in range(size):
        pivot = max(range(k, size), key=lambda i: abs(a[i][k]))
        value = a[pivot][k]
        if value == 0.0:
            return 0, -math.inf
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        log_abs += math.log(abs(value))
        if value < 0.0:
            sign = -sign
        inv = 1.0 / value
        row_k = a[k]
        for i in range(k + 1, size):
            factor = a[i][k] * inv
            if factor:
                row_i = a[i]
                for j in range(k + 1, size):
                    row_i[j] -= factor * row_k[j]
    return sign, log_abs


def _indicator(wavenumber: float, rayleigh: float, order: int, wall: str) -> int:
    sign, _ = determinant_sign_log(linear_stability_matrix(wavenumber, rayleigh, order, wall))
    return sign


# --------------------------------------------------------------------------- #
# Ra bisection scan (the case-pool criterion)
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class OnsetScan:
    """Result of a critical-Rayleigh-number search."""

    walls: str
    order: int
    rayleigh_critical: float
    wavenumber_critical: float
    wavenumber_scan: tuple[tuple[float, float], ...]
    reference: float
    reference_source: str
    tolerance_relative: float
    verdict: str

    @property
    def relative_error(self) -> float:
        """Signed relative error of ``Ra_c`` against the reference value."""
        return (self.rayleigh_critical - self.reference) / self.reference

    @property
    def accepted(self) -> bool:
        return abs(self.relative_error) <= self.tolerance_relative

    def as_dict(self) -> dict[str, Any]:
        return {
            "walls": self.walls,
            "order": self.order,
            "rayleigh_critical": self.rayleigh_critical,
            "wavenumber_critical": self.wavenumber_critical,
            "relative_error": self.relative_error,
            "reference": self.reference,
            "reference_source": self.reference_source,
            "tolerance_relative": self.tolerance_relative,
            "accepted": self.accepted,
            "verdict": self.verdict,
            "wavenumber_scan": [[a, ra] for a, ra in self.wavenumber_scan],
        }


def critical_rayleigh(
    wavenumber: float,
    *,
    order: int = 32,
    walls: str = "rigid",
    rayleigh_max: float = 2.0e5,
    rayleigh_min: float = 1.0,
    growth: float = 1.25,
    tol: float = 1e-6,
    max_iterations: int = 80,
) -> float | None:
    """Marginal Rayleigh number for one wavenumber, by bisection on ``det M``.

    Scanning ``Ra`` upward from ``rayleigh_min``, the first sign change of
    ``det M(a, Ra)`` is the marginal state: below it the conduction profile admits
    no non-trivial *stationary* perturbation, above it one appears
    (exchange of stability).  Returns ``None`` when no sign change occurs below
    ``rayleigh_max`` (the wavenumber is then stable in the scanned range).
    """
    wall = _normalise_walls(walls)
    lo = rayleigh_min
    sign_lo = _indicator(wavenumber, lo, order, wall)
    if sign_lo == 0:
        return lo
    ra = lo
    sign_prev = sign_lo
    hi = None
    while True:
        ra *= growth
        if ra > rayleigh_max:
            return None
        sign_new = _indicator(wavenumber, ra, order, wall)
        if sign_new != sign_prev:
            hi = ra
            break
        sign_prev = sign_new
    for _ in range(max_iterations):
        if (hi - lo) <= tol * 0.5 * (hi + lo):
            break
        mid = 0.5 * (lo + hi)
        if _indicator(wavenumber, mid, order, wall) == sign_prev:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def neutral_curve(
    wavenumbers: Iterable[float],
    *,
    order: int = 32,
    walls: str = "rigid",
    rayleigh_max: float = 2.0e5,
    tol: float = 1e-6,
) -> list[tuple[float, float]]:
    """Marginal ``Ra`` for each wavenumber (wavenumbers without a root are dropped)."""
    curve: list[tuple[float, float]] = []
    for a in wavenumbers:
        value = critical_rayleigh(a, order=order, walls=walls, rayleigh_max=rayleigh_max, tol=tol)
        if value is not None:
            curve.append((float(a), value))
    return curve


def _parabolic_vertex(points: Sequence[tuple[float, float]]) -> float | None:
    """Vertex of the parabola through three ``(a, Ra)`` samples near a minimum."""
    if len(points) != 3:
        return None
    (a0, f0), (a1, f1), (a2, f2) = points
    if a0 in (a1, a2) or a1 == a2:
        return None
    slope1 = (f1 - f0) / (a1 - a0)
    slope2 = (f2 - f1) / (a2 - a1)
    quad = (slope2 - slope1) / (a2 - a0)
    if quad == 0.0:
        return None
    linear = slope1 - quad * (a0 + a1)
    vertex = -linear / (2.0 * quad)
    lo, hi = min(a0, a2), max(a0, a2)
    if not (lo < vertex < hi):
        return None
    return vertex


def _golden_minimum(f: Callable[[float], float], lo: float, hi: float, iters: int = 26) -> float:
    inv_phi = (math.sqrt(5.0) - 1.0) / 2.0
    c = hi - inv_phi * (hi - lo)
    d = lo + inv_phi * (hi - lo)
    fc, fd = f(c), f(d)
    for _ in range(iters):
        if fc < fd:
            hi, d, fd = d, c, fc
            c = hi - inv_phi * (hi - lo)
            fc = f(c)
        else:
            lo, c, fc = c, d, fd
            d = lo + inv_phi * (hi - lo)
            fd = f(d)
    return 0.5 * (lo + hi)


def onset(
    *,
    order: int = 32,
    walls: str = "rigid",
    wavenumber_range: tuple[float, float] = (1.0, 5.0),
    wavenumber_points: int = 17,
    rayleigh_max: float = 2.0e5,
    tol: float = 1e-6,
    refine_iterations: int = 16,
    reference: float | None = None,
) -> OnsetScan:
    """Critical Rayleigh number and wavenumber of the conduction state.

    A wavenumber scan locates the minimiser of the neutral curve, then a golden
    section refine on ``a`` finishes the job (the marginal curve is unimodal
    around its minimum for both plate pairs).
    """
    wall = _normalise_walls(walls)
    lo, hi = wavenumber_range
    if not (0.0 < lo < hi):
        raise ValueError("wavenumber_range must satisfy 0 < lo < hi")
    if wavenumber_points < 3:
        raise ValueError("wavenumber_points must be >= 3")
    grid = [lo + (hi - lo) * i / (wavenumber_points - 1) for i in range(wavenumber_points)]
    curve = neutral_curve(grid, order=order, walls=wall, rayleigh_max=rayleigh_max, tol=tol)
    if not curve:
        raise ValueError("no marginal state found: raise rayleigh_max or widen wavenumber_range")
    best = min(range(len(curve)), key=lambda i: curve[i][1])
    a_best = curve[best][0]
    span = (hi - lo) / (wavenumber_points - 1)
    a_lo = max(lo, a_best - span)
    a_hi = min(hi, a_best + span)

    def marginal(a: float) -> float:
        value = critical_rayleigh(a, order=order, walls=wall, rayleigh_max=rayleigh_max, tol=tol)
        return value if value is not None else float("inf")

    # The marginal curve is flat in ``a`` near its minimum, so a vertex from the
    # three lowest grid points ("parabola through samples") pins ``a_c`` far
    # better than bisection-limited objective values alone.
    candidates: list[tuple[float, float]] = [(a_best, curve[best][1])]
    if 0 < best < len(curve) - 1:
        vertex = _parabolic_vertex(curve[best - 1 : best + 2])
        if vertex is not None:
            value = marginal(vertex)
            if math.isfinite(value):
                candidates.append((vertex, value))
                left = max((p for p in curve if p[0] < vertex), default=None, key=lambda p: p[0])
                right = min((p for p in curve if p[0] > vertex), default=None, key=lambda p: p[0])
                if left is not None and right is not None:
                    vertex2 = _parabolic_vertex([left, (vertex, value), right])
                    if vertex2 is not None:
                        value2 = marginal(vertex2)
                        if math.isfinite(value2):
                            candidates.append((vertex2, value2))
    if a_hi > a_lo:
        golden = _golden_minimum(marginal, a_lo, a_hi, iters=refine_iterations)
        value = marginal(golden)
        if math.isfinite(value):
            candidates.append((golden, value))
    # the curve is flat in ``a``: only trust a candidate that improves Ra_c by more
    # than the bisection noise floor
    a_refined, ra_refined = min(candidates, key=lambda pair: pair[1])
    best_a, best_ra = candidates[0]
    if abs(ra_refined - best_ra) <= tol * best_ra:
        a_refined, ra_refined = best_a, best_ra

    if reference is None:
        reference = (
            CHANDRASEKHAR_1961["rayleigh_critical"]
            if wall == _RIGID
            else FREE_FREE_ANALYTIC["rayleigh_critical"]
        )
    source = CHANDRASEKHAR_1961["source"] if wall == _RIGID else FREE_FREE_ANALYTIC["source"]
    tol_rel = float(RB_ONSET_ACCEPTANCE["tolerance_relative"])
    rel = (ra_refined - float(reference)) / float(reference)
    return OnsetScan(
        walls=wall,
        order=order,
        rayleigh_critical=ra_refined,
        wavenumber_critical=a_refined,
        wavenumber_scan=tuple(curve),
        reference=float(reference),
        reference_source=source,
        tolerance_relative=tol_rel,
        verdict="PASS" if abs(rel) <= tol_rel else "FAIL",
    )


# --------------------------------------------------------------------------- #
# Critical eigenfunction (mode shape) -- what a solver compares against
# --------------------------------------------------------------------------- #


def _lu_solve(matrix: Sequence[Sequence[float]], rhs: Sequence[float]) -> list[float]:
    size = len(matrix)
    a = [list(row) for row in matrix]
    b = list(rhs)
    for k in range(size):
        pivot = max(range(k, size), key=lambda i: abs(a[i][k]))
        if a[pivot][k] == 0.0:
            raise ValueError("singular matrix in linear solve")
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            b[k], b[pivot] = b[pivot], b[k]
        inv = 1.0 / a[k][k]
        for i in range(k + 1, size):
            factor = a[i][k] * inv
            if factor:
                row_i = a[i]
                for j in range(k, size):
                    row_i[j] -= factor * a[k][j]
                b[i] -= factor * b[k]
    x = [0.0] * size
    for i in range(size - 1, -1, -1):
        total = b[i] - sum(a[i][j] * x[j] for j in range(i + 1, size))
        x[i] = total / a[i][i]
    return x


def _residual(matrix: Sequence[Sequence[float]], vec: Sequence[float]) -> float:
    """``||M x||_inf / (||M||_inf ||x||_inf)`` -- a dimensionless residual."""
    x_scale = max(abs(v) for v in vec) or 1.0
    m_scale = max(sum(abs(c) for c in row) for row in matrix) or 1.0
    worst = 0.0
    for row in matrix:
        worst = max(worst, abs(sum(c * v for c, v in zip(row, vec, strict=True))))
    return worst / (m_scale * x_scale)


@dataclass(frozen=True)
class OnsetMode:
    """Critical eigenfunction of the stationary problem at ``(a, Ra)``."""

    wavenumber: float
    rayleigh: float
    order: int
    walls: str
    nodes: tuple[float, ...]
    psi: tuple[float, ...]
    theta: tuple[float, ...]
    residual: float

    def stream_function(self, x: float, z: float) -> float:
        """``psi(x, z) = Psi(z) cos(a x)`` -- rolls with axes along ``y``."""
        return chebyshev_interpolate(self.nodes, self.psi, z) * math.cos(self.wavenumber * x)

    def temperature_perturbation(self, x: float, z: float) -> float:
        """``theta(x, z) = Theta(z) cos(a x)``, normalised to unit peak amplitude."""
        return chebyshev_interpolate(self.nodes, self.theta, z) * math.cos(self.wavenumber * x)

    def as_dict(self) -> dict[str, Any]:
        return {
            "wavenumber": self.wavenumber,
            "rayleigh": self.rayleigh,
            "order": self.order,
            "walls": self.walls,
            "z": list(self.nodes),
            "psi": list(self.psi),
            "theta": list(self.theta),
            "residual": self.residual,
        }


def onset_mode(
    wavenumber: float,
    rayleigh: float,
    *,
    order: int = 32,
    walls: str = "rigid",
    iterations: int = 40,
) -> OnsetMode:
    """Null vector of the bordered matrix by inverse iteration.

    At the marginal state the smallest singular value is (near) zero, so repeated
    ``M^-1 v`` converges to the critical mode.  Dependency-free: only the LU solve
    above is used.  ``residual`` reports ``||M x||_inf / ||x||_inf`` after
    normalisation -- a solver author can compare this mode shape *and* check that
    the returned profile really is a marginal solution.
    """
    wall = _normalise_walls(walls)
    matrix = linear_stability_matrix(wavenumber, rayleigh, order, wall)
    size = len(matrix)
    vec = [1.0 + 0.05 * (i % 7) for i in range(size)]
    for _ in range(iterations):
        solution = _lu_solve(matrix, vec)
        norm = max(abs(v) for v in solution) or 1.0
        vec = [v / norm for v in solution]
    theta_part = vec[order + 1 :]
    peak = max(theta_part, key=abs) if theta_part else 0.0
    sign = 1.0 if peak >= 0.0 else -1.0
    amplitude = max(abs(v) for v in theta_part) or 1.0
    psi = [sign * v / amplitude for v in vec[: order + 1]]
    theta = [sign * v / amplitude for v in theta_part]
    return OnsetMode(
        wavenumber=wavenumber,
        rayleigh=rayleigh,
        order=order,
        walls=wall,
        nodes=_chebyshev(order).nodes,
        psi=tuple(psi),
        theta=tuple(theta),
        residual=_residual(matrix, vec),
    )


# --------------------------------------------------------------------------- #
# Solver-side harness: growth-rate probe and Ra bisection
# --------------------------------------------------------------------------- #


def growth_rate(times: Sequence[float], amplitudes: Sequence[float], *, discard: int = 0) -> float:
    """Exponential growth rate of a perturbation amplitude history.

    Least-squares slope of ``log(amplitude)`` against time -- the standard
    ``sigma`` probe: start a solver from a small perturbation, fit the velocity or
    energy amplitude history; negative slope means the base state is stable.
    ``discard`` drops startup samples.
    """
    if len(times) != len(amplitudes):
        raise ValueError("times and amplitudes must have the same length")
    pairs = [
        (t, math.log(v))
        for t, v in list(zip(times, amplitudes, strict=True))[discard:]
        if v > 0.0 and math.isfinite(v)
    ]
    if len(pairs) < 2:
        raise ValueError("need at least two positive samples")
    n = len(pairs)
    mean_t = sum(p[0] for p in pairs) / n
    mean_y = sum(p[1] for p in pairs) / n
    sxx = sum((p[0] - mean_t) ** 2 for p in pairs)
    if sxx == 0.0:
        raise ValueError("times must not all be identical")
    sxy = sum((p[0] - mean_t) * (p[1] - mean_y) for p in pairs)
    return sxy / sxx


@dataclass(frozen=True)
class OnsetVerdict:
    """Outcome of :func:`bisect_onset`."""

    rayleigh_critical: float
    iterations: int
    history: tuple[tuple[float, float], ...]
    stable_side: float
    unstable_side: float
    tolerance_relative: float
    reference: float | None = None

    @property
    def relative_error(self) -> float | None:
        if self.reference is None:
            return None
        return (self.rayleigh_critical - self.reference) / self.reference

    @property
    def accepted(self) -> bool | None:
        error = self.relative_error
        return None if error is None else abs(error) <= self.tolerance_relative

    def as_dict(self) -> dict[str, Any]:
        return {
            "rayleigh_critical": self.rayleigh_critical,
            "iterations": self.iterations,
            "history": [[ra, sigma] for ra, sigma in self.history],
            "stable_side": self.stable_side,
            "unstable_side": self.unstable_side,
            "tolerance_relative": self.tolerance_relative,
            "reference": self.reference,
            "relative_error": self.relative_error,
            "accepted": self.accepted,
        }


def bisect_onset(
    probe: Callable[[float], float],
    *,
    rayleigh_stable: float,
    rayleigh_unstable: float,
    tol: float = 1e-3,
    max_iterations: int = 60,
    reference: float | None = None,
    tolerance_relative: float | None = None,
) -> OnsetVerdict:
    """Bisect ``Ra`` with a solver-supplied growth-rate probe.

    ``probe(Ra)`` must return the measured growth rate ``sigma`` (positive =
    unstable).  The bracketing pair is validated before the loop: a probe that
    does not change sign reports a broken run, not an onset, and raises
    ``ValueError`` instead of fabricating a number.
    """
    if not (0.0 < rayleigh_stable < rayleigh_unstable):
        raise ValueError("need 0 < rayleigh_stable < rayleigh_unstable")
    sigma_stable = probe(rayleigh_stable)
    sigma_unstable = probe(rayleigh_unstable)
    if not (sigma_stable < 0.0 < sigma_unstable):
        raise ValueError(
            f"probe does not bracket an onset: sigma({rayleigh_stable})={sigma_stable}, "
            f"sigma({rayleigh_unstable})={sigma_unstable}"
        )
    tol_rel = (
        float(tolerance_relative)
        if tolerance_relative is not None
        else float(RB_ONSET_ACCEPTANCE["tolerance_relative"])
    )
    lo, hi = rayleigh_stable, rayleigh_unstable
    history: list[tuple[float, float]] = [(lo, sigma_stable), (hi, sigma_unstable)]
    for _ in range(max_iterations):
        mid = 0.5 * (lo + hi)
        sigma = probe(mid)
        history.append((mid, sigma))
        if sigma < 0.0:
            lo = mid
        else:
            hi = mid
        if (hi - lo) <= tol * 0.5 * (hi + lo):
            break
    return OnsetVerdict(
        rayleigh_critical=0.5 * (lo + hi),
        iterations=len(history),
        history=tuple(history),
        stable_side=lo,
        unstable_side=hi,
        tolerance_relative=tol_rel,
        reference=reference,
    )


def nusselt_number(temperature: Sequence[float], dz: float) -> float:
    """Nusselt number from a wall-normal temperature profile (hot wall -> cold wall).

    The profile must include the wall values, i.e. be nodal rather than
    cell-centred: ``temperature[0]`` sits on the hot plate (``z = 0``).

    ``Nu = 1`` is pure conduction.  The hot-wall gradient is taken with a
    second-order one-sided difference and normalised by the conduction gradient
    ``dT / d = 1``.
    """
    if len(temperature) < 3:
        raise ValueError("need at least three samples")
    if dz <= 0.0:
        raise ValueError("dz must be positive")
    # Nu = -(dT/dz)|_hot, normalised by the conduction gradient dT/d = 1
    return (3.0 * temperature[0] - 4.0 * temperature[1] + temperature[2]) / (2.0 * dz)


# --------------------------------------------------------------------------- #
# Cell specification for a 2-D solver
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class RbCell:
    """Two-dimensional Rayleigh-Benard cell a solver can consume directly.

    Domain ``[0, width] x [0, 1]``: bottom plate hot (``T = 1``, ``z = 0``), top
    plate cold (``T = 0``), both isothermal; side walls periodic in ``x`` so that
    one convection roll pair fits ``width``.
    """

    wavenumber: float
    width: float
    height: float
    walls: str
    resolutions: tuple[tuple[int, int], ...]
    boundary: dict[str, Any]
    probe_points: tuple[dict[str, float], ...]
    monitor: dict[str, Any]

    @property
    def roll_wavelength(self) -> float:
        return 2.0 * math.pi / self.wavenumber

    def cell_centres(self, nx: int, nz: int) -> tuple[list[float], list[float]]:
        """Uniform cell-centre coordinates for a structured ``nx x nz`` grid."""
        if nx < 1 or nz < 1:
            raise ValueError("resolution must be positive")
        xs = [(i + 0.5) * self.width / nx for i in range(nx)]
        zs = [(k + 0.5) * self.height / nz for k in range(nz)]
        return xs, zs

    def base_state(self, z: float) -> float:
        """Conduction temperature ``T = 1 - z`` (bottom plate hot)."""
        return 1.0 - z

    def as_dict(self) -> dict[str, Any]:
        return {
            "wavenumber": self.wavenumber,
            "width": self.width,
            "height": self.height,
            "roll_wavelength": self.roll_wavelength,
            "walls": self.walls,
            "resolutions": [list(r) for r in self.resolutions],
            "boundary": self.boundary,
            "probe_points": [dict(p) for p in self.probe_points],
            "monitor": self.monitor,
            "base_state": "T(z) = 1 - z, u = w = 0",
            "nondimensionalisation": NONDIMENSIONALISATION,
        }


def rb_cell(
    *,
    wavenumber: float | None = None,
    walls: str = "rigid",
    resolutions: Sequence[tuple[int, int]] = ((64, 32), (128, 64), (256, 128)),
) -> RbCell:
    """Build the solver-facing cell specification at the critical wavenumber.

    ``wavenumber`` defaults to the reference ``a_c`` of the requested plate pair
    (Chandrasekhar for rigid, exact ``pi/sqrt(2)`` for free), so the cell width is
    the wavelength of one roll pair: ``width = 2 pi / a_c``.
    """
    wall = _normalise_walls(walls)
    if wavenumber is None:
        wavenumber = (
            CHANDRASEKHAR_1961["wavenumber_critical"]
            if wall == _RIGID
            else FREE_FREE_ANALYTIC["wavenumber_critical"]
        )
    if wavenumber <= 0.0:
        raise ValueError("wavenumber must be positive")
    width = 2.0 * math.pi / wavenumber
    kinematic = "no-slip" if wall == _RIGID else "free-slip"
    boundary = {
        "bottom": {
            "locator": "z = 0",
            "type": "wall",
            "kinematic": kinematic,
            "thermal": "T = 1 (isothermal, hot)",
            "tag": "hot_wall",
        },
        "top": {
            "locator": "z = 1",
            "type": "wall",
            "kinematic": kinematic,
            "thermal": "T = 0 (isothermal, cold)",
            "tag": "cold_wall",
        },
        "left": {
            "locator": "x = 0",
            "type": "periodic",
            "pairs_with": "right",
            "tag": "periodic_west",
        },
        "right": {
            "locator": f"x = {width:.10f}",
            "type": "periodic",
            "pairs_with": "left",
            "tag": "periodic_east",
        },
        "gravity": "g = -g z_hat (destabilising, bottom plate hot)",
    }
    probe_points = (
        {"name": "roll_core", "x": 0.5 * width, "z": 0.5},
        {"name": "upwelling", "x": 0.0, "z": 0.5},
        {"name": "quarter_depth", "x": 0.5 * width, "z": 0.25},
    )
    monitor = {
        "perturbation": "w_rms = sqrt(mean((w - w_conduction)^2)) with w_conduction = 0",
        "nusselt": "Nu = -(dT/dz)|_hot / (dT/d)",
        "stability_criterion": (
            "sigma = d log(w_rms)/dt over the exponential window; sigma < 0 stable, "
            "sigma > 0 unstable -- bisect Ra with tensorcad.rb_onset.bisect_onset"
        ),
        "onset_reference": {
            "rigid": CHANDRASEKHAR_1961["rayleigh_critical"],
            "free": FREE_FREE_ANALYTIC["rayleigh_critical"],
        }[wall],
        "conduction_limit": 1.0,
    }
    return RbCell(
        wavenumber=wavenumber,
        width=width,
        height=1.0,
        walls=wall,
        resolutions=tuple((int(nx), int(nz)) for nx, nz in resolutions),
        boundary=boundary,
        probe_points=tuple(probe_points),
        monitor=monitor,
    )


def write_bundle(
    path: str | Path,
    *,
    scan: OnsetScan | None = None,
    mode: OnsetMode | None = None,
    cell: RbCell | None = None,
    extra: dict[str, Any] | None = None,
) -> Path:
    """Write a JSON bundle for solver consumption, with a sha256 digest.

    The digest covers the canonical JSON of the payload (sorted keys), so a
    consumer can detect a hand-edited reference value.
    """
    payload: dict[str, Any] = {
        "case": CASE_SPEC,
        "nondimensionalisation": NONDIMENSIONALISATION,
        "acceptance": RB_ONSET_ACCEPTANCE,
    }
    if scan is not None:
        payload["onset"] = scan.as_dict()
    if mode is not None:
        payload["critical_mode"] = mode.as_dict()
    if cell is not None:
        payload["cell"] = cell.as_dict()
    if extra:
        payload.update(extra)
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    payload["sha256"] = hashlib.sha256(canonical.encode("utf-8")).hexdigest()
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return target
