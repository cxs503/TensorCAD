"""Coefficient-calibrated hull: a parametric surface that *realises* Cb, Cm, Cwp.

Why this module exists
======================
:mod:`tensorcad.hull` lofts a hull from a fixed breadth law and a fixed vertical
profile, so the block coefficient it produces is whatever the exponents happen to
give (about 0.47 for the shipped defaults). A caller who asks for ``Cb = 0.82``
gets a mesh whose volume is *not* ``0.82 * L * B * T``: the geometry cannot honour
the parameter vector. The GAP-CAD-C2 rebuild campaign measured exactly that --
mismatched block coefficient was the single largest failure mode of the parameter
vector -> mesh pipeline.

This module closes the gap by constructing the surface from the *coefficients*
instead of hoping the coefficients come out:

* breadth law ``b(x) = B/2 * sin(pi*x/L)**s``            -> ``Cwp``
* local fullness ``c(x) = Cm * sin(pi*x/L)**r``          -> ``Cb`` (with ``Cm`` at midship)
* section profile ``y(x,z) = b(x) * zeta**(p(x))``       -> local sectional area

Each of ``s``, ``r`` and the per-station ``p`` is obtained by *bisection on the
discrete functional the mesh actually integrates*, not on its continuum limit.
That matters: a hull that is exact in the continuum is still 2-3 % off once the
loft is sampled at 81 x 21, which is precisely the error the rebuild gate rejects.
Solving the discrete functional makes the sampled mesh hit the coefficients to
round-off.

Definitions used throughout (matching :mod:`tensorcad.hull`):

``x``  longitudinal, ``0`` at the aft end, ``L`` at the forward end;
``y``  transverse, ``0`` on the centreplane;
``z``  vertical, ``0`` at the design waterline, ``-draft`` at the keel.

Coefficient definitions realised here:

``Cb  = V / (L*B*T)``                     displacement volume
``Cm  = A_midship / (B*T)``               maximum sectional area
``Cwp = A_waterplane / (L*B)``            design waterplane area

Feasibility
===========
The construction needs ``Cb/Cm <= Cwp`` (i.e. ``Cp <= Cwp``). Every real ship
satisfies it -- the waterplane is fuller than the prismatic coefficient -- and
:meth:`CalibratedHullParameters.problems` refuses the parameter set when it does
not, rather than silently returning a mesh that misses its coefficients.

Known limitations (deliberate, documented rather than hidden)
=============================================================
* The family has **no bulbous bow, no transom, no flare and no hollow waterlines**:
  it is a *coefficient-faithful* hull, not a ship's lines. It is therefore the right
  input for mesh-quality and hydrostatic self-consistency work, and the wrong input
  for resistance or seakeeping claims.
* It does carry two real hull-form features, both off by default so that the
  round-3 family is reproduced exactly: a **parallel middle body**
  (``parallel_body_fraction``) and a **rising keel toward the ends**, i.e. the
  forefoot/run cut-away (``end_draft_ratio``). Both are needed to reproduce a real
  hull's *wetted surface* at fixed coefficients: without the keel rise the family
  gives every station the full design draft, so the ends are full-height blades and
  the surface comes out high (up to +14.7 % against the Holtrop & Mennen regression
  for fine hulls). With the rise the same coefficients are realised by a hollowed
  run and entrance, which is what a real hull looks like.
* ``LCB`` is fixed by the symmetric sine law (``LCB = L/2`` for ``s = r``); only
  symmetric fore/aft families are generated.
* Only ``Cb``, ``Cm`` and ``Cwp`` are honoured. Anything else a caller thinks a
  parameter vector carries (LCB, Cp, L/R ratios beyond ``L, B, T``) is derived.
"""
from __future__ import annotations

import functools
import math
from bisect import bisect_left
from dataclasses import dataclass
from typing import Any

from . import mesh as meshlib
from .geometry import Mesh, MeshBuilder
from .hydrostatics import compute_hydrostatics, max_section_area

UNITS = "m"
COORDINATE_FRAME = (
    "x=longitudinal, from the aft end (x=0) toward the forward end (x=length); "
    "y=transverse, positive to starboard; z=vertical, positive up; "
    "z=0 is the design waterline, the keel lies at z=-draft"
)

#: Exponent of the vertical level clustering, ``zeta_j = (j/(nz-1))**LEVEL_CLUSTER``.
#: ``2.0`` puts the first level at ``1/(nz-1)**2`` above the keel, so the keel point
#: costs only ``0.5 * zeta_1`` of the section area (0.125 % at ``nz = 21``) and a
#: nearly rectangular section (``Cm -> 1``) stays representable.
LEVEL_CLUSTER = 2.0

#: Bounds on the per-station profile exponent solved by bisection/secant.
PROFILE_EXPONENT_MINIMUM = 1.0e-4
PROFILE_EXPONENT_MAXIMUM = 1.0e2

#: Table resolution used to invert ``R(p)`` (area ratio as a function of the
#: profile exponent). ``R`` is universal in ``p`` for a given level count, so the
#: table is built once per ``nz`` and reused for every station of every hull.
_INVERSE_TABLE_SIZE = 2001

#: Below this requested fullness the section is treated as a vanishing point
#: (both ``b`` and the area go to zero at the ends, so the profile exponent is moot).
FULLNESS_FLOOR = 1.0e-12

#: Exponent concentrating the keel rise toward the ends. A real hull keeps a flat
#: keel over the parallel middle body and lifts it through the entrance and the run,
#: so the draft ratio is ``eps + (1-eps) * taper**KEEL_RISE_EXPONENT``: at the
#: middle of a taper body (``taper = 0.5``) the keel has risen only 12.5 % of the way
#: to its end value, which is the shape of a forefoot, not a straight wedge.
KEEL_RISE_EXPONENT = 3.0


# --------------------------------------------------------------------------- #
# Discrete functionals -- the ones the mesh actually integrates
# --------------------------------------------------------------------------- #
def _trapezoid_mean(values: list[float] | tuple[float, ...]) -> float:
    """Mean of ``values`` on a uniform grid, using the trapezoid rule.

    ``(sum - 0.5*(first+last)) / (n-1)``: the trapezoid weights read ``1/2, 1, ...,
    1, 1/2`` normalised by ``n-1``, which is exactly what a lofted mesh of ``n``
    stations integrates.
    """
    n = len(values)
    if n < 2:
        raise ValueError("need at least two samples")
    return (math.fsum(values) - 0.5 * (values[0] + values[-1])) / (n - 1)


def taper_profile(nx: int, parallel_body_fraction: float = 0.0) -> tuple[float, ...]:
    """Longitudinal taper shape ``g_i`` in ``[0, 1]``: the planform of the hull.

    ``g = 1`` over the parallel middle body (the middle ``parallel_body_fraction`` of
    the length) and ``cos(pi*u/2)`` through the entrance and the run, with ``u`` the
    normalised distance from the end of the parallel body to the end of the ship. At
    ``parallel_body_fraction = 0`` this is exactly ``sin(pi*x/L)``, so the round-3
    family is reproduced bit for bit.

    The end stations are pinned to zero *exactly*: ``sin(pi)`` is 1.2e-16 and
    ``1.2e-16 ** 0.07`` is 0.08, so leaving it in would both fatten the stem and --
    worse -- keep the end ring from degenerating, which leaves the loft open and the
    divergence-theorem volume short by a whole end cap (measured: 14 % on a fine
    hull). Both the solve and the loft read this one function, so they cannot drift.
    """
    if not (0.0 <= parallel_body_fraction < 1.0):
        raise ValueError(
            f"parallel_body_fraction must lie in [0, 1), got {parallel_body_fraction!r}"
        )
    last = nx - 1
    half = 0.5 * parallel_body_fraction
    out: list[float] = []
    for i in range(nx):
        if i == 0 or i == last:
            out.append(0.0)
            continue
        xi = i / last
        distance = abs(xi - 0.5)
        if distance <= half:
            out.append(1.0)
        else:
            u = (distance - half) / (0.5 - half)
            out.append(math.cos(0.5 * math.pi * u))
    return tuple(out)


def keel_ratio_profile(
    taper: Sequence[float], end_draft_ratio: float = 1.0
) -> tuple[float, ...]:
    """Local draft / design draft: ``1`` amidships, ``end_draft_ratio`` at the ends.

    This is the forefoot and the run: the keel is flat over the parallel middle body
    and rises through the entrance and the run, with the rise concentrated toward
    the ends by :data:`KEEL_RISE_EXPONENT`. ``end_draft_ratio = 1`` keeps every
    station at the design draft (the round-3 family, no keel rise).
    """
    if not (0.0 < end_draft_ratio <= 1.0):
        raise ValueError(f"end_draft_ratio must lie in (0, 1], got {end_draft_ratio!r}")
    if end_draft_ratio == 1.0:
        return tuple(1.0 for _ in taper)
    return tuple(
        end_draft_ratio + (1.0 - end_draft_ratio) * (g**KEEL_RISE_EXPONENT) for g in taper
    )


def power_mean(samples: Sequence[float], beta: float) -> float:
    """Discrete ``mean`` of ``samples ** beta`` over the station grid.

    The trapezoid weights are the ones the loft itself integrates, so a coefficient
    solved against this mean is realised by the *sampled* mesh, not merely in the
    continuum limit.
    """
    return _trapezoid_mean([max(value, 0.0) ** beta for value in samples])


def weighted_power_mean(
    samples: Sequence[float], beta: float, weights: Sequence[float]
) -> float:
    """``mean(weights * samples ** beta)`` -- the keel-rise form of :func:`power_mean`.

    The sectional area at a station is ``2 * (B/2 * g**s) * (T * h) * Cm * g**r``, so
    ``Cb = Cm * mean(h * g**(s+r))``: the draft ratio ``h`` enters the block
    coefficient as a weight, which is why it needs its own inverse.
    """
    if len(samples) != len(weights):
        raise ValueError("samples and weights must be the same length")
    return _trapezoid_mean(
        [max(value, 0.0) ** beta * weight for value, weight in zip(samples, weights)]
    )


def sine_power_mean(beta: float, nx: int) -> float:
    """Discrete ``mean_xi sin(pi*xi)**beta`` over ``nx`` stations.

    Kept as the round-3 entry point; equal to :func:`power_mean` of
    :func:`taper_profile` at zero parallel body.
    """
    return power_mean(taper_profile(nx, 0.0), beta)


def _solve_power(
    target: float,
    samples: Sequence[float],
    weights: Sequence[float] | None = None,
) -> float:
    """Invert :func:`power_mean` / :func:`weighted_power_mean`.

    The mean falls monotonically from ``max(samples)`` (``beta -> 0``) to ``0``
    (``beta -> inf``), so a bisection in ``log beta`` converges for any target in
    between. ``target`` must be strictly below the ``beta -> 0`` limit, which the
    end stations pin at the value of the interior samples: if every station were
    ``1`` there would be no taper to solve for.
    """
    if not (0.0 < target < 1.0):
        raise ValueError(f"power-mean target must be in (0, 1), got {target!r}")
    weights = list(weights) if weights is not None else [1.0] * len(samples)

    def mean(beta: float) -> float:
        return weighted_power_mean(samples, beta, weights)

    low, high = -12.0, 12.0  # log10(beta)
    if mean(10.0**low) < target:  # pragma: no cover - extreme input
        raise ValueError(f"target {target!r} above the representable power mean")
    if mean(10.0**high) > target:  # pragma: no cover - extreme input
        raise ValueError(f"target {target!r} below the representable power mean")
    for _ in range(200):
        mid = 0.5 * (low + high)
        if mean(10.0**mid) > target:
            low = mid
        else:
            high = mid
        if high - low < 1.0e-14:
            break
    return 10.0 ** (0.5 * (low + high))


def _solve_sine_power(target: float, nx: int) -> float:
    """Round-3 entry point: the sine-law exponent whose discrete mean is ``target``."""
    return _solve_power(target, taper_profile(nx, 0.0))


# --------------------------------------------------------------------------- #
# Section profile: y_j = b * zeta_j**p, area ratio R(p)
# --------------------------------------------------------------------------- #
@functools.lru_cache(maxsize=None)
def _level_zetas(nz: int) -> tuple[float, ...]:
    """Vertical levels from the keel (``zeta = 0``) to the waterline (``zeta = 1``)."""
    last = nz - 1
    return tuple((j / last) ** LEVEL_CLUSTER for j in range(nz))


def section_area_ratio(p: float, nz: int) -> float:
    """Area of the sampled half-section divided by the bounding rectangle ``b * T``.

    ``y_j = b * zeta_j**p`` with ``zeta_j`` from :func:`_level_zetas`; the polygon
    the mesh builds is the trapezoid sum over the panels, and ``R`` is that sum
    itself -- no continuum formula is used anywhere in the calibration.
    """
    zetas = _level_zetas(nz)
    total = 0.0
    for lower, upper in zip(zetas, zetas[1:]):
        total += 0.5 * (lower**p + upper**p) * (upper - lower)
    return total


@functools.lru_cache(maxsize=None)
def _ratio_table(nz: int) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Monotone ``(log10 p, R(p))`` table, decreasing in ``p``."""
    logs: list[float] = []
    ratios: list[float] = []
    low = math.log10(PROFILE_EXPONENT_MINIMUM)
    high = math.log10(PROFILE_EXPONENT_MAXIMUM)
    for k in range(_INVERSE_TABLE_SIZE):
        log_p = low + (high - low) * k / (_INVERSE_TABLE_SIZE - 1)
        logs.append(log_p)
        ratios.append(section_area_ratio(10.0**log_p, nz))
    return tuple(logs), tuple(ratios)


def solve_profile_exponent(fullness: float, nz: int) -> float:
    """Profile exponent whose sampled sectional area is ``fullness * b * T``.

    ``R`` decreases with ``p``, so the inverse is single valued. The table gives a
    first guess and two secant steps refine it with exact evaluations, which keeps
    the realised area within ``~1e-12`` of the target without a bisection per
    station (81 stations x 500 hulls would otherwise dominate the runtime).
    """
    if fullness <= FULLNESS_FLOOR:
        return PROFILE_EXPONENT_MAXIMUM
    zetas = _level_zetas(nz)
    ceiling = 1.0 - 0.5 * zetas[1]  # p -> 0: rectangle minus the keel triangle
    if fullness >= ceiling:
        # A section this full is the rectangle limit; no finite exponent reaches it.
        return PROFILE_EXPONENT_MINIMUM
    logs, ratios = _ratio_table(nz)
    # ratios[] is decreasing: find the bracket containing `fullness`.
    lo, hi = 0, len(ratios) - 1
    if fullness >= ratios[0]:
        return 10.0 ** logs[0]
    if fullness <= ratios[-1]:
        return 10.0 ** logs[-1]
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if ratios[mid] > fullness:
            lo = mid
        else:
            hi = mid
    # Linear interpolation in log p between the bracketing table entries.
    span = ratios[lo] - ratios[hi]
    frac = (ratios[lo] - fullness) / span if span > 0 else 0.0
    log_p = logs[lo] + frac * (logs[hi] - logs[lo])
    p = 10.0**log_p
    # Two secant refinements on the exact (non-interpolated) functional.
    r_prev, p_prev = ratios[lo], 10.0 ** logs[lo]
    for _ in range(2):
        r_now = section_area_ratio(p, nz)
        if r_now == r_prev:
            break
        p_next = p - (r_now - fullness) * (p - p_prev) / (r_now - r_prev)
        if not math.isfinite(p_next) or p_next <= 0.0:
            break
        p_prev, r_prev = p, r_now
        p = p_next
    return min(max(p, PROFILE_EXPONENT_MINIMUM), PROFILE_EXPONENT_MAXIMUM)


# --------------------------------------------------------------------------- #
# Parameters
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class CalibratedHullParameters:
    """Principal dimensions plus the three coefficients the surface must realise."""

    length: float
    beam: float
    draft: float
    block_coefficient: float
    midship_coefficient: float
    waterplane_coefficient: float
    nx: int = 81
    nz: int = 21
    water_density: float = 1025.0
    #: Middle fraction of the length held at full breadth and full section area.
    #: ``0`` (default) is the round-3 sine planform with no parallel body.
    parallel_body_fraction: float = 0.0
    #: Draft at the extreme ends divided by the design draft -- the forefoot and the
    #: run. ``1`` (default) keeps every station at the design draft, i.e. no rise.
    end_draft_ratio: float = 1.0

    @property
    def prismatic_coefficient(self) -> float:
        return self.block_coefficient / self.midship_coefficient

    def problems(self) -> list[str]:
        """Everything that stops this parameter set from being realisable."""
        out: list[str] = []
        for name, value in (
            ("length", self.length),
            ("beam", self.beam),
            ("draft", self.draft),
        ):
            if not math.isfinite(value) or value <= 0.0:
                out.append(f"{name} must be a finite positive number, got {value!r}")
        for name, value in (
            ("block_coefficient", self.block_coefficient),
            ("midship_coefficient", self.midship_coefficient),
            ("waterplane_coefficient", self.waterplane_coefficient),
        ):
            if not math.isfinite(value) or not (0.0 < value < 1.0):
                out.append(f"{name} must lie in (0, 1), got {value!r}")
        if out:
            return out
        if self.nx < 3 or self.nz < 3:
            out.append(f"nx and nz must be >= 3, got nx={self.nx}, nz={self.nz}")
        cp = self.prismatic_coefficient
        if cp > self.waterplane_coefficient:
            out.append(
                "inconsistent coefficients: Cp = Cb/Cm = "
                f"{cp:.4f} exceeds Cwp = {self.waterplane_coefficient:.4f}; "
                "the waterplane cannot be finer than the prismatic coefficient"
            )
        if not (0.0 < cp < 1.0):
            out.append(f"prismatic coefficient Cb/Cm must lie in (0, 1), got {cp!r}")
        if not (0.0 <= self.parallel_body_fraction < 1.0):
            out.append(
                "parallel_body_fraction must lie in [0, 1), got "
                f"{self.parallel_body_fraction!r}"
            )
        if not (0.0 < self.end_draft_ratio <= 1.0):
            out.append(
                f"end_draft_ratio must lie in (0, 1], got {self.end_draft_ratio!r}"
            )
        if out:
            return out
        # A keel rise hollows out the ends, so the fullness law has to be *keener*
        # than the breadth law to still hold Cb: when the required fullness exponent
        # would go negative the parameter set is geometrically impossible (the ends
        # are asked to carry more volume than they have box), and it is refused here
        # rather than silently returned as a hull that misses its block coefficient.
        try:
            if self.fullness_exponents()[1] < 0.0:
                out.append(
                    "inconsistent shape request: end_draft_ratio "
                    f"{self.end_draft_ratio:.3f} with parallel_body_fraction "
                    f"{self.parallel_body_fraction:.3f} hollows the ends so far that "
                    "no positive fullness exponent realises Cb; use a larger "
                    "end_draft_ratio or a fuller hull"
                )
        except ValueError as error:
            out.append(
                f"inconsistent shape request: end_draft_ratio "
                f"{self.end_draft_ratio:.3f} with parallel_body_fraction "
                f"{self.parallel_body_fraction:.3f} leaves the ends too little volume "
                f"for Cb = {self.block_coefficient:.3f} ({error}); use a larger "
                "end_draft_ratio or a fuller hull"
            )
        return out

    def taper(self) -> tuple[float, ...]:
        """The planform shape sampled on this hull's station grid."""
        return taper_profile(self.nx, self.parallel_body_fraction)

    def keel_ratios(self) -> tuple[float, ...]:
        """Draft / design draft at each station."""
        return keel_ratio_profile(self.taper(), self.end_draft_ratio)

    def fullness_exponents(self) -> tuple[float, float]:
        """``(s, r)``: breadth exponent from ``Cwp``, fullness exponent from ``Cb``.

        Both are solved on the discrete functional the loft integrates:
        ``Cwp = mean(g**s)`` and ``Cb = Cm * mean(h * g**(s+r))``.
        """
        taper = self.taper()
        keel = self.keel_ratios()
        s = _solve_power(self.waterplane_coefficient, taper)
        q = _solve_power(self.prismatic_coefficient, taper, keel)
        return s, q - s

    def as_dict(self) -> dict[str, Any]:
        return {
            "length": self.length,
            "beam": self.beam,
            "draft": self.draft,
            "block_coefficient": self.block_coefficient,
            "midship_coefficient": self.midship_coefficient,
            "waterplane_coefficient": self.waterplane_coefficient,
            "parallel_body_fraction": self.parallel_body_fraction,
            "end_draft_ratio": self.end_draft_ratio,
            "prismatic_coefficient": self.prismatic_coefficient,
            "nx": self.nx,
            "nz": self.nz,
        }


@dataclass(frozen=True)
class CalibratedHull:
    """The generated surface plus the calibration it used."""

    mesh: Mesh
    parameters: CalibratedHullParameters
    #: Exponent of the breadth law solved from ``Cwp``.
    breadth_exponent: float
    #: Exponent of the fullness law solved from ``Cb`` (``Cm`` at midship).
    fullness_exponent: float
    #: Per-station half-breadth and profile exponent actually used.
    half_breadths: tuple[float, ...]
    profile_exponents: tuple[float, ...]
    #: Per-station draft / design draft actually used (empty = flat keel).
    keel_ratios: tuple[float, ...] = ()

    def calibration_report(self) -> dict[str, Any]:
        """Requested vs realised coefficients, integrated from the mesh itself."""
        realised = realised_coefficients(
            self.mesh,
            length=self.parameters.length,
            beam=self.parameters.beam,
            draft=self.parameters.draft,
        )
        requested = {
            "Cb": self.parameters.block_coefficient,
            "Cm": self.parameters.midship_coefficient,
            "Cwp": self.parameters.waterplane_coefficient,
        }
        return {
            "requested": requested,
            "realised": realised["coefficients"],
            "relative_error": relative_errors(realised["coefficients"], requested),
            "breadth_exponent": self.breadth_exponent,
            "fullness_exponent": self.fullness_exponent,
            "settings": {
                "nx": self.parameters.nx,
                "nz": self.parameters.nz,
                "level_cluster": LEVEL_CLUSTER,
                "parallel_body_fraction": self.parameters.parallel_body_fraction,
                "end_draft_ratio": self.parameters.end_draft_ratio,
                "keel_rise_exponent": KEEL_RISE_EXPONENT,
            },
        }


# --------------------------------------------------------------------------- #
# Realised coefficients -- always integrated from the mesh, never from the inputs
# --------------------------------------------------------------------------- #
def realised_coefficients(mesh: Mesh, *, length: float, beam: float, draft: float) -> dict[str, Any]:
    """Integrate ``Cb``, ``Cm`` and ``Cwp`` off the triangles.

    Volume comes from the divergence theorem, ``Cm`` from the largest sampled
    station (:func:`tensorcad.hydrostatics.max_section_area`) and ``Cwp`` from the
    waterplane flux in :func:`tensorcad.hydrostatics.compute_hydrostatics`, so the
    three numbers are *measurements of the mesh*, independent of the laws that
    built it. That is the whole point: the calibration is verified, not assumed.
    """
    if length <= 0.0 or beam <= 0.0 or draft <= 0.0:
        raise ValueError("length, beam and draft must be positive")
    hydro = compute_hydrostatics(mesh, 0.0, water_density=1025.0)
    volume = meshlib.total_signed_volume(mesh)
    cb = volume / (length * beam * draft)
    cm = max_section_area(mesh) / (beam * draft)
    cwp = hydro.waterplane_area / (length * beam)
    return {
        "coefficients": {"Cb": cb, "Cm": cm, "Cwp": cwp},
        "volume_m3": volume,
        "wetted_surface_m2": hydro.wetted_surface_area,
        "waterplane_area_m2": hydro.waterplane_area,
        "max_section_area_m2": max_section_area(mesh),
        "relative_error": {},
    }


def relative_errors(realised: dict[str, float], requested: dict[str, float]) -> dict[str, float]:
    """Symmetric relative error per coefficient, keyed like the inputs."""
    out: dict[str, float] = {}
    for key, want in requested.items():
        if key not in realised:
            continue
        got = realised[key]
        if want == 0.0:
            out[key] = abs(got)
        else:
            out[key] = abs(got - want) / abs(want)
    return out


# --------------------------------------------------------------------------- #
# Generation
# --------------------------------------------------------------------------- #
def _section_point(
    x: float, half_breadth: float, profile_exponent: float, draft: float, nz: int,
    j: int,
) -> tuple[float, float, float]:
    """Ring position ``j`` of one station, walking starboard -> keel -> port.

    The keel vertex is shared by both flanks (it appears once in the ring), which
    is what makes the loft closed without a welding pass -- the same construction
    :mod:`tensorcad.hull` uses, and the reason its Euler characteristic comes out
    at 2.
    """
    zetas = _level_zetas(nz)
    keelward = j <= nz - 1
    level = (nz - 1 - j) if keelward else (j - (nz - 1))
    zeta = zetas[level]
    z = -draft * (1.0 - zeta)
    breadth = half_breadth * (zeta**profile_exponent)
    return (x, breadth if keelward else -breadth, z)


def build_calibrated_hull(parameters: CalibratedHullParameters) -> CalibratedHull:
    """Loft a hull that realises ``Cb``, ``Cm`` and ``Cwp`` as asked."""
    problems = parameters.problems()
    if problems:
        raise ValueError("; ".join(problems))
    p = parameters
    nx, nz = p.nx, p.nz

    # 1) Breadth law from Cwp and 2) fullness law from Cb, both on the discrete
    #    functional the loft integrates: Cwp = mean(g**s), Cb = Cm*mean(h*g**(s+r)),
    #    with ``g`` the planform taper (parallel middle body) and ``h`` the local
    #    draft ratio (forefoot / run).
    taper = p.taper()
    keel_ratios = p.keel_ratios()
    breadth_exponent, fullness_exponent = p.fullness_exponents()
    if fullness_exponent < 0.0:  # pragma: no cover - problems() rejects this request
        raise ValueError("no realisable fullness law for this shape request")

    half_breadths: list[float] = []
    fullnesses: list[float] = []
    for g in taper:
        half_breadths.append(0.5 * p.beam * g**breadth_exponent)
        fullnesses.append(p.midship_coefficient * g**fullness_exponent)

    # 3) Per-station section profile so the sampled area is exactly fullness*B*T.
    profile_exponents = [solve_profile_exponent(f, nz) for f in fullnesses]

    ring_size = 2 * nz - 1
    builder = MeshBuilder(units=UNITS, coordinate_frame=COORDINATE_FRAME)
    rings: list[list[int]] = []
    for i in range(nx):
        x = p.length * i / (nx - 1)
        station_draft = p.draft * keel_ratios[i]
        rings.append([
            builder.add_vertex(
                _section_point(
                    x, half_breadths[i], profile_exponents[i], station_draft, nz, j
                )
            )
            for j in range(ring_size)
        ])

    for i in range(nx - 1):
        for j in range(ring_size - 1):
            builder.add_quad(rings[i][j], rings[i + 1][j], rings[i + 1][j + 1], rings[i][j + 1])
        # Deck: the strip between the two waterline ends of consecutive rings.
        builder.add_quad(
            rings[i][ring_size - 1],
            rings[i + 1][ring_size - 1],
            rings[i + 1][0],
            rings[i][0],
        )

    mesh = meshlib.orient_outward(builder.build()).mesh
    return CalibratedHull(
        mesh=mesh,
        parameters=p,
        breadth_exponent=breadth_exponent,
        fullness_exponent=fullness_exponent,
        half_breadths=tuple(half_breadths),
        profile_exponents=tuple(profile_exponents),
        keel_ratios=tuple(keel_ratios),
    )


def sectional_area_curve(hull: CalibratedHull) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """``(x, A(x))`` of the generated hull, straight from the per-station profile.

    ``A_i = 2 * b_i * T * R(p_i)``: the sampled polygon area, not its analytic
    limit, so the curve integrates to the mesh volume by construction.
    """
    p = hull.parameters
    nz = p.nz
    xs = tuple(p.length * i / (p.nx - 1) for i in range(p.nx))
    drafts = hull.keel_ratios or tuple(1.0 for _ in hull.half_breadths)
    areas = tuple(
        2.0 * b * p.draft * h * section_area_ratio(pe, nz)
        for b, h, pe in zip(hull.half_breadths, drafts, hull.profile_exponents)
    )
    return xs, areas