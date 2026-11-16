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
* local fullness ``c(x) = Cm * sin(pi*x/L)**r(x)``     -> ``Cb`` (with ``Cm`` at midship),
  where ``r(x) = r + skew`` aft of midship and ``r - skew`` forward of it
                                                        -> ``LCB``
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
* ``LCB`` **is** controllable: ``lcb_offset`` splits the fullness exponent into an
  aft value ``r + skew`` and a fore value ``r - skew`` at the midship station, so the
  fore/aft volume centroid moves while ``Cb``, ``Cm`` and ``Cwp`` stay exact
  (:func:`solve_lcb_skew`). The default ``lcb_offset = 0`` reproduces the symmetric
  round-3 family bit for bit. What is still *not* honoured is a hollow/knuckled
  waterline or a transom: this is a two-parameter fore/aft family, not a lines plan.
* Only ``Cb``, ``Cm``, ``Cwp`` and (optionally) ``LCB`` are honoured. Anything else a
  caller thinks a parameter vector carries (``Cp``, ``LCF``, ``L/R`` ratios beyond
  ``L, B, T``) is derived.
"""
from __future__ import annotations

import functools
import math
from bisect import bisect_left
from dataclasses import dataclass
from typing import Any, Sequence

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

#: Largest ``|LCB offset|`` from mid-length, as a fraction of L, that the family will
#: attempt. Real ships sit inside +/-3.5 % L (KVLCC2 is +3.48 %, KCS is -1.48 %); 5 %
#: is a generous bound that still refuses a request no hull would make.
LCB_OFFSET_LIMIT = 0.05

#: Bound on the fore/aft tilt used by the outer LCB bisection. A tilt of 1 would
#: collapse one half's exponents to zero (a full-breadth slab), so the search stays
#: strictly inside and an unreachable request fails loudly instead of degenerating.
TILT_LIMIT = 0.9

#: Smallest exponent the tilted family will use. A non-positive breadth or fullness
#: exponent makes an end section blunter than midship, which breaks ``Cm`` as the
#: maximum sectional area, so the bracket is closed rather than silently violated.
MIN_FULLNESS_EXPONENT = 1.0e-3


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
# Fore/aft asymmetry: LCB control at fixed Cb, Cm, Cwp
# --------------------------------------------------------------------------- #
class LcbRequestUnreachable(ValueError):
    """The requested LCB cannot be reached without breaking Cb, Cm or Cwp."""


def section_side_signs(nx: int) -> tuple[float, ...]:
    """``+1`` aft of midship, ``-1`` forward of it, ``0`` at the midship station.

    This is the multiplier of the fore/aft tilt: under a tilt ``tau`` the breadth
    exponent of station ``i`` is ``s * (1 + side_i * tau)`` and the fullness exponent
    is ``r * (1 + side_i * tau)``. A positive tilt therefore stretches both exponents
    on the run and shrinks them on the entrance -- slimmer afterbody, fuller
    forebody -- which moves the centre of buoyancy forward. The midship station keeps
    the untilted exponents, so ``Cm`` stays realised there.
    """
    if nx < 2:
        raise ValueError(f"need at least two stations, got nx={nx}")
    mid = 0.5 * (nx - 1)
    centre = (nx - 1) // 2 if nx % 2 == 1 else None
    return tuple(
        0.0 if i == centre else (1.0 if i < mid else -1.0) for i in range(nx)
    )


def asymmetric_mean(
    taper: Sequence[float],
    exponent: float,
    tilt: float,
    weights: Sequence[float] | None = None,
) -> float:
    """``mean(weights * g**(exponent * (1 + side*tilt)))`` -- the tilted mean.

    ``tilt = 0`` reproduces :func:`weighted_power_mean` exactly, so the symmetric
    round-3 family is a special case rather than a re-derivation.
    """
    sides = section_side_signs(len(taper))
    ws = list(weights) if weights is not None else [1.0] * len(taper)
    if len(ws) != len(taper):
        raise ValueError("taper and weights must be the same length")
    return _trapezoid_mean(
        [
            weight * max(g, 0.0) ** (exponent * (1.0 + side * tilt))
            for g, side, weight in zip(taper, sides, ws)
        ]
    )


def solve_asymmetric_exponents(
    taper: Sequence[float],
    keel: Sequence[float],
    waterplane_coefficient: float,
    prismatic_coefficient: float,
    tilt: float,
) -> tuple[float, float]:
    """``(s, r)`` reproducing ``Cwp`` and ``Cp`` at a given fore/aft ``tilt``.

    Two nested bisections on the discrete functionals the loft integrates:

    ``Cwp = mean(g**(s*(1+side*tilt)))``   and
    ``Cp  = mean(h * g**((s+r)*(1+side*tilt)))``

    Both means fall monotonically with their exponent, so each solve is a bisection in
    ``log`` exponent. A target outside the reachable band raises
    :class:`LcbRequestUnreachable` instead of being clamped -- a silently clamped hull
    is exactly the failure mode this module exists to remove.
    """
    if not (-1.0 < tilt < 1.0):
        raise LcbRequestUnreachable(f"tilt must lie in (-1, 1), got {tilt!r}")
    sides = section_side_signs(len(taper))
    factors = [1.0 + side * tilt for side in sides]
    if any(float(f) <= 0.0 for f in factors):  # pragma: no cover - guarded above
        raise LcbRequestUnreachable(f"tilt {tilt!r} flips an exponent's sign")

    def solve(target: float, base: Sequence[float], what: str) -> float:
        """Bisection on the exponent whose tilted mean equals ``target``.

        ``base`` is everything already known at this station (the draft ratio, and the
        breadth law once it is solved), so each call adds just ``g**(exponent*factor)``
        -- the term being solved for. Splitting it this way is what keeps the two
        solves independent; folding ``s`` into the ``r`` solve would silently drop the
        breadth factor and miss ``Cb`` by percent-level.
        """

        def mean_of(exponent: float) -> float:
            return _trapezoid_mean(
                [
                    value * max(g, 0.0) ** (exponent * factor)
                    for g, value, factor in zip(taper, base, factors)
                ]
            )

        low, high = math.log10(MIN_FULLNESS_EXPONENT), math.log10(
            PROFILE_EXPONENT_MAXIMUM
        )
        if not (mean_of(10.0**low) > target > mean_of(10.0**high)):
            raise LcbRequestUnreachable(
                f"{what} = {target!r} is outside the band "
                f"({mean_of(10.0**high):.6f}, {mean_of(10.0**low):.6f}) reachable at "
                f"tilt = {tilt:.4f}"
            )
        for _ in range(200):
            mid = 0.5 * (low + high)
            if mean_of(10.0**mid) > target:
                low = mid
            else:
                high = mid
            if high - low < 1.0e-14:
                break
        return 10.0 ** (0.5 * (low + high))

    breadth_exponent = solve(waterplane_coefficient, [1.0] * len(taper), "Cwp")
    # The block coefficient sees the breadth law the waterplane solve just fixed, so
    # the draft ratio and ``g**s_i`` both enter as the base of the second solve.
    block_base = [
        abs(h) * max(g, 0.0) ** (breadth_exponent * factor)
        for g, h, factor in zip(taper, keel, factors)
    ]
    fullness_exponent = solve(prismatic_coefficient, block_base, "Cp")
    return breadth_exponent, fullness_exponent


def station_volume_weights(
    taper: Sequence[float],
    keel: Sequence[float],
    breadth_exponent: float,
    fullness_exponent: float,
    tilt: float,
) -> tuple[float, ...]:
    """``A_i / (B*T*Cm)``: the sectional areas the tilted loft integrates.

    ``A_i = (B*g_i**s_i) * (T*h_i) * Cm*g_i**r_i`` with ``s_i = s*(1+side*tilt)``,
    so the normalised station areas are ``h_i * g_i**((s+r)*(1+side*tilt))``. They
    carry the whole fore/aft distribution, hence the whole LCB, with no scale factor
    left to get wrong.
    """
    sides = section_side_signs(len(taper))
    if len(keel) != len(taper):
        raise ValueError("taper and keel must be the same length")
    total = breadth_exponent + fullness_exponent
    return tuple(
        h * g ** (total * (1.0 + side * tilt))
        for g, h, side in zip(taper, keel, sides)
    )


def longitudinal_offset(weights: Sequence[float]) -> float:
    """First moment of a sampled station-area curve, as a fraction of L, ``+`` forward.

    ``offset = sum(w_i * (x_i - L/2) * A_i) / sum(w_i * A_i)`` with trapezoid
    weights: the same integral the lofted mesh carries, so the value is the
    construction's own LCB, not a separate approximation of it.
    """
    n = len(weights)
    if n < 2:
        raise ValueError("need at least two stations")
    trapezoid = [1.0] * n
    trapezoid[0] = trapezoid[-1] = 0.5
    mid = 0.5 * (n - 1)
    moment = math.fsum(
        w * ((i - mid) / (n - 1)) * value
        for i, (w, value) in enumerate(zip(trapezoid, weights))
    )
    total = math.fsum(w * value for w, value in zip(trapezoid, weights))
    if total <= 0.0:
        raise ValueError("the station-area curve carries no volume")
    return moment / total


def lcb_reach(
    taper: Sequence[float],
    keel: Sequence[float],
    waterplane_coefficient: float,
    prismatic_coefficient: float,
    limit: float = TILT_LIMIT,
) -> float:
    """Largest ``|LCB offset|`` this coefficient set can hold, as a fraction of L.

    The tilted solve refuses a coefficient set once ``Cp`` leaves the band the tilted
    exponents can reach, so the family has a *saturation* offset: past it the centroid
    cannot be moved further without missing ``Cb``. That number is a property of the
    geometry, not of the API, and is what tells a caller whether a requested LCB is
    inside the envelope before anything is built. Fast: one bisection on the tilt, each
    step two exponent solves.
    """
    def offset_at(tilt: float) -> float | None:
        try:
            breadth, fullness = solve_asymmetric_exponents(
                taper, keel, waterplane_coefficient, prismatic_coefficient, tilt
            )
        except LcbRequestUnreachable:
            return None
        return longitudinal_offset(
            station_volume_weights(taper, keel, breadth, fullness, tilt)
        )

    if not (0.0 < limit < 1.0):
        raise ValueError(f"limit must lie in (0, 1), got {limit!r}")
    lo, hi = 0.0, limit
    if offset_at(hi) is not None:
        return abs(offset_at(hi) or 0.0)
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if offset_at(mid) is not None:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1.0e-12:
            break
    reached = offset_at(lo)
    return 0.0 if reached is None else abs(reached)


def solve_lcb_tilt(
    target_offset: float,
    taper: Sequence[float],
    keel: Sequence[float],
    waterplane_coefficient: float,
    prismatic_coefficient: float,
    limit: float = TILT_LIMIT,
) -> float:
    """Fore/aft tilt whose realised LCB offset is ``target_offset`` (fraction of L).

    Outer bisection on the tilt, inner solves of ``(s, r)`` for ``Cwp`` and ``Cp``: the
    realised offset rises monotonically with the tilt, so the inversion is
    single-valued. ``target_offset = 0`` returns ``0.0`` without touching the search,
    which is what keeps the symmetric round-3 family bit-for-bit reproducible.
    """
    if not math.isfinite(target_offset) or abs(target_offset) >= LCB_OFFSET_LIMIT:
        raise LcbRequestUnreachable(
            f"LCB offset {target_offset!r} is outside (-{LCB_OFFSET_LIMIT}, "
            f"{LCB_OFFSET_LIMIT}) of L"
        )
    if target_offset == 0.0:
        return 0.0
    sign = 1.0 if target_offset > 0.0 else -1.0
    magnitude = abs(target_offset)

    def offset_at(tilt: float) -> float | None:
        try:
            breadth, fullness = solve_asymmetric_exponents(
                taper, keel, waterplane_coefficient, prismatic_coefficient, tilt
            )
        except LcbRequestUnreachable:
            return None
        return longitudinal_offset(
            station_volume_weights(taper, keel, breadth, fullness, tilt)
        )

    reach = limit
    edge: float | None = None
    while reach > 1.0e-6:
        edge = offset_at(sign * reach)
        if edge is not None:
            break
        reach *= 0.5
    if edge is None:  # pragma: no cover - every coefficient set admits some tilt
        raise LcbRequestUnreachable("no admissible tilt for this coefficient set")
    if abs(edge) < magnitude:
        raise LcbRequestUnreachable(
            f"LCB offset {target_offset:+.4f} L is unreachable for Cb/Cm/Cwp = "
            f"({prismatic_coefficient!r}): the family saturates at {abs(edge):+.4f} L "
            "before the deck runs out of breadth"
        )
    lo, hi = 0.0, reach
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        value = offset_at(sign * mid)
        if value is None or abs(value) < magnitude:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1.0e-15:
            break
    return sign * 0.5 * (lo + hi)


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
    #: Longitudinal centre of buoyancy, as a signed fraction of ``length`` measured
    #: forward from mid-length (``+3.48 %`` is KVLCC2, ``-1.48 %`` is KCS). ``0``
    #: (default) is the symmetric round-3 family. Realised by splitting the fullness
    #: exponent fore/aft at midship; ``Cb``, ``Cm`` and ``Cwp`` stay exact.
    lcb_offset: float = 0.0

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
        if not math.isfinite(self.lcb_offset) or abs(self.lcb_offset) >= LCB_OFFSET_LIMIT:
            out.append(
                f"lcb_offset must be a finite fraction of L with |offset| < "
                f"{LCB_OFFSET_LIMIT}, got {self.lcb_offset!r}"
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

    def lcb_tilt(self) -> float:
        """Fore/aft fullness skew that realises :attr:`lcb_offset` (``0`` when unasked).

        Solved on the same discrete functional the loft integrates, so the realised
        LCB is the requested one rather than a continuum estimate of it. Raises
        :class:`LcbRequestUnreachable` when the coefficients leave no room to move the
        centroid far enough.
        """
        if self.lcb_offset == 0.0:
            return 0.0
        return solve_lcb_tilt(
            self.lcb_offset,
            self.taper(),
            self.keel_ratios(),
            self.waterplane_coefficient,
            self.prismatic_coefficient,
        )

    def tilted_exponents(self, tilt: float) -> tuple[float, float]:
        """``(s, r)`` for a given fore/aft tilt: breadth from ``Cwp``, fullness from ``Cb``.

        ``tilt = 0`` reproduces :meth:`fullness_exponents` to round-off, and then the
        tilt is carried per station by :func:`station_volume_weights`.
        """
        return solve_asymmetric_exponents(
            self.taper(),
            self.keel_ratios(),
            self.waterplane_coefficient,
            self.prismatic_coefficient,
            tilt,
        )

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
            "lcb_offset": self.lcb_offset,
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
    #: Fore/aft tilt actually applied (``0`` = the symmetric family).
    lcb_tilt: float = 0.0
    #: Per-station ``A_i / (B*T*Cm)`` the loft integrates, in station order.
    station_weights: tuple[float, ...] = ()

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
        requested_offset = self.parameters.lcb_offset
        realised_offset = realised["lcb_fraction"]
        return {
            "requested": requested,
            "realised": realised["coefficients"],
            "relative_error": relative_errors(realised["coefficients"], requested),
            "lcb": {
                "requested_offset": requested_offset,
                "realised_offset_mesh": realised_offset,
                # The construction's own moment: an independent route to the same
                # number, so a bug in either one shows up as a disagreement.
                "realised_offset_sections": (
                    longitudinal_offset(self.station_weights)
                    if self.station_weights
                    else 0.0
                ),
                "tilt": self.lcb_tilt,
            },
            "breadth_exponent": self.breadth_exponent,
            "fullness_exponent": self.fullness_exponent,
            "settings": {
                "nx": self.parameters.nx,
                "nz": self.parameters.nz,
                "level_cluster": LEVEL_CLUSTER,
                "parallel_body_fraction": self.parameters.parallel_body_fraction,
                "end_draft_ratio": self.parameters.end_draft_ratio,
                "keel_rise_exponent": KEEL_RISE_EXPONENT,
                "lcb_offset": requested_offset,
            },
        }


# --------------------------------------------------------------------------- #
# Realised coefficients -- always integrated from the mesh, never from the inputs
# --------------------------------------------------------------------------- #
def realised_lcb_fraction(mesh: Mesh) -> float:
    """Volume centroid of a closed mesh, as a signed fraction of its length.

    Divergence theorem on the triangle soup: every face ``(a, b, c)`` closes a signed
    tetrahedron ``(0, a, b, c)`` of volume ``a . (b x c) / 6`` whose centroid is
    ``(a + b + c) / 4``. The volume-weighted mean of those centroids is the centre of
    buoyancy of the solid, and referring it to the mesh's own longitudinal mid-point
    keeps the number independent of where the generator put ``x = 0``. Positive is
    forward of mid-length.

    Pure Python and dependency-free on purpose: this is the *measurement* side of the
    LCB contract, kept apart from the construction that has to satisfy it.
    """
    vertices = mesh.vertices
    faces = mesh.faces
    if not faces:
        raise ValueError("mesh has no faces")
    volume = 0.0
    moment = 0.0
    for i, j, k in faces:
        ax, ay, az = vertices[i]
        bx, by, bz = vertices[j]
        cx, cy, cz = vertices[k]
        signed = (
            ax * (by * cz - bz * cy)
            - ay * (bx * cz - bz * cx)
            + az * (bx * cy - by * cx)
        ) / 6.0
        volume += signed
        moment += signed * (ax + bx + cx) / 4.0
    if abs(volume) < 1.0e-12:
        raise ValueError("mesh encloses no volume")
    xs = [vertex[0] for vertex in vertices]
    x_min, x_max = min(xs), max(xs)
    if x_max - x_min <= 0.0:
        raise ValueError("mesh has zero longitudinal extent")
    return (moment / volume - 0.5 * (x_min + x_max)) / (x_max - x_min)


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
        # Measured, not derived from the request: the LCB contract gets the same
        # treatment as the three coefficients.
        "lcb_fraction": realised_lcb_fraction(mesh),
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
    sides = section_side_signs(nx)
    if p.lcb_offset == 0.0:
        # Untouched round-3 path: one symmetric fullness exponent, solved exactly as
        # before, so the historical family is reproduced bit for bit.
        tilt = 0.0
        breadth_exponent, fullness_exponent = p.fullness_exponents()
    else:
        tilt = p.lcb_tilt()
        breadth_exponent, fullness_exponent = p.tilted_exponents(tilt)
    if fullness_exponent < 0.0:  # pragma: no cover - problems() rejects this request
        raise ValueError("no realisable fullness law for this shape request")

    half_breadths: list[float] = []
    fullnesses: list[float] = []
    for g, side in zip(taper, sides):
        # The tilt stretches both exponents on the run and shrinks them on the
        # entrance for ``tilt > 0``; the midship station keeps ``(s, r)`` untouched,
        # so ``Cm`` is still realised there.
        half_breadths.append(0.5 * p.beam * g ** (breadth_exponent * (1.0 + side * tilt)))
        fullnesses.append(
            p.midship_coefficient * g ** (fullness_exponent * (1.0 + side * tilt))
        )

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
        lcb_tilt=tilt,
        station_weights=station_volume_weights(
            taper, keel_ratios, breadth_exponent, fullness_exponent, tilt
        ),
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