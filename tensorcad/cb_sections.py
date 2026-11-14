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
* The family has no parallel middle body, no bulbous bow, no transom, no flare and
  no hollow waterlines: it is a *coefficient-faithful* hull, not a ship's lines.
  It is therefore the right input for mesh-quality and hydrostatic self-consistency
  work, and the wrong input for resistance or seakeeping claims.
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


def sine_power_mean(beta: float, nx: int) -> float:
    """Discrete ``mean_xi sin(pi*xi)**beta`` over ``nx`` stations -- the quantity a
    sine-law breadth distribution contributes to ``Cwp`` (and, with ``Cm`` in
    front, to ``Cb``).
    """
    last = nx - 1
    values = [
        # The end stations are pinned to zero exactly as the generator pins them
        # (``sin(pi)`` is 1.2e-16, and ``1.2e-16 ** 0.07`` is 0.08 -- leaving it in
        # would make the solve and the loft disagree by the whole end panel).
        0.0 if i == 0 or i == last else math.sin(math.pi * i / last) ** beta
        for i in range(nx)
    ]
    return _trapezoid_mean(values)


def _solve_sine_power(target: float, nx: int) -> float:
    """Invert :func:`sine_power_mean`: the exponent whose discrete mean is ``target``.

    The mean falls monotonically from ``1`` (``beta -> 0``) to ``0``
    (``beta -> inf``), so a plain bisection in ``log beta`` converges for any
    ``0 < target < 1``. The endpoints of the station grid contribute
    ``sin(0) = sin(pi) = 0``, hence the ``< 1`` ceiling.
    """
    if not (0.0 < target < 1.0):
        raise ValueError(f"sine-power mean target must be in (0, 1), got {target!r}")
    low, high = -12.0, 12.0  # log10(beta)
    # Guard: make sure the bracket really straddles the target.
    if sine_power_mean(10.0**low, nx) < target:  # pragma: no cover - extreme input
        raise ValueError(f"target {target!r} above the representable sine-power mean")
    if sine_power_mean(10.0**high, nx) > target:  # pragma: no cover - extreme input
        raise ValueError(f"target {target!r} below the representable sine-power mean")
    for _ in range(200):
        mid = 0.5 * (low + high)
        if sine_power_mean(10.0**mid, nx) > target:
            low = mid
        else:
            high = mid
        if high - low < 1.0e-14:
            break
    return 10.0 ** (0.5 * (low + high))


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
        return out

    def as_dict(self) -> dict[str, Any]:
        return {
            "length": self.length,
            "beam": self.beam,
            "draft": self.draft,
            "block_coefficient": self.block_coefficient,
            "midship_coefficient": self.midship_coefficient,
            "waterplane_coefficient": self.waterplane_coefficient,
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
            "settings": {"nx": self.parameters.nx, "nz": self.parameters.nz,
                         "level_cluster": LEVEL_CLUSTER},
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

    # 1) Breadth law from Cwp: the trapezoid mean of sin^s over the station grid.
    breadth_exponent = _solve_sine_power(p.waterplane_coefficient, nx)

    # 2) Fullness law from Cb: Cb = Cm * mean(sin**(s+r)).
    target = p.prismatic_coefficient
    fullness_exponent = _solve_sine_power(target, nx) - breadth_exponent
    if fullness_exponent < 0.0:  # pragma: no cover - problems() rejects Cp > Cwp
        raise ValueError("Cp exceeds Cwp: no realisable fullness law")

    last = nx - 1
    half_breadths: list[float] = []
    fullnesses: list[float] = []
    for i in range(nx):
        # ``math.sin(math.pi)`` is 1.2e-16, not zero: the end stations are forced
        # to a point explicitly, otherwise the loft leaves an open slot at the stem
        # and the post (the same trap ``tensorcad.hull`` documents).
        sine = 0.0 if i == 0 or i == last else math.sin(math.pi * i / last)
        half_breadths.append(0.5 * p.beam * sine**breadth_exponent)
        fullnesses.append(p.midship_coefficient * sine**fullness_exponent)

    # 3) Per-station section profile so the sampled area is exactly fullness*B*T.
    profile_exponents = [solve_profile_exponent(f, nz) for f in fullnesses]

    ring_size = 2 * nz - 1
    builder = MeshBuilder(units=UNITS, coordinate_frame=COORDINATE_FRAME)
    rings: list[list[int]] = []
    for i in range(nx):
        x = p.length * i / last
        rings.append([
            builder.add_vertex(
                _section_point(x, half_breadths[i], profile_exponents[i], p.draft, nz, j)
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
    )


def sectional_area_curve(hull: CalibratedHull) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """``(x, A(x))`` of the generated hull, straight from the per-station profile.

    ``A_i = 2 * b_i * T * R(p_i)``: the sampled polygon area, not its analytic
    limit, so the curve integrates to the mesh volume by construction.
    """
    p = hull.parameters
    nz = p.nz
    xs = tuple(p.length * i / (p.nx - 1) for i in range(p.nx))
    areas = tuple(
        2.0 * b * p.draft * section_area_ratio(pe, nz)
        for b, pe in zip(hull.half_breadths, hull.profile_exponents)
    )
    return xs, areas