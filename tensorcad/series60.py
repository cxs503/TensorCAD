"""Series 60 hull family: the systematic displacement-hull benchmark.

The Series 60 (Todd 1951) is the most widely cited parametric hull family in naval
architecture. Five parent forms cover Cb = 0.60, 0.65, 0.70, 0.75, 0.80 with
closed-form hydrostatic coefficients for verification.

Geometry: y(xi, zeta) = (B/2) * (1 - |2*xi-1|^n) * zeta**(1/m)
where n, m are derived from target Cb via published Cp values.

Closed forms: Cp=n/(n+1), Cm=m/(m+1), Cb=Cp*Cm, VCB/T=(m+1)/(2m+1)

**The single longitudinal exponent couples Cp and Cwp, and that coupling is wrong.**
With one exponent ``n`` the deck width and the section area are both proportional to
``1 - |2*xi-1|^n``, so the waterplane coefficient comes out *equal* to the prismatic
coefficient (``Cwp = Cp``) and the longitudinal centre of buoyancy sits at ``L/2``
exactly.  No real hull behaves that way -- a Series 60 parent at ``Cb=0.60`` carries
``Cp=0.625`` with ``Cwp=0.700`` -- so in that respect the family was not a Series 60
replica; ``docs/series60_fidelity.md`` records the measurement (-10.7 % on ``Cwp`` at
``Cb=0.60``; verdict: not a replica below ``Cb`` about 0.78).

``Series60Parameters.cwp`` and ``.lcb`` lift the coupling: the section fullness
``Cm(xi)`` becomes free to vary along the hull while the deck width keeps its own
taper, so ``Cwp`` (deck plan) and ``Cb`` (volume) are independent targets, and each
half of the hull carries its own exponent.  Both default to ``None``, which
reproduces the legacy single-exponent hull bit for bit.
"""
from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field
from typing import Any

from . import mesh as meshlib
from .geometry import Mesh, MeshBuilder
from .hull import COORDINATE_FRAME, UNITS
from .hydrostatics import Hydrostatics, compute_hydrostatics
from .validation import (
    WARNING,
    Check,
    DimensionTarget,
    ValidationReport,
    check_dimensions,
    inspect_geometry,
    merge,
)

_PARENT_FORMS: dict[float, tuple[float, float]] = {
    # (Cp, n) where n = Cp/(1-Cp), m = Cb/(Cp-Cb)
    # Cm = Cb/Cp increases with Cb (U-sections get fuller)
    0.60: (0.625, 5.0 / 3.0),  # Cp=0.625, n=1.667, m=24.0
    0.65: (0.670, 0.670/(1-0.670)),  # Cp=0.670, n=2.030, m=32.5
    0.70: (0.715, 0.715/(1-0.715)),  # Cp=0.715, n=2.509, m=46.7
    0.75: (0.760, 0.760/(1-0.760)),  # Cp=0.760, n=3.167, m=75.0
    0.80: (0.805, 0.805/(1-0.805)),  # Cp=0.805, n=4.128, m=160.0
}

SAMPLING_ERROR_WARNING = 0.02
MAX_STATIONS = 4001
MAX_LEVELS = 513


def _interpolate_exponents(cb: float) -> tuple[float, float]:
    if not math.isfinite(cb) or not 0.55 <= cb <= 0.85:
        raise ValueError(f"Cb must lie in [0.55, 0.85], got {cb!r}")
    sorted_cb = sorted(_PARENT_FORMS.keys())
    if cb <= sorted_cb[0]:
        cp, n = _PARENT_FORMS[sorted_cb[0]]
    elif cb >= sorted_cb[-1]:
        cp, n = _PARENT_FORMS[sorted_cb[-1]]
    else:
        for lo, hi in zip(sorted_cb, sorted_cb[1:]):
            if lo <= cb <= hi:
                frac = (cb - lo) / (hi - lo)
                cp = _PARENT_FORMS[lo][0] + frac * (_PARENT_FORMS[hi][0] - _PARENT_FORMS[lo][0])
                n = _PARENT_FORMS[lo][1] + frac * (_PARENT_FORMS[hi][1] - _PARENT_FORMS[lo][1])
                break
    if cp <= cb:
        raise ValueError(f"Cp ({cp:.6f}) must exceed Cb ({cb:.6f})")
    m = cb / (cp - cb)
    return n, m


def analytic_volume(L, B, T, n, m):
    return L * B * T * n * m / ((n + 1.0) * (m + 1.0))

def analytic_waterplane_area(L, B, n):
    return L * B * n / (n + 1.0)

def analytic_maximum_section_area(B, T, m):
    return B * T * m / (m + 1.0)

def analytic_vcb_above_keel(T, m):
    """Vertical centre of buoyancy above the keel for the ``zeta**(1/m)`` section.

    ``int zeta*y dzeta / int y dzeta`` with ``y ~ zeta**(1/m)`` gives
    ``(m+1)/(2m+1)``, not ``(m+1)/(2(m+2))``: the previous expression was off by
    6.5 % at m=24 (0.481 T against 0.510 T).  Both are inside the "0 < VCB < T"
    sanity band, which is why only an independent recomputation catches it.
    """
    return T * (m + 1.0) / (2.0 * m + 1.0)

def analytic_vcb_z(T, m):
    return analytic_vcb_above_keel(T, m) - T

def analytic_cb(n, m):
    return n * m / ((n + 1.0) * (m + 1.0))

def analytic_cp(n):
    return n / (n + 1.0)

def analytic_cm(m):
    return m / (m + 1.0)

def analytic_cwp(n):
    return n / (n + 1.0)


#: Gauss-Legendre order used for the one integral with no tidy closed form (``VCB``
#: once the section fullness varies along the hull).  The integrand is smooth in the
#: substituted variable ``t = |2*xi-1|``, and doubling the order moves the result by
#: less than 1e-15 relative, so the quadrature is not the error term.
_GL_ORDER = 96

_GL_CACHE: dict[int, tuple[tuple[float, ...], tuple[float, ...]]] = {}


def gauss_legendre_unit(order: int = _GL_ORDER) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Nodes and weights of the ``order``-point Gauss-Legendre rule on ``[0, 1]``.

    Newton iteration on the Legendre polynomial, standard library only: the core
    package declares no dependencies, so ``numpy.polynomial`` is not available here.
    """
    cached = _GL_CACHE.get(order)
    if cached is not None:
        return cached
    if order < 1:
        raise ValueError(f"order must be >= 1, got {order!r}")
    nodes: list[float] = []
    weights: list[float] = []
    for i in range(1, order + 1):
        x = math.cos(math.pi * (i - 0.25) / (order + 0.5))
        deriv = 1.0
        for _ in range(100):
            p_prev, p_cur = 1.0, x
            for k in range(2, order + 1):
                p_prev, p_cur = p_cur, ((2 * k - 1) * x * p_cur - (k - 1) * p_prev) / k
            deriv = order * (x * p_cur - p_prev) / (x * x - 1.0)
            step = p_cur / deriv
            x -= step
            if abs(step) <= 1e-16:
                break
        nodes.append(0.5 * (x + 1.0))
        weights.append(1.0 / ((1.0 - x * x) * deriv * deriv))
    rule = (tuple(nodes), tuple(weights))
    _GL_CACHE[order] = rule
    return rule


def _quad_unit(fn, order: int = _GL_ORDER) -> float:
    """Integrate ``fn`` over ``[0, 1]``."""
    nodes, weights = gauss_legendre_unit(order)
    return math.fsum(w * fn(x) for x, w in zip(nodes, weights))


def _waterplane_factor(n: float) -> float:
    """``int_0^1 (1 - |2*xi-1|^n) dxi`` over one half of the hull."""
    return n / (n + 1.0)


def _waterplane_square_factor(n: float) -> float:
    """``int_0^1 (1 - |2*xi-1|^n)^2 dxi`` over one half of the hull."""
    return 1.0 - 2.0 / (n + 1.0) + 1.0 / (2.0 * n + 1.0)


def _aft_moment_factor(n: float) -> float:
    """``int_0^{1/2} xi (1 - |2*xi-1|^n) dxi``: the aft half's first moment."""
    return 0.25 * (0.5 - 1.0 / (n + 1.0) + 1.0 / (n + 2.0))


def _fore_moment_factor(n: float) -> float:
    """``int_{1/2}^1 xi (1 - |2*xi-1|^n) dxi``: the fore half's first moment."""
    return 0.25 * (1.5 - 1.0 / (n + 1.0) - 1.0 / (n + 2.0))


def _waterplane_square_aft_factor(n: float) -> float:
    """``int_0^{1/2} xi (1 - |2*xi-1|^n)^2 dxi``: aft half's second moment integral."""
    return 0.25 * ((1.0 - 2.0 / (n + 1.0) + 1.0 / (2.0 * n + 1.0))
                   - (0.5 - 2.0 / (n + 2.0) + 1.0 / (2.0 * n + 2.0)))


def _waterplane_square_fore_factor(n: float) -> float:
    """``int_{1/2}^1 xi (1 - |2*xi-1|^n)^2 dxi``: fore half's second moment integral."""
    return 0.25 * ((1.0 - 2.0 / (n + 1.0) + 1.0 / (2.0 * n + 1.0))
                   + (0.5 - 2.0 / (n + 2.0) + 1.0 / (2.0 * n + 2.0)))


def _resolve_halves(n_aft: float, n_fore: float, cp: float) -> tuple[float, float, float, float]:
    """``(Cwp, LCB/L, int f^2, gain)`` for a pair of half exponents.

    ``LCB`` is the centroid of the *volume* distribution ``deck_taper * fullness``,
    not of the deck plan: the two only coincide on the legacy hull.  With
    ``Cm(xi) = Cm_mid (1 + gain (1 - f))`` the volume moment is
    ``Cm_mid [(1+gain) int xi f - gain int xi f^2]`` and the volume integral is
    ``Cm_mid * Cp`` by construction, so the ratio is the expression below.  Getting
    this wrong (using the waterplane moment) shows up as an LCB that converges to
    a nonzero value under refinement, which is how it was found.
    """
    cwp = 0.5 * (_waterplane_factor(n_aft) + _waterplane_factor(n_fore))
    square = 0.5 * (_waterplane_square_factor(n_aft) + _waterplane_square_factor(n_fore))
    first = _aft_moment_factor(n_aft) + _fore_moment_factor(n_fore)
    first_square = (_waterplane_square_aft_factor(n_aft)
                    + _waterplane_square_fore_factor(n_fore))
    gain = (cp - cwp) / (cwp - square)
    moment = (1.0 + gain) * first - gain * first_square
    return cwp, moment / cp, square, gain


def _halves(n_aft: float, n_fore: float) -> tuple[float, float, float]:
    """``(Cwp, int f^2, LCF/L)`` for a hull whose halves carry different exponents."""
    cwp = 0.5 * (_waterplane_factor(n_aft) + _waterplane_factor(n_fore))
    square = 0.5 * (_waterplane_square_factor(n_aft) + _waterplane_square_factor(n_fore))
    first = _aft_moment_factor(n_aft) + _fore_moment_factor(n_fore)
    return cwp, square, first / cwp


def _solve_halves(cwp_target: float, lcb_target: float, cp: float) -> tuple[float, float]:
    """Solve the two half exponents that put ``Cwp`` and ``LCB`` on their targets.

    Two unknowns, two residuals, Newton with a numerical Jacobian and a damped step.
    Both residuals are smooth in ``(n_aft, n_fore)`` and monotone in the useful
    neighbourhood (a fuller aft half moves buoyancy aft), so the iteration converges
    in a handful of steps from the symmetric seed.
    """
    n = cwp_target / (1.0 - cwp_target)
    n_aft = n_fore = n
    for _ in range(80):
        cwp, lcb, _square, _gain = _resolve_halves(n_aft, n_fore, cp)
        r1 = cwp - cwp_target
        r2 = lcb - lcb_target
        if abs(r1) < 1e-15 and abs(r2) < 1e-14:
            return n_aft, n_fore
        h = 1e-6
        cwp_a, lcb_a, _, _ = _resolve_halves(n_aft + h, n_fore, cp)
        cwp_b, lcb_b, _, _ = _resolve_halves(n_aft - h, n_fore, cp)
        cwp_c, lcb_c, _, _ = _resolve_halves(n_aft, n_fore + h, cp)
        cwp_d, lcb_d, _, _ = _resolve_halves(n_aft, n_fore - h, cp)
        j11 = (cwp_a - cwp_b) / (2.0 * h)
        j21 = (lcb_a - lcb_b) / (2.0 * h)
        j12 = (cwp_c - cwp_d) / (2.0 * h)
        j22 = (lcb_c - lcb_d) / (2.0 * h)
        det = j11 * j22 - j12 * j21
        if abs(det) < 1e-18:
            raise ValueError(
                f"cannot separate Cwp={cwp_target!r} from LCB/L={lcb_target!r}: the two "
                "targets need opposite changes in the half exponents"
            )
        dn_aft = (-r1 * j22 + r2 * j12) / det
        dn_fore = (r1 * j21 - r2 * j11) / det
        biggest = max(abs(dn_aft), abs(dn_fore))
        if biggest > 0.5:
            scale = 0.5 / biggest
            dn_aft *= scale
            dn_fore *= scale
        n_aft = min(max(n_aft + dn_aft, 1e-3), 1.0e4)
        n_fore = min(max(n_fore + dn_fore, 1e-3), 1.0e4)
    raise ValueError(
        f"no exponent pair reproduces Cwp={cwp_target!r} with LCB/L={lcb_target!r} "
        f"(best {n_aft!r}, {n_fore!r})"
    )


@dataclass(frozen=True)
class Series60SectionLaw:
    """Resolved station law: deck taper and section fullness held independently.

    ``deck_taper(xi)`` is the half-breadth at the design waterline over ``B/2``, so
    ``Cwp = int deck_taper dxi``.  ``fullness(xi)`` is the section coefficient of the
    *same* station, so ``Cb = int deck_taper * fullness dxi``.  A single shared
    exponent makes the two proportional and ``Cwp = Cp`` falls out; here the fullness
    is free to fall away from midship,

    ``Cm(xi) = Cm_mid * (1 + gain * (1 - deck_taper(xi)))``,

    with ``gain`` chosen so the volume lands on ``Cb`` exactly:

    ``Cb = Cm_mid * (Cwp + gain * (Cwp - int f^2))``.

    ``gain = 0`` is the legacy hull.
    """

    n_aft: float
    n_fore: float
    cm_mid: float
    gain: float
    cb: float
    cwp: float
    lcb_ratio: float
    lcf_ratio: float
    vcb_ratio: float
    is_legacy: bool

    def side_exponent(self, xi: float) -> float:
        return self.n_aft if xi < 0.5 else self.n_fore

    def deck_taper(self, xi: float) -> float:
        if xi <= 0.0 or xi >= 1.0:
            return 0.0
        return 1.0 - math.pow(abs(2.0 * xi - 1.0), self.side_exponent(xi))

    def fullness(self, xi: float) -> float:
        return self.cm_mid * (1.0 + self.gain * (1.0 - self.deck_taper(xi)))

    def exponent(self, xi: float) -> float:
        """Vertical exponent ``m(xi)`` with ``m/(m+1) = fullness(xi)``."""
        q = self.fullness(xi)
        return q / (1.0 - q)

    def to_dict(self) -> dict[str, Any]:
        return {
            "n_aft": self.n_aft,
            "n_fore": self.n_fore,
            "cm_mid": self.cm_mid,
            "section_taper_gain": self.gain,
            "block_coefficient": self.cb,
            "waterplane_coefficient": self.cwp,
            "lcb_over_length": self.lcb_ratio,
            "lcf_over_length": self.lcf_ratio,
            "vcb_over_draft": self.vcb_ratio,
            "legacy_single_exponent": self.is_legacy,
        }


def _half_vcb_integral(n: float, cm_mid: float, gain: float) -> float:
    """One half's share of ``int f * Cm(xi) * VCB/T``, exact in ``t = |2*xi-1|``."""

    def integrand(t: float) -> float:
        f = 1.0 - math.pow(t, n)
        q = cm_mid * (1.0 + gain * (1.0 - f))
        m = q / (1.0 - q)
        return f * q * ((m + 1.0) / (2.0 * m + 1.0))

    return _quad_unit(integrand)


def series60_section_law(p: "Series60Parameters") -> Series60SectionLaw:
    """Resolve ``p`` into the station law shared by the loft and the closed forms.

    ``p.cwp`` and ``p.lcb`` are optional.  Omitting both keeps the legacy
    single-exponent hull, whose ``Cwp`` is whatever ``Cp`` happens to be.
    """
    n_anchor, m = _interpolate_exponents(p.cb)
    cm_mid = m / (m + 1.0)
    cp = p.cb / cm_mid
    if p.cwp is None and p.lcb is None:
        return Series60SectionLaw(
            n_aft=n_anchor, n_fore=n_anchor, cm_mid=cm_mid, gain=0.0, cb=p.cb,
            cwp=_waterplane_factor(n_anchor), lcb_ratio=0.5, lcf_ratio=0.5,
            vcb_ratio=(m + 1.0) / (2.0 * m + 1.0), is_legacy=True,
        )
    if p.cwp is None:
        raise ValueError(
            "lcb needs cwp: LCB is bought by trading waterplane fullness between the aft "
            "and fore halves, so the waterplane coefficient has to be given too"
        )
    if not 0.5 < p.cwp < 0.98:
        raise ValueError(f"cwp must lie in (0.5, 0.98), got {p.cwp!r}")
    if p.lcb is None:
        n_aft = n_fore = p.cwp / (1.0 - p.cwp)
    else:
        if not 0.30 < p.lcb < 0.70:
            raise ValueError(f"lcb (LCB/L, 0.5 = midship) must lie in (0.30, 0.70), got {p.lcb!r}")
        n_aft, n_fore = _solve_halves(float(p.cwp), float(p.lcb), cp)
    cwp, lcb, square, gain = _resolve_halves(n_aft, n_fore, cp)
    if gain <= -1.0:
        raise ValueError(
            f"cwp={p.cwp!r} with cb={p.cb!r} needs a negative section fullness "
            f"(gain={gain:.4f} <= -1): the waterplane is too full for this block "
            f"coefficient at midship coefficient {cm_mid:.4f}"
        )
    if cm_mid * (1.0 + gain) >= 1.0:
        raise ValueError(
            f"cwp={p.cwp!r} with cb={p.cb!r} needs a section coefficient of "
            f"{cm_mid * (1.0 + gain):.4f} at the ends: the waterplane is too fine for "
            "this block coefficient, raise cwp"
        )
    cb = cm_mid * (cwp + gain * (cwp - square))
    if abs(cb - p.cb) > 1e-12 * max(1.0, p.cb):
        raise AssertionError(f"volume closure failed: {cb!r} against {p.cb!r}")
    vcb = 0.5 * (_half_vcb_integral(n_aft, cm_mid, gain)
                 + _half_vcb_integral(n_fore, cm_mid, gain)) / p.cb
    _cwp_check, _square_check, lcf = _halves(n_aft, n_fore)
    return Series60SectionLaw(
        n_aft=n_aft, n_fore=n_fore, cm_mid=cm_mid, gain=gain, cb=p.cb, cwp=cwp,
        lcb_ratio=lcb, lcf_ratio=lcf, vcb_ratio=vcb, is_legacy=False,
    )



@dataclass(frozen=True)
class Series60Parameters:
    length: float = 10.0
    beam: float = 1.5
    draft: float = 0.6
    nx: int = 41
    nz: int = 13
    water_density: float = 1025.0
    cb: float = 0.60
    #: Waterplane coefficient.  ``None`` keeps the legacy hull, where ``Cwp == Cp``
    #: by construction; a value decouples the deck plan from the section area curve.
    cwp: float | None = None
    #: Longitudinal centre of buoyancy as a fraction of ``L`` measured from ``x = 0``
    #: (0.5 = midship).  Real ships carry it aft of midship -- about 0.475 for a
    #: Series 60 parent at ``Cb = 0.60`` -- and needs ``cwp``.
    lcb: float | None = None

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)

    @property
    def exponents(self) -> tuple[float, float]:
        return _interpolate_exponents(self.cb)

    @property
    def shape(self) -> Series60SectionLaw:
        """Resolved station law: deck taper, section fullness, VCB."""
        return series60_section_law(self)

    @property
    def n(self) -> float:
        return self.exponents[0]

    @property
    def m(self) -> float:
        return self.exponents[1]

    def problems(self) -> tuple[list[str], list[str]]:
        errors: list[str] = []
        for name, mx in (("length", 1.0e4), ("beam", 1.0e3), ("draft", 1.0e3)):
            v = float(getattr(self, name))
            if not math.isfinite(v) or v <= 0.0 or v > mx:
                errors.append(f"{name} must be finite and in (0, {mx:g}], got {v!r}")
        if not 3 <= int(self.nx) <= MAX_STATIONS:
            errors.append(f"nx must be in [3, {MAX_STATIONS}], got {self.nx!r}")
        if not 3 <= int(self.nz) <= MAX_LEVELS:
            errors.append(f"nz must be in [3, {MAX_LEVELS}], got {self.nz!r}")
        if not math.isfinite(self.water_density) or not 300.0 <= self.water_density <= 2000.0:
            errors.append(f"water_density must be in [300, 2000], got {self.water_density!r}")
        if not math.isfinite(self.cb) or not 0.55 <= self.cb <= 0.85:
            errors.append(f"cb must lie in [0.55, 0.85], got {self.cb!r}")
        for name in ("cwp", "lcb"):
            value = getattr(self, name)
            if value is not None and not math.isfinite(value):
                errors.append(f"{name} must be finite or None, got {value!r}")
        if not errors:
            try:
                series60_section_law(self)
            except ValueError as exc:
                errors.append(str(exc))
        warnings: list[str] = []
        if not errors:
            lb = self.length / self.beam
            bt = self.beam / self.draft
            if not 1.5 <= lb <= 15.0:
                warnings.append(f"length_to_beam = {lb:.3g} outside [1.5, 15]")
            if not 1.0 <= bt <= 12.0:
                warnings.append(f"beam_to_draft = {bt:.3g} outside [1.0, 12]")
        return errors, warnings

    def validate(self) -> None:
        errors, _ = self.problems()
        if errors:
            raise ValueError("; ".join(errors))

    @property
    def grid_is_symmetric(self) -> bool:
        return self.nx % 2 == 1


def _longitudinal_taper(p: Series60Parameters, i: int,
                        law: Series60SectionLaw | None = None) -> float:
    """Deck half-breadth over ``B/2`` at station ``i`` (0 at the two ends)."""
    if i == 0 or i == p.nx - 1:
        return 0.0
    law = law if law is not None else series60_section_law(p)
    return law.deck_taper(i / (p.nx - 1))


def _vertical_taper_at(nz: int, level: int, m: float) -> float:
    """``zeta ** (1/m)`` for a station whose own fullness fixes ``m``."""
    if level >= nz - 1:
        return 1.0
    if level <= 0:
        return 0.0
    return math.pow(level / (nz - 1), 1.0 / m)


def _vertical_taper(p: Series60Parameters, level: int) -> float:
    """Section fullness at ``level`` (0 at the keel, ``nz-1`` at the waterline).

    ``zeta ** (1/m)``, not ``zeta ** m``.  The closed forms in this module --
    ``Cm = m/(m+1)`` and ``Cb = Cp*Cm``, which are what pins the exponents to the
    published Series 60 parent forms -- follow from ``int_0^1 zeta**(1/m) dzeta =
    m/(m+1)``.  With the exponent inverted the lofted sections collapse to a sliver
    over the bottom 90 % of the draft and the mesh integrates 21 times less volume
    than its own analytic reference; the error does not shrink under refinement
    (measured -95.8 % at 241x101), which is how it was found during the geometry
    audit (docs/series60.md records the measurement).
    """
    if level >= p.nz - 1:
        return 1.0
    if level <= 0:
        return 0.0
    return math.pow(level / (p.nz - 1), 1.0 / p.m)


def _section_point(p: Series60Parameters, i: int, j: int,
                   law: Series60SectionLaw | None = None) -> tuple[float, float, float]:
    law = law if law is not None else series60_section_law(p)
    xi = i / (p.nx - 1)
    x = p.length * xi
    taper = law.deck_taper(xi)
    if law.is_legacy:
        def vertical(level: int) -> float:
            return _vertical_taper(p, level)  # legacy path, and its negative control
    else:
        station_m = law.exponent(xi)

        def vertical(level: int) -> float:
            return _vertical_taper_at(p.nz, level, station_m)
    sign = 1.0
    if j <= p.nz - 1:
        level = p.nz - 1 - j
    else:
        level = j - (p.nz - 1)
        sign = -1.0
    b = 0.5 * p.beam * taper * vertical(level)
    z = -p.draft * (1.0 - level / (p.nz - 1))
    return (x, sign * b, z)


def build_series60_mesh(p: Series60Parameters) -> tuple[Mesh, tuple[tuple[float, float], ...]]:
    p.validate()
    law = series60_section_law(p)
    ring_size = 2 * p.nz - 1
    builder = MeshBuilder(units=UNITS, coordinate_frame=COORDINATE_FRAME)
    rings: list[list[int]] = []
    for i in range(p.nx):
        rings.append([builder.add_vertex(_section_point(p, i, j, law)) for j in range(ring_size)])
    for i in range(p.nx - 1):
        for j in range(ring_size - 1):
            builder.add_quad(rings[i][j], rings[i + 1][j], rings[i + 1][j + 1], rings[i][j + 1])
        builder.add_quad(rings[i][ring_size - 1], rings[i + 1][ring_size - 1],
                         rings[i + 1][0], rings[i][0])
    mesh = meshlib.orient_outward(builder.build()).mesh
    offsets = tuple((p.length * i / (p.nx - 1), 0.5 * p.beam * _longitudinal_taper(p, i, law))
                    for i in range(p.nx))
    return mesh, offsets


def _dimension_target(p: Series60Parameters) -> DimensionTarget:
    return DimensionTarget(length=p.length, beam=p.beam, draft=p.draft)


def _sampling_check(p: Series60Parameters, hydro: Hydrostatics) -> Check:
    ref_vol = analytic_series60_hydrostatics(p)["displacement_volume"]
    rel = (hydro.displacement_volume - ref_vol) / ref_vol if ref_vol > 0 else float("nan")
    ok = abs(rel) <= SAMPLING_ERROR_WARNING
    return Check("series60_volume_convergence", ok,
                 f"volume {hydro.displacement_volume:.6f} vs analytic {ref_vol:.6f} = {rel:+.4%}",
                 severity=WARNING)


@dataclass(frozen=True)
class Series60Model:
    parameters: Series60Parameters
    mesh: Mesh
    offsets: tuple[tuple[float, float], ...]
    validation: ValidationReport = field(default_factory=ValidationReport)
    hydrostatics: Hydrostatics | None = None

    @property
    def vertices(self): return self.mesh.vertices
    @property
    def faces(self): return self.mesh.faces
    @property
    def units(self): return self.mesh.units
    @property
    def coordinate_frame(self): return self.mesh.coordinate_frame

    def metadata(self, *, with_discretisation: bool = False) -> dict[str, Any]:
        n, m = self.parameters.exponents
        data: dict[str, Any] = {
            "schema_version": "2.0", "model_type": "series_60_parametric_displacement_hull",
            "parameters": self.parameters.as_dict(), "exponents": {"n": n, "m": m},
            "units": self.units, "coordinate_frame": self.coordinate_frame,
            "station_law": self.parameters.shape.to_dict(),
            "generator": {"kind": "series_60_lofted_surface",
                          "longitudinal_taper": f"1 - |2*xi - 1|^{n:.4f}"
                                                if self.parameters.shape.is_legacy
                                                else "1 - |2*xi - 1|^n(side)",
                          "vertical_taper": f"zeta^{m:.4f}"
                                            if self.parameters.shape.is_legacy
                                            else "zeta^(1/m(xi))",
                          "topology": "one closed section ring per station, capped at WL",
                          "backend": "tensorcad.geometry (standard library)"},
            "mesh": self.mesh.metadata(), "offset_count": len(self.offsets),
            "validation": self.validation.to_dict(),
            "engineering_status": "geometry and hydrostatics only",
        }
        if self.hydrostatics is not None:
            data["hydrostatics"] = self.hydrostatics.to_dict()
        return data


def generate_series60(parameters: Series60Parameters | None = None,
                      **overrides: Any) -> Series60Model:
    if parameters is None:
        p = Series60Parameters(**overrides)
    elif overrides:
        p = Series60Parameters(**{**parameters.as_dict(), **overrides})
    else:
        p = parameters
    p.validate()
    mesh, offsets = build_series60_mesh(p)
    report = merge(inspect_geometry(mesh), check_dimensions(mesh, _dimension_target(p)))
    if report.errors:
        raise ValueError("generated Series 60 hull failed validation: " + "; ".join(report.errors))
    hydrostatics = compute_hydrostatics(mesh, 0.0, water_density=p.water_density, closed=True)
    report = merge(list(report.checks), [_sampling_check(p, hydrostatics)])
    return Series60Model(parameters=p, mesh=mesh, offsets=offsets,
                         validation=report, hydrostatics=hydrostatics)


def analytic_series60_hydrostatics(p: Series60Parameters) -> dict[str, float]:
    """Closed forms of the hull that ``p`` actually describes.

    Every quantity comes from the resolved station law, so the same function answers
    for the legacy single-exponent hull and for one whose ``cwp``/``lcb`` were set.
    """
    law = p.shape
    return {
        "displacement_volume": p.length * p.beam * p.draft * law.cb,
        "waterplane_area": p.length * p.beam * law.cwp,
        "maximum_section_area": p.beam * p.draft * law.cm_mid,
        "vertical_centre_buoyancy": p.draft * law.vcb_ratio,
        "vertical_centre_buoyancy_z": p.draft * (law.vcb_ratio - 1.0),
        "block_coefficient": law.cb,
        "prismatic_coefficient": law.cb / law.cm_mid,
        "midship_coefficient": law.cm_mid,
        "waterplane_coefficient": law.cwp,
        "longitudinal_centre_buoyancy": p.length * law.lcb_ratio,
        "longitudinal_centre_flotation": p.length * law.lcf_ratio,
    }


def series60_analytic_relative_errors(p: Series60Parameters, hydro: Hydrostatics) -> dict[str, float]:
    ref = analytic_series60_hydrostatics(p)
    errors: dict[str, float] = {}
    for key in ("displacement_volume", "waterplane_area", "maximum_section_area",
                "vertical_centre_buoyancy", "longitudinal_centre_buoyancy",
                "longitudinal_centre_flotation", "block_coefficient",
                "prismatic_coefficient", "midship_coefficient", "waterplane_coefficient"):
        rv = ref.get(key)
        mv = getattr(hydro, key, None) if hydro else None
        if rv is not None and mv is not None and abs(rv) > 1e-15:
            errors[key] = (mv - rv) / rv
        else:
            errors[key] = float("nan")
    return errors


def estimate_series60_discretisation_error(p: Series60Parameters,
                                           max_triangles: int = 2_000_000) -> dict[str, Any]:
    ref = analytic_series60_hydrostatics(p)
    ref_vol, ref_awp = ref["displacement_volume"], ref["waterplane_area"]
    results: list[dict[str, Any]] = []
    nx, nz = p.nx, p.nz
    while True:
        pp = Series60Parameters(length=p.length, beam=p.beam, draft=p.draft,
                                nx=nx, nz=nz, water_density=p.water_density, cb=p.cb,
                                cwp=p.cwp, lcb=p.lcb)
        model = generate_series60(pp)
        tc = len(model.mesh.faces)
        ve = (model.hydrostatics.displacement_volume - ref_vol) / ref_vol if ref_vol else float("nan")
        ae = (model.hydrostatics.waterplane_area - ref_awp) / ref_awp if ref_awp else float("nan")
        results.append({"nx": nx, "nz": nz, "triangles": tc,
                        "volume_m3": model.hydrostatics.displacement_volume,
                        "volume_error": ve, "waterplane_area_m2": model.hydrostatics.waterplane_area,
                        "waterplane_error": ae})
        if tc > max_triangles or nx > MAX_STATIONS // 2:
            break
        nx = 2 * nx - 1
        nz = 2 * nz - 1
    return {"reference": ref, "convergence": results}
