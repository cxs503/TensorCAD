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

#: The five Todd (1951) parent forms, ascending: the family's published ladder.
PARENT_CB: tuple[float, ...] = tuple(sorted(_PARENT_FORMS))

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


def _aft_half_integrals(n: float, transom: float) -> tuple[float, float, float, float]:
    """``(int f, int f^2, int xi f, int xi f^2)`` over the aft half ``xi in [0, 1/2]``.

    The transom blends the base taper towards a full-width stern,

    ``f(xi) = (1 - tau) * (1 - |2*xi-1|^n) + tau``,

    which is the only blend that keeps the midship station full (``f(1/2) = 1``) while
    lifting the aft end from 0 to ``tau``: ``f(0) = tau``, and since the immersed
    section area at a station is ``(B/2) f T Cm``, ``tau`` *is* the immersed transom
    area over the midship section area.  Each moment is affine or quadratic in
    ``tau``, so the closed forms survive; ``tau = 0`` returns the legacy integrals
    bit for bit (``(1 - 0) * h + 0 * 1/2`` reproduces ``h`` exactly).
    """
    h = 0.5 * _waterplane_factor(n)
    q = 0.5 * _waterplane_square_factor(n)
    m = _aft_moment_factor(n)
    mq = _waterplane_square_aft_factor(n)
    if transom <= 0.0:
        return h, q, m, mq
    keep = 1.0 - transom
    return (
        keep * h + transom * 0.5,
        keep * keep * q + 2.0 * keep * transom * h + transom * transom * 0.5,
        keep * m + transom * 0.125,
        keep * keep * mq + 2.0 * keep * transom * m + transom * transom * 0.125,
    )


def _solve_single_exponent(cwp_target: float, transom: float) -> float:
    """Longitudinal exponent whose blended waterplane lands on ``cwp_target``.

    The legacy hull has ``gain = 0``, so ``Cb = Cm * Cwp``: once a transom fills the
    stern the exponent has to be re-solved, otherwise the aft-end area raises the
    block coefficient above the requested one.  ``Cwp`` falls monotonically with
    ``n``, so a bracketed bisection is enough.
    """
    def cwp(n: float) -> float:
        h_aft, _q, _m, _mq = _aft_half_integrals(n, transom)
        return h_aft + 0.5 * _waterplane_factor(n)

    lo, hi = 1e-3, 1.0e4
    # ``Cwp`` rises with the exponent (``n/(n+1)`` plus the transom's stern share), so
    # the bracket is lo -> small waterplane, hi -> full waterplane.
    if not cwp(lo) < cwp_target < cwp(hi):
        raise ValueError(
            f"Cwp={cwp_target!r} is out of reach for transom_ratio={transom!r}: "
            f"the blended family spans ({cwp(lo)!r}, {cwp(hi)!r})"
        )
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        if cwp(mid) < cwp_target:
            lo = mid
        else:
            hi = mid
        if hi / lo - 1.0 < 1e-15:
            break
    return math.sqrt(lo * hi)


def _resolve_halves(n_aft: float, n_fore: float, cp: float,
                    transom: float = 0.0) -> tuple[float, float, float, float]:
    """``(Cwp, LCB/L, int f^2, gain)`` for a pair of half exponents.

    ``LCB`` is the centroid of the *volume* distribution ``deck_taper * fullness``,
    not of the deck plan: the two only coincide on the legacy hull.  With
    ``Cm(xi) = Cm_mid (1 + gain (1 - f))`` the volume moment is
    ``Cm_mid [(1+gain) int xi f - gain int xi f^2]`` and the volume integral is
    ``Cm_mid * Cp`` by construction, so the ratio is the expression below.  Getting
    this wrong (using the waterplane moment) shows up as an LCB that converges to
    a nonzero value under refinement, which is how it was found.
    """
    h_aft, q_aft, m_aft, mq_aft = _aft_half_integrals(n_aft, transom)
    cwp = h_aft + 0.5 * _waterplane_factor(n_fore)
    square = q_aft + 0.5 * _waterplane_square_factor(n_fore)
    first = m_aft + _fore_moment_factor(n_fore)
    first_square = mq_aft + _waterplane_square_fore_factor(n_fore)
    if cwp == square:
        raise ValueError(
            f"cwp={cwp!r} equals int(f^2)={square!r}: the section fullness cannot carry "
            "the requested block coefficient"
        )
    gain = (cp - cwp) / (cwp - square)
    moment = (1.0 + gain) * first - gain * first_square
    return cwp, moment / cp, square, gain


def _halves(n_aft: float, n_fore: float,
            transom: float = 0.0) -> tuple[float, float, float]:
    """``(Cwp, int f^2, LCF/L)`` for a hull whose halves carry different exponents."""
    h_aft, q_aft, m_aft, _mq_aft = _aft_half_integrals(n_aft, transom)
    cwp = h_aft + 0.5 * _waterplane_factor(n_fore)
    square = q_aft + 0.5 * _waterplane_square_factor(n_fore)
    first = m_aft + _fore_moment_factor(n_fore)
    return cwp, square, first / cwp


def _solve_halves(cwp_target: float, lcb_target: float, cp: float,
                  transom: float = 0.0) -> tuple[float, float]:
    """Solve the two half exponents that put ``Cwp`` and ``LCB`` on their targets.

    Two unknowns, two residuals, Newton with a numerical Jacobian and a damped step.
    Both residuals are smooth in ``(n_aft, n_fore)`` and monotone in the useful
    neighbourhood (a fuller aft half moves buoyancy aft), so the iteration converges
    in a handful of steps from the symmetric seed.
    """
    n = cwp_target / (1.0 - cwp_target)
    n_aft = n_fore = n
    for _ in range(80):
        cwp, lcb, _square, _gain = _resolve_halves(n_aft, n_fore, cp, transom)
        r1 = cwp - cwp_target
        r2 = lcb - lcb_target
        if abs(r1) < 1e-15 and abs(r2) < 1e-14:
            return n_aft, n_fore
        h = 1e-6
        cwp_a, lcb_a, _, _ = _resolve_halves(n_aft + h, n_fore, cp, transom)
        cwp_b, lcb_b, _, _ = _resolve_halves(n_aft - h, n_fore, cp, transom)
        cwp_c, lcb_c, _, _ = _resolve_halves(n_aft, n_fore + h, cp, transom)
        cwp_d, lcb_d, _, _ = _resolve_halves(n_aft, n_fore - h, cp, transom)
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
    #: Immersed transom area over the midship section area; ``0`` is a pointed stern.
    transom_ratio: float = 0.0

    def side_exponent(self, xi: float) -> float:
        return self.n_aft if xi < 0.5 else self.n_fore

    def deck_taper(self, xi: float) -> float:
        """Deck half-breadth over ``B/2``; ``deck_taper(0) == transom_ratio``."""
        if xi < 0.0 or xi > 1.0:
            return 0.0
        base = 1.0 - math.pow(abs(2.0 * xi - 1.0), self.side_exponent(xi))
        if self.transom_ratio > 0.0 and xi < 0.5:
            return (1.0 - self.transom_ratio) * base + self.transom_ratio
        return base

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
            "transom_ratio": self.transom_ratio,
            "legacy_single_exponent": self.is_legacy,
        }


def _half_vcb_integral(n: float, cm_mid: float, gain: float,
                       transom: float = 0.0) -> float:
    """One half's share of ``int f * Cm(xi) * VCB/T``, exact in ``t = |2*xi-1|``."""

    def integrand(t: float) -> float:
        f = 1.0 - math.pow(t, n)
        if transom > 0.0:
            f = (1.0 - transom) * f + transom
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
    transom = float(getattr(p, "transom_ratio", 0.0) or 0.0)
    if p.cwp is None and p.lcb is None:
        if transom > 0.0:
            # The transom adds stern area, so the single exponent is re-solved to hold
            # Cb = Cm * Cwp instead of letting the filled stern push Cb up.
            n_anchor = _solve_single_exponent(cp, transom)
            cwp, _square, lcf = _halves(n_anchor, n_anchor, transom)
            lcb = _resolve_halves(n_anchor, n_anchor, cp, transom)[1]
            return Series60SectionLaw(
                n_aft=n_anchor, n_fore=n_anchor, cm_mid=cm_mid, gain=0.0, cb=p.cb,
                cwp=cwp, lcb_ratio=lcb, lcf_ratio=lcf,
                vcb_ratio=(m + 1.0) / (2.0 * m + 1.0), is_legacy=False,
                transom_ratio=transom,
            )
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
        if transom > 0.0:
            n_aft = n_fore = _solve_single_exponent(float(p.cwp), transom)
        else:
            n_aft = n_fore = p.cwp / (1.0 - p.cwp)
    else:
        if not 0.30 < p.lcb < 0.70:
            raise ValueError(f"lcb (LCB/L, 0.5 = midship) must lie in (0.30, 0.70), got {p.lcb!r}")
        n_aft, n_fore = _solve_halves(float(p.cwp), float(p.lcb), cp, transom)
    cwp, lcb, square, gain = _resolve_halves(n_aft, n_fore, cp, transom)
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
    vcb = 0.5 * (_half_vcb_integral(n_aft, cm_mid, gain, transom)
                 + _half_vcb_integral(n_fore, cm_mid, gain)) / p.cb
    _cwp_check, _square_check, lcf = _halves(n_aft, n_fore, transom)
    return Series60SectionLaw(
        n_aft=n_aft, n_fore=n_fore, cm_mid=cm_mid, gain=gain, cb=p.cb, cwp=cwp,
        lcb_ratio=lcb, lcf_ratio=lcf, vcb_ratio=vcb, is_legacy=False,
        transom_ratio=transom,
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
    #: Aft-end deck-breadth ratio ``f(0) = tau``: the half-breadth at the stern over
    #: ``B/2``, i.e. how wide the immersed flat transom is.  ``0`` (default) is the
    #: pointed stern every earlier model carried.  This is a *breadth* ratio; the
    #: Holtrop-style area ratio ``A_T / A_m`` is reported separately as
    #: ``transom_area_ratio_achieved``, because the section fullness at the stern is
    #: not the midship fullness in the ``cwp``/``lcb`` family (it is 0.74 of it at
    #: ``Cb = 0.70``), so the two ratios differ by that factor.
    transom_ratio: float = 0.0
    #: Bulb semi-axis ahead of the forward perpendicular over ``L``.  ``0`` = no bulb.
    bulb_length_ratio: float = 0.0
    #: Transverse bulb area at the forward perpendicular over the midship section
    #: area -- the ``A_BT`` of Holtrop & Mennen (1984).
    bulb_area_ratio: float = 0.0
    #: Bulb section centre height above the keel over the draft.
    bulb_centre_height_ratio: float = 0.35

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
        if not math.isfinite(self.transom_ratio) or not 0.0 <= self.transom_ratio < 0.45:
            errors.append(f"transom_ratio must lie in [0, 0.45), got {self.transom_ratio!r}")
        for name, limit in (("bulb_length_ratio", 0.25), ("bulb_area_ratio", 0.5),
                            ("bulb_centre_height_ratio", 2.0)):
            value = float(getattr(self, name))
            if not math.isfinite(value) or not 0.0 <= value <= limit:
                errors.append(f"{name} must lie in [0, {limit}], got {value!r}")
        if (self.bulb_length_ratio > 0.0) != (self.bulb_area_ratio > 0.0):
            errors.append(
                "a bulb needs both bulb_length_ratio and bulb_area_ratio: one long "
                "enough to see, one wide enough to matter"
            )
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
    """Deck half-breadth over ``B/2`` at station ``i``.

    Zero at the bow; ``transom_ratio`` at the stern, so a transom hull carries a real
    (non-degenerate) aft ring and its offsets table ends at ``tau * B/2``.
    """
    if i == p.nx - 1:
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
    if law.transom_ratio > 0.0:
        # The aft ring is no longer degenerate, so the stern needs a real face: a fan
        # over the ring polygon closed by the deck edge (starboard waterline ->
        # keel -> port waterline -> back along the deck). Every edge of the fan is
        # used by exactly one lofted face, which is what keeps the stern watertight.
        ring = rings[0]
        points = [_section_point(p, 0, j, law) for j in range(ring_size)]
        centre = builder.add_vertex((0.0, 0.0, sum(pt[2] for pt in points) / ring_size))
        for j in range(ring_size - 1):
            builder.add_triangle(centre, ring[j + 1], ring[j])
        builder.add_triangle(centre, ring[0], ring[ring_size - 1])
    mesh = meshlib.orient_outward(builder.build()).mesh
    offsets = tuple((p.length * i / (p.nx - 1), 0.5 * p.beam * _longitudinal_taper(p, i, law))
                    for i in range(p.nx))
    return mesh, offsets


def replace_check(check: Check, prefix: str) -> Check:
    """Same check, namespaced, so the bulb's topology is not confused with the hull's."""
    return Check(prefix + check.name, check.passed, check.detail)


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
    #: Closed mesh of the bulb appendage, when one was requested.
    bulb_mesh: Mesh | None = None

    @property
    def bulb(self) -> "Series60Bulb | None":
        return series60_bulb(self.parameters)

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
            "transom": {"transom_ratio": self.parameters.shape.transom_ratio,
                        "definition": "aft-end deck-breadth ratio f(0) = b_T/(B/2)",
                        "immersed_area_m2": transom_area(self.parameters),
                        "immersed_area_ratio_of_mid_section": (
                            transom_area(self.parameters)
                            / (self.parameters.beam * self.parameters.draft
                               * self.parameters.shape.cm_mid))},
            "appendages": {
                "bulb": self.bulb.to_dict() if self.bulb is not None else None,
                "bodies": 2 if self.bulb is not None else 1,
                "note": "the bulb is a separate closed solid, not boolean-unioned with "
                        "the hull; totals add its closed-form displacement",
            },
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
    checks = list(report.checks) + [_sampling_check(p, hydrostatics)]
    bulb_mesh = build_bulb_mesh(p)
    if bulb_mesh is not None:
        checks.extend(replace_check(check, "bulb_") for check in inspect_geometry(bulb_mesh))
        bulb = series60_bulb(p)
        checks.append(Check(
            "bulb_waterline_crossing", True,
            "" if bulb is None or not bulb.pierces_waterline else
            f"the bulb crosses the design waterline: only its submerged part "
            f"({bulb.volume:.6g} m^3 of {bulb.volume_if_submerged:.6g} m^3) is displaced, "
            f"and its waterplane patch ({bulb.waterplane_area:.6g} m^2) is added to the "
            f"hull waterplane",
            severity=WARNING))
    report = merge(checks, [])
    return Series60Model(parameters=p, mesh=mesh, offsets=offsets,
                         validation=report, hydrostatics=hydrostatics,
                         bulb_mesh=bulb_mesh)


def transom_area(p: Series60Parameters, law: Series60SectionLaw | None = None) -> float:
    """Immersed area of the flat stern face [m^2].

    The stern station's own section area: ``(B/2) f(0) T Cm(0)``.  Dividing it by the
    midship section area ``B T Cm_mid`` is the Holtrop-style transom area ratio, which
    differs from ``transom_ratio`` by the ratio of the two section coefficients.
    """
    law = law if law is not None else p.shape
    if law.transom_ratio <= 0.0:
        return 0.0
    # Port and starboard: the section area is ``B f T Cm``, not half of it.
    return p.beam * law.deck_taper(0.0) * p.draft * law.fullness(0.0)


def analytic_series60_hydrostatics(p: Series60Parameters) -> dict[str, float]:
    """Closed forms of the hull that ``p`` actually describes.

    Every quantity comes from the resolved station law, so the same function answers
    for the legacy single-exponent hull, for one whose ``cwp``/``lcb`` were set, and
    for one carrying a transom and a bulb.  The ``bulb_*`` and ``total_*`` keys are the
    appendage accounting: the hull terms never include the bulb, so nothing is counted
    twice.
    """
    law = p.shape
    bulb = series60_bulb(p)
    hull_volume = p.length * p.beam * p.draft * law.cb
    stern_area = transom_area(p, law)
    bulb_volume = bulb.volume if bulb else 0.0
    total_volume = hull_volume + bulb_volume
    hull_lcb = p.length * law.lcb_ratio
    hull_vcb = p.draft * law.vcb_ratio
    moment_x = hull_volume * hull_lcb
    moment_z = hull_volume * hull_vcb
    if bulb is not None and bulb.volume > 0.0:
        moment_x += bulb.volume * (p.length + bulb.centroid_x_ahead)
        moment_z += bulb.volume * bulb.centroid_z_above_keel
    reference: dict[str, float] = {
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
        "transom_area": stern_area,
        "transom_area_ratio": (stern_area / (p.beam * p.draft * law.cm_mid)
                               if law.cm_mid else 0.0),
        "bulb_volume": bulb_volume,
        "bulb_wetted_area": bulb.wetted_area if bulb else 0.0,
        "bulb_fp_area": bulb.fp_area if bulb else 0.0,
        "bulb_fp_area_submerged": bulb.fp_area_submerged if bulb else 0.0,
        "bulb_waterplane_area": bulb.waterplane_area if bulb else 0.0,
        "bulb_centroid_x_ahead": bulb.centroid_x_ahead if bulb else 0.0,
        "total_displacement_volume": total_volume,
        "total_longitudinal_centre_buoyancy": moment_x / total_volume,
        "total_vertical_centre_buoyancy": moment_z / total_volume,
        "total_waterplane_area": p.length * p.beam * law.cwp
        + (bulb.waterplane_area if bulb else 0.0),
    }
    return reference


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

# ---------------------------------------------------------------- bulbous bow

#: Quadrature order for the one bulb integral with no elementary form (its wetted
#: surface).  The integrand is smooth away from the waterline crossing circle, and
#: the order is doubled in the audit to publish the quadrature's own error term.
BULB_QUADRATURE_ORDER = 192


@dataclass(frozen=True)
class Series60Bulb:
    """A bulbous bow ahead of the forward perpendicular, as an exact solid.

    Shape: **half** of a spheroid of revolution lying forward of ``x = L``,

    ``(x - L)^2/a^2 + (y^2 + (z - z_b)^2)/b^2 <= 1,  x >= L``.

    Half, not whole, for a reason that matters to every number below: the part of a
    full bulb behind the FP would sit *inside* the hull whose displacement is already
    integrated, so the two volumes would double count.  The flat face at the FP is an
    ellipse of semi-axes ``b, b``, so its area ``pi b^2`` **is** the ``A_BT``
    ("transverse area of the bulbous bow at the forward perpendicular") that Holtrop
    & Mennen (1984) use, and ``bulb_area_ratio`` is therefore a real input to a
    resistance calculation rather than a shape knob.

    ``volume``, ``fp_area``, ``waterplane_area`` and ``centroid_x_ahead`` are closed
    forms; ``wetted_area`` is Gauss-Legendre and reports its own quadrature error.
    """

    length_ratio: float
    area_ratio: float
    centre_height_ratio: float
    #: Semi-axis along ``x``, ahead of the FP [m].
    semi_axis: float
    #: Radius of the transverse circle [m].
    radius: float
    #: Bulb section centre above the keel [m].
    centre_z: float
    #: ``T - centre_z``: how far the waterline sits above the bulb axis [m].
    waterline_offset: float
    volume: float
    volume_if_submerged: float
    wetted_area: float
    wetted_quadrature_error: float
    fp_area: float
    fp_area_submerged: float
    waterplane_area: float
    centroid_x_ahead: float
    centroid_z_above_keel: float
    pierces_waterline: bool

    def to_dict(self) -> dict[str, Any]:
        return {
            "length_ratio": self.length_ratio,
            "area_ratio": self.area_ratio,
            "centre_height_ratio": self.centre_height_ratio,
            "semi_axis_m": self.semi_axis,
            "radius_m": self.radius,
            "centre_z_above_keel_m": self.centre_z,
            "waterline_offset_m": self.waterline_offset,
            "volume_m3": self.volume,
            "volume_if_fully_submerged_m3": self.volume_if_submerged,
            "wetted_area_m2": self.wetted_area,
            "wetted_quadrature_error_m2": self.wetted_quadrature_error,
            "fp_area_m2": self.fp_area,
            "fp_area_submerged_m2": self.fp_area_submerged,
            "waterplane_area_m2": self.waterplane_area,
            "centroid_x_ahead_m": self.centroid_x_ahead,
            "centroid_z_above_keel_m": self.centroid_z_above_keel,
            "pierces_waterline": self.pierces_waterline,
        }


def bulb_volume(length: float, radius: float, waterline_offset: float) -> float:
    """Submerged volume of the half-spheroid, exact in the waterline offset.

    Integrating the half-ellipse section ``(pi a b / 2) (1 - (z/b)^2)`` from the
    bottom of the bulb up to the waterline gives

    ``V(h) = (pi a b / 2) [ (h + b) - (h^3 + b^3) / (3 b^2) ]``,

    with ``h = T - z_b`` measured from the bulb axis.  ``h = b`` is the fully
    submerged case ``2/3 pi a b^2`` and ``h = 0`` is the half-of-half
    ``pi a b^2 / 3``; both are pinned by tests, because a formula that is wrong only
    at the waterline is exactly the one nobody notices.
    """
    a = _positive(length, "length")
    b = _positive(radius, "radius")
    if waterline_offset >= b:
        return (2.0 / 3.0) * math.pi * a * b * b
    if waterline_offset <= -b:
        return 0.0
    h = waterline_offset
    return 0.5 * math.pi * a * b * ((h + b) - (h ** 3 + b ** 3) / (3.0 * b * b))


def bulb_centroid_x_ahead(length: float, radius: float, waterline_offset: float) -> float:
    """Closed-form ``x`` centroid of the submerged bulb, ahead of the FP.

    Each slice at height ``z`` is a half-ellipse of semi-axes ``a s`` and ``b s``
    (``s = sqrt(1 - (z/b)^2)``) whose own centroid is ``4 a s / (3 pi)`` ahead of the
    FP, so the integral is ``(2/3) a^2 b^2 K(h/b)`` with
    ``K(s) = (1/8)[s (5 - 2 s^2) sqrt(1 - s^2) + 3 asin s] + 3 pi / 16``.
    Fully submerged it collapses to ``3a/8`` -- the hemisphere result -- which is
    what the test pins.
    """
    a = _positive(length, "length")
    b = _positive(radius, "radius")
    volume = bulb_volume(a, b, waterline_offset)
    if volume <= 0.0:
        return float("nan")
    s = min(max(waterline_offset / b, -1.0), 1.0)
    root = math.sqrt(max(0.0, 1.0 - s * s))
    k = 0.125 * (s * (5.0 - 2.0 * s * s) * root + 3.0 * math.asin(s)) + 3.0 * math.pi / 16.0
    return (2.0 / 3.0) * a * a * b * b * k / volume


def _bulb_wetted_area(a: float, b: float, waterline_offset: float,
                      order: int = BULB_QUADRATURE_ORDER) -> float:
    """Wetted surface of the submerged half-spheroid, excluding the FP face.

    ``dA = sqrt(r^2 a^2 + b^4 u^2) du dtheta`` with ``r = b sqrt(1 - u^2)``; the arc
    of the ring that is under water has measure ``2 (pi - acos(h / r))``, clamped
    because a ring wholly below the waterline contributes the full circle and one
    wholly above contributes nothing.
    """
    nodes, weights = gauss_legendre_unit(order)
    total = 0.0
    for u, w in zip(nodes, weights):
        r2 = b * b * (1.0 - u * u)
        if r2 <= 0.0:
            continue
        r = math.sqrt(r2)
        ratio = waterline_offset / r
        if ratio >= 1.0:
            arc = 2.0 * math.pi
        elif ratio <= -1.0:
            continue
        else:
            arc = 2.0 * (math.pi - math.acos(ratio))
        # |dS/du x dS/dtheta| = sqrt((r r')^2 + (a r)^2) = sqrt((b^2 u)^2 + (a r)^2);
        # the axial stretch ``a`` is already inside it, so no extra factor of ``a``.
        ring = math.sqrt(r2 * a * a + b ** 4 * u * u)
        total += w * arc * ring
    return total


def _positive(value: float, name: str) -> float:
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be a finite positive number, got {value!r}")
    return float(value)


def series60_bulb(p: Series60Parameters) -> Series60Bulb | None:
    """Resolve ``p`` into the bulb solid, or ``None`` when no bulb was asked for."""
    if p.bulb_length_ratio <= 0.0 or p.bulb_area_ratio <= 0.0:
        return None
    law = p.shape
    mid_section = p.beam * p.draft * law.cm_mid
    fp_area = p.bulb_area_ratio * mid_section
    b = math.sqrt(fp_area / math.pi)
    a = p.bulb_length_ratio * p.length
    centre_z = p.bulb_centre_height_ratio * p.draft
    offset = p.draft - centre_z
    volume = bulb_volume(a, b, offset)
    wetted = _bulb_wetted_area(a, b, offset)
    wetted_refined = _bulb_wetted_area(a, b, offset, order=2 * BULB_QUADRATURE_ORDER)
    if offset >= b:
        fp_submerged = fp_area
        waterplane = 0.0
    elif offset <= -b:
        fp_submerged = 0.0
        waterplane = 0.0
    else:
        s = offset / b
        # Circular segment of the FP ellipse below the waterline, exact.
        theta = math.acos(max(-1.0, min(1.0, -s)))
        fp_submerged = b * b * (theta - math.sin(theta) * math.cos(theta))
        waterplane = 0.5 * math.pi * a * b * (1.0 - s * s)
    return Series60Bulb(
        length_ratio=p.bulb_length_ratio,
        area_ratio=p.bulb_area_ratio,
        centre_height_ratio=p.bulb_centre_height_ratio,
        semi_axis=a,
        radius=b,
        centre_z=centre_z,
        waterline_offset=offset,
        volume=volume,
        volume_if_submerged=(2.0 / 3.0) * math.pi * a * b * b,
        wetted_area=wetted,
        wetted_quadrature_error=wetted_refined - wetted,
        fp_area=fp_area,
        fp_area_submerged=fp_submerged,
        waterplane_area=waterplane,
        centroid_x_ahead=bulb_centroid_x_ahead(a, b, offset),
        centroid_z_above_keel=centre_z,
        pierces_waterline=volume > 0.0 and offset < b,
    )


def build_bulb_mesh(p: Series60Parameters, *, n_theta: int = 25,
                    n_u: int = 13) -> Mesh | None:
    """Closed triangle mesh of the bulb solid, for export alongside the hull.

    The solid is capped at the FP with a fan, so it is a closed surface on its own:
    the hull is *not* merged with it (that would need a boolean union this module does
    not have), and the metadata says so rather than letting a multi-body STL pass for
    one watertight body.
    """
    bulb = series60_bulb(p)
    if bulb is None:
        return None
    if n_theta < 5 or n_u < 3:
        raise ValueError("n_theta >= 5 and n_u >= 3 are needed for a closed bulb")
    builder = MeshBuilder(units=UNITS, coordinate_frame=COORDINATE_FRAME)
    rings: list[list[int]] = []
    for k in range(n_u - 1):   # the last station is the apex, added separately
        u = k / (n_u - 1)
        radius = bulb.radius * math.sqrt(max(0.0, 1.0 - u * u))
        ring = []
        for t in range(n_theta):
            angle = 2.0 * math.pi * t / n_theta
            ring.append(builder.add_vertex((
                p.length + bulb.semi_axis * u,
                radius * math.sin(angle),
                bulb.centre_z + radius * math.cos(angle),
            )))
        rings.append(ring)
    centre = builder.add_vertex((p.length, 0.0, bulb.centre_z))
    for t in range(n_theta):
        builder.add_triangle(centre, rings[0][(t + 1) % n_theta], rings[0][t])
    for k in range(n_u - 2):
        for t in range(n_theta):  # loft between consecutive non-degenerate rings
            nxt = (t + 1) % n_theta
            builder.add_quad(rings[k][t], rings[k][nxt], rings[k + 1][nxt], rings[k + 1][t])
    apex = builder.add_vertex((p.length + bulb.semi_axis, 0.0, bulb.centre_z))
    for t in range(n_theta):
        builder.add_triangle(apex, rings[-1][t], rings[-1][(t + 1) % n_theta])
    return meshlib.orient_outward(builder.build()).mesh


#: Waterline offset tolerance under which the bulb counts as fully submerged.
BULB_SUBMERGED_TOLERANCE = 1e-12


def series60_appended_hydrostatics(p: Series60Parameters,
                                   hydro: Hydrostatics) -> dict[str, Any]:
    """Total hydrostatics of hull + bulb + transom, each term attributed.

    The hull numbers come from the mesh (``hydro``); the bulb terms from the closed
    forms above.  LCB of the combination is the volume-weighted mean, taken from the
    stern (``x = 0``) so it is directly comparable with the hull's own centre.
    """
    law = p.shape
    bulb = series60_bulb(p)
    hull_volume = hydro.displacement_volume
    total_volume = hull_volume
    moment_x = hull_volume * hydro.longitudinal_centre_buoyancy
    moment_z = hull_volume * hydro.vertical_centre_buoyancy
    wetted = hydro.wetted_surface_area
    waterplane = hydro.waterplane_area
    if bulb is not None and bulb.volume > 0.0:
        total_volume += bulb.volume
        moment_x += bulb.volume * (p.length + bulb.centroid_x_ahead)
        moment_z += bulb.volume * bulb.centroid_z_above_keel
        wetted += bulb.wetted_area
        waterplane += bulb.waterplane_area
    transom_area = 0.0
    if law.transom_ratio > 0.0:
        # Immersed transom area: the stern station's own section area.
        transom_area = (p.beam * law.deck_taper(0.0) * p.draft
                        * law.fullness(0.0))
    return {
        "hull_displacement_volume_m3": hull_volume,
        "bulb_displacement_volume_m3": bulb.volume if bulb else 0.0,
        "total_displacement_volume_m3": total_volume,
        "total_displacement_mass_kg": total_volume * hydro.water_density,
        "total_block_coefficient": total_volume / (p.length * p.beam * p.draft),
        "hull_lcb_m": hydro.longitudinal_centre_buoyancy,
        "total_lcb_m": moment_x / total_volume,
        "total_lcb_over_length": moment_x / total_volume / p.length,
        "total_vcb_above_keel_m": moment_z / total_volume,
        "total_wetted_surface_area_m2": wetted,
        "total_waterplane_area_m2": waterplane,
        "transom_area_m2": transom_area,
        "transom_area_ratio_achieved": (transom_area / (p.beam * p.draft * law.cm_mid)
                                        if law.cm_mid else 0.0),
        "bulb_fp_area_m2": bulb.fp_area if bulb else 0.0,
        "bulb_fp_area_submerged_m2": bulb.fp_area_submerged if bulb else 0.0,
        "bulb_wetted_area_m2": bulb.wetted_area if bulb else 0.0,
        "bulb_pierces_waterline": bool(bulb.pierces_waterline) if bulb else False,
    }
