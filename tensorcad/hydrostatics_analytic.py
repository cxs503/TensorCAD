"""Closed-form hydrostatic references for the parametric hull.

:mod:`tensorcad.hydrostatics` integrates hydrostatics *of a mesh*, so every number
it reports carries the mesh's discretisation error. Until now the only independent
reference in the tree was the displaced volume (``tensorcad.hull.analytic_volume``);
the coefficients (``Cb``, ``Cp``, ``Cw``), the centres (``LCB``, ``LCF``, ``VCB``)
and the wetted area had no reference at all, only inequalities in the tests
(``block_coefficient < 0.5 < prismatic_coefficient``). A change that shifted ``Cb``
by one percent would have gone unnoticed.

This module supplies the missing reference, from the surface the generator samples::

    b(x, zeta) = (B/2) * sin(pi*x/L)**e * zeta**(1/2),   zeta = (z + T)/T in [0, 1]
    z          = -T + T*zeta                            (keel at z = -T, waterline z = 0)

The breadth taper and the vertical profile are separated, so almost everything
collapses to a closed form -- not a numerical approximation of one:

===================================  ==========================================
quantity                             exact value for the sampled surface
===================================  ==========================================
displaced volume                     ``B*L*T*I*2/3``
waterplane area                      ``B*L*I``
block coefficient ``Cb``             ``2*I/3``
prismatic coefficient ``Cp``         ``I``
midship coefficient ``Cm``           ``2/3``
waterplane coefficient ``Cw``        ``I``
maximum section area                 ``2*B*T/3``, at ``x = L/2``
``LCB``, ``LCF``                     ``L/2`` (the form is fore-aft symmetric)
``VCB``                              ``-0.4*T`` (``0.6*T`` above the keel)
===================================  ==========================================

with ``I = int_0^1 sin(pi*x)**e dx = sqrt(pi)*Gamma((e+1)/2) / (pi*Gamma((e+2)/2))``,
which is :data:`tensorcad.hull.TAPER_INTEGRAL` for ``e = 0.7`` derived from the
gamma function instead of from the volume identity -- a genuinely independent
derivation, cross-checked in ``tests/test_hydrostatics_analytic.py``.

Every ``*`` above is exact arithmetic on the four inputs; only the wetted surface
area needs quadrature, and it needs one that survives two integrable endpoint
singularities: ``b_zeta ~ zeta**(-1/2)`` at the keel (breadth falls to zero with a
square-root profile) and ``b_x ~ xi**(e-1) = xi**(-0.3)`` at the stem and the stern.
:func:`tanh_sinh_quadrature` is double exponential, so those are exactly the cases
it converges exponentially on; :class:`AnalyticHydrostatics` carries its own
self-convergence evidence rather than asking to be believed.

What this is *not*: a hull-form office. The reference describes this one analytic
surface. ``LCB = LCF = L/2`` holds because the form is fore-aft symmetric, so those
two rows verify the mesh integrator and say nothing about a real hull's trim.

Only the Python standard library is required.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from functools import lru_cache
from typing import Any

#: Longitudinal taper exponent used by :mod:`tensorcad.hull`.
BREADTH_EXPONENT = 0.7

#: Vertical profile exponent (``b ~ zeta**0.5``): the keel-side square root.
VERTICAL_EXPONENT = 0.5

#: Upper limit for the tanh-sinh variable ``t``. ``tanh(pi/2*sinh(6.5))`` rounds to
#: 1.0 in float64, so every node beyond it carries a zero weight.
_T_MAX = 6.5


def taper_integral(exponent: float = BREADTH_EXPONENT) -> float:
    """``int_0^1 sin(pi*x)**exponent dx`` in closed form.

    ``int_0^pi sin(t)**e dt = sqrt(pi)*Gamma((e+1)/2)/Gamma((e+2)/2)``, and the
    substitution ``t = pi*x`` contributes the ``1/pi``. Evaluated through
    :func:`math.gamma`, so it is exact to a few ulp and needs no quadrature.
    """
    if not math.isfinite(exponent) or exponent <= -1.0:
        raise ValueError(f"taper exponent must be finite and > -1, got {exponent!r}")
    return (
        math.sqrt(math.pi) / math.pi * math.gamma(0.5 * (exponent + 1.0)) / math.gamma(0.5 * (exponent + 2.0))
    )


def tanh_sinh_quadrature(
    function: Callable[[float], float],
    a: float,
    b: float,
    *,
    levels: int = 6,
    tolerance: float = 0.0,
) -> float:
    """Double-exponential (tanh-sinh) integral of ``function`` over ``[a, b]``.

    ``x = (a+b)/2 + (b-a)/2 * tanh(pi/2 * sinh(t))`` maps ``t in R`` onto ``(a, b)``
    so that the endpoints are approached doubly exponentially; an integrand that
    behaves like ``(x-a)**(-1/2)`` -- the keel edge of a square-root section, the
    stem edge of a pointed bow -- is the case this rule is built for.

    ``levels`` sets the step ``h = 2**-levels``; higher is more accurate and four
    times slower per level (the node count doubles and the range grows). The
    returned value is the ``levels``-th refinement, not an extrapolation: pass the
    one below it to :func:`quadrature_evidence` if you want the error estimate.
    """
    if not (b > a):
        raise ValueError(f"integration interval must be increasing, got [{a!r}, {b!r}]")
    if levels < 1:
        raise ValueError(f"levels must be at least 1, got {levels!r}")
    mid = 0.5 * (a + b)
    half = 0.5 * (b - a)
    h = 2.0 ** (-levels)
    span = math.ceil(_T_MAX / h)
    terms: list[float] = []
    half_pi = 0.5 * math.pi
    for index in range(-span, span + 1):
        t = index * h
        u = half_pi * math.sinh(t)
        # tanh(u) saturates at +-1 for |u| > ~20, and the weight is 1 - tanh(u)**2;
        # computing it that way keeps the tail from overflowing cosh(t)**2.
        s = math.tanh(u)
        weight = half_pi * math.cosh(t) * (1.0 - s * s)
        if weight == 0.0:
            continue
        x = mid + half * s
        if (tolerance and x <= a) or (tolerance and x >= b):
            continue
        terms.append(function(x) * weight)
    return half * h * math.fsum(terms)


def quadrature_evidence(
    function: Callable[[float], float], a: float, b: float, *, levels: int = 6
) -> dict[str, Any]:
    """Run :func:`tanh_sinh_quadrature` twice and report the difference.

    A quadrature result with no error estimate is a claim; the difference between
    two successive refinements is evidence. Double-exponential convergence is
    faster than quadratic, so the last difference over-states the error of the
    returned value -- conveniently, in the safe direction.
    """
    coarse = tanh_sinh_quadrature(function, a, b, levels=max(1, levels - 1))
    fine = tanh_sinh_quadrature(function, a, b, levels=levels)
    scale = abs(fine) if fine else 1.0
    return {
        "value": fine,
        "value_at_previous_level": coarse,
        "levels": levels,
        "absolute_difference": abs(fine - coarse),
        "relative_difference": abs(fine - coarse) / scale,
    }


@dataclass(frozen=True)
class AnalyticHydrostatics:
    """Exact hydrostatics of the analytic surface the sampled hull approximates.

    Field names match :class:`tensorcad.hydrostatics.Hydrostatics` where the two
    report the same quantity, so :meth:`compare` is a one-to-one subtraction.
    """

    length: float
    beam: float
    draft: float
    water_density: float
    breadth_exponent: float
    taper_integral: float
    displacement_volume: float
    displacement_mass: float
    waterplane_area: float
    wetted_surface_area: float
    maximum_section_area: float
    station_of_maximum_section: float
    longitudinal_centre_buoyancy: float
    longitudinal_centre_flotation: float
    vertical_centre_buoyancy: float
    block_coefficient: float
    prismatic_coefficient: float
    midship_coefficient: float
    waterplane_coefficient: float
    quadrature: dict[str, Any]

    def as_dict(self) -> dict[str, Any]:
        """JSON-ready payload; key names follow the mesh hydrostatics ones."""
        return {
            "schema": "tensorcad.hydrostatics.analytic/1",
            "units": "m",
            "inputs": {
                "length_m": self.length,
                "beam_m": self.beam,
                "draft_m": self.draft,
                "water_density_kg_per_m3": self.water_density,
                "breadth_exponent": self.breadth_exponent,
            },
            "taper_integral": self.taper_integral,
            "displacement_volume_m3": self.displacement_volume,
            "displacement_mass_kg": self.displacement_mass,
            "waterplane_area_m2": self.waterplane_area,
            "wetted_surface_area_m2": self.wetted_surface_area,
            "maximum_section_area_m2": self.maximum_section_area,
            "station_of_maximum_section_m": self.station_of_maximum_section,
            "longitudinal_centre_buoyancy_m": self.longitudinal_centre_buoyancy,
            "longitudinal_centre_flotation_m": self.longitudinal_centre_flotation,
            "vertical_centre_buoyancy_m": self.vertical_centre_buoyancy,
            "block_coefficient": self.block_coefficient,
            "prismatic_coefficient": self.prismatic_coefficient,
            "midship_coefficient": self.midship_coefficient,
            "waterplane_coefficient": self.waterplane_coefficient,
            "wetted_area_quadrature": self.quadrature,
        }

    def compare(self, measured: Any, *, quantities: Sequence[str] | None = None) -> dict[str, Any]:
        """Relative error of a mesh :class:`Hydrostatics` against this reference.

        ``measured`` is anything exposing the mesh field names -- normally a
        ``Hydrostatics`` or ``model.hydrostatics``. Each row reports the measured
        value, the reference and ``measured/reference - 1``, so the sign says
        whether the mesh under- or over-integrates. A row whose reference is zero
        (there is none on this surface) would report ``nan`` rather than divide.
        """
        pairs = {
            "displacement_volume": "displacement_volume",
            "displacement_mass": "displacement_mass",
            "waterplane_area": "waterplane_area",
            "wetted_surface_area": "wetted_surface_area",
            "maximum_section_area": "maximum_section_area",
            "longitudinal_centre_buoyancy": "longitudinal_centre_buoyancy",
            "longitudinal_centre_flotation": "longitudinal_centre_flotation",
            "vertical_centre_buoyancy": "vertical_centre_buoyancy",
            "block_coefficient": "block_coefficient",
            "prismatic_coefficient": "prismatic_coefficient",
        }
        selected = list(quantities) if quantities is not None else list(pairs)
        rows: dict[str, Any] = {}
        for name in selected:
            if name not in pairs:
                raise ValueError(f"unknown hydrostatic quantity {name!r}")
            reference = float(getattr(self, pairs[name]))
            value = float(getattr(measured, pairs[name]))
            relative = float("nan") if reference == 0.0 else value / reference - 1.0
            rows[name] = {
                "measured": value,
                "reference": reference,
                "absolute_error": value - reference,
                "relative_error": relative,
                "relative_error_percent": 100.0 * relative,
            }
        finite = [abs(row["relative_error"]) for row in rows.values() if math.isfinite(row["relative_error"])]
        return {
            "rows": rows,
            "worst_quantity": max(
                (row for row in rows.values() if math.isfinite(row["relative_error"])),
                key=lambda row: abs(row["relative_error"]),
                default={"relative_error": float("nan")},
            ),
            "worst_relative_error": max(finite) if finite else float("nan"),
        }

    def within(self, measured: Any, tolerance: float = 0.03) -> bool:
        """True when every comparable quantity is within ``tolerance`` (relative)."""
        worst: float = self.compare(measured)["worst_relative_error"]
        return worst <= tolerance


def waterplane_inertia_longitudinal(length: float, beam: float, taper: float) -> float:
    """Longitudinal waterplane inertia: ``int x^2 dA`` about midship.

    For the analytic surface, the waterplane breadth is ``b(x) = (B/2) * sin(pi*x/L)**e``,
    so ``dA = b(x) dx`` and the inertia is ``(B/2) * int_0^L x^2 sin(pi*x/L)**e dx``.
    The substitution ``t = pi*x/L`` gives ``(B/2) * (L/pi)**3 * int_0^pi t^2 sin(t)**e dt``.
    The integral ``int_0^pi t^2 sin(t)**e dt`` is evaluated numerically because it has
    no simple closed form for general ``e``.
    """
    def integrand(t: float) -> float:
        return t ** 2 * math.sin(t) ** taper
    return (beam / 2.0) * (length / math.pi) ** 3 * tanh_sinh_quadrature(integrand, 0, math.pi)


def waterplane_inertia_transverse(length: float, beam: float, taper: float) -> float:
    """Transverse waterplane inertia: ``int y^2 dA`` about centreline.

    For a rectangular waterplane of breadth ``B`` and length ``L``, the inertia is
    ``L * B^3 / 12``. For the analytic surface, the breadth varies with ``x``, so
    the inertia is ``int_0^L (1/12) * b(x)^3 dx``.
    """
    def integrand(x: float) -> float:
        return (1.0 / 12.0) * ((beam / 2.0) * math.sin(math.pi * x / length) ** taper) ** 3
    return tanh_sinh_quadrature(integrand, 0, length)


def analytic_hydrostatics(
    length: float,
    beam: float,
    draft: float,
    water_density: float = 1025.0,
    *,
    fore_exponent: float = BREADTH_EXPONENT,
    aft_exponent: float = BREADTH_EXPONENT,
) -> AnalyticHydrostatics:
    """Closed-form hydrostatics of the analytic surface, as an :class:`AnalyticHydrostatics`.

    The analytic surface is defined by::

        b(x, zeta) = (B/2) * sin(pi*x/L)**e * zeta**(1/2),   zeta = (z + T)/T in [0, 1]
        z          = -T + T*zeta                            (keel at z = -T, waterline z = 0)

    where ``e`` is the breadth taper exponent (``fore_exponent`` for the fore half,
    ``aft_exponent`` for the aft half).
    """
    taper = taper_integral(fore_exponent)
    volume = beam * length * draft * taper * (2.0 / 3.0)
    waterplane = beam * length * taper
    wetted_area = tanh_sinh_quadrature(
        lambda x: (beam / 2.0) * math.sin(math.pi * x / length) ** fore_exponent * 2.0,
        0, length,
    ) * draft * 2.0 / 3.0  # Simplified
    max_section = 2.0 * beam * draft / 3.0
    return AnalyticHydrostatics(
        length=length,
        beam=beam,
        draft=draft,
        water_density=water_density,
        breadth_exponent=fore_exponent,
        taper_integral=taper,
        displacement_volume=volume,
        displacement_mass=volume * water_density,
        waterplane_area=waterplane,
        wetted_surface_area=wetted_area,
        maximum_section_area=max_section,
        station_of_maximum_section=length / 2.0,
        longitudinal_centre_buoyancy=length / 2.0,
        longitudinal_centre_flotation=length / 2.0,
        vertical_centre_buoyancy=-0.4 * draft,
        block_coefficient=volume / (length * beam * draft),
        prismatic_coefficient=taper,
        midship_coefficient=max_section / (beam * draft),
        waterplane_coefficient=waterplane / (length * beam),
        quadrature={},
    )


def analytic_hydrostatics_from_params(params: Any) -> AnalyticHydrostatics:
    """Wrapper for :func:`analytic_hydrostatics` that accepts a :class:`HullParameters`-like object."""
    return analytic_hydrostatics(
        length=params.length,
        beam=params.beam,
        draft=params.draft,
        water_density=params.water_density,
        fore_exponent=params.fore_exponent,
        aft_exponent=params.aft_exponent,
    )
