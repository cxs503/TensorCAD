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


def _resolve_parameters(parameters: Any) -> tuple[float, float, float, float]:
    """Accept ``HullParameters``, a ``(L, B, T)`` triple, or ``None`` for defaults."""
    if parameters is None:
        from .hull import HullParameters

        parameters = HullParameters()
    if hasattr(parameters, "length") and hasattr(parameters, "beam") and hasattr(parameters, "draft"):
        density = float(getattr(parameters, "water_density", 1025.0))
        return (
            float(parameters.length),
            float(parameters.beam),
            float(parameters.draft),
            density,
        )
    if isinstance(parameters, (tuple, list)) and len(parameters) == 3:
        return (float(parameters[0]), float(parameters[1]), float(parameters[2]), 1025.0)
    raise TypeError(
        "analytic_hydrostatics takes HullParameters, an (L, B, T) triple or None, "
        f"got {type(parameters).__name__}"
    )


def _check_positive(name: str, value: float) -> float:
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be a finite positive number, got {value!r}")
    return float(value)


@lru_cache(maxsize=64)
def _wetted_area(
    length: float, beam: float, draft: float, exponent: float, levels: int
) -> tuple[float, float]:
    """Wetted surface of the half-body, doubled, plus its self-convergence error.

    Parametrise the flank as ``r(x, zeta) = (x, b(x, zeta), -T + T*zeta)``; the area
    element is ``sqrt(T**2*(1 + b_x**2) + b_zeta**2) dx dzeta``. Substituting
    ``zeta = u**2`` -- which removes the ``zeta**(-1/2)`` keel singularity of
    ``b_zeta`` into a smooth integrand -- leaves

        ``int_0^1 sqrt(4*T**2*u**2*(1 + (gamma*u)**2) + beta**2) du``

    with ``beta = (B/2)*sin(pi*x/L)**e`` and ``gamma = d beta / dx``. The outer
    integral over ``x`` still carries the ``xi**(-0.3)`` stem/stern edge, which is
    what the double-exponential rule is for.
    """
    half_beam = 0.5 * beam
    inner_levels = max(1, levels - 2)

    def flank(x: float) -> float:
        s = math.sin(math.pi * x / length)
        if s <= 0.0:  # only reachable at the endpoints, which tanh-sinh never hits
            return 0.0
        beta = half_beam * s**exponent
        gamma = (
            half_beam * (math.pi * exponent / length) * math.cos(math.pi * x / length) * s ** (exponent - 1.0)
        )
        four_t_squared = 4.0 * draft * draft
        beta_squared = beta * beta

        def inner(u: float) -> float:
            slope = gamma * u
            return math.sqrt(four_t_squared * u * u * (1.0 + slope * slope) + beta_squared)

        return tanh_sinh_quadrature(inner, 0.0, 1.0, levels=inner_levels)

    coarse = 2.0 * tanh_sinh_quadrature(flank, 0.0, length, levels=max(1, levels - 1))
    fine = 2.0 * tanh_sinh_quadrature(flank, 0.0, length, levels=levels)
    return fine, abs(fine - coarse)


def analytic_hydrostatics(
    parameters: Any = None,
    *,
    breadth_exponent: float = BREADTH_EXPONENT,
    quadrature_levels: int = 6,
    water_density: float | None = None,
) -> AnalyticHydrostatics:
    """Exact hydrostatics of the analytic surface the sampled hull approximates.

    ``parameters`` is a :class:`tensorcad.hull.HullParameters` (or anything with
    ``length``/``beam``/``draft``), an ``(L, B, T)`` triple, or ``None`` for the
    defaults. ``nx``/``nz`` are ignored: the reference describes the surface, not
    the sampling of it, which is the whole point of having it. ``water_density``
    overrides the density that would otherwise come from ``parameters``.
    """
    length, beam, draft, density = _resolve_parameters(parameters)
    if water_density is not None:
        density = float(water_density)
    _check_positive("length", length)
    _check_positive("beam", beam)
    _check_positive("draft", draft)
    _check_positive("water_density", density)

    taper = taper_integral(breadth_exponent)
    volume = beam * length * draft * taper * (2.0 / 3.0)
    waterplane = beam * length * taper
    section_area = (2.0 / 3.0) * beam * draft
    wetted, wetted_difference = _wetted_area(length, beam, draft, breadth_exponent, int(quadrature_levels))
    return AnalyticHydrostatics(
        length=length,
        beam=beam,
        draft=draft,
        water_density=density,
        breadth_exponent=breadth_exponent,
        taper_integral=taper,
        displacement_volume=volume,
        displacement_mass=volume * density,
        waterplane_area=waterplane,
        wetted_surface_area=wetted,
        maximum_section_area=section_area,
        station_of_maximum_section=0.5 * length,
        # sin(pi*x/L) is symmetric about L/2 and the section shape does not depend
        # on x, so both centres sit exactly at midship -- no quadrature needed.
        longitudinal_centre_buoyancy=0.5 * length,
        longitudinal_centre_flotation=0.5 * length,
        # z-weighted mean of zeta**(1/2) on [-T, 0]: (-4T/15) / (2/3) = -0.4*T.
        vertical_centre_buoyancy=-0.4 * draft,
        block_coefficient=taper * (2.0 / 3.0),
        prismatic_coefficient=taper,
        midship_coefficient=2.0 / 3.0,
        waterplane_coefficient=taper,
        quadrature={
            "wetted_area_levels": int(quadrature_levels),
            "wetted_area_absolute_difference": wetted_difference,
            "wetted_area_relative_difference": wetted_difference / wetted if wetted else float("nan"),
            "method": "tanh-sinh, nested (zeta = u**2 substitution for the keel edge)",
        },
    )


def waterplane_inertia_transverse(
    length: float, beam: float, *, breadth_exponent: float = BREADTH_EXPONENT
) -> float:
    """``I_T = int y**2 dA`` of the reference waterplane, in closed form.

    The waterplane of the analytic surface is ``|y| <= (B/2) sin(pi x/L)**e``, so

        I_T = int_0^L (2/3) (B/2 sin(pi x/L)**e)**3 dx
            = (B**3 L / 12) * int_0^1 sin(pi u)**(3e) du
            = (B**3 L / 12) * taper_integral(3e),

    and ``taper_integral`` is a Beta integral evaluated through ``math.gamma``, so this
    needs no quadrature at all. Note the exponent is ``3e``, not ``e``: the *half
    breadth* tapers as ``sin**e`` and the moment integrates its cube. Getting that
    wrong is a factor-of-two-sized error in ``BM`` that no mesh convergence study would
    reveal, because both the mesh and the reference would be wrong in the same way --
    which is why this is derived here and cross-checked against the mesh in
    ``tests/test_stability.py`` rather than being copied from the mesh integrator.
    """
    _check_positive("length", length)
    _check_positive("beam", beam)
    return beam**3 * length / 12.0 * taper_integral(3.0 * breadth_exponent)


def waterplane_inertia_longitudinal(
    length: float,
    beam: float,
    *,
    breadth_exponent: float = BREADTH_EXPONENT,
    quadrature_levels: int = 8,
) -> tuple[float, float]:
    """``I_L`` about midship, plus the quadrature's self-convergence evidence.

    Substituting ``u = x/L - 1/2`` (so ``sin(pi x/L) = cos(pi u)``) and using the
    fore-aft symmetry of the reference form,

        I_L = int_0^L (x - L/2)**2 * 2 y dx = 2 B L**3 int_0^{1/2} u**2 cos(pi u)**e du.

    ``e = 0.7 < 1`` leaves a derivative singularity at ``u = 1/2`` (the waterline closes
    to a point there) which is integrable but which no fixed-order rule handles well;
    :func:`tanh_sinh_quadrature` is built for exactly this and its two successive
    refinements are returned as the second element. Unlike ``I_T`` there is no elementary
    closed form, so the returned error estimate -- not a promise -- is what the tests
    hold this to.
    """
    _check_positive("length", length)
    _check_positive("beam", beam)
    if quadrature_levels < 2:
        raise ValueError(f"quadrature_levels must be at least 2, got {quadrature_levels!r}")

    def integrand(u: float) -> float:
        return u * u * math.cos(math.pi * u) ** breadth_exponent

    fine = tanh_sinh_quadrature(integrand, 0.0, 0.5, levels=int(quadrature_levels))
    coarse = tanh_sinh_quadrature(integrand, 0.0, 0.5, levels=int(quadrature_levels) - 1)
    scale = 2.0 * beam * length**3
    return scale * fine, abs(fine - coarse) * scale
