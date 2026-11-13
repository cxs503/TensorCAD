"""Wigley parabolic hull form: the hull every solver benchmark starts from.

``tensorcad.hull`` lofts one family -- ``sin(pi x/L)**0.7 * sqrt(z/T)`` -- which
is fore-aft symmetric, has no transom, and (the reason this module exists) is
awkward to quote against literature: the ``sqrt`` keel taper is not the surface
anybody publishes resistance or diffraction data on. The Wigley form is:

    ``y(x, w) = (B/2) * (1 - |2x/L - 1|**n) * (1 - (w/T)**m)``

with ``w`` the depth below the design waterline (``w = 0`` at the waterline,
``w = T`` at the keel) and the classic Wigley being ``n = m = 2``. It is the
standard geometry source for wave-resistance and seakeeping comparisons, and its
hydrostatic coefficients are exact rationals:

    ``Cb  = n m / ((n+1)(m+1))``      ``4/9`` at ``n = m = 2``
    ``Cp  = Cwp = n / (n+1)``         ``2/3``
    ``Cm  = m / (m+1)``               ``2/3``
    ``VCB = T (m+1) / (2 (m+2))``     above the keel, ``5T/8``
    ``Awp = L B n / (n+1)``           ``2LB/3``
    ``V   = Cb L B T``                ``4LBT/9``

Everything closed-form here is derived by hand in the docstrings and re-derived
numerically by the mesh integrator, so the two are independent: the module is
simultaneously a hull-form generator and a verification case. ``examples/wigley``
refines the grid and shows the integrated quantities converge to the closed
forms; ``tests/test_wigley.py`` pins the numbers.

Only the Python standard library is required.
"""
from __future__ import annotations

import math
from collections.abc import Callable
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

#: Quantities both the sampled mesh and the closed forms report, in table order.
ANALYTIC_COMPARABLE_QUANTITIES: tuple[str, ...] = (
    "displacement_volume",
    "waterplane_area",
    "wetted_surface_area",
    "maximum_section_area",
    "longitudinal_centre_buoyancy",
    "longitudinal_centre_flotation",
    "vertical_centre_buoyancy",
    "block_coefficient",
    "prismatic_coefficient",
)

#: The geometry resistance/diffraction benchmarks actually quote: L/B = 10,
#: B/T = 1.6 (so T = L/16). Recorded here so a case can be reproduced by name.
REFERENCE_CASE: dict[str, float | str] = {
    "name": "Wigley (ITTC reference proportions)",
    "length": 10.0,
    "beam": 1.0,
    "draft": 0.625,
    "length_to_beam": 10.0,
    "beam_to_draft": 1.6,
    "longitudinal_exponent": 2.0,
    "vertical_exponent": 2.0,
}

#: Above this relative deviation from the closed-form volume the sampled grid is
#: reported as too coarse to trust (same policy as ``tensorcad.hull``).
SAMPLING_ERROR_WARNING = 0.02

MAX_STATIONS = 4001
MAX_LEVELS = 513


@dataclass(frozen=True)
class WigleyParameters:
    """Principal dimensions, sampling resolution, shape exponents and fluid.

    Defaults are the ITTC reference proportions (``L/B = 10``, ``B/T = 1.6``) on
    the grid resolution ``tensorcad.hull`` uses for its default case, so the two
    families can be compared at equal sampling cost.

    Construction never fails; call :meth:`validate` (the generator does) to get
    every problem at once. Exponents ``>= 1`` keep the waterline corner and the
    keel angle finite; ``n = m = 2`` is the classic Wigley, smaller values fatten
    the shoulders, larger values sharpen them towards a wedge.
    """

    length: float = 10.0
    beam: float = 1.0
    draft: float = 0.625
    nx: int = 41
    nz: int = 13
    water_density: float = 1025.0
    longitudinal_exponent: float = 2.0
    vertical_exponent: float = 2.0

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)

    def problems(self) -> tuple[list[str], list[str]]:
        """Return ``(errors, warnings)``; warnings are advisory, never fatal."""
        errors: list[str] = []
        for name, maximum in (("length", 1.0e4), ("beam", 1.0e3), ("draft", 1.0e3)):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0 or value > maximum:
                errors.append(f"{name} must be finite and in (0, {maximum:g}], got {value!r}")
        if not 3 <= int(self.nx) <= MAX_STATIONS:
            errors.append(f"nx must be in [3, {MAX_STATIONS}], got {self.nx!r}")
        if not 3 <= int(self.nz) <= MAX_LEVELS:
            errors.append(f"nz must be in [3, {MAX_LEVELS}], got {self.nz!r}")
        if not math.isfinite(self.water_density) or not 300.0 <= self.water_density <= 2000.0:
            errors.append(f"water_density must be in [300, 2000] kg/m^3, got {self.water_density!r}")
        for name in ("longitudinal_exponent", "vertical_exponent"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value < 1.0 or value > 8.0:
                errors.append(f"{name} must be in [1, 8], got {value!r}")

        warnings: list[str] = []
        if not errors:
            for name, value, low, high in (
                ("length_to_beam", self.length / self.beam, 1.5, 15.0),
                ("beam_to_draft", self.beam / self.draft, 1.0, 12.0),
            ):
                if not low <= value <= high:
                    warnings.append(f"{name} = {value:.3g} is outside the usual [{low:g}, {high:g}]")
        return errors, warnings

    def validate(self) -> None:
        """Raise :class:`ValueError` naming every problem found, not just the first."""
        errors, _ = self.problems()
        if errors:
            raise ValueError("; ".join(errors))

    @property
    def grid_is_symmetric(self) -> bool:
        """True when a station sits exactly on midship (odd ``nx``).

        With an even ``nx`` no station lands on the widest section, so the realised
        beam falls short of the requested beam. The deficit is reported, never
        silently absorbed (see :func:`check_wigley_consistency`).
        """
        return self.nx % 2 == 1


def reference_parameters(**overrides: Any) -> WigleyParameters:
    """The ITTC reference proportions as a parameter set (overridable)."""
    base: dict[str, Any] = {
        "length": float(REFERENCE_CASE["length"]),
        "beam": float(REFERENCE_CASE["beam"]),
        "draft": float(REFERENCE_CASE["draft"]),
    }
    base.update(overrides)
    return WigleyParameters(**base)


def longitudinal_taper(parameters: WigleyParameters, index: int) -> float:
    """``1 - |2x/L - 1|**n`` at station ``index``, forced to zero at both ends.

    ``abs(u)**n`` with ``u = +-1`` is exactly zero, but the endpoint branches are
    kept so the stem and the post close on the centreline by construction rather
    than by floating-point luck.
    """
    if index == 0 or index == parameters.nx - 1:
        return 0.0
    u = 2.0 * index / (parameters.nx - 1) - 1.0
    return 1.0 - math.pow(abs(u), parameters.longitudinal_exponent)


def vertical_taper(parameters: WigleyParameters, level: int) -> float:
    """``1 - (w/T)**m`` at vertical ``level`` (0 = keel, ``nz-1`` = waterline).

    ``level`` counts up from the keel, so the depth fraction below the waterline is
    ``w/T = 1 - level/(nz-1)``. The two end levels are returned exactly (``0`` at
    the keel, ``1`` at the waterline) so the keel ring collapses cleanly.
    """
    if level >= parameters.nz - 1:
        return 1.0
    if level <= 0:
        return 0.0
    depth_fraction = 1.0 - level / (parameters.nz - 1)
    return 1.0 - math.pow(depth_fraction, parameters.vertical_exponent)


def half_breadth(parameters: WigleyParameters, x: float, depth: float) -> float:
    """Analytic half-breadth at ``(x, depth)``; ``depth`` is measured down from the WL.

    This *is* the surface the mesh samples rather than an approximation of it, so
    it is the right thing to test station offsets against.
    """
    u = min(abs(2.0 * x / parameters.length - 1.0), 1.0)
    taper = max(1.0 - math.pow(u, parameters.longitudinal_exponent), 0.0)
    if depth <= 0.0:
        section = 1.0
    elif depth >= parameters.draft:
        section = 0.0
    else:
        section = max(1.0 - math.pow(depth / parameters.draft, parameters.vertical_exponent), 0.0)
    return 0.5 * parameters.beam * taper * section


def _section_point(parameters: WigleyParameters, i: int, j: int) -> tuple[float, float, float]:
    """Ring position ``j`` at station ``i`` (same walk as the ``hull`` family).

    ``j`` goes down the starboard flank from the design waterline to the keel, then
    up the port flank. The keel point is shared by both flanks, which is what makes
    the lofted surface closed without a welding pass.
    """
    p = parameters
    x = p.length * i / (p.nx - 1)
    keelward = j <= p.nz - 1
    level = (p.nz - 1 - j) if keelward else (j - (p.nz - 1))
    # Written as (1 - fraction) so the waterline is exactly z = 0 (see hull.py).
    z = -p.draft * (1.0 - level / (p.nz - 1))
    breadth = 0.5 * p.beam * longitudinal_taper(p, i) * vertical_taper(p, level)
    return (x, breadth if keelward else -breadth, z)


# ---------------------------------------------------------------- closed forms


def analytic_volume(
    length: float,
    beam: float,
    draft: float,
    longitudinal_exponent: float = 2.0,
    vertical_exponent: float = 2.0,
) -> float:
    """Exact displaced volume ``L B T n m / ((n+1)(m+1))``.

    Derivation: ``V = 2 int_0^L int_0^T (B/2) taper(x) section(w) dw dx`` and the
    two factors separate -- ``int_0^L taper dx = L n/(n+1)`` (substitute
    ``u = 2x/L - 1``, ``int_{-1}^{1} (1 - |u|**n) du = 2 - 2/(n+1)``) and
    ``int_0^T section dw = T m/(m+1)``. At ``n = m = 2`` this is ``4LBT/9``.
    """
    n, m = longitudinal_exponent, vertical_exponent
    return length * beam * draft * n * m / ((n + 1.0) * (m + 1.0))


def analytic_waterplane_area(
    length: float, beam: float, longitudinal_exponent: float = 2.0
) -> float:
    """Exact waterplane area ``2 int_0^L (B/2) taper dx = L B n/(n+1)``."""
    n = longitudinal_exponent
    return length * beam * n / (n + 1.0)


def analytic_waterplane_inertia(
    length: float, beam: float, longitudinal_exponent: float = 2.0
) -> float:
    """Exact transverse waterplane inertia ``I_T = int y^2 dA`` about the centreline.

    The waterplane of the family is ``|y| <= b(x)`` with ``b = (B/2) taper(x)``, so
    ``I_T = int_0^L (2/3) b^3 dx``. Substituting ``u = 2x/L - 1`` (``dx = L du / 2``)
    and expanding the cube,

        I_T = (B^3 L / 12) * K(n),
        K(n) = 1 - 3/(n+1) + 3/(2n+1) - 1/(3n+1),

    because ``int_0^1 (1 - u^n)^3 du = 1 - 3/(n+1) + 3/(2n+1) - 1/(3n+1)``. At
    ``n = 2``: ``K = 16/35`` and ``I_T = 4 L B^3 / 105``.

    Note which exponent enters: the *half breadth* tapers as ``(1-|u|^n)``, the
    moment integrates its **cube**, so the Beta integral carries ``3n`` and not
    ``n``. The mesh never sees that distinction -- it integrates whatever polygon it
    is given -- which is why this reference is derived here and not read back from
    the integrator.
    """
    n = float(longitudinal_exponent)
    if not (math.isfinite(n) and n > 0.0):
        raise ValueError(f"longitudinal_exponent must be finite and > 0, got {n!r}")
    k = 1.0 - 3.0 / (n + 1.0) + 3.0 / (2.0 * n + 1.0) - 1.0 / (3.0 * n + 1.0)
    return (beam**3 * length / 12.0) * k


def analytic_waterplane_inertia_longitudinal(
    length: float, beam: float, longitudinal_exponent: float = 2.0
) -> float:
    """Exact longitudinal waterplane inertia about the centre of flotation.

    The waterplane is symmetric fore-aft, so the centre of flotation is at ``L/2``
    and ``I_L = int (x - L/2)^2 2 b(x) dx``. With ``u = 2x/L - 1`` and
    ``s = x - L/2 = (L/2) u``,

        I_L = (B L^3 / 4) * int_0^1 u^2 (1 - u^n) du
            = (B L^3 / 4) * (1/3 - 1/(n+3)),

    because ``int_0^1 u^2 (1-u^n) du = 1/3 - 1/(n+3)``. At ``n = 2`` the bracket is
    ``2/15`` and ``I_L = B L^3 / 30`` -- a third of the box value ``B L^3/12``,
    which is what a hull that tapers to a point at both ends should give.

    The parallel-axis term is **already removed**: the mesh routes report ``I_L``
    about the centre of flotation (see :func:`tensorcad.stability.waterplane_moments`),
    and comparing a centroidal reference with an about-origin one is the classic
    way to manufacture a disagreement that is not a defect.
    """
    n = float(longitudinal_exponent)
    if not (math.isfinite(n) and n > 0.0):
        raise ValueError(f"longitudinal_exponent must be finite and > 0, got {n!r}")
    k = 1.0 / 3.0 - 1.0 / (n + 3.0)
    return (beam * length**3 / 4.0) * k


def analytic_metacentric_radius_transverse(
    length: float,
    beam: float,
    draft: float,
    longitudinal_exponent: float = 2.0,
    vertical_exponent: float = 2.0,
) -> float:
    """Exact ``BM_T = I_T / nabla``: ``B^2 (n+1)(m+1) K(n) / (12 T n m)``.

    ``nabla = L B T n m / ((n+1)(m+1))`` is :func:`analytic_volume`, so the length
    cancels and only the breadth-to-draft ratio survives. At ``n = m = 2`` this is
    ``3 B^2 / (35 T)`` -- for the ITTC reference proportions ``B = 1 m``,
    ``T = 0.625 m`` that is ``0.137142857 m``.
    """
    n, m = float(longitudinal_exponent), float(vertical_exponent)
    volume = analytic_volume(length, beam, draft, n, m)
    return analytic_waterplane_inertia(length, beam, n) / volume


def analytic_metacentric_radius_longitudinal(
    length: float,
    beam: float,
    draft: float,
    longitudinal_exponent: float = 2.0,
    vertical_exponent: float = 2.0,
) -> float:
    """Exact ``BM_L = I_L / nabla``: ``L^2 (n+1)(m+1) (1/3 - 1/(n+3)) / (4 T n m)``.

    At ``n = m = 2`` this is ``3 L^2 / (40 T)``.
    """
    n, m = float(longitudinal_exponent), float(vertical_exponent)
    volume = analytic_volume(length, beam, draft, n, m)
    return analytic_waterplane_inertia_longitudinal(length, beam, n) / volume


def analytic_maximum_section_area(
    beam: float, draft: float, vertical_exponent: float = 2.0
) -> float:
    """Exact midship section area ``2 int_0^T (B/2) section dw = B T m/(m+1)``."""
    m = vertical_exponent
    return beam * draft * m / (m + 1.0)


def analytic_vertical_centre_buoyancy(
    draft: float, vertical_exponent: float = 2.0
) -> float:
    """Exact VCB **above the keel**: ``T - T(m+1)/(2(m+2))``.

    ``int_0^T w section(w) dw / int_0^T section(w) dw = T (m+1)/(2(m+2))`` is the
    depth of the centroid below the waterline; subtracting it from the draft
    converts to the keel-referenced height, which is what a stability booklet
    quotes. At ``m = 2`` that is ``5T/8``.
    """
    m = vertical_exponent
    return draft - draft * (m + 1.0) / (2.0 * (m + 2.0))


def analytic_vertical_centre_buoyancy_z(
    draft: float, vertical_exponent: float = 2.0
) -> float:
    """The same VCB in the package frame (``z = 0`` at the waterline), i.e. negative."""
    return analytic_vertical_centre_buoyancy(draft, vertical_exponent) - draft


def analytic_block_coefficient(
    longitudinal_exponent: float = 2.0, vertical_exponent: float = 2.0
) -> float:
    """``Cb = n m / ((n+1)(m+1))`` (``4/9`` for the classic Wigley)."""
    n, m = longitudinal_exponent, vertical_exponent
    return n * m / ((n + 1.0) * (m + 1.0))


def analytic_midship_coefficient(vertical_exponent: float = 2.0) -> float:
    """``Cm = m/(m+1)`` (``2/3`` for the classic Wigley)."""
    m = vertical_exponent
    return m / (m + 1.0)


def analytic_waterplane_coefficient(longitudinal_exponent: float = 2.0) -> float:
    """``Cwp = n/(n+1)`` (``2/3`` for the classic Wigley)."""
    n = longitudinal_exponent
    return n / (n + 1.0)


def analytic_prismatic_coefficient(
    longitudinal_exponent: float = 2.0, vertical_exponent: float = 2.0
) -> float:
    """``Cp = Cb/Cm = n/(n+1)`` (``2/3`` for the classic Wigley)."""
    n, m = longitudinal_exponent, vertical_exponent
    return analytic_block_coefficient(n, m) / analytic_midship_coefficient(m)


def _legendre_nodes_weights(order: int) -> tuple[list[float], list[float]]:
    """Gauss-Legendre nodes/weights on ``[-1, 1]`` by Newton iteration.

    Written out rather than imported so the quadrature path also keeps the
    package's standard-library-only promise.
    """
    nodes: list[float] = []
    weights: list[float] = []
    for i in range(1, order + 1):
        x = math.cos(math.pi * (i - 0.25) / (order + 0.5))
        for _ in range(100):
            p_prev, p = 1.0, x
            for k in range(2, order + 1):
                p_prev, p = p, ((2 * k - 1) * x * p - (k - 1) * p_prev) / k
            derivative = order * (x * p - p_prev) / (x * x - 1.0)
            step = p / derivative
            x -= step
            if abs(step) <= 1e-16:
                break
        p_prev, p = 1.0, x
        for k in range(2, order + 1):
            p_prev, p = p, ((2 * k - 1) * x * p - (k - 1) * p_prev) / k
        derivative = order * (x * p - p_prev) / (x * x - 1.0)
        nodes.append(x)
        weights.append(2.0 / ((1.0 - x * x) * derivative * derivative))
    return nodes, weights


def _composite_gauss(
    function: Callable[[float], float], low: float, high: float, *, panels: int, order: int
) -> float:
    """Composite Gauss-Legendre integral of ``function`` over ``[low, high]``."""
    nodes, weights = _legendre_nodes_weights(order)
    total = 0.0
    for p in range(panels):
        a = low + (high - low) * p / panels
        b = low + (high - low) * (p + 1) / panels
        half = 0.5 * (b - a)
        mid = 0.5 * (a + b)
        total += half * math.fsum(w * function(mid + half * x) for x, w in zip(nodes, weights, strict=False))
    return total


def _wetted_integrand(parameters: WigleyParameters, x: float, depth: float) -> float:
    """``sqrt(1 + (dy/dx)^2 + (dy/dw)^2)`` at ``(x, depth)`` on the starboard flank."""
    p = parameters
    n, m = p.longitudinal_exponent, p.vertical_exponent
    u = 2.0 * x / p.length - 1.0
    taper = 1.0 - math.pow(min(abs(u), 1.0), n)
    # d(taper)/dx = -n |u|^(n-1) sign(u) * 2/L
    dtaper = -n * math.pow(min(abs(u), 1.0), n - 1.0) * (1.0 if u >= 0.0 else -1.0) * (2.0 / p.length)
    if depth <= 0.0:
        section = 1.0
        dsection = -1.0 / p.draft if m == 1.0 else 0.0
    elif depth >= p.draft:
        section = 0.0
        dsection = -m / p.draft
    else:
        section = 1.0 - math.pow(depth / p.draft, m)
        dsection = -m * math.pow(depth / p.draft, m - 1.0) / p.draft
    half = 0.5 * p.beam
    dy_dx = half * dtaper * section
    dy_dw = half * taper * dsection
    return math.sqrt(1.0 + dy_dx * dy_dx + dy_dw * dy_dw)


def analytic_wetted_surface_area(
    length: float,
    beam: float,
    draft: float,
    longitudinal_exponent: float = 2.0,
    vertical_exponent: float = 2.0,
    *,
    panels: int = 4,
    order: int = 24,
) -> float:
    """Wetted surface by composite Gauss-Legendre on the analytic flank.

    ``S = 2 int_0^L int_0^T sqrt(1 + yx^2 + yw^2) dw dx`` has no elementary
    primitive, so unlike the other quantities this is a *quadrature* reference
    rather than a closed form -- but an independent one: it evaluates the analytic
    surface, not the mesh. The integrand is smooth for ``n, m >= 1`` (the Wigley
    has a finite keel angle and a finite waterline rake), so Gauss-Legendre
    converges geometrically; :func:`wetted_area_evidence` reports the
    order-to-order difference so the claim stays checkable rather than asserted.

    Use an **even** panel count: the second derivative of ``|2x/L - 1|**n`` jumps at
    midship for non-integer ``n``, so an even count puts a panel boundary exactly on
    the kink. Measured on the reference case the default (4 panels x 24 nodes) is
    exact to 4e-16 at ``n = m = 2`` and to 8e-10 at ``n = 1.5, m = 3``, in 78 ms.
    """
    parameters = WigleyParameters(
        length=length,
        beam=beam,
        draft=draft,
        longitudinal_exponent=longitudinal_exponent,
        vertical_exponent=vertical_exponent,
    )

    def inner(x: float) -> float:
        return _composite_gauss(
            lambda w: _wetted_integrand(parameters, x, w), 0.0, draft, panels=panels, order=order
        )

    return 2.0 * _composite_gauss(inner, 0.0, length, panels=panels, order=order)


def _resolve(parameters: Any, overrides: dict[str, Any]) -> WigleyParameters:
    """Accept a parameter object, an ``(L, B, T)`` triple, or ``None``."""
    if parameters is None:
        base = reference_parameters()
    elif isinstance(parameters, WigleyParameters):
        base = parameters
    elif isinstance(parameters, (tuple, list)) and len(parameters) == 3:
        base = reference_parameters(
            length=float(parameters[0]), beam=float(parameters[1]), draft=float(parameters[2])
        )
    else:
        base = WigleyParameters(
            length=float(parameters.length),
            beam=float(parameters.beam),
            draft=float(parameters.draft),
            water_density=float(getattr(parameters, "water_density", 1025.0)),
            longitudinal_exponent=float(getattr(parameters, "longitudinal_exponent", 2.0)),
            vertical_exponent=float(getattr(parameters, "vertical_exponent", 2.0)),
        )
    if overrides:
        values = base.as_dict()
        values.update(overrides)
        base = WigleyParameters(**values)
    return base


@dataclass(frozen=True)
class WigleyAnalytic:
    """Exact hydrostatics of the Wigley surface.

    Field names match :class:`tensorcad.hydrostatics.Hydrostatics` where the two
    report the same quantity, so :meth:`compare` is a one-to-one subtraction.
    """

    length: float
    beam: float
    draft: float
    water_density: float
    longitudinal_exponent: float
    vertical_exponent: float
    displacement_volume: float
    displacement_mass: float
    waterplane_area: float
    wetted_surface_area: float
    maximum_section_area: float
    longitudinal_centre_buoyancy: float
    longitudinal_centre_flotation: float
    vertical_centre_buoyancy: float
    block_coefficient: float
    prismatic_coefficient: float
    midship_coefficient: float
    waterplane_coefficient: float
    basis: str

    def as_dict(self) -> dict[str, Any]:
        """JSON-ready payload with the same key names as the mesh hydrostatics."""
        return {
            "schema": "tensorcad.wigley.analytic/1",
            "units": "m",
            "form": "y = (B/2) (1 - |2x/L - 1|^n) (1 - (w/T)^m)",
            "inputs": {
                "length_m": self.length,
                "beam_m": self.beam,
                "draft_m": self.draft,
                "water_density_kg_per_m3": self.water_density,
                "longitudinal_exponent": self.longitudinal_exponent,
                "vertical_exponent": self.vertical_exponent,
            },
            "displacement_volume_m3": self.displacement_volume,
            "displacement_mass_kg": self.displacement_mass,
            "waterplane_area_m2": self.waterplane_area,
            "wetted_surface_area_m2": self.wetted_surface_area,
            "maximum_section_area_m2": self.maximum_section_area,
            "longitudinal_centre_buoyancy_m": self.longitudinal_centre_buoyancy,
            "longitudinal_centre_flotation_m": self.longitudinal_centre_flotation,
            "vertical_centre_buoyancy_m": self.vertical_centre_buoyancy,
            "block_coefficient": self.block_coefficient,
            "prismatic_coefficient": self.prismatic_coefficient,
            "midship_coefficient": self.midship_coefficient,
            "waterplane_coefficient": self.waterplane_coefficient,
            "basis": self.basis,
        }

    def compare(
        self, measured: Any, *, quantities: tuple[str, ...] | None = None
    ) -> dict[str, Any]:
        """Signed relative error of a mesh :class:`Hydrostatics` against this reference.

        Signed, not absolute: a sampled loft sits *inside* the analytic surface, so
        volume, areas and section area are expected to come out low at every
        resolution. A sign flip would be a defect, not a refinement level.
        """
        names = quantities if quantities is not None else ANALYTIC_COMPARABLE_QUANTITIES
        errors: dict[str, float] = {}
        values: dict[str, dict[str, float]] = {}
        for name in names:
            exact = float(getattr(self, name))
            value = float(getattr(measured, name))
            errors[name] = (value - exact) / abs(exact) if exact else float("nan")
            values[name] = {"measured": value, "reference": exact}
        finite = [abs(v) for v in errors.values() if math.isfinite(v)]
        return {
            "quantities": values,
            "relative_error": errors,
            "max_abs_relative_error": max(finite) if finite else float("nan"),
            "reference_basis": self.basis,
        }


def wigley_analytic_hydrostatics(
    parameters: Any = None, **overrides: Any
) -> WigleyAnalytic:
    """Closed-form hydrostatics of the analytic surface the mesh samples.

    ``parameters`` is a :class:`WigleyParameters` (or anything carrying
    ``length``/``beam``/``draft``/exponents), ``None`` for the ITTC reference
    case. Sampling counts are ignored on purpose: the reference describes the
    surface, not the sampling of it.
    """
    p = _resolve(parameters, overrides)
    n, m = p.longitudinal_exponent, p.vertical_exponent
    volume = analytic_volume(p.length, p.beam, p.draft, n, m)
    return WigleyAnalytic(
        length=p.length,
        beam=p.beam,
        draft=p.draft,
        water_density=p.water_density,
        longitudinal_exponent=n,
        vertical_exponent=m,
        displacement_volume=volume,
        displacement_mass=volume * p.water_density,
        waterplane_area=analytic_waterplane_area(p.length, p.beam, n),
        wetted_surface_area=analytic_wetted_surface_area(p.length, p.beam, p.draft, n, m),
        maximum_section_area=analytic_maximum_section_area(p.beam, p.draft, m),
        longitudinal_centre_buoyancy=0.5 * p.length,
        longitudinal_centre_flotation=0.5 * p.length,
        vertical_centre_buoyancy=analytic_vertical_centre_buoyancy_z(p.draft, m),
        block_coefficient=analytic_block_coefficient(n, m),
        prismatic_coefficient=analytic_prismatic_coefficient(n, m),
        midship_coefficient=analytic_midship_coefficient(m),
        waterplane_coefficient=analytic_waterplane_coefficient(n),
        basis="closed forms of the Wigley algebraic surface; wetted area by composite Gauss-Legendre",
    )


#: Quantities :class:`WigleyStability` shares with
#: :class:`tensorcad.stability.InitialStability`, so :meth:`WigleyStability.compare` is
#: a one-to-one subtraction and not a name-matching exercise.
STABILITY_COMPARABLE_QUANTITIES: tuple[str, ...] = (
    "transverse_waterplane_inertia",
    "longitudinal_waterplane_inertia",
    "kb",
    "bm_transverse",
    "bm_longitudinal",
    "km_transverse",
    "km_longitudinal",
    "gm_transverse",
    "gm_longitudinal",
)

#: A measured stability result reports its waterplane moments one level down, which is
#: why the two inertia names need a route to the nested record.
_WATERPLANE_ROUTE: dict[str, str] = {
    "transverse_waterplane_inertia": "transverse_second_moment",
    "longitudinal_waterplane_inertia": "longitudinal_second_moment",
}


@dataclass(frozen=True)
class WigleyStability:
    """Exact initial stability of the Wigley surface: ``I_T``, ``I_L``, ``BM``, ``KM``, ``GM``.

    Field names match :class:`tensorcad.stability.InitialStability` (and
    :class:`tensorcad.hydrostatics.Hydrostatics` for the two inertias) so
    :meth:`compare` is a subtraction, not a translation. ``KG`` is an input for the same
    reason it is an input there: it is a loading condition, not a property of the shape.

    Heights are above the keel, inertias about the centreline / centre of flotation, and
    ``GM = KB + BM - KG`` -- no free-surface correction, no trim, no heel. This is the
    intact, upright reference, which is exactly what a mesh route should be measured
    against before it is trusted on anything else.
    """

    length: float
    beam: float
    draft: float
    water_density: float
    longitudinal_exponent: float
    vertical_exponent: float
    kg: float
    displacement_volume: float
    transverse_waterplane_inertia: float
    longitudinal_waterplane_inertia: float
    kb: float
    bm_transverse: float
    bm_longitudinal: float
    km_transverse: float
    km_longitudinal: float
    gm_transverse: float
    gm_longitudinal: float
    basis: str

    def as_dict(self) -> dict[str, Any]:
        """JSON-ready payload; keys mirror the mesh stability record."""
        return {
            "schema": "tensorcad.wigley.stability/1",
            "units": "m",
            "form": "y = (B/2) (1 - |2x/L - 1|^n) (1 - (w/T)^m)",
            "inputs": {
                "length_m": self.length,
                "beam_m": self.beam,
                "draft_m": self.draft,
                "water_density_kg_per_m3": self.water_density,
                "longitudinal_exponent": self.longitudinal_exponent,
                "vertical_exponent": self.vertical_exponent,
                "kg_above_keel_m": self.kg,
            },
            "displacement_volume_m3": self.displacement_volume,
            "transverse_waterplane_inertia_m4": self.transverse_waterplane_inertia,
            "longitudinal_waterplane_inertia_m4": self.longitudinal_waterplane_inertia,
            "kb_m": self.kb,
            "bm_transverse_m": self.bm_transverse,
            "bm_longitudinal_m": self.bm_longitudinal,
            "km_transverse_m": self.km_transverse,
            "km_longitudinal_m": self.km_longitudinal,
            "gm_transverse_m": self.gm_transverse,
            "gm_longitudinal_m": self.gm_longitudinal,
            "basis": self.basis,
        }

    def compare(
        self, measured: Any, *, quantities: tuple[str, ...] | None = None
    ) -> dict[str, Any]:
        """Signed relative error of a measured stability result against this reference.

        ``measured`` is a :class:`tensorcad.stability.InitialStability` (or anything with
        the same attribute names). Signed, because a sampled loft sits inside the analytic
        surface: the waterplane *chords* cut the corner off a convex taper, so ``I_T`` is
        expected to come out low at every resolution. A sign flip is a defect, not a
        refinement level.
        """
        names = quantities if quantities is not None else STABILITY_COMPARABLE_QUANTITIES
        errors: dict[str, float] = {}
        values: dict[str, dict[str, float]] = {}
        for name in names:
            exact = float(getattr(self, name))
            value = _measured_value(measured, name)
            errors[name] = (value - exact) / abs(exact) if exact else float("nan")
            values[name] = {"measured": value, "reference": exact}
        finite = [abs(v) for v in errors.values() if math.isfinite(v)]
        return {
            "quantities": values,
            "relative_error": errors,
            "max_abs_relative_error": max(finite) if finite else float("nan"),
            "reference_basis": self.basis,
        }


def _measured_value(measured: Any, name: str) -> float:
    """Read ``name`` off a stability result, following one level into its waterplane."""
    if hasattr(measured, name):
        return float(getattr(measured, name))
    nested = getattr(measured, "waterplane", None)
    leaf = _WATERPLANE_ROUTE.get(name)
    if nested is not None and leaf is not None and hasattr(nested, leaf):
        return float(getattr(nested, leaf))
    raise AttributeError(f"measured result has no {name!r} (nor a waterplane.{leaf})")


def analytic_initial_stability(parameters: Any = None, **overrides: Any) -> WigleyStability:
    """Closed-form ``KM``/``GM`` of the Wigley surface for one loading condition.

    ``parameters`` is a :class:`WigleyParameters` (or anything carrying
    ``length``/``beam``/``draft``/exponents) and ``kg`` is passed in ``overrides`` or
    read off the parameters object; it is **required** for the same reason
    :func:`tensorcad.stability.initial_stability` requires it -- a default of zero would
    silently report ``GM = KM``.

    Sampling counts are ignored: the reference describes the surface, not the sampling
    of it. Pass ``kg`` above the keel.
    """
    kg = overrides.get("kg", getattr(parameters, "kg", None))
    geometry_overrides = {k: v for k, v in overrides.items() if k != "kg"}
    p = _resolve(parameters, geometry_overrides)
    if kg is None:
        raise ValueError("kg (height of the centre of gravity above the keel) is required")
    kg = float(kg)
    if not math.isfinite(kg):
        raise ValueError(f"kg must be finite, got {kg!r}")
    n, m = p.longitudinal_exponent, p.vertical_exponent
    volume = analytic_volume(p.length, p.beam, p.draft, n, m)
    inertia_t = analytic_waterplane_inertia(p.length, p.beam, n)
    inertia_l = analytic_waterplane_inertia_longitudinal(p.length, p.beam, n)
    kb = analytic_vertical_centre_buoyancy(p.draft, m)
    bm_t = inertia_t / volume
    bm_l = inertia_l / volume
    km_t = kb + bm_t
    km_l = kb + bm_l
    return WigleyStability(
        length=p.length,
        beam=p.beam,
        draft=p.draft,
        water_density=p.water_density,
        longitudinal_exponent=n,
        vertical_exponent=m,
        kg=kg,
        displacement_volume=volume,
        transverse_waterplane_inertia=inertia_t,
        longitudinal_waterplane_inertia=inertia_l,
        kb=kb,
        bm_transverse=bm_t,
        bm_longitudinal=bm_l,
        km_transverse=km_t,
        km_longitudinal=km_l,
        gm_transverse=km_t - kg,
        gm_longitudinal=km_l - kg,
        basis="closed forms: I_T=(B^3 L/12)K(n), I_L=(B L^3/4)(1/3-1/(n+3)), KB=VCB, BM=I/nabla",
    )


def wetted_area_evidence(parameters: Any = None) -> dict[str, Any]:
    """Quadrature convergence of the only non-closed-form reference quantity."""
    p = _resolve(parameters, {})
    # Even panel counts throughout: a boundary then lands on midship, which is
    # where the integrand kinks for non-integer exponents.
    settings = ((2, 12), (4, 24), (8, 48))
    values = {
        f"{panels}x{order}": analytic_wetted_surface_area(
            p.length,
            p.beam,
            p.draft,
            p.longitudinal_exponent,
            p.vertical_exponent,
            panels=panels,
            order=order,
        )
        for panels, order in settings
    }
    ordered = [values[f"{panels}x{order}"] for panels, order in settings]
    return {
        "method": "composite Gauss-Legendre on the analytic flank (not on the mesh)",
        "levels": {
            f"{panels}x{order}": value
            for (panels, order), value in zip(settings, ordered, strict=False)
        },
        "relative_change_coarse_to_fine": abs(ordered[1] - ordered[0]) / abs(ordered[1]),
        "relative_change_fine_to_finer": abs(ordered[2] - ordered[1]) / abs(ordered[2]),
    }


# -------------------------------------------------------------------- geometry


@dataclass(frozen=True)
class WigleyModel:
    """A generated Wigley hull surface plus everything needed to reproduce it."""

    parameters: WigleyParameters
    mesh: Mesh
    offsets: tuple[tuple[float, float], ...]
    validation: ValidationReport = field(default_factory=ValidationReport)
    hydrostatics: Hydrostatics | None = None

    @property
    def vertices(self) -> tuple[tuple[float, float, float], ...]:
        return self.mesh.vertices

    @property
    def faces(self) -> tuple[tuple[int, int, int], ...]:
        return self.mesh.faces

    @property
    def units(self) -> str:
        return self.mesh.units

    @property
    def coordinate_frame(self) -> str:
        return self.mesh.coordinate_frame

    def metadata(self) -> dict[str, Any]:
        data: dict[str, Any] = {
            "schema_version": "1.0",
            "model_type": "closed_sampled_wigley_displacement_hull",
            "hull_form": "wigley",
            "parameters": self.parameters.as_dict(),
            "units": self.units,
            "coordinate_frame": self.coordinate_frame,
            "generator": {
                "kind": "analytic_lofted_surface",
                "surface": "y = (B/2) (1 - |2x/L - 1|^n) (1 - (w/T)^m)",
                "longitudinal_taper": f"1 - |2x/L - 1|^{self.parameters.longitudinal_exponent:g}",
                "vertical_taper": f"1 - (w/T)^{self.parameters.vertical_exponent:g}",
                "topology": (
                    "one closed section ring per station, lofted, capped at the design waterline"
                ),
                "backend": "tensorcad.geometry (standard library)",
            },
            "mesh": self.mesh.metadata(),
            "offset_count": len(self.offsets),
            "validation": self.validation.to_dict(),
            "engineering_status": (
                "geometry and hydrostatic integration of a standard benchmark hull form; "
                "appendages, transom, propeller, stability, resistance and structural "
                "checks are not implemented"
            ),
        }
        if self.hydrostatics is not None:
            data["hydrostatics"] = self.hydrostatics.to_dict()
        data["analytic"] = wigley_analytic_hydrostatics(self.parameters).as_dict()
        return data


def build_mesh(parameters: WigleyParameters) -> tuple[Mesh, tuple[tuple[float, float], ...]]:
    """Loft and cap the Wigley surface; returns the mesh and the waterline offsets.

    Deterministic: identical parameters produce identical vertex and face tuples,
    so a solver run can be replayed from the parameters alone. The deck cap at the
    design waterline closes the solid, which is what lets the volume be integrated
    from the triangles at all.
    """
    parameters.validate()
    p = parameters
    ring_size = 2 * p.nz - 1
    builder = MeshBuilder(units=UNITS, coordinate_frame=COORDINATE_FRAME)

    rings: list[list[int]] = []
    for i in range(p.nx):
        rings.append([builder.add_vertex(_section_point(p, i, j)) for j in range(ring_size)])

    for i in range(p.nx - 1):
        for j in range(ring_size - 1):
            builder.add_quad(rings[i][j], rings[i + 1][j], rings[i + 1][j + 1], rings[i][j + 1])
        # Deck cap: the strip between the two waterline ends of the rings.
        builder.add_quad(
            rings[i][ring_size - 1],
            rings[i + 1][ring_size - 1],
            rings[i + 1][0],
            rings[i][0],
        )

    mesh = meshlib.orient_outward(builder.build()).mesh
    offsets = tuple(
        (
            p.length * i / (p.nx - 1),
            mesh.vertices[rings[i][0]][1] - mesh.vertices[rings[i][ring_size - 1]][1],
        )
        for i in range(p.nx)
    )
    return mesh, offsets


def check_wigley_consistency(
    parameters: WigleyParameters, mesh: Mesh, offsets: tuple[tuple[float, float], ...]
) -> list[Check]:
    """Checks that tie the generated surface back to the requested parameters."""
    p = parameters
    checks: list[Check] = [
        Check(
            "offset_per_station",
            len(offsets) == p.nx,
            f"{len(offsets)} offset rows for {p.nx} stations",
        )
    ]
    ends_ok = bool(offsets) and offsets[0][1] == 0.0 and offsets[-1][1] == 0.0
    checks.append(
        Check(
            "offsets_close_at_both_ends",
            ends_ok,
            "" if ends_ok else "breadth must be exactly zero at the stem and the post",
        )
    )
    realised = max((width for _, width in offsets), default=0.0)
    achievable = max(longitudinal_taper(p, i) for i in range(p.nx))
    expected = p.beam * achievable
    close = math.isclose(realised, expected, rel_tol=1e-12, abs_tol=1e-12)
    checks.append(
        Check(
            "beam_matches_grid_achievable_maximum",
            close,
            f"realised maximum breadth {realised:.9g} vs grid-achievable {expected:.9g}",
        )
    )
    if achievable < 1.0:
        checks.append(
            Check(
                "requested_beam_not_achievable_on_this_grid",
                False,
                f"no station lands on midship with nx={p.nx}, so the realised beam "
                f"{realised:.9g} is {100 * (1 - realised / p.beam):.3f}% under the requested "
                f"{p.beam:.9g}; use an odd nx",
                severity=WARNING,
            )
        )
    # The loft only interpolates along x, so at the sampled stations the offsets
    # must sit on the analytic surface exactly, not merely close to it.
    worst = max(
        (abs(width - 2.0 * half_breadth(p, x, 0.0)) for x, width in offsets), default=0.0
    )
    checks.append(
        Check(
            "waterline_offsets_on_analytic_surface",
            worst <= 1e-12 * max(p.beam, 1.0),
            f"largest deviation of a waterline offset from the analytic breadth is {worst:.3e} m",
        )
    )
    for name, value, low, high in (
        ("length_to_beam", p.length / p.beam, 1.5, 15.0),
        ("beam_to_draft", p.beam / p.draft, 1.0, 12.0),
    ):
        inside = low <= value <= high
        checks.append(
            Check(
                f"ratio_{name}",
                inside,
                f"{name} = {value:.4g}" if inside else f"{name} = {value:.4g} outside [{low:g}, {high:g}]",
                severity=WARNING,
            )
        )
    return checks


def _dimension_target(parameters: WigleyParameters) -> DimensionTarget:
    return DimensionTarget(
        length=parameters.length,
        beam=parameters.beam * max(longitudinal_taper(parameters, i) for i in range(parameters.nx)),
        draft=parameters.draft,
        waterline_z=0.0,
    )


def _sampling_check(parameters: WigleyParameters, hydrostatics: Hydrostatics) -> Check:
    """Warn when the grid, rather than the form, is what makes the numbers wrong."""
    exact = analytic_volume(
        parameters.length,
        parameters.beam,
        parameters.draft,
        parameters.longitudinal_exponent,
        parameters.vertical_exponent,
    )
    measured = hydrostatics.displacement_volume
    error = (measured / exact - 1.0) if exact else float("nan")
    ok = math.isfinite(error) and abs(error) <= SAMPLING_ERROR_WARNING
    detail = (
        f"grid {parameters.nx}x{parameters.nz} integrates {measured:.6g} m^3 against a closed "
        f"form {exact:.6g} m^3 ({100.0 * error:+.3f} %)"
    )
    if not ok:
        detail += (
            f"; the discretisation error exceeds {100.0 * SAMPLING_ERROR_WARNING:.0f} % -- raise "
            "nx/nz, or read docs/wigley.md for the measured convergence"
        )
    return Check("sampling_error_within_tolerance", ok, detail, severity=WARNING)


def generate_wigley(parameters: Any = None, **overrides: Any) -> WigleyModel:
    """Generate a closed, validated Wigley hull from principal dimensions.

    Raises :class:`ValueError` when the parameters are invalid or when the
    generated surface fails any geometric validity check -- an unusable mesh is
    never handed back silently.
    """
    p = _resolve(parameters, overrides)
    mesh, offsets = build_mesh(p)
    report = merge(
        inspect_geometry(mesh),
        check_dimensions(mesh, _dimension_target(p)),
        check_wigley_consistency(p, mesh, offsets),
    )
    errors = report.errors
    if errors:
        raise ValueError("generated Wigley hull failed validation: " + "; ".join(errors))
    hydrostatics = compute_hydrostatics(mesh, 0.0, water_density=p.water_density, closed=True)
    report = merge(list(report.checks), [_sampling_check(p, hydrostatics)])
    return WigleyModel(
        parameters=p, mesh=mesh, offsets=offsets, validation=report, hydrostatics=hydrostatics
    )


def analytic_relative_errors(parameters: Any = None, **overrides: Any) -> dict[str, Any]:
    """Mesh hydrostatics against the closed forms, per quantity, signed.

    Also reports the grid the numbers came from and the analytic mass, so a stored
    JSON is self-describing. The sign convention is documented on
    :meth:`WigleyAnalytic.compare`.
    """
    p = _resolve(parameters, overrides)
    mesh, _ = build_mesh(p)
    measured = compute_hydrostatics(mesh, 0.0, water_density=p.water_density, closed=True)
    reference = wigley_analytic_hydrostatics(p)
    comparison = reference.compare(measured)
    return {
        "status": "measured",
        "method": "mesh integrals against closed forms of the analytic Wigley surface",
        "grid": {"nx": p.nx, "nz": p.nz},
        "parameters": p.as_dict(),
        "measured": {
            name: comparison["quantities"][name]["measured"]
            for name in ANALYTIC_COMPARABLE_QUANTITIES
        },
        "reference": {
            name: comparison["quantities"][name]["reference"]
            for name in ANALYTIC_COMPARABLE_QUANTITIES
        },
        "relative_error": comparison["relative_error"],
        "max_abs_relative_error": comparison["max_abs_relative_error"],
        "reference_basis": comparison["reference_basis"],
    }


def convergence_study(
    resolutions: tuple[tuple[int, int], ...] = ((21, 9), (41, 13), (81, 21), (161, 41)),
    parameters: Any = None,
    **overrides: Any,
) -> dict[str, Any]:
    """Refine the grid and report how each integrated quantity approaches the closed form.

    The spacing measure is ``h = 1/(nz-1)`` (vertical spacing as a fraction of the
    draft): the loft error is dominated by the chords cutting inside the section
    curve, so it is the vertical count that sets the rate. Reported so a caller can
    see the order instead of trusting the coarsest grid.
    """
    p = _resolve(parameters, overrides)
    reference = wigley_analytic_hydrostatics(p)
    rows: list[dict[str, Any]] = []
    for nx, nz in resolutions:
        case = _resolve(p, {"nx": int(nx), "nz": int(nz)})
        mesh, _ = build_mesh(case)
        hydrostatics = compute_hydrostatics(mesh, 0.0, water_density=case.water_density, closed=True)
        comparison = reference.compare(hydrostatics)
        measured = {
            "measured": {
                name: comparison["quantities"][name]["measured"]
                for name in ANALYTIC_COMPARABLE_QUANTITIES
            },
            "relative_error": comparison["relative_error"],
            "max_abs_relative_error": comparison["max_abs_relative_error"],
        }
        rows.append(
            {
                "nx": int(nx),
                "nz": int(nz),
                "h": 1.0 / (int(nz) - 1),
                "triangles": len(mesh.faces),
                **measured["measured"],
                "relative_error": measured["relative_error"],
                "max_abs_relative_error": measured["max_abs_relative_error"],
            }
        )
    orders: dict[str, float] = {}
    for name in ("displacement_volume", "waterplane_area", "wetted_surface_area"):
        first, last = rows[0]["relative_error"][name], rows[-1]["relative_error"][name]
        if first and last and first * last > 0.0:
            span = rows[0]["h"] / rows[-1]["h"]
            orders[name] = math.log(abs(first) / abs(last)) / math.log(span)
    return {
        "case": p.as_dict(),
        "reference": reference.as_dict(),
        "rows": rows,
        "observed_order": orders,
        "note": (
            "order = log(|e_coarse|/|e_fine|)/log(h_coarse/h_fine) fitted across the whole "
            "sweep; the loft error is O(h^2) once h is small"
        ),
    }


def export_wigley(
    model: WigleyModel,
    output_dir: str | Any,
    *,
    stl: bool = False,
) -> dict[str, str]:
    """Write ``wigley.obj``, ``wigley_offsets.csv`` and ``wigley.json``.

    Reuses the package writers so the artifacts are byte-compatible with the
    ``hull`` family's exports; the metadata document additionally carries the
    closed-form reference under ``analytic``.
    """
    from pathlib import Path

    from .export import write_metadata, write_obj, write_offsets_csv, write_stl

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    paths = {
        "mesh": write_obj(model.mesh, out / "wigley.obj"),
        "offsets": write_offsets_csv(model.offsets, out / "wigley_offsets.csv"),
        "metadata": write_metadata(model.metadata(), out / "wigley.json"),
    }
    if stl:
        paths["stl"] = write_stl(model.mesh, out / "wigley.stl")
    return {key: str(value) for key, value in paths.items()}


# --------------------------------------------------- curvature energy (fairness)
@dataclass(frozen=True)
class WigleyFlankEnergy:
    """Smoothing energy of the analytic Wigley flank, by quadrature on the surface.

    Every number here is an integral of an **analytic** differential-geometry
    integrand over the analytic surface ``y = h(x, w)`` -- not of a mesh estimate --
    so it is the reference the mesh path in :mod:`tensorcad.fairness_energy` is
    measured against. Both flanks are included: the starboard one is integrated and
    doubled, which is exact because ``h`` enters the metric squared.

    The integrands are ``k1**2 + k2**2``, ``H**2`` and ``K``, all dimensionless
    (``dA`` carries ``1/m**2``), so a geometrically similar hull must return the
    same three numbers at any model scale. That is the first claim the tests make.

    Why ``gaussian`` is *not* ``4 pi``. Gauss-Bonnet charges ``4 pi`` to a closed,
    boundary-free, genus-0 surface, and the lofted Wigley solid is closed -- but it
    is closed by a flat deck cap **and** by both flanks meeting along the keel and
    stem lines, which leaves four conical corner points (stem/stern x keel,
    stem/stern x waterline). On the default 41 x 13 mesh each of those corners
    carries a 180-degree angle defect; the smooth flank therefore carries
    ``4 pi`` *minus* what the corners hold, and here that residue is a small
    positive number (``~0.09``). The mesh path in
    :mod:`tensorcad.fairness_energy` reproduces the same split -- its defect sum is
    ``4 pi`` identically and the corner vertices hold almost all of it -- so the two
    paths are comparable on the smooth part only, which is exactly the comparison
    the energy study performs.
    """

    length: float
    beam: float
    draft: float
    longitudinal_exponent: float
    vertical_exponent: float
    wetted_area: float
    bending: float
    willmore: float
    gaussian: float
    corner_charge: float
    order: int
    panels: int
    quadrature_relative_change: float

    @property
    def bending_per_unit_area(self) -> float:
        return self.bending / self.wetted_area if self.wetted_area > 0.0 else float("nan")

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["schema"] = "tensorcad.wigley.flank_energy/1"
        payload["bending_per_unit_area"] = self.bending_per_unit_area
        payload["willmore_inequality_slack"] = self.willmore - 4.0 * math.pi
        payload["units"] = {
            "bending": "dimensionless",
            "willmore": "dimensionless",
            "gaussian": "dimensionless",
            "wetted_area": "m^2",
        }
        return payload


def _flank_derivatives(
    parameters: WigleyParameters, xi: float, depth: float
) -> tuple[float, float, float, float, float, float]:
    """``(h, hx, hw, hxx, hxw, hww)`` of the analytic flank at ``(xi, depth)``.

    ``xi = 2x/L - 1`` in ``[-1, 1]`` and ``depth`` measured down from the waterline,
    so the surface is the graph ``y = h(x, w)`` over the rectangle
    ``[0, L] x [0, T]``. The taper ``1 - |xi|**n`` and the section ``1 - (w/T)**m``
    are polynomials for the classic ``n = m = 2``; ``|xi|**n`` is ``C1`` across
    ``xi = 0`` for every admissible ``n > 1`` and the Gauss nodes avoid the kink when
    an even panel count puts a panel boundary there.
    """
    p = parameters
    n, m = p.longitudinal_exponent, p.vertical_exponent
    scale = 0.5 * p.beam
    magnitude = min(abs(xi), 1.0)
    sign = 1.0 if xi >= 0.0 else -1.0

    taper = 1.0 - math.pow(magnitude, n)
    dtaper = -n * math.pow(magnitude, n - 1.0) * sign
    d2taper = -n * (n - 1.0) * math.pow(magnitude, n - 2.0)

    ratio = min(max(depth / p.draft, 0.0), 1.0)
    section = 1.0 - math.pow(ratio, m)
    dsection = -m * math.pow(ratio, m - 1.0) / p.draft
    d2section = -m * (m - 1.0) * math.pow(ratio, m - 2.0) / (p.draft * p.draft)

    inverse_l = 2.0 / p.length
    return (
        scale * taper * section,
        scale * dtaper * section * inverse_l,
        scale * taper * dsection,
        scale * d2taper * section * inverse_l * inverse_l,
        scale * dtaper * dsection * inverse_l,
        scale * taper * d2section,
    )


def _flank_integrand(
    parameters: WigleyParameters, xi: float, depth: float
) -> tuple[float, float, float, float]:
    """``(bending, willmore, gaussian, jacobian)`` at one point of the flank.

    The mean and Gaussian curvatures of the graph ``y = h(x, w)`` are the standard
    first/second-fundamental-form expressions::

        H = ((1 + hw**2) hxx - 2 hx hw hxw + (1 + hx**2) hww) / (2 J**3)
        K = (hxx hww - hxw**2) / J**4,        J = sqrt(1 + hx**2 + hw**2)

    and the bending density follows from the identity
    ``k1**2 + k2**2 = 4 H**2 - 2 K``, so the three numbers cannot drift apart.
    """
    _, hx, hw, hxx, hxw, hww = _flank_derivatives(parameters, xi, depth)
    jacobian = math.sqrt(1.0 + hx * hx + hw * hw)
    mean = ((1.0 + hw * hw) * hxx - 2.0 * hx * hw * hxw + (1.0 + hx * hx) * hww) / (
        2.0 * jacobian**3
    )
    gaussian = (hxx * hww - hxw * hxw) / (jacobian**4)
    return 4.0 * mean * mean - 2.0 * gaussian, mean * mean, gaussian, jacobian


def analytic_flank_energy(
    parameters: Any = None,
    *,
    order: int = 64,
    panels: int = 4,
    **overrides: Any,
) -> WigleyFlankEnergy:
    """Smoothing energy of the analytic Wigley flank: ``int (k1**2+k2**2) dA`` etc.

    Tensor-product Gauss-Legendre over ``xi in [-1, 1]`` (``panels`` composite
    panels, so an even count puts a panel boundary on the midship kink of
    ``|xi|**n``) and ``w in [0, T]``. ``order`` is the node count per direction; the
    default ``64 x 4`` is converged to the round-off floor on the classic
    ``n = m = 2`` Wigley, and
    :attr:`WigleyFlankEnergy.quadrature_relative_change` reports the order-to-order
    difference so the claim stays checkable rather than asserted.

    The reason this exists at all: the mesh path's energy can only be judged
    against something that is not also a mesh estimate. Nothing in the literature
    quotes a Wigley bending energy, so it is derived here -- and the derivation is
    checked two ways, by quadrature doubling and against the mesh at three
    resolutions (``examples/fairness_energy``).
    """
    p = _resolve(parameters, dict(overrides))
    p.validate()
    if order < 4:
        raise ValueError(f"order must be at least 4, got {order!r}")
    if panels < 1 or panels % 2:
        raise ValueError(f"panels must be a positive even number, got {panels!r}")

    def integrate(node_count: int) -> tuple[float, float, float, float]:
        nodes, weights = _legendre_nodes_weights(node_count)
        depth_nodes = [0.5 * p.draft * (node + 1.0) for node in nodes]
        depth_weights = [0.5 * p.draft * weight for weight in weights]
        bending = willmore = gaussian = area = 0.0
        bending_terms: list[float] = []
        willmore_terms: list[float] = []
        gaussian_terms: list[float] = []
        area_terms: list[float] = []
        for panel in range(panels):
            low = -1.0 + 2.0 * panel / panels
            high = -1.0 + 2.0 * (panel + 1) / panels
            half = 0.5 * (high - low)
            mid = 0.5 * (high + low)
            # d xi -> dx: x = L (xi + 1) / 2, so dx = (L/2) d xi
            xi_weights = [0.5 * p.length * half * weight for weight in weights]
            for node, wx in zip(nodes, xi_weights, strict=True):
                xi = mid + half * node
                for depth, wd in zip(depth_nodes, depth_weights, strict=True):
                    b, wf, g, jac = _flank_integrand(p, xi, depth)
                    weight = wx * wd * jac
                    bending_terms.append(weight * b)
                    willmore_terms.append(weight * wf)
                    gaussian_terms.append(weight * g)
                    area_terms.append(weight)
        bending = math.fsum(bending_terms)
        willmore = math.fsum(willmore_terms)
        gaussian = math.fsum(gaussian_terms)
        area = math.fsum(area_terms)
        # Both flanks are congruent: h enters the metric and the curvatures
        # quadratically, so doubling the starboard integral is exact.
        return 2.0 * bending, 2.0 * willmore, 2.0 * gaussian, 2.0 * area

    bending, willmore, gaussian, area = integrate(order)
    fine = integrate(2 * order)
    scale = max(abs(fine[0]), abs(fine[1]), abs(fine[2]), 1e-30)
    quadrature_relative_change = max(
        abs(bending - fine[0]), abs(willmore - fine[1]), abs(gaussian - fine[2])
    ) / scale

    return WigleyFlankEnergy(
        length=p.length,
        beam=p.beam,
        draft=p.draft,
        longitudinal_exponent=p.longitudinal_exponent,
        vertical_exponent=p.vertical_exponent,
        wetted_area=area,
        bending=bending,
        willmore=willmore,
        gaussian=gaussian,
        corner_charge=4.0 * math.pi - gaussian,
        order=order,
        panels=panels,
        quadrature_relative_change=quadrature_relative_change,
    )
