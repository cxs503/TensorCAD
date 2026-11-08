"""Prismatic planing-craft hull family: constant deadrise, chine planform power law.

This is the geometry Savitsky's planing method is written for (Savitsky 1964,
"Hydrodynamic design of planing hulls", *Marine Technology* 1(1)): a bottom of
constant deadrise :math:`\\beta`, a chine whose half-breadth grows along the length,
and a transom. Savitsky needs, as input, the deadrise, the chine beam and the wetted
length/beam ratio of that surface; what he does *not* provide is the hull itself, so a
planing study that wants more than one parent form has to loft it. That is the gap this
module fills: it turns ``(L, B, T, beta, q)`` into a watertight submerged hull, and it
publishes the closed-form hydrostatics of the same parameters so the loft can be
checked against something that was not computed from the loft.

Geometry
--------
Bow at ``x = 0``, transom at ``x = L``, waterline ``z = 0``, keel ``z = -T``,
centreplane ``y = 0``. The chine half-breadth is a power law

    b(x) = (B/2) * (x/L)**q,            q > 0

so ``q = 1`` is a straight chine (a true prism with constant deadrise), ``q < 1``
hollows the bow and fills the stern, which is how real planing bottoms are drawn. The
section at ``x`` is a V of deadrise ``beta`` from the keel up to the chine at
``z_chine(x) = -T + b(x) tan(beta)`` and a vertical topside above it, so the wetted
half-breadth is

    y(x, z) = min( (z + T)/tan(beta),  b(x) )

Two regimes therefore exist inside one hull: near the bow the chine sits *below* the
waterline and the topside is wetted; near a sufficiently wide transom the chine rises
*above* the waterline and the waterline cuts the V-bottom. The wetted length

    x0 = L * (T / ((B/2) tan(beta)))**(1/q)

is where the two meet (clipped to ``[0, L]``). Both regimes are integrated in closed
form, and the mesh puts a grid line exactly on ``z_chine`` and on ``z = 0`` so neither
the chine kink nor the waterline is smeared by the sampling.

What the mesh is
----------------
The **wetted** hull only: the V-bottom and the wetted strip of topside, closed by the
planar waterline cap. That is the convention `tensorcad.hydrostatics` states, and it is
why the reported beam of a chine-dry hull is ``2 T / tan(beta)`` and not ``B`` -- the
dry part of the chine simply is not in the mesh. Freeboard, sheer, strakes, steps and
spray rails are not modelled.

Limitations, stated plainly: this is geometry and hydrostatics, not a resistance
prediction. Savitsky's lift, drag and running-trim relations are a *dynamic* problem
(they need speed, weight and a trim equilibrium) and are out of scope here; what the
module supplies is the geometry those relations consume. Nothing here claims seakeeping,
class-rule compliance or structural adequacy.
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
    ERROR,
    WARNING,
    Check,
    DimensionTarget,
    ValidationReport,
    check_dimensions,
    inspect_geometry,
    merge,
)

SAMPLING_ERROR_WARNING = 0.01
WETTED_AREA_WARNING = 0.02
MAX_STATIONS = 2001
MAX_LEVELS = 513
CHINE_SNAP_TOLERANCE = 1e-12

#: 10-point Gauss-Legendre nodes and weights on ``[-1, 1]`` (exact to degree 19).
_GL10_X = (
    -0.9739065285171717, -0.8650633666889845, -0.6794095682990244,
    -0.4333953941292472, -0.1488743389816312, 0.1488743389816312,
    0.4333953941292472, 0.6794095682990244, 0.8650633666889845,
    0.9739065285171717,
)
_GL10_W = (
    0.0666713443086881, 0.1494513491505806, 0.2190863625159820,
    0.2692667193099963, 0.2955242247147529, 0.2955242247147529,
    0.2692667193099963, 0.2190863625159820, 0.1494513491505806,
    0.0666713443086881,
)
TOPSIDE_QUADRATURE_RTOL = 1e-12
TOPSIDE_QUADRATURE_MAX_DEPTH = 24




@dataclass(frozen=True)
class PlaningHullParameters:
    """Parameters of a constant-deadrise planing bottom.

    ``draft`` is the static waterline above the keel. The family is defined for any
    draft; when ``draft > (beam/2) tan(deadrise)`` the chine is immersed at the
    transom and the topsides are wetted there, otherwise the waterline cuts the
    V-bottom (the usual light-ship case for a planing craft).
    """

    length: float = 10.0
    beam: float = 3.0
    draft: float = 0.55
    deadrise_deg: float = 18.0
    planform_exponent: float = 0.6
    nx: int = 81
    nz: int = 17
    water_density: float = 1025.0

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)

    @property
    def tan_beta(self) -> float:
        return math.tan(math.radians(self.deadrise_deg))

    @property
    def half_beam_transom(self) -> float:
        return 0.5 * self.beam

    @property
    def chine_height_transom(self) -> float:
        """Height of the transom chine above the keel, ``(B/2) tan(beta)``."""
        return self.half_beam_transom * self.tan_beta

    @property
    def chine_immersed_at_transom(self) -> bool:
        return self.draft > self.chine_height_transom

    def problems(self) -> tuple[list[str], list[str]]:
        errors: list[str] = []
        for name, lo, hi in (("length", 1e-3, 1e4), ("beam", 1e-3, 1e3), ("draft", 1e-3, 1e3)):
            value = float(getattr(self, name))
            if not math.isfinite(value) or not lo <= value <= hi:
                errors.append(f"{name} must be finite and in [{lo:g}, {hi:g}], got {value!r}")
        if not math.isfinite(self.deadrise_deg) or not 0.0 < self.deadrise_deg < 60.0:
            errors.append(f"deadrise_deg must be in (0, 60), got {self.deadrise_deg!r}")
        if not math.isfinite(self.planform_exponent) or not 0.05 <= self.planform_exponent <= 4.0:
            errors.append(
                f"planform_exponent must be in [0.05, 4.0], got {self.planform_exponent!r}"
            )
        if not 3 <= int(self.nx) <= MAX_STATIONS:
            errors.append(f"nx must be in [3, {MAX_STATIONS}], got {self.nx!r}")
        if not 3 <= int(self.nz) <= MAX_LEVELS:
            errors.append(f"nz must be in [3, {MAX_LEVELS}], got {self.nz!r}")
        if not math.isfinite(self.water_density) or not 300.0 <= self.water_density <= 2000.0:
            errors.append(f"water_density must be in [300, 2000], got {self.water_density!r}")
        warnings: list[str] = []
        if not errors:
            if self.draft > 0.5 * self.length:
                warnings.append(f"draft/length = {self.draft / self.length:.3g} looks like a barge")
            if not self.chine_immersed_at_transom and self.draft > self.half_beam_transom:
                warnings.append(
                    "chine is dry at the transom but the draft exceeds the transom "
                    "half-breadth; is the beam a chine beam?"
                )
        return errors, warnings

    def validate(self) -> None:
        errors, _ = self.problems()
        if errors:
            raise ValueError("; ".join(errors))


# --------------------------------------------------------------------------------------
# Geometry primitives, single source of truth for both the mesh and the closed forms
# --------------------------------------------------------------------------------------
def chine_half_beam(p: PlaningHullParameters, x: float) -> float:
    """``b(x) = (B/2) (x/L)**q``, clipped to the hull length."""
    if x <= 0.0:
        return 0.0
    if x >= p.length:
        return p.half_beam_transom
    return p.half_beam_transom * math.pow(x / p.length, p.planform_exponent)


def wetted_length(p: PlaningHullParameters) -> float:
    """``x0``: bow to the station where the chine reaches the waterline."""
    ratio = p.draft / (p.half_beam_transom * p.tan_beta)
    if ratio >= 1.0:
        return p.length
    if ratio <= 0.0:
        return 0.0
    return p.length * math.pow(ratio, 1.0 / p.planform_exponent)


def waterline_half_breadth(p: PlaningHullParameters, x: float) -> float:
    return min(chine_half_beam(p, x), p.draft / p.tan_beta)


def section_area(p: PlaningHullParameters, x: float) -> float:
    """Submerged area of the section at ``x``: V below the chine, wall above it."""
    tau = p.tan_beta
    b = chine_half_beam(p, x)
    chine_depth = b * tau
    if chine_depth >= p.draft:
        return p.draft * p.draft / tau
    return b * chine_depth + 2.0 * b * (p.draft - chine_depth)


def _power_integral(coefficient: float, exponent: float, lo: float, hi: float, power: int) -> float:
    """``coefficient**power * int_lo^hi x**(power*exponent) dx``."""
    order = power * exponent + 1.0
    return coefficient**power * (hi**order - lo**order) / order


def analytic_volume(p: PlaningHullParameters) -> float:
    """Volume below the waterline, both regimes integrated in closed form."""
    tau = p.tan_beta
    x0 = wetted_length(p)
    c = p.half_beam_transom / p.length**p.planform_exponent
    # Bow side (x < x0): chine below the waterline, V plus wetted topside.
    wet = 2.0 * p.draft * _power_integral(c, p.planform_exponent, 0.0, x0, 1)
    wet -= tau * _power_integral(c, p.planform_exponent, 0.0, x0, 2)
    # Stern side (x > x0): chine above the waterline, the waterline cuts the V.
    dry = p.draft * p.draft / tau * (p.length - x0)
    return wet + dry


def analytic_waterplane_area(p: PlaningHullParameters) -> float:
    x0 = wetted_length(p)
    c = p.half_beam_transom / p.length**p.planform_exponent
    chine_side = 2.0 * _power_integral(c, p.planform_exponent, 0.0, x0, 1)
    v_side = 2.0 * (p.draft / p.tan_beta) * (p.length - x0)
    return chine_side + v_side


def _moment_x_b(c: float, q: float, lo: float, hi: float, coefficient: float,
                power: int = 1) -> float:
    """``coefficient * int_lo^hi x * (c x**q)**power dx``."""
    order = power * q + 2.0
    return coefficient * c**power * (hi**order - lo**order) / order


def analytic_longitudinal_centre_buoyancy(p: PlaningHullParameters) -> float:
    """``LCB`` measured from the bow (``x = 0``)."""
    tau = p.tan_beta
    x0 = wetted_length(p)
    c = p.half_beam_transom / p.length**p.planform_exponent
    moment = _moment_x_b(c, p.planform_exponent, 0.0, x0, 2.0 * p.draft, 1)
    moment -= _moment_x_b(c, p.planform_exponent, 0.0, x0, tau, 2)
    moment += p.draft * p.draft / tau * 0.5 * (p.length**2 - x0**2)
    return moment / analytic_volume(p)


def analytic_longitudinal_centre_flotation(p: PlaningHullParameters) -> float:
    """``LCF`` from the bow, on the waterplane the volume was integrated over."""
    x0 = wetted_length(p)
    c = p.half_beam_transom / p.length**p.planform_exponent
    moment = _moment_x_b(c, p.planform_exponent, 0.0, x0, 1.0, 1)
    moment += (p.draft / p.tan_beta) * 0.5 * (p.length**2 - x0**2)
    return moment / (analytic_waterplane_area(p) / 2.0)


def analytic_vertical_centre_buoyancy(p: PlaningHullParameters) -> float:
    """``VCB`` in hull coordinates (negative below the waterline)."""
    tau = p.tan_beta
    x0 = wetted_length(p)
    c = p.half_beam_transom / p.length**p.planform_exponent
    moment = -p.draft**2 * _power_integral(c, p.planform_exponent, 0.0, x0, 1)
    moment += p.draft * tau * _power_integral(c, p.planform_exponent, 0.0, x0, 2)
    moment -= (1.0 / 3.0) * tau**2 * _power_integral(c, p.planform_exponent, 0.0, x0, 3)
    moment += -p.draft**3 / (3.0 * tau) * (p.length - x0)
    return moment / analytic_volume(p)


def analytic_wetted_area(p: PlaningHullParameters,
                         rtol: float = TOPSIDE_QUADRATURE_RTOL,
                         *, samples: int | None = None) -> float:
    """Wetted surface: V-bottom panels + wetted topside strip + transom.

    The V-bottom is closed form: its surface element is ``b sec(beta) dx`` independent
    of the chine slope, because the ruling lines are straight from keel to chine. The
    wetted topside is **not** its vertical projection: the side follows the chine curve,
    so its element is ``(T - b tan(beta)) sqrt(1 + b'(x)**2) dx``. That integral has no
    elementary closed form once ``q < 1`` (``b'`` diverges at a knife bow), so it is
    evaluated by adaptive Gauss-Legendre quadrature in the substituted variable
    ``x = x0 t**e`` with ``e = 2/q`` (``q < 2``) or ``e = 1`` (``q >= 2``), which removes
    the endpoint singularity without creating one at ``t = 0``; see
    :func:`_wetted_topside_area`. The transom is closed form.

    ``samples`` is the fixed Simpson interval count this integral used before the rule
    became adaptive; it is accepted for call compatibility and ignored.
    """
    tau = p.tan_beta
    x0 = wetted_length(p)
    bottom = analytic_waterplane_area(p) * math.sqrt(1.0 + tau * tau)
    return bottom + _wetted_topside_area(p, x0, rtol) + section_area(p, p.length)


def _gauss10(f, lo: float, hi: float) -> float:
    centre = 0.5 * (lo + hi)
    half = 0.5 * (hi - lo)
    total = 0.0
    for node, weight in zip(_GL10_X, _GL10_W):
        total += weight * f(centre + half * node)
    return half * total


def _adaptive_gauss10(f, lo: float, hi: float, rtol: float, depth: int = 0) -> float:
    """Adaptive 10-point Gauss-Legendre on ``[lo, hi]``, relative tolerance ``rtol``.

    Bisection is driven by the difference between one panel and its two halves, which
    for a 10-point rule overestimates the error by orders of magnitude and so stops
    early on the smooth parts of the integrand while still resolving a boundary layer
    (``q = 0.05`` puts a ``t**39`` shape against ``t = 1``).
    """
    whole = _gauss10(f, lo, hi)
    middle = 0.5 * (lo + hi)
    split = _gauss10(f, lo, middle) + _gauss10(f, middle, hi)
    if depth >= TOPSIDE_QUADRATURE_MAX_DEPTH:
        return split
    if abs(split - whole) <= rtol * abs(split):
        return split
    return _adaptive_gauss10(f, lo, middle, rtol, depth + 1) + _adaptive_gauss10(
        f, middle, hi, rtol, depth + 1
    )


def _wetted_topside_area(p: PlaningHullParameters, x0: float,
                         rtol: float = TOPSIDE_QUADRATURE_RTOL) -> float:
    """``2 int_0^x0 (T - b tan(beta)) sqrt(1 + b'^2) dx``, adaptively quadratured.

    Evaluated in the substituted variable ``x = x0 t**e`` so that neither the integrand
    nor its factors may leave double range. Three numerical traps are avoided:

    * ``b'(x) = c q x**(q-1)`` diverges at a knife bow (``q < 1``). ``e = 2/q`` turns
      ``b`` into ``b(x0) t**2`` and the arc element into
      ``sqrt((dx/dt)**2 + (db/dt)**2)`` with ``db/dt = 2 b(x0) t``, i.e. a smooth
      integrand. For ``q >= 2`` that substitution would *introduce* a ``t**(2/q-1)``
      singularity at ``t = 0``, so ``e = 1`` there instead (``b'`` already vanishes at
      the stem for ``q > 1``).
    * ``sqrt(1 + b'**2) * dx`` is formed as ``hypot(dx/dt, db/dt)`` rather than as
      ``sqrt(1 + slope**2) * dx/dt``. At ``q = 0.05``, ``beta >= 30 deg`` and a light
      draft the slope reaches ``~1e21``; squaring it overflowed to ``inf`` and the
      whole wetted area came back ``inf``.
    * the integral is adaptive instead of a fixed 4000-interval Simpson rule: with a
      small ``q`` the arc element varies by ``e**(2/q)`` across ``[0, 1]``, so the
      uniform grid was 0.06 % short at ``q = 0.05`` and 0.85 % short at ``q >= 2``
      (see ``docs/planing_numerical_robustness.md``).
    """
    if x0 <= 0.0:
        return 0.0
    q = p.planform_exponent
    tau = p.tan_beta
    exponent = 2.0 / q if q < 2.0 else 1.0
    b_x0 = chine_half_beam(p, x0)
    power = exponent * q

    def integrand(t: float) -> float:
        height = p.draft - tau * b_x0 * t**power
        if height <= 0.0:
            return 0.0
        dx_dt = x0 * exponent * (1.0 if exponent == 1.0 else t ** (exponent - 1.0))
        db_dt = b_x0 * power * (1.0 if power == 1.0 else t ** (power - 1.0))
        return 2.0 * height * math.hypot(dx_dt, db_dt)

    if not rtol > 0.0:
        raise ValueError(f"rtol must be positive, got {rtol!r}")
    return _adaptive_gauss10(integrand, 0.0, 1.0, rtol)


def analytic_planing_hydrostatics(p: PlaningHullParameters) -> dict[str, float]:
    """Closed-form hydrostatics of the family, independent of the mesh."""
    volume = analytic_volume(p)
    awp = analytic_waterplane_area(p)
    transom_area = section_area(p, p.length)
    waterline_beam = 2.0 * waterline_half_breadth(p, p.length)
    return {
        "displacement_volume": volume,
        "waterplane_area": awp,
        "wetted_surface_area": analytic_wetted_area(p),
        "transom_section_area": transom_area,
        "maximum_section_area": transom_area,
        "waterline_beam": waterline_beam,
        "wetted_length": wetted_length(p),
        "longitudinal_centre_buoyancy": analytic_longitudinal_centre_buoyancy(p),
        "vertical_centre_buoyancy": analytic_vertical_centre_buoyancy(p),
        "longitudinal_centre_flotation": analytic_longitudinal_centre_flotation(p),
        "vertical_centre_buoyancy_above_keel": analytic_vertical_centre_buoyancy(p) + p.draft,
        "block_coefficient": volume / (p.length * p.beam * p.draft),
        "prismatic_coefficient": volume / (p.length * transom_area) if transom_area else float("nan"),
        "waterplane_coefficient": awp / (p.length * p.beam),
        "chine_height_transom": p.chine_height_transom,
        "chine_immersed_at_transom": float(p.chine_immersed_at_transom),
    }


# --------------------------------------------------------------------------------------
# Mesh
# --------------------------------------------------------------------------------------
def _station_levels(p: PlaningHullParameters, x: float) -> list[float]:
    """``nz`` z-levels from the keel to the waterline, with ``z_chine`` exactly on one.

    The chine is a kink of the section, and a kink that lands between two grid lines is
    rounded; snapping the nearest interior level onto it keeps the loft exact in ``z``
    for both regimes, and it keeps the ring size fixed so the loft stays uniform.
    """
    levels = [-p.draft + p.draft * j / (p.nz - 1) for j in range(p.nz)]
    levels[0] = -p.draft
    levels[-1] = 0.0
    chine_z = -p.draft + chine_half_beam(p, x) * p.tan_beta
    if -p.draft < chine_z < 0.0:
        nearest = min(range(1, p.nz - 1), key=lambda j: abs(levels[j] - chine_z))
        if abs(levels[nearest] - chine_z) > CHINE_SNAP_TOLERANCE:
            levels[nearest] = chine_z
            levels.sort()
    return levels


def _ring_points(p: PlaningHullParameters, x: float) -> list[tuple[float, float, float]]:
    b = chine_half_beam(p, x)
    tau = p.tan_beta
    half = [(min((z + p.draft) / tau, b), z) for z in _station_levels(p, x)]
    ring = [(x, y, z) for (y, z) in reversed(half)]
    ring += [(x, -y, z) for (y, z) in half[1:]]
    return ring


def build_planing_mesh(p: PlaningHullParameters) -> Mesh:
    """Watertight wetted hull: lofted sections, waterline cap, triangulated transom."""
    p.validate()
    builder = MeshBuilder(units=UNITS, coordinate_frame=COORDINATE_FRAME)
    rings: list[list[int]] = []
    for i in range(p.nx):
        x = p.length * i / (p.nx - 1)
        rings.append([builder.add_vertex(point) for point in _ring_points(p, x)])
    size = len(rings[0])
    for i in range(p.nx - 1):
        for j in range(size - 1):
            builder.add_quad(rings[i][j], rings[i + 1][j], rings[i + 1][j + 1], rings[i][j + 1])
        builder.add_quad(rings[i][size - 1], rings[i + 1][size - 1], rings[i + 1][0], rings[i][0])
    # Transom: the last ring is a convex polygon in the plane x = L (the wetted section),
    # so a fan from an interior point triangulates it without degenerate facets.
    centre = builder.add_vertex((p.length, 0.0, -p.draft / 3.0))
    last = rings[-1]
    for j in range(size):
        builder.add_triangle(centre, last[j], last[(j + 1) % size])
    return meshlib.orient_outward(builder.build()).mesh


# --------------------------------------------------------------------------------------
# Model
# --------------------------------------------------------------------------------------
def _dimension_target(p: PlaningHullParameters) -> DimensionTarget:
    return DimensionTarget(
        length=p.length,
        beam=2.0 * waterline_half_breadth(p, p.length),
        draft=p.draft,
    )


def _relative(measured: float, reference: float) -> float:
    if not reference:
        return float("nan")
    return (measured - reference) / reference


def _hydrostatic_checks(p: PlaningHullParameters, hydro: Hydrostatics) -> list[Check]:
    reference = analytic_planing_hydrostatics(p)
    checks: list[Check] = []
    for name, key, measured, tolerance in (
        ("planing_volume", "displacement_volume", hydro.displacement_volume, SAMPLING_ERROR_WARNING),
        ("planing_waterplane", "waterplane_area", hydro.waterplane_area, SAMPLING_ERROR_WARNING),
        ("planing_wetted_area", "wetted_surface_area", hydro.wetted_surface_area, WETTED_AREA_WARNING),
        ("planing_lcb", "longitudinal_centre_buoyancy", hydro.longitudinal_centre_buoyancy, 0.02),
        ("planing_vcb", "vertical_centre_buoyancy", hydro.vertical_centre_buoyancy, 0.02),
    ):
        ref = reference[key]
        if not math.isfinite(ref) or ref == 0.0:
            checks.append(
                Check(
                    name,
                    False,
                    f"closed form is not a usable reference ({ref!r}); the mesh value "
                    f"{measured:.6f} cannot be validated against it",
                    severity=ERROR,
                )
            )
            continue
        if name == "planing_lcb":
            relative = _relative(measured, ref) if ref else _relative(measured, 1.0)
        elif name == "planing_vcb":
            relative = _relative(measured, ref) if ref else _relative(measured, 1.0)
        else:
            relative = _relative(measured, ref)
        ok = math.isfinite(relative) and abs(relative) <= tolerance
        checks.append(
            Check(
                name,
                ok,
                f"mesh {measured:.6f} vs closed form {ref:.6f} = {relative:+.4%} "
                f"(tolerance {tolerance:.1%})",
                severity=WARNING,
            )
        )
    return checks


@dataclass(frozen=True)
class PlaningHullModel:
    parameters: PlaningHullParameters
    mesh: Mesh
    validation: ValidationReport = field(default_factory=ValidationReport)
    hydrostatics: Hydrostatics | None = None

    @property
    def vertices(self):  # noqa: ANN201 - mirrors the other hull family models
        return self.mesh.vertices

    @property
    def faces(self):  # noqa: ANN201
        return self.mesh.faces

    @property
    def units(self):  # noqa: ANN201
        return self.mesh.units

    @property
    def coordinate_frame(self):  # noqa: ANN201
        return self.mesh.coordinate_frame

    def metadata(self) -> dict[str, Any]:
        data: dict[str, Any] = {
            "schema_version": "2.0",
            "model_type": "constant_deadrise_planing_hull",
            "parameters": self.parameters.as_dict(),
            "units": self.units,
            "coordinate_frame": self.coordinate_frame,
            "generator": {
                "kind": "prismatic_planing_loft",
                "chine_half_beam": "(B/2) (x/L)^q",
                "section": "V of deadrise beta to the chine, vertical topside above",
                "wetted_length": wetted_length(self.parameters),
                "topology": "lofted wetted sections, waterline cap, triangulated transom",
                "backend": "tensorcad.geometry (standard library)",
            },
            "mesh": self.mesh.metadata(),
            "validation": self.validation.to_dict(),
            "engineering_status": "geometry and hydrostatics only; no resistance, "
            "no seakeeping, no rule compliance",
        }
        if self.hydrostatics is not None:
            data["hydrostatics"] = self.hydrostatics.to_dict()
        data["analytic_reference"] = analytic_planing_hydrostatics(self.parameters)
        return data


def generate_planing_hull(parameters: PlaningHullParameters | None = None,
                          **overrides: Any) -> PlaningHullModel:
    if parameters is None:
        p = PlaningHullParameters(**overrides)
    elif overrides:
        p = PlaningHullParameters(**{**parameters.as_dict(), **overrides})
    else:
        p = parameters
    p.validate()
    mesh = build_planing_mesh(p)
    report = merge(inspect_geometry(mesh), check_dimensions(mesh, _dimension_target(p)))
    if report.errors:
        raise ValueError("generated planing hull failed validation: " + "; ".join(report.errors))
    hydrostatics = compute_hydrostatics(mesh, 0.0, water_density=p.water_density, closed=True)
    report = merge(list(report.checks), _hydrostatic_checks(p, hydrostatics))
    if report.errors:
        # A non-finite or zero closed form is not sampling noise: it means the reference
        # the mesh is validated against is broken, so the model must not be handed back
        # with ``inf`` in its metadata (see docs/planing_numerical_robustness.md).
        raise ValueError("generated planing hull failed validation: " + "; ".join(report.errors))
    return PlaningHullModel(parameters=p, mesh=mesh, validation=report,
                            hydrostatics=hydrostatics)


def planing_relative_errors(p: PlaningHullParameters,
                            model: PlaningHullModel | None = None) -> dict[str, float]:
    """Mesh-minus-closed-form relative errors of every compared quantity."""
    if model is None:
        model = generate_planing_hull(p)
    if model.hydrostatics is None:
        raise ValueError("model has no hydrostatics")
    reference = analytic_planing_hydrostatics(p)
    hydro = model.hydrostatics
    pairs = {
        "displacement_volume": (hydro.displacement_volume, reference["displacement_volume"]),
        "waterplane_area": (hydro.waterplane_area, reference["waterplane_area"]),
        "wetted_surface_area": (hydro.wetted_surface_area, reference["wetted_surface_area"]),
        "longitudinal_centre_buoyancy": (
            hydro.longitudinal_centre_buoyancy,
            reference["longitudinal_centre_buoyancy"],
        ),
        "vertical_centre_buoyancy": (
            hydro.vertical_centre_buoyancy,
            reference["vertical_centre_buoyancy"],
        ),
        "longitudinal_centre_flotation": (
            hydro.longitudinal_centre_flotation,
            analytic_longitudinal_centre_flotation(p),
        ),
    }
    return {name: _relative(measured, ref) for name, (measured, ref) in pairs.items()}


def simpson_section_integrals(p: PlaningHullParameters, intervals: int = 4096) -> dict[str, float]:
    """Independent quadrature of the section integrals, as a second reference route.

    The closed forms are polynomial integrals on paper; this evaluates the *sections*
    (the same :func:`section_area` the loft is built from) with composite Simpson, so a
    mistake in the algebra shows up as a disagreement between the two routes instead of
    being validated by the mesh that inherited it.
    """
    if intervals % 2:
        intervals += 1
    h = p.length / intervals

    def integrate(function) -> float:  # noqa: ANN001 - local helper
        total = function(0.0) + function(p.length)
        for i in range(1, intervals):
            weight = 4.0 if i % 2 else 2.0
            total += weight * function(i * h)
        return total * h / 3.0

    volume = integrate(lambda x: section_area(p, x))
    lcb = integrate(lambda x: x * section_area(p, x)) / volume

    def vertical_moment(x: float) -> float:
        # Section first moment about the waterline, from the same half-breadth rule the
        # mesh uses: integrate y(z) over z with weight z, doubled for both sides.
        tau = p.tan_beta
        b = chine_half_beam(p, x)
        levels = _station_levels(p, x)
        half = [min((z + p.draft) / tau, b) for z in levels]
        moment = 0.0
        for j in range(len(levels) - 1):
            z0, z1 = levels[j], levels[j + 1]
            y0, y1 = half[j], half[j + 1]
            trapezoid = 0.5 * (y0 + y1) * (z1 - z0)
            centroid = (z0 * (2 * y0 + y1) + z1 * (y0 + 2 * y1)) / (3.0 * (y0 + y1)) \
                if (y0 + y1) else 0.0
            moment += 2.0 * trapezoid * centroid
        return moment

    vcb = integrate(vertical_moment) / volume
    return {"displacement_volume": volume, "longitudinal_centre_buoyancy": lcb,
            "vertical_centre_buoyancy": vcb}


def estimate_planing_discretisation_error(p: PlaningHullParameters,
                                          max_triangles: int = 400_000) -> dict[str, Any]:
    reference = analytic_planing_hydrostatics(p)
    results: list[dict[str, Any]] = []
    nx, nz = p.nx, p.nz
    while True:
        trial = PlaningHullParameters(
            length=p.length, beam=p.beam, draft=p.draft,
            deadrise_deg=p.deadrise_deg, planform_exponent=p.planform_exponent,
            nx=nx, nz=nz, water_density=p.water_density,
        )
        model = generate_planing_hull(trial)
        hydro = model.hydrostatics
        assert hydro is not None
        results.append({
            "nx": nx,
            "nz": nz,
            "triangles": len(model.mesh.faces),
            "volume_m3": hydro.displacement_volume,
            "volume_error": _relative(hydro.displacement_volume, reference["displacement_volume"]),
            "waterplane_area_m2": hydro.waterplane_area,
            "waterplane_error": _relative(hydro.waterplane_area, reference["waterplane_area"]),
            "wetted_area_m2": hydro.wetted_surface_area,
            "wetted_area_error": _relative(hydro.wetted_surface_area,
                                           reference["wetted_surface_area"]),
        })
        if len(model.mesh.faces) > max_triangles or nx > MAX_STATIONS // 2:
            break
        nx = 2 * nx - 1
        nz = 2 * nz - 1
    return {"reference": reference, "convergence": results}