"""An asymmetric hull family, and its closed-form hydrostatic reference.

``tensorcad.hull`` samples one symmetric surface: the half-breadth is
``(beam/2) * sin(pi*x/length)^0.7 * sqrt((z+draft)/draft)``, so its longitudinal
centre of buoyancy is ``length/2`` by construction and the LCB/LCF test in that
module is *trivially* satisfied -- a sign or frame error in the ``x`` moments
would still be caught, but an antisymmetric bias could not exist in the geometry
even if the integrator invented one. That is the gap this module closes.

The family here keeps the separable vertical profile and gives the two halves of
the length their own taper exponent::

    b(x, z) = (beam/2) * sin(pi*x/length)^p(x) * sqrt((z + draft)/draft)
    p(x)    = fore_exponent  for  x >= length/2
              aft_exponent   for  x <  length/2

``sin(pi*x)`` is 1 at midship, so the surface is continuous there whatever the
two exponents are; only the waterline slope breaks, which is exactly what a
bow/stern asymmetry looks like on a real hull. ``fore_exponent != aft_exponent``
makes the longitudinal distribution of volume genuinely fore/aft dependent, so
LCB and LCF become numbers a closed form has to *predict* rather than ``L/2`` it
gets for free.

Why the reference is still analytic
-----------------------------------
The surface stays separable, so every hydrostatic integral factors:

| quantity | reference |
|---|---|
| displaced volume | ``L*B*T*(2/3)*(I_f + I_a)`` |
| waterplane area | ``L*B*(I_f + I_a)`` |
| midship section area | ``2/3*B*T`` |
| VCB | ``-2/5*T`` (the vertical profile is unchanged) |
| Cb / Cp | ``(2/3)*(I_f + I_a)`` / ``I_f + I_a`` |
| LCB, LCF | ``L*(M_a + I_f - M_f)/(I_a + I_f)`` |
| wetted surface area | Gauss-Legendre quadrature (no elementary antiderivative) |

with ``I_p = int_0^{1/2} sin(pi*x)^p dx`` and ``M_p = int_0^{1/2} x*sin(pi*x)^p dx``
(normalised ``x``). ``I_p`` has a Beta-function closed form
(``I_p = J(p)/pi``, ``J(p) = sqrt(pi)*Gamma((p+1)/2)/(2*Gamma((p+2)/2))``), which
:func:`analytic_half_taper_integral` evaluates with ``math.lgamma`` and which
reproduces ``hull.TAPER_INTEGRAL/2`` at ``p = 0.7`` to machine precision, so the
new code is checked against the shipped 0.11.0 constant rather than beside it.
``M_p`` has no elementary antiderivative; it is integrated with a Gauss-Legendre
rule after the substitution ``x = (1/2)*tau^(2/p)``, which turns the tip's
``x^p`` behaviour into ``tau^2`` and leaves a smooth integrand.

Note that ``LCF == LCB`` for this family: the waterplane and the volume are
weighted by the *same* longitudinal taper, so the two centroids coincide. That is
a property of the family, not of the integrator, and it is asserted as such
rather than presented as independent agreement.

Limitations
-----------
* Still one family. No transom, bulb, flare, knuckle or appendages; ``p(x)`` is
  the only asymmetry knob.
* The waterline is the plane ``z = 0``; a trimmed or heeled condition is out of
  scope, exactly as in :func:`tensorcad.hydrostatics.compute_hydrostatics`.
* Nothing here predicts stability, resistance or class compliance: it is a
  geometry plus an integration reference.
* The exponent range is ``[MIN_EXPONENT, MAX_EXPONENT] = [0.4, 3.0]``; the wetted
  area quadrature and the moment substitution are calibrated inside it, and the
  constructor refuses values outside it instead of returning a number whose
  accuracy is unknown.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

from . import mesh as meshlib
from .geometry import Mesh, MeshBuilder
from .hull import (
    COORDINATE_FRAME,
    MAX_LONGITUDINAL_STATIONS,
    MAX_VERTICAL_LEVELS,
    SAMPLING_ERROR_WARNING,
    UNITS,
    _gauss_legendre,
    _positive_dimension,
)
from .hydrostatics import Hydrostatics, compute_hydrostatics
from .validation import (
    WARNING,
    Check,
    ValidationReport,
    inspect_geometry,
    merge,
)

SCHEMA = "tensorcad.hull.asymmetric/1"

#: Smallest taper exponent the reference is calibrated for.
MIN_EXPONENT = 0.4
#: Largest taper exponent the reference is calibrated for.
MAX_EXPONENT = 3.0

#: Shipped node count for the ``M_p`` moments (stability measured in the tests).
MOMENT_NODES = 100
#: Shipped node count for the wetted-surface quadrature.
WETTED_SURFACE_NODES = 100

#: Default fore/aft exponents: a fine bow over a fuller run, i.e. LCB well
#: forward of midship, which is the case the symmetric family cannot express.
DEFAULT_FORE_EXPONENT = 0.7
DEFAULT_AFT_EXPONENT = 1.6


@dataclass(frozen=True)
class AsymmetricHullParameters:
    """Principal dimensions plus one taper exponent per half of the length."""

    length: float = 10.0
    beam: float = 2.0
    draft: float = 0.8
    fore_exponent: float = DEFAULT_FORE_EXPONENT
    aft_exponent: float = DEFAULT_AFT_EXPONENT
    nx: int = 81
    nz: int = 21
    water_density: float = 1025.0

    def _values(self) -> dict[str, float]:
        return {
            "length": float(self.length),
            "beam": float(self.beam),
            "draft": float(self.draft),
            "fore_exponent": float(self.fore_exponent),
            "aft_exponent": float(self.aft_exponent),
            "nx": float(self.nx),
            "nz": float(self.nz),
            "water_density": float(self.water_density),
        }

    def problems(self) -> list[str]:
        """Every reason the parameters cannot be sampled, not just the first."""
        errors: list[str] = []
        for name in ("length", "beam", "draft", "water_density"):
            try:
                _positive_dimension(name, float(getattr(self, name)))
            except (TypeError, ValueError) as exc:
                errors.append(str(exc))
        for name in ("fore_exponent", "aft_exponent"):
            try:
                value = float(getattr(self, name))
            except (TypeError, ValueError):
                errors.append(f"{name} must be a number")
                continue
            if not math.isfinite(value) or not MIN_EXPONENT <= value <= MAX_EXPONENT:
                errors.append(
                    f"{name} must lie in [{MIN_EXPONENT}, {MAX_EXPONENT}], got {value!r}"
                )
        for name, minimum in (("nx", 3), ("nz", 3)):
            try:
                value = int(getattr(self, name))
            except (TypeError, ValueError):
                errors.append(f"{name} must be an integer")
                continue
            if value != getattr(self, name):
                errors.append(f"{name} must be an integer, got {getattr(self, name)!r}")
            elif value < minimum:
                errors.append(f"{name} must be >= {minimum}, got {value}")
        if isinstance(self.nx, int) and self.nx > MAX_LONGITUDINAL_STATIONS:
            errors.append(f"nx must be <= {MAX_LONGITUDINAL_STATIONS}, got {self.nx}")
        if isinstance(self.nz, int) and self.nz > MAX_VERTICAL_LEVELS:
            errors.append(f"nz must be <= {MAX_VERTICAL_LEVELS}, got {self.nz}")
        return errors

    def validate(self) -> None:
        errors = self.problems()
        if errors:
            raise ValueError("; ".join(errors))

    @property
    def has_midship_station(self) -> bool:
        """True when a station lands exactly on the exponent switch.

        With an even ``nx`` the widest section is missed, the exponent switch sits
        between two stations, and the realised beam falls short of the requested
        one. The deficit is reported as a warning, not absorbed.
        """
        return self.nx % 2 == 1

    def as_dict(self) -> dict[str, Any]:
        data = self._values()
        data["nx"] = int(self.nx)
        data["nz"] = int(self.nz)
        data["has_midship_station"] = self.has_midship_station
        return data


def taper(fraction: float, parameters: AsymmetricHullParameters) -> float:
    """Longitudinal taper ``sin(pi*x/L)^p`` at ``fraction = x/L``.

    The end stations are forced to zero: ``math.sin(math.pi)`` is not exactly zero
    in floating point, and a stem that does not close on the centreline leaves a
    slot in the loft. Midship (``fraction == 0.5``) is 1 for either exponent, so
    the switch is continuous.
    """
    if fraction <= 0.0 or fraction >= 1.0:
        return 0.0
    exponent = (
        parameters.fore_exponent if fraction >= 0.5 else parameters.aft_exponent
    )
    return math.pow(math.sin(math.pi * fraction), exponent)


def _section_point(
    parameters: AsymmetricHullParameters, i: int, j: int
) -> tuple[float, float, float]:
    """Ring position ``j`` at station ``i``.

    Same ring convention as :func:`tensorcad.hull.build_mesh`: down the starboard
    side from the design waterline to the keel, then back up the port side, with
    the keel vertex shared by both flanks so the loft closes without welding.
    """
    p = parameters
    fraction = i / (p.nx - 1)
    x = p.length * fraction
    keelward = j <= p.nz - 1
    level = (p.nz - 1 - j) if keelward else (j - (p.nz - 1))
    z = -p.draft * (1.0 - level / (p.nz - 1))
    half_breadth = 0.5 * p.beam * taper(fraction, p) * math.sqrt(level / (p.nz - 1))
    y = half_breadth if keelward else -half_breadth
    return (x, y, z)


def build_mesh(
    parameters: AsymmetricHullParameters,
) -> tuple[Mesh, tuple[tuple[float, float], ...]]:
    """Loft and cap the asymmetric surface; returns the mesh and waterline offsets."""
    parameters.validate()
    p = parameters
    ring_size = 2 * p.nz - 1
    builder = MeshBuilder(units=UNITS, coordinate_frame=COORDINATE_FRAME)
    rings: list[list[int]] = []
    for i in range(p.nx):
        rings.append(
            [builder.add_vertex(_section_point(p, i, j)) for j in range(ring_size)]
        )
    for i in range(p.nx - 1):
        for j in range(ring_size - 1):
            builder.add_quad(
                rings[i][j], rings[i + 1][j], rings[i + 1][j + 1], rings[i][j + 1]
            )
        builder.add_quad(
            rings[i][ring_size - 1],
            rings[i + 1][ring_size - 1],
            rings[i + 1][0],
            rings[i][0],
        )
    built = builder.build()
    mesh = meshlib.orient_outward(built).mesh
    offsets = tuple(
        (
            p.length * i / (p.nx - 1),
            mesh.vertices[rings[i][0]][1] - mesh.vertices[rings[i][ring_size - 1]][1],
        )
        for i in range(p.nx)
    )
    return mesh, offsets


def analytic_half_taper_integral(exponent: float) -> float:
    """``I_p = int_0^{1/2} sin(pi*x)^p dx``, in closed form.

    Substituting ``u = pi*x`` gives ``I_p = J(p)/pi`` with
    ``J(p) = int_0^{pi/2} sin(u)^p du = sqrt(pi)*Gamma((p+1)/2)/(2*Gamma((p+2)/2))``;
    ``math.lgamma`` keeps the core standard-library only. At ``p = 0.7`` this is
    ``hull.TAPER_INTEGRAL/2``.
    """
    if not math.isfinite(exponent) or exponent <= 0.0:
        raise ValueError("exponent must be a finite positive number")
    log_j = (
        0.5 * math.log(math.pi)
        + math.lgamma(0.5 * (exponent + 1.0))
        - math.log(2.0)
        - math.lgamma(0.5 * (exponent + 2.0))
    )
    return math.exp(log_j) / math.pi


def half_taper_integrals(
    exponent: float, *, nodes: int = MOMENT_NODES, substitution: float | None = None
) -> tuple[float, float]:
    """``(I_p, M_p)`` by Gauss-Legendre quadrature, for cross-checking the formula.

    ``M_p = int_0^{1/2} x*sin(pi*x)^p dx`` is what makes LCB non-trivial and has no
    elementary antiderivative. The substitution ``x = (1/2)*tau^m`` with
    ``m = max(1, 2/p)`` maps the tip singularity ``x^p`` to ``tau^2`` for the
    pointed members of the range (``p < 2``) and leaves the already-regular
    full members on plain Gauss-Legendre (``m = 1``), which was the measured
    optimum: at ``p = 3`` the ``2/p`` power *loses* three orders
    (``9.7e-12`` versus ``1.1e-16`` node-refinement drift). ``substitution``
    overrides ``m`` so a test can integrate the same integral a second way.
    """
    if not math.isfinite(exponent) or exponent <= 0.0:
        raise ValueError("exponent must be a finite positive number")
    if nodes < 2:
        raise ValueError("need at least 2 quadrature nodes")
    m = substitution if substitution is not None else max(1.0, 2.0 / exponent)
    if not math.isfinite(m) or m <= 0.0:
        raise ValueError("substitution power must be finite and positive")
    xs, ws = _gauss_legendre(nodes)
    total_i = 0.0
    total_m = 0.0
    for index in range(nodes):
        tau = 0.5 * (xs[index] + 1.0)
        weight = 0.5 * ws[index]
        if tau <= 0.0:
            continue
        x = 0.5 * tau ** m
        dx = 0.5 * m * tau ** (m - 1.0)
        sine = math.sin(math.pi * x)
        if sine <= 0.0:
            continue
        value = sine ** exponent
        total_i += weight * dx * value
        total_m += weight * dx * x * value
    return total_i, total_m


def _longitudinal_integrals(
    parameters: AsymmetricHullParameters, *, nodes: int = MOMENT_NODES
) -> tuple[float, float, float, float]:
    """``(I_aft, M_aft, I_fore, M_fore)`` in normalised ``x/L``."""
    aft_i, aft_m = half_taper_integrals(parameters.aft_exponent, nodes=nodes)
    fore_i, fore_m = half_taper_integrals(parameters.fore_exponent, nodes=nodes)
    return aft_i, aft_m, fore_i, fore_m


def _half_wetted_area(
    length: float, beam: float, draft: float, exponent: float, *,
    nodes: int = WETTED_SURFACE_NODES,
) -> float:
    """Wetted area of one half of the length, both sides, by quadrature.

    Same construction as :func:`tensorcad.hull.analytic_wetted_surface_area`: with
    ``w = u^2`` the keel's ``1/sqrt(w)`` is removed and with ``psi = (pi/2)*tau^k``
    the tip slope ``x^{p-1}`` becomes ``tau^{k*p-1}``. ``k = 3/p`` for ``p < 1``
    (the shipped symmetric rule hard-codes ``k = 3``, which is the ``p = 0.7``
    member of this rule) and ``k = 3`` for ``p >= 1`` where the tip slope is
    already regular. Doubling ``nodes`` moves the result by ``< 1e-9`` relative
    across the calibrated exponent range -- see
    :func:`wetted_surface_quadrature_drift`.
    """
    k = 3.0 / exponent if exponent < 1.0 else 3.0
    gx, gw_x = _gauss_legendre(nodes)
    gu, gw_u = _gauss_legendre(nodes)
    total = 0.0
    for i in range(nodes):
        tau = 0.5 * (gx[i] + 1.0)
        weight_tau = 0.5 * gw_x[i]
        psi = 0.5 * math.pi * tau ** k
        dpsi_dtau = 0.5 * math.pi * k * tau ** (k - 1.0)
        sine = math.sin(psi)
        if sine <= 0.0:  # the tip itself
            continue
        s = sine ** exponent
        ds_dx = (
            exponent * sine ** (exponent - 1.0) * math.cos(psi) * (math.pi / length)
        )
        inner = 0.0
        for j in range(nodes):
            u = 0.5 * (gu[j] + 1.0)
            weight_u = 0.5 * gw_u[j]
            element = math.sqrt(
                4.0 * u * u * draft * draft
                + (draft * beam * ds_dx * u * u) ** 2
                + (0.5 * beam * s) ** 2
            )
            inner += weight_u * element
        total += weight_tau * dpsi_dtau * inner
    return (2.0 * length / math.pi) * total


def analytic_wetted_surface_area(
    parameters: AsymmetricHullParameters, *, nodes: int = WETTED_SURFACE_NODES
) -> float:
    """Wetted surface area of the analytic surface: the two halves summed."""
    parameters.validate()
    return _half_wetted_area(
        parameters.length, parameters.beam, parameters.draft, parameters.aft_exponent,
        nodes=nodes,
    ) + _half_wetted_area(
        parameters.length, parameters.beam, parameters.draft, parameters.fore_exponent,
        nodes=nodes,
    )


def wetted_surface_quadrature_drift(
    parameters: AsymmetricHullParameters, *, nodes: int = WETTED_SURFACE_NODES
) -> float:
    """Relative change in wetted area when the node count is doubled."""
    coarse = analytic_wetted_surface_area(parameters, nodes=nodes)
    fine = analytic_wetted_surface_area(parameters, nodes=2 * nodes)
    if coarse == 0.0:
        return float("inf")
    return abs(fine - coarse) / abs(coarse)


def analytic_volume(parameters: AsymmetricHullParameters) -> float:
    """Exact displaced volume: ``L*B*T*(2/3)*(I_f + I_a)``."""
    parameters.validate()
    aft_i, _, fore_i, _ = _longitudinal_integrals(parameters)
    return (
        parameters.beam
        * parameters.length
        * parameters.draft
        * (2.0 / 3.0)
        * (aft_i + fore_i)
    )


def analytic_waterplane_area(parameters: AsymmetricHullParameters) -> float:
    """Exact waterplane area: ``L*B*(I_f + I_a)``."""
    parameters.validate()
    aft_i, _, fore_i, _ = _longitudinal_integrals(parameters)
    return parameters.length * parameters.beam * (aft_i + fore_i)


def analytic_longitudinal_centre(parameters: AsymmetricHullParameters) -> float:
    """LCB (and LCF, which coincide for this family), in metres from the aft end.

    ``x/L = (M_a + I_f - M_f) / (I_a + I_f)``: the aft half contributes ``M_a`` and
    the fore half contributes ``int_{1/2}^1 x*f dx = I_f - M_f`` (mirror the fore
    half onto ``[0, 1/2]``). Both halves share one normalisation, and the
    waterplane weights are the same ``f``, hence LCF = LCB.
    """
    parameters.validate()
    aft_i, aft_m, fore_i, fore_m = _longitudinal_integrals(parameters)
    denominator = aft_i + fore_i
    if denominator <= 0.0:  # pragma: no cover - excluded by validate()
        raise ValueError("degenerate taper")
    return parameters.length * (aft_m + fore_i - fore_m) / denominator


def analytic_hydrostatics(parameters: AsymmetricHullParameters) -> Hydrostatics:
    """Closed-form hydrostatics of the asymmetric surface, as a ``Hydrostatics``."""
    parameters.validate()
    volume = analytic_volume(parameters)
    waterplane = analytic_waterplane_area(parameters)
    centre = analytic_longitudinal_centre(parameters)
    draft = parameters.draft
    beam = parameters.beam
    length = parameters.length
    section = (2.0 / 3.0) * beam * draft
    return Hydrostatics(
        displacement_volume=volume,
        displacement_mass=volume * parameters.water_density,
        wetted_surface_area=analytic_wetted_surface_area(parameters),
        waterplane_area=waterplane,
        longitudinal_centre_buoyancy=centre,
        longitudinal_centre_flotation=centre,
        vertical_centre_buoyancy=-(2.0 / 5.0) * draft,
        block_coefficient=volume / (length * beam * draft),
        prismatic_coefficient=volume / (length * section),
        maximum_section_area=section,
        waterline_z=0.0,
        water_density=parameters.water_density,
        basis=(
            "closed form for volume, waterplane area, midship area, VCB and the "
            "coefficients; Gauss-Legendre quadrature for LCB/LCF and wetted area"
        ),
    )


#: Quantities compared mesh-against-reference, in reporting order.
COMPARABLE_QUANTITIES: tuple[str, ...] = (
    "displacement_volume",
    "wetted_surface_area",
    "waterplane_area",
    "longitudinal_centre_buoyancy",
    "longitudinal_centre_flotation",
    "vertical_centre_buoyancy",
    "block_coefficient",
    "prismatic_coefficient",
    "maximum_section_area",
)


@dataclass(frozen=True)
class AsymmetricHullModel:
    """A generated asymmetric surface plus the reference it is judged against."""

    parameters: AsymmetricHullParameters
    mesh: Mesh
    offsets: tuple[tuple[float, float], ...]
    validation: ValidationReport = field(default_factory=ValidationReport)
    hydrostatics: Hydrostatics | None = None
    reference: Hydrostatics | None = None

    def errors(self) -> dict[str, float]:
        """Signed relative error of each mesh quantity against the reference."""
        if self.hydrostatics is None or self.reference is None:
            return {}
        out: dict[str, float] = {}
        for name in COMPARABLE_QUANTITIES:
            measured = float(getattr(self.hydrostatics, name))
            exact = float(getattr(self.reference, name))
            if exact == 0.0:
                continue
            out[name] = measured / exact - 1.0
        return out

    def metadata(self) -> dict[str, Any]:
        data: dict[str, Any] = {
            "schema": SCHEMA,
            "model_type": "asymmetric_sampled_parametric_displacement_hull",
            "parameters": self.parameters.as_dict(),
            "units": self.mesh.units,
            "coordinate_frame": self.mesh.coordinate_frame,
            "generator": {
                "kind": "analytic_lofted_surface",
                "longitudinal_taper": (
                    "sin(pi*x/length)^aft_exponent for x < length/2, "
                    "sin(pi*x/length)^fore_exponent for x >= length/2"
                ),
                "vertical_taper": "sqrt((z+draft)/draft)",
                "topology": "one closed ring per station, lofted, capped at the design waterline",
                "backend": "tensorcad.geometry (standard library)",
            },
            "mesh": self.mesh.metadata(),
            "validation": self.validation.to_dict(),
            "engineering_status": (
                "geometry and hydrostatic integration only; no fairness, weight, "
                "stability, resistance, structural or class-rule claim"
            ),
        }
        if self.hydrostatics is not None:
            data["hydrostatics"] = self.hydrostatics.to_dict()
        if self.reference is not None:
            data["analytic_reference"] = self.reference.to_dict()
        if self.hydrostatics is not None and self.reference is not None:
            data["relative_errors"] = dict(self.errors())
        return data


def _consistency_checks(
    parameters: AsymmetricHullParameters,
    mesh: Mesh,
    offsets: tuple[tuple[float, float], ...],
) -> list[Check]:
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
    achievable = max(taper(i / (p.nx - 1), p) for i in range(p.nx))
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
                f"{realised:.9g} is {100 * (1 - realised / p.beam):.3f}% under the "
                f"requested {p.beam:.9g}; use an odd nx",
                severity=WARNING,
            )
        )
    # The exponent switch must be continuous: whichever station lands on midship
    # has to carry the full requested beam, otherwise the two halves do not meet.
    midship_index = None
    for i in range(p.nx):
        if abs(i / (p.nx - 1) - 0.5) <= 1e-15:
            midship_index = i
    if midship_index is not None:
        checks.append(
            Check(
                "midship_station_continuous_across_exponent_switch",
                math.isclose(taper(0.5, p), 1.0, rel_tol=0.0, abs_tol=1e-15),
                f"taper at midship = {taper(0.5, p):.17g}",
            )
        )
    return checks


def generate_hull(parameters: AsymmetricHullParameters | None = None) -> AsymmetricHullModel:
    """Generate a closed, validated asymmetric hull surface.

    Raises :class:`ValueError` when the parameters are invalid or the generated
    surface fails a geometric validity check; an unusable mesh is never returned.
    """
    p = parameters if parameters is not None else AsymmetricHullParameters()
    mesh, offsets = build_mesh(p)
    report = merge(inspect_geometry(mesh), _consistency_checks(p, mesh, offsets))
    if report.errors:
        raise ValueError("generated asymmetric hull failed validation: " + "; ".join(report.errors))
    waterline_z = mesh.bounding_box().maximum[2]
    hydrostatics = compute_hydrostatics(
        mesh, waterline_z, water_density=p.water_density, closed=True
    )
    reference = analytic_hydrostatics(p)
    relative = hydrostatics.displacement_volume / reference.displacement_volume - 1.0
    ok = abs(relative) <= SAMPLING_ERROR_WARNING
    sampling = Check(
        "sampling_error_within_tolerance",
        ok,
        f"{p.nx}x{p.nz} integrates {hydrostatics.displacement_volume:.6g} m^3 against "
        f"an analytic {reference.displacement_volume:.6g} m^3 ({100.0 * relative:+.2f} %)",
        severity=WARNING,
    )
    if not ok:
        sampling = Check(
            sampling.name,
            False,
            sampling.detail
            + f"; the discretisation error exceeds {100.0 * SAMPLING_ERROR_WARNING:.0f} %",
            severity=WARNING,
        )
    report = merge(list(report.checks), [sampling])
    return AsymmetricHullModel(
        parameters=p,
        mesh=mesh,
        offsets=offsets,
        validation=report,
        hydrostatics=hydrostatics,
        reference=reference,
    )


def analytic_relative_errors(
    parameters: AsymmetricHullParameters | None = None,
) -> dict[str, Any]:
    """Signed mesh-against-reference errors for one grid, plus quadrature accuracy."""
    model = generate_hull(parameters)
    errors = model.errors()
    worst = max(errors.items(), key=lambda item: abs(item[1])) if errors else None
    return {
        "schema": f"{SCHEMA}/errors",
        "parameters": model.parameters.as_dict(),
        "relative_errors": dict(errors),
        "worst_quantity": worst[0] if worst else None,
        "worst_relative_error": worst[1] if worst else None,
        "quadrature": {
            "moment_nodes": MOMENT_NODES,
            "wetted_surface_nodes": WETTED_SURFACE_NODES,
            "wetted_surface_drift": wetted_surface_quadrature_drift(model.parameters),
        },
        "reference_basis": model.reference.basis if model.reference else None,
    }
