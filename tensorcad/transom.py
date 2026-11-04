"""Transom-stern hull family: the aft body truncated by a flat transom plate.

The parametric hull in :mod:`tensorcad.hull` tapers to a point at both ends
(``sin(pi*x/length)**0.7``), which makes it closed but geometrically useless aft:
no real displacement hull ends in a knife edge at the waterline. This module adds
the missing stern form -- an **immersed transom**, the flat, roughly trapezoidal
plate that most ships above Froude ~0.25 actually carry -- by truncating the body
at a station aft of midship and capping the opening with a planar face.

Two parameters turn the family into a design space rather than a single shape:

``transom_fraction`` ``tau``
    Longitudinal station of the transom, as a fraction of the length from the
    (untruncated) aft end. ``tau = 0`` reproduces the pointed hull of
    :mod:`tensorcad.hull` exactly -- same stations, same vertices, same topology
    -- so this module is a strict generalisation, not a fork.
``breadth_exponent`` ``p``
    Exponent of the sine taper. ``p < 1`` keeps the ends finely pointed,
    ``p > 1`` fattens the shoulders. It is the knob that moves block coefficient
    across the practical range.

Why the closed forms survive the truncation
-------------------------------------------
Every hydrostatic reference in :mod:`tensorcad.hull` is an integral of the
analytic surface. Truncating the body replaces a complete integral of
``sin(pi*u)**p`` with an *incomplete* one, which has no elementary closed form
for general ``p``. Rather than drop the analytic reference (the whole point of
which is to be independent of the mesh), the incomplete integrals are evaluated
by a substitution that removes the integrand's only singularity:

``K(a) = int_0^a sin(pi*s)**p ds``, and the integrand splits into a closed-form
singular factor and an analytic remainder: ``sin(pi*s)**p = (pi*s)**p *
(sinc(pi*s))**p``, where the second factor is even in ``s`` and therefore a power
series in ``s**2``. Term by term the integral of ``s**p`` is a power, so ``K``
is a rapidly convergent series (:func:`sinc_power_coefficients`), evaluated to
round-off in double precision for every exponent in the supported range.

The series is self-validated against an independent closed form: at ``tau = 0``
the complete integral must equal the gamma-function form in
``tensorcad.hydrostatics_analytic.taper_integral``,
and ``tests/test_transom.py`` asserts agreement at the 1e-14 level for p from
0.25 to 3. A reference that agrees with an independent closed form to round-off
is a reference a mesh can be measured against; one that only agrees with itself
is not.

All quantities are per the usual naval-architecture frame of :mod:`tensorcad.hull`:
``x`` from the untruncated aft end toward the stem, ``z = 0`` at the design
waterline, the keel at ``z = -draft``.
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass
from typing import Any

from . import mesh as meshlib
from .geometry import Mesh, MeshBuilder
from .hull import (
    BREADTH_EXPONENT,
    COORDINATE_FRAME,
    UNITS,
    analytic_vertical_centre_buoyancy,
)
from .hydrostatics import Hydrostatics, compute_hydrostatics
from .hydrostatics_analytic import taper_integral
from .validation import (
    WARNING,
    Check,
    DimensionTarget,
    ValidationReport,
    check_dimensions,
    inspect_geometry,
    merge,
)

#: Transom stations must stay strictly aft of midship: at and beyond ``0.5`` the
#: truncation would cut into the parallel middle body and change the maximum
#: section area, which the whole hull family is parameterised around.
MAX_TRANSOM_FRACTION = 0.5

#: Exponent range of the sine taper. Below ~0.3 the ends are needles that the
#: sampler cannot resolve; above ~2.0 the "pointed" ends are not pointed at all.
MIN_BREADTH_EXPONENT = 0.25
MAX_BREADTH_EXPONENT = 3.0

#: Terms kept in the ``sinc**p`` power series that gives the incomplete taper
#: integrals. Convergence is factorial, so 24 terms carry the aft half of the
#: integral (|w| <= 1/4) past double precision; the order is fixed, never
#: adaptive, so results are bit-reproducible across runs and platforms.
SERIES_ORDER = 24

#: Relative deviation from the analytic displacement at which the sampled grid is
#: reported as too coarse to trust (same policy as :mod:`tensorcad.hull`).
SAMPLING_ERROR_WARNING = 0.02

_SERIES_CACHE: dict[tuple[float, int], tuple[float, ...]] = {}


def sinc_power_coefficients(exponent: float, order: int = SERIES_ORDER) -> tuple[float, ...]:
    """Power-series coefficients of ``(sinc(pi*sqrt(w)))**p`` in ``w``.

    ``sin(pi*s)**p = (pi*s)**p * (sinc(pi*s))**p``, and ``sinc(pi*s)`` is even in
    ``s``, so the second factor is a power series in ``w = s**2`` starting at 1.
    That split is what makes the incomplete integrals exactly integrable: the
    ``s**p`` factor carries the end singularity in closed form, and the remaining
    series has a finite radius (its first zero is at ``s = 1``).

    The coefficients come from ``exp(p * log(sinc))`` by the two standard
    O(n**2) recurrences -- ``log`` from ``C' S = S'``, ``exp`` from
    ``E' = p A' E`` -- so no factorial or series library is needed.
    """
    p = float(exponent)
    if not (math.isfinite(p) and p > 0.0):
        raise ValueError(f"exponent must be a finite positive number, got {exponent!r}")
    key = (p, order)
    cached = _SERIES_CACHE.get(key)
    if cached is not None:
        return cached
    if order < 2:
        raise ValueError("order must be >= 2")
    # sinc(pi*sqrt(w)) = sum_k (-1)^k pi^(2k) w^k / (2k+1)!
    sinc = [(-1.0) ** k * math.pi ** (2 * k) / math.factorial(2 * k + 1) for k in range(order)]
    log_sinc = [0.0] * order
    for n in range(1, order):
        acc = 0.0
        for k in range(1, n):
            acc += k * log_sinc[k] * sinc[n - k]
        log_sinc[n] = sinc[n] - acc / n
    exp_p = [0.0] * order
    exp_p[0] = 1.0
    for n in range(1, order):
        acc = 0.0
        for k in range(1, n + 1):
            acc += k * log_sinc[k] * exp_p[n - k]
        exp_p[n] = p * acc / n
    result = tuple(exp_p)
    _SERIES_CACHE[key] = result
    return result


def _rising_taper_moment(limit: float, exponent: float, moment: int = 0) -> float:
    """``int_0^limit s**moment * sin(pi*s)**exponent ds`` for ``0 <= limit <= 1/2``.

    Term by term the integrand is a power of ``s``, so the integral is the
    series ``pi**p * sum_k c_k * limit**(moment+p+2k+1) / (moment+p+2k+1)`` with
    ``c_k`` from :func:`sinc_power_coefficients`. Staying inside ``[0, 1/2]`` keeps
    ``w = s**2 <= 1/4`` far inside the series radius; the forward half of the
    complete integral is recovered by symmetry in
    :func:`truncated_taper_integral`.
    """
    if not (0.0 <= limit <= 0.5):
        raise ValueError(f"limit must lie in [0, 0.5], got {limit!r}")
    if limit == 0.0:
        return 0.0
    p = float(exponent)
    total = 0.0
    for k, coefficient in enumerate(sinc_power_coefficients(p)):
        power = moment + p + 2 * k + 1
        total += coefficient * limit**power / power
    return math.pi**p * total


def truncated_taper_integral(tau: float, exponent: float = BREADTH_EXPONENT) -> float:
    """``int_tau^1 sin(pi*u)**p du``: the aperture left after the aft cut.

    Split at midship and mirrored, ``J(tau) = 2*K(1/2) - K(tau)`` with ``K`` the
    substituted rising integral, so only the aft singularity ever needs handling.
    ``tau = 0`` returns :func:`tensorcad.hydrostatics_analytic.taper_integral` -- the
    independent gamma-function closed form -- to round-off.
    """
    _check_fraction("tau", tau)
    top = 2.0 * _rising_taper_moment(0.5, exponent, 0) - _rising_taper_moment(tau, exponent, 0)
    if top < 0.0:  # pragma: no cover - defensive: the integral is positive by construction
        raise AssertionError(f"negative taper integral {top!r}")
    return top


def truncated_taper_moment(tau: float, exponent: float = BREADTH_EXPONENT) -> float:
    """``int_tau^1 u * sin(pi*u)**p du``: the first moment that locates LCB."""
    _check_fraction("tau", tau)
    return _rising_taper_moment(0.5, exponent, 0) - _rising_taper_moment(tau, exponent, 1)


def _check_fraction(name: str, value: float) -> None:
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite, got {value!r}")
    if not (0.0 <= value < MAX_TRANSOM_FRACTION):
        raise ValueError(f"{name} must lie in [0, {MAX_TRANSOM_FRACTION:g}), got {value!r}")


def transom_area_ratio(tau: float, exponent: float = BREADTH_EXPONENT) -> float:
    """Immersed transom area as a fraction of the midship section area.

    Every section is ``(2/3) * beam * draft * f(x)`` with ``f`` the sine taper, so
    the ratio is simply ``sin(pi*tau)**p``: the transom plate is the section at
    ``x = 0`` and midship carries the unit taper. For the default ``p = 0.7`` a
    realistic few-percent transom therefore needs ``tau`` of order 0.01-0.02 -- the
    sine run is so fine aft that any larger cut opens an unphysically large plate.
    Fuller bodies (larger ``p``) allow a visibly larger transom at the same
    immersed-area ratio, which is the design trade this pair of parameters exposes.
    """
    _check_fraction("tau", tau)
    return math.pow(math.sin(math.pi * tau), float(exponent))


def transom_fraction_for_area_ratio(
    area_ratio: float, exponent: float = BREADTH_EXPONENT
) -> float:
    """Inverse of :func:`transom_area_ratio`: the tau that gives a stated ratio.

    Closed form, ``tau = asin(ratio**(1/p)) / pi``, monotone on ``[0, 1)``, so a
    design loop can dial an immersed transom area directly instead of bisecting.
    """
    if not (math.isfinite(area_ratio) and 0.0 < area_ratio < 1.0):
        raise ValueError(f"area_ratio must lie in (0, 1), got {area_ratio!r}")
    tau = math.asin(area_ratio ** (1.0 / float(exponent))) / math.pi
    if not (0.0 < tau < MAX_TRANSOM_FRACTION):
        raise ValueError(
            f"area_ratio {area_ratio!r} at exponent {exponent!r} implies tau = {tau!r}, "
            f"outside (0, {MAX_TRANSOM_FRACTION:g})"
        )
    return tau


@dataclass(frozen=True)
class TransomHullParameters:
    """Principal dimensions plus the two parameters that define the stern form.

    ``length`` is the overall length from the transom to the stem (LOA).
    ``transom_fraction = 0`` is the pointed aft body of :mod:`tensorcad.hull`; any
    positive value below 0.5 truncates the body and adds the transom plate, and
    the parent body the taper is defined on is then ``LOA / (1 - tau)`` long.
    """

    length: float = 10.0
    beam: float = 2.0
    draft: float = 0.8
    nx: int = 41
    nz: int = 13
    water_density: float = 1025.0
    transom_fraction: float = 0.2
    breadth_exponent: float = BREADTH_EXPONENT

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)

    def problems(self) -> list[str]:
        """Every problem with this parameter set, not just the first."""
        errors: list[str] = []
        for name in ("length", "beam", "draft"):
            value = getattr(self, name)
            if not math.isfinite(value) or value <= 0.0:
                errors.append(f"{name} must be a finite positive number, got {value!r}")
        if not (isinstance(self.nx, int) and self.nx >= 3):
            errors.append(f"nx must be an integer >= 3, got {self.nx!r}")
        if not (isinstance(self.nz, int) and self.nz >= 3):
            errors.append(f"nz must be an integer >= 3, got {self.nz!r}")
        if not (math.isfinite(self.water_density) and self.water_density > 0.0):
            errors.append(f"water_density must be a finite positive number, got {self.water_density!r}")
        if not (math.isfinite(self.transom_fraction) and 0.0 <= self.transom_fraction < MAX_TRANSOM_FRACTION):
            errors.append(
                f"transom_fraction must lie in [0, {MAX_TRANSOM_FRACTION:g}) -- at or beyond "
                f"midship the cut would remove the maximum section; got {self.transom_fraction!r}"
            )
        if not (
            math.isfinite(self.breadth_exponent)
            and MIN_BREADTH_EXPONENT <= self.breadth_exponent <= MAX_BREADTH_EXPONENT
        ):
            errors.append(
                f"breadth_exponent must lie in [{MIN_BREADTH_EXPONENT:g}, {MAX_BREADTH_EXPONENT:g}], "
                f"got {self.breadth_exponent!r}"
            )
        return errors

    def validate(self) -> None:
        errors = self.problems()
        if errors:
            raise ValueError("; ".join(errors))


@dataclass(frozen=True)
class TransomModel:
    """A generated transom-stern hull plus the references it was checked against."""

    parameters: TransomHullParameters
    mesh: Mesh
    offsets: tuple[tuple[float, float], ...]
    hydrostatics: Hydrostatics
    analytic: dict[str, float]
    validation: ValidationReport
    transom_station_x: float
    transom_immersed_area_analytic: float
    transom_breadth_at_waterline: float

    @property
    def units(self) -> str:
        return self.mesh.units

    @property
    def coordinate_frame(self) -> str:
        return self.mesh.coordinate_frame

    def relative_errors(self) -> dict[str, float]:
        """Mesh integration against the analytic references, as relative errors."""
        exact = self.analytic
        measured = self.hydrostatics
        pairs = {
            "displacement_volume": (measured.displacement_volume, exact["displacement_volume"]),
            "waterplane_area": (measured.waterplane_area, exact["waterplane_area"]),
            "longitudinal_centre_buoyancy": (
                measured.longitudinal_centre_buoyancy,
                exact["longitudinal_centre_buoyancy"],
            ),
            "vertical_centre_buoyancy": (
                measured.vertical_centre_buoyancy,
                exact["vertical_centre_buoyancy"],
            ),
            "block_coefficient": (measured.block_coefficient, exact["block_coefficient"]),
        }
        out: dict[str, float] = {}
        for name, (got, want) in pairs.items():
            out[name] = (got - want) / want if abs(want) > 1e-15 else float("nan")
        return out

    def metadata(self) -> dict[str, Any]:
        data: dict[str, Any] = {
            "schema_version": "1.0",
            "model_type": "closed_sampled_transom_stern_hull",
            "parameters": self.parameters.as_dict(),
            "units": self.units,
            "coordinate_frame": self.coordinate_frame,
            "generator": {
                "kind": "analytic_lofted_surface_with_planar_transom_cap",
                "longitudinal_taper": f"sin(pi*x/length)^{self.parameters.breadth_exponent:g}",
                "vertical_taper": "sqrt((z+draft)/draft)",
                "truncation": (
                    f"aft body cut at x = {self.transom_station_x:.9g} m "
                    f"(transom_fraction = {self.parameters.transom_fraction:.6g}), "
                    "opening closed by a planar fan"
                ),
                "topology": "one closed section ring per station, lofted, capped at the design waterline",
                "backend": "tensorcad.geometry (standard library)",
            },
            "transom": {
                "station_x_m": self.transom_station_x,
                "breadth_at_waterline_m": self.transom_breadth_at_waterline,
                "immersed_area_analytic_m2": self.transom_immersed_area_analytic,
            },
            "analytic_reference": dict(self.analytic),
            "mesh": self.mesh.metadata(),
            "validation": self.validation.to_dict(),
            "engineering_status": (
                "geometry and hydrostatic integration only; no resistance, stability, "
                "structure or classification-rule checks"
            ),
        }
        data["hydrostatics"] = self.hydrostatics.to_dict()
        data["relative_error"] = self.relative_errors()
        return data


def parent_length(parameters: TransomHullParameters) -> float:
    """Length of the untruncated pointed hull this transom body was cut from.

    ``length`` in :class:`TransomHullParameters` is the overall length from the
    transom to the stem (LOA), which is what a naval architect measures. The sine
    taper, however, is defined on the *untruncated* body, so the reference length
    is ``LOA / (1 - tau)``.
    """
    p = parameters
    return p.length / (1.0 - p.transom_fraction)


def _station_x(parameters: TransomHullParameters, index: int) -> float:
    """Axial position of station ``index`` in the mesh frame; the transom is at x=0."""
    return parameters.length * index / (parameters.nx - 1)


def _taper_at(parameters: TransomHullParameters, x: float) -> float:
    """Sine taper at mesh position ``x``, forced to exactly zero at the stem.

    ``math.sin(math.pi)`` is not zero in binary floating point, so the stem station
    is pinned: a residual breadth there would leave a slot at the stem.
    """
    p = parameters
    u = p.transom_fraction + x / parent_length(p)
    if u >= 1.0:
        return 0.0
    if u <= 0.0:
        return 0.0
    return math.pow(math.sin(math.pi * u), p.breadth_exponent)


def _ring_size(parameters: TransomHullParameters) -> int:
    return 2 * parameters.nz - 1


def _ring_point(parameters: TransomHullParameters, x: float, taper: float, j: int) -> tuple[float, float, float]:
    """Ring position ``j`` of the section at ``x``.

    ``j`` walks the outline once: down the starboard flank from the design
    waterline to the keel, then up the port flank. The keel point belongs to the
    ring once, which is what makes the lofted surface closed without welding.
    """
    p = parameters
    keelward = j <= p.nz - 1
    level = (p.nz - 1 - j) if keelward else (j - (p.nz - 1))
    z = -p.draft * (1.0 - level / (p.nz - 1))
    half_breadth = 0.5 * p.beam * taper * math.sqrt(level / (p.nz - 1))
    y = half_breadth if keelward else -half_breadth
    return (x, y, z)


def build_transom_mesh(
    parameters: TransomHullParameters,
) -> tuple[Mesh, tuple[tuple[float, float], ...]]:
    """Loft the truncated body, cap the waterline, close the stern with a plate.

    Deterministic: identical parameters give identical vertex and face tuples.
    With ``transom_fraction = 0`` the transom plate is skipped (the aft ring is
    already degenerate on the centreline) and the result is the pointed hull.
    """
    parameters.validate()
    p = parameters
    ring_size = _ring_size(p)
    builder = MeshBuilder(units=UNITS, coordinate_frame=COORDINATE_FRAME)

    xs = [_station_x(p, i) for i in range(p.nx)]
    tapers = [_taper_at(p, x) for x in xs]
    rings: list[list[int]] = []
    for x, taper in zip(xs, tapers):
        rings.append([builder.add_vertex(_ring_point(p, x, taper, j)) for j in range(ring_size)])

    for i in range(p.nx - 1):
        for j in range(ring_size - 1):
            builder.add_quad(rings[i][j], rings[i + 1][j], rings[i + 1][j + 1], rings[i][j + 1])
        builder.add_quad(
            rings[i][ring_size - 1],
            rings[i + 1][ring_size - 1],
            rings[i + 1][0],
            rings[i][0],
        )

    if p.transom_fraction > 0.0:
        # The section is convex (a parabolic keel profile under the waterline), so a
        # fan from the starboard waterline corner covers it with non-degenerate
        # triangles and introduces no new vertices.
        for j in range(1, ring_size - 1):
            builder.add_triangle(rings[0][0], rings[0][j], rings[0][j + 1])

    mesh = meshlib.orient_outward(builder.build()).mesh
    offsets = tuple((xs[i], mesh.vertices[rings[i][0]][1] - mesh.vertices[rings[i][ring_size - 1]][1])
                    for i in range(p.nx))
    return mesh, offsets


def analytic_transom_hydrostatics(parameters: TransomHullParameters) -> dict[str, float]:
    """Closed-form (to round-off) hydrostatics of the analytic transom surface.

    With ``J = int_tau^1 sin(pi*u)**p du`` and its first moment ``M``::

        volume            = (2/3) * B * draft * L_parent * J
        waterplane area   = B * L_parent * J
        LCB / LCF (fwd of transom) = L_parent * (M/J - tau)
        VCB               = -2/5 * draft                 (sqrt profile, tau-independent)
        block coefficient = (2/3) * J / (1 - tau)         (over the LOA)
        waterplane coeff. = J / (1 - tau)

    ``L_parent = LOA / (1 - tau)`` is the untruncated reference length, and the
    coefficients are taken over the LOA so they agree with what
    :func:`tensorcad.hydrostatics.compute_hydrostatics` reads off the mesh.
    """
    p = parameters
    p.validate()
    exponent = p.breadth_exponent
    tau = p.transom_fraction
    j_integral = truncated_taper_integral(tau, exponent)
    m_integral = truncated_taper_moment(tau, exponent)
    if not (math.isfinite(j_integral) and j_integral > 0.0):
        raise AssertionError(f"degenerate taper integral {j_integral!r}")
    reference_length = parent_length(p)
    lcb_from_parent_aft = reference_length * m_integral / j_integral
    section_scale = (2.0 / 3.0) * p.beam * p.draft
    volume = section_scale * reference_length * j_integral
    waterplane = p.beam * reference_length * j_integral
    transom_taper = math.pow(math.sin(math.pi * tau), exponent) if tau > 0.0 else 0.0
    return {
        "parent_length": reference_length,
        "displacement_volume": volume,
        "displacement_mass": volume * p.water_density,
        "waterplane_area": waterplane,
        "longitudinal_centre_buoyancy": lcb_from_parent_aft - tau * reference_length,
        "longitudinal_centre_flotation": lcb_from_parent_aft - tau * reference_length,
        "vertical_centre_buoyancy": analytic_vertical_centre_buoyancy(p.draft),
        "block_coefficient": (2.0 / 3.0) * j_integral / (1.0 - tau),
        "waterplane_coefficient": j_integral / (1.0 - tau),
        "maximum_section_area": section_scale,
        "transom_station_x": 0.0,
        "transom_breadth_at_waterline": p.beam * transom_taper,
        "transom_immersed_area": section_scale * transom_taper,
    }


def _dimension_target(parameters: TransomHullParameters) -> DimensionTarget:
    p = parameters
    achievable = max(_taper_at(p, _station_x(p, i)) for i in range(p.nx))
    return DimensionTarget(
        length=p.length,
        beam=p.beam * achievable,
        draft=p.draft,
        waterline_z=0.0,
    )


def _transom_checks(parameters: TransomHullParameters, mesh: Mesh, analytic: dict[str, float]) -> list[Check]:
    """Checks that only exist because the stern is truncated."""
    p = parameters
    checks: list[Check] = []
    if p.transom_fraction <= 0.0:
        return [
            Check(
                "transom_absent",
                True,
                "transom_fraction = 0: pointed aft body, no transom plate, no immersed transom area",
            )
        ]
    measured_area = 0.0
    station = analytic["transom_station_x"]
    for face in mesh.faces:
        a, b, c = (mesh.vertices[index] for index in face)
        if all(abs(vertex[0] - station) <= 1e-12 for vertex in (a, b, c)):
            measured_area += 0.5 * abs(
                (b[1] - a[1]) * (c[2] - a[2]) - (c[1] - a[1]) * (b[2] - a[2])
            )
    expected = analytic["transom_immersed_area"]
    relative = (measured_area - expected) / expected if expected > 0.0 else float("nan")
    checks.append(
        Check(
            "transom_plate_area",
            math.isfinite(relative) and abs(relative) <= 0.05,
            f"transom plate {measured_area:.6g} m^2 vs analytic section {expected:.6g} m^2 "
            f"({100.0 * relative:+.2f} %; a piecewise-linear section under-reads the parabolic profile)",
            severity=WARNING,
        )
    )
    checks.append(
        Check(
            "transom_aft_of_midship",
            p.transom_fraction < MAX_TRANSOM_FRACTION,
            f"transom station x = {station:.6g} m is aft of midship x = {0.5 * p.length:.6g} m",
        )
    )
    return checks


def generate_transom_hull(parameters: TransomHullParameters | None = None) -> TransomModel:
    """Generate and validate a closed transom-stern hull.

    Raises :class:`ValueError` for invalid parameters or when the surface fails a
    geometric validity check -- an unwatertight or non-manifold stern is never
    handed back silently.
    """
    p = parameters if parameters is not None else TransomHullParameters()
    mesh, offsets = build_transom_mesh(p)
    analytic = analytic_transom_hydrostatics(p)
    waterline_z = mesh.bounding_box().maximum[2]
    hydrostatics = compute_hydrostatics(mesh, waterline_z, water_density=p.water_density, closed=True)
    report = merge(
        inspect_geometry(mesh),
        check_dimensions(mesh, _dimension_target(p)),
        _transom_checks(p, mesh, analytic),
        [
            Check(
                "analytic_displacement_match",
                abs(hydrostatics.displacement_volume - analytic["displacement_volume"])
                <= SAMPLING_ERROR_WARNING * analytic["displacement_volume"],
                f"{p.nx}x{p.nz} grid integrates {hydrostatics.displacement_volume:.6g} m^3 against an "
                f"analytic {analytic['displacement_volume']:.6g} m^3 "
                f"({100.0 * (hydrostatics.displacement_volume / analytic['displacement_volume'] - 1.0):+.2f} %)",
                severity=WARNING,
            )
        ],
    )
    errors = report.errors
    if errors:
        raise ValueError("generated transom hull failed validation: " + "; ".join(errors))
    return TransomModel(
        parameters=p,
        mesh=mesh,
        offsets=offsets,
        hydrostatics=hydrostatics,
        analytic=analytic,
        validation=report,
        transom_station_x=analytic["transom_station_x"],
        transom_immersed_area_analytic=analytic["transom_immersed_area"],
        transom_breadth_at_waterline=analytic["transom_breadth_at_waterline"],
    )


def estimate_transom_discretisation_error(
    parameters: TransomHullParameters, *, growth: float = 2.0
) -> dict[str, Any]:
    """Refine the grid once and report the relative change of the integrals.

    Same caveat as the pointed-hull version: this measures *self*consistency, and
    the analytic references in :func:`analytic_transom_hydrostatics` are the independent
    check. The two are reported side by side so a reader can see both.
    """
    if not (math.isfinite(growth) and growth > 1.0):
        raise ValueError("growth must be finite and > 1.0")
    p = parameters
    refined = TransomHullParameters(
        length=p.length,
        beam=p.beam,
        draft=p.draft,
        nx=min(4001, math.ceil((p.nx - 1) * growth) + 1),
        nz=min(513, math.ceil((p.nz - 1) * growth) + 1),
        water_density=p.water_density,
        transom_fraction=p.transom_fraction,
        breadth_exponent=p.breadth_exponent,
    )
    coarse_hull = generate_transom_hull(p).hydrostatics
    fine_hull = generate_transom_hull(refined).hydrostatics
    exact = analytic_transom_hydrostatics(p)
    out: dict[str, Any] = {
        "status": "measured",
        "coarse_grid": {"nx": p.nx, "nz": p.nz},
        "refined_grid": {"nx": refined.nx, "nz": refined.nz},
        "analytic": {
            "displacement_volume": exact["displacement_volume"],
            "waterplane_area": exact["waterplane_area"],
            "longitudinal_centre_buoyancy": exact["longitudinal_centre_buoyancy"],
        },
    }
    for name in ("displacement_volume", "waterplane_area", "longitudinal_centre_buoyancy"):
        coarse = float(getattr(coarse_hull, name))
        fine = float(getattr(fine_hull, name))
        out[name] = {
            "coarse": coarse,
            "refined": fine,
            "refinement_change": (coarse - fine) / fine if abs(fine) > 1e-15 else float("nan"),
            "coarse_error_vs_analytic": (coarse - exact[name]) / exact[name],
            "refined_error_vs_analytic": (fine - exact[name]) / exact[name],
        }
    return out