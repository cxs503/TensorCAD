"""Port/starboard asymmetric transom-stern hull family.

:mod:`tensorcad.transom` gives the aft body a flat transom plate but keeps the two
sides mirror images of each other; :mod:`tensorcad.hull`'s ``fore_exponent`` /
``aft_exponent`` asymmetry is *longitudinal*. The half of the backlog item
"asymmetric hull + transom" that stayed open (see ``docs/transom_stern.md`` §6) is
the transverse one: a hull whose **starboard side is wider than its port side**.
This module closes it, and it closes it the way every other family here is closed
-- with closed-form hydrostatics that a mesh can be measured against, not with a
mesh and a story.

How the asymmetry is built
--------------------------
The section keeps the separable profile of the transom family, with one half
breadth per side::

    y_s(x, z) = +(beam/2) * (1 + skew) * f(x) * sqrt(1 + z/draft)
    y_p(x, z) = -(beam/2) * (1 - skew) * f(x) * sqrt(1 + z/draft)
    f(x)      = sin(pi * (tau + x/L_parent))**p,     L_parent = LOA/(1 - tau)

``beam`` is therefore the **overall waterline breadth** (starboard half plus port
half), not one side's: ``y_s(0) + |y_p(0)| = beam`` at every station, so the
bounding box is unchanged by the asymmetry. ``skew = 0`` reproduces
:mod:`tensorcad.transom` bit for bit -- same stations, same vertices, same
topology -- which is the strict-generalisation test the other families use.

What the asymmetry does, and does not, move
-------------------------------------------
Because the two half breadths add to the same ``beam * f(x)``, everything that is
an integral of the *sum* is skew-independent, to round-off:

* displacement volume, waterplane area, LCB, LCF, VCB, ``Cb``, ``Cp``, ``Cwp``;
* the immersed transom area and the maximum section area.

What does move is the **transverse** first moment, which integrates the
*difference*:

    y_B = (1/V) int y dV = (3/8) * skew * beam * J2/J1
    y_F = (1/A_w) int y dA = (1/2) * skew * beam * J2/J1

with ``J1 = int_tau^1 f du`` and ``J2 = int_tau^1 f**2 du`` (the exponents double
inside the second moment, which is why ``J2 = truncated_taper_integral(tau, 2p)``
and nothing new is needed for it). Their ratio ``y_F/y_B = 4/3`` is exact and
independent of ``tau``, ``p`` and the grid -- a shape invariant, not a fitted
number. The transverse waterplane inertia is

    I_T (about the centreline) = (B_s**3 + B_p**3)/24 * L_parent * J3,
    J3 = truncated_taper_integral(tau, 3p)

and the inertia a stability calculation wants is about the *floating* centre,
``I_T - A_w * y_F**2``; both are reported because confusing them is the standard
way an asymmetric waterplane gets a wrong metacentric height.

Scope (what this module deliberately does not do)
-------------------------------------------------
* **Upright equilibrium only.** The body is generated and integrated at one
  attitude. It says nothing about the heel at which this hull would actually
  float (the asymmetric waterplane has a heel it *prefers*), about `GZ`, or about
  resistance.
* The asymmetry is a **breadth asymmetry with a common taper**. Port and starboard
  do not get their own taper exponent, their own deadrise or their own flare.
* ``tau > 0`` is required: with no transom plate there is nothing left of this
  family that :mod:`tensorcad.hull` does not already provide, and the degenerate
  aft ring would only add zero-area faces to check.

The frame is the package's: ``x`` from the transom toward the stem, ``z = 0`` at
the design waterline, the keel at ``z = -draft``, metres, kg, s.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

from . import mesh as meshlib
from .geometry import Mesh, MeshBuilder
from .hull import BREADTH_EXPONENT, COORDINATE_FRAME, UNITS, analytic_vertical_centre_buoyancy
from .hydrostatics import Hydrostatics, compute_hydrostatics
from .transom import (
    MAX_TRANSOM_FRACTION,
    MAX_BREADTH_EXPONENT,
    MIN_BREADTH_EXPONENT,
    truncated_taper_integral,
    truncated_taper_moment,
)
from .validation import (
    WARNING,
    Check,
    DimensionTarget,
    ValidationReport,
    check_dimensions,
    inspect_geometry,
    merge,
)

__all__ = [
    "MAX_SKEW",
    "AsymmetricTransomModel",
    "AsymmetricTransomParameters",
    "analytic_asymmetric_hydrostatics",
    "build_asymmetric_mesh",
    "generate_asymmetric_transom_hull",
    "mesh_transverse_centroid",
    "simpson_reference",
]

#: Largest breadth skew the family accepts. ``|skew| = 1`` would collapse one side
#: to zero breadth (a wall-sided half hull); 0.75 keeps the narrow side at a
#: quarter of the wide one, which is already more asymmetry than any real hull.
MAX_SKEW = 0.75

#: Relative deviation from the analytic displacement at which the sampled grid is
#: reported as too coarse to trust (same policy as :mod:`tensorcad.transom`).
SAMPLING_ERROR_WARNING = 0.02


def _check_skew(value: float) -> None:
    if not (math.isfinite(value) and abs(value) <= MAX_SKEW):
        raise ValueError(f"skew must be finite and satisfy |skew| <= {MAX_SKEW:g}, got {value!r}")


@dataclass(frozen=True)
class AsymmetricTransomParameters:
    """A transom-stern body whose port and starboard breadths differ."""

    length: float = 10.0
    beam: float = 2.0
    draft: float = 0.8
    nx: int = 41
    nz: int = 13
    water_density: float = 1025.0
    transom_fraction: float = 0.02
    breadth_exponent: float = BREADTH_EXPONENT
    #: ``B_starboard = beam*(1+skew)``, ``B_port = beam*(1-skew)``; the sum is
    #: ``2*beam`` and the overall waterline breadth is ``beam`` at every station.
    skew: float = 0.0

    def problems(self) -> list[str]:
        errors: list[str] = []
        if not (math.isfinite(self.length) and self.length > 0.0):
            errors.append(f"length must be a finite positive number, got {self.length!r}")
        if not (math.isfinite(self.beam) and self.beam > 0.0):
            errors.append(f"beam must be a finite positive number, got {self.beam!r}")
        if not (math.isfinite(self.draft) and self.draft > 0.0):
            errors.append(f"draft must be a finite positive number, got {self.draft!r}")
        if not (isinstance(self.nx, int) and self.nx >= 3):
            errors.append(f"nx must be an integer >= 3, got {self.nx!r}")
        if not (isinstance(self.nz, int) and self.nz >= 3):
            errors.append(f"nz must be an integer >= 3, got {self.nz!r}")
        if not (math.isfinite(self.water_density) and self.water_density > 0.0):
            errors.append(f"water_density must be a finite positive number, got {self.water_density!r}")
        if not (
            math.isfinite(self.transom_fraction)
            and 0.0 < self.transom_fraction < MAX_TRANSOM_FRACTION
        ):
            errors.append(
                f"transom_fraction must lie in (0, {MAX_TRANSOM_FRACTION:g}) -- this family is "
                f"defined by an actual transom plate; use tensorcad.hull for a pointed stern; "
                f"got {self.transom_fraction!r}"
            )
        if not (
            math.isfinite(self.breadth_exponent)
            and MIN_BREADTH_EXPONENT <= self.breadth_exponent <= MAX_BREADTH_EXPONENT
        ):
            errors.append(
                f"breadth_exponent must lie in [{MIN_BREADTH_EXPONENT:g}, "
                f"{MAX_BREADTH_EXPONENT:g}], got {self.breadth_exponent!r}"
            )
        if not (math.isfinite(self.skew) and abs(self.skew) <= MAX_SKEW):
            errors.append(f"skew must be finite and satisfy |skew| <= {MAX_SKEW:g}, got {self.skew!r}")
        return errors

    def validate(self) -> None:
        errors = self.problems()
        if errors:
            raise ValueError("; ".join(errors))

    def half_breadths(self) -> tuple[float, float]:
        """``(starboard, port)`` half breadths at the design waterline, midship."""
        return 0.5 * self.beam * (1.0 + self.skew), 0.5 * self.beam * (1.0 - self.skew)


@dataclass(frozen=True)
class AsymmetricTransomModel:
    """A generated asymmetric transom hull plus the references it was checked against."""

    parameters: AsymmetricTransomParameters
    mesh: Mesh
    hydrostatics: Hydrostatics
    analytic: dict[str, float]
    validation: ValidationReport
    #: Mesh-measured transverse centre of buoyancy (signed tetrahedra), metres.
    transverse_centre_buoyancy: float

    @property
    def units(self) -> str:
        return self.mesh.units

    @property
    def coordinate_frame(self) -> str:
        return self.mesh.coordinate_frame

    def relative_errors(self) -> dict[str, float]:
        """Mesh minus analytic, relative, for the three quantities that matter."""
        a = self.analytic
        return {
            "displacement_volume": self.hydrostatics.displacement_volume
            / a["displacement_volume"]
            - 1.0,
            "waterplane_area": self.hydrostatics.waterplane_area / a["waterplane_area"] - 1.0,
            "longitudinal_centre_buoyancy": (
                self.hydrostatics.longitudinal_centre_buoyancy
                - a["longitudinal_centre_buoyancy"]
            )
            / a["displacement_volume"],
            "transverse_centre_buoyancy": (
                self.transverse_centre_buoyancy - a["transverse_centre_buoyancy"]
            )
            / a["beam"],
        }

    def to_dict(self) -> dict[str, Any]:
        return {
            "parameters": {
                "length_m": self.parameters.length,
                "beam_m": self.parameters.beam,
                "draft_m": self.parameters.draft,
                "nx": self.parameters.nx,
                "nz": self.parameters.nz,
                "transom_fraction": self.parameters.transom_fraction,
                "breadth_exponent": self.parameters.breadth_exponent,
                "skew": self.parameters.skew,
                "starboard_half_breadth_m": self.parameters.half_breadths()[0],
                "port_half_breadth_m": self.parameters.half_breadths()[1],
            },
            "hydrostatics": self.hydrostatics.to_dict(),
            "analytic": {key: _round(value) for key, value in self.analytic.items()},
            "transverse_centre_buoyancy_m": _round(self.transverse_centre_buoyancy),
            "relative_errors": {key: _round(value) for key, value in self.relative_errors().items()},
            "validation": [check.to_dict() for check in self.validation.checks],
            "validation_passed": self.validation.passed,
        }


def _round(value: float, digits: int = 12) -> float:
    return float(f"{value:.{digits}g}")


def parent_length(parameters: AsymmetricTransomParameters) -> float:
    """Length of the untruncated pointed body this transom was cut from."""
    p = parameters
    return p.length / (1.0 - p.transom_fraction)


def _taper_at(parameters: AsymmetricTransomParameters, x: float) -> float:
    p = parameters
    u = p.transom_fraction + x / parent_length(p)
    if u >= 1.0:
        return 0.0
    if u <= 0.0:
        return 0.0
    return math.pow(math.sin(math.pi * u), p.breadth_exponent)


def _ring_size(parameters: AsymmetricTransomParameters) -> int:
    return 2 * parameters.nz - 1


def _ring_point(
    parameters: AsymmetricTransomParameters, x: float, taper: float, j: int
) -> tuple[float, float, float]:
    """Ring point ``j``: starboard flank down to the keel, then the port flank up.

    The keel point belongs to the ring once (``j = nz - 1``), which is what makes
    the lofted surface closed without welding two halves.
    """
    p = parameters
    starboard_half, port_half = p.half_breadths()
    keelward = j <= p.nz - 1
    level = (p.nz - 1 - j) if keelward else (j - (p.nz - 1))
    z = -p.draft * (1.0 - level / (p.nz - 1))
    # sqrt(level/(nz-1)) == sqrt(1 + z/draft): the parabolic keel profile.
    profile = math.sqrt(level / (p.nz - 1))
    y = (starboard_half if keelward else -port_half) * taper * profile
    return (x, y, z)


def build_asymmetric_mesh(
    parameters: AsymmetricTransomParameters,
) -> tuple[Mesh, tuple[tuple[float, float], ...]]:
    """Loft the body, cap the waterline, close the stern with the transom plate.

    ``offsets`` pairs each station's ``x`` with its **overall** waterline breadth
    ``y_starboard - y_port``, so an offsets table for an asymmetric hull is a
    single number per station (the two half breadths are recoverable from
    :meth:`AsymmetricTransomParameters.half_breadths`). Deterministic: identical
    parameters give identical vertex and face tuples.
    """
    parameters.validate()
    p = parameters
    ring_size = _ring_size(p)
    builder = MeshBuilder(units=UNITS, coordinate_frame=COORDINATE_FRAME)
    xs = [p.length * i / (p.nx - 1) for i in range(p.nx)]
    rings: list[list[int]] = []
    for x in xs:
        taper = _taper_at(p, x)
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
    # Transom plate: the section is convex (a parabolic keel profile under the
    # waterline), so a fan from the starboard waterline corner covers it with
    # non-degenerate triangles. The fan vertex is on the outline at any skew.
    for j in range(1, ring_size - 1):
        builder.add_triangle(rings[0][0], rings[0][j], rings[0][j + 1])

    mesh = meshlib.orient_outward(builder.build()).mesh
    offsets = tuple(
        (xs[i], mesh.vertices[rings[i][0]][1] - mesh.vertices[rings[i][ring_size - 1]][1])
        for i in range(p.nx)
    )
    return mesh, offsets


def analytic_asymmetric_hydrostatics(parameters: AsymmetricTransomParameters) -> dict[str, float]:
    """Closed-form (to round-off) hydrostatics of the analytic asymmetric surface.

    With ``J1 = int_tau^1 f du``, ``M = int_tau^1 u f du``, ``J2 = int_tau^1 f**2 du``
    and ``J3 = int_tau^1 f**3 du`` (``f = sin(pi u)**p``, so ``Jk`` is the same
    incomplete integral at exponent ``k*p``)::

        volume            = (2/3) * beam * draft * L_parent * J1
        waterplane area   = beam * L_parent * J1
        LCB / LCF         = L_parent * (M/J1 - tau)      (forward of the transom)
        VCB               = -2/5 * draft
        y_B               = (3/8) * skew * beam * J2/J1   (toward starboard for +skew)
        y_F               = (1/2) * skew * beam * J2/J1
        I_T (centreline)  = (B_s**3 + B_p**3)/24 * L_parent * J3
        I_T (centroid)    = I_T - waterplane_area * y_F**2

    Every quantity that is an integral of the *sum* of the half breadths carries
    no ``skew`` at all; that is the statement the tests make about it.
    """
    p = parameters
    p.validate()
    exponent = p.breadth_exponent
    tau = p.transom_fraction
    j1 = truncated_taper_integral(tau, exponent)
    m1 = truncated_taper_moment(tau, exponent)
    j2 = truncated_taper_integral(tau, 2.0 * exponent)
    j3 = truncated_taper_integral(tau, 3.0 * exponent)
    for name, value in (("J1", j1), ("J2", j2), ("J3", j3)):
        if not (math.isfinite(value) and value > 0.0):  # pragma: no cover - defensive
            raise AssertionError(f"degenerate taper integral {name} = {value!r}")
    reference_length = parent_length(p)
    section_scale = (2.0 / 3.0) * p.beam * p.draft
    volume = section_scale * reference_length * j1
    waterplane = p.beam * reference_length * j1
    starboard_half, port_half = p.half_breadths()
    transverse_inertia = (
        (8.0 * starboard_half**3 + 8.0 * port_half**3) / 24.0 * reference_length * j3
    )
    y_b = (3.0 / 8.0) * p.skew * p.beam * j2 / j1
    y_f = (1.0 / 2.0) * p.skew * p.beam * j2 / j1
    taper_at_transom = _taper_at(p, 0.0)
    return {
        "parent_length": reference_length,
        "beam": p.beam,
        "starboard_half_breadth": starboard_half,
        "port_half_breadth": port_half,
        "displacement_volume": volume,
        "displacement_mass": volume * p.water_density,
        "waterplane_area": waterplane,
        "longitudinal_centre_buoyancy": reference_length * m1 / j1 - tau * reference_length,
        "longitudinal_centre_flotation": reference_length * m1 / j1 - tau * reference_length,
        "vertical_centre_buoyancy": analytic_vertical_centre_buoyancy(p.draft),
        "transverse_centre_buoyancy": y_b,
        "transverse_centre_flotation": y_f,
        "transverse_waterplane_inertia": transverse_inertia,
        "transverse_waterplane_inertia_centroid": transverse_inertia - waterplane * y_f * y_f,
        "block_coefficient": (2.0 / 3.0) * j1 / (1.0 - tau),
        "waterplane_coefficient": j1 / (1.0 - tau),
        "maximum_section_area": section_scale,
        "transom_area": section_scale * taper_at_transom,
        "transom_starboard_half_breadth": starboard_half * taper_at_transom,
        "transom_port_half_breadth": port_half * taper_at_transom,
        "taper_integral_j1": j1,
        "taper_integral_j2": j2,
        "taper_integral_j3": j3,
    }


def mesh_transverse_centroid(mesh: Mesh) -> float:
    """``(1/V) int y dV`` by signed tetrahedra from the origin.

    Independent of the parametric closed forms and of the station quadrature: each
    triangle with the origin forms a tetrahedron whose volume and ``y`` first
    moment are evaluated exactly (the centroid of a tetrahedron is the mean of its
    four vertices). The origin is not required to be inside the body -- the signs
    of the tetrahedra are what carry the orientation.
    """
    volume = 0.0
    moment = 0.0
    for face in mesh.faces:
        a, b, c = (mesh.vertices[index] for index in face)
        cross = (
            b[1] * c[2] - b[2] * c[1],
            b[2] * c[0] - b[0] * c[2],
            b[0] * c[1] - b[1] * c[0],
        )
        tet = (a[0] * cross[0] + a[1] * cross[1] + a[2] * cross[2]) / 6.0
        volume += tet
        moment += tet * (a[1] + b[1] + c[1]) / 4.0
    if abs(volume) < 1e-300:  # pragma: no cover - defensive
        raise ValueError("degenerate mesh: zero signed volume")
    return moment / volume


def simpson_reference(
    parameters: AsymmetricTransomParameters, panels: int = 2000
) -> dict[str, float]:
    """An independent route to the same numbers: composite Simpson in ``x``.

    The closed forms reduce the surface to incomplete taper integrals; this route
    never touches them -- it evaluates the analytic section expressions directly on
    a fine station grid. Agreement between the two is agreement between a series
    and a quadrature, i.e. a check on both.
    """
    if not (isinstance(panels, int) and panels >= 2 and panels % 2 == 0):
        raise ValueError(f"panels must be an even integer >= 2, got {panels!r}")
    p = parameters
    p.validate()
    starboard_half, port_half = p.half_breadths()
    reference_length = parent_length(p)
    h = p.length / panels
    volume = 0.0
    y_moment = 0.0
    squared = 0.0
    for k in range(panels + 1):
        x = k * h
        taper = _taper_at(p, x)
        weight = 1.0 if k in (0, panels) else (4.0 if k % 2 else 2.0)
        section_area = (starboard_half + port_half) * p.draft * taper * (2.0 / 3.0)
        section_y_moment = (starboard_half**2 - port_half**2) * p.draft * taper * taper / 4.0
        volume += weight * section_area
        y_moment += weight * section_y_moment
        squared += weight * (starboard_half**3 + port_half**3) * taper**3 / 3.0
    volume *= h / 3.0
    y_moment *= h / 3.0
    transom_area = (starboard_half + port_half) * p.draft * _taper_at(p, 0.0) * (2.0 / 3.0)
    # The transom plate at x=0 is one end of the composite rule; the pointed stem
    # contributes nothing, so the Simpson total is already the body's volume.
    del reference_length
    return {
        "displacement_volume": volume,
        "transverse_centre_buoyancy": y_moment / volume,
        "transverse_waterplane_inertia": squared * h / 3.0,
        "transom_area": transom_area,
    }


def _dimension_target(parameters: AsymmetricTransomParameters) -> DimensionTarget:
    p = parameters
    achievable = max(_taper_at(p, p.length * i / (p.nx - 1)) for i in range(p.nx))
    # The bounding box spans y_starboard + y_port = beam * taper, so the *overall*
    # breadth is still `beam`; the asymmetry slides the box, it does not widen it.
    return DimensionTarget(
        length=p.length, beam=p.beam * achievable, draft=p.draft, waterline_z=0.0
    )


def _asymmetry_checks(
    parameters: AsymmetricTransomParameters, mesh: Mesh, analytic: dict[str, float]
) -> list[Check]:
    p = parameters
    checks: list[Check] = []
    measured_plate = 0.0
    for face in mesh.faces:
        a, b, c = (mesh.vertices[index] for index in face)
        if all(abs(vertex[0]) <= 1e-12 for vertex in (a, b, c)):
            measured_plate += 0.5 * abs(
                (b[1] - a[1]) * (c[2] - a[2]) - (c[1] - a[1]) * (b[2] - a[2])
            )
    expected = analytic["transom_area"]
    relative = (measured_plate - expected) / expected if expected > 0.0 else float("nan")
    checks.append(
        Check(
            "transom_plate_area",
            math.isfinite(relative) and abs(relative) <= 0.05,
            f"transom plate {measured_plate:.6g} m^2 vs analytic section {expected:.6g} m^2 "
            f"({100.0 * relative:+.2f} %; the piecewise-linear section under-reads the parabola)",
            severity=WARNING,
        )
    )
    box = mesh.bounding_box()
    measured_span = box.maximum[1] - box.minimum[1]
    expected_span = analytic["beam"] * max(
        _taper_at(p, p.length * i / (p.nx - 1)) for i in range(p.nx)
    )
    checks.append(
        Check(
            "breadth_span_preserved",
            abs(measured_span - expected_span) <= 1e-12 * max(1.0, expected_span),
            f"bounding-box breadth {measured_span:.6g} m vs declared overall breadth "
            f"{expected_span:.6g} m (the asymmetry must not widen the box)",
        )
    )
    starboard_half, port_half = p.half_breadths()
    checks.append(
        Check(
            "sides_match_declared_skew",
            abs(starboard_half - port_half - p.skew * p.beam) <= 1e-12 * max(1.0, p.beam),
            f"starboard half breadth {starboard_half:.6g} m, port {port_half:.6g} m, "
            f"difference {starboard_half - port_half:+.6g} m vs beam*skew "
            f"{p.skew * p.beam:+.6g} m",
        )
    )
    return checks


def generate_asymmetric_transom_hull(
    parameters: AsymmetricTransomParameters | None = None,
) -> AsymmetricTransomModel:
    """Generate and validate an asymmetric transom hull.

    Raises :class:`ValueError` on invalid parameters or when the surface fails a
    geometric validity check -- an unwatertight or non-manifold stern is never
    handed back silently.
    """
    p = parameters if parameters is not None else AsymmetricTransomParameters()
    mesh, _offsets = build_asymmetric_mesh(p)
    analytic = analytic_asymmetric_hydrostatics(p)
    hydrostatics = compute_hydrostatics(
        mesh, mesh.bounding_box().maximum[2], water_density=p.water_density, closed=True
    )
    transverse = mesh_transverse_centroid(mesh)
    pct_displacement = 100.0 * (hydrostatics.displacement_volume / analytic["displacement_volume"] - 1.0)
    report = merge(
        inspect_geometry(mesh),
        check_dimensions(mesh, _dimension_target(p)),
        _asymmetry_checks(p, mesh, analytic),
        [
            Check(
                "analytic_displacement_match",
                abs(hydrostatics.displacement_volume - analytic["displacement_volume"])
                <= SAMPLING_ERROR_WARNING * analytic["displacement_volume"],
                f"{p.nx}x{p.nz} grid integrates {hydrostatics.displacement_volume:.6g} m^3 against "
                f"an analytic {analytic['displacement_volume']:.6g} m^3 ({pct_displacement:+.2f} %)",
                severity=WARNING,
            ),
            Check(
                "transverse_centre_buoyancy_direction",
                (p.skew == 0.0 and abs(transverse) <= 1e-4 * p.beam)
                or (p.skew * transverse > 0.0),
                f"measured y_B = {transverse:+.6g} m for skew {p.skew:+.3f} "
                f"(the wider side must carry the buoyancy centre)",
            ),
        ],
    )
    errors = report.errors
    if errors:
        raise ValueError("generated asymmetric transom hull failed validation: " + "; ".join(errors))
    return AsymmetricTransomModel(
        parameters=p,
        mesh=mesh,
        hydrostatics=hydrostatics,
        analytic=analytic,
        validation=report,
        transverse_centre_buoyancy=transverse,
    )