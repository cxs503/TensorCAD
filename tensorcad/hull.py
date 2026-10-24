"""A deterministic parametric hull surface, closed into a usable solid.

This first TensorCAD model generates a symmetric displacement-hull surface from
principal dimensions, lofts it into a triangle mesh, caps it at the design
waterline, orients it outward, and integrates its hydrostatics. The result is a
watertight, edge-manifold, positively-oriented closed surface -- which is the
minimum a mesh needs before a volume mesher or solver pre-processor will accept
it. It is still *not* a validated hull form: the breadth distribution is analytic
and intentionally simple, and nothing here checks fairness, stability, resistance
or class rules.

The mesh is built from one closed section ring per station rather than from two
mirrored sheets welded afterwards. Mirrored sheets plus tolerance welding were the
source of two real defects: inverted winding on one side (a mesh whose signed
volume came out exactly zero) and tolerance-dependent merging at the keel. Sharing
vertices by construction means the topology is reproducible bit-for-bit and the
Euler characteristic comes out at 2 without any repair step.

Only the Python standard library is required.
"""
from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field
from typing import Any

from . import mesh as meshlib
from .geometry import Mesh, MeshBuilder
from .hydrostatics import Hydrostatics, compute_hydrostatics
from .parameters import (
    FLOAT,
    INT,
    ParameterSpec,
    RatioSpec,
    check_ratios,
    validate,
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

UNITS = "m"
COORDINATE_FRAME = (
    "x=longitudinal, from the aft end (x=0) toward the forward end (x=length); "
    "y=transverse, positive to starboard; z=vertical, positive up; "
    "z=0 is the design waterline, the keel lies at z=-draft"
)

#: Closed-form integral of the longitudinal taper over the length,
#: ``int_0^1 sin(pi*x)^0.7 dx``, used to check sampled volumes analytically.
TAPER_INTEGRAL = (2.0 / math.pi) * math.sqrt(math.pi) * math.gamma(0.85) / (2.0 * math.gamma(1.35))

#: Exponent on the sine taper. ``< 1`` keeps the ends sharply pointed, ``> 1``
#: fattens the run of shoulders; the keel-side ``sqrt`` profile is fixed.
BREADTH_EXPONENT = 0.7

MAX_LONGITUDINAL_STATIONS = 4001
MAX_VERTICAL_LEVELS = 513
#: Above this mesh size the refinement used for the convergence estimate is
#: skipped rather than doubling an already expensive generation.
MAX_MESH_SIZE_FOR_CONVERGENCE = 400_000

PARAMETER_SPECS: tuple[ParameterSpec, ...] = (
    ParameterSpec("length", FLOAT, minimum=0.0, maximum=1.0e4, unit=UNITS,
                  description="longitudinal extent, aft end to forward end",
                  exclusive_minimum=True),
    ParameterSpec("beam", FLOAT, minimum=0.0, maximum=1.0e3, unit=UNITS,
                  description="maximum moulded breadth at the design waterline",
                  exclusive_minimum=True),
    ParameterSpec("draft", FLOAT, minimum=0.0, maximum=1.0e3, unit=UNITS,
                  description="keel to design waterline", exclusive_minimum=True),
    ParameterSpec("nx", INT, minimum=3.0, maximum=float(MAX_LONGITUDINAL_STATIONS),
                  unit="stations", description="longitudinal station count, including both ends"),
    ParameterSpec("nz", INT, minimum=3.0, maximum=float(MAX_VERTICAL_LEVELS),
                  unit="levels", description="vertical levels from keel to waterline, inclusive"),
    ParameterSpec("water_density", FLOAT, minimum=300.0, maximum=2000.0, unit="kg/m^3",
                  description="density used to convert displaced volume to mass"),
)

RATIO_SPECS: tuple[RatioSpec, ...] = (
    RatioSpec("length_to_beam", "length", "beam", minimum=1.5, maximum=15.0,
              guidance="displacement hulls are usually 2.5-8"),
    RatioSpec("beam_to_draft", "beam", "draft", minimum=1.0, maximum=12.0,
              guidance="a shallow-draft hull has large beam/draft"),
    RatioSpec("draft_to_length", "draft", "length", minimum=0.01, maximum=1.0,
              guidance="beyond ~0.4 the form is a barge, not a hull"),
)


@dataclass(frozen=True)
class HullParameters:
    """Principal dimensions, sampling resolution, and the fluid for displacement.

    Construction never fails: call :meth:`validate` (or use a generator, which
    calls it for you) to get every problem at once. ``water_density`` is kept here
    rather than passed to a reporting call so that a stored parameter set alone
    reproduces the numbers in the metadata (1025 is seawater, 1000 is fresh).
    """

    length: float = 10.0
    beam: float = 2.0
    draft: float = 0.8
    nx: int = 41
    nz: int = 13
    water_density: float = 1025.0

    def as_dict(self) -> dict[str, Any]:
        data: dict[str, Any] = asdict(self)
        return data

    def _values(self) -> dict[str, Any]:
        return {spec.name: getattr(self, spec.name) for spec in PARAMETER_SPECS}

    def problems(self) -> tuple[list[str], list[str]]:
        """Return ``(errors, warnings)`` for this parameter set."""
        values = self._values()
        errors = validate(PARAMETER_SPECS, values)
        warnings: list[str] = []
        if not errors:
            warnings = check_ratios(RATIO_SPECS, {k: float(v) for k, v in values.items()})
        return errors, warnings

    def validate(self) -> None:
        """Raise :class:`ValueError` naming every problem found, not just the first."""
        errors, _ = self.problems()
        if errors:
            raise ValueError("; ".join(errors))

    @property
    def grid_is_symmetric(self) -> bool:
        """True when the station grid lands exactly on the widest section.

        With an even ``nx`` no station sits at midship, so the lofted surface
        misses the maximum of the taper and the realised beam falls short of the
        requested beam. The deficit is reported, not silently absorbed.
        """
        return self.nx % 2 == 1


def _longitudinal_taper(parameters: HullParameters, index: int) -> float:
    """Half-breadth taper at station ``index``, forced to zero at both ends.

    ``math.sin(math.pi)`` is not exactly zero in floating point, so the end
    stations are set to zero explicitly: the ends must close on the centreline or
    the loft leaves a slot at the stem and post.
    """
    if index == 0 or index == parameters.nx - 1:
        return 0.0
    # math.pow rather than ** because ``float ** float`` is typed Any in typeshed.
    return math.pow(math.sin(math.pi * index / (parameters.nx - 1)), BREADTH_EXPONENT)


def _vertical_taper(parameters: HullParameters, level: int) -> float:
    """Half-breadth scale at vertical ``level`` (0 = keel, nz-1 = waterline).

    ``sqrt`` gives the hard chine-like knuckle at the keel that keeps the sampled
    surface single-valued in z; it also makes the section derivative singular
    there, which is why volume convergence is slower than the mesh density.
    """
    return math.sqrt(level / (parameters.nz - 1))


@dataclass(frozen=True)
class HullModel:
    """A generated hull surface plus everything needed to reproduce it."""

    parameters: HullParameters
    mesh: Mesh
    offsets: tuple[tuple[float, float], ...]
    validation: ValidationReport = field(default_factory=ValidationReport)
    hydrostatics: Hydrostatics | None = None

    # Back-compat accessors: callers of the 0.1 API read these two directly.
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

    def metadata(self, *, with_discretisation: bool = False) -> dict[str, Any]:
        data: dict[str, Any] = {
            "schema_version": "2.0",
            "model_type": "closed_sampled_parametric_displacement_hull",
            "parameters": self.parameters.as_dict(),
            "units": self.units,
            "coordinate_frame": self.coordinate_frame,
            "generator": {
                "kind": "analytic_lofted_surface",
                "longitudinal_taper": f"sin(pi*x/length)^{BREADTH_EXPONENT}",
                "vertical_taper": "sqrt((z+draft)/draft)",
                "topology": "one closed section ring per station, lofted, capped at the design waterline",
                "backend": "tensorcad.geometry (standard library)",
            },
            "mesh": self.mesh.metadata(),
            "offset_count": len(self.offsets),
            "validation": self.validation.to_dict(),
            "engineering_status": (
                "geometry and hydrostatic integration only; fairness, lightship weight, "
                "stability, resistance, structural and classification-rule checks are not implemented"
            ),
        }
        if self.hydrostatics is not None:
            data["hydrostatics"] = self.hydrostatics.to_dict()
        if with_discretisation:
            data["discretisation"] = estimate_discretisation_error(self.parameters)
        return data


def _section_point(parameters: HullParameters, i: int, j: int) -> tuple[float, float, float]:
    """Ring position ``j`` at station ``i``.

    ``j`` walks the section outline once: down the starboard side from the design
    waterline to the keel, then up the port side back to the waterline. The keel
    point belongs to the ring once, so both flanks share it instead of
    duplicating it -- that is what makes the lofted surface closed without any
    welding step.
    """
    p = parameters
    x = p.length * i / (p.nx - 1)
    keelward = j <= p.nz - 1
    level = (p.nz - 1 - j) if keelward else (j - (p.nz - 1))
    # Written as (1 - fraction) so the top level is exactly z = 0: -draft +
    # draft * 1.0 does not cancel in binary floating point and would leave the
    # waterline vertices a few ulp apart, which breaks coplanarity checks.
    z = -p.draft * (1.0 - level / (p.nz - 1))
    half_breadth = 0.5 * p.beam * _longitudinal_taper(p, i) * _vertical_taper(p, level)
    y = half_breadth if keelward else -half_breadth
    return (x, y, z)


def build_mesh(parameters: HullParameters) -> tuple[Mesh, tuple[tuple[float, float], ...]]:
    """Loft and cap the hull surface; returns the mesh and waterline offsets.

    The generator is deterministic: identical parameters produce identical
    vertex and face tuples, so a simulation run can be replayed from its
    parameters alone.
    """
    parameters.validate()
    p = parameters
    ring_size = 2 * p.nz - 1
    builder = MeshBuilder(units=UNITS, coordinate_frame=COORDINATE_FRAME)

    rings: list[list[int]] = []
    for i in range(p.nx):
        rings.append([builder.add_vertex(_section_point(p, i, j)) for j in range(ring_size)])

    # Shell: loft consecutive rings. Consecutive ring indices are consecutive
    # points on the same section outline, so every quad is a real panel.
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


def check_consistency(
    parameters: HullParameters, mesh: Mesh, offsets: tuple[tuple[float, float], ...]
) -> list[Check]:
    """Checks that tie the generated surface back to the requested parameters."""
    p = parameters
    checks: list[Check] = []
    checks.append(
        Check(
            "offset_per_station",
            len(offsets) == p.nx,
            f"{len(offsets)} offset rows for {p.nx} stations",
        )
    )
    ends_ok = bool(offsets) and offsets[0][1] == 0.0 and offsets[-1][1] == 0.0
    checks.append(
        Check("offsets_close_at_both_ends", ends_ok,
              "" if ends_ok else "breadth must be exactly zero at the stem and the post")
    )
    bad = [index for index, (x, width) in enumerate(offsets) if not math.isfinite(width) or width < 0]
    checks.append(
        Check("offsets_finite_non_negative", not bad,
              "" if not bad else f"{len(bad)} invalid offset breadths")
    )
    ordered = all(offsets[k][0] <= offsets[k + 1][0] for k in range(len(offsets) - 1)) if offsets else False
    checks.append(
        Check(
            "offsets_station_order",
            ordered,
            "" if ordered else "stations are not monotone in x",
        )
    )

    realised = max((width for _, width in offsets), default=0.0)
    achievable = max(_longitudinal_taper(p, i) for i in range(p.nx))
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

    values = {name: float(value) for name, value in p._values().items()}
    for spec in RATIO_SPECS:
        problem = spec.check(values)
        checks.append(
            Check(
                f"ratio_{spec.name}",
                problem is None,
                problem or f"{spec.numerator}/{spec.denominator} within [{spec.minimum:g}, {spec.maximum:g}]",
                severity=WARNING,
            )
        )
    return checks


def _dimension_target(parameters: HullParameters) -> DimensionTarget:
    return DimensionTarget(
        length=parameters.length,
        beam=parameters.beam * max(_longitudinal_taper(parameters, i) for i in range(parameters.nx)),
        draft=parameters.draft,
        waterline_z=0.0,
    )


def generate_hull(parameters: HullParameters | None = None) -> HullModel:
    """Generate a closed, validated hull surface from principal dimensions.

    Raises :class:`ValueError` when the parameters are invalid or when the
    generated surface fails any geometric validity check -- an unusable mesh is
    never handed back silently.
    """
    p = parameters if parameters is not None else HullParameters()
    mesh, offsets = build_mesh(p)
    report = merge(
        inspect_geometry(mesh),
        check_dimensions(mesh, _dimension_target(p)),
        check_consistency(p, mesh, offsets),
    )
    errors = report.errors
    if errors:
        raise ValueError("generated hull failed validation: " + "; ".join(errors))
    waterline_z = mesh.bounding_box().maximum[2]
    hydrostatics = compute_hydrostatics(
        mesh, waterline_z, water_density=p.water_density, closed=True
    )
    return HullModel(parameters=p, mesh=mesh, offsets=offsets, validation=report, hydrostatics=hydrostatics)


def validate_hull(model: HullModel) -> list[str]:
    """Back-compat helper: the list of error strings for an already built model."""
    return model.validation.errors


def estimate_discretisation_error(parameters: HullParameters) -> dict[str, Any]:
    """Compare the requested resolution against a uniformly refined one.

    A sampled surface is only an approximation of the analytic form, and a
    downstream solver inherits that error. Regenerating at twice the density and
    reporting the relative change in integrated quantities is cheap, honest
    evidence of convergence -- far better than implying that mesh numbers are
    exact.
    """
    p = parameters
    if p.nx * p.nz > MAX_MESH_SIZE_FOR_CONVERGENCE:
        return {
            "status": "skipped",
            "reason": f"mesh size {p.nx * p.nz} exceeds {MAX_MESH_SIZE_FOR_CONVERGENCE} vertices per side",
        }
    refined = HullParameters(p.length, p.beam, p.draft, 2 * p.nx - 1, 2 * p.nz - 1)
    coarse_mesh, _ = build_mesh(p)
    fine_mesh, _ = build_mesh(refined)
    wl_coarse = coarse_mesh.bounding_box().maximum[2]
    wl_fine = fine_mesh.bounding_box().maximum[2]
    coarse_values = compute_hydrostatics(coarse_mesh, wl_coarse, closed=True)
    fine_values = compute_hydrostatics(fine_mesh, wl_fine, closed=True)

    def relative(name: str) -> float:
        fine = float(getattr(fine_values, name))
        coarse = float(getattr(coarse_values, name))
        if not math.isfinite(fine) or abs(fine) < 1e-15:
            return float("nan")
        return abs(coarse - fine) / abs(fine)

    return {
        "status": "measured",
        "method": "Richardson-style comparison against a uniformly refined grid",
        "coarse_grid": {"nx": p.nx, "nz": p.nz},
        "refined_grid": {"nx": refined.nx, "nz": refined.nz},
        "relative_change": {
            "displacement_volume": relative("displacement_volume"),
            "wetted_surface_area": relative("wetted_surface_area"),
            "waterplane_area": relative("waterplane_area"),
        },
        "note": "relative change is an error estimate for the coarse grid, not a proof of convergence",
    }
