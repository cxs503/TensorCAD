"""Parametric PC6 ice-class hull form: raked stem, flared topsides, transom stern.

The campaign design point is L = 200 m, B = 24 m, T = 11 m. The body is described
by two analytic ratios and lofted into a closed, watertight surface; the same body
is available as an implicit field so a Cartesian solver never has to consume the
triangulation.

Coordinate frame (identical to :mod:`tensorcad.hull`)::

    x = longitudinal, x = 0 at the aft perpendicular (transom), x = L at the stem
        at the design waterline
    y = transverse, positive to starboard
    z = vertical, positive up, z = 0 at the design waterline, keel at z = -T

The solid is closed at the design waterline: its "deck" is the waterplane. That is
the body a resistance or seakeeping solver meshes, and it is stated explicitly
because it is *not* a complete above-water hull -- topside flare above the DWL is
outside this contract.

Geometry
--------
Stations are parameterised by ``s`` in [0, 1], heights by ``u`` in [0, 1]::

    z(u)    = -T * (1 - u)
    x(s, u) = s * x_stem(z),      x_stem(z) = L + stem_rake * z / T
    y(s, u) = +-0.5 * B * g(s) * sec(u)

``g`` is the longitudinal half-breadth ratio (transom ratio ``t0`` at ``s = 0``,
unity at the midship fraction ``s_m``, zero at the stem) and ``sec`` is the section
ratio (flat of bottom ``f_b`` at ``u = 0``, wall-sided above the bilge height, DWL
breadth of unity at ``u = 1`` with a cubic flare term between). ``g(1) = 0`` while
``sec`` stays positive, so the starboard and port sheets meet exactly on the raked
stem curve and share those vertices by construction: the loft closes without
welding, and no tolerance decides the topology.

Named limitations (see ``docs/pc6.md``): the form is single-valued in ``y`` at
every station, so it cannot represent a bulbous bow or a reentrant stern; there is
no skeg, no bilge keel, no appendage; the entrance/run exponents meet with a
kink in slope at the maximum section (C0, not C1).

Only the Python standard library is required.
"""
from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass, replace
from typing import Any

from . import mesh as meshlib
from .geometry import Mesh, MeshBuilder
from .hydrostatics import Hydrostatics, compute_hydrostatics
from .intersect import find_self_intersections
from .validation import Check, ValidationReport, inspect_geometry, vertex_link_defects

UNITS = "m"
COORDINATE_FRAME = (
    "x=longitudinal from the aft perpendicular (x=0) to the stem (x=L); "
    "y=transverse, starboard positive; z=vertical, positive up; "
    "z=0 is the design waterline and the keel lies at z=-draft"
)
SEAWATER_DENSITY = 1025.0
GRAVITY = 9.80665

#: Relative tolerance on the sampled-vs-analytic displacement. The default grid
#: lands two orders of magnitude inside it (``tests/test_pc6.py`` pins the value).
VOLUME_TOLERANCE = 0.005
#: Relative tolerance on the bounding box against L, B, T. The station and level
#: lists insert the maximum section and the bilge break exactly, so the sampled
#: solid reproduces the three principal dimensions to rounding.
DIMENSION_TOLERANCE = 1e-9

MAX_LONGITUDINAL_STATIONS = 801
MAX_VERTICAL_LEVELS = 257


@dataclass(frozen=True)
class PC6Parameters:
    """Principal dimensions and form coefficients of the PC6 hull.

    Construction never fails: :meth:`validate` returns every problem at once and
    :func:`generate` calls it.
    """

    length: float = 200.0
    beam: float = 24.0
    draft: float = 11.0
    #: Horizontal setback of the stem at the keel relative to the stem at the DWL
    #: (a forward-raked stem: the bow leans forward going up).
    stem_rake: float = 8.0
    #: Transom half-breadth as a fraction of the maximum half-breadth.
    transom_ratio: float = 0.88
    #: Half-breadth at the keel as a fraction of the DWL half-breadth: the flat of
    #: bottom. Zero collapses the section onto a V and the keel becomes a line.
    bilge_ratio: float = 0.95
    #: Fraction of the draft over which the bilge turns into wall-sided topsides.
    bilge_height: float = 0.08
    #: Cubic flare coefficient of the topsides; 0 is wall-sided down to the DWL.
    flare: float = 0.15
    #: Exponent of the aft run (> 1 gives a fair, slack run).
    aft_exponent: float = 1.1
    #: Exponent of the fore entrance (< 1 gives a full entrance, > 1 a fine one).
    fore_exponent: float = 0.22
    #: Longitudinal position of the maximum section as a fraction of L.
    midship_fraction: float = 0.62
    longitudinal_stations: int = 97
    vertical_levels: int = 19
    water_density: float = SEAWATER_DENSITY

    @property
    def transom_breadth(self) -> float:
        return self.transom_ratio * self.beam

    def validate(self) -> list[str]:
        """Every problem with this parameter set, as human-readable strings."""
        errors: list[str] = []
        for name, value in (
            ("length", self.length),
            ("beam", self.beam),
            ("draft", self.draft),
            ("water_density", self.water_density),
        ):
            if not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0.0:
                errors.append(f"{name} must be a finite positive number, got {value!r}")
        for name, value in (
            ("midship_fraction", self.midship_fraction),
            ("bilge_height", self.bilge_height),
            ("transom_ratio", self.transom_ratio),
        ):
            if not math.isfinite(value) or not 0.0 < value < 1.0:
                errors.append(f"{name} must lie strictly inside (0, 1), got {value!r}")
        if not math.isfinite(self.bilge_ratio) or not 0.0 <= self.bilge_ratio < 1.0:
            errors.append(
                f"bilge_ratio must lie in [0, 1), got {self.bilge_ratio!r}: a negative flat "
                "of bottom folds the section through the centreline"
            )
        if not math.isfinite(self.flare) or not -1.0 < self.flare <= 3.0:
            errors.append(f"flare must lie in (-1, 3], got {self.flare!r}")
        if not math.isfinite(self.stem_rake) or self.stem_rake < 0.0:
            errors.append(f"stem_rake must be finite and non-negative, got {self.stem_rake!r}")
        if (
            math.isfinite(self.length)
            and self.length > 0.0
            and math.isfinite(self.stem_rake)
            and self.stem_rake > 0.35 * self.length
        ):
            errors.append(
                f"stem_rake={self.stem_rake!r} m exceeds 35 % of length={self.length!r} m: the "
                "stem curve passes aft of the transom and the loft folds through itself"
            )
        for name, value in (
            ("aft_exponent", self.aft_exponent),
            ("fore_exponent", self.fore_exponent),
        ):
            if not math.isfinite(value) or not 0.05 <= value <= 8.0:
                errors.append(f"{name} must lie in [0.05, 8], got {value!r}")
        for name, value, ceiling in (
            ("longitudinal_stations", self.longitudinal_stations, MAX_LONGITUDINAL_STATIONS),
            ("vertical_levels", self.vertical_levels, MAX_VERTICAL_LEVELS),
        ):
            if not isinstance(value, int) or not 5 <= value <= ceiling:
                errors.append(f"{name} must be an integer in [5, {ceiling}], got {value!r}")
        if not errors:
            beam_to_draft = self.beam / self.draft
            if not 2.0 <= beam_to_draft <= 8.0:
                errors.append(f"beam/draft={beam_to_draft:.3f} is outside [2, 8] for a PC6 hull form")
            slenderness = self.length / self.beam
            if not 5.0 <= slenderness <= 10.0:
                errors.append(f"length/beam={slenderness:.3f} is outside [5, 10] for a PC6 hull form")
        return errors

    def require_valid(self) -> None:
        errors = self.validate()
        if errors:
            raise ValueError("invalid PC6 parameters: " + "; ".join(errors))

    def to_dict(self) -> dict[str, Any]:
        return {
            "length": {"value": self.length, "unit": UNITS,
                       "description": "length between perpendiculars L"},
            "beam": {"value": self.beam, "unit": UNITS,
                     "description": "moulded breadth at the design waterline B"},
            "draft": {"value": self.draft, "unit": UNITS,
                      "description": "keel to design waterline T"},
            "stem_rake": {"value": self.stem_rake, "unit": UNITS,
                          "description": "horizontal setback of the stem at the keel (forward rake)"},
            "transom_ratio": {"value": self.transom_ratio, "unit": "1",
                              "description": "transom half-breadth / maximum half-breadth"},
            "bilge_ratio": {"value": self.bilge_ratio, "unit": "1",
                            "description": "half-breadth at the keel / DWL half-breadth"},
            "bilge_height": {"value": self.bilge_height, "unit": "1",
                             "description": "draft fraction over which the bilge turns"},
            "flare": {"value": self.flare, "unit": "1",
                      "description": "cubic topside flare coefficient"},
            "aft_exponent": {"value": self.aft_exponent, "unit": "1",
                             "description": "run fullness exponent"},
            "fore_exponent": {"value": self.fore_exponent, "unit": "1",
                              "description": "entrance fullness exponent"},
            "midship_fraction": {"value": self.midship_fraction, "unit": "1",
                                 "description": "position of the maximum section as a fraction of L"},
            "longitudinal_stations": {"value": self.longitudinal_stations, "unit": "count",
                                      "description": "requested station count (breaks inserted on top)"},
            "vertical_levels": {"value": self.vertical_levels, "unit": "count",
                                "description": "requested level count (breaks inserted on top)"},
            "water_density": {"value": self.water_density, "unit": "kg/m^3",
                              "description": "density turning displacement volume into mass"},
        }


# ---------------------------------------------------------------------------
# analytic surface description
# ---------------------------------------------------------------------------
def longitudinal_ratio(parameters: PC6Parameters, s: float) -> float:
    """``g(s)``: half-breadth ratio along the length, unity at the maximum section."""
    if s <= 0.0:
        return parameters.transom_ratio
    if s >= 1.0:
        return 0.0
    midship = parameters.midship_fraction
    if s <= midship:
        fraction = s / midship
        return float(
            parameters.transom_ratio
            + (1.0 - parameters.transom_ratio) * float(fraction) ** parameters.aft_exponent
        )
    fraction = (s - midship) / (1.0 - midship)
    return float((1.0 - fraction) ** parameters.fore_exponent)


def section_ratio(parameters: PC6Parameters, u: float) -> float:
    """``sec(u)``: half-breadth ratio over the draft, unity at the design waterline."""
    if u <= 0.0:
        base = parameters.bilge_ratio
    elif u >= parameters.bilge_height:
        base = 1.0
    else:
        base = parameters.bilge_ratio + (1.0 - parameters.bilge_ratio) * math.sqrt(
            u / parameters.bilge_height
        )
    return base * (1.0 + parameters.flare * u ** 3) / (1.0 + parameters.flare)


def stem_x(parameters: PC6Parameters, z: float) -> float:
    """Longitudinal position of the raked stem at height ``z``."""
    return parameters.length + parameters.stem_rake * (z / parameters.draft)


def half_breadth(parameters: PC6Parameters, s: float, u: float) -> float:
    """Starboard half-breadth at longitudinal parameter ``s`` and height ``u``."""
    return 0.5 * parameters.beam * longitudinal_ratio(parameters, s) * section_ratio(parameters, u)


# ---------------------------------------------------------------------------
# analytic integrals: the independent reference the sampled mesh is judged against
# ---------------------------------------------------------------------------
def _poly_mul(a: Sequence[float], b: Sequence[float]) -> list[float]:
    out = [0.0] * (len(a) + len(b) - 1)
    for i, left in enumerate(a):
        if left == 0.0:
            continue
        for j, right in enumerate(b):
            out[i + j] += left * right
    return out


def _power_poly(a: float, b: float, power: int) -> list[float]:
    result = [1.0]
    for _ in range(power):
        result = _poly_mul(result, [a, b])
    return result


def _integrate_section(parameters: PC6Parameters, poly: Sequence[float]) -> float:
    """``int_0^1 sec(u) * poly(u) du`` in exact power terms.

    ``sec`` is piecewise ``flat + (1-flat) * sqrt(u/u_bilge)`` times a cubic flare
    factor, so the integral is a finite sum of powers -- no quadrature, and no code
    path shared with the mesh integration.
    """
    flared = [0.0] * (len(poly) + 3)
    for index, coefficient in enumerate(poly):
        flared[index] += coefficient
        flared[index + 3] += parameters.flare * coefficient
    bilge = parameters.bilge_height
    flat = parameters.bilge_ratio
    total = 0.0
    for index, coefficient in enumerate(flared):
        if coefficient == 0.0:
            continue
        power = index + 1
        below = coefficient * bilge ** power / power
        root = coefficient * bilge ** (index + 1.5) / (index + 1.5)
        above = coefficient * (1.0 - bilge ** power) / power
        total += flat * below + (1.0 - flat) * root / math.sqrt(bilge) + above
    return total / (1.0 + parameters.flare)


def longitudinal_moments(parameters: PC6Parameters) -> tuple[float, float]:
    """``int g ds`` and ``int s g ds`` over the length."""
    midship = parameters.midship_fraction
    transom = parameters.transom_ratio
    aft = parameters.aft_exponent
    fore = parameters.fore_exponent
    area = midship * (transom + (1.0 - transom) / (aft + 1.0)) + (1.0 - midship) / (fore + 1.0)
    moment = (
        transom * midship * midship / 2.0
        + (1.0 - transom) * midship * midship / (aft + 2.0)
        + midship * (1.0 - midship) / (fore + 1.0)
        + (1.0 - midship) ** 2 / ((fore + 1.0) * (fore + 2.0))
    )
    return area, moment


@dataclass(frozen=True)
class AnalyticForm:
    """Closed-form integrals of the analytic surface the mesh generator samples."""

    volume: float
    waterplane_area: float
    lcb: float
    lcf: float
    vcb: float
    midship_section_area: float
    block_coefficient: float
    midship_coefficient: float
    prismatic_coefficient: float
    waterplane_coefficient: float

    def to_dict(self) -> dict[str, float]:
        return {
            "displacement_volume_m3": self.volume,
            "waterplane_area_m2": self.waterplane_area,
            "lcb_from_aft_m": self.lcb,
            "lcf_from_aft_m": self.lcf,
            "vcb_from_dwl_m": self.vcb,
            "midship_section_area_m2": self.midship_section_area,
            "block_coefficient": self.block_coefficient,
            "midship_coefficient": self.midship_coefficient,
            "prismatic_coefficient": self.prismatic_coefficient,
            "waterplane_coefficient": self.waterplane_coefficient,
        }


def analytic_form(parameters: PC6Parameters) -> AnalyticForm:
    """Exact volume, areas and centres of the surface :func:`build_mesh` samples.

    ``V = B T (int g ds) (int sec (A + B u) du)`` with ``A`` the stem position at the
    keel and ``B`` the rake: the rake enters only through the section integral, which
    is why the displacement does not collapse when the stem is raked.
    """
    p = parameters
    area, moment = longitudinal_moments(p)
    keel_x = p.length - p.stem_rake
    rake = p.stem_rake
    i0 = _integrate_section(p, [1.0])
    i1 = _integrate_section(p, [keel_x, rake])
    i2 = _integrate_section(p, _power_poly(keel_x, rake, 2))
    iz = _integrate_section(p, [-p.draft * keel_x, -p.draft * (rake - keel_x), p.draft * rake])
    volume = p.beam * p.draft * area * i1
    waterplane = p.beam * p.length * area
    midship_area = p.beam * p.draft * i0
    return AnalyticForm(
        volume=volume,
        waterplane_area=waterplane,
        lcb=(moment / area) * (i2 / i1),
        lcf=p.length * moment / area,
        vcb=iz / i1,
        midship_section_area=midship_area,
        block_coefficient=volume / (p.length * p.beam * p.draft),
        midship_coefficient=midship_area / (p.beam * p.draft),
        prismatic_coefficient=volume / (p.length * midship_area),
        waterplane_coefficient=waterplane / (p.length * p.beam),
    )


# ---------------------------------------------------------------------------
# mesh generation
# ---------------------------------------------------------------------------
def station_parameters(parameters: PC6Parameters) -> tuple[float, ...]:
    """Station positions in ``s``, with the maximum section inserted exactly.

    Sampling the maximum section exactly is what lets the generated solid reach
    ``g = 1`` and therefore reproduce the requested beam instead of undershooting it
    by the curvature of ``g`` between two stations.
    """
    count = parameters.longitudinal_stations
    values = [index / (count - 1) for index in range(count)]
    values.append(parameters.midship_fraction)
    return tuple(sorted(set(values)))


def height_parameters(parameters: PC6Parameters) -> tuple[float, ...]:
    """Height positions in ``u``, with the bilge break inserted exactly."""
    count = parameters.vertical_levels
    values = [index / (count - 1) for index in range(count)]
    values.append(parameters.bilge_height)
    return tuple(sorted(set(values)))


def build_mesh(parameters: PC6Parameters, *, validate: bool = True) -> Mesh:
    """Loft the closed hull solid: two shell sheets plus bottom, deck and transom.

    The sheets share the vertices on the raked stem curve (half-breadth exactly
    zero), which closes the loft at the bow without welding. ``validate=False``
    skips the parameter check: that is the escape hatch the negative controls in the
    test suite use to prove the *geometric* detector works, never a public path.
    """
    p = parameters
    if validate:
        p.require_valid()
    stations = station_parameters(p)
    levels = height_parameters(p)
    builder = MeshBuilder(units=UNITS, coordinate_frame=COORDINATE_FRAME)
    starboard: list[list[int]] = []
    port: list[list[int]] = []
    for s in stations:
        star_row: list[int] = []
        port_row: list[int] = []
        for u in levels:
            z = -p.draft * (1.0 - u)
            x = s * stem_x(p, z)
            breadth = half_breadth(p, s, u)
            star_row.append(builder.add_vertex((x, breadth, z)))
            # add_vertex deduplicates by exact coordinate, so a zero half-breadth
            # gives both sheets the same index: a crease, not a hole.
            port_row.append(builder.add_vertex((x, -breadth, z)))
        starboard.append(star_row)
        port.append(port_row)

    last_station = len(stations) - 1
    last_level = len(levels) - 1
    for i in range(last_station):
        for j in range(last_level):
            builder.add_quad(
                starboard[i][j], starboard[i][j + 1], starboard[i + 1][j + 1], starboard[i + 1][j]
            )
            builder.add_quad(
                port[i][j], port[i + 1][j], port[i + 1][j + 1], port[i][j + 1]
            )
        builder.add_quad(
            starboard[i][last_level],
            port[i][last_level],
            port[i + 1][last_level],
            starboard[i + 1][last_level],
        )
        builder.add_quad(starboard[i][0], starboard[i + 1][0], port[i + 1][0], port[i][0])
    for j in range(last_level):
        builder.add_quad(starboard[0][j], port[0][j], port[0][j + 1], starboard[0][j + 1])
    return builder.build()


# ---------------------------------------------------------------------------
# implicit field
# ---------------------------------------------------------------------------
def field_value(parameters: PC6Parameters, x: float, y: float, z: float) -> float:
    """Implicit margin ``f(x, y, z)``: negative inside the solid, zero on the shell.

    The sign is exact -- it is the analytic membership test of the same surface
    :func:`build_mesh` samples. The magnitude is a *margin*, not a distance: it is
    minus the smallest slack against the five bounding constraints (deck, keel,
    transom, raked stem, side), so it is a lower bound on the true distance. That is
    enough for inside/outside queries, voxel classification and level-set
    initialisation, and it is said here because a consumer must not read it as a
    signed distance function.
    """
    p = parameters
    margins = [
        -z,  # below the design waterline
        z + p.draft,  # above the keel
        x,  # forward of the transom
        stem_x(p, z) - x,  # aft of the raked stem
    ]
    local_stem = stem_x(p, z)
    if local_stem > 0.0:
        s = x / local_stem
        if 0.0 <= s <= 1.0:
            margins.append(half_breadth(p, s, min(1.0, max(0.0, (z + p.draft) / p.draft))) - abs(y))
        else:
            margins.append(-abs(y) - 1.0)
    else:  # pragma: no cover - excluded by parameter validation
        margins.append(-abs(y) - 1.0)
    return -min(margins)


def is_inside(parameters: PC6Parameters, point: Sequence[float]) -> bool:
    """Membership in the hull solid -- exactly the sign of :func:`field_value`."""
    return field_value(parameters, point[0], point[1], point[2]) < 0.0


def sampled_field(
    parameters: PC6Parameters, *, nx: int = 33, ny: int = 9, nz: int = 7
) -> dict[str, Any]:
    """A coarse sampled copy of the implicit field for the JSON contract.

    The analytic definition travels with the numbers, so a consumer can evaluate the
    field anywhere; the samples exist so a solver that only wants to initialise (or a
    reviewer diffing two revisions) never has to reimplement the formula.
    """
    p = parameters
    bounds = {
        "x": (-0.02 * p.length, 1.02 * p.length),
        "y": (-0.55 * p.beam, 0.55 * p.beam),
        "z": (-1.02 * p.draft, 0.02 * p.draft),
    }
    values: list[float] = []
    for iz in range(nz):
        z = bounds["z"][0] + (bounds["z"][1] - bounds["z"][0]) * iz / (nz - 1)
        for iy in range(ny):
            y = bounds["y"][0] + (bounds["y"][1] - bounds["y"][0]) * iy / (ny - 1)
            for ix in range(nx):
                x = bounds["x"][0] + (bounds["x"][1] - bounds["x"][0]) * ix / (nx - 1)
                values.append(round(field_value(p, x, y, z), 4))
    return {
        "grid": [nx, ny, nz],
        "index_order": "z-major, then y, then x (C order with iz slowest)",
        "bounds": {axis: [round(low, 6), round(high, 6)] for axis, (low, high) in bounds.items()},
        "values": values,
        "values_units": UNITS,
        "definition": {
            "kind": "constraint_margin",
            "sign_convention": "negative inside the solid, zero on the shell, positive outside",
            "formula": "f = -min(-z, z+T, x, x_stem(z)-x, b(s,u)-|y|)",
            "x_stem": "x_stem(z) = L + stem_rake * z / T",
            "section": "b = 0.5*B*longitudinal_ratio(x/x_stem)*section_ratio((z+T)/T)",
            "magnitude": "lower bound on the distance to the nearest boundary, not a signed distance",
        },
    }


# ---------------------------------------------------------------------------
# topology, model, generation
# ---------------------------------------------------------------------------
def _shell_components(mesh: Mesh) -> int:
    """Connected components of the vertex graph -- one per closed shell."""
    parent = list(range(len(mesh.vertices)))

    def find(item: int) -> int:
        while parent[item] != item:
            parent[item] = parent[parent[item]]
            item = parent[item]
        return item

    for face in mesh.faces:
        roots = [find(vertex) for vertex in face]
        for root in roots[1:]:
            if root != roots[0]:
                parent[root] = roots[0]
    return len({find(index) for index in range(len(mesh.vertices))})


def topology(mesh: Mesh) -> dict[str, Any]:
    """Topology numbers a volume mesher cares about, computed from the mesh alone."""
    edges = meshlib.edge_face_counts(mesh)
    return {
        "vertices": len(mesh.vertices),
        "faces": len(mesh.faces),
        "edges": len(edges),
        "watertight": meshlib.is_watertight(mesh),
        "edge_manifold": meshlib.is_edge_manifold(mesh),
        "boundary_edges": len(meshlib.boundary_edges(mesh)),
        "euler_characteristic": len(mesh.vertices) - len(edges) + len(mesh.faces),
        "components": _shell_components(mesh),
        "zero_area_faces": len(meshlib.zero_area_faces(mesh)),
        "vertex_link_defects": len(vertex_link_defects(mesh)),
        "signed_volume_m3": meshlib.total_signed_volume(mesh),
        "surface_area_m2": meshlib.surface_area(mesh),
    }


@dataclass(frozen=True)
class PC6Model:
    """A generated hull plus every number needed to justify handing it on."""

    parameters: PC6Parameters
    mesh: Mesh
    form: AnalyticForm
    hydrostatics: Hydrostatics
    validation: ValidationReport
    report: dict[str, Any]


def _relative(measured: float, reference: float) -> float:
    if reference == 0.0:
        return math.inf if measured != 0.0 else 0.0
    return (measured - reference) / reference


def generate(
    parameters: PC6Parameters | None = None,
    *,
    check_self_intersection: bool = True,
    volume_tolerance: float = VOLUME_TOLERANCE,
) -> PC6Model:
    """Generate, validate and gate the PC6 hull.

    Raises :class:`ValueError` when the parameters are invalid, when the surface
    fails a topology check, when the sampled displacement misses the analytic value
    by more than ``volume_tolerance`` (default :data:`VOLUME_TOLERANCE`; raise it
    for a deliberately coarse grid), or when the mesh intersects itself. An unusable
    hull is never returned silently.
    """
    p = parameters if parameters is not None else PC6Parameters()
    p.require_valid()
    mesh = build_mesh(p)
    form = analytic_form(p)
    hydrostatics = compute_hydrostatics(mesh, 0.0, water_density=p.water_density, closed=True)
    topo = topology(mesh)
    box = mesh.bounding_box()

    checks: list[Check] = list(inspect_geometry(mesh))
    checks.append(
        Check(
            "watertight",
            bool(topo["watertight"]),
            "" if topo["watertight"] else f"{topo['boundary_edges']} boundary edges",
        )
    )
    checks.append(
        Check(
            "edge_manifold",
            bool(topo["edge_manifold"]),
            "" if topo["edge_manifold"] else "an edge is shared by more than two faces",
        )
    )
    checks.append(
        Check(
            "single_component",
            topo["components"] == 1,
            "" if topo["components"] == 1 else f"{topo['components']} shells",
        )
    )
    checks.append(
        Check(
            "euler_characteristic_genus_zero",
            topo["euler_characteristic"] == 2,
            "" if topo["euler_characteristic"] == 2 else f"chi={topo['euler_characteristic']}",
        )
    )
    checks.append(
        Check(
            "no_vertex_link_defects",
            topo["vertex_link_defects"] == 0,
            f"{topo['vertex_link_defects']} pinched vertices",
        )
    )
    volume_error = _relative(hydrostatics.displacement_volume, form.volume)
    checks.append(
        Check(
            "displacement_matches_analytic",
            abs(volume_error) <= volume_tolerance,
            f"mesh/analytic - 1 = {volume_error:+.3e}, tolerance {volume_tolerance:g}",
        )
    )
    dimension_errors = {
        "length": _relative(box.sizes[0], p.length),
        "beam": _relative(box.sizes[1], p.beam),
        "draft": _relative(box.sizes[2], p.draft),
    }
    checks.append(
        Check(
            "bounding_box_matches_principal_dimensions",
            all(abs(value) <= DIMENSION_TOLERANCE for value in dimension_errors.values()),
            "; ".join(f"{name}: {value:+.3e}" for name, value in dimension_errors.items()),
        )
    )
    intersections: tuple[tuple[int, int], ...] = ()
    if check_self_intersection:
        intersections = find_self_intersections(mesh)
        checks.append(
            Check(
                "no_self_intersection",
                not intersections,
                f"{len(intersections)} intersecting face pairs, first {list(intersections[:3])}",
            )
        )
    report = ValidationReport(tuple(checks))
    if report.errors:
        raise ValueError("generated PC6 hull failed validation: " + "; ".join(report.errors))

    summary: dict[str, Any] = {
        "topology": topo,
        "bounding_box": box.to_dict(),
        "principal_dimensions": {
            "length": box.sizes[0],
            "beam": box.sizes[1],
            "draft": box.sizes[2],
        },
        "dimension_relative_errors": dimension_errors,
        "sampling": {
            "stations": len(station_parameters(p)),
            "levels": len(height_parameters(p)),
        },
        "hydrostatics": hydrostatics.to_dict(),
        "analytic": form.to_dict(),
        "reference_volume_error": volume_error,
        "self_intersections": [list(pair) for pair in intersections],
    }
    return PC6Model(
        parameters=p,
        mesh=mesh,
        form=form,
        hydrostatics=hydrostatics,
        validation=report,
        report=summary,
    )


def beam_sweep(
    beams: Sequence[float] = (24.0, 26.0, 28.0),
    base: PC6Parameters | None = None,
    *,
    volume_tolerance: float = VOLUME_TOLERANCE,
) -> list[PC6Model]:
    """Regenerate the hull for each beam -- the parametric-scaling evidence.

    Every entry goes through the same validated path, so a beam that breaks
    watertightness or the analytic displacement check raises instead of being
    recorded as a pass.
    """
    if not beams:
        raise ValueError("beam_sweep needs at least one beam")
    start = base if base is not None else PC6Parameters()
    return [
        generate(replace(start, beam=float(beam)), volume_tolerance=volume_tolerance)
        for beam in beams
    ]


def inside_by_winding_number(mesh: Mesh, point: Sequence[float]) -> bool:
    """Independent inside test: the solid angle subtended by the mesh.

    Used to check the implicit field against the mesh rather than against itself.
    A closed and consistently oriented surface subtends ``4*pi`` at an interior
    point and ``0`` outside, so this shares no code with the analytic membership
    test -- which is the whole point of having it.
    """
    px, py, pz = float(point[0]), float(point[1]), float(point[2])
    total = 0.0
    for face in mesh.faces:
        ax, ay, az = mesh.vertices[face[0]]
        bx, by, bz = mesh.vertices[face[1]]
        cx, cy, cz = mesh.vertices[face[2]]
        ux, uy, uz = ax - px, ay - py, az - pz
        vx, vy, vz = bx - px, by - py, bz - pz
        wx, wy, wz = cx - px, cy - py, cz - pz
        lu = math.sqrt(ux * ux + uy * uy + uz * uz)
        lv = math.sqrt(vx * vx + vy * vy + vz * vz)
        lw = math.sqrt(wx * wx + wy * wy + wz * wz)
        numerator = ux * (vy * wz - vz * wy) - uy * (vx * wz - vz * wx) + uz * (vx * wy - vy * wx)
        denominator = (
            lu * lv * lw
            + (ux * vx + uy * vy + uz * vz) * lw
            + (ux * wx + uy * wy + uz * wz) * lv
            + (vx * wx + vy * wy + vz * wz) * lu
        )
        total += 2.0 * math.atan2(numerator, denominator)
    return abs(total) > 2.0 * math.pi


__all__ = [
    "COORDINATE_FRAME",
    "DIMENSION_TOLERANCE",
    "GRAVITY",
    "SEAWATER_DENSITY",
    "UNITS",
    "VOLUME_TOLERANCE",
    "AnalyticForm",
    "PC6Model",
    "PC6Parameters",
    "analytic_form",
    "beam_sweep",
    "build_mesh",
    "field_value",
    "generate",
    "half_breadth",
    "height_parameters",
    "inside_by_winding_number",
    "is_inside",
    "longitudinal_moments",
    "longitudinal_ratio",
    "sampled_field",
    "section_ratio",
    "station_parameters",
    "stem_x",
    "topology",
]
