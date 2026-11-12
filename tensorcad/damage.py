"""Deterministic damage stability: free-flooding compartments by lost buoyancy.

Everything else in the stability stack answers "what does the intact ship do".
:mod:`tensorcad.intact_stability` evaluates the IS Code 2008 part A criteria on a
measured ``GZ`` curve and both it and :mod:`tensorcad.gz` say so on their pages:
*damage stability is not implemented*. This module is the deterministic half of
that gap -- the classic "one compartment open to the sea" problem a damage control
book is worked out of.

The bookkeeping, written down rather than assumed
-------------------------------------------------
Let ``V_ext`` be the volume of hull below the *external* waterline, ``V_c`` the part
of it inside the flooded boxes, and ``mu`` the compartment permeability (the
fraction of the space water can occupy: 1.0 empty hold, 0.6 cargo hold, 0.95
machinery space, 0 for a tank already full). Flooding adds water *inside* the hull,
so the ship's weight is unchanged and buoyancy is lost where the water sits:

* buoyancy balance ``V_ext - mu V_c = Delta / rho``;
* moment balance: the centroid of ``V_ext - mu V_c`` is the ship's centre of
  gravity, because the external buoyancy acts at one centroid, the flooded water's
  weight at another, and their difference carries the ship.

So the problem is: integrate the *net* volume ``V_ext - mu V_c`` and its first
moments below a waterline plane ``z = z0 + tan(trim)(x - x_ref) + tan(heel) y``, and
solve the three equations ``V = Delta/rho``, ``LCB = x_G``, ``TCB = y_G`` for the
three unknowns ``(z0, tan(trim), tan(heel))``. That is
:func:`flooded_equilibrium`. It is not linearised and the heel is an *output*, so
an asymmetric compartment finds its own list.

What is exact and what is discretised
-------------------------------------
A wall-sided box has a rectangular section, so nothing is discretised there and the
closed form is reproduced to round-off (``docs/damage_stability.md`` measures it).
Along the hull ``A(x)`` jumps at the compartment ends, so the integration grid makes
each compartment boundary a panel edge -- no Gauss point ever sits on a break -- and
each panel is integrated with three-point Gauss-Legendre, exact for the
piecewise-smooth integrand of a polygonal mesh. The remaining discretisation is the
mesh itself.

Metacentric height after damage
-------------------------------
``GM = KB + I_T/(Delta/rho) - KG`` still holds, with ``KB`` the centroid of the net
buoyancy and ``I_T`` the second moment of the *intact* waterplane: the external
waterplane minus ``mu`` times the flooded part. That subtraction is the
lost-buoyancy method's whole answer to the free-surface question -- the water in a
compartment open to the sea keeps its surface parallel to the sea, so there is no
free-surface correction to add, and adding one would count it twice. A compartment
that is *not* open to the sea (an internally ruptured tank) is the added-weight
problem with a free surface and is out of scope. Given ``KG``, :func:`damage_gm`
evaluates it; a negative ``GM`` is reported as a negative number, not clamped.

The righting arm of the damaged ship is here too (:func:`flooded_righting_arm`,
:func:`flooded_righting_arm_curve`): the heel is prescribed instead of solved, so the
three balances become two, and ``GZ`` follows from the net buoyancy centroid. See
``docs/damage_gz.md``.

Still not implemented, and still said so: probabilistic subdivision (SOLAS part B-1,
index A vs R), progressive flooding through internal openings, cross-flooding and the
intermediate stages of flooding.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Iterable, Sequence

from .bonjean import section_loops
from .geometry import Mesh

__all__ = [
    "SCHEMA_FLOODED",
    "DamageError",
    "Compartment",
    "DamageWaterline",
    "FloodedArm",
    "FloodedArmCurve",
    "FloodedState",
    "buoyancy_curve",
    "damage_gm",
    "flooded_equilibrium",
    "flooded_righting_arm",
    "flooded_righting_arm_curve",
    "intact_equilibrium",
]

SCHEMA_FLOODED = "tensorcad.damage.flooded_state/1"

#: Equilibrium residuals: volume to 1e-12 of the displacement, the two centroids to
#: 1e-10 of the hull size. Tighter than any measurement this package makes, loose
#: enough that Newton cannot stall on round-off.
VOLUME_TOLERANCE = 1e-12
MOMENT_TOLERANCE = 1e-10

#: Max Newton iterations before the case is reported as not converged.
MAX_ITERATIONS = 80

#: Three-point Gauss-Legendre per panel: exact for degree 5, which covers ``A(x)``
#: and ``x A(x)`` on a mesh panel.
_GAUSS_X = (-math.sqrt(0.6), 0.0, math.sqrt(0.6))
_GAUSS_W = (5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0)


class DamageError(ValueError):
    """A damage case that cannot be posed, or an equilibrium that was not found."""


def _round(value: float, digits: int = 9) -> float | None:
    if not math.isfinite(value):
        return None
    return round(value, digits)


def _finite(value: Any, name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise DamageError(f"{name} must be a number, got {value!r}") from None
    if not math.isfinite(number):
        raise DamageError(f"{name} must be finite, got {value!r}")
    return number


@dataclass(frozen=True)
class Compartment:
    """A box-shaped flooded space, in the mesh's own frame.

    ``x0``/``x1`` are mandatory; ``y0``/``y1`` default to the hull's full breadth,
    ``z0`` to the keel and ``z1`` to above the deck, so the shortest case -- a
    full-depth, full-breadth slice -- is ``Compartment(x0, x1)``.

    ``permeability`` is ``mu``, the fraction of the space the sea can fill. It is an
    input and is never guessed; :class:`FloodedState` reports
    ``lost_buoyancy_volume = mu * V_c`` so the number used is on the record.
    """

    x0: float
    x1: float
    y0: float | None = None
    y1: float | None = None
    z0: float | None = None
    z1: float | None = None
    permeability: float = 1.0
    name: str = "compartment"

    def __post_init__(self) -> None:
        x0 = _finite(self.x0, "x0")
        x1 = _finite(self.x1, "x1")
        if x1 <= x0:
            raise DamageError(f"compartment {self.name!r}: x1={x1!r} must exceed x0={x0!r}")
        permeability = _finite(self.permeability, "permeability")
        if not 0.0 <= permeability <= 1.0:
            raise DamageError(
                f"compartment {self.name!r}: permeability must lie in [0, 1], "
                f"got {permeability!r}"
            )
        for name in ("y0", "y1", "z0", "z1"):
            value = getattr(self, name)
            if value is not None:
                _finite(value, name)

    @property
    def length(self) -> float:
        return self.x1 - self.x0

    def resolved(self, mesh: Mesh) -> "Compartment":
        """This compartment with the ``None`` bounds filled in from the hull."""
        box = mesh.bounding_box()
        y0 = box.minimum[1] if self.y0 is None else float(self.y0)
        y1 = box.maximum[1] if self.y1 is None else float(self.y1)
        z0 = box.minimum[2] if self.z0 is None else float(self.z0)
        z1 = (box.maximum[2] + box.sizes[2]) if self.z1 is None else float(self.z1)
        if y1 <= y0:
            raise DamageError(f"compartment {self.name!r}: y1={y1!r} must exceed y0={y0!r}")
        if z1 <= z0:
            raise DamageError(f"compartment {self.name!r}: z1={z1!r} must exceed z0={z0!r}")
        return Compartment(
            x0=self.x0,
            x1=self.x1,
            y0=y0,
            y1=y1,
            z0=z0,
            z1=z1,
            permeability=self.permeability,
            name=self.name,
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "x0_m": _round(self.x0),
            "x1_m": _round(self.x1),
            "y0_m": None if self.y0 is None else _round(self.y0),
            "y1_m": None if self.y1 is None else _round(self.y1),
            "z0_m": None if self.z0 is None else _round(self.z0),
            "z1_m": None if self.z1 is None else _round(self.z1),
            "permeability": _round(self.permeability, 12),
        }


@dataclass(frozen=True)
class DamageWaterline:
    """The plane ``z = base + tan(trim) (x - x_ref) + tan(heel) y``."""

    base_z: float
    trim_tangent: float = 0.0
    heel_tangent: float = 0.0
    reference_x: float = 0.0

    def at(self, x: float, y: float = 0.0) -> float:
        return self.base_z + self.trim_tangent * (x - self.reference_x) + self.heel_tangent * y

    @property
    def trim_angle_deg(self) -> float:
        return math.degrees(math.atan(self.trim_tangent))

    @property
    def heel_angle_deg(self) -> float:
        return math.degrees(math.atan(self.heel_tangent))

    def to_dict(self) -> dict[str, Any]:
        return {
            "base_z_m": _round(self.base_z),
            "trim_tangent": _round(self.trim_tangent, 12),
            "heel_tangent": _round(self.heel_tangent, 12),
            "trim_angle_deg": _round(self.trim_angle_deg, 9),
            "heel_angle_deg": _round(self.heel_angle_deg, 9),
            "reference_x_m": _round(self.reference_x),
        }


# --------------------------------------------------------------------------- #
# polygon algebra: half-plane clip, box clip, area and first moments
# --------------------------------------------------------------------------- #
def _clip_halfplane(
    polygon: Sequence[tuple[float, float]], a: float, b: float, c: float
) -> list[tuple[float, float]]:
    """Sutherland-Hodgman clip of a ``(y, z)`` polygon against ``a y + b z <= c``."""
    count = len(polygon)
    if count < 3:
        return []
    kept: list[tuple[float, float]] = []
    for index in range(count):
        y0, z0 = polygon[index]
        y1, z1 = polygon[(index + 1) % count]
        f0 = a * y0 + b * z0 - c
        f1 = a * y1 + b * z1 - c
        inside0 = f0 <= 0.0
        inside1 = f1 <= 0.0
        if inside0:
            kept.append((y0, z0))
        if inside0 != inside1:
            span = f0 - f1
            t = 0.0 if span == 0.0 else f0 / span
            kept.append((y0 + t * (y1 - y0), z0 + t * (z1 - z0)))
    return kept


def _clip_box(
    polygon: Sequence[tuple[float, float]], y0: float, y1: float, z0: float, z1: float
) -> list[tuple[float, float]]:
    """Intersection of a ``(y, z)`` polygon with the axis-aligned rectangle."""
    current = list(polygon)
    for a, b, c in ((1.0, 0.0, y1), (-1.0, 0.0, -y0), (0.0, 1.0, z1), (0.0, -1.0, -z0)):
        if len(current) < 3:
            return []
        current = _clip_halfplane(current, a, b, c)
    return current


def _area_moments(polygon: Sequence[tuple[float, float]]) -> tuple[float, float, float]:
    """``(area, INT y dA, INT z dA)`` of a ``(y, z)`` polygon, area unsigned.

    The shoelace pairing gives the signed area and both first moments together; the
    traversal sign cancels inside ``moment / signed_area`` (the centroid), so the
    physical moment is the absolute area times that centroid. A plane cut of a closed
    mesh can produce loops traversed either way -- both ways at once for a hull with
    a hole -- and this is what makes that irrelevant.
    """
    count = len(polygon)
    if count < 3:
        return 0.0, 0.0, 0.0
    twice = 0.0
    six_y = 0.0
    six_z = 0.0
    for index in range(count):
        py, pz = polygon[index]
        qy, qz = polygon[(index + 1) % count]
        cross = py * qz - qy * pz
        twice += cross
        six_y += (py + qy) * cross
        six_z += (pz + qz) * cross
    signed = 0.5 * twice
    if signed == 0.0:
        return 0.0, 0.0, 0.0
    area = abs(signed)
    centroid_y = (six_y / 6.0) / signed
    centroid_z = (six_z / 6.0) / signed
    return area, area * centroid_y, area * centroid_z


# --------------------------------------------------------------------------- #
# the sampled hull: sections cut once, integrated with Gauss panels
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class _Integrals:
    external_volume: float
    external_x: float
    external_y: float
    external_z: float
    flooded_volume: float
    flooded_x: float
    flooded_y: float
    flooded_z: float

    @property
    def volume(self) -> float:
        return self.external_volume - self.flooded_volume

    @property
    def lcb(self) -> float:
        return (self.external_x - self.flooded_x) / self.volume

    @property
    def tcb(self) -> float:
        return (self.external_y - self.flooded_y) / self.volume

    @property
    def z_moment(self) -> float:
        """``INT z dV`` of the net buoyancy, unnormalised."""
        return self.external_z - self.flooded_z


@dataclass(frozen=True)
class _Waterplane:
    area: float
    centre_y: float
    transverse_second_moment: float
    flooded_area: float


class _HullSections:
    """Cross-sections of one mesh at fixed stations, cut once and then reused.

    ``section_loops`` depends on the station only, never on the waterline, so the
    expensive part (walking every face) happens once at set-up. Each Newton
    iteration then only clips cached polygons.
    """

    def __init__(self, mesh: Mesh, stations: Sequence[float], weights: Sequence[float]) -> None:
        self.mesh = mesh
        self.stations = tuple(stations)
        self.weights = tuple(weights)
        self.loops = tuple(section_loops(mesh, x) for x in self.stations)
        box = mesh.bounding_box()
        self.keel_z = box.minimum[2]
        self.length = box.sizes[0]
        self.beam = box.sizes[1]
        self.scale = max(abs(value) for value in box.maximum) or 1.0
        self.reference_x = 0.5 * (box.minimum[0] + box.maximum[0])

    def integrate(
        self, waterline: DamageWaterline, compartments: Sequence[Compartment]
    ) -> _Integrals:
        """Net buoyancy integrals of the flooded hull below ``waterline``."""
        ext_v = ext_x = ext_y = ext_z = 0.0
        cmp_v = cmp_x = cmp_y = cmp_z = 0.0
        slope = waterline.heel_tangent
        for station, weight, loops in zip(self.stations, self.weights, self.loops, strict=True):
            level = waterline.at(station)
            area_ext = m_y = m_z = 0.0
            area_cmp = cm_y = cm_z = 0.0
            for loop in loops:
                clipped = _clip_halfplane(loop, -slope, 1.0, level)
                if len(clipped) < 3:
                    continue
                area, iy, iz = _area_moments(clipped)
                if area == 0.0:
                    continue
                area_ext += area
                m_y += iy
                m_z += iz
                for compartment in compartments:
                    if not compartment.x0 <= station <= compartment.x1:
                        continue
                    part = _clip_box(
                        clipped,
                        compartment.y0 or 0.0,
                        compartment.y1 or 0.0,
                        compartment.z0 or 0.0,
                        compartment.z1 or 0.0,
                    )
                    if len(part) < 3:
                        continue
                    p_area, p_iy, p_iz = _area_moments(part)
                    if p_area == 0.0:
                        continue
                    mu = compartment.permeability
                    area_cmp += mu * p_area
                    cm_y += mu * p_iy
                    cm_z += mu * p_iz
            ext_v += weight * area_ext
            ext_x += weight * area_ext * station
            ext_y += weight * m_y
            ext_z += weight * m_z
            cmp_v += weight * area_cmp
            cmp_x += weight * area_cmp * station
            cmp_y += weight * cm_y
            cmp_z += weight * cm_z
        return _Integrals(
            external_volume=ext_v,
            external_x=ext_x,
            external_y=ext_y,
            external_z=ext_z,
            flooded_volume=cmp_v,
            flooded_x=cmp_x,
            flooded_y=cmp_y,
            flooded_z=cmp_z,
        )

    def waterplane(
        self, waterline: DamageWaterline, compartments: Sequence[Compartment]
    ) -> _Waterplane:
        """Intact waterplane area and transverse second moment about its own axis.

        The waterplane chords are read off the *clipped* section: an edge of the
        submerged polygon lying on the waterline plane is a chord of the waterplane,
        so the polygon machinery already in hand produces the plane section without
        a second mesh cut in an arbitrary orientation. ``INT y dA`` and ``INT y^2 dA``
        come from the chord width and centre, and ``I_T`` is about the waterplane's
        own centroidal axis (``INT y^2 dA - A ybar^2``), which is the axis ``BM`` uses.
        """
        tolerance = 1e-9 * self.scale
        area = m1 = m2 = 0.0
        cmp_area = cmp_m1 = cmp_m2 = 0.0
        slope = waterline.heel_tangent
        for station, weight, loops in zip(self.stations, self.weights, self.loops, strict=True):
            level = waterline.at(station)
            for loop in loops:
                clipped = _clip_halfplane(loop, -slope, 1.0, level)
                count = len(clipped)
                if count < 3:
                    continue
                for index in range(count):
                    y0, z0 = clipped[index]
                    y1, z1 = clipped[(index + 1) % count]
                    if abs(z0 - (level + slope * y0)) > tolerance:
                        continue
                    if abs(z1 - (level + slope * y1)) > tolerance:
                        continue
                    width = abs(y1 - y0)
                    if width == 0.0:
                        continue
                    centre = 0.5 * (y0 + y1)
                    area += weight * width
                    m1 += weight * width * centre
                    m2 += weight * width * (centre * centre + width * width / 12.0)
                    for compartment in compartments:
                        if not compartment.x0 <= station <= compartment.x1:
                            continue
                        low = max(min(y0, y1), compartment.y0 or 0.0)
                        high = min(max(y0, y1), compartment.y1 or 0.0)
                        if high <= low:
                            continue
                        mu = compartment.permeability
                        for band_low, band_high in _level_band(
                            waterline,
                            station,
                            low,
                            high,
                            compartment.z0 or 0.0,
                            compartment.z1 or 0.0,
                        ):
                            part = band_high - band_low
                            if part <= 0.0:
                                continue
                            part_centre = 0.5 * (band_low + band_high)
                            cmp_area += mu * weight * part
                            cmp_m1 += mu * weight * part * part_centre
                            cmp_m2 += mu * weight * part * (
                                part_centre * part_centre + part * part / 12.0
                            )
        intact_area = area - cmp_area
        if intact_area <= 0.0:
            raise DamageError(
                "the intact waterplane has no area -- the flooding removed all of it "
                f"(external {area!r} m^2, flooded {cmp_area!r} m^2)"
            )
        centre_y = (m1 - cmp_m1) / intact_area
        second = m2 - cmp_m2
        inertia = max(0.0, second - intact_area * centre_y * centre_y)
        return _Waterplane(
            area=intact_area,
            centre_y=centre_y,
            transverse_second_moment=inertia,
            flooded_area=cmp_area,
        )


def _level_band(
    waterline: DamageWaterline,
    station: float,
    y_low: float,
    y_high: float,
    z_low: float,
    z_high: float,
) -> tuple[tuple[float, float], ...]:
    """Sub-interval(s) of ``[y_low, y_high]`` where the waterline is inside the box.

    The waterline is linear in ``y`` at fixed ``x``, so the set where
    ``z_low <= z <= z_high`` is one interval or none.
    """
    base = waterline.at(station)
    slope = waterline.heel_tangent
    if slope == 0.0:
        if z_low <= base <= z_high:
            return ((y_low, y_high),)
        return ()
    lo = (z_low - base) / slope
    hi = (z_high - base) / slope
    if lo > hi:
        lo, hi = hi, lo
    low = max(y_low, lo)
    high = min(y_high, hi)
    if high <= low:
        return ()
    return ((low, high),)


def _panels(
    mesh: Mesh, compartments: Sequence[Compartment], panels: int
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Gauss stations and weights, every compartment boundary on a panel edge.

    A compartment end is a jump in ``A(x)``. Sampling across it would smear the step
    over a whole panel, so the breaks become panel edges and no Gauss point sits on
    one.
    """
    if panels < 1:
        raise DamageError(f"panels must be at least 1, got {panels!r}")
    box = mesh.bounding_box()
    x_min = box.minimum[0]
    x_max = box.maximum[0]
    if x_max <= x_min:
        raise DamageError("the mesh has no length along x")
    breaks = {x_min, x_max}
    for compartment in compartments:
        for edge in (compartment.x0, compartment.x1):
            if x_min < edge < x_max:
                breaks.add(edge)
    edges = sorted(breaks)
    stations: list[float] = []
    weights: list[float] = []
    for index in range(len(edges) - 1):
        a = edges[index]
        b = edges[index + 1]
        if b - a <= 0.0:
            continue
        mid = 0.5 * (a + b)
        half = 0.5 * (b - a)
        for offset, weight in zip(_GAUSS_X, _GAUSS_W, strict=True):
            stations.append(mid + half * offset)
            weights.append(half * weight)
    return tuple(stations), tuple(weights)


def _validated_compartments(
    mesh: Mesh, compartments: Iterable[Compartment]
) -> tuple[Compartment, ...]:
    resolved = tuple(compartment.resolved(mesh) for compartment in compartments)
    for i, first in enumerate(resolved):
        for second in resolved[i + 1 :]:
            if (
                first.x0 < second.x1
                and second.x0 < first.x1
                and (first.y0 or 0.0) < (second.y1 or 0.0)
                and (second.y0 or 0.0) < (first.y1 or 0.0)
                and (first.z0 or 0.0) < (second.z1 or 0.0)
                and (second.z0 or 0.0) < (first.z1 or 0.0)
            ):
                raise DamageError(
                    f"compartments {first.name!r} and {second.name!r} overlap: their "
                    "flooded volumes would be counted twice"
                )
    return resolved


@dataclass(frozen=True)
class FloodedState:
    """The damaged ship's equilibrium, with every number its derivation used."""

    waterline: DamageWaterline
    displacement_volume: float
    volume_net: float
    lcb: float
    tcb: float
    kb: float
    keel_z: float
    draft_at_reference: float
    draft_min: float
    draft_max: float
    trim_angle_deg: float
    heel_angle_deg: float
    lost_buoyancy_volume: float
    waterplane_area: float
    flooded_waterplane_area: float
    transverse_second_moment: float
    bm_transverse: float
    km_transverse: float
    gm_transverse: float | None
    kg: float | None
    intact_base_z: float | None
    sinkage: float | None
    compartments: tuple[Compartment, ...]
    panels: int
    stations: int
    volume_residual: float
    lcg_residual_m: float
    tcb_residual_m: float
    iterations: int
    converged: bool
    schema: str = SCHEMA_FLOODED

    @property
    def damaged(self) -> bool:
        return bool(self.compartments)

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": self.schema,
            "waterline": self.waterline.to_dict(),
            "displacement_volume_m3": _round(self.displacement_volume),
            "volume_net_m3": _round(self.volume_net),
            "lcb_m": _round(self.lcb),
            "tcb_m": _round(self.tcb),
            "kb_m": _round(self.kb),
            "keel_z_m": _round(self.keel_z),
            "draft_at_reference_m": _round(self.draft_at_reference),
            "draft_min_m": _round(self.draft_min),
            "draft_max_m": _round(self.draft_max),
            "trim_angle_deg": _round(self.trim_angle_deg, 9),
            "heel_angle_deg": _round(self.heel_angle_deg, 9),
            "lost_buoyancy_volume_m3": _round(self.lost_buoyancy_volume),
            "waterplane_area_m2": _round(self.waterplane_area),
            "flooded_waterplane_area_m2": _round(self.flooded_waterplane_area),
            "transverse_second_moment_m4": _round(self.transverse_second_moment),
            "bm_transverse_m": _round(self.bm_transverse),
            "km_transverse_m": _round(self.km_transverse),
            "gm_transverse_m": None if self.gm_transverse is None else _round(self.gm_transverse),
            "kg_m": None if self.kg is None else _round(self.kg),
            "intact_base_z_m": None if self.intact_base_z is None else _round(self.intact_base_z),
            "sinkage_m": None if self.sinkage is None else _round(self.sinkage),
            "compartments": [compartment.to_dict() for compartment in self.compartments],
            "panels": self.panels,
            "stations": self.stations,
            "volume_residual": _round(self.volume_residual, 15),
            "lcg_residual_m": _round(self.lcg_residual_m, 12),
            "tcb_residual_m": _round(self.tcb_residual_m, 12),
            "iterations": self.iterations,
            "converged": self.converged,
        }


@dataclass(frozen=True)
class BuoyancyCurve:
    """``A(x)`` of the flooded hull: what the damage took out of the buoyancy curve."""

    stations: tuple[float, ...]
    external_areas: tuple[float, ...]
    net_areas: tuple[float, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "stations_m": [_round(value) for value in self.stations],
            "external_areas_m2": [_round(value) for value in self.external_areas],
            "net_areas_m2": [_round(value) for value in self.net_areas],
        }


def _solve3(matrix: list[list[float]], rhs: list[float]) -> list[float]:
    """Gaussian elimination with partial pivoting on a 3x3 system."""
    rows = [list(matrix[i]) + [rhs[i]] for i in range(3)]
    for column in range(3):
        pivot = max(range(column, 3), key=lambda row: abs(rows[row][column]))
        if abs(rows[pivot][column]) < 1e-300:
            raise DamageError(
                "the equilibrium Jacobian is singular: the damage leaves the waterplane "
                "unswept in a direction the solver cannot see"
            )
        rows[column], rows[pivot] = rows[pivot], rows[column]
        for row in range(column + 1, 3):
            factor = rows[row][column] / rows[column][column]
            if factor == 0.0:
                continue
            for entry in range(column, 4):
                rows[row][entry] -= factor * rows[column][entry]
    solution = [0.0, 0.0, 0.0]
    for row in (2, 1, 0):
        total = rows[row][3]
        for column in range(row + 1, 3):
            total -= rows[row][column] * solution[column]
        solution[row] = total / rows[row][row]
    return solution


def _residuals(
    sections: _HullSections,
    waterline: DamageWaterline,
    compartments: Sequence[Compartment],
    target_volume: float,
    lcg: float,
    tcb: float | None,
) -> tuple[list[float], _Integrals]:
    integrals = sections.integrate(waterline, compartments)
    if integrals.volume <= 0.0:
        raise DamageError(
            "the net buoyant volume vanished: the flooded compartments removed more "
            "buoyancy than the hull had at this waterline"
        )
    residuals = [
        (integrals.volume - target_volume) / target_volume,
        (integrals.lcb - lcg) / sections.length,
    ]
    if tcb is not None:
        residuals.append((integrals.tcb - tcb) / sections.beam)
    return residuals, integrals


def _norm(residuals: Sequence[float]) -> float:
    return max(abs(value) for value in residuals)


def _converged(residuals: Sequence[float], volume: float) -> bool:
    if abs(residuals[0]) > VOLUME_TOLERANCE:
        return False
    return all(abs(value) <= MOMENT_TOLERANCE for value in residuals[1:])


def _newton(
    sections: _HullSections,
    compartments: Sequence[Compartment],
    target_volume: float,
    lcg: float,
    tcb: float | None,
    initial: DamageWaterline,
) -> tuple[DamageWaterline, _Integrals, list[float], int, bool]:
    """Damped Newton on ``(z0, tan(trim), tan(heel))`` against the three balances."""
    variables = [initial.base_z, initial.trim_tangent, initial.heel_tangent]
    steps = [1e-6 * sections.scale, 1e-7, 1e-7]
    residuals, integrals = _residuals(
        sections, initial, compartments, target_volume, lcg, tcb
    )
    best = (_norm(residuals), list(variables), residuals, integrals)
    iterations = 0
    for iterations in range(1, MAX_ITERATIONS + 1):
        if _converged(residuals, target_volume):
            return (
                DamageWaterline(
                    base_z=variables[0], trim_tangent=variables[1], heel_tangent=variables[2],
                    reference_x=sections.reference_x,
                ),
                integrals,
                residuals,
                iterations - 1,
                True,
            )
        jacobian: list[list[float]] = [[0.0, 0.0, 0.0] for _ in range(3)]
        for column in range(3):
            trial = list(variables)
            trial[column] += steps[column]
            probe = DamageWaterline(
                base_z=trial[0], trim_tangent=trial[1], heel_tangent=trial[2],
                reference_x=sections.reference_x,
            )
            try:
                probe_residuals, _ = _residuals(
                    sections, probe, compartments, target_volume, lcg, tcb
                )
            except DamageError:
                # a finite-difference probe that leaves the valid region: shorten it
                steps[column] *= 0.5
                probe_residuals = residuals
            for row in range(len(residuals)):
                jacobian[row][column] = (probe_residuals[row] - residuals[row]) / steps[column]
        delta = _solve3(jacobian, [-value for value in residuals])
        # backtrack so a Newton step can never make the residuals worse
        scale = 1.0
        improvements = False
        for _ in range(20):
            trial = [variables[i] + scale * delta[i] for i in range(3)]
            probe = DamageWaterline(
                base_z=trial[0], trim_tangent=trial[1], heel_tangent=trial[2],
                reference_x=sections.reference_x,
            )
            try:
                probe_residuals, probe_integrals = _residuals(
                    sections, probe, compartments, target_volume, lcg, tcb
                )
            except DamageError:
                scale *= 0.5
                continue
            if _norm(probe_residuals) < _norm(residuals):
                variables = trial
                residuals = probe_residuals
                integrals = probe_integrals
                improvements = True
                break
            scale *= 0.5
        if not improvements:
            break
        if _norm(residuals) < best[0]:
            best = (_norm(residuals), list(variables), residuals, integrals)
    score, variables, residuals, integrals = best
    waterline = DamageWaterline(
        base_z=variables[0], trim_tangent=variables[1], heel_tangent=variables[2],
        reference_x=sections.reference_x,
    )
    return waterline, integrals, residuals, iterations, _converged(residuals, target_volume)


def _level_waterline(
    sections: _HullSections, target_volume: float
) -> tuple[float, _Integrals]:
    """Bisect the level waterline ``z0`` that displaces ``target_volume`` of hull."""
    box = sections.mesh.bounding_box()
    low = box.minimum[2]
    high = box.maximum[2]

    def volume(z0: float) -> float:
        return sections.integrate(
            DamageWaterline(base_z=z0, reference_x=sections.reference_x), ()
        ).external_volume

    low_volume = volume(low + 1e-9 * (high - low))
    high_volume = volume(high)
    if high_volume < target_volume:
        raise DamageError(
            f"displacement_volume={target_volume!r} m^3 exceeds the hull's own volume "
            f"{high_volume!r} m^3 -- the ship would be a submarine"
        )
    if low_volume > target_volume:
        raise DamageError(
            f"displacement_volume={target_volume!r} m^3 is below the volume at the keel"
        )
    for _ in range(200):
        mid = 0.5 * (low + high)
        if volume(mid) < target_volume:
            low = mid
        else:
            high = mid
        if high - low <= 1e-15 * (abs(high) + abs(low) + 1.0):
            break
    base = 0.5 * (low + high)
    waterline = DamageWaterline(base_z=base, reference_x=sections.reference_x)
    integrals = sections.integrate(waterline, ())
    return base, integrals


def _build_state(
    sections: _HullSections,
    waterline: DamageWaterline,
    integrals: _Integrals,
    residuals: Sequence[float],
    iterations: int,
    converged: bool,
    compartments: Sequence[Compartment],
    panels: int,
    target_volume: float,
    kg: float | None,
    intact_base: float | None,
) -> FloodedState:
    waterplane = sections.waterplane(waterline, compartments)
    volume = integrals.volume
    kb = integrals.z_moment / volume - sections.keel_z
    bm = waterplane.transverse_second_moment / volume
    box = sections.mesh.bounding_box()
    drafts = [
        waterline.at(box.minimum[0], box.minimum[1]) - sections.keel_z,
        waterline.at(box.minimum[0], box.maximum[1]) - sections.keel_z,
        waterline.at(box.maximum[0], box.minimum[1]) - sections.keel_z,
        waterline.at(box.maximum[0], box.maximum[1]) - sections.keel_z,
    ]
    gm = None if kg is None else kb + bm - kg
    return FloodedState(
        waterline=waterline,
        displacement_volume=target_volume,
        volume_net=volume,
        lcb=integrals.lcb,
        tcb=integrals.tcb,
        kb=kb,
        keel_z=sections.keel_z,
        draft_at_reference=waterline.base_z - sections.keel_z,
        draft_min=min(drafts),
        draft_max=max(drafts),
        trim_angle_deg=waterline.trim_angle_deg,
        heel_angle_deg=waterline.heel_angle_deg,
        lost_buoyancy_volume=integrals.flooded_volume,
        waterplane_area=waterplane.area,
        flooded_waterplane_area=waterplane.flooded_area,
        transverse_second_moment=waterplane.transverse_second_moment,
        bm_transverse=bm,
        km_transverse=kb + bm,
        gm_transverse=gm,
        kg=kg,
        intact_base_z=intact_base,
        sinkage=None if intact_base is None else waterline.base_z - intact_base,
        compartments=tuple(compartments),
        panels=panels,
        stations=len(sections.stations),
        volume_residual=residuals[0] * target_volume,
        lcg_residual_m=residuals[1] * sections.length,
        tcb_residual_m=(residuals[2] * sections.beam) if len(residuals) > 2 else 0.0,
        iterations=iterations,
        converged=converged,
    )


def intact_equilibrium(
    mesh: Mesh,
    displacement_volume: float,
    *,
    panels: int = 8,
) -> FloodedState:
    """The undamaged ship on an even keel at ``displacement_volume``.

    The level waterline is a monotone bisection on ``V(z0)``, not a Newton solve, so
    the intact answer is independent of the damaged solver and can be used to check it
    (a zero-permeability compartment has to reproduce this state exactly).
    """
    target = _finite(displacement_volume, "displacement_volume")
    if target <= 0.0:
        raise DamageError(f"displacement_volume must be positive, got {target!r}")
    stations, weights = _panels(mesh, (), panels)
    sections = _HullSections(mesh, stations, weights)
    base, integrals = _level_waterline(sections, target)
    waterline = DamageWaterline(base_z=base, reference_x=sections.reference_x)
    residuals = [(integrals.volume - target) / target, 0.0, 0.0]
    return _build_state(
        sections,
        waterline,
        integrals,
        residuals,
        iterations=0,
        converged=abs(integrals.volume - target) <= VOLUME_TOLERANCE * target,
        compartments=(),
        panels=panels,
        target_volume=target,
        kg=None,
        intact_base=None,
    )


def damage_gm(state: FloodedState, kg: float) -> float:
    """``GM = KB + I_T/(Delta/rho) - KG`` of the damaged state, no free-surface term.

    The intact waterplane in ``state`` is already the *net* one, so the loss of
    waterplane caused by flooding is inside ``I_T`` and ``BM``; adding a free-surface
    correction on top would subtract the same effect twice.
    """
    weight = _finite(kg, "kg")
    return state.kb + state.transverse_second_moment / state.volume_net - weight


def flooded_equilibrium(
    mesh: Mesh,
    compartments: Sequence[Compartment] | Compartment,
    displacement_volume: float,
    *,
    lcg: float | None = None,
    tcb: float | None = None,
    kg: float | None = None,
    panels: int = 12,
    raise_on_failure: bool = True,
) -> FloodedState:
    """Solve the free-flooding damage case for sinkage, trim and heel at once.

    ``displacement_volume`` is ``Delta/rho`` of the **undamaged** ship: flooding adds
    water inside the hull, not weight to it. ``lcg``/``tcb`` default to the intact
    ship's own centres of buoyancy at that displacement, i.e. "the ship is intact
    trimmed, then opened up"; pass them explicitly to place the weight elsewhere.
    ``kg`` is only needed for ``GM``; without it the ``GM`` fields are ``None`` rather
    than guessed.

    Raises :class:`DamageError` when the case cannot be posed (a compartment off the
    hull, overlapping compartments, flooding that removes the whole waterplane) and,
    unless ``raise_on_failure=False``, when Newton does not converge -- a
    ``converged=False`` state is never quietly returned as an answer.
    """
    target = _finite(displacement_volume, "displacement_volume")
    if target <= 0.0:
        raise DamageError(f"displacement_volume must be positive, got {target!r}")
    if isinstance(compartments, Compartment):
        compartments = (compartments,)
    resolved = _validated_compartments(mesh, compartments)
    box = mesh.bounding_box()
    for compartment in resolved:
        if compartment.x1 <= box.minimum[0] or compartment.x0 >= box.maximum[0]:
            raise DamageError(
                f"compartment {compartment.name!r} lies outside the hull along x: "
                f"[{compartment.x0!r}, {compartment.x1!r}] against "
                f"[{box.minimum[0]!r}, {box.maximum[0]!r}]"
            )
    stations, weights = _panels(mesh, resolved, panels)
    sections = _HullSections(mesh, stations, weights)

    intact_base, intact_integrals = _level_waterline(sections, target)
    x_g = intact_integrals.lcb if lcg is None else _finite(lcg, "lcg")
    y_g = intact_integrals.tcb if tcb is None else _finite(tcb, "tcb")
    if any(compartment.permeability > 0.0 for compartment in resolved):
        start = DamageWaterline(base_z=intact_base, reference_x=sections.reference_x)
        waterline, integrals, residuals, iterations, converged = _newton(
            sections, resolved, target, x_g, y_g, start
        )
    else:
        # no flooding: the intact state already satisfies every balance, and it does so
        # exactly, so the solver must not perturb it
        waterline = DamageWaterline(base_z=intact_base, reference_x=sections.reference_x)
        integrals = intact_integrals
        residuals = [
            (integrals.volume - target) / target,
            (integrals.lcb - x_g) / sections.length,
            (integrals.tcb - y_g) / sections.beam,
        ]
        iterations = 0
        converged = _converged(residuals, target)
    state = _build_state(
        sections,
        waterline,
        integrals,
        residuals,
        iterations=iterations,
        converged=converged,
        compartments=resolved,
        panels=panels,
        target_volume=target,
        kg=None if kg is None else _finite(kg, "kg"),
        intact_base=intact_base,
    )
    if raise_on_failure and not state.converged:
        raise DamageError(
            "the damaged equilibrium did not converge: "
            f"volume residual {state.volume_residual!r} m^3, "
            f"LCG residual {state.lcg_residual_m!r} m, "
            f"TCB residual {state.tcb_residual_m!r} m after {state.iterations} iterations"
        )
    return state


def buoyancy_curve(
    mesh: Mesh,
    compartments: Sequence[Compartment] | Compartment = (),
    *,
    waterline: DamageWaterline | None = None,
    displacement_volume: float | None = None,
    stations: int = 41,
) -> BuoyancyCurve:
    """``A(x)`` of the external hull and of the net (intact) buoyancy at a waterline.

    With neither ``waterline`` nor ``displacement_volume`` the intact level waterline
    of the hull's own volume is used, so the curve is defined for any mesh. This is
    the figure a damage control book draws: buoyancy curve before and after, with the
    compartment's worth of area cut out of it.
    """
    if isinstance(compartments, Compartment):
        compartments = (compartments,)
    resolved = _validated_compartments(mesh, compartments)
    if stations < 3:
        raise DamageError(f"stations must be at least 3, got {stations!r}")
    box = mesh.bounding_box()
    grid = tuple(
        box.minimum[0] + (box.maximum[0] - box.minimum[0]) * index / (stations - 1)
        for index in range(stations)
    )
    if waterline is None:
        if displacement_volume is None:
            raise DamageError(
                "buoyancy_curve needs either a waterline or the displacement volume that "
                "places one: a hull with freeboard has no default waterline to draw at"
            )
        target = _finite(displacement_volume, "displacement_volume")
        if target <= 0.0:
            raise DamageError(f"displacement_volume must be positive, got {target!r}")
        solver_stations, solver_weights = _panels(mesh, resolved, 8)
        solver = _HullSections(mesh, solver_stations, solver_weights)
        _base, _integrals = _level_waterline(solver, target)
        waterline = DamageWaterline(base_z=_base, reference_x=solver.reference_x)
    external: list[float] = []
    net: list[float] = []
    for station in grid:
        loops = section_loops(mesh, station)
        level = waterline.at(station)
        slope = waterline.heel_tangent
        area_ext = 0.0
        area_cmp = 0.0
        for loop in loops:
            clipped = _clip_halfplane(loop, -slope, 1.0, level)
            if len(clipped) < 3:
                continue
            area, _iy, _iz = _area_moments(clipped)
            area_ext += area
            for compartment in resolved:
                if not compartment.x0 <= station <= compartment.x1:
                    continue
                part = _clip_box(
                    clipped,
                    compartment.y0 or 0.0,
                    compartment.y1 or 0.0,
                    compartment.z0 or 0.0,
                    compartment.z1 or 0.0,
                )
                if len(part) < 3:
                    continue
                p_area, _piy, _piz = _area_moments(part)
                area_cmp += compartment.permeability * p_area
        external.append(area_ext)
        net.append(area_ext - area_cmp)
    return BuoyancyCurve(
        stations=grid, external_areas=tuple(external), net_areas=tuple(net)
    )
# --------------------------------------------------------------------------- #
# the righting arm of the flooded ship: the waterline is free, the heel is not
# --------------------------------------------------------------------------- #
SCHEMA_ARM = "tensorcad.damage.flooded_arm/1"

#: Heel ladder limit for the damaged arm. Not a geometric limit of the integrator --
#: that would take any angle -- but the point past which the lost-buoyancy model
#: stops describing the case: a compartment open to the sea loses its water back
#: through the opening once the waterline in the compartment rises above the damage,
#: and on the hulls this package generates that has happened well before 60 degrees.
MAX_ARM_HEEL_DEG = 60.0


@dataclass(frozen=True)
class FloodedArm:
    """One point of the *damaged* ship's righting-arm curve.

    ``righting_arm`` is ``(y_B - y_G) cos(phi) + (KB - KG) sin(phi)``, both points
    taken in the world frame of the heeled body and every vertical coordinate measured
    above the keel. That is the same transverse convention :func:`tensorcad.gz.righting_arm`
    uses for the intact ship (``GZ = y_B - y_G`` there), written out in the two body-frame
    coordinates this module integrates so the arm can be formed without rotating the mesh.

    ``volume``, ``lcb``, ``tcb``, ``kb`` are the **net** (external minus flooded)
    quantities: the ship floats on the difference, and the arm is the offset of that
    difference's centroid from ``G``.
    """

    heel_deg: float
    righting_arm: float
    kg: float
    lcg: float
    tcb_gravity: float
    keel_z: float
    volume: float
    lcb: float
    tcb: float
    kb: float
    lost_buoyancy_volume: float
    draft_at_reference: float
    draft_min: float
    draft_max: float
    waterline: DamageWaterline
    volume_residual: float
    lcb_residual_m: float
    iterations: int
    converged: bool
    schema: str = SCHEMA_ARM

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": self.schema,
            "heel_deg": _round(self.heel_deg, 9),
            "righting_arm_m": _round(self.righting_arm, 12),
            "kg_m": _round(self.kg),
            "lcg_m": _round(self.lcg),
            "tcb_gravity_m": _round(self.tcb_gravity, 12),
            "keel_z_m": _round(self.keel_z, 12),
            "volume_m3": _round(self.volume),
            "lcb_m": _round(self.lcb),
            "tcb_m": _round(self.tcb),
            "kb_m": _round(self.kb),
            "lost_buoyancy_volume_m3": _round(self.lost_buoyancy_volume),
            "draft_at_reference_m": _round(self.draft_at_reference),
            "draft_min_m": _round(self.draft_min),
            "draft_max_m": _round(self.draft_max),
            "waterline": self.waterline.to_dict(),
            "volume_residual_m3": _round(self.volume_residual, 15),
            "lcb_residual_m": _round(self.lcb_residual_m, 12),
            "iterations": self.iterations,
            "converged": self.converged,
        }


@dataclass(frozen=True)
class FloodedArmCurve:
    """The damaged ship's ``GZ(phi)`` with the rule quantities read off the curve.

    ``positive_range_deg``, ``vanishing_angle_deg``, ``max_arm``, ``heel_at_max_arm_deg``
    and ``area_to_max_arm`` are measurements of the curve this module produced -- they
    are *not* a comparison against a statutory criterion, which is a separate question
    with its own thresholds and its own loading conditions.
    """

    points: tuple[FloodedArm, ...]
    kg: float
    lcg: float
    displacement_volume: float
    basis: str = "constant displacement, sinkage and trim solved at every heel"

    @property
    def peak(self) -> FloodedArm:
        return max(self.points, key=lambda point: point.righting_arm)

    @property
    def max_arm(self) -> float:
        return self.peak.righting_arm

    @property
    def heel_at_max_arm_deg(self) -> float:
        return self.peak.heel_deg

    @property
    def vanishing_angle_deg(self) -> float | None:
        """First heel past the arm's peak where ``GZ`` reaches zero, or ``None``.

        ``None`` means the ladder ended while the arm was still positive: the range is
        *at least* what the ladder covered, and saying so is the honest answer.
        """
        previous = None
        for point in self.points:
            if previous is not None and previous.heel_deg >= self.heel_at_max_arm_deg:
                if point.righting_arm <= 0.0 < previous.righting_arm:
                    span = previous.righting_arm - point.righting_arm
                    fraction = previous.righting_arm / span if span else 0.0
                    return previous.heel_deg + fraction * (point.heel_deg - previous.heel_deg)
            previous = point
        return None

    @property
    def positive_range_deg(self) -> float:
        vanishing = self.vanishing_angle_deg
        if vanishing is None:
            return self.points[-1].heel_deg
        return vanishing

    @property
    def area_to_max_arm(self) -> float:
        """``INT GZ dphi`` from upright to the peak, in m rad (trapezoid on the ladder)."""
        total = 0.0
        for first, second in zip(self.points, self.points[1:], strict=False):
            if first.heel_deg >= self.heel_at_max_arm_deg:
                break
            phi_a = math.radians(first.heel_deg)
            phi_b = math.radians(second.heel_deg)
            arm_a = first.righting_arm
            arm_b = second.righting_arm
            if second.heel_deg > self.heel_at_max_arm_deg:
                fraction = (self.heel_at_max_arm_deg - first.heel_deg) / (
                    second.heel_deg - first.heel_deg
                )
                phi_b = math.radians(first.heel_deg + fraction * (second.heel_deg - first.heel_deg))
                arm_b = arm_a + fraction * (arm_b - arm_a)
            total += 0.5 * (arm_a + arm_b) * (phi_b - phi_a)
        return total

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA_ARM,
            "kg_m": _round(self.kg),
            "lcg_m": _round(self.lcg),
            "displacement_volume_m3": _round(self.displacement_volume),
            "righting_arms": [point.to_dict() for point in self.points],
            "max_arm_m": _round(self.max_arm, 12),
            "heel_at_max_arm_deg": _round(self.heel_at_max_arm_deg, 9),
            "vanishing_angle_deg": (
                None
                if self.vanishing_angle_deg is None
                else _round(self.vanishing_angle_deg, 9)
            ),
            "positive_range_deg": _round(self.positive_range_deg, 9),
            "area_to_max_arm_m_rad": _round(self.area_to_max_arm, 12),
            "basis": self.basis,
        }


def _solve2(matrix: list[list[float]], rhs: list[float]) -> list[float]:
    """Gaussian elimination with partial pivoting on a 2x2 system."""
    a, b = matrix[0]
    c, d = matrix[1]
    determinant = a * d - b * c
    if abs(determinant) < 1e-300:
        raise DamageError(
            "the heel-fixed equilibrium Jacobian is singular: sinkage and trim cannot "
            "both be recovered at this heel"
        )
    return [
        (rhs[0] * d - b * rhs[1]) / determinant,
        (a * rhs[1] - rhs[0] * c) / determinant,
    ]


def _newton_heel_fixed(
    sections: _HullSections,
    compartments: Sequence[Compartment],
    target_volume: float,
    lcg: float,
    heel_tangent: float,
    initial: DamageWaterline,
) -> tuple[DamageWaterline, _Integrals, list[float], int, bool]:
    """Damped Newton on ``(z0, tan(trim))`` with the heel held at ``heel_tangent``.

    The transverse balance is dropped, not solved: at a prescribed heel the ship is
    held there by an external moment and only volume and longitudinal position are
    free, so the third unknown disappears instead of being pinned by an equation. That
    is why this is not :func:`_newton` with a column deleted -- the residual *count*
    changes and with it the system that has to be eliminated.
    """
    variables = [initial.base_z, initial.trim_tangent]
    steps = [1e-6 * sections.scale, 1e-8]

    def waterline_of(values: Sequence[float]) -> DamageWaterline:
        return DamageWaterline(
            base_z=values[0],
            trim_tangent=values[1],
            heel_tangent=heel_tangent,
            reference_x=sections.reference_x,
        )

    residuals, integrals = _residuals(
        sections, waterline_of(variables), compartments, target_volume, lcg, None
    )
    best = (_norm(residuals), list(variables), residuals, integrals)
    iterations = 0
    for iterations in range(1, MAX_ITERATIONS + 1):
        if _converged(residuals, target_volume):
            return (
                waterline_of(variables), integrals, residuals, iterations - 1, True
            )
        jacobian: list[list[float]] = [[0.0, 0.0] for _ in range(2)]
        for column in range(2):
            trial = list(variables)
            trial[column] += steps[column]
            try:
                probe_residuals, _ = _residuals(
                    sections, waterline_of(trial), compartments, target_volume, lcg, None
                )
            except DamageError:
                steps[column] *= 0.5
                probe_residuals = residuals
            for row in range(2):
                jacobian[row][column] = (
                    probe_residuals[row] - residuals[row]
                ) / steps[column]
        delta = _solve2(jacobian, [-value for value in residuals])
        improved = False
        scale = 1.0
        for _ in range(20):
            trial = [variables[index] + scale * delta[index] for index in range(2)]
            try:
                probe_residuals, probe_integrals = _residuals(
                    sections, waterline_of(trial), compartments, target_volume, lcg, None
                )
            except DamageError:
                scale *= 0.5
                continue
            if _norm(probe_residuals) < _norm(residuals):
                variables = trial
                residuals = probe_residuals
                integrals = probe_integrals
                improved = True
                break
            scale *= 0.5
        if not improved:
            break
        if _norm(residuals) < best[0]:
            best = (_norm(residuals), list(variables), residuals, integrals)
    score, variables, residuals, integrals = best
    return (
        waterline_of(variables),
        integrals,
        residuals,
        iterations,
        _converged(residuals, target_volume),
    )


def flooded_righting_arm(
    mesh: Mesh,
    compartments: Sequence[Compartment] | Compartment,
    displacement_volume: float,
    *,
    heel_deg: float,
    kg: float,
    lcg: float | None = None,
    tcb: float | None = None,
    panels: int = 12,
    start: DamageWaterline | None = None,
    raise_on_failure: bool = True,
) -> FloodedArm:
    """``GZ`` of the free-flooding damaged ship at a *prescribed* heel.

    One point of :func:`flooded_righting_arm_curve`, and the function to call when a
    single angle is what is wanted. ``kg`` is above the keel, as everywhere else in the
    stability stack, and is required: an arm without a centre of gravity is not an arm.
    ``lcg``/``tcb`` default to the intact ship's own centres at this displacement, i.e.
    the weight stays where the intact ship carried it.

    ``start`` warm-starts the sinkage/trim solve from a neighbouring heel; a curve
    should pass the previous point's waterline, which keeps the Newton history short
    and the answers continuous in ``heel_deg``.
    """
    target = _finite(displacement_volume, "displacement_volume")
    if target <= 0.0:
        raise DamageError(f"displacement_volume must be positive, got {target!r}")
    heel = _finite(heel_deg, "heel_deg")
    if abs(heel) > MAX_ARM_HEEL_DEG:
        raise DamageError(
            f"heel_deg {heel!r} is outside [-{MAX_ARM_HEEL_DEG:g}, {MAX_ARM_HEEL_DEG:g}] "
            "degrees: past that the lost-buoyancy bookkeeping no longer describes a "
            "compartment open to the sea"
        )
    weight = _finite(kg, "kg")
    if isinstance(compartments, Compartment):
        compartments = (compartments,)
    resolved = _validated_compartments(mesh, compartments)
    stations, weights = _panels(mesh, resolved, panels)
    sections = _HullSections(mesh, stations, weights)
    intact_base, intact_integrals = _level_waterline(sections, target)
    x_g = intact_integrals.lcb if lcg is None else _finite(lcg, "lcg")
    y_g = intact_integrals.tcb if tcb is None else _finite(tcb, "tcb")
    initial = (
        DamageWaterline(base_z=intact_base, reference_x=sections.reference_x)
        if start is None
        else start
    )
    tangent = math.tan(math.radians(heel))
    waterline, integrals, residuals, iterations, converged = _newton_heel_fixed(
        sections, resolved, target, x_g, tangent, initial
    )
    volume = integrals.volume
    if volume <= 0.0:
        raise DamageError(
            "the net buoyant volume vanished at this heel: the flooded compartments "
            "removed more buoyancy than the hull had"
        )
    kb = integrals.z_moment / volume - sections.keel_z
    phi = math.radians(heel)
    arm = (integrals.tcb - y_g) * math.cos(phi) + (kb - weight) * math.sin(phi)
    box = mesh.bounding_box()
    drafts = [
        waterline.at(box.minimum[0], box.minimum[1]) - sections.keel_z,
        waterline.at(box.minimum[0], box.maximum[1]) - sections.keel_z,
        waterline.at(box.maximum[0], box.minimum[1]) - sections.keel_z,
        waterline.at(box.maximum[0], box.maximum[1]) - sections.keel_z,
    ]
    point = FloodedArm(
        heel_deg=heel,
        righting_arm=arm,
        kg=weight,
        lcg=x_g,
        tcb_gravity=y_g,
        keel_z=sections.keel_z,
        volume=volume,
        lcb=integrals.lcb,
        tcb=integrals.tcb,
        kb=kb,
        lost_buoyancy_volume=integrals.flooded_volume,
        draft_at_reference=waterline.base_z - sections.keel_z,
        draft_min=min(drafts),
        draft_max=max(drafts),
        waterline=waterline,
        volume_residual=residuals[0] * target,
        lcb_residual_m=residuals[1] * sections.length,
        iterations=iterations,
        converged=converged,
    )
    if raise_on_failure and not point.converged:
        raise DamageError(
            f"the flooded equilibrium did not converge at heel {heel!r} deg: volume "
            f"residual {point.volume_residual!r} m^3, LCB residual "
            f"{point.lcb_residual_m!r} m after {point.iterations} iterations"
        )
    return point


def flooded_righting_arm_curve(
    mesh: Mesh,
    compartments: Sequence[Compartment] | Compartment,
    displacement_volume: float,
    *,
    kg: float,
    heels_deg: Sequence[float] | None = None,
    lcg: float | None = None,
    tcb: float | None = None,
    panels: int = 12,
    raise_on_failure: bool = True,
) -> FloodedArmCurve:
    """The damaged ship's ``GZ`` over a heel ladder, each point re-solved.

    ``heels_deg`` defaults to 0..60 in 5-degree steps. Every point re-solves sinkage
    and trim at constant displacement, so the curve is a locus of floating attitudes
    rather than a rotation of the upright one -- which is the whole reason the damaged
    curve differs from the intact ship's.
    """
    if heels_deg is None:
        heels = tuple(float(value) for value in range(0, 61, 5))
    else:
        heels = tuple(sorted({float(value) for value in heels_deg}))
    if not heels:
        raise DamageError("heels_deg is empty: a curve needs at least one heel")
    weight = _finite(kg, "kg")
    points: list[FloodedArm] = []
    start: DamageWaterline | None = None
    for heel in heels:
        point = flooded_righting_arm(
            mesh,
            compartments,
            displacement_volume,
            heel_deg=heel,
            kg=weight,
            lcg=lcg,
            tcb=tcb,
            panels=panels,
            start=start,
            raise_on_failure=raise_on_failure,
        )
        start = point.waterline
        points.append(point)
        if lcg is None:
            lcg = point.lcg
    return FloodedArmCurve(
        points=tuple(points),
        kg=weight,
        lcg=float(lcg),
        displacement_volume=float(displacement_volume),
    )