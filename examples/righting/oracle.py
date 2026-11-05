"""An independent righting-arm integrator, kept as a verification oracle.

This is the *third* code path for ``GZ(phi)`` in the tree, and it is deliberately
not part of the public API. :mod:`tensorcad.gz` is canonical (it rotates the mesh,
solves the waterline with Newton's method and takes volume and first moments from
the divergence theorem); :mod:`tensorcad.righting_arm` is a second, bisection-based
implementation kept for cross-checking. Both of those integrate the *closed
surface*. This one integrates the *heeled sections* instead:

* the mesh is cut with the plane ``x = station`` once per station
  (:func:`tensorcad.bonjean.section_loops`), giving closed ``(y, z)`` polygons;
* each polygon is clipped against the tilted half-plane
  ``s(y, z) = z cos(phi) - y sin(phi) <= offset`` with Sutherland-Hodgman, and its
  area and first moments come from Green's theorem, which is *exact* for the
  polygon the clip produced;
* the station areas and moments are integrated along the hull with composite
  Simpson, which gives the heeled volume and the transverse moments in one pass;
* the offset that keeps the displacement equal to the upright displacement is
  found by bisection, and

      ``GZ = (y_B - y_G) cos(phi) + (z_B - z_G) sin(phi)``

  with ``z`` measured from the keel (the mesh's lowest vertex) and ``z_G = keel + kg``.

The point of a third path is that its error sources are *different* from the
canonical ones: here the residual error is the longitudinal quadrature of the
station areas plus the polygonisation of the sections, and there is no cap, no
Newton iteration and no surface triangulation sum to get wrong. Agreement between
the two therefore catches coding errors that a single implementation cannot, and
``examples/righting/run_example.py`` measures exactly that.

Scope, inherited from the method: sections must be single-valued in ``z`` per
``y`` per station, so no overhangs and no immersed deck below the modelled body;
``kg`` is a load input above the keel, never inferred; free surface, damage and
statutory criteria are all outside. Only the standard library is used.
"""

from __future__ import annotations

import itertools
import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import Any

from tensorcad.bonjean import section_loops
from tensorcad.geometry import Mesh

#: Default heel stations for a curve: the 5-degree table a cross curve is quoted at.
DEFAULT_HEEL_ANGLES_DEG: tuple[float, ...] = tuple(float(angle) for angle in range(0, 91, 5))

#: Longitudinal intervals for the Simpson sum; even, because Simpson needs it.
DEFAULT_STATION_INTERVALS = 64

#: Below this the quadrature of a hull with pointed ends is not worth reporting.
MINIMUM_STATION_INTERVALS = 8

#: Bisection stops when the bracket is this fraction of the model scale.
OFFSET_TOLERANCE = 1e-12

MAX_BISECTION_STEPS = 200


class OracleError(ValueError):
    """Raised when this integrator cannot define a curve for the geometry given."""


@dataclass(frozen=True)
class OraclePoint:
    """One heel angle: the arm, the displacement it floats at, and where ``B`` went."""

    heel_deg: float
    gz: float
    submerged_volume: float
    transverse_centre_buoyancy: float
    vertical_centre_buoyancy: float
    offset: float
    deck_immersion_fraction: float
    converged: bool = True

    @property
    def deck_edge_immersed(self) -> bool:
        """True when any deck station's edge is under the heeled waterplane."""
        return self.deck_immersion_fraction > 0.0

    def to_dict(self) -> dict[str, Any]:
        return {
            "heel_deg": _round(self.heel_deg),
            "gz_m": _round(self.gz),
            "submerged_volume_m3": _round(self.submerged_volume),
            "y_b_m": _round(self.transverse_centre_buoyancy),
            "z_b_m": _round(self.vertical_centre_buoyancy),
            "offset_m": _round(self.offset),
            "deck_immersion_fraction": _round(self.deck_immersion_fraction),
            "converged": self.converged,
        }


@dataclass(frozen=True)
class OracleCurve:
    """A sampled righting-arm curve from the section integration."""

    kg: float
    y_g: float
    keel_z: float
    waterline_z: float
    displacement_volume: float
    points: tuple[OraclePoint, ...]
    station_intervals: int
    heel_angles_deg: tuple[float, ...]
    basis: str = "heeled section areas (Green), Simpson along the hull"

    @property
    def gz(self) -> tuple[float, ...]:
        return tuple(point.gz for point in self.points)

    @property
    def gz_max(self) -> float:
        """Largest sampled arm; ``nan`` when the curve never rights the body."""
        values = [point.gz for point in self.points if point.converged]
        if not values:
            return float("nan")
        best = max(values)
        return best if best > 0.0 else float("nan")

    @property
    def angle_gz_max_deg(self) -> float:
        """Heel of the largest sampled arm, refined by a parabola through its neighbours."""
        converged = [point for point in self.points if point.converged]
        if not converged:
            return float("nan")
        best = max(range(len(converged)), key=lambda index: converged[index].gz)
        if converged[best].gz <= 0.0:
            return float("nan")
        if best == 0 or best == len(converged) - 1:
            return converged[best].heel_deg
        left, middle, right = converged[best - 1], converged[best], converged[best + 1]
        step = middle.heel_deg - left.heel_deg
        if step <= 0.0 or right.heel_deg - middle.heel_deg != step:
            return middle.heel_deg
        denominator = left.gz - 2.0 * middle.gz + right.gz
        if denominator == 0.0:
            return middle.heel_deg
        return middle.heel_deg + 0.5 * step * (left.gz - right.gz) / denominator

    @property
    def max_volume_error(self) -> float:
        """Worst relative drift of the displacement over the sampled heels."""
        worst = 0.0
        for point in self.points:
            if not point.converged:
                continue
            worst = max(worst, abs(point.submerged_volume - self.displacement_volume))
        return worst / self.displacement_volume if self.displacement_volume else float("nan")

    def to_dict(self) -> dict[str, Any]:
        return {
            "kg_m": _round(self.kg),
            "y_g_m": _round(self.y_g),
            "keel_z_m": _round(self.keel_z),
            "waterline_z_m": _round(self.waterline_z),
            "displacement_volume_m3": _round(self.displacement_volume),
            "gz_max_m": _round(self.gz_max),
            "angle_gz_max_deg": _round(self.angle_gz_max_deg),
            "max_volume_error": _round(self.max_volume_error),
            "station_intervals": self.station_intervals,
            "heel_angles_deg": [_round(angle) for angle in self.heel_angles_deg],
            "points": [point.to_dict() for point in self.points],
            "basis": self.basis,
        }


def _round(value: float, digits: int = 9) -> float | None:
    """Round for reporting; JSON has no ``nan``, so non-finite values become ``None``."""
    return float(f"{value:.{digits}f}") if math.isfinite(value) else None


def _polygon_moments(polygon: Sequence[tuple[float, float]]) -> tuple[float, float, float]:
    """``(area, INT y dA, INT z dA)`` of a simple polygon, by Green's theorem.

    The shoelace sum is signed, so a loop traversed the other way still contributes a
    positive area with a first moment in the same frame.
    """
    signed = 0.0
    sum_y = 0.0
    sum_z = 0.0
    count = len(polygon)
    for index in range(count):
        y0, z0 = polygon[index]
        y1, z1 = polygon[(index + 1) % count]
        cross = y0 * z1 - y1 * z0
        signed += cross
        sum_y += cross * (y0 + y1)
        sum_z += cross * (z0 + z1)
    signed *= 0.5
    if signed == 0.0:
        return (0.0, 0.0, 0.0)
    sign = 1.0 if signed > 0.0 else -1.0
    return (abs(signed), sign * sum_y / 6.0, sign * sum_z / 6.0)


def _clip_below_waterplane(
    polygon: Sequence[tuple[float, float]], sin_heel: float, cos_heel: float, offset: float
) -> list[tuple[float, float]]:
    """Sutherland-Hodgman clip against ``z cos(phi) - y sin(phi) <= offset``.

    ``phi = 0`` reduces to the horizontal clip of :mod:`tensorcad.bonjean`; the same
    routine serves both because the half-plane is written as a signed value rather
    than as a ``z`` threshold.
    """
    if not polygon:
        return []
    clipped: list[tuple[float, float]] = []
    count = len(polygon)
    for index in range(count):
        y0, z0 = polygon[index]
        y1, z1 = polygon[(index + 1) % count]
        s0 = z0 * cos_heel - y0 * sin_heel - offset
        s1 = z1 * cos_heel - y1 * sin_heel - offset
        inside0 = s0 <= 0.0
        if inside0:
            clipped.append((y0, z0))
        if inside0 != (s1 <= 0.0):
            span = s0 - s1
            t = 0.0 if span == 0.0 else s0 / span
            clipped.append((y0 + t * (y1 - y0), z0 + t * (z1 - z0)))
    return clipped


@dataclass(frozen=True)
class _SectionGrid:
    """The station outlines, cut once and reused for every heel angle."""

    stations: tuple[float, ...]
    loops: tuple[tuple[tuple[tuple[float, float], ...], ...], ...]
    deck_points: tuple[tuple[tuple[float, float], ...], ...]
    spacing: float
    scale: float
    intervals: int
    keel_z: float


def section_grid(mesh: Mesh, intervals: int = DEFAULT_STATION_INTERVALS) -> _SectionGrid:
    """Cut ``mesh`` into ``intervals + 1`` station outlines spanning its length."""
    if intervals < MINIMUM_STATION_INTERVALS:
        raise OracleError(
            f"station_intervals must be at least {MINIMUM_STATION_INTERVALS}, got {intervals!r}"
        )
    if intervals % 2:
        intervals += 1  # Simpson needs an even count
    box = mesh.bounding_box()
    x0, x1 = box.minimum[0], box.maximum[0]
    if not (math.isfinite(x0) and math.isfinite(x1)) or x1 <= x0:
        raise OracleError(f"the mesh has no longitudinal extent: x in [{x0!r}, {x1!r}]")
    spacing = (x1 - x0) / intervals
    stations = tuple(x0 + spacing * index for index in range(intervals + 1))
    loops = tuple(section_loops(mesh, station) for station in stations)
    scale = max(abs(value) for value in box.maximum + box.minimum) or 1.0
    # Only a station that reaches the top of the body has a deck edge that can go
    # under; the pointed ends of a hull have none and must not be counted as safe.
    tolerance = 1e-9 * scale
    top = box.maximum[2]
    deck_points = tuple(
        tuple(point for loop in station_loops for point in loop if point[1] >= top - tolerance)
        for station_loops in loops
    )
    return _SectionGrid(
        stations=stations,
        loops=loops,
        deck_points=deck_points,
        spacing=spacing,
        scale=scale,
        intervals=intervals,
        keel_z=min(vertex[2] for vertex in mesh.vertices),
    )


def _simpson(values: Sequence[float], spacing: float) -> float:
    """Composite Simpson over an even number of intervals."""
    total = values[0] + values[-1]
    for index in range(1, len(values) - 1):
        total += (4.0 if index % 2 else 2.0) * values[index]
    return total * spacing / 3.0


def _heeled_integrals(
    grid: _SectionGrid, heel_rad: float, offset: float
) -> tuple[float, float, float]:
    """``(volume, INT A y_c dx, INT A z_c dx)`` for one heel and waterplane offset."""
    sin_heel = math.sin(heel_rad)
    cos_heel = math.cos(heel_rad)
    areas: list[float] = []
    moment_y: list[float] = []
    moment_z: list[float] = []
    for loops in grid.loops:
        area = 0.0
        my = 0.0
        mz = 0.0
        for loop in loops:
            clipped = _clip_below_waterplane(loop, sin_heel, cos_heel, offset)
            if len(clipped) < 3:
                continue
            loop_area, loop_my, loop_mz = _polygon_moments(clipped)
            area += loop_area
            my += loop_my
            mz += loop_mz
        areas.append(area)
        moment_y.append(my)
        moment_z.append(mz)
    return (
        _simpson(areas, grid.spacing),
        _simpson(moment_y, grid.spacing),
        _simpson(moment_z, grid.spacing),
    )


def _vertical_extent(grid: _SectionGrid, heel_rad: float) -> tuple[float, float]:
    """Bracket for the waterplane offset: the extreme ``s`` over all stations."""
    sin_heel = math.sin(heel_rad)
    cos_heel = math.cos(heel_rad)
    low = math.inf
    high = -math.inf
    for loops in grid.loops:
        for loop in loops:
            for y, z in loop:
                value = z * cos_heel - y * sin_heel
                low = min(low, value)
                high = max(high, value)
    return low, high


def _deck_immersion(grid: _SectionGrid, heel_rad: float, offset: float) -> float:
    """Share of deck stations whose deck edge is under the heeled waterplane."""
    sin_heel = math.sin(heel_rad)
    cos_heel = math.cos(heel_rad)
    total = 0
    awash = 0
    for points in grid.deck_points:
        if not points:
            continue
        total += 1
        if min(z * cos_heel - y * sin_heel for y, z in points) < offset:
            awash += 1
    return 0.0 if total == 0 else awash / total


@dataclass(frozen=True)
class _Immersed:
    """The heeled body at equilibrium: displacement, centre of buoyancy, waterplane."""

    volume: float
    y_b: float
    z_b: float
    offset: float


def _solve_offset(grid: _SectionGrid, heel_rad: float, displacement: float) -> float:
    """Offset of the tilted waterplane that displaces ``displacement`` (bisection).

    ``s`` is monotone in the offset -- raising the plane only adds area -- so the
    bracket over the body's bounding box always contains the root and no initial
    guess is involved.
    """
    low, high = _vertical_extent(grid, heel_rad)
    if high - low <= 0.0:  # pragma: no cover - a body with no vertical extent
        raise OracleError("the mesh has no vertical extent in the heeled frame")
    if _heeled_integrals(grid, heel_rad, high)[0] < displacement:
        raise OracleError(
            f"the body cannot displace {displacement!r} m3 at any waterline in this "
            "attitude, so it would be fully immersed -- add freeboard to the mesh"
        )
    lo, hi = low, high
    offset = high
    tolerance = OFFSET_TOLERANCE * grid.scale
    for _ in range(MAX_BISECTION_STEPS):
        offset = 0.5 * (lo + hi)
        volume = _heeled_integrals(grid, heel_rad, offset)[0]
        if volume < displacement:
            lo = offset
        else:
            hi = offset
        if hi - lo <= tolerance or abs(volume - displacement) <= 1e-12 * displacement:
            break
    return offset


def _immersed(grid: _SectionGrid, heel_rad: float, displacement: float) -> _Immersed:
    """Integrate the body at one heel angle, holding the displacement constant."""
    offset = _solve_offset(grid, heel_rad, displacement)
    volume, moment_y, moment_z = _heeled_integrals(grid, heel_rad, offset)
    if not math.isfinite(volume) or volume <= 0.0:  # pragma: no cover - guarded above
        raise OracleError("the heeled body displaced no volume")
    keel_z = grid.keel_z
    return _Immersed(volume, moment_y / volume, moment_z / volume - keel_z, offset)


def immersed_body(
    mesh: Mesh,
    heel_deg: float,
    *,
    waterline_z: float = 0.0,
    station_intervals: int = DEFAULT_STATION_INTERVALS,
    displacement_volume: float | None = None,
) -> _Immersed:
    """Volume, ``y_B``, ``z_B`` (both from the keel) and waterplane offset of one heel.

    ``displacement_volume`` defaults to the upright volume at ``waterline_z``, which is
    what holds a floating body's displacement constant along a curve. With
    ``heel_deg = 0`` the volume returned is the upright displacement, so it can be
    compared with the closed-surface integral of :func:`tensorcad.gz.heeled_hydrostatics`:
    two different integrations of the same body.
    """
    if not math.isfinite(heel_deg):
        raise OracleError(f"heel_deg must be finite, got {heel_deg!r}")
    if not math.isfinite(waterline_z):
        raise OracleError(f"waterline_z must be finite, got {waterline_z!r}")
    grid = section_grid(mesh, station_intervals)
    if displacement_volume is None:
        displacement_volume = _heeled_integrals(grid, 0.0, waterline_z)[0]
        if not math.isfinite(displacement_volume) or displacement_volume <= 0.0:
            raise OracleError(
                f"the waterline {waterline_z!r} m cuts no immersed volume; a heeled body "
                "needs one to float"
            )
    try:
        return _immersed(grid, math.radians(heel_deg), displacement_volume)
    except OracleError as error:
        raise OracleError(f"heel {heel_deg!r} deg: {error}") from error


def righting_arm_curve(
    mesh: Mesh,
    *,
    kg: float,
    heel_angles_deg: Iterable[float] | None = None,
    waterline_z: float = 0.0,
    station_intervals: int = DEFAULT_STATION_INTERVALS,
    y_g: float = 0.0,
) -> OracleCurve:
    """Integrate ``GZ(phi)`` at constant displacement over ``heel_angles_deg``.

    ``kg`` is the centre of gravity **above the keel** (the mesh's lowest vertex) and is
    required, because a default would invite the silent zero the parameter exists to
    avoid. ``waterline_z`` is the *upright* waterline in mesh coordinates.
    """
    if not all(math.isfinite(value) for value in (kg, y_g, waterline_z)):
        raise OracleError("kg, y_g and waterline_z must be finite")
    angles = DEFAULT_HEEL_ANGLES_DEG if heel_angles_deg is None else tuple(float(a) for a in heel_angles_deg)
    if not angles:
        raise OracleError("at least one heel angle is required")
    if any(not math.isfinite(angle) for angle in angles):
        raise OracleError("heel angles must be finite")
    if any(angle < 0.0 or angle > 180.0 for angle in angles):
        raise OracleError("heel angles must lie in [0, 180] degrees")
    if any(right <= left for left, right in itertools.pairwise(angles)):
        raise OracleError("heel angles must be strictly increasing")
    if 0.0 not in angles:
        raise OracleError("heel angles must include 0 degrees (the upright point is the datum)")

    keel_z = min(vertex[2] for vertex in mesh.vertices)
    deck_z = max(vertex[2] for vertex in mesh.vertices)
    if deck_z - waterline_z <= 0.0:
        raise OracleError(
            f"freeboard {deck_z - waterline_z!r} m at waterline {waterline_z!r} m: the "
            "waterline is at or above the top of the body, so every heel floods the deck"
        )
    grid = section_grid(mesh, station_intervals)
    displacement = _heeled_integrals(grid, 0.0, waterline_z)[0]
    if not math.isfinite(displacement) or displacement <= 0.0:
        raise OracleError(
            f"the waterline {waterline_z!r} m cuts no immersed volume; a righting arm "
            "needs a floating body"
        )

    points: list[OraclePoint] = []
    for heel_deg in angles:
        heel_rad = math.radians(heel_deg)
        try:
            heeled = _immersed(grid, heel_rad, displacement)
        except OracleError:
            points.append(
                OraclePoint(
                    heel_deg=heel_deg,
                    gz=float("nan"),
                    submerged_volume=float("nan"),
                    transverse_centre_buoyancy=float("nan"),
                    vertical_centre_buoyancy=float("nan"),
                    offset=float("nan"),
                    deck_immersion_fraction=float("nan"),
                    converged=False,
                )
            )
            continue
        gz = (heeled.y_b - y_g) * math.cos(heel_rad) + (heeled.z_b - kg) * math.sin(heel_rad)
        points.append(
            OraclePoint(
                heel_deg=heel_deg,
                gz=gz,
                submerged_volume=heeled.volume,
                transverse_centre_buoyancy=heeled.y_b,
                vertical_centre_buoyancy=heeled.z_b,
                offset=heeled.offset,
                deck_immersion_fraction=_deck_immersion(grid, heel_rad, heeled.offset),
            )
        )

    return OracleCurve(
        kg=kg,
        y_g=y_g,
        keel_z=keel_z,
        waterline_z=waterline_z,
        displacement_volume=displacement,
        points=tuple(points),
        station_intervals=grid.intervals,
        heel_angles_deg=angles,
    )
