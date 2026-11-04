"""Large-angle righting arms: the ``GZ`` curve, from heeled section areas.

:mod:`tensorcad.stability` answers "how far can G be raised before upright stops
being stable" and says, in its own docstring, that it does not answer the next
question: *how much righting moment is left at 30 degrees of heel?* This module
answers that one. ``GM`` is the slope of ``GZ(phi)`` at ``phi = 0``; the curve
itself is a different integration, because at finite heel the displaced volume is
bounded by a *tilted* waterplane and the centre of buoyancy moves sideways as
well as down.

**How the arm is measured.** The body is heeled instead of the water, which is
the same computation seen from a rotating frame. Heel the body ``phi`` to
starboard about the longitudinal axis and the still-water surface, expressed in
*body* coordinates, is the plane

    ``s(y, z) = z cos(phi) - y sin(phi) = offset``

with the immersed side being ``s <= offset``. (The sign is the one that puts the
starboard side *deeper* for ``phi > 0``: a point at fixed ``z`` and positive ``y``
has a smaller ``s``, so it is immersed first.) For a given heel the offset that
holds the displacement equal to the upright displacement is found by bisection --
``s`` is monotone in the offset, so the bracket ``[s_min, s_max]`` over the
bounding box always contains the root, and no initial guess is involved.

Each station of the mesh is cut with the plane ``x = station`` exactly once
(:func:`tensorcad.bonjean.section_loops`), giving closed polygons in ``(y, z)``.
Those polygons are then clipped against the tilted half-plane by
Sutherland-Hodgman and their area and first moments are taken by Green's theorem:

    ``A      = 1/2 SUM (y_i z_{i+1} - y_{i+1} z_i)``
    ``INT y dA = 1/6 SUM (y_i + y_{i+1})(y_i z_{i+1} - y_{i+1} z_i)``
    ``INT z dA = 1/6 SUM (z_i + z_{i+1})(y_i z_{i+1} - y_{i+1} z_i)``

which is *exact* for the polygon the clip produced, so the only approximation
left is the one already in the mesh: the section is a polyline and the station
grid is finite. Both are measured, not assumed -- ``examples/righting`` refines
the station grid and ``tests/test_righting.py`` pins the wall-sided case, where
the arm is a closed form.

The station areas are then integrated along the hull with composite Simpson,
which gives the heeled volume and the two transverse moments of that volume in
one pass:

    ``V  = INT A dx``   ``y_B = (INT A y_c dx) / V``   ``z_B = (INT A z_c dx) / V``

and the righting arm is the lever from G to B about the heeled axis,

    ``GZ = (y_B - y_G) cos(phi) + (z_B - z_G) sin(phi)``,

``z`` measured from the keel (the mesh's lowest vertex), ``z_G = keel + kg``,
``y_G = 0`` on the centreline unless a loading condition says otherwise. Positive
``GZ`` means the pair ``(displacement, GZ)`` is a *righting* moment.

**What this is not.** The body is treated as a rigid watertight solid with no
free surface, no flooding and no tankage: the free-surface correction
``sum(rho i_x) / displacement`` is the caller's, exactly as in
:mod:`tensorcad.stability`. Sections are assumed to be single-valued in ``z`` per
``y`` (no overhangs, no immersed deck), which is what a monohull hull below its
deck edge is; a small craft heeled past its deck edge, or a hull with a wide
transom that immerses, needs a general immersed-volume integration that this
module does not have. Damage stability, weathertight-integrity assessment and any
statutory criterion (IS Code 2008, SOLAS) are outside the module: an angle of
vanishing stability computed here is geometry, not compliance.

Only the Python standard library is required.
"""
from __future__ import annotations

import itertools
import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import Any

from .bonjean import section_loops
from .geometry import Mesh

SCHEMA_POINT = "tensorcad.righting.point/1"
SCHEMA_CURVE = "tensorcad.righting.curve/1"

#: Default heel stations: 0 to 90 degrees in 5-degree steps, which is the
#: resolution a cross-curve table is normally quoted at.
DEFAULT_HEEL_ANGLES_DEG: tuple[float, ...] = tuple(float(angle) for angle in range(0, 91, 5))

#: Number of longitudinal intervals used by :func:`righting_arm_curve` when the
#: caller does not ask for one. Simpson's rule needs an even count, and the
#: default hull lofts 40 intervals, so the station grid lands on the mesh's own
#: rings at this setting only by coincidence -- the convergence study in
#: ``examples/righting`` is what says the number is enough.
DEFAULT_STATION_INTERVALS = 80

#: Smallest number of longitudinal intervals accepted. Below this the Simpson
#: integration of a hull with pointed ends is not worth reporting.
MINIMUM_STATION_INTERVALS = 8

#: Bisection stops when the bracket is narrower than this fraction of the model
#: scale, which is far below the discretisation error of any mesh the package
#: generates.
OFFSET_TOLERANCE = 1e-12

MAX_BISECTION_STEPS = 200


class RightingError(ValueError):
    """Raised when a righting-arm curve cannot be defined for the given geometry."""


@dataclass(frozen=True)
class RightingArmPoint:
    """One heel angle of the curve.

    ``offset`` is the coordinate of the heeled waterplane in body coordinates
    (see the module docstring): it is the quantity the bisection solves for, kept
    because it is the only way to reproduce one station of the table by hand.
    ``sinkage`` is ``offset - waterline_z``, positive when the body has to go down
    to keep the displacement. ``converged`` is ``False`` when the target
    displacement is not reachable at this heel (the whole body would have to be
    immersed); such a point carries ``nan`` for the geometric quantities.

    ``deck_immersion_fraction`` is the share of the deck stations where the
    heeled waterplane has risen above the deck edge, i.e. where the low side's
    deck edge is under water (the classic ``(B/2) tan(phi) > freeboard``
    condition, applied per station instead of to the midship section). It stays
    ``0.0`` while the waterplane still cuts every deck edge, and it is the only
    honest warning this module can give: past that point the geometry is still
    integrated correctly (the body is closed), but the answer describes a
    watertight box, not a ship taking water over the deck edge.
    """

    heel_deg: float
    gz: float
    submerged_volume: float
    offset: float
    sinkage: float
    transverse_centre_buoyancy: float
    vertical_centre_buoyancy: float
    deck_immersion_fraction: float = 0.0
    converged: bool = True
    schema: str = SCHEMA_POINT

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": self.schema,
            "heel_deg": _round(self.heel_deg),
            "heel_rad": _round(math.radians(self.heel_deg)),
            "gz_m": _round(self.gz),
            "submerged_volume_m3": _round(self.submerged_volume),
            "waterplane_offset_m": _round(self.offset),
            "sinkage_m": _round(self.sinkage),
            "transverse_centre_buoyancy_m": _round(self.transverse_centre_buoyancy),
            "vertical_centre_buoyancy_m": _round(self.vertical_centre_buoyancy),
            "deck_immersion_fraction": _round(self.deck_immersion_fraction),
            "converged": self.converged,
        }

    @property
    def deck_edge_immersed(self) -> bool:
        """True once the low-side deck edge is under the heeled waterplane somewhere."""
        return self.deck_immersion_fraction > 0.0


@dataclass(frozen=True)
class RightingArmCurve:
    """A ``GZ`` curve plus the two angles a naval architect reads off it.

    ``gm_transverse`` is *not* recomputed here; it is whatever the caller passed
    in from :func:`tensorcad.stability.initial_stability` (or ``nan`` when it was
    not given), so that the small-angle slope and the curve are two independent
    numbers that can be compared -- which is exactly what
    ``tests/test_righting.py`` does.
    """

    kg: float
    y_g: float
    keel_z: float
    waterline_z: float
    freeboard: float
    displacement_volume: float
    gm_transverse: float
    points: tuple[RightingArmPoint, ...]
    station_intervals: int
    heel_angles_deg: tuple[float, ...]
    basis: str = "heeled section areas (Green), Simpson along the hull"
    schema: str = SCHEMA_CURVE

    def gz_at(self, heel_deg: float) -> float:
        """Linear interpolation of ``GZ`` between the sampled heel angles.

        Linear, deliberately: the curve is piecewise smooth and a cubic through
        five-degree samples overshoots near the peak, which is the one place the
        value matters. Values outside the sampled range raise -- extrapolating a
        righting arm past the last computed heel would be an invention.
        """
        if not self.points:
            raise RightingError("the curve has no points")
        angles = [point.heel_deg for point in self.points]
        if heel_deg < angles[0] or heel_deg > angles[-1]:
            raise RightingError(
                f"heel {heel_deg!r} deg is outside the sampled range "
                f"[{angles[0]!r}, {angles[-1]!r}]"
            )
        for index in range(1, len(angles)):
            if heel_deg <= angles[index]:
                left, right = self.points[index - 1], self.points[index]
                span = right.heel_deg - left.heel_deg
                if span == 0.0:  # pragma: no cover - duplicate angles are rejected
                    return left.gz
                weight = (heel_deg - left.heel_deg) / span
                return left.gz + weight * (right.gz - left.gz)
        return self.points[-1].gz  # pragma: no cover - covered by the range check

    @property
    def gz_max(self) -> float:
        """Largest sampled righting arm (``nan`` if nothing is positive)."""
        values = [point.gz for point in self.points if point.converged]
        if not values:
            return float("nan")
        best = max(values)
        return best if best > 0.0 else float("nan")

    @property
    def angle_gz_max_deg(self) -> float:
        """Heel angle of the largest sampled arm, refined to the parabola peak.

        The three points around the maximum are fitted with a parabola; when the
        maximum sits on the first or last sampled point there is nothing to fit
        and the sampled angle is returned.
        """
        converged = [point for point in self.points if point.converged]
        if len(converged) < 3:
            return float("nan")
        best = max(range(len(converged)), key=lambda index: converged[index].gz)
        if best == 0 or best == len(converged) - 1 or converged[best].gz <= 0.0:
            if converged[best].gz <= 0.0:
                return float("nan")
            return converged[best].heel_deg
        left, middle, right = converged[best - 1], converged[best], converged[best + 1]
        h = middle.heel_deg - left.heel_deg
        if h <= 0.0 or right.heel_deg - middle.heel_deg != h:
            return middle.heel_deg
        denominator = left.gz - 2.0 * middle.gz + right.gz
        if denominator == 0.0:
            return middle.heel_deg
        shift = 0.5 * h * (left.gz - right.gz) / denominator
        return middle.heel_deg + shift

    @property
    def angle_vanishing_stability_deg(self) -> float:
        """Heel where ``GZ`` falls back to zero after being positive (``nan`` if not)."""
        converged = [point for point in self.points if point.converged]
        for index in range(1, len(converged)):
            previous, current = converged[index - 1], converged[index]
            if previous.gz > 0.0 and current.gz <= 0.0:
                span = current.gz - previous.gz
                if span == 0.0:  # pragma: no cover - a flat crossing
                    return current.heel_deg
                weight = previous.gz / -span
                return previous.heel_deg + weight * (current.heel_deg - previous.heel_deg)
        return float("nan")

    @property
    def deck_immersion_onset_deg(self) -> float:
        """First sampled heel angle with an immersed deck edge (``nan`` if never).

        Below this angle a wall-sided hull's ``GZ`` is the closed form
        :func:`wall_sided_righting_arm`; above it the closed form is wrong by
        construction and the curve here keeps integrating the closed body. Report
        the onset rather than truncating the curve, because truncating would hide
        where the model stopped describing a ship.
        """
        for point in self.points:
            if point.converged and point.deck_edge_immersed:
                return point.heel_deg
        return float("nan")

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": self.schema,
            "kg_m": _round(self.kg),
            "y_g_m": _round(self.y_g),
            "keel_z_m": _round(self.keel_z),
            "waterline_z_m": _round(self.waterline_z),
            "freeboard_m": _round(self.freeboard),
            "displacement_volume_m3": _round(self.displacement_volume),
            "gm_transverse_input_m": _round(self.gm_transverse),
            "gz_max_m": _round(self.gz_max),
            "angle_gz_max_deg": _round(self.angle_gz_max_deg),
            "angle_vanishing_stability_deg": _round(self.angle_vanishing_stability_deg),
            "deck_immersion_onset_deg": _round(self.deck_immersion_onset_deg),
            "station_intervals": self.station_intervals,
            "heel_angles_deg": [_round(angle) for angle in self.heel_angles_deg],
            "points": [point.to_dict() for point in self.points],
            "basis": self.basis,
        }


def _round(value: float, digits: int = 9) -> float | None:
    """Round for reporting, mapping non-finite values to ``None`` (JSON has no nan)."""
    if not math.isfinite(value):
        return None
    return round(value, digits)


def _polygon_moments(polygon: Sequence[tuple[float, float]]) -> tuple[float, float, float]:
    """``(area, INT y dA, INT z dA)`` of a closed polygon, all three signed alike.

    The area is returned as a positive number and the moments are multiplied by
    the sign of the shoelace sum, so a loop traversed the other way (which
    happens for a hole, or for a body whose winding the caller rebuilt) still
    contributes a positive area with a first moment in the same frame.
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

    ``phi = 0`` reduces to the horizontal clip :mod:`tensorcad.bonjean` uses, and
    the same routine serves both because the half-plane is written as a signed
    value rather than as a ``z`` threshold.
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
        inside1 = s1 <= 0.0
        if inside0:
            clipped.append((y0, z0))
        if inside0 != inside1:
            span = s0 - s1
            t = 0.0 if span == 0.0 else s0 / span
            clipped.append((y0 + t * (y1 - y0), z0 + t * (z1 - z0)))
    return clipped


@dataclass(frozen=True)
class _SectionGrid:
    """The hull's station outlines, cut once and reused for every heel angle."""

    stations: tuple[float, ...]
    loops: tuple[tuple[tuple[tuple[float, float], ...], ...], ...]
    deck_points: tuple[tuple[tuple[float, float], ...], ...]
    spacing: float
    scale: float
    intervals: int


def _section_grid(mesh: Mesh, intervals: int) -> _SectionGrid:
    if intervals < MINIMUM_STATION_INTERVALS:
        raise RightingError(
            f"station_intervals must be at least {MINIMUM_STATION_INTERVALS}, got {intervals!r}"
        )
    if intervals % 2:
        intervals += 1  # Simpson needs an even count of intervals
    box = mesh.bounding_box()
    x0, x1 = box.minimum[0], box.maximum[0]
    if not (math.isfinite(x0) and math.isfinite(x1)) or x1 <= x0:
        raise RightingError(f"the mesh has no longitudinal extent: x in [{x0!r}, {x1!r}]")
    spacing = (x1 - x0) / intervals
    stations = tuple(x0 + spacing * index for index in range(intervals + 1))
    loops = tuple(section_loops(mesh, station) for station in stations)
    scale = max(abs(value) for value in box.maximum + box.minimum) or 1.0
    # A "deck station" is one whose section reaches the top of the body: only
    # those have a deck edge that can go under. A station that stops lower (the
    # pointed ends of a hull, the flanks of a sphere) has no deck and is skipped
    # by the immersion flag rather than counted as safe.
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


@dataclass(frozen=True)
class ImmersedBody:
    """The submerged part of a closed body below one (possibly heeled) waterplane.

    This is the quantity the whole module is built on, exposed because it is also
    the cross-check: with ``heel_deg = 0`` and the waterplane at the top of a mesh
    capped at its design waterline, :attr:`volume` must equal the closed-surface
    integral :func:`tensorcad.hydrostatics.compute_hydrostatics` reports -- two
    different integrations of the same body -- and
    :attr:`vertical_centre_buoyancy` must equal its ``VCB``.
    """

    volume: float
    transverse_centre_buoyancy: float
    vertical_centre_buoyancy: float
    heel_deg: float
    waterline_z: float
    station_intervals: int
    basis: str = "heeled section areas (Green), Simpson along the hull"

    @property
    def transverse_centre_buoyancy_above_keel(self) -> float:
        """``y_B`` is report-only; this exists so the name says what the datum is."""
        return self.transverse_centre_buoyancy

    def to_dict(self) -> dict[str, Any]:
        return {
            "volume_m3": _round(self.volume),
            "transverse_centre_buoyancy_m": _round(self.transverse_centre_buoyancy),
            "vertical_centre_buoyancy_m": _round(self.vertical_centre_buoyancy),
            "heel_deg": _round(self.heel_deg),
            "waterline_z_m": _round(self.waterline_z),
            "station_intervals": self.station_intervals,
            "basis": self.basis,
        }


def submerged_body(
    mesh: Mesh,
    *,
    heel_deg: float = 0.0,
    waterline_z: float = 0.0,
    station_intervals: int = DEFAULT_STATION_INTERVALS,
) -> ImmersedBody:
    """Integrate the part of ``mesh`` below the waterplane ``s = waterline_z``.

    ``heel_deg`` tilts the waterplane about the longitudinal axis (positive heels
    the starboard side down); ``waterline_z`` is where that plane sits when
    ``heel_deg`` is zero. The heeled case is the same integral
    :func:`righting_arm_curve` solves for a constant displacement, offered here
    without the displacement constraint so that a reviewer can watch the volume
    change with the plane.
    """
    if not math.isfinite(heel_deg):
        raise RightingError(f"heel_deg must be finite, got {heel_deg!r}")
    if not math.isfinite(waterline_z):
        raise RightingError(f"waterline_z must be finite, got {waterline_z!r}")
    grid = _section_grid(mesh, station_intervals)
    volume, moment_y, moment_z = _heeled_integrals(grid, math.radians(heel_deg), waterline_z)
    if volume <= 0.0:
        return ImmersedBody(
            volume=volume,
            transverse_centre_buoyancy=float("nan"),
            vertical_centre_buoyancy=float("nan"),
            heel_deg=heel_deg,
            waterline_z=waterline_z,
            station_intervals=grid.intervals,
        )
    return ImmersedBody(
        volume=volume,
        transverse_centre_buoyancy=moment_y / volume,
        vertical_centre_buoyancy=moment_z / volume,
        heel_deg=heel_deg,
        waterline_z=waterline_z,
        station_intervals=grid.intervals,
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
    """Share of deck stations whose deck edge is under the heeled waterplane.

    The heeled waterplane is ``s = offset`` with ``s = z cos(phi) - y sin(phi)``,
    so a deck point is under water when its ``s`` is below the offset. The lowest
    deck point of a station is the deck edge on the low (immersed) side, which is
    the classic ``(B/2) tan(phi) > freeboard`` limit once a flat deck is assumed.
    """
    sin_heel = math.sin(heel_rad)
    cos_heel = math.cos(heel_rad)
    total = 0
    awash = 0
    for points in grid.deck_points:
        if not points:
            continue
        total += 1
        lowest = min(z * cos_heel - y * sin_heel for y, z in points)
        if lowest < offset:
            awash += 1
    if total == 0:
        return 0.0
    return awash / total


def righting_arm_curve(
    mesh: Mesh,
    *,
    kg: float,
    waterline_z: float = 0.0,
    heel_angles_deg: Iterable[float] | None = None,
    station_intervals: int = DEFAULT_STATION_INTERVALS,
    y_g: float = 0.0,
    gm_transverse: float = float("nan"),
) -> RightingArmCurve:
    """Integrate the righting-arm curve of a floating body at constant displacement.

    ``kg`` is the centre of gravity **above the keel** (the mesh's lowest vertex)
    and is a required keyword for the same reason as in
    :func:`tensorcad.stability.initial_stability`: a default would invite the
    silent zero the parameter exists to avoid. ``waterline_z`` is the *upright*
    waterplane; the freeboard available to heel with is the body above it, so a
    hull capped at its design waterline should be floated at a negative
    ``waterline_z`` if a curve with a vanishing angle is wanted.

    The upright displacement is integrated with the same station grid as the
    heeled ones, so ``GZ(0)`` comes out as round-off zero rather than as the
    difference of two discretisation errors; the price is that the displacement
    quoted here is the section integral, which differs from
    :func:`tensorcad.hydrostatics.compute_hydrostatics` by the usual grid error.

    ``gm_transverse`` is copied into the report only -- pass the value from
    :func:`tensorcad.stability.initial_stability` and the two independent
    integrations can be checked against each other.
    """
    if not math.isfinite(kg):
        raise RightingError(f"kg must be finite, got {kg!r}")
    if not math.isfinite(y_g):
        raise RightingError(f"y_g must be finite, got {y_g!r}")
    if not math.isfinite(waterline_z):
        raise RightingError(f"waterline_z must be finite, got {waterline_z!r}")

    angles = DEFAULT_HEEL_ANGLES_DEG if heel_angles_deg is None else tuple(float(a) for a in heel_angles_deg)
    if not angles:
        raise RightingError("at least one heel angle is required")
    if any(not math.isfinite(angle) for angle in angles):
        raise RightingError("heel angles must be finite")
    if any(angle < 0.0 or angle > 180.0 for angle in angles):
        raise RightingError("heel angles must lie in [0, 180] degrees")
    if any(right <= left for left, right in itertools.pairwise(angles)):
        raise RightingError("heel angles must be strictly increasing")
    if 0.0 not in angles:
        raise RightingError("heel angles must include 0 degrees (the upright point is the datum)")

    keel_z = min(vertex[2] for vertex in mesh.vertices)
    deck_z = max(vertex[2] for vertex in mesh.vertices)
    freeboard = deck_z - waterline_z
    if freeboard <= 0.0:
        raise RightingError(
            f"freeboard {freeboard!r} m at waterline {waterline_z!r} m: the waterline is at or "
            f"above the top of the body ({deck_z!r} m), so every heel angle floods the deck. "
            "Float the body lower (negative waterline_z) or supply a mesh that includes freeboard."
        )
    grid = _section_grid(mesh, station_intervals)
    upright_volume, _, _ = _heeled_integrals(grid, 0.0, waterline_z)
    if not math.isfinite(upright_volume) or upright_volume <= 0.0:
        raise RightingError(
            f"the waterline {waterline_z!r} m cuts no immersed volume; a righting arm needs a floating body"
        )

    z_g = keel_z + kg
    points: list[RightingArmPoint] = []
    for heel_deg in angles:
        heel_rad = math.radians(heel_deg)
        low, high = _vertical_extent(grid, heel_rad)
        if high - low <= 0.0:  # pragma: no cover - a body with no vertical extent
            points.append(
                RightingArmPoint(
                    heel_deg=heel_deg,
                    gz=float("nan"),
                    submerged_volume=float("nan"),
                    offset=float("nan"),
                    sinkage=float("nan"),
                    transverse_centre_buoyancy=float("nan"),
                    vertical_centre_buoyancy=float("nan"),
                    converged=False,
                )
            )
            continue
        tolerance = OFFSET_TOLERANCE * grid.scale
        volume_high, _, _ = _heeled_integrals(grid, heel_rad, high)
        if volume_high < upright_volume:
            points.append(
                RightingArmPoint(
                    heel_deg=heel_deg,
                    gz=float("nan"),
                    submerged_volume=volume_high,
                    offset=high,
                    sinkage=high - waterline_z,
                    transverse_centre_buoyancy=float("nan"),
                    vertical_centre_buoyancy=float("nan"),
                    converged=False,
                )
            )
            continue
        lo, hi = low, high
        offset = high
        for _ in range(MAX_BISECTION_STEPS):
            offset = 0.5 * (lo + hi)
            volume, _, _ = _heeled_integrals(grid, heel_rad, offset)
            if volume < upright_volume:
                lo = offset
            else:
                hi = offset
            if hi - lo <= tolerance or abs(volume - upright_volume) <= 1e-12 * upright_volume:
                break
        volume, moment_y, moment_z = _heeled_integrals(grid, heel_rad, offset)
        if not math.isfinite(volume) or volume <= 0.0:  # pragma: no cover - guarded above
            raise RightingError(f"heel {heel_deg!r} deg produced no immersed volume")
        y_b = moment_y / volume
        z_b = moment_z / volume
        immersion = _deck_immersion(grid, heel_rad, offset)
        gz = (y_b - y_g) * math.cos(heel_rad) + (z_b - z_g) * math.sin(heel_rad)
        points.append(
            RightingArmPoint(
                heel_deg=heel_deg,
                gz=gz,
                submerged_volume=volume,
                offset=offset,
                sinkage=offset - waterline_z,
                transverse_centre_buoyancy=y_b,
                vertical_centre_buoyancy=z_b,
                deck_immersion_fraction=immersion,
            )
        )

    return RightingArmCurve(
        kg=kg,
        y_g=y_g,
        keel_z=keel_z,
        waterline_z=waterline_z,
        freeboard=freeboard,
        displacement_volume=upright_volume,
        gm_transverse=gm_transverse,
        points=tuple(points),
        station_intervals=grid.intervals,
        heel_angles_deg=angles,
    )


def wall_sided_righting_arm(heel_deg: float, *, gm: float, bm: float) -> float:
    """The closed form the numerical curve is checked against.

    For a body whose sides stay vertical through the heeled waterplane (a box, a
    barge, a hull up to its deck edge) the exact righting arm is

        ``GZ = sin(phi) (GM + BM tan^2(phi) / 2)``

    with ``BM = I_T / V`` and ``GM = KB + BM - KG``. It is *not* an
    approximation in ``phi``: it is what the waterplane widening ``B/cos(phi)``
    and the accompanying trapezoid centroid shift add up to exactly, and it fails
    only when the waterplane leaves the vertical sides. ``examples/righting`` and
    ``tests/test_righting.py`` use it as the reference for the box barge; nothing
    in this module depends on it.
    """
    if not math.isfinite(heel_deg):
        raise RightingError(f"heel_deg must be finite, got {heel_deg!r}")
    phi = math.radians(heel_deg)
    return math.sin(phi) * (gm + 0.5 * bm * math.tan(phi) ** 2)
