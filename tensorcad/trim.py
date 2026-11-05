"""Trim, and the two particulars that set it: TPC and MCT.

:mod:`tensorcad.hydrostatics` reports where buoyancy *is* (volume, LCB, LCF, KB);
:mod:`tensorcad.stability` reports whether it is *upright* (KM, GM, free surface);
:mod:`tensorcad.bonjean` reports the stations a bending calculation needs. None of
them answers the question a loading computer asks first, and the one that decides
whether a displacement is even valid: **given that a weight has moved along the
hull, how far does the ship trim, and what do a centimetre of immersion and a
centimetre of trim cost?** The classic constants are

* ``TPC``, tonnes per centimetre of immersion -- ``rho * A_wp / 1e5`` with ``rho``
  in kg/m^3 and ``A_wp`` in m^2 (the 1e5 is 1000 kg/t times 100 cm/m);
* ``MCT``, the moment that changes the trim by one centimetre --
  ``rho_t * I_L / (100 * L)`` in tonne-metres per centimetre, where ``rho_t`` is
  in t/m^3 and ``I_L`` is the waterplane's second moment about the transverse
  axis through its own centroid, the axis the trimming moment is taken about.

Both are pure derivatives of quantities this package already integrates, so the
interesting part is not the arithmetic, it is the *evidence*. This module therefore
computes every longitudinal quantity twice, by routes that share no code:

* the **section route** -- integrate the immersed area of each station,
  :func:`tensorcad.bonjean.section_properties`, along the hull with Simpson's
  rule, and divide the first moment by the area to get LCB. This is the route a
  naval architect does by hand and it never looks at a face normal;
* the **waterplane route** -- the plane section's own moments
  (:func:`tensorcad.stability.waterplane_moments`), which are what BM, KM and
  therefore MCT are built from.

The two are compared in the docs and in the tests, and the residual is the
discretisation error of the hull, not of either route.

Trim is where the discretisation stops being bookkeeping. The waterline of a
trimmed ship is the plane ``z(x) = z_ref + (x - x_ref) * tan(theta)``, and the
volume below it is -- exactly -- the integral of the station areas *clipped at
their own station's waterline height*. No solid clipping is needed: at a fixed
station, "below a tilted plane" is a half plane in ``(y, z)``, which is what
:func:`tensorcad.bonjean.section_properties` already clips. So
:func:`trimmed_hydrostatics` evaluates the same section machinery at a per-station
waterline, :func:`waterline_for_displacement` bisects the mean height until the
displaced volume is preserved, and :func:`equilibrium_trim` bisects the angle
until the LCB has moved by ``w * d / displacement``. Every one of those solves is
a conservation law, not a fitted constant -- which is why a box-shaped hull comes
out bit-exact (its station area is *linear* in ``x``, so Simpson is exact) and a
curved hull's residual converges as the mesh is refined.

Sign convention, stated once because trim conventions are a classic source of
sign errors: ``x`` increases from the ``-x`` end to the ``+x`` end of the bounding
box, ``trim_angle > 0`` raises the waterline toward ``+x`` (so the ``+x`` end
immerses deeper), and a weight shifted by ``shift > 0`` moves toward ``+x``.
``draft_fore`` is measured at the ``+x`` end (the stem, for the hull families in
:mod:`tensorcad.hull` and :mod:`tensorcad.wigley`, whose stations are ordered the
same way), ``draft_aft`` at the ``-x`` end.

Standard library only, like the rest of the core.
"""
from __future__ import annotations

import bisect
import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

from .bonjean import STATION_NUDGE, SectionProperties, section_properties
from .geometry import Mesh
from .stability import StabilityError, WaterplaneMoments, waterplane_moments

#: Density of sea water, kg/m^3, used when the caller does not say otherwise.
SEAWATER_DENSITY = 1025.0

#: Standard gravity, m/s^2, reported so that mass and force moments can be converted.
GRAVITY = 9.80665

SCHEMA_PARTICULARS = "tensorcad.hydrostatic_particulars/1"
SCHEMA_TRIMMED = "tensorcad.trimmed_hydrostatics/1"
SCHEMA_TRIM_SOLUTION = "tensorcad.trim_solution/1"


class TrimError(ValueError):
    """A trim or particular was asked for that does not have an answer."""


def _round(value: float, digits: int = 9) -> float:
    if not isinstance(value, float) or not math.isfinite(value):
        return value
    return round(value, digits)


def _stations(mesh: Mesh, station_count: int) -> tuple[float, ...]:
    """``station_count`` evenly spaced stations, made odd so Simpson's rule applies."""
    box = mesh.bounding_box()
    if station_count < 3:
        raise TrimError(f"station_count must be at least 3, got {station_count!r}")
    count = station_count + 1 if station_count % 2 == 0 else station_count
    x_min, x_max = box.minimum[0], box.maximum[0]
    return tuple(x_min + (x_max - x_min) * index / (count - 1) for index in range(count))


def _integrate(stations: Sequence[float], values: Sequence[float]) -> float:
    """Simpson when the stations are uniform and odd in number, trapezoid otherwise."""
    count = len(values)
    if count < 2:
        return 0.0
    spacing = (stations[-1] - stations[0]) / (count - 1)
    uniform = all(
        abs((stations[index + 1] - stations[index]) - spacing) <= 1e-9 * max(1.0, abs(spacing))
        for index in range(count - 1)
    )
    if uniform and count % 2 == 1:
        total = values[0] + values[-1]
        for index in range(1, count - 1):
            total += (4.0 if index % 2 == 1 else 2.0) * values[index]
        return total * spacing / 3.0
    total = 0.0
    for index in range(count - 1):
        total += 0.5 * (values[index] + values[index + 1]) * (stations[index + 1] - stations[index])
    return total


class _SectionCutter:
    """Station areas at arbitrary waterline heights, one face pass per bucket.

    :func:`tensorcad.bonjean.section_properties` is the reference route, and it is
    what a caller should use for a one-off station. A trim solve, however, needs a
    few hundred *sweeps* (every station, and the waterline height changes at every
    iteration), and re-scanning every face for every station would make a loading
    calculation cost seconds. This class indexes the faces by the stations they can
    straddle -- a face spans one station interval, so it belongs to two or three
    buckets -- and then reproduces the plane cut, the clip at the waterline and the
    shoelace area for one station only. It must agree with the reference route to
    round-off, and ``test_the_indexed_cutter_matches_the_reference_route`` holds it
    to that; the two are independent implementations of the same integral, so
    agreement is evidence rather than a tautology.
    """

    __slots__ = ("_buckets", "_nudge", "mesh", "scale", "stations")

    def __init__(self, mesh: Mesh, stations: Sequence[float]) -> None:
        self.mesh = mesh
        self.stations = tuple(stations)
        box = mesh.bounding_box()
        self.scale = max(abs(value) for value in box.maximum) or 1.0
        self._nudge = STATION_NUDGE * self.scale
        buckets: list[list[int]] = [[] for _ in self.stations]
        vertices = mesh.vertices
        tolerance = 2.0 * self._nudge
        for index, face in enumerate(mesh.faces):
            xs = (vertices[face[0]][0], vertices[face[1]][0], vertices[face[2]][0])
            low, high = min(xs), max(xs)
            start = bisect.bisect_left(self.stations, low - tolerance) - 1
            stop = bisect.bisect_right(self.stations, high + tolerance) + 1
            for station_index in range(max(start, 0), min(stop, len(self.stations))):
                buckets[station_index].append(index)
        self._buckets = buckets

    def _plane(self, station_index: int) -> float:
        box = self.mesh.bounding_box()
        station = self.stations[station_index]
        offset = self._nudge
        if station + offset > box.maximum[0]:
            return station - offset
        return station + offset

    def area(self, station_index: int, waterline_z: float) -> float:
        """Clipped section area at station ``station_index`` below ``waterline_z``."""
        x0 = self._plane(station_index)
        vertices = self.mesh.vertices
        points: dict[tuple[int, int], tuple[float, float]] = {}
        adjacency: dict[tuple[int, int], list[tuple[int, int]]] = {}
        for face_index in self._buckets[station_index]:
            face = self.mesh.faces[face_index]
            crossings: list[tuple[int, int]] = []
            for position in range(3):
                a = face[position]
                b = face[(position + 1) % 3]
                da = vertices[a][0] - x0
                db = vertices[b][0] - x0
                if (da > 0.0) == (db > 0.0):
                    continue
                key = (a, b) if a < b else (b, a)
                if key not in points:
                    va = vertices[key[0]]
                    vb = vertices[key[1]]
                    da = va[0] - x0
                    db = vb[0] - x0
                    t = da / (da - db)
                    points[key] = (va[1] + t * (vb[1] - va[1]), va[2] + t * (vb[2] - va[2]))
                crossings.append(key)
            if len(crossings) != 2:
                continue
            adjacency.setdefault(crossings[0], []).append(crossings[1])
            adjacency.setdefault(crossings[1], []).append(crossings[0])
        if not adjacency:
            return 0.0
        unused: set[tuple[tuple[int, int], tuple[int, int]]] = set()
        for key, neighbours in adjacency.items():
            for neighbour in neighbours:
                unused.add((key, neighbour) if key < neighbour else (neighbour, key))
        total = 0.0
        while unused:
            first, second = next(iter(unused))
            unused.discard((first, second))
            chain = [first, second]
            previous, current = first, second
            while current != chain[0]:
                following = None
                for candidate in adjacency.get(current, []):
                    if candidate == previous:
                        continue
                    edge = (current, candidate) if current < candidate else (candidate, current)
                    if edge in unused:
                        following = candidate
                        break
                if following is None:
                    break
                edge = (current, following) if current < following else (following, current)
                unused.discard(edge)
                chain.append(following)
                previous, current = current, following
            polygon = [points[key] for key in chain]

            clipped = _clip_below(polygon, waterline_z)
            if len(clipped) < 3:
                continue
            total += abs(_signed_area(clipped))
        return total


def _clip_below(polygon: Sequence[tuple[float, float]], waterline_z: float) -> list[tuple[float, float]]:
    """Sutherland-Hodgman clip of a ``(y, z)`` polygon to ``z <= waterline_z``."""
    clipped: list[tuple[float, float]] = []
    count = len(polygon)
    for index in range(count):
        current = polygon[index]
        following = polygon[(index + 1) % count]
        current_in = current[1] <= waterline_z
        following_in = following[1] <= waterline_z
        if current_in:
            clipped.append(current)
        if current_in != following_in:
            span = following[1] - current[1]
            t = 0.0 if span == 0.0 else (waterline_z - current[1]) / span
            clipped.append((current[0] + t * (following[0] - current[0]), waterline_z))
    return clipped


def _signed_area(polygon: Sequence[tuple[float, float]]) -> float:
    total = 0.0
    count = len(polygon)
    for index in range(count):
        y0, z0 = polygon[index]
        y1, z1 = polygon[(index + 1) % count]
        total += y0 * z1 - y1 * z0
    return 0.5 * total


@dataclass(frozen=True)
class HydrostaticParticulars:
    """The particulars of one level waterline, including TPC and MCT.

    ``volume``, ``lcb`` and ``kb`` come from the section route (Simpson over station
    areas); ``waterplane_area``, ``lcf`` and the second moments come from the plane
    section's own moments. ``draft`` is measured from the mesh's lowest point, which
    is the keel for every hull in this package.
    """

    waterline_z: float
    keel_z: float
    draft: float
    length: float
    volume: float
    displacement_mass: float
    lcb: float
    lcf: float
    kb: float
    waterplane_area: float
    transverse_second_moment: float
    longitudinal_second_moment: float
    bm_transverse: float
    bm_longitudinal: float
    km_transverse: float
    km_longitudinal: float
    tpc: float
    mct: float
    sections: int
    water_density: float
    waterline_nudged: bool = False
    basis: str = "section areas (Simpson) + plane-section moments"
    schema: str = SCHEMA_PARTICULARS

    @property
    def lcb_relative_to_lcf(self) -> float:
        """``LCB - LCF``: the lever the longitudinal weight distribution trims about."""
        return self.lcb - self.lcf

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": self.schema,
            "waterline_z_m": _round(self.waterline_z),
            "keel_z_m": _round(self.keel_z),
            "draft_m": _round(self.draft),
            "length_m": _round(self.length),
            "displacement_volume_m3": _round(self.volume),
            "displacement_mass_kg": _round(self.displacement_mass, 3),
            "lcb_m": _round(self.lcb),
            "lcf_m": _round(self.lcf),
            "lcb_minus_lcf_m": _round(self.lcb_relative_to_lcf),
            "kb_m": _round(self.kb),
            "waterplane_area_m2": _round(self.waterplane_area),
            "transverse_second_moment_m4": _round(self.transverse_second_moment),
            "longitudinal_second_moment_m4": _round(self.longitudinal_second_moment),
            "bm_transverse_m": _round(self.bm_transverse),
            "bm_longitudinal_m": _round(self.bm_longitudinal),
            "km_transverse_m": _round(self.km_transverse),
            "km_longitudinal_m": _round(self.km_longitudinal),
            "tpc_t_per_cm": _round(self.tpc),
            "mct_t_m_per_cm": _round(self.mct),
            "sections": self.sections,
            "water_density_kg_per_m3": _round(self.water_density, 3),
            "waterline_nudged": self.waterline_nudged,
            "basis": self.basis,
        }


def hydrostatic_particulars(
    mesh: Mesh,
    waterline_z: float = 0.0,
    *,
    water_density: float = SEAWATER_DENSITY,
    station_count: int = 41,
) -> HydrostaticParticulars:
    """Particulars of the level waterline ``waterline_z``, TPC and MCT included.

    ``waterline_z`` is a height in the mesh's own frame, exactly as in
    :func:`tensorcad.stability.waterplane_moments`; it has to cross the body, so a
    waterline below the keel or above the deck is an error rather than a zero.
    """
    if not math.isfinite(waterline_z):
        raise TrimError(f"waterline_z must be finite, got {waterline_z!r}")
    if not math.isfinite(water_density) or water_density <= 0.0:
        raise TrimError(f"water_density must be positive, got {water_density!r}")
    box = mesh.bounding_box()
    keel_z = box.minimum[2]
    draft = waterline_z - keel_z
    if draft <= 0.0:
        raise TrimError(
            f"waterline_z={waterline_z!r} is at or below the keel ({keel_z!r}); nothing is immersed"
        )
    stations = _stations(mesh, station_count)
    sections: list[SectionProperties] = [section_properties(mesh, x, waterline_z) for x in stations]
    areas = [section.area for section in sections]
    volume = _integrate(stations, areas)
    if volume <= 0.0:
        raise TrimError(f"the waterline at z={waterline_z!r} displaces no volume")
    lcb = _integrate(stations, [area * x for area, x in zip(areas, stations, strict=True)]) / volume
    centroid_z = _integrate(
        stations,
        [area * section.centroid_z for area, section in zip(areas, sections, strict=True)],
    )
    kb = centroid_z / volume - keel_z

    nudged = False
    try:
        waterplane: WaterplaneMoments = waterplane_moments(mesh, waterline_z)
    except StabilityError:
        # The generated hulls put vertex rings on evenly spaced levels, so a "round"
        # waterline such as -T/2 lands exactly on a ring where the plane cut is
        # degenerate. ``bonjean`` nudges its station planes for the same reason; the
        # offset here is 1e-6 of the mesh scale, which moves the waterplane area by
        # O(1e-6) relative -- three orders below the discretisation error it is used
        # to measure -- and the report says so rather than hiding it.
        offset = 1e-6 * max(abs(value) for value in mesh.bounding_box().maximum)
        try:
            waterplane = waterplane_moments(mesh, waterline_z - offset)
        except StabilityError:
            raise
        nudged = True
    if waterplane.area <= 0.0:
        raise TrimError(f"the waterplane at z={waterline_z!r} has no area")
    length = box.sizes[0]
    bm_t = waterplane.transverse_second_moment / volume
    bm_l = waterplane.longitudinal_second_moment / volume
    density_t = water_density / 1000.0  # t/m^3
    return HydrostaticParticulars(
        waterline_z=waterline_z,
        keel_z=keel_z,
        draft=draft,
        length=length,
        volume=volume,
        displacement_mass=water_density * volume,
        lcb=lcb,
        lcf=waterplane.longitudinal_centre_flotation,
        kb=kb,
        waterplane_area=waterplane.area,
        transverse_second_moment=waterplane.transverse_second_moment,
        longitudinal_second_moment=waterplane.longitudinal_second_moment,
        bm_transverse=bm_t,
        bm_longitudinal=bm_l,
        km_transverse=kb + bm_t,
        km_longitudinal=kb + bm_l,
        tpc=density_t * waterplane.area / 100.0,
        mct=density_t * waterplane.longitudinal_second_moment / (100.0 * length),
        sections=len(stations),
        water_density=water_density,
        waterline_nudged=nudged,
    )


def hydrostatic_table(
    mesh: Mesh,
    waterlines: Sequence[float],
    *,
    water_density: float = SEAWATER_DENSITY,
    station_count: int = 41,
) -> tuple[HydrostaticParticulars, ...]:
    """The hydrostatic table: one :class:`HydrostaticParticulars` per waterline.

    This is the deliverable a loading manual prints -- draft against displacement,
    LCB, LCF, KB, KM and the two centimetre constants -- and asking for it is the
    cheapest way to see whether TPC and MCT vary the way the geometry says they
    should (TPC tracks the waterplane area, MCT the waterplane's length cubed).
    """
    rows = tuple(
        hydrostatic_particulars(
            mesh, waterline, water_density=water_density, station_count=station_count
        )
        for waterline in waterlines
    )
    if not rows:
        raise TrimError("hydrostatic_table needs at least one waterline")
    return rows


@dataclass(frozen=True)
class TrimmedHydrostatics:
    """A hull floating at a trim angle, with the displaced volume held fixed.

    ``waterline_at_reference`` is the height of the trimmed waterline at
    ``reference_x``; the plane itself is ``z(x) = waterline_at_reference +
    (x - reference_x) * tan_trim``. ``volume`` and ``lcb`` are the section route of
    :func:`hydrostatic_particulars`, evaluated at each station's own waterline.
    """

    trim_angle_rad: float
    trim_angle_deg: float
    tan_trim: float
    reference_x: float
    waterline_at_reference: float
    volume: float
    lcb: float
    draft_fore: float
    draft_aft: float
    keel_z: float
    sections: int
    basis: str = "per-station clipped section areas (Simpson)"
    schema: str = SCHEMA_TRIMMED

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": self.schema,
            "trim_angle_rad": _round(self.trim_angle_rad),
            "trim_angle_deg": _round(self.trim_angle_deg, 9),
            "tan_trim": _round(self.tan_trim),
            "reference_x_m": _round(self.reference_x),
            "waterline_at_reference_m": _round(self.waterline_at_reference),
            "displacement_volume_m3": _round(self.volume),
            "lcb_m": _round(self.lcb),
            "draft_fore_m": _round(self.draft_fore),
            "draft_aft_m": _round(self.draft_aft),
            "keel_z_m": _round(self.keel_z),
            "sections": self.sections,
            "basis": self.basis,
        }


def _reference_x(mesh: Mesh) -> float:
    box = mesh.bounding_box()
    return 0.5 * (box.minimum[0] + box.maximum[0])


def _trimmed_volume_and_lcb(
    cutter: _SectionCutter,
    stations: Sequence[float],
    *,
    tan_trim: float,
    waterline_at_reference: float,
    reference_x: float,
) -> tuple[float, float]:
    """Volume and LCB below the tilted waterline, by the section route."""
    areas = [
        cutter.area(index, waterline_at_reference + (x - reference_x) * tan_trim)
        for index, x in enumerate(stations)
    ]
    volume = _integrate(stations, areas)
    if volume <= 0.0:
        return 0.0, float("nan")
    lcb = _integrate(stations, [area * x for area, x in zip(areas, stations, strict=True)]) / volume
    return volume, lcb


def trimmed_hydrostatics(
    mesh: Mesh,
    *,
    trim_angle: float = 0.0,
    waterline_at_reference: float = 0.0,
    reference_x: float | None = None,
    preserved: float | None = None,
    station_count: int = 41,
    tolerance: float = 1e-10,
    max_iterations: int = 200,
) -> TrimmedHydrostatics:
    """The trimmed condition at ``trim_angle`` (radians).

    ``preserved`` is the displaced volume to hold: when it is given the mean height
    of the waterline is solved for, which is what floats a ship at a fixed
    displacement with a prescribed trim; when it is ``None`` the waterline is taken
    as given at ``reference_x``.
    """
    if not math.isfinite(trim_angle):
        raise TrimError(f"trim_angle must be finite, got {trim_angle!r}")
    stations = _stations(mesh, station_count)
    x_ref = _reference_x(mesh) if reference_x is None else float(reference_x)
    tan_trim = math.tan(trim_angle)
    height = waterline_at_reference
    if preserved is not None:
        height = waterline_for_displacement(
            mesh,
            trim_angle=trim_angle,
            displaced_volume=float(preserved),
            reference_x=x_ref,
            station_count=station_count,
            tolerance=tolerance,
            max_iterations=max_iterations,
        )
    volume, lcb = _trimmed_volume_and_lcb(
        _SectionCutter(mesh, stations),
        stations,
        tan_trim=tan_trim,
        waterline_at_reference=height,
        reference_x=x_ref,
    )
    if volume <= 0.0:
        raise TrimError(
            f"trim_angle={trim_angle!r} with waterline {height!r} displaces no volume"
        )
    box = mesh.bounding_box()
    keel_z = box.minimum[2]
    return TrimmedHydrostatics(
        trim_angle_rad=trim_angle,
        trim_angle_deg=math.degrees(trim_angle),
        tan_trim=tan_trim,
        reference_x=x_ref,
        waterline_at_reference=height,
        volume=volume,
        lcb=lcb,
        draft_fore=height + (box.maximum[0] - x_ref) * tan_trim - keel_z,
        draft_aft=height + (box.minimum[0] - x_ref) * tan_trim - keel_z,
        keel_z=keel_z,
        sections=len(stations),
    )


def _false_position(a: float, fa: float, b: float, fb: float) -> float:
    """Regula falsi step; a flat residual (no root between the brackets) bisects."""
    if fb == fa:  # pragma: no cover - a flat residual cannot bracket a root
        return 0.5 * (a + b)
    return b - fb * (b - a) / (fb - fa)


def waterline_for_displacement(
    mesh: Mesh,
    *,
    trim_angle: float,
    displaced_volume: float,
    reference_x: float | None = None,
    station_count: int = 41,
    tolerance: float = 1e-10,
    volume_tolerance: float = 1e-12,
    max_iterations: int = 200,
    initial_height: float | None = None,
    cutter: _SectionCutter | None = None,
) -> float:
    """Mean waterline height at ``reference_x`` that displaces ``displaced_volume``.

    Displaced volume is strictly increasing in the waterline height on a closed
    body -- raising the plane can only add area to every station -- so the root is
    unique and bracketed by the bounding box widened by ``L * |tan|`` on each side,
    which covers the tilted plane's excursion along the hull. The root is found by
    a safeguarded secant (Illinois) iteration rather than plain bisection: volume is
    smooth in the height, so the secant lands in a handful of evaluations, and any
    step that leaves the bracket is replaced by a bisection step, so the iteration
    cannot diverge. ``initial_height`` warm-starts the search from a neighbouring
    trim angle, which is what makes :func:`equilibrium_trim` affordable: the inner
    solve then costs two or three section sweeps instead of forty.
    """
    if not math.isfinite(displaced_volume) or displaced_volume <= 0.0:
        raise TrimError(f"displaced_volume must be positive, got {displaced_volume!r}")
    stations = cutter.stations if cutter is not None else _stations(mesh, station_count)
    index = cutter if cutter is not None else _SectionCutter(mesh, stations)
    box = mesh.bounding_box()
    x_ref = _reference_x(mesh) if reference_x is None else float(reference_x)
    tan_trim = math.tan(trim_angle)
    slack = box.sizes[0] * abs(tan_trim) + box.sizes[2]
    low = box.minimum[2] - slack
    high = box.maximum[2] + slack

    def _residual(height: float) -> float:
        return _trimmed_volume_and_lcb(
            index,
            stations,
            tan_trim=tan_trim,
            waterline_at_reference=height,
            reference_x=x_ref,
        )[0] - displaced_volume

    a, fa = low, _residual(low)
    b, fb = high, _residual(high)
    if fb < 0.0:
        raise TrimError(
            f"even the highest waterline displaces less than {displaced_volume!r} m^3"
        )
    if fa > 0.0:  # pragma: no cover - the lowest waterline is below the keel
        raise TrimError("the lowest waterline already displaces too much")
    if initial_height is not None and math.isfinite(initial_height) and a < initial_height < b:
        fm = _residual(float(initial_height))
        if fm < 0.0:
            a, fa = float(initial_height), fm
        else:
            b, fb = float(initial_height), fm
    side = 0
    for _ in range(max_iterations):
        candidate = _false_position(a, fa, b, fb)
        if not (min(a, b) < candidate < max(a, b)):
            candidate = 0.5 * (a + b)
            side = 0
        fc = _residual(candidate)
        if abs(fc) <= volume_tolerance * displaced_volume:
            return candidate
        if fc < 0.0:
            a, fa = candidate, fc
            if side == -1:
                fb *= 0.5  # Illinois: keep the bracket from stalling on one side
            side = -1
        else:
            b, fb = candidate, fc
            if side == 1:
                fa *= 0.5
            side = 1
        if abs(b - a) <= tolerance * max(1.0, abs(b)):
            break
    return 0.5 * (a + b)


@dataclass(frozen=True)
class TrimSolution:
    """The equilibrium reached when ``weight`` is shifted ``shift`` along the hull.

    ``lcb_target`` is the LCB that conserves the longitudinal moment,
    ``lcb_0 + (weight / displacement) * shift``; ``lcb`` is the LCB the trimmed
    hull actually has there. ``small_angle_deg`` is the closed-form answer
    ``atan(weight * shift / (rho * I_L))``, reported beside the solved angle so the
    caller can see how much of the answer is the discretisation and how much is the
    small-angle assumption.
    """

    weight: float
    shift: float
    displacement_mass: float
    lcb_0: float
    lcb_target: float
    lcb: float
    trim_angle_rad: float
    trim_angle_deg: float
    small_angle_deg: float
    change_of_trim_cm: float
    mct: float
    waterline_at_reference: float
    displaced_volume: float
    iterations: int
    residual_m: float
    basis: str = "LCB conservation on a per-station tilted waterline"
    schema: str = SCHEMA_TRIM_SOLUTION

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": self.schema,
            "weight_kg": _round(self.weight, 6),
            "shift_m": _round(self.shift),
            "displacement_mass_kg": _round(self.displacement_mass, 3),
            "lcb_0_m": _round(self.lcb_0),
            "lcb_target_m": _round(self.lcb_target),
            "lcb_m": _round(self.lcb),
            "trim_angle_rad": _round(self.trim_angle_rad),
            "trim_angle_deg": _round(self.trim_angle_deg, 9),
            "small_angle_deg": _round(self.small_angle_deg, 9),
            "change_of_trim_cm": _round(self.change_of_trim_cm),
            "mct_t_m_per_cm": _round(self.mct),
            "waterline_at_reference_m": _round(self.waterline_at_reference),
            "displaced_volume_m3": _round(self.displaced_volume),
            "iterations": self.iterations,
            "residual_m": _round(self.residual_m),
            "basis": self.basis,
        }


def small_angle_trim(
    particulars: HydrostaticParticulars, weight: float, shift: float
) -> tuple[float, float, float]:
    """Closed-form trim from a weight shift: ``(tan_trim, angle_deg, change_of_trim_cm)``.

    The trimming moment is ``weight * shift`` and the ship answers it with
    ``rho * g * I_L * tan(theta)``, so ``tan(theta) = weight * shift / (rho * I_L)``
    -- the displacement cancels, which is why a weight shift trims a light ship and
    a loaded one by the same angle as long as both keep the same waterplane. The
    same answer through the loading-computer route is ``w * d / MCT`` centimetres,
    and both are reported: they must agree, and
    ``test_the_two_centimetre_routes_agree`` holds them to it.
    """
    if not math.isfinite(weight) or weight <= 0.0:
        raise TrimError(f"weight must be positive, got {weight!r}")
    if not math.isfinite(shift):
        raise TrimError(f"shift must be finite, got {shift!r}")
    if particulars.longitudinal_second_moment <= 0.0:
        raise TrimError("the waterplane has no longitudinal second moment to trim about")
    tan_trim = weight * shift / (particulars.water_density * particulars.longitudinal_second_moment)
    angle = math.atan(tan_trim)
    change_cm = 100.0 * particulars.length * tan_trim
    return tan_trim, math.degrees(angle), change_cm


def equilibrium_trim(
    mesh: Mesh,
    *,
    weight: float,
    shift: float,
    waterline_z: float = 0.0,
    water_density: float = SEAWATER_DENSITY,
    station_count: int = 41,
    angle_tolerance: float = 1e-14,
    lcb_tolerance: float = 1e-9,
    max_trim_deg: float = 15.0,
    max_iterations: int = 120,
) -> TrimSolution:
    """Solve the trim angle that conserves longitudinal moment, and report it.

    The hull floats at ``waterline_z`` upright; ``weight`` kilograms then move
    ``shift`` metres toward ``+x``. Longitudinally the ship is in equilibrium when
    its LCB has moved to ``lcb_0 + (weight / displacement) * shift``, because the
    centre of gravity has moved by ``weight * shift / displacement`` and buoyancy
    must sit under it. So the solve is: bisect the trim angle, at every angle hold
    the displaced volume fixed (that is :func:`waterline_for_displacement`, drafted
    from the previous angle so it converges in a few iterations), and read off the
    LCB. Both are strictly monotone in the angle, so the root is unique; it is found
    with a safeguarded secant (Illinois) on the tangent, which lands in about ten
    iterations because the LCB is smooth in the angle. The answer is exact to the
    mesh's own discretisation -- there is no linearisation anywhere in it, which is
    what lets a box hull reproduce the wallsided closed form at 8 degrees of trim
    rather than merely near zero.
    """
    upright = hydrostatic_particulars(
        mesh, waterline_z, water_density=water_density, station_count=station_count
    )
    if not math.isfinite(weight) or weight < 0.0:
        raise TrimError(f"weight must be non-negative, got {weight!r}")
    if not math.isfinite(shift):
        raise TrimError(f"shift must be finite, got {shift!r}")
    lcb_target = upright.lcb + (weight / upright.displacement_mass) * shift
    x_ref = _reference_x(mesh)
    stations = _stations(mesh, station_count)
    cutter = _SectionCutter(mesh, stations)
    tangent_max = math.tan(math.radians(max_trim_deg))
    warm: dict[str, float] = {"height": upright.waterline_z}

    def _residual(tan_trim: float) -> tuple[float, float]:
        """LCB minus its target at this trim, and the waterline that holds the volume."""
        height = waterline_for_displacement(
            mesh,
            trim_angle=math.atan(tan_trim),
            displaced_volume=upright.volume,
            reference_x=x_ref,
            station_count=station_count,
            initial_height=warm["height"],
            cutter=cutter,
        )
        warm["height"] = height
        _, lcb = _trimmed_volume_and_lcb(
            cutter,
            stations,
            tan_trim=tan_trim,
            waterline_at_reference=height,
            reference_x=x_ref,
        )
        return lcb - lcb_target, height

    sign = 1.0 if lcb_target >= upright.lcb else -1.0
    a, fa = 0.0, upright.lcb - lcb_target
    if abs(fa) <= lcb_tolerance:
        tan_trim = 0.0
        height = upright.waterline_z
        iterations = 1
    else:
        b = sign * tangent_max
        fb = _residual(b)[0]
        if fb * sign < 0.0:  # pragma: no cover - geometry cannot deliver the shift
            raise TrimError(
                f"no equilibrium within {max_trim_deg!r} degrees of trim: the hull displaces "
                f"{upright.volume!r} m^3 at this waterline and its sections above it cannot "
                f"carry an LCB shift of {lcb_target - upright.lcb!r} m"
            )
        a, fa = 0.0, upright.lcb - lcb_target
        iterations = 1
        height = upright.waterline_z
        side = 0
        tan_trim = b
        for _ in range(max_iterations):
            candidate = _false_position(a, fa, b, fb)
            if not (min(a, b) < candidate < max(a, b)):
                candidate = 0.5 * (a + b)
                side = 0
            fc, height_c = _residual(candidate)
            iterations += 1
            tan_trim, height = candidate, height_c
            if abs(fc) <= lcb_tolerance:
                break
            if fc * fa < 0.0:
                b, fb = candidate, fc
                if side == -1:
                    fa *= 0.5
                side = -1
            else:
                a, fa = candidate, fc
                if side == 1:
                    fb *= 0.5
                side = 1
            if abs(b - a) <= angle_tolerance * max(1.0, abs(b)):
                break
    height = warm["height"]
    _, lcb = _trimmed_volume_and_lcb(
        cutter,
        stations,
        tan_trim=tan_trim,
        waterline_at_reference=height,
        reference_x=x_ref,
    )
    trim_angle = math.atan(tan_trim)
    _, small_deg, change_cm = small_angle_trim(upright, weight, shift)
    return TrimSolution(
        weight=weight,
        shift=shift,
        displacement_mass=upright.displacement_mass,
        lcb_0=upright.lcb,
        lcb_target=lcb_target,
        lcb=lcb,
        trim_angle_rad=trim_angle,
        trim_angle_deg=math.degrees(trim_angle),
        small_angle_deg=small_deg,
        change_of_trim_cm=change_cm,
        mct=upright.mct,
        waterline_at_reference=height,
        displaced_volume=upright.volume,
        iterations=iterations,
        residual_m=lcb - lcb_target,
    )


@dataclass(frozen=True)
class BoxParticulars:
    """Closed-form particulars of a box-shaped hull, for holding the mesh route to account.

    Units follow :class:`HydrostaticParticulars`: metres, m^3, kg, t/cm, t*m/cm.
    Keel at ``-draft``, deck above the waterline, ``length`` along ``x`` and
    ``beam`` along ``y``, centred so that ``LCF = LCB = L/2``.
    """

    length: float
    beam: float
    draft: float
    water_density: float
    volume: float
    displacement_mass: float
    lcb: float
    lcf: float
    kb: float
    waterplane_area: float
    transverse_second_moment: float
    longitudinal_second_moment: float
    bm_transverse: float
    bm_longitudinal: float
    km_transverse: float
    km_longitudinal: float
    tpc: float
    mct: float
    basis: str = "closed forms for a rectangular hull"
    schema: str = "tensorcad.box_particulars/1"

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": self.schema,
            "length_m": _round(self.length),
            "beam_m": _round(self.beam),
            "draft_m": _round(self.draft),
            "displacement_volume_m3": _round(self.volume),
            "displacement_mass_kg": _round(self.displacement_mass, 3),
            "lcb_m": _round(self.lcb),
            "lcf_m": _round(self.lcf),
            "kb_m": _round(self.kb),
            "waterplane_area_m2": _round(self.waterplane_area),
            "transverse_second_moment_m4": _round(self.transverse_second_moment),
            "longitudinal_second_moment_m4": _round(self.longitudinal_second_moment),
            "bm_transverse_m": _round(self.bm_transverse),
            "bm_longitudinal_m": _round(self.bm_longitudinal),
            "km_transverse_m": _round(self.km_transverse),
            "km_longitudinal_m": _round(self.km_longitudinal),
            "tpc_t_per_cm": _round(self.tpc),
            "mct_t_m_per_cm": _round(self.mct),
            "water_density_kg_per_m3": _round(self.water_density, 3),
            "basis": self.basis,
        }


def box_analytic_particulars(
    length: float,
    beam: float,
    draft: float,
    *,
    water_density: float = SEAWATER_DENSITY,
) -> BoxParticulars:
    """Exact particulars of a box-shaped hull: ``V = LBT``, ``KB = T/2``, ``I_L = B L^3/12``.

    Every quantity is an elementary integral of a rectangle, so this is a genuine
    reference and not a fit: it is what the mesh route has to reproduce, and for a
    box it should reproduce it to round-off, because the station area is linear in
    ``x`` and Simpson's rule is exact on a quadratic.
    """
    for name, value in (("length", length), ("beam", beam), ("draft", draft)):
        if not math.isfinite(value) or value <= 0.0:
            raise TrimError(f"{name} must be positive, got {value!r}")
    if not math.isfinite(water_density) or water_density <= 0.0:
        raise TrimError(f"water_density must be positive, got {water_density!r}")
    volume = length * beam * draft
    area = length * beam
    i_t = length * beam**3 / 12.0
    i_l = beam * length**3 / 12.0
    bm_t = i_t / volume
    bm_l = i_l / volume
    kb = 0.5 * draft
    density_t = water_density / 1000.0
    return BoxParticulars(
        length=length,
        beam=beam,
        draft=draft,
        water_density=water_density,
        volume=volume,
        displacement_mass=water_density * volume,
        lcb=0.5 * length,
        lcf=0.5 * length,
        kb=kb,
        waterplane_area=area,
        transverse_second_moment=i_t,
        longitudinal_second_moment=i_l,
        bm_transverse=bm_t,
        bm_longitudinal=bm_l,
        km_transverse=kb + bm_t,
        km_longitudinal=kb + bm_l,
        tpc=density_t * area / 100.0,
        mct=density_t * i_l / (100.0 * length),
    )


def box_trim_tangent(
    length: float, beam: float, weight: float, shift: float, *, water_density: float = SEAWATER_DENSITY
) -> float:
    """Exact ``tan(trim)`` for a box after a weight shift -- valid at *any* angle.

    The waterplane of a rectangular hull does not change shape as it trims, so the
    wallsided relation ``tan(theta) = w d / (rho I_L)`` is exact rather than
    linearised, and the trimmed LCB is ``tan(theta) * L^2 / 12`` about midship with
    no small-angle term. That makes the box the one case where the mesh solve can
    be checked against a closed form without a "for small angles" caveat, and
    ``docs/trim.md`` uses angles up to 8 degrees to show it.
    """
    for name, value in (("length", length), ("beam", beam), ("weight", weight)):
        if not math.isfinite(value) or value <= 0.0:
            raise TrimError(f"{name} must be positive, got {value!r}")
    if not math.isfinite(shift):
        raise TrimError(f"shift must be finite, got {shift!r}")
    if not math.isfinite(water_density) or water_density <= 0.0:
        raise TrimError(f"water_density must be positive, got {water_density!r}")
    i_l = beam * length**3 / 12.0
    return weight * shift / (water_density * i_l)


__all__ = [
    "GRAVITY",
    "SCHEMA_PARTICULARS",
    "SCHEMA_TRIMMED",
    "SCHEMA_TRIM_SOLUTION",
    "SEAWATER_DENSITY",
    "BoxParticulars",
    "HydrostaticParticulars",
    "TrimError",
    "TrimSolution",
    "TrimmedHydrostatics",
    "box_analytic_particulars",
    "box_trim_tangent",
    "equilibrium_trim",
    "hydrostatic_particulars",
    "hydrostatic_table",
    "small_angle_trim",
    "trimmed_hydrostatics",
    "waterline_for_displacement",
]
