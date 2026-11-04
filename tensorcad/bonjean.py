"""Bonjean curves: a station's immersed area as a function of the waterline.

Everything else in this package answers a *whole-hull* question. The hydrostatics
module integrates over the wetted surface and reports one displaced volume, one
waterplane, one centre; :func:`tensorcad.hydrostatics.section_area_at` returns the
area of the station the mesh happens to have a vertex ring on, at the full depth
of that ring. Nothing in the tree can answer the question a naval architect asks
next: *at station x, how much area is under the waterline z_w?* That table --
area against waterline, one curve per station -- is the Bonjean curve, and it is
the input to damage-stability (lost buoyancy of one flooded compartment) and to
still-water shear force and bending moment (the load distribution along the hull
is the difference between the weight curve and the buoyancy curve ``gamma * A(x)``
at the design waterline).

So this module adds:

* :func:`section_loops` -- the outline(s) where the plane ``x = station`` cuts a
  closed mesh, in ``(y, z)``;
* :func:`section_properties` -- area, perimeter, centroid and extent of that
  outline below a chosen waterline, by clipping the polygon and taking its area;
* :func:`bonjean_curve` / :func:`bonjean_table` -- the same, repeated over
  waterlines and stations, plus :meth:`BonjeanTable.integrated_volume`, which
  recovers the displacement by integrating the areas along the hull -- a check
  that uses none of the surface integrals :mod:`tensorcad.hydrostatics` uses;
* :func:`analytic_section_area` -- the closed form for the parametric hull
  family, ``(2/3) * beam * sin(pi x / L)^0.7 * zeta^(3/2) / sqrt(draft)``, which
  the sampled mesh never sees.

The cut is taken at ``station + 1e-9 * scale`` rather than exactly on the station,
because the hull loft puts vertex rings exactly on the station planes: a plane
through a vertex is a degenerate cut (a face can have two vertices on it), and
nudging the plane by a millionth of a millimetre converts the degeneracy into an
ordinary crossing whose area differs from the exact station by ``O(1e-9)`` of the
station spacing. That is smaller than the discretisation error of any mesh this
package produces by nine orders of magnitude, and it is stated here rather than
hidden.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from .geometry import Mesh
from .hull import BREADTH_EXPONENT

#: Relative offset applied to the cutting plane so that no vertex lies on it.
STATION_NUDGE = 1e-9


def _nudge(mesh: Mesh, station: float) -> float:
    """Shift the cutting plane off any vertex ring, staying inside the mesh.

    The offset is ``1e-9`` of the mesh scale. At the ends of a body the inward
    direction is the one that exists: nudging ``x = length`` outwards would cut
    empty space and report a zero section for the stem, so the sign flips there
    and the end station keeps its section.
    """
    box = mesh.bounding_box()
    scale = max(abs(value) for value in box.maximum) or 1.0
    offset = STATION_NUDGE * scale
    if station + offset > box.maximum[0]:
        return station - offset
    return station + offset


def section_loops(mesh: Mesh, station: float) -> tuple[tuple[tuple[float, float], ...], ...]:
    """Cut a closed mesh with the plane ``x = station``; return ``(y, z)`` loops.

    A closed triangle mesh is cut into closed polylines: every triangle that
    straddles the plane contributes exactly one segment, and the segments chain
    end to end because the mesh is a manifold. Vertices exactly on the plane are
    avoided by :data:`STATION_NUDGE`, so each crossing is an ordinary two-point
    intersection and the segments can be matched by their shared edge key --
    identical floats, not "close enough" floats, which is why the chaining cannot
    drift apart on a 200 k-triangle mesh.

    Returns the loops without a repeated closing point. A station outside the
    mesh, or on a stem, returns ``()``.
    """
    if not math.isfinite(station):
        raise ValueError(f"station must be finite, got {station!r}")
    x0 = _nudge(mesh, station)
    vertices = mesh.vertices

    points: dict[tuple[int, int], tuple[float, float]] = {}
    segments: list[tuple[tuple[int, int], tuple[int, int]]] = []

    for face in mesh.faces:
        crossings: list[tuple[int, int]] = []
        for index in range(3):
            a = face[index]
            b = face[(index + 1) % 3]
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
        if len(crossings) == 2:
            segments.append((crossings[0], crossings[1]))
        elif len(crossings) > 2:  # pragma: no cover - impossfor a triangle
            raise ValueError("a triangle crossed the station plane more than twice")

    if not segments:
        return ()

    adjacency: dict[tuple[int, int], list[tuple[int, int]]] = {}
    for first, second in segments:
        adjacency.setdefault(first, []).append(second)
        adjacency.setdefault(second, []).append(first)

    unused = {tuple(sorted(segment)) for segment in segments}
    loops: list[tuple[tuple[float, float], ...]] = []
    while unused:
        first, second = next(iter(unused))
        unused.discard((first, second))
        chain = [first, second]
        previous, current = first, second
        while current != chain[0]:
            following = None
            for candidate in adjacency[current]:
                if candidate == previous:
                    continue
                if tuple(sorted((current, candidate))) in unused:
                    following = candidate
                    break
            if following is None:
                break  # open chain: a non-manifold mesh, reported as an open loop
            unused.discard(tuple(sorted((current, following))))
            chain.append(following)
            previous, current = current, following
        # Drop the repeated closing vertex if the walk came back to the start.
        if len(chain) > 2 and chain[-1] == chain[0]:
            chain = chain[:-1]
        loops.append(tuple(points[key] for key in chain))
    return tuple(loops)


def _clip_below(
    polygon: tuple[tuple[float, float], ...], waterline_z: float
) -> list[tuple[float, float]]:
    """Sutherland-Hodgman clip of a polygon against the half plane ``z <= waterline_z``."""
    if not polygon:
        return []
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


def _signed_ratio(value: float, reference: float) -> float:
    """``value / reference`` with the convention that a zero reference is NaN."""
    if reference == 0.0 or not math.isfinite(reference):
        return float("nan")
    return value / reference


def _signed_area(polygon: list[tuple[float, float]]) -> float:
    total = 0.0
    count = len(polygon)
    for index in range(count):
        y0, z0 = polygon[index]
        y1, z1 = polygon[(index + 1) % count]
        total += y0 * z1 - y1 * z0
    return 0.5 * total


def _perimeter(polygon: tuple[tuple[float, float], ...]) -> float:
    count = len(polygon)
    if count < 2:
        return 0.0
    total = 0.0
    for index in range(count):
        y0, z0 = polygon[index]
        y1, z1 = polygon[(index + 1) % count]
        total += math.hypot(y1 - y0, z1 - z0)
    return total


@dataclass(frozen=True)
class SectionProperties:
    """One station's immersed cross-section: area, perimeter, centroid, extent."""

    station: float
    waterline_z: float
    area: float
    perimeter: float
    centroid_y: float
    centroid_z: float
    minimum_z: float
    maximum_z: float
    loops: int
    basis: str = "polygon clipped at the waterline"

    def to_dict(self) -> dict[str, float | int | str]:
        return {
            "station_m": _round(self.station),
            "waterline_z_m": _round(self.waterline_z),
            "area_m2": _round(self.area),
            "perimeter_m": _round(self.perimeter),
            "centroid_y_m": _round(self.centroid_y),
            "centroid_z_m": _round(self.centroid_z),
            "minimum_z_m": _round(self.minimum_z),
            "maximum_z_m": _round(self.maximum_z),
            "loops": self.loops,
            "basis": self.basis,
        }


def _round(value: float, digits: int = 9) -> float:
    if not isinstance(value, float) or not math.isfinite(value):
        return value
    return round(value, digits)


EMPTY_SECTION = SectionProperties(
    station=float("nan"),
    waterline_z=float("nan"),
    area=0.0,
    perimeter=0.0,
    centroid_y=float("nan"),
    centroid_z=float("nan"),
    minimum_z=float("nan"),
    maximum_z=float("nan"),
    loops=0,
    basis="no material at this station",
)


def section_properties(
    mesh: Mesh, station: float, waterline_z: float | None = None
) -> SectionProperties:
    """Area and geometry of the section at ``station`` below ``waterline_z``.

    ``waterline_z`` defaults to the top of the section, so the call returns the
    whole section -- the quantity :func:`tensorcad.hydrostatics.section_area_at`
    reports from the vertex profile, reached here by clipping a cut polygon
    instead, which lets the waterline move.
    """
    loops = section_loops(mesh, station)
    if not loops:
        return EMPTY_SECTION
    tops = [z for loop in loops for _, z in loop]
    if waterline_z is None:
        waterline_z = max(tops)
    if not math.isfinite(waterline_z):
        raise ValueError(f"waterline_z must be finite, got {waterline_z!r}")

    weighted_area = 0.0
    weighted_y = 0.0
    weighted_z = 0.0
    signed_total = 0.0
    perimeter = 0.0
    minimum_z = math.inf
    maximum_z = -math.inf
    kept = 0
    for loop in loops:
        clipped = _clip_below(loop, waterline_z)
        if len(clipped) < 3:
            continue
        signed = _signed_area(clipped)
        if signed == 0.0:
            continue
        # First moment of the clipped polygon, exactly, by the same shoelace
        # pairing that gives its area: for a polygon with vertices (y_i, z_i),
        # INT y dA = (1/6) SUM (y_i + y_{i+1}) (y_i z_{i+1} - y_{i+1} z_i), and the
        # same with y and z exchanged. The traversal sign cancels between the
        # moment and the area, so the centroid is (moment / 6) / signed area.
        sum_y = 0.0
        sum_z = 0.0
        for index in range(len(clipped)):
            y0, z0 = clipped[index]
            y1, z1 = clipped[(index + 1) % len(clipped)]
            cross = y0 * z1 - y1 * z0
            sum_y += cross * (y0 + y1)
            sum_z += cross * (z0 + z1)
        weighted_area += abs(signed)
        signed_total += signed
        weighted_y += sum_y / 6.0
        weighted_z += sum_z / 6.0
        perimeter += _perimeter(loop)
        minimum_z = min(minimum_z, min(z for _, z in loop))
        maximum_z = max(maximum_z, max(z for _, z in loop))
        kept += 1

    if kept == 0:
        return SectionProperties(
            station=station,
            waterline_z=waterline_z,
            area=0.0,
            perimeter=perimeter,
            centroid_y=float("nan"),
            centroid_z=float("nan"),
            minimum_z=minimum_z if math.isfinite(minimum_z) else float("nan"),
            maximum_z=maximum_z if math.isfinite(maximum_z) else float("nan"),
            loops=len(loops),
            basis="waterline below the section",
        )

    return SectionProperties(
        station=station,
        waterline_z=waterline_z,
        area=weighted_area,
        perimeter=perimeter,
        centroid_y=_signed_ratio(weighted_y, signed_total),
        centroid_z=_signed_ratio(weighted_z, signed_total),
        minimum_z=minimum_z,
        maximum_z=maximum_z,
        loops=len(loops),
    )


@dataclass(frozen=True)
class BonjeanCurve:
    """One station's area against waterline: ``areas[i]`` at ``waterlines[i]``."""

    station: float
    waterlines: tuple[float, ...]
    areas: tuple[float, ...]
    section: SectionProperties

    def area_at(self, waterline_z: float) -> float:
        """Linear interpolation between the tabulated waterlines."""
        if not self.waterlines:
            return 0.0
        if waterline_z <= self.waterlines[0]:
            return self.areas[0] if waterline_z == self.waterlines[0] else 0.0
        if waterline_z >= self.waterlines[-1]:
            return self.areas[-1]
        for index in range(1, len(self.waterlines)):
            upper = self.waterlines[index]
            if waterline_z <= upper:
                lower = self.waterlines[index - 1]
                span = upper - lower
                fraction = 0.0 if span == 0.0 else (waterline_z - lower) / span
                return self.areas[index - 1] + fraction * (self.areas[index] - self.areas[index - 1])
        return self.areas[-1]

    def to_dict(self) -> dict[str, object]:
        return {
            "station_m": _round(self.station),
            "waterlines_m": [_round(value) for value in self.waterlines],
            "areas_m2": [_round(value) for value in self.areas],
            "section": self.section.to_dict(),
        }


def default_waterlines(section: SectionProperties, count: int = 13) -> tuple[float, ...]:
    """Evenly spaced waterlines from the bottom of the section to its top."""
    if count < 2:
        raise ValueError("count must be at least 2")
    bottom = section.minimum_z
    top = section.maximum_z
    if not math.isfinite(bottom) or not math.isfinite(top):
        return ()
    return tuple(bottom + (top - bottom) * index / (count - 1) for index in range(count))


def bonjean_curve(
    mesh: Mesh,
    station: float,
    *,
    waterlines: tuple[float, ...] | None = None,
    count: int = 13,
) -> BonjeanCurve:
    """The Bonjean curve at one station, sampled on ``waterlines`` (or ``count``)."""
    whole = section_properties(mesh, station)
    if not math.isfinite(whole.maximum_z):
        return BonjeanCurve(station=station, waterlines=(), areas=(), section=whole)
    if waterlines is None:
        waterlines = default_waterlines(whole, count)
    ordered = tuple(sorted(float(value) for value in waterlines))
    areas = tuple(section_properties(mesh, station, waterline).area for waterline in ordered)
    return BonjeanCurve(station=station, waterlines=ordered, areas=areas, section=whole)


@dataclass(frozen=True)
class BonjeanTable:
    """A full Bonjean table: stations down the rows, waterlines across the columns."""

    stations: tuple[float, ...]
    waterlines: tuple[float, ...]
    areas: tuple[tuple[float, ...], ...]
    curves: tuple[BonjeanCurve, ...] = ()

    def area_at(self, station: float, waterline_z: float) -> float:
        """Interpolate the table in both directions (bilinear on the grid)."""
        if not self.stations or not self.waterlines:
            return 0.0
        rows = [curve.area_at(waterline_z) for curve in self.curves]
        if station <= self.stations[0]:
            return rows[0]
        if station >= self.stations[-1]:
            return rows[-1]
        for index in range(1, len(self.stations)):
            if station <= self.stations[index]:
                lower = self.stations[index - 1]
                span = self.stations[index] - lower
                fraction = 0.0 if span == 0.0 else (station - lower) / span
                return rows[index - 1] + fraction * (rows[index] - rows[index - 1])
        return rows[-1]

    def integrated_volume(
        self, waterline_z: float | None = None, *, method: str = "simpson"
    ) -> float:
        """Displacement recovered by integrating the station areas along ``x``.

        ``(1/3) * h * (f_0 + 4 f_1 + 2 f_2 + ... + f_n)`` when the stations are
        evenly spaced and there is an odd number of them, otherwise the trapezoid
        rule. Both are stated rather than silently switching.
        """
        if not self.stations:
            return 0.0
        if waterline_z is None:
            waterline_z = self.waterlines[-1] if self.waterlines else 0.0
        values = [curve.area_at(waterline_z) for curve in self.curves]
        return _integrate(self.stations, values, method=method)

    def to_dict(self) -> dict[str, object]:
        return {
            "stations_m": [_round(value) for value in self.stations],
            "waterlines_m": [_round(value) for value in self.waterlines],
            "areas_m2": [[_round(value) for value in row] for row in self.areas],
        }


def _integrate(stations: tuple[float, ...], values: list[float], *, method: str) -> float:
    count = len(values)
    if count < 2:
        return 0.0
    spacing = (stations[-1] - stations[0]) / (count - 1)
    uniform = all(
        abs((stations[index + 1] - stations[index]) - spacing) <= 1e-9 * max(1.0, abs(spacing))
        for index in range(count - 1)
    )
    if method not in {"simpson", "trapezoid"}:
        raise ValueError(f"unknown integration method {method!r}")
    if method == "simpson" and uniform and count % 2 == 1:
        total = values[0] + values[-1]
        for index in range(1, count - 1):
            total += (4.0 if index % 2 == 1 else 2.0) * values[index]
        return total * spacing / 3.0
    total = 0.0
    for index in range(count - 1):
        total += 0.5 * (values[index] + values[index + 1]) * (stations[index + 1] - stations[index])
    return total


def bonjean_table(
    mesh: Mesh,
    *,
    stations: tuple[float, ...] | None = None,
    waterlines: tuple[float, ...] | None = None,
    station_count: int = 41,
    waterline_count: int = 13,
) -> BonjeanTable:
    """The Bonjean table over ``stations`` and ``waterlines``.

    Defaults follow the hull generator: the stations are the design ones (evenly
    spaced along the bounding box in ``x``) and the waterlines run from the keel
    to the top of the sections.
    """
    box = mesh.bounding_box()
    if stations is None:
        if station_count < 2:
            raise ValueError("station_count must be at least 2")
        length = box.sizes[0]
        x_min = box.minimum[0]
        stations = tuple(x_min + length * index / (station_count - 1) for index in range(station_count))
    stations = tuple(sorted(float(value) for value in stations))

    curves: list[BonjeanCurve] = []
    for station in stations:
        curves.append(bonjean_curve(mesh, station, waterlines=waterlines, count=waterline_count))

    if waterlines is None:
        top = max((curve.waterlines[-1] for curve in curves if curve.waterlines), default=float("nan"))
        bottom = min((curve.waterlines[0] for curve in curves if curve.waterlines), default=float("nan"))
        if math.isfinite(top) and math.isfinite(bottom):
            waterlines = tuple(
                bottom + (top - bottom) * index / (waterline_count - 1)
                for index in range(waterline_count)
            )
        else:
            waterlines = ()
    waterlines = tuple(sorted(float(value) for value in waterlines))

    areas: list[tuple[float, ...]] = []
    for station in stations:
        row = tuple(section_properties(mesh, station, waterline).area for waterline in waterlines)
        areas.append(row)

    return BonjeanTable(
        stations=stations,
        waterlines=waterlines,
        areas=tuple(areas),
        curves=tuple(curves),
    )


def taper_exponent_at(
    length: float,
    station: float,
    fore_exponent: float = BREADTH_EXPONENT,
    aft_exponent: float = BREADTH_EXPONENT,
) -> float:
    """Taper exponent the generator applies at ``station``.

    The split family (``HullParameters.fore_exponent`` / ``aft_exponent``) indexes
    the exponent by station *number*, ``2 * i >= nx - 1`` meaning fore, which is
    exactly the half-open split ``x >= length / 2`` once ``x = L * i / (nx - 1)``.
    Midship therefore takes the fore exponent, and ``sin(pi/2) = 1`` makes the
    choice immaterial there. For the shipped family both arguments default to
    ``BREADTH_EXPONENT`` and the answer is that constant.
    """
    for name, value in (("length", length), ("fore_exponent", fore_exponent), ("aft_exponent", aft_exponent)):
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be a finite positive number, got {value!r}")
    return fore_exponent if 2.0 * station >= length else aft_exponent


def analytic_section_area(
    length: float,
    beam: float,
    draft: float,
    station: float,
    waterline_z: float,
    *,
    fore_exponent: float = BREADTH_EXPONENT,
    aft_exponent: float = BREADTH_EXPONENT,
) -> float:
    """Closed form for the parametric hull family: area under ``waterline_z``.

    The generator's half breadth is ``(beam/2) * sin(pi x / L)^p * sqrt(zeta/draft)``
    for ``zeta = z + draft`` measured up from the keel, so the section area is the
    integral of ``2 * half_breadth`` over the immersed depth::

        A(x, zeta) = (2/3) * beam * sin(pi x / L)^p * zeta^(3/2) / sqrt(draft)

    At midship and ``zeta = draft`` this is ``(2/3) * beam * draft``, the
    :func:`tensorcad.hull.analytic_maximum_section_area` the hydrostatics module
    already compares against -- but the mesh-facing function here varies with the
    station and with the waterline, which is what a Bonjean curve needs and what
    nothing in the tree had.

    ``fore_exponent`` / ``aft_exponent`` are the per-half sine exponents of the
    asymmetric family (:mod:`tensorcad.hull` 0.17.0). Both default to
    ``BREADTH_EXPONENT``, so the shipped surface keeps this reference bit for bit;
    passing the two exponents the mesh was generated with extends the same closed
    form to the asymmetric family through :func:`taper_exponent_at`. The vertical
    taper is ``sqrt`` for every member of this family, which is why the exponent
    split does not touch the ``zeta^(3/2)`` factor.
    """
    for name, value in (("length", length), ("beam", beam), ("draft", draft)):
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be a finite positive number, got {value!r}")
    if not 0.0 <= station <= length:
        raise ValueError(f"station must lie in [0, {length}], got {station!r}")
    zeta = waterline_z + draft
    if not math.isfinite(zeta):
        raise ValueError(f"waterline_z must be finite, got {waterline_z!r}")
    if zeta <= 0.0:
        return 0.0
    exponent = taper_exponent_at(length, station, fore_exponent, aft_exponent)
    taper = math.pow(math.sin(math.pi * station / length), exponent)
    return (2.0 / 3.0) * beam * taper * math.pow(zeta, 1.5) / math.sqrt(draft)


__all__ = [
    "EMPTY_SECTION",
    "STATION_NUDGE",
    "BonjeanCurve",
    "BonjeanTable",
    "SectionProperties",
    "analytic_section_area",
    "bonjean_curve",
    "bonjean_table",
    "default_waterlines",
    "section_loops",
    "section_properties",
    "taper_exponent_at",
]
