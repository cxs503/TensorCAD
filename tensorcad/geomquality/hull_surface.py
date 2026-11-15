r"""Physical validation of an imported hull surface against published particulars.

Why this module exists
---------------------
The pipeline of :mod:`tensorcad.geomquality.trim_mesh` ends in a *topology* gate
(watertight / manifold / single shell).  For the publisher files that actually
ship hull geometry -- ``KCS_hull_SVA.igs`` (55 patches, entity census
``{128: 55}``) and ``kvlcc2.igs`` -- that gate can never pass: they are
half-breadth **open** surface models.  Measured on the KCS file: ``y in [0, 16.11]``
(= B/2), no centreplane cap, no deck cap, and no trimmed entities at all, so the
only honest verdict for ``watertight`` is *false*, no matter how the tessellator
behaves.

A green topology gate would not prove the geometry right anyway.  What those
workshop files *are* is a correct description of the moulded hull, and that is
checkable without any closure: the outermost breadth of the surface at a station
``x`` and height ``z`` is the moulded half-breadth, so sectional areas, the
displacement volume, ``Cb``/``Cm`` and the wetted area all follow from the
surface as it stands.  Vendors publish exactly those numbers (KCS: nabla =
52 030 m3, S = 9 424 m2; KVLCC2: nabla = 312 622 m3, S = 27 194 m2), which makes
them an *external* benchmark instead of the self-consistent round trips that
earlier rounds of the IGES anchor could only offer.

Method
------
Pure standard library (the package charter keeps ``dependencies = []``).

* **Envelope breadth.**  A ray along ``+-y`` at ``(x, z)`` is intersected with
  every triangle whose ``x``-span brackets ``x`` (x-slabs keep this cheap).  Since
  the ray is axis aligned, each intersection reduces to a point-in-triangle test
  in the ``x-z`` plane plus linear interpolation of ``y``; the largest ``|y|``
  over all hits is the moulded breadth.  Overlapping patches are harmless -- the
  envelope takes the outermost surface, it does not sum the soup.
* **Sectional area** ``A(x) = int_0^T b(x, z) dz`` by the trapezoid rule on a
  uniform ``z`` grid; **volume** ``nabla = k * int A(x) dx`` with ``k = 2`` for a
  half-breadth model and ``k = 1`` for a full one (detected from the model).
* **Wetted area**: triangle area below the waterline, clipped by Sutherland-
  Hodgman at ``z = T``, doubled for a half-breadth model.

The result carries its own discretisation knob (``stations``/``levels``) so the
convergence of the integrals can be shown rather than assumed.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Iterable, Sequence

Vec3 = tuple[float, float, float]

#: Published particulars of the two workshop hulls this repository can import.
#: ``source`` is the workshop/vendor documentation the numbers are quoted from;
#: ``volume`` is the moulded displacement volume at the design draft, ``s_wetted``
#: the wetted surface at the same draft (full hull, both sides).
PARTICULARS: dict[str, dict[str, Any]] = {
    "KCS": {
        "lpp": 230.0,
        "beam": 32.2,
        "draft": 10.8,
        "volume": 52030.0,
        "s_wetted": 9424.0,
        "cb": 0.651,
        "source": "KRISO container ship (KCS) hull-form particulars as distributed "
                  "with the SIMMAN / Tokyo-2015 workshop geometry (Lpp 230.0 m, "
                  "B 32.2 m, T 10.8 m, Cb 0.651, nabla 52 030 m3, S 9 424 m2)",
    },
    "KVLCC2": {
        "lpp": 320.0,
        "beam": 58.0,
        "draft": 20.8,
        "volume": 312622.0,
        "s_wetted": 27194.0,
        "cb": 0.8098,
        "source": "KRISO tanker (KVLCC2) particulars as distributed with the "
                  "SIMMAN / Tokyo-2015 workshop geometry (Lpp 320.0 m, B 58.0 m, "
                  "T 20.8 m, Cb 0.8098, nabla 312 622 m3, S 27 194 m2)",
    },
}


# --------------------------------------------------------------------------- #
# geometry primitives
# --------------------------------------------------------------------------- #
@dataclass
class _Slabs:
    """Triangles bucketed by their ``x`` span, so a station only looks at a few."""

    triangles: list[tuple[Vec3, Vec3, Vec3]]
    x_lo: float
    width: float
    buckets: list[list[int]] = field(default_factory=list)

    @classmethod
    def build(cls, vertices: Sequence[Vec3], faces: Iterable[Sequence[int]],
              count: int = 256) -> "_Slabs":
        triangles = [tuple(vertices[i] for i in face) for face in faces]
        if not triangles:
            return cls(triangles, 0.0, 1.0, [[] for _ in range(count)])
        x_min = min(min(p[0] for p in tri) for tri in triangles)
        x_max = max(max(p[0] for p in tri) for tri in triangles)
        width = max((x_max - x_min) / count, 1e-12)
        buckets: list[list[int]] = [[] for _ in range(count)]
        for index, tri in enumerate(triangles):
            first = int((min(p[0] for p in tri) - x_min) / width)
            last = int((max(p[0] for p in tri) - x_min) / width)
            for bucket in range(max(first, 0), min(last, count - 1) + 1):
                buckets[bucket].append(index)
        return cls(triangles, x_min, width, buckets)

    def near(self, station: float) -> list[int]:
        bucket = int((station - self.x_lo) / self.width)
        if bucket < 0 or bucket >= len(self.buckets):
            return []
        return self.buckets[bucket]


def envelope_breadth(slabs: _Slabs, station: float, height: float) -> float:
    """Moulded breadth (``max |y|``) of the surface at ``x = station, z = height``.

    Zero when the station/height lies outside the model or the surface has a hole
    there -- "no surface here", never an interpolated guess.
    """
    best = 0.0
    for index in slabs.near(station):
        a, b, c = slabs.triangles[index]
        det = ((b[0] - a[0]) * (c[2] - a[2]) - (c[0] - a[0]) * (b[2] - a[2]))
        if abs(det) < 1e-14:
            continue  # triangle is edge-on in the x-z plane: no area to hit
        px = station - a[0]
        pz = height - a[2]
        u = (px * (c[2] - a[2]) - pz * (c[0] - a[0])) / det
        v = ((b[0] - a[0]) * pz - (b[2] - a[2]) * px) / det
        if u < -1e-12 or v < -1e-12 or u + v > 1.0 + 1e-12:
            continue
        y = a[1] + u * (b[1] - a[1]) + v * (c[1] - a[1])
        if abs(y) > best:
            best = abs(y)
    return best


def sectional_area(slabs: _Slabs, station: float, waterline: float,
                   levels: int = 48) -> float:
    """Area of the immersed part of the section at ``station`` (one side only)."""
    if levels < 1:
        raise ValueError("levels must be >= 1")
    step = waterline / levels
    total = 0.0
    previous = envelope_breadth(slabs, station, 0.0)
    for index in range(1, levels + 1):
        current = envelope_breadth(slabs, station, index * step)
        total += 0.5 * (previous + current) * step
        previous = current
    return total


def _clip_area_below(triangle: Sequence[Vec3], waterline: float) -> float:
    """Area of the part of a 3-D triangle with ``z <= waterline``."""
    polygon = list(triangle)
    if max(p[2] for p in polygon) <= waterline:
        return _triangle_area(polygon[0], polygon[1], polygon[2])
    if min(p[2] for p in polygon) >= waterline:
        return 0.0
    clipped: list[Vec3] = []
    for index in range(len(polygon)):
        start = polygon[index]
        end = polygon[(index + 1) % len(polygon)]
        inside_start = start[2] <= waterline
        inside_end = end[2] <= waterline
        if inside_start:
            clipped.append(start)
        if inside_start != inside_end:
            t = (waterline - start[2]) / (end[2] - start[2])
            clipped.append((start[0] + t * (end[0] - start[0]),
                            start[1] + t * (end[1] - start[1]),
                            waterline))
    if len(clipped) < 3:
        return 0.0
    area = 0.0
    for index in range(1, len(clipped) - 1):
        area += _triangle_area(clipped[0], clipped[index], clipped[index + 1])
    return area


def _triangle_area(a: Vec3, b: Vec3, c: Vec3) -> float:
    ux, uy, uz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
    vx, vy, vz = c[0] - a[0], c[1] - a[1], c[2] - a[2]
    cx = uy * vz - uz * vy
    cy = uz * vx - ux * vz
    cz = ux * vy - uy * vx
    return 0.5 * math.sqrt(cx * cx + cy * cy + cz * cz)


# --------------------------------------------------------------------------- #
# model classification
# --------------------------------------------------------------------------- #
def classify_model(vertices: Sequence[Vec3]) -> dict[str, Any]:
    """Is this a half-breadth model, a full hull, or neither?

    A publisher hull file with all ``y >= -tol`` and a maximum breadth well above
    zero is a half-breadth model: it describes one side of the ship and the other
    side is symmetry, so integrals double and the centreplane can never be sealed.
    """
    if not vertices:
        return {"half_breadth": False, "y_min": 0.0, "y_max": 0.0, "breadth": 0.0}
    ys = [p[1] for p in vertices]
    y_min, y_max = min(ys), max(ys)
    breadth = y_max - y_min
    tolerance = 1e-6 * max(abs(y_max), abs(y_min), 1.0)
    half = y_min >= -tolerance and y_max > 0.0
    return {
        "half_breadth": half,
        "y_min": y_min,
        "y_max": y_max,
        "breadth": breadth,
        "symmetry_factor": 2 if half else 1,
    }


# --------------------------------------------------------------------------- #
# whole-hull characteristics
# --------------------------------------------------------------------------- #
@dataclass
class HullCharacteristics:
    waterline: float
    stations: tuple[float, ...]
    sectional_areas: tuple[float, ...]
    volume: float
    s_wetted: float
    lwl: float
    beam: float
    mid_ship_area: float
    cb: float
    cm: float
    symmetry_factor: int
    station_count: int
    level_count: int

    def to_dict(self) -> dict[str, Any]:
        return {
            "waterline": self.waterline,
            "volume": self.volume,
            "s_wetted": self.s_wetted,
            "lwl": self.lwl,
            "beam": self.beam,
            "mid_ship_area": self.mid_ship_area,
            "cb": self.cb,
            "cm": self.cm,
            "symmetry_factor": self.symmetry_factor,
            "station_count": self.station_count,
            "level_count": self.level_count,
            "x_min": self.stations[0] if self.stations else None,
            "x_max": self.stations[-1] if self.stations else None,
        }


def hull_characteristics(vertices: Sequence[Vec3], faces: Sequence[Sequence[int]], *,
                         waterline: float, stations: int = 61, levels: int = 48,
                         slab_count: int = 256) -> HullCharacteristics:
    """Displacement volume, ``Cb``/``Cm``, ``Lwl``, beam and wetted area.

    Every quantity is measured on the surface as imported: no closure surface is
    invented, so this works on the open half-breadth workshop files.
    """
    if waterline <= 0.0:
        raise ValueError("waterline must be positive")
    slabs = _Slabs.build(vertices, faces, slab_count)
    if not slabs.triangles:
        raise ValueError("mesh has no triangles")
    xs = [min(p[0] for p in tri) for tri in slabs.triangles]
    x_lo, x_hi = min(xs), max(max(p[0] for p in tri) for tri in slabs.triangles)
    if x_hi <= x_lo:
        raise ValueError("mesh has zero x extent")
    grid = [x_lo + (x_hi - x_lo) * (index + 0.5) / stations for index in range(stations)]
    areas = [sectional_area(slabs, station, waterline, levels) for station in grid]
    # Midpoint rule over the full x extent: each cell contributes its centre
    # station's area times the cell width, so a constant-area body integrates exactly.
    step = (x_hi - x_lo) / stations
    half_volume = sum(areas) * step
    characteristics = classify_model(vertices)
    factor = characteristics["symmetry_factor"]

    wetted = 0.0
    symmetry_tolerance = 1e-6 * max(abs(p[1]) for tri in slabs.triangles for p in tri)
    for triangle in slabs.triangles:
        if characteristics["half_breadth"] and all(abs(p[1]) <= symmetry_tolerance
                                                  for p in triangle):
            # A face lying in the centreplane is a cut through the interior, not
            # wetted surface: counting it (and doubling) would invent area.
            continue
        wetted += _clip_area_below(triangle, waterline)

    beam = 0.0
    for station in grid:
        for level in range(1, levels + 1):
            beam = max(beam, envelope_breadth(slabs, station, waterline * level / levels))
    # Half-breadth model: the printed breadth is already the full beam.
    beam = 2.0 * beam if characteristics["half_breadth"] else beam

    wet_stations = [station for station, area in zip(grid, areas) if area > 2e-3 * max(areas or [0.0])]
    lwl = (wet_stations[-1] - wet_stations[0] + step) if len(wet_stations) >= 2 else 0.0
    volume = factor * half_volume
    mid_ship = max(areas) * (2.0 if characteristics["half_breadth"] else 1.0)
    cb = volume / (lwl * beam * waterline) if lwl > 0 and beam > 0 else 0.0
    cm = mid_ship / (beam * waterline) if beam > 0 else 0.0
    return HullCharacteristics(
        waterline=waterline,
        stations=tuple(grid),
        sectional_areas=tuple(areas),
        volume=volume,
        s_wetted=factor * wetted,
        lwl=lwl,
        beam=beam,
        mid_ship_area=mid_ship,
        cb=cb,
        cm=cm,
        symmetry_factor=factor,
        station_count=stations,
        level_count=levels,
    )


# --------------------------------------------------------------------------- #
# benchmark against published numbers
# --------------------------------------------------------------------------- #
def validate_against_particulars(characteristics: HullCharacteristics,
                                 particulars: dict[str, Any], *,
                                 tolerance_percent: float = 3.0) -> dict[str, Any]:
    """Per-quantity error of a measured hull against the published particulars.

    Returns one row per published quantity with reference, measured, relative
    error in per cent and the pass verdict against ``tolerance_percent``; the
    ``passed`` key is the conjunction, so a caller cannot accidentally read a
    green subset as a green whole.
    """
    rows: dict[str, Any] = {}
    for key, measured in (("volume", characteristics.volume),
                          ("s_wetted", characteristics.s_wetted),
                          ("cb", characteristics.cb)):
        reference = particulars.get(key)
        if reference in (None, 0):
            continue
        if key == "cb" and particulars.get("lpp") and particulars.get("beam"):
            # Published Cb is defined on Lpp, not on the length the file happens to
            # span (the bulb projects ahead of the FP), so the measured Cb is put on
            # the same basis before it is compared.
            measured = characteristics.volume / (float(particulars["lpp"])
                                                * characteristics.beam
                                                * characteristics.waterline)
        error = (measured - float(reference)) / float(reference) * 100.0
        rows[key] = {
            "reference": float(reference),
            "measured": measured,
            "error_percent": error,
            "passed": abs(error) <= tolerance_percent,
        }
    return {
        "source": particulars.get("source", ""),
        "tolerance_percent": tolerance_percent,
        "quantities": rows,
        "passed": bool(rows) and all(row["passed"] for row in rows.values()),
        "characteristics": characteristics.to_dict(),
    }