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
  Hodgman at ``z = T``, doubled for a half-breadth model.  Publisher files skin
  the hull with *overlapping* patches, so the raw sum counts the same skin twice:
  :func:`wetted_area` therefore removes hidden triangles first (centroid ray
  along the outward normal, plus a coplanar-coincidence test for sheets that lie
  on top of each other) and reports how much area it discarded.  ``dedup=False``
  reproduces the old raw sum, which stays available as the measured defect.

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


def envelope_breadth(slabs: _Slabs, station: float, height: float, *,
                     full: bool = False) -> float:
    """Breadth of the surface at ``x = station, z = height``.

    ``full=False`` (default) returns the moulded **half**-breadth ``max |y|``: the
    one-sided sheet that the publisher half-breadth models carry, which the caller
    doubles with ``symmetry_factor``.

    ``full=True`` returns the **full** breadth ``y_max - y_min`` of the hits, which is
    what a model carrying both sides needs.  Reporting ``max |y|`` there is exactly half
    the breadth, and it halves A(x), the displacement and the beam while leaving Cb
    untouched (both its numerator and its denominator shrink by two) -- see
    ``docs/gap_cad_t1_full_model_halving.md``.

    Zero when the station/height lies outside the model or the surface has a hole
    there -- "no surface here", never an interpolated guess.
    """
    best = 0.0
    best_hi = float("-inf")
    best_lo = float("inf")
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
        if y > best_hi:
            best_hi = y
        if y < best_lo:
            best_lo = y
    if full:
        return (best_hi - best_lo) if best_hi >= best_lo else 0.0
    return best


def sectional_area(slabs: _Slabs, station: float, waterline: float,
                   levels: int = 48, *, full: bool = False) -> float:
    """Area of the immersed part of the section at ``station``.

    ``full=False`` integrates the one-sided (half-breadth) breadth, so the result still
    has to be doubled by ``symmetry_factor``; ``full=True`` integrates ``y_max - y_min``
    and is the whole section area of a both-sides model.
    """
    if levels < 1:
        raise ValueError("levels must be >= 1")
    step = waterline / levels
    total = 0.0
    previous = envelope_breadth(slabs, station, 0.0, full=full)
    for index in range(1, levels + 1):
        current = envelope_breadth(slabs, station, index * step, full=full)
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
# hidden-surface removal on a patch soup
# --------------------------------------------------------------------------- #
def _sub(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _cross(u: Vec3, v: Vec3) -> Vec3:
    return (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])


def _dot(u: Vec3, v: Vec3) -> float:
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]


def _unit_normal(triangle: Sequence[Vec3]) -> Vec3:
    """Unit normal of a triangle, or ``(0, 0, 0)`` when it is degenerate."""
    normal = _cross(_sub(triangle[1], triangle[0]), _sub(triangle[2], triangle[0]))
    length = math.sqrt(_dot(normal, normal))
    if length == 0.0:
        return (0.0, 0.0, 0.0)
    return (normal[0] / length, normal[1] / length, normal[2] / length)


def _point_in_triangle(a: Vec3, b: Vec3, c: Vec3, point: Vec3, normal: Vec3,
                       tolerance: float) -> bool:
    """Is ``point`` inside (or within ``tolerance`` of) the triangle ``a, b, c``?"""
    for start, end in ((a, b), (b, c), (c, a)):
        edge = _sub(end, start)
        if _dot(_cross(edge, _sub(point, start)), normal) < -tolerance * math.sqrt(_dot(edge, edge)):
            return False
    return True


class _TriGrid:
    """Uniform 3-D cell index, so a ray only meets the triangles near it."""

    def __init__(self, triangles: Sequence[Sequence[Vec3]], cell: float) -> None:
        self.triangles = triangles
        self.cell = cell
        self.cells: dict[tuple[int, int, int], list[int]] = {}
        for index, triangle in enumerate(triangles):
            for key in self._span(triangle):
                self.cells.setdefault(key, []).append(index)

    def _span(self, triangle: Sequence[Vec3]) -> list[tuple[int, int, int]]:
        lo = [min(p[axis] for p in triangle) for axis in range(3)]
        hi = [max(p[axis] for p in triangle) for axis in range(3)]
        return [(i, j, k)
                for i in range(int(math.floor(lo[0] / self.cell)), int(math.floor(hi[0] / self.cell)) + 1)
                for j in range(int(math.floor(lo[1] / self.cell)), int(math.floor(hi[1] / self.cell)) + 1)
                for k in range(int(math.floor(lo[2] / self.cell)), int(math.floor(hi[2] / self.cell)) + 1)]

    def at(self, point: Sequence[float]) -> list[int]:
        key = (int(math.floor(point[0] / self.cell)),
               int(math.floor(point[1] / self.cell)),
               int(math.floor(point[2] / self.cell)))
        return self.cells.get(key, [])


def _march(candidates: list[int], grid: _TriGrid, origin: Vec3, direction: Vec3,
           stride: float, reach: float, skip: int, tolerance: float,
           source: int | None = None) -> bool:
    """Walk the ray ``origin + t*direction`` to ``reach``; True if a sheet is hit."""
    steps = int(math.ceil(reach / stride))
    visited: set[int] = set()
    for step in range(steps + 1):
        point = (origin[0] + direction[0] * step * stride,
                 origin[1] + direction[1] * step * stride,
                 origin[2] + direction[2] * step * stride)
        candidates.extend(grid.at(point))
    for index in dict.fromkeys(candidates):
        if index == skip:
            continue
        a, b, c = grid.triangles[index]
        normal = _unit_normal((a, b, c))
        if normal == (0.0, 0.0, 0.0):
            continue
        # Coincident sheets do not cross the ray at a positive distance -- the
        # centroid we cast from lies *on* them.  That is still an occlusion, and
        # it is how two patches that tile the same skin get counted once.
        # Tie-break by index: two *exactly* coincident sheets would otherwise hide
        # each other and both vanish from the sum.  Only a lower-numbered sheet may
        # claim a coincident one, so each coincident cluster keeps exactly one copy.
        if (index < source and abs(_dot(normal, _sub(origin, a))) <= tolerance
                and _point_in_triangle(a, b, c, origin, normal, tolerance)):
            return True
        denominator = _dot(normal, direction)
        if abs(denominator) < 1e-9:
            continue
        distance = _dot(normal, _sub(a, origin)) / denominator
        if distance <= 1e-6:
            continue
        hit = (origin[0] + distance * direction[0],
               origin[1] + distance * direction[1],
               origin[2] + distance * direction[2])
        if _point_in_triangle(a, b, c, hit, normal, tolerance):
            return True
    return False


def wetted_area(vertices: Sequence[Vec3], faces: Iterable[Sequence[int]], *,
                waterline: float, dedup: bool = True,
                cell_size: float | None = None,
                tolerance: float | None = None) -> tuple[float, float]:
    """Wetted surface at ``waterline``: ``(area, discarded)``, both full-hull.

    ``area`` is the skin of the hull below the waterline, doubled for a
    half-breadth model; ``discarded`` is the area of the triangles that turned out
    to be hidden under another sheet (``0.0`` when ``dedup`` is off).

    Publisher IGES hull files describe the skin with untrimmed patches that
    *overlap*: summing triangle areas then counts the shared skin twice (on
    ``KCS_hull_SVA.igs`` this is +19.4 % against the published 9 424 m²).  With
    ``dedup`` each triangle casts a ray from its centroid along its outward
    normal; a hit means another sheet covers it and it is dropped from the sum.
    Two coincident sheets have no crossing point, so a coplanar test at the same
    tolerance catches those as well.

    The outward direction is taken from the bounding-box centroid, which is exact
    for a star-shaped (convex-ish) hull skin and the reason this is not a general
    solid modeller: a deeply concave pocket can be hidden in the wrong direction.
    ``cell_size`` and ``tolerance`` default to a fraction of the model size.
    """
    if waterline <= 0.0:
        raise ValueError("waterline must be positive")
    triangles = [tuple(vertices[i] for i in face) for face in faces]
    if not triangles:
        raise ValueError("mesh has no triangles")
    characteristics = classify_model(vertices)
    factor = characteristics["symmetry_factor"]
    span = max(max(p[axis] for t in triangles for p in t) - min(p[axis] for t in triangles for p in t)
               for axis in range(3))
    if span <= 0.0:
        raise ValueError("mesh has zero extent")
    cell = cell_size if cell_size is not None else max(span / 64.0, 1e-9)
    tol = tolerance if tolerance is not None else 1e-4 * span
    symmetry_tolerance = 1e-6 * max(abs(p[1]) for t in triangles for p in t)

    def is_centreplane(triangle: Sequence[Vec3]) -> bool:
        return characteristics["half_breadth"] and all(
            abs(p[1]) <= symmetry_tolerance for p in triangle)

    if not dedup:
        total = 0.0
        for triangle in triangles:
            if not is_centreplane(triangle):
                total += _clip_area_below(triangle, waterline)
        return factor * total, 0.0

    points = [p for t in triangles for p in t]
    centre = tuple(sum(p[axis] for p in points) / len(points) for axis in range(3))
    grid = _TriGrid(triangles, cell)
    stride = 0.5 * cell
    kept = 0.0
    discarded = 0.0
    for index, triangle in enumerate(triangles):
        if is_centreplane(triangle):
            continue
        area = _clip_area_below(triangle, waterline)
        if area == 0.0:
            continue
        normal = _unit_normal(triangle)
        if normal == (0.0, 0.0, 0.0):
            kept += area
            continue
        centroid = tuple(sum(p[axis] for p in triangle) / 3.0 for axis in range(3))
        if _dot(normal, _sub(centroid, centre)) < 0.0:
            normal = (-normal[0], -normal[1], -normal[2])
        offset = 1e-3 * tol
        origin = (centroid[0] + offset * normal[0],
                  centroid[1] + offset * normal[1],
                  centroid[2] + offset * normal[2])
        if _march([], grid, origin, normal, stride, reach=4.0 * cell,
                  skip=index, tolerance=tol, source=index):
            discarded += area
        else:
            kept += area
    return factor * kept, factor * discarded


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
    #: Raw triangle sum, i.e. what the soup measures before hidden sheets are
    #: removed; kept so the de-duplication stays auditable.
    s_wetted_naive: float = 0.0
    #: Area dropped as hidden under an overlapping sheet.
    s_wetted_discarded: float = 0.0

    def to_dict(self) -> dict[str, Any]:
        return {
            "waterline": self.waterline,
            "volume": self.volume,
            "s_wetted": self.s_wetted,
            "s_wetted_naive": self.s_wetted_naive,
            "s_wetted_discarded": self.s_wetted_discarded,
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
                         slab_count: int = 256, dedup: bool = True,
                         full_breadth: bool | None = None) -> HullCharacteristics:
    """Displacement volume, ``Cb``/``Cm``, ``Lwl``, beam and wetted area.

    Every quantity is measured on the surface as imported: no closure surface is
    invented, so this works on the open half-breadth workshop files.  ``dedup``
    passes through to :func:`wetted_area`; the raw sum is reported alongside in
    ``s_wetted_naive`` so the effect of the de-duplication is always visible.

    ``full_breadth`` selects how a section's breadth is read: ``None`` (default) infers
    it from :func:`classify_model` -- one-sided ``max |y|`` for a half-breadth sheet,
    ``y_max - y_min`` for a model carrying both sides.  Forcing ``False`` on a
    both-sides model reproduces the pre-fix behaviour (displacement and beam exactly
    halved, Cb unchanged); see ``docs/gap_cad_t1_full_model_halving.md``.
    """
    if waterline <= 0.0:
        raise ValueError("waterline must be positive")
    slabs = _Slabs.build(vertices, faces, slab_count)
    if not slabs.triangles:
        raise ValueError("mesh has no triangles")
    characteristics = classify_model(vertices)
    if full_breadth is None:
        full_breadth = not characteristics["half_breadth"]
    xs = [min(p[0] for p in tri) for tri in slabs.triangles]
    x_lo, x_hi = min(xs), max(max(p[0] for p in tri) for tri in slabs.triangles)
    if x_hi <= x_lo:
        raise ValueError("mesh has zero x extent")
    grid = [x_lo + (x_hi - x_lo) * (index + 0.5) / stations for index in range(stations)]
    areas = [sectional_area(slabs, station, waterline, levels, full=full_breadth)
             for station in grid]
    # Midpoint rule over the full x extent: each cell contributes its centre
    # station's area times the cell width, so a constant-area body integrates exactly.
    step = (x_hi - x_lo) / stations
    half_volume = sum(areas) * step
    factor = characteristics["symmetry_factor"]

    # A face lying in the centreplane is a cut through the interior, not wetted
    # surface: counting it (and doubling) would invent area -- `wetted_area`
    # skips those, and removes the sheets that hide under an overlapping patch.
    wetted, discarded = wetted_area(vertices, faces, waterline=waterline, dedup=dedup)
    wetted_naive = wetted + discarded if dedup else wetted
    if dedup:
        wetted_naive = wetted_area(vertices, faces, waterline=waterline, dedup=False)[0]

    beam = 0.0
    for station in grid:
        for level in range(1, levels + 1):
            beam = max(beam, envelope_breadth(slabs, station, waterline * level / levels,
                                              full=full_breadth))
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
        s_wetted=wetted,
        lwl=lwl,
        beam=beam,
        mid_ship_area=mid_ship,
        cb=cb,
        cm=cm,
        symmetry_factor=factor,
        station_count=stations,
        level_count=levels,
        s_wetted_naive=wetted_naive,
        s_wetted_discarded=discarded,
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