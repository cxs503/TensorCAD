"""Exact hydrostatics of a rectangular barge at **combined heel and trim**.

Why this module exists
----------------------
Every hydrostatic reference in this package is one-dimensional in attitude.
:mod:`tensorcad.gz` / :mod:`tensorcad.box_gz` heel a section about a fixed axis and never
trim; :mod:`tensorcad.trim` assumes the hull is upright; :func:`tensorcad.hydrostatics.
compute_hydrostatics` cuts the mesh with a *horizontal* plane.  A barge that is both heeled
and trimmed therefore has no exact anchor, although it is the attitude a ship actually
floats at in a turn and one a damage case can force.  The mesh route can be pointed at it
(rotate the mesh, then cut with a horizontal plane), but a mesh route verified only against
other mesh routes cannot see a modelling error they share -- the argument of
``docs/box_large_heel.md``.

What is exact here
------------------
The submerged body of a box is ``box intersected with a half-space``, i.e. a convex
polyhedron.  Clipping the six box faces by the waterline plane and closing the cut with the
waterplane polygon gives that polyhedron *exactly*: no discretisation and no case analysis.
Bilge emergence, deck-edge immersion and a cut that passes through a single corner are all
one code path, which is what makes this usable as a reference beyond the wall-sided window
where the one-line formula :func:`tensorcad.gz.box_righting_arms_analytic` refuses to
answer.

* volume and buoyancy centre -- divergence theorem over the triangulated polyhedron;
* waterplane area, centroid and second moments -- Green's theorem over the cap polygon,
  measured *in the waterplane* (the plane is tilted, so the in-plane area is what ``BM``
  wants);
* the righting arm ``GZ`` -- the lever of the buoyancy centre about the centre of gravity
  along the horizontal transverse direction *of that waterplane*, so both attitudes enter
  one expression instead of two formulas glued together.

Conventions
-----------
``x`` forward, ``y`` starboard, ``z`` up from the keel; the box occupies
``[-L/2, L/2] x [-B/2, B/2] x [0, depth]``.  ``heel_deg > 0`` puts the starboard side down
(the waterline in body coordinates rises with ``y``); ``trim_deg > 0`` puts the bow down.
The waterline plane is ``n . p = c`` with the up-normal

    n = normalize(-sin(tau), -sin(phi) cos(tau), cos(phi) cos(tau))

so upright it is ``(0, 0, 1)`` and the one-axis limits match the documented conventions of
:mod:`tensorcad.gz`.  The submerged side is ``n . p <= c``.

What this is not
----------------
A box.  It says nothing about a hull with a shaped section -- for that the mesh routes are
what to use, and this module is the attitude reference they can be checked against.  No
free-surface correction, no sinkage/trim *solve* (the attitude is an input and the
displacement is what is held), no windage, no stability criterion.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable, Sequence

Vec3 = tuple[float, float, float]

#: Points closer than this fraction of the model size are the same point when the cut
#: polygon is assembled from per-face intersection points.
_RELATIVE_EPS = 1e-12


class BoxInclinedError(ValueError):
    """The requested box, attitude or displacement cannot be realised."""


# ----------------------------------------------------------------------------------
# vector helpers (standard library only: the zero-dependency core policy of this package)
# ----------------------------------------------------------------------------------
def _sub(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _dot(a: Vec3, b: Vec3) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a: Vec3, b: Vec3) -> Vec3:
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _norm(a: Vec3) -> float:
    return math.sqrt(_dot(a, a))


def _unit(a: Vec3) -> Vec3:
    length = _norm(a)
    if length == 0.0:
        raise BoxInclinedError("a zero vector has no direction")
    return (a[0] / length, a[1] / length, a[2] / length)


def _project_into_plane(vector: Vec3, n: Vec3) -> Vec3:
    return _sub(vector, tuple(n[k] * _dot(vector, n) for k in range(3)))


def _finite(name: str, value: float) -> float:
    if not isinstance(value, (int, float)) or not math.isfinite(value):
        raise BoxInclinedError(f"{name} must be a finite number, got {value!r}")
    return float(value)


def _positive(name: str, value: float) -> float:
    value = _finite(name, value)
    if value <= 0.0:
        raise BoxInclinedError(f"{name} must be positive, got {value!r}")
    return value


# ----------------------------------------------------------------------------------
# the waterplane
# ----------------------------------------------------------------------------------
def waterline_normal(heel_deg: float, trim_deg: float) -> Vec3:
    """Unit up-normal of the waterline, expressed in body coordinates.

    Upright this is ``(0, 0, 1)``; heel alone gives ``(0, -sin phi, cos phi)`` and trim
    alone ``(-sin tau, 0, cos tau)`` -- the compositions ``R_y(-tau) R_x(phi)`` applied to
    the earth up-vector, written out so there is no matrix convention to get wrong.
    """
    phi = math.radians(_finite("heel_deg", heel_deg))
    tau = math.radians(_finite("trim_deg", trim_deg))
    return _unit((-math.sin(tau), -math.sin(phi) * math.cos(tau), math.cos(phi) * math.cos(tau)))


# ----------------------------------------------------------------------------------
# clipping
# ----------------------------------------------------------------------------------
def _clip_polygon(
    poly: Sequence[Vec3], n: Vec3, c: float, eps: float
) -> tuple[list[Vec3], list[Vec3]]:
    """Keep ``n . p <= c`` on ``poly``; also return the points that landed on the plane."""
    kept: list[Vec3] = []
    cut: list[Vec3] = []
    count = len(poly)
    for index in range(count):
        a = poly[index]
        b = poly[(index + 1) % count]
        da = _dot(n, a) - c
        db = _dot(n, b) - c
        if da <= eps:
            kept.append(a)
        crosses = (da < -eps and db > eps) or (db < -eps and da > eps)
        if crosses:
            t = da / (da - db)
            point = (a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]), a[2] + t * (b[2] - a[2]))
            kept.append(point)
            cut.append(point)
        elif abs(da) <= eps:
            cut.append(a)
    return kept, cut


def _dedupe(points: Iterable[Vec3], eps: float) -> list[Vec3]:
    unique: list[Vec3] = []
    for point in points:
        if not any(_norm(_sub(point, other)) <= eps for other in unique):
            unique.append(point)
    return unique


def _order_in_plane(points: Sequence[Vec3], n: Vec3) -> list[Vec3]:
    """Order a convex planar point set counter-clockwise about its centroid."""
    if len(points) < 3:
        return list(points)
    centroid = (
        sum(p[0] for p in points) / len(points),
        sum(p[1] for p in points) / len(points),
        sum(p[2] for p in points) / len(points),
    )
    seed = (1.0, 0.0, 0.0) if abs(n[0]) < 0.9 else (0.0, 1.0, 0.0)
    u = _unit(_project_into_plane(seed, n))
    v = _cross(n, u)

    def angle(point: Vec3) -> float:
        d = _sub(point, centroid)
        return math.atan2(_dot(d, v), _dot(d, u))

    return sorted(points, key=angle)


def _box_faces(length: float, beam: float, depth: float) -> list[list[Vec3]]:
    """The six faces of the box, each wound outwardly."""
    half_l, half_b = length / 2.0, beam / 2.0
    x0, x1, y0, y1, z0, z1 = -half_l, half_l, -half_b, half_b, 0.0, depth
    return [
        [(x0, y0, z0), (x0, y1, z0), (x1, y1, z0), (x1, y0, z0)],  # keel, normal -z
        [(x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1)],  # deck,  normal +z
        [(x0, y0, z0), (x1, y0, z0), (x1, y0, z1), (x0, y0, z1)],  # port,  normal -y
        [(x1, y1, z0), (x0, y1, z0), (x0, y1, z1), (x1, y1, z1)],  # stbd,  normal +y
        [(x0, y1, z0), (x0, y0, z0), (x0, y0, z1), (x0, y1, z1)],  # transom, normal -x
        [(x1, y0, z0), (x1, y1, z0), (x1, y1, z1), (x1, y0, z1)],  # stem,   normal +x
    ]


def _polyhedron_moments(faces: Sequence[Sequence[Vec3]]) -> tuple[float, Vec3]:
    """Signed volume and first moments, by the divergence theorem about the origin."""
    volume = 0.0
    mx = my = mz = 0.0
    for poly in faces:
        anchor = poly[0]
        for index in range(1, len(poly) - 1):
            b, c = poly[index], poly[index + 1]
            tetra = _dot(anchor, _cross(b, c)) / 6.0
            volume += tetra
            mx += tetra * (anchor[0] + b[0] + c[0]) / 4.0
            my += tetra * (anchor[1] + b[1] + c[1]) / 4.0
            mz += tetra * (anchor[2] + b[2] + c[2]) / 4.0
    return volume, (mx, my, mz)


def _triangle_area(a: Vec3, b: Vec3, c: Vec3) -> float:
    return 0.5 * _norm(_cross(_sub(b, a), _sub(c, a)))


def _polygon_props(
    points2d: Sequence[tuple[float, float]],
) -> tuple[float, float, float, float, float]:
    """``(area, cx, cy, ixx_c, iyy_c)`` of a planar polygon by Green's theorem.

    ``ixx_c`` is ``integral (y - cy)^2 dA``.  Every term carries the same oriented
    ``cross_i`` factor, so the winding only changes the sign, which the caller normalises.
    """
    count = len(points2d)
    area = 0.0
    cx = cy = 0.0
    ixx_o = iyy_o = 0.0
    for index in range(count):
        x0, y0 = points2d[index]
        x1, y1 = points2d[(index + 1) % count]
        cross = x0 * y1 - x1 * y0
        area += cross
        cx += (x0 + x1) * cross
        cy += (y0 + y1) * cross
        ixx_o += (y0 * y0 + y0 * y1 + y1 * y1) * cross
        iyy_o += (x0 * x0 + x0 * x1 + x1 * x1) * cross
    area *= 0.5
    if area == 0.0:
        return 0.0, 0.0, 0.0, 0.0, 0.0
    cx /= 6.0 * area
    cy /= 6.0 * area
    ixx_o /= 12.0
    iyy_o /= 12.0
    return area, cx, cy, ixx_o - area * cy * cy, iyy_o - area * cx * cx


# ----------------------------------------------------------------------------------
# the exact reference
# ----------------------------------------------------------------------------------
@dataclass(frozen=True)
class InclinedBoxHydrostatics:
    """Exact hydrostatics of a rectangular barge floating at a heel *and* a trim."""

    heel_deg: float
    trim_deg: float
    length: float
    beam: float
    depth: float
    #: waterline plane ``n . p = waterline_offset``, in body coordinates
    waterline_normal: Vec3
    waterline_offset: float
    volume_m3: float
    #: buoyancy centre in body coordinates, i.e. ``(LCB, TCB, VCB)``
    centre_of_buoyancy: Vec3
    waterplane_area_m2: float
    waterplane_centroid: Vec3
    #: in-plane second moments about the waterplane centroid: across the width, along it
    i_transverse_m4: float
    i_longitudinal_m4: float
    wetted_area_m2: float
    waterplane_polygon: tuple[Vec3, ...]
    vertices: int
    faces: int
    #: volume integrated a second time from the clipped polyhedron (closure check)
    polyhedron_volume_m3: float

    @property
    def lcb(self) -> float:
        return self.centre_of_buoyancy[0]

    @property
    def tcb(self) -> float:
        return self.centre_of_buoyancy[1]

    @property
    def vcb(self) -> float:
        return self.centre_of_buoyancy[2]

    @property
    def bm_transverse(self) -> float | None:
        return None if self.volume_m3 <= 0.0 else self.i_transverse_m4 / self.volume_m3

    @property
    def bm_longitudinal(self) -> float | None:
        return None if self.volume_m3 <= 0.0 else self.i_longitudinal_m4 / self.volume_m3

    def gz(self, kg: float) -> float:
        """Righting lever about ``G = (0, 0, kg)``, positive when restoring.

        Measured along the horizontal transverse direction *of this waterplane*: the body
        ``x`` axis projected into the plane, crossed with the plane normal.  A trimmed
        waterplane therefore measures the arm across its own width, not the body's.
        """
        n = self.waterline_normal
        forward = (1.0, 0.0, 0.0)
        theta = math.radians(_finite("heel_deg", self.heel_deg))
        transverse = _unit(_cross(n, _unit(_project_into_plane(forward, n))))
        if math.cos(theta) < 0.0:  # heeled past 90 deg: keep the lever's sense continuous
            transverse = (-transverse[0], -transverse[1], -transverse[2])
        lever = _sub(self.centre_of_buoyancy, (0.0, 0.0, _finite("kg", kg)))
        return _dot(transverse, lever)

    def to_dict(self) -> dict[str, object]:
        return {
            "heel_deg": round(self.heel_deg, 12),
            "trim_deg": round(self.trim_deg, 12),
            "length": self.length,
            "beam": self.beam,
            "depth": self.depth,
            "waterline_normal": [round(v, 15) for v in self.waterline_normal],
            "waterline_offset": round(self.waterline_offset, 15),
            "volume_m3": round(self.volume_m3, 15),
            "polyhedron_volume_m3": round(self.polyhedron_volume_m3, 15),
            "centre_of_buoyancy": [round(v, 15) for v in self.centre_of_buoyancy],
            "lcb": round(self.lcb, 15),
            "tcb": round(self.tcb, 15),
            "vcb": round(self.vcb, 15),
            "waterplane_area_m2": round(self.waterplane_area_m2, 15),
            "waterplane_centroid": [round(v, 15) for v in self.waterplane_centroid],
            "i_transverse_m4": round(self.i_transverse_m4, 15),
            "i_longitudinal_m4": round(self.i_longitudinal_m4, 15),
            "bm_transverse": None if self.bm_transverse is None else round(self.bm_transverse, 15),
            "bm_longitudinal": (
                None if self.bm_longitudinal is None else round(self.bm_longitudinal, 15)
            ),
            "wetted_area_m2": round(self.wetted_area_m2, 15),
            "vertices": self.vertices,
            "faces": self.faces,
            "waterplane_polygon": [[round(c, 12) for c in p] for p in self.waterplane_polygon],
        }


def submerged_box(
    length: float,
    beam: float,
    depth: float,
    *,
    heel_deg: float = 0.0,
    trim_deg: float = 0.0,
    waterline_offset: float,
) -> InclinedBoxHydrostatics:
    """Clip the box by ``n . p <= waterline_offset`` and integrate the result exactly."""
    length = _positive("length", length)
    beam = _positive("beam", beam)
    depth = _positive("depth", depth)
    n = waterline_normal(heel_deg, trim_deg)
    c = _finite("waterline_offset", waterline_offset)
    eps = _RELATIVE_EPS * max(length, beam, depth)

    faces: list[list[Vec3]] = []
    cap_points: list[Vec3] = []
    wetted_area = 0.0
    for face in _box_faces(length, beam, depth):
        kept, cut = _clip_polygon(face, n, c, eps)
        cap_points.extend(cut)
        if len(kept) >= 3:
            faces.append(kept)
            anchor = kept[0]
            for index in range(1, len(kept) - 1):
                wetted_area += _triangle_area(anchor, kept[index], kept[index + 1])

    cap = _order_in_plane(_dedupe(cap_points, eps), n)
    if len(cap) >= 3:
        normal = _cross(_sub(cap[1], cap[0]), _sub(cap[2], cap[0]))
        faces.append(cap if _dot(normal, n) > 0.0 else list(reversed(cap)))

    volume, moments = _polyhedron_moments(faces)
    if volume < 0.0:
        faces = [list(reversed(poly)) for poly in faces]
        volume, moments = _polyhedron_moments(faces)
    if volume <= 0.0:
        raise BoxInclinedError(
            f"the waterline offset {waterline_offset!r} leaves no submerged volume "
            f"(heel {heel_deg!r} deg, trim {trim_deg!r} deg)"
        )

    area = i_trans = i_long = 0.0
    centroid = (0.0, 0.0, 0.0)
    if len(cap) >= 3:
        seed = (1.0, 0.0, 0.0) if abs(n[0]) < 0.9 else (0.0, 1.0, 0.0)
        u = _unit(_project_into_plane(seed, n))
        v = _cross(n, u)
        coords = [(_dot(p, u), _dot(p, v)) for p in cap]
        area, cu, cv, ixx, iyy = _polygon_props(coords)
        area, ixx, iyy = abs(area), abs(ixx), abs(iyy)
        centroid = tuple(cu * u[k] + cv * v[k] for k in range(3))  # type: ignore[assignment]
        i_trans, i_long = ixx, iyy

    vertices = len({tuple(round(co, 12) for co in p) for poly in faces for p in poly})
    return InclinedBoxHydrostatics(
        heel_deg=float(heel_deg),
        trim_deg=float(trim_deg),
        length=length,
        beam=beam,
        depth=depth,
        waterline_normal=n,
        waterline_offset=c,
        volume_m3=volume,
        centre_of_buoyancy=(moments[0] / volume, moments[1] / volume, moments[2] / volume),
        waterplane_area_m2=area,
        waterplane_centroid=centroid,
        i_transverse_m4=i_trans,
        i_longitudinal_m4=i_long,
        wetted_area_m2=wetted_area,
        waterplane_polygon=tuple(cap),
        vertices=vertices,
        faces=len(faces),
        polyhedron_volume_m3=volume,
    )


def box_inclined_hydrostatics(
    length: float,
    beam: float,
    depth: float,
    *,
    volume_m3: float,
    heel_deg: float = 0.0,
    trim_deg: float = 0.0,
    tolerance: float = 1e-15,
) -> InclinedBoxHydrostatics:
    """The exact hydrostatics of the box floating at ``volume_m3`` with that attitude.

    The waterline offset is solved by bisection on ``volume(c)``, which is strictly
    monotone because the submerged region only grows with ``c``; no derivative and no
    interpolation are involved, so the returned geometry is exact to round-off.
    """
    length = _positive("length", length)
    beam = _positive("beam", beam)
    depth = _positive("depth", depth)
    target = _positive("volume_m3", volume_m3)
    if target > length * beam * depth * (1.0 + 1e-12):
        raise BoxInclinedError(
            f"volume {target!r} m^3 exceeds the box itself ({length * beam * depth!r} m^3)"
        )
    n = waterline_normal(heel_deg, trim_deg)
    scale = max(length, beam, depth)
    projections = [
        _dot(n, corner)
        for corner in (
            (x, y, z)
            for x in (-length / 2.0, length / 2.0)
            for y in (-beam / 2.0, beam / 2.0)
            for z in (0.0, depth)
        )
    ]
    lo = min(projections) - scale
    hi = max(projections) + scale
    while hi - lo > tolerance * scale:
        mid = 0.5 * (lo + hi)
        try:
            probe = submerged_box(
                length, beam, depth, heel_deg=heel_deg, trim_deg=trim_deg, waterline_offset=mid
            ).volume_m3
        except BoxInclinedError:
            probe = 0.0
        if probe < target:
            lo = mid
        else:
            hi = mid
    return submerged_box(
        length, beam, depth, heel_deg=heel_deg, trim_deg=trim_deg, waterline_offset=0.5 * (lo + hi)
    )


def box_inclined_gz(
    length: float,
    beam: float,
    depth: float,
    *,
    kg: float,
    volume_m3: float,
    heel_deg: float,
    trim_deg: float = 0.0,
) -> float:
    """Righting lever of the box at ``(heel_deg, trim_deg)`` and fixed displacement."""
    hydro = box_inclined_hydrostatics(
        length, beam, depth, volume_m3=volume_m3, heel_deg=heel_deg, trim_deg=trim_deg
    )
    return hydro.gz(kg)


def box_inclined_curve(
    length: float,
    beam: float,
    depth: float,
    *,
    kg: float,
    volume_m3: float,
    angles_deg: Sequence[float],
    trim_deg: float = 0.0,
) -> tuple[float, ...]:
    """``GZ`` at a list of heels, one exact solve per angle."""
    return tuple(
        box_inclined_gz(
            length, beam, depth, kg=kg, volume_m3=volume_m3, heel_deg=angle, trim_deg=trim_deg
        )
        for angle in angles_deg
    )


__all__ = [
    "BoxInclinedError",
    "InclinedBoxHydrostatics",
    "box_inclined_curve",
    "box_inclined_gz",
    "box_inclined_hydrostatics",
    "submerged_box",
    "waterline_normal",
]