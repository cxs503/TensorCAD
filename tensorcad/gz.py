"""Large-angle stability: the righting-arm curve ``GZ(phi)``.

:mod:`tensorcad.stability` answers the small-angle question -- where is ``G`` allowed to
sit before the upright position stops being stable -- from the waterplane moments, and
says so: ``GM`` is the *slope* of the righting-arm curve at ``phi = 0`` and nothing more.
This module crosses that line and integrates the **heeled displacement** instead, which
is a different sum: at heel the waterplane moves, the displaced volume has to be held
constant by finding the new waterline, and the centre of buoyancy walks out sideways.
That sideways walk is the whole of the righting moment.

The righting arm of a body heeled by ``phi`` about its longitudinal axis is

    GZ(phi) = KN(phi) - KG sin(phi),      KN(phi) = y_B cos(phi) + z_B sin(phi),

where ``(y_B, z_B)`` is the centre of buoyancy of the *heeled* displaced volume measured
in the body frame from the keel datum, and ``KN`` is its lever about the keel line. The
positive sense used throughout is: ``phi > 0`` heels the ``+y`` (starboard) side **down**,
the buoyancy centre therefore moves towards ``+y``, and ``GZ > 0`` means the resulting
moment is a *righting* one.

**How the heeled volume is measured.** The mesh is rotated about the ``x`` axis, then cut
by the horizontal plane ``z = c`` and the part below it is integrated. ``c`` is not a free
parameter -- a floating body keeps its displacement -- so it is solved for each heel angle
with Newton's method, using the exact derivative ``dV/dc = A(c)`` (the waterplane area,
which the same pass accumulates as the projected area of the wetted triangles; over the
closed surface that sum vanishes, so the wetted part's contribution is minus the cap's,
and its magnitude *is* the cap area).

Volume and first moments come from the divergence theorem with the apex at the origin,
placed **on the waterline plane** -- so the flat cap that closes the wetted surface
contributes exactly zero and the wetted triangles alone give the answer. For a polyhedral
mesh this is exact; the only error left is double-precision round-off. ``displacement_volume``
is reported per point and is not a fit: every point re-solves the waterline so that
``V(phi) = V(0)`` holds to the quoted tolerance.

**What is verified, and against what.** Two shapes have an exact righting arm with no
discretisation error at all, and both are used as references in the tests:

* a **box barge** with vertical sides. The submerged section is a trapezoid with parallel
  sides ``h1 = T - (B/2)tan(phi)`` and ``h2 = T + (B/2)tan(phi)``, so integrating
  ``int z dA`` and ``int y dA`` over it and nothing else gives

      KB(phi) = T/2 + (B^2/(24 T)) tan^2(phi)  = KB0 + BM0/2 * tan^2(phi)
      y_B(phi) = (B^2/(12 T)) tan(phi)         = BM0 tan(phi)
      KN(phi)  = sin(phi) (KB0 + BM0 + BM0/2 tan^2(phi))
      GZ(phi)  = sin(phi) (GM0 + BM0/2 tan^2(phi))     with BM0 = B^2/(12T),

  which is the classical *wall-sided* formula -- derived here rather than quoted, so the
  test against it exercises the integrator instead of a copied constant. It stops being
  exact when the waterline reaches the bilge on the rising side (``tan(phi) = 2T/B``) or
  the deck edge on the dipping side (``tan(phi) = 2(D-T)/B``);
  :func:`box_wall_sided_limit_deg` returns that bound, and the tests assert the formula
  *disagrees* beyond it so the comparison cannot pass by being vacuous.
* a **cylinder of circular section floating at half depth** ("a log"). Any line through
  the centre of a centrally symmetric section splits it into equal areas, so the waterline
  stays through the section centre at every heel; the submerged centroid then sits a
  distance ``d`` below that centre along the waterline normal, and rotating it back into the
  body frame gives ``KN = R sin(phi) + s``, with ``R`` the centre's height above the keel and
  ``s`` the centroid's offset *along* the cut line. For a circle ``d`` is constant and ``s``
  vanishes, so ``KN = R sin(phi)`` and ``GZ(phi) = (R - KG) sin(phi)`` for every heel up to
  90 degrees: the round bilge is stable only with ``G`` below the section centre, which is
  the same statement as ``GM = KB + BM - KG = R - KG`` -- because ``KB = R - 4R/(3 pi)`` and
  ``BM = 4R/(3 pi)`` (``I_T = L(2R)^3/12``) sum to exactly ``R``. A regular even-sided polygon
  is centrally symmetric too, so the same closed form holds for the polygonal mesh; ``s``
  vanishes exactly on cut directions that are symmetry axes of the polygon.

**What this is not.** ``KG`` is the only loading input -- a height above the *keel*, as
:func:`tensorcad.stability.initial_stability` takes it, because ``KN`` is levered about the
keel line -- and there are no tank contents, so
:func:`free_surface_correction` offers the flat reduction ``sum(i_x)/V`` that a designer
applies by hand, and nothing more. :func:`with_freeboard` exists because a hull generated
to its design waterline simply *cannot* be heeled -- on the dipping side the waterline
leaves the modelled surface immediately -- so a closed freeboard band has to be added
first. Damage stability, deck-edge immersion criteria, grain heeling, statutory maxima and
class rules are all outside this module, and nothing here should be read as a check
against a statutory criterion.

Only the Python standard library is used.

**A second integration route, kept alongside this one.** This module carries two
independent large-angle implementations, merged rather than replaced. The note below
documents the ``submerged_solid`` / ``solve_waterline`` / ``gz_curve`` route (its own
crossing-triangle clip, its own sinkage solve, and the IMO intact-stability screening);
the note above documents the ``heeled_hydrostatics`` / ``righting_arm_curve`` route.
They share the ``GZ = KN - KG sin(phi)`` definition and cross-check each other through
``tests/test_gz.py``, which exercises both.

Large-angle stability: the ``GZ`` curve from heeled plane cuts.

:mod:`tensorcad.stability` stops at the upright position: it reports ``GM``, the slope
of the righting-arm curve at zero heel. ``GZ(phi)`` itself -- the quantity an intact
stability booklet, an IMO area criterion and a ``KG`` limit all actually need -- requires
something that module deliberately does not do: a **heeled** waterplane, and, because the
displaced volume must stay constant as the hull goes over, a sinkage solve at every angle.

This module supplies both and nothing else. It does not add a new integration scheme: the
volume of the submerged part of the surface, the transverse position of its centroid and
the waterplane area are all read off a cut that reuses the machinery in
:mod:`tensorcad.stability` (the same mesh-identity-keyed section segments, the same
Green's-theorem loop moments), so a hull that is watertight and Euler-2 for hydrostatics
is watertight and Euler-2 here too.

**What is measured.** For a heel angle ``phi`` the body is rotated about a longitudinal
axis through the keel (the mesh's lowest point, ``y = 0``) and cut with the horizontal
plane ``z = waterline_z`` in the *rotated* frame; the waterline is then moved until the
submerged volume equals the upright displaced volume. Writing ``KN(phi)`` for the
transverse distance from the centreline to the vertical through the centre of buoyancy,

    GZ(phi) = KN(phi) - KG * sin(phi)          (G on the centreline, at height KG)

which is the standard cross-curve form: ``KN`` is a property of the shape alone, ``KG`` is
a loading input, and the ``sin`` factor is the transverse offset of G itself (G, being on
the centreline, swings out by ``KG sin(phi)`` when the hull heels). Positive heel means
starboard down and positive ``GZ`` means righting, and both conventions are asserted
rather than assumed -- see the mirrored-hull test and ``examples/gz``.

**How the cut is integrated.** A triangle is classified by the signed height of its three
vertices; wholly-below triangles are kept, wholly-above ones dropped, and crossing ones
clipped to a polygon and fan-triangulated. Faces lying *in* the waterplane are dropped
rather than counted -- their boundary is recovered from the neighbouring faces -- and the
waterplane cap is then added back from the section loops, with the largest loop taken as
the outer boundary and any loop inside it reversed as a hole. Volume and centroid are the
tetrahedra ``(origin, a, b, c)`` summed with their signs, so the answer is exact for the
polygon soup the mesh actually is; the only approximation is the one already baked into
the mesh. The fan triangulation of a non-convex loop is safe here because the tetrahedron
sum of a fan reproduces the boundary integral exactly whether or not the triangles
overlap -- their overlap cancels.

**What it does not do** (all of it stated in ``docs/gz.md`` too): free trim at heel (the
waterline is kept parallel to the still-water line, so the reported arm is a
*constant-trim* cross curve -- the classical one, and slightly conservative for a hull
that would trim aft as it heels); free-surface corrections; grain/heeling-moment allowances;
damage or floodable-length calculations; and any class-rule approval. The IMO check here
is a *screening* of the 2008 IS Code Part A 2.2 numerical criteria against a measured
curve, not a stability booklet.
"""

from __future__ import annotations

import itertools
import math
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

from . import mesh as meshlib
from .geometry import Mesh
from .mesh import total_signed_volume
from .stability import (
    PLANE_TOLERANCE,
    StabilityError,
    _bounding_scale,
    waterplane_moments,
)

SCHEMA_HEELED = "tensorcad.stability.heeled/1"
SCHEMA_RIGHTING_ARM = "tensorcad.stability.righting_arm/1"

#: Heel angles a righting-arm curve is sampled at unless the caller says otherwise.
DEFAULT_ANGLES_DEG: tuple[float, ...] = tuple(float(angle) for angle in range(0, 91, 5))


class GzError(ValueError):
    """Raised when a heeled integration cannot be set up honestly."""


@dataclass(frozen=True)
class HeeledHydrostatics:
    """Displaced volume and centre of buoyancy of a body heeled by one angle.

    ``volume_error`` is ``(V(phi) - V(0)) / V(0)`` after the waterline solve. A curve whose
    points sit at different displacements is not a righting-arm curve, so the residual is
    reported per point rather than assumed away.
    """

    heel_angle_deg: float
    displacement_volume: float
    target_volume: float
    volume_error: float
    waterline_offset: float
    waterplane_area: float
    wetted_surface_area: float
    centre_of_buoyancy_y: float
    centre_of_buoyancy_z: float
    kn: float
    iterations: int

    def to_dict(self) -> dict[str, object]:
        return {
            "heel_angle_deg": _round(self.heel_angle_deg),
            "displacement_volume_m3": _round(self.displacement_volume),
            "target_volume_m3": _round(self.target_volume),
            "volume_error": _round(self.volume_error, 15),
            "waterline_offset_m": _round(self.waterline_offset),
            "waterplane_area_m2": _round(self.waterplane_area),
            "wetted_surface_area_m2": _round(self.wetted_surface_area),
            "centre_of_buoyancy_y_m": _round(self.centre_of_buoyancy_y),
            "centre_of_buoyancy_z_m": _round(self.centre_of_buoyancy_z),
            "kn_m": _round(self.kn),
            "iterations": self.iterations,
        }


@dataclass(frozen=True)
class RightingArmCurve:
    """``GZ`` sampled over heel angles, with the summary numbers a curve is read for."""

    angles_deg: tuple[float, ...]
    gz: tuple[float, ...]
    kn: tuple[float, ...]
    points: tuple[HeeledHydrostatics, ...]
    kg: float
    free_surface_correction: float
    waterline_z: float
    keel_z: float
    displacement_volume: float

    @property
    def max_volume_error(self) -> float:
        """Worst ``|V(phi)/V(0) - 1|`` over the curve: the displacement residual."""
        return max(abs(point.volume_error) for point in self.points)

    @property
    def max_gz(self) -> float:
        return max(self.gz)

    @property
    def angle_of_max_gz_deg(self) -> float:
        return self.angles_deg[self.gz.index(max(self.gz))]

    @property
    def range_of_stability_deg(self) -> float | None:
        """First angle past the peak where ``GZ`` returns to zero, interpolated.

        ``None`` when the sampled range ends while the body is still stable: that is
        \"beyond the range sampled\", not \"infinite\", and the caller should widen it.
        """
        peak = self.gz.index(max(self.gz))
        for index in range(peak + 1, len(self.gz)):
            if self.gz[index] <= 0.0:
                previous_angle = self.angles_deg[index - 1]
                previous_gz = self.gz[index - 1]
                span = previous_gz - self.gz[index]
                if span <= 0.0:  # pragma: no cover - degenerate flat segment
                    return previous_angle
                return previous_angle + (self.angles_deg[index] - previous_angle) * (previous_gz / span)
        return None

    def area_under_curve(self, limit_deg: float | None = None) -> float:
        """Dynamical stability ``int GZ dphi`` in m*rad, trapezoidal on the sample grid.

        Truncated at ``limit_deg`` when given (interpolating ``GZ`` there), which is the
        form the number is quoted in -- an area to 30 or 40 degrees, not to vanishing GZ.
        """
        limit = self.angles_deg[-1] if limit_deg is None else limit_deg
        if limit <= self.angles_deg[0]:
            return 0.0
        total = 0.0
        for index in range(len(self.angles_deg) - 1):
            left, right = self.angles_deg[index], self.angles_deg[index + 1]
            if left >= limit:
                break
            gz_left, gz_right = self.gz[index], self.gz[index + 1]
            if right > limit:
                fraction = (limit - left) / (right - left)
                gz_right = gz_left + fraction * (gz_right - gz_left)
            total += 0.5 * (gz_left + gz_right) * math.radians(min(right, limit) - left)
        return total

    def to_dict(self) -> dict[str, object]:
        limit = 30.0
        return {
            "schema": SCHEMA_RIGHTING_ARM,
            "heel_angles_deg": list(self.angles_deg),
            "gz_m": [_round(value) for value in self.gz],
            "kn_m": [_round(value) for value in self.kn],
            "kg_m": _round(self.kg),
            "keel_z_m": _round(self.keel_z),
            "free_surface_correction_m": _round(self.free_surface_correction),
            "waterline_z_m": _round(self.waterline_z),
            "displacement_volume_m3": _round(self.displacement_volume),
            "max_gz_m": _round(self.max_gz),
            "angle_of_max_gz_deg": _round(self.angle_of_max_gz_deg),
            "range_of_stability_deg": (
                None if self.range_of_stability_deg is None else _round(self.range_of_stability_deg)
            ),
            "area_under_curve_m_rad": _round(self.area_under_curve(limit)),
            "area_under_curve_limit_deg": limit,
            "max_volume_error": self.max_volume_error,
            "points": [point.to_dict() for point in self.points],
        }


@dataclass(frozen=True)
class FreeSurfaceCorrection:
    """The flat ``sum(i_x) / V`` reduction of a solid righting-arm curve."""

    inertias: tuple[float, ...]
    total_inertia: float
    displacement_volume: float
    loss_m: float

    def to_dict(self) -> dict[str, object]:
        return {
            "tank_inertias_m4": [_round(value) for value in self.inertias],
            "total_inertia_m4": _round(self.total_inertia),
            "displacement_volume_m3": _round(self.displacement_volume),
            "loss_m": _round(self.loss_m),
        }


@dataclass(frozen=True)
class _Integral:
    """Volume, first moments and areas of the part of a surface below a plane."""

    volume: float
    moment: tuple[float, float, float]
    waterplane_area: float
    wetted_area: float


def _round(value: float, digits: int = 9) -> float | None:
    if not math.isfinite(value):
        return None
    return round(value, digits)


def _rotate_vertices(
    vertices: Sequence[tuple[float, float, float]], heel_rad: float
) -> list[tuple[float, float, float]]:
    """Heel about ``+x``: the ``+y`` side goes down for positive angles."""
    cos_heel, sin_heel = math.cos(heel_rad), math.sin(heel_rad)
    return [(x, y * cos_heel + z * sin_heel, -y * sin_heel + z * cos_heel) for (x, y, z) in vertices]


def _clip_below(
    a: tuple[float, float, float],
    b: tuple[float, float, float],
    c: tuple[float, float, float],
    level: float,
) -> list[tuple[float, float, float]]:
    """Sutherland-Hodgman clip of one triangle against ``z <= level``, order preserved."""
    polygon: list[tuple[float, float, float]] = []
    for start, end in ((a, b), (b, c), (c, a)):
        depth_start = start[2] - level
        depth_end = end[2] - level
        if depth_start <= 0.0:
            polygon.append(start)
        if (depth_start < 0.0 < depth_end) or (depth_end < 0.0 < depth_start):
            fraction = depth_start / (depth_start - depth_end)
            polygon.append(
                (
                    start[0] + fraction * (end[0] - start[0]),
                    start[1] + fraction * (end[1] - start[1]),
                    level,
                )
            )
    return polygon


def _integrate_below(
    rotated: Sequence[tuple[float, float, float]],
    faces: Sequence[tuple[int, int, int]],
    level: float,
) -> _Integral:
    """Integrate the region below ``z = level``, apex on the plane (so the cap is free)."""
    volumes: list[float] = []
    moments_x: list[float] = []
    moments_y: list[float] = []
    moments_z: list[float] = []
    projected: list[float] = []
    wetted: list[float] = []
    for i, j, k in faces:
        first, second, third = rotated[i], rotated[j], rotated[k]
        if first[2] == level and second[2] == level and third[2] == level:
            # A face lying *in* the waterline plane adds no volume and no moment, but it
            # would corrupt the projected-area sum that gives the waterplane area: at
            # phi = 0 that face *is* the waterplane, and counting it twice cancels it.
            continue
        polygon = _clip_below(first, second, third, level)
        for index in range(1, len(polygon) - 1):
            p, q, r = polygon[0], polygon[index], polygon[index + 1]
            # waterline-relative coordinates put the apex on the closing cap
            ax, ay, az = p[0], p[1], p[2] - level
            bx, by, bz = q[0], q[1], q[2] - level
            cx, cy, cz = r[0], r[1], r[2] - level
            volume = (ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx)) / 6.0
            volumes.append(volume)
            if volume:
                moments_x.append(volume * (ax + bx + cx) / 4.0)
                moments_y.append(volume * (ay + by + cy) / 4.0)
                moments_z.append(volume * (az + bz + cz) / 4.0)
            ux, uy = q[0] - p[0], q[1] - p[1]
            vx, vy = r[0] - p[0], r[1] - p[1]
            projected.append(0.5 * (ux * vy - uy * vx))
            wetted.append(0.5 * math.sqrt(_norm2(p, q, r)))
    return _Integral(
        volume=math.fsum(volumes),
        moment=(math.fsum(moments_x), math.fsum(moments_y), math.fsum(moments_z)),
        waterplane_area=abs(math.fsum(projected)),
        wetted_area=math.fsum(wetted),
    )


def _norm2(
    p: tuple[float, float, float], q: tuple[float, float, float], r: tuple[float, float, float]
) -> float:
    """Squared norm of the triangle's area vector (twice the area, as a vector)."""
    ux, uy, uz = q[0] - p[0], q[1] - p[1], q[2] - p[2]
    vx, vy, vz = r[0] - p[0], r[1] - p[1], r[2] - p[2]
    cx = uy * vz - uz * vy
    cy = uz * vx - ux * vz
    cz = ux * vy - uy * vx
    return cx * cx + cy * cy + cz * cz


def _solve_level(
    rotated: Sequence[tuple[float, float, float]],
    faces: Sequence[tuple[int, int, int]],
    target: float,
    start: float,
    lower: float,
    upper: float,
    tolerance: float,
) -> tuple[float, _Integral, int]:
    """Find the waterline height that restores the displaced volume, by Newton/bisection."""
    level = min(max(start, lower), upper)
    integral = _integrate_below(rotated, faces, level)
    for iteration in range(1, 81):
        residual = integral.volume - target
        if abs(residual) <= tolerance:
            return level, integral, iteration
        if residual > 0.0:
            upper = level
        else:
            lower = level
        step = -residual / integral.waterplane_area if integral.waterplane_area > 0.0 else 0.0
        candidate = level + step
        if not lower < candidate < upper:
            candidate = 0.5 * (lower + upper)
        if abs(candidate - level) <= 1e-16 * max(1.0, abs(level)):
            return level, integral, iteration
        level = candidate
        integral = _integrate_below(rotated, faces, level)
    return level, integral, 80


def _shifted(mesh: Mesh, keel_z: float) -> list[tuple[float, float, float]]:
    """Vertex list with ``z = 0`` on the keel, so levers are measured as naval architects do."""
    return [(x, y, z - keel_z) for (x, y, z) in mesh.vertices]


def _upright_volume(
    vertices: Sequence[tuple[float, float, float]],
    faces: Sequence[tuple[int, int, int]],
    level: float,
) -> float:
    return _integrate_below(vertices, faces, level).volume


def _keel_of(mesh: Mesh, keel_z: float | None) -> float:
    if keel_z is None:
        return mesh.bounding_box().minimum[2]
    if not math.isfinite(keel_z):
        raise GzError(f"keel_z must be finite: {keel_z!r}")
    return float(keel_z)


def heeled_hydrostatics(
    mesh: Mesh,
    heel_deg: float,
    *,
    waterline_z: float = 0.0,
    keel_z: float | None = None,
    displacement_volume: float | None = None,
    tolerance: float = 1e-13,
) -> HeeledHydrostatics:
    """Integrate the displaced volume of ``mesh`` heeled by ``heel_deg`` degrees.

    ``displacement_volume`` defaults to the upright displaced volume at ``waterline_z``;
    pass a previously computed value to keep every point of a curve on one displacement
    instead of re-deriving it per angle. The mesh must be closed and must extend above the
    waterline over the whole range of heel asked for -- :func:`with_freeboard` adds that
    band to a hull that stops at its design waterline.
    """
    if not math.isfinite(heel_deg) or abs(heel_deg) >= 180.0:
        raise GzError(f"heel angle must be finite and within (-180, 180) degrees: {heel_deg!r}")
    if not math.isfinite(waterline_z):
        raise GzError(f"waterline_z must be finite: {waterline_z!r}")
    if not meshlib.is_watertight(mesh):
        raise GzError(
            "a heeled displacement needs a closed mesh: an open surface has no inside, "
            "so the displaced volume is not defined (see mesh.weld/orient_outward)"
        )
    keel = _keel_of(mesh, keel_z)
    if waterline_z <= keel:
        raise GzError(
            f"the waterline ({waterline_z!r}) must sit above the keel datum ({keel!r}): "
            "otherwise there is no displaced volume to heel"
        )
    vertices = _shifted(mesh, keel)
    base_level = waterline_z - keel
    target = (
        _upright_volume(vertices, mesh.faces, base_level)
        if displacement_volume is None
        else float(displacement_volume)
    )
    if not math.isfinite(target) or target <= 0.0:
        raise GzError(f"the displaced volume at the waterline must be positive: {target!r}")
    if tolerance <= 0.0:
        raise GzError(f"tolerance must be positive: {tolerance!r}")

    heel_rad = math.radians(heel_deg)
    rotated = _rotate_vertices(vertices, heel_rad)
    heights = [point[2] for point in rotated]
    level, integral, iterations = _solve_level(
        rotated,
        mesh.faces,
        target,
        base_level,
        min(heights),
        max(heights),
        tolerance * target,
    )
    volume = integral.volume
    if volume <= 0.0:  # pragma: no cover - a closed mesh above the plane only
        raise GzError("the solved waterline leaves no displaced volume")
    # centroid in waterline-relative coordinates, then in the body frame from the keel
    y_world = integral.moment[1] / volume
    z_world = integral.moment[2] / volume + level
    cos_heel, sin_heel = math.cos(heel_rad), math.sin(heel_rad)
    # rotate the world centroid back into the (keel-datum) body frame
    centre_y = y_world * cos_heel - z_world * sin_heel
    centre_z = y_world * sin_heel + z_world * cos_heel
    return HeeledHydrostatics(
        heel_angle_deg=heel_deg,
        displacement_volume=volume,
        target_volume=target,
        volume_error=(volume - target) / target,
        waterline_offset=level - base_level,
        waterplane_area=integral.waterplane_area,
        wetted_surface_area=integral.wetted_area,
        centre_of_buoyancy_y=centre_y,
        centre_of_buoyancy_z=centre_z,
        kn=y_world,
        iterations=iterations,
    )


def righting_arm_curve(
    mesh: Mesh,
    kg: float,
    angles_deg: Sequence[float] = DEFAULT_ANGLES_DEG,
    *,
    waterline_z: float = 0.0,
    keel_z: float | None = None,
    displacement_volume: float | None = None,
    free_surface_i_x: float = 0.0,
    tolerance: float = 1e-13,
) -> RightingArmCurve:
    """Sample ``GZ`` over ``angles_deg`` for a centre of gravity ``kg`` above the keel.

    ``free_surface_i_x`` is the sum of the free-surface second moments ``i_x`` of any slack
    tanks; it enters as the flat reduction ``sum(i_x)/V`` of the whole curve, which is what
    that correction is.
    """
    if not math.isfinite(kg):
        raise GzError(f"kg must be finite: {kg!r}")
    if not math.isfinite(free_surface_i_x) or free_surface_i_x < 0.0:
        raise GzError(f"free_surface_i_x must be finite and non-negative: {free_surface_i_x!r}")
    angles = tuple(float(angle) for angle in angles_deg)
    if not angles:
        raise GzError("a curve needs at least one heel angle")
    if len(set(angles)) != len(angles):
        raise GzError(f"heel angles must be distinct: {angles!r}")

    keel = _keel_of(mesh, keel_z)
    if kg < keel:
        raise GzError(f"kg ({kg!r}) is below the keel datum ({keel!r})")
    target = (
        _upright_volume(_shifted(mesh, keel), mesh.faces, waterline_z - keel)
        if displacement_volume is None
        else float(displacement_volume)
    )
    points = tuple(
        heeled_hydrostatics(
            mesh,
            angle,
            waterline_z=waterline_z,
            keel_z=keel,
            displacement_volume=target,
            tolerance=tolerance,
        )
        for angle in angles
    )
    correction = free_surface_i_x / target
    # ``kg`` is measured above the keel and the levers ``kn`` are taken about the keel
    # line (the rotation axis), so the arm is simply ``kn - kg sin(phi)``.
    gz = tuple(point.kn - kg * math.sin(math.radians(point.heel_angle_deg)) - correction for point in points)
    return RightingArmCurve(
        angles_deg=angles,
        gz=gz,
        kn=tuple(point.kn for point in points),
        points=points,
        kg=kg,
        free_surface_correction=correction,
        waterline_z=waterline_z,
        keel_z=keel,
        displacement_volume=target,
    )


def wall_sided_righting_arm(angle_deg: float, *, gm: float, bm: float) -> float:
    """The classical wall-sided righting arm ``sin(phi) (GM + BM/2 tan^2 phi)``."""
    if not math.isfinite(angle_deg) or abs(angle_deg) >= 90.0:
        raise GzError(f"the wall-sided form needs |heel| < 90 degrees: {angle_deg!r}")
    heel_rad = math.radians(angle_deg)
    return math.sin(heel_rad) * (gm + 0.5 * bm * math.tan(heel_rad) ** 2)


def box_wall_sided_limit_deg(breadth: float, draft: float, depth: float | None = None) -> float:
    """Heel at which a box stops being wall-sided: bilge emergence or deck-edge immersion."""
    if breadth <= 0.0 or draft <= 0.0:
        raise GzError("breadth and draft must be positive")
    limits = [math.degrees(math.atan(2.0 * draft / breadth))]
    if depth is not None:
        if depth <= draft:
            raise GzError(f"depth {depth!r} must exceed draft {draft!r}")
        limits.append(math.degrees(math.atan(2.0 * (depth - draft) / breadth)))
    return min(limits)


def box_righting_arms_analytic(
    breadth: float,
    draft: float,
    kg: float,
    angles_deg: Sequence[float],
    *,
    depth: float | None = None,
) -> tuple[float, ...]:
    """Wall-sided reference for a box barge: ``KB0 = T/2``, ``BM0 = B^2/(12T)``."""
    if breadth <= 0.0 or draft <= 0.0:
        raise GzError("breadth and draft must be positive")
    bm = breadth * breadth / (12.0 * draft)
    gm = 0.5 * draft + bm - kg
    if depth is not None:
        limit = box_wall_sided_limit_deg(breadth, draft, depth)
        for angle in angles_deg:
            if abs(angle) > limit:
                raise GzError(
                    f"heel {angle!r} exceeds the wall-sided limit {limit:.3f} degrees: the "
                    "waterline has reached the bilge or the deck edge, so the closed form "
                    "no longer applies"
                )
    return tuple(wall_sided_righting_arm(angle, gm=gm, bm=bm) for angle in angles_deg)


def rectangular_tank_inertia(length: float, breadth: float) -> float:
    """Second moment ``i_x = l b^3 / 12`` of a rectangular free surface about its long axis."""
    if length <= 0.0 or breadth <= 0.0:
        raise GzError("tank length and breadth must be positive")
    return length * breadth**3 / 12.0


def free_surface_correction(
    tank_inertias: Sequence[float], displacement_volume: float
) -> FreeSurfaceCorrection:
    """``sum(i_x) / V``: the flat amount a slack tank takes off every point of ``GZ``."""
    if displacement_volume <= 0.0:
        raise GzError(f"displacement_volume must be positive: {displacement_volume!r}")
    inertias = tuple(float(value) for value in tank_inertias)
    for value in inertias:
        if not math.isfinite(value) or value < 0.0:
            raise GzError(f"a free-surface inertia must be finite and non-negative: {value!r}")
    total = math.fsum(inertias)
    return FreeSurfaceCorrection(
        inertias=inertias,
        total_inertia=total,
        displacement_volume=displacement_volume,
        loss_m=total / displacement_volume,
    )


def with_freeboard(mesh: Mesh, freeboard: float, *, waterline_z: float = 0.0) -> Mesh:
    """Close a hull that stops at its waterline with a vertical freeboard band and a deck.

    A hull from :mod:`tensorcad.hull` or :mod:`tensorcad.wigley` is modelled only up to the
    design waterline, so heeling it lets the waterline leave the surface on the dipping side
    and the displaced volume silently stops meaning anything. This removes the waterline cap,
    extrudes the resulting boundary ring upwards by ``freeboard`` with vertical walls, and
    closes it with a fan. The upright displacement at ``waterline_z`` is unchanged, which the
    tests assert to round-off; the heel range that becomes valid is bounded by deck-edge
    immersion, :func:`deck_edge_immersion_deg`.
    """
    if not math.isfinite(freeboard) or freeboard <= 0.0:
        raise GzError(f"freeboard must be finite and positive: {freeboard!r}")
    if not mesh.faces:
        raise GzError("cannot add freeboard to an empty mesh")
    box = mesh.bounding_box()
    tolerance = 1e-9 * max(box.sizes) if max(box.sizes) > 0.0 else 1e-12
    cap = [
        face
        for face in mesh.faces
        if all(abs(mesh.vertices[index][2] - waterline_z) <= tolerance for index in face)
    ]
    if not cap:
        raise GzError(
            "no faces lie on the waterline: with_freeboard needs the waterline cap a "
            "generated hull is closed with"
        )
    cap_set = set(cap)
    kept = [face for face in mesh.faces if face not in cap_set]
    counts: dict[tuple[int, int], int] = defaultdict(int)
    for face in kept:
        for first, second in ((face[0], face[1]), (face[1], face[2]), (face[2], face[0])):
            counts[(min(first, second), max(first, second))] += 1
    boundary = [edge for edge, count in counts.items() if count == 1]
    if not boundary:
        raise GzError("the waterline cap left no open boundary to extrude")
    neighbours: dict[int, list[int]] = defaultdict(list)
    for first, second in boundary:
        neighbours[first].append(second)
        neighbours[second].append(first)
    start = boundary[0][0]
    loop = [start]
    previous: int | None = None
    current = start
    while True:
        options = [vertex for vertex in neighbours[current] if vertex != previous]
        if not options:
            raise GzError("the waterline boundary is not a closed ring")
        following = options[0]
        if following == start:
            break
        if following in loop:
            raise GzError("the waterline boundary revisits a vertex (not a simple ring)")
        loop.append(following)
        previous, current = current, following
    if len(loop) != len(boundary):
        raise GzError("the waterline has more than one open ring; only one hull ring is supported")

    top = waterline_z + freeboard
    vertices = list(mesh.vertices)
    extruded = []
    for index in loop:
        x, y, _ = mesh.vertices[index]
        extruded.append(len(vertices))
        vertices.append((x, y, top))
    faces = list(kept)
    for position, index in enumerate(loop):
        lower = index
        lower_next = loop[(position + 1) % len(loop)]
        upper = extruded[position]
        upper_next = extruded[(position + 1) % len(loop)]
        faces.append((lower, lower_next, upper_next))
        faces.append((lower, upper_next, upper))
    mean_x = math.fsum(vertices[index][0] for index in loop) / len(loop)
    mean_y = math.fsum(vertices[index][1] for index in loop) / len(loop)
    apex = len(vertices)
    vertices.append((mean_x, mean_y, top))
    for position in range(len(loop)):
        faces.append((apex, extruded[position], extruded[(position + 1) % len(loop)]))

    built = Mesh(
        vertices=tuple(vertices),
        faces=tuple(faces),
        units=mesh.units,
        coordinate_frame=mesh.coordinate_frame,
    )
    oriented = meshlib.orient_outward(built)
    return oriented.mesh


def deck_edge_immersion_deg(half_breadth_at_waterline: float, freeboard: float) -> float:
    """Heel at which the deck edge reaches the waterline, ``atan(freeboard / (B/2))``."""
    if half_breadth_at_waterline <= 0.0:
        raise GzError("the waterline half-breadth must be positive")
    if freeboard <= 0.0:
        raise GzError("the freeboard must be positive")
    return math.degrees(math.atan(freeboard / half_breadth_at_waterline))

# ---------------------------------------------------------------------------
# the second route (``submerged_solid`` .. ``imo_intact_stability``)
# ---------------------------------------------------------------------------
#: The two routes raise their own error type (``GzError`` for the cross-curve route,
#: ``GZError`` for the ``gz_curve`` route below); both mean the same thing -- a righting
#: arm that cannot be defined honestly -- and both are exported.
SCHEMA_SUBMERGED = "tensorcad.gz.submerged/1"
SCHEMA_GZ_CURVE = "tensorcad.gz.curve/1"
SCHEMA_IMO = "tensorcad.gz.imo_intact/1"

#: A heel grid that carries every angle the IMO criteria are evaluated at (30 and 40
#: degrees included, so ``area(30, 40)`` needs no extrapolation).
DEFAULT_HEEL_ANGLES_DEG: tuple[float, ...] = (
    0.0,
    5.0,
    10.0,
    15.0,
    20.0,
    25.0,
    30.0,
    35.0,
    40.0,
    45.0,
    50.0,
    60.0,
    70.0,
    80.0,
    90.0,
)

Vec3 = tuple[float, float, float]


class GZError(ValueError):
    """Raised when a righting arm cannot be defined for the geometry asked for."""


def _gz_round(value: float, digits: int = 9) -> float:
    return round(value + 0.0, digits)


# ---------------------------------------------------------------------------
# the submerged solid at one waterline
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class SubmergedSolid:
    """Volume, centre and waterplane area of the part of a mesh below ``waterline_z``.

    ``centroid`` is the centre of buoyancy in the frame the mesh was given in, so for a
    heeled mesh its ``y`` component is ``KN``; its ``z`` component is the height of B in
    the same frame (useful as an independent check against ``hydrostatics`` at zero heel).
    """

    volume: float
    centroid: Vec3
    waterplane_area: float
    waterline_z: float
    section_segments: int
    faces_below: int
    faces_clipped: int
    faces_coplanar: int
    schema: str = SCHEMA_SUBMERGED
    basis: str = "clipped surface plus waterplane cap, tetrahedra from the origin"

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": self.schema,
            "volume_m3": _gz_round(self.volume),
            "centroid_m": tuple(_gz_round(value) for value in self.centroid),
            "waterplane_area_m2": _gz_round(self.waterplane_area),
            "waterline_z_m": _gz_round(self.waterline_z),
            "section_segments": self.section_segments,
            "faces_below": self.faces_below,
            "faces_clipped": self.faces_clipped,
            "faces_coplanar": self.faces_coplanar,
            "basis": self.basis,
        }


def _clip_polygon_below(points: list[Vec3], waterline_z: float, tolerance: float) -> list[Vec3]:
    """Sutherland-Hodgman clip of one polygon against ``z <= waterline_z``."""
    clipped: list[Vec3] = []
    count = len(points)
    for index in range(count):
        current = points[index]
        following = points[(index + 1) % count]
        d_current = current[2] - waterline_z
        d_following = following[2] - waterline_z
        inside_current = d_current <= tolerance
        inside_following = d_following <= tolerance
        if inside_current:
            clipped.append(current)
        if inside_current != inside_following:
            t = d_current / (d_current - d_following)
            clipped.append(
                (
                    current[0] + t * (following[0] - current[0]),
                    current[1] + t * (following[1] - current[1]),
                    waterline_z,
                )
            )
    return clipped


def submerged_solid(mesh: Mesh, waterline_z: float = 0.0) -> SubmergedSolid:
    """Volume, centre of buoyancy and waterplane area below ``waterline_z``.

    The surface is clipped rather than sampled, so at a waterline that coincides with the
    mesh's own cap (a hull generated at its design waterline) the answer is the mesh's
    exact volume instead of one triangle's worth short.

    The waterplane cap is *not* triangulated, and the segments are not chained into loops.
    Each crossing triangle's outward normal gives the section direction ``z_hat x n``,
    which walks the material counter-clockwise seen from above, so the signed sum over the
    oriented segments is the section area with holes already subtracted. The cap's
    contribution to the volume integrals is then closed form: ``2*wl*A`` to the six-volume
    sum and ``12*wl^2*A`` to the ``z`` moment, nothing to ``x`` and ``y`` because the cap
    normal has no transverse component. This is exact for a concave section too, and it
    does not depend on chaining segments into loops -- which a facetted surface can defeat,
    because two faces meeting at a shallow angle key their crossing points differently even
    when the points coincide.
    """
    if not math.isfinite(waterline_z):
        raise GZError(f"waterline_z must be finite, got {waterline_z!r}")
    tolerance = PLANE_TOLERANCE * _bounding_scale(mesh)
    triangles: list[tuple[Vec3, Vec3, Vec3]] = []
    coplanar: list[tuple[Vec3, Vec3, Vec3]] = []
    faces_below = 0
    faces_clipped = 0
    section_area = 0.0
    section_first_x = 0.0
    section_first_y = 0.0
    section_segments = 0
    oriented: list[tuple[Vec3, Vec3]] = []

    for face in mesh.faces:
        points = [mesh.vertices[index] for index in face]
        heights = [point[2] - waterline_z for point in points]
        codes = [
            0 if abs(height) <= tolerance else (1 if height > 0.0 else -1) for height in heights
        ]
        if codes[0] == 0 and codes[1] == 0 and codes[2] == 0:
            coplanar.append((points[0], points[1], points[2]))
            continue
        has_above = any(code > 0 for code in codes)
        has_below = any(code < 0 for code in codes)
        if has_above and has_below:
            clipped = _clip_polygon_below(points, waterline_z, tolerance)
            faces_clipped += 1
            for corner in range(1, len(clipped) - 1):
                triangles.append((clipped[0], clipped[corner], clipped[corner + 1]))
            crossings: list[Vec3] = []
            for first, second in ((0, 1), (1, 2), (2, 0)):
                if abs(heights[first]) <= tolerance:
                    crossings.append(points[first])
                if heights[first] * heights[second] < 0.0:
                    weight = heights[first] / (heights[first] - heights[second])
                    crossings.append(
                        (
                            points[first][0] + weight * (points[second][0] - points[first][0]),
                            points[first][1] + weight * (points[second][1] - points[first][1]),
                            waterline_z,
                        )
                    )
            unique: list[Vec3] = []
            for point in crossings:
                if all(
                    abs(point[0] - other[0]) > tolerance or abs(point[1] - other[1]) > tolerance
                    for other in unique
                ):
                    unique.append(point)
            if len(unique) != 2:  # pragma: no cover - a triangle meets a plane in <= 2 points
                continue
            edge = (points[1][0] - points[0][0], points[1][1] - points[0][1], points[1][2] - points[0][2])
            other = (
                points[2][0] - points[0][0],
                points[2][1] - points[0][1],
                points[2][2] - points[0][2],
            )
            normal = (
                edge[1] * other[2] - edge[2] * other[1],
                edge[2] * other[0] - edge[0] * other[2],
                edge[0] * other[1] - edge[1] * other[0],
            )
            towards = (-normal[1], normal[0])
            head, tail = unique
            if (tail[0] - head[0]) * towards[0] + (tail[1] - head[1]) * towards[1] < 0.0:
                head, tail = tail, head
            oriented.append((head, tail))
            section_area += 0.5 * (head[0] * tail[1] - tail[0] * head[1])
            section_first_x += (
                (tail[1] - head[1])
                * (head[0] * head[0] + head[0] * tail[0] + tail[0] * tail[0])
                / 6.0
            )
            section_first_y -= (
                (tail[0] - head[0])
                * (head[1] * head[1] + head[1] * tail[1] + tail[1] * tail[1])
                / 6.0
            )
            section_segments += 1
            continue
        if has_above:
            continue  # wholly above the waterline
        faces_below += 1
        triangles.append((points[0], points[1], points[2]))

    if section_segments == 0 and coplanar:
        # The waterline coincides with the mesh's own cap: those faces *are* the section,
        # already outward-facing, so they close the solid with no extra area.
        triangles.extend(coplanar)
        for cap_a, cap_b, cap_c in coplanar:
            edge = (cap_b[0] - cap_a[0], cap_b[1] - cap_a[1], cap_b[2] - cap_a[2])
            other = (cap_c[0] - cap_a[0], cap_c[1] - cap_a[1], cap_c[2] - cap_a[2])
            section_area += 0.5 * (edge[0] * other[1] - edge[1] * other[0])

    six_volume = 0.0
    moment = [0.0, 0.0, 0.0]
    for a, b, c in triangles:
        cross = (
            b[1] * c[2] - b[2] * c[1],
            b[2] * c[0] - b[0] * c[2],
            b[0] * c[1] - b[1] * c[0],
        )
        six = a[0] * cross[0] + a[1] * cross[1] + a[2] * cross[2]
        six_volume += six
        for axis in range(3):
            moment[axis] += six * (a[axis] + b[axis] + c[axis])

    if section_segments:
        # A section is a closed curve: every section point is shared by exactly two segments
        # (collinear interior points included). An odd count means the plane ran along mesh
        # edges, where the cut is degenerate and the area it would report belongs to nothing.
        hair = max(tolerance, 1e-15)
        counts: dict[tuple[int, int], int] = {}
        for head, tail in oriented:
            for point in (head, tail):
                key = (round(point[0] / hair), round(point[1] / hair))
                counts[key] = counts.get(key, 0) + 1
        if section_area == 0.0 or any(count % 2 for count in counts.values()):
            raise GZError(
                f"the waterline z={waterline_z!r} m meets the surface along mesh edges: the "
                f"section does not close, so it has no area and no displaced volume can be "
                f"assigned (nudge the heel angle or the waterline height off the edges)"
            )
    if section_segments:
        # The origin-cone over the section closes the open side surface: its six-volume and
        # first moments are what the tetra sum is missing.
        six_volume += 2.0 * waterline_z * section_area
        moment[0] += 6.0 * waterline_z * section_first_x
        moment[1] += 6.0 * waterline_z * section_first_y
        moment[2] += 6.0 * waterline_z * waterline_z * section_area

    if six_volume < 0.0:
        six_volume = -six_volume
        moment = [-value for value in moment]
        section_area = -section_area
    if six_volume == 0.0:
        centroid: Vec3 = (0.0, 0.0, 0.0)
        volume = 0.0
    else:
        volume = six_volume / 6.0
        scale = 4.0 * six_volume
        centroid = (moment[0] / scale, moment[1] / scale, moment[2] / scale)

    return SubmergedSolid(
        volume=volume,
        centroid=centroid,
        waterplane_area=section_area,
        waterline_z=waterline_z,
        section_segments=section_segments,
        faces_below=faces_below,
        faces_clipped=faces_clipped,
        faces_coplanar=len(coplanar),
    )


def solve_waterline(
    mesh: Mesh,
    displacement_volume: float,
    *,
    initial_z: float | None = None,
    tolerance: float = 1e-12,
    max_iterations: int = 200,
) -> tuple[float, SubmergedSolid]:
    """Move the waterline until the submerged volume equals ``displacement_volume``.

    Newton on ``V(z) - V_target`` with the waterplane area as ``dV/dz``, safeguarded by
    bisection so it cannot leave the bracket. This is the sinkage solve a heeled
    righting-arm needs and an upright hydrostatic call does not.
    """
    if not math.isfinite(displacement_volume) or displacement_volume <= 0.0:
        raise GZError(
            f"displacement volume {displacement_volume!r} m^3 is not positive; a righting "
            f"arm needs a floating body"
        )
    heights = [vertex[2] for vertex in mesh.vertices]
    if not heights:
        raise GZError("cannot float an empty mesh")
    low, high = min(heights), max(heights)
    span = max(high - low, 1.0)
    margin = 1e-9 * span
    lower, upper = low - margin, high + margin
    total = abs(total_signed_volume(mesh))
    if displacement_volume > total * (1.0 + 1e-12):
        raise GZError(
            f"displacement volume {displacement_volume!r} m^3 is not below the mesh's own "
            f"volume {total!r} m^3; no waterline exists at or above the model"
        )
    z = initial_z if initial_z is not None else 0.5 * (lower + upper)
    solid = submerged_solid(mesh, z)
    for _ in range(max_iterations):
        difference = solid.volume - displacement_volume
        if abs(difference) <= tolerance * displacement_volume:
            return z, solid
        if difference < 0.0:
            lower = z
        else:
            upper = z
        step = solid.waterplane_area
        candidate = z - difference / step if step > 0.0 else math.nan
        if not (lower < candidate < upper):
            candidate = 0.5 * (lower + upper)
        if not math.isfinite(candidate) or abs(candidate - z) <= 1e-15 * span:
            return z, solid
        z = candidate
        solid = submerged_solid(mesh, z)
    raise GZError(
        f"waterline solve did not converge to {displacement_volume!r} m^3 in "
        f"{max_iterations} iterations"
    )


# ---------------------------------------------------------------------------
# heeling
# ---------------------------------------------------------------------------
def heel_mesh(mesh: Mesh, heel_angle_deg: float, *, pivot_z: float | None = None) -> Mesh:
    """Rotate the model about a longitudinal axis through the keel centreline.

    Positive ``heel_angle_deg`` puts the **starboard** side down, i.e. the point ``(y>0,
    z=pivot_z)`` is driven below the pivot and the point ``(0, pivot_z + h)`` swings to
    ``y = +h sin(phi)``. That is the convention under which ``GZ = KN - KG sin(phi)`` comes
    out positive for a stable hull, and the mirrored-hull test in ``tests/test_gz.py``
    pins the sign.
    """
    if not math.isfinite(heel_angle_deg):
        raise GZError(f"heel_angle_deg must be finite, got {heel_angle_deg!r}")
    angle = math.radians(heel_angle_deg)
    cosine, sine = math.cos(angle), math.sin(angle)
    pivot = min(vertex[2] for vertex in mesh.vertices) if pivot_z is None else pivot_z
    vertices = tuple(
        (x, y * cosine + (z - pivot) * sine, -y * sine + (z - pivot) * cosine + pivot)
        for x, y, z in mesh.vertices
    )
    return Mesh(
        vertices=vertices,
        faces=mesh.faces,
        units=mesh.units,
        coordinate_frame=mesh.coordinate_frame,
    )


@dataclass(frozen=True)
class HeelPoint:
    """One point of the cross curve: the shape's ``KN`` and the loading's ``GZ``."""

    heel_angle_deg: float
    kn: float
    kg: float
    gz: float
    displacement_volume: float
    waterline_z: float
    centre_buoyancy: Vec3
    waterplane_area: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "heel_angle_deg": _gz_round(self.heel_angle_deg, 6),
            "kn_m": _gz_round(self.kn),
            "kg_m": _gz_round(self.kg),
            "gz_m": _gz_round(self.gz),
            "displacement_volume_m3": _gz_round(self.displacement_volume),
            "waterline_z_m": _gz_round(self.waterline_z),
            "centre_buoyancy_m": tuple(_gz_round(value) for value in self.centre_buoyancy),
            "waterplane_area_m2": _gz_round(self.waterplane_area),
        }


@dataclass(frozen=True)
class GZCurve:
    """A constant-trim righting-arm curve plus the loading it was cut for."""

    points: tuple[HeelPoint, ...]
    kg: float
    keel_z: float
    displacement_volume: float
    displacement_mass: float
    water_density: float
    gm_upright: float
    upright_waterline_z: float
    schema: str = SCHEMA_GZ_CURVE
    basis: str = "constant-trim cross curves, clipped sections, GZ = KN - KG sin(phi)"

    @property
    def angles(self) -> tuple[float, ...]:
        return tuple(point.heel_angle_deg for point in self.points)

    def gz_at(self, heel_angle_deg: float) -> float:
        """``GZ`` interpolated linearly between the sampled angles (no extrapolation)."""
        angles = self.angles
        if not angles:
            raise GZError("the curve has no points")
        if heel_angle_deg < angles[0] or heel_angle_deg > angles[-1]:
            raise GZError(
                f"heel angle {heel_angle_deg!r} deg is outside the sampled range "
                f"[{angles[0]!r}, {angles[-1]!r}]; add it to heel_angles_deg"
            )
        for index in range(len(angles) - 1):
            if angles[index] <= heel_angle_deg <= angles[index + 1]:
                left, right = self.points[index], self.points[index + 1]
                if right.heel_angle_deg == left.heel_angle_deg:
                    return left.gz
                weight = (heel_angle_deg - left.heel_angle_deg) / (
                    right.heel_angle_deg - left.heel_angle_deg
                )
                return left.gz + weight * (right.gz - left.gz)
        return self.points[-1].gz

    def area(self, start_deg: float = 0.0, stop_deg: float | None = None) -> float:
        """Area under the curve between two heel angles, in **m·rad**.

        The curve is piecewise linear between the sampled angles, so the trapezoid rule on
        the segments -- split exactly at ``start_deg`` and ``stop_deg`` -- is the integral,
        not an approximation of it. Returns a signed value: area below the heel axis is
        negative, which is what an unstable range looks like.
        """
        angles = self.angles
        if not angles:
            raise GZError("the curve has no points")
        stop = angles[-1] if stop_deg is None else stop_deg
        if start_deg < angles[0] or stop > angles[-1]:
            raise GZError(
                f"the interval [{start_deg!r}, {stop!r}] deg is outside the sampled range "
                f"[{angles[0]!r}, {angles[-1]!r}]; add the endpoints to heel_angles_deg"
            )
        breaks = sorted({start_deg, stop} | {angle for angle in angles if start_deg < angle < stop})
        total = 0.0
        for left, right in itertools.pairwise(breaks):
            total += 0.5 * (self.gz_at(left) + self.gz_at(right)) * math.radians(right - left)
        return total

    def max_gz(self) -> tuple[float, float]:
        """``(max GZ, angle at which it occurs)`` over the sampled angles."""
        best = max(self.points, key=lambda point: point.gz)
        return best.gz, best.heel_angle_deg

    def to_dict(self) -> dict[str, Any]:
        maximum, maximum_angle = self.max_gz()
        return {
            "schema": self.schema,
            "kg_m": _gz_round(self.kg),
            "keel_z_m": _gz_round(self.keel_z),
            "displacement_volume_m3": _gz_round(self.displacement_volume),
            "displacement_mass_kg": _gz_round(self.displacement_mass, 6),
            "water_density_kg_m3": self.water_density,
            "gm_upright_m": _gz_round(self.gm_upright),
            "upright_waterline_z_m": _gz_round(self.upright_waterline_z),
            "max_gz_m": _gz_round(maximum),
            "angle_at_max_gz_deg": _gz_round(maximum_angle, 6),
            # the IMO areas are only reported when the curve actually spans them: an area
            # extrapolated past the last sampled angle would be a number nobody measured
            "area_0_30_m_rad": _gz_round(self.area(0.0, 30.0))
            if self.angles[0] <= 0.0 <= self.angles[-1] and self.angles[-1] >= 30.0
            else None,
            "area_0_40_m_rad": _gz_round(self.area(0.0, 40.0))
            if self.angles[0] <= 0.0 <= self.angles[-1] and self.angles[-1] >= 40.0
            else None,
            "points": [point.to_dict() for point in self.points],
            "basis": self.basis,
        }


def gz_curve(
    mesh: Mesh,
    *,
    kg: float,
    heel_angles_deg: tuple[float, ...] = DEFAULT_HEEL_ANGLES_DEG,
    waterline_z: float = 0.0,
    water_density: float = 1025.0,
    displacement_volume: float | None = None,
    pivot_z: float | None = None,
    closed: bool | None = None,
) -> GZCurve:
    """The righting-arm curve of one loading condition.

    ``kg`` is the height of the centre of gravity above the keel, exactly as in
    :func:`tensorcad.stability.initial_stability`, and it is a required keyword for the
    same reason: a default of zero would silently report the bare shape's cross curve as
    a ship's stability. The displaced volume is either read off the upright waterline
    (``waterline_z``) or given explicitly, in which case the upright floatation is solved
    for first -- so the curve is always at one constant displacement.
    """
    if not math.isfinite(kg):
        raise GZError(f"kg must be finite, got {kg!r}")
    angles = tuple(float(angle) for angle in heel_angles_deg)
    if not angles:
        raise GZError("heel_angles_deg is empty")
    if any(not math.isfinite(angle) for angle in angles):
        raise GZError("heel_angles_deg contains a non-finite angle")
    if any(abs(angle) > 180.0 for angle in angles):
        raise GZError("heel angles are measured in degrees and must lie in [-180, 180]")
    if list(angles) != sorted(angles):
        raise GZError("heel_angles_deg must be increasing; the interpolation assumes it")

    keel_z = min(vertex[2] for vertex in mesh.vertices)
    pivot = keel_z if pivot_z is None else pivot_z

    if displacement_volume is None:
        upright_solid = submerged_solid(mesh, waterline_z)
        if upright_solid.waterplane_area <= 0.0:
            raise GZError(
                f"the waterline z={waterline_z!r} m gives no section; the model is above or "
                f"below the water there"
            )
        target = upright_solid.volume
        upright_waterline = waterline_z
    else:
        target = float(displacement_volume)
        upright_waterline, upright_solid = solve_waterline(mesh, target, initial_z=waterline_z)

    # GM is assembled here rather than taken from ``initial_stability`` because that
    # function integrates the *whole* closed mesh: correct for a hull capped at its design
    # waterline (where the two agree to the last bit, see the hull test), but it would use
    # the deck volume of a model that extends above the waterline. The two pieces are the
    # same ones that module uses -- the clipped volume and the waterplane's second moment.
    section = waterplane_moments(mesh, upright_waterline)
    kb = upright_solid.centroid[2] - keel_z
    bm = section.transverse_second_moment / target
    gm_upright = kb + bm - kg

    points: list[HeelPoint] = []
    for angle in angles:
        heeled = heel_mesh(mesh, angle, pivot_z=pivot)
        solved_z, solid = solve_waterline(heeled, target, initial_z=None)
        if solid.waterplane_area <= 0.0:
            raise GZError(
                f"at {angle!r} deg the waterline finds no section; the model is submerged "
                f"or clear of the water at that attitude"
            )
        kn = solid.centroid[1]
        points.append(
            HeelPoint(
                heel_angle_deg=angle,
                kn=kn,
                kg=kg,
                gz=kn - kg * math.sin(math.radians(angle)),
                displacement_volume=solid.volume,
                waterline_z=solved_z,
                centre_buoyancy=solid.centroid,
                waterplane_area=solid.waterplane_area,
            )
        )

    return GZCurve(
        points=tuple(points),
        kg=kg,
        keel_z=keel_z,
        displacement_volume=target,
        displacement_mass=target * water_density,
        water_density=water_density,
        gm_upright=gm_upright,
        upright_waterline_z=upright_waterline,
    )


# ---------------------------------------------------------------------------
# IMO intact-stability screening
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class CriterionCheck:
    """One numbered criterion, its measured value and the margin over the threshold."""

    key: str
    requirement: str
    measured: float
    unit: str
    satisfied: bool
    margin: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "key": self.key,
            "requirement": self.requirement,
            "measured": _gz_round(self.measured, 6),
            "unit": self.unit,
            "satisfied": self.satisfied,
            "margin": _gz_round(self.margin, 6),
        }


@dataclass(frozen=True)
class IMOIntactStability:
    """Six-criterion screening of a measured curve (see :func:`imo_intact_stability`)."""

    checks: tuple[CriterionCheck, ...]
    passed: bool
    gm_upright: float
    schema: str = SCHEMA_IMO
    note: str = (
        "screening against the numerical criteria of the 2008 IS Code Part A 2.2 only "
        "(no weather criterion, no free-surface or grain correction, no class approval); "
        "the curve is constant-trim and comes from a measured mesh, not a stability booklet"
    )

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": self.schema,
            "passed": self.passed,
            "gm_upright_m": _gz_round(self.gm_upright),
            "checks": [check.to_dict() for check in self.checks],
            "note": self.note,
        }


def imo_intact_stability(curve: GZCurve, *, gm_upright: float | None = None) -> IMOIntactStability:
    """Screen ``curve`` against the 2008 IS Code Part A 2.2 numerical criteria.

    The six items are: area under the curve 0-30 deg >= 0.055 m·rad, 0-40 deg >= 0.090,
    30-40 deg >= 0.030, ``GZ`` at 30 deg >= 0.20 m, maximum ``GZ`` >= 0.20 m reached at not
    less than 25 deg, and an upright ``GM`` >= 0.15 m. The curve must span 0-40 deg (the
    default grid does) because the areas are taken over the measured curve rather than
    extrapolated.
    """
    angles = curve.angles
    if not angles or angles[0] > 0.0 or angles[-1] < 40.0:
        raise GZError(
            f"the IMO criteria need a curve from 0 to at least 40 deg of heel; this one "
            f"spans [{angles[0] if angles else None}, {angles[-1] if angles else None}]"
        )
    gm = curve.gm_upright if gm_upright is None else gm_upright
    maximum, angle_at_maximum = curve.max_gz()
    area_0_30 = curve.area(0.0, 30.0)
    area_0_40 = curve.area(0.0, 40.0)
    area_30_40 = curve.area(30.0, 40.0)
    gz_30 = curve.gz_at(30.0)
    checks = (
        CriterionCheck(
            "area_0_30", ">= 0.055 m rad", area_0_30, "m rad", area_0_30 >= 0.055, area_0_30 - 0.055
        ),
        CriterionCheck(
            "area_0_40", ">= 0.090 m rad", area_0_40, "m rad", area_0_40 >= 0.090, area_0_40 - 0.090
        ),
        CriterionCheck(
            "area_30_40", ">= 0.030 m rad", area_30_40, "m rad", area_30_40 >= 0.030,
            area_30_40 - 0.030,
        ),
        CriterionCheck("gz_at_30", ">= 0.20 m", gz_30, "m", gz_30 >= 0.20, gz_30 - 0.20),
        CriterionCheck("max_gz", ">= 0.20 m", maximum, "m", maximum >= 0.20, maximum - 0.20),
        CriterionCheck(
            "angle_at_max_gz",
            ">= 25 deg",
            angle_at_maximum,
            "deg",
            angle_at_maximum >= 25.0,
            angle_at_maximum - 25.0,
        ),
        CriterionCheck("gm_upright", ">= 0.15 m", gm, "m", gm >= 0.15, gm - 0.15),
    )
    return IMOIntactStability(
        checks=checks,
        passed=all(check.satisfied for check in checks),
        gm_upright=gm,
    )



# ``StabilityError`` is re-exported so callers that catch the initial-stability error also
# catch the heeled-waterline one without importing two modules.
StabilityError = StabilityError

__all__ = [
    "DEFAULT_ANGLES_DEG",
    "DEFAULT_HEEL_ANGLES_DEG",
    "SCHEMA_GZ_CURVE",
    "SCHEMA_HEELED",
    "SCHEMA_IMO",
    "SCHEMA_RIGHTING_ARM",
    "SCHEMA_SUBMERGED",
    "CriterionCheck",
    "FreeSurfaceCorrection",
    "GZCurve",
    "GZError",
    "GzError",
    "HeelPoint",
    "HeeledHydrostatics",
    "IMOIntactStability",
    "RightingArmCurve",
    "StabilityError",
    "SubmergedSolid",
    "_Integral",
    "box_righting_arms_analytic",
    "box_wall_sided_limit_deg",
    "deck_edge_immersion_deg",
    "free_surface_correction",
    "gz_curve",
    "heel_mesh",
    "heeled_hydrostatics",
    "imo_intact_stability",
    "rectangular_tank_inertia",
    "righting_arm_curve",
    "solve_waterline",
    "submerged_solid",
    "wall_sided_righting_arm",
    "with_freeboard",
]
