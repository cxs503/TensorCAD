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
keel line -- and there is no tank model inside the integrator.
:func:`free_surface_correction` supplies the classical ``sum(i_x)/V`` loss of ``GM``, applied
the way that loss belongs on a curve -- ``-loss * sin(phi)``, i.e. a raise of ``KG``, so an
upright symmetric ship still has ``GZ(0) = 0`` -- and
:mod:`tensorcad.free_surface` supplies the *exact* heel-dependent shift of a prismatic slack
tank for callers who want the finite-angle answer instead of the first-order one. :func:`with_freeboard` exists because a hull generated
to its design waterline simply *cannot* be heeled -- on the dipping side the waterline
leaves the modelled surface immediately -- so a closed freeboard band has to be added
first. Damage stability, deck-edge immersion criteria, grain heeling, statutory maxima and
class rules are all outside this module, and nothing here should be read as a check
against a statutory criterion.

Only the Python standard library is used.
"""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass

from . import mesh as meshlib
from .free_surface import SlackTank, free_surface_gz_shift_curve, free_surface_loss_m
from .geometry import Mesh

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
    free_surface_model: str = "kg_raise"

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
            "free_surface_model": self.free_surface_model,
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
    """The ``sum(i_x) / V`` loss of ``GM`` a set of slack tanks costs.

    ``loss_m`` is the metres of metacentric height, i.e. the raise of ``KG``: a righting-arm
    curve loses ``loss_m * sin(phi)`` at heel ``phi``, not ``loss_m`` at every ordinate.
    """

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
    free_surface_model: str = "kg_raise",
    slack_tanks: Sequence[SlackTank] = (),
    tolerance: float = 1e-13,
) -> RightingArmCurve:
    """Sample ``GZ`` over ``angles_deg`` for a centre of gravity ``kg`` above the keel.

    ``free_surface_i_x`` is the sum of the free-surface second moments ``i_x`` of any slack
    tanks, i.e. their loss of ``GM`` is ``sum(i_x)/V``. How that loss enters the curve is the
    caller's choice, and it is not cosmetic:

    * ``free_surface_model="kg_raise"`` (default) -- ``GZ - sum(i_x)/V * sin(phi)``. This is
      what raising ``KG`` by the loss does to the curve, it is the form the classical
      correction is derived in, and it keeps ``GZ(0) = 0`` for an upright symmetric ship.
    * ``free_surface_model="exact"`` -- the finite-angle shift of the actual tanks, from
      :func:`tensorcad.free_surface.free_surface_gz_shift_curve`, which needs the tank
      geometry in ``slack_tanks`` and no ``free_surface_i_x``. Exact for a prismatic
      rectangular tank; the first-order model above is its small-angle limit.
    * ``free_surface_model="flat"`` -- the historical ``GZ - sum(i_x)/V``, retained only so a
      curve published by an earlier release can be reproduced. It is wrong at zero heel: an
      upright symmetric ship gets ``GZ(0) = -sum(i_x)/V`` instead of zero, and every
      intercepted area moves by ``sum(i_x)/V`` times the heel interval rather than by
      ``sum(i_x)/V`` times ``(1 - cos)``.
    """
    if not math.isfinite(kg):
        raise GzError(f"kg must be finite: {kg!r}")
    if not math.isfinite(free_surface_i_x) or free_surface_i_x < 0.0:
        raise GzError(f"free_surface_i_x must be finite and non-negative: {free_surface_i_x!r}")
    if free_surface_model not in ("kg_raise", "flat", "exact"):
        raise GzError(
            f"free_surface_model must be 'kg_raise', 'flat' or 'exact': {free_surface_model!r}"
        )
    tanks = tuple(slack_tanks)
    if free_surface_model == "exact":
        if not tanks:
            raise GzError("free_surface_model='exact' needs at least one slack tank")
        if free_surface_i_x:
            raise GzError(
                "pass the tanks, not free_surface_i_x, when free_surface_model='exact': the "
                "exact shift already carries their inertias"
            )
    elif tanks:
        raise GzError(
            "slack tanks are only used by free_surface_model='exact'; for the first-order "
            "model pass their summed free_surface_i_x"
        )
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
    if free_surface_model == "exact":
        shift = free_surface_gz_shift_curve(tanks, angles, target)
        correction = free_surface_loss_m(tanks, target)
        gz = tuple(
            point.kn - kg * math.sin(math.radians(point.heel_angle_deg)) + offset
            for point, offset in zip(points, shift, strict=True)
        )
    else:
        correction = free_surface_i_x / target
        if free_surface_model == "flat":
            gz = tuple(
                point.kn - kg * math.sin(math.radians(point.heel_angle_deg)) - correction
                for point in points
            )
        else:
            gz = tuple(
                point.kn
                - (kg + correction) * math.sin(math.radians(point.heel_angle_deg))
                for point in points
            )
    # ``kg`` is measured above the keel and the levers ``kn`` are taken about the keel
    # line (the rotation axis), so the arm is ``kn - kg sin(phi)`` before the tanks.
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
        free_surface_model=free_surface_model,
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
    """``sum(i_x) / V``: the loss of ``GM`` a slack tank costs, as a raise of ``KG``.

    The number is a *metacentric* quantity. :func:`righting_arm_curve` applies it as
    ``-loss * sin(phi)``, which is what a raise of ``KG`` by ``loss`` does to a curve; the
    historical flat subtraction ``GZ - loss`` is still available as
    ``free_surface_model="flat"`` and is wrong at zero heel (it gives an upright symmetric
    ship a non-zero arm), which is why it is no longer the default.
    """
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
