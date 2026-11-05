"""Large-angle righting arms: the heeled *displacement* integration behind ``GZ``.

:mod:`tensorcad.stability` answers the upright question -- ``GM = KB + BM - KG``, the
slope of the righting-arm curve at zero heel. This module answers the curve itself:
hold the displaced volume fixed, heel the hull, and measure how far the line of
buoyancy has moved from the centre of gravity.

**What "GZ" is, precisely.** In the body frame (``x`` forward, ``y`` to starboard,
``z`` up) heeling to starboard by ``phi`` puts world-up at

    w = (0, -sin(phi), cos(phi)),

so the submerged region is ``{p : w . p <= c}`` for some plane offset ``c``, and the
horizontal direction (in the body frame, perpendicular to ``w``) is
``e_h = (0, w_z, -w_y)``. Then

    GZ(phi) = (B - G) . e_h = (y_B - y_G) cos(phi) + (z_B - z_G) sin(phi),

and the righting moment is ``displacement * GZ``. The sign convention is chosen so
that heeling to starboard is ``phi > 0`` and a *positive* ``GZ`` is a *righting* arm;
at small angles this reduces to ``GZ ~ GM sin(phi)``, which is exactly the derivate
check the tests make against :func:`tensorcad.stability.initial_stability`.

**The offset is solved, not assumed.** A GZ curve is drawn at constant displacement:
for every heel angle the waterline has to rise or fall until the submerged volume
equals the upright displacement. That is a one-dimensional root find on ``c`` with a
monotone function, and it is solved here by a Newton iteration whose derivative is
exact -- ``dV/dc`` is the waterplane area measured at that offset, so no finite
difference enters. The residual volume error is reported on every point instead of
being assumed small.

**The integral is exact for the volume the mesh actually encloses.** The surface is
cut with the tilted waterplane by clipping each triangle against the plane and
keeping the part below; the submerged solid is closed by the planar cut, and a
tetrahedron decomposition about a point *on that plane* makes the cut contribute
exactly zero (a tetrahedron with all four corners coplanar has no volume). So the
volume and the centroid of buoyancy are obtained from the submerged part of the
hull surface alone, without ever triangulating the waterplane -- no fan triangulation
of a possibly non-convex section, and no hole bookkeeping. The cut area is reported
as the flux of that same surface, which is exact for the same reason.

**What this is not.** Solid-hull, still-water, free-trim heel about a longitudinal
axis: no free-surface correction, no damage or floodable compartments, no trim, no
wind or lifting surfaces, no stability criteria. The numbers this returns are raw
arms; whether the vessel "passes" is a rule question that needs KG, the tankage and
the statutory criterion, none of which a geometry library knows. A mesh capped at the
design waterline also has no freeboard: floating it at ``waterline_z < 0`` leaves
deck to immerse, and past that angle the arms belong to a box, not to a ship.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

from . import mesh as meshlib
from .geometry import Mesh, Vec3

VOLUME_TOLERANCE = 1e-12
GEOMETRY_TOLERANCE = 1e-12


class RightingArmError(ValueError):
    """Raised when a righting arm cannot be defined for the arguments given."""


# --------------------------------------------------------------------- vector bits


def _unit(vector: Sequence[float]) -> Vec3:
    norm = math.sqrt(vector[0] ** 2 + vector[1] ** 2 + vector[2] ** 2)
    if not math.isfinite(norm) or norm <= 0.0:
        raise RightingArmError(f"cannot normalise {tuple(vector)!r}")
    return (vector[0] / norm, vector[1] / norm, vector[2] / norm)


def _determinant(u: Sequence[float], v: Sequence[float], w: Sequence[float]) -> float:
    return (
        u[0] * (v[1] * w[2] - v[2] * w[1])
        - u[1] * (v[0] * w[2] - v[2] * w[0])
        + u[2] * (v[0] * w[1] - v[1] * w[0])
    )


def heel_normal(heel_rad: float) -> Vec3:
    """World-up expressed in the body frame for a starboard heel of ``heel_rad``.

    ``heel_rad > 0`` dips the starboard side, so the submerged region is
    ``{p : heel_normal . p <= offset}`` and :func:`righting_arm` returns a righting
    arm that is positive when the hull is stable at that angle.
    """
    if not math.isfinite(heel_rad):
        raise RightingArmError(f"heel angle must be finite, got {heel_rad!r}")
    return (0.0, -math.sin(heel_rad), math.cos(heel_rad))


def _horizontal_direction(normal: Vec3) -> Vec3:
    """In-plane horizontal direction with the GZ sign convention (see the module doc)."""
    return (0.0, normal[2], -normal[1])


# ------------------------------------------------------------------ plane clipping


@dataclass(frozen=True)
class SubmergedProperties:
    """Volume and buoyancy centroid of the part of a closed surface below a plane."""

    volume: float
    centre_of_buoyancy: Vec3
    cut_area: float
    plane_normal: Vec3
    plane_offset: float
    faces_below: int
    faces_cut: int

    def to_dict(self) -> dict[str, object]:
        return {
            "volume_m3": _round(self.volume),
            "centre_of_buoyancy_m": [_round(value) for value in self.centre_of_buoyancy],
            "cut_area_m2": _round(self.cut_area),
            "plane_normal": [_round(value, 12) for value in self.plane_normal],
            "plane_offset_m": _round(self.plane_offset),
            "faces_below": self.faces_below,
            "faces_cut": self.faces_cut,
        }


def _round(value: float, digits: int = 12) -> float:
    if not math.isfinite(value):
        return value
    return round(value, digits)


def _clip_polygon(
    points: Sequence[Vec3], heights: Sequence[float], tolerance: float
) -> list[Vec3]:
    """Sutherland-Hodgman clip of an ordered polygon against ``height <= 0``.

    Winding is preserved, so a fan triangulation of the result keeps the original
    triangle's outward orientation.
    """
    clipped: list[Vec3] = []
    count = len(points)
    for index in range(count):
        point = points[index]
        height = heights[index]
        next_point = points[(index + 1) % count]
        next_height = heights[(index + 1) % count]
        inside = height <= tolerance
        next_inside = next_height <= tolerance
        if inside:
            clipped.append(point)
        if inside != next_inside:
            fraction = height / (height - next_height)
            clipped.append(
                (
                    point[0] + fraction * (next_point[0] - point[0]),
                    point[1] + fraction * (next_point[1] - point[1]),
                    point[2] + fraction * (next_point[2] - point[2]),
                )
            )
    return clipped


def submerged_properties(
    mesh: Mesh,
    normal: Sequence[float] = (0.0, 0.0, 1.0),
    plane_offset: float = 0.0,
    *,
    tolerance: float = GEOMETRY_TOLERANCE,
    closed: bool | None = None,
) -> SubmergedProperties:
    """Integrate volume, buoyancy centroid and cut area below ``normal . p = offset``.

    The cut face is not triangulated: the tetrahedra are taken about a point on the
    waterplane, so the coplanar cut contributes nothing and the submerged part of the
    hull surface alone is a sufficient closed boundary for the volume integral.
    """
    unit = _unit(normal)
    if not math.isfinite(plane_offset):
        raise RightingArmError(f"plane_offset must be finite, got {plane_offset!r}")
    if closed is None:
        closed = meshlib.is_watertight(mesh)
    if not closed:
        raise RightingArmError("a righting arm needs a closed surface; the mesh is not watertight")

    vertices = mesh.vertices
    heights = [
        unit[0] * vertex[0] + unit[1] * vertex[1] + unit[2] * vertex[2] - plane_offset
        for vertex in vertices
    ]
    if not heights:
        raise RightingArmError("cannot integrate an empty mesh")
    scale = max(abs(value) for value in heights)
    point_tolerance = tolerance * scale if scale > 0.0 else tolerance
    # A point in the plane: with |n| = 1, n . (offset * n) = offset.
    origin = (unit[0] * plane_offset, unit[1] * plane_offset, unit[2] * plane_offset)

    volumes: list[float] = []
    moment_x: list[float] = []
    moment_y: list[float] = []
    moment_z: list[float] = []
    flux_x: list[float] = []
    flux_y: list[float] = []
    flux_z: list[float] = []
    faces_below = 0
    faces_cut = 0

    for face in mesh.faces:
        a, b, c = face
        ha, hb, hc = heights[a], heights[b], heights[c]
        if ha > point_tolerance and hb > point_tolerance and hc > point_tolerance:
            continue
        if ha <= point_tolerance and hb <= point_tolerance and hc <= point_tolerance:
            faces_below += 1
            polygon: Sequence[Vec3] = (vertices[a], vertices[b], vertices[c])
        else:
            faces_cut += 1
            polygon = _clip_polygon(
                (vertices[a], vertices[b], vertices[c]), (ha, hb, hc), point_tolerance
            )
            if len(polygon) < 3:
                continue
        head = polygon[0]
        hx, hy, hz = head[0] - origin[0], head[1] - origin[1], head[2] - origin[2]
        for index in range(1, len(polygon) - 1):
            p, q = polygon[index], polygon[index + 1]
            px, py, pz = p[0] - origin[0], p[1] - origin[1], p[2] - origin[2]
            qx, qy, qz = q[0] - origin[0], q[1] - origin[1], q[2] - origin[2]
            volume = _determinant((hx, hy, hz), (px, py, pz), (qx, qy, qz)) / 6.0
            volumes.append(volume)
            # Tetrahedron (origin, head, p, q) centroid = (origin + head + p + q) / 4.
            moment_x.append(volume * 0.25 * (origin[0] + head[0] + p[0] + q[0]))
            moment_y.append(volume * 0.25 * (origin[1] + head[1] + p[1] + q[1]))
            moment_z.append(volume * 0.25 * (origin[2] + head[2] + p[2] + q[2]))
            ux, uy, uz = p[0] - head[0], p[1] - head[1], p[2] - head[2]
            vx, vy, vz = q[0] - head[0], q[1] - head[1], q[2] - head[2]
            flux_x.append(0.5 * (uy * vz - uz * vy))
            flux_y.append(0.5 * (uz * vx - ux * vz))
            flux_z.append(0.5 * (ux * vy - uy * vx))

    volume = math.fsum(volumes)
    if volume <= 0.0 or not math.isfinite(volume):
        raise RightingArmError(
            f"the plane {tuple(unit)!r} . p = {plane_offset!r} bounds no positive volume"
        )
    centroid: Vec3 = (
        math.fsum(moment_x) / volume,
        math.fsum(moment_y) / volume,
        math.fsum(moment_z) / volume,
    )
    # Flux of the closed submerged surface is zero, so the submerged hull surface
    # carries the cut face's flux with the opposite sign: its magnitude is the cut area.
    cut_area = math.sqrt(
        math.fsum(flux_x) ** 2 + math.fsum(flux_y) ** 2 + math.fsum(flux_z) ** 2
    )
    return SubmergedProperties(
        volume=volume,
        centre_of_buoyancy=centroid,
        cut_area=cut_area,
        plane_normal=unit,
        plane_offset=plane_offset,
        faces_below=faces_below,
        faces_cut=faces_cut,
    )


# ---------------------------------------------------------------- constant volume


@dataclass(frozen=True)
class PlaneSolution:
    """Offset that holds the displaced volume fixed at a given heel normal."""

    offset: float
    volume: float
    target_volume: float
    cut_area: float
    iterations: int
    residual: float

    @property
    def relative_residual(self) -> float:
        return abs(self.volume - self.target_volume) / abs(self.target_volume)


def solve_plane_offset(
    mesh: Mesh,
    normal: Sequence[float],
    target_volume: float,
    *,
    tolerance: float = VOLUME_TOLERANCE,
    max_iterations: int = 60,
    initial_offset: float | None = None,
    closed: bool | None = None,
) -> PlaneSolution:
    """Find the waterline offset that displaces ``target_volume`` at this heel.

    Newton's method on ``V(c) - target`` with the exact derivative ``dV/dc`` equal to
    the cut area, safeguarded by bisection on the bracket ``min(n . p) <= c <= max(n . p)``
    so the iteration cannot leave the flooded/emerged limits. Monotone in ``c`` by
    construction, so the solution is unique.
    """
    if not math.isfinite(target_volume) or target_volume <= 0.0:
        raise RightingArmError(f"target_volume must be finite and positive, got {target_volume!r}")
    if closed is None:
        closed = meshlib.is_watertight(mesh)
    if not closed:
        raise RightingArmError("a righting arm needs a closed surface; the mesh is not watertight")
    unit = _unit(normal)
    heights = [
        unit[0] * vertex[0] + unit[1] * vertex[1] + unit[2] * vertex[2]
        for vertex in mesh.vertices
    ]
    if not heights:
        raise RightingArmError("cannot integrate an empty mesh")
    low = min(heights)
    high = max(heights)
    total = abs(meshlib.total_signed_volume(mesh))
    if total <= 0.0:
        raise RightingArmError("the mesh encloses no volume")
    if target_volume >= total:
        raise RightingArmError(
            f"target volume {target_volume!r} m^3 is not less than the {total!r} m^3 the mesh encloses"
        )
    if initial_offset is None:
        guess = low + (high - low) * (target_volume / total)
    else:
        guess = min(max(float(initial_offset), low), high)
    offset = guess
    best = submerged_properties(mesh, unit, offset, closed=True)
    iterations = 1
    while iterations < max_iterations:
        error = best.volume - target_volume
        if abs(error) <= tolerance * target_volume:
            break
        if error < 0.0:
            low = offset
        else:
            high = offset
        derivative = best.cut_area
        candidate = offset - error / derivative if derivative > 0.0 else float("nan")
        if not (low < candidate < high) or abs(candidate - offset) <= 1e-16 * max(1.0, abs(offset)):
            candidate = 0.5 * (low + high)
        if candidate == offset:
            break
        offset = candidate
        best = submerged_properties(mesh, unit, offset, closed=True)
        iterations += 1
    return PlaneSolution(
        offset=offset,
        volume=best.volume,
        target_volume=target_volume,
        cut_area=best.cut_area,
        iterations=iterations,
        residual=best.volume - target_volume,
    )


def upright_displacement(mesh: Mesh, waterline_z: float = 0.0, *, closed: bool | None = None) -> float:
    """Displaced volume under the level waterline ``z = waterline_z``, by the same cut."""
    return submerged_properties(mesh, (0.0, 0.0, 1.0), waterline_z, closed=closed).volume


# ------------------------------------------------------------------------ GZ curve


@dataclass(frozen=True)
class HeeledFloat:
    """One point of a righting-arm curve: the heeled floating body at fixed volume."""

    heel_deg: float
    heel_rad: float
    displacement_volume: float
    plane_offset: float
    centre_of_buoyancy: Vec3
    waterplane_area: float
    kg: float
    keel_z: float
    righting_arm: float
    volume_residual: float
    iterations: int

    @property
    def relative_volume_residual(self) -> float:
        return abs(self.volume_residual) / self.displacement_volume

    def to_dict(self) -> dict[str, object]:
        return {
            "heel_deg": _round(self.heel_deg),
            "displacement_volume_m3": _round(self.displacement_volume),
            "plane_offset_m": _round(self.plane_offset),
            "centre_of_buoyancy_m": [_round(value) for value in self.centre_of_buoyancy],
            "waterplane_area_m2": _round(self.waterplane_area),
            "kg_m": _round(self.kg),
            "gz_m": _round(self.righting_arm),
            "volume_residual_m3": _round(self.volume_residual, 15),
            "iterations": self.iterations,
        }


@dataclass(frozen=True)
class RightingArmCurve:
    """``GZ`` against heel at constant displacement."""

    points: tuple[HeeledFloat, ...]
    kg: float
    displacement_volume: float
    keel_z: float
    y_g: float = 0.0
    basis: str = "heeled displacement integration, solid hull, free trim"

    @property
    def angles_deg(self) -> tuple[float, ...]:
        return tuple(point.heel_deg for point in self.points)

    @property
    def arms(self) -> tuple[float, ...]:
        return tuple(point.righting_arm for point in self.points)

    @property
    def max_arm(self) -> float:
        return max(self.arms)

    @property
    def angle_of_max_arm_deg(self) -> float:
        best = max(self.points, key=lambda point: point.righting_arm)
        return best.heel_deg

    @property
    def zero_crossing_deg(self) -> float | None:
        """First heel angle where the arm changes sign, interpolated linearly."""
        previous: HeeledFloat | None = None
        for point in self.points:
            if previous is not None and previous.righting_arm > 0.0 >= point.righting_arm:
                span = previous.righting_arm - point.righting_arm
                if span <= 0.0:
                    return point.heel_deg
                fraction = previous.righting_arm / span
                return previous.heel_deg + fraction * (point.heel_deg - previous.heel_deg)
            previous = point
        return None

    def arm_at(self, heel_deg: float) -> float:
        return _interpolate(self.angles_deg, self.arms, heel_deg)

    def area(self, limit_deg: float = 30.0) -> float:
        """Area under the curve from ``0`` to ``limit_deg`` in metre-radians.

        Trapezoidal on the computed points (the arm is a smooth function of heel, so
        the quadrature error is set by the sample spacing, which is reported rather
        than hidden); callers writing a rule check can supply their own spacing.
        """
        if limit_deg <= 0.0:
            raise RightingArmError(f"limit_deg must be positive, got {limit_deg!r}")
        angles = list(self.angles_deg)
        arms = list(self.arms)
        if not any(abs(angle) <= 1e-12 for angle in angles):
            raise RightingArmError("the curve must contain the upright point to integrate from 0")
        if limit_deg > angles[-1]:
            raise RightingArmError(
                f"limit {limit_deg!r} deg is outside the sampled curve (last point {angles[-1]!r} deg)"
            )
        total = 0.0
        for index in range(len(angles) - 1):
            start, end = angles[index], angles[index + 1]
            if start >= limit_deg:
                break
            stop = min(end, limit_deg)
            start_arm = arms[index]
            end_arm = (
                arms[index + 1] if end <= limit_deg else _interpolate(self.angles_deg, self.arms, stop)
            )
            total += 0.5 * (start_arm + end_arm) * math.radians(stop - start)
        return total

    def to_dict(self) -> dict[str, object]:
        return {
            "kg_m": _round(self.kg),
            "y_g_m": _round(self.y_g),
            "keel_z_m": _round(self.keel_z),
            "displacement_volume_m3": _round(self.displacement_volume),
            "max_arm_m": _round(self.max_arm),
            "angle_of_max_arm_deg": _round(self.angle_of_max_arm_deg),
            "zero_crossing_deg": (
                None if self.zero_crossing_deg is None else _round(self.zero_crossing_deg)
            ),
            "area_to_30_deg_m_rad": _round(self.area(30.0)) if self.angles_deg[-1] >= 30.0 else None,
            "area_to_40_deg_m_rad": _round(self.area(40.0)) if self.angles_deg[-1] >= 40.0 else None,
            "basis": self.basis,
            "points": [point.to_dict() for point in self.points],
        }


def _interpolate(angles: Sequence[float], arms: Sequence[float], heel_deg: float) -> float:
    if heel_deg <= angles[0]:
        return arms[0]
    if heel_deg >= angles[-1]:
        return arms[-1]
    for index in range(len(angles) - 1):
        start, end = angles[index], angles[index + 1]
        if start <= heel_deg <= end:
            span = end - start
            if span <= 0.0:
                return arms[index]
            fraction = (heel_deg - start) / span
            return arms[index] + fraction * (arms[index + 1] - arms[index])
    return arms[-1]


def righting_arm(
    mesh: Mesh,
    heel_deg: float = 30.0,
    *,
    kg: float,
    waterline_z: float = 0.0,
    displacement_volume: float | None = None,
    y_g: float = 0.0,
    tolerance: float = VOLUME_TOLERANCE,
    closed: bool | None = None,
) -> HeeledFloat:
    """Righting arm at one heel angle, at the upright displacement.

    ``kg`` is measured above the keel, matching
    :func:`tensorcad.stability.initial_stability`; ``kg = 0`` puts G on the keel and
    makes ``GZ`` the horizontal offset of the buoyancy centroid. ``waterline_z``
    fixes the reference displacement when ``displacement_volume`` is not given.
    """
    if not math.isfinite(kg):
        raise RightingArmError(f"kg must be finite, got {kg!r}")
    if not math.isfinite(heel_deg):
        raise RightingArmError(f"heel_deg must be finite, got {heel_deg!r}")
    if closed is None:
        closed = meshlib.is_watertight(mesh)
    if not closed:
        raise RightingArmError("a righting arm needs a closed surface; the mesh is not watertight")
    if displacement_volume is None:
        displacement_volume = upright_displacement(mesh, waterline_z, closed=True)
    if displacement_volume <= 0.0:
        raise RightingArmError(
            f"displacement volume {displacement_volume!r} m^3 is not positive; "
            "a righting arm needs a floating body"
        )
    phi = math.radians(heel_deg)
    normal = heel_normal(phi)
    solution = solve_plane_offset(
        mesh, normal, displacement_volume, tolerance=tolerance, closed=True
    )
    properties = submerged_properties(mesh, normal, solution.offset, closed=True)
    keel_z = min(vertex[2] for vertex in mesh.vertices)
    buoyancy = properties.centre_of_buoyancy
    horizontal = _horizontal_direction(normal)
    arm = (buoyancy[1] - y_g) * horizontal[1] + (buoyancy[2] - (keel_z + kg)) * horizontal[2]
    return HeeledFloat(
        heel_deg=heel_deg,
        heel_rad=phi,
        displacement_volume=displacement_volume,
        plane_offset=solution.offset,
        centre_of_buoyancy=buoyancy,
        waterplane_area=properties.cut_area,
        kg=kg,
        keel_z=keel_z,
        righting_arm=arm,
        volume_residual=solution.residual,
        iterations=solution.iterations,
    )


def righting_arm_curve(
    mesh: Mesh,
    angles_deg: Iterable[float] = (0.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0),
    *,
    kg: float,
    waterline_z: float = 0.0,
    displacement_volume: float | None = None,
    y_g: float = 0.0,
    tolerance: float = VOLUME_TOLERANCE,
    closed: bool | None = None,
) -> RightingArmCurve:
    """Sample the righting-arm curve, all points at the same displacement."""
    if closed is None:
        closed = meshlib.is_watertight(mesh)
    if not closed:
        raise RightingArmError("a righting arm needs a closed surface; the mesh is not watertight")
    angles = tuple(sorted(float(angle) for angle in angles_deg))
    if not angles:
        raise RightingArmError("at least one heel angle is required")
    if displacement_volume is None:
        displacement_volume = upright_displacement(mesh, waterline_z, closed=True)
    points = tuple(
        righting_arm(
            mesh,
            angle,
            kg=kg,
            waterline_z=waterline_z,
            displacement_volume=displacement_volume,
            y_g=y_g,
            tolerance=tolerance,
            closed=True,
        )
        for angle in angles
    )
    return RightingArmCurve(
        points=points,
        kg=kg,
        displacement_volume=displacement_volume,
        keel_z=points[0].keel_z,
        y_g=y_g,
    )


# ------------------------------------------------------------------- references


def wall_sided_righting_arm(heel_rad: float, *, gm: float, bm: float) -> float:
    """Closed form for a vessel with vertical sides: ``GZ = sin(phi) (GM + BM tan^2(phi) / 2)``.

    Exact for a rectangular section until a deck edge immerses, and the classical
    first check on any large-angle code: it fixes the whole curve from two upright
    quantities, so an implementation that agrees with it over the wall-sided range
    cannot have the sign, the offset solve or the moment arm wrong.
    """
    tangent = math.tan(heel_rad)
    return math.sin(heel_rad) * (gm + 0.5 * bm * tangent * tangent)


def bilge_emergence_deg(*, beam: float, draft: float) -> float:
    """Heel angle at which the high-side bilge of a box leaves the water.

    Above it one side of the section is dry, the section stops being wall-sided, and
    :func:`wall_sided_righting_arm` stops being exact -- which is the range boundary
    the tests use rather than an assumption buried in a comment.
    """
    if beam <= 0.0 or draft <= 0.0:
        raise RightingArmError("beam and draft must be positive")
    return math.degrees(math.atan(2.0 * draft / beam))


def deck_edge_immersion_deg(*, beam: float, depth: float) -> float:
    """Heel angle at which the low-side deck edge of a box reaches the waterline.

    ``depth`` is the full height above the keel; the wall-sided form above is exact
    only below this angle, which for the reference box (``B/T`` and freeboard as
    used in the tests) lands at ``45 deg``.
    """
    if beam <= 0.0 or depth <= 0.0:
        raise RightingArmError("beam and depth must be positive")
    return math.degrees(math.atan(2.0 * depth / beam))


def box_mesh(length: float, beam: float, height: float, *, z0: float = 0.0) -> Mesh:
    """A closed rectangular box: the one hull with a closed form at every heel angle."""
    from .geometry import MeshBuilder

    for name, value in (("length", length), ("beam", beam), ("height", height)):
        if not math.isfinite(value) or value <= 0.0:
            raise RightingArmError(f"{name} must be finite and positive, got {value!r}")
    half_length = 0.5 * length
    half_beam = 0.5 * beam
    top = z0 + height
    builder = MeshBuilder(units="m", coordinate_frame="x forward, y starboard, z up")
    corners = {}
    for xi, x in enumerate((-half_length, half_length)):
        for yi, y in enumerate((-half_beam, half_beam)):
            for zi, z in enumerate((z0, top)):
                corners[(xi, yi, zi)] = builder.add_vertex((x, y, z))
    quads = (
        ((0, 0, 0), (0, 1, 0), (0, 1, 1), (0, 0, 1)),  # aft face
        ((1, 1, 0), (1, 0, 0), (1, 0, 1), (1, 1, 1)),  # forward face
        ((1, 0, 0), (0, 0, 0), (0, 0, 1), (1, 0, 1)),  # port side
        ((0, 1, 0), (1, 1, 0), (1, 1, 1), (0, 1, 1)),  # starboard side
        ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)),  # bottom
        ((0, 1, 1), (1, 1, 1), (1, 0, 1), (0, 0, 1)),  # top
    )
    for quad in quads:
        builder.add_quad(*(corners[key] for key in quad))
    return meshlib.orient_outward(builder.build()).mesh


def _polygon_area_centroid(polygon: Sequence[tuple[float, float]]) -> tuple[float, float, float]:
    """Signed area and centroid of a 2D polygon (Green's theorem)."""
    twice_area = 0.0
    cx = 0.0
    cy = 0.0
    count = len(polygon)
    for index in range(count):
        x0, y0 = polygon[index]
        x1, y1 = polygon[(index + 1) % count]
        cross = x0 * y1 - x1 * y0
        twice_area += cross
        cx += (x0 + x1) * cross
        cy += (y0 + y1) * cross
    area = 0.5 * twice_area
    if abs(area) <= 0.0:
        return 0.0, float("nan"), float("nan")
    return area, cx / (3.0 * twice_area), cy / (3.0 * twice_area)


def analytic_box_righting_arm(
    length: float,
    beam: float,
    height: float,
    draft: float,
    heel_deg: float,
    *,
    kg: float,
    y_g: float = 0.0,
) -> float:
    """Independent 2D reference for a box: clip the section, integrate, then arm.

    The box is prismatic, so its heeled displacement integral is a *section* problem:
    clip the rectangle by the tilted waterline, get the section area and centroid,
    and solve for the waterline that recovers the upright area ``beam * draft``. This
    shares no code with the 3D tetrahedron integration in :func:`righting_arm`, which
    is what makes it a check rather than a restatement. Valid at *every* heel angle,
    including past deck-edge immersion.
    """
    if draft <= 0.0 or draft >= height:
        raise RightingArmError(f"draft must lie between 0 and the box height, got {draft!r}")
    phi = math.radians(heel_deg)
    half_beam = 0.5 * beam
    # Section frame: y across, z up, upright waterline at z = 0, keel at z = -draft.
    rectangle = (
        (-half_beam, -draft),
        (half_beam, -draft),
        (half_beam, height - draft),
        (-half_beam, height - draft),
    )
    target = beam * draft

    def clipped(offset: float) -> list[tuple[float, float]]:
        # Waterline: -y sin(phi) + z cos(phi) = offset; keep the part below it.
        heights = [-y * math.sin(phi) + z * math.cos(phi) - offset for y, z in rectangle]
        out: list[tuple[float, float]] = []
        count = len(rectangle)
        for index in range(count):
            point = rectangle[index]
            height = heights[index]
            next_point = rectangle[(index + 1) % count]
            next_height = heights[(index + 1) % count]
            inside = height <= 0.0
            if inside:
                out.append(point)
            if inside != (next_height <= 0.0):
                fraction = height / (height - next_height)
                out.append(
                    (
                        point[0] + fraction * (next_point[0] - point[0]),
                        point[1] + fraction * (next_point[1] - point[1]),
                    )
                )
        return out

    def area(offset: float) -> float:
        polygon = clipped(offset)
        if len(polygon) < 3:
            return 0.0
        return abs(_polygon_area_centroid(polygon)[0])

    low, high = -height - beam, height + beam
    for _ in range(200):
        middle = 0.5 * (low + high)
        if area(middle) < target:
            low = middle
        else:
            high = middle
        if high - low <= 1e-15 * max(1.0, height):
            break
    offset = 0.5 * (low + high)
    polygon = clipped(offset)
    section_area, y_b, z_b = _polygon_area_centroid(polygon)
    if section_area <= 0.0 or not math.isfinite(y_b):
        raise RightingArmError("the box section is empty at this heel angle")
    # The centroid of the submerged prism is the section centroid, independent of x.
    # kg enters above the keel; this section frame has z = 0 at the upright waterline.
    horizontal_y, horizontal_z = math.cos(phi), math.sin(phi)
    return (y_b - y_g) * horizontal_y + (z_b + draft - kg) * horizontal_z
