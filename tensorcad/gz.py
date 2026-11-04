"""Large-angle stability: the righting-arm curve ``GZ(phi)``.

:mod:`tensorcad.stability` stops at the slope of the righting-arm curve at zero heel.
It says so in its own docstring: ``GM`` is the slope, and "stability at large angles is
a different integration (a heeled displacement integration, not a waterplane moment)".
This module is that integration.

**What is computed.** For a closed, outward-oriented surface mesh floating at a fixed
displacement, heel the waterplane to ``phi`` about the longitudinal axis and find the
plane that still displaces the upright volume -- the ship is *not* re-displaced as it
heels (no flooding, no sinkage), so the plane offset is an unknown that has to be
solved for, and it is. The centre of buoyancy of the heeled region is then a volume
centroid in the ship's own frame, and the righting arm is its projection on the
heeled waterplane:

    up        n(phi) = (0, -sin phi, cos phi)      (world up in the ship frame)
    in-plane  e(phi) = (0,  cos phi, sin phi)
    GZ        = (B(phi) - G) . e = (y_B - y_G) cos phi + (z_B - z_G) sin phi

with ``G`` at ``(y_G, z_keel + KG)``: ``KG`` is measured above the keel (the mesh's
lowest point), the same datum :func:`tensorcad.stability.initial_stability` uses, and
``y_G`` defaults to the centreline. Positive ``GZ`` rights the ship.

**How the heeled displacement is integrated.** The immersed region is a half-space
intersection of a polytope, so it has an *exact* volume and centroid: every triangle is
clipped against ``n.p <= d`` (Sutherland-Hodgman on a triangle, at most four vertices
out), and the clipped surface is summed as signed tetrahedra from a point ``o = d*n``
that lies *on* the cut plane. That reference point is the whole trick -- the cap that
would close the clipped surface lies in a plane through ``o``, so every cap tetrahedron
has zero volume and the cap never has to be built. The result is exact for the
polyhedron the mesh actually is, in the same spirit as the Green's-theorem waterplane
integral in :mod:`tensorcad.stability`: no sampling, no tolerance, and the only error
is the discretisation the mesh already carries.

The offset ``d`` is solved by a safeguarded secant iteration on
``V(d) - V_upright``, which is monotone in ``d``, and the residual is reported in
``relative_volume_error`` rather than assumed away.

**What it is not.** No free-surface correction (that is the caller's tankage), no
flooding or damage stability (a flooded compartment changes the body, not the heel
angle), no trim coupling (the waterplane is rotated about ``x`` only), no wind or
launching curves, and no statutory criterion -- :func:`gz_curve` reports the geometry
of the curve, not whether it passes IMO A.749 / IS Code 2008. A hull whose mesh is
capped at the design waterline has no freeboard to heel into; the parametric family in
:mod:`tensorcad.hull` is exactly that, so large-angle work needs a hull modelled to the
deck edge (the barge in ``examples/gz`` is one).

Only the Python standard library is used.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Sequence

from .geometry import Mesh

SCHEMA_HEELED = "tensorcad.stability.heeled/1"
SCHEMA_ARM = "tensorcad.stability.righting_arm/1"
SCHEMA_CURVE = "tensorcad.stability.gz_curve/1"

#: Relative tolerance on the displacement balance of the heeled plane solve, in
#: displaced-volume units. Secant reaches ``1e-15`` on the analytic bodies in
#: ``examples/gz`` well inside the iteration cap, so this is a floor and not a budget.
VOLUME_TOLERANCE = 1e-13

#: Cap on the secant iterations. Six is typical, eight is common on a coarse mesh
#: where ``V(d)`` has kinks at the tangent planes; 80 means "this is not converging".
MAX_SOLVE_ITERATIONS = 80


class HeelError(ValueError):
    """Raised when a heeled displacement or a righting arm is not defined."""


# ---------------------------------------------------------------------------
# the frame
# ---------------------------------------------------------------------------
def up_direction(heel_angle: float) -> tuple[float, float, float]:
    """World up expressed in the ship frame for a heel of ``heel_angle`` radians.

    Positive heel puts starboard (``+y``) down, so up acquires a *negative* ``y``
    component -- the sign that makes ``GZ`` positive on a stable ship.
    """
    return (0.0, -math.sin(heel_angle), math.cos(heel_angle))


def waterplane_direction(heel_angle: float) -> tuple[float, float, float]:
    """The transverse, in-waterplane direction ``e`` that ``GZ`` is projected on."""
    return (0.0, math.cos(heel_angle), math.sin(heel_angle))


def _check_angle(heel_angle: float) -> None:
    if not math.isfinite(heel_angle):
        raise HeelError("heel angle must be finite")
    if abs(heel_angle) > 0.5 * math.pi:
        raise HeelError(
            f"heel angle {math.degrees(heel_angle):.6g} deg is past 90 deg; the "
            f"righting-arm convention here is transverse heel only"
        )


def _check_loading(mesh: Mesh, kg: float, yg: float) -> None:
    if not math.isfinite(kg):
        raise HeelError("KG must be finite")
    if not math.isfinite(yg):
        raise HeelError("y_G must be finite")
    if not mesh.faces:
        raise HeelError("the mesh has no faces, so there is nothing to heel")


def keel_z(mesh: Mesh) -> float:
    """The ``KG`` datum: the lowest point of the mesh."""
    return mesh.bounding_box().minimum[2]


def _scale(mesh: Mesh) -> float:
    sizes = mesh.bounding_box().sizes
    scale = max(abs(size) for size in sizes)
    return scale if scale > 0.0 else 1.0


# ---------------------------------------------------------------------------
# exact immersed moments of a polytope below a plane
# ---------------------------------------------------------------------------
def _clip_to_halfspace(
    points: Sequence[tuple[float, float, float]], distances: Sequence[float]
) -> list[tuple[float, float, float]]:
    """Keep the part of a convex polygon where ``distance >= 0`` (``n.p <= d``)."""
    out: list[tuple[float, float, float]] = []
    count = len(points)
    for index in range(count):
        following = (index + 1) % count
        here, there = points[index], points[following]
        s_here, s_there = distances[index], distances[following]
        if s_here >= 0.0:
            out.append(here)
        if (s_here > 0.0 > s_there) or (s_here < 0.0 < s_there):
            t = s_here / (s_here - s_there)
            out.append(
                (
                    here[0] + t * (there[0] - here[0]),
                    here[1] + t * (there[1] - here[1]),
                    here[2] + t * (there[2] - here[2]),
                )
            )
    return out


def _immersed_volume_moments(
    mesh: Mesh, normal: tuple[float, float, float], offset: float
) -> tuple[float, tuple[float, float, float], int]:
    """Volume, centroid and number of emerged faces of ``{p in mesh : n.p <= offset}``.

    Signed tetrahedra are taken from ``o = offset * n``, a point on the cut plane, so
    the unclosed clipped surface integrates the closed region on its own: the missing
    cap is coplanar with ``o`` and contributes zero.

    The third return value counts the faces with at least one vertex *above* the
    waterplane. Zero means the plane does not cut the surface at all: either nothing
    is immersed (which raises) or the whole body is below it, and the whole-body
    case is a real trap -- a mesh capped at the design waterline has no freeboard, so
    every heel "conserves" the displaced volume by submerging the entire model and
    reports the constant arm ``KB - KG`` instead of a righting arm.
    """
    nx, ny, nz = normal
    ox, oy, oz = offset * nx, offset * ny, offset * nz
    vertices = mesh.vertices
    volume = 0.0
    mx = my = mz = 0.0
    emerged = 0
    for face in mesh.faces:
        points = (vertices[face[0]], vertices[face[1]], vertices[face[2]])
        distances = [offset - (nx * p[0] + ny * p[1] + nz * p[2]) for p in points]
        if distances[0] >= 0.0 and distances[1] >= 0.0 and distances[2] >= 0.0:
            polygon = list(points)
        elif distances[0] < 0.0 and distances[1] < 0.0 and distances[2] < 0.0:
            continue
        else:
            if distances[0] < 0.0 or distances[1] < 0.0 or distances[2] < 0.0:
                emerged += 1
            polygon = _clip_to_halfspace(points, distances)
        if len(polygon) < 3:
            continue
        ax, ay, az = polygon[0][0] - ox, polygon[0][1] - oy, polygon[0][2] - oz
        for k in range(1, len(polygon) - 1):
            bx, by, bz = polygon[k][0] - ox, polygon[k][1] - oy, polygon[k][2] - oz
            cx, cy, cz = polygon[k + 1][0] - ox, polygon[k + 1][1] - oy, polygon[k + 1][2] - oz
            det = (
                ax * (by * cz - bz * cy)
                - ay * (bx * cz - bz * cx)
                + az * (bx * cy - by * cx)
            )
            volume += det
            mx += det * (ax + bx + cx)
            my += det * (ay + by + cy)
            mz += det * (az + bz + cz)
    volume /= 6.0
    if volume == 0.0:
        raise HeelError(
            "the waterplane displaces no volume; the heeled plane does not cut the body"
        )
    return (
        volume,
        (
            ox + mx / (24.0 * volume),
            oy + my / (24.0 * volume),
            oz + mz / (24.0 * volume),
        ),
        emerged,
    )


def upright_displacement(
    mesh: Mesh, waterline_z: float = 0.0
) -> tuple[float, float, float, float]:
    """``(volume, x_B, y_B, z_B)`` below ``z = waterline_z`` with the same integrator."""
    volume, centre, _ = _immersed_volume_moments(mesh, (0.0, 0.0, 1.0), waterline_z)
    return volume, centre[0], centre[1], centre[2]


# ---------------------------------------------------------------------------
# the plane solve
# ---------------------------------------------------------------------------
def _solve_plane_offset(
    mesh: Mesh,
    normal: tuple[float, float, float],
    target_volume: float,
    initial_offset: float,
    scale: float,
    tolerance: float,
) -> tuple[float, int]:
    """Secant iteration for the offset whose immersed volume is ``target_volume``."""
    if target_volume <= 0.0:
        raise HeelError("target displacement volume must be positive")

    def residual(offset: float) -> float:
        volume, _, _ = _immersed_volume_moments(mesh, normal, offset)
        return volume - target_volume

    previous_offset = initial_offset
    previous = residual(initial_offset)
    if abs(previous) <= tolerance * target_volume:
        return initial_offset, 0

    delta = 0.05 * scale
    offset = initial_offset + delta
    current = residual(offset)
    iterations = 1
    while iterations < MAX_SOLVE_ITERATIONS:
        if abs(current) <= tolerance * target_volume:
            return offset, iterations
        denominator = current - previous
        if denominator == 0.0:  # pragma: no cover - needs a degenerate V(d)
            step = 0.5 * delta
        else:
            step = -current * (offset - previous_offset) / denominator
        # The secant is not bracketed; keep the walk inside a sane neighbourhood of
        # the upright waterplane so a kink cannot send it off the model.
        limit = 0.5 * scale
        step = max(-limit, min(limit, step))
        if step == 0.0:  # pragma: no cover - only on an exact non-monotone V(d)
            step = 0.25 * delta
        previous_offset, previous = offset, current
        offset = offset + step
        current = residual(offset)
        iterations += 1
    raise HeelError(
        f"the heeled plane solve did not reach {tolerance:g} relative in "
        f"{MAX_SOLVE_ITERATIONS} iterations (residual {current / target_volume:.3e})"
    )


@dataclass(frozen=True)
class HeeledDisplacement:
    """The heeled displacement: how much, where, and on which plane."""

    heel_angle_deg: float
    displacement_volume: float
    target_volume: float
    centre_of_buoyancy: tuple[float, float, float]
    waterplane_offset: float
    iterations: int
    keel_z: float
    yg: float
    schema: str = SCHEMA_HEELED

    @property
    def relative_volume_error(self) -> float:
        if self.target_volume == 0.0:  # pragma: no cover - guarded upstream
            return 0.0
        return (self.displacement_volume - self.target_volume) / self.target_volume

    @property
    def kb(self) -> float:
        """Height of the heeled ``B`` above the keel."""
        return self.centre_of_buoyancy[2] - self.keel_z

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": self.schema,
            "heel_angle_deg": _round(self.heel_angle_deg),
            "displacement_volume_m3": _round(self.displacement_volume),
            "target_volume_m3": _round(self.target_volume),
            "relative_volume_error": _round(self.relative_volume_error, 15),
            "centre_of_buoyancy_m": tuple(_round(value) for value in self.centre_of_buoyancy),
            "waterplane_offset_m": _round(self.waterplane_offset),
            "keel_z_m": _round(self.keel_z),
            "iterations": self.iterations,
        }


def heeled_displacement(
    mesh: Mesh,
    heel_angle: float,
    *,
    waterline_z: float = 0.0,
    displacement_volume: float | None = None,
    tolerance: float = VOLUME_TOLERANCE,
) -> HeeledDisplacement:
    """Displacement and centre of buoyancy of ``mesh`` heeled by ``heel_angle``.

    The displaced volume is held at its upright value (``displacement_volume``
    overrides it for a different loading condition) and the waterplane offset is
    solved for, which is the difference between "the ship heels" and "the tank
    tilts".
    """
    _check_angle(heel_angle)
    if not mesh.faces:
        raise HeelError("the mesh has no faces, so there is nothing to heel")
    if not math.isfinite(waterline_z):
        raise HeelError("waterline_z must be finite")

    scale = _scale(mesh)
    if displacement_volume is None:
        target, _, _, _ = upright_displacement(mesh, waterline_z)
    else:
        if not math.isfinite(displacement_volume) or displacement_volume <= 0.0:
            raise HeelError("displacement_volume must be a finite positive number")
        target = float(displacement_volume)

    normal = up_direction(heel_angle)
    guess = waterline_z * normal[2]
    offset, iterations = _solve_plane_offset(
        mesh, normal, target, guess, scale, tolerance
    )
    volume, centre, emerged = _immersed_volume_moments(mesh, normal, offset)
    if heel_angle != 0.0 and emerged == 0:
        raise HeelError(
            f"at {math.degrees(heel_angle):.6g} deg the heeled waterplane does not cut "
            f"the surface: the mesh has no freeboard above it, so the displaced volume "
            f"is only conserved by submerging the whole body. Model the hull to the "
            f"deck edge (a mesh capped at the design waterline cannot be heeled)."
        )
    return HeeledDisplacement(
        heel_angle_deg=math.degrees(heel_angle),
        displacement_volume=volume,
        target_volume=target,
        centre_of_buoyancy=centre,
        waterplane_offset=offset,
        iterations=iterations,
        keel_z=keel_z(mesh),
        yg=0.0,
    )


# ---------------------------------------------------------------------------
# the righting arm
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class RightingArm:
    """One point of the ``GZ`` curve, with the buoyancy it was computed from."""

    heel_angle_deg: float
    righting_arm: float
    displacement_volume: float
    centre_of_buoyancy: tuple[float, float, float]
    kg: float
    yg: float
    kb: float
    waterplane_offset: float
    relative_volume_error: float
    iterations: int
    basis: str = "heeled displacement, exact polytope integration"
    schema: str = SCHEMA_ARM

    @property
    def gz_over_sin(self) -> float:
        """``GZ / sin(phi)``: the metacentric height at small heel, a curve later."""
        angle = math.radians(self.heel_angle_deg)
        sine = math.sin(angle)
        if sine == 0.0:
            return float("nan")
        return self.righting_arm / sine

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": self.schema,
            "heel_angle_deg": _round(self.heel_angle_deg),
            "righting_arm_m": _round(self.righting_arm, 12),
            "displacement_volume_m3": _round(self.displacement_volume),
            "centre_of_buoyancy_m": tuple(_round(value) for value in self.centre_of_buoyancy),
            "kb_m": _round(self.kb),
            "kg_m": _round(self.kg),
            "yg_m": _round(self.yg),
            "waterplane_offset_m": _round(self.waterplane_offset),
            "relative_volume_error": _round(self.relative_volume_error, 15),
            "iterations": self.iterations,
            "basis": self.basis,
        }


def righting_arm(
    mesh: Mesh,
    heel_angle: float,
    kg: float,
    *,
    waterline_z: float = 0.0,
    yg: float = 0.0,
    displacement_volume: float | None = None,
    tolerance: float = VOLUME_TOLERANCE,
) -> RightingArm:
    """The righting arm ``GZ`` at ``heel_angle`` radians for a loading ``KG`` above the keel."""
    _check_loading(mesh, kg, yg)
    heeled = heeled_displacement(
        mesh,
        heel_angle,
        waterline_z=waterline_z,
        displacement_volume=displacement_volume,
        tolerance=tolerance,
    )
    x_b, y_b, z_b = heeled.centre_of_buoyancy
    z_g = heeled.keel_z + kg
    cos = math.cos(heel_angle)
    sin = math.sin(heel_angle)
    arm = (y_b - yg) * cos + (z_b - z_g) * sin
    return RightingArm(
        heel_angle_deg=heeled.heel_angle_deg,
        righting_arm=arm,
        displacement_volume=heeled.displacement_volume,
        centre_of_buoyancy=(x_b, y_b, z_b),
        kg=kg,
        yg=yg,
        kb=heeled.kb,
        waterplane_offset=heeled.waterplane_offset,
        relative_volume_error=heeled.relative_volume_error,
        iterations=heeled.iterations,
    )


# ---------------------------------------------------------------------------
# the curve
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class GZCurve:
    """A righting-arm curve and the four numbers a designer reads off it."""

    angles_deg: tuple[float, ...]
    arms: tuple[float, ...]
    kg: float
    keel_z: float
    maximum_arm: float
    angle_of_maximum_deg: float
    vanishing_angle_deg: float | None
    loll_angle_deg: float | None
    dynamical_stability_m_rad: float
    integrated_to_deg: float
    integration_rule: str
    gm_from_curve: float | None
    worst_relative_volume_error: float
    basis: str = "heeled displacement at constant displacement"
    schema: str = SCHEMA_CURVE

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": self.schema,
            "kg_m": _round(self.kg),
            "keel_z_m": _round(self.keel_z),
            "angles_deg": tuple(_round(value, 9) for value in self.angles_deg),
            "righting_arms_m": tuple(_round(value, 12) for value in self.arms),
            "maximum_arm_m": _round(self.maximum_arm, 12),
            "angle_of_maximum_deg": _round(self.angle_of_maximum_deg),
            "vanishing_angle_deg": (
                None if self.vanishing_angle_deg is None else _round(self.vanishing_angle_deg)
            ),
            "loll_angle_deg": (
                None if self.loll_angle_deg is None else _round(self.loll_angle_deg)
            ),
            "dynamical_stability_m_rad": _round(self.dynamical_stability_m_rad, 12),
            "integrated_to_deg": _round(self.integrated_to_deg),
            "integration_rule": self.integration_rule,
            "gm_from_curve_m": (
                None if self.gm_from_curve is None else _round(self.gm_from_curve, 12)
            ),
            "worst_relative_volume_error": _round(self.worst_relative_volume_error, 15),
            "basis": self.basis,
        }


def _integrate_arm(
    angles: Sequence[float],
    arms: Sequence[float],
    tail: tuple[float, float] | None = None,
) -> tuple[float, str]:
    """Composite Simpson on evenly spaced angles, trapezoid otherwise.

    ``tail`` is a final ``(angle, arm)`` point that does not fall on the grid -- the
    angle of vanishing stability is a bisected root, so its step is not the sampled
    one and the last interval is closed with a trapezoid.
    """
    # The arm is in metres and the abscissa is heel: metres times *radians*.
    angles = tuple(math.radians(angle) for angle in angles)
    arms = tuple(arms)
    if tail is not None:
        tail = (math.radians(tail[0]), tail[1])
    if len(angles) < 2:
        if tail is None:
            return 0.0, "none"
        return 0.5 * (arms[0] + tail[1]) * (tail[0] - angles[0]), "trapezoid"
    step = angles[1] - angles[0]
    uniform = all(
        abs((angles[i + 1] - angles[i]) - step) <= 1e-12 * max(1.0, abs(step))
        for i in range(len(angles) - 1)
    )
    if uniform and len(angles) >= 3 and len(angles) % 2 == 1:
        total = arms[0] + arms[-1]
        total += 4.0 * sum(arms[1:-1:2])
        total += 2.0 * sum(arms[2:-1:2])
        value = total * step / 3.0
        if tail is None:
            return value, "composite Simpson"
        return value + 0.5 * (arms[-1] + tail[1]) * (tail[0] - angles[-1]), (
            "composite Simpson + trapezoid tail"
        )
    total = 0.0
    for i in range(len(angles) - 1):
        total += 0.5 * (arms[i] + arms[i + 1]) * (angles[i + 1] - angles[i])
    rule = "trapezoid"
    if tail is not None:
        total += 0.5 * (arms[-1] + tail[1]) * (tail[0] - angles[-1])
        rule = "trapezoid + tail"
    return total, rule


def _crossing(
    mesh: Mesh,
    kg: float,
    waterline_z: float,
    yg: float,
    previous_angle: float,
    previous_arm: float,
    next_angle: float,
    next_arm: float,
) -> float:
    """Bisect the zero crossing of ``GZ`` between two angles of opposite sign."""
    low = math.radians(previous_angle)
    high = math.radians(next_angle)
    low_value = previous_arm
    for _ in range(200):
        mid = 0.5 * (low + high)
        arm = righting_arm(mesh, mid, kg, waterline_z=waterline_z, yg=yg).righting_arm
        if (arm < 0.0) == (next_arm < 0.0):
            high = mid
        else:
            low = mid
            low_value = arm
        if high - low <= 1e-13 * max(1.0, abs(high)):
            break
    del low_value
    return math.degrees(0.5 * (low + high))


def gz_curve(
    mesh: Mesh,
    angles_deg: Sequence[float] = tuple(float(a) for a in range(0, 91, 5)),
    kg: float = 0.0,
    *,
    waterline_z: float = 0.0,
    yg: float = 0.0,
    displacement_volume: float | None = None,
    tolerance: float = VOLUME_TOLERANCE,
) -> GZCurve:
    """Righting-arm curve over ``angles_deg`` (degrees, must be non-negative and sorted).

    Beside the arms it reports the peak, the angle of vanishing stability (the first
    zero crossing after the peak, by bisection; ``None`` if the curve is still positive
    at the last angle), the angle of loll (the first crossing *before* the peak, which
    exists when ``GM < 0``), the dynamical stability as the integral of ``GZ`` over heel in
    metre-radians, and ``gm_from_curve`` -- ``GZ / sin(phi)`` at the smallest non-zero
    angle, which converges to the metacentric height from
    :func:`tensorcad.stability.initial_stability` as that angle shrinks.
    """
    angles = tuple(float(angle) for angle in angles_deg)
    if not angles:
        raise HeelError("gz_curve needs at least one heel angle")
    for angle in angles:
        if not math.isfinite(angle) or angle < 0.0:
            raise HeelError("heel angles must be finite and non-negative (degrees)")
        _check_angle(math.radians(angle))
    points = [righting_arm(mesh, math.radians(angle), kg, waterline_z=waterline_z,
                           yg=yg, displacement_volume=displacement_volume,
                           tolerance=tolerance) for angle in angles]
    arms = tuple(point.righting_arm for point in points)
    peak_index = max(range(len(arms)), key=lambda index: arms[index])
    maximum = arms[peak_index]

    # A crossing has to be a crossing, not a round-off wobble around GZ(0) = 0.
    tiny = 1e-11 * max(1e-12, max(abs(arm) for arm in arms))

    vanishing: float | None = None
    for index in range(peak_index, len(arms) - 1):
        if arms[index] > tiny and arms[index + 1] < -tiny:
            vanishing = _crossing(
                mesh, kg, waterline_z, yg,
                angles[index], arms[index], angles[index + 1], arms[index + 1],
            )
            break
    # An upright-unstable ship (GM < 0) reaches positive GZ at the angle of loll,
    # which is the first crossing between zero heel and the peak.
    loll: float | None = None
    for index in range(0, peak_index):
        if (arms[index] < -tiny and arms[index + 1] > tiny) or (
            arms[index] > tiny and arms[index + 1] < -tiny
        ):
            loll = _crossing(
                mesh, kg, waterline_z, yg,
                angles[index], arms[index], angles[index + 1], arms[index + 1],
            )
            break
    stop = len(arms) - 1
    if vanishing is not None:
        # The convention: dynamical stability is the area under GZ up to the angle of
        # vanishing stability, not over the whole sampled range (which is negative
        # afterwards and would silently cancel the useful number).
        stop = max(index for index in range(len(angles)) if angles[index] <= vanishing)
    tail = (vanishing, 0.0) if vanishing is not None and angles[stop] < vanishing else None
    integral, rule = _integrate_arm(angles[: stop + 1], arms[: stop + 1], tail)

    gm_from_curve: float | None = None
    for index, angle in enumerate(angles):
        if angle > 0.0:
            gm_from_curve = arms[index] / math.sin(math.radians(angle))
            break

    return GZCurve(
        angles_deg=angles,
        arms=arms,
        kg=kg,
        keel_z=keel_z(mesh),
        maximum_arm=maximum,
        angle_of_maximum_deg=angles[peak_index],
        vanishing_angle_deg=vanishing,
        loll_angle_deg=loll,
        dynamical_stability_m_rad=integral,
        integrated_to_deg=angles[stop] if tail is None else tail[0],
        integration_rule=rule,
        gm_from_curve=gm_from_curve,
        worst_relative_volume_error=max(
            abs(point.relative_volume_error) for point in points
        ),
    )


# ---------------------------------------------------------------------------
# closed-form and independent references
# ---------------------------------------------------------------------------
def wall_sided_metacentric_height(beam: float, draft: float, kg: float) -> float:
    """``GM = KB + BM - KG`` of a rectangular section, ``KB = T/2``, ``BM = B^2/(12T)``."""
    _check_box(beam, draft)
    if not math.isfinite(kg):
        raise HeelError("KG must be finite")
    return 0.5 * draft + beam * beam / (12.0 * draft) - kg


def wall_sided_range(beam: float, draft: float) -> float:
    """Largest heel (radians) at which a rectangle's waterplane still cuts both sides.

    Beyond ``atan(2T/B)`` the waterline reaches the keel of the high side and the
    section is no longer wall-sided (a bilge turns, or a box's keel emerges).
    """
    _check_box(beam, draft)
    return math.atan(2.0 * draft / beam)


def wall_sided_righting_arm(
    beam: float, draft: float, kg: float, heel_angle: float, *, yg: float = 0.0
) -> float:
    """Exact ``GZ`` of a wall-sided body: ``GZ = sin(phi) (GM + BM tan^2(phi) / 2)``.

    This is not a small-angle truncation. Holding the displaced area at ``B*T`` gives a
    waterplane offset ``d = T cos(phi)`` exactly, the immersed trapezoid has
    ``y_B = B^2 tan(phi) / (12 T)`` and ``z_B = T/2 + B^2 tan^2(phi)/(24 T)`` in closed
    form, and projecting ``B - G`` on the heeled waterplane collapses to the formula
    above -- which is also why the classical "wall-sided" expression is used as a
    large-angle check rather than an approximation.
    """
    _check_angle(heel_angle)
    if not math.isfinite(yg):
        raise HeelError("y_G must be finite")
    limit = wall_sided_range(beam, draft)
    if abs(heel_angle) >= limit:
        raise HeelError(
            f"heel {math.degrees(heel_angle):.6g} deg is outside the wall-sided range "
            f"of {math.degrees(limit):.6g} deg for B/T = {beam / draft:.6g}"
        )
    bm = beam * beam / (12.0 * draft)
    gm = wall_sided_metacentric_height(beam, draft, kg)
    tangent = math.tan(heel_angle)
    return math.sin(heel_angle) * (gm + 0.5 * bm * tangent * tangent) - yg * math.cos(heel_angle)


def wall_sided_dynamical_stability(
    beam: float, draft: float, kg: float, heel_angle: float, *, yg: float = 0.0
) -> float:
    """Exact ``int_0^phi GZ dtheta`` (metre-radians) for a wall-sided body.

    ``int sin(t) dt = 1 - cos(t)`` and ``int sin^3(t)/cos^2(t) dt = sec(t) + cos(t) - 2``,
    so the wall-sided arm integrates in closed form and a quadrature of a *measured*
    curve can be checked against an exact area rather than against itself.
    """
    _check_angle(heel_angle)
    limit = wall_sided_range(beam, draft)
    if abs(heel_angle) >= limit:
        raise HeelError(
            f"heel {math.degrees(heel_angle):.6g} deg is outside the wall-sided range "
            f"of {math.degrees(limit):.6g} deg for B/T = {beam / draft:.6g}"
        )
    bm = beam * beam / (12.0 * draft)
    gm = wall_sided_metacentric_height(beam, draft, kg)
    cosine = math.cos(heel_angle)
    return (
        gm * (1.0 - cosine)
        + 0.5 * bm * (1.0 / cosine + cosine - 2.0)
        - yg * math.sin(heel_angle)
    )


def _check_box(beam: float, draft: float) -> None:
    for name, value in (("beam", beam), ("draft", draft)):
        if not math.isfinite(value) or value <= 0.0:
            raise HeelError(f"{name} must be a finite positive number")


# ---------------------------------------------------------------------------
# an independent reference: a barge section, clipped in two dimensions
# ---------------------------------------------------------------------------
def _clip_polygon_2d(
    polygon: Sequence[tuple[float, float]], normal: tuple[float, float], offset: float
) -> list[tuple[float, float]]:
    out: list[tuple[float, float]] = []
    count = len(polygon)
    distances = [offset - (normal[0] * y + normal[1] * z) for y, z in polygon]
    for index in range(count):
        following = (index + 1) % count
        s_here, s_there = distances[index], distances[following]
        if s_here >= 0.0:
            out.append(polygon[index])
        if (s_here > 0.0 > s_there) or (s_here < 0.0 < s_there):
            t = s_here / (s_here - s_there)
            here, there = polygon[index], polygon[following]
            out.append((here[0] + t * (there[0] - here[0]), here[1] + t * (there[1] - here[1])))
    return out


def _polygon_area_centroid(
    polygon: Sequence[tuple[float, float]]
) -> tuple[float, float, float]:
    area = 0.0
    cy = 0.0
    cz = 0.0
    count = len(polygon)
    for index in range(count):
        y0, z0 = polygon[index]
        y1, z1 = polygon[(index + 1) % count]
        cross = y0 * z1 - y1 * z0
        area += cross
        cy += (y0 + y1) * cross
        cz += (z0 + z1) * cross
    area *= 0.5
    if area == 0.0:  # pragma: no cover - only a degenerate cut
        raise HeelError("the heeled section has zero area")
    return area, cy / (6.0 * area), cz / (6.0 * area)


@dataclass(frozen=True)
class BargeSection:
    """A heeled rectangular section: area, centroid, and the plane it was cut with."""

    area: float
    y_centroid: float
    z_centroid: float
    waterplane_offset: float
    iterations: int
    keel_z: float
    deck_z: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "area_m2": _round(self.area),
            "y_centroid_m": _round(self.y_centroid),
            "z_centroid_m": _round(self.z_centroid),
            "waterplane_offset_m": _round(self.waterplane_offset),
            "iterations": self.iterations,
            "keel_z_m": _round(self.keel_z),
            "deck_z_m": _round(self.deck_z),
        }


def barge_section(
    beam: float,
    draft: float,
    depth: float,
    heel_angle: float,
    *,
    waterline_z: float = 0.0,
    tolerance: float = 1e-14,
) -> BargeSection:
    """Immerse a rectangular section to a fixed area and report its centroid.

    A barge is a prism, so its three-dimensional centre of buoyancy *is* the section
    centroid and its displacement is the area times the length. This path clips a
    quadrilateral and integrates by the shoelace formula -- a different algorithm from
    the tetrahedron flux used above, on the same polytope -- which is what makes it
    useful as a reference rather than a restatement.
    """
    _check_angle(heel_angle)
    _check_box(beam, draft)
    if not math.isfinite(depth) or depth <= draft:
        raise HeelError(
            "depth must be finite and greater than the draft: a section with no "
            "freeboard cannot be heeled at constant displacement"
        )
    if not math.isfinite(waterline_z):
        raise HeelError("waterline_z must be finite")

    half = 0.5 * beam
    keel = waterline_z - draft
    deck = keel + depth
    rectangle = [(-half, keel), (half, keel), (half, deck), (-half, deck)]
    target = beam * draft
    normal = (-math.sin(heel_angle), math.cos(heel_angle))

    def area(offset: float) -> tuple[float, tuple[float, float, float]]:
        return _polygon_area_centroid(_clip_polygon_2d(rectangle, normal, offset))

    guess = waterline_z * normal[1]
    if abs(area(guess)[0] - target) <= tolerance * target:
        offset, iterations = guess, 0
    else:
        previous_offset, previous = guess, area(guess)[0] - target
        offset = guess + 0.05 * depth
        current = area(offset)[0] - target
        iterations = 1
        while abs(current) > tolerance * target:
            if iterations >= MAX_SOLVE_ITERATIONS:  # pragma: no cover
                raise HeelError("the barge section solve did not converge")
            denominator = current - previous
            step = -current * (offset - previous_offset) / denominator if denominator else 0.05 * depth
            step = max(-0.5 * depth, min(0.5 * depth, step))
            previous_offset, previous = offset, current
            offset += step
            current = area(offset)[0] - target
            iterations += 1

    section_area, y_c, z_c = area(offset)
    return BargeSection(
        area=section_area,
        y_centroid=y_c,
        z_centroid=z_c,
        waterplane_offset=offset,
        iterations=iterations,
        keel_z=keel,
        deck_z=deck,
    )


def barge_righting_arm(
    beam: float,
    draft: float,
    depth: float,
    kg: float,
    heel_angle: float,
    *,
    waterline_z: float = 0.0,
    yg: float = 0.0,
) -> float:
    """``GZ`` of a rectangular barge from the independent section integration."""
    if not math.isfinite(kg):
        raise HeelError("KG must be finite")
    section = barge_section(beam, draft, depth, heel_angle, waterline_z=waterline_z)
    z_g = section.keel_z + kg
    return (section.y_centroid - yg) * math.cos(heel_angle) + (
        section.z_centroid - z_g
    ) * math.sin(heel_angle)


def barge_deck_edge_angle(beam: float, draft: float, depth: float) -> float:
    """Heel (radians) at which the low side's deck edge reaches the waterplane."""
    _check_box(beam, draft)
    freeboard = depth - draft
    if freeboard <= 0.0:
        return 0.0
    return math.atan(2.0 * freeboard / beam)


def _round(value: float, digits: int = 9) -> float:
    if not isinstance(value, float) or not math.isfinite(value):
        return value
    return round(value, digits)