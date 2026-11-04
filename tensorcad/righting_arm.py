"""Large-angle stability: the righting-arm curve ``GZ(phi)``.

``stability.py`` answers the upright question -- ``GM = KB + BM - KG`` -- and says in its
own docstring that it cannot answer the heeled one, because the two are different
integrations: ``GM`` is a *waterplane moment* taken at one waterline, ``GZ`` is a
*heeled displacement* of the whole body. This module is that second integration.

**Definition.** The body is rotated about a longitudinal axis through the centreline of
the model frame (``y = z = 0`` for a hull generated at its design waterline) by the heel
angle ``phi``; the waterline stays horizontal and the immersion is *solved* so that the
displaced volume equals the loading displacement; ``B`` is the centroid of that immersed
volume and ``G`` is the loading point at height ``KG`` above the keel. Then

    GZ(phi) = y_B(phi) - y_G(phi)

is the horizontal lever between the vertical through ``B`` and the vertical through
``G`` -- the righting arm, positive when it tends to bring the body upright. Heel is
positive to starboard (``+y``), and for ``GM > 0`` the curve leaves the origin with
slope ``GM``, which is the identity ``stability.py`` states without being able to check
it.

**Why the draft has to be solved.** A body that keeps its displacement while it heels
cannot keep its waterline: the low side goes deeper and the high side shallower, and the
displaced volume is preserved only by sinking or rising. ``GZ`` is therefore not a
property of a mesh rotated about a fixed waterline; the waterline is an unknown. This
module solves it by bisection on the immersed volume to a relative residual of ``1e-13``
and *reports the residual per angle* instead of assuming it away. The lever itself does
not care where the waterline ends up -- a vertical shift moves neither ``y_B`` nor
``y_G`` -- so the solve only has to make the volume right.

**How the immersed volume is integrated.** The closed surface is clipped with the
horizontal plane, and the moments of the part below it are read off the *clipped surface
alone*, with no cap. That is what the vector field ``F = (0, 0, z - w)`` buys:

    volume    = int int int dV        = oint (z - w) n_z dA
    int x dV  = int int int x dV      = oint x (z - w) n_z dA
    int y dV                          = oint y (z - w) n_z dA
    int z dV                          = w * volume + oint (z - w)^2 / 2 n_z dA

Every field is z-directed and every integrand carries the factor ``(z - w)``, so the flux
through the waterplane cap -- whose normal is ``z`` and on which ``z = w`` -- is
identically zero, and the cap never has to be built. (Putting ``x`` on the x-component
instead, ``F = (x (z - w), 0, 0)``, is the trap: its divergence is ``z - w`` and its flux
through the cap does not vanish, so the moment silently comes out as ``-w * volume``.) Clipping and then chaining the cut segments
into loops to triangulate the cap is the other route -- and it is the right one for
``stability.py``, which *needs* the waterplane area -- but here the cap is not wanted,
and leaving it out also removes the cases that a cap makes delicate: a face lying exactly
in the clipping plane, or the plane running exactly along a mesh edge, contributes
exactly zero to all four integrals above instead of doubling the cap.

**Exactness of each clipped polygon.** Clipping a triangle to a half-space gives a convex
polygon of at most four vertices. It is fan-triangulated and each triangle is integrated
with the three-point edge-midpoint rule, which is exact for the degree-two integrands
above (``z**2``, ``x*z``, ``y*z``) as well as for the linear ones. Nothing here
quadratures an approximation of the surface: the only discretisation is the one already
baked into the mesh, and ``examples/righting_arm`` measures it rather than assuming it
small.

**What this is not.** No free-surface correction (a slack tank is the caller's problem,
exactly as in ``stability.py``), no trim coupling (the body is heeled about a fixed axis;
a real heel also changes trim), no damage or flooding case, and no statutory verdict:
``IS Code 2008`` needs wind heeling arms, roll amplitudes and tank free surfaces, none of
which is here. :meth:`RightingArmCurve.range_of_stability` reports where the computed
``GZ`` first returns to zero; it is not a compliance statement.

**Geometry caveat that is easy to get wrong.** ``generate_hull`` and ``generate_wigley``
cap the surface *at the design waterline*, and a surface capped there has **zero
freeboard**: its displacement is the whole body volume, so no heeled waterline can hold
that displacement -- the only way to displace the entire body is to submerge it, and the
arm of a fully submerged body is identically zero. The shipped hull family therefore has
no large-angle curve at all, and asking for one raises instead of returning a row of
zeros. For a hull with topsides, generate the deeper body and float it at the lower
waterline (``HullParameters(draft=2 * t)`` with ``waterline_z=-t``), which is what
``examples/righting_arm`` does; the topsides are then the analytic continuation of the
same section law, and are reported as such rather than as a designed deck.

Only the Python standard library is used.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import Any

from .geometry import Mesh
from .mesh import is_watertight

#: Standard gravity, used only to turn a mass displacement into a righting moment.
GRAVITY = 9.80665

SCHEMA_IMMERSED = "tensorcad.righting_arm.immersed/1"
SCHEMA_ARM = "tensorcad.righting_arm.arm/1"
SCHEMA_CURVE = "tensorcad.righting_arm.curve/1"
SCHEMA_WALL_SIDED = "tensorcad.righting_arm.wall_sided/1"

#: Bisection stops when the displaced volume is this close to the target, relative to
#: the target. A tighter value buys nothing: the lever is a centroid coordinate, and
#: ``y_B`` is insensitive to the waterline to first order near a plane ship.
VOLUME_TOLERANCE = 1e-13

MAX_ITERATIONS = 200


class RightingArmError(ValueError):
    """Raised when a righting arm cannot be defined for the given geometry/loading."""


def _round(value: float, digits: int = 9) -> float:
    if not math.isfinite(value):
        return value
    return round(value, digits)


# ---------------------------------------------------------------------------
# immersed volume below a horizontal plane
# ---------------------------------------------------------------------------
def _clip_triangle_below(
    a: tuple[float, float, float],
    b: tuple[float, float, float],
    c: tuple[float, float, float],
    waterline_z: float,
) -> list[tuple[float, float, float]]:
    """Sutherland-Hodgman clip of one triangle to the half-space ``z <= waterline_z``.

    The vertex order of the input is preserved, so the polygon that comes back is
    wound the same way as the face it came from and its area vector still points along
    the outward normal. Vertices exactly on the plane are kept and no interpolation is
    attempted across an edge that lies in the plane, which is what keeps a face lying
    exactly in the clipping plane from producing a spurious sliver.
    """
    polygon: list[tuple[float, float, float]] = []
    triangle = (a, b, c)
    for index in range(3):
        p = triangle[index]
        q = triangle[(index + 1) % 3]
        dp = p[2] - waterline_z
        dq = q[2] - waterline_z
        if dp <= 0.0:
            polygon.append(p)
        if (dp < 0.0 < dq) or (dq < 0.0 < dp):
            t = dp / (dp - dq)
            polygon.append((p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1]), waterline_z))
    return polygon


def _polygon_moments(
    polygon: Sequence[tuple[float, float, float]],
) -> tuple[float, tuple[float, float, float], tuple[float, float, float, float, float, float]] | None:
    """Area vector and planar integrals of a convex, planar, outward-wound polygon.

    Returns ``(area, unit_normal, (Ix, Iy, Iz, Izz, Ixz, Iyz))`` where ``I`` are the
    area integrals ``int f dA`` of the six integrands the divergence theorem needs, or
    ``None`` when the polygon is degenerate. The integrals are exact for degree-two
    integrands: a fan from the first vertex splits the polygon into triangles, and the
    edge-midpoint rule is exact for quadratics on each triangle.
    """
    count = len(polygon)
    if count < 3:
        return None
    ax = ay = az = 0.0
    for index in range(count):
        x1, y1, z1 = polygon[index]
        x2, y2, z2 = polygon[(index + 1) % count]
        ax += y1 * z2 - z1 * y2
        ay += z1 * x2 - x1 * z2
        az += x1 * y2 - y1 * x2
    ax *= 0.5
    ay *= 0.5
    az *= 0.5
    area = math.sqrt(ax * ax + ay * ay + az * az)
    if not area > 0.0:
        return None

    ox, oy, oz = polygon[0]
    i_x = i_y = i_z = i_zz = i_xz = i_yz = 0.0
    for index in range(1, count - 1):
        px, py, pz = polygon[index]
        qx, qy, qz = polygon[index + 1]
        ux, uy, uz = px - ox, py - oy, pz - oz
        vx, vy, vz = qx - ox, qy - oy, qz - oz
        cx = uy * vz - uz * vy
        cy = uz * vx - ux * vz
        cz = ux * vy - uy * vx
        area_triangle = 0.5 * math.sqrt(cx * cx + cy * cy + cz * cz)
        if area_triangle == 0.0:
            continue
        weight = area_triangle / 3.0
        # edge midpoints (o,p), (p,q), (q,o); unrolled because this is the hot loop
        mx = 0.5 * (ox + px)
        my = 0.5 * (oy + py)
        mz = 0.5 * (oz + pz)
        i_x += weight * mx
        i_y += weight * my
        i_z += weight * mz
        i_zz += weight * mz * mz
        i_xz += weight * mx * mz
        i_yz += weight * my * mz
        mx = 0.5 * (px + qx)
        my = 0.5 * (py + qy)
        mz = 0.5 * (pz + qz)
        i_x += weight * mx
        i_y += weight * my
        i_z += weight * mz
        i_zz += weight * mz * mz
        i_xz += weight * mx * mz
        i_yz += weight * my * mz
        mx = 0.5 * (qx + ox)
        my = 0.5 * (qy + oy)
        mz = 0.5 * (qz + oz)
        i_x += weight * mx
        i_y += weight * my
        i_z += weight * mz
        i_zz += weight * mz * mz
        i_xz += weight * mx * mz
        i_yz += weight * my * mz
    return area, (ax / area, ay / area, az / area), (i_x, i_y, i_z, i_zz, i_xz, i_yz)


def _face_coordinates(mesh: Mesh) -> tuple[tuple[Any, Any, Any], ...]:
    """Face vertex triples, resolved once so a bisection does not re-index per step."""
    vertices = mesh.vertices
    return tuple((vertices[a], vertices[b], vertices[c]) for a, b, c in mesh.faces)


def _immersed_integrals(
    faces: Sequence[tuple[Any, Any, Any]], waterline_z: float
) -> tuple[float, float, float, float, int]:
    """Volume and first moments below ``waterline_z``, from the clipped surface alone.

    Returns ``(volume, int x dV, int y dV, int z dV, clipped_face_count)``. All four
    quantities are fluxes of z-directed fields ``(0, 0, g(z))`` with ``g`` carrying the
    factor ``(z - w)``, so the waterplane cap -- whose normal is ``z`` and on which
    ``z = w`` -- contributes exactly zero and is absent by construction.
    """
    volume = 0.0
    moment_x = 0.0
    moment_y = 0.0
    moment_z = 0.0
    clipped = 0
    for face in faces:
        polygon = _clip_triangle_below(face[0], face[1], face[2], waterline_z)
        if len(polygon) < 3:
            continue
        measured = _polygon_moments(polygon)
        if measured is None:
            continue
        area, (nx, ny, nz), (i_x, i_y, i_z, i_zz, i_xz, i_yz) = measured
        clipped += 1
        volume += nz * (i_z - waterline_z * area)
        moment_x += nz * (i_xz - waterline_z * i_x)
        moment_y += nz * (i_yz - waterline_z * i_y)
        moment_z += nz * 0.5 * (i_zz - 2.0 * waterline_z * i_z + waterline_z * waterline_z * area)
    moment_z += waterline_z * volume
    return volume, moment_x, moment_y, moment_z, clipped


@dataclass(frozen=True)
class ImmersedVolume:
    """The part of a closed body below a horizontal plane, integrated exactly.

    ``centre_of_buoyancy`` is the centroid of the immersed volume in the mesh frame, so
    it is the point the buoyancy force acts through. ``clipped_faces`` counts the faces
    the plane actually cut (a face entirely below it is used whole), which is the
    number to look at when a volume comes out odd.
    """

    volume: float
    centre_of_buoyancy: tuple[float, float, float]
    waterline_z: float
    clipped_faces: int
    basis: str = "clipped closed surface, divergence theorem with F = (0, 0, z - w)"
    schema: str = SCHEMA_IMMERSED

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": self.schema,
            "volume_m3": _round(self.volume),
            "centre_of_buoyancy_m": [_round(value) for value in self.centre_of_buoyancy],
            "waterline_z_m": _round(self.waterline_z),
            "clipped_faces": self.clipped_faces,
            "basis": self.basis,
        }


def immersed_volume(mesh: Mesh, waterline_z: float = 0.0, *, closed: bool | None = None) -> ImmersedVolume:
    """Volume and buoyancy centre of ``mesh`` below the plane ``z = waterline_z``.

    The mesh must be closed: a righting arm taken from an open surface is a number
    about nothing, so an open surface raises instead of returning a vector of NaNs.
    ``waterline_z`` may cut the surface anywhere between the keel and the top, which is
    what makes this usable for a *heeled* body -- :func:`tensorcad.hydrostatics.
    compute_hydrostatics` expects a surface already cut at its waterline.
    """
    if not math.isfinite(waterline_z):
        raise RightingArmError(f"waterline_z must be finite, got {waterline_z!r}")
    if closed is None:
        closed = is_watertight(mesh)
    if not closed:
        raise RightingArmError(
            "the surface is not closed (it has boundary edges): an immersed volume needs a "
            "watertight surface, otherwise the integral is the volume of nothing"
        )
    volume, moment_x, moment_y, moment_z, clipped = _immersed_integrals(
        _face_coordinates(mesh), waterline_z
    )
    if not math.isfinite(volume) or volume <= 0.0:
        raise RightingArmError(
            f"immersed volume below z={waterline_z!r} is {volume!r} m^3: the waterline is "
            "outside the body or the faces are wound inward"
        )
    centroid = (moment_x / volume, moment_y / volume, moment_z / volume)
    return ImmersedVolume(
        volume=volume,
        centre_of_buoyancy=centroid,
        waterline_z=waterline_z,
        clipped_faces=clipped,
    )


# ---------------------------------------------------------------------------
# heeled geometry
# ---------------------------------------------------------------------------
def heeled_mesh(mesh: Mesh, heel_deg: float) -> Mesh:
    """Rotate ``mesh`` by ``heel_deg`` about the longitudinal axis through ``y = z = 0``.

    Positive heel puts the starboard side (``+y``) down, which is the convention the
    arm signs in this module are written against:

        y' = y cos(phi) + z sin(phi)
        z' = -y sin(phi) + z cos(phi)
    """
    if not math.isfinite(heel_deg):
        raise RightingArmError(f"heel_deg must be finite, got {heel_deg!r}")
    angle = math.radians(heel_deg)
    cos = math.cos(angle)
    sin = math.sin(angle)
    vertices = tuple((x, y * cos + z * sin, -y * sin + z * cos) for x, y, z in mesh.vertices)
    return Mesh(vertices=vertices, faces=mesh.faces, units=mesh.units, coordinate_frame=mesh.coordinate_frame)


def _solve_waterline(
    faces: Sequence[tuple[Any, Any, Any]],
    target_volume: float,
    lower: float,
    upper: float,
    *,
    tolerance: float = VOLUME_TOLERANCE,
) -> tuple[float, float, int]:
    """Bisect the waterline height until the immersed volume matches ``target_volume``.

    ``volume(waterline)`` is monotone in the waterline height for a closed surface --
    raising the plane can only add to the region below it -- so bisection cannot fail
    for a target inside ``(0, total volume)``, which is why a bracketed root finder is
    preferred to a Newton iteration over a function whose slope jumps wherever the plane
    passes a vertex.
    """
    low_volume, _, _, _, _ = _immersed_integrals(faces, lower)
    high_volume, _, _, _, _ = _immersed_integrals(faces, upper)
    if target_volume > high_volume * (1.0 - 1e-12):
        raise RightingArmError(
            f"the requested displacement {target_volume!r} m^3 is the whole volume of the "
            f"body ({high_volume!r} m^3), so the body has no freeboard: a closed surface "
            "capped at its waterline cannot heel at constant displacement (any heel would "
            "have to hold the entire body below the surface, and the arm is then "
            "identically zero). Generate a deeper body and float it at a lower "
            "waterline_z, e.g. HullParameters(draft=2 * t) with waterline_z=-t"
        )
    if target_volume <= 0.0 or target_volume < low_volume:
        raise RightingArmError(
            f"target volume {target_volume!r} m^3 is outside the range the body can "
            f"displace in this frame ({low_volume!r} to {high_volume!r} m^3)"
        )
    middle = 0.5 * (lower + upper)
    volume = 0.0
    residual = math.inf
    iterations = 0
    for iterations in range(1, MAX_ITERATIONS + 1):
        middle = 0.5 * (lower + upper)
        volume = _immersed_integrals(faces, middle)[0]
        residual = volume - target_volume
        if abs(residual) <= tolerance * max(target_volume, 1e-30):
            return middle, residual, iterations
        if residual > 0.0:
            upper = middle
        else:
            lower = middle
    return middle, residual, iterations


# ---------------------------------------------------------------------------
# the curve
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class RightingArm:
    """One point of a righting-arm curve: a heel angle and the arm measured there."""

    heel_deg: float
    gz: float
    waterline_z: float
    displaced_volume: float
    target_volume: float
    volume_residual: float
    centre_of_buoyancy: tuple[float, float, float]
    centre_of_gravity: tuple[float, float, float]
    kg: float
    water_density: float
    iterations: int
    schema: str = SCHEMA_ARM

    @property
    def heel_rad(self) -> float:
        return math.radians(self.heel_deg)

    @property
    def displaced_mass(self) -> float:
        return self.displaced_volume * self.water_density

    @property
    def righting_moment(self) -> float:
        """``displacement * g * GZ`` in kN m, the moment the hull actually applies."""
        return self.displaced_mass * GRAVITY * self.gz / 1000.0

    @property
    def relative_volume_error(self) -> float:
        if self.target_volume == 0.0:
            return math.inf
        return self.volume_residual / self.target_volume

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": self.schema,
            "heel_deg": _round(self.heel_deg),
            "gz_m": _round(self.gz),
            "waterline_z_m": _round(self.waterline_z),
            "displaced_volume_m3": _round(self.displaced_volume),
            "target_volume_m3": _round(self.target_volume),
            "relative_volume_error": _round(self.relative_volume_error, 15),
            "centre_of_buoyancy_m": [_round(value) for value in self.centre_of_buoyancy],
            "centre_of_gravity_m": [_round(value) for value in self.centre_of_gravity],
            "kg_m": _round(self.kg),
            "righting_moment_kn_m": _round(self.righting_moment, 6),
            "bisection_iterations": self.iterations,
        }


@dataclass(frozen=True)
class RightingArmCurve:
    """A ``GZ`` curve: the arms at the requested heels plus the loading behind them.

    The aggregate properties (:meth:`max_gz`, :meth:`angle_of_max_gz`,
    :meth:`range_of_stability`, :meth:`area`) are read off the *computed grid*, not off a
    fitted curve; a grid that is too coarse to resolve the maximum says so by returning
    an angle equal to a grid point, which is why the grid is part of the record.
    """

    arms: tuple[RightingArm, ...]
    kg: float
    target_volume: float
    water_density: float
    waterline_z: float
    keel_z: float
    basis: str = "heeled displacement, exact clipped-surface moments"
    schema: str = SCHEMA_CURVE

    def __post_init__(self) -> None:
        if not self.arms:
            raise RightingArmError("a righting-arm curve needs at least one heel angle")

    @property
    def heel_deg(self) -> tuple[float, ...]:
        return tuple(arm.heel_deg for arm in self.arms)

    @property
    def gz(self) -> tuple[float, ...]:
        return tuple(arm.gz for arm in self.arms)

    def arm_at(self, heel_deg: float, *, tolerance: float = 1e-9) -> RightingArm:
        """The computed arm at a grid angle (no interpolation: use :meth:`gz_at`)."""
        for arm in self.arms:
            if abs(arm.heel_deg - heel_deg) <= tolerance:
                return arm
        raise RightingArmError(
            f"heel {heel_deg!r} deg is not on this curve's grid {self.heel_deg!r}"
        )

    def gz_at(self, heel_deg: float) -> float:
        """Linearly interpolated arm; outside the grid the nearest end is held."""
        points = sorted(zip(self.heel_deg, self.gz))
        if heel_deg <= points[0][0]:
            return points[0][1]
        if heel_deg >= points[-1][0]:
            return points[-1][1]
        for (x0, y0), (x1, y1) in zip(points, points[1:]):
            if x0 <= heel_deg <= x1:
                if x1 == x0:
                    return y1
                weight = (heel_deg - x0) / (x1 - x0)
                return y0 + weight * (y1 - y0)
        return points[-1][1]

    @property
    def max_gz(self) -> float:
        return max(self.gz)

    @property
    def angle_of_max_gz(self) -> float:
        arms = self.arms
        best = max(arms, key=lambda arm: arm.gz)
        return best.heel_deg

    def range_of_stability(self) -> float | None:
        """Heel at which ``GZ`` first returns to zero after the maximum, or ``None``.

        Linearly interpolated between the two grid points that bracket the crossing, so
        the answer is grid-limited, not exact. ``None`` means ``GZ`` never came back to
        zero inside the computed grid -- the range is at least the largest heel
        computed, and this method will not guess beyond it.
        """
        arms = self.arms
        peak = max(range(len(arms)), key=lambda index: arms[index].gz)
        if arms[peak].gz <= 0.0:
            return None
        for lower, upper in zip(arms[peak:], arms[peak + 1 :]):
            if upper.gz <= 0.0 < lower.gz:
                if lower.gz == upper.gz:
                    return upper.heel_deg
                weight = lower.gz / (lower.gz - upper.gz)
                return lower.heel_deg + weight * (upper.heel_deg - lower.heel_deg)
        return None

    def area(self, to_deg: float | None = None) -> float:
        """Trapezoidal area under ``GZ`` from the first heel to ``to_deg``, in m rad.

        This is the "dynamic stability" a statutory criterion integrates; it is a
        property of the computed grid (trapezoid on the grid steps), so a coarse grid
        under-resolves it. Angles beyond the last grid point are not extrapolated: the
        integration stops at the last point at or below ``to_deg``.
        """
        arms = sorted(self.arms, key=lambda arm: arm.heel_deg)
        total = 0.0
        for lower, upper in zip(arms, arms[1:]):
            if to_deg is not None and upper.heel_deg > to_deg:
                break
            total += 0.5 * (lower.gz + upper.gz) * math.radians(upper.heel_deg - lower.heel_deg)
        return total

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": self.schema,
            "kg_m": _round(self.kg),
            "keel_z_m": _round(self.keel_z),
            "design_waterline_z_m": _round(self.waterline_z),
            "target_volume_m3": _round(self.target_volume),
            "water_density_kg_per_m3": _round(self.water_density),
            "max_gz_m": _round(self.max_gz),
            "angle_of_max_gz_deg": _round(self.angle_of_max_gz),
            "range_of_stability_deg": None
            if self.range_of_stability() is None
            else _round(self.range_of_stability()),
            "area_under_gz_m_rad": _round(self.area()),
            "max_relative_volume_error": _round(
                max(abs(arm.relative_volume_error) for arm in self.arms), 15
            ),
            "arms": [arm.to_dict() for arm in self.arms],
            "basis": self.basis,
        }


def righting_arms(
    mesh: Mesh,
    *,
    kg: float,
    heel_angles: Iterable[float] = (0.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0),
    waterline_z: float = 0.0,
    displacement: float | None = None,
    water_density: float = 1025.0,
    keel_z: float | None = None,
    closed: bool | None = None,
    tolerance: float = VOLUME_TOLERANCE,
) -> RightingArmCurve:
    """Compute ``GZ(phi)`` for a closed body at a fixed displacement.

    Parameters
    ----------
    mesh:
        Closed body. It may be capped anywhere; the heeled waterline is solved inside,
        so the mesh does *not* have to be pre-cut at its waterline (unlike
        :func:`tensorcad.hydrostatics.compute_hydrostatics`).
    kg:
        Height of the centre of gravity **above the keel** (the mesh's lowest vertex),
        a required loading input for the same reason ``stability.py`` requires it: a
        default of zero would silently answer a different question.
    heel_angles:
        Heel angles in degrees, positive to starboard. Duplicates are rejected.
    waterline_z:
        Height of the *design* waterline in the mesh frame (``0`` for a hull generated
        at its waterline). The heeled waterline floats around this by the solve.
    displacement:
        Displaced volume to hold constant, in m^3. Defaults to the immersed volume at
        ``waterline_z``, i.e. a body floating at its design waterline.
    water_density:
        Only used for the reported mass and righting moment; the geometry is density
        free.
    keel_z:
        Height of the keel datum, defaulting to the mesh's lowest vertex.

    The volume residual of the draft solve is carried on every :class:`RightingArm`
    rather than being asserted, so a curve whose bisection failed to converge cannot be
    mistaken for a converged one.
    """
    if not math.isfinite(kg):
        raise RightingArmError(f"kg must be finite, got {kg!r}")
    angles = list(heel_angles)
    if not angles:
        raise RightingArmError("heel_angles must not be empty")
    for angle in angles:
        if not math.isfinite(angle):
            raise RightingArmError(f"heel angle must be finite, got {angle!r}")
        if abs(angle) > 180.0:
            raise RightingArmError(f"heel angle {angle!r} deg is outside [-180, 180]")
    if len(set(angles)) != len(angles):
        raise RightingArmError(f"heel_angles contains duplicates: {angles!r}")
    if closed is None:
        closed = is_watertight(mesh)
    if not closed:
        raise RightingArmError(
            "the surface is not closed (it has boundary edges): a heeled displacement "
            "needs a watertight body"
        )

    if keel_z is None:
        keel_z = min(vertex[2] for vertex in mesh.vertices)
    if not math.isfinite(keel_z):
        raise RightingArmError(f"keel_z must be finite, got {keel_z!r}")

    if displacement is None:
        # Not ``compute_hydrostatics``: that one integrates the *whole* closed surface and
        # is only correct for a mesh already cut at its waterline. This module's own
        # integrator is the general one, so the default displacement comes from it.
        displacement = immersed_volume(mesh, waterline_z, closed=True).volume
    if not math.isfinite(displacement) or displacement <= 0.0:
        raise RightingArmError(
            f"displacement {displacement!r} m^3 is not positive: a righting arm needs a "
            "floating body"
        )

    gravity_height = keel_z + kg
    arms: list[RightingArm] = []
    for angle in angles:
        rotated = heeled_mesh(mesh, angle)
        faces = _face_coordinates(rotated)
        box = rotated.bounding_box()
        lower = box.minimum[2]
        upper = box.maximum[2]
        waterline, residual, iterations = _solve_waterline(
            faces, displacement, lower, upper, tolerance=tolerance
        )
        volume, moment_x, moment_y, moment_z, _ = _immersed_integrals(faces, waterline)
        centroid = (moment_x / volume, moment_y / volume, moment_z / volume)
        angle_rad = math.radians(angle)
        # G is fixed in the body: at (0, 0, keel + KG) in the mesh frame, so after the
        # heel it sits at y = (keel + KG) sin(phi), z = (keel + KG) cos(phi).
        centre_of_gravity = (0.0, gravity_height * math.sin(angle_rad), gravity_height * math.cos(angle_rad))
        arms.append(
            RightingArm(
                heel_deg=angle,
                gz=centroid[1] - centre_of_gravity[1],
                waterline_z=waterline,
                displaced_volume=volume,
                target_volume=displacement,
                volume_residual=residual,
                centre_of_buoyancy=centroid,
                centre_of_gravity=centre_of_gravity,
                kg=kg,
                water_density=water_density,
                iterations=iterations,
            )
        )
    return RightingArmCurve(
        arms=tuple(arms),
        kg=kg,
        target_volume=displacement,
        water_density=water_density,
        waterline_z=waterline_z,
        keel_z=keel_z,
    )


def righting_arm(
    mesh: Mesh,
    *,
    kg: float,
    heel_deg: float,
    **kwargs: Any,
) -> RightingArm:
    """The arm at one heel angle; a one-point curve's only element."""
    curve = righting_arms(mesh, kg=kg, heel_angles=(heel_deg,), **kwargs)
    return curve.arms[0]


# ---------------------------------------------------------------------------
# closed-form reference: the wall-sided formula
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class WallSidedValidity:
    """Where the wall-sided formula is exact, and which edge leaves it first.

    The formula assumes the waterline still spans the full beam, so it is exact for a
    rectangular body only while neither the bilge has emerged on the high side nor the
    deck edge has immersed on the low side. Both limits are ``2 * height / beam``
    tangents, and the smaller one binds.
    """

    bilge_emergence_deg: float
    deck_immersion_deg: float
    binding_deg: float
    binding: str
    schema: str = SCHEMA_WALL_SIDED

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": self.schema,
            "bilge_emergence_deg": _round(self.bilge_emergence_deg),
            "deck_immersion_deg": _round(self.deck_immersion_deg),
            "binding_deg": _round(self.binding_deg),
            "binding": self.binding,
        }


def wall_sided_validity(*, beam: float, draft: float, depth: float) -> WallSidedValidity:
    """Heel limits of the wall-sided idealisation of a rectangular section.

    ``draft`` is the immersed depth and ``depth`` the full depth (so ``depth - draft`` is
    the freeboard). Raises when the draft does not fit in the depth, because then there
    is no upright waterline to speak of.
    """
    for name, value in (("beam", beam), ("draft", draft), ("depth", depth)):
        if not math.isfinite(value) or value <= 0.0:
            raise RightingArmError(f"{name} must be positive and finite, got {value!r}")
    if draft > depth:
        raise RightingArmError(f"draft {draft!r} m exceeds depth {depth!r} m")
    bilge = math.degrees(math.atan2(2.0 * draft, beam))
    freeboard = depth - draft
    deck = math.degrees(math.atan2(2.0 * freeboard, beam)) if freeboard > 0.0 else 0.0
    if bilge <= deck:
        return WallSidedValidity(bilge, deck, bilge, "bilge emergence")
    return WallSidedValidity(bilge, deck, deck, "deck immersion")


def wall_sided_righting_arm(heel_deg: float, *, kb: float, bm: float, kg: float) -> float:
    """``GZ = sin(phi) * (GM + BM * tan(phi)**2 / 2)``, the wall-sided formula.

    Exact for a body whose waterline keeps the full beam while it heels -- a rectangular
    barge with freeboard, inside :func:`wall_sided_validity` -- and the classic closed
    form to check a numerical righting arm against. ``kb`` and ``bm`` are the upright
    ``KB`` and ``BM`` measured from the keel, ``kg`` the loading height above the keel.
    """
    angle = math.radians(heel_deg)
    tangent = math.tan(angle)
    return math.sin(angle) * ((kb + bm - kg) + 0.5 * bm * tangent * tangent)


def box_reference_arm(
    heel_deg: float,
    *,
    length: float,
    beam: float,
    draft: float,
    kg: float,
) -> float:
    """Wall-sided arm of the rectangular barge ``L x B`` floating at ``draft``."""
    volume = length * beam * draft
    kb = 0.5 * draft
    bm = (length * beam**3 / 12.0) / volume
    return wall_sided_righting_arm(heel_deg, kb=kb, bm=bm, kg=kg)


def zero_crossings(heel_range_deg: Sequence[float], gz: Sequence[float]) -> tuple[float, ...]:
    """All heels at which ``gz`` crosses zero, linearly interpolated between samples.

    Both directions are returned in increasing heel order, because the two mean
    different things on a righting-arm curve: a downward crossing is the end of the
    range of stability, and an upward one is where a body with ``GM < 0`` comes to rest
    (the loll angle of the wall-sided idealisation). Reporting the crossings instead of
    the "vanishing angle" the caller has in mind keeps the caller from silently getting
    the wrong one.
    """
    if len(heel_range_deg) != len(gz):
        raise RightingArmError("heel_range_deg and gz must have the same length")
    crossings: list[float] = []
    for (x0, y0), (x1, y1) in zip(zip(heel_range_deg, gz), zip(heel_range_deg[1:], gz[1:])):
        if y0 == 0.0:
            crossings.append(x0)
            continue
        if (y0 < 0.0 < y1) or (y1 < 0.0 < y0):
            crossings.append(x0 + (x1 - x0) * (y0 / (y0 - y1)))
    if len(gz) and gz[-1] == 0.0:
        crossings.append(heel_range_deg[-1])
    return tuple(crossings)


def wall_sided_vanishing_angle(*, kb: float, bm: float, kg: float) -> float | None:
    """Heel at which the wall-sided arm returns to zero, ``tan(phi) = sqrt(-2 GM / BM)``.

    Exists only for ``GM < 0`` (a body that is unstable upright falls to this heel and
    stays there, the "angle of loll" of the wall-sided idealisation); returns ``None``
    when ``GM >= 0``, in which case the wall-sided curve never returns to zero.
    """
    gm = kb + bm - kg
    if gm >= 0.0:
        return None
    return math.degrees(math.atan(math.sqrt(-2.0 * gm / bm)))