"""Large-angle stability: the righting-arm curve ``GZ(phi)``.

:mod:`tensorcad.stability` reports ``GM``, which is the *slope* of the righting-arm
curve at the upright: ``GM = dGZ/dphi`` at ``phi = 0``. That is one number about one
point. The curve itself is what a stability booklet, an intact-stability check and a
damage-stability calculation consume: the maximum righting arm and the heel at which it
occurs, the range of positive stability, and the area under the curve up to 30 and 40
degrees. This module adds it.

Geometry and convention
-----------------------
The body is heeled by ``phi`` about the longitudinal (``x``) axis, the waterline stays
horizontal, and the displacement volume is held at its upright value -- heel at constant
displacement, which is what a stability calculation is. The rotation is

    y' =  y cos(phi) + z sin(phi)
    z' = -y sin(phi) + z cos(phi)

so **positive heel puts starboard (``+y``) down**. ``phi`` is in degrees at the API
boundary and radians inside.

The centre of gravity sits on the centreplane at ``keel_z + kg`` (no list, no trim), so
after the rotation its transverse ordinate is ``z_G sin(phi)``. The righting arm is the
transverse distance between the line of action of buoyancy (the vertical through the
centre of buoyancy) and ``G``:

    GZ(phi) = y_B(phi) - z_G sin(phi),      z_G = keel_z + kg

``GZ > 0`` is a righting arm, ``GZ < 0`` a heeling arm. At ``phi = 0`` this collapses to
``y_B`` of the upright body, which is zero for a symmetric hull and *not* zero for an
asymmetric one -- the module reports the list moment instead of hiding it.

Why the waterline is solved for, not assumed
-------------------------------------------
Holding the displacement volume fixed means the waterline height is an unknown: the
immersed volume of the heeled body must equal the upright displacement. It is found by a
bracketed root solve on the volume below a plane, not by an iterative sinkage correction
and not by assuming the waterline passes through the origin.

The immersed volume and centre of buoyancy come from the divergence theorem over the
*clipped* hull surface, with the apex of the tetrahedra placed **on the waterplane**::

    V   = (1/3) oint x.n dA
    M_i = (1/2) oint x_i^2 n_i dA

Every waterplane cap triangle passes through the apex, so its tetrahedron is degenerate:
the cap contributes exactly zero volume and zero first moment, and it never has to be
built. ``y_B`` therefore needs no cap polygon, no winding fix-up and no "close the loop"
step, and it is exact for the triangulation given -- the only error is the discretisation
of the hull surface itself.

What is not modelled
--------------------
* No free-surface correction (``GM`` is the solid value, labelled as such).
* No trim: heel only, about the longitudinal axis.
* The deck is whatever the mesh closes with. The parametric hull of
  :mod:`tensorcad.hull` is closed by a flat lid *on the design waterline*, i.e. it has
  **zero freeboard**, so a heeled condition immediately submerges part of the lid. The
  box and circular-cylinder references in ``examples/gz`` carry freeboard and are the
  analytic checks; the hull row there is a convergence and consistency check.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

from .geometry import Mesh

#: Default heel ladder for :func:`gz_curve`: the 5 degree stations a stability booklet
#: uses, run out to 60 degrees so the vanishing angle is bracketed for ordinary hulls.
DEFAULT_ANGLES_DEG: tuple[float, ...] = (0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 60.0)

#: Relative tolerance on the displacement balance (immersed volume vs target volume).
VOLUME_TOLERANCE = 1e-12

#: Cap on the Illinois iterations of the waterline solve. It converges in ~15 for
#: ordinary hulls; hitting the cap means the target volume is out of reach.
MAX_WATERLINE_ITERATIONS = 200

Point = tuple[float, float, float]


class GZError(ValueError):
    """Raised when a righting arm cannot be defined for the given input."""


# ---------------------------------------------------------------------------
# rotation
# ---------------------------------------------------------------------------
def _rotate_about_x(point: Point, angle_rad: float) -> Point:
    """Rotate ``point`` about the ``x`` axis: positive angle puts ``+y`` down."""
    cos_a = math.cos(angle_rad)
    sin_a = math.sin(angle_rad)
    x, y, z = point
    return (x, y * cos_a + z * sin_a, -y * sin_a + z * cos_a)


def heel_mesh(mesh: Mesh, heel_deg: float) -> Mesh:
    """The same body heeled by ``heel_deg`` about the ``x`` axis through the origin.

    Positive ``heel_deg`` puts positive ``y`` down, matching :func:`righting_arm`. Faces
    and vertex indices are untouched, so topology checks on the result are checks on the
    input.
    """
    if not isinstance(heel_deg, (int, float)) or not math.isfinite(heel_deg):
        raise GZError(f"heel angle must be a finite number, got {heel_deg!r}")
    angle = math.radians(heel_deg)
    return Mesh(
        vertices=tuple(_rotate_about_x(vertex, angle) for vertex in mesh.vertices),
        faces=mesh.faces,
        units=mesh.units,
        coordinate_frame=mesh.coordinate_frame,
    )


# ---------------------------------------------------------------------------
# volume and first moments of the part of a closed mesh below a plane
# ---------------------------------------------------------------------------
def _clip_below(polygon: Sequence[Point], waterline_z: float) -> list[Point]:
    """Sutherland-Hodgman clip of a convex polygon to ``z <= waterline_z``."""
    clipped: list[Point] = []
    count = len(polygon)
    for index in range(count):
        current = polygon[index]
        following = polygon[(index + 1) % count]
        current_in = current[2] <= waterline_z
        following_in = following[2] <= waterline_z
        if current_in:
            clipped.append(current)
        if current_in != following_in:
            span = following[2] - current[2]
            if span == 0.0:  # pragma: no cover - excluded by current_in != following_in
                continue
            t = (waterline_z - current[2]) / span
            clipped.append(
                (
                    current[0] + t * (following[0] - current[0]),
                    current[1] + t * (following[1] - current[1]),
                    waterline_z,
                )
            )
    return clipped


def _below_plane_moments(
    vertices: Sequence[Point], faces: Sequence[Sequence[int]], waterline_z: float
) -> tuple[float, float, float]:
    """``(V, M_y, M_z)`` of the material below ``z = waterline_z``.

    Uses ``V = (1/3) oint x.n dA`` and ``M_i = (1/2) oint x_i^2 n_i dA`` with the apex of
    each tetrahedron on the cutting plane, so the (unbuilt) cap cancels identically.
    """
    volume = 0.0
    moment_y = 0.0
    moment_z = 0.0
    for face in faces:
        p0 = vertices[face[0]]
        p1 = vertices[face[1]]
        p2 = vertices[face[2]]
        if p0[2] <= waterline_z and p1[2] <= waterline_z and p2[2] <= waterline_z:
            polygon: Sequence[Point] = (p0, p1, p2)
        elif p0[2] > waterline_z and p1[2] > waterline_z and p2[2] > waterline_z:
            continue
        else:
            polygon = _clip_below((p0, p1, p2), waterline_z)
            if len(polygon) < 3:
                continue
        head = polygon[0]
        for index in range(1, len(polygon) - 1):
            b = polygon[index]
            c = polygon[index + 1]
            ax = head[0]
            ay = head[1]
            az = head[2] - waterline_z
            bx = b[0]
            by = b[1]
            bz = b[2] - waterline_z
            cx = c[0]
            cy = c[1]
            cz = c[2] - waterline_z
            six_volume = ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx)
            tetra = six_volume / 6.0
            volume += tetra
            # apex ordinates are (0, 0, waterline_z): y_O = 0, z_O = waterline_z
            moment_y += tetra * (head[1] + b[1] + c[1]) / 4.0
            moment_z += tetra * (waterline_z + head[2] + b[2] + c[2]) / 4.0
    return volume, moment_y, moment_z


def immersed_moments(mesh: Mesh, waterline_z: float = 0.0) -> tuple[float, tuple[float, float]]:
    """Immersed volume and transverse/vertical centre of buoyancy below a plane.

    ``waterline_z`` is a world ordinate and the body is taken as modelled. A mesh that is
    not closed has no well-defined immersed volume -- the integral silently omits the
    waterplane cap -- so check it with :func:`tensorcad.mesh.is_watertight` first.
    """
    if not math.isfinite(waterline_z):
        raise GZError(f"waterline_z must be finite, got {waterline_z!r}")
    volume, moment_y, moment_z = _below_plane_moments(mesh.vertices, mesh.faces, waterline_z)
    if volume <= 0.0:
        raise GZError(
            f"immersed volume {volume!r} m^3 is not positive at waterline z={waterline_z!r}; "
            "the waterline is below the body or the mesh is inside-out"
        )
    return volume, (moment_y / volume, moment_z / volume)


def _total_volume(vertices: Sequence[Point], faces: Sequence[Sequence[int]]) -> float:
    volume = 0.0
    for face in faces:
        p0 = vertices[face[0]]
        p1 = vertices[face[1]]
        p2 = vertices[face[2]]
        volume += (
            p0[0] * (p1[1] * p2[2] - p1[2] * p2[1])
            - p0[1] * (p1[0] * p2[2] - p1[2] * p2[0])
            + p0[2] * (p1[0] * p2[1] - p1[1] * p2[0])
        ) / 6.0
    return volume


# ---------------------------------------------------------------------------
# the constant-displacement waterline
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class HeeledCondition:
    """The equilibrium of a heeled body at constant displacement."""

    heel_deg: float
    waterline_z: float
    volume: float
    centre_buoyancy: tuple[float, float]
    target_volume: float
    iterations: int

    @property
    def volume_error(self) -> float:
        """Relative error of the solved displacement volume against the target."""
        return self.volume / self.target_volume - 1.0

    def to_dict(self) -> dict[str, float]:
        return {
            "heel_deg": self.heel_deg,
            "waterline_z": self.waterline_z,
            "volume": self.volume,
            "target_volume": self.target_volume,
            "volume_error": self.volume_error,
            "centre_buoyancy_y": self.centre_buoyancy[0],
            "centre_buoyancy_z": self.centre_buoyancy[1],
            "iterations": self.iterations,
        }


def heel_state(
    mesh: Mesh,
    heel_deg: float,
    *,
    displacement_volume: float | None = None,
    waterline_z: float = 0.0,
    tolerance: float = VOLUME_TOLERANCE,
) -> HeeledCondition:
    """Heel ``mesh`` by ``heel_deg`` and solve for the constant-displacement waterline.

    ``displacement_volume`` defaults to the upright displacement, i.e. the immersed volume
    of ``mesh`` below ``waterline_z``. Pass it explicitly to study another loading
    condition: the heeled body then floats at that volume.
    """
    if not math.isfinite(tolerance) or tolerance <= 0.0:
        raise GZError(f"tolerance must be a positive finite number, got {tolerance!r}")
    upright_volume, _ = immersed_moments(mesh, waterline_z)
    target = upright_volume if displacement_volume is None else float(displacement_volume)
    if not math.isfinite(target) or target <= 0.0:
        raise GZError(f"displacement volume must be positive and finite, got {target!r}")

    heeled = heel_mesh(mesh, heel_deg)
    vertices = heeled.vertices
    faces = heeled.faces
    total = _total_volume(vertices, faces)
    if total <= 0.0:
        raise GZError(
            "mesh is inside-out (signed volume <= 0): the divergence integral would report "
            "the complement of the body; fix the winding first"
        )
    if target > total * (1.0 + 1e-12):
        raise GZError(
            f"displacement volume {target!r} m^3 exceeds the body volume {total!r} m^3; "
            "a body cannot displace more than its own volume"
        )
    if target >= total * (1.0 - 1e-9):
        raise GZError(
            f"displacement volume {target!r} m^3 is the whole body volume {total!r} m^3: "
            "the body floats exactly at full submergence, where the waterplane area vanishes and "
            "a righting arm is undefined. Give the body freeboard (a waterline below the lid) or "
            "pass a smaller displacement"
        )

    low = min(vertex[2] for vertex in vertices)
    high = max(vertex[2] for vertex in vertices)
    span = max(high - low, 1e-12)
    lo = low - span
    hi = high + 1e-6 * span
    f_lo = _below_plane_moments(vertices, faces, lo)[0] - target
    f_hi = _below_plane_moments(vertices, faces, hi)[0] - target
    if f_lo > 0.0 or f_hi < 0.0:  # pragma: no cover - the bracket follows from the bounds
        raise GZError("failed to bracket the equilibrium waterline")
    iterations = 0
    waterline = 0.5 * (lo + hi)
    volume = 0.0
    moment_y = 0.0
    moment_z = 0.0
    converged = False
    # Illinois (regula falsi with one-sided halving) on a monotone smooth root: a handful
    # of volume evaluations instead of ~50 bisections.
    while iterations < MAX_WATERLINE_ITERATIONS:
        iterations += 1
        denominator = f_hi - f_lo
        candidate = hi - f_hi * (hi - lo) / denominator if denominator != 0.0 else 0.5 * (lo + hi)
        waterline = candidate if lo < candidate < hi else 0.5 * (lo + hi)
        volume, moment_y, moment_z = _below_plane_moments(vertices, faces, waterline)
        residual = volume - target
        if abs(residual) <= tolerance * target:
            converged = True
            break
        if residual < 0.0:
            lo = waterline
            f_lo = residual
            f_hi *= 0.5
        else:
            hi = waterline
            f_hi = residual
            f_lo *= 0.5
    if not converged:  # pragma: no cover - only on pathological geometry
        raise GZError(
            f"the waterline solve did not reach {tolerance!r} in {MAX_WATERLINE_ITERATIONS} iterations"
        )
    return HeeledCondition(
        heel_deg=float(heel_deg),
        waterline_z=waterline,
        volume=volume,
        centre_buoyancy=(moment_y / volume, moment_z / volume),
        target_volume=target,
        iterations=iterations,
    )


def _gravity_z(mesh: Mesh, kg: float) -> float:
    if not isinstance(kg, (int, float)) or not math.isfinite(kg):
        raise GZError(f"kg must be a finite number, got {kg!r}")
    return min(vertex[2] for vertex in mesh.vertices) + float(kg)


def righting_arm(
    mesh: Mesh,
    heel_deg: float,
    *,
    kg: float,
    displacement_volume: float | None = None,
    waterline_z: float = 0.0,
) -> float:
    """The righting arm ``GZ`` in metres at ``heel_deg``; positive means righting.

    ``kg`` is the height of the centre of gravity above the keel (the mesh's lowest
    vertex), as in :func:`tensorcad.stability.initial_stability`: a loading input, hence a
    required keyword. ``G`` is taken on the centreplane.
    """
    z_g = _gravity_z(mesh, kg)
    state = heel_state(
        mesh,
        heel_deg,
        displacement_volume=displacement_volume,
        waterline_z=waterline_z,
    )
    return state.centre_buoyancy[0] - z_g * math.sin(math.radians(heel_deg))

@dataclass(frozen=True)
class GZCurve:
    """A sampled righting-arm curve plus the quantities read off it.

    ``conditions`` keeps the per-angle waterline solve (waterline height, volume error,
    centre of buoyancy) so the curve can be audited rather than only plotted: a curve is
    only as good as the displacement balance at each station, and that balance is reported.
    """

    angles_deg: tuple[float, ...]
    arms: tuple[float, ...]
    kg: float
    keel_z: float
    displacement_volume: float
    conditions: tuple[HeeledCondition, ...] = ()

    def __post_init__(self) -> None:
        if len(self.angles_deg) != len(self.arms):
            raise GZError("angles and arms must have the same length")
        for left, right in zip(self.angles_deg, self.angles_deg[1:]):
            if not right > left:
                raise GZError("heel angles must be strictly increasing")
        for angle in self.angles_deg:
            if not math.isfinite(angle):
                raise GZError(f"heel angles must be finite, got {angle!r}")

    @property
    def z_g(self) -> float:
        """Height of the centre of gravity above the keel datum used here."""
        return self.keel_z + self.kg

    @property
    def max_arm(self) -> tuple[float, float]:
        """``(heel_deg, GZ)`` of the largest sampled arm."""
        index = max(range(len(self.arms)), key=lambda i: self.arms[i])
        return (self.angles_deg[index], self.arms[index])

    @property
    def max_volume_error(self) -> float:
        """Largest relative displacement error over the sampled stations."""
        return max((abs(condition.volume_error) for condition in self.conditions), default=0.0)

    def arm_at(self, angle_deg: float) -> float:
        """Linearly interpolated arm; refuses to extrapolate past the sampled range."""
        if not self.arms:
            raise GZError("the curve is empty")
        if angle_deg < self.angles_deg[0] or angle_deg > self.angles_deg[-1]:
            raise GZError(
                f"heel {angle_deg!r} deg is outside the sampled range "
                f"[{self.angles_deg[0]!r}, {self.angles_deg[-1]!r}]"
            )
        for index in range(len(self.angles_deg) - 1):
            left = self.angles_deg[index]
            right = self.angles_deg[index + 1]
            if left <= angle_deg <= right:
                span = right - left
                weight = 0.0 if span == 0.0 else (angle_deg - left) / span
                return self.arms[index] + weight * (self.arms[index + 1] - self.arms[index])
        return self.arms[-1]  # pragma: no cover - excluded by the range check

    def vanishing_angle(self, *, tolerance: float = 0.0) -> float | None:
        """The heel where the arm comes back down to ``-tolerance``, interpolated.

        ``None`` when the curve is still above ``tolerance`` at the last station -- the
        honest answer, and the reason the default ladder runs to 60 degrees.
        """
        for index in range(1, len(self.arms)):
            left = self.arms[index - 1]
            right = self.arms[index]
            if left > tolerance and right <= tolerance:
                span = left - right
                weight = 0.0 if span == 0.0 else (left - tolerance) / span
                return self.angles_deg[index - 1] + weight * (
                    self.angles_deg[index] - self.angles_deg[index - 1]
                )
        return None

    def area_under_curve(self, upto_deg: float | None = None) -> float:
        """Trapezoid area under the arm curve in m*rad up to ``upto_deg``.

        The trapezoid rule on a 5 degree ladder carries an O(h^2) error; the samples that
        produced it are part of the curve, so the spacing is auditable. A limit between
        two stations clips the last panel linearly.
        """
        if len(self.arms) < 2:
            raise GZError("an area needs at least two samples")
        limit = self.angles_deg[-1] if upto_deg is None else float(upto_deg)
        if limit < self.angles_deg[0]:
            raise GZError(f"upto_deg={limit!r} is below the first sampled station")
        area = 0.0
        for index in range(len(self.angles_deg) - 1):
            left = self.angles_deg[index]
            right = self.angles_deg[index + 1]
            if limit <= left:
                break
            end = min(limit, right)
            arm_left = self.arms[index]
            arm_right = self.arms[index + 1]
            if end < right:
                arm_right = arm_left + (arm_right - arm_left) * (end - left) / (right - left)
            area += 0.5 * (arm_left + arm_right) * math.radians(end - left)
        return area

    def gm_from_slope(self, angle_deg: float | None = None) -> float:
        """``GM`` as the slope of the curve at the upright: ``GZ(phi) / sin(phi)``.

        Evaluated at the smallest strictly positive station by default; the difference to
        the waterplane integral of :func:`tensorcad.stability.initial_stability` is
        O(phi^2) and is measured in ``examples/gz``.
        """
        if angle_deg is None:
            candidates = [index for index, angle in enumerate(self.angles_deg) if angle > 0.0]
            if not candidates:
                raise GZError("a slope needs at least one strictly positive heel angle")
            angle_deg = self.angles_deg[candidates[0]]
        angle = float(angle_deg)
        if angle == 0.0:
            raise GZError("the slope at zero heel needs a non-zero station")
        return self.arm_at(angle) / math.sin(math.radians(angle))

    def to_dict(self) -> dict[str, object]:
        max_angle, max_value = self.max_arm
        vanishing = self.vanishing_angle()
        return {
            "kg": self.kg,
            "keel_z": self.keel_z,
            "z_g": self.z_g,
            "displacement_volume": self.displacement_volume,
            "angles_deg": list(self.angles_deg),
            "arms": list(self.arms),
            "max_arm": max_value,
            "angle_of_max_arm": max_angle,
            "vanishing_angle": vanishing,
            "area_to_30_deg": self.area_under_curve(30.0) if self.angles_deg[-1] >= 30.0 else None,
            "area_to_40_deg": self.area_under_curve(40.0) if self.angles_deg[-1] >= 40.0 else None,
            "gm_from_slope": self.gm_from_slope(),
            "max_volume_error": self.max_volume_error,
            "conditions": [condition.to_dict() for condition in self.conditions],
        }


def gz_curve(
    mesh: Mesh,
    angles_deg: Sequence[float] = DEFAULT_ANGLES_DEG,
    *,
    kg: float,
    displacement_volume: float | None = None,
    waterline_z: float = 0.0,
    tolerance: float = VOLUME_TOLERANCE,
) -> GZCurve:
    """Sample the righting-arm curve of a closed mesh over ``angles_deg``.

    Each station is an independent constant-displacement waterline solve, so the curve
    carries its own convergence evidence (``GZCurve.max_volume_error``) instead of being
    a smooth interpolation of two end states.
    """
    angles = tuple(float(angle) for angle in angles_deg)
    if not angles:
        raise GZError("at least one heel angle is required")
    z_g = _gravity_z(mesh, kg)
    keel_z = min(vertex[2] for vertex in mesh.vertices)
    arms: list[float] = []
    conditions: list[HeeledCondition] = []
    displacement: float | None = displacement_volume
    for angle in angles:
        state = heel_state(
            mesh,
            angle,
            displacement_volume=displacement,
            waterline_z=waterline_z,
            tolerance=tolerance,
        )
        if displacement is None:
            displacement = state.target_volume
        arms.append(state.centre_buoyancy[0] - z_g * math.sin(math.radians(angle)))
        conditions.append(state)
    return GZCurve(
        angles_deg=angles,
        arms=tuple(arms),
        kg=float(kg),
        keel_z=keel_z,
        displacement_volume=float(displacement if displacement is not None else 0.0),
        conditions=tuple(conditions),
    )
