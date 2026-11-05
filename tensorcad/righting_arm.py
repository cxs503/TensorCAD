"""Righting-arm curve ``GZ(phi)``: what a metacentric height only approximates.

:mod:`tensorcad.stability` measures ``GM``, the *slope* of the righting-arm curve at
``phi = 0``, and says in as many words that a curve for finite heel is a different
integration: at heel the displaced volume is no longer the upright volume below a
horizontal waterline, and the buoyancy centre has to be found where the *heeled* hull
displaces the same volume. This module is that integration.

Two things make it exact rather than sampled
-------------------------------------------
* **The immersed volume and its centroid come from the divergence theorem on the hull
  surface, with no cap and no station sum.** For the region ``z <= h`` of a closed
  surface, ``V = oint (z-h) n_z dA``, ``int y dV = oint (1/2) y^2 n_y dA`` and
  ``int z dV = oint (1/2)(z-h)^2 n_z dA + h V``: the waterplane cap contributes nothing
  to any of them (its normal is ``(0,0,1)`` and ``z = h`` on it), so the cut never has
  to be stitched into loops and the result does not depend on topology. Every integrand
  is a polynomial of degree <= 2 over a planar triangle, so each contribution is
  **exact** for the polygonal hull; what remains is the mesh's own discretisation.
  A face lying *in* the waterplane contributes exactly zero, which is why a hull capped
  at its design waterline needs no special case.
* **Equilibrium is solved, not assumed.** The waterline height ``h`` is the root of
  ``V(h) = displacement`` for the heeled attitude; ``V(h)`` is strictly increasing in
  ``h`` (raising the plane only adds volume), so a bracket-and-bisect solve cannot fail
  once the upright displacement is inside ``(0, V_total)``. The residual is reported
  instead of being left to be trusted.

The transverse righting arm is the horizontal distance between the vertical through the
centre of gravity and the vertical through the centre of buoyancy, both in the world
frame of the heeled body:

    GZ(phi) = y_B(phi) - y_G(phi),     y_G = (z_keel + KG) sin(phi)

with ``phi > 0`` heeled to starboard (the body is rotated by ``R_x(-phi)``, so that the
masthead moves to ``+y``). The sign is the usual one: ``GZ > 0`` is a restoring arm.

Sign and magnitude are pinned by a box, where the closed form is exact
----------------------------------------------------------------------
For a wall-sided barge of breadth ``B``, draft ``T`` and ``BM = B^2/(12 T)``,

    GZ(phi) = sin(phi) (GM + (1/2) BM tan^2(phi)),    GM = KB + BM - KG

which is *not* an approximation for a box: it follows from ``y_B = BM tan(phi)`` and
``z_B = KB + (1/2) BM tan^2(phi)`` at constant displacement, and holds until the
waterline reaches the deck edge (``tan phi = 2 (D - T) / B``) or uncovers the bilge
(``tan phi = 2 T / B``). :func:`wall_sided_righting_arm` returns that form and
:func:`wall_sided_validity_deg` returns where it stops being exact, so the comparison in
``examples/righting_arm`` is a verification and not a fit.

**What this is not.** No free-surface correction (a slack tank's ``sum(rho i_x) / disp``
must still be subtracted by the caller), no damage or flooded-compartment case, no trim
or sinkage coupling (the heel is prescribed, the waterline height is solved -- a free
*trim* solve is a second root and is not here), no ``KG`` default and no class-rule
criterion. ``KG`` is a loading input, never inferred.

Only the Python standard library is used.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from itertools import pairwise
from typing import Any

from .geometry import Mesh, Vec3
from .mesh import is_watertight

SCHEMA_HEELED = "tensorcad.stability.heeled/1"
SCHEMA_RIGHTING_ARM = "tensorcad.stability.gz/1"
SCHEMA_CURVE = "tensorcad.stability.gzcurve/1"

#: Heels a full curve is normally reported on: dense where ``GZ`` is nearly linear.
DEFAULT_HEELS_DEG = (0.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0)

#: Beyond this the transverse convention (a heel about the longitudinal axis) stops
#: describing the problem and the caller wants a different analysis, not a bigger table.
MAX_HEEL_DEG = 90.0

_DEFAULT_RELATIVE_TOLERANCE = 1e-12
_MAX_ITERATIONS = 200


class RightingArmError(ValueError):
    """Raised when a righting arm cannot be computed from the inputs given."""


def _round(value: float, digits: int = 9) -> float:
    return float(f"{value:.{digits}f}")


# ---------------------------------------------------------------------------
# surface integrals of a closed mesh below a plane
# ---------------------------------------------------------------------------
def _triangle_terms(
    a: Vec3, b: Vec3, c: Vec3
) -> tuple[float, float, float, float, float] | None:
    """Moments of one planar triangle, as the five numbers the integral needs.

    Returns ``(A n_z, A n_z mean_z, Q_z, Q_x, Q_y)`` with ``Q_a = (1/2) A n_a mean2_a``:
    the linear coefficient the volume needs and the three ``1/2 int a^2 n_a dA``
    moments, evaluated once per triangle. ``None`` for a degenerate (zero-area)
    triangle, which contributes nothing.
    """
    e1 = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
    e2 = (c[0] - a[0], c[1] - a[1], c[2] - a[2])
    nx = e1[1] * e2[2] - e1[2] * e2[1]
    ny = e1[2] * e2[0] - e1[0] * e2[2]
    nz = e1[0] * e2[1] - e1[1] * e2[0]
    norm = math.sqrt(nx * nx + ny * ny + nz * nz)
    if norm <= 0.0 or not math.isfinite(norm):
        return None
    area = 0.5 * norm
    ux, uy, uz = nx / norm, ny / norm, nz / norm
    mean_z = (a[2] + b[2] + c[2]) / 3.0
    mean2_x = (a[0] * a[0] + b[0] * b[0] + c[0] * c[0]
               + a[0] * b[0] + b[0] * c[0] + c[0] * a[0]) / 6.0
    mean2_y = (a[1] * a[1] + b[1] * b[1] + c[1] * c[1]
               + a[1] * b[1] + b[1] * c[1] + c[1] * a[1]) / 6.0
    mean2_z = (a[2] * a[2] + b[2] * b[2] + c[2] * c[2]
               + a[2] * b[2] + b[2] * c[2] + c[2] * a[2]) / 6.0
    return (
        area * uz,
        area * uz * mean_z,
        0.5 * area * uz * mean2_z,
        0.5 * area * ux * mean2_x,
        0.5 * area * uy * mean2_y,
    )


def _clip_below(a: Vec3, b: Vec3, c: Vec3, height: float) -> list[Vec3]:
    """Sutherland-Hodgman clip of a triangle to ``z <= height``.

    A vertex exactly on the plane is kept once and no intersection point is emitted for
    an edge that merely touches the plane, so the polygon cannot carry duplicate points.
    """
    polygon = (a, b, c)
    out: list[Vec3] = []
    for i in range(3):
        current = polygon[i]
        following = polygon[(i + 1) % 3]
        dc = current[2] - height
        df = following[2] - height
        if dc <= 0.0:
            out.append(current)
        if (dc < 0.0 < df) or (df < 0.0 < dc):
            t = dc / (dc - df)
            out.append((
                current[0] + t * (following[0] - current[0]),
                current[1] + t * (following[1] - current[1]),
                height,
            ))
    return out


class _HeeledSection:
    """A mesh rotated to a prescribed heel, ready to be cut at any waterline height.

    The coordinates are shifted by the mesh centre first: the divergence-theorem
    integrals are translation invariant, but a hull sitting far from the origin would
    lose digits to cancellation, and this costs one subtraction.
    """

    __slots__ = (
        "faces",
        "pivot_z",
        "shift",
        "terms",
        "vertices",
        "z_max",
        "z_min",
        "z_span",
    )

    def __init__(self, mesh: Mesh, heel_deg: float, pivot_z: float = 0.0) -> None:
        # The heel is a rotation about the longitudinal axis *through the upright
        # waterline*, not about the coordinate origin: a hull whose z origin is not its
        # waterline would otherwise be swung through a large arc (which leaves GZ
        # invariant but makes the reported waterline and buoyancy centre meaningless).
        phi = math.radians(heel_deg)
        cos_phi, sin_phi = math.cos(phi), math.sin(phi)
        rotated = [
            (x, y * cos_phi + (z - pivot_z) * sin_phi, -y * sin_phi + (z - pivot_z) * cos_phi)
            for (x, y, z) in mesh.vertices
        ]
        self.pivot_z = pivot_z
        xs = [p[0] for p in rotated]
        ys = [p[1] for p in rotated]
        zs = [p[2] for p in rotated]
        shift = (0.5 * (min(xs) + max(xs)), 0.5 * (min(ys) + max(ys)),
                 0.5 * (min(zs) + max(zs)))
        vertices = [(p[0] - shift[0], p[1] - shift[1], p[2] - shift[2]) for p in rotated]
        self.shift = shift
        self.vertices = vertices
        self.faces = mesh.faces
        # Brackets are kept in the *unshifted* frame, which is what integrals() takes:
        # storing them shifted would silently break for a mesh far from the origin.
        self.z_min = min(zs)
        self.z_max = max(zs)
        self.z_span = max(self.z_max - self.z_min, 0.0)
        self.terms = [
            _triangle_terms(vertices[i], vertices[j], vertices[k])
            for (i, j, k) in mesh.faces
        ]

    def integrals(self, height: float) -> tuple[float, float, float, float]:
        """``(volume, int x dV, int y dV, int z dV)`` for the part below ``height``.

        ``height`` is measured from the upright waterline (the pivot), in metres.
        """
        h = height - self.shift[2]
        vertices = self.vertices
        volume = 0.0
        moment_x = 0.0
        moment_y = 0.0
        moment_z = 0.0
        for face, terms in zip(self.faces, self.terms, strict=True):
            if terms is None:
                continue
            a = vertices[face[0]]
            b = vertices[face[1]]
            c = vertices[face[2]]
            if a[2] > h and b[2] > h and c[2] > h:
                continue
            az, az_mean_z, qz, qx, qy = terms
            if a[2] <= h and b[2] <= h and c[2] <= h:
                # The whole triangle is inside: the coefficients are already summed,
                # so the cap-free forms below are the closed-form contributions.
                volume += az_mean_z - h * az
                moment_x += qx
                moment_y += qy
                moment_z += qz - 0.5 * h * h * az
                continue
            polygon = _clip_below(a, b, c, h)
            if len(polygon) < 3:
                continue
            first = polygon[0]
            for index in range(1, len(polygon) - 1):
                sub = _triangle_terms(first, polygon[index], polygon[index + 1])
                if sub is None:
                    continue
                s_az, s_az_mean_z, s_qz, s_qx, s_qy = sub
                face_volume = s_az_mean_z - h * s_az
                volume += face_volume
                moment_x += s_qx
                moment_y += s_qy
                moment_z += s_qz - h * s_az_mean_z + 0.5 * h * h * s_az + h * face_volume
        return volume, moment_x, moment_y, moment_z


# ---------------------------------------------------------------------------
# heeled displacement at constant volume
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class HeeledDisplacement:
    """Displaced volume and buoyancy centre of a hull held at a prescribed heel."""

    heel_deg: float
    waterline_z: float
    volume: float
    centre_of_buoyancy: Vec3
    iterations: int
    volume_residual: float
    basis: str = "divergence theorem on the heeled hull, waterline solved for volume"

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA_HEELED,
            "heel_deg": _round(self.heel_deg),
            "waterline_z": _round(self.waterline_z, 12),
            "volume": _round(self.volume, 12),
            "centre_of_buoyancy": [_round(value, 12) for value in self.centre_of_buoyancy],
            "iterations": self.iterations,
            "volume_residual": self.volume_residual,
            "basis": self.basis,
        }


def _require_closed(mesh: Mesh, closed: bool | None) -> None:
    if closed is False:
        raise RightingArmError(
            "the immersed-volume integral is a closed-surface integral; an open mesh "
            "would silently integrate a part of a hull"
        )
    if closed is None and not is_watertight(mesh):
        raise RightingArmError(
            "the mesh is not watertight (every edge must be shared by exactly two "
            "faces), so the displaced volume below the waterline is not defined; "
            "pass closed=True only if you know the surface closes"
        )


def _upright_volume(mesh: Mesh, waterline_z: float) -> float:
    """Volume the hull displaces at the upright waterline, integrated here.

    Deliberately *not* :func:`tensorcad.hydrostatics.compute_hydrostatics`: that
    function integrates the whole closed mesh, which is the displacement only for a
    hull capped at its waterline. A mesh carrying freeboard would hand back its total
    volume instead of the immersed part, and every heeled point after it would be
    solved for the wrong displacement. The two agree exactly on a capped hull, which is
    a test rather than an assumption.
    """
    volume = _HeeledSection(mesh, 0.0, waterline_z).integrals(0.0)[0]
    if not math.isfinite(volume) or volume <= 0.0:
        raise RightingArmError(
            f"the upright displacement at waterline_z={waterline_z!r} is {volume!r} m^3; "
            "a heeled equilibrium needs a floating body"
        )
    return volume


def _solve_waterline(
    section: _HeeledSection,
    target: float,
    *,
    max_iterations: int,
    relative_tolerance: float,
) -> tuple[float, int]:
    """Bisect the waterline height until the immersed volume equals ``target``."""
    total = section.integrals(section.z_max)[0]
    if target > total * (1.0 + 1e-12):
        raise RightingArmError(
            f"requested displacement {target!r} m^3 exceeds the hull's whole volume "
            f"{total!r} m^3: the body would be fully submerged, not floating"
        )
    low = section.z_min
    high = section.z_max
    tolerance = relative_tolerance * max(section.z_span, 1.0)
    height = high
    iterations = 0
    while iterations < max_iterations:
        iterations += 1
        height = 0.5 * (low + high)
        volume = section.integrals(height)[0]
        if volume < target:
            low = height
        else:
            high = height
        if high - low <= tolerance:
            break
    return 0.5 * (low + high), iterations


def heeled_displacement(
    mesh: Mesh,
    *,
    heel_deg: float,
    displacement_volume: float | None = None,
    waterline_z: float = 0.0,
    closed: bool | None = None,
    relative_tolerance: float = _DEFAULT_RELATIVE_TOLERANCE,
    max_iterations: int = _MAX_ITERATIONS,
) -> HeeledDisplacement:
    """Displaced volume and buoyancy centre with the hull held at ``heel_deg``.

    ``displacement_volume`` defaults to the **upright** displacement at ``waterline_z``
    (the loading the upright hull carries), which is then held constant as the body
    heels -- that is what a ship does, and it is the reason a heeled waterline is not
    simply the upright one. ``waterline_z`` only fixes that upright loading.
    """
    if not math.isfinite(heel_deg):
        raise RightingArmError(f"heel_deg must be finite, got {heel_deg!r}")
    if abs(heel_deg) > MAX_HEEL_DEG:
        raise RightingArmError(
            f"heel_deg {heel_deg!r} is outside (-{MAX_HEEL_DEG:g}, {MAX_HEEL_DEG:g}]; "
            "past that the transverse heel convention no longer describes the problem"
        )
    if not math.isfinite(waterline_z):
        raise RightingArmError(f"waterline_z must be finite, got {waterline_z!r}")
    if not math.isfinite(relative_tolerance) or relative_tolerance <= 0.0:
        raise RightingArmError(
            f"relative_tolerance must be positive and finite, got {relative_tolerance!r}"
        )
    _require_closed(mesh, closed)
    target = (
        _upright_volume(mesh, waterline_z)
        if displacement_volume is None
        else float(displacement_volume)
    )
    if not math.isfinite(target) or target <= 0.0:
        raise RightingArmError(
            f"displacement_volume must be positive and finite, got {target!r} m^3"
        )
    section = _HeeledSection(mesh, heel_deg, waterline_z)
    height, iterations = _solve_waterline(
        section, target, max_iterations=max_iterations, relative_tolerance=relative_tolerance
    )
    volume, moment_x, moment_y, moment_z = section.integrals(height)
    if volume <= 0.0 or not math.isfinite(volume):
        raise RightingArmError(
            f"the heeled hull displaces {volume!r} m^3 at waterline {height!r}: "
            "the solve did not bracket a floating attitude"
        )
    centre = (
        moment_x / volume + section.shift[0],
        moment_y / volume + section.shift[1],
        moment_z / volume + section.shift[2] + waterline_z,
    )
    return HeeledDisplacement(
        heel_deg=heel_deg,
        waterline_z=height + waterline_z,
        volume=volume,
        centre_of_buoyancy=centre,
        iterations=iterations,
        volume_residual=(volume - target) / target,
    )


# ---------------------------------------------------------------------------
# the righting arm
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class RightingArm:
    """One point of a righting-arm curve."""

    heel_deg: float
    righting_arm: float
    kg: float
    lcg: float
    keel_z: float
    displacement_volume: float
    waterline_z: float
    centre_of_buoyancy: Vec3
    centre_of_gravity: Vec3
    iterations: int
    volume_residual: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA_RIGHTING_ARM,
            "heel_deg": _round(self.heel_deg),
            "righting_arm": _round(self.righting_arm, 12),
            "kg": _round(self.kg),
            "lcg": _round(self.lcg),
            "keel_z": _round(self.keel_z),
            "displacement_volume": _round(self.displacement_volume, 12),
            "waterline_z": _round(self.waterline_z, 12),
            "centre_of_buoyancy": [_round(value, 12) for value in self.centre_of_buoyancy],
            "centre_of_gravity": [_round(value, 12) for value in self.centre_of_gravity],
            "iterations": self.iterations,
            "volume_residual": self.volume_residual,
        }


def righting_arm(
    mesh: Mesh,
    *,
    heel_deg: float,
    kg: float,
    lcg: float = 0.0,
    waterline_z: float = 0.0,
    displacement_volume: float | None = None,
    closed: bool | None = None,
    relative_tolerance: float = _DEFAULT_RELATIVE_TOLERANCE,
) -> RightingArm:
    """Transverse righting arm ``GZ`` at one heel for one loading condition.

    ``kg`` is the height of the centre of gravity **above the keel** (the mesh's lowest
    vertex) and ``lcg`` its longitudinal position; both are loading inputs. The arm is
    ``GZ = y_B - y_G`` with both points in the world frame of the heeled body, so
    ``GZ(0) = 0`` for a hull whose centre of gravity is on the centreline.
    """
    if not math.isfinite(kg):
        raise RightingArmError(f"kg must be finite, got {kg!r}")
    if not math.isfinite(lcg):
        raise RightingArmError(f"lcg must be finite, got {lcg!r}")
    state = heeled_displacement(
        mesh,
        heel_deg=heel_deg,
        displacement_volume=displacement_volume,
        waterline_z=waterline_z,
        closed=closed,
        relative_tolerance=relative_tolerance,
    )
    phi = math.radians(heel_deg)
    keel_z = min(vertex[2] for vertex in mesh.vertices)
    # Both G and B are expressed in the pivot frame (z measured from the upright
    # waterline): the arm is invariant under that choice, the two points are not.
    gravity_body = (lcg, 0.0, keel_z + kg - waterline_z)
    gravity = (
        gravity_body[0],
        gravity_body[1] * math.cos(phi) + gravity_body[2] * math.sin(phi),
        -gravity_body[1] * math.sin(phi) + gravity_body[2] * math.cos(phi),
    )
    arm = state.centre_of_buoyancy[1] - gravity[1]
    return RightingArm(
        heel_deg=heel_deg,
        righting_arm=arm,
        kg=kg,
        lcg=lcg,
        keel_z=keel_z,
        displacement_volume=state.volume,
        waterline_z=state.waterline_z,
        centre_of_buoyancy=state.centre_of_buoyancy,
        centre_of_gravity=gravity,
        iterations=state.iterations,
        volume_residual=state.volume_residual,
    )


@dataclass(frozen=True)
class RightingArmCurve:
    """``GZ`` sampled over a heel ladder, with the curve quantities a rule asks for."""

    points: tuple[RightingArm, ...]
    kg: float
    lcg: float
    displacement_volume: float
    basis: str = "constant displacement, waterline solved at every heel"

    def gz_at(self, heel_deg: float) -> float:
        """Linearly interpolated arm at ``heel_deg`` (linear in the sampled ladder)."""
        heels = [point.heel_deg for point in self.points]
        arms = [point.righting_arm for point in self.points]
        if heel_deg < heels[0] or heel_deg > heels[-1]:
            raise RightingArmError(
                f"heel {heel_deg!r} deg is outside the sampled curve "
                f"[{heels[0]!r}, {heels[-1]!r}]; sample it or use righting_arm()"
            )
        for index in range(len(heels) - 1):
            if heels[index] <= heel_deg <= heels[index + 1]:
                span = heels[index + 1] - heels[index]
                if span == 0.0:
                    return arms[index]
                weight = (heel_deg - heels[index]) / span
                return arms[index] + weight * (arms[index + 1] - arms[index])
        return arms[-1]

    def peak(self) -> tuple[float, float]:
        """``(heel_deg, GZ)`` at the maximum of the sampled curve."""
        best = max(self.points, key=lambda point: point.righting_arm)
        return best.heel_deg, best.righting_arm

    def positive_range_deg(self) -> float:
        """Heel at which the interpolated curve first returns to zero, from upright.

        A conventional *range of stability* needs the curve sampled past its vanishing
        point; when the ladder stops while ``GZ`` is still positive this returns the last
        sampled heel and the caller must not read it as a vanishing angle -- that is why
        the value is reported together with :meth:`vanishing_crossed`.
        """
        return self._vanishing()[0]

    def vanishing_crossed(self) -> bool:
        """True when the sampled ladder contains a change of sign of ``GZ``."""
        return self._vanishing()[1]

    def _vanishing(self) -> tuple[float, bool]:
        # Arms below this are the upright point's round-off, not a vanishing angle.
        zero = 1e-12
        previous = self.points[0]
        if previous.righting_arm < -zero:
            return previous.heel_deg, True
        for point in self.points[1:]:
            if point.righting_arm <= zero:
                span = point.heel_deg - previous.heel_deg
                if span == 0.0 or previous.righting_arm == point.righting_arm:
                    return previous.heel_deg, True
                weight = previous.righting_arm / (previous.righting_arm - point.righting_arm)
                return previous.heel_deg + weight * span, True
            previous = point
        return self.points[-1].heel_deg, False

    def area_under_curve(self) -> float:
        """Dynamical stability ``int GZ dphi`` in ``m rad`` over the sampled ladder."""
        total = 0.0
        for previous, point in pairwise(self.points):
            span = math.radians(point.heel_deg - previous.heel_deg)
            total += 0.5 * (previous.righting_arm + point.righting_arm) * span
        return total

    def to_dict(self) -> dict[str, Any]:
        crossed = self._vanishing()[1]
        peak_heel, peak_arm = self.peak()
        return {
            "schema": SCHEMA_CURVE,
            "kg": _round(self.kg),
            "lcg": _round(self.lcg),
            "displacement_volume": _round(self.displacement_volume, 12),
            "heels_deg": [_round(point.heel_deg) for point in self.points],
            "righting_arms": [_round(point.righting_arm, 12) for point in self.points],
            "peak": {"heel_deg": _round(peak_heel), "righting_arm": _round(peak_arm, 12)},
            "positive_range_deg": _round(self.positive_range_deg()),
            "vanishing_crossed": crossed,
            "area_under_curve_m_rad": _round(self.area_under_curve(), 12),
            "basis": self.basis,
            "points": [point.to_dict() for point in self.points],
        }


def righting_arm_curve(
    mesh: Mesh,
    *,
    kg: float,
    heels_deg: Sequence[float] = DEFAULT_HEELS_DEG,
    lcg: float = 0.0,
    waterline_z: float = 0.0,
    displacement_volume: float | None = None,
    closed: bool | None = None,
    relative_tolerance: float = _DEFAULT_RELATIVE_TOLERANCE,
) -> RightingArmCurve:
    """``GZ`` over a heel ladder at one constant displacement."""
    heels = tuple(float(heel) for heel in heels_deg)
    if len(heels) < 2:
        raise RightingArmError(
            f"a curve needs at least two heels, got {len(heels)}"
        )
    if any(later <= earlier for earlier, later in pairwise(heels)):
        raise RightingArmError(f"heels_deg must be strictly increasing, got {heels!r}")
    points = tuple(
        righting_arm(
            mesh,
            heel_deg=heel,
            kg=kg,
            lcg=lcg,
            waterline_z=waterline_z,
            displacement_volume=displacement_volume,
            closed=closed,
            relative_tolerance=relative_tolerance,
        )
        for heel in heels
    )
    return RightingArmCurve(
        points=points,
        kg=kg,
        lcg=lcg,
        displacement_volume=points[0].displacement_volume,
    )


# ---------------------------------------------------------------------------
# the closed form a box pins the sign and magnitude with
# ---------------------------------------------------------------------------
def wall_sided_righting_arm(
    *, breadth: float, draft: float, kg: float, heel_deg: float
) -> float:
    """Exact ``GZ`` of a wall-sided (box-shaped) barge at constant displacement.

    ``GZ = sin(phi) (GM + (1/2) BM tan^2(phi))`` with ``KB = T/2``,
    ``BM = B^2/(12 T)`` and ``GM = KB + BM - KG``. Exact while the waterline stays on
    the vertical sides; :func:`wall_sided_validity_deg` says where that stops.
    """
    if not (math.isfinite(breadth) and math.isfinite(draft) and math.isfinite(kg)):
        raise RightingArmError("breadth, draft and kg must all be finite")
    if breadth <= 0.0 or draft <= 0.0:
        raise RightingArmError(
            f"breadth {breadth!r} and draft {draft!r} must both be positive"
        )
    if not math.isfinite(heel_deg) or abs(heel_deg) > MAX_HEEL_DEG:
        raise RightingArmError(f"heel_deg {heel_deg!r} is outside (-90, 90]")
    phi = math.radians(heel_deg)
    kb = 0.5 * draft
    bm = breadth * breadth / (12.0 * draft)
    gm = kb + bm - kg
    tangent = math.tan(phi)
    return math.sin(phi) * (gm + 0.5 * bm * tangent * tangent)


def wall_sided_validity_deg(
    *, breadth: float, draft: float, depth: float | None = None
) -> float:
    """Heel (degrees) beyond which :func:`wall_sided_righting_arm` is no longer exact.

    The waterline must stay on the vertical sides: the bilge is uncovered when
    ``tan(phi) = 2 T / B`` and the deck edge immerses when ``tan(phi) = 2 (D - T) / B``.
    The smaller of the two governs; without ``depth`` only the bilge limit is known.
    """
    if not (math.isfinite(breadth) and math.isfinite(draft)):
        raise RightingArmError("breadth and draft must be finite")
    if breadth <= 0.0 or draft <= 0.0:
        raise RightingArmError(
            f"breadth {breadth!r} and draft {draft!r} must both be positive"
        )
    limits = [math.degrees(math.atan2(2.0 * draft, breadth))]
    if depth is not None:
        if not math.isfinite(depth) or depth <= draft:
            raise RightingArmError(
                f"depth {depth!r} must be finite and above the draft {draft!r}"
            )
        limits.append(math.degrees(math.atan2(2.0 * (depth - draft), breadth)))
    return min(limits)
