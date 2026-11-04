"""Large-angle stability: righting arms integrated from the submerged volume.

:mod:`tensorcad.stability` answers the small-angle question (``KM``, ``GM``) from the
upright waterplane. That is enough to say whether a loading condition is stable and
nothing more: a ``GM`` cannot produce a ``GZ`` curve, and the deck-edge immersion and
the loss of waterplane that decide whether a hull survives 30-60 degrees of heel are
second-order in the small-angle expansion, i.e. absent from it. This module computes
the righting arm directly:

    ``GZ(phi) = KN(phi) - KG sin(phi)``

with ``KN`` the transverse lever of the centre of buoyancy about the keel, at
*constant displacement* and level trim. For every heel angle the immersed region is
therefore the hull clipped by a **tilted** waterplane whose sinkage is solved so that
the displaced volume is unchanged -- exactly the geometry
:mod:`tensorcad.hydrostatics` declines to handle ("no face straddles the waterline").

Exactness
---------
The immersed region is not sampled on a grid. A triangle mesh cut by a plane of
constant ``z`` produces a cross-section whose boundary vertices move **linearly** in
``z``, so its area and first moments are **quadratic in ``z``** between two
consecutive vertex heights. Three samples per band (at ``1/4, 1/2, 3/4`` of the band,
chosen so no sample ever lands on a vertex height) determine that quadratic exactly,
so

* volume, first moments and waterplane integrals below *any* plane are exact to
  floating point rather than to a discretisation tolerance;
* a heel angle needs no march in ``phi``, and the sinkage solve is a bisection on a
  monotone function instead of a re-mesh.

Two consequences carry the verification: a box-shaped hull reproduces the textbook
wall-sided righting-arm formula ``sin(phi)(GM + BM tan^2(phi)/2)`` to ~1e-15, and
``dGZ/dphi`` at zero equals the ``GM`` of
:func:`tensorcad.stability.initial_stability` -- references that share no code with
this integrator.

Conventions
-----------
* Body frame as in the rest of the package: ``+x`` forward, ``+y`` to starboard,
  ``+z`` up, upright waterline at ``z = 0``.
* ``heel_deg > 0`` heels **to starboard** (the starboard side goes down): the body is
  rotated about the ``x`` axis by ``-heel`` so the water stays level.
* ``kg`` is the height of the centre of gravity **above the keel** -- a loading input,
  never a shape property, therefore a required keyword.
* Fixed trim: the hull is held at level trim at every heel angle. A trim excursion
  changes ``GZ`` by an amount this module does not model.

What this is not
----------------
No free-surface correction (tank contents are solid here), no wind or launching
heels, no damage or flooded-compartment cases, no IMO/IS-Code compliance statement.
``GZ`` is a geometric function of one hull, one displacement and one ``KG``; whether
a ship satisfies a criterion is a regulatory question this package does not answer.
"""
from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, replace
from itertools import pairwise
from typing import Any

from . import mesh as meshlib
from .geometry import Mesh, MeshBuilder, Vec3

__all__ = [
    "GzCurve",
    "HeeledHydrostatics",
    "LargeAngleError",
    "RightingArm",
    "box_dynamical_stability",
    "box_hull_mesh",
    "box_righting_arm",
    "dynamical_stability",
    "extend_freeboard",
    "gz_curve",
    "heeled_hydrostatics",
    "righting_arm",
    "wall_sided_limit_deg",
]

#: Sea water density used when a caller does not override it (kg/m^3).
WATER_DENSITY = 1025.0


class LargeAngleError(ValueError):
    """Raised when a large-angle request cannot be answered honestly."""


def _round(value: float, digits: int = 9) -> float:
    if not math.isfinite(value):
        return value
    return round(value, digits)


# --------------------------------------------------------------------------- #
# public value types
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class HeeledHydrostatics:
    """Properties of the immersed volume at one heel angle, level trim."""

    heel_deg: float
    heel_rad: float
    volume: float
    waterline_z: float
    centre_of_buoyancy: Vec3
    displacement_mass: float
    water_density: float
    target_volume: float
    volume_error: float
    keel_z: float
    basis: str = "exact section integrals, tilted waterplane, level trim"

    @property
    def kn(self) -> float:
        """Transverse lever of the centre of buoyancy about the keel (m).

        Measured from the upright keel line, the convention ``KN`` is tabulated in:
        ``KN(0) = 0`` and ``GZ = KN - KG sin(phi)``.
        """
        return self.centre_of_buoyancy[1] - self.keel_z * math.sin(self.heel_rad)

    def to_dict(self) -> dict[str, Any]:
        return {
            "heel_deg": _round(self.heel_deg),
            "heel_rad": _round(self.heel_rad),
            "volume_m3": _round(self.volume),
            "target_volume_m3": _round(self.target_volume),
            "relative_volume_error": _round(self.volume_error, 12),
            "waterline_z_m": _round(self.waterline_z),
            "centre_of_buoyancy_m": tuple(_round(value) for value in self.centre_of_buoyancy),
            "kn_m": _round(self.kn),
            "displacement_mass_kg": _round(self.displacement_mass, 3),
            "water_density_kg_per_m3": _round(self.water_density, 3),
            "keel_z_m": _round(self.keel_z),
            "basis": self.basis,
        }


@dataclass(frozen=True)
class RightingArm:
    """One point of a ``GZ`` curve: the lever for one heel angle and one ``KG``."""

    heel_deg: float
    gz: float
    kn: float
    kg: float
    hydrostatics: HeeledHydrostatics

    @property
    def heel_rad(self) -> float:
        return self.hydrostatics.heel_rad

    @property
    def volume(self) -> float:
        return self.hydrostatics.volume

    def to_dict(self) -> dict[str, Any]:
        return {
            "heel_deg": _round(self.heel_deg),
            "gz_m": _round(self.gz),
            "kn_m": _round(self.kn),
            "kg_m": _round(self.kg),
            "volume_m3": _round(self.volume),
            "waterline_z_m": _round(self.hydrostatics.waterline_z),
            "centre_of_buoyancy_m": tuple(
                _round(value) for value in self.hydrostatics.centre_of_buoyancy
            ),
        }


@dataclass(frozen=True)
class GzCurve:
    """A ``GZ`` curve, with the loading condition it was computed for."""

    points: tuple[RightingArm, ...]
    kg: float
    displacement_volume: float
    water_density: float

    def __len__(self) -> int:
        return len(self.points)

    @property
    def angles_deg(self) -> tuple[float, ...]:
        return tuple(point.heel_deg for point in self.points)

    @property
    def arms(self) -> tuple[float, ...]:
        return tuple(point.gz for point in self.points)

    @property
    def peak(self) -> RightingArm:
        """Point of maximum righting arm (the curve's stability peak)."""
        if not self.points:
            raise LargeAngleError("an empty GZ curve has no peak")
        return max(self.points, key=lambda point: point.gz)

    @property
    def vanishing_angle_deg(self) -> float | None:
        """Angle at which ``GZ`` first returns to zero, linearly interpolated.

        ``None`` means the computed sector did not bracket the vanishing angle (or
        the curve never became positive); the caller is told so instead of being
        handed the last computed angle as if it were the answer.
        """
        previous: RightingArm | None = None
        for point in self.points:
            if point.heel_deg <= 0.0:
                previous = point
                continue
            if point.gz <= 0.0 and previous is not None and previous.gz > 0.0:
                span = previous.gz - point.gz
                if span <= 0.0:
                    return point.heel_deg
                return previous.heel_deg + (previous.gz / span) * (point.heel_deg - previous.heel_deg)
            previous = point
        return None

    def dynamical_stability(self, upto_deg: float | None = None) -> float:
        """Area under the curve in ``m*rad`` (trapezoid rule on the points given)."""
        return dynamical_stability(self.points, upto_deg=upto_deg)

    def at(self, heel_deg: float) -> RightingArm:
        """Exact curve point at ``heel_deg``; the angle must be one of the computed ones."""
        for point in self.points:
            if abs(point.heel_deg - heel_deg) <= 1e-9:
                return point
        raise LargeAngleError(f"heel_deg {heel_deg!r} is not a computed angle of this curve")

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": "tensorcad.stability.gz_curve/1",
            "kg_m": _round(self.kg),
            "displacement_volume_m3": _round(self.displacement_volume),
            "water_density_kg_per_m3": _round(self.water_density, 3),
            "points": [point.to_dict() for point in self.points],
            "peak_gz_m": _round(self.peak.gz),
            "peak_angle_deg": _round(self.peak.heel_deg),
            "vanishing_angle_deg": (
                None if self.vanishing_angle_deg is None else _round(self.vanishing_angle_deg)
            ),
            "dynamical_stability_m_rad_upto_40deg": _round(self.dynamical_stability(40.0)),
        }


# --------------------------------------------------------------------------- #
# exact section integrals
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class _Band:
    """One ``z`` band between consecutive vertex heights, with quadratic fits."""

    z0: float
    dz: float
    area: tuple[float, float, float]
    moment_x: tuple[float, float, float]
    moment_y: tuple[float, float, float]


@dataclass(frozen=True)
class _Ladder:
    bands: tuple[_Band, ...]
    z_min: float
    z_max: float


def _fit_quadratic(f1: float, f2: float, f3: float) -> tuple[float, float, float]:
    """Coefficients of the quadratic with values ``f1, f2, f3`` at ``u = 1/4, 1/2, 3/4``.

    Those three nodes integrate a quadratic (indeed a cubic, Milne's rule) exactly and
    never coincide with a band edge, where a cross-section is degenerate.
    """
    c2 = 8.0 * (f1 + f3 - 2.0 * f2)
    c1 = 2.0 * (f3 - f1) - c2
    c0 = f2 - 0.5 * c1 - 0.25 * c2
    return (c0, c1, c2)


def _integral(poly: tuple[float, float, float], s: float = 1.0) -> float:
    """``int_0^s poly(u) du``."""
    c0, c1, c2 = poly
    return c0 * s + 0.5 * c1 * s * s + (c2 / 3.0) * s * s * s


def _integral_times_u(poly: tuple[float, float, float], s: float = 1.0) -> float:
    """``int_0^s u * poly(u) du``."""
    c0, c1, c2 = poly
    return 0.5 * c0 * s * s + (c1 / 3.0) * s ** 3 + 0.25 * c2 * s ** 4


def _face_normal_xy(vertices: Sequence[Vec3], face: tuple[int, int, int]) -> tuple[float, float]:
    """``z_hat x n`` for the face's outward normal ``n``, projected into the xy plane.

    That vector is tangent to the cross-section boundary and points the way a
    counter-clockwise (positive-area) loop traverses it, so orienting each cut segment
    along it makes the shoelace sums below come out positive for an outward-oriented
    surface -- without linking loops, which is what lets a section split into several
    closed curves (or develop a hole) without special cases.
    """
    a, b, c = (vertices[index] for index in face)
    ux, uy, uz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
    vx, vy, vz = c[0] - a[0], c[1] - a[1], c[2] - a[2]
    nx = uy * vz - uz * vy
    ny = uz * vx - ux * vz
    return (-ny, nx)


def _slice_segment(
    face: tuple[int, int, int],
    vertices: Sequence[Vec3],
    z: float,
    normal_xy: tuple[float, float],
) -> tuple[tuple[float, float], tuple[float, float]] | None:
    """Directed cut segment of one triangle at height ``z`` (``None`` if it misses).

    The direction is ``z_hat x n`` for the face's outward normal ``n``, the way a
    positive-area loop traverses the boundary, so the shoelace sums below come out
    positive without ever linking loops. A section that splits into several loops (or
    develops a hole) is handled by the same sign rule.
    """
    below = [vertices[index][2] < z for index in face]
    if all(below) or not any(below):
        return None
    crossings: list[tuple[float, float]] = []
    for i, j in ((0, 1), (1, 2), (2, 0)):
        if below[i] == below[j]:
            continue
        first = vertices[face[i]]
        second = vertices[face[j]]
        t = (z - first[2]) / (second[2] - first[2])
        crossings.append(
            (first[0] + t * (second[0] - first[0]), first[1] + t * (second[1] - first[1]))
        )
    if len(crossings) != 2:
        return None
    first_xy, second_xy = crossings
    direction = (second_xy[0] - first_xy[0]) * normal_xy[0] + (
        second_xy[1] - first_xy[1]
    ) * normal_xy[1]
    if direction < 0.0:
        return second_xy, first_xy
    return first_xy, second_xy


def _build_ladder(mesh: Mesh) -> _Ladder:
    """Precompute the exact cross-section ladder of a closed triangle mesh.

    Face records are swept in ``z`` so each of the three samples per band visits only
    the triangles that actually straddle it: the cost is proportional to the faces
    per band, not to the whole mesh times the number of bands.
    """
    vertices = mesh.vertices
    if not vertices:
        raise LargeAngleError("cannot integrate an empty mesh")
    levels = sorted({vertex[2] for vertex in vertices})
    if len(levels) < 2:
        raise LargeAngleError(
            "the mesh is flat in z; a waterplane integral needs a mesh with vertical extent"
        )

    faces = mesh.faces
    normals = [_face_normal_xy(vertices, face) for face in faces]
    z_low = [min(vertices[index][2] for index in face) for face in faces]
    z_high = [max(vertices[index][2] for index in face) for face in faces]

    band_count = len(levels) - 1
    schedule: list[tuple[float, int, int]] = []
    for band_index in range(band_count):
        z0 = levels[band_index]
        dz = levels[band_index + 1] - z0
        if dz <= 0.0:
            continue
        for slot, fraction in enumerate((0.25, 0.5, 0.75)):
            schedule.append((z0 + fraction * dz, band_index, slot))
    if not schedule:
        raise LargeAngleError("no vertical band with positive height; nothing to integrate")
    schedule.sort()

    samples: dict[int, list[list[float]]] = {
        band_index: [[0.0, 0.0, 0.0] for _ in range(3)] for band_index in range(band_count)
    }
    order = sorted(range(len(faces)), key=lambda index: z_low[index])
    cursor = 0
    active: set[int] = set()
    for z, band_index, slot in schedule:
        while cursor < len(order) and z_low[order[cursor]] <= z:
            active.add(order[cursor])
            cursor += 1
        for face_index in [index for index in active if z_high[index] < z]:
            active.discard(face_index)
        area = moment_x = moment_y = 0.0
        for face_index in active:
            segment = _slice_segment(faces[face_index], vertices, z, normals[face_index])
            if segment is None:
                continue
            (px, py), (qx, qy) = segment
            cross = px * qy - qx * py
            area += 0.5 * cross
            moment_x += cross * (px + qx) / 6.0
            moment_y += cross * (py + qy) / 6.0
        samples[band_index][slot] = [area, moment_x, moment_y]

    bands: list[_Band] = []
    for band_index in range(band_count):
        z0 = levels[band_index]
        dz = levels[band_index + 1] - z0
        if dz <= 0.0:
            continue
        values = samples[band_index]
        bands.append(
            _Band(
                z0=z0,
                dz=dz,
                area=_fit_quadratic(values[0][0], values[1][0], values[2][0]),
                moment_x=_fit_quadratic(values[0][1], values[1][1], values[2][1]),
                moment_y=_fit_quadratic(values[0][2], values[1][2], values[2][2]),
            )
        )
    return _Ladder(bands=tuple(bands), z_min=levels[0], z_max=levels[-1])


def _integrals_below(ladder: _Ladder, waterline_z: float) -> tuple[float, float, float, float]:
    """``(volume, int x dV, int y dV, int z dV)`` of the solid below ``waterline_z``."""
    volume = first_x = first_y = first_z = 0.0
    for band in ladder.bands:
        if band.z0 >= waterline_z:
            break
        s = 1.0
        if band.z0 + band.dz > waterline_z:
            s = (waterline_z - band.z0) / band.dz
            if s <= 0.0:
                break
        area_integral = _integral(band.area, s)
        volume += band.dz * area_integral
        first_x += band.dz * _integral(band.moment_x, s)
        first_y += band.dz * _integral(band.moment_y, s)
        first_z += band.dz * (band.z0 * area_integral + band.dz * _integral_times_u(band.area, s))
    return volume, first_x, first_y, first_z


def _solve_waterline(ladder: _Ladder, target_volume: float) -> tuple[float, float]:
    """Bisect the waterplane height that displaces ``target_volume``."""
    if target_volume <= 0.0:
        raise LargeAngleError(f"target displacement volume must be positive, got {target_volume!r}")
    low, high = ladder.z_min, ladder.z_max
    span = high - low
    tolerance = 4.0 * math.ulp(max(abs(low), abs(high), span)) + 1e-14 * max(span, 1.0)
    volume_low, *_ = _integrals_below(ladder, low)
    volume_high, *_ = _integrals_below(ladder, high)
    if not (volume_low <= target_volume <= volume_high):
        raise LargeAngleError(
            f"target volume {target_volume!r} m^3 is outside the mesh "
            f"({volume_low!r} .. {volume_high!r})"
        )
    for _ in range(200):
        middle = 0.5 * (low + high)
        volume_middle, *_ = _integrals_below(ladder, middle)
        if volume_middle < target_volume:
            low = middle
        else:
            high = middle
        if high - low <= tolerance:
            break
    waterline_z = 0.5 * (low + high)
    volume, *_ = _integrals_below(ladder, waterline_z)
    return waterline_z, abs(volume - target_volume) / target_volume


# --------------------------------------------------------------------------- #
# rotation and hydrostatics at a heel angle
# --------------------------------------------------------------------------- #
def _finite(name: str, value: float) -> float:
    try:
        numeric = float(value)
    except (TypeError, ValueError) as error:
        raise LargeAngleError(f"{name} must be a number, got {value!r}") from error
    if not math.isfinite(numeric):
        raise LargeAngleError(f"{name} must be finite, got {value!r}")
    return numeric


def _positive(name: str, value: float) -> float:
    numeric = _finite(name, value)
    if numeric <= 0.0:
        raise LargeAngleError(f"{name} must be positive, got {value!r}")
    return numeric


def _non_negative(name: str, value: float) -> float:
    numeric = _finite(name, value)
    if numeric < 0.0:
        raise LargeAngleError(f"{name} must not be negative, got {value!r}")
    return numeric


def _heeled_mesh(mesh: Mesh, heel_rad: float) -> Mesh:
    """Rotate the body by ``-heel`` about the x axis so the water stays level."""
    cosine = math.cos(heel_rad)
    sine = math.sin(heel_rad)
    vertices = tuple(
        (x, y * cosine + z * sine, -y * sine + z * cosine) for (x, y, z) in mesh.vertices
    )
    return replace(mesh, vertices=vertices)


def _checked_ladder(mesh: Mesh) -> _Ladder:
    if not meshlib.is_watertight(mesh):
        raise LargeAngleError(
            "large-angle stability needs a closed, outward-oriented surface: an open mesh has no "
            "displacement to conserve, so GZ would be meaningless"
        )
    if meshlib.total_signed_volume(mesh) <= 0.0:
        raise LargeAngleError(
            "the mesh encloses no positive volume; large-angle GZ needs outward-oriented faces"
        )
    return _build_ladder(mesh)


def heeled_hydrostatics(
    mesh: Mesh,
    heel_deg: float,
    *,
    displacement_volume: float | None = None,
    waterline_z: float | None = None,
    water_density: float = WATER_DENSITY,
    keel_z: float | None = None,
) -> HeeledHydrostatics:
    """Immersed volume and centre of buoyancy at one heel angle, level trim.

    Either ``displacement_volume`` (the usual case: hold the upright displacement and
    let the sinkage follow) or ``waterline_z`` (hold the waterplane and let the
    displacement follow, which is what a fixed-carriage model test does) may be given
    -- not both. With neither, the upright displacement at ``z = 0`` is the target:
    the hull as it floats at rest.
    """
    heel_rad = math.radians(_finite("heel_deg", heel_deg))
    if abs(math.degrees(heel_rad)) > 180.0:
        raise LargeAngleError(f"heel_deg must be within +-180, got {heel_deg!r}")
    density = _positive("water_density", water_density)
    if displacement_volume is not None and waterline_z is not None:
        raise LargeAngleError(
            "give either displacement_volume or waterline_z, not both: the two define each other"
        )
    upright_ladder = _checked_ladder(mesh)
    keel = min(vertex[2] for vertex in mesh.vertices) if keel_z is None else _finite("keel_z", keel_z)

    if waterline_z is not None:
        held = _finite("waterline_z", waterline_z)
        ladder = _build_ladder(_heeled_mesh(mesh, heel_rad))
        volume, first_x, first_y, first_z = _integrals_below(ladder, held)
        if volume <= 0.0:
            raise LargeAngleError(
                f"waterline_z {held!r} leaves no immersed volume at heel {heel_deg!r} deg"
            )
        return HeeledHydrostatics(
            heel_deg=heel_deg,
            heel_rad=heel_rad,
            volume=volume,
            waterline_z=held,
            centre_of_buoyancy=(first_x / volume, first_y / volume, first_z / volume),
            displacement_mass=volume * density,
            water_density=density,
            target_volume=volume,
            volume_error=0.0,
            keel_z=keel,
        )

    if displacement_volume is None:
        target, *_ = _integrals_below(upright_ladder, 0.0)
        if target <= 0.0:
            raise LargeAngleError(
                "the upright waterline z=0 is at or below the keel, so the hull displaces nothing; "
                "pass displacement_volume or waterline_z explicitly"
            )
    else:
        target = _positive("displacement_volume", displacement_volume)

    whole = meshlib.total_signed_volume(mesh)
    if abs(heel_rad) > 0.0 and abs(target - whole) <= 1e-9 * whole:
        raise LargeAngleError(
            f"the mesh is capped at its own top surface: it displaces its whole volume "
            f"({whole!r} m^3), so no tilted waterplane can hold that displacement and a heeled "
            "righting arm is undefined. This is a hull *without freeboard* -- wrap the mesh in "
            "extend_freeboard() (or pass a displacement below the full volume) first"
        )

    ladder = _build_ladder(_heeled_mesh(mesh, heel_rad))
    solved_waterline, relative_error = _solve_waterline(ladder, target)
    volume, first_x, first_y, first_z = _integrals_below(ladder, solved_waterline)
    return HeeledHydrostatics(
        heel_deg=heel_deg,
        heel_rad=heel_rad,
        volume=volume,
        waterline_z=solved_waterline,
        centre_of_buoyancy=(first_x / volume, first_y / volume, first_z / volume),
        displacement_mass=volume * density,
        water_density=density,
        target_volume=target,
        volume_error=relative_error,
        keel_z=keel,
    )


def righting_arm(
    mesh: Mesh,
    heel_deg: float,
    *,
    kg: float,
    displacement_volume: float | None = None,
    waterline_z: float | None = None,
    water_density: float = WATER_DENSITY,
    keel_z: float | None = None,
) -> RightingArm:
    """Righting arm ``GZ = KN - KG sin(phi)`` at one heel angle.

    ``kg`` is measured **above the keel** and is required: a default would make a
    silent ``KG = 0`` look like a valid loading condition.
    """
    kg_value = _finite("kg", kg)
    hydrostatics = heeled_hydrostatics(
        mesh,
        heel_deg,
        displacement_volume=displacement_volume,
        waterline_z=waterline_z,
        water_density=water_density,
        keel_z=keel_z,
    )
    return RightingArm(
        heel_deg=heel_deg,
        gz=hydrostatics.kn - kg_value * math.sin(hydrostatics.heel_rad),
        kn=hydrostatics.kn,
        kg=kg_value,
        hydrostatics=hydrostatics,
    )


def gz_curve(
    mesh: Mesh,
    angles_deg: Iterable[float] | None = None,
    *,
    kg: float,
    displacement_volume: float | None = None,
    waterline_z: float | None = None,
    water_density: float = WATER_DENSITY,
    up_to_deg: float = 60.0,
    step_deg: float = 5.0,
) -> GzCurve:
    """``GZ`` curve from ``0`` to ``up_to_deg``, or at the ``angles_deg`` given.

    Every point holds the same displacement, so the sinkage differs from point to
    point: that is the physical content of a righting-arm curve and the reason it
    cannot be read off a single waterplane.
    """
    kg_value = _finite("kg", kg)
    if angles_deg is None:
        step = _positive("step_deg", step_deg)
        limit = _positive("up_to_deg", up_to_deg)
        count = math.floor(limit / step + 1e-9)
        angles = [index * step for index in range(count + 1)]
        if limit - angles[-1] > 1e-9:
            angles.append(limit)
    else:
        angles = [_finite("angles_deg entry", angle) for angle in angles_deg]
        if not angles:
            raise LargeAngleError("angles_deg must contain at least one angle")

    reference = heeled_hydrostatics(
        mesh,
        0.0,
        displacement_volume=displacement_volume,
        waterline_z=waterline_z,
        water_density=water_density,
    )
    points = tuple(
        righting_arm(
            mesh,
            angle,
            kg=kg_value,
            displacement_volume=reference.volume,
            water_density=reference.water_density,
            keel_z=reference.keel_z,
        )
        for angle in angles
    )
    return GzCurve(
        points=points,
        kg=kg_value,
        displacement_volume=reference.volume,
        water_density=reference.water_density,
    )


def dynamical_stability(points: Sequence[RightingArm], upto_deg: float | None = None) -> float:
    """Area under a ``GZ`` curve in ``m*rad``, trapezoid rule on the points given.

    The trapezoid rule is what the naval-architecture literature uses on 5-10 degree
    tabulations; ``tests/test_large_angle.py`` measures its ``O(h^2)`` convergence
    against the closed-form integral of the wall-sided formula instead of claiming
    exactness.
    """
    ordered = sorted(points, key=lambda point: point.heel_deg)
    if len(ordered) < 2:
        raise LargeAngleError("dynamical stability needs at least two curve points")
    limit = math.inf if upto_deg is None else _finite("upto_deg", upto_deg)
    area = 0.0
    for left, right in pairwise(ordered):
        if left.heel_deg >= limit:
            break
        span = right.heel_deg - left.heel_deg
        if span <= 0.0:
            continue
        fraction = min(1.0, (limit - left.heel_deg) / span)
        if fraction <= 0.0:
            continue
        arm_right = left.gz + fraction * (right.gz - left.gz)
        area += 0.5 * (left.gz + arm_right) * math.radians(span) * fraction
    return area


# --------------------------------------------------------------------------- #
# reference cases
# --------------------------------------------------------------------------- #
def box_hull_mesh(
    length: float,
    breadth: float,
    draft: float,
    freeboard: float = 0.0,
) -> Mesh:
    """Closed rectangular barge: sides from ``z = -draft`` up to ``z = +freeboard``.

    ``freeboard > 0`` matters. The wall-sided reference formula is only valid while
    the deck edge stays clear of the water, and a box capped exactly at the waterline
    (``freeboard = 0``) submerges a deck wedge from the first degree of heel -- a real
    geometric effect, not a numerical one, and the reason the tests sweep the
    freeboard.
    """
    half_length = 0.5 * _positive("length", length)
    half_breadth = 0.5 * _positive("breadth", breadth)
    bottom = -_positive("draft", draft)
    top = _non_negative("freeboard", freeboard)
    corners = [
        (-half_length, -half_breadth, bottom),
        (half_length, -half_breadth, bottom),
        (half_length, half_breadth, bottom),
        (-half_length, half_breadth, bottom),
        (-half_length, -half_breadth, top),
        (half_length, -half_breadth, top),
        (half_length, half_breadth, top),
        (-half_length, half_breadth, top),
    ]
    faces = [
        (0, 2, 1),
        (0, 3, 2),  # bottom, normal -z
        (4, 5, 6),
        (4, 6, 7),  # top, normal +z
        (0, 1, 5),
        (0, 5, 4),  # -y side
        (1, 2, 6),
        (1, 6, 5),  # +x end
        (2, 3, 7),
        (2, 7, 6),  # +y side
        (3, 0, 4),
        (3, 4, 7),  # -x end
    ]
    mesh = Mesh(vertices=tuple(corners), faces=tuple(faces))
    if meshlib.total_signed_volume(mesh) <= 0.0:  # pragma: no cover - sanity guard
        mesh = Mesh(
            vertices=tuple(corners),
            faces=tuple((a, c, b) for (a, b, c) in faces),
        )
    return mesh


def extend_freeboard(mesh: Mesh, height: float, *, deck_z: float | None = None) -> Mesh:
    """Raise a mesh's flat top surface by ``height``, adding vertical sides and a new deck.

    Every generator in :mod:`tensorcad.hull` and :mod:`tensorcad.wigley` caps the
    surface *at the design waterline*, which is right for a wetted-surface
    calculation and fatal for large-angle stability: such a body displaces its entire
    volume upright, so no tilted waterplane can hold that displacement and ``GZ`` does
    not exist (``heeled_hydrostatics`` says so instead of returning a number). This
    helper turns such a mesh into the hull-with-freeboard the calculation needs, by

    1. deleting the flat deck faces (all three vertices at ``deck_z``, default the
       highest ``z`` in the mesh),
    2. extruding the deck's boundary loop vertically by ``height``, and
    3. re-capping with the deleted deck's own triangulation, lifted.

    The extension is **waterline neutral**: no vertex below ``deck_z`` moves, so the
    immersed volume and centre of buoyancy reported by
    ``heeled_hydrostatics(extended, 0.0, waterline_z=0.0)`` reproduce
    :func:`tensorcad.hydrostatics.compute_hydrostatics` on the original capped mesh to
    floating point (asserted in ``tests/test_large_angle.py``). ``GM`` therefore does
    not change either, which is what makes this safe as a stability input rather than
    a new geometry model.

    One boundary to keep in mind: :func:`tensorcad.hydrostatics.compute_hydrostatics`
    and :func:`tensorcad.stability.initial_stability` integrate a surface that is
    *already clipped at the waterline* and treat a closed mesh as the whole wetted
    solid, so on a freeboarded mesh they report the entire solid instead of the
    displacement. Feed those two the capped mesh; feed this module the extended one.
    """
    lift = _positive("height", height)
    if not mesh.vertices:
        raise LargeAngleError("cannot extend an empty mesh")
    z_top = max(vertex[2] for vertex in mesh.vertices) if deck_z is None else _finite("deck_z", deck_z)
    scale = max(1.0, abs(z_top))
    tolerance = 1e-12 * scale
    vertices = mesh.vertices
    deck_indices: list[int] = []
    for index, face in enumerate(mesh.faces):
        if all(abs(vertices[vertex][2] - z_top) <= tolerance for vertex in face):
            deck_indices.append(index)
    if not deck_indices:
        raise LargeAngleError(
            f"no flat deck found at z = {z_top!r}; extend_freeboard needs a cap whose triangles "
            "all lie in one horizontal plane"
        )
    deck_set = set(deck_indices)
    deck_faces = [mesh.faces[index] for index in deck_indices]

    edge_counts: dict[tuple[int, int], int] = {}
    for a, b, c in deck_faces:
        for edge in ((a, b), (b, c), (c, a)):
            key = (min(edge), max(edge))
            edge_counts[key] = edge_counts.get(key, 0) + 1
    boundary = [edge for edge, count in edge_counts.items() if count == 1]
    if not boundary:
        raise LargeAngleError(
            "the deck triangles have no boundary loop; the cap is not an open surface to extrude"
        )

    builder = MeshBuilder(units=mesh.units, coordinate_frame=mesh.coordinate_frame)
    lifted: dict[int, int] = {}

    def lifted_index(index: int) -> int:
        cached = lifted.get(index)
        if cached is None:
            x, y, z = vertices[index]
            cached = builder.add_vertex((x, y, z + lift))
            lifted[index] = cached
        return cached

    for index, face in enumerate(mesh.faces):
        if index in deck_set:
            continue
        builder.add_triangle(*(builder.add_vertex(vertices[vertex]) for vertex in face))
    for a, b in boundary:
        if vertices[a][:2] == vertices[b][:2]:
            continue  # collapsed waterline end: no side wall to add
        a_lift, b_lift = lifted_index(a), lifted_index(b)
        a_index, b_index = builder.add_vertex(vertices[a]), builder.add_vertex(vertices[b])
        builder.add_triangle(a_index, b_index, b_lift)
        builder.add_triangle(a_index, b_lift, a_lift)
    for a, b, c in deck_faces:
        builder.add_triangle(lifted_index(a), lifted_index(b), lifted_index(c))

    extended = builder.build()
    oriented = meshlib.orient_outward(extended)
    result = oriented.mesh
    if not meshlib.is_watertight(result):
        raise LargeAngleError(
            "the freeboard extension did not produce a closed surface; the input mesh's deck was "
            "not a simple capped surface"
        )
    return result


def box_righting_arm(*, breadth: float, draft: float, kg: float, heel_deg: float) -> float:
    """Wall-sided reference ``GZ = sin(phi)(GM + BM tan^2(phi)/2)`` for a box barge.

    ``GM = KB + BM - KG`` with ``KB = T/2`` and ``BM = B^2/(12 T)``. The formula is
    exact for a rectangular section at constant displacement -- the immersed area of a
    heeled rectangle is a trapezoid whose mean depth is unchanged -- and it is the
    external reference :func:`righting_arm` is tested against. It stops being exact at
    the deck-edge immersion angle, see :func:`wall_sided_limit_deg`.
    """
    beam = _positive("breadth", breadth)
    draught = _positive("draft", draft)
    kg_value = _finite("kg", kg)
    heel_rad = math.radians(_finite("heel_deg", heel_deg))
    bm = beam * beam / (12.0 * draught)
    gm = 0.5 * draught + bm - kg_value
    tangent = math.tan(heel_rad)
    return math.sin(heel_rad) * (gm + 0.5 * bm * tangent * tangent)


def wall_sided_limit_deg(*, breadth: float, freeboard: float) -> float:
    """Heel angle at which the deck edge of a box barge reaches the water (degrees).

    Beyond it the wall-sided formula over-predicts the righting arm, because the
    immersed section loses its outer waterplane: the geometry, not the algebra, is
    what changes.
    """
    beam = _positive("breadth", breadth)
    height = _non_negative("freeboard", freeboard)
    return math.degrees(math.atan(2.0 * height / beam))


def box_dynamical_stability(
    *, breadth: float, draft: float, kg: float, upto_deg: float
) -> float:
    """Closed-form ``int_0^Phi GZ dphi`` for the wall-sided box, in ``m*rad``.

    ``GM (1 - cos Phi) + (BM/2) (cos Phi + sec Phi - 2)``. It integrates
    ``sin(phi)(GM + BM tan^2(phi)/2)`` exactly, so it grades the trapezoid rule in
    :func:`dynamical_stability` without sharing a single line of code with it.
    """
    beam = _positive("breadth", breadth)
    draught = _positive("draft", draft)
    kg_value = _finite("kg", kg)
    heel_rad = math.radians(_finite("upto_deg", upto_deg))
    bm = beam * beam / (12.0 * draught)
    gm = 0.5 * draught + bm - kg_value
    cosine = math.cos(heel_rad)
    return gm * (1.0 - cosine) + 0.5 * bm * (cosine + 1.0 / cosine - 2.0)
