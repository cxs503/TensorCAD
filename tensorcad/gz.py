"""Large-angle stability: the righting-arm curve ``GZ(phi)``.

:mod:`tensorcad.stability` ends its own docstring with a warning that it measures the
*slope* of the righting-arm curve at the upright and nothing else::

    GM is the slope of that curve at phi = 0, and stability at large angles is a
    different integration (a heeled displacement integration, not a waterplane moment).

This module is that integration. A hull heeled to starboard by ``phi`` floats on the
plane whose up-normal, written in the ship's own axes, is

    n = (0, -sin phi, cos phi)          (unit, exact)

so the submerged set is the half-space ``p . n <= d`` for a single scalar ``d``. For a
fixed loading the displacement does not change with heel, so ``d`` is not free: it is
the root of

    V(d) = displacement volume          V monotone increasing in d,

which is solved by bisection on the same clip. With ``B`` the centroid of the submerged
polyhedron and ``G = (0, 0, keel + KG)``, the righting arm is the transverse lever

    m  = (0, cos phi, sin phi)          (unit, perpendicular to n)
    GZ = (B - G) . m = y_B cos phi + (z_B - z_G) sin phi

and ``GZ > 0`` means the pair (buoyancy up through B, weight down through G) rotates the
ship back towards upright. The sign convention is pinned by the small-angle limit: with
``y_B ~ BM tan phi`` and ``z_B ~ KB`` this gives ``GZ ~ GM sin phi``, so the curve's
slope at the origin *is* the ``GM`` of :mod:`tensorcad.stability` -- which is exactly the
cross-check :func:`righting_arm_curve` is validated against in ``examples/gz``.

**Why the volume is integrated the way it is.** The half-space clip turns every triangle
into a convex polygon of at most four vertices, and the clipped surface plus the planar
cap is a closed boundary. Taking the apex of the tetrahedral fan *in the cutting plane*
makes the cap contribute exactly zero volume and zero first moment, so

    V   = sum_faces V_tet(O, a, b, c)
    M   = sum_faces V_tet * (O + a + b + c) / 4

need no cap at all and stay exact for the polyhedron the mesh actually is. The only
discretisation in the whole chain is the mesh's own faceting -- a box reproduces its
closed form to round-off, which is what ``tests/test_gz.py`` pins.

**A second, independent route.** :func:`cross_curve_righting_arm` gets the same ``GZ``
from the classical cross-curves: cut the hull at each station (via
:func:`tensorcad.bonjean.section_loops`, a 2-D plane cut), clip each station outline with
the *inclined* waterline, integrate area and centroid along ``x``. It shares nothing with
the 3-D route except the mesh, so agreement between them is evidence about the
integration rather than about one implementation's habits.

**What this is not.** There is no free-surface correction (a slack tank's
``sum(rho i_x) / displacement`` must be subtracted by the caller who knows the tankage),
no wind heeling arm, no damage/permability case, and no class-rule criterion: the
``area under the curve to 30 deg / 40 deg`` that
:class:`RightingArmCurve.area_under_curve` reports is the *dynamical stability* used by
those criteria, not a pass/fail verdict. Longitudinal righting (trim) is not integrated.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import Any

from . import mesh as meshlib
from .bonjean import section_loops
from .geometry import Mesh

#: Default heel stations of a righting-arm curve, in degrees.
DEFAULT_HEEL_ANGLES: tuple[float, ...] = (0.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 75.0, 90.0)

#: Relative displacement tolerance accepted by the waterline solve.
VOLUME_TOLERANCE = 1e-10

#: Hard cap on bisection steps of the waterline solve (each step halves the bracket).
MAX_WATERLINE_STEPS = 200


class GZError(ValueError):
    """Raised when a righting arm cannot be defined for the given input."""


# ---------------------------------------------------------------------------
# results
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class HeelPoint:
    """One point of the righting-arm curve: the floating body at one heel angle."""

    heel_degrees: float
    submerged_volume: float
    target_volume: float
    transverse_centre_buoyancy: float
    centre_buoyancy_above_keel: float
    waterline_offset: float
    righting_arm: float
    iterations: int

    @property
    def volume_error(self) -> float:
        """Signed relative error of the displacement balance, ``V - V0`` over ``V0``."""
        return (self.submerged_volume - self.target_volume) / self.target_volume

    def to_dict(self) -> dict[str, Any]:
        return {
            "heel_degrees": _round(self.heel_degrees, 9),
            "submerged_volume_m3": _round(self.submerged_volume, 12),
            "target_volume_m3": _round(self.target_volume, 12),
            "volume_relative_error": _round(self.volume_error, 15),
            "transverse_centre_buoyancy_m": _round(self.transverse_centre_buoyancy, 12),
            "centre_buoyancy_above_keel_m": _round(self.centre_buoyancy_above_keel, 12),
            "waterline_offset_m": _round(self.waterline_offset, 12),
            "righting_arm_m": _round(self.righting_arm, 12),
            "waterline_iterations": self.iterations,
        }


@dataclass(frozen=True)
class RightingArmCurve:
    """``GZ`` against heel for one loading condition, plus its summary numbers."""

    points: tuple[HeelPoint, ...]
    kg: float
    displacement_volume: float
    water_density: float
    method: str
    basis: str

    @property
    def heel_angles(self) -> tuple[float, ...]:
        return tuple(point.heel_degrees for point in self.points)

    @property
    def righting_arms(self) -> tuple[float, ...]:
        return tuple(point.righting_arm for point in self.points)

    @property
    def maximum_righting_arm(self) -> float:
        return max(point.righting_arm for point in self.points)

    @property
    def angle_of_maximum(self) -> float:
        """Heel angle of the largest righting arm, at the resolution of the sample grid."""
        return max(self.points, key=lambda point: point.righting_arm).heel_degrees

    @property
    def range_of_stability(self) -> float | None:
        """Heel angle at which ``GZ`` returns to zero, by linear interpolation.

        ``None`` when the curve is still positive at the last sampled angle, or never
        becomes positive at all: the first case means the sampling stopped too early,
        the second that the loading is unstable and there is no range to report. The
        interpolation is between two sampled points and is therefore first-order
        accurate in the angle step; it is a diagnostic, not a certificate.
        """
        # An upright body's arm is zero to round-off, so "still positive" has to mean
        # positive at the curve's own scale; otherwise a 1e-15 residue at phi = 0 would
        # be interpolated into a range of stability of 1e-15 degrees.
        scale = max(abs(point.righting_arm) for point in self.points)
        floor = 1e-9 * scale
        previous = self.points[0]
        if previous.righting_arm <= floor:
            return None
        for point in self.points[1:]:
            if point.righting_arm <= 0.0:
                span = point.righting_arm - previous.righting_arm
                if span == 0.0:
                    return point.heel_degrees
                fraction = -previous.righting_arm / span
                return previous.heel_degrees + fraction * (point.heel_degrees - previous.heel_degrees)
            previous = point
        return None

    def area_under_curve(self, limit_degrees: float) -> float:
        """Dynamical stability up to ``limit_degrees``: ``integral GZ dphi`` in rad*m.

        Trapezoidal on the sampled grid, truncated at the limit by interpolating the last
        partial panel, so ``area_under_curve(30.0)`` matches the 30 deg criterion's
        quantity rather than the nearest sampled angle.
        """
        if limit_degrees <= 0.0:
            raise GZError(f"limit_degrees must be positive, got {limit_degrees!r}")
        total = 0.0
        for first, second in zip(self.points, self.points[1:], strict=False):
            a0, a1 = first.heel_degrees, second.heel_degrees
            if a0 >= limit_degrees:
                break
            if a1 > limit_degrees:
                fraction = (limit_degrees - a0) / (a1 - a0)
                gz1 = first.righting_arm + fraction * (second.righting_arm - first.righting_arm)
                total += 0.5 * (first.righting_arm + gz1) * math.radians(limit_degrees - a0)
                break
            total += 0.5 * (first.righting_arm + second.righting_arm) * math.radians(a1 - a0)
        return total

    def to_dict(self) -> dict[str, Any]:
        return {
            "method": self.method,
            "basis": self.basis,
            "kg_above_keel_m": _round(self.kg, 12),
            "displacement_volume_m3": _round(self.displacement_volume, 12),
            "water_density_kg_m3": self.water_density,
            "maximum_righting_arm_m": _round(self.maximum_righting_arm, 12),
            "angle_of_maximum_degrees": _round(self.angle_of_maximum, 9),
            "range_of_stability_degrees": (
                None if self.range_of_stability is None else _round(self.range_of_stability, 9)
            ),
            "area_under_curve_0_30_m_rad": _round(self.area_under_curve(30.0), 12),
            "area_under_curve_0_40_m_rad": _round(self.area_under_curve(40.0), 12),
            "points": [point.to_dict() for point in self.points],
        }


@dataclass(frozen=True)
class SubmergedProperties:
    """Volume and centroid of the part of a closed mesh below an arbitrary plane."""

    volume: float
    centroid: tuple[float, float, float]
    clipped_faces: int
    kept_faces: int

    @property
    def transverse_centre(self) -> float:
        return self.centroid[1]


def _round(value: float, digits: int = 9) -> float:
    return float(f"{value:.{digits}g}")


# ---------------------------------------------------------------------------
# half-space integration
# ---------------------------------------------------------------------------
def _unit_triple(normal: Sequence[float]) -> tuple[float, float, float]:
    """Normalise a 3-vector into a fixed-length tuple (a generator would not typecheck)."""
    length = math.sqrt(normal[0] ** 2 + normal[1] ** 2 + normal[2] ** 2)
    return (normal[0] / length, normal[1] / length, normal[2] / length)


def _plane_point(normal: Sequence[float], offset: float) -> tuple[float, float, float]:
    """A point in the plane ``p . normal = offset``: the foot of the normal."""
    return (normal[0] * offset, normal[1] * offset, normal[2] * offset)


def _clip_triangle(
    a: Sequence[float],
    b: Sequence[float],
    c: Sequence[float],
    normal: Sequence[float],
    offset: float,
    heights: tuple[float, float, float],
) -> list[tuple[float, float, float]]:
    """Clip one triangle to ``p . normal <= offset``; return the polygon's vertices."""
    corners = (a, b, c)
    polygon: list[tuple[float, float, float]] = []
    for index in range(3):
        current = corners[index]
        nxt = corners[(index + 1) % 3]
        h_current = heights[index]
        h_next = heights[(index + 1) % 3]
        inside_current = h_current <= 0.0
        if inside_current:
            polygon.append((float(current[0]), float(current[1]), float(current[2])))
        if inside_current != (h_next <= 0.0):
            # Signed heights are distances along the unit normal, so the crossing
            # parameter is exact in that coordinate and needs no rescaling.
            t = h_current / (h_current - h_next)
            polygon.append(
                (
                    float(current[0]) + t * (nxt[0] - current[0]),
                    float(current[1]) + t * (nxt[1] - current[1]),
                    float(current[2]) + t * (nxt[2] - current[2]),
                )
            )
    return polygon


def _fan_volume_moment(
    polygon: Sequence[Sequence[float]],
    apex: tuple[float, float, float],
) -> tuple[float, float, float, float]:
    """Volume and first moments of the cone ``apex`` x ``polygon`` (fan of triangles)."""
    volume = 0.0
    moment_x = 0.0
    moment_y = 0.0
    moment_z = 0.0
    p0 = polygon[0]
    for index in range(1, len(polygon) - 1):
        p1 = polygon[index]
        p2 = polygon[index + 1]
        ax, ay, az = p0[0] - apex[0], p0[1] - apex[1], p0[2] - apex[2]
        bx, by, bz = p1[0] - apex[0], p1[1] - apex[1], p1[2] - apex[2]
        cx, cy, cz = p2[0] - apex[0], p2[1] - apex[1], p2[2] - apex[2]
        double = ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx)
        piece = double / 6.0
        volume += piece
        if piece == 0.0:
            continue
        # The tetrahedron's centroid is (apex + p0 + p1 + p2) / 4.
        moment_x += piece * 0.25 * (apex[0] + p0[0] + p1[0] + p2[0])
        moment_y += piece * 0.25 * (apex[1] + p0[1] + p1[1] + p2[1])
        moment_z += piece * 0.25 * (apex[2] + p0[2] + p1[2] + p2[2])
    return volume, moment_x, moment_y, moment_z


def _submerged(
    mesh: Mesh,
    unit: tuple[float, float, float],
    scaled_offset: float,
) -> SubmergedProperties:
    """Clip-and-integrate without the positivity contract; ``solve_waterline`` needs 0."""
    apex = _plane_point(unit, scaled_offset)
    vertices = mesh.vertices
    height = [
        vertex[0] * unit[0] + vertex[1] * unit[1] + vertex[2] * unit[2] - scaled_offset
        for vertex in vertices
    ]

    volume = 0.0
    moment_x = 0.0
    moment_y = 0.0
    moment_z = 0.0
    kept = 0
    clipped = 0
    polygon: Sequence[Sequence[float]]
    for face in mesh.faces:
        ha, hb, hc = height[face[0]], height[face[1]], height[face[2]]
        if ha <= 0.0 and hb <= 0.0 and hc <= 0.0:
            kept += 1
            polygon = (vertices[face[0]], vertices[face[1]], vertices[face[2]])
        elif ha > 0.0 and hb > 0.0 and hc > 0.0:
            continue
        else:
            clipped += 1
            polygon = _clip_triangle(
                vertices[face[0]],
                vertices[face[1]],
                vertices[face[2]],
                unit,
                scaled_offset,
                (ha, hb, hc),
            )
            if len(polygon) < 3:
                continue
        piece, px, py, pz = _fan_volume_moment(polygon, apex)
        volume += piece
        moment_x += px
        moment_y += py
        moment_z += pz

    if volume == 0.0:
        centroid = (float("nan"), float("nan"), float("nan"))
    else:
        centroid = (moment_x / volume, moment_y / volume, moment_z / volume)
    return SubmergedProperties(
        volume=volume,
        centroid=centroid,
        clipped_faces=clipped,
        kept_faces=kept,
    )


def submerged_properties(
    mesh: Mesh,
    normal: Sequence[float],
    offset: float,
) -> SubmergedProperties:
    """Volume and centroid of ``{p : p . normal <= offset}`` inside a closed mesh.

    ``normal`` need not be normalised; it is scaled internally. The apex of the
    tetrahedral fan sits in the cutting plane, so the cap that closes the clipped
    surface contributes zero volume and zero first moment and is never built --
    which is what makes the result exact for the clipped polyhedron instead of
    approximate at the waterline.

    The mesh is expected to be closed and *outward*-oriented. An inverted or open
    surface is rejected rather than returning a plausible-looking negative volume:
    every number downstream of here would be wrong in a way that looks like physics.
    """
    length = math.sqrt(sum(component * component for component in normal))
    if not math.isfinite(length) or length == 0.0:
        raise GZError(f"normal must be a finite non-zero vector, got {tuple(normal)!r}")
    if not math.isfinite(offset):
        raise GZError(f"offset must be finite, got {offset!r}")
    unit = _unit_triple(normal)
    properties = _submerged(mesh, unit, offset / length)
    if properties.volume <= 0.0 or not math.isfinite(properties.volume):
        raise GZError(
            f"the half-space p.n <= {offset!r} bounds no volume inside the mesh "
            f"(integrated {properties.volume!r} m^3): the plane misses the body, "
            "or the surface is open or inward-oriented"
        )
    return properties


def _normal_bracket(mesh: Mesh, normal: Sequence[float]) -> tuple[float, float]:
    projections = [
        vertex[0] * normal[0] + vertex[1] * normal[1] + vertex[2] * normal[2] for vertex in mesh.vertices
    ]
    return min(projections), max(projections)


def solve_waterline(
    mesh: Mesh,
    normal: Sequence[float],
    target_volume: float,
    *,
    tolerance: float = VOLUME_TOLERANCE,
    max_steps: int = MAX_WATERLINE_STEPS,
) -> tuple[float, SubmergedProperties, int]:
    """Find the plane offset whose submerged volume equals ``target_volume``.

    ``V(d)`` is strictly increasing in ``d`` for a closed body, so bisection cannot
    be trapped by a local flat spot the way a Newton step on a nearly
    vertical-sided hull can (``dV/dd`` is the waterplane area, which collapses to
    zero at a deck edge -- precisely where a large-angle curve spends its time).
    Returns the offset, the properties at it, and the number of steps used.
    """
    if not math.isfinite(target_volume) or target_volume <= 0.0:
        raise GZError(f"target_volume must be a finite positive volume, got {target_volume!r}")
    length = math.sqrt(sum(component * component for component in normal))
    if not math.isfinite(length) or length == 0.0:
        raise GZError(f"normal must be a finite non-zero vector, got {tuple(normal)!r}")
    low, high = _normal_bracket(mesh, normal)
    low -= abs(low) * 1e-15 + 1e-12
    high += abs(high) * 1e-15 + 1e-12

    unit = _unit_triple(normal)
    low_properties = _submerged(mesh, unit, low)
    if low_properties.volume >= target_volume:
        raise GZError(
            f"the target volume {target_volume!r} m^3 is not above the smallest bracketed volume "
            f"{low_properties.volume!r} m^3"
        )
    total = _total_volume(mesh)
    if target_volume > total * (1.0 + 1e-12):
        raise GZError(
            f"the target volume {target_volume!r} m^3 is not below the mesh volume {total!r} m^3: "
            "the body cannot float on the requested waterline"
        )

    middle = 0.5 * (low + high)
    properties = low_properties
    for step in range(1, max_steps + 1):
        middle = 0.5 * (low + high)
        properties = _submerged(mesh, unit, middle)
        error = properties.volume - target_volume
        if abs(error) <= tolerance * target_volume:
            return middle, properties, step
        if error < 0.0:
            low = middle
        else:
            high = middle
    return middle, properties, max_steps


def _total_volume(mesh: Mesh) -> float:
    """Enclosed volume of the whole mesh, by the same clipped-path integral at ``+inf``."""
    high = _normal_bracket(mesh, (0.0, 0.0, 1.0))[1]
    return submerged_properties(mesh, (0.0, 0.0, 1.0), high + abs(high) + 1.0).volume


# ---------------------------------------------------------------------------
# righting arms
# ---------------------------------------------------------------------------
def _heel_geometry(
    heel_degrees: float,
) -> tuple[tuple[float, float, float], tuple[float, float, float], float]:
    """Up-normal, transverse lever direction and the rotation sine/cosine."""
    if not math.isfinite(heel_degrees):
        raise GZError(f"heel_degrees must be finite, got {heel_degrees!r}")
    if abs(heel_degrees) > 180.0:
        raise GZError(f"heel_degrees must lie in [-180, 180], got {heel_degrees!r}")
    phi = math.radians(heel_degrees)
    sine = math.sin(phi)
    cosine = math.cos(phi)
    normal = (0.0, -sine, cosine)
    lever = (0.0, cosine, sine)
    return normal, lever, phi


def _resolve_displacement(
    mesh: Mesh,
    waterline_z: float,
    displacement_volume: float | None,
    water_density: float,
) -> float:
    if displacement_volume is not None:
        if not math.isfinite(displacement_volume) or displacement_volume <= 0.0:
            raise GZError(
                f"displacement_volume must be a finite positive volume, got {displacement_volume!r}"
            )
        return displacement_volume
    # The displaced volume is the part of the hull *below* the waterline, so it comes
    # from the same half-space clip the heeled states use -- not from
    # ``compute_hydrostatics``, whose volume is the whole closed surface and therefore
    # the whole hull however high the waterline sits.
    upright = submerged_properties(mesh, (0.0, 0.0, 1.0), waterline_z)
    return upright.volume


def keel_height(mesh: Mesh) -> float:
    """``z`` of the lowest vertex -- the datum every righting arm here is measured from."""
    return min(vertex[2] for vertex in mesh.vertices)


def righting_arm(
    mesh: Mesh,
    heel_degrees: float,
    *,
    kg: float,
    displacement_volume: float | None = None,
    waterline_z: float = 0.0,
    water_density: float = 1025.0,
    tolerance: float = VOLUME_TOLERANCE,
) -> HeelPoint:
    """One point of the righting-arm curve, from the 3-D heeled displacement.

    ``kg`` is the centre of gravity height **above the keel** (the mesh's lowest
    vertex), matching :func:`tensorcad.stability.initial_stability`; it is a loading
    input and has no default. ``displacement_volume`` defaults to the upright
    displacement at ``waterline_z``, so ``righting_arm(mesh, 0.0, kg=..)`` is the
    upright state and its ``GZ`` is zero for a symmetric hull to round-off.
    """
    if not math.isfinite(kg):
        raise GZError(f"kg must be finite, got {kg!r}")
    normal, lever, _ = _heel_geometry(heel_degrees)
    target = _resolve_displacement(mesh, waterline_z, displacement_volume, water_density)
    offset, properties, steps = solve_waterline(mesh, normal, target, tolerance=tolerance)
    keel = keel_height(mesh)
    centroid = properties.centroid
    gravity = (0.0, 0.0, keel + kg)
    arm = (
        (centroid[0] - gravity[0]) * lever[0]
        + (centroid[1] - gravity[1]) * lever[1]
        + (centroid[2] - gravity[2]) * lever[2]
    )
    return HeelPoint(
        heel_degrees=heel_degrees,
        submerged_volume=properties.volume,
        target_volume=target,
        transverse_centre_buoyancy=centroid[1],
        centre_buoyancy_above_keel=centroid[2] - keel,
        waterline_offset=offset,
        righting_arm=arm,
        iterations=steps,
    )


def righting_arm_curve(
    mesh: Mesh,
    angles: Iterable[float] | None = None,
    *,
    kg: float,
    displacement_volume: float | None = None,
    waterline_z: float = 0.0,
    water_density: float = 1025.0,
    closed: bool | None = None,
) -> RightingArmCurve:
    """``GZ`` over a list of heel angles (default :data:`DEFAULT_HEEL_ANGLES`)."""
    grid = tuple(DEFAULT_HEEL_ANGLES if angles is None else angles)
    if not grid:
        raise GZError("at least one heel angle is required")
    is_closed = meshlib.is_watertight(mesh) if closed is None else closed
    target = _resolve_displacement(mesh, waterline_z, displacement_volume, water_density)
    points = tuple(
        righting_arm(
            mesh,
            angle,
            kg=kg,
            displacement_volume=target,
            waterline_z=waterline_z,
            water_density=water_density,
        )
        for angle in grid
    )
    basis = (
        "3-D clip: exact for the faceted hull, no free-surface correction"
        if is_closed
        else "OPEN surface: volume is not trustworthy, righting arms are not evidence"
    )
    return RightingArmCurve(
        points=points,
        kg=kg,
        displacement_volume=target,
        water_density=water_density,
        method="heeled displacement (half-space clip + bisection on volume)",
        basis=basis,
    )


# ---------------------------------------------------------------------------
# second route: cross curves from station cuts
# ---------------------------------------------------------------------------
def _clip_inclined(
    polygon: Sequence[tuple[float, float]],
    sine: float,
    cosine: float,
    offset: float,
) -> list[tuple[float, float]]:
    """Keep the part of a ``(y, z)`` polygon where ``z cos phi - y sin phi <= offset``."""
    out: list[tuple[float, float]] = []
    count = len(polygon)
    for index in range(count):
        y0, z0 = polygon[index]
        y1, z1 = polygon[(index + 1) % count]
        f0 = z0 * cosine - y0 * sine - offset
        f1 = z1 * cosine - y1 * sine - offset
        inside0 = f0 <= 0.0
        if inside0:
            out.append((y0, z0))
        if inside0 != (f1 <= 0.0):
            t = f0 / (f0 - f1)
            out.append((y0 + t * (y1 - y0), z0 + t * (z1 - z0)))
    return out


def _polygon_area_centroid(
    polygon: Sequence[tuple[float, float]],
) -> tuple[float, float, float] | None:
    """Area and centroid of a simple ``(y, z)`` polygon, by the shoelace pairing."""
    if len(polygon) < 3:
        return None
    area2 = 0.0
    moment_y6 = 0.0
    moment_z6 = 0.0
    count = len(polygon)
    for index in range(count):
        y0, z0 = polygon[index]
        y1, z1 = polygon[(index + 1) % count]
        cross = y0 * z1 - y1 * z0
        area2 += cross
        moment_y6 += (y0 + y1) * cross
        moment_z6 += (z0 + z1) * cross
    if area2 == 0.0:
        return None
    signed = 0.5 * area2
    # The traversal sign cancels between the moment and the area, so the centroid is
    # (moment / 6) / signed area whichever way the loop runs.
    return abs(signed), moment_y6 / (6.0 * signed), moment_z6 / (6.0 * signed)


def _integrate_along_x(stations: Sequence[float], values: Sequence[float]) -> float:
    """Integrate ``values`` over ``stations``: Simpson when the grid allows it, else trapezoid."""
    if len(stations) < 2:
        return 0.0
    if len(stations) % 2 == 1:
        total = 0.0
        for index in range(0, len(stations) - 2, 2):
            width = stations[index + 2] - stations[index]
            total += width * (values[index] + 4.0 * values[index + 1] + values[index + 2]) / 6.0
        return total
    total = 0.0
    for index in range(len(stations) - 1):
        width = stations[index + 1] - stations[index]
        total += 0.5 * width * (values[index] + values[index + 1])
    return total


def station_offsets(mesh: Mesh, count: int | None = None) -> tuple[float, ...]:
    """Station positions along ``x``: the mesh's own vertex rings, or ``count`` even slices.

    The default uses every distinct vertex coordinate, because those are exactly the
    rings the loft put there and the station cut is then exact rather than
    interpolated between rings.
    """
    if count is None:
        values = sorted({vertex[0] for vertex in mesh.vertices})
        if len(values) >= 3:
            return tuple(values)
    box = mesh.bounding_box()
    x0 = box.minimum[0]
    x1 = box.maximum[0]
    slices = 121 if count is None else count
    if slices < 3:
        raise GZError(f"count must be at least 3 slices, got {count!r}")
    return tuple(x0 + (x1 - x0) * index / (slices - 1) for index in range(slices))


def cross_curve_righting_arm(
    mesh: Mesh,
    heel_degrees: float,
    *,
    kg: float,
    displacement_volume: float | None = None,
    waterline_z: float = 0.0,
    water_density: float = 1025.0,
    stations: Sequence[float] | None = None,
    tolerance: float = VOLUME_TOLERANCE,
    outlines: Sequence[Sequence[tuple[tuple[float, float], ...]]] | None = None,
) -> HeelPoint:
    """``GZ`` from cross curves: inclined-waterline areas integrated along the hull.

    The hull is cut at each station into ``(y, z)`` outlines (a 2-D operation), each
    outline is clipped with the *inclined* waterline, and area and centroid are
    integrated along ``x``. Nothing here uses the 3-D surface integration of
    :func:`righting_arm`, so a comparison between the two is a check on the
    integration and not a self-consistency exercise.
    """
    if not math.isfinite(kg):
        raise GZError(f"kg must be finite, got {kg!r}")
    normal, lever, _ = _heel_geometry(heel_degrees)
    sine = math.sin(math.radians(heel_degrees))
    cosine = math.cos(math.radians(heel_degrees))
    target = _resolve_displacement(mesh, waterline_z, displacement_volume, water_density)
    grid = tuple(station_offsets(mesh) if stations is None else stations)
    cuts = tuple(section_loops(mesh, station) for station in grid) if outlines is None else outlines
    keel = keel_height(mesh)

    def state(offset: float) -> tuple[float, float, float]:
        areas: list[float] = []
        moment_y: list[float] = []
        moment_z: list[float] = []
        for loops in cuts:
            area = 0.0
            my = 0.0
            mz = 0.0
            for loop in loops:
                clipped = _clip_inclined(loop, sine, cosine, offset)
                if not clipped:
                    continue
                measured = _polygon_area_centroid(clipped)
                if measured is None:
                    continue
                piece, cy, cz = measured
                area += piece
                my += piece * cy
                mz += piece * cz
            areas.append(area)
            moment_y.append(my)
            moment_z.append(mz)
        volume = _integrate_along_x(grid, areas)
        if volume <= 0.0 or not math.isfinite(volume):
            raise GZError(
                f"the inclined waterline at offset {offset!r} bounds no volume; "
                "the plane misses the hull"
            )
        centre_y = _integrate_along_x(grid, moment_y) / volume
        centre_z = _integrate_along_x(grid, moment_z) / volume
        return volume, centre_y, centre_z

    low = min(
        vertex[1] * normal[1] + vertex[2] * normal[2] for vertex in mesh.vertices
    )
    high = max(
        vertex[1] * normal[1] + vertex[2] * normal[2] for vertex in mesh.vertices
    )
    low -= abs(low) * 1e-9 + 1e-9
    high += abs(high) * 1e-9 + 1e-9
    middle = 0.5 * (low + high)
    steps = MAX_WATERLINE_STEPS
    volume, centre_y, centre_z = state(middle)
    for step in range(1, MAX_WATERLINE_STEPS + 1):
        middle = 0.5 * (low + high)
        volume, centre_y, centre_z = state(middle)
        error = volume - target
        if abs(error) <= tolerance * target:
            steps = step
            break
        if error < 0.0:
            low = middle
        else:
            high = middle
    gravity = (0.0, 0.0, keel + kg)
    arm = centre_y * lever[1] + (centre_z - gravity[2]) * lever[2]
    return HeelPoint(
        heel_degrees=heel_degrees,
        submerged_volume=volume,
        target_volume=target,
        transverse_centre_buoyancy=centre_y,
        centre_buoyancy_above_keel=centre_z - keel,
        waterline_offset=middle,
        righting_arm=arm,
        iterations=steps,
    )


def cross_curve_righting_arm_curve(
    mesh: Mesh,
    angles: Iterable[float] | None = None,
    *,
    kg: float,
    displacement_volume: float | None = None,
    waterline_z: float = 0.0,
    water_density: float = 1025.0,
    stations: Sequence[float] | None = None,
) -> RightingArmCurve:
    """The cross-curve route over a list of heel angles (see :func:`cross_curve_righting_arm`)."""
    grid = tuple(DEFAULT_HEEL_ANGLES if angles is None else angles)
    if not grid:
        raise GZError("at least one heel angle is required")
    target = _resolve_displacement(mesh, waterline_z, displacement_volume, water_density)
    grid_stations = tuple(station_offsets(mesh) if stations is None else stations)
    outlines = tuple(section_loops(mesh, station) for station in grid_stations)
    points = tuple(
        cross_curve_righting_arm(
            mesh,
            angle,
            kg=kg,
            displacement_volume=target,
            waterline_z=waterline_z,
            water_density=water_density,
            stations=grid_stations,
            outlines=outlines,
        )
        for angle in grid
    )
    return RightingArmCurve(
        points=points,
        kg=kg,
        displacement_volume=target,
        water_density=water_density,
        method="cross curves (station cuts + inclined waterline, integrated along x)",
        basis="2-D station integration: independent of the 3-D surface route",
    )


# ---------------------------------------------------------------------------
# literature reference
# ---------------------------------------------------------------------------
def wall_sided_righting_arm(heel_degrees: float, gm: float, bm: float) -> float:
    """Closed-form ``GZ`` of a wall-sided body, ``sin(phi) * (GM + BM tan^2(phi) / 2)``.

    This is the textbook large-angle result for a hull whose sides are vertical at
    the waterline and which does not immerse a deck edge or emerge a bilge (Biran,
    *Ship Hydrostatics and Stability*, and Rawson & Tupper, *Basic Ship Theory*).
    It exists here as the reference a *box* righting arm must reproduce, not as an
    approximation to a real hull: a ship with flare, tumblehome or a bilge radius
    leaves the formula's assumptions immediately, and the cross- and 3-D
    integrations are the honest answer for those.
    """
    phi = math.radians(heel_degrees)
    return math.sin(phi) * (gm + 0.5 * bm * math.tan(phi) ** 2)
