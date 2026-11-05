"""Large-angle stability: the ``GZ`` curve from heeled plane cuts.

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
from dataclasses import dataclass
from typing import Any

from .geometry import Mesh
from .mesh import total_signed_volume
from .stability import (
    PLANE_TOLERANCE,
    StabilityError,
    _bounding_scale,
    waterplane_moments,
)

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


def _round(value: float, digits: int = 9) -> float:
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
            "volume_m3": _round(self.volume),
            "centroid_m": tuple(_round(value) for value in self.centroid),
            "waterplane_area_m2": _round(self.waterplane_area),
            "waterline_z_m": _round(self.waterline_z),
            "section_segments": self.section_segments,
            "faces_below": self.faces_below,
            "faces_clipped": self.faces_clipped,
            "faces_coplanar": self.faces_coplanar,
            "basis": self.basis,
        }


def _clip_below(points: list[Vec3], waterline_z: float, tolerance: float) -> list[Vec3]:
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
            clipped = _clip_below(points, waterline_z, tolerance)
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
            "heel_angle_deg": _round(self.heel_angle_deg, 6),
            "kn_m": _round(self.kn),
            "kg_m": _round(self.kg),
            "gz_m": _round(self.gz),
            "displacement_volume_m3": _round(self.displacement_volume),
            "waterline_z_m": _round(self.waterline_z),
            "centre_buoyancy_m": tuple(_round(value) for value in self.centre_buoyancy),
            "waterplane_area_m2": _round(self.waterplane_area),
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
            "kg_m": _round(self.kg),
            "keel_z_m": _round(self.keel_z),
            "displacement_volume_m3": _round(self.displacement_volume),
            "displacement_mass_kg": _round(self.displacement_mass, 6),
            "water_density_kg_m3": self.water_density,
            "gm_upright_m": _round(self.gm_upright),
            "upright_waterline_z_m": _round(self.upright_waterline_z),
            "max_gz_m": _round(maximum),
            "angle_at_max_gz_deg": _round(maximum_angle, 6),
            # the IMO areas are only reported when the curve actually spans them: an area
            # extrapolated past the last sampled angle would be a number nobody measured
            "area_0_30_m_rad": _round(self.area(0.0, 30.0))
            if self.angles[0] <= 0.0 <= self.angles[-1] and self.angles[-1] >= 30.0
            else None,
            "area_0_40_m_rad": _round(self.area(0.0, 40.0))
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
            "measured": _round(self.measured, 6),
            "unit": self.unit,
            "satisfied": self.satisfied,
            "margin": _round(self.margin, 6),
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
            "gm_upright_m": _round(self.gm_upright),
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


__all__ = [
    "DEFAULT_HEEL_ANGLES_DEG",
    "SCHEMA_GZ_CURVE",
    "SCHEMA_IMO",
    "SCHEMA_SUBMERGED",
    "CriterionCheck",
    "GZCurve",
    "GZError",
    "HeelPoint",
    "IMOIntactStability",
    "SubmergedSolid",
    "gz_curve",
    "heel_mesh",
    "imo_intact_stability",
    "solve_waterline",
    "submerged_solid",
]

# ``StabilityError`` is re-exported so callers that catch the initial-stability error also
# catch the heeled-waterline one without importing two modules.
StabilityError = StabilityError
