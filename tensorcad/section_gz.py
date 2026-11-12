"""Exact righting arms for a *prismatic* body: one section, every heel angle.

:mod:`tensorcad.gz` integrates a righting-arm curve off a closed triangle mesh and checks
it against two closed forms -- a box barge (the wall-sided formula) and a circular section
floating at *half* submergence (``GZ = (R - KG) sin(phi)``).  Both are exact, but both are
special: the box is exact only while its sides stay wall-sided, and the circular section is
degenerate, because a circle is invariant under the heel rotation, so its heeled hydrostatics
never change shape and the arm is a pure sine at *any* draft.  Neither can say anything about
a section that is neither wall-sided nor rotationally symmetric -- a V-shaped forefoot, a
hard-chine pontoon, an elliptical strut -- which is exactly where ``GZ(phi)`` stops being a
sine and where a designer needs to know the shape of the curve.

This module supplies the missing reference: an **exact** righting arm for a prismatic body
whose transverse section is known in closed form.  A prism's heeled displacement is a purely
two-dimensional problem -- the volume is the section area times the length, and the buoyancy
centre is the section centroid -- so the arm can be computed without any discretisation at
all, for any heel angle, at any draft.

Two section families are implemented:

* :func:`polygon_prism_righting_arms` -- any *convex polygon* section.  At each heel the
  section is rotated, cut by the (world-horizontal) waterline, and the submerged part's area
  and centroid come from a half-plane clip of the polygon and Green's theorem.  Both are
  exact for the polygon, so the oracle has no discretisation error and the waterline level is
  the only solved quantity (bisection on the area, monotone for a convex section).
* :func:`ellipse_prism_righting_arms` -- an *elliptical* section, closed form.  The heeled
  unit disk cut by a line at signed distance ``d`` has area
  ``pi/2 + asin(d) + d sqrt(1 - d**2)`` and first moment ``-(2/3) (1 - d**2)**(3/2)`` about
  its centre; an affine map carries that to the ellipse, so the exact arm needs no polygon at
  all and no discretisation to converge away.  This is the reference a *smooth* non-circular
  section has never had here.

**Sign convention.**  Identical to :mod:`tensorcad.gz`, so the two can be compared point by
point: ``phi > 0`` heels the ``+y`` side down, ``KN`` is the buoyancy centre's lever about the
keel line (``KN = y_B``), and ``GZ = KN - KG sin(phi)`` with ``KG`` measured above the keel.
The section is given in the body frame the mesh would be given in, and the keel is its lowest
point; ``waterline_z`` is measured in that same frame, exactly as
:func:`tensorcad.gz.righting_arm_curve` measures it.

**What this is not.**  Prismatic means *no variation along x*: no rocker, no flare, no
entrance or run, so the arm is that of an infinitely (or uniformly) long body of that section.
It validates the *integrator*, not a hull.  Nor does it model tanks, damage, deck-edge
immersion, or any statutory criterion.  Sections must be convex: a concave section can have
two waterline levels at one displacement, and :func:`_validate_section` rejects it rather than
silently returning the first root bisection happens to bracket.

Only the Python standard library is used.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

from .geometry import Mesh, MeshBuilder
from .gz import GzError

SCHEMA_SECTION_ARM = "tensorcad.stability.section_arm/1"
SCHEMA_SECTION_CURVE = "tensorcad.stability.section_arm_curve/1"

Point2 = tuple[float, float]

#: Bisection budget for the waterline level.  200 halvings take a double interval to its
#: neighbouring floats from anywhere, so the level is solved to round-off, not to a tolerance
#: the caller has to trust.
_LEVEL_ITERATIONS = 200


@dataclass(frozen=True)
class SectionArm:
    """One point of an exact righting-arm curve, with the quantities behind it."""

    heel_angle_deg: float
    righting_arm: float
    kn: float
    submerged_area: float
    waterline_level: float
    centre_y: float
    centre_z: float

    def to_dict(self) -> dict[str, float]:
        return {
            "heel_angle_deg": self.heel_angle_deg,
            "righting_arm": self.righting_arm,
            "kn": self.kn,
            "submerged_area": self.submerged_area,
            "waterline_level": self.waterline_level,
            "centre_y": self.centre_y,
            "centre_z": self.centre_z,
        }


@dataclass(frozen=True)
class SectionArmCurve:
    """Exact righting arms of a prismatic body, one :class:`SectionArm` per heel angle."""

    arms: tuple[SectionArm, ...]
    kg: float
    displaced_area: float
    length: float
    section: str

    @property
    def gz(self) -> tuple[float, ...]:
        return tuple(arm.righting_arm for arm in self.arms)

    @property
    def angles_deg(self) -> tuple[float, ...]:
        return tuple(arm.heel_angle_deg for arm in self.arms)

    def angle_of(self, angle_deg: float) -> SectionArm:
        for arm in self.arms:
            if arm.heel_angle_deg == angle_deg:
                return arm
        raise KeyError(f"heel angle {angle_deg!r} is not on this curve")

    def max_arm(self) -> tuple[float, float] | None:
        """``(GZ_max, heel of GZ_max)`` over the sampled angles, or ``None`` if never righting."""
        if not self.arms:
            return None
        best = max(self.arms, key=lambda arm: arm.righting_arm)
        if best.righting_arm <= 0.0:
            return None
        return best.righting_arm, best.heel_angle_deg

    def vanishing_stability_angle_deg(self) -> float | None:
        """First sign change of ``GZ`` after its maximum, by linear interpolation.

        ``None`` when the sampled curve never returns to zero -- which is the honest answer
        for a short angle list, and the reason this returns an angle *and* the caller can see
        the sampled points it was interpolated from.
        """
        angles, arms = self.angles_deg, self.arms
        if len(angles) < 2:
            return None
        start = 0
        values = [arm.righting_arm for arm in arms]
        if max(values) <= 0.0:
            return None
        start = values.index(max(values))
        for index in range(start, len(values) - 1):
            if values[index] > 0.0 >= values[index + 1]:
                left, right = values[index], values[index + 1]
                fraction = left / (left - right)
                return angles[index] + fraction * (angles[index + 1] - angles[index])
        return None

    def area_under_curve(self, limit_deg: float | None = None) -> float:
        """Trapezoidal dynamic-stability area, ``int GZ dphi`` in metre-radians."""
        total = 0.0
        for left, right in zip(self.arms, self.arms[1:], strict=False):
            if limit_deg is not None and right.heel_angle_deg > limit_deg:
                break
            total += (
                0.5
                * (left.righting_arm + right.righting_arm)
                * math.radians(right.heel_angle_deg - left.heel_angle_deg)
            )
        return total

    def to_dict(self) -> dict[str, object]:
        return {
            "schema": SCHEMA_SECTION_CURVE,
            "section": self.section,
            "kg": self.kg,
            "displaced_area": self.displaced_area,
            "length": self.length,
            "arms": [arm.to_dict() for arm in self.arms],
        }


def _validate_section(section: Sequence[Point2]) -> tuple[Point2, ...]:
    """Return the section as a counter-clockwise tuple, rejecting what the maths cannot take."""
    points = tuple((float(y), float(z)) for (y, z) in section)
    if len(points) < 3:
        raise GzError(f"a section needs at least 3 vertices, got {len(points)}")
    for index, (y, z) in enumerate(points):
        if not (math.isfinite(y) and math.isfinite(z)):
            raise GzError(f"section vertex {index} is not finite: {(y, z)!r}")
    if len(set(points)) != len(points):
        raise GzError("the section repeats a vertex, so it is not a simple polygon")
    area = _polygon_area(points)
    if abs(area) <= 1e-15:
        raise GzError("the section has zero area")
    ordered = points if area > 0.0 else tuple(reversed(points))
    # convexity: every vertex turns the same way (collinear vertices are allowed, a
    # straight run of edges is a legal polygon boundary and clipping still works)
    sign = 0
    for index in range(len(ordered)):
        a, b, c = ordered[index], ordered[(index + 1) % len(ordered)], ordered[(index + 2) % len(ordered)]
        cross = (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0])
        if abs(cross) <= 1e-14 * max(1.0, abs(area)):
            continue
        current = 1 if cross > 0.0 else -1
        if sign == 0:
            sign = current
        elif current != sign:
            raise GzError(
                "the section is not convex: a concave section can float at two waterline "
                "levels with the same displacement, so the level solve is not unique "
                "(use the mesh path in tensorcad.gz instead)"
            )
    if sign == 0:
        raise GzError("the section is degenerate (all vertices are collinear)")
    return ordered


def _polygon_area(points: Sequence[Point2]) -> float:
    """Signed area by the shoelace rule; exact for the polygon."""
    total = 0.0
    for index, (y, z) in enumerate(points):
        next_y, next_z = points[(index + 1) % len(points)]
        total += y * next_z - next_y * z
    return 0.5 * total


def _area_and_centroid(points: Sequence[Point2]) -> tuple[float, float, float]:
    """``(area, centroid_y, centroid_z)`` of a simple polygon, via Green's theorem.

    Exact for the polygon: this is the same identity :mod:`tensorcad.stability` uses on the
    waterplane, applied to a 2D section instead of a 3D surface.
    """
    count = len(points)
    if count < 3:
        return 0.0, 0.0, 0.0
    twice_area = 0.0
    moment_y = 0.0
    moment_z = 0.0
    for index in range(count):
        y0, z0 = points[index]
        y1, z1 = points[(index + 1) % count]
        cross = y0 * z1 - y1 * z0
        twice_area += cross
        moment_y += (y0 + y1) * cross
        moment_z += (z0 + z1) * cross
    if twice_area == 0.0:
        return 0.0, 0.0, 0.0
    area = 0.5 * twice_area
    return area, moment_y / (3.0 * twice_area), moment_z / (3.0 * twice_area)


def _clip_below(section: Sequence[Point2], level: float) -> list[Point2]:
    """Sutherland-Hodgman clip of a convex section against ``z <= level`` (world frame)."""
    clipped: list[Point2] = []
    for index, start in enumerate(section):
        end = section[(index + 1) % len(section)]
        depth_start = start[1] - level
        depth_end = end[1] - level
        if depth_start <= 0.0:
            clipped.append(start)
        if (depth_start < 0.0 < depth_end) or (depth_end < 0.0 < depth_start):
            fraction = depth_start / (depth_start - depth_end)
            clipped.append((start[0] + fraction * (end[0] - start[0]), level))
    return clipped


def section_properties(section: Sequence[Point2], level: float) -> tuple[float, float, float]:
    """``(area, centroid_y, centroid_z)`` of the part of ``section`` below ``level``.

    The section is taken in the body frame; ``level`` is a *world* height, so callers that
    have heeled the section must pass the already-rotated points.
    """
    return _area_and_centroid(_clip_below(section, level))


def _bisect_level(area_of_level, target: float, low: float, high: float) -> float:
    """Solve ``area_of_level(level) = target`` on a bracket where the area is increasing."""
    low_area = area_of_level(low)
    high_area = area_of_level(high)
    if not (low_area <= target <= high_area):
        raise GzError(
            f"the displaced area {target!r} m^2 is outside what the section holds between "
            f"levels {low!r} and {high!r} ({low_area!r} .. {high_area!r} m^2)"
        )
    for _ in range(_LEVEL_ITERATIONS):
        middle = 0.5 * (low + high)
        if middle <= low or middle >= high:  # the bracket has collapsed to neighbouring floats
            break
        if area_of_level(middle) < target:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def _righting_arms(
    area_of_level_factory,
    *,
    angles_deg: Sequence[float],
    kg: float,
    displaced_area: float,
    length: float,
    section: str,
) -> SectionArmCurve:
    """Shared heel sweep: solve the level per angle, then ``GZ = KN - KG sin(phi)``."""
    if not math.isfinite(kg):
        raise GzError(f"KG must be finite: {kg!r}")
    if not math.isfinite(displaced_area) or displaced_area <= 0.0:
        raise GzError(f"the displaced area must be positive and finite: {displaced_area!r}")
    if not math.isfinite(length) or length <= 0.0:
        raise GzError(f"the prism length must be positive and finite: {length!r}")
    arms: list[SectionArm] = []
    for angle in angles_deg:
        angle = float(angle)
        if not math.isfinite(angle) or abs(angle) >= 180.0:
            raise GzError(f"heel angle must be finite and within (-180, 180) degrees: {angle!r}")
        area_of_level, rotated_properties, bounds = area_of_level_factory(math.radians(angle))
        level = _bisect_level(area_of_level, displaced_area, bounds[0], bounds[1])
        area, centre_y, centre_z = rotated_properties(level)
        kn = centre_y
        arms.append(
            SectionArm(
                heel_angle_deg=angle,
                righting_arm=kn - kg * math.sin(math.radians(angle)),
                kn=kn,
                submerged_area=area,
                waterline_level=level,
                centre_y=centre_y,
                centre_z=centre_z,
            )
        )
    return SectionArmCurve(
        arms=tuple(arms),
        kg=kg,
        displaced_area=displaced_area,
        length=length,
        section=section,
    )


def polygon_prism_righting_arms(
    section: Sequence[Point2],
    angles_deg: Sequence[float],
    *,
    kg: float,
    waterline_z: float = 0.0,
    length: float = 1.0,
    displaced_area: float | None = None,
) -> SectionArmCurve:
    """Exact righting arms of a prism with a convex polygon section.

    ``section`` is the transverse section in the body frame ``(y, z)``; the keel is its
    lowest point and ``waterline_z`` is the upright waterline height above that same datum
    (so a section drawn with its keel at ``z = 0`` takes ``waterline_z = draft``).
    ``displaced_area`` defaults to the area below ``waterline_z`` upright; pass it explicitly
    to keep every point of the curve on one displacement.
    """
    ordered = _validate_section(section)
    # the keel-shifted frame that ``tensorcad.gz`` rotates in: the lowest point of the
    # section goes to z = 0 and the heel rotation turns the section about the origin, so
    # ``KG`` is measured above the keel exactly as the mesh path measures it
    keel = min(z for (_, z) in ordered)
    shifted = tuple((y, z - keel) for (y, z) in ordered)
    base_level = waterline_z - keel
    height = max(z for (_, z) in shifted)
    if base_level <= 0.0:
        raise GzError(f"the waterline ({waterline_z!r}) must sit above the keel ({keel!r})")
    if base_level >= height:
        raise GzError(
            f"the waterline ({waterline_z!r}) must sit below the section's top "
            f"({keel + height!r}): a fully submerged section has no waterplane to heel about"
        )
    upright_area, _, _ = section_properties(shifted, base_level)
    target = upright_area if displaced_area is None else float(displaced_area)

    def factory(heel_rad: float):
        cos_heel, sin_heel = math.cos(heel_rad), math.sin(heel_rad)
        rotated = tuple(
            (y * cos_heel + z * sin_heel, -y * sin_heel + z * cos_heel) for (y, z) in shifted
        )
        heights = [z for (_, z) in rotated]
        return (
            lambda level: _area_and_centroid(_clip_below(rotated, level))[0],
            lambda level: _area_and_centroid(_clip_below(rotated, level)),
            (min(heights), max(heights)),
        )

    return _righting_arms(
        factory,
        angles_deg=angles_deg,
        kg=kg,
        displaced_area=target,
        length=length,
        section=f"polygon[{len(ordered)}]",
    )


def ellipse_prism_righting_arms(
    semi_beam: float,
    semi_height: float,
    angles_deg: Sequence[float],
    *,
    kg: float,
    waterline_z: float,
    length: float = 1.0,
    displaced_area: float | None = None,
) -> SectionArmCurve:
    """Exact righting arms of a prism with an elliptical section, in closed form.

    The section is the ellipse ``(y/a)**2 + ((z - b)/b)**2 = 1``: semi-beam ``a`` and
    semi-height ``b``, its bottom on ``z = 0`` so that ``waterline_z`` *is* the draft.  No
    polygon is involved: the heeled cut is a line against the unit disk (reached by the affine
    map ``u = y/a``, ``v = (z - b)/b``), whose segment area and first moment are elementary,
    and the affine map carries the centroid back.

    This is the reference for a smooth, *non*-circular section.  A circle cannot provide one:
    it is invariant under the heel rotation, so its arm is a sine at every draft and tests
    nothing about the shape of a real curve.
    """
    a, b = float(semi_beam), float(semi_height)
    if not (a > 0.0 and b > 0.0):
        raise GzError(f"the semi-axes must be positive: a={semi_beam!r}, b={semi_height!r}")
    draft = waterline_z
    if not (0.0 < draft < 2.0 * b):
        raise GzError(
            f"the waterline ({waterline_z!r}) must sit strictly inside the section height "
            f"(0 .. {2.0 * b!r} m above the keel)"
        )

    def disk_area(d: float) -> float:
        d = min(1.0, max(-1.0, d))
        return 0.5 * math.pi + math.asin(d) + d * math.sqrt(max(0.0, 1.0 - d * d))

    def disk_moment(d: float) -> float:
        d = min(1.0, max(-1.0, d))
        return -(2.0 / 3.0) * max(0.0, 1.0 - d * d) ** 1.5

    def factory(heel_rad: float):
        cos_heel, sin_heel = math.cos(heel_rad), math.sin(heel_rad)
        # world z of a section point, in disk coordinates:
        #   z_w = -a sin(phi) u + b cos(phi) v + b cos(phi)
        alpha = -a * sin_heel
        beta = b * cos_heel
        norm = math.hypot(alpha, beta)
        offset = b * cos_heel  # world z of the section centre
        centre_y = b * sin_heel  # world y of the section centre

        def level_of(d: float) -> float:
            return offset + norm * d

        def area_of_level(level: float) -> float:
            return a * b * disk_area((level - offset) / norm)

        def rotated_properties(level: float) -> tuple[float, float, float]:
            d = (level - offset) / norm
            area_disk = disk_area(d)
            moment = disk_moment(d)
            if area_disk <= 0.0:
                return 0.0, centre_y, offset
            # centroid in disk coordinates: distance |moment/area| along the cut normal
            scale = moment / area_disk
            u = scale * alpha / norm
            v = scale * beta / norm
            delta_y, delta_z = a * u, b * v
            world_y = centre_y + delta_y * cos_heel + delta_z * sin_heel
            world_z = offset - delta_y * sin_heel + delta_z * cos_heel
            return a * b * area_disk, world_y, world_z

        return area_of_level, rotated_properties, (level_of(-1.0), level_of(1.0))

    total_area = math.pi * a * b
    upright = a * b * disk_area((draft - b) / b)
    target = upright if displaced_area is None else float(displaced_area)
    if target > total_area:
        raise GzError(
            f"the displaced area {target!r} m^2 exceeds the section area {total_area!r} m^2"
        )
    return _righting_arms(
        factory,
        angles_deg=angles_deg,
        kg=kg,
        displaced_area=target,
        length=length,
        section=f"ellipse[a={a},b={b},draft={draft}]",
    )


def regular_polygon_section(sides: int, radius: float, *, rotation_deg: float = 0.0) -> tuple[Point2, ...]:
    """A regular polygon section with its bottom vertex on ``z = 0`` (a keel at the origin).

    ``rotation_deg`` turns the polygon about its own centre before the keel is taken, which is
    how a hard-chine section is *not* aligned with the waterplane -- the case where a
    symmetry-aligned reference would pass vacuously.
    """
    if sides < 3:
        raise GzError(f"a polygon needs at least 3 sides, got {sides}")
    if not (radius > 0.0):
        raise GzError(f"the radius must be positive: {radius!r}")
    turn = math.radians(rotation_deg)
    points = [
        (
            radius * math.cos(2.0 * math.pi * index / sides - 0.5 * math.pi + turn),
            radius * math.sin(2.0 * math.pi * index / sides - 0.5 * math.pi + turn),
        )
        for index in range(sides)
    ]
    keel = min(z for (_, z) in points)
    return tuple((y, z - keel) for (y, z) in points)


def ellipse_section(semi_beam: float, semi_height: float, sides: int, *, rotation_deg: float = 0.0) -> tuple[Point2, ...]:
    """``sides`` points on the ellipse of :func:`ellipse_prism_righting_arms`, keel at ``z = 0``.

    An *inscribed* polygon, so the mesh built from it under-fills the exact section and the
    chord error has a known sign.
    """
    if sides < 3:
        raise GzError(f"a polygon needs at least 3 sides, got {sides}")
    turn = math.radians(rotation_deg)
    points = [
        (
            semi_beam * math.cos(2.0 * math.pi * index / sides + turn),
            semi_height + semi_height * math.sin(2.0 * math.pi * index / sides + turn),
        )
        for index in range(sides)
    ]
    keel = min(z for (_, z) in points)
    return tuple((y, z - keel) for (y, z) in points)


def prism_mesh(section: Sequence[Point2], length: float, *, x0: float = 0.0) -> Mesh:
    """Extrude a convex section along ``+x`` into a closed, outward-oriented prism mesh.

    The mesh is what :func:`tensorcad.gz.righting_arm_curve` consumes, so it exists to make
    the comparison between the mesh integrator and the exact section oracle possible at all.
    Winding is constructed rather than repaired: sides are wound so their normals point away
    from the section interior, and the two caps are fanned with the section's own orientation.
    """
    ordered = _validate_section(section)
    if not (length > 0.0):
        raise GzError(f"the prism length must be positive: {length!r}")
    builder = MeshBuilder()
    near = [(x0, y, z) for (y, z) in ordered]
    far = [(x0 + length, y, z) for (y, z) in ordered]
    near_index = [builder.add_vertex(point) for point in near]
    far_index = [builder.add_vertex(point) for point in far]
    count = len(ordered)
    for index in range(count):
        nxt = (index + 1) % count
        a, b = near_index[index], near_index[nxt]
        c, d = far_index[nxt], far_index[index]
        builder.add_triangle(a, b, c)
        builder.add_triangle(a, c, d)
    for index in range(1, count - 1):
        builder.add_triangle(far_index[0], far_index[index], far_index[index + 1])
        builder.add_triangle(near_index[0], near_index[index + 1], near_index[index])
    return builder.build()