"""Large-angle stability: the righting-arm curve ``GZ(phi)`` from a heeled displacement.

``tensorcad.stability`` stops where it says it stops: ``GM`` is the *slope* of
``GZ(phi)`` at ``phi = 0``, a waterplane moment, and it says nothing about the curve
away from upright. Everything a naval architect uses after that -- the maximum righting
arm, the angle at which the arm vanishes (the range of stability), the area under the
curve that a dynamic heel has to spend -- needs a different integral: for each heel the
*solid* displacement is re-equilibrated and its centre of buoyancy found.

That is what this module does. A body heeled by ``phi`` keeps gravity vertical, so the
honest way to write the problem is to leave the water flat and tilt the body:

    submerged region = { p : p_z - tan(phi) * p_y < wl }

in the *body* frame, where ``phi > 0`` puts the starboard side down and ``wl`` is the
waterline's intercept at ``y = 0``. For each heel angle ``wl`` is solved so that the
submerged volume equals the displacement; the centroid ``B`` of that volume then gives

    GZ = y_B cos(phi) + (z_B - z_G) sin(phi)

the horizontal distance from ``G`` to the vertical through ``B``, measured in the upright
transverse direction. At small angles this collapses to ``GZ = GM sin(phi)``, which is
asserted as a test rather than assumed.

**How the volume is measured.** The tilted body is cut by the plane and the submerged
part is treated as a closed polyhedron: the clipped triangles of the wetted surface plus
the planar cut. Volume and first moments are summed as signed tetrahedra from an apex
that lies *in the cut plane* -- so the cut faces, which all lie in that plane, span
zero volume and never have to be triangulated. The cut itself is a per-triangle
Sutherland-Hodgman clip against one scalar field, exact for the polyhedron the mesh
actually is: a box reproduces its closed-form righting arm to round-off, which is the
first claim ``examples/righting_arm`` checks.

**What this is not.** Free-surface corrections (that is a tank property, not a hull
property), damage or flooded-compartment stability, trim change with heel (a cross curve
holds trim constant, and this module holds it constant too), and any statutory criterion
(IS Code / SOLAS intact-stability rules are outside the package). A mesh cut off at its
design waterline has no freeboard above it: past the heel at which the waterline reaches
that rim the answer describes a body that does not exist, so each point carries
``deck_submerged`` and the curve reports ``physical_limit_deg``. Feed a body that
actually models the freeboard (a box with depth, a hull lofted above the waterline) and
that limit moves up with the geometry.

Only the Python standard library is used.
"""

from __future__ import annotations

import itertools
import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass

from . import mesh as meshlib
from .geometry import Mesh

SCHEMA_RIGHTING_ARM = "tensorcad.righting_arm/1"
SCHEMA_RIGHTING_ARM_CURVE = "tensorcad.righting_arm.curve/1"

#: Heel angles measured by :func:`righting_arm_curve` when none are given. The two
#: hundredth-degree points are what makes ``initial_slope`` metacentric-grade: the
#: Richardson extrapolation there removes the ``phi^2`` term and leaves an ``O(phi^4)``
#: residual of order 1e-9 m/rad, where a pair of five-degree points would leave 1e-6.
DEFAULT_HEEL_ANGLES_DEG: tuple[float, ...] = (
    0.0, 0.05, 0.1, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0
)

#: Relative tolerance on the displacement the waterline has to reproduce.
DEFAULT_VOLUME_TOLERANCE = 1e-12


class RightingArmError(ValueError):
    """Raised when a righting arm cannot be measured from the geometry given."""


def _round(value: float, digits: int = 9) -> float:
    if not math.isfinite(value):
        return value
    return round(value, digits)


def _bounding_scale(mesh: Mesh) -> float:
    sizes = mesh.bounding_box().sizes
    scale = max(abs(size) for size in sizes)
    return scale if scale > 0.0 else 1.0


def _clip_below(
    polygon: Sequence[tuple[float, float, float]], tan_heel: float, waterline: float
) -> list[tuple[float, float, float]]:
    """Clip a convex polygon to ``z - tan(phi) * y <= waterline`` (Sutherland-Hodgman).

    One plane, one scalar field, so the clipped polygon is convex and has at most one
    edge more than the input. Interpolated points are computed the same way for both
    triangles that share the cut edge, from the same two vertex coordinates, so the two
    halves of a cut face meet without a gap.
    """

    def height(point: tuple[float, float, float]) -> float:
        return point[2] - tan_heel * point[1] - waterline

    output: list[tuple[float, float, float]] = []
    count = len(polygon)
    for index in range(count):
        start = polygon[index]
        end = polygon[(index + 1) % count]
        h_start = height(start)
        h_end = height(end)
        start_in = h_start <= 0.0
        if start_in:
            output.append(start)
        if start_in != (h_end <= 0.0):
            t = h_start / (h_start - h_end)
            output.append(
                (
                    start[0] + t * (end[0] - start[0]),
                    start[1] + t * (end[1] - start[1]),
                    start[2] + t * (end[2] - start[2]),
                )
            )
    return output


def _submerged_volume_and_centroid(
    mesh: Mesh,
    tan_heel: float,
    waterline: float,
    apex: tuple[float, float, float],
) -> tuple[float, tuple[float, float, float]]:
    """Volume and centroid of the part of ``mesh`` below the tilted waterline.

    Signed tetrahedra are taken from ``apex``, which lies *in* the waterline plane, so
    the planar cut contributes exactly zero and needs no triangulation. The sign of the
    sum is the mesh's outward orientation, so an inward-facing mesh shows up as a
    negative volume instead of a plausible-looking answer.
    """
    vertices = mesh.vertices
    ax, ay, az = apex
    volumes: list[float] = []
    moments_x: list[float] = []
    moments_y: list[float] = []
    moments_z: list[float] = []

    for face in mesh.faces:
        lower, middle, upper = (vertices[index] for index in face)
        h_lower = lower[2] - tan_heel * lower[1] - waterline
        h_middle = middle[2] - tan_heel * middle[1] - waterline
        h_upper = upper[2] - tan_heel * upper[1] - waterline
        if h_lower > 0.0 or h_middle > 0.0 or h_upper > 0.0:
            # The fast path is the whole face below the plane, which is nearly every
            # face of a hull at a moderate heel: skip the Sutherland-Hodgman machinery.
            clipped = _clip_below([lower, middle, upper], tan_heel, waterline)
            if len(clipped) < 3:
                continue
        else:
            clipped = [lower, middle, upper]
        first = clipped[0]
        for k in range(1, len(clipped) - 1):
            b = clipped[k]
            c = clipped[k + 1]
            ux, uy, uz = first[0] - ax, first[1] - ay, first[2] - az
            vx, vy, vz = b[0] - ax, b[1] - ay, b[2] - az
            wx, wy, wz = c[0] - ax, c[1] - ay, c[2] - az
            cross_x = vy * wz - vz * wy
            cross_y = vz * wx - vx * wz
            cross_z = vx * wy - vy * wx
            tetra = (ux * cross_x + uy * cross_y + uz * cross_z) / 6.0
            if tetra == 0.0:
                continue
            volumes.append(tetra)
            quarter = 0.25
            moments_x.append(tetra * quarter * (ax + first[0] + b[0] + c[0]))
            moments_y.append(tetra * quarter * (ay + first[1] + b[1] + c[1]))
            moments_z.append(tetra * quarter * (az + first[2] + b[2] + c[2]))

    volume = math.fsum(volumes)
    if volume <= 0.0 or not math.isfinite(volume):
        return volume, (float("nan"), float("nan"), float("nan"))
    return volume, (
        math.fsum(moments_x) / volume,
        math.fsum(moments_y) / volume,
        math.fsum(moments_z) / volume,
    )


def _apex_for(mesh: Mesh, waterline: float) -> tuple[float, float, float]:
    """A point in the waterline plane, inside the body's horizontal extent.

    ``y = 0`` is the centreline for every geometry this package generates, and the
    mid-length abscissa keeps the apex inside the section so the tetrahedra are not
    thin and the signed sum does not lose digits to cancellation.
    """
    bounds = mesh.bounding_box()
    (x_min, _, _), (x_max, _, _) = (bounds.minimum, bounds.maximum)
    return (0.5 * (x_min + x_max), 0.0, waterline)


def heeled_displacement(
    mesh: Mesh,
    heel_angle_deg: float,
    waterline: float,
) -> tuple[float, tuple[float, float, float]]:
    """Submerged volume and centre of buoyancy for an explicit (heel, waterline).

    ``waterline`` is the intercept of the tilted plane at ``y = 0`` in the body frame.
    This is the root-finding kernel exposed for callers that already know the waterline
    (and for the tests that compare a solved waterline against a closed form).
    """
    if not math.isfinite(heel_angle_deg):
        raise RightingArmError(f"heel angle must be finite, got {heel_angle_deg!r}")
    if not math.isfinite(waterline):
        raise RightingArmError(f"waterline must be finite, got {waterline!r}")
    tan_heel = math.tan(math.radians(heel_angle_deg))
    return _submerged_volume_and_centroid(mesh, tan_heel, waterline, _apex_for(mesh, waterline))


@dataclass(frozen=True)
class HeeledEquilibrium:
    """The waterline that makes a heeled body displace the required volume."""

    heel_angle_deg: float
    waterline_z: float
    displacement_volume: float
    centre_of_buoyancy: tuple[float, float, float]
    iterations: int
    residual: float
    basis: str = "clipped-polyhedron volume, signed tetrahedra"

    def to_dict(self) -> dict[str, object]:
        return {
            "heel_angle_deg": _round(self.heel_angle_deg),
            "waterline_z_m": _round(self.waterline_z),
            "displacement_volume_m3": _round(self.displacement_volume),
            "centre_of_buoyancy_m": [_round(value) for value in self.centre_of_buoyancy],
            "iterations": self.iterations,
            "residual_m3": _round(self.residual),
            "basis": self.basis,
        }


def _waterline_bracket(mesh: Mesh, tan_heel: float) -> tuple[float, float, float]:
    """The waterline range that leaves the body empty, and the one that submerges it.

    Both bounds are exact functions of the *vertices*, not of a bisection estimate:
    ``h = z - tan(phi) y - wl`` attains its minimum and maximum at vertices, so
    ``wl_min = min(h + wl)`` is the last waterline at which nothing is submerged and
    ``wl_max = max(h + wl)`` is the first at which everything is. Bracketing on those
    two numbers is what keeps the solver away from the degenerate solution below.
    """
    heights = [vertex[2] - tan_heel * vertex[1] for vertex in mesh.vertices]
    empty = min(heights)
    full = max(heights)
    return empty, full, full - empty


def _solve_waterline(
    mesh: Mesh,
    tan_heel: float,
    target_volume: float,
    *,
    tolerance: float = DEFAULT_VOLUME_TOLERANCE,
) -> tuple[float, float, tuple[float, float, float], int]:
    """Bisect-and-secant the waterline that displaces ``target_volume``.

    ``volume(waterline)`` increases monotonically up to the waterline at which the body
    is fully submerged and is *constant* above it, so the equation has a whole flat
    interval of solutions there. The floating equilibrium is the least of them, and the
    bracket is the exact empty-to-full vertex range, which makes the root unique and
    strictly bracketed. That matters: a mesh capped at its design waterline encloses
    exactly the upright displacement, and accepting any root in the flat region returns
    a fully submerged body whose centre of buoyancy cannot move -- a plausible-looking
    ``GZ`` curve with no hydrostatics in it. Illinois regula falsi guides the search and
    bisection takes over whenever the secant leaves the bracket or stops halving it, so
    the flat volume of a box (kinks at the bilge) cannot stall it either.
    """
    scale = _bounding_scale(mesh)
    empty, full, span = _waterline_bracket(mesh, tan_heel)
    if span <= 0.0:  # pragma: no cover - a flat body along the heel axis
        raise RightingArmError("the model has no height in the heeled frame")

    volume_full, _ = _submerged_volume_and_centroid(mesh, tan_heel, full, _apex_for(mesh, full))
    if volume_full <= 0.0:
        raise RightingArmError(
            "the model's outward orientation is inverted (the submerged volume is negative)"
        )
    if target_volume >= volume_full * (1.0 - 1e-13):
        raise RightingArmError(
            f"a body enclosing {volume_full!r} m^3 cannot float displacing "
            f"{target_volume!r} m^3 at heel {math.degrees(math.atan(tan_heel)):.6f} deg: "
            f"it would have to be fully submerged, where its centre of buoyancy cannot move"
        )

    low, f_low = empty, -target_volume
    high, f_high = full, volume_full - target_volume
    side = 0
    iterations = 0
    waterline = 0.5 * (low + high)
    centre_best: tuple[float, float, float] = (float("nan"), float("nan"), float("nan"))
    f_best = f_low
    while iterations < 400:
        iterations += 1
        if f_high - f_low > 0.0:
            candidate = (low * f_high - high * f_low) / (f_high - f_low)
        else:  # pragma: no cover - equal residuals collapse to bisection
            candidate = 0.5 * (low + high)
        if not (low < candidate < high) or 4.0 * abs(candidate - waterline) > (high - low):
            candidate = 0.5 * (low + high)
        waterline = candidate
        volume, centre = _submerged_volume_and_centroid(
            mesh, tan_heel, waterline, _apex_for(mesh, waterline)
        )
        f_value = volume - target_volume
        f_best = f_value
        centre_best = centre
        if f_value < 0.0:
            low, f_low = waterline, f_value
            if side == -1:
                f_high *= 0.5
            side = -1
        else:
            high, f_high = waterline, f_value
            if side == 1:
                f_low *= 0.5
            side = 1
        if abs(f_value) <= tolerance * target_volume or high - low <= 1e-14 * max(scale, 1.0):
            break
    return waterline, f_best, centre_best, iterations


def heeled_equilibrium(
    mesh: Mesh,
    heel_angle_deg: float,
    displacement_volume: float,
    *,
    tolerance: float = DEFAULT_VOLUME_TOLERANCE,
) -> HeeledEquilibrium:
    """Re-equilibrate a body heeled by ``heel_angle_deg`` to a fixed displacement.

    Trim is held constant (that is what a cross curve is), so the only unknown is the
    waterline's height ``wl``.
    """
    if not math.isfinite(displacement_volume) or displacement_volume <= 0.0:
        raise RightingArmError(
            f"displacement volume must be finite and positive, got {displacement_volume!r}"
        )
    if not math.isfinite(heel_angle_deg):
        raise RightingArmError(f"heel angle must be finite, got {heel_angle_deg!r}")
    tan_heel = math.tan(math.radians(heel_angle_deg))
    waterline, residual, centre, iterations = _solve_waterline(
        mesh, tan_heel, displacement_volume, tolerance=tolerance
    )
    return HeeledEquilibrium(
        heel_angle_deg=heel_angle_deg,
        waterline_z=waterline,
        displacement_volume=displacement_volume,
        centre_of_buoyancy=centre,
        iterations=iterations,
        residual=residual,
    )


def _keel_z(mesh: Mesh) -> float:
    return min(vertex[2] for vertex in mesh.vertices)


def _deck_edge_wetted(mesh: Mesh, tan_heel: float, waterline: float) -> bool:
    """True when the tilted waterline rises above the model's own top rim.

    A hull generated at its design waterline stops there; the water above it is not
    modelled, so a heel that puts the rim under water is answered by a shape the caller
    did not supply. The test is exact for the polyhedral rim: a rim vertex below the
    plane means part of the (artificial) top cap is now wetted surface.
    """
    top = max(vertex[2] for vertex in mesh.vertices)
    scale = _bounding_scale(mesh)
    tolerance = 1e-12 * scale
    for vertex in mesh.vertices:
        if vertex[2] >= top - tolerance and vertex[2] - tan_heel * vertex[1] - waterline < -tolerance:
            return True
    return False


@dataclass(frozen=True)
class RightingArm:
    """One point of the righting-arm curve: a heeled, re-equilibrated floating body."""

    heel_angle_deg: float
    gz: float
    waterline_z: float
    centre_of_buoyancy: tuple[float, float, float]
    centre_of_gravity_z: float
    displacement_volume: float
    deck_submerged: bool
    iterations: int
    residual: float
    basis: str = "clipped-polyhedron volume, signed tetrahedra"

    @property
    def centre_of_gravity_offset(self) -> float:
        """``z_G`` above the keel, the datum the arm is measured from."""
        return self.centre_of_gravity_z

    def to_dict(self) -> dict[str, object]:
        return {
            "heel_angle_deg": _round(self.heel_angle_deg),
            "gz_m": _round(self.gz),
            "waterline_z_m": _round(self.waterline_z),
            "centre_of_buoyancy_m": [_round(value) for value in self.centre_of_buoyancy],
            "centre_of_gravity_z_m": _round(self.centre_of_gravity_z),
            "displacement_volume_m3": _round(self.displacement_volume),
            "deck_submerged": self.deck_submerged,
            "iterations": self.iterations,
            "residual_m3": _round(self.residual),
            "basis": self.basis,
        }


def righting_arm(
    mesh: Mesh,
    *,
    heel_angle_deg: float,
    kg: float,
    displacement_volume: float,
    gy: float = 0.0,
    tolerance: float = DEFAULT_VOLUME_TOLERANCE,
) -> RightingArm:
    """The righting arm at one heel angle.

    ``kg`` is the height of the centre of gravity above the keel (the model's lowest
    vertex), matching :func:`tensorcad.stability.initial_stability`; ``gy`` is its
    transverse offset from the centreline and is a required part of the answer for a
    list to port, so it is exposed even though every geometry this package generates has
    ``gy = 0``. ``displacement_volume`` is the volume the body has to displace, which
    fixes the sinkage at every heel.
    """
    if not math.isfinite(kg):
        raise RightingArmError(f"kg must be finite, got {kg!r}")
    if not math.isfinite(gy):
        raise RightingArmError(f"gy must be finite, got {gy!r}")
    equilibrium = heeled_equilibrium(
        mesh, heel_angle_deg, displacement_volume, tolerance=tolerance
    )
    keel = _keel_z(mesh)
    z_g = keel + kg
    heel = math.radians(heel_angle_deg)
    cos_heel, sin_heel = math.cos(heel), math.sin(heel)
    _, y_b, z_b = equilibrium.centre_of_buoyancy
    gz = (y_b - gy) * cos_heel + (z_b - z_g) * sin_heel
    tan_heel = math.tan(heel)
    return RightingArm(
        heel_angle_deg=heel_angle_deg,
        gz=gz,
        waterline_z=equilibrium.waterline_z,
        centre_of_buoyancy=equilibrium.centre_of_buoyancy,
        centre_of_gravity_z=z_g,
        displacement_volume=displacement_volume,
        deck_submerged=_deck_edge_wetted(mesh, tan_heel, equilibrium.waterline_z),
        iterations=equilibrium.iterations,
        residual=equilibrium.residual,
    )


@dataclass(frozen=True)
class RightingArmCurve:
    """A righting-arm curve plus the three numbers read off it."""

    points: tuple[RightingArm, ...]
    kg: float
    gy: float
    displacement_volume: float
    keel_z: float
    maximum_gz: float
    angle_of_maximum_gz_deg: float
    vanishing_angle_deg: float | None
    zero_crossings_deg: tuple[float, ...]
    physical_limit_deg: float | None
    initial_slope: float
    basis: str = "heeled displacement integration at constant trim"

    def at(self, heel_angle_deg: float) -> RightingArm:
        for point in self.points:
            if math.isclose(point.heel_angle_deg, heel_angle_deg, rel_tol=0.0, abs_tol=1e-12):
                return point
        raise KeyError(f"heel angle {heel_angle_deg!r} deg is not on this curve")

    def dynamical_stability(self, upper_deg: float | None = None) -> float:
        """Area under the curve (the dynamic heel it can absorb), trapezoidal on the points.

        Restricted to ``upper_deg`` when given; the trapezoid rule is the definition of
        the reported number, not an approximation of a smooth curve the caller supplied.
        """
        limit = math.inf if upper_deg is None else upper_deg
        area = 0.0
        for left, right in itertools.pairwise(self.points):
            if left.heel_angle_deg >= limit:
                break
            span = min(right.heel_angle_deg, limit) - left.heel_angle_deg
            if span <= 0.0:
                continue
            fraction = span / (right.heel_angle_deg - left.heel_angle_deg)
            right_gz = left.gz + fraction * (right.gz - left.gz)
            area += 0.5 * (left.gz + right_gz) * math.radians(span)
        return area

    def to_dict(self) -> dict[str, object]:
        return {
            "schema": SCHEMA_RIGHTING_ARM_CURVE,
            "kg_m_above_keel": _round(self.kg),
            "gy_m": _round(self.gy),
            "keel_z_m": _round(self.keel_z),
            "displacement_volume_m3": _round(self.displacement_volume),
            "maximum_gz_m": _round(self.maximum_gz),
            "angle_of_maximum_gz_deg": _round(self.angle_of_maximum_gz_deg),
            "vanishing_angle_deg": (
                None if self.vanishing_angle_deg is None else _round(self.vanishing_angle_deg)
            ),
            "zero_crossings_deg": [_round(angle) for angle in self.zero_crossings_deg],
            "physical_limit_deg": (
                None if self.physical_limit_deg is None else _round(self.physical_limit_deg)
            ),
            "initial_slope_m_per_rad": _round(self.initial_slope),
            "basis": self.basis,
            "points": [point.to_dict() for point in self.points],
        }


def righting_arm_curve(
    mesh: Mesh,
    *,
    kg: float,
    displacement_volume: float,
    angles_deg: Sequence[float] = DEFAULT_HEEL_ANGLES_DEG,
    gy: float = 0.0,
    tolerance: float = DEFAULT_VOLUME_TOLERANCE,
) -> RightingArmCurve:
    """Sweep the heel angles and read the curve's maxima and zeros off the result.

    ``vanishing_angle_deg`` is the first heel past the maximum at which the arm changes
    sign (the range of stability); it is ``None`` when the curve stays positive over the
    angles asked for, which is an honest "not in the range scanned", not "infinite".
    ``physical_limit_deg`` is the first heel at which the model's own top rim goes under
    water, reported because past it the curve describes a body with no freeboard.
    """
    if not angles_deg:
        raise RightingArmError("at least one heel angle is required")
    angles = sorted(float(angle) for angle in angles_deg)
    points = tuple(
        righting_arm(
            mesh,
            heel_angle_deg=angle,
            kg=kg,
            displacement_volume=displacement_volume,
            gy=gy,
            tolerance=tolerance,
        )
        for angle in angles
    )
    maximum = max(points, key=lambda point: point.gz)
    crossings: list[float] = []
    for left, right in itertools.pairwise(points):
        if left.gz * right.gz < 0.0 or (left.gz != 0.0 and right.gz == 0.0):
            crossings.append(
                _refine_zero(
                    mesh,
                    left.heel_angle_deg,
                    right.heel_angle_deg,
                    kg=kg,
                    gy=gy,
                    displacement_volume=displacement_volume,
                    tolerance=tolerance,
                )
            )
    # The range of stability ends where a positive arm first goes non-positive; a curve
    # that starts non-positive (an unstable upright body, GM < 0) has no such crossing --
    # its zero is the angle of loll, which is reported in ``zero_crossings_deg`` but is
    # not a vanishing angle and is not silently relabelled as one.
    vanishing = next(
        (angle for angle in crossings if angle > maximum.heel_angle_deg - 1e-12),
        None,
    ) if maximum.gz > 0.0 else None
    physical_limit = _physical_limit(mesh, points, kg=kg, gy=gy, displacement_volume=displacement_volume,
                                     tolerance=tolerance)
    initial_slope = _initial_slope(points)
    return RightingArmCurve(
        points=points,
        kg=kg,
        gy=gy,
        displacement_volume=displacement_volume,
        keel_z=_keel_z(mesh),
        maximum_gz=maximum.gz,
        angle_of_maximum_gz_deg=maximum.heel_angle_deg,
        vanishing_angle_deg=vanishing,
        zero_crossings_deg=tuple(crossings),
        physical_limit_deg=physical_limit,
        initial_slope=initial_slope,
    )


def _refine_root(
    function: Callable[[float], float],
    low_deg: float,
    high_deg: float,
    *,
    tolerance_deg: float = 1e-7,
    max_iterations: int = 60,
) -> float:
    """Bisect-and-secant a scalar root between two bracketing heel angles.

    The reported landmarks are roots of the curve, not grid points: interpolating on a
    five-degree grid would put a zero several hundredths of a degree off and make a
    comparison with a closed form look like a method error when it is a sampling error.
    Illinois keeps the secant fast without letting one end of the bracket pin it, which
    matters because every evaluation is a full heeled equilibrium.
    """
    low, high = low_deg, high_deg
    f_low, f_high = function(low), function(high)
    if f_low == 0.0:
        return low
    if f_high == 0.0:
        return high
    if (f_low < 0.0) == (f_high < 0.0):  # pragma: no cover - the caller brackets the root
        raise RightingArmError(
            f"no sign change on [{low_deg!r}, {high_deg!r}] deg ({f_low!r}, {f_high!r})"
        )
    side = 0
    for _ in range(max_iterations):
        if high - low <= tolerance_deg:
            break
        candidate = (low * f_high - high * f_low) / (f_high - f_low)
        if not (low < candidate < high):
            candidate = 0.5 * (low + high)
        value = function(candidate)
        if value == 0.0:
            return candidate
        if (value < 0.0) == (f_low < 0.0):
            low, f_low = candidate, value
            if side == -1:
                f_high *= 0.5
            side = -1
        else:
            high, f_high = candidate, value
            if side == 1:
                f_low *= 0.5
            side = 1
    return 0.5 * (low + high)


def _refine_zero(
    mesh: Mesh,
    low_deg: float,
    high_deg: float,
    *,
    kg: float,
    gy: float,
    displacement_volume: float,
    tolerance: float,
) -> float:
    """The heel at which the arm changes sign, inside a bracketing pair."""

    def arm(angle: float) -> float:
        return righting_arm(
            mesh,
            heel_angle_deg=angle,
            kg=kg,
            gy=gy,
            displacement_volume=displacement_volume,
            tolerance=tolerance,
        ).gz

    return _refine_root(arm, low_deg, high_deg)


def _physical_limit(
    mesh: Mesh,
    points: Sequence[RightingArm],
    *,
    kg: float,
    gy: float,
    displacement_volume: float,
    tolerance: float,
) -> float | None:
    """First heel at which the model's own top rim goes under water, bisected.

    ``None`` means the rim stayed clear at every angle scanned -- a statement about the
    range asked for, not about the geometry.
    """
    clean: RightingArm | None = None
    wet: RightingArm | None = None
    for point in points:
        if not point.deck_submerged:
            clean = point
        elif clean is not None:
            wet = point
            break
    if clean is None:
        return points[0].heel_angle_deg if points and points[0].deck_submerged else None
    if wet is None:
        return None

    def clearance(angle: float) -> float:
        point = righting_arm(
            mesh,
            heel_angle_deg=angle,
            kg=kg,
            gy=gy,
            displacement_volume=displacement_volume,
            tolerance=tolerance,
        )
        return -1.0 if point.deck_submerged else 1.0

    return _refine_root(clearance, clean.heel_angle_deg, wet.heel_angle_deg)


def _initial_slope(points: Sequence[RightingArm]) -> float:
    """``GZ / sin(phi)`` in the limit ``phi -> 0``, by two-point Richardson extrapolation.

    The arm is odd in the heel (``GZ = GM sin(phi) + O(phi^3)``), so the ratio
    ``GZ / sin(phi)`` at two small non-zero heels extrapolates to the slope at the
    origin with the ``phi^2`` term removed. Using the widest pair of the smallest two
    positive angles keeps the extrapolation away from the ``0 -> small`` round-off of a
    directly solved ``GZ(0) = 0`` point.
    """
    positive = [point for point in points if point.heel_angle_deg > 0.0 and abs(point.gz) > 0.0]
    positive.sort(key=lambda point: point.heel_angle_deg)
    if len(positive) >= 2:
        coarse, fine = positive[1], positive[0]
        ratio_coarse = coarse.gz / math.sin(math.radians(coarse.heel_angle_deg))
        ratio_fine = fine.gz / math.sin(math.radians(fine.heel_angle_deg))
        observed = (coarse.heel_angle_deg / fine.heel_angle_deg) ** 2
        if observed > 1.0 + 1e-6:
            return (ratio_fine * observed - ratio_coarse) / (observed - 1.0)
        return ratio_fine
    if positive:
        return positive[0].gz / math.sin(math.radians(positive[0].heel_angle_deg))
    if len(points) >= 2:  # pragma: no cover - needs a curve with no positive heel at all
        raise RightingArmError("cannot estimate the initial slope without a positive heel angle")
    raise RightingArmError("a curve needs at least two heel angles to have a slope")


def wall_sided_righting_arm(
    heel_angle_deg: float,
    *,
    beam: float,
    depth: float,
    draft: float,
    kg: float,
    gy: float = 0.0,
) -> float:
    """Closed-form righting arm of a rectangular, wall-sided section.

        GZ = sin(phi) * (GM + 0.5 * BM * tan(phi)^2),      GM = KB + BM - KG

    with ``KB = draft / 2`` and ``BM = beam^2 / (12 * draft)``. The closed form holds
    only while the waterline still meets both vertical sides, that is while
    ``tan(phi) < 2 * min(draft, depth - draft) / beam``: past that either the deck edge
    is under water or the bilge has come out, and the section is no longer wall-sided.
    The function raises rather than returning a formula its assumptions do not cover,
    which is what lets the tests compare it against an integration at full precision.
    """
    if beam <= 0.0 or depth <= 0.0 or draft <= 0.0:
        raise RightingArmError("beam, depth and draft must be positive")
    if not 0.0 < draft < depth:
        raise RightingArmError(f"draft {draft!r} must lie between 0 and depth {depth!r}")
    tan_heel = math.tan(math.radians(heel_angle_deg))
    limit = 2.0 * min(draft, depth - draft) / beam
    if abs(tan_heel) >= limit:
        raise RightingArmError(
            f"heel {heel_angle_deg!r} deg is outside the wall-sided range (|tan| < {limit!r})"
        )
    kb = 0.5 * draft
    bm = beam * beam / (12.0 * draft)
    gm = kb + bm - kg
    return math.sin(math.radians(heel_angle_deg)) * (
        gm + 0.5 * bm * tan_heel * tan_heel
    ) + gy * math.cos(math.radians(heel_angle_deg))


def wall_sided_vanishing_angle(*, beam: float, depth: float, draft: float, kg: float) -> float | None:
    """The heel at which the wall-sided arm vanishes, when it vanishes inside the range.

    ``sin(phi) (GM + BM tan^2(phi) / 2) = 0`` has the non-trivial root
    ``tan(phi) = sqrt(-2 GM / BM)``, which exists only for ``GM < 0`` -- an unstable
    upright body -- and is returned only when that root is still wall-sided.
    """
    kb = 0.5 * draft
    bm = beam * beam / (12.0 * draft)
    gm = kb + bm - kg
    if gm >= 0.0:
        return None
    tan_vanishing = math.sqrt(-2.0 * gm / bm)
    limit = 2.0 * min(draft, depth - draft) / beam
    if tan_vanishing >= limit:
        return None
    return math.degrees(math.atan(tan_vanishing))


def is_watertight(mesh: Mesh) -> bool:
    """Re-export the watertightness test the volume integral depends on."""
    return meshlib.is_watertight(mesh)
