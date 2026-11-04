"""Large-angle stability: the righting arm ``GZ(phi)`` and the cross curves ``KN``.

``tensorcad.stability`` answers the first question a naval architect asks of a hull -- how
far can G be raised before the upright position stops being stable -- but it answers it
*at* the upright position: ``GM`` is the slope of ``GZ(phi)`` at ``phi = 0``, and a ship
whose ``GM`` is positive can still capsize, because the righting arm grows, peaks and then
vanishes again. This module computes the curve, not just its first derivative.

**What ``GZ`` is.** Heel the ship to starboard by ``phi`` and let it find the waterline
that displaces the loading condition's volume. The buoyancy force acts upwards through
``B``, the weight downwards through ``G``, and the righting arm is the horizontal distance
between those two parallel lines,

    GZ(phi) = cos(phi) * (y_B - y_G) + sin(phi) * (z_B - z_G),

with ``(y_B, z_B)`` the centroid of the displaced volume **in the heeled body frame** and
``(y_G, z_G)`` the centre of gravity. Two consequences are worth stating because they are
what the tests use:

* ``GZ(0) = y_B - y_G``, which is ``0`` for a hull symmetric about its centreline -- the
  upright righting arm is zero, and all upright stability lives in the slope.
* ``dGZ/dphi (0) = GM``. :func:`gz_initial_slope` measures that derivative and the tests
  compare it against ``tensorcad.stability.initial_stability``, which gets ``GM`` from two
  *different* integrals (``KB`` from the solid, ``BM`` from the waterplane moments).
  Agreement is therefore not a tautology.

**How the submerged volume is measured.** Not by sampling and not by a strip sum: the
closed surface is clipped against the heeled waterplane triangle by triangle, and the
volume and its first moment are read off the clipped surface with the divergence theorem,

    V      = (1/3) * closed_integral (p . n) dS,
    V * c  = (1/4) * closed_integral p (p . n) dS,

evaluated with the origin placed **on** the cutting plane. That choice is the whole trick:
the cap that closes the clipped solid then has ``p . n = 0`` over its whole area, so it
contributes exactly nothing to either integral and never has to be assembled into polygon
loops. The only error left is the one already in the mesh -- a curved hull's sides are
chords -- and ``examples/heeling`` measures it on a refinement ladder instead of assuming
it is small.

**Constant displacement, not constant waterline.** A ``GZ`` curve is defined at fixed
displacement, so the waterline offset is solved for at every heel angle (monotone
bisection on ``V(offset)``, monotone because a plane sweeping through a solid encloses a
non-decreasing volume). This is also why a hull mesh capped *at* the design waterline
cannot produce a curve: it has no freeboard, every heel angle would need the whole solid
under water, and :func:`waterline_offset` says so instead of returning nonsense.

**The reference the tests lean on.** A wall-sided form -- vertical sides, constant
waterplane breadth, waterline between the bilge and the deck -- has the exact closed form

    GZ(phi) = sin(phi) * (GM + BM * tan^2(phi) / 2),

so a rectangular barge has a ``GZ`` known analytically at *every* heel angle, not just
near zero. :func:`wall_sided_righting_arm` exposes it and the tests require agreement to
round-off.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Sequence

from . import mesh as meshlib
from .geometry import Mesh, MeshBuilder, Vec3
from .stability import waterplane_moments

__all__ = [
    "HeelError",
    "HeeledFloat",
    "RightingArmCurve",
    "barge_mesh",
    "cross_curve",
    "cross_curves",
    "gz_curve",
    "gz_initial_slope",
    "heel_normal",
    "heeled_float",
    "submerged_properties",
    "wall_sided_righting_arm",
    "waterline_offset",
]

Vec3Like = tuple[float, float, float]

# Bisection stops when the bracket on the waterline offset is this fraction of the body's
# extent along the waterline normal. At 1e-14 a box barge's GZ is reproduced to round-off;
# loosening it does not buy measurable speed.
OFFSET_TOLERANCE = 1e-14
MAX_BISECTIONS = 200
# A target volume this close to the whole solid has no heeled equilibrium: the plane has
# left the body and every higher plane displaces the same volume.
FULL_VOLUME_TOLERANCE = 1e-12
# A waterline offset whose displaced volume matches the target this closely is accepted.
# 1e-15 relative is one double-precision step away from the best any float solve can do.
VOLUME_TOLERANCE = 1e-15


class HeelError(ValueError):
    """Raised when a heeled equilibrium does not exist for the requested inputs."""


def _round(value: float, digits: int = 9) -> float:
    if not math.isfinite(value):
        return value
    return round(value, digits)


def _dot(a: Sequence[float], b: Sequence[float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a: Sequence[float], b: Sequence[float]) -> Vec3Like:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def heel_normal(heel_degrees: float) -> Vec3Like:
    """Upward waterline normal in the *heeled body frame*, for a starboard heel.

    Positive heel puts the starboard (``+y``) side down, so the outward normal of the
    water surface, written in body coordinates, tilts towards ``-y``:
    ``n = (0, -sin(phi), cos(phi))``. At ``phi = 90`` the starboard side is straight down
    and ``n = (0, -1, 0)``, as it must be.
    """
    if not math.isfinite(heel_degrees):
        raise HeelError(f"heel angle must be finite, got {heel_degrees!r}")
    if abs(heel_degrees) > 90.0:
        raise HeelError(f"heel angle must lie in [-90, 90] degrees, got {heel_degrees!r}")
    radians = math.radians(heel_degrees)
    return (0.0, -math.sin(radians), math.cos(radians))


def _clip_polygon(
    polygon: Sequence[Vec3Like], normal: Vec3Like
) -> tuple[list[Vec3Like], tuple[Vec3Like, Vec3Like] | None]:
    """Keep the part of ``polygon`` with ``normal . p <= 0``.

    Returns the clipped polygon in order, plus -- when the plane actually cut it -- the
    directed segment the cut contributes to the cap boundary. The direction matters: every
    face hands over ``exit -> entry``, so the segments already carry one consistent
    winding around the cap and its area can be summed without assembling loops.
    """
    kept: list[Vec3Like] = []
    exit_point: Vec3Like | None = None
    entry_point: Vec3Like | None = None
    count = len(polygon)
    for index in range(count):
        a = polygon[index]
        b = polygon[(index + 1) % count]
        sa = _dot(normal, a)
        sb = _dot(normal, b)
        if sa <= 0.0:
            kept.append(a)
        if (sa < 0.0 < sb) or (sa > 0.0 > sb):
            t = sa / (sa - sb)
            crossing = (
                a[0] + t * (b[0] - a[0]),
                a[1] + t * (b[1] - a[1]),
                a[2] + t * (b[2] - a[2]),
            )
            kept.append(crossing)
            if sa < 0.0:  # leaving the submerged half-space
                exit_point = crossing
            else:  # re-entering it
                entry_point = crossing
    segment = None
    if exit_point is not None and entry_point is not None:
        segment = (exit_point, entry_point)
    return kept, segment


def _triangle_moments(triangle: Sequence[Vec3Like]) -> tuple[float, Vec3Like, tuple[Vec3Like, Vec3Like, Vec3Like]]:
    """Area, unit normal and ``S_ij = integral x_i x_j dA`` of one triangle.

    ``S_ij = (A/12) (sum_k v_ki v_kj + (sum_k v_ki)(sum_k v_kj))`` is the exact three-point
    vertex rule for a quadratic integrand over a triangle -- no quadrature, no sampling.
    """
    a, b, c = triangle[0], triangle[1], triangle[2]
    ab = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
    ac = (c[0] - a[0], c[1] - a[1], c[2] - a[2])
    nx, ny, nz = _cross(ab, ac)
    norm = math.sqrt(nx * nx + ny * ny + nz * nz)
    if norm <= 0.0:
        zero: Vec3Like = (0.0, 0.0, 0.0)
        return 0.0, zero, (zero, zero, zero)
    area = 0.5 * norm
    normal: Vec3Like = (nx / norm, ny / norm, nz / norm)
    total = (a[0] + b[0] + c[0], a[1] + b[1] + c[1], a[2] + b[2] + c[2])
    rows: list[Vec3Like] = []
    for i in range(3):
        row = []
        for j in range(3):
            vertex_sum = a[i] * a[j] + b[i] * b[j] + c[i] * c[j] + total[i] * total[j]
            row.append(area * vertex_sum / 12.0)
        rows.append((row[0], row[1], row[2]))
    return area, normal, (rows[0], rows[1], rows[2])


def _require_closed(mesh: Mesh, closed: bool | None) -> None:
    if closed is False:
        return
    if not meshlib.is_watertight(mesh):
        raise HeelError(
            "a heeled displacement needs a closed solid: the surface has boundary edges, "
            "so the volume below the waterplane is not defined"
        )


def submerged_properties(
    mesh: Mesh,
    normal: Vec3Like,
    offset: float,
    *,
    closed: bool | None = None,
) -> tuple[float, Vec3Like, float]:
    """Volume, centroid and waterplane area of ``normal . p <= offset``.

    Exact for the polygon mesh: every triangle is clipped against the plane and the
    resulting closed surface is integrated with the divergence theorem, with the origin
    moved onto the plane so the closing cap drops out of both integrals (``p . n = 0``
    on it). The centroid is returned in the mesh's own coordinates.
    """
    if not math.isfinite(offset):
        raise HeelError(f"waterline offset must be finite, got {offset!r}")
    _require_closed(mesh, closed)
    nx, ny, nz = normal
    norm = math.sqrt(nx * nx + ny * ny + nz * nz)
    if norm <= 0.0:
        raise HeelError("waterline normal must be a non-zero vector")
    unit: Vec3Like = (nx / norm, ny / norm, nz / norm)
    origin: Vec3Like = (unit[0] * offset, unit[1] * offset, unit[2] * offset)

    volume_terms: list[float] = []
    moments: list[Vec3Like] = []
    cap_x = cap_y = cap_z = 0.0

    for face in mesh.faces:
        va, vb, vc = mesh.vertices[face[0]], mesh.vertices[face[1]], mesh.vertices[face[2]]
        polygon = (
            (va[0] - origin[0], va[1] - origin[1], va[2] - origin[2]),
            (vb[0] - origin[0], vb[1] - origin[1], vb[2] - origin[2]),
            (vc[0] - origin[0], vc[1] - origin[1], vc[2] - origin[2]),
        )
        kept, segment = _clip_polygon(polygon, unit)
        if segment is not None:
            seg = _cross(segment[0], segment[1])
            cap_x += seg[0]
            cap_y += seg[1]
            cap_z += seg[2]
        if len(kept) < 3:
            continue
        first = kept[0]
        for index in range(1, len(kept) - 1):
            triangle = (first, kept[index], kept[index + 1])
            area, tri_normal, tensor = _triangle_moments(triangle)
            if area <= 0.0:
                continue
            centroid = (
                (triangle[0][0] + triangle[1][0] + triangle[2][0]) / 3.0,
                (triangle[0][1] + triangle[1][1] + triangle[2][1]) / 3.0,
                (triangle[0][2] + triangle[1][2] + triangle[2][2]) / 3.0,
            )
            # V = (1/3) closed_integral (p . n) dS: p . n is linear, so the centroid rule
            # is exact over the clipped triangle.
            volume_terms.append(area * _dot(centroid, tri_normal) / 3.0)
            # V * c = (1/4) closed_integral p (p . n) dS: p is linear and p . n is linear.
            projected = (
                tensor[0][0] * tri_normal[0] + tensor[0][1] * tri_normal[1] + tensor[0][2] * tri_normal[2],
                tensor[1][0] * tri_normal[0] + tensor[1][1] * tri_normal[1] + tensor[1][2] * tri_normal[2],
                tensor[2][0] * tri_normal[0] + tensor[2][1] * tri_normal[1] + tensor[2][2] * tri_normal[2],
            )
            moments.append((projected[0] / 4.0, projected[1] / 4.0, projected[2] / 4.0))

    volume = math.fsum(volume_terms)
    waterplane_area = 0.5 * math.sqrt(cap_x * cap_x + cap_y * cap_y + cap_z * cap_z)
    if abs(volume) <= 1e-15:
        return 0.0, (float("nan"), float("nan"), float("nan")), waterplane_area
    centroid: Vec3 = (
        origin[0] + math.fsum(moment[0] for moment in moments) / volume,
        origin[1] + math.fsum(moment[1] for moment in moments) / volume,
        origin[2] + math.fsum(moment[2] for moment in moments) / volume,
    )
    return volume, centroid, waterplane_area


def _extent(mesh: Mesh, normal: Vec3Like) -> tuple[float, float]:
    projections = [_dot(normal, vertex) for vertex in mesh.vertices]
    return min(projections), max(projections)


def waterline_offset(
    mesh: Mesh,
    normal: Vec3Like,
    target_volume: float,
    *,
    closed: bool | None = None,
) -> float:
    """Offset whose waterplane displaces exactly ``target_volume``.

    ``V(offset)`` is non-decreasing, so the root is bracketed by the projections of the
    body's own vertices and found by bisection: no initial guess, and no convergence
    failure mode that depends on the body's shape.
    """
    if not math.isfinite(target_volume):
        raise HeelError(f"target volume must be finite, got {target_volume!r}")
    if target_volume <= 0.0:
        raise HeelError(f"target volume must be positive, got {target_volume!r}")
    _require_closed(mesh, closed)
    low, high = _extent(mesh, normal)
    solid_volume = abs(meshlib.total_signed_volume(mesh))
    if target_volume >= solid_volume * (1.0 - FULL_VOLUME_TOLERANCE):
        raise HeelError(
            f"target volume {target_volume!r} m^3 leaves no freeboard on a solid of "
            f"{solid_volume!r} m^3: every heel angle would need the whole body under water. "
            "A GZ curve needs geometry above the waterline, and a mesh capped at the design "
            "waterline has none (float it deeper, or loft the hull to a deck edge)."
        )
    span = high - low
    if span <= 0.0:
        raise HeelError("degenerate body: zero extent along the waterline normal")
    tolerance = OFFSET_TOLERANCE * span
    # Illinois (regula falsi with a halved stale end) rather than plain bisection:
    # V(offset) is smooth and nearly linear over the band that matters, so the secant
    # iterate lands in ~12 evaluations where halving the bracket would take ~47. The
    # bracket itself is never widened, so the method cannot diverge.
    low_volume = 0.0 - target_volume
    high_volume = solid_volume - target_volume
    side = 0
    root: float | None = None
    for _ in range(MAX_BISECTIONS):
        denominator = high_volume - low_volume
        if denominator == 0.0:
            break
        middle = (low * high_volume - high * low_volume) / denominator
        if not (low < middle < high):
            middle = 0.5 * (low + high)
        volume, _, _ = submerged_properties(mesh, normal, middle, closed=closed)
        residual = volume - target_volume
        if abs(residual) <= VOLUME_TOLERANCE * target_volume:
            # The iterate itself is the answer; the bracket may still be wide, so returning
            # its midpoint here would throw away exactly the accuracy that was just won.
            root = middle
            break
        if residual < 0.0:
            low, low_volume = middle, residual
            if side == -1:
                high_volume *= 0.5
            side = -1
        else:
            high, high_volume = middle, residual
            if side == 1:
                low_volume *= 0.5
            side = 1
        if high - low <= tolerance:
            break
    if root is not None:
        return root
    return 0.5 * (low + high)


def _finite(name: str, value: float) -> float:
    if not math.isfinite(value):
        raise HeelError(f"{name} must be finite, got {value!r}")
    return value


@dataclass(frozen=True)
class HeeledFloat:
    """One heel angle at a fixed displacement: the geometry behind a ``GZ`` value."""

    heel_degrees: float
    heel_radians: float
    displacement_volume: float
    waterline_offset: float
    centre_of_buoyancy: Vec3
    waterplane_area: float
    kn: float
    keel_z: float = 0.0

    def to_dict(self) -> dict[str, Any]:
        return {
            "heel_degrees": _round(self.heel_degrees, 9),
            "keel_z": _round(self.keel_z, 12),
            "heel_radians": _round(self.heel_radians, 12),
            "displacement_volume": _round(self.displacement_volume, 12),
            "waterline_offset": _round(self.waterline_offset, 12),
            "centre_of_buoyancy": [_round(component, 9) for component in self.centre_of_buoyancy],
            "waterplane_area": _round(self.waterplane_area, 9),
            "kn": _round(self.kn, 12),
        }


def heeled_float(
    mesh: Mesh,
    heel_degrees: float,
    *,
    displacement_volume: float,
    closed: bool | None = None,
) -> HeeledFloat:
    """Equilibrium geometry at ``heel_degrees`` for a given displaced volume."""
    normal = heel_normal(heel_degrees)
    offset = waterline_offset(mesh, normal, displacement_volume, closed=closed)
    volume, centroid, area = submerged_properties(mesh, normal, offset, closed=closed)
    radians = math.radians(heel_degrees)
    kn = math.cos(radians) * centroid[1] + math.sin(radians) * centroid[2]
    return HeeledFloat(
        heel_degrees=heel_degrees,
        heel_radians=radians,
        displacement_volume=volume,
        waterline_offset=offset,
        centre_of_buoyancy=centroid,
        waterplane_area=area,
        kn=kn,
        keel_z=min(vertex[2] for vertex in mesh.vertices),
    )


def wall_sided_righting_arm(heel_radians: float, *, kb: float, bm: float, kg: float) -> float:
    """Exact ``GZ`` for a wall-sided form: ``sin(phi) (GM + BM tan^2(phi) / 2)``.

    Valid only while the waterplane stays on both vertical sides -- above the bilge and
    below the deck. For a rectangular barge that is the range in which the closed form
    holds *exactly*, which is why the tests use it as the reference: it is not an
    approximation of a barge, it is the barge.
    """
    gm = kb + bm - kg
    return math.sin(heel_radians) * (gm + 0.5 * bm * math.tan(heel_radians) ** 2)


def barge_mesh(length: float, beam: float, depth: float) -> Mesh:
    """A wall-sided reference body: an axis-aligned parallelepiped, keel on ``z = 0``.

    Twelve triangles, no discretisation of anything, oriented outward with the bottom
    face wound so that ``total_signed_volume`` is positive. It exists because the
    wall-sided closed form above is *exact* for it, which makes it the one body on which
    ``GZ`` can be checked at every heel angle rather than only near upright.
    """
    for name, value in (("length", length), ("beam", beam), ("depth", depth)):
        if not math.isfinite(value) or value <= 0.0:
            raise HeelError(f"{name} must be a finite positive number, got {value!r}")
    half_length = 0.5 * length
    half_beam = 0.5 * beam
    builder = MeshBuilder()
    corners = {
        "000": builder.add_vertex((-half_length, -half_beam, 0.0)),
        "100": builder.add_vertex((half_length, -half_beam, 0.0)),
        "110": builder.add_vertex((half_length, half_beam, 0.0)),
        "010": builder.add_vertex((-half_length, half_beam, 0.0)),
        "001": builder.add_vertex((-half_length, -half_beam, depth)),
        "101": builder.add_vertex((half_length, -half_beam, depth)),
        "111": builder.add_vertex((half_length, half_beam, depth)),
        "011": builder.add_vertex((-half_length, half_beam, depth)),
    }
    faces = (
        ("000", "010", "110", "100"),  # bottom, normal -z
        ("001", "101", "111", "011"),  # deck, normal +z
        ("000", "100", "101", "001"),  # -y side
        ("010", "011", "111", "110"),  # +y side
        ("000", "001", "011", "010"),  # -x transom
        ("100", "110", "111", "101"),  # +x bow
    )
    for a, b, c, d in faces:
        builder.add_quad(corners[a], corners[b], corners[c], corners[d])
    return builder.build()


@dataclass(frozen=True)
class RightingArmCurve:
    """``GZ(phi)`` at constant displacement, with the ``KN`` decomposition kept visible.

    ``GZ = KN - KG sin(phi) - y_G cos(phi)`` holds pointwise, which is the identity that
    makes cross curves useful: ``KN`` depends on the hull and the displacement alone, so
    one table serves every loading condition.
    """

    kg: float
    yg: float
    displacement_volume: float
    points: tuple[HeeledFloat, ...]
    keel_z: float = 0.0

    @property
    def angles(self) -> tuple[float, ...]:
        return tuple(point.heel_degrees for point in self.points)

    @property
    def kn(self) -> tuple[float, ...]:
        return tuple(point.kn for point in self.points)

    @property
    def gravity_z(self) -> float:
        """``G`` in the mesh's own z coordinate: ``kg`` is a height above the keel."""
        return self.keel_z + self.kg

    @property
    def gz(self) -> tuple[float, ...]:
        return tuple(self._gz(point) for point in self.points)

    def _gz(self, point: HeeledFloat) -> float:
        return (
            point.kn
            - self.gravity_z * math.sin(point.heel_radians)
            - self.yg * math.cos(point.heel_radians)
        )

    def gz_at(self, heel_degrees: float) -> float:
        """``GZ`` at an existing sample; interpolating is a design choice left to callers."""
        for point in self.points:
            if abs(point.heel_degrees - heel_degrees) <= 1e-12:
                return self._gz(point)
        raise HeelError(f"heel angle {heel_degrees!r} is not one of the sampled angles {self.angles!r}")

    def maximum(self) -> tuple[float, float]:
        """``(angle, GZ)`` at the largest sampled righting arm (degrees, metres)."""
        values = self.gz
        if not values:
            raise HeelError("empty curve")
        best = max(range(len(values)), key=lambda index: values[index])
        return self.points[best].heel_degrees, values[best]

    def first_vanishing_angle(self) -> float | None:
        """Heel angle where ``GZ`` returns to zero after its peak, by linear interpolation.

        ``None`` when the sampled range does not reach it -- reported rather than guessed,
        because a curve that is still positive at its last sample has no vanishing angle
        inside the sampled range.
        """
        values = self.gz
        peak = max(range(len(values)), key=lambda index: values[index])
        for index in range(peak, len(values) - 1):
            a, b = values[index], values[index + 1]
            if a > 0.0 >= b:
                x0, x1 = self.points[index].heel_degrees, self.points[index + 1].heel_degrees
                return x0 + (x1 - x0) * (a / (a - b))
        return None

    def to_dict(self) -> dict[str, Any]:
        return {
            "kg": _round(self.kg, 9),
            "yg": _round(self.yg, 9),
            "keel_z": _round(self.keel_z, 12),
            "displacement_volume": _round(self.displacement_volume, 12),
            "angles_degrees": [_round(angle, 9) for angle in self.angles],
            "gz": [_round(value, 12) for value in self.gz],
            "kn": [_round(value, 12) for value in self.kn],
            "points": [point.to_dict() for point in self.points],
        }


def gz_curve(
    mesh: Mesh,
    heels_degrees: Sequence[float],
    *,
    kg: float,
    displacement_volume: float | None = None,
    waterline_z: float = 0.0,
    yg: float = 0.0,
    closed: bool | None = None,
) -> RightingArmCurve:
    """``GZ`` at each heel angle, all at the *same* displaced volume.

    ``displacement_volume`` defaults to the volume the body displaces at ``waterline_z``,
    i.e. the upright loading condition. The volume is then held fixed across heel -- the
    definition of a ``GZ`` curve -- so the returned ``displacement_volume`` is the target
    and each point reports the volume it actually solved (they agree to the bisection
    tolerance; ``tests/test_heeling.py`` asserts it).
    """
    _finite("kg", kg)
    _finite("yg", yg)
    _require_closed(mesh, closed)
    if displacement_volume is None:
        # The upright displacement is the volume below ``waterline_z``, measured by the
        # same clipping kernel that will be used at heel -- not the enclosed volume of the
        # solid, which would be wrong for any body that carries freeboard.
        displacement_volume = submerged_properties(mesh, (0.0, 0.0, 1.0), waterline_z, closed=closed)[0]
    _finite("displacement_volume", displacement_volume)
    if displacement_volume <= 0.0:
        raise HeelError(f"displacement volume must be positive, got {displacement_volume!r}")
    points = tuple(
        heeled_float(mesh, heel, displacement_volume=displacement_volume, closed=closed)
        for heel in heels_degrees
    )
    keel_z = min(vertex[2] for vertex in mesh.vertices)
    return RightingArmCurve(
        kg=kg, yg=yg, displacement_volume=displacement_volume, points=points, keel_z=keel_z
    )


def gz_initial_slope(
    mesh: Mesh,
    *,
    kg: float,
    displacement_volume: float | None = None,
    waterline_z: float = 0.0,
    yg: float = 0.0,
    step_degrees: float = 0.25,
    closed: bool | None = None,
) -> dict[str, float]:
    """``dGZ/dphi`` at upright by a central difference of the curve itself.

    The reference is ``GM = KB + BM - KG`` assembled from *upright* geometry only:
    ``BM = I_T / volume`` with ``I_T`` from
    :func:`tensorcad.stability.waterplane_moments` -- which cuts the waterplane and
    integrates it round its own loops with Green's theorem -- and ``KB``/volume from the
    clipped solid. None of that touches a heeled equilibrium, so agreement between the two
    numbers is a real check on the clip, the centroid, the waterline solve and the arm, not
    a restatement of one formula twice.
    """
    _finite("step_degrees", step_degrees)
    if step_degrees <= 0.0 or step_degrees > 5.0:
        raise HeelError(f"step_degrees must lie in (0, 5], got {step_degrees!r}")
    if displacement_volume is None:
        displacement_volume = submerged_properties(mesh, (0.0, 0.0, 1.0), waterline_z, closed=closed)[0]
    upright = submerged_properties(mesh, (0.0, 0.0, 1.0), waterline_z, closed=closed)
    keel_z = min(vertex[2] for vertex in mesh.vertices)
    kb = upright[1][2] - keel_z
    bm = waterplane_moments(mesh, waterline_z).transverse_second_moment / upright[0]
    # ``kg`` is a height above the keel, so GM = KB + BM - KG needs no datum fix-up; the
    # heeled curve below gets the same datum (G at keel_z + kg).
    gm = kb + bm - kg
    curve = gz_curve(
        mesh,
        (-step_degrees, step_degrees),
        kg=kg,
        displacement_volume=displacement_volume,
        yg=yg,
        closed=closed,
    )
    minus, plus = curve.gz  # the first sample is the negative heel
    slope = (plus - minus) / (2.0 * math.radians(step_degrees))
    return {
        "slope": slope,
        "gm": gm,
        "kb": kb,
        "bm": bm,
        "difference": slope - gm,
        "step_degrees": step_degrees,
    }


def cross_curve(
    mesh: Mesh,
    heel_degrees: float,
    *,
    displacement_volume: float,
    closed: bool | None = None,
) -> HeeledFloat:
    """``KN`` at one heel angle and one displacement: hull-only, no ``KG`` in sight."""
    return heeled_float(mesh, heel_degrees, displacement_volume=displacement_volume, closed=closed)


def cross_curves(
    mesh: Mesh,
    heels_degrees: Sequence[float],
    volumes: Sequence[float],
    *,
    closed: bool | None = None,
) -> tuple[tuple[float, ...], ...]:
    """The cross-curve table: rows indexed by ``volumes``, columns by ``heels_degrees``."""
    table: list[tuple[float, ...]] = []
    for volume in volumes:
        _finite("displacement_volume", volume)
        row = tuple(
            heeled_float(mesh, heel, displacement_volume=volume, closed=closed).kn
            for heel in heels_degrees
        )
        table.append(row)
    return tuple(table)