"""Righting arms at finite heel: the ``GZ`` curve of a closed body.

:mod:`tensorcad.stability` answers the upright question -- ``GM = KB + BM - KG``, the
slope of the righting-arm curve at zero heel -- and says in as many words that it does
**not** give the curve itself: *"Nor does it give a righting-arm curve ``GZ(phi)`` for
finite heel ... stability at large angles is a different integration (a heeled
displacement integration, not a waterplane moment)."* This module is that integration.

The definition is the one a loading computer uses. A body of displacement ``V0`` is
heeled by ``phi`` about the longitudinal axis at **constant displacement**; the water
surface is a plane whose normal, written in the body frame, is

    u(phi) = (0, -sin(phi), cos(phi))

(positive ``phi`` = starboard side down), so the submerged side is ``u . r <= c``. For
each ``phi`` the waterline offset ``c`` is solved so that the submerged volume is exactly
``V0`` -- that is the physical constraint, and it is the piece a naive "cut at the
upright waterline and rotate the answer" gets wrong. The centre of buoyancy ``B`` comes
out of the same integral, and the righting arm is the lever between ``B`` and ``G``
measured perpendicular to the water surface:

    GZ(phi) = y_B cos(phi) + (z_B - z_G) sin(phi),      z_G = z_keel + KG

(with ``B`` and ``G`` in the body frame, ``y`` to starboard). At small heel this reduces
to ``GZ = GM sin(phi)``, which is the consistency check against ``initial_stability``:
the two modules must agree on the slope even though they compute different objects.

**How the submerged body is measured.** Each triangle is clipped against the half-space
``u . r <= c`` (Sutherland-Hodgman on the three signed vertex heights; no geometric
tolerance anywhere -- a vertex exactly on the plane is kept once, a crossing is
interpolated), and the clipped polygon is split into tetrahedra with a **single apex on
the water surface**, ``o = c u``. The tetrahedra of the waterplane lid are then
degenerate and contribute exactly zero, so the signed tetrahedron sum over the clipped
hull triangles alone is the exact submerged volume and first moment:

    V = (1/6) sum_faces det[a-o, b-o, c-o],      int r dV = sum_tet V_tet (o+a+b+c)/4

which is an identity for a closed surface, not a sampler. It is *exact* for the
polyhedron the mesh actually is, which is why the box and sphere checks close to
round-off rather than to a mesh error: the only remaining error is the one-dimensional
root-find on ``c``. A heeled waterline that leaves the deck, a non-convex section and an
asymmetric body all need no special case.

**What this is not.** Nothing here is a statutory intact-stability check. There is no
free-surface correction (that is a tank property, ``sum(rho i_x)/displacement``, for the
caller who knows the tankage), no damage stability, no wind or launching heeling moment,
and no weathertight-closure assumption: the closed body *is* the hull, so its deck is
whatever face the mesh closes on. A hull mesh capped at the design waterline has zero
freeboard and therefore no heeled equilibrium at all (``HeelError``) -- computing a
ship's ``GZ`` needs the body closed at the deck or bulwark, not at the waterline. ``KG``
is an input loading condition, exactly as in
:func:`tensorcad.stability.initial_stability`, and ``G`` is assumed on the centreline.

Only the Python standard library is used.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

from .geometry import Mesh
from .stability import waterplane_moments

SCHEMA_SUBMERGED = "tensorcad.heel.submerged/1"
SCHEMA_GZ = "tensorcad.heel.gz_curve/1"

#: Heel angles sampled by :func:`gz_curve` when the caller does not pass a list.
DEFAULT_ANGLES_DEG: tuple[float, ...] = tuple(float(angle) for angle in range(0, 65, 5))


class HeelError(ValueError):
    """Raised when a heeled equilibrium at the requested displacement does not exist."""


@dataclass(frozen=True)
class SubmergedBody:
    """The immersed part of a closed body below one tilted water surface.

    ``waterline_offset`` is ``c = u . r`` on the surface in the body frame, so the water
    plane is ``{r : u . r = c}``; for ``heel_deg = 0`` and a body floating at
    ``waterline_z`` it is just ``waterline_z``. ``centre_of_buoyancy`` is the volume
    centroid of the immersed part, in the same frame as the mesh.
    """

    heel_deg: float
    waterplane_normal: tuple[float, float, float]
    waterline_offset: float
    displacement_volume: float
    centre_of_buoyancy: tuple[float, float, float]
    residual_volume: float
    basis: str = "clipped closed surface, signed tetrahedron sum"
    schema: str = SCHEMA_SUBMERGED

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": self.schema,
            "heel_deg": _round(self.heel_deg),
            "waterplane_normal": [_round(value) for value in self.waterplane_normal],
            "waterline_offset_m": _round(self.waterline_offset),
            "displacement_volume_m3": _round(self.displacement_volume),
            "centre_of_buoyancy_m": [_round(value) for value in self.centre_of_buoyancy],
            "residual_volume_m3": _round(self.residual_volume),
            "basis": self.basis,
        }


@dataclass(frozen=True)
class HeeledSection:
    """One point of a righting-arm curve."""

    heel_deg: float
    waterline_offset: float
    displacement_volume: float
    centre_of_buoyancy: tuple[float, float, float]
    gz: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "heel_deg": _round(self.heel_deg),
            "waterline_offset_m": _round(self.waterline_offset),
            "displacement_volume_m3": _round(self.displacement_volume),
            "centre_of_buoyancy_m": [_round(value) for value in self.centre_of_buoyancy],
            "gz_m": _round(self.gz),
        }


@dataclass(frozen=True)
class GZCurve:
    """Righting arms at constant displacement, with the shape metrics of the curve."""

    kg: float
    keel_z: float
    displacement_volume: float
    waterline_z: float
    sections: tuple[HeeledSection, ...]
    gm_transverse: float
    tolerance: float
    basis: str = "heeled displacement integration at constant displacement"
    schema: str = SCHEMA_GZ

    @property
    def angles_deg(self) -> tuple[float, ...]:
        return tuple(section.heel_deg for section in self.sections)

    @property
    def gz(self) -> tuple[float, ...]:
        return tuple(section.gz for section in self.sections)

    @property
    def max_gz(self) -> tuple[float, float]:
        """``(heel_deg, gz_m)`` at the largest righting arm on the sampled angles."""
        best = max(self.sections, key=lambda section: section.gz)
        return best.heel_deg, best.gz

    @property
    def vanishing_angle_deg(self) -> float | None:
        """First zero crossing after the maximum, linearly interpolated.

        ``None`` when the curve has not returned to zero by the last sampled angle --
        the honest answer for a truncated range, not a guess at the capsizing angle.
        """
        peak = max(range(len(self.sections)), key=lambda index: self.sections[index].gz)
        for index in range(peak, len(self.sections) - 1):
            before = self.sections[index]
            after = self.sections[index + 1]
            if before.gz > 0.0 >= after.gz:
                span = before.gz - after.gz
                if span <= 0.0:
                    return after.heel_deg
                return before.heel_deg + (after.heel_deg - before.heel_deg) * (before.gz / span)
        return None

    def at(self, heel_deg: float) -> float:
        """Linearly interpolated ``GZ`` at ``heel_deg`` inside the sampled range."""
        angles = self.angles_deg
        if not angles:
            raise HeelError("empty GZ curve")
        if heel_deg < angles[0] or heel_deg > angles[-1]:
            raise HeelError(
                f"heel {heel_deg} deg is outside the sampled range "
                f"{angles[0]}..{angles[-1]} deg"
            )
        for index in range(len(angles) - 1):
            lo, hi = angles[index], angles[index + 1]
            if lo <= heel_deg <= hi:
                before = self.sections[index]
                after = self.sections[index + 1]
                if hi == lo:
                    return before.gz
                weight = (heel_deg - lo) / (hi - lo)
                return before.gz + weight * (after.gz - before.gz)
        return self.sections[-1].gz

    def to_dict(self) -> dict[str, Any]:
        max_angle, max_gz = self.max_gz
        vanishing = self.vanishing_angle_deg
        return {
            "schema": self.schema,
            "kg_m": _round(self.kg),
            "keel_z_m": _round(self.keel_z),
            "displacement_volume_m3": _round(self.displacement_volume),
            "waterline_z_m": _round(self.waterline_z),
            "gm_transverse_m": _round(self.gm_transverse),
            "max_gz_m": _round(max_gz),
            "angle_of_max_gz_deg": _round(max_angle),
            "vanishing_angle_deg": None if vanishing is None else _round(vanishing),
            "tolerance_m3": _round(self.tolerance),
            "sections": [section.to_dict() for section in self.sections],
            "basis": self.basis,
        }


def _round(value: float, digits: int = 9) -> float:
    rounded = round(value, digits)
    return 0.0 if rounded == 0.0 else rounded


def waterplane_normal(heel_deg: float) -> tuple[float, float, float]:
    """The body-frame "up" of a body heeled by ``heel_deg`` (positive = starboard down)."""
    if not math.isfinite(heel_deg):
        raise HeelError(f"heel angle must be finite, got {heel_deg!r}")
    phi = math.radians(heel_deg)
    return (0.0, -math.sin(phi), math.cos(phi))


def _clip(
    triangle: tuple[tuple[float, float, float], ...],
    signed: tuple[float, float, float],
) -> list[tuple[float, float, float]]:
    """Keep the part of one triangle with ``signed >= 0`` (distance above the plane)."""
    polygon: list[tuple[float, float, float]] = []
    for index in range(3):
        a = triangle[index]
        b = triangle[(index + 1) % 3]
        sa = signed[index]
        sb = signed[(index + 1) % 3]
        if sa >= 0.0:
            polygon.append(a)
        if (sa > 0.0 > sb) or (sa < 0.0 < sb):
            t = sa / (sa - sb)
            polygon.append(
                (
                    a[0] + t * (b[0] - a[0]),
                    a[1] + t * (b[1] - a[1]),
                    a[2] + t * (b[2] - a[2]),
                )
            )
    return polygon


def _submerged_terms(
    mesh: Mesh, normal: tuple[float, float, float], offset: float
) -> tuple[float, tuple[float, float, float]]:
    """Exact volume and first moment of the part of ``mesh`` below ``u . r = offset``."""
    ux, uy, uz = normal
    ax, ay, az = ux * offset, uy * offset, uz * offset
    volume = 0.0
    mx = my = mz = 0.0
    vertices = mesh.vertices
    for i0, i1, i2 in mesh.faces:
        p0, p1, p2 = vertices[i0], vertices[i1], vertices[i2]
        s0 = offset - (ux * p0[0] + uy * p0[1] + uz * p0[2])
        s1 = offset - (ux * p1[0] + uy * p1[1] + uz * p1[2])
        s2 = offset - (ux * p2[0] + uy * p2[1] + uz * p2[2])
        if s0 < 0.0 and s1 < 0.0 and s2 < 0.0:
            continue
        polygon = _clip((p0, p1, p2), (s0, s1, s2))
        if len(polygon) < 3:
            continue
        x0 = polygon[0][0] - ax
        y0 = polygon[0][1] - ay
        z0 = polygon[0][2] - az
        for index in range(1, len(polygon) - 1):
            x1 = polygon[index][0] - ax
            y1 = polygon[index][1] - ay
            z1 = polygon[index][2] - az
            x2 = polygon[index + 1][0] - ax
            y2 = polygon[index + 1][1] - ay
            z2 = polygon[index + 1][2] - az
            cross_x = y1 * z2 - z1 * y2
            cross_y = z1 * x2 - x1 * z2
            cross_z = x1 * y2 - y1 * x2
            six = x0 * cross_x + y0 * cross_y + z0 * cross_z
            volume += six
            mx += six * (x0 + x1 + x2)
            my += six * (y0 + y1 + y2)
            mz += six * (z0 + z1 + z2)
    volume *= 1.0 / 6.0
    scale = 1.0 / 24.0
    return volume, (
        mx * scale + ax * volume,
        my * scale + ay * volume,
        mz * scale + az * volume,
    )


def _normal_extremes(mesh: Mesh, normal: tuple[float, float, float]) -> tuple[float, float]:
    if not mesh.vertices:
        raise HeelError("cannot heel an empty mesh")
    ux, uy, uz = normal
    heights = [ux * v[0] + uy * v[1] + uz * v[2] for v in mesh.vertices]
    return min(heights), max(heights)


def _solve_offset(
    mesh: Mesh,
    normal: tuple[float, float, float],
    target: float,
    tolerance: float,
) -> tuple[float, float]:
    """Solve ``V(c) = target`` for the waterline offset by bracketed regula falsi.

    ``V(c)`` is monotone increasing in ``c`` (raising the water surface cannot remove
    submerged volume), so a bracket always exists; Illinois damping keeps the iteration
    from stalling against one end on a nearly piecewise-linear function. Returns
    ``(offset, residual)``.
    """
    lo, hi = _normal_extremes(mesh, normal)
    span = hi - lo
    vol_lo = _submerged_terms(mesh, normal, lo)[0]
    vol_hi = _submerged_terms(mesh, normal, hi)[0]
    total = vol_hi - vol_lo
    if total <= 0.0:
        raise HeelError(
            "the mesh encloses no volume below its own bounding planes; it is probably "
            "not closed or not outward-oriented"
        )
    if not math.isfinite(target) or target <= 0.0:
        raise HeelError(f"displacement volume must be positive, got {target!r} m^3")
    if target > total * (1.0 + 1e-12) or total - target <= tolerance:
        raise HeelError(
            f"requested displacement {target:.9g} m^3 is the body's whole volume "
            f"{total:.9g} m^3 (or more): the water surface would not cut the body at all. "
            "A hull mesh capped at the design waterline has zero freeboard and therefore "
            "no heeled equilibrium; close the body at the deck before asking for GZ"
        )

    a, b = lo, hi
    fa = vol_lo - target
    fb = vol_hi - target
    best, best_residual = hi, abs(fb)
    for _ in range(200):
        denominator = fb - fa
        c = b if denominator == 0.0 else (a * fb - b * fa) / denominator
        if not a < c < b:
            c = 0.5 * (a + b)
        fc = _submerged_terms(mesh, normal, c)[0] - target
        if abs(fc) < best_residual:
            best, best_residual = c, abs(fc)
        if abs(fc) <= tolerance or b - a <= 1e-15 * abs(span):
            return c, fc
        if fc < 0.0:
            a, fa = c, fc
            fb *= 0.5
        else:
            b, fb = c, fc
            fa *= 0.5
    return best, best_residual


def heeled_displacement(
    mesh: Mesh,
    *,
    heel_deg: float,
    displacement_volume: float,
    tolerance: float | None = None,
) -> SubmergedBody:
    """The immersed part of ``mesh`` below the water surface holding ``displacement_volume``.

    The displacement is held fixed and the waterline offset is solved for, which is what
    a righting-arm definition requires. Raises :class:`HeelError` when no such surface
    exists -- a displacement above the body's own volume, or a body with no freeboard --
    rather than returning a silently clipped volume.
    """
    normal = waterplane_normal(heel_deg)
    if tolerance is None:
        tolerance = 1e-12 * max(displacement_volume, 1.0)
    if not math.isfinite(tolerance) or tolerance <= 0.0:
        raise HeelError(f"tolerance must be positive and finite, got {tolerance!r}")
    offset, residual = _solve_offset(mesh, normal, displacement_volume, tolerance)
    volume, moment = _submerged_terms(mesh, normal, offset)
    if abs(residual) > 10.0 * tolerance:
        raise HeelError(
            f"heeled waterline did not converge at {heel_deg} deg: residual "
            f"{residual:.3e} m^3 exceeds {10.0 * tolerance:.3e} m^3"
        )
    return SubmergedBody(
        heel_deg=heel_deg,
        waterplane_normal=normal,
        waterline_offset=offset,
        displacement_volume=volume,
        centre_of_buoyancy=(
            moment[0] / volume,
            moment[1] / volume,
            moment[2] / volume,
        ),
        residual_volume=residual,
    )


def righting_arm(
    body: SubmergedBody,
    *,
    kg: float,
    keel_z: float,
) -> float:
    """``GZ`` of one solved heeled body for a centre of gravity ``KG`` above the keel."""
    if not math.isfinite(kg):
        raise HeelError(f"kg must be finite, got {kg!r}")
    phi = math.radians(body.heel_deg)
    y_b, z_b = body.centre_of_buoyancy[1], body.centre_of_buoyancy[2]
    return y_b * math.cos(phi) + (z_b - (keel_z + kg)) * math.sin(phi)


def gz_curve(
    mesh: Mesh,
    *,
    kg: float,
    angles_deg: Any = None,
    displacement_volume: float | None = None,
    waterline_z: float = 0.0,
    tolerance: float | None = None,
) -> GZCurve:
    """Righting arms of a closed body over a list of heel angles.

    ``kg`` is measured above the keel (the mesh's lowest vertex) and is a required
    keyword for the same reason as in :func:`tensorcad.stability.initial_stability`: a
    default would silently produce ``GZ`` for a loading condition nobody chose.
    ``displacement_volume`` defaults to the volume of the body below the horizontal
    plane ``z = waterline_z`` -- the *submerged* volume, measured here rather than taken
    from :func:`tensorcad.hydrostatics.compute_hydrostatics`, whose ``displacement`` is
    the volume of the whole closed body (correct only when the mesh is capped at the
    waterline, i.e. when the body has no freeboard and no heeled equilibrium exists).

    ``gm_transverse`` on the result is ``KB + BM - KG`` computed from *this* module's
    upright submerged body and :func:`tensorcad.stability.waterplane_moments`, so
    ``dGZ/dphi`` at zero heel can be compared against a number that came from a plane
    section rather than from the heeled integration.
    """
    angles = DEFAULT_ANGLES_DEG if angles_deg is None else tuple(float(a) for a in angles_deg)
    if not angles:
        raise HeelError("at least one heel angle is required")
    if not math.isfinite(kg):
        raise HeelError(f"kg must be finite, got {kg!r}")
    keel_z = min(vertex[2] for vertex in mesh.vertices)
    if displacement_volume is None:
        target = _submerged_terms(mesh, (0.0, 0.0, 1.0), waterline_z)[0]
    else:
        target = displacement_volume
    if not math.isfinite(target) or target <= 0.0:
        raise HeelError(
            f"displacement volume must be positive, got {target!r} m^3; a waterline at "
            f"z={waterline_z!r} that does not submerge the body has no righting arm"
        )
    if tolerance is None:
        tolerance = 1e-12 * max(target, 1.0)
    sections = []
    for angle in sorted(angles):
        body = heeled_displacement(
            mesh, heel_deg=angle, displacement_volume=target, tolerance=tolerance
        )
        sections.append(
            HeeledSection(
                heel_deg=angle,
                waterline_offset=body.waterline_offset,
                displacement_volume=body.displacement_volume,
                centre_of_buoyancy=body.centre_of_buoyancy,
                gz=righting_arm(body, kg=kg, keel_z=keel_z),
            )
        )
    upright = heeled_displacement(
        mesh, heel_deg=0.0, displacement_volume=target, tolerance=tolerance
    )
    kb = upright.centre_of_buoyancy[2] - keel_z
    try:
        section = waterplane_moments(mesh, upright.waterline_offset)
        bm = section.transverse_second_moment / target
    except Exception:  # noqa: BLE001 - GM is reported, not required, for a heeled body
        bm = float("nan")
    return GZCurve(
        kg=kg,
        keel_z=keel_z,
        displacement_volume=target,
        waterline_z=waterline_z,
        sections=tuple(sections),
        gm_transverse=kb + bm - kg,
        tolerance=tolerance,
    )


def wall_sided_gz(*, kb: float, bm: float, kg: float, heel_deg: float) -> float:
    """The closed-form wall-sided righting arm, for verification.

        GZ(phi) = sin(phi) * (GM + BM tan^2(phi) / 2),      GM = KB + BM - KG

    Exact for a body whose sides are vertical at the waterline, up to the angle at which
    the deck edge immerses -- for a rectangular barge of beam ``B`` and freeboard ``f``
    that is ``atan(2 f / B)``. Beyond it the formula is *not* a reference and this module's
    integration is the one to trust; ``examples/gz_curve`` measures where the two part
    company instead of assuming the range.
    """
    for name, value in (("kb", kb), ("bm", bm), ("kg", kg), ("heel_deg", heel_deg)):
        if not math.isfinite(value):
            raise HeelError(f"{name} must be finite, got {value!r}")
    if bm < 0.0:
        raise HeelError(f"bm must be non-negative, got {bm!r}")
    phi = math.radians(heel_deg)
    gm = kb + bm - kg
    return math.sin(phi) * (gm + 0.5 * bm * math.tan(phi) ** 2)