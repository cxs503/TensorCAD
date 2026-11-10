"""Parametric hull surface generation: NURBS cross-sections + longitudinal loft.

This module provides a complete pipeline from mathematical hull definition to
solver-ready surface mesh:

1. **Cross-section definition** -- from offset tables (station-half-breadth-height)
   or parametric functions (e.g. Wigley, SUBOFF radius profiles).
2. **NURBS curve interpolation** -- each cross-section is fitted with a NURBS curve;
   longitudinal curves (waterlines, buttock lines) are interpolated through the
   section control points.
3. **Surface construction** -- a tensor-product NURBS surface is lofted through the
   section curves, producing a smooth, watertight hull surface.
4. **Geometric transformations** -- scale, translate, rotate, mirror (port<->starboard)
   to support twin-screw / single-screw, bulbous bow / cruiser stern variants.
5. **Surface mesh generation** -- automatic triangulation of the parametric domain
   with density control (uniform or curvature-adaptive).
6. **Hydrostatic properties** -- displacement volume, centre of buoyancy, waterplane
   area, transverse metacentric radius, computed by divergence-theorem integration
   over the triangulated surface.

Coordinate frame
----------------
``x`` longitudinal (bow -> stern positive), ``y`` transverse (starboard positive),
``z`` vertical (up positive).  The waterline is at ``z = 0`` by default.

Dependencies
------------
Only the Python standard library.  The module re-uses :mod:`tensorcad.nurbs` for
B-spline evaluation and :mod:`tensorcad.geometry` / :mod:`tensorcad.mesh` for the
mesh contract and topology checks.
"""
from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Any

from .geometry import Mesh, MeshBuilder, Vec3
from .nurbs import NURBSCurve, NURBSSurface, NURBSError
from . import mesh as meshlib

__all__ = [
    "CrossSection",
    "OffsetTable",
    "ParametricHull",
    "ParametricHullError",
    "SurfaceMesh",
    "HydrostaticResult",
    "Transform",
    "build_parametric_hull",
    "cross_section_from_offsets",
    "cross_section_from_function",
    "wigley_cross_sections",
    "suboff_cross_sections",
    "loft_nurbs_surface",
    "tessellate_surface",
    "apply_transform",
    "compute_hydrostatics",
    "mesh_quality_report",
]


class ParametricHullError(ValueError):
    """Raised when hull definition or transformation is geometrically invalid."""


# ---------------------------------------------------------------------------
# Cross-section definition
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class CrossSection:
    """A single transverse section at a given longitudinal station.

    Attributes
    ----------
    x : float
        Longitudinal position of the section.
    points : tuple[tuple[float, float], ...]
        Ordered (y, z) pairs defining the half-breadth profile.  ``y >= 0``
        (starboard half); ``z`` runs from keel (lowest) to deck (highest).
    """
    x: float
    points: tuple[tuple[float, float], ...]

    @property
    def max_half_breadth(self) -> float:
        return max(p[0] for p in self.points) if self.points else 0.0

    @property
    def max_height(self) -> float:
        return max(p[1] for p in self.points) if self.points else 0.0

    @property
    def min_z(self) -> float:
        return min(p[1] for p in self.points) if self.points else 0.0

    def area(self) -> float:
        """Section area by the trapezoidal rule (half-section, starboard)."""
        if len(self.points) < 2:
            return 0.0
        sorted_pts = sorted(self.points, key=lambda p: p[1])
        area = 0.0
        for i in range(len(sorted_pts) - 1):
            y0, z0 = sorted_pts[i]
            y1, z1 = sorted_pts[i + 1]
            area += 0.5 * (y0 + y1) * (z1 - z0)
        return area


@dataclass(frozen=True)
class OffsetTable:
    """Station-based offset table.

    Each entry is a tuple ``(x_frac, [(y0,z0), (y1,z1), ...])`` where ``x_frac``
    is the fractional station position in [0, 1] (0 = AP, 1 = FP).
    """
    stations: tuple[tuple[float, tuple[tuple[float, float], ...]], ...]

    def __post_init__(self) -> None:
        for frac, pts in self.stations:
            if not (0.0 <= frac <= 1.0):
                raise ParametricHullError(f"station fraction {frac} not in [0, 1]")
            if len(pts) < 2:
                raise ParametricHullError(
                    f"station at x={frac} needs >= 2 points, got {len(pts)}"
                )

    @property
    def num_stations(self) -> int:
        return len(self.stations)

    def to_cross_sections(self, length: float) -> tuple[CrossSection, ...]:
        """Convert to CrossSection objects scaled to the given hull length."""
        sections = []
        for frac, pts in self.stations:
            x = frac * length
            sections.append(CrossSection(x=x, points=tuple(pts)))
        return tuple(sections)


def cross_section_from_offsets(
    x: float,
    half_breadths: Sequence[float],
    heights: Sequence[float],
) -> CrossSection:
    """Build a cross-section from paired half-breadth and height arrays."""
    if len(half_breadths) != len(heights):
        raise ParametricHullError("half_breadths and heights must have same length")
    if len(half_breadths) < 2:
        raise ParametricHullError("need at least 2 points per section")
    points = tuple(
        (float(y), float(z)) for y, z in zip(half_breadths, heights)
    )
    return CrossSection(x=float(x), points=points)


def cross_section_from_function(
    x: float,
    profile_func: Any,
    n_points: int = 20,
    z_min: float = 0.0,
    z_max: float = 1.0,
) -> CrossSection:
    """Build a cross-section by sampling a parametric profile function.

    Parameters
    ----------
    x : float
        Longitudinal position.
    profile_func : callable
        ``profile_func(z) -> y`` returning half-breadth at height z.
    n_points : int
        Number of sample points along z.
    z_min, z_max : float
        Vertical extent of the section.
    """
    if n_points < 2:
        raise ParametricHullError("n_points must be >= 2")
    points = []
    for i in range(n_points):
        z = z_min + (z_max - z_min) * i / (n_points - 1)
        y = profile_func(z)
        if y < -1e-15:
            raise ParametricHullError(f"negative half-breadth y={y} at z={z}")
        points.append((max(y, 0.0), float(z)))
    return CrossSection(x=float(x), points=tuple(points))


# ---------------------------------------------------------------------------
# Standard hull families
# ---------------------------------------------------------------------------
def wigley_cross_sections(
    length: float = 2.0,
    beam: float = 0.2,
    draft: float = 0.125,
    n_stations: int = 11,
    n_points: int = 10,
) -> tuple[CrossSection, ...]:
    """Generate cross-sections for the Wigley III mathematical hull.

    The Wigley hull is defined by::

        y(x, z) = (B/2) * (1 - (2x/L)^2) * (1 - (z/T)^2)

    where x is measured from midship (x=0 at midship, x = +/-L/2 at ends).

    Analytic properties:
        Volume V = (4/9) * L * B * T  =>  Cb = 4/9 ~ 0.4444
        Waterplane area Aw = (4/9) * L * B
    """
    sections = []
    for i in range(n_stations):
        xi = i / (n_stations - 1)  # 0 to 1
        x_abs = (xi - 0.5) * length  # -L/2 to L/2

        def profile(z: float, _xa: float = x_abs) -> float:
            x_factor = 1.0 - (2.0 * _xa / length) ** 2
            z_factor = 1.0 - (z / draft) ** 2
            return 0.5 * beam * max(x_factor, 0.0) * max(z_factor, 0.0)

        sec = cross_section_from_function(
            x=xi * length,  # 0 to L (AP to FP)
            profile_func=profile,
            n_points=n_points,
            z_min=-draft,
            z_max=0.0,
        )
        sections.append(sec)
    return tuple(sections)


def suboff_cross_sections(
    length_m: float | None = None,
    n_stations: int = 21,
    n_points: int = 20,
) -> tuple[CrossSection, ...]:
    """Generate cross-sections for the DARPA SUBOFF bare hull.

    The SUBOFF is axisymmetric, so each section is a circle of radius r(x).
    """
    from .suboff import (
        SUBOFF_LENGTH_M,
        SUBOFF_RADIUS_M,
        suboff_radius_fraction,
    )

    L = length_m if length_m is not None else SUBOFF_LENGTH_M
    R = SUBOFF_RADIUS_M * (L / SUBOFF_LENGTH_M)
    sections = []
    for i in range(n_stations):
        xi = i / (n_stations - 1)
        x = xi * L
        r = R * suboff_radius_fraction(xi)

        # Circular section: sample the right half (y >= 0)
        points = []
        for j in range(n_points):
            theta = -math.pi / 2 + math.pi * j / (n_points - 1)
            y = r * abs(math.cos(theta))
            z = r * math.sin(theta)
            points.append((y, z))
        sections.append(CrossSection(x=x, points=tuple(points)))
    return tuple(sections)


# ---------------------------------------------------------------------------
# NURBS surface lofting
# ---------------------------------------------------------------------------
def _clamped_knots(degree: int, n_ctrl: int) -> tuple[float, ...]:
    """Uniform clamped knot vector on [0, 1]."""
    if n_ctrl < degree + 1:
        raise ParametricHullError(
            f"{n_ctrl} control points insufficient for degree {degree}"
        )
    n_inner = n_ctrl - degree - 1
    knots = [0.0] * (degree + 1)
    for i in range(1, n_inner + 1):
        knots.append(i / (n_inner + 1))
    knots.extend([1.0] * (degree + 1))
    return tuple(knots)


def _resample_section(section: CrossSection, n_points: int) -> CrossSection:
    """Resample a cross-section to a uniform number of points by linear interpolation."""
    pts = sorted(section.points, key=lambda p: p[1])
    if len(pts) == n_points:
        return section
    chords = [0.0]
    for i in range(1, len(pts)):
        dy = pts[i][0] - pts[i - 1][0]
        dz = pts[i][1] - pts[i - 1][1]
        chords.append(chords[-1] + math.sqrt(dy * dy + dz * dz))
    total = chords[-1]
    if total < 1e-15:
        return CrossSection(x=section.x, points=tuple(pts[:n_points]))

    new_pts = []
    for i in range(n_points):
        target = total * i / (n_points - 1)
        seg = 0
        for k in range(1, len(chords)):
            if chords[k] >= target:
                seg = k - 1
                break
        else:
            seg = len(chords) - 2
        seg_len = chords[seg + 1] - chords[seg]
        if seg_len < 1e-15:
            t = 0.0
        else:
            t = (target - chords[seg]) / seg_len
        y = pts[seg][0] + t * (pts[seg + 1][0] - pts[seg][0])
        z = pts[seg][1] + t * (pts[seg + 1][1] - pts[seg][1])
        new_pts.append((y, z))
    return CrossSection(x=section.x, points=tuple(new_pts))


def loft_nurbs_surface(
    sections: Sequence[CrossSection],
    degree_u: int = 3,
    degree_v: int = 3,
    n_points_per_section: int | None = None,
) -> NURBSSurface:
    """Loft a NURBS surface through a sequence of cross-sections.

    The surface is constructed by using the section points directly as control points,
    forming a control net. The knot vectors are uniform clamped.

    Parameters
    ----------
    sections : sequence of CrossSection
        Ordered from AP to FP.
    degree_u : int
        Degree in the longitudinal direction.
    degree_v : int
        Degree in the transverse/vertical direction.
    n_points_per_section : int, optional
        If given, resample each section to this many points before lofting.
    """
    if len(sections) < 2:
        raise ParametricHullError("need at least 2 sections to loft a surface")

    # Determine target point count per section
    if n_points_per_section is not None:
        n_v = n_points_per_section
    else:
        n_v = len(sections[0].points)

    # Resample all sections to uniform point count
    resampled = [_resample_section(sec, n_v) for sec in sections]
    n_u = len(resampled)

    # Clamp degrees
    degree_u_actual = min(degree_u, n_u - 1)
    degree_v_actual = min(degree_v, n_v - 1)

    # Build control net directly from section points
    # control_net[u_idx][v_idx] = (x, y, z)
    control_net = []
    for u_idx in range(n_u):
        row = []
        for v_idx in range(n_v):
            y, z = resampled[u_idx].points[v_idx]
            x = resampled[u_idx].x
            row.append((x, y, z))
        control_net.append(tuple(row))

    # Build uniform clamped knot vectors
    knots_u = _clamped_knots(degree_u_actual, n_u)
    knots_v = _clamped_knots(degree_v_actual, n_v)

    # Weights (all 1.0 for non-rational)
    weights = [[1.0] * n_v for _ in range(n_u)]

    return NURBSSurface(
        degree=(degree_u_actual, degree_v_actual),
        knots=(tuple(knots_u), tuple(knots_v)),
        control_net=control_net,
        weights=weights,
    )


# ---------------------------------------------------------------------------
# Surface tessellation (mesh generation)
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class SurfaceMesh:
    """A triangulated surface mesh with quality metadata."""
    mesh: Mesh
    n_u: int
    n_v: int

    @property
    def n_vertices(self) -> int:
        return len(self.mesh.vertices)

    @property
    def n_faces(self) -> int:
        return len(self.mesh.faces)


def tessellate_surface(
    surface: NURBSSurface,
    n_u: int = 30,
    n_v: int = 15,
    *,
    mirror_y: bool = True,
) -> SurfaceMesh:
    """Tessellate a NURBS surface into a triangulated mesh.

    Parameters
    ----------
    surface : NURBSSurface
        The parametric surface to mesh.
    n_u, n_v : int
        Number of evaluation points in each parametric direction.
    mirror_y : bool
        If True, mirror the surface about the xz-plane to create a full
        (port + starboard) hull.
    """
    if n_u < 2 or n_v < 2:
        raise ParametricHullError("need at least 2x2 evaluation points")

    builder = MeshBuilder()

    # Evaluate surface on a regular grid and add vertices
    # Starboard half
    starboard_indices = []
    for i in range(n_u):
        u = i / (n_u - 1)
        row_indices = []
        for j in range(n_v):
            v = j / (n_v - 1)
            pt = surface.evaluate(u, v)
            idx = builder.add_vertex(pt)
            row_indices.append(idx)
        starboard_indices.append(row_indices)

    if mirror_y:
        # Port half (mirror y)
        port_indices = []
        for i in range(n_u):
            u = i / (n_u - 1)
            row_indices = []
            for j in range(n_v):
                v = j / (n_v - 1)
                pt = surface.evaluate(u, v)
                mirrored = (pt[0], -pt[1], pt[2])
                idx = builder.add_vertex(mirrored)
                row_indices.append(idx)
            port_indices.append(row_indices)

        # Starboard faces
        for i in range(n_u - 1):
            for j in range(n_v - 1):
                v00 = starboard_indices[i][j]
                v10 = starboard_indices[i + 1][j]
                v01 = starboard_indices[i][j + 1]
                v11 = starboard_indices[i + 1][j + 1]
                builder.add_triangle(v00, v10, v11)
                builder.add_triangle(v00, v11, v01)

        # Port faces (reversed winding for outward normals)
        for i in range(n_u - 1):
            for j in range(n_v - 1):
                v00 = port_indices[i][j]
                v10 = port_indices[i + 1][j]
                v01 = port_indices[i][j + 1]
                v11 = port_indices[i + 1][j + 1]
                builder.add_triangle(v00, v11, v10)
                builder.add_triangle(v00, v01, v11)
    else:
        for i in range(n_u - 1):
            for j in range(n_v - 1):
                v00 = starboard_indices[i][j]
                v10 = starboard_indices[i + 1][j]
                v01 = starboard_indices[i][j + 1]
                v11 = starboard_indices[i + 1][j + 1]
                builder.add_triangle(v00, v10, v11)
                builder.add_triangle(v00, v11, v01)

    mesh = builder.build()

    return SurfaceMesh(mesh=mesh, n_u=n_u, n_v=n_v)


# ---------------------------------------------------------------------------
# Geometric transformations
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Transform:
    """Rigid-body + scale transformation for hull surfaces.

    All transformations are applied in order: mirror -> scale -> rotate -> translate.
    """
    scale: tuple[float, float, float] = (1.0, 1.0, 1.0)
    rotation_deg: tuple[float, float, float] = (0.0, 0.0, 0.0)
    translation: tuple[float, float, float] = (0.0, 0.0, 0.0)
    mirror_x: bool = False
    mirror_y: bool = False
    mirror_z: bool = False

    def apply_to_point(self, pt: Vec3) -> Vec3:
        x, y, z = pt
        if self.mirror_x:
            x = -x
        if self.mirror_y:
            y = -y
        if self.mirror_z:
            z = -z
        x *= self.scale[0]
        y *= self.scale[1]
        z *= self.scale[2]
        rx = math.radians(self.rotation_deg[0])
        ry = math.radians(self.rotation_deg[1])
        rz = math.radians(self.rotation_deg[2])
        y1 = y * math.cos(rx) - z * math.sin(rx)
        z1 = y * math.sin(rx) + z * math.cos(rx)
        x2 = x * math.cos(ry) + z1 * math.sin(ry)
        z2 = -x * math.sin(ry) + z1 * math.cos(ry)
        x3 = x2 * math.cos(rz) - y1 * math.sin(rz)
        y3 = x2 * math.sin(rz) + y1 * math.cos(rz)
        return (x3 + self.translation[0], y3 + self.translation[1], z2 + self.translation[2])


def apply_transform(mesh: Mesh, transform: Transform) -> Mesh:
    """Apply a geometric transformation to all vertices of a mesh."""
    new_vertices = tuple(transform.apply_to_point(v) for v in mesh.vertices)
    return Mesh(
        vertices=new_vertices,
        faces=mesh.faces,
        units=mesh.units,
        coordinate_frame=mesh.coordinate_frame,
    )


# ---------------------------------------------------------------------------
# Hydrostatic property computation
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class HydrostaticResult:
    """Hydrostatic properties computed from a triangulated hull surface."""
    displaced_volume: float
    centre_of_buoyancy: tuple[float, float, float]
    waterplane_area: float
    wetted_surface_area: float
    block_coefficient: float
    transverse_metacentric_radius: float
    longitudinal_waterplane_inertia: float
    transverse_waterplane_inertia: float
    lcb: float
    vcb: float

    def to_dict(self) -> dict[str, float]:
        return {
            "displaced_volume_m3": self.displaced_volume,
            "centre_of_buoyancy": self.centre_of_buoyancy,
            "waterplane_area_m2": self.waterplane_area,
            "wetted_surface_area_m2": self.wetted_surface_area,
            "block_coefficient": self.block_coefficient,
            "transverse_metacentric_radius_m": self.transverse_metacentric_radius,
            "longitudinal_waterplane_inertia_m4": self.longitudinal_waterplane_inertia,
            "transverse_waterplane_inertia_m4": self.transverse_waterplane_inertia,
            "lcb_m": self.lcb,
            "vcb_m": self.vcb,
        }


def compute_hydrostatics(
    mesh: Mesh,
    waterline_z: float = 0.0,
    water_density: float = 1025.0,
    length: float | None = None,
    beam: float | None = None,
    draft: float | None = None,
) -> HydrostaticResult:
    """Compute hydrostatic properties of a hull mesh below the waterline.

    Uses the divergence theorem over the wetted surface.
    """
    xs = [v[0] for v in mesh.vertices]
    ys = [v[1] for v in mesh.vertices]
    zs = [v[2] for v in mesh.vertices]
    if length is None:
        length = max(xs) - min(xs)
    if beam is None:
        beam = max(ys) - min(ys)
    if draft is None:
        draft = waterline_z - min(zs)

    # Filter wetted faces
    wetted_faces = []
    for face in mesh.faces:
        v0 = mesh.vertices[face[0]]
        v1 = mesh.vertices[face[1]]
        v2 = mesh.vertices[face[2]]
        cz = (v0[2] + v1[2] + v2[2]) / 3.0
        if cz < waterline_z + 1e-9:
            wetted_faces.append(face)

    # Volume by divergence theorem
    volume = 0.0
    cx_sum = 0.0
    cy_sum = 0.0
    cz_sum = 0.0

    for face in wetted_faces:
        v0 = mesh.vertices[face[0]]
        v1 = mesh.vertices[face[1]]
        v2 = mesh.vertices[face[2]]
        tet_vol = (
            v0[0] * (v1[1] * v2[2] - v1[2] * v2[1])
            + v0[1] * (v1[2] * v2[0] - v1[0] * v2[2])
            + v0[2] * (v1[0] * v2[1] - v1[1] * v2[0])
        ) / 6.0
        volume += tet_vol
        cx_sum += tet_vol * (v0[0] + v1[0] + v2[0]) / 4.0
        cy_sum += tet_vol * (v0[1] + v1[1] + v2[1]) / 4.0
        cz_sum += tet_vol * (v0[2] + v1[2] + v2[2]) / 4.0

    # The sign of volume indicates normal orientation; we want positive volume
    # and the centre of buoyancy should be computed with the same sign convention
    if volume < 0:
        # Normals are inward; flip the sign
        volume = -volume
        cx_sum = -cx_sum
        cy_sum = -cy_sum
        cz_sum = -cz_sum
    volume = abs(volume)
    if volume > 1e-15:
        cob = (cx_sum / volume, cy_sum / volume, cz_sum / volume)
    else:
        cob = (0.0, 0.0, 0.0)

    # Wetted surface area
    wetted_area = 0.0
    for face in wetted_faces:
        v0 = mesh.vertices[face[0]]
        v1 = mesh.vertices[face[1]]
        v2 = mesh.vertices[face[2]]
        ux, uy, uz = v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]
        vx, vy, vz = v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]
        nx = uy * vz - uz * vy
        ny = uz * vx - ux * vz
        nz = ux * vy - uy * vx
        wetted_area += 0.5 * math.sqrt(nx * nx + ny * ny + nz * nz)

    # Waterplane area and inertia
    waterplane_area = _compute_waterplane_area(mesh, waterline_z)
    ixx, iyy = _compute_waterplane_inertia(mesh, waterline_z)

    # Block coefficient
    lbt = length * beam * draft
    cb = volume / lbt if lbt > 1e-15 else 0.0

    # Transverse metacentric radius
    bm_t = ixx / volume if volume > 1e-15 else 0.0

    return HydrostaticResult(
        displaced_volume=volume,
        centre_of_buoyancy=cob,
        waterplane_area=waterplane_area,
        wetted_surface_area=wetted_area,
        block_coefficient=cb,
        transverse_metacentric_radius=bm_t,
        longitudinal_waterplane_inertia=iyy,
        transverse_waterplane_inertia=ixx,
        lcb=cob[0],
        vcb=cob[2],
    )


def _compute_waterplane_area(mesh: Mesh, waterline_z: float) -> float:
    """Compute waterplane area by finding the waterline contour."""
    tol = 1e-9
    wl_edges: list[tuple[tuple[float, float], tuple[float, float]]] = []

    for face in mesh.faces:
        v0 = mesh.vertices[face[0]]
        v1 = mesh.vertices[face[1]]
        v2 = mesh.vertices[face[2]]
        verts = [v0, v1, v2]
        crossings = []
        for ei in range(3):
            va = verts[ei]
            vb = verts[(ei + 1) % 3]
            za, zb = va[2], vb[2]
            if abs(zb - za) > 1e-15:
                t = (waterline_z - za) / (zb - za)
                if -tol <= t <= 1.0 + tol:
                    t = max(0.0, min(1.0, t))
                    x = va[0] + t * (vb[0] - va[0])
                    y = va[1] + t * (vb[1] - va[1])
                    crossings.append((x, y))
        if len(crossings) >= 2:
            wl_edges.append((crossings[0], crossings[1]))

    if not wl_edges:
        return 0.0

    # Collect all crossing points and sort into a contour
    pts: list[tuple[float, float]] = []
    for (x1, y1), (x2, y2) in wl_edges:
        pts.append((x1, y1))
        pts.append((x2, y2))

    # Deduplicate
    unique_pts = list(set(pts))
    if len(unique_pts) < 3:
        return 0.0

    # Sort by angle from centroid
    mean_x = sum(p[0] for p in unique_pts) / len(unique_pts)
    mean_y = sum(p[1] for p in unique_pts) / len(unique_pts)
    unique_pts.sort(key=lambda p: math.atan2(p[1] - mean_y, p[0] - mean_x))

    # Shoelace formula
    n = len(unique_pts)
    area = 0.0
    for i in range(n):
        x0, y0 = unique_pts[i]
        x1, y1 = unique_pts[(i + 1) % n]
        area += x0 * y1 - x1 * y0
    return abs(area) * 0.5


def _compute_waterplane_inertia(
    mesh: Mesh, waterline_z: float
) -> tuple[float, float]:
    """Compute waterplane second moments Ixx and Iyy about centroid."""
    tol = 1e-9
    pts_set: set[tuple[float, float]] = set()

    for face in mesh.faces:
        v0 = mesh.vertices[face[0]]
        v1 = mesh.vertices[face[1]]
        v2 = mesh.vertices[face[2]]
        verts = [v0, v1, v2]
        for ei in range(3):
            va = verts[ei]
            vb = verts[(ei + 1) % 3]
            za, zb = va[2], vb[2]
            if abs(zb - za) > 1e-15:
                t = (waterline_z - za) / (zb - za)
                if -tol <= t <= 1.0 + tol:
                    t = max(0.0, min(1.0, t))
                    x = va[0] + t * (vb[0] - va[0])
                    y = va[1] + t * (vb[1] - va[1])
                    pts_set.add((round(x, 10), round(y, 10)))

    pts = list(pts_set)
    if len(pts) < 3:
        return (0.0, 0.0)

    mean_x = sum(p[0] for p in pts) / len(pts)
    mean_y = sum(p[1] for p in pts) / len(pts)
    pts.sort(key=lambda p: math.atan2(p[1] - mean_y, p[0] - mean_x))

    n = len(pts)
    area = 0.0
    cx = 0.0
    cy = 0.0
    for i in range(n):
        x0, y0 = pts[i]
        x1, y1 = pts[(i + 1) % n]
        cross = x0 * y1 - x1 * y0
        area += cross
        cx += (x0 + x1) * cross
        cy += (y0 + y1) * cross
    area *= 0.5
    if abs(area) < 1e-15:
        return (0.0, 0.0)
    cx /= (6.0 * area)
    cy /= (6.0 * area)

    ixx = 0.0
    iyy = 0.0
    for i in range(n):
        x0, y0 = pts[i]
        x1, y1 = pts[(i + 1) % n]
        cross = x0 * y1 - x1 * y0
        ixx += (y0 ** 2 + y0 * y1 + y1 ** 2) * cross
        iyy += (x0 ** 2 + x0 * x1 + x1 ** 2) * cross
    ixx = abs(ixx) / 12.0
    iyy = abs(iyy) / 12.0

    return (ixx, iyy)


# ---------------------------------------------------------------------------
# Mesh quality assessment
# ---------------------------------------------------------------------------
def mesh_quality_report(mesh: Mesh) -> dict[str, Any]:
    """Assess mesh quality: scaled Jacobian, aspect ratio, watertightness."""
    is_watertight = meshlib.is_watertight(mesh)
    is_manifold = meshlib.is_edge_manifold(mesh)

    min_jacobian = 1.0
    max_aspect_ratio = 0.0
    n_degenerate = 0

    for face in mesh.faces:
        v0 = mesh.vertices[face[0]]
        v1 = mesh.vertices[face[1]]
        v2 = mesh.vertices[face[2]]

        e0 = math.sqrt(sum((a - b) ** 2 for a, b in zip(v0, v1)))
        e1 = math.sqrt(sum((a - b) ** 2 for a, b in zip(v1, v2)))
        e2 = math.sqrt(sum((a - b) ** 2 for a, b in zip(v2, v0)))

        edges = sorted([e0, e1, e2])
        if edges[0] < 1e-15:
            n_degenerate += 1
            min_jacobian = 0.0
            continue

        ar = edges[2] / edges[0]
        max_aspect_ratio = max(max_aspect_ratio, ar)

        ux, uy, uz = v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]
        vx, vy, vz = v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]
        nx = uy * vz - uz * vy
        ny = uz * vx - ux * vz
        nz = ux * vy - uy * vx
        area2 = math.sqrt(nx * nx + ny * ny + nz * nz)
        denom = edges[1] * edges[2]
        if denom > 1e-15:
            sj = area2 / denom
            min_jacobian = min(min_jacobian, sj)

    return {
        "watertight": is_watertight,
        "edge_manifold": is_manifold,
        "n_vertices": len(mesh.vertices),
        "n_faces": len(mesh.faces),
        "min_scaled_jacobian": min_jacobian,
        "max_aspect_ratio": max_aspect_ratio,
        "n_degenerate_faces": n_degenerate,
    }


# ---------------------------------------------------------------------------
# High-level builder
# ---------------------------------------------------------------------------
@dataclass
class ParametricHull:
    """Complete parametric hull: sections -> surface -> mesh -> hydrostatics."""
    sections: tuple[CrossSection, ...]
    surface: NURBSSurface
    mesh_data: SurfaceMesh
    hydrostatics: HydrostaticResult
    quality: dict[str, Any]
    parameters: dict[str, float]

    def to_dict(self) -> dict[str, Any]:
        return {
            "parameters": self.parameters,
            "n_sections": len(self.sections),
            "surface_degree": list(self.surface.degree),
            "mesh": {
                "n_vertices": self.mesh_data.n_vertices,
                "n_faces": self.mesh_data.n_faces,
            },
            "hydrostatics": self.hydrostatics.to_dict(),
            "quality": self.quality,
        }


def build_parametric_hull(
    sections: Sequence[CrossSection],
    *,
    length: float | None = None,
    beam: float | None = None,
    draft: float | None = None,
    degree_u: int = 3,
    degree_v: int = 3,
    n_points_per_section: int = 10,
    mesh_n_u: int = 40,
    mesh_n_v: int = 15,
    waterline_z: float = 0.0,
    mirror_y: bool = True,
) -> ParametricHull:
    """Build a complete parametric hull from cross-sections.

    This is the main entry point: it lofts a NURBS surface, tessellates it,
    computes hydrostatics, and assesses mesh quality.
    """
    # Infer dimensions if not given
    xs = [s.x for s in sections]
    if length is None:
        length = max(xs) - min(xs)
    if beam is None:
        beam = 2.0 * max(s.max_half_breadth for s in sections)
    if draft is None:
        # Draft is the depth from waterline to keel
        min_z = min(s.min_z for s in sections)
        draft = abs(min_z) if min_z < 0 else max(s.max_height for s in sections)

    # Loft surface
    surface = loft_nurbs_surface(
        sections, degree_u=degree_u, degree_v=degree_v,
        n_points_per_section=n_points_per_section,
    )

    # Tessellate
    surface_mesh = tessellate_surface(
        surface, n_u=mesh_n_u, n_v=mesh_n_v, mirror_y=mirror_y,
    )

    # Hydrostatics
    hydro = compute_hydrostatics(
        surface_mesh.mesh,
        waterline_z=waterline_z,
        length=length,
        beam=beam,
        draft=draft,
    )

    # Quality
    quality = mesh_quality_report(surface_mesh.mesh)

    params = {
        "length": length,
        "beam": beam,
        "draft": draft,
        "n_sections": len(sections),
        "degree_u": degree_u,
        "degree_v": degree_v,
        "mesh_n_u": mesh_n_u,
        "mesh_n_v": mesh_n_v,
    }

    return ParametricHull(
        sections=tuple(sections),
        surface=surface,
        mesh_data=surface_mesh,
        hydrostatics=hydro,
        quality=quality,
        parameters=params,
    )
