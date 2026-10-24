"""Hydrostatic integrals evaluated directly on a triangulated surface.

These are exact integrals *of the mesh*, computed with the divergence theorem, so
they inherit the mesh's discretisation error -- which is precisely why
:func:`tensorcad.hull.estimate_discretisation_error` reports it alongside them.
They are not a hull-form office calculation and nothing here claims stability,
resistance or class-rule compliance.

Conventions and assumptions
---------------------------
* +z is up and the waterline is a plane of constant ``z``.
* The mesh is the wetted surface closed by a planar waterline cap, and no face
  straddles the waterline; faces are classified as wetted or coplanar by their
  centroid. Slicing a mesh at an arbitrary waterline is out of scope.
* Section areas assume a symmetric section that is single-valued in ``z``
  (a trapezoidal integral of the half-breadth), which is what the generated
  sampled hull provides. A bulbous or reentrant section needs a real loft.
* Volume and centre integrals need a closed, outward-oriented surface; when the
  mesh is open the result is reported with ``basis="open surface"`` so a caller
  cannot mistake it for a trustworthy displacement.
"""
from __future__ import annotations

import itertools
import math
from collections import defaultdict
from dataclasses import dataclass

from . import mesh as meshlib
from .geometry import Mesh

WATERLINE_TOLERANCE = 1e-9


@dataclass(frozen=True)
class Hydrostatics:
    """Integrated hydrostatic properties of a surface below the waterline."""

    displacement_volume: float
    displacement_mass: float
    wetted_surface_area: float
    waterplane_area: float
    longitudinal_centre_buoyancy: float
    longitudinal_centre_flotation: float
    vertical_centre_buoyancy: float
    block_coefficient: float
    prismatic_coefficient: float
    maximum_section_area: float
    waterline_z: float
    water_density: float
    basis: str = "closed surface integrals"

    def to_dict(self) -> dict[str, object]:
        return {
            "displacement_volume_m3": _round(self.displacement_volume),
            "displacement_mass_kg": _round(self.displacement_mass, 3),
            "wetted_surface_area_m2": _round(self.wetted_surface_area),
            "waterplane_area_m2": _round(self.waterplane_area),
            "longitudinal_centre_buoyancy_m": _round(self.longitudinal_centre_buoyancy),
            "longitudinal_centre_flotation_m": _round(self.longitudinal_centre_flotation),
            "vertical_centre_buoyancy_m": _round(self.vertical_centre_buoyancy),
            "block_coefficient": _round(self.block_coefficient, 6),
            "prismatic_coefficient": _round(self.prismatic_coefficient, 6),
            "maximum_section_area_m2": _round(self.maximum_section_area),
            "waterline_z_m": _round(self.waterline_z),
            "water_density_kg_per_m3": _round(self.water_density, 3),
            "basis": self.basis,
        }


def _round(value: float, digits: int = 9) -> float:
    if not math.isfinite(value):
        return value
    return round(value, digits)


def compute_hydrostatics(
    mesh: Mesh,
    waterline_z: float = 0.0,
    *,
    water_density: float = 1025.0,
    closed: bool | None = None,
) -> Hydrostatics:
    """Integrate displacement, areas and centres from the surface triangles.

    ``water_density`` defaults to seawater (1025 kg/m^3); pass 1000.0 for fresh
    water. ``closed`` overrides the watertightness assumption for callers that
    have already validated it.
    """
    if water_density <= 0 or not math.isfinite(water_density):
        raise ValueError("water_density must be a finite positive number")
    if closed is None:
        closed = meshlib.is_watertight(mesh)

    volume_terms: list[float] = []
    moment_x: list[float] = []
    moment_z: list[float] = []
    wetted_areas: list[float] = []
    flotation_area: list[float] = []
    flotation_moment_x: list[float] = []

    for face in mesh.faces:
        tetra = _tetra_volume(mesh, face)
        centroid = meshlib.face_centroid(mesh, face)
        area = meshlib.triangle_area(mesh, face)
        volume_terms.append(tetra)
        # The tetrahedron is (origin, a, b, c), so its centroid is (a+b+c)/4: the
        # face centroid (a+b+c)/3 scaled by 3/4, not divided by 4.
        moment_x.append(tetra * 0.75 * centroid[0])
        moment_z.append(tetra * 0.75 * centroid[2])
        if centroid[2] < waterline_z - WATERLINE_TOLERANCE:
            wetted_areas.append(area)
            # Projected area of the face on the waterplane: summing it over the
            # wetted surface gives the waterplane area, because the flux of a
            # closed surface integrates to zero.
            normal = meshlib.face_normal(mesh, face)
            projected = 0.5 * normal[2]
            flotation_area.append(projected)
            flotation_moment_x.append(centroid[0] * projected)

    volume = math.fsum(volume_terms)
    wetted = math.fsum(wetted_areas)
    signed_flotation = math.fsum(flotation_area)
    waterplane_area = abs(signed_flotation)

    basis = "closed surface integrals" if closed else "open surface (volume not trustworthy)"
    if abs(volume) <= 1e-15:
        return Hydrostatics(
            displacement_volume=volume,
            displacement_mass=volume * water_density,
            wetted_surface_area=wetted,
            waterplane_area=waterplane_area,
            longitudinal_centre_buoyancy=float("nan"),
            longitudinal_centre_flotation=_safe_ratio(flotation_moment_x, signed_flotation),
            vertical_centre_buoyancy=float("nan"),
            block_coefficient=float("nan"),
            prismatic_coefficient=float("nan"),
            maximum_section_area=float("nan"),
            waterline_z=waterline_z,
            water_density=water_density,
            basis=basis,
        )

    box = mesh.bounding_box()
    length, beam = box.sizes[0], box.sizes[1]
    depth = box.sizes[2]
    maximum_section = max_section_area(mesh)
    block = volume / (length * beam * depth) if length * beam * depth else float("nan")
    prismatic = (
        volume / (length * maximum_section)
        if length and maximum_section and math.isfinite(maximum_section)
        else float("nan")
    )
    return Hydrostatics(
        displacement_volume=volume,
        displacement_mass=volume * water_density,
        wetted_surface_area=wetted,
        waterplane_area=waterplane_area,
        longitudinal_centre_buoyancy=math.fsum(moment_x) / volume,
        longitudinal_centre_flotation=_safe_ratio(flotation_moment_x, signed_flotation),
        vertical_centre_buoyancy=math.fsum(moment_z) / volume,
        block_coefficient=block,
        prismatic_coefficient=prismatic,
        maximum_section_area=maximum_section,
        waterline_z=waterline_z,
        water_density=water_density,
        basis=basis,
    )


def _safe_ratio(weights: list[float], total: float) -> float:
    if abs(total) <= 1e-15:
        return float("nan")
    return math.fsum(weights) / total


def _tetra_volume(mesh: Mesh, face: tuple[int, int, int]) -> float:
    """Signed volume of the tetrahedron spanned by the origin and the face."""
    ax, ay, az = mesh.vertices[face[0]]
    bx, by, bz = mesh.vertices[face[1]]
    cx, cy, cz = mesh.vertices[face[2]]
    return (
        ax * (by * cz - bz * cy)
        + ay * (bz * cx - bx * cz)
        + az * (bx * cy - by * cx)
    ) / 6.0


def stations(mesh: Mesh) -> tuple[float, ...]:
    """Distinct longitudinal positions present in the mesh, smallest first."""
    return tuple(sorted({round(point[0], 9) for point in mesh.vertices}))


def section_area_at(mesh: Mesh, x: float, tolerance: float | None = None) -> float:
    """Wetted cross-section area at a station, from the half-breadth profile."""
    if tolerance is None:
        tolerance = max(1e-9, abs(x) * 1e-9)
    profile: dict[float, float] = defaultdict(float)
    for px, py, pz in mesh.vertices:
        if abs(px - x) <= tolerance:
            profile[pz] = max(profile[pz], abs(py))
    levels = sorted(profile)
    if len(levels) < 2:
        return float("nan")
    area = 0.0
    for lower, upper in itertools.pairwise(levels):
        area += 0.5 * (profile[lower] + profile[upper]) * (upper - lower)
    return 2.0 * area


def max_section_area(mesh: Mesh) -> float:
    """Largest wetted cross-section area over the stations the mesh actually has."""
    areas = [section_area_at(mesh, x) for x in stations(mesh)]
    finite = [area for area in areas if math.isfinite(area)]
    return max(finite) if finite else float("nan")
