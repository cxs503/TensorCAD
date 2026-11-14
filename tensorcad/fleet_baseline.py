"""Baseline hull fleet: from published IGES geometry to a watertight solid.

The platform's geometry library carried six *placeholder* hulls (304 triangles
each, two byte-identical).  The gate in :mod:`tensorcad.fleet_gate` can only
*reject* such files; replacing them needs a producer that turns the workshop's
published geometry into the closed, oriented, volume-verified solid a solver
pre-processor accepts.  This module is that producer, and it is the whole chain:

``IGES patches`` -> ``point samples`` -> ``half-breadth field b(x, z)`` ->
``waterline slices`` -> ``closed triangle mesh`` -> ``hydrostatics`` ->
``verdicts against the publisher's own numbers``

Why a field and not a patch-by-patch tessellation
-------------------------------------------------
The workshop files model the hull as dozens of independent B-spline patches
(55 for KCS, 38 for KVLCC2) with no guarantee that neighbouring patch boundaries
share vertices.  Tessellating each patch on its own parametric grid therefore
leaves T-vertices and cracks along every seam, i.e. a mesh that is *not*
watertight -- exactly the property the deliverable is graded on.  Collapsing the
patch set to the half-breadth field ``b(x, z) = max |y|`` first removes seams
entirely: topology is then built from structured waterline slices, so every edge
is shared by exactly two triangles *by construction* and no repair/welding step
can silently distort the geometry.

The cost is explicit and measured, not hidden: ``b(x, z)`` is single-valued in
``y`` (true for a displacement hull, bulbs included) and is reconstructed from
samples on an ``(x, z)`` grid whose ``max`` binning biases the surface outwards
by O(cell^2 / radius).  That bias is bounded by the reported
``mesh_vs_field_volume_percent`` cross-check (mesh divergence theorem versus an
independent trapezoid integration of the same field).

Units and frame
---------------
Input IGES is in metres (global section unit flag 6), right-handed, ``x`` aft to
forward, ``z`` up from the keel.  The mesh keeps that frame and is capped at the
design waterline, so it is a *hull below the waterline*: displacement volume,
waterplane area and LCB are comparable with the publisher's table, and nothing
above the waterline is claimed.

Only the Python standard library is required.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import struct
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

from . import mesh as meshlib
from .geometry import Mesh, MeshBuilder
from .iges import BSplineSurface, read_iges, unit_scale

SCHEMA = "tensorcad.fleet_baseline/1.0"

#: Default reconstruction resolution: field cell sizes in metres, vertices per
#: waterline slice, and the target point density on the IGES surface.
DEFAULT_DX = 0.5
DEFAULT_DZ = 0.1
DEFAULT_RING_POINTS = 240
DEFAULT_SAMPLES_PER_METRE = 8.0  # ~0.125 m spacing on the surface

#: Volume acceptance against the publisher's displacement (the gate's tolerance).
VOLUME_TOLERANCE = 0.01

VERIFIED = "VERIFIED"
WITHHELD = "WITHHELD"


# --------------------------------------------------------------------------- #
# sampling
# --------------------------------------------------------------------------- #
def _net_lengths(surface: BSplineSurface) -> tuple[float, float]:
    """Polyline length of the control net in each parametric direction.

    This bounds the patch's own extent (B-spline convex-hull property) and gives
    per-direction grid counts instead of one isotropic count, which matters for
    the workshop patches: they are long and thin (a hull side strip is ~30 m in
    one direction and ~2 m in the other), so an isotropic grid either wastes
    millions of evaluations or under-samples one direction.
    """
    count_u, count_v = surface.count_u, surface.count_v
    control = surface.control
    length_u = 0.0
    for j in range(count_v):
        row = control[j * count_u:(j + 1) * count_u]
        length_u += sum(math.dist(row[i], row[i + 1]) for i in range(count_u - 1))
    length_u /= max(count_v, 1)
    length_v = 0.0
    for i in range(count_u):
        length_v += sum(
            math.dist(control[j * count_u + i], control[(j + 1) * count_u + i])
            for j in range(count_v - 1)
        )
    length_v /= max(count_u, 1)
    return length_u, length_v


def sample_surface(
    surface: BSplineSurface,
    spacing: float = 1.0 / DEFAULT_SAMPLES_PER_METRE,
    max_axis: int = 200,
    min_axis: int = 4,
) -> list[tuple[float, float, float]]:
    """Sample ``surface`` on an (u, v) grid fine enough for ``spacing`` metres."""
    spacing = max(spacing, 1e-6)
    length_u, length_v = _net_lengths(surface)
    count_u = int(math.ceil(length_u / spacing)) + 1
    count_v = int(math.ceil(length_v / spacing)) + 1
    count_u = max(min_axis, min(max_axis, count_u))
    count_v = max(min_axis, min(max_axis, count_v))
    return surface.sample(count_u, count_v)


def sample_patches(
    surfaces: Sequence[BSplineSurface],
    spacing: float = 1.0 / DEFAULT_SAMPLES_PER_METRE,
) -> list[tuple[float, float, float]]:
    points: list[tuple[float, float, float]] = []
    for surface in surfaces:
        points.extend(sample_surface(surface, spacing))
    return points


# --------------------------------------------------------------------------- #
# half-breadth field
# --------------------------------------------------------------------------- #
@dataclass
class HalfBreadthField:
    """``b[i][k] = max |y|`` of the hull surface at ``x = x0 + i*dx``, ``z = k*dz``."""

    x0: float
    dx: float
    dz: float
    b: list[list[float]]
    filled_nodes: int = 0
    empty_nodes: int = 0
    interpolated_nodes: int = 0

    @property
    def nx(self) -> int:
        return len(self.b)

    @property
    def nz(self) -> int:
        return len(self.b[0]) if self.b else 0

    def x(self, i: int) -> float:
        return self.x0 + i * self.dx

    def z(self, k: int) -> float:
        return k * self.dz

    def cover(self) -> float:
        total = self.filled_nodes + self.empty_nodes
        return 0.0 if total == 0 else self.filled_nodes / total

    def to_dict(self) -> dict[str, Any]:
        levels_with_data = [
            k for k in range(self.nz) if any(self.b[i][k] > 0.0 for i in range(self.nx))
        ]
        return {
            "x0": self.x0,
            "dx": self.dx,
            "dz": self.dz,
            "nx": self.nx,
            "nz": self.nz,
            "filled_nodes": self.filled_nodes,
            "empty_nodes": self.empty_nodes,
            "interpolated_nodes": self.interpolated_nodes,
            "cover": round(self.cover(), 6),
            "max_half_breadth_m": round(max(max(row) for row in self.b), 6),
            "levels_with_data": len(levels_with_data),
            "max_z_with_data_m": round(self.z(levels_with_data[-1]) if levels_with_data else 0.0, 6),
        }


def field_from_points(
    points: Iterable[tuple[float, float, float]],
    draft: float,
    dx: float = DEFAULT_DX,
    dz: float = DEFAULT_DZ,
    min_x: float | None = None,
    max_x: float | None = None,
) -> HalfBreadthField:
    """Accumulate the half-breadth field from surface samples below ``draft``."""
    kept = [(x, y, z) for x, y, z in points if 0.0 <= z <= draft and abs(y) > 1e-9]
    if not kept:
        raise ValueError("no surface samples below the design waterline")
    x_lo = min_x if min_x is not None else min(x for x, _, _ in kept)
    x_hi = max_x if max_x is not None else max(x for x, _, _ in kept)
    nx = int(math.floor((x_hi - x_lo) / dx)) + 1
    nz = int(math.floor(draft / dz)) + 1
    grid = [[0.0] * nz for _ in range(nx)]
    seen = [[False] * nz for _ in range(nx)]
    for x, y, z in kept:
        i = int(round((x - x_lo) / dx))
        k = int(round(z / dz))
        if not (0 <= i < nx and 0 <= k < nz):
            continue
        magnitude = abs(y)
        if magnitude > grid[i][k]:
            grid[i][k] = magnitude
        seen[i][k] = True
    filled = sum(1 for i in range(nx) for k in range(nz) if seen[i][k])
    empty = nx * nz - filled
    interpolated = _fill_interior_holes(grid, seen)
    return HalfBreadthField(
        x0=x_lo,
        dx=dx,
        dz=dz,
        b=grid,
        filled_nodes=filled,
        empty_nodes=empty,
        interpolated_nodes=interpolated,
    )


def _fill_interior_holes(grid: list[list[float]], seen: list[list[bool]]) -> int:
    """Linearly bridge empty nodes that have filled neighbours on both sides in ``x``.

    Such a node is *inside* the hull footprint and was only missed by the
    sampling grid; leaving it at zero would punch a spurious hole through the
    surface.  Nodes at the ends of a run are genuine hull boundaries and stay
    zero.
    """
    filled = 0
    nx = len(grid)
    for k in range(len(grid[0])):
        i = 0
        while i < nx:
            if seen[i][k]:
                i += 1
                continue
            start = i
            while i < nx and not seen[i][k]:
                i += 1
            end = i
            if start == 0 or end >= nx:
                continue
            left, right = start - 1, end
            span = right - left
            for j in range(start, end):
                weight = (j - left) / span
                grid[j][k] = grid[left][k] * (1.0 - weight) + grid[right][k] * weight
                filled += 1
    return filled


# --------------------------------------------------------------------------- #
# waterline slices -> closed mesh
# --------------------------------------------------------------------------- #
def slice_profile(field: HalfBreadthField, k: int) -> tuple[list[float], list[float]] | None:
    """Starboard profile ``(x, b)`` of the waterline at level ``k``.

    End points are interpolated to ``b = 0`` so a slice closes where the hull
    actually ends at that immersion level instead of snapping to the grid.
    """
    values = [field.b[i][k] for i in range(field.nx)]
    nonzero = [i for i, value in enumerate(values) if value > 0.0]
    if not nonzero:
        return None
    first, last = nonzero[0], nonzero[-1]
    xs: list[float] = []
    bs: list[float] = []
    if first > 0:
        left = values[first - 1]
        weight = values[first] / (values[first] - left) if values[first] != left else 0.0
        xs.append(field.x(first) - weight * field.dx)
        bs.append(0.0)
    for i in range(first, last + 1):
        xs.append(field.x(i))
        bs.append(values[i])
    if last < field.nx - 1:
        right = values[last + 1]
        weight = values[last] / (values[last] - right) if values[last] != right else 0.0
        xs.append(field.x(last) + weight * field.dx)
        bs.append(0.0)
    return xs, bs


def closed_ring(
    xs: Sequence[float], bs: Sequence[float], z: float, count: int
) -> list[tuple[float, float, float]]:
    """Resample the closed waterline loop (starboard + mirrored port) to ``count`` points."""
    starboard = list(zip(xs, bs))
    loop = list(starboard) + [(x, -y) for x, y in reversed(starboard[1:-1])]
    lengths = [
        math.hypot(loop[(i + 1) % len(loop)][0] - loop[i][0], loop[(i + 1) % len(loop)][1] - loop[i][1])
        for i in range(len(loop))
    ]
    total = sum(lengths)
    if total <= 0.0:
        return [(loop[0][0], loop[0][1], z)] * count
    step = total / count
    ring: list[tuple[float, float, float]] = []
    segment = 0
    consumed = 0.0
    for index in range(count):
        target = index * step
        while segment < len(lengths) - 1 and consumed + lengths[segment] < target:
            consumed += lengths[segment]
            segment += 1
        a, b = loop[segment], loop[(segment + 1) % len(loop)]
        length = lengths[segment]
        if length <= 0.0:
            x, y = a
        else:
            ratio = (target - consumed) / length
            x = a[0] + (b[0] - a[0]) * ratio
            y = a[1] + (b[1] - a[1]) * ratio
        ring.append((x, y, z))
    return ring


def ring_area(ring: Sequence[tuple[float, float, float]]) -> float:
    area = 0.0
    for index in range(len(ring)):
        x0, y0, _ = ring[index]
        x1, y1, _ = ring[(index + 1) % len(ring)]
        area += x0 * y1 - x1 * y0
    return 0.5 * abs(area)


def build_solid(
    field: HalfBreadthField,
    ring_points: int = DEFAULT_RING_POINTS,
    draft: float | None = None,
) -> tuple[Mesh, dict[str, Any]]:
    """Loft the waterline slices of ``field`` into a closed, outward-oriented mesh."""
    draft = draft if draft is not None else field.z(field.nz - 1)
    rings: list[list[tuple[float, float, float]]] = []
    slice_areas: list[float] = []
    for k in range(field.nz):
        if field.z(k) > draft + 1e-9:
            break
        profile = slice_profile(field, k)
        if profile is None:
            continue
        ring = closed_ring(profile[0], profile[1], field.z(k), ring_points)
        area = ring_area(ring)
        if area <= 1e-9 and rings:
            continue
        rings.append(ring)
        slice_areas.append(area)
    if len(rings) < 2:
        raise ValueError("field produced fewer than two usable waterline slices")

    builder = MeshBuilder(
        units="m",
        coordinate_frame="x aft->fwd, y positive to starboard, z up from keel (m)",
    )
    indices = [[builder.add_vertex(point) for point in ring] for ring in rings]
    for lower, upper in zip(indices, indices[1:]):
        for index in range(ring_points):
            nxt = (index + 1) % ring_points
            builder.add_quad(lower[index], lower[nxt], upper[nxt], upper[index])
    caps = 0
    for position, ring, area in ((0, indices[0], slice_areas[0]), (-1, indices[-1], slice_areas[-1])):
        if area <= 1e-6:
            continue  # pinched end already closes the surface
        points = rings[position]
        centre = (
            sum(p[0] for p in points) / ring_points,
            sum(p[1] for p in points) / ring_points,
            points[0][2],
        )
        centre_index = builder.add_vertex(centre)
        for index in range(ring_points):
            nxt = (index + 1) % ring_points
            builder.add_triangle(ring[index], ring[nxt], centre_index)
        caps += 1
    raw = builder.build()
    oriented = meshlib.orient_outward(raw)
    info = {
        "slices": len(rings),
        "slices_capped": caps,
        "first_slice_area_m2": slice_areas[0],
        "last_slice_area_m2": slice_areas[-1],
        "max_slice_area_m2": max(slice_areas),
        "degenerate_skipped": builder.degenerate_skipped,
        "orientation_flipped": getattr(oriented, "flipped", None),
        "z_min": rings[0][0][2],
        "z_max": rings[-1][0][2],
    }
    return oriented.mesh, info


# --------------------------------------------------------------------------- #
# measurement
# --------------------------------------------------------------------------- #
def measure_mesh(mesh: Mesh, draft: float) -> dict[str, Any]:
    """Volume, wetted area and topology facts of a closed hull-below-waterline mesh."""
    volume = meshlib.total_signed_volume(mesh)
    bbox = mesh.bounding_box()
    return {
        "volume_m3": volume,
        "wetted_surface_area_m2": meshlib.surface_area(mesh),
        "vertex_count": len(mesh.vertices),
        "triangle_count": len(mesh.faces),
        "zero_area_faces": len(meshlib.zero_area_faces(mesh)),
        "is_watertight": meshlib.is_watertight(mesh),
        "is_edge_manifold": meshlib.is_edge_manifold(mesh),
        "boundary_edges": len(meshlib.boundary_edges(mesh)),
        "bbox_min": list(bbox.minimum),
        "bbox_max": list(bbox.maximum),
        "extent_m": [hi - lo for lo, hi in zip(bbox.minimum, bbox.maximum)],
        "draft": draft,
    }


def field_integrals(field: HalfBreadthField) -> dict[str, Any]:
    """Independent piecewise-trapezoid hydrostatics straight from the field.

    This is the cross-check on the mesh: it integrates the *same* field without
    going through triangles, so an orientation flip or a topology bug shows up as
    a disagreement between the two volumes.
    """
    nz, dx, dz = field.nz, field.dx, field.dz
    section_areas: list[float] = []
    moments: list[float] = []
    for i in range(field.nx):
        column = [field.b[i][k] for k in range(nz)]
        area = 2.0 * dz * (0.5 * column[0] + sum(column[1:-1]) + 0.5 * column[-1])
        section_areas.append(area)
        moments.append(field.x(i) * area)
    volume = 0.0
    moment = 0.0
    for i in range(field.nx - 1):
        volume += 0.5 * (section_areas[i] + section_areas[i + 1]) * dx
        moment += 0.5 * (moments[i] + moments[i + 1]) * dx
    top = nz - 1
    waterplane = 2.0 * sum(
        0.5 * (field.b[i][top] + field.b[i + 1][top]) * dx for i in range(field.nx - 1)
    )
    return {
        "volume_m3": volume,
        "lcb_x_m": moment / volume if volume else float("nan"),
        "max_section_area_m2": max(section_areas),
        "waterplane_area_m2": waterplane,
        "half_breadth_max_m": max(max(row) for row in field.b),
    }


# --------------------------------------------------------------------------- #
# one ship
# --------------------------------------------------------------------------- #
@dataclass
class ShipSpec:
    """A baseline hull plus the provenance needed to grade it."""

    name: str
    draft: float
    length_pp: float
    beam: float
    published_volume_m3: float | None = None
    published_cb: float | None = None
    published_lcb_percent: float | None = None
    #: model-space -> full-scale factor of the *geometry file*.  Workshop archives
    #: ship some hulls at model scale (ONRT 1:48.9355 as published, JBC in a
    #: model frame) and some at full scale (KVLCC2).  ``None`` means "infer it by
    #: normalising the point cloud's x-extent to the published L_PP"; that inference
    #: is recorded in the report and turns the volume gate into a block-coefficient
    #: (shape) gate, so it is never silent.
    geometry_scale: float | None = 1.0
    scale_basis: str = ""
    published_source: str = ""
    source_url: str = ""
    publisher: str = ""
    license: str = ""
    iges_sha256: str = ""
    iges_path: str = ""
    notes: str = ""

    def to_dict(self) -> dict[str, Any]:
        return dict(self.__dict__)


def partition_surfaces(
    surfaces: Sequence[BSplineSurface],
) -> tuple[list[BSplineSurface], list[dict[str, Any]]]:
    """Split patches into ``(mesh_these, quarantined_with_reason)``.

    A patch whose control-point layout is ambiguous (two candidate layouts within a
    factor of four in control-net roughness) may be a mis-read net, so it is not meshed:
    guessing would put a fictitious surface into the half-breadth field, and the field's
    interior-hole fill already covers a small missing patch from its neighbours. The
    census travels with the hull report so the deliverable never hides a hole.
    """
    usable: list[BSplineSurface] = []
    quarantined: list[dict[str, Any]] = []
    for index, surface in enumerate(surfaces):
        if surface.layout_ambiguous:
            quarantined.append(
                {
                    "patch_index": index,
                    "layout_chosen": surface.layout,
                    "reason": "ambiguous control-point layout (two candidates within 4x)",
                }
            )
            continue
        usable.append(surface)
    return usable, quarantined


def build_ship(
    spec: ShipSpec,
    dx: float = DEFAULT_DX,
    dz: float = DEFAULT_DZ,
    ring_points: int = DEFAULT_RING_POINTS,
    samples_per_metre: float = DEFAULT_SAMPLES_PER_METRE,
) -> tuple[Mesh, dict[str, Any]]:
    """Run the whole IGES -> mesh -> hydrostatics chain for one hull."""
    model = read_iges(spec.iges_path)
    surfaces = model.surfaces()
    if not surfaces:
        raise ValueError(f"{spec.name}: no B-spline surfaces in {spec.iges_path}")
    # Scale the samples to metres before any hydrostatics: units flag (metres,
    # millimetres, ...) times the model-space scale.  The read of the flag is now
    # Hollerith-positional (see iges.hollerith_fields); the previous content scan read
    # the JBC file -- millimetres throughout -- as metres, i.e. a 1000x hull.
    factor = unit_scale(model.unit_flag) * model.model_scale
    usable, quarantined = partition_surfaces(surfaces)
    if not usable:
        raise ValueError(f"{spec.name}: every surface was quarantined")
    raw_points = list(
        sample_patches(usable, spacing=1.0 / max(samples_per_metre, 1e-6))
    )
    geometry_scale = spec.geometry_scale
    scale_basis = spec.scale_basis or "geometry file is full scale"
    if geometry_scale is None:
        # No published model scale for this hull: normalise the point cloud's
        # x-extent to the publisher's L_PP.  The factor is recorded, and the
        # volume gate below degrades to a block-coefficient (shape) gate because
        # one length has been fitted.
        xs_raw = [point[0] for point in raw_points]
        extent_raw = (max(xs_raw) - min(xs_raw)) * factor if xs_raw else float("nan")
        if not (extent_raw and extent_raw == extent_raw and extent_raw > 0.0):
            raise ValueError(f"{spec.name}: cannot infer a scale from an empty point cloud")
        geometry_scale = spec.length_pp / extent_raw
        scale_basis = (
            f"inferred (no published model scale): published L_PP {spec.length_pp:.3f} m / "
            f"geometric x-extent {extent_raw:.4f} m = {geometry_scale:.5f}"
        )
    factor = factor * geometry_scale
    points = [(x * factor, y * factor, z * factor) for x, y, z in raw_points]
    # Frame guard.  Workshop files do not share one origin: KCS/KVLCC2/JBC put the
    # keel at z=0, but the KVLCC2 patch set has the *waterline* at z=0 (keel at
    # z=-20.8 m).  Detect that layout from the data (every point at or below a plane
    # that sits within one draft of the deepest point) and shift the hull up so the
    # keel is z=0, which is the frame the waterline cut expects.  A translation does
    # not change the volume, only the placement of the cut.
    z_min = min(z for _, _, z in points)
    z_max = max(z for _, _, z in points)
    frame_shift_z = 0.0
    if z_min < -0.5 * spec.draft and -0.05 * spec.draft <= z_max <= 0.05 * spec.draft:
        frame_shift_z = -z_max
        points = [(x, y, z + frame_shift_z) for x, y, z in points]
    field = field_from_points(points, spec.draft, dx=dx, dz=dz)
    mesh, mesh_info = build_solid(field, ring_points=ring_points, draft=spec.draft)
    mesh_stats = measure_mesh(mesh, spec.draft)
    integrals = field_integrals(field)

    top_profile = slice_profile(field, field.nz - 1)
    length_wl = (top_profile[0][-1] - top_profile[0][0]) if top_profile else float("nan")

    volume = mesh_stats["volume_m3"]
    cb = (
        volume / (spec.length_pp * spec.beam * spec.draft)
        if spec.length_pp and spec.beam and spec.draft
        else float("nan")
    )
    cm = (
        integrals["max_section_area_m2"] / (spec.beam * spec.draft)
        if spec.beam and spec.draft
        else float("nan")
    )
    cwp = (
        integrals["waterplane_area_m2"] / (length_wl * spec.beam)
        if length_wl and spec.beam
        else float("nan")
    )
    lcb_percent = (
        (integrals["lcb_x_m"] - 0.5 * spec.length_pp) / spec.length_pp * 100.0
        if spec.length_pp
        else float("nan")
    )

    volume_error = None
    volume_ok: bool | None = None
    volume_gate_basis = (
        "shape (Cb): scale was fitted to the published L_PP, so a volume match is a "
        "block-coefficient match, not an independent check of the absolute scale"
        if spec.geometry_scale is None
        else "absolute displacement at the publisher's own scale"
    )
    if spec.published_volume_m3:
        volume_error = (volume - spec.published_volume_m3) / spec.published_volume_m3 * 100.0
        volume_ok = abs(volume_error) <= 100.0 * VOLUME_TOLERANCE

    # Scale guard: a model-scale file or a mis-read units flag shows up as a geometric
    # extent that is not the published length. This is checked *before* the volume gate so
    # that "1 m of model" cannot masquerade as a 99.99% volume error.
    extent_x = mesh_stats.get("extent_m", [float("nan")] * 3)[0]
    scale_ratio = extent_x / spec.length_pp if spec.length_pp and extent_x else float("nan")
    scale_consistent = bool(
        spec.length_pp and extent_x and 0.85 <= scale_ratio <= 1.20
    )

    topology_ok = (
        mesh_stats["is_watertight"]
        and mesh_stats["is_edge_manifold"]
        and mesh_stats["zero_area_faces"] == 0
        and volume > 0.0
    )
    verdict = (
        VERIFIED if (topology_ok and volume_ok is True and scale_consistent) else WITHHELD
    )

    reasons: list[str] = []
    if not mesh_stats["is_watertight"]:
        reasons.append("mesh is not watertight")
    if not mesh_stats["is_edge_manifold"]:
        reasons.append("mesh is not edge-manifold")
    if mesh_stats["zero_area_faces"]:
        reasons.append(f"{mesh_stats['zero_area_faces']} zero-area faces")
    if volume_ok is False:
        reasons.append(
            f"volume {volume:.0f} m3 is {volume_error:+.2f}% off the published "
            f"{spec.published_volume_m3:.0f} m3 (tolerance {VOLUME_TOLERANCE:.0%})"
        )
    if spec.published_volume_m3 is None:
        reasons.append("no published displacement available for this hull")
    if not scale_consistent:
        reasons.append(
            f"geometry extent {extent_x:.3f} m vs published L_PP {spec.length_pp:.1f} m "
            f"(ratio {scale_ratio:.4f}): model-scale bytes or a units error, so a volume "
            f"comparison against the full-scale displacement is not meaningful"
            + (
                f"; scaling by {1.0 / scale_ratio:.3f}^3 would give "
                f"{volume * (1.0 / scale_ratio) ** 3:.1f} m3"
                if scale_ratio and scale_ratio == scale_ratio and scale_ratio != 0
                else ""
            )
        )

    report = {
        "schema": SCHEMA,
        "volume_gate_basis": volume_gate_basis,
        "spec": spec.to_dict(),
        "source": {
            "iges_stats": model.stats(),
            "units": model.units,
            "unit_flag": model.unit_flag,
            "unit_name": model.unit_name,
            "unit_source": model.unit_source,
            "model_space_scale": model.model_scale,
            "geometry_scale_to_full": geometry_scale,
            "scale_basis": scale_basis,
            "scale_inferred": spec.geometry_scale is None,
            "frame_shift_z": frame_shift_z,
            "unit_scale_to_metres": factor,
            "unit_warnings": list(model.unit_warnings),
            "patches_total": len(surfaces),
            "patches_used": len(usable),
            "patches_quarantined": quarantined,
            "surface_count": len(surfaces),
            "sample_count": len(points),
            "sampling_spacing_m": 1.0 / max(samples_per_metre, 1e-6),
        },
        "field": field.to_dict(),
        "mesh_info": mesh_info,
        "mesh": mesh_stats,
        "integrals": integrals,
        "derived": {
            "volume_m3": volume,
            "field_volume_m3": integrals["volume_m3"],
            "mesh_vs_field_volume_percent": (
                (volume - integrals["volume_m3"]) / integrals["volume_m3"] * 100.0
                if integrals["volume_m3"]
                else float("nan")
            ),
            "cb": cb,
            "cm": cm,
            "cwp": cwp,
            "lcb_percent_of_lpp": lcb_percent,
            "length_waterline_m": length_wl,
            "published_volume_rel_error_percent": volume_error,
        },
        "verdict": verdict,
        "reasons": reasons,
        "checks": {
            "watertight": mesh_stats["is_watertight"],
            "edge_manifold": mesh_stats["is_edge_manifold"],
            "positive_volume": volume > 0.0,
            "volume_within_tolerance": volume_ok,
            "scale_consistent_with_published_lpp": scale_consistent,
            "extent_over_lpp": scale_ratio if scale_ratio == scale_ratio else None,
            "lwl_vs_lpp_percent": (
                (length_wl - spec.length_pp) / spec.length_pp * 100.0
                if length_wl and spec.length_pp
                else None
            ),
        },
    }
    return mesh, report


# --------------------------------------------------------------------------- #
# STL export (the placeholder replacement format)
# --------------------------------------------------------------------------- #
def write_binary_stl(mesh: Mesh, path: str | Path, name: str = "hull") -> dict[str, Any]:
    """Write a binary STL and return its size/digest."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    header = f"tensorcad fleet_baseline {name}: watertight hull below WL, metres".encode()[:80]
    payload = bytearray(header + b" " * (80 - len(header)))
    payload += struct.pack("<I", len(mesh.faces))
    for face in mesh.faces:
        a, b, c = (mesh.vertices[index] for index in face)
        payload += struct.pack("<12fH", *meshlib.face_normal(mesh, face), *a, *b, *c, 0)
    target.write_bytes(bytes(payload))
    return {
        "path": str(target),
        "bytes": len(payload),
        "triangles": len(mesh.faces),
        "sha256": hashlib.sha256(payload).hexdigest(),
    }


# --------------------------------------------------------------------------- #
# fleet driver
# --------------------------------------------------------------------------- #
def specs_from_json(path: str | Path, *, asset_registry: object | None = None) -> list[ShipSpec]:
    """Read a fleet spec file, resolving each geometry pointer through the asset library.

    Spec files written before the durable store existed point at mission work trees
    (``/root/tp-agent/work/m18488/...``); those paths die with the mission. Pointers are
    therefore resolved by :func:`tensorcad.assetlib.consumer.resolve_pointer`, which
    prefers an ``asset://`` URI, then the path itself, then the recorded SHA-256, so a
    benchmark can be re-run from the store without editing the spec. A pointer that
    resolves nowhere is passed through unchanged: the reader reports the absence, which
    is more honest than failing here with a different message.
    """
    from tensorcad.assetlib.consumer import resolve_pointer

    payload = json.loads(Path(path).read_text())
    ships = payload["ships"] if isinstance(payload, dict) else payload
    specs = []
    for entry in ships:
        iges_path = resolve_pointer(
            entry["iges_path"],
            entry.get("iges_sha256", ""),
            registry=asset_registry,  # type: ignore[arg-type]
        )
        specs.append(
            ShipSpec(
                name=entry["name"],
                draft=float(entry["draft"]),
                length_pp=float(entry["length_pp"]),
                beam=float(entry["beam"]),
                published_volume_m3=entry.get("published_volume_m3"),
                published_cb=entry.get("published_cb"),
                published_lcb_percent=entry.get("published_lcb_percent"),
                published_source=entry.get("published_source", ""),
                source_url=entry.get("source_url", ""),
                publisher=entry.get("publisher", ""),
                license=entry.get("license", ""),
                iges_sha256=entry.get("iges_sha256", ""),
                iges_path=iges_path,
                notes=entry.get("notes", ""),
            )
        )
    return specs


def build_fleet(
    specs: Sequence[ShipSpec],
    out_dir: str | Path,
    *,
    dx: float = DEFAULT_DX,
    dz: float = DEFAULT_DZ,
    ring_points: int = DEFAULT_RING_POINTS,
    samples_per_metre: float = DEFAULT_SAMPLES_PER_METRE,
    stl: bool = True,
) -> dict[str, Any]:
    """Build every hull in ``specs`` and write per-ship manifests.

    Failures are recorded per ship and do not abort the fleet: a hull that cannot
    be rebuilt is reported as ``WITHHELD`` with its exception, because silently
    dropping it is how the placeholder problem started.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    results: dict[str, Any] = {}
    summary: dict[str, Any] = {
        "schema": SCHEMA,
        "hulls": {},
        "verified": 0,
        "withheld": 0,
        "errors": 0,
    }
    for spec in specs:
        try:
            mesh, report = build_ship(
                spec,
                dx=dx,
                dz=dz,
                ring_points=ring_points,
                samples_per_metre=samples_per_metre,
            )
        except Exception as error:  # noqa: BLE001 - recorded, not swallowed
            summary["hulls"][spec.name] = {
                "verdict": WITHHELD,
                "reasons": [f"{type(error).__name__}: {error}"],
            }
            summary["errors"] += 1
            summary["withheld"] += 1
            results[spec.name] = summary["hulls"][spec.name]
            continue
        if stl:
            report["stl"] = write_binary_stl(mesh, out / "meshes" / f"{spec.name}_Hull.stl", spec.name)
        (out / f"{spec.name}_manifest.json").write_text(json.dumps(report, indent=2, sort_keys=False))
        summary["hulls"][spec.name] = {
            "verdict": report["verdict"],
            "volume_m3": report["derived"]["volume_m3"],
            "published_volume_m3": spec.published_volume_m3,
            "volume_rel_error_percent": report["derived"]["published_volume_rel_error_percent"],
            "watertight": report["checks"]["watertight"],
            "edge_manifold": report["checks"]["edge_manifold"],
            "triangles": report["mesh"]["triangle_count"],
            "reasons": report["reasons"],
            "manifest": f"{spec.name}_manifest.json",
        }
        if report["verdict"] == VERIFIED:
            summary["verified"] += 1
        else:
            summary["withheld"] += 1
        results[spec.name] = report
    (out / "fleet_summary.json").write_text(json.dumps(summary, indent=2))
    return summary


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m tensorcad.fleet_baseline",
        description="Build watertight baseline hull meshes from published IGES geometry.",
    )
    parser.add_argument("--specs", required=True, help="fleet spec JSON (ships[] with iges_path)")
    parser.add_argument("--out", default="artifacts/baseline_fleet", help="output directory")
    parser.add_argument("--dx", type=float, default=DEFAULT_DX)
    parser.add_argument("--dz", type=float, default=DEFAULT_DZ)
    parser.add_argument("--ring-points", type=int, default=DEFAULT_RING_POINTS)
    parser.add_argument("--samples-per-metre", type=float, default=DEFAULT_SAMPLES_PER_METRE)
    parser.add_argument("--no-stl", action="store_true")
    args = parser.parse_args(argv)

    specs = specs_from_json(args.specs)
    summary = build_fleet(
        specs,
        args.out,
        dx=args.dx,
        dz=args.dz,
        ring_points=args.ring_points,
        samples_per_metre=args.samples_per_metre,
        stl=not args.no_stl,
    )
    print(json.dumps(summary["hulls"], indent=2))
    return 0 if summary["errors"] == 0 else 1


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())