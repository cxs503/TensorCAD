"""Extract a watertight triangle mesh from a continuous implicit field.

The neural field is the canonical geometry; a mesh is a *derived* artifact for
viewers, solvers and exchange formats such as STL/OBJ. This module therefore
does one job and does it strictly: sample the field on a padded grid, march
tetrahedra over a Kuhn (Freudenthal) subdivision, and return a ``tensorcad.Mesh``
that the existing :mod:`tensorcad.validation` checks can accept or reject.

Why marching *tetrahedra* rather than marching cubes: the six-tetrahedron
subdivision of a cube is unambiguous, so the surface cannot develop the cracks
that force a marching-cubes implementation to carry case tables and asymptotic
deciders. Every tet face is shared by exactly two cells with the same diagonal,
so identical intersection points are keyed by the same integer pair and
vertices are shared by construction rather than by tolerance welding.

Face winding is fixed *locally*: each triangle knows which of its source tet
vertices lies inside and which lies outside, so the normal is compared against
the inside->outside direction and flipped when it points the wrong way. That is
one vectorised pass instead of a breadth-first repair over the whole mesh, which
keeps extraction near-linear in the number of triangles.

Only numpy is required, never torch: the caller passes a callable mapping an
``(n, 3)`` point array to an ``(n,)`` value array.
"""
from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import numpy.typing as npt

from .. import mesh as meshlib
from ..geometry import Mesh, MeshBuilder

#: Array annotations are written as ``npt.NDArray[Any]`` on purpose: older NumPy
#: stubs require the type arguments and newer ones make them optional, so a bare
#: ``np.ndarray`` or a dtype-pinned alias type-checks differently depending on
#: which NumPy the environment happens to resolve. Shape and dtype are checked at
#: runtime by the tests; the alias only has to mean the same thing everywhere.
Array = npt.NDArray[Any]

#: Corner bit layout of a cell: index ``4*x + 2*y + 1*z`` with each bit 0/1.
_CORNER_BITS = np.array(
    [[(index >> 2) & 1, (index >> 1) & 1, index & 1] for index in range(8)], dtype=np.int64
)
_CORNER_OFFSETS = _CORNER_BITS[:, 0] * 1 + _CORNER_BITS[:, 1] + _CORNER_BITS[:, 2]

#: Kuhn subdivision: the six tetrahedra sharing the 0-7 main diagonal.
_TETRAHEDRA = np.array(
    [[0, 4, 6, 7], [0, 4, 5, 7], [0, 2, 6, 7], [0, 2, 3, 7], [0, 1, 5, 7], [0, 1, 3, 7]],
    dtype=np.int64,
)


class ExtractionError(ValueError):
    """Raised when the sampled field cannot produce a usable surface."""


@dataclass(frozen=True)
class ExtractionReport:
    """Measured quality of an extracted mesh, written into the object metadata."""

    resolution: int | tuple[int, int, int]
    vertices: int
    triangles: int
    watertight: bool
    edge_manifold: bool
    euler_characteristic: int
    signed_volume: float
    surface_area: float
    smallest_triangle_area: float
    bounds_min: tuple[float, float, float]
    bounds_max: tuple[float, float, float]

    def to_dict(self) -> dict[str, object]:
        return {
            "resolution_cells": list(self.resolution)
            if isinstance(self.resolution, tuple)
            else self.resolution,
            "vertex_count": self.vertices,
            "triangle_count": self.triangles,
            "watertight": self.watertight,
            "edge_manifold": self.edge_manifold,
            "euler_characteristic": self.euler_characteristic,
            "signed_volume_m3": self.signed_volume,
            "surface_area_m2": self.surface_area,
            "smallest_triangle_area_m2": self.smallest_triangle_area,
            "bounds_min_m": list(self.bounds_min),
            "bounds_max_m": list(self.bounds_max),
        }


def _evaluate(
    field: Callable[[Array], Array], points: Array, chunk: int = 1 << 18
) -> Array:
    values = np.empty(len(points), dtype=np.float64)
    for start in range(0, len(points), chunk):
        block = points[start : start + chunk]
        block_values = np.asarray(field(block), dtype=np.float64).reshape(-1)
        if len(block_values) != len(block):
            raise ExtractionError(f"field returned {len(block_values)} values for {len(block)} points")
        values[start : start + chunk] = block_values
    if not np.all(np.isfinite(values)):
        raise ExtractionError("field produced non-finite values on the sampling grid")
    return values


def cells_per_axis(resolution: int | Sequence[int]) -> tuple[int, int, int]:
    """Normalise a scalar or per-axis resolution into three cell counts.

    An axisymmetric or slender body does not want a cube of cells: a 8.6:1 hull
    sampled at 64^3 wastes 98 % of the grid on empty space outside the body and
    still leaves only 64 cells along the length. Per-axis resolution keeps the
    cell size near-isotropic at the same total cost.
    """
    if isinstance(resolution, (Sequence, np.ndarray)):
        cells = tuple(int(value) for value in resolution)
        if len(cells) != 3:
            raise ExtractionError(f"resolution needs 1 or 3 cell counts, got {len(cells)}")
        values: tuple[int, int, int] = (cells[0], cells[1], cells[2])
    else:
        scalar = int(resolution)
        values = (scalar, scalar, scalar)
    if min(values) < 2:
        raise ExtractionError(f"resolution must be at least 2 cells per axis, got {values}")
    return values


def sample_grid(
    field: Callable[[Array], Array],
    minimum: Sequence[float],
    maximum: Sequence[float],
    resolution: int | Sequence[int],
    margin: int = 3,
    *,
    isovalue: float = 0.0,
) -> tuple[Array, Array, tuple[int, int, int]]:
    """Sample ``field`` on a padded grid; returns ``(points, values, counts)``.

    ``resolution`` is cells per axis, either one number or three.
    """
    low = np.asarray(minimum, dtype=np.float64).reshape(3)
    high = np.asarray(maximum, dtype=np.float64).reshape(3)
    if np.any(high <= low):
        raise ExtractionError("bounds must have positive extent on every axis")
    cells = np.asarray(cells_per_axis(resolution), dtype=np.float64)
    if margin < 1:
        raise ExtractionError("margin must be at least 1 cell")
    spacing = (high - low) / cells
    origin = low - spacing * margin
    counts: tuple[int, int, int] = (
        int(cells[0]) + 1 + 2 * margin,
        int(cells[1]) + 1 + 2 * margin,
        int(cells[2]) + 1 + 2 * margin,
    )
    axes = [origin[axis] + spacing[axis] * np.arange(counts[axis]) for axis in range(3)]
    points = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3)
    values = _evaluate(field, points)
    if not np.any(values < isovalue):
        raise ExtractionError("the field has no inside (negative) samples in the given bounds")
    if np.all(values < isovalue):
        raise ExtractionError("the field is inside everywhere in the given bounds; enlarge them")
    # Treat an exact isovalue sample as inside so interpolation never divides by
    # a zero difference and never lands exactly on a grid point.
    values = np.where(np.abs(values) < 1e-12, -1e-12, values)
    return points, values, counts


@dataclass
class _Slab:
    edge_lo: Array
    edge_hi: Array
    triangles: Array
    inside_ref: Array
    outside_ref: Array


def _march_slab(
    values: Array,
    corners: Array,
    strides: tuple[int, int, int],
    start: int,
    stop: int,
    isovalue: float,
) -> _Slab:
    """March the cells between grid planes ``start`` and ``stop``."""
    slab = corners[start:stop, :-1, :-1]
    offsets = _CORNER_BITS[:, 0] * strides[0] + _CORNER_BITS[:, 1] * strides[1] + _CORNER_BITS[:, 2]
    cell_corners = (slab[..., None] + offsets).reshape(-1, 8).astype(np.int64)
    cell_values = values[cell_corners]
    del slab

    tet_corners = cell_corners[:, _TETRAHEDRA].reshape(-1, 4)
    tet_values = cell_values[:, _TETRAHEDRA].reshape(-1, 4)
    del cell_corners, cell_values

    inside = tet_values < isovalue
    inside_count = inside.sum(axis=1)
    active = np.nonzero((inside_count > 0) & (inside_count < 4))[0]
    if len(active) == 0:
        empty = np.empty(0, dtype=np.int64)
        return _Slab(empty, empty, np.empty((0, 3), dtype=np.int64), empty, empty)

    counts_active = inside_count[active]
    ordering = np.argsort(inside[active], axis=1, kind="stable")
    active_corners = tet_corners[active]

    edge_lo: list[Array] = []
    edge_hi: list[Array] = []
    triangles: list[Array] = []
    inside_refs: list[Array] = []
    outside_refs: list[Array] = []
    edge_total = 0

    def emit(
        group: Array,
        slot_a: Array,
        slot_b: Array,
        inside_slot: Array,
        outside_slot: Array,
        local_triangles: tuple[tuple[int, int, int], ...],
    ) -> None:
        nonlocal edge_total
        row_index = group[:, None]
        first = active_corners[row_index, slot_a]
        second = active_corners[row_index, slot_b]
        edge_lo.append(np.minimum(first, second).reshape(-1))
        edge_hi.append(np.maximum(first, second).reshape(-1))
        per_row = slot_a.shape[1]
        base = edge_total + np.arange(len(group))[:, None] * per_row + np.arange(per_row)[None, :]
        in_ref = active_corners[group, inside_slot]
        out_ref = active_corners[group, outside_slot]
        for triangle in local_triangles:
            triangles.append(base[:, list(triangle)])
        # Triangle blocks are concatenated group-by-group (all rows of the first
        # triangle, then all rows of the second), so the references must be
        # stacked the same way -- not interleaved per row.
        inside_refs.append(np.concatenate([in_ref] * len(local_triangles)))
        outside_refs.append(np.concatenate([out_ref] * len(local_triangles)))
        edge_total += len(group) * per_row

    single = np.nonzero(counts_active == 1)[0]
    if len(single) > 0:
        row = ordering[single]
        emit(
            single,
            np.repeat(row[:, 3:4], 3, axis=1),
            row[:, :3],
            row[:, 3],
            row[:, 0],
            ((0, 1, 2),),
        )

    triple = np.nonzero(counts_active == 3)[0]
    if len(triple) > 0:
        row = ordering[triple]
        emit(
            triple,
            row[:, 1:4],
            np.repeat(row[:, 0:1], 3, axis=1),
            row[:, 1],
            row[:, 0],
            ((0, 1, 2),),
        )

    pair = np.nonzero(counts_active == 2)[0]
    if len(pair) > 0:
        row = ordering[pair]
        i0, i1 = row[:, 2:3], row[:, 3:4]
        o0, o1 = row[:, 0:1], row[:, 1:2]
        emit(
            pair,
            np.concatenate((i0, i0, i1, i1), axis=1),
            np.concatenate((o0, o1, o1, o0), axis=1),
            row[:, 2],
            row[:, 0],
            ((0, 1, 2), (0, 2, 3)),
        )

    if not edge_lo:
        empty = np.empty(0, dtype=np.int64)
        return _Slab(empty, empty, np.empty((0, 3), dtype=np.int64), empty, empty)
    return _Slab(
        np.concatenate(edge_lo),
        np.concatenate(edge_hi),
        np.concatenate(triangles).reshape(-1, 3),
        np.concatenate(inside_refs),
        np.concatenate(outside_refs),
    )


def extract_surface(
    field: Callable[[Array], Array],
    minimum: Sequence[float],
    maximum: Sequence[float],
    resolution: int | Sequence[int],
    *,
    isovalue: float = 0.0,
    margin: int = 3,
    units: str = "m",
    coordinate_frame: str = "object_local",
    cell_chunk: int = 60_000,
) -> Mesh:
    """March the zero level set of ``field`` and return a closed triangle mesh.

    ``minimum``/``maximum`` bound the region of interest; sampling is extended by
    ``margin`` cells on every side so the surface closes inside the grid instead
    of hitting its boundary. ``resolution`` is cells per axis (one number or three).
    """
    points, values, counts = sample_grid(field, minimum, maximum, resolution, margin, isovalue=isovalue)
    corners = np.arange(len(values)).reshape(counts)
    strides = (counts[1] * counts[2], counts[2], 1)

    slabs: list[_Slab] = []
    edge_total = 0
    triangle_total = 0
    cells_per_slab = max(1, cell_chunk // (counts[1] * counts[2]))
    for start in range(0, counts[0] - 1, cells_per_slab):
        stop = min(start + cells_per_slab, counts[0] - 1)
        slab = _march_slab(values, corners, strides, start, stop, isovalue)
        if len(slab.edge_lo) == 0:
            continue
        slabs.append(
            _Slab(
                slab.edge_lo,
                slab.edge_hi,
                slab.triangles + edge_total,
                slab.inside_ref,
                slab.outside_ref,
            )
        )
        edge_total += len(slab.edge_lo)
        triangle_total += len(slab.triangles)

    if not slabs:
        raise ExtractionError("marching tetrahedra produced no triangles")

    edge_lo = np.concatenate([slab.edge_lo for slab in slabs])
    edge_hi = np.concatenate([slab.edge_hi for slab in slabs])
    triangles = np.concatenate([slab.triangles for slab in slabs])
    inside_ref = np.concatenate([slab.inside_ref for slab in slabs])
    outside_ref = np.concatenate([slab.outside_ref for slab in slabs])
    assert len(triangles) == triangle_total == len(inside_ref) == len(outside_ref)

    keys = edge_lo.astype(np.int64) * np.int64(len(values)) + edge_hi
    unique_keys, inverse = np.unique(keys, return_inverse=True)
    unique_lo = unique_keys // len(values)
    unique_hi = unique_keys % len(values)

    endpoint_a = points[unique_lo]
    endpoint_b = points[unique_hi]
    value_a = values[unique_lo]
    value_b = values[unique_hi]
    denominator = value_b - value_a
    safe = np.abs(denominator) >= 1e-300
    fraction = np.where(safe, (isovalue - value_a) / np.where(safe, denominator, 1.0), 0.5)
    # Keep crossings strictly inside their edge: a crossing exactly at a grid
    # sample would make two different edges produce the same point, which collapses
    # a triangle and can open the surface when the degenerate face is dropped.
    # (This does not remove slivers: two crossings can still land very close in
    # space when the surface passes near a grid vertex. See docs/neural-geometry-format.md.)
    fraction = np.clip(fraction, 1e-6, 1.0 - 1e-6)
    intersections = endpoint_a + fraction[:, None] * (endpoint_b - endpoint_a)

    faces = inverse[triangles].astype(np.int64)
    # Orient each triangle away from the solid using the inside->outside pairing
    # recorded at emission time; no global repair pass is needed afterwards.
    p0, p1, p2 = intersections[faces[:, 0]], intersections[faces[:, 1]], intersections[faces[:, 2]]
    normals = np.cross(p1 - p0, p2 - p0)
    direction = points[outside_ref] - points[inside_ref]
    flip = np.einsum("ij,ij->i", normals, direction) < 0.0
    faces[flip] = faces[flip][:, [0, 2, 1]]

    builder = MeshBuilder(units=units, coordinate_frame=coordinate_frame)
    for x, y, z in intersections:
        builder.add_vertex((float(x), float(y), float(z)))
    for a, b, c in faces:
        builder.add_triangle(int(a), int(b), int(c))
    mesh = meshlib.remove_degenerate_faces(builder.build())
    if not mesh.faces:
        raise ExtractionError("every extracted triangle was degenerate")
    if meshlib.total_signed_volume(mesh) < 0.0:  # pragma: no cover - defensive
        mesh = meshlib.orient_outward(mesh).mesh
    return mesh


def summarise(
    mesh: Mesh,
    resolution: int | Sequence[int],
    *,
    expected_units: str | None = None,
) -> ExtractionReport:
    """Measure an extracted mesh with the same checks used for CAD output."""
    cells = cells_per_axis(resolution)
    edges = meshlib.edge_face_counts(mesh)
    areas = meshlib.triangle_areas(mesh)
    box = mesh.bounding_box()
    if expected_units is not None and mesh.units != expected_units:
        raise ExtractionError(f"mesh units {mesh.units!r} do not match the object's {expected_units!r}")
    return ExtractionReport(
        # Isotropic callers keep the scalar they passed; anisotropic ones get the
        # per-axis tuple, so an existing reader of ``resolution`` is unaffected.
        resolution=cells if isinstance(resolution, (tuple, list)) else cells[0],
        vertices=len(mesh.vertices),
        triangles=len(mesh.faces),
        watertight=meshlib.is_watertight(mesh),
        edge_manifold=meshlib.is_edge_manifold(mesh),
        euler_characteristic=len(mesh.vertices) - len(edges) + len(mesh.faces),
        signed_volume=meshlib.total_signed_volume(mesh),
        surface_area=meshlib.surface_area(mesh),
        smallest_triangle_area=min(areas) if areas else 0.0,
        bounds_min=box.minimum,
        bounds_max=box.maximum,
    )
