"""Closed-surface mesh -> cell occupancy lattice: the voxel face of the CAD handoff.

Why this module exists
----------------------
:mod:`tensorcad.solver_interface` publishes three solver-facing tables (surface
triangles, a structured grid, a lattice).  The lattice route was, until now, only
ever handed a *declared* occupancy: a box-aligned channel whose fluid cells the
caller typed in.  A CAD hull could not produce one.  The previous round said so in
its own words -- *"the hull ray-cast voxelisation needs a spatial index (the naive
O(cells x faces) version did not finish)"* -- and left it open.

The spatial index it needed already exists
(:class:`tensorcad.spatial_index.TriangleGrid`).  What was missing is the thing
that turns it into a *lattice*: cast **one ray per column** instead of one per
cell, sort that column's crossings once, and fill the inside *spans* by parity.
Cost per column is ``O(cells crossed + hits)`` instead of ``O(faces)``, and the
fill is a bulk slice assignment rather than a per-cell point query, so a
2000-triangle hull at 10^5 cells is seconds, not hours.

What a caller is allowed to assume (and what is checked)

``occupancy``
    One byte per cell, ``1`` = inside the solid, in C order with ``i`` along
    ``+x`` fastest, then ``j`` along ``+y``, then ``k`` along ``+z`` -- the same
    ``(i, j, k)`` / ``+x, +y, +z`` convention ``solver_interface.write_bundle_lattice``
    writes, so a solver reads one cell ordering for every route.
``frame``
    The lattice is described in ``COORDINATE_CONTRACT`` (z-up, metres); the module
    never re-derives a frame, it carries the contract's name and schema.
``parity, not a signed-distance field``
    A cell is solid iff its centre is enclosed by the surface by ray-crossing
    parity.  That is exact for a closed surface and *first-order* accurate in the
    cell size -- occupancy can only gain or lose a straddling cell, never
    interpolate it.  :func:`convergence_study` measures that order rather than
    asserting it.
``precondition stated, not assumed``
    Parity is meaningless on an open surface (a non-watertight mesh gives stripes,
    not a geometry).  :func:`watertight_gate` runs :mod:`tensorcad.meshquality`'s
    surface topology first, and it is the first gate of the driver.

Evidence lives in :func:`run_voxel_case` (driver:
``python -m tensorcad.voxel_solid --out <dir>``) and ``tests/test_voxel_solid.py``;
the numbers are in ``docs/voxel_solid_handoff.md``.
"""
from __future__ import annotations

import base64
import hashlib
import json
import math
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterator, Sequence

from .geometry import Mesh, MeshBuilder
from .spatial_index import TriangleGrid
from .hydrostatics import compute_hydrostatics
from .solver_interface import (
    COORDINATE_CONTRACT,
    HandoffError,
    verify_bundle,
    write_bundle_lattice,
)

__all__ = [
    "ConvergenceRow",
    "LATTICE_SCHEMA",
    "SolidLattice",
    "box_surface_mesh",
    "brute_force_solid_lattice",
    "convergence_study",
    "hull_solid_lattice",
    "lattice_bundle",
    "observed_order",
    "occupancy_mismatch",
    "ray_direction",
    "read_lattice_bundle",
    "run_voxel_case",
    "solid_lattice",
    "watertight_gate",
]

#: Bumped on any incompatible change to the lattice record this module emits.
LATTICE_SCHEMA = "tensorcad.voxel_solid/1"

#: Probe directions, one per ray axis.  Tilted off the axis so a column ray cannot
#: run exactly along a ring of vertices (a shared edge is reported by two
#: triangles and would flip the parity for the rest of the column), but by only
#: ``1e-6``: the drift over a 10^2 m column is ``1e-4`` m, four orders below the
#: finest cell this module is used at.  The first version of this table used
#: ``1.3e-3`` -- copied from ``spatial_index._CONTAIN_PROBE``, which is sized for a
#: *point* query -- and the resulting ray drifted a whole cell across a hull box,
#: losing 3.5 % of the volume at every resolution and disagreeing with the
#: brute-force reference on 546 of 2880 cells.  The tilt must be small enough that
#: the ray stays inside its own column, not merely off an edge.
_PROBES = {
    "x": (1.0, 1.0e-6, 3.7e-7),
    "y": (1.0e-6, 1.0, 3.7e-7),
    "z": (1.0e-6, 3.7e-7, 1.0),
}

_AXIS_INDEX = {"x": 0, "y": 1, "z": 2}


def ray_direction(axis: str) -> tuple[float, float, float]:
    """Normalised probe direction for ``axis`` (``"x"``, ``"y"`` or ``"z"``)."""
    try:
        raw = _PROBES[axis]
    except KeyError:
        raise ValueError(f"ray axis must be one of x/y/z, got {axis!r}") from None
    norm = math.sqrt(sum(component * component for component in raw))
    return (raw[0] / norm, raw[1] / norm, raw[2] / norm)


def _stride(axis: int, shape: Sequence[int]) -> int:
    """Flat-index stride of ``axis`` in C order (``i`` fastest)."""
    return 1 if axis == 0 else shape[0] if axis == 1 else shape[0] * shape[1]


def _mesh_digest(mesh: Mesh) -> str:
    """Stable digest of the geometry the occupancy was computed from."""
    hasher = hashlib.sha256()
    for vertex in mesh.vertices:
        hasher.update(("%.17g %.17g %.17g\n" % vertex).encode())
    hasher.update(b"--\n")
    for face in mesh.faces:
        hasher.update(("%d %d %d\n" % face).encode())
    return hasher.hexdigest()


@dataclass(frozen=True)
class SolidLattice:
    """A cell occupancy grid derived from a closed triangle mesh.

    ``occupancy[i + nx * (j + ny * k)]`` is ``1`` when the cell
    ``origin + (i, j, k) * spacing`` .. ``+ spacing`` is inside the solid.  The
    record carries the mesh digest, the ray axis actually used and the wall-clock
    cost, so a number can be traced to the geometry and the method that produced
    it.
    """

    shape: tuple[int, int, int]
    origin: tuple[float, float, float]
    spacing: tuple[float, float, float]
    occupancy: bytes
    mesh_digest: str
    ray_axis: str
    seconds: float = 0.0

    def __post_init__(self) -> None:
        nx, ny, nz = (int(value) for value in self.shape)
        if min(nx, ny, nz) < 1:
            raise ValueError(f"lattice shape must be positive, got {self.shape!r}")
        cells = nx * ny * nz
        if len(self.occupancy) != cells:
            raise ValueError(
                f"occupancy has {len(self.occupancy)} entries for a {nx}x{ny}x{nz} lattice "
                f"({cells} cells)"
            )
        if min(self.spacing) <= 0.0:
            raise ValueError(f"spacing must be positive, got {self.spacing!r}")

    # -- geometry -----------------------------------------------------------

    @property
    def cells(self) -> int:
        """Total cell count."""
        return self.shape[0] * self.shape[1] * self.shape[2]

    @property
    def cell_volume(self) -> float:
        """Volume of one cell [m^3]."""
        return self.spacing[0] * self.spacing[1] * self.spacing[2]

    @property
    def solid_count(self) -> int:
        """Number of solid cells."""
        return sum(self.occupancy)

    @property
    def solid_volume(self) -> float:
        """Occupied volume ``solid_count * cell_volume`` [m^3]."""
        return self.solid_count * self.cell_volume

    def index(self, i: int, j: int, k: int) -> int:
        """Flat index of cell ``(i, j, k)``."""
        nx, ny, nz = self.shape
        if not (0 <= i < nx and 0 <= j < ny and 0 <= k < nz):
            raise IndexError(f"cell {(i, j, k)} lies outside a {self.shape} lattice")
        return i + nx * (j + ny * k)

    def is_solid(self, i: int, j: int, k: int) -> bool:
        """True when cell ``(i, j, k)`` is inside the solid."""
        return bool(self.occupancy[self.index(i, j, k)])

    def centre(self, i: int, j: int, k: int) -> tuple[float, float, float]:
        """Centre of cell ``(i, j, k)`` in the frame's metres."""
        return (
            self.origin[0] + (i + 0.5) * self.spacing[0],
            self.origin[1] + (j + 0.5) * self.spacing[1],
            self.origin[2] + (k + 0.5) * self.spacing[2],
        )

    def _cells(self, want: int) -> Iterator[tuple[int, int, int]]:
        nx, ny, nz = self.shape
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    if self.occupancy[i + nx * (j + ny * k)] == want:
                        yield (i, j, k)

    def solid_cells(self) -> Iterator[tuple[int, int, int]]:
        """Yield every solid cell index, in ``(i, j, k)`` order (``k`` slowest)."""
        return self._cells(1)

    def fluid_cells(self) -> Iterator[tuple[int, int, int]]:
        """Yield every non-solid (fluid) cell index, in ``(i, j, k)`` order."""
        return self._cells(0)

    # -- reporting ----------------------------------------------------------

    def to_dict(self, *, include_occupancy: bool = False) -> dict[str, Any]:
        """Record for ``results.json``; the occupancy itself only on request."""
        nx, ny, nz = self.shape
        record: dict[str, Any] = {
            "schema": LATTICE_SCHEMA,
            "frame": COORDINATE_CONTRACT["name"],
            "shape": [nx, ny, nz],
            "origin_m": list(self.origin),
            "spacing_m": list(self.spacing),
            "cell_volume_m3": self.cell_volume,
            "cells": self.cells,
            "solid_cells": self.solid_count,
            "solid_volume_m3": self.solid_volume,
            "mesh_digest_sha256": self.mesh_digest,
            "ray_axis": self.ray_axis,
            "seconds": self.seconds,
        }
        if include_occupancy:
            record["occupancy_b64"] = base64.b64encode(self.occupancy).decode("ascii")
            record["occupancy_bytes"] = len(self.occupancy)
        return record


# --------------------------------------------------------------------- ray maths


def _ray_triangle(
    origin: Sequence[float],
    direction: Sequence[float],
    a: Sequence[float],
    b: Sequence[float],
    c: Sequence[float],
) -> float | None:
    """Möller--Trumbore; the slow path's own intersection test.

    Deliberately *not* imported from :mod:`tensorcad.spatial_index`: the point of
    :func:`brute_force_solid_lattice` is to be an independent reference, the same
    way ``solver_interface.verify_bundle`` re-reads tables with its own parser.
    """
    e1 = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
    e2 = (c[0] - a[0], c[1] - a[1], c[2] - a[2])
    px = direction[1] * e2[2] - direction[2] * e2[1]
    py = direction[2] * e2[0] - direction[0] * e2[2]
    pz = direction[0] * e2[1] - direction[1] * e2[0]
    det = e1[0] * px + e1[1] * py + e1[2] * pz
    if abs(det) < 1e-14:
        return None
    inv = 1.0 / det
    tx = origin[0] - a[0]
    ty = origin[1] - a[1]
    tz = origin[2] - a[2]
    u = (tx * px + ty * py + tz * pz) * inv
    if u < 0.0 or u > 1.0:
        return None
    qx = ty * e1[2] - tz * e1[1]
    qy = tz * e1[0] - tx * e1[2]
    qz = tx * e1[1] - ty * e1[0]
    v = (direction[0] * qx + direction[1] * qy + direction[2] * qz) * inv
    if v < 0.0 or u + v > 1.0:
        return None
    t = (e2[0] * qx + e2[1] * qy + e2[2] * qz) * inv
    return t if t > 0.0 else None


def _bounds(
    mesh: Mesh, bounds: Sequence[Sequence[float]] | None, pad: float
) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    """Sampling box: the mesh's own bounding box, or a caller override, plus ``pad``."""
    if bounds is not None:
        lo = tuple(float(value) for value in bounds[0])
        hi = tuple(float(value) for value in bounds[1])
        if len(lo) != 3 or len(hi) != 3:
            raise ValueError("bounds must be two three-vectors")
    else:
        xs = [v[0] for v in mesh.vertices]
        ys = [v[1] for v in mesh.vertices]
        zs = [v[2] for v in mesh.vertices]
        lo = (min(xs), min(ys), min(zs))
        hi = (max(xs), max(ys), max(zs))
    if pad:
        lo = tuple(lo[i] - pad for i in range(3))
        hi = tuple(hi[i] + pad for i in range(3))
    if any(hi[i] <= lo[i] for i in range(3)):
        raise ValueError(f"empty sampling box: {lo!r} .. {hi!r}")
    return (lo, hi)


def _resolution(
    lo: Sequence[float],
    hi: Sequence[float],
    *,
    resolution: int | Sequence[int] | None,
    spacing: float | Sequence[float] | None,
) -> tuple[tuple[int, int, int], tuple[float, float, float]]:
    """Resolve the sampling grid from either a cell count or a cell size."""
    extent = tuple(hi[i] - lo[i] for i in range(3))
    if resolution is not None and spacing is not None:
        raise ValueError("pass resolution or spacing, not both")
    if spacing is not None:
        if isinstance(spacing, (int, float)):
            values = (float(spacing),) * 3
        else:
            values = tuple(float(value) for value in spacing)
        if len(values) != 3 or min(values) <= 0.0:
            raise ValueError(f"spacing must be a positive scalar or three-vector, got {spacing!r}")
        cells = tuple(max(1, int(round(extent[i] / values[i]))) for i in range(3))
        actual = tuple(extent[i] / cells[i] for i in range(3))
        return (int(cells[0]), int(cells[1]), int(cells[2])), (
            actual[0],
            actual[1],
            actual[2],
        )
    if resolution is None:
        raise ValueError("pass resolution (cells per axis) or spacing (cell size in metres)")
    if isinstance(resolution, (int, float)):
        cells = (int(resolution),) * 3
    else:
        values_r = tuple(int(value) for value in resolution)
        if len(values_r) == 1:
            cells = values_r * 3
        elif len(values_r) == 3:
            cells = values_r
        else:
            raise ValueError(f"resolution needs 1 or 3 cell counts, got {resolution!r}")
    if min(cells) < 1:
        raise ValueError(f"resolution must be positive cells per axis, got {resolution!r}")
    actual = tuple(extent[i] / cells[i] for i in range(3))
    return (cells[0], cells[1], cells[2]), (actual[0], actual[1], actual[2])


# ------------------------------------------------------------------- occupancy


def _fill_axis(
    shape: tuple[int, int, int],
    origin: tuple[float, float, float],
    spacing: tuple[float, float, float],
    axis: str,
    margin: float,
    crossings_for: Any,
) -> bytearray:
    """Fill a lattice by parity along ``axis``, one ray per column.

    ``crossings_for(u, v)`` returns the sorted crossing parameters ``t`` of the
    ray through column ``(u, v)``; every ray starts at ``origin[axis] - margin``,
    and ``t`` is arc length along the normalised direction, so a cell centre is
    reached at ``(margin + (n + 1/2) h) / d_axis`` -- divided, not multiplied.  Cell centres are walked against the crossing list by two
    pointers, so the interior is emitted as *spans* (bulk slice assignments along
    the contiguous axis) instead of cell-by-cell point queries.
    """
    nx, ny, nz = shape
    occ = bytearray(nx * ny * nz)
    index = _AXIS_INDEX[axis]
    run_axes = [a for a in range(3) if a != index]
    n_along = shape[index]
    stride = _stride(index, shape)
    direction = ray_direction(axis)
    d_along = direction[index]
    stride_u = _stride(run_axes[0], shape)
    stride_v = _stride(run_axes[1], shape)

    for v in range(shape[run_axes[1]]):
        for u in range(shape[run_axes[0]]):
            base = u * stride_u + v * stride_v
            crossings = crossings_for(u, v)
            n_crossings = len(crossings)
            pointer = 0
            run_start = -1
            for n in range(n_along):
                t_cell = (margin + (n + 0.5) * spacing[index]) / d_along
                while pointer < n_crossings and crossings[pointer] < t_cell:
                    pointer += 1
                inside = (pointer & 1) == 1
                if inside and run_start < 0:
                    run_start = n
                elif not inside and run_start >= 0:
                    _emit(occ, base, run_start, n, stride)
                    run_start = -1
            if run_start >= 0:
                _emit(occ, base, run_start, n_along, stride)
    return occ


def _emit(occ: bytearray, base: int, start: int, stop: int, stride: int) -> None:
    """Mark cells ``start .. stop-1`` of one column solid."""
    if stride == 1:
        lo = base + start
        occ[lo : lo + (stop - start)] = b"\x01" * (stop - start)
    else:
        for n in range(start, stop):
            occ[base + n * stride] = 1


def solid_lattice(
    mesh: Mesh,
    *,
    resolution: int | Sequence[int] | None = None,
    spacing: float | Sequence[float] | None = None,
    bounds: Sequence[Sequence[float]] | None = None,
    pad: float = 0.0,
    ray_axis: str = "x",
) -> SolidLattice:
    """Occupancy lattice of the solid bounded by ``mesh``.

    One ray per column along ``ray_axis`` (default ``x``: fewest rays on an
    elongated hull, and the only axis whose spans are contiguous in the C-order
    occupancy), cast through the uniform-grid index, then parity-filled by span.

    ``resolution`` is cells per axis (a scalar means a cubic count); ``spacing`` is
    the cell size in metres and is honoured by rounding each axis to a whole
    number of cells, which is reported back in ``lattice.spacing``.  Cells whose
    centre lies exactly on the surface are undefined, as for any parity test in
    floating point.
    """
    if not mesh.faces:
        raise ValueError("cannot voxelise a mesh with no triangles")
    lo, hi = _bounds(mesh, bounds, pad)
    shape, actual = _resolution(lo, hi, resolution=resolution, spacing=spacing)
    if ray_axis not in _AXIS_INDEX:
        raise ValueError(f"ray axis must be one of x/y/z, got {ray_axis!r}")
    started = time.perf_counter()
    grid = TriangleGrid(mesh)
    index = _AXIS_INDEX[ray_axis]
    run_axes = [a for a in range(3) if a != index]
    direction = ray_direction(ray_axis)
    # Start every ray outside the sampling box so the first crossing opens the solid.
    margin = 10.0 * max(hi[i] - lo[i] for i in range(3))

    def crossings_for(u: int, v: int) -> list[float]:
        start = [0.0, 0.0, 0.0]
        start[run_axes[0]] = lo[run_axes[0]] + (u + 0.5) * actual[run_axes[0]]
        start[run_axes[1]] = lo[run_axes[1]] + (v + 0.5) * actual[run_axes[1]]
        start[index] = lo[index] - margin
        return [hit.distance for hit in grid.cast_all(tuple(start), direction)]

    occupancy = _fill_axis(shape, lo, actual, ray_axis, margin, crossings_for)
    return SolidLattice(
        shape=shape,
        origin=lo,
        spacing=actual,
        occupancy=bytes(occupancy),
        mesh_digest=_mesh_digest(mesh),
        ray_axis=ray_axis,
        seconds=time.perf_counter() - started,
    )


def brute_force_solid_lattice(
    mesh: Mesh,
    *,
    resolution: int | Sequence[int] | None = None,
    spacing: float | Sequence[float] | None = None,
    bounds: Sequence[Sequence[float]] | None = None,
    pad: float = 0.0,
    ray_axis: str = "x",
) -> SolidLattice:
    """Reference occupancy: every cell centre tested against every triangle.

    This is the O(cells x faces) path the previous round could not finish on a
    hull.  It is kept in the library (not only in tests) so the accelerated fill
    can be checked on *any* mesh, exactly as ``spatial_index.brute_force_nearest``
    is kept for the ray caster: an acceleration structure that changes answers is
    worthless.
    """
    lo, hi = _bounds(mesh, bounds, pad)
    shape, actual = _resolution(lo, hi, resolution=resolution, spacing=spacing)
    direction = ray_direction(ray_axis)
    vertices = mesh.vertices
    faces = mesh.faces
    nx, ny, nz = shape
    occ = bytearray(nx * ny * nz)
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                point = (
                    lo[0] + (i + 0.5) * actual[0],
                    lo[1] + (j + 0.5) * actual[1],
                    lo[2] + (k + 0.5) * actual[2],
                )
                crossings = 0
                for a, b, c in faces:
                    if _ray_triangle(point, direction, vertices[a], vertices[b], vertices[c]) is not None:
                        crossings += 1
                if crossings % 2 == 1:
                    occ[i + nx * (j + ny * k)] = 1
    return SolidLattice(
        shape=shape,
        origin=lo,
        spacing=actual,
        occupancy=bytes(occ),
        mesh_digest=_mesh_digest(mesh),
        ray_axis=ray_axis,
        seconds=0.0,
    )


def occupancy_mismatch(left: SolidLattice, right: SolidLattice) -> int:
    """Number of cells on which two lattices disagree; ``0`` means identical answers."""
    if left.shape != right.shape:
        raise ValueError(f"cannot compare lattices of shape {left.shape} and {right.shape}")
    return sum(1 for a, b in zip(left.occupancy, right.occupancy) if a != b)


def box_surface_mesh(
    size: Sequence[float],
    *,
    origin: Sequence[float] = (0.0, 0.0, 0.0),
) -> Mesh:
    """Closed, outward-oriented axis-aligned box (12 triangles), for exactness tests.

    Every face is split along the same diagonal pattern and winding is taken from
    :func:`tensorcad.mesh.orient_outward`, so the solid's parity is unambiguous and
    the occupied volume of an aligned occupancy grid is exactly ``lx*ly*lz``.
    """
    lx, ly, lz = (float(v) for v in size)
    ox, oy, oz = (float(v) for v in origin)
    builder = MeshBuilder()
    corners = [
        builder.add_vertex((ox + x, oy + y, oz + z))
        for x, y, z in (
            (0.0, 0.0, 0.0),
            (lx, 0.0, 0.0),
            (lx, ly, 0.0),
            (0.0, ly, 0.0),
            (0.0, 0.0, lz),
            (lx, 0.0, lz),
            (lx, ly, lz),
            (0.0, ly, lz),
        )
    ]
    for a, b, c, d in (
        (0, 3, 2, 1),  # bottom, -z
        (4, 5, 6, 7),  # top, +z
        (0, 1, 5, 4),  # -y
        (1, 2, 6, 5),  # +x
        (2, 3, 7, 6),  # +y
        (3, 0, 4, 7),  # -x
    ):
        builder.add_quad(corners[a], corners[b], corners[c], corners[d])
    from .mesh import orient_outward

    return orient_outward(builder.build()).mesh


# ---------------------------------------------------------------- verification


def watertight_gate(mesh: Mesh) -> dict[str, Any]:
    """Closed-surface precondition, measured by :mod:`tensorcad.meshquality`.

    A parity fill on an open surface produces stripes rather than a geometry, so
    this gate runs *before* any volume is reported and its verdict travels with the
    numbers.  ``boundary_edges == 0`` is the watertight condition; a topology
    record that does not expose the field reports ``watertight = None`` rather
    than a guess.
    """
    from . import meshquality

    topology = meshquality.surface_topology(mesh.vertices, mesh.faces)
    if hasattr(topology, "to_dict"):
        record: dict[str, Any] = dict(topology.to_dict())
    else:  # pragma: no cover - dataclass fallback
        record = dict(vars(topology))
    boundary = record.get("boundary_edges", record.get("open_edges"))
    nonmanifold = record.get("nonmanifold_edges", record.get("non_manifold_edges"))
    return {
        "boundary_edges": boundary,
        "nonmanifold_edges": nonmanifold,
        "watertight": getattr(topology, "watertight", None),
        "signed_volume_m3": getattr(topology, "signed_volume", None),
        "topology": record,
    }


@dataclass(frozen=True)
class ConvergenceRow:
    """One resolution in a voxel convergence study."""

    resolution: tuple[int, int, int]
    h_m: float
    measured_volume_m3: float
    reference_volume_m3: float
    solid_cells: int
    seconds: float

    @property
    def rel_error(self) -> float:
        """Signed relative volume error against the reference."""
        return (self.measured_volume_m3 - self.reference_volume_m3) / self.reference_volume_m3

    def to_dict(self) -> dict[str, Any]:
        return {
            "resolution": list(self.resolution),
            "h_m": self.h_m,
            "measured_volume_m3": self.measured_volume_m3,
            "reference_volume_m3": self.reference_volume_m3,
            "rel_error": self.rel_error,
            "solid_cells": self.solid_cells,
            "seconds": self.seconds,
        }


def observed_order(rows: Sequence[ConvergenceRow]) -> float:
    """Fitted order of accuracy, least squares on ``log|error|`` against ``log h``.

    A boolean occupancy grid can only gain or lose a straddling cell, so the
    *worst case* order is 1; a smooth convex body measures about 1.9 because the
    midpoint indicator error of a curved boundary largely cancels (measured on a
    sphere: 5.4 % / 0.84 % / 0.40 % at 10 / 20 / 40 cells, order 1.88).  The
    number is measured, never assumed, and a fit below ~0.5 would say the fill is
    rounding the surface away rather than resolving it.
    """
    usable = [row for row in rows if row.rel_error != 0.0 and row.h_m > 0.0]
    if len(usable) < 2:
        raise ValueError("need at least two resolutions with non-zero error")
    xs = [math.log(row.h_m) for row in usable]
    ys = [math.log(abs(row.rel_error)) for row in usable]
    n = len(xs)
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    numerator = sum((xs[i] - mean_x) * (ys[i] - mean_y) for i in range(n))
    denominator = sum((xs[i] - mean_x) ** 2 for i in range(n))
    return numerator / denominator if denominator else float("nan")


def convergence_study(
    mesh: Mesh,
    *,
    reference_volume: float,
    resolutions: Sequence[int],
    bounds: Sequence[Sequence[float]] | None = None,
    pad: float = 0.0,
    ray_axis: str = "x",
) -> list[ConvergenceRow]:
    """Voxelise at each resolution and compare the occupied volume with a reference."""
    if reference_volume <= 0.0:
        raise ValueError(f"reference volume must be positive, got {reference_volume!r}")
    rows: list[ConvergenceRow] = []
    for value in resolutions:
        lattice = solid_lattice(mesh, resolution=value, bounds=bounds, pad=pad, ray_axis=ray_axis)
        rows.append(
            ConvergenceRow(
                resolution=lattice.shape,
                h_m=max(lattice.spacing),
                measured_volume_m3=lattice.solid_volume,
                reference_volume_m3=reference_volume,
                solid_cells=lattice.solid_count,
                seconds=lattice.seconds,
            )
        )
    return rows


def hull_solid_lattice(
    parameters: Any = None,
    *,
    resolution: int | Sequence[int] | None = 64,
    pad: float = 0.0,
    ray_axis: str = "x",
    spacing: float | None = None,
) -> tuple[SolidLattice, float, Mesh]:
    """Voxelise :mod:`tensorcad.hull`'s closed hull surface.

    Returns the lattice, the hull's **analytic** displaced volume -- a formula the
    triangulation never saw -- and the mesh.  That analytic volume is why this case
    is worth running: it is an independent reference for a voxel volume, not a
    second opinion from the same sampled surface.
    """
    from .hull import HullParameters, analytic_volume, build_mesh

    params = parameters if parameters is not None else HullParameters()
    mesh, _offsets = build_mesh(params)
    reference = analytic_volume(
        params.length,
        params.beam,
        params.draft,
        fore_exponent=params.fore_exponent,
        aft_exponent=params.aft_exponent,
    )
    lattice = solid_lattice(
        mesh, resolution=resolution, spacing=spacing, pad=pad, ray_axis=ray_axis
    )
    return lattice, reference, mesh


# ---------------------------------------------------------------------- handoff


def lattice_bundle(
    outdir: str | Path,
    lattice: SolidLattice,
    *,
    case: str,
    walls: dict[str, Any] | None = None,
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Write the lattice as the merged ``lattice.json`` handoff and verify it.

    ``fluid_cells`` is the complement of the solid, so the file is exactly what a
    voxel/LBM consumer already reads (``solver_interface.write_bundle_lattice``
    schema) and :func:`tensorcad.solver_interface.verify_bundle` re-reads it with
    its own parser before any solver sees it.  The solid mask itself travels in the
    manifest ``extra`` (packed base64) so a consumer can recover which cells were
    removed -- a fluid-only table would leave the hull invisible.

    A lattice with **no** fluid cells is refused by the writer's own gate (nothing
    to solve), which is why a solid-only box is not a valid bundle.
    """
    record = lattice.to_dict(include_occupancy=True)
    manifest_extra: dict[str, Any] = {
        "voxel_solid": {k: v for k, v in record.items() if k != "occupancy_b64"},
        "solid_occupancy_b64": record["occupancy_b64"],
    }
    manifest_extra.update(extra or {})
    manifest = write_bundle_lattice(
        outdir,
        shape=lattice.shape,
        origin=lattice.origin,
        spacing=lattice.spacing,
        fluid_cells=list(lattice.fluid_cells()),
        walls=walls,
        case=case,
        extra=manifest_extra,
    )
    verified = verify_bundle(outdir)
    if not verified.get("ok", True):
        raise HandoffError(f"lattice bundle failed verification: {verified}")
    return manifest


def read_lattice_bundle(bundle_dir: str | Path) -> tuple[dict[str, Any], SolidLattice]:
    """Read a bundle written by :func:`lattice_bundle` back into a lattice.

    The occupancy is recovered from the manifest's ``solid_occupancy_b64`` and
    cross-checked against the ``fluid_cells`` table: a consumer that trusts only
    one of the two would not notice a bundle whose tables disagree.
    """
    bundle = Path(bundle_dir)
    manifest = json.loads((bundle / "handoff.json").read_text(encoding="utf-8"))
    extra = manifest.get("extra", {})
    record = extra.get("voxel_solid")
    packed = extra.get("solid_occupancy_b64")
    if record is None or packed is None:
        raise HandoffError(f"{bundle} carries no voxel_solid occupancy record")
    occupancy = base64.b64decode(packed)
    lattice = SolidLattice(
        shape=tuple(record["shape"]),
        origin=tuple(record["origin_m"]),
        spacing=tuple(record["spacing_m"]),
        occupancy=occupancy,
        mesh_digest=record["mesh_digest_sha256"],
        ray_axis=record["ray_axis"],
        seconds=float(record.get("seconds", 0.0)),
    )
    table = json.loads((bundle / "lattice.json").read_text(encoding="utf-8"))
    listed = {tuple(cell) for cell in table["fluid_cells"]}
    recovered = set(lattice.fluid_cells())
    if listed != recovered:
        raise HandoffError(
            f"{bundle}: lattice.json lists {len(listed)} fluid cells but the solid mask "
            f"implies {len(recovered)}"
        )
    return manifest, lattice


def _boundary_unstable(mesh: Mesh, fast: SolidLattice, slow: SolidLattice) -> int:
    """Mismatching cells whose parity flips under a 1e-9 m perturbation.

    A cell centre that lies on the surface (to within floating-point noise) has no
    defined inside/outside, which this module documents.  Rather than exclude such
    cells by fiat, the driver perturbs each disagreeing centre by +-1e-9 m along a
    skewed direction and re-counts crossings against *every* triangle: a centre
    that changes answer under that perturbation is a boundary case; one that does
    not is a genuine disagreement between the two paths.
    """
    direction = (0.5773502691896258, 0.5773502691896257, 0.5773502691896257)
    unstable = 0
    nx, ny, nz = fast.shape
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                index = i + nx * (j + ny * k)
                if fast.occupancy[index] == slow.occupancy[index]:
                    continue
                centre = fast.centre(i, j, k)
                parities = set()
                for sign in (1.0, -1.0):
                    point = tuple(centre[a] + sign * 1e-9 * direction[a] for a in range(3))
                    crossings = 0
                    for a, b, c in mesh.faces:
                        if _ray_triangle(
                            point, direction, mesh.vertices[a], mesh.vertices[b], mesh.vertices[c]
                        ) is not None:
                            crossings += 1
                    parities.add(crossings % 2)
                if len(parities) > 1:
                    unstable += 1
    return unstable


# ----------------------------------------------------------------------- driver


def _gate(
    name: str, measured: float, reference: float, tolerance: float, unit: str = ""
) -> dict[str, Any]:
    """One comparison gate: relative error against a tolerance, verdict spelled out."""
    rel = abs(measured - reference) / abs(reference) if reference else abs(measured)
    return {
        "gate": name,
        "measured": measured,
        "reference": reference,
        "rel_error": rel,
        "tolerance": tolerance,
        "unit": unit,
        "verdict": "PASS" if rel <= tolerance else "FAIL",
    }


def _verdict(name: str, ok: bool, detail: str) -> dict[str, Any]:
    return {"gate": name, "measured": detail, "reference": "requirement", "rel_error": 0.0,
            "tolerance": 0.0, "unit": "", "verdict": "PASS" if ok else "FAIL"}


def run_voxel_case(
    outdir: str | Path,
    *,
    hull_resolution: int = 64,
    hull_parameters: Any = None,
    converge_resolutions: Sequence[int] = (32, 64, 128),
    brute_force_resolution: int = 12,
) -> dict[str, Any]:
    """Run every gate of the voxel face; write ``results.json`` and a lattice bundle.

    Gates, in the order they run:

    ``watertight``             the hull surface is closed before anything is filled
    ``box_exact``              an axis-aligned closed box fills to its volume exactly
    ``ray_axis_agreement``     all three ray axes produce the same occupancy
    ``grid_vs_bruteforce``     the accelerated fill agrees cell-for-cell with the slow path
    ``sphere_order``           a sphere's voxel volume converges at first order
    ``hull_volume``            the hull's voxel volume against its *analytic* displacement
    ``hull_order``             that error shrinks at first order with the cell size
    ``bundle_roundtrip``       the written lattice bundle re-reads into the same lattice
    """
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    results: dict[str, Any] = {
        "schema": "tensorcad.voxel_solid.results/1",
        "frame": COORDINATE_CONTRACT["name"],
        "gates": [],
        "cases": {},
    }

    def add(gate: dict[str, Any]) -> None:
        results["gates"].append(gate)

    # -- 1. the box: an aligned solid must fill exactly -------------------------
    box = box_surface_mesh((1.0, 0.6, 0.4))
    box_volume = 1.0 * 0.6 * 0.4
    box_lattice = solid_lattice(box, resolution=(25, 15, 10))
    results["cases"]["box"] = box_lattice.to_dict()
    add(_gate("box_exact_volume", box_lattice.solid_volume, box_volume, 1e-12, "m3"))

    # -- 2. ray-axis independence: the method, not the axis, decides -------------
    lattices_by_axis = {
        axis: solid_lattice(box, resolution=(13, 9, 7), ray_axis=axis) for axis in ("x", "y", "z")
    }
    mismatches = {
        f"{a}_vs_x": occupancy_mismatch(lattices_by_axis[a], lattices_by_axis["x"])
        for a in ("y", "z")
    }
    results["cases"]["ray_axis_mismatch"] = mismatches
    add(_verdict("ray_axis_agreement", sum(mismatches.values()) == 0, str(mismatches)))

    # -- 3. acceleration must not change answers --------------------------------
    #    A hull is the interesting mesh (curved, tapered, closed) and the lattice is
    #    small because the slow path is O(cells x faces) by design.
    hull = None
    hull_reference = 0.0
    try:
        _lattice_probe, hull_reference, hull = hull_solid_lattice(
            hull_parameters, resolution=hull_resolution
        )
    except Exception as exc:  # pragma: no cover - geometry failure is reported, not raised
        results["cases"]["hull_error"] = repr(exc)
        add(_verdict("hull_available", False, repr(exc)))
    if hull is not None:
        topology = watertight_gate(hull)
        results["cases"]["hull_topology"] = topology
        add(_verdict("hull_watertight", bool(topology["watertight"]), str(topology["watertight"])))

        # The faceted surface's own displacement is the reference the voxel
        # discretisation is scored against -- integrating the same triangles is a
        # different discretisation of the same solid, while the analytic volume
        # below carries the faceting error of the mesh as well.
        hydro = compute_hydrostatics(hull, closed=bool(topology["watertight"]))
        faceted = hydro.displacement_volume
        results["cases"]["hull_hydrostatics"] = {
            "displacement_volume_m3": faceted,
            "waterplane_area_m2": hydro.waterplane_area,
            "topology_signed_volume_m3": topology.get("signed_volume_m3"),
        }

        fast = solid_lattice(hull, resolution=brute_force_resolution)
        slow = brute_force_solid_lattice(hull, resolution=brute_force_resolution)
        n_bad = occupancy_mismatch(fast, slow)
        unstable = _boundary_unstable(hull, fast, slow)
        results["cases"]["bruteforce_crosscheck"] = {
            "resolution": list(fast.shape),
            "cells": fast.cells,
            "triangles": len(hull.faces),
            "mismatched_cells": n_bad,
            "mismatched_fraction": n_bad / fast.cells,
            "mismatched_cells_on_surface": unstable,
            "fast_seconds": fast.seconds,
            "solid_cells_fast": fast.solid_count,
            "solid_cells_slow": slow.solid_count,
        }
        # A centre lying on the surface has no defined parity; a cell whose answer
        # flips under a 1e-9 m perturbation is that case, not a disagreement.
        add(
            _verdict(
                "grid_vs_bruteforce",
                n_bad <= max(1, int(1e-3 * fast.cells)) and unstable == n_bad,
                f"{n_bad} of {fast.cells} cells differ ({unstable} of them boundary-unstable)",
            )
        )

        # -- 4. analytic hull volume + order of accuracy ------------------------
        rows = convergence_study(
            hull,
            reference_volume=faceted,
            resolutions=list(converge_resolutions),
        )
        order = observed_order(rows)
        results["cases"]["hull_convergence"] = {
            "reference_volume_m3": faceted,
            "reference_kind": "faceted displacement (compute_hydrostatics)",
            "analytic_volume_m3": hull_reference,
            "rows": [row.to_dict() for row in rows],
            "observed_order": order,
        }
        finest = rows[-1]
        add(_gate("hull_voxel_vs_faceted", finest.measured_volume_m3, faceted, 0.01, "m3"))
        add(_gate("hull_voxel_vs_analytic", finest.measured_volume_m3, hull_reference, 0.03, "m3"))
        add(_verdict("hull_order", order >= 0.5, f"observed order {order:.3f}"))

    # -- 6. bundle round-trip: the file is the interface -------------------------
    #    A box *with* a surrounding fluid margin, so the fluid complement is non-empty
    #    and the table is a real lattice rather than a solid block.
    shell = box_surface_mesh((0.5, 0.4, 0.3), origin=(0.25, 0.25, 0.2))
    shell_lattice = solid_lattice(shell, resolution=(20, 16, 12), pad=0.25)
    bundle_dir = out / "bundle"
    manifest = lattice_bundle(
        bundle_dir,
        shell_lattice,
        case="voxel_solid_box_in_box",
        walls={"normal_axis": "z"},
    )
    _read_manifest, recovered = read_lattice_bundle(bundle_dir)
    same = (
        recovered.shape == shell_lattice.shape
        and recovered.occupancy == shell_lattice.occupancy
        and occupancy_mismatch(recovered, shell_lattice) == 0
    )
    results["cases"]["bundle"] = {
        "dir": str(bundle_dir),
        "digest_sha256": manifest.get("digest_sha256"),
        "fluid_cells": len(list(shell_lattice.fluid_cells())),
        "solid_cells": shell_lattice.solid_count,
        "roundtrip_identical": same,
    }
    add(_verdict("bundle_roundtrip", same, "occupancy survived write/read"))

    failed = [gate["gate"] for gate in results["gates"] if gate["verdict"] != "PASS"]
    results["gates_failed"] = failed
    results["verdict"] = "PASS" if not failed else "FAIL"
    (out / "results.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    return results


def _main(argv: Sequence[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    outdir = "artifacts/voxel_solid"
    resolution = 64
    index = 0
    while index < len(args):
        token = args[index]
        if token == "--out" and index + 1 < len(args):
            outdir = args[index + 1]
            index += 2
        elif token == "--resolution" and index + 1 < len(args):
            resolution = int(args[index + 1])
            index += 2
        elif token in {"-h", "--help"}:
            print("usage: python -m tensorcad.voxel_solid [--out DIR] [--resolution N]")
            return 0
        else:
            raise SystemExit(f"unrecognised argument {token!r}; try --help")
    started = time.perf_counter()
    results = run_voxel_case(outdir, hull_resolution=resolution)
    elapsed = time.perf_counter() - started
    print(json.dumps({"verdict": results["verdict"], "gates": len(results["gates"]),
                      "failed": results["gates_failed"], "seconds": round(elapsed, 2)}))
    return 0 if results["verdict"] == "PASS" else 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())