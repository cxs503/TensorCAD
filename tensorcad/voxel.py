"""SDF voxel-grid export: the solver-facing form of a :class:`GeometryField`.

``to_mesh`` answers "give me a surface". A voxel/LBM/FVM solver, an L2 geometry
encoder and the scan_runner SDF stack instead want the *volume* sampling: one
signed distance per cell, on a regular grid, with the cell centres, the spacing
and the sign convention written down. This module is that export, for every field
in the protocol (four analytic primitives, the SUBOFF hull and appendages, and a
learned ``.ngo`` field), so a solver never has to know which representation it is
reading.

Two decisions are deliberate and are the ones a downstream solver must know:

* **Cell-centred, not vertex-centred.** ``sdf[i, j, k]`` is the signed distance at
  ``origin + (i + 1/2, j + 1/2, k + 1/2) * spacing`` -- the value of the *cell*, which
  is what a finite-volume cell wants. A vertex-centred grid would have ``shape + 1``
  samples and half a cell of position ambiguity at the boundary.
* **Negative inside.** The protocol's convention (``inside_sign="negative"``) travels
  in the handoff record rather than being re-derived by the consumer.

Volume comes in two flavours, because they converge at different orders and the
difference is itself evidence (:func:`VoxelGrid.volume` vs
:func:`VoxelGrid.occupied_volume`).

NumPy is imported inside the functions, never at module import: ``import tensorcad``
keeps working in a standard-library-only environment, and only the caller who asks
for a grid pays for NumPy (same rule as ``to_mesh``).
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass, field as dataclass_field
from typing import Any, Sequence

from .geometry_field import FieldError, GeometryField

#: Schema tag for :meth:`VoxelGrid.as_handoff`.
VOXEL_HANDOFF_SCHEMA = "tensorcad.geometry.voxel_handoff/1"

#: Sign convention of ``sdf``; identical to :data:`GeometryField.inside_sign`.
SIGN_CONVENTION = "negative_inside"

#: Cells per sampling block. 2**18 float64 triples is ~6 MB of query points, which
#: keeps a 256^3 export inside a few hundred MB without a caller-set parameter.
DEFAULT_CHUNK_CELLS = 1 << 18


def _numpy() -> Any:
    """Import NumPy on demand, with an error that names the fix."""
    try:
        import numpy as np
    except ImportError as exc:  # pragma: no cover - exercised only without NumPy
        raise FieldError(
            "voxel export needs NumPy (as does to_mesh); the analytic fields themselves "
            "do not. Install it with `pip install numpy`."
        ) from exc
    return np


def _cells_per_axis(resolution: int | Sequence[int]) -> tuple[int, int, int]:
    """Normalise a scalar or per-axis resolution into three cell counts."""
    if isinstance(resolution, (str, bytes)):
        raise FieldError(f"resolution must be a cell count, got {resolution!r}")
    if hasattr(resolution, "__len__"):
        values = tuple(int(value) for value in resolution)  # type: ignore[union-attr]
        if len(values) != 3:
            raise FieldError(f"resolution needs 1 or 3 cell counts, got {len(values)}")
        cells = values
    else:
        cells = (int(resolution),) * 3
    if min(cells) < 2:
        raise FieldError(f"resolution must be at least 2 cells per axis, got {cells}")
    return cells  # type: ignore[return-value]


def _resolve_bounds(
    field: GeometryField, bounds: Sequence[Sequence[float]] | None
) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    """Return the sampling box: the field's own bounds, or a caller override."""
    if bounds is None:
        box = field.bounds()
        if box is None:
            raise FieldError(
                f"{field.kind} is unbounded and has no grid to sample; pass bounds=((x0, y0, z0), (x1, y1, z1))"
            )
        return (
            tuple(float(value) for value in box.minimum),  # type: ignore[return-value]
            tuple(float(value) for value in box.maximum),  # type: ignore[return-value]
        )
    corners = list(bounds)
    if len(corners) != 2:
        raise FieldError(f"bounds needs two corners (min, max), got {len(corners)}")
    resolved = []
    for name, corner in zip(("minimum", "maximum"), corners, strict=True):
        values = list(corner)
        if len(values) != 3:
            raise FieldError(f"bounds {name} needs 3 coordinates, got {len(values)}")
        try:
            resolved.append(tuple(float(value) for value in values))
        except (TypeError, ValueError) as exc:
            raise FieldError(f"bounds {name} is not numeric: {values!r}") from exc
    minimum, maximum = resolved
    if not all(maximum[axis] > minimum[axis] for axis in range(3)):
        raise FieldError(f"bounds are degenerate or inverted: {minimum!r} -> {maximum!r}")
    if not all(abs(value) < float("inf") for value in (*minimum, *maximum)):
        raise FieldError(f"bounds are not finite: {minimum!r} -> {maximum!r}")
    return minimum, maximum  # type: ignore[return-value]


@dataclass(frozen=True, eq=False)
class VoxelGrid:
    """A signed-distance grid plus the metadata a solver needs to read it.

    ``eq=False`` on purpose: two grids are the same when their bytes are the same, and
    that question is answered by :meth:`sdf_digest`, not by an elementwise ``==``
    whose result is an array.
    """

    sdf: Any
    origin: tuple[float, float, float]
    spacing: tuple[float, float, float]
    provenance: dict[str, Any] = dataclass_field(default_factory=dict)

    # -- geometry ---------------------------------------------------------
    @property
    def shape(self) -> tuple[int, int, int]:
        return tuple(int(value) for value in self.sdf.shape)  # type: ignore[return-value]

    @property
    def cells(self) -> int:
        nx, ny, nz = self.shape
        return nx * ny * nz

    @property
    def cell_volume(self) -> float:
        hx, hy, hz = self.spacing
        return hx * hy * hz

    @property
    def characteristic_cell_size(self) -> float:
        """Isotropic-equivalent cell size ``(hx·hy·hz)^(1/3)``.

        The partial-fill rule below is a plane-cut fraction, which is exact for a face
        cut and curvature-biased otherwise; on an anisotropic grid the bias scales with
        this length, so the estimator uses it rather than, say, ``min(h)``.
        """
        return self.cell_volume ** (1.0 / 3.0)

    @property
    def bounds_min(self) -> tuple[float, float, float]:
        return self.origin

    @property
    def bounds_max(self) -> tuple[float, float, float]:
        return tuple(  # type: ignore[return-value]
            self.origin[axis] + self.shape[axis] * self.spacing[axis] for axis in range(3)
        )

    def axis(self, index: int) -> Any:
        """Cell-centre coordinates along ``index`` (0, 1, 2)."""
        if index not in (0, 1, 2):
            raise FieldError(f"axis must be 0, 1 or 2, got {index!r}")
        np = _numpy()
        return self.origin[index] + (np.arange(self.shape[index]) + 0.5) * self.spacing[index]

    # -- measures ---------------------------------------------------------
    def volume_fraction(self) -> Any:
        """Per-cell solid fraction, linear (plane-cut) estimate.

        ``f = clip(0.5 - s / h, 0, 1)``: a cell whose centre sits on the surface is half
        full, one ``h/2`` inside is full, one ``h/2`` outside is empty. Exact for a plane
        cutting a cell parallel to a face, second-order for a curved surface -- which is
        why :meth:`volume` beats :meth:`occupied_volume` (see docs/voxel_export.md).
        """
        np = _numpy()
        return np.clip(0.5 - self.sdf / self.characteristic_cell_size, 0.0, 1.0)

    def volume(self) -> float:
        """Volume by the partial-fill rule (second-order in the cell size)."""
        np = _numpy()
        return float(np.sum(self.volume_fraction())) * self.cell_volume

    def occupied_volume(self) -> float:
        """Volume by hard cell counting (first-order in the cell size).

        Kept because it is what a boolean occupancy grid actually represents: an LBM/FVM
        solid-fraction array set to 0/1. The gap between the two volumes measures the
        discretisation error a solver inherits, so both are reported.
        """
        np = _numpy()
        return float(np.count_nonzero(self.sdf < 0.0)) * self.cell_volume

    def inside_count(self) -> int:
        np = _numpy()
        return int(np.count_nonzero(self.sdf < 0.0))

    # -- handoff ----------------------------------------------------------
    def sdf_digest(self) -> str:
        """SHA-256 over the grid's float64 bits, so a solver can name what it read."""
        np = _numpy()
        payload = np.ascontiguousarray(self.sdf, dtype=np.float64).tobytes()
        return hashlib.sha256(payload).hexdigest()

    def as_handoff(self) -> dict[str, Any]:
        """JSON-ready record for a downstream solver (TensorLBM/TensorFVM/TensorFEM)."""
        return {
            "schema": VOXEL_HANDOFF_SCHEMA,
            "quantity": "signed_distance",
            "inside_sign": SIGN_CONVENTION,
            "sampling": "cell_centred",
            "units": self.provenance.get("units", "m"),
            "coordinate_frame": self.provenance.get("coordinate_frame", "object_local"),
            "shape": list(self.shape),
            "cells": self.cells,
            "spacing_m": [float(value) for value in self.spacing],
            "origin_m": [float(value) for value in self.origin],
            "bounds_min_m": [float(value) for value in self.bounds_min],
            "bounds_max_m": [float(value) for value in self.bounds_max],
            "cell_volume_m3": self.cell_volume,
            "volume_m3": self.volume(),
            "occupied_volume_m3": self.occupied_volume(),
            "sdf_dtype": str(self.sdf.dtype),
            "sdf_sha256": self.sdf_digest(),
            "provenance": dict(self.provenance),
        }


def to_voxel_grid(
    field: GeometryField,
    resolution: int | Sequence[int] = 48,
    *,
    padding: int = 0,
    bounds: Sequence[Sequence[float]] | None = None,
    dtype: str = "float64",
    chunk_cells: int = DEFAULT_CHUNK_CELLS,
) -> VoxelGrid:
    """Sample ``field`` onto a regular cell-centred signed-distance grid.

    ``resolution`` is cells per axis, one number or three (a slender body wants the
    anisotropic form, exactly as for ``to_mesh``). ``padding`` adds that many cells of
    margin on every side *without* changing the spacing -- the grid stays aligned with
    the field's own bounds, so the same field at the same resolution always produces the
    same cells. ``bounds`` overrides the sampling box (a far field for a solver, or a
    window onto part of a body).

    Sampling goes through ``field.numpy_callable()`` when the field offers one and falls
    back to the point-list contract otherwise, in blocks of ``chunk_cells`` so a 256^3
    export does not allocate one giant point array.
    """
    np = _numpy()
    cells = _cells_per_axis(resolution)
    minimum, maximum = _resolve_bounds(field, bounds)

    pad = int(padding)
    if pad < 0:
        raise FieldError(f"padding must be a non-negative number of cells, got {padding!r}")
    if int(chunk_cells) < 1:
        raise FieldError(f"chunk_cells must be positive, got {chunk_cells!r}")
    try:
        grid_dtype = np.dtype(dtype)
    except TypeError as exc:
        raise FieldError(f"dtype {dtype!r} is not a NumPy dtype") from exc
    if grid_dtype.kind != "f":
        raise FieldError(f"dtype must be a floating type, got {grid_dtype}")

    sizes = tuple(maximum[axis] - minimum[axis] for axis in range(3))
    spacing = tuple(sizes[axis] / cells[axis] for axis in range(3))
    counts = tuple(cells[axis] + 2 * pad for axis in range(3))
    origin = tuple(minimum[axis] - pad * spacing[axis] for axis in range(3))

    axes = [
        origin[axis] + (np.arange(counts[axis], dtype=np.float64) + 0.5) * spacing[axis]
        for axis in range(3)
    ]

    sdf = np.empty(counts, dtype=grid_dtype)
    evaluate = field.numpy_callable()

    def sample(points: Any) -> Any:
        if evaluate is not None:
            return np.asarray(evaluate(points), dtype=np.float64)
        return np.asarray(field.signed_distance(np.asarray(points).tolist()), dtype=np.float64)

    ny, nz = counts[1], counts[2]
    plane = ny * nz
    block = max(1, min(counts[0], int(chunk_cells) // plane))
    yy, zz = np.meshgrid(axes[1], axes[2], indexing="ij")
    plane_y = yy.ravel()
    plane_z = zz.ravel()
    for start in range(0, counts[0], block):
        stop = min(start + block, counts[0])
        rows = stop - start
        points = np.empty((rows * plane, 3), dtype=np.float64)
        points[:, 0] = np.repeat(axes[0][start:stop], plane)
        points[:, 1] = np.tile(plane_y, rows)
        points[:, 2] = np.tile(plane_z, rows)
        values = sample(points)
        if values.size != rows * plane:
            raise FieldError(
                f"{field.kind} returned {values.size} distances for {rows * plane} points"
            )
        sdf[start:stop] = values.reshape(rows, ny, nz)

    provenance = {
        "producer": "TensorCAD",
        "representation": "voxel_grid",
        "field_kind": field.kind,
        "field_parameters": dict(field.parameters()),
        "field_provenance": dict(field.provenance()),
        "units": field.units,
        "coordinate_frame": field.coordinate_frame,
        "field_bounds_min_m": [float(value) for value in minimum],
        "field_bounds_max_m": [float(value) for value in maximum],
        "resolution_cells": list(cells),
        "padding_cells": pad,
        "sampled_on": "cell_centres",
    }
    return VoxelGrid(sdf=sdf, origin=origin, spacing=spacing, provenance=provenance)


__all__ = [
    "DEFAULT_CHUNK_CELLS",
    "SIGN_CONVENTION",
    "VOXEL_HANDOFF_SCHEMA",
    "VoxelGrid",
    "to_voxel_grid",
]