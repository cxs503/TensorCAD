# SDF voxel-grid export (0.10)

`to_mesh` answers "give me a surface"; a finite-volume, LBM or L2-encoder consumer wants
the *volume*. `to_voxel_grid` samples any `GeometryField` — the four analytic primitives,
the SUBOFF hull and its appendages, and a learned `.ngo` field — onto a regular
cell-centred signed-distance grid, and returns a `VoxelGrid` that carries its own
provenance and a JSON handoff record. It closes the gap logged in
[`proposals/tensorcad-sdf-voxel-export.txt`](../proposals/tensorcad-sdf-voxel-export.txt)
and in the audit backlog ("no SDF voxel grid export: field → mesh only").

```python
from tensorcad import SuboffHullField, to_voxel_grid

grid = to_voxel_grid(SuboffHullField(), (218, 26, 26))
grid.volume()                      # 0.698982 m^3  (published 0.699211, -0.033 %)
grid.occupied_volume()             # 0.707470 m^3  (boolean occupancy, +1.18 %)
grid.as_handoff()["sdf_sha256"]    # names the exact grid that was read
grid.sdf                           # (218, 26, 26) float64, negative inside
```

## Conventions the consumer has to know

| Item | Value |
|---|---|
| Sample position | cell centres: `sdf[i, j, k]` is the distance at `origin + (i+½, j+½, k+½)·spacing` |
| Sign | negative inside (the protocol's `inside_sign`, carried in the handoff record) |
| `spacing` | `box_size / cells` per axis, so `shape = cells + 2·padding` at unchanged spacing |
| `origin` | `bounds_min − padding·spacing`; `bounds[0]` of the field unless `bounds=` overrides |
| Default box | `field.bounds()` — no padding, so the grid is reproducible from the field alone |
| dtype | `float64` default; `float32` halves the bytes and moves the volume by < 1e-5 relative |
| Units / frame | taken from the field (`m`, `object_local`) |
| Missing capability | there is no `surface_area()` — see "What this does not provide" |

Padding *grows the box without moving the grid*: interior samples are bit-identical to the
unpadded grid, so a solver can add a far field without invalidating a cached digest.

## Accuracy: volume against closed forms

Two estimators are offered because they converge at different orders and the difference is
what a solver inherits if it only has a boolean occupancy grid:

* `volume()` — partial-fill rule `f = clip(½ − s/h, 0, 1)` per cell, summed.
* `occupied_volume()` — hard cell counting (`sdf < 0`), i.e. a 0/1 solid-fraction array.

Sphere, r = 0.5 m (exact 0.523598776 m³), box ±0.6 m so the surface is not grid-coincident:

| cells | cells | `volume()` [m³] | err % | `occupied_volume()` [m³] | err % |
|---|---|---|---|---|---|
| 8 | 512 | 0.526265575 | +0.5093 | 0.540000000 | +3.1324 |
| 16 | 4,096 | 0.526280060 | +0.5121 | 0.509625000 | -2.6688 |
| 32 | 32,768 | 0.524775777 | +0.2248 | 0.524812500 | +0.2318 |
| 64 | 262,144 | 0.523791244 | +0.0368 | 0.524601563 | +0.1915 |

Least-squares order of the two error columns: partial-fill **1.26**, hard count **1.56**
(the hard-count column changes sign between resolutions — it is aliasing, not a clean order,
which is the argument for the partial-fill rule). At 64³ the partial-fill volume is
+0.0368 % and hard counting is +0.1915 %.

## Anisotropic cells, at equal cost

Ellipsoid, semi-axes (0.5, 0.25, 0.125) m, exact 0.065449847 m³:

| cells (nx, ny, nz) | spacing [m] | `volume()` [m³] | err % |
|---|---|---|---|
| (16, 16, 16) | (0.06250, 0.03125, 0.01562) | 0.066369974 | +1.4059 |
| (32, 16, 8) | (0.03125, 0.03125, 0.03125) | 0.065884237 | +0.6637 |
| (64, 32, 16) | (0.01562, 0.01562, 0.01562) | 0.065515047 | +0.0996 |

Same cell count, same body: an isotropic 32³ grid gives **+0.2426 %**,
while (64, 32, 16) — the same 32,768 cells, distributed along the axes — gives
**+0.0996 %**. That is the whole reason `to_voxel_grid` takes three cell counts, as
`to_mesh` already does.

## The published hull

SUBOFF bare hull (AFF-1), published displacement volume 0.699211 m³ (Groves, Huang & Chang
1989, via the offset table). This is the same external reference the mesh path answers to
([docs/suboff.md](suboff.md)), applied here to the voxel path:

| cells (nx, ny, nz) | cell size [m] | `volume()` [m³] | err % | hard-count err % | occupancy | time [s] |
|---|---|---|---|---|---|---|
| (109, 13, 13) | 0.03937 | 0.700481556 | +0.1817 | +2.5078 | 63.76 % | 0.22 |
| (218, 26, 26) | 0.01969 | 0.698982292 | -0.0327 | +1.1812 | 62.93 % | 2.69 |
| (436, 52, 52) | 0.00984 | 0.699410899 | +0.0286 | +0.1704 | 62.30 % | 21.09 |

The partial-fill error falls from +0.1817 % to +0.0286 %;
the boolean-occupancy error starts at +2.508 % because a cell whose
centre is inside counts as fully solid. Sampling runs at 84 k cells/s, 55 k cells/s, 56 k cells/s on the vectorised path.

The grid is *not* an approximation of the hull in the way a mesh is: every sample is the
field's own analytic distance (zero sampling error); the error above is entirely the
integral of that field over cells.

## Handoff record

`VoxelGrid.as_handoff()` is a JSON document (`tensorcad.geometry.voxel_handoff/1`), so a
consumer can record exactly which grid it read:

```json
{
  "schema": "tensorcad.geometry.voxel_handoff/1",
  "quantity": "signed_distance",
  "inside_sign": "negative_inside",
  "sampling": "cell_centred",
  "units": "m",
  "coordinate_frame": "object_local",
  "shape": [
    218,
    26,
    26
  ],
  "cells": 147368,
  "spacing_m": [
    0.019982110557798166,
    0.01953846075692308,
    0.01953846075692308
  ],
  "origin_m": [
    0.0,
    -0.25399998984000005,
    -0.25399998984000005
  ],
  "cell_volume_m3": 7.628199654518588e-06,
  "volume_m3": 0.6989823039274436,
  "occupied_volume_m3": 0.707469748758672,
  "sdf_dtype": "float32",
  "sdf_sha256": "83e4d51d483184be396a07dbf73e669b0bfdc814a41cb8110de40e1f28f94f21"
}
```

`sdf_sha256` is taken over the contiguous float64 bits of the grid, so it names the data, not
the file. `grid.sdf.tobytes()` (or `numpy.save`) is the transport.

## What this does not provide

| Not provided | Why / what to use instead |
|---|---|
| `surface_area()` | Measured and **removed**: lattice delta quadrature gave −0.159 % on the sphere at 64³ but **−8.4 %** on the SUBOFF hull at (218, 26, 26), because the deeply diagonal stern taper is where a lattice delta sum loses area. `to_mesh(field, cells).report.surface_area` is −0.066 % there, so the mesh path keeps that job. A regression test pins the absence. |
| Boolean/CSG, adaptive refinement, narrow-band storage | Out of scope for this export; the field is what you edit, the grid is what you hand over. |
| Marching cubes / isosurface extraction | `to_mesh` (marching tetrahedra) already exists and reports watertightness. |
| Compressed/streamed formats (`.npy` is your call, no HDF5) | Keeping the module standard-library-plus-NumPy; a consumer writes its own container. |
| Time-to-first-cell for a 512³ grid | Not measured here; memory is `cells × 8 B` (`float64`) or `× 4 B` (`float32`). |

Sharp-tip bodies (SUBOFF bow `ρ ~ ξ^0.476`, elliptical stern pole) are sampled, not
smoothed: no cell centre is special-cased, so the volume error at the tips is the same
partial-fill rule that applies everywhere else.

## Zero-dependency rule

NumPy is imported inside the functions, never at module import. `import tensorcad` and
`import tensorcad.voxel` still work with NumPy blocked; only the call to `to_voxel_grid`
raises (a `FieldError` that names the fix). A test enforces both halves of that.

## Reproduction

```bash
python examples/geometry_core/run_voxel_export.py --output artifacts/voxel_export
python -m pytest -q tests/test_voxel.py
```

The tables above are generated from `artifacts/voxel_export/voxel_export.json`, which the
example script rewrites on every run; the script exits non-zero if a claim stops holding.

## Provenance

| Item | Value |
|---|---|
| Implementation | `tensorcad/voxel.py` (0.10.0), re-exported as `tensorcad.to_voxel_grid` / `tensorcad.VoxelGrid` |
| Handoff schema | `tensorcad.geometry.voxel_handoff/1` |
| Proposal | `proposals/tensorcad-sdf-voxel-export.txt` |
| Reference volume | Groves, Huang & Chang (1989), DTRC/SHD-1298-01 (offset table), as transcribed in `tensorcad/suboff.py` |
| Tests | `tests/test_voxel.py` (20 tests: convention, exactness, convergence, sign, errors, no-NumPy import, example run) |
