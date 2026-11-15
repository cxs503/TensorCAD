# CAD → voxel face: occupancy from a closed surface, and what it costs

Anchor: `GAP-CAD-B3-SOLVER-INTERFACE` (CAD→求解器消费面). This closes the item the
previous round wrote down about itself: *"the hull ray-cast voxelisation needs a
spatial index (the naive O(cells × faces) version did not finish)"*. The third
solver face (lattice / voxel) now has a CAD producer, and the producer is
verified against a reference, not against itself.

## 1. What was missing

`solver_interface.write_bundle_lattice` published the lattice table, and
`solver_consumers.lbm_channel_from_bundle` consumed one — but the only lattices
ever written had their fluid cells **declared by the caller** (a box-aligned
channel). Nothing turned a CAD surface into an occupancy grid, so the LBM route
could not see a hull. Two things were absent, not one:

1. **the fill**: per-cell parity is O(cells × faces) — 2.1×10⁶ cells × 1998
   triangles ≈ 4×10⁹ triangle tests for a 128³ lattice, which is the round that
   "did not finish";
2. **the reference**: an occupancy grid that is only compared with itself proves
   nothing. The hull module already ships `analytic_volume`, a closed form the
   triangulation never saw.

## 2. Method

`tensorcad/voxel_solid.py`. One ray per **column** through
`spatial_index.TriangleGrid` (uniform grid + 3-D DDA), crossings sorted once per
column, interior emitted as **spans** by a two-pointer walk on the cell centre
parameters. Cost per column is O(cells crossed + hits) instead of O(faces), and
the fill is a bulk slice assignment.

Cell ordering is the contract's: `occupancy[i + nx*(j + ny*k)]`, `i` along `+x`
fastest, `(i, j, k)` ⇒ `+x, +y, +z` — the same ordering
`write_bundle_lattice` writes, in the `z-up` frame of `COORDINATE_CONTRACT`.
Precondition: the surface must be closed. `watertight_gate` runs
`meshquality.surface_topology` first and its verdict is the first gate.

## 3. A defect found on the way (and fixed)

The first probe table was copied from `spatial_index._CONTAIN_PROBE`
(`(1.0, 0.0013, 0.0007)`), which is sized for a **point** query. Over a 100 m
column the 1.3×10⁻³ tilt drifts 0.13 m in `y` and 0.077 m in `z` — for a hull
lattice with `h_z = 0.08 m` the "column" walks out of its own column. Measured
consequences: a **−3.5 % volume bias at every resolution** (32/48/64/96 cells:
−2.81 %/−2.00 %/−2.85 %/−2.81 %, i.e. no convergence at all) and **546 of 2880
cells disagreeing with the brute-force reference**. One ellipsoid test showed the
same plateau. Fixed by tilting only 1×10⁻⁶ (drift 10⁻⁴ m over 100 m) and by
dividing instead of multiplying by the direction's axial component when mapping a
cell centre to a ray parameter. Both facts are pinned by tests
(`test_ray_direction_is_unit_and_axial` asserts the tilt, the cross-check test
asserts cell-for-cell agreement).

## 4. Verification (driver: `python -m tensorcad.voxel_solid --out artifacts/voxel_solid`)

Geometry: `hull.build_mesh(HullParameters(nx=41, nz=13))` — 1001 nodes, 1998
triangles, `L/B/T = 10.0/2.0/0.8 m`, waterline `z = 0`, keel `z = −0.8 m`.

| gate | measured | reference | error | tolerance | verdict |
|---|---|---|---|---|---|
| box fills exactly (25×15×10) | 0.24000000000000005 m³ | 0.24 m³ | 2×10⁻¹⁶ | 1e-12 | PASS |
| ray axis y vs x / z vs x | 0 / 0 cells differ | 0 | — | 0 | PASS |
| hull watertight | boundary_edges 0, non-manifold 0 | closed | — | — | PASS |
| grid vs brute force (12³ = 1728 cells, 1998 tris) | 0 differ | 0 | — | ≤1e-3·cells | PASS |
| hull voxel vs faceted displacement (128³) | 7.447571 m³ | 7.446554 m³ | +1.37e-4 | 1 % | PASS |
| hull voxel vs **analytic** volume (128³) | 7.447571 m³ | 7.512698 m³ | −8.67e-3 | 3 % | PASS |
| convergence order, log-log fit | 1.959 | ≥ 0.5 | — | — | PASS |
| lattice bundle round-trip | occupancy identical | identical | 0 | 0 | PASS |

`gates_failed = []` — 8 gates, 16 s of wall clock for the whole driver.

### Convergence, hull against its own faceted displacement

| resolution | h [m] | voxel volume [m³] | rel. error | seconds |
|---:|---:|---:|---:|---:|
| 32³ | 0.3125 | 7.461914 | +2.06e-3 | 0.49 |
| 64³ | 0.1563 | 7.449585 | +4.07e-4 | 1.81 |
| 128³ | 0.0781 | 7.447571 | +1.37e-4 | 7.57 |

Sphere (`ellipsoid_mesh((0.5,0.5,0.5))`, reference 4/3·π·r³): 5.42 %/0.84 %/0.40 %/
0.38 % at 10/20/40/56 cells, fitted order **1.88**. The worst case for an
indicator-function quadrature is first order; a smooth convex body measures about
1.9 because the curved-boundary error largely cancels. The hull measures 1.96 for
the same reason: its taper is smooth, and the residual is the flatter
`sqrt(z/draft)` region near the keel. Figure:
`docs/voxel_solid_convergence.png` (mid-ship section, deck-plan slice, log-log
convergence with both references).

### Cost

| lattice | cells | fill [s] | per-cell parity would need |
|---|---:|---:|---:|
| 1728 (12³) | 1.7e3 | 0.09 | 3.5e6 triangle tests |
| 262144 (64³) | 2.6e5 | 1.81 | 5.2e8 |
| 2097152 (128³) | 2.1e6 | 7.57 | 4.2e9 |

The right column is the arithmetic of the path that did not finish; the middle is
the same problem solved by casting 16 k (64³) or 65 k (128³) column rays instead
of 2.1 M point rays.

### Handoff

`lattice_bundle` writes `lattice.json` (fluid cells = complement of the solid)
through `write_bundle_lattice`, so a voxel/LBM consumer reads the same schema it
already reads; `verify_bundle` re-reads it with its own parser; the solid mask
itself travels packed in the manifest `extra`
(`solid_occupancy_b64` + the `voxel_solid` record with shape/origin/spacing/
digest/frame), because a fluid-only table would leave the hull invisible.
`read_lattice_bundle` recovers the lattice and refuses a bundle whose two tables
disagree. Artifact: `artifacts/voxel_solid/bundle/`.

## 5. Points of caution

1. **Parity needs a closed surface.** `watertight_gate` is a gate, not advice: on
   an open surface the fill produces stripes (`test_open_surface_is_reported_not_voxelised_silently`
   measures the difference). Cells whose centre lies *exactly* on the surface are
   undefined, as for any floating-point parity test; the driver perturbs every
   disagreeing centre by ±1e-9 m and reports how many flip (`mismatched_cells_on_surface`).
2. **Occupancy is a staircase.** The volume is exact for aligned solids and
   first-order (in the worst case) otherwise; there is no sub-cell fraction. A
   solver consuming this lattice inherits a staircase boundary of one cell.
3. **The analytic gap is the mesh, not the lattice.** Voxel vs analytic is
   −0.87 %, of which the *faceted* surface itself is −0.88 % (nx=41, nz=13
   sampling). The lattice reproduces its input surface to 1.4e-4; a product
   decision about hull faceting is a separate question.
4. **The driver's brute-force cross-check is at 12³/1728 cells** because the slow
   path is O(cells × faces) by design. Agreement is cell-for-cell there; the fast
   path is the same code at every resolution, but this is not a proof at 128³.
5. **This is the geometry face of the LBM route, not an LBM result.** The 2-D
   channel case in `docs/solver_consumers.md` is still the LBM consumption
   evidence; nothing here claims a flow, a resistance or a drag. `tensorlbm.stream`
   still transports against its own `d2q9.C` (documented in that doc and in
   `solver_consumers.stream_probe`), so a 3-D hull LBM case would need its own
   streaming step; that is the next layer, not this round.
6. **`tensordg` remains unavailable** in this environment, so the DG curved-wall
   face is still open (as the previous round stated).

## 6. Error verdict against the acceptance rule (≤15 % unacceptable; is it <3 %?)

| key quantity | reference | error | <3 %? |
|---|---|---:|---:|
| hull voxel volume (128³) | analytic displacement `analytic_volume()` | 8.7e-03 | yes |
| hull voxel volume (64³) | same | 8.4e-03 | yes |
| hull voxel volume (128³) | faceted displacement `compute_hydrostatics()` | 1.4e-04 | yes |
| box occupancy | lx·ly·lz | 2e-16 | yes |
| fast fill vs brute force | cell-for-cell | 0 cells | yes |

What these numbers do **not** say: they verify a *discretisation* (surface →
cells) against closed forms on the same geometry. They say nothing about whether a
voxel mesh is the right mesh for a hull flow — that is a solver-side question and
belongs to the LBM/FVM routes, not to this handoff.