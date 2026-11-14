# CAD → solver consumption: the files are the interface

Anchor: `GAP-CAD-B3-SOLVER-INTERFACE` (CAD→求解器消费面). This document covers the
**consumption** half of the handoff contract: the geometry library writes a bundle,
:func:`tensorcad.solver_interface.verify_bundle` re-reads it, and
:mod:`tensorcad.solver_consumers` then parses the table **from disk** and runs a real
solver on it.

## 1. Why this half was missing

`tensorcad/solver_interface.py` (merged earlier in this repository) could write a
z-up bundle and verify it, but nothing consumed one: no code read `surface.json`
back and handed it to a panel code. The earlier attempt at this anchor that lived
only in a worktree passed **in-memory** arrays to the solvers while its report said
"the files are the interface", and its `verify_bundle` never ran. Both halves are
now in one place and the claim is machine-checked.

## 2. Routes (three tables, one frame)

| route | table | consumer | case | closed-form reference |
|---|---|---|---|---|
| `surface` | `surface.json` (triangles) | TensorBEM `TriangleMesh` + `compute_added_mass_matrix_3d` | submerged sphere, added mass | `2/3 π ρ a³` (Lamb) |
| `grid` | `grid.json` (structured quads) | TensorFEM `plate.q4_mindlin_stiffness` + assembly | simply-supported plate, sine load | `q0 a⁴/(4 π⁴ D)` (Navier) |
| `lattice` | `lattice.json` (cell occupancy) | TensorLBM `equilibrium`/`collide_bgk`/`macroscopic` | body-force channel | `F H²/(8 ν)` |

Frame, for every route: **z-up**, right-handed, metres, `x` aft→fore, `y` starboard,
`z` up, waterline at `z = 0` (`COORDINATE_CONTRACT` in `solver_interface.py`). The
manifest declares it, `verify_bundle` re-reads it, and each consumer refuses a
bundle whose route/frame is not the one it needs.

New since this round: two writers and two independent readers for the routes that
had none — `write_bundle_grid`/`read_grid_table` (plate and shell solvers consume a
2-D mesh, not triangles) and `write_bundle_lattice`/`read_lattice_table` (LBM and
voxel solvers consume an occupancy mask, not a surface). Both are wired into
`verify_bundle`, so a moved node or a dropped cell fails the digest before any
solver runs.

## 3. Refusals (the part that makes it a contract)

* wrong frame (`y-up` manifest) → `HandoffError`
* tampered table (a node moved by `1e-6`) → digest mismatch → `HandoffError`
* surface with a hole → `bad_edges = 3` → refused
* inward-facing surface (negative signed volume) → refused at write time
* quads not wound about `+z` → negative projected area → refused
* occupancy cell outside the lattice box → refused
* a solver library that is not importable → `status="unavailable"` with the
  exception text (never a zero, never a pass)

## 4. Results (regenerated: `python -m tensorcad.solver_consumers --out artifacts/solver_consumers`)

### sphere_added_mass_16x32 — status `ran` (solver `tensorbem`, repo 9abf785c19ec3cb922aa3c9d22d2950bc8416a69)

| gate | measured | reference | error | tolerance | verdict |
|---|---|---|---|---|---|
| enclosed_volume | 0.5152427176 m3 | 0.5235987756 m3 | -1.5959 % | 2 % | PASS |
| added_mass_heave | 258.4479396 kg | 261.7993878 kg | -1.2802 % | 5 % | PASS |
| added_mass_surge | 259.8540989 kg | 261.7993878 kg | -0.7430 % | 5 % | PASS |

### sphere_added_mass_24x48 — status `ran` (solver `tensorbem`, repo 9abf785c19ec3cb922aa3c9d22d2950bc8416a69)

| gate | measured | reference | error | tolerance | verdict |
|---|---|---|---|---|---|
| enclosed_volume | 0.5198714295 m3 | 0.5235987756 m3 | -0.7119 % | 2 % | PASS |
| added_mass_heave | 260.2959742 kg | 261.7993878 kg | -0.5743 % | 5 % | PASS |
| added_mass_surge | 260.9331 kg | 261.7993878 kg | -0.3309 % | 5 % | PASS |

### plate_sine_8x8 — status `ran` (solver `tensorfem.plate`, repo 3dbc16d6336d18c587849f775afe88e5ae061f1e)

| gate | measured | reference | error | tolerance | verdict |
|---|---|---|---|---|---|
| centre_deflection | 0.0005098266717 m | 0.0005124778342 m | -0.5173 % | 5 % | PASS |
| plate_area | 4 m2 | 4 m2 | +0.0000 % | 1e-10 % | PASS |

### plate_sine_16x16 — status `ran` (solver `tensorfem.plate`, repo 3dbc16d6336d18c587849f775afe88e5ae061f1e)

| gate | measured | reference | error | tolerance | verdict |
|---|---|---|---|---|---|
| centre_deflection | 0.0005131808675 m | 0.0005124778342 m | +0.1372 % | 5 % | PASS |
| plate_area | 4 m2 | 4 m2 | +0.0000 % | 1e-10 % | PASS |

### channel_poiseuille_8 — status `ran` (solver `tensorlbm`, repo 1ae21b9def95cb82fcc525e7f55ed978b97cdb55)

| gate | measured | reference | error | tolerance | verdict |
|---|---|---|---|---|---|
| peak_velocity_lattice | 0.001551999824 lattice | 0.0016 lattice | -3.0000 % | 5 % | PASS |
| centreline_curvature | 0.0001999999853 lattice | 0.0002 lattice | -0.0000 % | 1 % | PASS |
| mass_conservation | 8.697043086e-12 lattice | 0 lattice | +0.0000 % | 1e-07 % | PASS |

### channel_poiseuille_16 — status `ran` (solver `tensorlbm`, repo 1ae21b9def95cb82fcc525e7f55ed978b97cdb55)

| gate | measured | reference | error | tolerance | verdict |
|---|---|---|---|---|---|
| peak_velocity_lattice | 0.006351999286 lattice | 0.0064 lattice | -0.7500 % | 5 % | PASS |
| centreline_curvature | 0.000199999993 lattice | 0.0002 lattice | -0.0000 % | 1 % | PASS |
| mass_conservation | 1.439275366e-10 lattice | 0 lattice | +0.0000 % | 1e-07 % | PASS |

## 5. CAD → voxelisation exactness (curved solid)

Ellipsoid semi-axes [0.5, 0.35, 0.25] m, analytic volume 0.183259571 m³:

| resolution | occupied volume [m³] | rel. error |
|---:|---:|---:|
| 16 | 0.185937500 | +1.4613 % |
| 32 | 0.184313965 | +0.5754 % |

First order in the cell size, as a boolean occupancy grid must be; the box-aligned
channel lattice is exact by construction (every cell declared fluid, count gate
0 %).

## 6. Cross-repository finding: `tensorlbm.stream` sign

`stream_probe()` places a unit population on the `+x` direction and streams once:
measured `dx_lattice = -1` where `+1` is required
(`measured`). The release measured here (repo
`1ae21b9def95cb82fcc525e7f55ed978b97cdb55`)
therefore transports against its own `d2q9.C` for a `(9, ny, nx)` table. The channel
case uses TensorLBM's equilibrium, BGK collision and moments and a **local** half-way
bounce-back streaming step; this is a hybrid and is labelled as one. Not fixed here —
TensorLBM owns that code.

## 7. Known limitations (stated, not papered over)

1. **Two-dimensional channel.** The LBM case is the x–z plane with one cell through
   `y`; no claim is made about drag prediction.
2. **Plain Guo source.** With the `(1 − 1/(2τ))` split forcing and the half-force
   velocity the measured steady amplitude came out at `(1 − 1/(2τ))F` — a 62.5 %
   shortfall at `τ = 0.8`, caught by the closed-form gate. The plain source is used
   instead; it leaves a constant `O(F)` offset, reported as
   `forcing_offset_lattice` and visible as the −3.0 % peak-velocity error at
   `nz = 8`. The *shape* is asserted through the discrete curvature (constant
   second difference `−F/ν`), which is independent of that offset.
3. **Host-bound solver checkouts.** The three solver libraries are editable installs
   from other worktrees; each `results.json` case records module path, version and
   the checkout's git commit. On a host without them the cases report
   `unavailable` and the local gates are the only evidence.
4. **Not a ship case.** A sphere, a square plate and a planar channel are standard
   verification geometries with closed-form answers, not hull flows. The hull →
   solver case belongs to the next layer of this anchor.
5. All three cases are hand-written drivers over the solver libraries' element /
   kernel APIs; none of the three libraries ships a bundle reader, so the readers
   here are the interface.
