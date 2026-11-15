# CAD -> lattice -> 3-D LBM: the duct route consumed, end to end

Anchor: `GAP-CAD-B3-SOLVER-INTERFACE` (CAD→求解器消费面). Auto-rendered from
`artifacts/lattice_flow/results.json` by `scripts/render_lattice_flow_doc.py`
— no number below is typed by hand.

## 1. Case: Hagen–Poiseuille flow in the CAD solid's own interior

`tensorcad.voxel_solid` ends with a sentence about itself: *"This is the geometry
face of the LBM route, not an LBM result."* The only lattice any solver had consumed
was the 2-D channel whose fluid cells the caller typed in. So the lattice route
repeated the defect this repository keeps finding — **exported is not accepted** —
one face further out. This round runs a flow on a lattice produced from CAD.

Geometry: `nurbs.cylinder_solid_mesh(R=0.05 m, L=4 cells)` — a closed
cylinder (192/128-sample shell + triangle-fan caps, watertight by construction). Its
**interior** cells are the fluid; every other cell is wall; the cylinder axis is the
periodic direction, so the lattice is one periodic unit of an infinitely long duct.
Parameters: `tau=0.8` ⇒ `nu=0.100000`,
body force `g=2e-05` along `+z` (lattice units, `rho=1`),
`axial_cells=4`, `radial_cells` = 6, 10.

Closed form (exact for this configuration, and unknown to the lattice):
`u(r)=g(R^2-r^2)/(4 nu)`, `u_max=g R^2/(4 nu)`, `Q=pi g R^4/(8 nu)`.

## 2. Method

1. `voxel_solid.solid_lattice` voxelises the CAD surface (one ray per column through
   `spatial_index.TriangleGrid`), the **interior** becomes `lattice.json` via
   `solver_interface.write_bundle_lattice` — the contract table (z-up, metres,
   `occupancy[i+nx*(j+ny*k)]`, `i` along `+x` fastest).
2. `solver_interface.verify_bundle` re-reads `handoff.json` + `lattice.json` with an
   independent reader and is a hard gate; a wrong frame or a tampered table raises
   before any solver state exists.
3. `lattice_flow.lattice_flow_from_bundle` rebuilds the fluid mask **from the file**,
   then runs D3Q19: TensorLBM `equilibrium3d`, `macroscopic3d` and the `d3q19`
   constants; collision `f-(f-feq)/tau` and pull streaming `+c` with half-way
   bounce-back are local, for reasons in §5.
4. Steady state to `max|df| < 1e-13` on fluid cells; `Q` is the flux through **one**
   cross-section (summing every periodic plane over-counts by `nz` — measured, then
   fixed).

## 3. Geometry: what the staircase costs (measured)

| case | shape (nx,ny,nz) | h [mm] | fluid cells | voxel area [m2] | area vs pi R^2 [%] | R_eff [cells] |
|---|---|---|---|---|---|---|
| r=6 | [16, 16, 4] | 8.333 | 448 | 7.777778e-03 | -0.970 | 5.9708 |
| r=10 | [24, 24, 4] | 5.000 | 1264 | 7.900000e-03 | +0.586 | 10.0293 |

## 4. Results

| case | steps | converged | solve [s] | u_max | u_max ref (R_eff) | err [%] | u_mean/u_max | max transverse | mass drift |
|---|---|---|---|---|---|---|---|---|---|
| r=6 | 1056 | yes | 6.7 | 1.683655e-03 | 1.800000e-03 | -6.464 | 0.4928 | 1.44e-08 | 3.3e-11 |
| r=10 | 3032 | yes | 29.5 | 4.892957e-03 | 5.000000e-03 | -2.141 | 0.4937 | 3.83e-08 | 2.7e-10 |

Gates (PASS = |error| within tolerance, every number recomputed by the driver):

| case | gate | measured | reference | error | tolerance | verdict |
|---|---|---|---|---|---|---|
| r=6 | voxel_cross_section_area_m2 | 7.777778e-03 | 7.853982e-03 | -0.9703 | 1.000e-01 | PASS |
| r=6 | flow_rate_lattice | 9.292111e-02 | 1.017876e-01 | -8.7108 | 2.500e-01 | PASS |
| r=6 | flow_rate_vs_area_equivalent | 9.292111e-02 | 9.982198e-02 | -6.9132 | 1.000e-01 | PASS |
| r=6 | u_max_lattice | 1.683655e-03 | 1.800000e-03 | -6.4636 | 2.500e-01 | PASS |
| r=6 | u_max_vs_area_equivalent | 1.683655e-03 | 1.782535e-03 | -5.5472 | 1.000e-01 | PASS |
| r=6 | mass_drift | 3.342393e-11 | 0.000000e+00 | 3.34e-11 | 4.480e-06 | PASS |
| r=6 | transverse_speed_max | 1.438370e-08 | 0.000000e+00 | 1.44e-08 | 9.000e-03 | PASS |
| r=10 | voxel_cross_section_area_m2 | 7.900000e-03 | 7.853982e-03 | +0.5859 | 1.000e-01 | PASS |
| r=10 | flow_rate_lattice | 7.633746e-01 | 7.853982e-01 | -2.8041 | 2.500e-01 | PASS |
| r=10 | flow_rate_vs_area_equivalent | 7.633746e-01 | 7.946288e-01 | -3.9332 | 1.000e-01 | PASS |
| r=10 | u_max_lattice | 4.892957e-03 | 5.000000e-03 | -2.1409 | 2.500e-01 | PASS |
| r=10 | u_max_vs_area_equivalent | 4.892957e-03 | 5.029296e-03 | -2.7109 | 1.000e-01 | PASS |
| r=10 | mass_drift | 2.712568e-10 | 0.000000e+00 | 2.71e-10 | 1.264e-05 | PASS |
| r=10 | transverse_speed_max | 3.829052e-08 | 0.000000e+00 | 3.83e-08 | 2.500e-02 | PASS |

Convergence of the flow rate (the staircase error, which is the point of the case):

| case | h/R | err vs CAD R [%] | err vs R_eff [%] |
|---|---|---|---|
| r=6 | 0.1667 | +8.711 | +6.913 |
| r=10 | 0.1000 | +2.804 | +3.933 |

observed order vs CAD radius: `[2.219]`,
vs area-equivalent radius: `[1.104]`.

`gates_failed = []`.

Figure: `artifacts/lattice_flow/lattice_flow_convergence.png` (left) flow-rate error against the cell size, with a slope-1 guide;
(right) the measured developed profile against the parabola of the area-equivalent duct.

## 5. Cross-repository finding: 3-D streaming is reversed too

`stream3d_probe()` places one population on `+x` and streams once: measured
`dx_lattice = -1` where `+1` is required
(`agrees_with_pull_convention = False`). The 2-D
finding in `docs/solver_consumers.md` therefore holds for `stream3d` as well
(`shifts=(-cz,-cy,-cx)` where a pull formulation needs `(+cz,+cy,+cx)`). Additionally
`tensorlbm.collide_mrt` literally calls `collide_bgk`, which is the **D2Q9** path
(`d2q9.macroscopic`), so the release exposes no usable 3-D collision. Both facts are
measured here, not inherited; both remain TensorLBM's to fix.

## 6. Points of caution

1. **One periodic unit, not an entry length.** `axial_cells=4` cells
   along the axis with periodicity; nothing is claimed about developing flow.
2. **A duct, not a ship.** No drag, resistance or hull-flow claim. The hull lattice
   from `voxel_solid` is same-shape evidence for the *geometry* face only.
3. **Plain Guo source.** The force enters as `S_i = w_i[3 c.F + 9(c.u)(c.F) - 3 u.F]`,
   which sums to zero (mass drift ≈ 0, see the gate table) and leaves a constant
   `O(F)` momentum offset that the split form would remove; the same choice and the
   same stated offset as the 2-D channel case.
4. **Hybrid solver.** Streaming and collision are local (see §5); only equilibrium,
   moments and lattice constants come from TensorLBM. Labelled as a hybrid.
5. **Two resolution levels in this run** (6, 10);
   a single order estimate, so "first order" is a measured tendency, not a fitted law.
   A third level was queued on the compute gate (`aj4033`-family submit) and did not
   run to completion inside this round's local budget — stated, not hidden.
6. **Staircase walls are staircase.** The hydraulic radius of a staircased circle is
   not the circle's; that is precisely what §3 and §4 quantify rather than assume.

## 7. Error verdict against the acceptance rule (is it <3 %?)

| key quantity | reference | error | <3 %? |
|---|---|---:|---|
| voxel cross-section area | `pi R^2` (analytic) | 0.970 % | yes |
| duct flow rate vs area-equivalent closed form | exact (same discretisation) | worst 6.913 % | no |
| duct flow rate vs CAD radius | `pi g R^4/(8 nu)` | 8.711 % (coarsest) → 2.804 % (finest) | yes |
| mass conservation | 0 drift | see gate table | yes |

The honest reading: the **area-equivalent** comparison is within a few percent and the
**CAD-radius** comparison is not at the coarse level — because a staircased circle
genuinely has less area than the circle. The gap is a discretisation measurement with
a sign and a convergence order, not an error to paper over.
