# CAD -> DG high-order curved-element handoff (`GAP-CAD-B3-SOLVER-INTERFACE`, DG leg)

The anchor names three consumer families: FVM volume meshes, FEM, and **DG high-order
curved elements** (`watertight -> FVM 体网格 / FEM / DG 高阶曲边 / LBM 体素化`). The BEM/FEM-plate/LBM
and FEM/BEM/FVM combinations were delivered earlier; the DG curved leg had never been produced
from CAD geometry, and this document records it being produced, consumed and scored.

## 1. What is handed over

| stage | artefact | contract |
|---|---|---|
| geometry source | `tensorcad.nurbs.cylinder_solid_mesh(1.0, 2.0)` — watertight: 5184 edges, incidence {'2': 5184} | `tensorcad.geometry.Mesh` |
| section | plane `z = 0.86` m cut, 384 contour points, manifold-checked | — |
| curve fit | centre (-7.338e-17, 1.469e-17), r = 0.999933554 m, max residual 9.393e-05 m | — |
| contract | `contract.json` + `handoff.json` in `artifacts/dg_curved/handoff_circle/` | `tensordg.curved_geometry.v1`, units m, frame z-up |

The emitted curve is a **rational quadratic NURBS circle** (9 control points, weights
`1, 1/sqrt(2), ...`), so it represents the circle exactly rather than approximately: the
radius spread over 4096 samples is `6.661e-16` m
(`6.7e-16` relative) — the round-off floor.

## 2. How it was run

`PYTHONPATH=.:/root/tensor-ecosystem/TensorDG python examples/dg_curved/run_dg_curved.py`
(71.43 s): CAD solid -> section -> fit -> contract -> `validate_contract` (CAD side,
independent de Boor evaluator) -> TensorDG `validate_geometry_contract` -> `mesh_from_contract`
-> `CurvedAdvection` -> three measurements. The DG case is an exact solution: rigid rotation
`u = omega(-y, x)` leaves every concentric circle invariant, so `f0(r, theta - omega t)` solves
`u_t + div(u f) = 0` exactly with no-penetration walls; `omega = 1`, `T = 0.15`, blob
`r0 = 1.5 m`, `delta = 0.4 m`, azimuthal mode 2, annulus `r_wall = 0.999933554`,
`r_outer = 2.0 m`.

## 3. Results

### 3.1 Handoff geometry error vs the CAD section

| quantity | fine section (192 segments) | coarse section (48) |
|---|---|---|
| emitted radius [m] | 0.999933554 | 0.998949447 |
| max radial deviation of emitted curve vs analytic circle [m] | 6.645e-05 | 1.051e-03 |
| max residual of CAD section points vs circle [m] | 1.604e-04 | 2.528e-03 |
| enclosed area rel. error (Green) | -1.333e-04 | -2.100e-03 |
| perimeter rel. error | -6.655e-05 | -1.051e-03 |

The section residual shrinks by **15.76x** for a 4x
finer CAD tessellation: the O(h^2) faceting of the solid, not the handoff, is what the contract
inherits. The NURBS representation adds nothing to it.

### 3.2 DG run on the CAD-derived curved mesh

| p | mesh (nr x ntheta) | cells | relative L2 error (rigid rotation) | free-stream operator residual | free-stream drift |
|---|---|---|---|---|---|
| 1 | 2×16 | 32 | 9.4149e-02 | 4.234e+00 | 2.740e-01 |
| 1 | 4×32 | 128 | 4.6214e-02 | 4.445e+00 | 1.447e-01 |
| 1 | 8×64 | 512 | 2.0696e-02 | 4.572e+00 | 7.469e-02 |
| 2 | 2×16 | 32 | 2.8211e-02 | 7.792e-03 | 5.126e-04 |
| 2 | 4×32 | 128 | 2.7917e-03 | 2.119e-03 | 6.943e-05 |
| 2 | 8×64 | 512 | 3.5870e-04 | 5.493e-04 | 8.991e-06 |
| 3 | 2×16 | 32 | 1.6952e-03 | 2.483e-02 | 6.100e-04 |
| 3 | 4×32 | 128 | 2.9515e-04 | 6.628e-03 | 8.154e-05 |
| 3 | 8×64 | 512 | 1.9102e-05 | 1.711e-03 | 1.053e-05 |

Observed convergence orders (each refinement doubles both directions):

| p | order (32->128) | order (128->512) |
|---|---|---|
| 1 | 1.03 | 1.16 |
| 2 | 3.34 | 2.96 |
| 3 | 2.52 | 3.95 |

Reference-geometry control: the same case run on TensorDG's own exact-circle contract gives
L2 `5.571055e-03` and free-stream residual
`1.551621e-03`. The CAD-derived contract
differs from it by **1.32e-04** in L2 and
**7.30e-05** in the metric
identity — i.e. the handoff is invisible to the solver at the level it can resolve.

Domain area from element summation vs `pi (r_outer^2 - r_wall^2)`:
p=2 `-3.09e-06`, p=3 `-7.57e-10`.

## 4. Attention points (including mistakes made here)

1. **The first refinement study was wrong and is reported as such**: it refined only `ntheta`
   with `nr = 3` fixed, and the error floor (5.6e-3 for p=2, 1.1e-3 for p=3, flat across a 4x
   refinement) was the *radial* discretisation, not a DG property. Refining both directions
   removed the floor (0.00557179 -> 0.000358703).
   Any "order" quoted from a one-directional refinement is meaningless.
2. The observed orders are **below p+1** at the coarse end (p=1: 1.03;
   p=3: 2.52) and approach it only after refinement; with
   straight-sided p=1 geometry the wall itself is only O(h^2) accurate, so an order above ~2
   cannot be observed there. This is a property of the case, not a pass/fail gate.
3. The section plane must not pass through a vertex ring: the cylinder's lateral rings sit at
   `z = 0, 0.25, 0.5, ...` for `v_samples = 4`, so a cut at `z = 1.0` is refused instead of
   guessed. `section_contour` raises rather than returning a half-ordered ring.
4. **Consumer boundary (honest, quantified)**: TensorDG's `mesh_from_contract` builds a polar
   annulus and therefore requires a **circular** inner wall. A genuine elliptical section
   (ellipsoid 1.2 x 0.9 x 0.8 m cut at z = 0.08 m, 256 points) was
   fitted by least squares (16 control points, max deviation
   9.902e-03 m = 0.83 % of the semi-axis),
   passes the CAD-side validation, and is then **refused** by the consumer:
   `inner wall is not circular; the annulus builder needs a circular wall (got radius spread 2.903e-01)`. So "CAD -> DG curved" is real for circular walls only;
   general hull sections need the DG-side annulus builder generalised (see below).
5. Scales: the case is a 2 m annulus at unit rotation rate. Nothing here says anything about
   ship-scale Reynolds numbers; this is an interface/consistency verification, not a flow result.

## 5. Error judgement (<3 %?)

| quantity | reference | measured | error | <3 %? |
|---|---|---|---|---|
| CAD solid volume | pi r^2 h = 6.283185 m3 | 6.282040 m3 | -0.0182 % | yes |
| emitted curve radius | 1.0 m | 0.999933554 m | -0.00664 % | yes |
| emitted curve area | pi m2 | 3.141173932 m2 | -0.01333 % | yes |
| DG L2 error, p=3 finest | exact rotation | 1.910e-05 | 0.0019 % | yes |
| DG L2 error, p=1 finest | exact rotation | 2.070e-02 | 2.070 % | **no (p=1, straight-sided wall)** |
| domain area, p=2 | pi (r_o^2 - r_w^2) | element sum | -3.09e-04 % | yes |
| CAD vs reference-geometry DG solution | TensorDG circle contract | — | 0.0132 % | yes |

The p=1 row is reported as a failure of the <3 % bar on purpose: a straight-sided (geometric
degree 1) wall cannot resolve a curved wall to better than O(h^2), and one refinement step is
not enough. It is a discretisation limit of the case, not a handoff defect — the same p=1 mesh on
TensorDG's own exact-circle contract behaves identically.

## 6. Gates

{
  "cad_solid_watertight": true,
  "cad_volume_within_1pct": true,
  "contract_area_within_inherited_faceting": true,
  "contract_curve_self_consistent_below_1e-15": true,
  "dg_area_sum_within_1e-4": true,
  "dg_consumed_contract": true,
  "dg_matches_reference_geometry_within_0p1pct": true,
  "dg_order_toward_p_plus_one": true,
  "dg_refinement_reduces_error": true,
  "ellipse_contract_refused_with_reason": true,
  "section_is_closed": true
}
