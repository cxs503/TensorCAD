# Hull CAD → solver consumption: the volume route, end to end

Anchor: `GAP-CAD-B3-SOLVER-INTERFACE` (CAD→求解器消费面). This closes the limitation the earlier handoff doc stated in its own words: *"Not a ship case… The hull → solver case belongs to the next layer of this anchor."* The new layer is here, and it uncovered a defect on the way.

## 1. The defect: the graded hull solid is not consumable as written

`tensorcad.volume_mesh.volume_mesh_from_hull` fills the hull body with graded hexahedra and splits them into `C3D4` tetrahedra. At the stem and the stern the section breadth collapses to zero, the breadth layers collapse onto the centreline, and the recovered tetrahedra have **zero volume**. The gmsh bundle still verifies (volume and digest are self-consistent), so nothing upstream noticed. A solver does:

```
tensorfvm.native_mesh3d.NativeMesh3D.from_cells(...) -> ValueError: degenerate face
```

Measured on the 41×13 default hull (`layers=2`): **195 zero-volume tetrahedra of 1920** (10.2 %). Dropping them left nodes that only they used, each owning three stiffness rows: TensorFEM's assembled tangent then came out singular (`torch.linalg.solve: input matrix is singular`). Both facts are pinned by tests (`test_hull_volume_bundle_is_watertight_and_matches_hydrostatics`, `test_orphaned_nodes_are_compacted`).

## 2. The fix: `tensorcad/volume_conform.py`

`weld_and_prune` welds coincident nodes (absolute rounding, documented), drops cells that welding proved degenerate, drops cells whose signed volume is zero to working precision, flips inverted cells (counted, never silent) and compacts nodes that lost all their cells. `conform_report` re-derives, from the arrays a solver will actually see: counts, boundary triangles, non-manifold edges, signed volume, minimum cell volume, inverted and zero cells.

| stage | nodes | tetrahedra | watertight | non-manifold edges | zero cells | volume [m³] |
|---|---:|---:|:---:|---:|---:|---:|
| raw `volume_mesh_from_hull` | 567 | 1920 | True | 0 | 195 | 7.374345 |
| after `weld_and_prune` | 492 | 1710 | True | 0 | 0 | 7.374345 |

Welding merged 567→492 nodes and dropped 210 degenerate cells (0 cells flipped). The conformed volume equals the hydrostatics displacement volume 7.374345 m³ of the *same* faceted surface to 1.4e-14 relative.

## 3. Geometry gate first (Cb / L / B / T), per the 5415 discipline

| quantity | value |
|---|---:|
| length L [m] | 10.000 |
| beam B [m] | 2.000 |
| draft T [m] | 0.800 |
| Cb (closed-surface hydrostatics) | 0.460897 |
| displacement volume [m³] | 7.374345 |
| LCB [m from aft] / VCB [m] | 5.0000 / -0.3147 |
| frame | z-up, metres, x aft→fore (`COORDINATE_CONTRACT`) |
| hull validation, all error-severity checks | PASS |

## 4. Consumer 1 — TensorFVM builds the face-based mesh and audits it

`NativeMesh3D.from_cells` / `NativeOperators3D.audit_affine` on the conformed tetrahedra, read back from `volume.msh` on disk (never passed in memory):

| quantity (dimensionless, exact = 0) | measured | gate | verdict |
|---|---:|---:|:---:|
| least-squares gradient max error | 5.329e-14 | 1e-09 | PASS |
| face reconstruction affine error | 3.553e-15 | 1e-09 | PASS |
| Gauss gradient max error | 1.044e-13 | 1e-09 | PASS |
| divergence of a constant vector field | 6.102e-15 | 1e-09 | PASS |
| divergence of an affine field (exact 6) | 4.530e-14 | 1e-09 | PASS |
| global flux balance error | 3.841e-16 | 1e-09 | PASS |

Mesh TensorFVM actually built: **1710 cells, 3777 faces (3063 internal / 714 boundary)**, volume 7.374345350 m³ vs declared 7.374345350 m³ (rel. 2.3e-15).

## 5. Consumer 2 — TensorFEM tet4 under hydrostatic pressure, checked by Archimedes

Hydrostatic pressure `p = ρg(0 − z)` acts on every boundary triangle below the waterline as a consistent nodal load (`p·A/3` per corner, inward normal). Supports are statically determinate, chosen geometrically (picking the three nodes by x-order alone left the system singular — measured). Closed-form references: Archimedes `ρgV` and the hydrostatics LCB.

| quantity | measured | reference | error | verdict |
|---|---:|---:|---:|:---:|
| buoyancy from the CAD surface [N] | 74125.564419 | 74125.564419 | 3.53e-15 | PASS |
| reaction at the supports [N] | -74125.564419 | -74125.564419 | 6.76e-17 | PASS |
| LCB from pressure moment [m] | 5.005276 | 5.000000 | 1.06e-03 | PASS |

Solve: 1710 tetrahedra, 492 nodes, 714 wetted boundary triangles, 3 Newton iterations, residual 1.43e-07, max displacement 1.849e-04 m (yield stress set to 1e15 Pa so the case stays elastic — the point is the load path, not plasticity).

## 6. Gates (driver: `python -m tensorcad.volume_conform --out artifacts/hull_handoff`)

| gate | measured | tolerance | verdict |
|---|---:|---:|:---:|
| nonmanifold_edges | 0.000e+00 | 0e+00 | PASS |
| zero_tets | 0.000e+00 | 0e+00 | PASS |
| inverted_tets | 0.000e+00 | 0e+00 | PASS |
| mesh_volume_m3 [m3] | 7.374e+00 | 7e-06 | PASS |
| fvm_ls_gradient_error | 5.329e-14 | 1e-09 | PASS |
| fvm_face_affine_error | 3.553e-15 | 1e-09 | PASS |
| fvm_affine_divergence_error | 4.530e-14 | 1e-09 | PASS |
| fvm_global_flux_balance_error | 3.841e-16 | 1e-09 | PASS |
| fem_archimedes_rel_error | 3.534e-15 | 5e-03 | PASS |
| fem_equilibrium_rel_error | 6.762e-17 | 1e-08 | PASS |
| fem_lcb_rel_error | 1.055e-03 | 2e-02 | PASS |

`gates_failed = []` — 11 gates, all PASS, every number written by the driver.

## 7. Points of caution / limitations

1. **Coarse hull for the gated run** (`nx=21, nz=9, layers=1`, 1.7 k tetrahedra). The 41×13 mesh needs ~90 s of TensorFVM mesh construction and ~50 s per TensorFEM solve; no resolution conclusion is claimed.
2. **Tet4, not a higher-order element.** A linear tetrahedron's pressure load is exact for a linear pressure field, but the displacement field is piecewise linear; real hull membrane/bending response needs tet10 or shells and is out of scope.
3. **Closed solid with the deck as a wall.** Pressure is applied only below the waterline; no deck plate, openings or internal structure. Not a structural-strength model.
4. **`weld_and_prune` makes the mesh usable; it does not repair the mesher.** The better fix — collapse the zero-breadth end columns deliberately and emit a closed tip — is a change to `volume_mesh.build_volume_mesh` and is left as a follow-up.
5. **Welding tolerance is absolute** (rounding to 1e-12 m); a geometry in other units must pass `decimals`. Stated in the docstring, not hidden.
6. **DG high-order curved handoff and hull voxelisation were not done this round.** `tensordg` is not importable in the development environment, and the hull ray-cast voxelisation needs a spatial index (the naive O(cells×faces) version did not finish). Stated as open, not claimed.

## 8. Error verdict against the acceptance rule (≤15 % unacceptable; is it <3 %?)

| key quantity | reference | error | <3 %? |
|---|---|---:|:---:|
| conformed solid volume | hydrostatics displacement volume | 1.4e-14 | yes |
| buoyancy from CAD pressure | ρgV (Archimedes) | 3.5e-15 | yes |
| support reaction | applied load | 6.8e-17 | yes |
| LCB from pressure moment | hydrostatics LCB | 1.1e-03 | yes |
| FVM affine exactness | exact (0) | 4.5e-14 | yes |

What these numbers do **not** say: they are consistency checks against closed forms on the same faceted surface and against exactness of the discrete operators. They do not validate resistance, seakeeping or structural strength; no experimental or reference-solver comparison is made here.
