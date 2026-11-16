# One solid from the unified geometry library → FVM / DG / LBM, one machine-checked gate

**Anchor**: `GAP-CAD-B3-SOLVER-INTERFACE` (CAD→求解器消费面：watertight→FVM 体网格 / FEM / DG 高阶曲边 / LBM 体素化)
**Module**: `tensorcad/solver_solid_gate.py` · **Driver**: `python -m tensorcad.solver_solid_gate`
**Tests**: `tests/test_solver_solid_gate.py` (22 tests)
**Artifacts**: `artifacts/solver_solid_gate/{gate.json,summary.md,fvm/,dg/,lattice/,dg_consumer/}`

## 1. Why this round exists: an audit of the claim, not a restatement of it

The anchor's acceptance is a *cross-route* sentence: "三求解器各一个端到端算例真跑通（坐标系契约统一
z-up）；输入=统一几何库输出". `tensorcad/solver_integration.py` claimed exactly that, and an
independent read of its code contradicted four parts of the claim:

| # | claim in `docs/solver_integration.md` | what the code did |
|---|---|---|
| A1 | "all three routes on the same CAD solid" | FVM used an ellipsoid `(a,b,c)`; DG was handed a **circle of radius a** (`_dg_route`, not a section of that ellipsoid, whose equator is an `a×b` ellipse); LBM voxelised a **cylinder radius a, length 2c**, whose volume `2πa²c` is `3a/(2b)` = **3×** the ellipsoid at the default `a=1, b=0.5` |
| A2 | "all signed volumes agree within 5 %" | `verify_frame_consistency` never compared a volume across routes; it compared each route's own metric with its own analytic reference (5 % / 1e-10 / 10 %) |
| A3 | "a tampered y-up manifest is refused by every route's verifier" | `frame_tamper_refused` ended in `except Exception: return True` — a *crash* counted as a refusal — and its DG branch tampered a bundle whose verifier had no frame field to refuse |
| A4 | "three routes ... really run" | every route returned `status="geometry_ready"`; the FVM leg never wrote a volume bundle and never built a solver mesh, the comment "Build a trivial tetrahedralisation" is followed by no tetrahedralisation |

Consequences that are now `defect_fix`-ed in code: `frame_tamper_refused` no longer
swallows exceptions (only `HandoffError` counts), and `docs/solver_integration.md`
carries a correction pointing here for the cross-route gate.

## 2. What this module does instead

**The solid comes from the unified geometry library, once.** `geometry_field.EllipsoidField`
is the CGC protocol object (`signed_distance` / `gradient` / `bounds` / `state`);
`geometry_field.to_mesh` is that library's own extraction path. The gate samples it
*once* and every route consumes that mesh or a section of it; the mesh SHA-256 travels
in every bundle's manifest, so "same solid" is a value that is compared, not prose.

Geometry gate first, in the order the 5415 rounds taught: the mesh must be closed and
consistently oriented (every directed edge used exactly once, its reverse exactly once),
outward-facing (positive signed volume), and within 5 % of `4/3 π a b c` — a bad fixture
is refused before a solver sees it.

| route | CAD side | solver side | reference |
|---|---|---|---|
| FVM | star-tetrahedralise the library mesh about its centroid → gmsh 2.2 `C3D4` + tagged wall (`write_bundle_volume`) → re-read by the independent parser | TensorFVM `NativeMesh3D.from_cells` + `NativeHybridScalar3D.solve` for `φ = 1 + 0.4x − 0.7y + 1.3z` (harmonic ⇒ exact) | `4/3 π a b c`, and the sampled surface volume |
| DG | `z = 0` section of the same solid as an **exact** rational-quadratic NURBS ellipse (`dg_handoff.ellipse_contract`, affine image of the tangent-intersection circle) → `write_dg_bundle` | TensorDG `curved_geometry.mesh_from_contract`, reading the bundle **from disk** (circular-section case, §4) | `4a·E(e)` from the θ-parametrised elliptic integral — a different parametrisation from the curve under test — and `πab` |
| LBM | `voxel_solid.solid_lattice` on the same mesh → occupancy lattice with true LBM semantics (fluid cells outside the solid, `write_bundle_lattice`) | voxel solver consumes `lattice.json` | `4/3 π a b c` (staircase tolerance) |

## 3. Results (`resolution 24`, semi-axes `(1.0, 0.5, 0.3)`, artifact `artifacts/solver_solid_gate/gate.json`)

Solid: 7794 nodes / 15584 triangles, digest `2f32851962c8fb1a…`, sampled volume
0.626048 m³ vs analytic 0.628319 m³ (**0.361 %**, fixture discretisation).

| route | measurement | vs reference | bundle re-read | consumer |
|---|---|---|---|---|
| FVM | 15584 tets, V = 0.626048 m³ | **1.6e-15** vs the sampled surface; 3.61e-03 vs analytic | PASS | **ran** (TensorFVM) |
| DG | ellipse a=1.0, b=0.5; P = 4.844224 m, A = 1.570796 m² | **7.3e-16** (perimeter), **7.1e-16** (area), 2.2e-16 m off the solid | PASS | ellipse section refused by TensorDG; circle section **ran** |
| LBM | [40, 20, 12] @ 0.05 m, V = 0.629250 m³ | **1.48e-03** | PASS | occupancy bundle |

Gate checks (all PASS; one WITHHELD, §4):

| check | measured | tolerance |
|---|---|---|
| frame_contract | `{z-up}` | z-up |
| same_solid_digest | one digest across all three manifests | 1 |
| fvm_tet_vs_surface_volume_rel_error | 1.596e-15 | 1e-12 |
| fvm_surface_vs_analytic_volume_rel_error | 3.614e-03 | 0.05 |
| fvm_patch_relative_linf / _l2 | 1.661e-08 / 6.475e-12 | 1e-06 / 1e-09 |
| fvm_global_relative_balance | 2.339e-08 | 1e-07 (fixture-bound, §5) |
| dg_perimeter / area / on_solid | 7.334e-16 / 7.068e-16 / 2.220e-16 m | 1e-09 / 1e-09 / 1e-12 |
| lattice_volume_rel_error | 1.482e-03 | 0.06 |
| lattice_fluid_plus_solid | 0 | 0 |
| cross_route_volume_rel_error | 5.097e-03 | 0.06 |
| dg_consumer_wall_radius_rel_error | 0.0 | 1e-09 |
| dg_consumer_ran | ran (CurvedMesh, degree 2, 4×16 elements) | ran |
| dg_consumer_accepts_non_circular_wall | **refused** | accepted → **WITHHELD** |
| frame_tamper_refused_all_routes | 3 / 3 | 3 |

## 4. The finding this round produced: the DG consumer only takes circles

TensorDG's `mesh_from_contract` builds an **annulus**: it samples the contract and
requires a constant wall radius. Handed the exact ellipse section, it refuses —
measured, verbatim:

```
ValueError: inner wall is not circular; the annulus builder needs a circular wall
(got radius spread 2.713e-01)
```

That is why the earlier integration handed it a circle and still called it "the same
solid". Two facts are now recorded separately instead of being conflated:

* **the CAD payload is exact** for a non-circular section (perimeter/area to 7e-16, every
  sampled point on the solid to 2e-16 m) — TensorCAD is not the limitation;
* **the consumer is** — until TensorDG's annulus builder accepts a general NURBS wall,
  a ship's non-circular waterline cannot reach it. Gate verdict for this check is
  `WITHHELD`, not PASS, so the debt stays visible.

The route that *can* be consumed is measured too: a spheroid `(0.5, 0.5, 0.3)` from the
same library pipeline has a circular equator, its bundle is written and verified, and
TensorDG builds a `CurvedMesh` from the file on disk with recovered wall radius
0.5 m (relative error 0.0, spread 2.2e-16).

## 5. Findings and limits, kept rather than smoothed

1. **Star tetrahedra are a hand-off fixture, not a production mesh.** They tile the
   surface exactly (1.6e-15) but are slivers: maximum non-orthogonality **88.7°**,
   minimum cell volume **1.5e-13 m³**. The plain orthogonal hybrid operator stays exact
   for a linear field (patch L∞ 1.7e-08, L2 6.5e-12) yet loses global balance at
   **2.3e-08**, against 3.8e-16 for the hull volume mesh in `docs/fvm_handoff.md`. The
   tolerance `1e-07` is declared as *fixture-bound* — it tracks cell quality, not the
   solver.
2. **`(a, c, b, centre)` is the right tet ordering.** Outward-facing triangles with
   `(a, b, c, centre)` give the *negative* of the cell volume in this repository's
   determinant convention (measured, then pinned by
   `test_an_inward_oriented_surface_is_refused`).
3. **Lattice error is staircase error**: 1.5e-03 at 0.05 m cells on a 1×0.5×0.3 m
   ellipsoid, i.e. a quarter of a cell layer, not a solver defect.
4. **Not claimed**: that the FVM solve here is a convergence study (it is an exactness
   patch test; the ladder is in `docs/fvm_handoff.md`), and that the hull-shaped solid
   traverses this gate (the hull volume route is `docs/hull_solver_handoff.md`).
5. **Environment**: `tensorfvm` resolves to another worktree's editable install; the
   resolved path is recorded in `gate.json` (`routes.fvm.consumer.solver_file`).
   `tensordg` is pinned via `TENSORDG_SRC` or the ecosystem checkout.

## 6. Reproduce

```
python -m tensorcad.solver_solid_gate --out artifacts/solver_solid_gate --resolution 24
python -m pytest tests/test_solver_solid_gate.py -q     # 22 passed in ~14 s
```

`--resolution 40` runs the same gate on a finer fixture (slower: the star decomposition
has ~43 k tetrahedra and TensorFVM's topology build is the cost).