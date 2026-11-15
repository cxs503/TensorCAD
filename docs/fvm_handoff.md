# Watertight CAD solid -> tetrahedral volume mesh -> FVM

Anchor: `GAP-CAD-B3-SOLVER-INTERFACE` (CAD→求解器消费面). This document covers the
route the anchor names first and that until now had **no number behind it**:
`watertight surface -> FVM volume mesh`. The three other routes already had a
consumer (`docs/solver_consumers.md`: surface→BEM, 2-D grid→FEM plate,
lattice→LBM); the volume route was written and re-read by the contract verifier but
never solved on -- the FVM leg of the earlier round was a 2-D planar channel, which
an independent board audit recorded as finding **F5** ("锚点『水密面→FVM 体网格』规范路径未演示").

Module: `tensorcad/fvm_handoff.py`. Tests: `tests/test_fvm_handoff.py`.
Artifacts: `artifacts/fvm_handoff/` (`results.json`, `summary.md`, per-case bundles).

## 1. The chain, link by link

| link | who does it | what is machine-checked |
|---|---|---|
| CAD solid | `tensorcad.volume_mesh.build_volume_mesh(SectionGrid)` | a hull section grid (waterline `z=0`, keel `z=-draft`, half breadth single-valued in `y`) |
| write | `tensorcad.solver_interface.write_bundle_volume` | gmsh 2.2 ASCII, one `C3D4` block + one physically tagged `wall` triangle block, `handoff.json` with counts/volume/digest/frame |
| verify | `tensorcad.solver_interface.verify_bundle` | independent re-read: volume, tetrahedron count, digest, frame name; **hard gate, runs first** |
| parse | `tensorcad.solver_interface.read_gmsh_volume` | second reader, not the writer's memory |
| mesh | `tensorfvm.native_mesh3d.NativeMesh3D.from_cells` | the *solver* rebuilds faces, orients them, rejects non-manifold/non-convex cells, audits closure and face projection |
| solve | `tensorfvm.native_hybrid_scalar3d.NativeHybridScalar3D` | conservative cell-centred FV for `-div(kappa grad phi) = s`, Dirichlet on the tagged wall |

Every case starts from the file on disk. Nothing in this module passes in-memory
arrays to a solver.

## 2. References (closed-form, computed here)

| case | field | why it is a reference |
|---|---|---|
| `patch` | `phi = 1 + 0.3x - 0.7y + 0.2z` | linear ⇒ harmonic ⇒ a consistent FV operator must reproduce it **exactly** on any conforming mesh. This is a real patch test, unlike the earlier round's `patch_test_relative_linf` headline, which the board audit showed was measuring a boundary-layer incompatibility (a Poisson-ratio control experiment moved it by 985×, finding **F3**) |
| `harmonic` | `sin(pi x/L) sin(pi y/B) sinh(gamma (z+T))`, `gamma = pi sqrt(1/L² + 1/B²)` | exact harmonic function on `[0,L]×[0,B]×[-T,0]` by separation of variables, so the Dirichlet problem has a non-trivial exact solution; `L=2.0, B=0.4, T=0.3 m` |
| `frame` | the `patch` field applied to the same table read as y-up | prices the z-up/y-up mix-up the 5415 rounds paid for, as a number |

`test_harmonic_field_really_satisfies_laplace` differentiates the second reference
numerically and asserts `|lap phi| / (|phi| pi²) < 1e-3`, so the reference cannot
silently stop being harmonic.

## 3. Results (regenerated: `python -m tensorcad.fvm_handoff --out artifacts/fvm_handoff`)

Geometry: `barge` = box-section hull solid `2.0 × 0.4 × 0.3 m` (volume 0.24 m³);
`hull` = curved Wigley-derived hull `10 × 1 × 0.6 m` with a flat of the bottom and a
transom. Frame for all of them: **z-up**, metres, waterline `z = 0`.

| case | cells | check | value | bound / reference | verdict |
|---|---|---|---|---|---|
| barge_2x0.4x0.3_9x5x2 | 768 | patch_relative_linf | 4.48575e-14 | <= 1e-09 | PASS |
| barge_2x0.4x0.3_9x5x2 | 768 | mesh_volume_vs_declared | 0.24 | ref 0.24 (rel -1.14e-14) | PASS |
| barge_2x0.4x0.3_9x5x2 | 768 | harmonic_relative_l2 | 0.0419248 | <= 0.15 | PASS |
| barge_2x0.4x0.3_9x5x2 | - | zup_relative_linf | 4.48575e-14 | <= 1e-09 | PASS |
| barge_2x0.4x0.3_17x9x3 | 4608 | patch_relative_linf | 2.07143e-14 | <= 1e-09 | PASS |
| barge_2x0.4x0.3_17x9x3 | 4608 | mesh_volume_vs_declared | 0.24 | ref 0.24 (rel +1.38e-14) | PASS |
| barge_2x0.4x0.3_17x9x3 | 4608 | harmonic_relative_l2 | 0.0298341 | <= 0.15 | PASS |
| hull_10x1x0.6_13x7x2 | 1728 | patch_relative_linf | 2.52382e-12 | <= 1e-09 | PASS |
| hull_10x1x0.6_13x7x2 | 1728 | mesh_volume_vs_declared | 2.97565 | ref 2.97565 (rel +7.16e-15) | PASS |

Integrity / ledger numbers that make the above believable rather than lucky:

| quantity | barge 768 | barge 4608 | hull 1728 |
|---|---|---|---|
| harmonic relative L∞ | 5.308e-02 | 3.592e-02 | – (patch only) |
| harmonic relative L2 | 4.192e-02 | 2.983e-02 | – |
| patch relative L∞ | 4.486e-14 | 2.071e-14 | 2.524e-12 |
| solver linear residual | 6.65e-16 / 5.77e-16 | 1.05e-15 / 8.98e-16 | 7.22e-16 |
| solver global balance | 4.73e-14 / 2.08e-15 | 3.45e-13 / 6.35e-15 | 7.70e-14 |
| max non-orthogonality | 87.96° | 85.93° | 89.66° |
| wall `bad_edges` / Euler χ | 0 / 2 | 0 / 2 | 0 / 2 |
| mesh volume vs declared | −1.14e-14 | +1.38e-14 | +7.16e-15 |
| frame price (y-up / z-up) | 23.76 % / 4.49e-14 = **5.3e12×** | – | – |

(The two numbers on a line are the `patch` and `harmonic` solves respectively.)

## 4. Error judgement against the anchor's 3 % bar, stated item by item

| key quantity | reference | error | < 3 %? |
|---|---|---|---|
| patch field, box, 4608 cells | exact linear field | 2.07e-14 | **yes** (round-off) |
| patch field, curved hull, 1728 cells | exact linear field | 2.52e-12 | **yes** (round-off) |
| enclosed volume, all cases | `handoff.json` declaration | ≤ 1.4e-14 | **yes** |
| harmonic field, box, 4608 cells, L2 | separation of variables | 2.98e-02 | **yes**, marginally |
| harmonic field, box, 4608 cells, L∞ | separation of variables | 3.59e-02 | **no** — 3.6 % |

Honest reading: the handoff itself contributes nothing measurable (patch tests at
round-off, volume to 1e-14); the 3 % bar is met on the L2 norm of the harmonic case
and missed on L∞ by 0.6 percentage points. Observed order between 768 and 4608 cells
is `p = ln(4.19/2.98)/ln((4608/768)^{1/3}) = 0.57`, i.e. **slow**. Two causes are
identifiable from the recorded mesh metrics rather than guessed: maximum
non-orthogonality is 86-90° (skew tets out of the boundary-layer growth of
`build_volume_mesh`), and the solver operator used here is the plain orthogonal
approximation. A non-orthogonal-corrected operator, or a CAD-side mesh with less
skew, is where the remaining error lives — not in the handoff. No claim of a
convergence *order* is made from two resolutions.

## 5. Findings recorded here rather than hidden

1. **A collapsed mesh column breaks the volume handoff.** A hull section grid whose
   half breadth reaches zero on the keel line or at a pointed stem is accepted by
   `tensorcad.volume_mesh` (its builder treats it as a collapsed column), but the
   tet split then contains topologically degenerate faces and TensorFVM's auditor
   refuses the mesh with `degenerate face`. Dropping those cells would open the
   surface, so the correct fix is geometric: `hull_section_grid` keeps a flat of the
   bottom (`FLAT_KEEL_FRACTION = 0.10`) and a transom
   (`TRANSOM_FRACTION = 0.15`) — both ordinary features of a real hull.
   `test_a_collapsed_keel_column_is_refused_rather_than_silently_dropped` pins the
   behaviour so the finding cannot quietly disappear.
2. **`NativeMesh3D.from_gmsh` needs the gmsh SDK, which has no loongarch64 wheel in
   this environment.** The file is therefore parsed by this repository's reader and
   the arrays enter TensorFVM through `from_cells` — the same constructor `from_gmsh`
   calls after reading. *Not* claimed: that the solver opened the gmsh file itself.
   Claimed: the numbers the solver used are the ones on disk (digest-verified) and the
   topology is the solver's own audit.
3. **`import tensorfvm` on this host resolves to another agent's worktree** through
   editable-install `.pth` files (three TensorFVM checkouts are on `sys.path`).
   `_pin_solver_path()` therefore pins `TENSORFVM_SRC` or the ecosystem checkout, and
   `solver_fingerprint` records the resolved path and git sha in every result. A
   number without that record is not reproducible.

## 6. Reproduce

```
PYTHONPATH=.:/srv/tensor/workspace/TensorFVM/src python -m tensorcad.fvm_handoff --out artifacts/fvm_handoff
PYTHONPATH=.:/srv/tensor/workspace/TensorFVM/src python -m pytest tests/test_fvm_handoff.py -q
```

`--quick` runs one barge and one hull resolution. The full ladder above is two barge
resolutions (768, 4608 cells) and one hull resolution (1728 cells); a third rung
(24576 cells) was measured during development (`harmonic` L2 2.25 %) but is not in
the default artifact set because it pushes the single-process run past ten minutes.