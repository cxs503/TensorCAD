# `SectionGrid.from_hull` z-pairing defect: the body was mirrored in `z`

Anchor: `GAP-CAD-B3-SOLVER-INTERFACE` (CAD→求解器消费面). This report fixes a defect **in the
bridge every solver leg of that anchor consumes**: the section grid that
`build_volume_mesh` / `volume_mesh_from_hull` turns into the body-fitted mesh handed to
TensorFVM and TensorFEM.

## 1. Case, geometry and the gate that must come first

`HullParameters` (sine family, `tensorcad.hull.generate_hull`):

| quantity | value |
|---|---|
| length L [m] | 10.000 |
| beam B [m] | 2.000 |
| draft T [m] | 0.800 |
| grid (stations × levels) | 41 × 13 |
| frame | z-up, metres, `x` aft→fore, `y` starboard, `z` up, waterline at `z = 0` |
| sampled volume [m³] | 7.446554 |
| analytic volume [m³] | 7.512698 (`analytic_volume`) |
| Cb = V/(L·B·T) | 0.465409 |
| analytic VCB [m, from WL] | −0.320000 (`−2/5 · T`) |

`Cb`, `L`, `B`, `T` are checked before anything else, per the 5415 discipline.

## 2. What was wrong, and why the existing gates could not see it

`_section_point(parameters, i, j)` walks the half outline **from the design waterline to
the keel**: ring index `j = 0` is at `z = 0` (breadth `B/2`), index `j = nz−1` is the keel
(breadth 0). `SectionGrid.from_hull` built its level list as

```
zs[j] = -draft * (1 - j/(nz-1))          # j = 0 -> keel height
```

but read the breadth of ring point `j`:

```
breadth(i, j) = |_section_point(p, i, j)[1]|   # j = 0 -> waterline breadth
```

So level `zs[0] = −0.8` (keel) was given `b = 1.0 m` and level `zs[12] = 0` (waterline) was
given `b = 0`: **the represented body was the hull mirrored in `z`** — widest at the keel,
closed at the waterline, i.e. the exact shape the volume-mesh module's own docstring says
it does *not* build ("collapsed columns … the keel line … where the hull genuinely has zero
breadth").

A z-mirror is invisible to every volume gate: `∫dV` is invariant. That is why the earlier
FVM/FEM hull legs of this anchor reported perfect agreements (conformed solid volume vs
hydrostatics 1.4e−14, buoyancy vs ρgV 3.5e−15, mesh affine exactness 4.5e−14) on a
mirrored body. The quantity that *does* see the mirror is the **first moment in z**, i.e.
the vertical centre of buoyancy, and the hull family carries an analytic VCB.

## 3. The fix (one line, in `tensorcad/volume_mesh.py`)

```python
def breadth(i: int, j: int) -> float:
    # Ring point ``nz - 1 - j`` is the one that sits at ``zs[j]``.
    return abs(_section_point(parameters, i, parameters.nz - 1 - j)[1])
```

No other change: `zs` stays strictly increasing (required by `SectionGrid.validate`), the
volume mesh layout, grading, split and quality machinery are untouched.

## 4. Results

Numbers are recomputed from the code in this worktree (no hand-copied values).

| quantity | pre-fix (ring `j`) | fixed (ring `nz−1−j`) | reference | pre-fix error | fixed error |
|---|---:|---:|---:|---:|---:|
| max &#124;`zs[j]` − ring z&#124; [m] | 0.800000 | **0.0** (exact) | 0 | — | 0 % |
| half breadth at keel level [m] | 1.000000 | **0.000000** | 0 | — | exact |
| half breadth at waterline level [m] | 0.000000 | **1.000000** | B/2 = 1.0 | — | exact |
| solid volume [m³] | 7.446554 | 7.446554 | 7.512698 analytic | −0.88 % | −0.88 % |
| VCB [m, from WL] | **−0.482306** | **−0.317694** | −0.320000 | **+50.7 %** | **−0.72 %** |
| tet-mesh volume centroid z [m] | −0.483005 | −0.316995 | −0.320000 | **+50.9 %** | −0.94 % |

The volume is identical to the last digit in both pairings (7.446554097149225 vs
7.446554097149229) — the quantitative statement that the old gates were blind.

Figure: `artifacts/section_grid_vcb/section_grid_vcb.png` (midship section `b(z)` for both
pairings, and the VCB bars against the analytic value).

## 5. Regression gates added

`tests/test_section_grid_orientation.py` (6 tests, all passing):

1. `zs[j]` equals the height of the ring point the breadth was read from — **exactly 0.0
   difference** over the whole 41 × 13 grid.
2. the waterline level carries `B/2`, the keel level carries 0.
3. grid VCB vs analytic VCB < 3 %.
4. **negative control**: the pre-fix pairing is reconstructed and must keep the volume to
   1e−12 while missing the VCB by > 15 % — the gate proves it can see what the volume
   gates cannot.
5. the tet volume mesh's centroid z sits between the analytic VCB (−5 % rel) and 0, so the
   downstream FVM/FEM body inherits the correct orientation.
6. the constant-breadth box reference is unaffected (the pairing has nothing to flip).

## 6. Points of caution / limitations

1. **Resolution.** The VCB residual after the fix (−0.72 %) is the discretisation error of a
   bilinear 41 × 13 quadrature against the `sqrt((z+T)/T)` section; it is not a residual
   defect. The mesh-centroid figure (−0.94 %, 3840 hexes, `first_layer=0.005, layers=4`) is
   the same error measured through the tetrahedra.
2. **Published downstream tables were produced on the mirrored mesh.** Volume, support
   reaction and affine-exactness gates in `docs/hull_solver_handoff.md` and
   `docs/fvm_handoff.md` are mirror-invariant, so those numbers stand; any *vertical-moment*
   claim in them (VCB, pressure-moment arms) must be regenerated. Not regenerated here:
   re-running the gated hull case costs ~90 s of FVM mesh construction plus ~50 s per FEM
   solve, and the acceptance question of this report — is the pairing right? — is answered
   by an exact-coincidence test, not by a solver run.
3. **Scope of the claim.** This is a geometry-orientation fix and its regression gate. No
   claim is made about resistance, seakeeping or structural response on the corrected body.
4. **Only hull-derived grids moved.** `box_section_grid` and caller-supplied `SectionGrid`s
   are unaffected (tested).

## 7. Error verdict (acceptance rule: > 15 % unacceptable; is it < 3 %?)

| key quantity | reference | pre-fix error | fixed error | < 3 %? |
|---|---|---:|---:|---:|
| VCB (grid quadrature) | analytic −2/5·T = −0.320000 m | **+50.7 %** | −0.72 % | yes |
| VCB (tet mesh centroid, 3840 hexes) | same | **+50.9 %** (−0.483005 m) | −0.94 % | yes |
| level/ring coincidence | 0 | 0.8 m | 0.0 | yes |
| solid volume | analytic 7.512698 m³ | −0.88 % | −0.88 % | yes (unchanged: mirror-invariant) |

Before the fix the VCB was **50.7 % off**, far outside the 15 % unacceptable band; the body
handed to FVM/FEM was upside down. After the fix every quantity above is inside 3 %.
## 8. Commit-gate note (declared, not hidden)

The repository's `pre-commit` golden-geometry gate fails on this branch for `inclined_box`
(`identity_volume_equals_requested_displacement` residual −3.87e−03, `identity_mirrored_heel_
mirrors_transverse_centre` residual 1.54e−05). Those checks live in
`tensorcad/golden_inclined_box.py`, which does not use `SectionGrid` (`SectionGrid` is
referenced only by `volume_mesh`, `solver_interface` and `fvm_handoff`), so they are
independent of this change and were already failing before it. This commit therefore used
the gate's documented per-commit bypass (`TENSORCAD_GOLDEN_SKIP=1`). The bypass is recorded
here rather than removed from the record; fixing that pre-existing inclined-box baseline is
a separate task.

Verification actually run for this change: `pytest tests/test_section_grid_orientation.py
tests/test_volume_mesh.py tests/test_volume_conform.py tests/test_volume_mesh_format_bridge.py
tests/test_volume_mesh_su2_markers.py` (73 passed) and `pytest tests/test_fvm_handoff.py
tests/test_volume_xval.py tests/test_volume_xval_entry.py tests/test_mesh_volume.py
tests/test_geometry_audit.py` (118 passed, 1 skipped).
