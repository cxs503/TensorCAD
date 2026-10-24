# Representation benchmark

One harness, one ray set, one closed form, four representations. The question is not
"which representation is better" — it is **which task each representation is good at, and by
how much**, measured rather than argued. This is the code behind
[issue #19](https://github.com/cxs503/TensorCAD/issues/19).

```bash
python examples/benchmarks/run_benchmark.py --quick --output docs/benchmarks/quick   # CI scale
python examples/benchmarks/run_benchmark.py           --output docs/benchmarks/full # reference
python examples/benchmarks/run_benchmark.py --neural-bundle results/sphere-dent.ngo \
                                                      --output docs/benchmarks/neural
```

Exit code 1 means a gate failed. `results.json` is the machine-readable record;
`report.md` is generated from the same data by `tensorcad.benchmark` — no number in this
document was typed by hand, and every table below comes from the runs committed in
`docs/benchmarks/`.

## What is measured

Two cases, both with an exact reference:

- **Case A — `case_a_sphere_dent`.** A sphere R = 0.5 m with a Gaussian dent of depth
  0.06 m and σ = 0.22 rad on +x. The reference radius is a formula, so "the right answer" is
  known to machine precision. Rows: `analytic_field` (the formula, triangulated for export),
  `mesh_local_edit` (move the vertices inside 2.5 σ, nothing else), `neural_field` (a trained
  `.ngo` object, queried as an implicit field), `ocp_step` (faceted STEP round trip).
- **Case B — `case_b_hull_handoff`.** A design target → converged hull → the gmsh artifact a
  solver receives. The reference is `hull.analytic_volume`, not the mesh. Rows:
  `parametric_hull`, `ocp_step`.

Every row is probed along **the same stratified ray set**: a third of the rays lie inside the
edit cone (a uniform set puts only ~4 of 512 rays there, which is too few to measure a small
edit), the rest are uniform on the sphere. No row gets an easier sample, and no row gets a
different tolerance:

| gate | threshold | applies to |
|---|---|---|
| `mesh_edit_hits_the_target_depth` | axis depth error ≤ **1 % of 0.06 m** = 6.0e-4 m | the mesh row |
| `the_depth_threshold_is_the_same_for_every_measured_row` | the same 6.0e-4 m, unchanged | **any** row that measured a depth, including the neural one |
| `mesh_edit_is_local` | worst protected-ray leak ≤ 5 % of the edit depth | the mesh row |
| `every_ray_hits_every_surface` | coverage ≥ 0.999 | the mesh and neural rows |
| `analytic_row_is_the_reference` | the harness reproduces the formula to 1e-12 m | the analytic row |
| `no_representation_went_unmeasured_silently` | every row is `measured`, `imported` or explicitly `not_supported` | all rows |
| `hull_mesh_matches_the_analytic_volume` / `design_loop_reaches_the_target` / `export_round_trip_is_exact` | 1 %, convergence, exact re-read | case B |

Three definitions matter more than the thresholds:

- **Leak is measured against the same representation, not against the formula.**
  The mesh row's `protected_displacement_*` is the distance between the edited mesh and the
  *unedited mesh on the same ray*. Comparing against the analytic surface instead would
  report the chord sag of the triangulation (≈1.9e-4 m on the shipped grid) as if it were
  edit leakage — the kind of mistake that flatters or damns a representation for the wrong
  reason. Total error against the formula is reported next to it, under a different name.
- **The edit boundary is a staircase, so it is measured, not hidden.** A triangle straddling
  the 2.5 σ support edge carries the edit up to one cell further out. The protected region
  therefore starts 1.5 cells beyond the support edge, and the band in between is reported as
  `transition_band_displacement_m`. It is not a leak; it is where the boundary is.
- **Missing is not zero.** A metric that could not be measured carries `status` ∈
  {`not_supported`, `not_run`, `failed`}, `value = null`, a reason, and fails any gate it is
  attached to. A row whose extra is not installed reports `not_supported` and the command
  that would install it.

## The discretisation floor, and why the grid is 80×160

A lat-long triangulation is inscribed, so a ray between vertices lands on a chord *below*
the true surface. That floor is the smallest error any row measured against the formula can
show, and it is reported as the metric `discretisation_floor_m` (max chord sag outside the
edit, on the unedited sphere, at the same grid).

Measured ladder at 128 stratified rays — the mesh row's own error as a function of its grid:

| grid | axis depth error | cone mean depth error | dent surface error mean | volume error | floor | gates |
|---|---|---|---|---|---|---|
| 12×24 | 1.02e-3 m | 6.3e-4 m | 8.0e-4 m | 2.9e-2 | 7.9e-3 m | **fail** (depth) |
| 24×48 | 8.3e-4 m | 1.1e-4 m | 6.8e-4 m | 7.2e-3 | 2.1e-3 m | **fail** (depth) |
| 40×80 | 4.8e-4 m | 2.1e-5 m | 2.1e-4 m | 2.5e-3 | 7.4e-4 m | pass — but the floor is *above* the 6.0e-4 m gate |
| 60×120 | 2.6e-4 m | 5.0e-6 m | 7.3e-5 m | 8.3e-4 | 3.4e-4 m | pass |
| **80×160** | **1.4e-4 m** | **3.8e-6 m** | **4.7e-5 m** | **3.8e-4** | **1.9e-4 m** | pass, floor ≈ ⅓ of the gate |
| 120×240 | 2.1e-5 m | 1.1e-6 m | 1.9e-5 m | 1.2e-5 | 8.5e-5 m | pass |

40×80 happens to pass, which is exactly the reason it is not the shipped grid: its floor is
larger than the gate, so a green result there is luck, not margin. 80×160 is the coarsest
grid on which the floor is comfortably below the gate. `--quick` reduces the **ray count**
(512 → 128), never the grid: the grid sets the thing being measured, so it is not a
convenience knob. Case A takes 0.8 s at `--quick` and 2.9 s at full sampling.

## Case A: the same edit, three representations (full sampling, 512 rays)

From `docs/benchmarks/neural/results.json`:

| metric | `analytic_field` | `mesh_local_edit` | `neural_field` |
|---|---|---|---|
| status | measured (reference) | measured | measured, 1.73 s |
| dent depth, cone mean | 0.047199 m (reference, by construction) | 0.047199 m | 0.045292 m |
| depth error, cone mean | 0 | **2.3e-6 m** (0.004 % of depth) | 1.9e-3 m (**3.2 %**, gate fails) |
| depth error, ray nearest the axis | 0 | 1.3e-4 m | 1.8e-3 m |
| peak depth vs 0.06 m target | 0.0597846 m (reference) | 0.0597846 m | 0.0591424 m |
| dent surface error, mean / max | 0 | 4.5e-5 / 1.5e-4 m | 2.1e-3 / 5.3e-3 m |
| protected leak (vs its own parent) | 0 | **0 m** (exactly) | 1.3e-3 m (imported; parent weights not shipped) |
| protected total error (vs formula) | 0 | 1.3e-4 m | 4.1e-3 m |
| transition band | 0 | 2.0e-3 m (staircase) | 1.1e-2 m |
| surface error over every ray, mean / max | 0 | 1.0e-4 / 1.2e-3 m | 3.4e-3 / 1.5e-2 m |
| volume error | 6.7e-4 (its own triangulation) | 3.8e-4 | not measurable: a field has no faces to integrate |
| watertight / Euler χ | 1 / 2 | 1 / 2 | n/a (no mesh) |
| ray coverage | 1 | 1.0 | 1.0 |

The trained bundle's own validation claimed `dent_surface_error_mean_m = 2.32e-3 m` on its
own 4096-ray set; this harness measures 2.12e-3 m on a different ray set. The two agree, so
the artifact reproduces its claim — that cross-check is why the imported numbers are kept
next to the measured ones instead of replacing them.

**What this says, and what it does not.** On a local edit with a closed-form target, moving
vertices on a triangulation is 13–840× more accurate here than the trained field, and its
leak outside the support is exactly zero because untouched vertices are untouched. The
network's cost is not only accuracy: 1.73 s to answer 512 rays (batched, 241 samples per ray
plus 48 bisection steps) against 0.05 s to build the analytic surface and 2.9 s to build and
intersect a 12,642-vertex mesh — the field sits between them and wins on neither axis. What
the field does have is a 236 KB artifact with no topology, a continuous guarantee between
samples, and an edit that is a weight update rather than a re-meshing; none of that is
scored by this benchmark, and none of it is claimed in the table above.

It is *also* true that this comparison is tilted in a way worth stating: the network was
trained against the same closed form the harness scores it with, so it cannot lose to a bad
target definition, while the mesh row is a single naive vertex move with no re-projection
and no fairing. Both effects are in the same direction — they exaggerate the network's
advantage — and it still loses on accuracy.

## Case B: a design target, a hull, and the artifact a solver reads

From `docs/benchmarks/full/results.json` (all gates green):

| metric | value | gate |
|---|---|---|
| mesh volume vs `analytic_volume` at (41, 13) | **0.880 %** | ≤ 1 % |
| design loop residual (displacement target) | 2.2e-16 | converged, 6 evaluations |
| design loop wall clock | 1.76 s | — |
| gmsh re-read: worst node displacement | **0 m** | exact |
| gmsh re-read: area difference | 6.8e-16 relative | exact |
| grid refinement (41,13) → (81,25) | +0.577 % volume | ≤ 5 % |
| watertight / χ / wetted area | 1 / 2 / 51.99 m² | — |
| STEP round trip (`ocp_step`) | `not_supported` here (no OCP) | — |

The `ocp_step` row is `not_supported` in every committed run because this machine has no
OpenCASCADE wheel. It is a real row, not a placeholder: with `pip install "tensorcad[ocp]"`
it exports the faceted solid and re-reads it, reporting facet count beside the round-trip
volume, because a faceted solid is not a feature B-rep and the row must not be read as one.

## Threats to validity

- **One seed, one artifact.** The neural row reflects a single trained bundle
  (`weights.sha256 db1df71e…`, reproducible bit-for-bit: re-running the example produced the
  same hash). No variance across seeds is reported, so no significance claim is made.
- **The parent is not shipped.** The `.ngo` format does not retain the parent revision's
  weights, so the neural leak can only be imported, never re-measured. Adding a parent
  payload would make that column measurable and is the single change that would most improve
  this benchmark.
- **Means over rays, not surface integrals.** The numbers are unweighted ray means; rows are
  comparable to each other, not to a flux over the surface.
- **One task per representation.** A local dent and a hull-volume round trip do not test
  topology changes, booleans, or feature history — exactly where a B-rep should win and where
  this harness has nothing to say yet.
- **The neural path needs two optional extras.** `tensorcad[neural]` for training and
  measurement; without it the row degrades to `imported` and says so.

## Reproducing

```bash
python -m pip install -e ".[test,dev]"          # core benchmark, no torch, no OCP
python examples/benchmarks/run_benchmark.py --quick --output /tmp/bench   # ~4 s, all gates

python -m pip install -e ".[neural]"            # the neural row becomes measurable
python examples/neural_geometry/edit_sphere_dent.py --output results/sphere-dent.ngo   # ~3 min CPU
python examples/benchmarks/run_benchmark.py --neural-bundle results/sphere-dent.ngo \
                                            --output /tmp/bench-neural
```

That last command **exits 1**, and that is the result: the learned field misses the mesh
row's depth gate by 3.2×, and the same threshold was applied to both rather than inventing a
looser one for the network.
