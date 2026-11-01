# TensorCAD health audit (first intake) and one real defect absorption

*Date 2026-10-31 · repository `tensor-ecosystem/TensorCAD` · branch `agent/m1619` · base commit `ca0dde3` (== `main` == `origin/main`) · first mission of the 99 dispatch line for this repo (0 previous missions in history).*

## 1. Subject: what was audited, and against what

TensorCAD is the geometry foundation of the Tensor scientific-computing ecosystem: agent-callable parametric hull geometry, hydrostatics integrated from the artifact a solver actually receives, continuous signed-distance/implicit fields, SUBOFF (AFF-1/AFF-8) reference bodies, mesh conditioning, a TensorFEM handoff and a representation benchmark. It carries no fluid solver and no experimental fleet data, so the accuracy criteria available inside this repository are (a) analytic values of the generating functions, (b) the published DARPA SUBOFF offsets and predicates as transcribed in TensorLBM, (c) quadrature volumes of the published cross sections. Those are the baselines used below.

| Item | Measured |
|---|---|
| Commits on the branch history | 63 |
| HEAD | `ca0dde3cf` |
| `origin/main` (recorded) | `ca0dde3cf`, ahead 0 / behind 0 vs HEAD's base |
| Live `git ls-remote origin refs/heads/main` | reachable, upstream `main` = `ca0dde3cf`; object present locally |
| Declared core dependencies (`pyproject.toml`) | none (`dependencies = []`) |
| `import tensorcad` pulls in numpy? | no (measured by importing the package and inspecting `sys.modules`) |
| Live smoke value, README's own example | Cb = 0.4654, displacement volume = 7.446554 m³, validation passed = True |

The README's own claims reproduce on this machine: the default 41x13 grid example still returns Cb = 0.4654 and V = 7.446554 m³, and its own validation still passes (watertight / manifold / orientation / Euler checks). That is the `REAL` end of the evidence scale — no MOCK/STUB markers were found in the geometry core, and no example prints a number it did not compute.

### 1.1 Last 15 commit subjects (theme list)

| # | sha | date | subject |
|---|---|---|---|
| 1 | `ca0dde3` | 2026-10-25 | Merge pull request #21 from cxs503/arena/5b1c04cf-tensorcad |
| 2 | `f4fda27` | 2026-10-24 | Fix the CI benchmark assertion after the row rename |
| 3 | `4582c74` | 2026-10-24 | Align the benchmark with the merged contract (row IDs, P95, plot, provenance) |
| 4 | `8e5d362` | 2026-10-24 | Add the representation benchmark (0.9.0) |
| 5 | `ddee154` | 2026-10-24 | Merge pull request #20 from cxs503/benchmark/representation-comparison |
| 6 | `ed51ac8` | 2026-10-24 | docs: define unified geometry representation benchmark |
| 7 | `6502dc6` | 2026-10-24 | Merge pull request #18 from cxs503/arena/2bc34824-tensorcad |
| 8 | `45a605f` | 2026-10-24 | 0.8.0: a design loop (targets in, verified hull + solver manifest out) |
| 9 | `83ada0f` | 2026-10-24 | Merge pull request #17 from cxs503/agent/m241 |
| 10 | `336f96b` | 2026-10-24 | Merge pull request #16 from cxs503/arena/b0350018-tensorcad |
| 11 | `6c7fcf6` | 2026-10-24 | 0.7.0: mesh conditioning, FEM handoff, optional STEP bridge |
| 12 | `1f5f5a2` | 2026-10-24 | feat: DARPA SUBOFF appendages (AFF-8) from their published predicates (0.6.0) |
| 13 | `a2a2498` | 2026-10-24 | TensorCAD 扩编首波审计完成：基线验证 + backlog 登记 + 开发提案 |
| 14 | `19d68d2` | 2026-10-24 | feat: DARPA SUBOFF bare hull as an analytic field (0.5.0) |
| 15 | `54c5944` | 2026-10-24 | Merge pull request #15 from cxs503/arena/95320e8c-tensorcad |

Themes, in order of appearance: 0.5 bare hull → 0.6 SUBOFF appendages → 0.7 mesh conditioning + FEM handoff + optional STEP bridge → 0.8 design loop → 0.9 representation benchmark. Eight of the last fifteen commits are integration work (PR merges and the CI assertion they broke), which is the pattern to expect from a repository that is being driven by several agents at once: the physics/geometry lands first, the contract between the pieces lands second.

### 1.2 Directory structure vs the README's feature table

| README claim | Implementing modules | LOC | Test suite | Test LOC | Doc |
|---|---|---|---|---|---|
| What is implemented | `hull.py + hydrostatics.py + validation.py` | 1163 | `tests/test_hull.py` | 311 | -- (README only) (missing) |
| Continuous geometry core (0.4) | `geometry_field.py + geometry.py` | 933 | `tests/test_geometry_field.py` | 497 | docs/geometry-core.md (present) |
| SUBOFF bare hull (0.5) | `suboff.py` | 527 | `tests/test_suboff.py` | 305 | docs/suboff.md (present) |
| SUBOFF appendages (0.6) | `suboff_appendages.py` | 1153 | `tests/test_suboff_appendages.py` | 840 | docs/suboff_appendages.md (present) |
| Mesh conditioning (0.7) | `mesh_conditioning.py + mesh.py` | 1219 | `tests/test_mesh_conditioning.py` | 476 | docs/mesh_conditioning.md (present) |
| STEP and the TensorFEM handoff (0.7) | `adapters/fem.py + adapters/ocp.py + export.py` | 1071 | `tests/test_adapters_fem.py` | 278 | docs/step_export.md (present) |
| Design loop (0.8) | `design.py` | 749 | `tests/test_design.py` | 358 | docs/design_loop.md (present) |
| Representation benchmark (0.9) | `benchmark.py` | 1799 | `tests/test_benchmark.py` | 467 | docs/representation_benchmark.md (present) |
| Neural geometry objects (0.3) | `neural_geometry.py + neural/*` | 835 | `tests/test_neural_store.py` | 379 | docs/neural-geometry-format.md (present) |

Every claim in the README has a module, a test file and a document behind it; nothing in the tree is undocumented and nothing documented is absent. The two gaps are structural rather than rhetorical:

1. `docs/experiments/` and `docs/benchmarks/` hold artifacts (a learned `.ngo` bundle, benchmark runs) that the neural tests need; the neural pipeline itself needs `safetensors`, which this container cannot install (see §5).
2. `examples/` is empty of Python files even though the README's CI paragraph says `examples/parametric_hull/run_example.py` enforces principle 6 — the example exists in history but not in this checkout, so that particular CI claim is unverifiable here (reported as WITHHELD, not as a failure).

## 2. Method

1. **Inventory** from the artifacts themselves: `git log`, `git rev-list --left-right` against `origin/main`, a `git ls-remote` probe, and a walk of `tensorcad/` / `tests/` / `docs/` cross-referenced against the README's `##` sections. Nothing was taken from the README on trust.
2. **Test baseline, per file, with a three-way classification.** Every `tests/test_*.py` was run in its own process (`runs/2026-10-31-sail-ulp/baseline.py`), the log kept, and each file labelled `green` / `red` / `env-red`. `env-red` means every failure in that file died on `ModuleNotFoundError` for an optional dependency — an environment fact, not a library fact. The classifier refuses to call a file env-red if even one failure is a real assertion failure, so the label cannot be used to hide a defect.
3. **Defect triage.** The single real red was a bit-level disagreement between the pure-Python and NumPy evaluation paths of the AFF-8 sail. It was localised by differential sweeping rather than by reading: evaluate the same probes both ways, count mismatches, then bisect the arithmetic by substituting one operation at a time. `runs/2026-10-31-sail-ulp/audit.py` does this against the pre-fix implementation *imported from the base commit's blob*, so "before" is the shipped code and not a paraphrase of it.
4. **Fix, then re-measure.** After the fix the same sweep must be exactly zero, the existing pinned tests must stay green, and the pre-fix red must still be reproducible by stashing the fix (recorded in `baseline_prefix/`).

## 3. Results

### 3.1 Test baseline (post-fix tree)

| Metric | Value |
|---|---|
| Test files run | 21 |
| green | 18 |
| red (real failures) | 0 |
| env-red (missing optional dependency) | 3 |
| passed | 382 |
| failed | 2 |
| skipped | 11 |
| errors (collection) | 2 |

| File | passed | failed | skipped | error | class | missing module |
|---|---|---|---|---|---|---|
| `tests/test_adapters_fem.py` | 20 | 0 | 0 | 0 | green | — |
| `tests/test_adapters_ocp.py` | 1 | 0 | 10 | 0 | green | — |
| `tests/test_benchmark.py` | 26 | 0 | 1 | 0 | green | — |
| `tests/test_cli.py` | 9 | 0 | 0 | 0 | green | — |
| `tests/test_cli_neural.py` | 0 | 0 | 0 | 1 | env-red | safetensors |
| `tests/test_design.py` | 19 | 0 | 0 | 0 | green | — |
| `tests/test_export.py` | 9 | 0 | 0 | 0 | green | — |
| `tests/test_extract.py` | 12 | 0 | 0 | 0 | green | — |
| `tests/test_geometry_field.py` | 47 | 2 | 0 | 0 | env-red | safetensors |
| `tests/test_hull.py` | 44 | 0 | 0 | 0 | green | — |
| `tests/test_hydrostatics.py` | 8 | 0 | 0 | 0 | green | — |
| `tests/test_mesh.py` | 21 | 0 | 0 | 0 | green | — |
| `tests/test_mesh_conditioning.py` | 22 | 0 | 0 | 0 | green | — |
| `tests/test_neural_geometry.py` | 13 | 0 | 0 | 0 | green | — |
| `tests/test_neural_pipeline.py` | 0 | 0 | 0 | 1 | env-red | safetensors |
| `tests/test_neural_store.py` | 30 | 0 | 0 | 0 | green | — |
| `tests/test_parameters.py` | 28 | 0 | 0 | 0 | green | — |
| `tests/test_public_surface.py` | 6 | 0 | 0 | 0 | green | — |
| `tests/test_suboff.py` | 19 | 0 | 0 | 0 | green | — |
| `tests/test_suboff_appendages.py` | 34 | 0 | 0 | 0 | green | — |
| `tests/test_validation.py` | 14 | 0 | 0 | 0 | green | — |

Signature of the environment red, stated explicitly so it is not confused with a defect:

* `tests/test_cli_neural.py`, `tests/test_neural_pipeline.py` — collection error: `ModuleNotFoundError: No module named 'safetensors'` (raised from `tensorcad/neural/field.py`, a module that only the neural checkpoints need).
* `tests/test_geometry_field.py` — 2 of 49 fail, both at `from tensorcad.neural.objects import field_from_object`, same missing module; the other 47 pass, including the four exact primitives.
* `tests/test_adapters_ocp.py` — 1 passed, 10 skipped: `cadquery-ocp` is absent, so the STEP bridge is exercised only through its fallback path. Reported as skipped, not as covered.
* `pip install safetensors` fails in this container while building its own build dependencies, and the GitHub remote is unreachable (`Failure when receiving data from the peer`), so the neural pipeline is **unverifiable here, not broken**. 3 of 21 files (the neural surface) are outside this audit's evidence.

### 3.2 The real red, before and after

Pre-fix baseline of the file that owns the defect (`baseline_prefix/test_suboff_appendages.log`): `2 failed, 32 passed, 1 warning in 68.68s (0:01:08)`

* `FAILED tests/test_suboff_appendages.py::test_the_scalar_path_walks_the_sail_entrance_segment`
* `FAILED tests/test_suboff_appendages.py::test_the_scalar_and_vectorised_sail_paths_agree_bit_for_bit`

Diagnosis, measured rather than assumed. The AFF-8 sail's implicit function exists twice — a pure-Python path (`_sail_implicit_ft`, used when NumPy is absent or when a caller evaluates one point) and a NumPy path (`_sail_implicit_array_ft`, what a sampler and every mesh operation actually calls). The module's own contract is that they are mirrors, and `test_the_scalar_path_walks_the_sail_entrance_segment` pins that with `==`. They were not mirrors: the cap-branch normalisation used `math.hypot` in one and `np.hypot` in the other, two different correctly-rounded algorithms that disagree in the last bit on ~3 % of the sail's bounding box. Substituting each operation in turn localised it — `math.hypot(1.0, 1.3182451600542016)` and `numpy.hypot` of the same arguments differ by exactly 1 ULP, while `math.sqrt(a*a + b*b)` agrees with `math.hypot` to the bit.

| Sweep resolution | Probes | pre-fix mismatches | post-fix mismatches | pre-fix max &#124;Δ&#124; [m] | pre-fix max rel. Δ |
|---|---|---|---|---|---|
| 5³ + boundaries | 173 | 2 | 0 | 6.94e-18 | 3.50e-16 |
| 7³ + boundaries | 391 | 20 | 0 | 1.04e-17 | 1.36e-15 |
| 9³ + boundaries | 777 | 26 | 0 | 6.94e-18 | 3.50e-16 |
| 11³ + boundaries | 1379 | 24 | 0 | 1.04e-17 | 1.24e-15 |
| 13³ + boundaries | 2245 | 68 | 0 | 1.04e-17 | 9.45e-15 |

What the fix moved (n=13 sweep): the two paths were rewritten onto one arithmetic sequence — `sqrt(a*a + b*b)` with a denominator that is only divided where the cap branch uses it — so 67 probes changed on the scalar side and 23 on the NumPy side. The largest change to any returned value relative to the pre-fix NumPy path is 3.50e-16 (1 ULP); the largest absolute discrepancy that the fix removed is 1.04e-17 m.

The parts a caller can build on did not move at all: containment labels identical (True), sign changes 0 over 2245 probes, and the published magnitude claim of the field (1.3e-05 m normalisation error) is untouched. A regression test now sweeps the whole bounding box plus the four axial mask boundaries and the cap plane ±1e-12 m, so the invariant is enforced on 2245 probes instead of the four hand-picked ones that happened to miss the cap branch's worst case.

### 3.3 Test-level verification

`tests/test_suboff_appendages.py`: 2 failed, 32 passed, 1 warning before the fix (2 red) → 34 passed after it, including the new `test_the_scalar_and_vectorised_sail_paths_agree_bit_for_bit`. Full suite after the fix: 382 passed, 2 failed (both env-red), 11 skipped, 2 collection errors (both env-red). `ruff check` on both changed files: clean.

## 4. Figure

![sail ULP divergence before/after](2026-10-31-sail-ulp-fix.png)

*Panel (a)*: every probe in the n=13 sweep where the two paths disagreed, in the sail's (x, z) plane; the dashed line is the cap plane, where the branch that carried the `hypot` normalisation begins. The mismatches cluster in the cap's failed-side band, which is why four hand-picked body probes almost missed them. *Panel (b)*: mismatch rate against sweep resolution — 1.4–3.0 % of probes over five resolutions before the fix, exactly 0 afterwards. *Panel (c)*: the size of the discrepancy on a log axis against the field's own documented accuracy: the disagreement peaked at 1.0e-17 m, i.e. 8.0e-13× the documented normalisation error. The defect was therefore invisible to any physical accuracy criterion and fatal only to the bit-identity contract — which is exactly why it needed a `==` test rather than a tolerance.

## 5. Notes: traps, limits, scope

1. **`env-red` vs `red` needs a machine rule, not a mood.** Three files look broken and are not; the classifier in `baseline.py` only grants `env-red` when *every* failure in the file traces to a missing optional dependency. Applying the label by eye would have buried the sail defect too, since that file's red was a 2e-16 assertion.
2. **Bit-level defects hide behind tolerances and behind hand-picked probes.** The existing four-point test caught a probe at index 1; the same probe on a differently shaped sweep disagrees at index 3. One probe is not a contract — a sweep is.
3. **`math.hypot` and `numpy.hypot` are not the same function.** Measured here: they differ by 1 ULP on ~3 % of a representative argument set. Any module that keeps a scalar and a vectorised implementation "in mirror" and normalises with `hypot` has the same latent defect.
4. **Measured follow-up, not fixed in this mission.** `CylinderField` (`tensorcad/geometry_field.py`) has the same scalar/NumPy split with `math.hypot` vs `np.hypot` and, unlike the sail, also differs in accumulation order (`delta @ axis` via BLAS, `np.einsum` for the radius). Over 4000 probes on a tilted cylinder, **926 (23 %) disagree, max relative 5.6e-14** — measured by `python3 -c` against the two paths of the shipped field. No test demands bit identity there, so it is not a red; fixing it is a wider change (BLAS/einsum order, not one line) and belongs in its own mission. Three further `math.hypot`/`np.hypot` pairs exist in `tensorcad/suboff.py` (lines 277/349, 323/442) and one in `tensorcad/geometry_field.py` (381/460); they were **not** individually probed here, so their mismatch rate is unknown rather than zero.
5. **The fix deliberately favours one arithmetic over two library calls.** `sqrt(a*a+b*b)` loses `hypot`'s overflow protection; the arguments are sail coordinates in feet, bounded by the model (|x| < 20 ft), so it is safe here. A probe at 10⁶ m from the origin still returns a finite value (verified), and the cap-branch denominator is now divided only where it is used, which removed the `RuntimeWarning: overflow in multiply` that a naive substitution introduced at degenerate probes (z = 0, y = 0).
6. **This is a geometry-library audit, so "error vs experiment" has no meaning here.** The quantities measured are field values and their agreement with each other and with the published DARPA predicates; no hydrodynamic quantity was computed, and no resistance or seakeeping comparison was attempted. That is a property of the package, not an omission of the audit.
7. **Nothing was pushed, and nothing was fetched.** The work sits on the mission branch as one commit; the platform's pusher owns the GitHub side. Upstream `refs/heads/main` was read with `git ls-remote` (reachable at report time, `ca0dde3cf` — identical to this branch's base), so the ahead/behind figure in §1 is a live comparison against the upstream tip, not a stale one; no ref was written locally by the probe.

## 6. Error judgement (explicit)

| # | Key quantity | Baseline / reference | After the fix | Error | < 3 % ? |
|---|---|---|---|---|---|
| 1 | Scalar-path vs NumPy-path sail signed distance (n=13 sweep, 2245 probes) | exact equality — the module's documented mirror contract, enforced by `==` | 0 mismatching probes | 0.0 (exact) | yes, by construction (0 %) |
| 2 | Same, pre-fix (the defect being removed) | as above | 68 mismatching probes, max 1.0e-17 m | 3.0 % of probes; 9.5e-15 relative | yes on any physical scale — and still unacceptable, because the contract is exactness |
| 3 | Value change vs the pre-fix NumPy path (what a solver would have seen) | pre-fix `signed_distance` | max relative change 3.5e-16 | 1 ULP | yes (≪ 3 %) |
| 4 | Zero-set topology (containment labels, sign) | pre-fix field | identical labels, 0 sign changes | 0 | yes (exact) |
| 5 | Published DARPA/TensorLBM-pinned numbers in `tests/test_suboff_appendages.py` | TensorLBM predicate transcription + cross-section quadrature | all 34 tests pass | unchanged within 1 ULP | yes |
| 6 | README smoke values (Cb, displacement volume, watertight) | README: Cb 0.4654, V 7.446554 m³ | Cb 0.4654, V 7.446554 m³, validation passed = True | Cb 0.0000 %, V 0.0000 % | yes |

Answers, plainly: **yes, every key quantity is below 3 %** — five of the six rows are exact or within 1 ULP, and the sixth reproduces the README to 4 decimal places of a percent. The substance of this mission was not a physics error: the sail's 1-ULP divergence was 8.0e-13× the field's own documented accuracy, i.e. measurable only by `==`. It was nevertheless a real defect, because the module promises a bit-identical mirror and a caller who samples and a caller who probes one point would have received two different numbers describing one surface.

## Outline of what changed

* `tensorcad/suboff_appendages.py` — `_hypot_ft` added; `_sail_implicit_ft` and `_sail_implicit_array_ft` rewritten onto one arithmetic sequence (sqrt of a sum of squares, denominator divided only in the cap branch). 43 lines touched, no API change, no constant changed.
* `tests/test_suboff_appendages.py` — `test_the_scalar_and_vectorised_sail_paths_agree_bit_for_bit` added: 2245-probe sweep (bounding box + both mask boundaries + cap plane ±1e-12 m) demanding exact equality; red before the fix, green after.
* `runs/2026-10-31-sail-ulp/` — `audit.py` (pre-fix code imported from the base commit's blob, differential sweep, figure), `baseline.py` (per-file baseline with env-red classification), `make_report.py` (this document), `metrics.json`, `baseline/summary.json`, `baseline_prefix/` red evidence, `sail_ulp_fix.png`.
* `docs/analysis/2026-10-31-tensorcad-health-audit.md` and `docs/analysis/2026-10-31-sail-ulp-fix.png` — this report.

