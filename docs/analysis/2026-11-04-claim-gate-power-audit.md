# Gate power of TensorCAD's documented claims

**Date**: 2026-11-04 · **Scope**: `tensorcad/{hydrostatics,bonjean,stability,validation,mesh}.py`
· **Harness**: `runs/2026-11-04-claim-gate-power/mutate.py` (rerun it; the page is generated
from its JSON by `report.py`)

## Why this audit

Every number in this repository is supposed to arrive with the check that produced it. That
claim answers a different question from the one this page asks. "There is a test" and "the
test would fail if the number were wrong" are not the same statement, and the difference is
exactly the size of the verification debt a mature module accumulates: a test written after
the implementation tends to re-assert whatever the implementation does.

So the claims were attacked instead of read. 14 single-defect mutants were
injected into the five modules that carry the package's documented numbers -- wrong constants,
a dropped integration term, a disabled refusal, and shifts of 1e-09 and 0.5 % placed
*inside* the accuracy a doc claims -- and the test files that are supposed to protect the
claim were run against each one. A mutant the suite does not notice is a documented number
with no teeth.

## Method

* **One defect at a time.** The source file is patched by exact-string replacement (the
  harness refuses a pattern that is not unique), one pytest run per mutant, and the file is
  restored from its original bytes with a SHA-256 check before the next mutant starts.
* **Clean-tree guard.** The harness aborts if `git diff -- tensorcad` is non-empty, because an
  interrupted run would leave one mutant in place and silently contaminate the next (this
  actually happened during development -- `STA-1` came back `pattern-not-unique` -- which is
  why the guard exists).
* **Targeted scope.** Each mutant runs only the test files that guard its claim, so a kill is
  attributable. Baselines (unmutated) are recorded for the same file sets:

| test files | baseline | wall time |
|---|---|---|
| `tests/test_analytic_hydrostatics.py | tests/test_hydrostatics.py` | 29 passed, 1 warning in 13.75s | 14.83 s |
| `tests/test_bonjean.py` | 39 passed, 1 warning in 6.00s | 6.95 s |
| `tests/test_stability.py` | 19 passed, 1 warning in 4.49s | 5.62 s |
| `tests/test_validation.py | tests/test_gate.py` | 26 passed, 1 warning in 2.92s | 4.34 s |
| `tests/test_mesh.py | tests/test_validation.py` | 35 passed, 1 warning in 0.31s | 1.18 s |

* **Four mutant kinds.** `value` (wrong constant or missing factor), `model` (a dropped or
  altered term in an integral), `gate` (a disabled refusal or a widened predicate),
  `tolerance` (a shift smaller than a documented accuracy band).

## Result

Gate power before the cure: **12/14 =
85.7 %**. After the two tests added at the end of this page:
**14/14 = 100.0 %.**

| mutant | kind | file | documented claim it should protect | before | after |
|---|---|---|---|---|---|
| `HYD-1-tetra-sixth` | value | `tensorcad/hydrostatics.py` | parametric hull displacement 7.446554 m^3 at 41x13 vs analytic 7.512698 (-0.88 %) | killed | killed |
| `HYD-2-moment-arm` | model | `tensorcad/hydrostatics.py` | VCB = -2/5*T and the LCB/VCB values reported by the mesh integrator | killed | killed |
| `HYD-3-waterplane-flux` | value | `tensorcad/hydrostatics.py` | waterplane area 14.061669 m^2 at 41x13, and the closed-surface flux identity | killed | killed |
| `HYD-4-volume-half-percent` | tolerance | `tensorcad/hydrostatics.py` | the documented 0.88 % / 0.136 % volume grid errors (a 0.5 % shift must be visible) | killed | killed |
| `BON-1-shoelace-half` | value | `tensorcad/bonjean.py` | barge section area is B*T with error 0, bit-exact | killed | killed |
| `BON-2-shoelace-1e-9` | tolerance | `tensorcad/bonjean.py` | the 'bit-exact' wording of the barge section area | killed | killed |
| `STA-1-parallel-axis-T` | model | `tensorcad/stability.py` | box I_T = 2.2e-16 and BM = I/nabla to 0.75 % | **SURVIVED** | killed |
| `STA-2-parallel-axis-L` | model | `tensorcad/stability.py` | box I_L bit-exact; hull I_L 76.41 m^4 vs quadrature 76.8708 (-0.60 %) | killed | killed |
| `STA-3-area-1e-9` | tolerance | `tensorcad/stability.py` | the 'bit-exact' wording on the box waterplane area | killed | killed |
| `VAL-1-watertight-gate` | gate | `tensorcad/validation.py` | an open surface is refused ('the surface is open') | killed | killed |
| `VAL-2-euler-gate` | gate | `tensorcad/validation.py` | V - E + F = 2 for a closed genus-0 shell | killed | killed |
| `VAL-3-manifold-gate` | gate | `tensorcad/validation.py` | an edge shared by three faces is refused (edge_manifold) | killed | killed |
| `MSH-1-watertight-predicate` | gate | `tensorcad/mesh.py` | is_watertight is True only for a closed surface | killed | killed |
| `MSH-2-manifold-predicate` | gate | `tensorcad/mesh.py` | is_edge_manifold refuses an edge shared by three faces | **SURVIVED** | killed |

![gate power](gate_power.png)

*Left*: one bar per mutant, green = the suite turned red (killed), red = survived. *Right*: the
same data aggregated by module, before and after the cure.

### What survived, and why

* `STA-1-parallel-axis-T` -- box I_T = 2.2e-16 and BM = I/nabla to 0.75 %
* `MSH-2-manifold-predicate` -- is_edge_manifold refuses an edge shared by three faces

* **`STA-1`** removes the parallel-axis transfer `A * y_c^2` from the transverse second moment
  of the waterplane. It survived because *every* stability fixture in the suite -- box, hemisphere,
  parametric hull -- has its waterplane centred on `y = 0`, where the dropped term is identically
  zero. The mutation was invisible for a structural reason, not a tolerance reason: the suite never
  looked at a waterplane whose centroid is off the centreline, which is precisely the case
  (asymmetric loading, a heeling hull) where the term matters.
* **`MSH-2`** widens `is_edge_manifold` from `count <= 2` to `count <= 3` per edge. It survived
  because across the whole suite `is_edge_manifold` is only ever asserted **True** -- twelve call
  sites, all of them "a good mesh is manifold" -- and the `edge_manifold` refusal in
  `validation.validate_mesh` is computed independently and stays correct. The predicate itself had
  no negative control.

### The cure

Two tests, appended by `runs/2026-11-04-claim-gate-power/close_survivors.py`:

* `tests/test_stability.py::test_waterplane_second_moment_is_taken_about_the_section_centroid` -- a
  `6 x 2 m` box slid sideways by 3 m so the raw integral about `y = 0` (`112 m^4`) and the
  centroidal `I_T` (`4 m^4`) differ by a factor of 28.
* `tests/test_mesh.py::test_a_three_faced_edge_is_not_manifold` -- three triangles sharing one edge,
  asserting the face count, `is_edge_manifold is False`, `is_watertight is False` and six boundary
  edges.

With them in place the audit is re-run from scratch: **14/14**
mutants killed.

### What did *not* survive, i.e. the claims that are already pinned

Two results are worth stating because they are the ones a reader might expect to be weak:

* The `tolerance` mutants (`BON-2`, `STA-3`, `HYD-4`) are killed. A 1e-09 relative error added
  to the barge section area breaks 10 tests, and a 0.5 %
  shift of the displaced volume breaks 8 tests -- so
  "bit-exact" and "0.88 %" in the docs are measurements with teeth, not adjectives.
* The `gate` mutants are killed at the level of the *predicate* for `is_watertight` but the two
  `validation.py` refusals are killed through a single test each (`VAL-2`, `VAL-3`): thin, but real.

## Limitations

* **The mutant set is hand-picked, not exhaustive.** 14 mutants over five modules
  is a probe of the *style* of each module's evidence, not a coverage metric. A surviving mutant
  proves a gap; a killed mutant does not prove the absence of one.
* **`tolerance`-kind mutants measure the band, not correctness.** They show that a documented
  number is pinned tighter than the injected shift; they say nothing about whether the number is
  right (that is what the analytic references are for).
* **Not covered here**: `hull.py` (the generator's own closed forms are shared constants -- a
  mutant of the exponent would have to be paired with the analytic side to be meaningful),
  the neural layer (the `safetensors` extra is not installed in this environment, so
  `tests/test_cli_neural.py` and the neural suites do not collect), `mesh_conditioning.py`
  (35 s per run makes it a separate, deliberate audit), and `export.py`/`gate.py` thresholds.
* **No timing or performance claim** is made; wall times in the baseline table are one machine,
  one run.

## Reproduce

```bash
python runs/2026-11-04-claim-gate-power/mutate.py --out results_before.json   # 12/14
python runs/2026-11-04-claim-gate-power/close_survivors.py                   # (once; tests are committed)
python runs/2026-11-04-claim-gate-power/mutate.py --out results_after.json    # 14/14
python runs/2026-11-04-claim-gate-power/report.py                            # this page + figure
```
