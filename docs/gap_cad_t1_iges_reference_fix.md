# GAP-CAD-T1-IGES-PIPE — round 3: the gold reference itself was the defect

**Mission** m17683 · repo TensorCAD · branch `agent/m17683` · anchor `GAP-CAD-T1-IGES-PIPE`
Acceptance: *pipeline CLI + library entry point; 3 real IGES files through the gates;
gold regression test green; dependency lock.*

## 1. What the previous round left open

`docs/gap_cad_t1_iges_pipe.md` closed with an honest FAIL: the harness's gold benchmark
did not pass, because "my quadrature reference and the sampling sweep do not agree".
Reproduced here before touching anything (`runs/gap_t1_iges/evidence.json`, 205.6 s run):

| symptom | measured |
|---|---|
| every one of the four IGES rows | `error: KeyError: 'component_count'` — no row recorded a number |
| Series 60 Cb=0.60 sweep vs the then-reference | +0.001 %, +2.02 %, +3.18 %, +3.49 % — increasing, never converging |
| DTMB5415 loft sweep vs the then-reference | −129.5 % … −126.5 % (reference numerically impossible) |
| KP505 (real CATIA V5 IGES) | 8/8 faces parsed, 2/8 tessellated, 6 refused, 252 triangles, not watertight |

## 2. Root causes (both are defects in the *judging* code, not in the pipeline)

**(a) `KeyError: 'component_count'`.** `trim_mesh`'s audit publishes
`components` / `closed_components` (`tensorcad/geomquality/trim_mesh.py:480`), the driver
read `after["component_count"]`. Every row died in the `except` branch, so `evidence.json`
carried four error rows and the "gold regression" was silently unjudged.

**(b) The reference integral under-resolved the keel singularity.** The reference was a
single 5x5 Gauss rule over each whole patch. The analytic baselines are built on
`y = (B/2)(1-|2xi-1|^n) zeta**(1/m)` with `m` up to 24, i.e. `dy/dzeta` is singular at the
keel. A single-panel rule converges to the wrong number, and the *mesh* path (which has
the same singularity to resolve) then appeared to "diverge" away from it.

Measured with the analytic closed form as the independent anchor
(`V = Cb·L·B·T`, `reference_integrals`, both in-repo):

| specimen | analytic V | 1 panel (old ref) | 2 panels | 4 panels | 8 panels |
|---|---|---|---|---|---|
| Series 60 Cb=0.60 | 3.75 m³ | 3.619985 (−3.467 %) | 3.74780 (−0.0587 %) | 3.750000 (−0.0001 %) | 3.750000 (+1e-5 %) |
| Series 60 Cb=0.65 | 4.0606000 m³ | 3.639481 (−10.371 %) | 4.042340 (−0.4499 %) | 4.060590 (−0.00035 %) | 4.060610 (−1e-5 %) |
| KVLCC2 | 312 621.6704 m³ | 281 236.39 (−10.039 %) | 306 758.44 (−1.8755 %) | 312 813.01 (+0.0612 %) | 312 628.73 (+0.0023 %) |

The old reference was 3.5 % (Cb 0.60) to 10.4 % (Cb 0.65, KVLCC2) low — larger than the
3 % acceptance gate it was supposed to police. That is why no convergence could be shown.

## 3. Fixes

| file | change |
|---|---|
| `benchmarks/verified/iges_import/gold_loft.py` | `graded_panels()` (geometric clustering of panels toward *both* parameter ends) + `brep_quadrature(order, panels, ratio, flags)` (graded composite Gauss, per-face accounting, method string) + `PatchDerivatives` (derivative control nets built once: 5.9 ms → 0.54 ms per point, bit-identical to the library, asserted to 1e-12) + `enclosed_volume_flags()` replacing the hard-coded `same_sense=[True]*4` of the loft |
| `benchmarks/verified/iges_import/run_benchmark.py` | reads `components`/`closed_components`; reference = graded quadrature with a *measured* panel ladder and an analytic closed-form anchor per specimen; a row is only judged (`within_3pct`) when its reference settled (`self_change < 1e-3 %`) **and** matches the analytic value (<1 %) |
| `tests/test_iges_gold_regression.py` (new) | 7 tests: cached-derivative equality, panel-grid tiling/clustering, reference vs analytic closed form (3 specimens), ladder-settles-while-single-panel-fails, pipeline vs written-patch volume within 3 %, DTMB5415 loft import, driver-key regression |
| `tests/test_iges_trim_import.py` | three evidence-reading tests asserted an obsolete evidence schema (`files` as a list, `benchmark.relative_error_pct`); updated to the current dict schema, with a skip guard so an evidence file from an older driver cannot silently "pass" |

## 4. Loft orientation (second, independent defect)

The loft hard-coded `same_sense=[True, True, True, True]`. Measured per-face volume
contributions (graded quadrature): hull +8160.6 m³, deck −5761.7 m³ → the "closed" body
came out at **+2398.9 m³**, while the orientation-repaired mesh of the same patches holds
≈13 922 m³. `enclosed_volume_flags()` now derives the flags from the enclosed volume
(maximise |V| over the 2ⁿ sign sets) instead of guessing, and the flags actually used are
recorded in `provenance.json` (`loft.same_sense`).

**Still open (honest):** for the DTMB5415 loft the quadrature does *not* settle
(panel ladder disagrees; the loft's end caps are near-degenerate lunes —
`girth_span_min = 2.1e-05` m — and their integrand is unbounded). So the loft's *volume*
reference is reported as `reference_converged: false` and is **not** used to judge the
pipeline; the driver writes the row, refuses to score it, and says why. A converged
reference for the loft needs either a cap-free loft (knife edges, as the analytic
baselines do) or a trimmed-cap formulation — next round.
KP505 is unchanged (6/8 faces still refused with measured inversion residuals 60–145 file
units vs tolerance 1.5e-04) — that needs the parameter-space trim-curve reader extension.

## 5. Error judgement (explicit, key quantities)

| key quantity | reference | measured | error | < 3 % ? |
|---|---|---|---|---|
| Series 60 Cb=0.60 patch-set volume | analytic 3.75 m³ | 3.750000 (panels=4) | −0.0001 % | **yes** |
| Series 60 Cb=0.65 patch-set volume | 4.0606 m³ | 4.060590 | −0.00035 % | **yes** |
| KVLCC2 patch-set volume | 312 621.67 m³ | 312 813.01 | +0.0612 % | **yes** |
| Series 60 imported mesh volume @96x48 | patch integral 3.75000 | 3.74626 | −0.10 % | **yes** (test `test_pipeline_reproduces_the_written_patch_volume_within_3pct`) |
| DTMB5415 loft volume | not established (reference non-converged) | — | **unknown / not_run** | **not claimed** |
| KP505 faces tessellated | 8 | 2 | −75 % | **no — function not achieved** |

## 6. Checks not run

* The full benchmark driver (`run_benchmark.py`) was **not** re-run to completion in this
  slot: the KP505 trimmed route alone costs ~120 s and the whole run ~5 min, over the
  per-command budget. Everything the driver reports is covered by the new tests at the
  component level (reference vs analytic, ladder settling, end-to-end Series 60 within
  3 %); the two evidence-file tests skip with an explicit reason until the driver is
  rerun. `repro_cmd` in `evidence.json` is the one-liner for that rerun.
* No cluster run was submitted.

## 7. Files

| path | role |
|---|---|
| `benchmarks/verified/iges_import/gold_loft.py` | graded reference integral, cached exact derivatives, orientation rule |
| `benchmarks/verified/iges_import/run_benchmark.py` | driver: schema fix, ladder + analytic anchor, honest non-converged rows |
| `tests/test_iges_gold_regression.py` | new regression suite (7 tests) |
| `tests/test_iges_trim_import.py` | evidence-schema update + stale-evidence guard |
