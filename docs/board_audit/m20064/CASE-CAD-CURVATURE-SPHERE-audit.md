# Board audit `CASE-CAD-CURVATURE-SPHERE` (mission 20002) -- closure verification

* audited commit `4558ef6d3d18e91090d777ab05987528a85fe72e` (benchmarks: verified case CASE-CAD-CURVATURE-SPHERE (sphere discrete mean curvature vs analytic 2/R))
* audited work tree `/root/tp-agent/work/m20002`, clean tree (no uncommitted change)
* auditor work tree `/root/tp-agent/work/m20064` (mission 20064), branch `agent/m20064`
* verdict: **L3, gate reopened** -- see §7

## 1. Contract

The mission claims a verified benchmark case: *the discrete mean curvature of a sphere,
`S = kappa_1 + kappa_2 = 2H`, against the analytic `2/R`, 2 % criterion*, delivered as
runner + `result.json` + gates + figures + a six-section document, reproducible by one
command. The board audit re-derives every one of those claims instead of reading the
mission's own report:

| claim | how the auditor checked it |
|---|---|
| the artifacts exist | sha256 + byte count of every committed file in the audited tree (§2) |
| the numbers are in tolerance | the ledgers are queried, then every gate is recomputed from `reference`/`measured` (§3) |
| the run reproduces | the committed runner is re-executed in a clean process; 590 numeric leaves are diffed against the committed `result.json` (§4) |
| the tests pin the claims | `pytest tests/test_case_cad_curvature_sphere.py` re-run (§5) |
| the number is a measurement, not a restatement of `2/R` | three anti-tautology probes that perturb the geometry / the reference (§6) |

## 2. Deliverables -- do the three pieces exist?

* `benchmarks/verified/sphere_mean_curvature/run_benchmark.py` -- True (34989 B, sha256 `fbc4ddc4ec4c0092...`)
* `benchmarks/verified/sphere_mean_curvature/result.json` -- True (31984 B, sha256 `3f19e9d4bf53193c...`)
* `benchmarks/verified/sphere_mean_curvature/gates_verified.json` -- True (6651 B, sha256 `6da7da79e9967415...`)
* `benchmarks/verified/sphere_mean_curvature/README.md` -- True (2958 B, sha256 `abdcb9811d8be759...`)
* `docs/sphere_mean_curvature.md` -- True (11739 B, sha256 `1dd3164c09c236bb...`)
* `tests/test_case_cad_curvature_sphere.py` -- True (4558 B, sha256 `ebb7186c2bd8ca7b...`)
* `benchmarks/verified/sphere_mean_curvature/figures/convergence.png` -- True (96726 B, sha256 `8254f56183ca3c8a...`)
* `benchmarks/verified/sphere_mean_curvature/figures/error_cloud.png` -- True (483422 B, sha256 `cc5310e963d2b95e...`)
* `benchmarks/verified/sphere_mean_curvature/figures/readings.png` -- True (39870 B, sha256 `4525849f3b59803d...`)

The document `docs/sphere_mean_curvature.md` carries the required six sections:
"1. Question and reference", "2. Method", "3. Parameters (all dimensionless groups)", "4. Reference vs measured", "5. Convergence and sensitivity", "6. Failure modes and artifacts".

## 3. The benchmark judge: what does the ledger say?

| query | rows |
|---|---|
| `tensor_benchmark_cases` where `anchor_slug = 'CASE-CAD-CURVATURE-SPHERE'` | **0** |
| `tensor_benchmark_cases` where `anchor_slug ilike '%CAD%'` | 9 |
| `tensor_benchmark_cases` total | 121 |
| `tensorbench_results` where `metric ilike '%curvature%'` | 0 |

Nine other TensorCAD anchors *are* registered (e.g. `GAP-CAD-D2-GOLDEN-SPHERE_ELLIPSOID`,
`GAP-CAD-T1-BASELINE-FLEET-KCS-DISPLACEMENT`), so registration is the established route and
this anchor simply did not take it. **The sole judge has no row to rule on**: the case's
`ref`/`tolerance` pair exists only inside its own `gates_verified.json`.

Re-derived gate table (auditor computes `|got-ref|/|ref|` or, for `G2` where the reference
is "zero error", the absolute value, and compares to the gate's own tolerance):

| gate | criterion | ref | got | error | tol | mission | audit |
|---|---|---|---|---|---|---|---|
| G1 | discrete mean curvature, area-weighted integral, 1 m sphere, marching  | 2 | 2.000237787 | +0.01189 % | 0.02 | pass | re-checked |
| G2 | pointwise discrete mean curvature, median relative error, structured U | 0 | 0.003938717571 | — | 0.02 | pass | re-checked |
| G3 | observed order of the pointwise estimate on the structured sphere | 2 | 1.946204219 | -2.69 % | 0.5 | pass | re-checked |
| G4 | field-path kappa_1+kappa_2, sphere R=0.1 m, 40 Fibonacci probes | 20 | 19.99999952 | -2.396e-06 % | 1e-06 | pass | re-checked |
| G4 | field-path kappa_1+kappa_2, sphere R=0.25 m, 40 Fibonacci probes | 8 | 7.999999808 | -2.396e-06 % | 1e-06 | pass | re-checked |
| G4 | field-path kappa_1+kappa_2, sphere R=0.5 m, 40 Fibonacci probes | 4 | 3.999999904 | -2.396e-06 % | 1e-06 | pass | re-checked |
| G4 | field-path kappa_1+kappa_2, sphere R=1 m, 40 Fibonacci probes | 2 | 1.999999952 | -2.396e-06 % | 1e-06 | pass | re-checked |
| G4 | field-path kappa_1+kappa_2, sphere R=10 m, 40 Fibonacci probes | 0.2 | 0.1999999952 | -2.396e-06 % | 1e-06 | pass | re-checked |
| G4 | field-path kappa_1+kappa_2, sphere R=2 m, 40 Fibonacci probes | 1 | 0.999999976 | -2.396e-06 % | 1e-06 | pass | re-checked |
| G4 | field-path kappa_1+kappa_2, sphere R=5 m, 40 Fibonacci probes | 0.4 | 0.3999999904 | -2.396e-06 % | 1e-06 | pass | re-checked |
| G5 | surface integral of kappa_1+kappa_2 over the mesh, level 48 | 25.132741 | 25.1242619 | -0.03374 % | 0.02 | pass | re-checked |
| G6 | discrete mean curvature, area-weighted integral, sphere R=0.1 m, level | 20 | 20.00979779 | +0.04899 % | 0.02 | pass | re-checked |
| G6 | discrete mean curvature, area-weighted integral, sphere R=0.5 m, level | 4 | 4.001959557 | +0.04899 % | 0.02 | pass | re-checked |
| G6 | discrete mean curvature, area-weighted integral, sphere R=1 m, level 2 | 2 | 2.000979779 | +0.04899 % | 0.02 | pass | re-checked |
| G6 | discrete mean curvature, area-weighted integral, sphere R=10 m, level  | 0.2 | 0.2000979779 | +0.04899 % | 0.02 | pass | re-checked |
| G6 | discrete mean curvature, area-weighted integral, sphere R=2 m, level 2 | 1 | 1.000489889 | +0.04899 % | 0.02 | pass | re-checked |

All 16 numeric gates pass under the auditor's own arithmetic; the one
`diagnostic` row (`G7`, raw pointwise reading on the marching-tets mesh) is a **fail** and
is reported as such by the mission -- it is not gated away.

## 4. Reproduction

```
$ cd /root/tp-agent/work/m20002 && python benchmarks/verified/sphere_mean_curvature/run_benchmark.py \
      --out /root/tp-agent/work/m20064/runs/m20064_sphere_curvature
case CASE-CAD-CURVATURE-SPHERE: 16 gates, 1 diagnostics
headline S vs 2/R: measured 2.000238 vs 2.000000 (+0.0119 %)  -> within 2 %
pointwise median: raw 50.12 %, conditioned 1.94 %
figures: ['figures/convergence.png', 'figures/readings.png', 'figures/error_cloud.png']
seconds: 46.7
```

**590 / 590 numeric leaves identical**
to the committed `result.json` (only wall-clock, environment and figure path fields excluded);
mismatches: 0. Fresh run: 46.7 s.
The `repro_cmd` in the committed artifact is therefore executable and deterministic.

## 5. Tests

`8 passed, 1 warning in 15.50s` -- 8 tests, including
test_the_analytic_reference_is_the_sum_of_the_principal_curvatures, test_the_integral_reading_is_inside_the_two_percent_criterion, test_the_pointwise_reading_needs_a_well_shaped_triangulation, test_the_pointwise_estimate_keeps_second_order_on_a_structured_sphere. They assert invariants (tolerance membership,
observed order, radius independence, and that the raw pointwise reading is *outside* 2 %),
not merely that code runs.

## 6. Is the headline a measurement? (anti-tautology probes)

| probe | result |
|---|---|
| `+-1 %` radial dent, same topology, same reference | integral moves from +0.0001 % to -0.0179 % -- response 1.797e-04 |
| mesh scaled x2, reference kept at `R=1` | measured integral = 1.000000706 = 0.500000353 x reference (expected 0.5) -> measured side does not depend on the reference |
| same sphere, reference at `R=2` | measured integral unchanged to 0.0e+00 1/m |
| **which reading discriminates?** | pointwise median error +0.0999 % -> +19.8627 % under the dent; the gated integral responds by only 1.80e-04 |

The measurement is genuinely mesh-driven (`mesh-driven`,
`reference-independent`). But probe 4 is the audit's central finding: **the gated reading is
the shape-insensitive one.** For a closed surface `int S dA` is nearly invariant, so a 2 %
gate on it cannot detect a 1 % shape error, while the pointwise reading -- the quantity a
fairing/bending-energy consumer actually uses -- moves by a factor of ~200 on the same
perturbation and sits at **+50.1 % raw** on the
marching-tets mesh (`G7`), i.e. far outside the platform's 15 % precision floor.

## 7. Findings, severity-ranked

| # | severity | finding | evidence |
|---|---|---|---|
| F1 | blocking (L4) | `tensor_benchmark_cases` has **0 rows** for `CASE-CAD-CURVATURE-SPHERE`: no `ref_source`/`metric`/`tolerance` in the ledger, so the benchmark judge cannot confirm L4 and no downstream mission can cite the criterion | §3; `evidence.json:ledger` |
| F2 | material | the gated metric is shape-insensitive; the sensitive pointwise reading is +50.1 % raw on unstructured geometry and only +1.94 % after `taubin_smooth(10) + project_to_field(2)` (the latter needs the implicit field, unavailable for a generic CAD B-rep) | §6; `mesh_dependence.json:which_reading_discriminates`; result `marching.series` |
| F3 | material (L3 caveat) | delivery commit is not on `origin/main`: `git branch --contains 4558ef6d` = ['agent/m20002']; no platform-side caller (`grep` over `TensorPlatform/backend`, `/api` = 0 hits); repo-internal mainline consumers only: `tensorcad/fairness_energy.py:170`, `:486`, `examples/fairness/run_example.py:174` | `evidence.json:integration_callers_*` |
| F4 | minor | unit convention mixing inside `gates_verified.json`: `error_percent` is in percent while `tolerance` is a fraction, and `G4`'s `1e-6` is relative while `G3`'s `0.5` is absolute; a machine reader comparing the two fields directly is off by 100x | `gates_verified.json` (see §3 sign/scale) |
| F5 | minor | the `conditioned` reading smooths the mesh by `taubin_drift = 0.0261` m at level 48 (2.6 % of R) before it reaches 1.94 %: the post-conditioning figure is a different geometry, not the CAD-like mesh | `result.json:marching.series[-1].conditioned` |
| S1 | strength | reproduction is exact (590/590 leaves), 8 passed, 1 warning in 15.50s, the failing pointwise reading is disclosed in README, doc and gates rather than gated away | §4, §5 |

## 8. Level and disposition

* **L1** (math/unit verification): yes -- analytic `2/R`, closed form `int S dA = 8 pi R`.
* **L2** (real code, input -> compute -> output): yes --
  `tensorcad/curvature.py:292 vertex_mean_curvature` (cotangent Laplacian,
  `:359` mixed-area/2 normalisation), radius-free `measure()` at
  `benchmarks/verified/sphere_mean_curvature/run_benchmark.py:178`,
  entry `run_benchmark.py:748 main` / `repro_cmd`.
* **L3** (integration): yes, with the caveat of F3 -- the audited estimator is on the
  library's mainline fairing-energy path (`tensorcad/fairness_energy.py:170`,
  `:486`, public export `tensorcad/__init__.py:751`); no platform-side caller.
* **L4** (benchmark + error analysis + reproducible): **withheld** -- the case work itself is
  L4-grade (16 gates, two convergence series, radius sweep, error analysis, exact
  reproduction), but F1 leaves the sole judge without a row, and F2 shows the gated metric
  is the weak one. Per the cumulative scheme L4 is not confirmable until the anchor is
  registered.
* **L5**: no -- no new research result; the estimator and its second-order behaviour are
  textbook (Meyer et al. 2003).

**Disposition: request changes** (gate `reopened`). Remediation, in order:
1. register the anchor in `tensor_benchmark_cases` with `ref_source` = do Carmo Sec. 3-3
   Example 2 + Meyer et al. Sec. 3.3, `metric` = "area-weighted integral of S=2H vs 2/R",
   `tolerance` = 2 %, and state in the same row that the pointwise reading is diagnostic
   only -- or gate the pointwise reading with a predeclared criterion once it is fixed;
2. either fix the pointwise estimator on unstructured meshes (mixed Voronoi areas,
   angle-based weights, or spline/normal-cycle reconstruction) to reach < 15 %, or narrow
   the case claim in `README`/doc to the integral reading;
3. merge `agent/m20002` to `main` so the L3 integration is real on the audited branch.
