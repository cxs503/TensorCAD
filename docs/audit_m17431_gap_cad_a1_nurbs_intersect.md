# Board audit — m17346 `GAP-CAD-A1-NURBS-INTERSECT` (TensorCAD)

Auditor mission **m17431** · verdict_by `board-audit` · subject commit **9f4697c** (branch
`agent/m17346`), merged as **#402 → 59f441e** · evidence snapshot with sha256 in
`artifacts/audit_m17431/evidence_m17346/` and `artifacts/audit_m17431/MANIFEST.json`.

This is a **cross-check (交差) audit**, not a redo: the subject's implementation was read,
its run artifacts re-executed, and its numbers re-derived by an independent route.

## 1. What was delivered (三件套 check)

| item | claim | audit result |
|---|---|---|
| implementation | `tensorcad/nurbs_intersect.py`, 1 176 lines, `nurbs_surface_intersection()` + `NURBSIntersectionResult` | **PRESENT, real** — marching squares on the oriented distance field `g(u,v) = (S1(u,v) − Q)·n2`, then a min-norm Gauss–Newton polish on the 4-parameter system `S1(u,v) = S2(s,t)` (3×4 Jacobian, Tikhonov damping), midpoint refinement, trim domains by the sign of `g`. No `NotImplementedError`, no placeholder body. |
| run artifacts | `artifacts/m17346/nurbs_intersect/{results.json, occt_xcheck.json, 3 × png}` | **PRESENT, non-empty** (results.json 101 669 B, figures 78–262 kB) |
| repro | `python examples/nurbs_intersect/run_example.py` (+ pytest) | **EXECUTABLE and deterministic** — see §3 |
| tests | 18 cases, `tests/test_nurbs_intersect.py` | **18 passed in 48.72 s** (re-run by auditor) |
| doc | `docs/nurbs_intersect.md`, rendered from results.json by the example script (no hand-copied number) | **PRESENT**, numbers regenerate |

## 2. Benchmark judgment (the only judge: `tensor_benchmark_cases`)

```sql
SELECT case_name, metric, tolerance FROM tensor_benchmark_cases
WHERE anchor_slug = 'GAP-CAD-A1-NURBS-INTERSECT';   -- 0 rows
```

**No benchmark case is registered for this anchor** (the table has 24 rows overall, none CAD
except `GAP-CAD-T1-BASELINE-FLEET`). The repo's own benchmark-record convention
(`benchmarks/verified/<case>/result.json` = committed ref-vs-got table, as used by
`parametric_hull` and `cad_highorder_geometry`) is also **absent** for this node. The
acceptance thresholds (radius ≤ 1e-8 R, lengths ≤ 3 %, cross-check ≤ 1e-4) live only in the
subject's own example script and doc — self-declared, not board-accepted.

Consequence for grading: **L4 (基准验证) cannot be certified by the board's own rule**, however
good the numbers are. And there is **no integration caller** (grep for
`nurbs_intersect` outside `tests/`, `examples/`, `tools/` returns nothing; `csg.py`,
`fillet.py`, `exchange.py` do not import it, and nothing in TensorPlatform calls it) →
**L3 (系统集成) not met either**. Verdict: **L2**.

## 3. Independent re-execution (reproducibility)

Environment: Python 3.11.6, loongarch64, clean copy of the subject tree at 9f4697c.

* `pytest tests/test_nurbs_intersect.py -q` → **18 passed, 48.72 s**.
* `python examples/nurbs_intersect/run_example.py --output artifacts/audit_ni` → exit 0,
  "all claims hold". The regenerated `results.json` was diffed leaf-by-leaf against the
  shipped one: **201 numeric leaves, 0 differences, 0 missing/extra keys** — the artifact is
  bit-reproducible, not a stale record.
* `tools/xcheck_occt_nurbs_intersect.py` → `occt_available: false`, `parity: WITHHELD`,
  exit 2. **Honest**: no invented OpenCASCADE numbers; the claim "parity with OCCT" is
  correctly *not* made.

Robustness nit (L0, non-blocking): `--doc` outside the repo raises
`ValueError: ... is not in the subpath of ...` (Path.relative_to in `main`) instead of
falling back to an absolute path. It bit the auditor on the first run; the subject's
documented invocation is unaffected.

## 4. Numbers re-derived by an independent route

Analytic closed forms recomputed here (not read from results.json):

| quantity | closed form | recomputed | subject value | rel. error |
|---|---|---|---|---|
| sphere R=2 × plane z=0.5 radius | `sqrt(R²−h²)` = 1.936491673104 | 1.936491673104 | 1.936491673104 | 5.6e-16 / R |
| same, length | `2πr` = 12.167336028 | 12.167336028 | 12.167258038 | 6.4e-06 |
| sphere 1 × 1 (d=1.2) radius, length | 0.800000, 5.026548246 | 0.800000, 5.026548246 | 0.800000, 5.026516027 | 4.4e-16 / R, 6.4e-06 |
| sphere 1 × 1.2 (d=1) radius, length | 0.960000, 6.031857895 | 0.960000, 6.031857895 | 0.960000, 6.031819232 | 5.6e-16 / R, 6.4e-06 |
| NURBS cylinder × 30° plane (ellipse) | `4aE(m)`, a=1/cos30, m=0.25 | 6.777917613 | 6.777877122 | 5.97e-06 |
| vertex residual on both patches | — | — | max **1.63e-15 m** at patch scale 8.5 m (1.9e-16 rel) | — |
| convergence order (cells 32→256) | O(h²) | 1.9992 / 1.9998 / 1.9999 | length errors 4.10e-04 → 6.41e-06 | order 2.00 confirmed |
| refinement (levels 0→3) | — | errors ×1/4 per level | 4.10e-04 → 6.41e-06 | order 2.00 |
| independent cross-check vs implicit-field marcher | different algorithm | 2.6e-10 / 2.8e-10 / 2.5e-10 | same | ✓ |
| trim area, plane disc (32 cells) | `πr²`=0.6597345 | 0.6597345 | 0.6594238 | 4.7e-04 |
| trim surface area, sphere band | 31.41592654 | — | 30.96347339 at (8,16) | 1.44e-02 (∈ 2 % test tolerance) |

Every headline number in `docs/nurbs_intersect.md` is backed by a leaf of `results.json`
that the auditor regenerated. **No overstated claim found in the report.**

## 5. Adversarial probes (`artifacts/audit_m17431/probe_audit.py`)

The subject's own tests are dominated by one side being a plane or a *rational* sphere; the
auditor therefore probed a generic **free-form bicubic clamped NURBS × tilted plane**, two
components, under-resolution, and scale invariance. References are obtained by *bisection of
the exact signed distance along dense v-isolines* (no marching squares, no Gauss–Newton, no
shared code).

| probe | result |
|---|---|
| free-form bicubic (6×6 net, deg 3) × tilted plane | 1 open curve; length **3.467680** vs independent reference **3.467622**, rel. err **1.66e-05**; every vertex on S2 to **6.9e-16** rel and on the plane to **5.5e-16** rel |
| same geometry at ×1e-3 and ×1e3 | identical relative errors (1.661974e-05 in all three) → tolerances are genuinely scale-relative, not tuned to one size |
| two-component cut (cubic patch crossing the plane twice) | **2 open curves**, lengths rel. err **2.3e-06** and **9.9e-06** vs reference; branch structure comes out of the cell topology as claimed |
| under-resolution (120×60 → 8×4 cells) | still **2 curves** at every level, errors degrade gracefully 9.9e-06 → 2.0e-03; **no silent branch loss observed** in this instance (this is one instance, not a robustness proof) |
| trim-domain polarity on the plane patch | 6/64 test points disagreed with the auditor's side classification; the auditor's classifier uses a coarse 30×30 nearest-point search, so this is **inconclusive, recorded as not-a-defect** (the subject's own analytic disc test asserts containment exactly and passes) |

The first probe run reported `curves = 0` and was traced to an **auditor error** (the probe
built the plane patch from a corner that was off the reference plane, so kernel and reference
solved different planes). Corrected; disclosed here because an audit that hides its own
false positive is worth less than one that does not.

## 6. Honest-failure coverage (checked, present)

Tangential contact returns no branch and raises `tangential_candidate` (`collapsed_branches`
counted); coincident patches raise instead of returning a cloud; a small finite patch
truncates the curve and reports it `open` (`closed_by == ""`); rank-deficient contact is
Tikhonov-damped and counted in `diagnostics["damped_steps"]`; `TrimLoop.area_reliable` is
False for the sphere latitude band that spans a seam (`loop_area_uv ≈ 1e-16` instead of a
meaningless shoelace area). These are the failure modes a B-rep kernel actually trips on, and
they are treated as first-class results, not hidden.

## 7. Verdict

| level | criterion | met |
|---|---|---|
| L1 | math/unit validation | yes — 18 tests, analytic anchors |
| L2 | real function, input → computation → output | **yes** — `nurbs_surface_intersection()`, 1 176 lines, re-executed, numbers re-derived |
| L3 | system integration (main solver / platform caller) | **no** — zero callers outside tests/examples |
| L4 | benchmark validation (complete case + error analysis + reproducible, registered) | **no** — no `tensor_benchmark_cases` row, no `benchmarks/verified/` record |

**Declared level: 2.** The numeric evidence (round-off vertex residuals, O(h²) convergence,
independent-method cross-check at 1e-10, scale invariance, bit-reproducible artifacts) is
plainly L4-*quality*; what is missing is the formal part — a registered benchmark row with a
board-owned tolerance and a committed `benchmarks/verified/nurbs_intersect/` record — plus any
caller that would make it L3. Remediation is therefore cheap and mechanical:

1. `INSERT` a `tensor_benchmark_cases` row: case = sphere×plane circle + NURBS cylinder×30°
   ellipse, metric = `radius_rel_error / length_rel_error`, tolerance = 1e-8 / 3 %.
2. Commit `benchmarks/verified/nurbs_intersect/` in the `parametric_hull` style
   (`run_benchmark.py` + ref-vs-got `result.json` + repro line).
3. Wire one caller (fillet or plate-boundary trimming in `csg.py`/`fillet.py`) to lift L3.

## 8. Checks not run

* OpenCASCADE parity — the OCP/pythonocc wheels are absent in this environment; the subject
  declares it WITHHELD and the auditor confirms the harness refuses to invent numbers.
* No stress test against degenerate multi-patch configurations (three patches meeting on a
  line, near-tangential grazing over a large region), no performance measurement, no
  robustness proof against pathological knot vectors — out of scope for this cross-check.