# m5479 geometry-audit gates artifact: restore + persistence note

**Mission**: m16119 (selfstart generalist seat)
**Date**: 2026-11-12
**Branch**: agent/m16119

## 1. Case

`tests/test_geometry_audit_m5479.py::TestC5GatesArtifact::test_gates_file_exists_and_valid`
was RED: it requires the committed evidence file
`artifacts/m5479_geometry_audit/gates_verified.json`, which was absent from the
worktree. The gap was first recorded by m6410 (2026-11-08) and re-flagged by the
exchange-verification report (id 15201, 2026-11-12): the file "从未入库" and the
merge gate "反复滤除" it, so the red test recurred across missions.

## 2. What was done

1. Ran the **committed** generator `tools/generate_gates_m5479.py` (no code change).
   It re-runs the `tensorcad.geometry_audit` pipeline and writes the machine-readable
   gate table.
2. Committed **only** the 2.5 KB evidence file `gates_verified.json`. The generator's
   462 KB byproduct `wigley_test.patchset.json` is a runtime byproduct (the
   `TestFullChainValidation` round-trip test regenerates it into `ARTIFACTS_DIR`), so it
   is deliberately **not** committed.
3. Re-ran the suite: `tests/test_geometry_audit_m5479.py` → **30 passed** (was 29 passed,
   1 failed).

## 3. Gate results (read from the regenerated artifact, zero hand-copy)

| Gate | Criterion | Measured | Threshold | Passed |
|------|-----------|----------|-----------|--------|
| C1 | seam-gap / flip / hole / duplicate recall | recall 1.000 (24/24) | >= 0.95 | yes |
| C2 | G2 curvature-break detection | bump detected, small bump not flagged | — | yes |
| C3 | offsets-table vs mesh hydrostatics | worst rel diff 7.8e-4 | < 1 % | yes |
| C4 | Wigley vs analytic closed forms | worst rel err 7.8e-4 | < 2 % | yes |
| C5 | gates_verified.json produced & valid | 5/5 gates, all ids present | — | yes |

All five gates `passed=true`, `all_passed=true`.

## 4. Error judgement

| Quantity | Reference | Measured | <3 %? |
|----------|-----------|----------|-------|
| C3 worst relative difference (vol/area/Cb) | independent Simpson route | 7.8e-4 (0.078 %) | **yes** |
| C4 worst relative error (Wigley hydrostatics) | analytic closed forms | 7.8e-4 (0.078 %) | **yes** |

Both key quantities are far below the 3 % line.

## 5. Caveats / open item

- The **merge gate keeps filtering this artifact** out of merged branches (see the
  "N 个文件被门禁滤除" subjects in recent merge commits). This mission restores it in the
  worktree and commits it; if the gate filters it again on merge, that is a platform
  gate-policy defect (evidence files under `artifacts/` should not be dropped) and should
  be raised as a platform issue, not papered over by weakening the test.
- The test was **not** weakened: it still asserts the file exists, the schema, the five
  gate ids, and `all_passed`.
