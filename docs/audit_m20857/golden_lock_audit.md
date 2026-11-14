# Anchor GAP-CAD-D2-GOLDEN-REGRESSION: the committed baseline did not lock three of the four golden meshes

**Mission**: m20857 (anchor `GAP-CAD-D2-GOLDEN-REGRESSION`, TensorCAD, branch `agent/m20857`)
**Audited tree**: `c772985` (merge #501) + this mission's fix
**Disposition**: request-changes → fixed in this mission (baseline re-blessed, severity raised, CI step added)

## 1. What the anchor claims, and what the tree actually did

The suite (`tensorcad/golden.py`, `.github/workflows/golden-geometry.yml`,
`.githooks/pre-commit`, `docs/golden_regression.md`) is the golden library required by the
anchor: four cases (`wigley`, `series60`, `sphere_ellipsoid`, `dtmb5415`), every reference
recomputed from a source named in `reference_source` (no hand-copied criteria), tolerances
derived from the two-grid change.

`docs/golden_regression.md` §"the mesh itself, not only the numbers integrated on it"
states that the baseline was re-blessed and that *all four* cases fingerprint their mesh.
The committed baseline contradicted it:

```
$ git show HEAD:docs/golden/tensorcad_geometry_golden.json | python -c "..."   # c772985
{'dtmb5415': 'EMPTY', 'series60': 'EMPTY', 'sphere_ellipsoid': 'EMPTY', 'wigley': '0ef07ca6'}
```

`geometry_sha256` was `""` for `series60`, `sphere_ellipsoid` and `dtmb5415`. Every master
revision of the file that `git log --follow` reports (`033bb62`, `d40e193`, `c34fc00`,
`9bb0560`) carries the same three empty strings, i.e. the documented re-bless never reached
`main`.

Consequences, in order of severity:

1. **The gate passed while 3 of 4 meshes were unlocked.** `geometry_fingerprint_not_blessed`
   was a `WARNING`, and `SuiteResult.passed` only looks at `ERROR` findings, so
   `python -m tensorcad.golden --check` exited 0 and the CI job's `assert data["passed"]`
   was satisfied. A re-tessellation, a vertex renumbering or a winding flip on the
   Series-60, SDF or 5415 path is invisible to every numeric metric (the area takes the
   norm, the sphere/ellipsoid volume the absolute value) — `tests/test_golden_regression.py::
   test_vertex_relabelling_is_caught_only_by_the_fingerprint` demonstrates exactly this —
   and the fingerprint that would have caught it was not compared.
2. **The suite's own regression test was red on HEAD.** `pytest tests/test_golden_regression.py`
   → `test_every_golden_case_fingerprints_its_mesh` FAILED:
   `assert 'f9da70a4…' == ''` (series60). The CI workflow never ran that file (the
   `golden` job runs the gate script only, `mutation-audit` runs
   `tests/test_golden_mutation.py`), so a red acceptance test never blocked a merge.
3. The failure mode is silent by construction: nothing in the pipeline distinguishes "this
   case locks its geometry" from "this case locks only its integrals".

## 2. Root cause

Two changes landed out of order: the per-case `geometry_sha256` implementations for
`series60` / `sphere_ellipsoid` / `dtmb5415` (and the read-back in `evaluate_case`) landed
without a re-blessed baseline reaching `main`; the missing-fingerprint finding was
classified `WARNING`, so the gate kept saying yes. A warning that means "this case is not
guarded" is not a warning — it is a gate that is not a gate.

## 3. Fix

| # | change | file |
|---|---|---|
| 1 | re-bless: the four `geometry_sha256` fields are filled (`git diff` = 3 fingerprints + 2 provenance fields, **no measurement, reference or tolerance moved**) | `docs/golden/tensorcad_geometry_golden.json` |
| 2 | `geometry_fingerprint_not_blessed` upgraded `WARNING` → `ERROR` (fail closed): a baseline that does not lock the mesh turns the gate red until someone runs `--bless` deliberately | `tensorcad/golden.py` |
| 3 | the unit test now asserts the error, not a pass-with-warning | `tests/test_golden_regression.py::test_baseline_without_fingerprint_is_an_error` |
| 4 | the `golden` CI job runs `tests/test_golden_regression.py`, so a baseline that loses its fingerprints cannot merge green | `.github/workflows/golden-geometry.yml` |
| 5 | docstrings/`docs/golden_regression.md` corrected to the enforced behaviour | `tensorcad/golden.py`, `docs/golden_regression.md` |

### The re-bless moved no number

```text
$ git diff --stat docs/golden/
 docs/golden/tensorcad_geometry_golden.json | 10 +++++-----
$ git diff docs/golden/ | grep '^[-+] ' | grep -v geometry_sha256 | grep -v -E 'git_commit|geometry_fingerprint'
(no output)
```

`holdout_rms_m` (dtmb5415) reproduced to `-2.65e-08` relative, i.e. five orders below the
`1e-9` drift lock, and the analytic cases reproduce bit-identically; the re-bless is
therefore a pure lock-strengthening, not a re-baseline of numbers (a re-bless that *had*
moved numbers would have had to justify each move).

### Verification

```text
$ python -m tensorcad.golden --check            # full two-grid, 22.7 s
passed=True, findings=[]                        # was: 3 warnings, exited 0 anyway
wigley 0ef07ca68299  series60 f9da70a467f6
sphere_ellipsoid f6ded09395f7  dtmb5415 30c05a782971
$ pytest tests/test_golden_regression.py tests/test_golden_mutation.py -q
# see the mission report for the exact counts (cluster job)
```

## 4. Negative control: does the raised severity actually bite?

`tests/test_golden_regression.py::test_baseline_without_fingerprint_is_an_error` empties one
fingerprint in a blessed temporary baseline and requires `result.passed is False` with
`geometry_fingerprint_not_blessed` among the errors. On the pre-fix engine the same input
passed with a warning — the test failed on the unmodified baseline and passes after the
fix, which is the invariant the fix targets.

## 5. Not run / residual

* The mutation-audit record (`docs/golden/mutation_audit.{json,md}`) was not regenerated:
  its findings are filtered to `ERROR` codes, and the control run (unmutated library vs the
  now-fingerprinted baseline) still produces no error, so the record is expected to be
  unchanged. `tests/test_golden_mutation.py` in the cluster job is the check that this
  holds; if it disagrees, the audit must be re-run rather than hand-edited.
* The committed audit still names its baseline path from another worktree (`m20561`); a
  cosmetic field, not touched here.
* Declared blind spot `series60/maximum_section_area_m2` (≥ 1 % blind per the audit
  envelope) is unchanged by this fix: the fingerprint locks the mesh, not the metric's
  sensitivity.
* The git hook still runs `--quick` (coarse only); with the severity change an un-blessed
  baseline now blocks a commit locally as well, but the hook only covers checkouts that ran
  `tools/golden/install_git_hooks.sh`. CI is the enforced layer.