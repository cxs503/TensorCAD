# The recorded criterion is a derivation, not a chosen number

Scope: `GAP-CAD-D2-GOLDEN-REGRESSION` (golden geometry regression suite + CI). This note
records two defects found by re-reviewing the gate that this mission ships, and the minimal
fixes with their evidence. Both are about the gate being *able to fail*, not about the hulls
it measures.

## D1 — in `--quick` mode the acceptance band was an input a hand could edit

`tensorcad/golden.py` derives every acceptance tolerance as twice the observed coarse-to-fine
change (a two-grid bound) and stores it in `docs/golden/tensorcad_geometry_golden.json`. That
is the right discipline, and CI (`--check`, full) re-derives it on every run.

`--quick` — the mode `.githooks/pre-commit` runs on every local commit — has no coarse level, so
it asserted the acceptance **against the stored value**, and nothing verified that the stored
value was the derivation. Reproduced on the committed baseline: multiplying the Wigley
block-coefficient band by 100 by hand made `python -m tensorcad.golden --check --quick` print
`PASS` with `tol = 4.33e-01`, i.e. the criterion in force was 43 % on a quantity whose honest
bound is 0.43 %. Every measurement, reference and mesh fingerprint stayed honest; only the
criterion moved. The coverage table still read green.

Fix (`tensorcad/golden.py`): `_check_recorded_criterion` re-derives each recorded band from the
inputs it was derived on — the blessed coarse/fine measurements for a two-grid bound, the case's
own declaration for a declared band — and refuses any value that does not reproduce
(`tolerance_not_derived`, relative 1e-9), in quick **and** full mode. Three adjacent escape
routes are closed with it: a band recorded for a metric no reference asserts
(`tolerance_unknown_metric`), a band whose recorded origin contradicts the case's
(`tolerance_origin_mismatch` — quick and full would otherwise assert different criteria on the
same tree), and a band without the measurements to derive it from. `DECLARED_BASIS_PREFIX` now
names the single marker that distinguishes a declared band from a two-grid one.

## D2 — the negative control did not load, so the CI step was red on `main`

`tests/test_golden_regression.py::_load_redteam` loaded `tools/golden/mutation_redteam.py` with
`importlib.util.module_from_spec` but never registered the module in `sys.modules`. Dataclasses
resolve their own module through `sys.modules` while the class body runs, so the load raised
`AttributeError: 'NoneType' object has no attribute '__dict__'` and
`test_redteam_harness_kills_an_identity_mutant` failed — on `HEAD`, before this mission's change
(verified by running the committed test file). The CI job "The suite's own tests"
(`.github/workflows/golden-geometry.yml`) therefore did not certify the gate. The working pattern
already existed two files away in `tests/test_golden_mutation.py::_load_audit`; the fix is that
one line plus its reason.

## Evidence

| evidence | file | how it is produced |
| --- | --- | --- |
| before/after, bound to commits and file hashes | `docs/golden/criterion_provenance.json` | `python tools/golden/criterion_provenance.py --before-rev 740c5f80 --json … --png …` |
| figure: recorded vs derived band, verdict matrix | `docs/golden/criterion_provenance.png` | same command |
| mutation red team, 8/8 killed | `docs/golden/redteam.md`, `docs/golden/redteam_quick.json` | `python tools/golden/mutation_redteam.py --json … --markdown …` |
| regression tests | `tests/test_golden_regression.py` | `python -m pytest tests/test_golden_regression.py -q` |

`criterion_provenance.json`: the *same* hand edit (Wigley Cb band ×100) is run against the
pre-fix tree (`before`, the whole gate-relevant tree at `740c5f80`, the parent of the fix
merge; `golden.py` sha256 `64ed9648…`) and against the fixed worktree (`after`, sha256
`88e0742b…`). Result: `before` quick gate **PASS** (the hole), `after` quick and full gate
**FAIL** with `tolerance_not_derived`, and both trees pass their honest-baseline control —
so the difference is attributable to the edit, not to a broken harness.

The `before` tree is a faithful snapshot, not a `golden.py` swap: the blessed baseline has
since grown a fifth case (`inclined_box`), and the pre-fix code fails closed on a baseline
it does not register (`baseline_case_unknown`), so a swap could no longer produce a passing
control.  Staging the baseline and the 5415 data from the same revision keeps the control
honest; the `after` tree additionally stages `tools/golden/inclined_box_oracle.py`, which the
inclined-box case loads by path.

Red team: the new mutant `baseline-tolerance-hand-inflated` (kind `criterion`) widens the blessed
Wigley Cb band 100× and no other character changes; it is killed by `tolerance_not_derived`, in a
run whose unmodified control passes.

## Notes and limits

- A change that moves the two-grid derivation now also trips `tolerance_not_derived` besides its
  own code (visible in the red-team table for the two `reference` mutants and the sphere-radius
  mutant). That is intended: the recorded band is stale until `--bless` records the new
  derivation, which is already the documented workflow for a reference change.
- The check binds the band to its derivation, not to the *quality* of that derivation: a
  measurement whose convergence is poor keeps a wide band honestly. The worst such case in the
  committed library remains `series60/maximum_section_area_m2` (declared tolerance 1.81e-02, the
  documented coherent blind spot of 2.72e-02 in `docs/golden/mutation_audit.md`). Tightening it
  needs a better-converging estimator for the maximum section area, not a tighter number.
- `docs/golden/mutation_audit.json/.md` (battery A/B of `tools/golden/mutation_audit.py`) was not
  re-run here: its battery A mutates library code and its policy is unchanged by this fix; the
  mutant class added here lives in the red team, which is the harness that sweeps the baseline.