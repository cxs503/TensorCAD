# Golden gate reach: a geometry change that arrives as a data edit

Anchor: `GAP-CAD-D2-GOLDEN-REGRESSION` (golden library + "any tensorcad geometry change
runs the gate"). This note records one hole found in the *reach* half of that acceptance
criterion, and the fix. Everything below is reproducible with the commands quoted.

## 1. The defect

The gate measured four cases; one of them (`dtmb5415`) is not generated from code alone.
It reads a shipped reduced model:

| file | role in the case |
| --- | --- |
| `docs/data/dtmb5415/design_space.npz` | the rank-40 basis the hull is rebuilt from (`DTMB5415Space`) |
| `docs/data/dtmb5415/design_space_metrics.json` | the fit that produced that basis |
| `docs/data/dtmb5415/holdout_fixture.npz` | **the reference source**: `per_config_rms_m[0]` is the `holdout_rms_m` reference |

None of them is under `tensorcad/`. The trigger set was:

* `staged_requires_golden()` — `tensorcad/**.py`, the suite's own test, the baseline,
  `docs/golden/**`;
* `geometry_files()` (the source fingerprint) — `tensorcad/**.py` + the suite's test;
* `.github/workflows/golden-geometry.yml` `paths:` — the same `tensorcad/**` + suite/tools.

So editing any of the three datasets changed the geometry being asserted **and** the
reference it is asserted against, while both the local hook and CI stayed quiet. Measured
on the pre-fix tree:

```
$ printf 'docs/data/dtmb5415/holdout_fixture.npz\n' > /tmp/s.txt
$ python -m tensorcad.golden --staged /tmp/s.txt --quick
golden: no staged geometry change, check skipped      # <- green pass over a changed dataset
```

That is the precise failure mode the golden suite exists to prevent: not a wrong number,
but a gate that was never asked.

## 2. The fix — declare the inputs, derive the reach

`GoldenCase.geometry_inputs()` is new: a case states the data files it reads. `dtmb5415`
returns the three files above; every code-only case returns nothing (the `tensorcad/` rule
already covers it).

Derived from that one declaration:

* `staged_requires_golden()` triggers on a declared input (the hook therefore runs);
* `geometry_files()` includes them, so `geometry_fingerprint()` is now a statement about
  the *inputs* of the measurement, not only about the code that reads them — a dataset
  edit moves neither a file name under `tensorcad/` nor a module, so a code-only
  fingerprint called it unchanged;
* the CI `paths:` filter names `docs/data/**`. A test (`_glob_matches`) holds that filter
  to the declared inputs, so the next data-driven case cannot be added without covering
  its paths.

Deliberately *not* done: triggering on everything under `docs/data/`. The trigger is the
declared set; a data file no case reads still skips the gate
(`test_staged_data_that_no_case_reads_is_still_skipped`).

## 3. A red self-test on `main`, found on the way

`tests/test_golden_regression.py::test_a_declared_tolerance_must_be_justified` was
**failing on the unmodified tree** (`git stash` baseline reproduced it), which makes the
golden CI job red for every geometry change — a gate that is always red is not a gate.

Root cause (test, not engine): the case declared the band-without-a-basis on
`coarse_only`, a key the probe's `reference()` does not assert, so the engine's
`declared_tolerance_unknown_metric` branch pre-empted `declared_tolerance_without_basis`.
The three ways a band can be worthless cannot all be expressed in one mapping when the
probe asserts one reference metric. The test now uses two declarations/two evaluations,
and the fail-closed claim is asserted through `run_suite` rather than on a bare
`CaseResult` — the old `assert result.passed is False` was vacuous, because
`evaluate_case` hands declaration errors back to the caller and the bare case result had
no accuracy error.

## 4. Verification

```
# engine + suite (submitted to the cluster queue; the file takes > 2 min)
python -m pytest tests/test_golden_regression.py tests/test_golden_mutation.py -q
# the gate itself, unchanged in outcome
python -m tensorcad.golden --check          # PASS, four cases
```

Regression proof (new tests against the pre-fix tree, code+workflow stashed):

```
test_every_declared_input_exists ......................... FAILED
test_staged_declared_input_triggers_the_gate ............. FAILED
test_geometry_files_include_the_declared_inputs .......... FAILED
test_source_fingerprint_moves_when_a_dataset_edit_lands .. FAILED
test_git_hook_runs_the_gate_on_a_dataset_edit ............ FAILED
test_ci_path_filter_watches_every_declared_input ......... FAILED
test_staged_data_that_no_case_reads_is_still_skipped ..... passed   (a guard, unchanged)
```

## 5. Limits

* The `data/` copy of the hull referenced by `tensorcad.dtmb5415.REFERENCE` lives outside
  the repository (absolute path, pinned `sha256`); it cannot be in a staged-path set, so
  it is pinned by hash only. A missing/changed copy is a runtime refusal of the acceptance
  script, not a gate trigger.
* Reach is still path-based: a dataset pulled from the network at run time is invisible.
* Truthfulness of the *numbers* is unchanged by this note; it is about whether the gate is
  asked at all.