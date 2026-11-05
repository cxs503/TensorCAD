# PR #199 CI red, round 2 — root cause and cure (mission m4548)

## 1. The red: what CI actually reported

PR #199 (branch `agent/m4306`, tip `c6f9659`) is the *cure* of the dirty merge of
`agent/m3391`. Round 1 fixed two files and CI stayed red:

| check | conclusion |
|---|---|
| `test (3.13)` | cancelled |
| `test (3.12)` | cancelled |
| `test (3.11)` | **failure** |
| `test (3.10)` | cancelled |

`.github/workflows/tests.yml` runs the 3.10–3.13 matrix with the default
`fail-fast: true`. Three jobs *cancelled* and one *failure* is the signature of a
fail-fast matrix: one job dies, GitHub cancels its siblings. Which version is
reported as the failure is a scheduling accident — the defect is version
independent, so the diagnosis has to come from the failing step, not from 3.11.

## 2. Root cause

Replaying the job locally (CPython 3.11.6, the reported failure) the first gate
that dies is `ruff check .`, before `mypy` or `pytest` ever run:

```
$ ruff check .            # ruff 0.15.0, resolved from `dev = ["ruff>=0.6", ...]`
tools/make_report_m3859.py:21:34: RUF100 Unused `noqa` directive (unused: `E402`)
tools/make_report_m3859.py:22:21: RUF100 Unused `noqa` directive (unused: `E402`)
tools/make_report_m3859.py:40:15: RUF059 Unpacked variable `rc` is never used
tools/make_report_m3859.py:132:15: UP031 Use format specifiers instead of percent format
tools/make_report_m3859.py:163:68: W292 No newline at end of file
tools/verify_merge_m3859.py:124:9: F841 Local variable `lost` is assigned to but never used
tools/verify_merge_m3859.py:178:21: W292 No newline at end of file
Found 7 errors.
```

Three facts make the attribution unambiguous:

1. **PR #199 never touched those two files.** Its diff against `main` is three
   lines (`zip(..., strict=True)`, two missing trailing newlines) in
   `tests/test_merge_cure_m3391.py` and `tools/verify_merge_m3391.py`.
2. **`origin/main` fails the same gate.** `git show origin/main:tools/verify_merge_m3859.py`
   already trips `F841` and `W292`, so the red is not introduced by the branch.
3. **No code changed, the linter did.** `ruff>=0.6` is unbounded in the `dev`
   extra, so CI now installs ruff 0.15.x, which grew rule families that the tree
   never satisfied: `RUF100` (unused `noqa`), `RUF059` (unused unpacked
   variable), plus the pre-existing `F841`/`W292`/`UP031` that round 1 missed in
   the two *other* m3859 scripts. Round 1 fixed the pair of files the PR touched;
   the audit scripts the merge governor itself depends on stayed dirty.

So: a green tree turned red on a ruff release day with zero commits in between.
That is a CI reproducibility defect, and it is what this round fixes.

## 3. The cure (minimal, behaviour preserving)

| file:line | rule | before | after |
|---|---|---|---|
| `tools/make_report_m3859.py:21,22` | RUF100 | `import matplotlib.pyplot as plt  # noqa: E402` | directive dropped (ruff never raises E402 there — the `noqa` was dead) |
| `tools/make_report_m3859.py:40` | RUF059 | `name, rc, summary = line.split("\t")` | `name, _rc, summary = ...` |
| `tools/make_report_m3859.py:132` | UP031 | `label="... %.2f° ..." % box["limit_deg"]` | `label=f"... {box['limit_deg']:.2f}° ..."` (byte-identical output, checked) |
| `tools/make_report_m3859.py:163` | W292 | no trailing newline | trailing newline |
| `tools/verify_merge_m3859.py:124` | F841 | `lost: list[str] = []` (never appended, never read) | line removed; the live counter `losses` is untouched |
| `tools/verify_merge_m3859.py:178` | W292 | no trailing newline | trailing newline |
| `pyproject.toml` | — | `ruff>=0.6` | `ruff>=0.6,<0.16` with a comment: the lint gate must not move on its own |

No `tensorcad/` source file, no test, no workflow file is modified. The two
`tools/` scripts are the merge-audit tooling that produced the round-1 verify
evidence; their behaviour is unchanged (see §4.3).

## 4. Verification (measured, this worktree, CPython 3.11.6 / loongarch64)

### 4.1 The failing gate, before → after

| gate | before | after |
|---|---|---|
| `ruff check .` (ruff 0.15.0) | `Found 7 errors.` | `All checks passed!` |
| `mypy` (1.18.2, strict, 34 modules) | success (unchanged) | `Success: no issues found in 34 source files` |

### 4.2 The rest of the `test` job, replayed step by step

The default `test` job deliberately runs **without** the `neural` extra (no
torch). This box has torch 2.7.1 installed, so the run was made torch-free with a
`sitecustomize.py` shim that answers `find_spec("torch")` with
`ModuleNotFoundError` — the exact exception a fresh runner raises, which is what
`pytest.importorskip` needs in order to skip instead of error.

| step (`tests.yml`) | exit |
|---|---|
| `pytest -q` full suite | **893 passed, 17 skipped** in 404 s |
| `tensorcad-neural --help` | 0 |
| `tensorcad-hull` artifacts | 0 |
| reproducibility `cmp` (same params, same bytes) | 0 |
| `examples/geometry_core/run_benchmark.py --resolution 40 --points 5000` | 0 |
| `examples/parametric_hull/run_example.py` | 0 |
| `examples/suboff/run_example.py` | 0 |
| `examples/suboff_appendages/run_example.py --coarse 218,26,38 --fine 320,38,56` | 0 |
| `examples/mesh_conditioning/run_example.py --resolution 110,14,20 --targets 30000,20000` | 0 |
| `examples/step_export/run_example.py --resolution 60,8,11` | 0 |
| `examples/design_loop/run_example.py` | 0 |
| `examples/benchmarks/run_benchmark.py --quick` (exit code is the verdict) | 0 |
| artifact assertions from `tests.yml` (aff8, conditioning, step, design, `results.json` schema + required rows) | all ok |

`pytest --cov` itself could not be timed on this host: the coverage C tracer is
absent for loongarch64 (`CoverageWarning: Couldn't import C tracer`), so the
pure-Python tracer makes the suite impractically slow. The gate is nevertheless
unaffected by this change: `[tool.coverage.run] source = ["tensorcad"]` and this
diff touches no file under `tensorcad/`, so the measured total is identical to
`main`'s.

### 4.3 Behaviour of the two edited tools

* `tools/verify_merge_m3859.py` was executed end to end after the edit
  (120 paths audited, `1 lost -> FAIL`) — same verdict and same exit-code
  contract as before; only the dead `lost` list is gone.
* `tools/make_report_m3859.py` — the `%` → f-string rewrite was checked on the
  actual value: both forms render `'closed form beyond its 35.00° limit'`
  (identical strings). `py_compile` clean for both files.

## 5. Notes, limits, applicability

* **Scope.** The lint debt was pre-existing on `main`; this PR is only the first
  to be judged by a ruff that can see it. Keeping the fix inside the two audit
  scripts is what "minimal" means here — no rule was disabled, no `# noqa` was
  added to silence the new families, and `[tool.ruff.lint] select` is untouched.
* **The pin.** `ruff<0.16` is a reproducibility guard, not a cure. It should be
  raised together with a `ruff check .` pass, deliberately, not on a release day.
* **What was not attempted.** `mypy` and the neural-geometry workflow were not
  touched; `neural`-extra tests are out of scope for PR #199's red (they live in
  the separate `neural-geometry-experiments` workflow, which installs torch).
* **Interesting for the next agent.** A torch-free CI replay on a box that *has*
  torch needs a `ModuleNotFoundError` from `find_spec`, not an `ImportError`
  raised inside the module body: pytest re-raises the latter, so the suite dies
  at collection instead of skipping.