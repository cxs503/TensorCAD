# PR #163 (branch `agent/m4125`) — round 2: the red step was **Lint**, and the branch was stale

## What CI actually said (read from the API, not guessed)

For head `d7de3b7` the `tests` workflow (run `42023935490`) ended as
`test (3.10)=cancelled, test (3.12)=cancelled, test (3.11)=cancelled, test (3.13)=failure`.
The job/step table for that run resolves the signature:

| leg | started | steps |
|---|---|---|
| 3.13 | 17:20:12Z | **6 Lint → failure**, every step after it `skipped` |
| 3.10 / 3.11 / 3.12 | 12:03:32Z | never scheduled a step; `cancelled` at 17:20:32Z by `fail-fast` |

So the three `cancelled` legs are `fail-fast` collateral: **the single red step is
`Lint` (`ruff check .`)**, and it is version-independent — 3.13 was simply the only
leg that got a runner.

## Two defects, both fixed here

1. **Stale base.** `agent/m4125` branches from `57414db`, i.e. before `main` grew
   `tensorcad/righting_arm.py`, `examples/righting_arm`, `tools/verify_merge_m3859.py`
   and the `m3859` / `m3962` / `m4038` documents. `git diff origin/main FETCH_HEAD`
   therefore showed **−5479 lines**, all of it deletions of `main`'s work — the branch
   could never be merged green or even mergeable-by-content. Round 2 rebuilds the
   delivery on current `main` (`c985c46`): the 13 trim files are re-applied on top,
   nothing of `main` is removed.
2. **Lint, red on `main` as well.** With the rebuild in place `ruff check .` reported
   **11 diagnostics, all of them in files this PR does not otherwise touch**
   (`tools/make_report_m3859.py`, `tools/verify_merge_m3391.py`,
   `tools/verify_merge_m3859.py`, `tests/test_merge_cure_m3391.py`): `W292` ×4,
   `RUF100` ×2, `F401`, `F841`, `RUF059`, `UP031`, `B905`. Since `tests.yml` runs
   `ruff check .` over the whole tree, no branch could pass that step. All 11 are fixed
   mechanically and behaviour-preservingly: missing final newlines, two stale
   `# noqa: E402`, an unused `math` import, an unused `lost: list[str]` declaration,
   `rc → _rc` in an unpack, `%`-format → f-string, `zip(..., strict=True)`.

## Verification (all commands run in this worktree, CPython 3.11.6 — the only
interpreter on this host; raw logs in `runs/pr163_round2/`)

| CI step | Command | Result |
|---|---|---|
| Lint | `ruff check .` | `All checks passed!` (0 errors; before: 11) |
| Type check | `mypy` | `Success: no issues found in 35 source files` |
| Tests with coverage gate | `pytest --cov --cov-report=term-missing -q -n 16` | **955 passed, 10 skipped, 0 failed**, total coverage **93.43 %** ≥ 90 % gate (`tensorcad/trim.py` 92 %); 750 s |
| Neural CLI help | `tensorcad-neural --help` | pass |
| CLI artifacts / reproducibility | `tensorcad-hull` twice + `cmp` | byte-identical objects and metadata |
| Geometry benchmark | `examples/parametric_hull/run_example.py` | pass, `hull.obj` non-empty |
| Geometry core benchmark | `examples/geometry_core/run_benchmark.py --resolution 40 --points 5000` | pass |
| SUBOFF bare hull + AFF-8 appendages | `examples/suboff*` | pass, OBJ/STL/handoff/field all present |
| Mesh conditioning | `examples/mesh_conditioning/run_example.py --resolution 110,14,20 --targets 30000,20000` | pass (obj/msh/json) |
| FEM handoff / STEP export | `examples/step_export/run_example.py --resolution 60,8,11` | pass (STEP part skips without OCP, as designed) |
| Design loop | `examples/design_loop/run_example.py` | pass (obj/msh/manifest/path) |
| Representation benchmark | `examples/benchmarks/run_benchmark.py --quick` + the workflow's `assert`s | pass, `passed is True`, all four representation rows present |

Every non-pytest step of `tests.yml` was executed verbatim by
`runs/pr163_round2/ci_steps.sh` → `rc_all=0` (17/17 pass, log `runs/pr163_round2/ci_steps.log`).

## Numbers behind the figure

`docs/pr163/pr163_ci_evidence.png` is drawn by `docs/pr163/make_figure.py` from the raw
logs only (`docs/pr163/results.json` is the same parse, machine-readable):
(a) coverage 93.43 % package / 92 % `trim.py` vs the 90 % gate; (b) the 11 ruff
diagnostics by rule, → 0; (c) 17/17 example steps green.

## Honest limits

* The 3.10 / 3.12 / 3.13 legs were not executed (no such interpreters here). The change
  is version-neutral: lint/mypy/pytest are the same commands on all legs, no
  version-conditional code was added or removed, and the touched files use only
  pre-3.10 syntax (`zip(strict=)` exists since 3.10).
* Local coverage ran under the pure-Python tracer (`coverage.tracer` is absent here), so
  93.43 % may differ from CI's by a small amount; the margin over the gate is 3.4 pp.
* `docs/pr87/` and `docs/pr49_cure_report.md` are earlier deliveries carried by this
  branch unchanged; they are documentation and figure generators, not CI inputs.