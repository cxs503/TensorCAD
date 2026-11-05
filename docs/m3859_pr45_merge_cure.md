# PR #45 cure — merging `agent/m3391` into `main` without losing either side

**Mission** m3859 · **repository** TensorCAD · **work branch** `agent/m3859`
**merge commit** `6153c6a` (parents `2ff80b2` = `origin/main`, `164fc63` = `agent/m3391`)
**old `agent/m3391` SHA** `164fc634c624e26092e6931ee554707a5747bdfe` -> **cure SHA**
`6153c6ad3187c4327dfbc05e077b499fb82fe875`; the cure is a fast-forward of `agent/m3391`
(no rewrite), carried by the platform automatic branch push because the push subcommand is
refused inside this sandbox by tool policy (see section 5).

| item | value |
|---|---|
| PR state before cure | dirty, `mergeable=false` |
| conflict files | **1** — `CHANGELOG.md` |
| files auto-merged but audited | 119 further paths (85 of them touched by both sides) |
| content lost | **0 paths, 0 unique lines** |
| tests merged vs main | **844 passed / 2 failed / 11 skipped** vs **818 / 2 / 11** |
| regressions | none (the 2 failures are the same pre-existing tests on both trees) |

## 1. Case

TensorCAD is the hull-geometry/hydrostatics kernel of the ship-performance stack: parametric
and imported hull surfaces, mesh conditioning, hydrostatics, stability. Two lines of work met
in PR #45:

* **main** (`0.21.0`, 2026-11-04) — waterplane second moments (`∫y²dA`, `∫(x−LCF)²dA`, `keel_z`)
  on the `Hydrostatics` record and on `analytic_hydrostatics`, plus the **transom-stern hull
  family** (`tensorcad/transom.py`, truncated-taper integrals, `A_t/A_m = sin(πτ)^p` inversion).
* **`agent/m3391`** (`164fc63`) — **large-angle stability**: `tensorcad/gz.py` (heeled
  displacement solve by Newton with the exact derivative `dV/dc = A(c)`), righting-arm curve
  `GZ = KN − KG sinφ`, `max_gz` / range / dynamical stability, `with_freeboard`, free-surface
  correction, verified against the box wall-sided trapezoid and the half-submerged log.

Both sides independently opened a `## 0.21.0` changelog section and bumped `pyproject.toml`
to `0.21.0` — which is exactly the shape that turns into a "clean" `git status` with one
side's notes quietly gone. That is the failure mode this cure had to exclude, not just the
`mergeable=false` flag.

Geometry/parameters that the merged tree must still reproduce (the branch's verification
hull): box barge 10 × 2 × 2 m with KG such that `GM₀ = 0.3167 m`, `BM₀ = 0.4167 m`; log of
radius 1 m in 64 segments; parametric hull L = 10 m, B = 2 m, T = 0.8 m, KG = 0.5 m,
freeboard 1 m.

## 2. Method

1. Worktree on `main` @ `2ff80b2`; `git merge 164fc63` (three-way, merge base `ca0dde3`).
2. Conflicts resolved by hand for semantic union; **no file was taken wholesale from one side**.
3. `tools/verify_merge_m3859.py` audits the *result*, not the exit code: for all 120 paths the
   two sides touched since the base it requires
   * single-side paths to be byte-identical to that side (binary files included — the check
     compares bytes, so a PNG cannot slip through a text-mode comparison), and
   * for paths touched by **both** sides, that the merge keeps every line that is *unique* to
     either side (multiset of merged lines ⊇ each side's unique lines). Lines that cannot
     survive because both sides wrote different statements for the same place are reported as
     `SUPERSEDED` with the line quoted, never as a silent pass.
4. Tests: every test file run in its own pytest process on **both** trees with identical flags
   (`-q -p no:randomly -p no:cacheprovider`, `OMP_NUM_THREADS=1`), 24 at a time, then compared
   file-by-file. The comparison is what the acceptance criterion is ("not below the main
   baseline"), so the baseline is a full run of `main` @ `2ff80b2`, not a cached log.
5. Evidence packaged by `tools/make_report_m3859.py` from the raw `RESULTS.tsv` logs and
   `examples/gz/outputs/results.json` — no number in this report is typed by hand.

Toolchain: Python 3.11.6, numpy 2.4.6, scipy 1.16.0, pytest 9.2.0, git 2.x, matplotlib (Agg).

## 3. Conflict resolution, file by file

| file | kind | resolution |
|---|---|---|
| `CHANGELOG.md` | **conflict (both sides added `## 0.21.0`)** | Union. Main's `## 0.21.0 — 2026-11-04` entry (waterplane inertias + transom family) stays first with its `### Added` list; the branch's entry becomes `### Large-angle stability: the righting-arm curve, and two exact benchmarks` **under the same 0.21.0 heading**, body verbatim. Both sides bumped the same unreleased version, so one version section with two subsections is the honest shape — no invented version number, no dropped bullets (39 main lines + 26 branch lines all present, verified line-by-line). Only one line changed level: the branch's `##` heading → `###`. |
| `tensorcad/__init__.py` | both sides edited the same block | git auto-merged; audited: 134 exported names, a strict superset of main's 120 and the branch's 121 (`from .gz import (...)` block plus 14 new `__all__` entries all present). |
| `tests/test_public_surface.py` | both sides edited `EXPECTED_API` | git auto-merged; the 26 public names in the merged file are a superset of both sides. Two `SUPERSEDED` flags are cosmetic: main had appended `# transom-stern family (0.21.0)` to the `"wigley_reference_parameters"` line, and its `"wigley", "transom", …` line replaced the branch's spacing of the same list — every name survives. |
| `tensorcad/gz.py`, `tests/test_gz.py`, `docs/gz.md`, `examples/gz/*` | added by the branch only | taken unchanged from `agent/m3391`, byte-identical (incl. `plot.png`). |
| `tensorcad/transom.py`, `tensorcad/hydrostatics*.py`, `tensorcad/hull.py`, `docs/transom_stern.md`, `tests/test_transom.py`, `tests/test_hydrostatics_inertia.py`, `artifacts/m3768/*`, `docs/asymmetric_hull/*`, `docs/m2873/*` | added/changed by main only | unchanged from `main`, byte-identical. `tensorcad/hull.py` is identical on both tips (the branch had already merged main), so the merge had nothing to choose there. |

## 4. Results

### 4.1 Content audit (`tools/merge_content_check.json`)

| metric | value |
|---|---|
| paths considered (union of both sides' changes since base) | 120 |
| touched by both sides | 86 |
| paths lost | 0 |
| unique lines lost | 0 |
| flags needing review | 1 `TRANSFORMED` (changelog heading, intentional) + 1 `SUPERSEDED` (comment/whitespace variant in `test_public_surface.py`, all names present) |

### 4.2 Tests (`docs/m3859/outputs/results.json`, raw logs in `artifacts/m3859/logs/`)

| tree | passed | failed | skipped |
|---|---|---|---|
| merged `6153c6a` | **844** | 2 | 11 |
| main `2ff80b2` | 818 | 2 | 11 |
| Δ | **+26** (the branch's `tests/test_gz.py`) | 0 | 0 |

Per-file comparison over 35 test files: 34 files identical to the baseline, 1 file
(`test_gz.py`) new, **0 files below baseline**. The two failures are
`test_geometry_field.py::test_the_learned_revision_speaks_the_same_protocol` and
`::test_the_learned_field_gradient_is_autograd_and_points_outward`, failing with the same
names on `main` — pre-existing, unrelated to this merge. `tests/test_cli_neural.py` and
`tests/test_neural_pipeline.py` cannot be collected on this host on either tree
(`ModuleNotFoundError: safetensors`), so both runs exclude them identically.

### 4.3 Numbers the merge had to preserve (branch evidence, re-pinned on the merged tree)

`tests/test_gz.py` (26 tests) passes on the merged tree, which re-asserts the branch's claims;
the values below come from `examples/gz/outputs/results.json`.

| quantity | merged-tree value | reference | error |
|---|---|---|---|
| box barge GZ vs wall-sided `sinφ(GM₀+½BM₀tan²φ)`, 0–35° | worst 6.94e-17 m | closed form | 6.9e-17 m |
| same form beyond its 38.66° limit (45°, 70°) | up to 1.8e-02 m deviation | form invalid there | expected |
| log in a circular section `(R−KG)sinφ`, 33 heels 0–90° | worst 1.11e-16 m | exact | 1.1e-16 m |
| `GM` from `initial_stability` vs from `dGZ/dφ` at 0 | 0.4216229965 vs 0.4216229958 m | each other | 1.8e-09 relative |

## 5. Traps met (and how they were handled)

* **The dangerous "clean" merge.** Only `CHANGELOG.md` conflicted; the other 85 jointly-touched
  paths merged silently. A green `git status` says nothing about content, so loss was measured
  with the unique-line criterion instead of trusted (this is the mechanical-overlay failure
  mode the mission warns about). A prior `cure/m3391` attempt in the tree used an
  `## Unreleased` heading for the branch's block; that is also defensible, but this cure keeps
  a single 0.21.0 section because both `pyproject.toml`s say 0.21.0.
* **`git show` in text mode corrupts binary comparisons** (`\r\n` → `\n`), which made every PNG
  look "lost" on the first pass. The audit now compares bytes.
* **A 128-core host at load ≈ 640** (other agents): a per-file pytest sweep hits the 100 s wall
  clock on `test_benchmark.py` on *both* trees, so it was re-run alone and in halves
  (3 + 29 tests) to get a real number instead of reporting a timeout as parity.
* **The push subcommand is refused in this sandbox** (`REJECTED: 命令被安全策略拦截`), so the
  commit could not be published from inside the mission. The cure is one merge commit on
  `agent/m3859` whose first parent is `origin/main` `2ff80b2` and second parent is the old
  `agent/m3391` tip `164fc63`, so landing it on `agent/m3391` is a fast-forward (no rewrite,
  no force) that the platform performs automatically.
* Scope limits: solver numerics were not re-derived here; this mission owns merge correctness
  and regression parity. Depth of evidence = one full test sweep per tree, no numerical
  re-validation of main's transom work.

## 6. Error judgement against the acceptance criterion

| criterion | measured | verdict |
|---|---|---|
| conflict files resolved without losing either side | 0 paths / 0 unique lines lost over 120 paths, 86 joint | **met** |
| tests not below the main baseline | 844 vs 818 passed, 2 vs 2 failed (same tests), 11 vs 11 skipped, 0 regressions over 35 files | **met** |
| PR mergeable again | merge commit `6153c6a` has `origin/main` `2ff80b2` as a direct parent | **met** (governor re-checks) |
| numeric error of the branch's claims | GZ vs closed forms ≤ 1.1e-16 m; GM cross-check 1.8e-09 relative | **met, ≪ 3 %** |

Not done / stated plainly: the two pre-existing `test_geometry_field` failures were left
untouched (out of scope, they fail identically on `main`); the neural test modules remain
uncollectable for a missing `safetensors` in this environment.