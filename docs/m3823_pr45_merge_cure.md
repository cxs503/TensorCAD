# PR #45 cure — `agent/m3391` ← `main`

**Verdict: merged clean, nothing lost, suite above the `main` baseline.**

PR #45 (branch `agent/m3391`, "Large-angle stability: righting-arm curves verified against
two exact closed forms") went dirty against `main`: the branch forked at `df31744`, while
`main` advanced to `2ff80b2` through the #21 → #39 merge-governor chain. GitHub reported
`mergeable=false`. The cure is a three-way merge of `main` **into** the PR branch, resolved
by hand where the text merge stopped, verified mechanically, and pushed back to the same
branch so `merge_governor` can re-check it on its next cycle.

## Revisions

| role | branch | SHA |
|---|---|---|
| merge base | — | `df31744cd5659c1fbb31b54607a13f833aa4a555` |
| ours (PR branch, before) | `agent/m3391` | `164fc634c624e26092e6931ee554707a5747bdfe` |
| theirs | `origin/main` | `2ff80b2b9c2285da440648b354382d9ff0cb90a8` |
| ours (PR branch, after) | `agent/m3391` | merge commit, recorded in `artifacts/m3823/push_record.json` |

The branch had **no** remote counterpart (`origin/agent/m3391` did not exist), so the cure
also publishes it. Old and new SHAs are recorded in `artifacts/m3823/push_record.json`.

## What the two sides changed

* **ours** — 11 files, 2080 added lines: `tensorcad/gz.py` (new: `heeled_hydrostatics`,
  `righting_arm_curve`, `GZ = KN − KG sin φ`, `with_freeboard`, free-surface correction),
  `tests/test_gz.py` (26 tests), `docs/gz.md`, `examples/gz/`, plus the public-surface,
  CHANGELOG and version bumps.
* **theirs** — 38 files, 3865 added lines: `tensorcad/transom.py` (transom-stern family),
  waterplane second moments on the `Hydrostatics` record, the `hull.py` analytic reference
  rework, `examples/transom_stern/`, `examples/asymmetric_hull/`, their tests, docs and
  tooling.
* **both** — 4 files: `CHANGELOG.md`, `pyproject.toml`, `tensorcad/__init__.py`,
  `tests/test_public_surface.py`. These are the only files where a text merge can silently
  drop a name, and they are the only files that got manual review.

## Conflict files and how each was resolved

| file | git outcome | resolution |
|---|---|---|
| `CHANGELOG.md` | **textual conflict** (`UU`) | Both sides had opened a `## 0.21.0` section. Main's section is the released, dated `## 0.21.0 — 2026-11-04` (transom family + waterplane inertias); ours is the GZ work that was never tagged. **Union: both sections kept, main's dated `0.21.0` first, ours relabelled `## Unreleased — Large-angle stability: …`** — the text of every line is preserved, only the heading number moves. Two identical version headings in one file is the defect, not the content. No line of either section was dropped (checked mechanically, see below). |
| `pyproject.toml` | auto (identical edits) | Both sides bumped `0.20.0 → 0.21.0`; the two hunks are **byte-identical**, so the merge is a single `0.21.0` and no intent is lost. No version bump beyond that: `main`'s release stays `0.21.0`, and `tests/test_public_surface.py`'s `__version__ == "0.21.0"` pin still holds. |
| `tensorcad/__init__.py` | auto (different regions) | Union of the public surface: ours adds the 14 GZ names (`RightingArmCurve`, `with_freeboard`, `heeled_hydrostatics`, …), theirs adds the 13 transom-stern names (`TransomModel`, `build_transom_mesh`, …). Merged `__all__` = **134** = `|ours ∪ theirs|` exactly (107 base + 14 + 13), with **zero** names lost and zero extra. Every name imports. |
| `tests/test_public_surface.py` | auto (adjacent hunks) | Both sides extended `EXPECTED_API`; the merged literal contains both groups (GZ block and the `wigley_reference_parameters` / transom block, main's trailing comment preserved). The version pin stays `0.21.0`, matching the merged `pyproject.toml`. |

The other **41** files were touched by one side only; each was required to come through
**byte-identical to that side's blob**, and all 41 did. That is the check that catches the
mechanical-overlay failure in its other direction (a whole-file copy from the wrong side,
or a stray whitespace hunk).

## How "nothing was silently lost" was established

`tools/verify_merge_m3823.py` (new) re-derives the merge base, then, for every file either
side touched, collects the lines that side **added** with respect to the base and requires
each one to appear in the merged tree under a whole-file, whitespace-normalised
containment test. The test is deliberately strict — a line that survives only as a prefix
of a longer line (the "list literal merged over the other side's entry" shape) counts as a
loss. Deletions get the mirror check, and the public surface (`__all__`) is compared as a
set against the union.

| check | result |
|---|---|
| added lines checked (ours / theirs) | 2080 / 3865 = **5945** |
| lines not found in the merged tree | **0** |
| allowlisted deviations | **1** — the CHANGELOG heading renumbering above, recorded with its reason in the tool |
| deleted lines checked | 73 (ours 0 + theirs 73); 21 still present after the merge, all but one in `tensorcad/hull.py`, a file **only main touched** — so those hits are the containment test matching identical strings elsewhere in main's own rewritten file, not union residues. The one in a shared file is main's own reformat `"wigley_reference_parameters",` → `"wigley_reference_parameters",  # transom-stern family (0.21.0)`, i.e. the name is *kept*, as the union requires. |
| single-side files byte-identical | 41 / 41 |
| public surface merged vs union | 134 / 134, lost 0, extra 0 |

Artefact: `artifacts/m3823/merge_content_check.json` (verdict `PASS`). Reproduce with
`python3 tools/verify_merge_m3823.py`; the tables and the figure are regenerated by
`python3 tools/report_m3823.py` from the JSON and the pytest logs, never typed by hand.

## Test evidence

Both runs are the same command in a pristine export of each tree, on the same machine,
with the two modules that need the absent `safetensors` package excluded:

```
python3 -m pytest -q --ignore=tests/test_cli_neural.py --ignore=tests/test_neural_pipeline.py
```

| tree | passed | failed | skipped | wall |
|---|---|---|---|---|
| `main` (`2ff80b2`) baseline | 818 | 2 | 11 | 1267 s |
| merged (`agent/m3391` + `main`) | **844** | 2 | 11 | 551 s |

* **+26 passed, no new failure** — exactly `tests/test_gz.py` (26 tests), the PR's own
  contribution.
* The **2 failures are identical and pre-existing on `main`**:
  `tests/test_geometry_field.py::test_the_learned_revision_speaks_the_same_protocol` and
  `::test_the_learned_field_gradient_is_autograd_and_points_outward` both die on
  `ModuleNotFoundError: No module named 'safetensors'` — an environment gap of this box,
  not a merge regression. The baseline exhibits the same two, at the same position in the
  run (`FF` at ~43 %).
* Targeted re-run of the merge-touching suites —
  `test_gz`, `test_transom`, `test_hydrostatics_inertia`, `test_public_surface`,
  `test_hull`, `test_analytic_hydrostatics`, `test_curvature`, `test_stability` —
  **251 passed, 0 failed** (`artifacts/m3823/logs/pytest_merged_targeted.log`). This is the
  integration that mattered: `gz.py` consumes the `hydrostatics`/`hull` code that `main`
  rewrote in the same window.
* `examples/gz/run_example.py` runs green on the merged tree and reports
  `worst |GZ − wall-sided| = 6.939e-17 m (0–38°)`, `1.110e-16 m` on the half-submerged log,
  `GZ'(0) = 0.421622996 m/rad` vs `GM = 0.421622996 m` — and it reproduces its committed
  outputs bit for bit (no diff in `examples/gz/outputs/`), so the merge did not perturb the
  numerics.

Logs: `artifacts/m3823/logs/`; parsed numbers: `docs/m3823/outputs/results.json`;
figure: `docs/m3823/outputs/merge_verification.png` (public-surface union, content-check
counts on a log axis, suite comparison); generated summary: `docs/m3823/report.md`.

## Notes, limits, applicability

* The cure is a **merge**, not a rebase: history keeps the PR's commits and adds main's, so
  the branch is a fast-forward candidate for `main` and `merge_governor` can re-verify it.
* The one editorial change — relabelling the PR's CHANGELOG section `Unreleased` — is a
  deliberate, allowlisted deviation, not a lost line; `pyproject.toml`/`__version__` stay
  `0.21.0` because that is what both sides independently bumped to, and a second bump would
  be an unrequested release.
* The 21 "deleted line still present" hits are a known false-positive class of the
  whole-file containment test (short, repeated lines such as `)` or `p = BREADTH_EXPONENT`).
  They are listed rather than filtered so the judgement is visible; none is in a file both
  sides edited except main's own comma reformat.
* Applicability: this cure resolves *text and content*. It cannot prove that the two
  features are semantically composable in ways the tests do not exercise — the coverage
  argued here is the 844-test suite plus the GZ example's own analytic claims.