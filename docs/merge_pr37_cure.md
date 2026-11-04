# PR #37 cure — semantic union of the transom-stern branch and `main`

Mission m3741 · repo **TensorCAD** · worktree `/root/tp-agent/work/m3741` · branch `agent/m3741`

## 1. What was broken

| item | value |
|---|---|
| PR | #37, head branch `agent/m2873` |
| remote head before the cure | `206698a` — *"m2873: rebuild on current main 6c568c31"* |
| `origin/main` at cure time | `8389c7e` — *"Merge #33 (m2575) via merge-governor"* (package 0.20.0, Wigley family) |
| merge base `206698a` vs `main` | `6c568c31` — the PR branch sat **7 commits behind** `main` |
| divergence | 7 commits on `main` / 1 commit on the branch |
| GitHub verdict | `mergeable = false` |

The branch was *not* dirty because of an evil hunk: it had been **rebuilt on an older
`main`** (a previous cure of "no common ancestor"), so every release note, the version
literal and the public-surface list on `main` collided with the branch's own additions.
The dangerous way to clear this is a mechanical overlay (`git checkout --ours/--theirs .`),
which is exactly how content was silently lost in this repository before. This cure
therefore treats the merge as a **union problem with a machine-checked proof**.

## 2. Method

```
git merge --no-ff --no-commit origin/agent/m2873      # onto origin/main
# 4 conflicts, 13 files merged cleanly
# each conflict resolved to the union of both sides, then:
python3 artifacts/merge_pr37_audit/audit_union.py      # path / line / token level proof
python3 artifacts/merge_pr37_audit/postmerge_smoke.py  # both families still compute
python3 -m pytest -q ...                               # full suite vs main baseline
```

The merge commit is `87542a3` with parents `8389c7e` (`origin/main`) and `206698a`
(`origin/agent/m2873`), so it **fast-forwards** the PR head instead of rewriting it.

## 3. Conflict files and their resolutions

| file | conflict | resolution (union) |
|---|---|---|
| `pyproject.toml` | `version = "0.20.0"` (main) vs `"0.19.0"` (branch) | single `0.21.0` — the merged tree ships both families, so neither parent's number is true |
| `tensorcad/__init__.py` | `__all__` line `"generate_wigley"` vs `"generate_transom_hull"`; `__version__` | **both names kept** (alphabetical: `generate_hull`, `generate_transom_hull`, `generate_wigley`); `__version__ = "0.21.0"` |
| `tests/test_public_surface.py` | `EXPECTED_API` = main's Wigley set vs the branch's 12 transom names; version assertion; layered-module list | `EXPECTED_API` = **Wigley set ∪ transom set** (branch block kept verbatim, comment records that it came from the 0.19.0 branch); module list keeps **both** `"wigley"` and `"transom"`; assertion follows `0.21.0` |
| `CHANGELOG.md` | main's `## 0.20.0` section vs the branch's `## 0.19.0` section | branch's three transom bullets kept **verbatim** and renumbered `## 0.21.0 — 2026-11-04`, placed **above** main's `0.20.0` entry; no bullet dropped, no duplication |

Thirteen files merged without help and were **not** touched by hand:
`tensorcad/transom.py`, `tests/test_transom.py`, `tools/report_m2873.py`,
`examples/transom_stern/run_example.py`, `docs/transom_stern.md`,
`docs/m2873/outputs/transom_stern/{convergence.csv,design_table.csv,reference_ladder.csv,*.png}`
(the branch's five artefacts) — i.e. the whole transom deliverable arrived intact, and
main's Wigley files (`tensorcad/wigley.py`, `tests/test_wigley.py`, …) were never at risk.

## 4. Proof that nothing was lost (`artifacts/merge_pr37_audit/union_audit.json`)

`audit_union.py` compares the merge against **both parents** on three levels:

* **path level** — `main_only_missing = []`, `branch_only_missing = []`, `deleted_vs_main = []`;
* **line level** — every line a parent added relative to the merge base is either present
  verbatim, present after whitespace/re-wrapping normalisation, or explicitly *superseded*
  by one of six auditable rules (a singleton version literal ×3, an assertion ×1, a
  release heading ×1 … 6 lines total, all printed in the JSON with their rule);
* **token level** — every identifier token of a re-wrapped line must survive in the file,
  which is what catches a list that lost one API name.

Verdict: **`UNION-COMPLETE`**, `lost_lines_total = 0`, 40 files audited, no conflict
marker anywhere in the tree.

## 5. Tests

Command (identical on both trees):

```
python3 -m pytest -q --ignore=tests/test_cli_neural.py --ignore=tests/test_neural_pipeline.py -p no:cacheprovider
```

| tree | collected | passed | failed | skipped |
|---|---|---|---|---|
| merged (`agent/m3741`, union tree) | 815 | **802** | 2 | 11 |
| main baseline, shared files (derived: merged − the branch's 96 new tests) | 719 | 706 | 2 | 11 |

* The two pre-existing collection errors (`tests/test_cli_neural.py`,
  `tests/test_neural_pipeline.py`) are `ModuleNotFoundError: No module named 'safetensors'`
  and occur **identically on `main`** — an image gap, not a merge effect.
* The 2 merged failures are `tests/test_geometry_field.py::test_the_learned_revision_speaks_the_same_protocol`
  and `::test_the_learned_field_gradient_is_autograd_and_points_outward`, both
  `safetensors` at `tensorcad/neural/field.py:26`. Run in isolation on a pristine
  `git archive origin/main` tree they fail the same way (47 passed, 2 failed) — the merged
  tree **adds no failure**.
* All 96 tests in `tests/test_transom.py` (the branch's own suite) pass on the merged tree.
* Honest gap: the *full* main-baseline run was started twice and both times was killed at
  the 120 s harness limit because the box was at load average 176–350 (other agents). The
  baseline number in the table is therefore derived from the merged run minus the branch's
  own test file, not measured directly; the two failures were verified on `main` in
  isolation instead.

`artifacts/merge_pr37_audit/postmerge_smoke.py` is the independent cross-check that the
merge is also *numerically* alive (both families, `postmerge_smoke.json`, all 8 checks pass):

| check | value | reference | rel. error | tol |
|---|---|---|---|---|
| Wigley `Cb` closed form | 0.444444 | 4/9 | 0.0 | 1e-12 |
| Wigley `V` vs `Cb·L·B·T` | 2.777778 m³ | 2.777778 m³ | 0.0 | 1e-12 |
| Wigley LCB symmetry | 5.0 m | L/2 | 0.0 | 1e-12 |
| transom `A_t/A_m` at τ = 0 | 0.0 | 0.0 | 0.0 | 1e-18 |
| transom τ = 0 mesh volume | 8.216624 m³ | 8.281656 m³ | **0.785 %** | 2 % |
| transom τ = 0.25 mesh volume | 8.243855 m³ | 8.308624 m³ | **0.780 %** | 2 % |
| area-ratio inversion round trip | 0.25 | 0.25 | 1.1e-16 | 1e-12 |
| topology τ = 0 / τ = 0.25 | watertight, edge-manifold, outward, no zero-area faces | same | — | — |

Both mesh-vs-closed-form errors are **< 3 %**; they are single-grid (41×13) discretisation
biases of the same sign as the branch's published convergence ladder
(`docs/m2873/outputs/transom_stern/convergence.csv`: −0.88 % → −0.31 % → −0.14 % under
refinement), not merge damage.

## 6. Push status / SHA record

| | SHA |
|---|---|
| old PR head (`origin/agent/m2873`) | `206698a` |
| merge base with `main` | `6c568c31` |
| `origin/main` | `8389c7e` |
| **new head** (merge, fast-forwardable over `206698a`) | **`87542a3`** |

`git push` from this sandbox is refused by the security policy (`REJECTED: 命令被安全策略拦截`,
as is any direct GitHub access), so the new head was **not** pushed by hand. It lives on
branch `agent/m3741`, whose history contains `206698a` as an ancestor, so
`git push origin agent/m3741:agent/m2873` is a plain fast-forward for whoever holds
credentials — no force, no history rewrite.

## 7. Reproduce

```
git merge --no-ff origin/agent/m2873          # expect exactly the 4 conflicts above
python3 artifacts/merge_pr37_audit/audit_union.py
python3 artifacts/merge_pr37_audit/postmerge_smoke.py
python3 artifacts/merge_pr37_audit/make_figures.py
```

Figures: `artifacts/merge_pr37_audit/fig1_union_and_parity.png` (union audit + test parity),
`artifacts/merge_pr37_audit/fig2_postmerge_numerics.png` (branch convergence data after the
merge + smoke errors against tolerance).