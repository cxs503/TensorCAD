# Curing PR #37 — `agent/m2873` × `main`

**Mission** m3682 · **repo** TensorCAD · **PR** #37 (transom-stern hull family)
**Branch cured in place** `agent/m2873`: `33243ef` (old) → `571af0f` (merge)
**Merge parents** `33243ef` (branch) + `8389c7e` (`origin/main`, "Merge #33 (m2575) via merge-governor")
**Merge base** `d75dd54`

The PR reported `mergeable=false`: the branch was 13 commits behind `main` and four
files conflicted. The cure is a **union merge on the branch**, plus a machine audit
that the union really is a union.

## 1. Why it was dirty

| | commits | files touched vs. base |
|---|---|---|
| `agent/m2873` | 1 ahead | 15 (`tensorcad/transom.py`, `tests/test_transom.py`, `docs/transom_stern.md`, `examples/transom_stern/`, `outputs/transom_stern/`, `tools/report_m2873.py`) |
| `origin/main` | 13 ahead | 65 (Wigley, Bonjean, stability, asymmetric hull, README/backlog) |

Both sides touched the same four "registry" files — `CHANGELOG.md`, `pyproject.toml`,
`tensorcad/__init__.py`, `tests/test_public_surface.py` — which is exactly the shape
where a mechanical `-X ours/-X theirs` overlay silently drops one side's content.

## 2. Per-file resolution

| file | conflict | resolution | why |
|---|---|---|---|
| `tensorcad/__init__.py` | 6 hunks (import block + 5 in `__all__` + `__version__`) | rebuilt from `main`'s file: main's 107-name `__all__` and import block kept **verbatim and in place**; the branch's `from .transom import (…)` block re-inserted before `from .validation import`; its 13 public names appended as a labelled group; `__version__ = "0.21.0"` | branch had also *reordered* `DESIGN_SCHEMA`/`SUBOFF_*`; taking both hunks would have duplicated those names. Appending keeps main's diff history readable and yields 120 unique names = 107 ∪ 13 |
| `pyproject.toml` | `version` | `0.21.0` | neither `0.14.0` (branch) nor `0.20.0` (main) is right for "main + new public API"; see §3 |
| `tests/test_public_surface.py` | `version` assert; module tuple | assert `"0.21.0"`; tuple lists **both** `"wigley"` (main) and `"transom"` (branch) on one line | the two sides added one module each; folding them onto one line keeps both |
| `CHANGELOG.md` | whole leading block | branch's transom section re-headered `## 0.14.0 → ## 0.21.0` and placed **above** main's `0.20.0 / 0.17.0 / 0.16.0 / 0.14.0 / 0.13.0` sections; every line of both sides retained | **main had already shipped a `0.14.0` section**, so keeping the branch's number would have produced two `0.14.0` entries. This was the one non-obvious trap of the merge |

Non-conflicting auto-merges worth naming: `main` deletes the tracked `.ruff_cache/`
files (cache, `.gitignore`d); that deletion is **kept deleted**, not resurrected from
the branch. All 12 branch files (`tensorcad/transom.py`, `tests/test_transom.py`,
`docs/transom_stern.md`, `examples/transom_stern/run_example.py`,
`outputs/transom_stern/*`, `tools/report_m2873.py`) merged in untouched, and all 65
main-side paths arrived untouched.

## 3. Version decision

`0.21.0` everywhere (package, pytest assertion, changelog header). Rationale: the merge
adds public API to a tree where `main` is the release line at `0.20.0`; a minor bump is
the smallest consistent choice, and it also removes the `0.14.0`-section collision
described above. The branch's `0.14.0` label was stale by construction (it predates
`main`'s Wigley/Bonjean/stability releases).

## 4. Union audit (anti-silent-loss)

`docs/merge_pr37/verify_union.py` re-derives both sides from git and checks the merge
commit against them, so the "nothing was dropped" claim is a test, not a promise:

* **file set** — every path either side added/modified vs. base `d75dd54` exists in the
  merge, unless a side deleted it on purpose (here: only `main`'s `.ruff_cache/`);
* **line set** — every non-trivial line either side *added* is present in the merged
  file, after normalisation; eight deliberate relabels are enumerated in
  `EXPECTED_RELABELS` and are themselves fatal if they stop matching (the five version
  strings, the two module-tuple lines folded into one, and one superseded comment
  marker).

Result (`docs/merge_pr37/union_audit.json`):

```
files_checked       76
files_missing       []          # both sides' added/modified paths all present
lines_dropped       []          # 0 dropped added-lines
lines_relabelled    8           # all listed above, each verified present in merge
public_surface      version 0.21.0, |__all__| = 120, duplicates []
ok                  true
```

## 5. Test evidence

Same command, same interpreter, same environment, `main` **before** the merge vs. the
merge commit:

```
python3 -m pytest tests/ -q --ignore=tests/test_cli_neural.py --ignore=tests/test_neural_pipeline.py
```

| tree | passed | failed | skipped | wall |
|---|---|---|---|---|
| `main` `8389c7e` (baseline) | 706 | 2 | 11 | 339.9 s |
| merge `571af0f` (cured branch) | **802** | 2 | 11 | 291.4 s |

* The **same two** failures exist on both sides —
  `tests/test_geometry_field.py::test_the_learned_revision_speaks_the_same_protocol` and
  `::test_the_learned_field_gradient_is_autograd_and_points_outward` — and both are
  `ModuleNotFoundError: No module named 'safetensors'` on this image. That is why
  `test_cli_neural.py` and `test_neural_pipeline.py` are excluded (they cannot even be
  collected for the same reason). **No regression, no new failure.**
* `+96` passing tests come from the branch's own `tests/test_transom.py` plus the
  public-surface assertions now covering both families.
* Logs: `runs/baseline_main.log`, `runs/merged_tests.log`.
* Branch acceptance still green after the merge: `examples/transom_stern/run_example.py`
  → **exit 0**, worst reference disagreement **1.11e-15**, converged mesh errors
  (τ=0.35: −0.770 / −0.274 / −0.125 / −0.044 % on 41×13 → 321×81), topology
  watertight/manifold/Euler-2 at every τ. Log: `runs/transom_example2.log`;
  figure: `docs/merge_pr37/merge_evidence.png`.

## 6. What this cure does *not* claim

* The two `test_geometry_field` failures are an image-level dependency gap
  (`safetensors` missing), not a code defect of either side; they are left as-is and
  reported as pre-existing.
* The merge does not re-verify the branch's numerics beyond the shipped example/tests —
  `merge_governor` re-runs the merge next cycle, and the physical claims (Cb/Cp/Cwp,
  transom area, convergence orders) remain those documented in `docs/transom_stern.md`.
* `runs/` and `outputs/` logs are evidence only; nothing was regenerated by hand.

## 7. Reproduction

```bash
git checkout -b cure agent/m2873
git merge --no-ff main            # 4 conflicts, resolved as in §2
PYTHONPATH=. python3 docs/merge_pr37/verify_union.py \
    --base d75dd54 --ours 33243ef --theirs 8389c7e --out docs/merge_pr37/union_audit.json
python3 -m pytest tests/ -q --ignore=tests/test_cli_neural.py --ignore=tests/test_neural_pipeline.py
PYTHONPATH=. python3 examples/transom_stern/run_example.py
PYTHONPATH=. python3 docs/merge_pr37/make_figure.py
```
## 8. Branch state and what the agent could not do

| ref | SHA | meaning |
|---|---|---|
| `agent/m2873` (GitHub, PR #37) | `33243ef` | still the pre-cure tip at the time of writing |
| `agent/m2873` (local, other worktree) | `33243ef` | untouched |
| cure merge commit | `571af0f` | parents `33243ef` + `8389c7e`, clean against `main` |
| this branch's tip | see §0 of the commit log | cure merge + this evidence set |
| local `agent/m2873-cure` | same as the tip | the cure, ready to be pointed at `agent/m2873` |

Publishing with `git` from an agent shell is blocked by harness policy (agents do not
push; the platform pusher does), so step 4 of the mission is **executed by the platform
pusher, not by this agent**: the pusher publishes this worktree's branch as
`agent/m<mid>` (`agent/m3682`), and `merge_governor` re-verifies PR #37's mergeability on
its next cycle. The cure itself is branch-agnostic — the local ref `agent/m2873-cure`
already contains `33243ef` as an ancestor, so repointing `agent/m2873` at it is a
fast-forward that lands exactly the tree audited above:

```
git branch -f agent/m2873 fcd8f84      # or, via the pusher: ref agent/m2873-cure -> agent/m2873
```

Honest status: **merge resolved and verified locally; publication of the `agent/m2873`
ref is pending the platform pusher.**

## 9. Error judgement (key quantities vs. baseline)

| quantity | reference | measured after merge | relative error | < 3 % ? |
|---|---|---|---|---|
| incomplete taper integral, `p ∈ {0.3…3}`, `τ ∈ {0…0.45}` | independent double-exponential rule / `Γ((1+p)/2)/(√π Γ(1+p/2))` at `τ=0` | `0.837355517709` … `0.424413181578` | worst **1.11e-15** | **yes** |
| volume `V`, `τ = 0 / 0.05 / 0.2 / 0.35`, 41×13 | closed form of the truncated family | 7.446554 … 8.113840 m³ | −0.880 … −0.770 % | **yes** |
| volume `V`, same, finest rung 321×81 | as above | 7.509115 … 8.173195 m³ | −0.048 … −0.044 % | **yes** |
| waterplane `A_wp`, 41×13 → 321×81 | as above | — | −0.175 … −0.002 % | **yes** |
| LCB, `τ > 0`, 41×13 | as above | 4.7910 / 4.3746 / 4.1076 m | +0.058 … 0.065 % | **yes** |
| topology (all `τ`) | watertight · edge-manifold · `V−E+F=2` · outward | 4/4 satisfied | exact | **yes** |
| core test suite | `main` baseline 706 passed / 2 failed / 11 skipped | 802 passed / 2 failed / 11 skipped | 0 regressions, +96 passing | **yes** |

The merge-invariant families (Wigley, Bonjean, stability) are not re-baselined here
beyond the suite: they pass identically on both sides, so the merge introduced no drift.
