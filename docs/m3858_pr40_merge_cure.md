# PR #40 cure — `agent/m3169` rebased onto main by semantic union (mission m3858)

**Result: the conflict is resolved losslessly; the cured branch merges into main with
zero conflicts. The network push of the branch was refused by this sandbox's command
policy, so `origin/agent/m3169` still points at `fcf15ad` — see §5 for the exact
one-command hand-off.**

## 1. Shapes of the two sides

| role | ref | SHA |
|---|---|---|
| merge base | `df31744` (PR #38 merge) | `df31744cd5659c1fbb31b54607a13f833aa4a555` |
| PR head ("ours") | `origin/agent/m3169` | `fcf15adccd0dcbae671868dcf51e064c5e398946` |
| target main ("theirs") | `origin/main` (PR #39 merge) | `2ff80b2b9c2285da440648b354382d9ff0cb90a8` |
| cure merge commit | `cure/m3858-m3169` | `c18a88bf171edfd4ac7cd6296d784ef82596a2fd` |

Since the merge base the branch adds 13 files (`+2524/−3`, the righting-arm curve
deliverable) and main adds 38 files (`+4425/−71`: waterplane inertias, the transom-stern
hull family, PR #37/#39 deliverables). The two-changed-file overlap is only **four**
paths, which is why the dirty PR is cheap to cure but easy to cure *wrongly*:
three of the four are registries where a one-sided pick silently deletes work.

## 2. Per-file resolution

| file | type | resolution |
|---|---|---|
| `CHANGELOG.md` | conflict | Both sides opened a `## 0.21.0` section. One section now carries the branch's righting-arm items **verbatim** and, under a `### Added` heading, main's waterplane-inertia / transom-stern / `examples/transom_stern` items **verbatim**. No bullet was reworded, merged or dropped. |
| `tests/test_public_surface.py` | conflict | The layered-module tuple listed `stability, righting_arm, bonjean` on the branch and `transom` on main. Resolution is the union of the ten module names (`stability, righting_arm, bonjean, wigley, hull, design, suboff, suboff_appendages, cli, transom`), reflowed to the file's line width; the `EXPECTED_API` set was auto-merged and already held both families. |
| `tensorcad/__init__.py` | auto (overlap) | Git merged both import blocks and both `__all__` additions. Verified: 441 tokens ≥ 418 (main) and ≥ 412 (branch); `RightingArm`, `wall_sided_righting_arm`, `TransomModel`, `transom_fraction_for_area_ratio` all importable, `__version__ == "0.21.0"`. |
| `pyproject.toml` | auto (overlap) | Both sides bumped `0.20.0 → 0.21.0` identically; no textual conflict, no semantic difference. |

Nothing was taken "wholesale": every other path was changed by exactly one side, and the
audit proves each of those landed with *that* side's blob (§3, `unexpected_edits = 0`).

## 3. Preservation audit (machine check, not self-report)

`artifacts/m3858/verify_preservation.py` audits a merge commit on three invariants and
exits non-zero if any fails:

1. **path presence** — every blob of either parent exists in the merge;
2. **token multiset containment** — for every file both parents touched, no token
   (identifier, quoted string, number, word) has fewer occurrences in the merge than on
   either side. Line containment is deliberately *not* used: a semantic union legally
   reflows lines, and a checker that forbids reflow pushes the agent back to one-sided
   picks;
3. **binary preservation** — a binary must be one parent's blob, byte for byte.

Result for `c18a88b` (`artifacts/m3858/preservation_audit.json`):

| overlap file | tokens ours | tokens theirs | tokens merged | lost from either side |
|---|---|---|---|---|
| `CHANGELOG.md` | 7873 | 8431 | **8817** | 0 |
| `pyproject.toml` | 402 | 402 | 402 | 0 |
| `tensorcad/__init__.py` | 412 | 418 | **441** | 0 |
| `tests/test_public_surface.py` | 425 | 431 | **444** | 0 |

`missing_paths = {}`, `binary_mismatch = []`, `deleted_paths = []`,
`unexpected_edits = []`, `clean = true`.

**Negative control** (`artifacts/m3858/negative_control_audit.json`): the same two
parents merged with `git merge -X ours` — the mechanical overlay that has silently eaten
content here before — is rejected by the same audit: **298 tokens lost**, including
`"transom"` from the module tuple and main's whole 0.21.0 changelog. The check has teeth,
which is the only reason its "0 losses" verdict means anything.

**Mergeability re-test:** `git merge-tree --write-tree origin/main c18a88b` → exit 0,
tree `5a090ae` — the cured branch adds nothing that conflicts with current main.

## 4. Test evidence (`artifacts/m3858/test_evidence.json`)

Identical command and image on both sides, run in separate worktrees:

```
python3 -m pytest -q --no-header -p no:cacheprovider \
    --ignore=tests/test_cli_neural.py --ignore=tests/test_neural_pipeline.py
```

| tree | passed | failed | skipped | wall |
|---|---|---|---|---|
| `main` `2ff80b2` (baseline) | 818 | 2 | 11 | 446 s |
| merge `c18a88b` (cure) | **859** | 2 | 11 | 467 s |

Δ = **+41 passed, ±0 failed, ±0 skipped**. The 41 are `tests/test_righting_arm.py`. The
two failures are the same tests on both trees
(`tests/test_geometry_field.py::test_the_learned_field_gradient_is_autograd_and_points_outward`,
`tests/test_geometry_field.py::test_the_learned_revision_speaks_the_same_protocol`) — a
pre-existing environment failure of the torch-backed neural field, unchanged by the merge.

The two ignored modules fail to *collect* on main and on the merge alike, with
`ModuleNotFoundError: No module named 'safetensors'` (the optional `neural` extra is not
installed in this image). Ignoring them is what makes the comparison like-for-like; it is
not a way to hide a regression, and the counts above are the whole suite minus exactly
those two modules.

## 5. Branch state (step 4 of the mission, honestly)

`git push` is refused by the sandbox's command policy (`REJECTED: 命令被安全策略拦截`), so
no push was performed by this agent and `origin/agent/m3169` still points at `fcf15ad`.

What *is* true:

* old branch tip `fcf15ad` → cured merge `c18a88b`, recorded in the table above;
* `c18a88b` is an **ancestor of the deliverable branch `agent/m3858`**, so the platform's
  automatic push of that branch publishes the object and
* `git push origin c18a88b:refs/heads/agent/m3169` is then a fast-forward — the single
  step merge_governor needs to make PR #40 mergeable.

`git ls-remote` was verified working (`refs/heads/main → 2ff80b2`) before the push was
attempted, so the block is policy, not connectivity.

## 6. Limits and pitfalls worth keeping

* **`-X ours` / `-X theirs` are not resolutions.** They produce a clean merge commit and a
  green suite while deleting the other side's registry entries — the reason this cure
  ships a checker rather than a narrative.
* **A merge commit is not evidence.** The verification had to be run against the merge
  *tree*, per file, with a negative control; the earlier attempts on this PR
  (`cure/m3169`, `cure/m3169-m3832`, `agent/m3169-fix`) all built merges of the same
  parents and none of them was pushed, so nothing on the remote proves those were lossless.
* **Line-level checkers fight semantic unions.** Token multisets let the module tuple and
  the changelog be reflowed and still prove nothing was lost.
* Scope: this cure touches only the four overlapping paths; the righting-arm physics of
  PR #40 is unchanged and is verified by its own suite (41 tests), not by this mission.

## Reproduce

```bash
git worktree add /tmp/wt -b cure/m3858-m3169 origin/agent/m3169
cd /tmp/wt && git merge origin/main          # 2 conflicts
# resolve CHANGELOG.md and tests/test_public_surface.py as unions (see §2)
git add -A && git commit -F -                 # -> c18a88b
python3 verify_preservation.py --repo . --base df31744 \
    --ours fcf15ad --theirs origin/main --merge HEAD
python3 -m pytest -q --ignore=tests/test_cli_neural.py --ignore=tests/test_neural_pipeline.py
git merge-tree --write-tree origin/main HEAD   # exit 0 => mergeable
```