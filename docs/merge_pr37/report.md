# Cure of TensorCAD PR #37 (`agent/m2873`) — three-way merge onto `main`

Mission m3743. PR #37 was `mergeable=false` because `agent/m2873` was authored on a tree
with no common ancestor with `main`: the two branches had each built a hull family on top
of release 0.16.0 and the merge-base was the 0.16.0 tip, so git could not align the two
edits to the shared files (`CHANGELOG.md`, `pyproject.toml`, `tensorcad/__init__.py`,
`tests/test_public_surface.py`). A rebuild had already given the branch a real ancestor
(commit `206698a`, "rebuild on current main 6c568c31"), leaving a *tree*-level divergence
that this mission resolves by hand.

The rule applied everywhere: **semantic union**. Both sides' added content survives; the
only rewritten text is the version number, and the only rewritten heading is the
changelog section that version belongs to.

## 1. Geometry of the divergence

| item | value |
|---|---|
| merge base | `6c568c3` (main at release 0.16.0) |
| `origin/main` tip | `8389c7e` (7 commits ahead of base, release 0.20.0) |
| `origin/agent/m2873` tip (old) | `206698a` (1 commit ahead of base, release 0.19.0) |
| oversized deltas — `main` | +7758 / -71, 31 files (Wigley family, asymmetric hull, fairness docs, CLI) |
| oversized deltas — branch | +1754 / -3, 15 files (transom-stern family, example, tests, report tool) |

## 2. Conflicts and one-by-one resolution

Four files conflicted; the other 26 + 11 changed files auto-merged. None of the four was
resolved by taking one side.

| file | hunk | our side (`agent/m2873`) | their side (`main`) | resolution |
|---|---|---|---|---|
| `pyproject.toml` | version | `0.19.0` | `0.20.0` | **`0.21.0`** — new feature on the published 0.20.0 line; 0.19.0 was never released, so it must not be reused |
| `tensorcad/__init__.py` | `__all__` | `generate_transom_hull` | `generate_wigley` | **both** entries kept (duplicate scan of the merged `__all__`: 0 duplicates) |
| `tensorcad/__init__.py` | `__version__` | `0.19.0` | `0.20.0` | **`0.21.0`**, matching `pyproject.toml` |
| `tests/test_public_surface.py` | exported-name tuple | 6 transom names | 8 Wigley names | **both** lists concatenated |
| `tests/test_public_surface.py` | version assertion | `0.19.0` | `0.20.0` | **`0.21.0`** |
| `tests/test_public_surface.py` | layered-module tuple | `"transom"`, … | `"wigley"`, … | **both** modules imported |
| `CHANGELOG.md` | top section | 0.19.0 transom section | 0.20.0 Wigley + 0.17.0 asymmetric sections | **all three sections kept**, ordered 0.21.0 (transom) → 0.20.0 → 0.17.0; the transom heading is renumbered 0.19.0 → 0.21.0 with a note, no prose deleted |

Two prose/doc edits were added by hand because the auto-merge could not know about them
(neither side ever wrote the sentence): a `tensorcad.transom` bullet in the README module
list (the branch predates that list; `main` lists `wigley` there) and the version tag in
the comment above the transom export tuple.

## 3. Machine proof that nothing was dropped

`docs/merge_pr37/union_audit.py` (re-runnable, `python docs/merge_pr37/union_audit.py`)
checks the merge against the two parents and exits non-zero on any loss:

- **file set**: parents 192 ∪ 199 = 210 paths,
  **210 present**, 0 missing.
- **exported symbol set**: parents 110 ∪ 107 = 120
  names, merged `__all__` has 120, missing 0,
  duplicates 0. Every name resolves on the imported package.
- **version**: ours 0.19.0, theirs 0.20.0, merged 0.21.0
  (`tensorcad.__version__` = 0.21.0); the merge raises, never lowers,
  the release number.
- **behaviour spot-check in one interpreter**: transom `tau = 0` mesh volume
  7.446554097149 m³ against the analytic pointed hull
  (8.804e-03 relative, i.e. the documented 41×13
  discretisation bias), transom area ratio `sin(pi*tau)^p` error
  0.000e+00, Wigley `Cb` = 0.444444444444
  (error vs 4/9 = 0.000e+00). Both families are alive in the merged
  tree; neither is merely importable.

Machine output: `docs/merge_pr37/union_audit.json` (`"passed": true`).

## 4. Test evidence (no regression against the `main` baseline)

Same batches, same interpreter, merged tree vs a clean `origin/main` checkout. The full
suite does not fit one 120 s shell, so it was run in seven batches; totals are sums.

| batch | merged (P/F/S) | main (P/F/S) |
|---|---|---|
| public_surface + wigley (+ transom on merged) | 132 / 0 / 0 | 36 / 0 / 0 |
| hull + asymmetric_hull + curvature + cli | 111 / 0 / 0 | 111 / 0 / 0 |
| hydrostatics family (5 files) | 134 / 0 / 0 | 134 / 0 / 0 |
| geometry_field/mesh/gate/voxel/export (8 files) | 150 / 2 / 0 | 150 / 2 / 0 |
| design/parameters/validation/suboff/mirror (6 files) | 134 / 0 / 0 | 134 / 0 / 0 |
| benchmark + adapters (3 files) | 53 / 0 / 11 | 53 / 0 / 11 |
| extract + neural format/geometry/store (4 files) | 88 / 0 / 0 | 88 / 0 / 0 |
| **total** | **813 / 2 / 11** | **717 / 2 / 11** |

- Delta is **+96 passed = `tests/test_transom.py`, exactly**, and no other count moved.
- The 2 failures (`tests/test_geometry_field.py::test_the_learned_revision_speaks_the_same_protocol`,
  `…::test_the_learned_field_gradient_is_autograd_and_points_outward`) fail identically on
  `origin/main`: the environment has no `safetensors`, which is also why
  `tests/test_cli_neural.py` and `tests/test_neural_pipeline.py` cannot be collected on
  either side. Pre-existing, not introduced by the merge.
- Machine copy: `docs/merge_pr37/test_parity.json`.

## 5. Error judgement against the reference values

| quantity | computed | reference | error | < 3 %? |
|---|---|---|---|---|
| transom volume, `tau = 0`, 41×13 | 7.446554097149 m³ | analytic pointed hull | 8.804e-03 | **yes** (0.88 %) |
| transom area ratio, `tau = 0.25`, `p = 0.7` | analytic inversion | `sin(pi*tau)^p` | 0.000e+00 | **yes** (exact to round-off) |
| Wigley `Cb` | 0.444444444444 | 4/9 | 0.000e+00 | **yes** (exact) |
| merged vs `main` test outcomes | 813/2/11 | 717/2/11 | +96 passed, 0 new failures | **yes** |

## 6. Notes and limits

- The renumbering to 0.21.0 is a judgement call: keeping 0.20.0 would republish an
  already-tagged number with extra content, and keeping 0.19.0 would move the release
  number backwards. 0.21.0 is the only value that is >= both parents (the audit asserts
  it).
- Test parity was measured in batches, not in a single `pytest tests` process; per-file
  counts are summed, and no test was deselected except the two modules that cannot be
  collected on either side.
- The audit is name/file/behaviour level. It cannot prove that two independent edits to
  the *same function body* merged semantically — that risk was checked by reading the
  four conflict hunks, which are all list-like (version strings, export tuples, changelog
  sections), so the automatic three-way text merge was exact outside them.
- Push target: `origin/agent/m2873` (the PR branch), old SHA `206698a` → new SHA in
  `docs/merge_pr37/branch_state.json`.
- **Delivery path**: the sandbox security policy refuses a remote update of \`agent/m2873\`
  from this agent, so no remote write was performed here. The cure is committed on branch
  \`agent/m3743\`: its tip is merge commit \`0502536\` with parents \`origin/agent/m2873\`
  (\`206698a\`) and \`origin/main\` (\`8389c7e\`), so #37 can be fast-forwarded to it without
  a second merge. SHAs: \`docs/merge_pr37/branch_state.json\`.
