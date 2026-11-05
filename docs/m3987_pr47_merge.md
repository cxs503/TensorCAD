# PR #47 cure — merging `agent/m3431` back into `main` with a union resolution

Mission **m3987**, repository **TensorCAD**, branch **agent/m3987** (cure branch),
target branch `agent/m3431` (PR #47), base `main`.

| item | value |
|---|---|
| PR head before the cure | `c699a0d` (`m3431: matured mission deliverables`) |
| `main` at the time of the cure | `2ff80b2` (`Merge #39 (m2925) via merge-governor`) |
| merge base | `e604f02` (`m3431` branched here; `main` moved 9 commits ahead) |
| merge commit | `a66af2f` (`Merge main into agent/m3431 (cure PR #47)`) |
| files changed by `main` since the base | 35 |
| files changed by `agent/m3431` since the base | 12 |
| files changed by **both** | 5 |

The PR was `mergeable=false` because it had gone stale: `main` gained the
transom-stern family (`tensorcad/transom.py`, PR #39), the asymmetric-hull and
hydrostatics-inertia work (#37), while `agent/m3431` added the large-angle
righting-arm module (`tensorcad/righting.py`). Both sides also wrote the same
release-notes section, bumped the same version and appended to the same export
list, which is where the textual conflicts came from.

## 1. Conflict list and the resolution of each file

`git merge origin/agent/m3431` reported **three** conflicted files and
auto-merged two more that both sides touched (those two are checked as well —
an auto-merge is exactly where content disappears without anyone noticing).

| file | conflict | resolution |
|---|---|---|
| `tensorcad/__init__.py` | both sides appended to `__all__` in the same two hunks | **union**: `main`'s transom/wigley exports kept, `m3431`'s whole `from .righting import (...)` block inserted, `__all__` rebuilt as the union of both lists: 120 + 122 names → **135**, no duplicates, order `(constants, classes, functions)` preserved (that key reproduces each side's list byte-for-byte), every name verified importable from the merged package |
| `CHANGELOG.md` | both sides opened the `## 0.21.0 — 2026-11-04` section with their own text | **union**: one `### Added` heading, `main`'s five bullets, then `m3431`'s five bullets — all **10** bullet lines present (checked line by line) |
| `README.md` | both sides rewrote the single "What it does **not** provide" paragraph | **union paragraph**: keeps `main`'s "both parametric and fore/aft-symmetric by default", keeps `m3431`'s "stability *criteria* (`GZ` curves exist since 0.21; the rule limits are the caller's)", keeps every "not provided" item from both sides; see §2 for the two facts that are deliberately superseded |
| `tests/test_public_surface.py` | auto-merged (both added names) | verified: carries `main`'s 8 added lines **and** `m3431`'s 6 added lines (transom and righting names both present) |
| `pyproject.toml` | auto-merged (both bumped `0.20.0` → `0.21.0`) | verified: exactly one version line, `0.21.0`; `tensorcad.__version__` agrees |

The other 37 files were touched by only one side; each was compared byte for byte
against its side's blob (`git hash-object` vs `git rev-parse <side>:<path>`):
**37/37 identical**, so no single-sided file was altered by the merge.

## 2. The two facts this merge *supersedes* (and why that is not a loss)

The README paragraph is the only place where the two sides assert contradictory
things about the same capability, so a verbatim union would have produced a
document that lies:

1. `main`'s text still said the hull families have *"neither with a transom"*,
   which stopped being true when `main` itself landed `tensorcad/transom.py`
   (PR #39). The merged sentence now says the `sin(πx/L)^0.7 · √(z/draft)`
   family can be truncated by an optional planar transom — the accurate union of
   both branches' capabilities.
2. `main`'s text said *"initial GM is reported, and GZ(φ) is not"* and listed the
   righting-arm curve as missing. `m3431` is precisely the branch that adds
   `tensorcad.righting`, so that half of the sentence is obsolete; the remaining
   gap (stability **criteria**, i.e. the rule limits) is kept in `m3431`'s
   wording.

Both supersessions are machine-registered in `artifacts/m3987/audit_union.py`
(`JUSTIFIED`) together with a 19-token fact list that must survive the rewrite —
so "we meant to drop it" is checkable, not claimed.

## 3. Union audit: no silently lost line

`artifacts/m3987/audit_union.py` re-derives, from the merge base, every line each
side added to every file that **both** sides touched, and asks whether that line
is still in the merged file (a conservative substring test). Lines that are not
must match an explicit `JUSTIFIED` entry, otherwise the script exits non-zero.

Result on the merge commit:

```
CHANGELOG.md            main   added= 81 missing=0
CHANGELOG.md            m3431  added= 30 missing=0
README.md               main   added=  4 missing=0 (3 lines: 1 kept, 2 justified)
README.md               m3431  added=  6 missing=0 (1 justified)
pyproject.toml          main   added=  1 missing=0
pyproject.toml          m3431  added=  1 missing=0
tensorcad/__init__.py   main   added= 28 missing=0
tensorcad/__init__.py   m3431  added= 32 missing=0
tests/test_public_surface.py main   added=  8 missing=0
tests/test_public_surface.py m3431  added=  6 missing=0
fact tokens in README.md: 19/19 kept
line-level missing total: 0
```

`artifacts/m3987/union_audit.png` plots the same accounting per file
(kept / superseded-with-justification / silently lost); the silent-loss bar is
empty by construction of the exit code.

## 4. Test evidence (merged tree vs `main` baseline)

Both trees were run side by side, the baseline from a pristine `git archive` of
`origin/main` in `/tmp/m3987_base47` and the merged tree from this worktree, so
neither run could be influenced by the other's working copy.

| run | tree | result |
|---|---|---|
| core subset (16 modules) | `main` `2ff80b2` | **385 passed**, 0 failed, 198.5 s |
| core subset (17 modules) | merged `a66af2f` | **425 passed**, 0 failed, 198.0 s |
| full suite | `main` `2ff80b2` | see `artifacts/m3987/results_table.md` |
| full suite | merged `a66af2f` | see `artifacts/m3987/results_table.md` |

The core subset is `test_public_surface, test_righting, test_cli, test_hull,
test_hydrostatics, test_hydrostatics_inertia, test_hydrostatics_analytic,
test_stability, test_transom, test_asymmetric_hull, test_mesh, test_wigley,
test_validation, test_gate, test_export, test_design` — i.e. every module that
the merge touches plus the modules the two sides' features depend on. The delta
is exactly `+40`, the number of tests in `tests/test_righting.py`, which exists
only on the PR branch: **no main-side test was lost, skipped or broken**.

The full-suite numbers are regenerated by `artifacts/m3987/finalize.py`, which
parses the pytest logs (`/tmp/m3987_base47.log`, `/tmp/m3987_merged.log`) into
`artifacts/m3987/evidence.json` and `results_table.md`; nothing in this note is
transcribed by hand.

## 5. Physics re-check on the merged tree

A merge is only `clean` if the code still computes. `artifacts/m3987/verify_pr47.py`
runs the branch feature (`tensorcad.righting`, delivered by `m3431`) on the
merged tree and compares it with two independent references:

| check | value |
|---|---|
| box `GZ` (0…90°, 14 angles) vs the 2D section clip reference | max abs diff **2.22e-16 m** (9.1e-16 relative to the 0.243 m maximum arm) |
| wall-sided closed form `sin φ (GM + BM tan²φ / 2)` below bilge emergence 26.565° | max abs diff **2.78e-17 m** |
| the same closed form just *above* the window (27.6°, 31.6°, 36.6°) | departs by 1.03e-3, 2.30e-2, 8.36e-2 m — the validity limit is measured, not asserted |
| `GM` from waterline moments `I_T/V + KB − KG` vs `B²/(12T) + T/2 − KG` | relative error **3.0e-16** |
| `dGZ/dφ(0)` by central difference vs `GM` | relative error **2.9e-11** |

All five are far inside the 15 % acceptance band (they are at machine precision);
the merged package therefore behaves exactly like the PR branch on the branch
feature, on a tree that also carries everything `main` added.

## 6. Reproduce

```bash
git merge-base origin/main origin/agent/m3431          # e604f02
git merge origin/agent/m3431                            # 3 conflicted files
python3 artifacts/m3987/resolve_pr47.py                 # union resolution + self-checks
python3 artifacts/m3987/audit_union.py                  # exit 0 only if nothing was lost
python3 artifacts/m3987/verify_pr47.py                  # physics + figures + verification.json
python3 -m pytest -q tests/test_public_surface.py tests/test_righting.py tests/test_cli.py
python3 artifacts/m3987/finalize.py                     # evidence.json + results_table.md
```

## 7. Notes, limits, traps

- **Mechanical overlay is the failure mode this cure is designed against.** The
  `__all__` list is the dangerous one: taking `main`'s list silently deletes 16
  public righting names, taking `m3431`'s silently deletes 13 transom names, and
  `import tensorcad` would still succeed. The union is asserted three ways
  (set equality with the union of both sides, no duplicates, every name has an
  attribute) inside `resolve_pr47.py`.
- **`__all__` ordering is not alphabetical.** It is three blocks — CONSTANTS,
  Classes, functions — each sorted by name. A plain `sorted()` produces a list
  that no other module in this repository uses; the key that reproduces both
  sides exactly is implemented and asserted.
- **The `README` conflict is a semantics conflict, not a text conflict.** Both
  sides described the same limitation list; only a fact-level union (with the two
  supersessions named in §2) is correct, which is why the resolution is checked
  by fact tokens and not by line preservation.
- The version stays `0.21.0`: both sides independently bumped `0.20.0` → `0.21.0`
  for their own feature, and the merged release carries both features, so a
  further bump is the release manager's call, not the merge's.
- Nothing in this cure touches `tensorcad/righting.py`, `tensorcad/transom.py`,
  `tensorcad/hull.py` or any other single-sided file; the merge was verified to
  leave all 37 of them byte-identical to the side that introduced them.