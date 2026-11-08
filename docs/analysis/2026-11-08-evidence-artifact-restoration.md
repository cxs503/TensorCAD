# Red on `main`: restoring the `box_gz_area` evidence chain, and what it proves

**Date** 2026-11-08 · **Mission** m6416 · **Scope** evidence plumbing, no physics changed
**Module** `tensorcad.box_gz` (areas under the exact multi-regime righting-arm curve)
**Baseline** `main` at `31948ca`

## 1. The case: two tests were red on a clean checkout, and it was known

`tests/test_box_gz_area.py` reads a committed evidence artifact and pins the numbers of
`docs/box_gz_areas.md` to it. On the checkout this mission started from, two of its tests fail:

```
FAILED tests/test_box_gz_area.py::test_the_committed_evidence_still_makes_its_claims
FAILED tests/test_box_gz_area.py::test_the_artifact_is_what_the_code_computes_today
FileNotFoundError: .../artifacts/m4873/box_gz_areas.json
```

This was not a surprise failure, it was **recorded debt**:
`docs/m5283_geometry_audit.md` (lines 133-139) reports it as *"pre-existing and unrelated"* --
the artifact "has never been committed (`git log -- artifacts/m4873` is empty), so it fails on a
clean checkout of the parent commit as well".

Why it happened. The feature was authored in commit `7747fd8`
(*box_gz_area: the IS Code areas under the exact curve, past the wall-sided limit*), which
carries **9 files**: `tensorcad/box_gz.py`, `tests/test_box_gz_area.py`, `docs/box_gz_areas.md`,
`tools/verify_box_gz_areas.py`, `CHANGELOG.md`, `tensorcad/__init__.py`,
`test_public_surface.py` -- **and** `artifacts/m4873/box_gz_areas.json` (+`.png`).

That commit is **not an ancestor of `main`** (`git merge-base --is-ancestor 7747fd8 HEAD` is
false). `main` received the node through merge `49df625` (#235), whose commit message ends with
*"N 个文件被门禁滤除"* -- the integration gate whitelists the code/doc/test paths and drops
`artifacts/`. So `main` got the code, the page and the tests, but not the evidence those tests
open. The code was correct; the *chain* was severed.

That is the generalisable defect: **an evidence directory outside the gate's whitelist is not
an evidence directory.** A feature can land green in review and permanently red on `main`, with
its page numbers unverifiable, because the harness output it depends on is filtered in transit.

## 2. What was done

The fix is plumbing plus a regeneration, deliberately not a physics change:

1. **Evidence moved under `docs/`**, following the convention the package already uses for
   `docs/bonjean/`, `docs/fairness/`, `docs/intact_stability/`: `artifacts/m4873/` ->
   `docs/box_gz_areas/`. `docs/` survives the gate, and the page's figure link becomes
   docs-relative (`box_gz_areas/box_gz_areas.png`) instead of climbing out to `artifacts/`.
2. `tools/verify_box_gz_areas.py`: `ARTIFACTS = ROOT / "docs" / "box_gz_areas"` (+ its two
   docstring/usage references).
3. `tests/test_box_gz_area.py`: `ARTIFACT` path follows.
4. `docs/box_gz_areas.md`: the image link follows.
5. Regenerated with the committed harness, not copied from the branch:

```bash
python3 tools/verify_box_gz_areas.py     # writes docs/box_gz_areas/{json,png}, syncs the page
```

## 3. Results

### 3.1 The harness' own claims, re-measured on `main`'s code

`python3 tools/verify_box_gz_areas.py`, one core, **32.5 s**:

| claim | measured | threshold | verdict |
|---|---|---|---|
| inside the wall-sided window the integral is the closed-form antiderivative | 1.55e-13 relative | <= 1e-11 | PASS |
| statutory's adaptive Simpson agrees past the wall-sided limit | 3.02e-12 relative | <= 1e-11 | PASS |
| the heeled-mesh curve agrees to engineering tolerance | 8.11e-06 relative | <= 1e-03 | PASS |
| the negative control breaks the threshold it is compared against | 1.53e-01 | >= 1e-02 | PASS |
| a millimetre of `KG` is visible in the area | 1.04e-03 relative | non-zero | PASS |
| the split is what makes the quadrature converge | unsplit/split error 6.6e6x, 1.5e4x, 2.2e5x | > 1 | PASS |

### 3.2 The page reproduces, and the artifact is bit-identical

Two independent statements, both machine-checked:

* **The page.** The harness rewrites the tables of `docs/box_gz_areas.md` from the numbers it
  just computed. After the run, `git diff docs/box_gz_areas.md` contains **one line**: the figure
  link. Every regenerated table matched the committed one byte for byte, so the 12-decimal
  values the page prints are what today's code computes.
* **The artifact.** The regenerated files are byte-identical to the ones produced when the node
  was authored:

| file | sha256 (regenerated 2026-11-08) | sha256 (`7747fd8`, original) | match |
|---|---|---|---|
| `box_gz_areas.json` | `9fa5acf0bd8431a68e6827fab5eda258213f3dde5c312d8c145a0e6fedbffa4e` | `9fa5acf0...ffa4e` | **identical** |
| `box_gz_areas.png` | `351a91ffab87b49d3f42b5f9d6c8ab705a8a7b2e1f095bf30847f6896492038c` | `351a91ff...2038c` | **identical** |

A structural comparison of the two JSONs, leaf by leaf: **1415 leaf keys, 1141 numeric fields,
0 differing**. The figure's PNG is byte-identical too, so the plot is deterministic down to
`matplotlib` output, not merely visually similar. Reproducibility across machines and days is
what makes "the artifact *is* the evidence" a defensible claim rather than a hope.

### 3.3 Test status

| selection | before | after |
|---|---|---|
| `test_box_gz_area.py` + `test_box_gz_exact.py` + `test_public_surface.py` + `test_statutory.py` | 2 failed, 192 passed, 2 skipped | **174 passed, 2 skipped** |
| GZ/statutory evidence chain: the above + `test_gz.py`, `test_righting_arm.py`, `test_intact_stability.py` | (2 failed) | **293 passed, 2 skipped in 29.8 s** |

## 4. Error judgement against the standards used in this package

| quantity | reference | error | `< 3 %`? |
|---|---|---|---|
| area under the wall-sided window | closed-form antiderivative `GM0(1-cos) + BM0/2(sec+cos-2)` | 1.55e-13 | yes (1.6e-11 %) |
| regulated areas past the wall-sided limit | `statutory.area_under_gz_curve` adaptive Simpson | 3.02e-12 | yes (3.0e-10 %) |
| regulated areas, heeled-mesh route | `gz.righting_arm_curve` on a closed-box mesh | 8.11e-06 | yes (8.1e-4 %) |
| regenerated artifact vs authored artifact | `7747fd8` copy | 0 (byte-identical) | yes |
| page tables vs regenerated tables | harness output | 0 (no text diff) | yes |

## 5. Notes, limits, and what this does not cover

* **No physics was touched.** `tensorcad/box_gz.py`, `statutory.py` and `gz.py` are unchanged;
  the diff is two paths, one link and the regenerated files. The numbers moved by zero.
* **`artifacts/` is still the wrong home for evidence.** This page fixes one node. Any other
  page whose only copy of its numbers lives under `artifacts/` is one gate filter away from the
  same failure -- `docs/box_large_heel.md` (2 figure links into `../artifacts/m4864/`) and
  `docs/planing_hull.md` (`artifacts/planing_hull/`) are visible in the same tree, and their
  files are absent from the checkout: their images render broken on GitHub today. They are
  *not* fixed here (their harnesses were not re-run in this mission) and are listed as debt, not
  as done. The durable rule: **harness output a page cites goes under `docs/<module>/`.**
* **The full suite was not re-run in this mission** (tool wall clock is capped well below the
  ~930 s the suite took at `docs/m5283`), so the claim made here is "293 tests covering the
  entire touched evidence chain are green", not "the whole suite is green". The two tests this
  mission removes were the only failures that audit recorded.
* **The gate itself is outside this repository.** The filter that dropped `artifacts/` lives in
  the integration tooling, and this page does not change it; it removes this node's dependence
  on it.

## 6. Reproduce

```bash
python3 tools/verify_box_gz_areas.py                  # 32 s, writes docs/box_gz_areas/, syncs the page
git diff --stat docs/box_gz_areas.md                  # expect: 1 file changed (the figure link only)
python -m pytest tests/test_box_gz_area.py -q         # expect: all green
sha256sum docs/box_gz_areas/box_gz_areas.json         # 9fa5acf0...ffa4e
```