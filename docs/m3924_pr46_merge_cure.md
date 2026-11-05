# PR #46 cure — semantic union of `main` and `agent/m3406` (m3924)

Repo `TensorCAD`. Pull request #46 (`agent/m3406`, mission m3406: the righting-arm
`GZ` route) came out `mergeable=false` against `main` (mission m3768/m2873/m2925:
asymmetric hull exponents, transom-stern family, waterplane inertias on the
`Hydrostatics` record). This document records the state before the cure, the conflict
inventory, the per-file resolution, and the numbers that show neither side lost
anything.

## 1. State before the cure

| ref | SHA | note |
|---|---|---|
| `origin/main` | `2ff80b2` | “Merge #39 (m2925) via merge-governor” |
| `origin/agent/m3406` (PR #46 head, old) | `0900b21` | “m3406: matured mission deliverables”, 1 commit unique to the branch |
| merge base | `e604f02` | “Merge #43 (m3313) via merge-governor” |
| merged head (new `agent/m3406`) | `1401fc0` | merge commit, parents `0900b21` + `2ff80b2` |

`main` was 10 commits ahead of the base and changed **35 files** (+4250/−68); the branch
was 1 commit ahead and changed **10 files** (+2001/−13). Only **two** paths were edited
on both sides — the public surface and its test — so a real 3-way merge can keep
everything provided those two are resolved as a semantic union rather than by picking a
side.

## 2. Conflict inventory and per-file resolution

`git merge --no-ff origin/main` produced **one textual conflict**
(`tensorcad/__init__.py`, a single hunk inside `__all__`) and one silent-but-safe
auto-merge (`tests/test_public_surface.py`, two additions in disjoint regions).

| file | conflict | resolution (union, nothing dropped) |
|---|---|---|
| `tensorcad/__init__.py` | 1 hunk: both sides inserted an alphabetically adjacent export at the same anchor (`sinc_power_coefficients` vs `submerged_moments`) | kept **both**: `"sinc_power_coefficients"` then `"submerged_moments"`, matching the file's alphabetical order. All other additions (13 from `main`, 6 from the branch) had already auto-merged: `transom` family, `parent_length`, `MAX_TRANSOM_FRACTION`, `truncated_taper_*`, `transom_*`, `estimate_transom_discretisation_error` on one side; `GZCurve`, `RightingArm`, `SubmergedMoments`, `gz_curve`, `righting_arm`, `submerged_moments` on the other. `__version__` stays `0.21.0` (main's bump). |
| `tests/test_public_surface.py` | none (auto) | verified, not assumed: the merged `EXPECTED_API` is exactly `EXPECTED_API(main) ∪ EXPECTED_API(branch)` = 126 symbols, equal to `set(tensorcad.__all__)`; the branch's `"GZCurve", "RightingArm", "SubmergedMoments", "gz_curve", "righting_arm", "submerged_moments"` line and main's transom block both survive, no assertion is weakened to a subset test. |
| 33 files touched only by `main` | none | byte-identical to `origin/main` after the merge (checked file by file, incl. binaries: `docs/transom_stern.md`, `docs/asymmetric_hull*`, `artifacts/m3768/*`, `tools/*`). |
| 8 files touched only by the branch | none | byte-identical to `0900b21` after the merge: `docs/stability.md`, `docs/stability/{gz_curves,gz_convergence}.png`, `docs/stability/{gz_report.md,gz_results.json}`, `examples/stability/run_gz.py`, `tensorcad/stability.py`, `tests/test_stability_gz.py`. |

Machine check (re-runnable): `python tools/verify_merge_m3406.py --full` →
`artifacts/m3924/merge_verification.json`, `all_ok: true`.

## 3. Why no deeper semantic conflict exists

The two sides are *additive on different physics*:

* `main` changed the buoyancy side — `tensorcad/hull.py` (asymmetric fore/aft
  exponents, analytic inertias), `tensorcad/hydrostatics.py` (three new fields on the
  `Hydrostatics` dataclass, all with `nan` defaults), added `tensorcad/transom.py`.
* the branch changed the stability side — `tensorcad/stability.py` (`waterplane_moments`
  by Green's theorem, `righting_arm`, `gz_curve`, `submerged_moments`).

The branch consumes `compute_hydrostatics` from `main`'s module. That call is unchanged
by the merge — signature `compute_hydrostatics(mesh, waterline_z=0.0, *, water_density=1025.0,
closed=None)` — and the three new fields are defaulted dataclass fields, so the branch's
existing call sites keep working and simply gain the inertias. Nothing on either side
relied on a name the other side renamed: `main`'s `analytic_hydrostatics` rewrite touched
`hull.py`/`hydrostatics.py` only, which the branch does not import.

## 4. Test evidence — merged tree vs `main` baseline

Both suites run with the same command and the same exclusions (two modules that need the
optional `safetensors` extra, absent in this environment — the exclusion is applied to
both sides, so the comparison is apples to apples):

```
python -m pytest tests -q -n 8 --ignore=tests/test_cli_neural.py \
    --ignore=tests/test_neural_pipeline.py
```

| tree | passed | failed | skipped | wall |
|---|---|---|---|---|
| `main` baseline (`2ff80b2`, extracted to a clean dir) | 818 | 2 | 11 | 298 s |
| merged (`1401fc0`) | **837** | **2** | 11 | 254 s |

* **+19 passed, zero new failures**, and the +19 is exactly the 19 tests of the
  branch's new `tests/test_stability_gz.py` (collected: 19) — no existing test was
  dropped or skipped, since `skipped` stays at 11 on both trees.
* The 2 failures are identical on both trees and are **pre-existing environment
  failures**, not merge damage: `tests/test_geometry_field.py::test_the_learned_revision_speaks_the_same_protocol`
  and `::test_the_learned_field_gradient_is_autograd_and_points_outward` fail with
  `ModuleNotFoundError: No module named 'safetensors'` while importing
  `tensorcad/neural/field.py` — the same extra that forces the two ignored modules.
* Logs committed: `artifacts/m3924/pytest_main_baseline.log`,
  `artifacts/m3924/pytest_merged.log`.

## 5. Both sides still compute their own numbers on the merged tree

Tests prove the contracts; regenerating the two sides' committed artifacts proves the
*numbers*. Both examples were re-run from the merged tree and compared to the committed
baselines.

**(a) Branch route — `GZ` (mission m3406).** `examples/stability/run_gz.py --output …`
on the merged 0.21.0 tree, compared leaf by leaf against the committed
`docs/stability/gz_results.json` (written at 0.20.0):

| quantity | value |
|---|---|
| numeric leaves compared (scalars + list entries) | 385 |
| max relative difference | **0.0** (identical to the last bit) |
| barge `KG` ladder — `GM` | 0.233333 / 0.183333 / 0.133333 m (unchanged) |
| barge `KG` ladder — vanishing angle | 110.669° / 90.000° / 69.331° (unchanged) |
| wall-sided deck-edge check | `GZ` vs `sin φ (GM + ½ BM tan²φ)` still round-off |

Regenerated copy: `artifacts/m3924/gz_results_merged_0.21.0.json`. The only field that
differs is the provenance stamp `"version": "0.20.0" → "0.21.0"`; the committed branch
artifact is left untouched so the branch's own history stays linear.

**(b) Main route — transom design space (missions m2873/m2925/m3768).**
`examples/transom_stern/run_example.py --output …` on the merged tree:

| file | vs `docs/m2873/outputs/transom_stern/` |
|---|---|
| `design_table.csv` (15 design points) | byte-identical |
| `convergence.csv` (tau ladder) | byte-identical |
| `reference_ladder.csv` | byte-identical |

Topology assertions also still hold on the merged tree: watertight, edge-manifold,
`V−E+F = 2` for `tau ∈ {0, 0.05, 0.2, 0.35}`.

## 6. Error judgement against the acceptance threshold

| key quantity | reference | merged vs reference | < 3 % ? |
|---|---|---|---|
| 385 GZ/`GM`/vanishing-angle numbers | branch committed baseline (0.20.0) | 0.0 % (bit-exact) | yes |
| `design_table.csv`, `convergence.csv`, `reference_ladder.csv` | main committed baseline | 0 % (byte-identical) | yes |
| public surface | `main.__all__ ∪ branch.__all__` | 126 vs 126, 0 missing / 0 extra / 0 duplicate | yes |
| test suite | `main` baseline (818 passed) | 837 passed, same 2 environmental failures | yes (no regression) |

Every quantity in this cure is compared against a *committed* reference, so the 3 %
threshold is not the binding constraint here — the observed error is zero except for the
two pre-existing environment failures, which are declared as **honest-FAIL** for the
optional `safetensors`/neural path on both trees (this work does not fix them; they exist
on `main` and are outside PR #46's scope).

## 7. Notes, limits, applicability

* **Mechanical overlay is the failure mode this cure avoids.** Copying one side's
  `__init__.py` over the other would have silently dropped 6 (or 13) exports and made the
  other side's feature unreachable from the package root while `tests/test_stability_gz.py`
  or `tests/test_transom.py` still passed — they import submodules directly. The union is
  asserted symbol by symbol in `tools/verify_merge_m3406.py`, and the *test* file's
  `EXPECTED_API` is compared against the union of the two sides' `EXPECTED_API`, so a
  drop on the test side cannot hide a drop on the API side.
* **Binary assets** (`png`, and the `main`-side `artifacts/m3768/*`) are compared as bytes,
  which is the only meaningful check for them; the *numbers behind* the branch's two PNGs
  are covered by the JSON reproduction in §5.
* **Version stamp drift**: the branch's `docs/stability/*` were generated at 0.20.0 and
  now ship under a 0.21.0 tree. Numbers are unaffected (proven above); regenerating them
  would be a provenance refresh, deliberately left out to keep this cure minimal.
* **Scope**: this document covers the merge only. Whether the *content* of PR #46 meets
  its own physics acceptance criteria is documented in `docs/stability.md` and
  `docs/stability/gz_report.md` (unchanged by this cure).

## 8. Reproduction

```
git checkout agent/m3406            # local cure branch, head 1401fc0
python tools/verify_merge_m3406.py --full --json artifacts/m3924/merge_verification.json
python tools/report_m3924.py       # regenerates artifacts/m3924/merge_verification.png
python -m pytest tests -q -n 8 --ignore=tests/test_cli_neural.py \
    --ignore=tests/test_neural_pipeline.py
```