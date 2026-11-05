# PR #45 cure — merging `agent/m3391` (large-angle GZ) onto `main` without losing a side

Mission m3781 · merge commit `27f7b56` · branch `agent/m3781` (pushed to `agent/m3391`)

## 1. Why the PR was dirty, and what "dirty" actually meant here

PR #45 adds `tensorcad/gz.py` (the righting-arm curve: `heeled_hydrostatics`,
`righting_arm_curve`, `with_freeboard`, free-surface correction) plus `docs/gz.md`,
`examples/gz/` and `tests/test_gz.py` — 2 311 added lines on 11 files. The branch forked
at `ca0dde3`; while it was open, `main` gained three whole capability sets:
waterplane second moments on the `Hydrostatics` record (`I_T`, `I_L`, `keel_z`), the
transom-stern hull family (`tensorcad/transom.py`, 650 lines), and the asymmetric
fore/aft hull exponents. Both sides therefore touched `CHANGELOG.md`, `pyproject.toml`,
`tensorcad/__init__.py` and `tests/test_public_surface.py`, and GitHub reported
`mergeable=false`.

Two preceding cures in this repo (PR #39, PR #40) died on *silent content loss* from
mechanical overlay merges, so this cure is judged by measurement, not intent:
`tools/verify_merge_m3391.py` writes `artifacts/m3781/merge_verification.json` and exits
non-zero if any of the four checks below fail.

## 2. Method

Work tree `/root/tp-agent/work/m3781` on `agent/m3781`, base `origin/main` = `2ff80b2`
(Merge #39), merged side `agent/m3391` = `164fc63`.

1. Three-way merge `git merge --no-ff agent/m3391`; one textual conflict, in
   `CHANGELOG.md`; `tensorcad/__init__.py` and `tests/test_public_surface.py` auto-merged.
2. Conflict resolved by **semantic union**, not by choosing a side: one `## 0.21.0`
   heading with two `### Added` blocks — main's waterplane/transom-stern block verbatim,
   then the branch's large-angle-stability block verbatim. The branch's own `h2` title was
   the single line whose level changed (`##` → `###`); its text is intact.
3. Auto-merged files were **audited, not trusted**: every identifier, import and declared
   public name of each parent must still be present on the merged tree.
4. Core tests run per file on the merge **and** on a clean `origin/main` work tree, so
   "not below baseline" is a measured comparison against the same machine load.

## 3. Conflict inventory and per-file resolution

| File | Kind | Resolution | Verification |
|---|---|---|---|
| `CHANGELOG.md` | textual conflict (both sides opened `0.21.0`) | one `0.21.0` heading + two `### Added` blocks, both bodies verbatim | 780 non-empty lines from `main` present, 726 from the branch present, 0 missing; 1 heading remapped (`##`→`###`), 0 unexplained |
| `tensorcad/__init__.py` | auto-merged | union of both import blocks: 14 `gz` names **and** 12 `transom` names | 242 identifiers from `main` + 244 from the branch, 0 missing |
| `tests/test_public_surface.py` | auto-merged | `EXPECTED_API` superset; `transom` kept in the layered-module list, `gz` names added | 179 identifiers per parent, 0 missing |
| `pyproject.toml` | both touched, same value (`0.21.0`) | no change needed | `version = "0.21.0"` on both parents and on the merge |

File-level union (no path dropped, no path invented): `main` 217 files, branch 200 files,
merge 224 files, `missing_from_main = []`, `missing_from_branch = []`, `invented = []`.

Public surface: `__all__` has 120 names on `main`, 121 on the branch, **134 on the
merge** — 107 shared plus 27 split names (`MAX_TRANSOM_FRACTION` … `with_freeboard`),
`main_only_lost = []`, `branch_only_lost = []`.

## 4. Results — both routes still compute (numbers, not claims)

Branch route, box barge 10 × 2 × 2 m, T = 0.8 m, KG = 0.5 m, wall-sided limit
38.6598°:

| heel [deg] | 0 | 5 | 10 | 20 | 30 | 35 |
|---|---|---|---|---|---|---|
| GZ mesh [m] | −0.000000 | 0.027738 | 0.056113 | 0.117746 | 0.193056 | 0.240220 |
| GZ closed form [m] | 0.000000 | 0.027738 | 0.056113 | 0.117746 | 0.193056 | 0.240220 |

`max |GZ_mesh − GZ_wall-sided| = 6.94 × 10⁻¹⁷ m`, max relative volume residual
`2.22 × 10⁻¹⁶` (displacement held at 16.000 m³ through the sweep), `max GZ = 0.240220 m`,
dynamical stability `0.0657612 m·rad` over 0–35°.

Main's routes on the same merged tree:

| check | 41 × 13 | 81 × 25 | reference |
|---|---|---|---|
| `I_T` rel. error vs closed form | −0.1386 % | −0.0345 % | `hydrostatics_analytic` |
| `I_L` rel. error vs closed form | −0.6038 % | −0.1979 % | tanh-sinh quadrature |
| transom (τ = 0.2, p = 1.0) | Cb = 0.479857, ∇ = 7.677707 m³, A_t,imm = 0.626971 m² | | sinc^p power series |
| Wigley reference | Cb = 0.444444 (= 4/9), ∇ = 2.777778 m³, Cm = 0.666667 | | closed forms |

## 5. Test evidence (per-file runs, same machine, merge vs clean `main`)

| tree | passed | failed | skipped |
|---|---|---|---|
| baseline `origin/main` (`2ff80b2`) | 722 | 2 | 11 |
| merged (`27f7b56`) | **748** | 2 | 11 |

748 − 722 = 26, exactly the new `tests/test_gz.py` file. **No file is below its
baseline** (`below_baseline = []`). The 2 failures are identical on both trees and are an
environment gap, not a merge effect: `tensorcad/neural/field.py` imports
`safetensors`, which is absent in this container, so
`test_geometry_field.py::test_the_learned_revision_speaks_the_same_protocol` and
`…::test_the_learned_field_gradient_is_autograd_and_points_outward` fail with
`ModuleNotFoundError: No module named 'safetensors'` on `main` as well. The two modules
`tests/test_cli_neural.py` and `tests/test_neural_pipeline.py` fail to *collect* for the
same reason on both trees and are excluded from the table.

Raw per-file logs: `artifacts/m3781/pytest_merged_per_file.log`,
`artifacts/m3781/pytest_baseline_per_file.log`, totals in
`artifacts/m3781/tests_per_file.json`. `tests/test_merge_cure_m3391.py` re-runs the audit
inside pytest (3 tests) so the next merge cannot quietly undo it.

## 6. Error judgement against the 3 % bar

| quantity | reference | merged-tree result | error | < 3 %? |
|---|---|---|---|---|
| GZ (8 heels 0–35°) | wall-sided closed form (derived in `gz.py`) | 0.240220 m at 35° | 6.9e-17 m abs (≈ 3e-16 rel) | **yes** |
| displacement across the sweep | upright volume 16.000 m³ | 16.000000 m³ | 2.2e-16 | **yes** |
| `I_T` / `I_L` on 41 × 13 | closed forms | — | 0.14 % / 0.60 % | **yes** |
| `I_T` / `I_L` on 81 × 25 | closed forms | — | 0.035 % / 0.20 % | **yes** |
| test parity | 722 baseline passes | 748 | ≥ 0 per file | **yes** |

No quantity required a reference that is missing, so no claim is withheld.

## 7. Caveats

- The GZ agreement is exact because a box barge's submerged section is a polygon the mesh
  represents exactly; it is a check of the *integration and the waterline solve*, not of
  discretisation error on a curved hull. Beyond 38.66° the closed form stops applying by
  construction (the waterline reaches the bilge), which is why the comparison stops at 35°.
- `safetensors` is missing in this image; the two neural tests fail identically on `main`
  and on the merge. Installing it (or adding `pytest.importorskip("safetensors")` in
  `tensorcad/neural/field.py`) would take the suite to 748 + 2 green — a separate fix, not
  part of this cure.
- Push: the platform pushes this work tree's branch to `agent/m3391`
  (`164fc63` → `27f7b56`); `merge_governor` re-checks automatically next cycle.