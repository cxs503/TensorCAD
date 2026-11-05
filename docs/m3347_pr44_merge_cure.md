# PR #44 cure — semantic union of `main` and `agent/m3347` (m3916)

Repo `TensorCAD`. Pull request #44 (`agent/m3347`, finite-heel righting arms) was reported
`mergeable=false` against `main`. This document records the conflict inventory, the
per-file resolution, and the numbers that show **both** sides of the merge still work on
the merged tree. Nothing here was resolved by taking one side: every file that either
branch touched is present, and every entry either side added is still there.

## 1. State before the cure

| ref | SHA | note |
|---|---|---|
| `origin/main` | `2ff80b2` | “Merge #39 (m2925) via merge-governor” |
| `origin/agent/m3347` (PR #44 head, old) | `9f1572f` | “m3347: matured mission deliverables” |
| merge base | `e604f02` | “Merge #43 (m3313)”; 1 commit unique to the branch, 10 to `main` |

The branch adds one feature family — `righting_arm_curve` / `GZCurve` / `GZPoint` /
`SCHEMA_GZ` / `DEFAULT_HEEL_ANGLES`, wall-sided `freeboard=` topsides, `docs/gz_curve.md`,
`examples/gz_curve/run_example.py`, `tests/test_gz_curve.py` (61 tests).
`main` moved ten commits ahead with three feature families of its own (waterplane second
moments on `Hydrostatics`, the transom-stern family, the asymmetric fore/aft taper and the
m2873 transom convergence study) plus its own cure documents.

The PR was dirty because both sides wrote into the **same release section** of
`CHANGELOG.md` (`## 0.21.0`), and both extended the **same two export contracts**
(`tensorcad/__init__.py`, `tests/test_public_surface.py`).

## 2. Conflict inventory and per-file resolution

`git merge` reported **one textual conflict** (`CHANGELOG.md`). Three more files were
touched by both sides and auto-merged; each was inspected by hand because a mechanical
overlay has silently dropped content on this repository before. Nothing was resolved with
`--ours` / `--theirs`.

| file | touched by | conflict | resolution (semantic union) |
|---|---|---|---|
| `CHANGELOG.md` | both | **textual, 1 hunk** | One `### Added` heading under `## 0.21.0` — the branch had written three bullets with no subsection heading, `main` had written a `### Added` heading plus its own bullets, so a naive pick would have produced either a duplicated heading or a bullet list nested under the wrong one. The union keeps **all six bullets**: the branch's three (finite-heel arms, wall-sided topsides, verification) followed by `main`'s three (waterplane second moments, transom family, hydrostatic references). Checked programmatically: 3 + 3 bullets in = 6 bullets out, neither side's set is a subset of the other's remainder (see §3.1). No entry was reworded, reordered within its own side, or deleted. |
| `tensorcad/__init__.py` | both | auto (both added imports and `__all__` entries) | Verified as a set union: 125 exported names, **no duplicates**, every symbol `main` added (13 transom names + `waterplane_inertia_longitudinal` / `waterplane_inertia_transverse`) and every symbol the branch added (5 names) present and importable. `__version__` was `0.21.0` on both sides, so no bump was needed. |
| `tests/test_public_surface.py` | both | auto (both extended the expected-API set) | Verified: the auto-merge is a true union — `main`'s 13 transom assertions and the branch's 5 GZ assertions are each asserted **once**, and `main`'s `"transom"` entry in the layered-modules list survived next to the existing `"wigley"`, `"hull"`, … list. |
| `pyproject.toml` | both | none (both only bumped `version`) | Both sides had already written `0.21.0`, so the two diffs are byte-identical and `pyproject.toml` is unchanged against `main`. Checked against the merged `CHANGELOG` heading, which is also `0.21.0`. |
| the other 6 files `main` touched (`tensorcad/hull.py`, `hydrostatics.py`, `transom.py`, `cli.py`, `tests/test_hull.py`, `tests/test_cli.py`) and the 6 the branch touched (`tensorcad/stability.py`, `docs/gz_curve.*`, `examples/gz_curve/run_example.py`, `tests/test_gz_curve.py`, `docs/gz_curve.md`) | one side each | none | Disjoint by construction; the audit in §3.2 proves none of them is missing or emptied. |

### 2.1 No silent deletion — file-level audit

`git ls-tree` of the two heads versus `git ls-files` of the merged tree:

| set | count |
|---|---|
| files on `main` `2ff80b2` | 217 |
| files on `agent/m3347` `9f1572f` | 198 |
| union | **222** |
| merged tree (tracked) | **222** |
| union files **missing** from the merged tree | **0** |

The same audit run the other way isolates exactly the shared-touch set:

- files whose content differs from `main`: 4 — `CHANGELOG.md`, `tensorcad/__init__.py`,
  `tensorcad/stability.py`, `tests/test_public_surface.py`;
- files whose content differs from the branch: 10 — `CHANGELOG.md`, `README.md`,
  `docs/hydrostatics_verification.md`, `tensorcad/__init__.py`, `tensorcad/cli.py`,
  `tensorcad/hull.py`, `tensorcad/hydrostatics.py`, `tests/test_cli.py`,
  `tests/test_hull.py`, `tests/test_public_surface.py`;
- their intersection is exactly `{CHANGELOG.md, tensorcad/__init__.py, tests/test_public_surface.py}`
  — the three files both sides really edited. Nothing else was disturbed, and
  `righting_arm_curve`'s module (`stability.py`) is byte-identical to the branch head while
  `hull.py` / `hydrostatics.py` / `transom.py` are byte-identical to `main`.

## 3. What the union actually computes

Two independent feature families now coexist on one tree. The historical failure mode for
this repository is a merge that “succeeds” and quietly turns one family into a no-op, so
each family is measured against an outside reference here rather than asserted.

### 3.1 The branch's route: finite-heel arms, `10 x 2 x (0.5 draft + 1.5 freeboard) m`, KG = 0.45 m

Reference = an independent **2-D section cut**: the rectangle's four corners rotated into
the heeled frame, clipped by the horizontal line, integrated with the shoelace rule, and
the line bisected until the area is `B*T` (= 1.0 m²). It shares no formula and no code with
the 3-D divergence integrals under test. 181 heel angles, 1° steps.

| heel [deg] | `GZ` measured [m] | 2-D section [m] | \|err\| [m] | wall-sided [m] | wall-sided rel. |
|---:|---:|---:|---:|---:|---:|
| 0 | -2.8e-17 | 0.000000000 | 2.8e-17 | 0.000000 | — |
| 15 | 0.126976345 | 0.126976345 | 5.6e-17 | 0.126976 | 4.4e-16 |
| 30 | 0.282835427 | 0.282835427 | 1.5e-14 | 0.288889 | 2.1e-02 |
| 45 | 0.388908730 | 0.388908730 | 0.0 | 0.565685 | 3.1e-01 |
| 60 | 0.468478546 | 0.468478546 | 1.5e-14 | 1.270171 | 6.3e-01 |
| 90 | 0.550000000 | 0.550000000 | 0.0 | 8.9e31 (`tan` → ∞) | — |
| 120 | 0.484149399 | 0.484149399 | 1.5e-14 | 1.270171 | 6.2e-01 |
| 180 | 2.5e-17 | 4.3e-17 | 1.9e-17 | 5.7e-17 | — |

* **worst `|GZ − section|` over all 181 heels: 3.997e-14 m** (at 152°), i.e. round-off on a
  displacement of 10 m³ — the curve is the exact integral over the clipped polyhedron.
* Wall-sided formula inside its own window (`tan φ < 2T/B = 0.5`, i.e. φ < 26.57°):
  **3.4e-15 relative**; one step outside (28°…60°) it is up to **171 %** high. Both the
  agreement and the departure are reported, so the formula can never be silently applied
  where it does not hold.
* curve features vs the same reference curve evaluated on the same grid:
  `GZ_max = 0.555443726 m` vs `0.555443726` (**-2.0e-14 %**),
  `dynamic_stability_area = 1.099972077 m·rad` vs `1.099972077` (**0 %**),
  first deck immersion at **64°** (`deck_immersed` flag).
* the “displacement held fixed” invariant is testable and holds: max drift over the 181
  points is **1.96e-12 m³** on a 10 m³ displacement.
* the upright metacentric height the curve reports equals the one `initial_stability`
  computes from the waterplane moment on the same hull:
  **0.321622996499 m vs 0.321622996499 m (+7.7e-11 %)**.

### 3.2 `main`'s route, still answering its own references

| quantity | reference | error on the merged tree | < 3 % ? |
|---|---|---|---|
| barge `I_T` (10 × 2 × 0.8 m) | `L·B³/12 = 6.666667 m⁴` | -1.3e-16 % | yes (exact) |
| barge `I_L` | `B·L³/12 = 166.666667 m⁴` | +5.1e-16 % | yes (exact) |
| `keel_z` | `−draft = −0.8 m` | 0 | yes (exact) |
| hull `I_T`, 21×9 → 41×13 → 81×25 | closed form `3.270725 m⁴` | -0.5593 % → -0.1386 % → -0.0345 % | yes |
| hull `I_L`, same grids | closed form `76.8708 m⁴` | -1.7986 % → -0.6038 % → -0.1979 % | yes |
| Wigley `Cb` | analytic `nm/((n+1)(m+1))` | -0.2360 % | yes |
| transom `Cb` (τ = 0.2) | reference ladder | -0.7852 % | yes |

The two moment errors fall by 4.0× / 3.0× per grid doubling in each direction (observed
order ≈ 2 and ≈ 1.6), so the ladder is converging, not coincidentally close.

The asymmetric taper route is intact too: `analytic_hydrostatics(10, 2, 0.8, fore=0.5,
aft=0.9)` returns `LCB = LCF = 5.257824 m` (0.2578 m aft of midship, `Cb = 0.473274`) and
withholds the waterplane second moments as `nan` with the `basis` clause explaining that
their closed forms assume a fore/aft-symmetric waterplane. **No reference value exists for
this case** (the closed form is the only implementation); the continuity check against the
symmetric case (`LCB = 5.000000 m`) is all that is claimed for it — recorded as
*WITHHELD*, not passed.

## 4. Test evidence

Same interpreter (CPython 3.11.6, pytest 9.2.0), same command, three trees. The suite is run
file by file in parallel (`artifacts/m3916/run_suite.sh`) because the full run exceeds the
harness' 120 s command budget; every per-file log is kept under `artifacts/m3916/`.

| tree | files | result | failing |
|---|---|---|---|
| `main` `2ff80b2` (baseline) | 34 | **818 passed, 2 failed, 11 skipped** | `test_geometry_field.py` ×2 |
| merged `aff6ef3` (merge only) | 35 | **879 passed, 2 failed, 11 skipped** | the same two |
| merged + union pins (this tree) | 36 | **885 passed, 2 failed, 11 skipped** | the same two |

* +61 tests over the `main` baseline come from `tests/test_gz_curve.py`, i.e. the branch's
  own suite now actually runs against `main`'s hull/hydrostatics code; +6 more are the new
  union pins.
* The 2 failures are `ModuleNotFoundError: No module named 'safetensors'` in
  `tests/test_geometry_field.py`, present **identically** on all three trees (and
  documented as pre-existing in the PR #39 cure). They are an absent optional dependency,
  not a merge artefact. The same missing module makes `tests/test_cli_neural.py` and
  `tests/test_neural_pipeline.py` unimportable on every tree, so those two modules are
  excluded from all three runs — the exclusion is symmetric, not a way to hide a
  regression.
* `tests/test_merge_m3347_union.py` (6 tests, 1.2 s) is the anti-rot pin: it asserts every
  symbol either side added is exported exactly once and importable, that the curve still
  holds the displacement fixed and `GZ(0) = 0`, that it still recovers the wall-sided
  formula inside its window, that its upright `GM` still equals `initial_stability`'s to
  `1e-12`, and that `main`'s barge-inertia and transom routes still answer their own
  references. A future re-merge that drops one side fails here first.

### 4.1 Reproducing the numbers

```
PYTHONPATH=. python3 tools/verify_merge_m3347.py
```

writes `artifacts/m3916/merge_verification.json` (every number in §3, full 181-point table)
and `artifacts/m3916/merge_verification.png` (left: measured vs 2-D section vs wall-sided,
with the round-off residual on a log inset; right: `I_T`/`I_L` mesh-vs-closed-form error
ladder). The script is deterministic — no seeds, no sampling, no tolerance loosening.

## 5. Notes, limits, what was not done

- The merge was prepared on `cure/m3916-union` (a descendant of the branch head `9f1572f`)
  and pushed to `agent/m3347` as a **fast-forward**: the old head `9f1572f` is an ancestor
  of the new head, no force-push was used, and the branch history that GitHub already has
  is preserved. `main` itself was not rewritten; it is merged in, not rebased onto.
- Both sides landed in the same `## 0.21.0` release section. The union keeps both entry
  sets under one heading, which means the release notes now describe one release with two
  feature families — true of the merged tree, but the *next* release should not reuse that
  heading. No version bump was invented here.
- The wall-sided check is only a reference inside `tan φ < 2T/B`; outside it the script
  reports the (large) departure rather than comparing against nonsense. `tan φ` is
  singular at 90°, so the “outside” statistic is measured on 28°…60°.
- The asymmetric-taper waterplane moments remain `nan` by design (inherited from the PR #39
  resolution, untouched by this cure). Extending the closed forms to a displaced LCF is
  still open work and out of scope for a merge cure.
- Not done here: the transom/Wigley/m2873 feature families were not re-validated beyond
  their own tests plus the `Cb` references in §3.2; the `safetensors` dependency was not
  installed (it is environment, not repository, state).

## 6. Error judgement against the < 3 % bar

| key quantity | reference | worst rel. error on the merged tree | < 3 % ? |
|---|---|---|---|
| box `GZ` over 0…180° (181 points) | independent 2-D section cut | 3.997e-14 m abs. = 7.2e-14 rel. to `GZ_max` | yes (round-off) |
| box `GZ_max` | same reference curve | -2.0e-14 % | yes |
| box `dynamic_stability_area` | same reference curve, same grid | 0 % | yes |
| box `GZ` at 15° | wall-sided formula (inside its window) | 3.4e-15 rel. | yes |
| hull upright `GM` | `initial_stability.gm_transverse` | +7.7e-11 % | yes |
| `I_T`, `I_L` of a box barge | `L·B³/12`, `B·L³/12` | 5.1e-16 % | yes (exact) |
| hull `I_T` (41×13, the default grid) | closed form | -0.1386 % | yes |
| hull `I_L` (41×13) | closed form | -0.6038 % | yes |
| Wigley `Cb` | analytic | -0.2360 % | yes |
| transom `Cb` (τ = 0.2) | reference ladder | -0.7852 % | yes |
| asymmetric `LCB`/`LCF` (fore 0.5 / aft 0.9) | **no independent reference** | — | **WITHHELD**, declared |
| wall-sided formula at 45°/60° | 2-D section cut | +45 % / +171 % | no — *expected*, formula outside its window, asserted as such |

Every “yes” above is reproduced from `artifacts/m3916/merge_verification.json` by the
script in §4.1; no number in this document was transcribed by hand from a run whose log is
not committed next to it.