# Cure of PR #39 (branch `agent/m2925`) — conflict resolution record

Mission m3755.  `mergeable=false` on GitHub; three files conflicted against
`main` and the resolution below is a **semantic union**, not an overlay: both
features are live and both keep their own evidence, and the one place where the
two features contradicted each other (a symmetric-only exponent on an asymmetric
hull) was fixed numerically rather than by picking a side.

| item | value |
|---|---|
| repository | TensorCAD |
| branch under cure | `agent/m2925` |
| old branch SHA | `e7e8b18` |
| merge base | `6c568c31` |
| merged `main` | `311f45a` |
| new branch SHA | the merge commit recorded in this branch (parents `e7e8b18`, `311f45a`) |

`agent/m2925` is one commit on top of the merge base (`0.21.0` waterplane second
moments); `main` is ten commits ahead of it (transom stern `0.21.0`, Wigley
`0.20.0`, asymmetric hull `0.17.0`, plus merge-governor merges).

## Conflicted files and how each was resolved

### 1. `tensorcad/hull.py` — two hunks, both semantic

*Hunk 1 — `analytic_hydrostatics` body.*  `main` grew a split sine-taper exponent
(`analytic_longitudinal_centre(length, fore_exponent=…, aft_exponent=…)` plus a
`basis` string that depends on whether the family is symmetric), `agent/m2925`
added the two waterplane second moments and a `basis` clause saying the
longitudinal one is quadrature.  Union: `main`'s exponent-aware centre **and**
m2925's inertia calls, with the basis string composed (main's family-dependent
text, then the inertia clause appended for both families).

The union exposed a real defect that a mechanical overlay would have shipped:
m2925 called `waterplane_inertia_transverse(..., breadth_exponent=BREADTH_EXPONENT)`,
i.e. it hard-coded the *symmetric* exponent.  On `main`'s asymmetric family that
returns a wrong number with no error signal (2.0 % low for
`(p_fore, p_aft) = (0.5, 0.9)`).  Fix: the two references in
`tensorcad/hydrostatics_analytic.py` now take `fore_exponent` / `aft_exponent`,

* `I_T = (B³L/12) · ½·(T(3·e_aft) + T(3·e_fore))`, `T` the gamma identity — still
  no quadrature; equal exponents short-circuit to the original expression, so the
  symmetric family is bit-for-bit unchanged;
* `I_L` per half, `B L³ (J(e_aft) + J(e_fore))`, plus the parallel-axis shift from
  midship to the LCF (`I_LCF = I_mid − A_wp (LCF − L/2)²`), which `main`'s
  `analytic_longitudinal_centre` supplies.  This is the first time the record's
  `∫(x−LCF)² dA` claim is true on the asymmetric family rather than only by
  symmetry.

*Hunk 2 — the `Hydrostatics(...)` construction.*  Union: m2925's
`transverse_waterplane_inertia` / `longitudinal_waterplane_inertia` / `keel_z`
fields **and** main's `basis=basis`.  Nothing dropped from either side.

The docstring table of closed forms was updated in the same commit (it had been
auto-merged to the symmetric-only formulae, which are no longer the whole story).
No exported name was removed: `hull.taper_integral`, `taper_integral_pair`,
`analytic_*` and `hydrostatics_analytic.waterplane_inertia_*` all keep their
previous signatures, the new arguments being keyword-only with defaults.

### 2. `CHANGELOG.md` — union of two `0.21.0` entry sets

Both sides added bullets under `## 0.21.0 — 2026-11-04 / ### Added`.  Kept both
sets in full (waterplane second moments first, then transom stern), and added one
bullet for the split-taper extension above so the changelog describes what the
code does.  Main's `0.20.0` / `0.17.0` sections are untouched.

### 3. `README.md` — one "what it does not provide" paragraph

m2925 added the stability limits (`righting-arm curve`, `free-surface or damage
stability`, `initial GM is reported, and GZ(φ) is not`); main rewrote the hull-form
clause (`two families exist -- sin^0.7·√(z/T) and the parabolic Wigley -- both
fore-aft symmetric, neither with a transom`).  Union: main's hull-form clause with
m2925's stability clauses restored in the same sentence.  Neither side's list of
missing capabilities lost an item.

### 4. `tests/test_asymmetric_hull.py` — one literal, not a behaviour

main asserted the symmetric `basis` string equals the pre-merge text.  The merged
`basis` must advertise the quadrature-derived longitudinal moment (that is what
m2925's half of the record reports), so the literal now includes that clause.
This is the only test edit and it strengthens rather than relaxes the assertion.

## Auto-merged files that were checked, not trusted

`tensorcad/hydrostatics.py` (the flux-identity record), `hydrostatics_analytic.py`
(main's LCB work), `tests/test_public_surface.py`, `pyproject.toml` /
`tensorcad/__init__.py` (both sides bumped to `0.21.0` — no version conflict), and
the whole `docs/` + `examples/` payload of both branches (transom stern, Wigley,
asymmetric hull, `tests/test_transom.py`, `tests/test_wigley.py`,
`tests/test_asymmetric_hull.py`, `tests/test_hydrostatics_inertia.py`).  Each side's
files are present in the merged tree; nothing was dropped by the merge driver.

## Test evidence (merged tree, local, Python 3.11)

| suite | result |
|---|---|
| `test_hydrostatics_inertia.py` + `test_asymmetric_hull.py` + `test_hydrostatics_analytic.py` + `test_analytic_hydrostatics.py` | 116 passed |
| `test_hull.py` + `test_hydrostatics.py` + `test_public_surface.py` + `test_parameters.py` + `test_validation.py` + `test_stability.py` + `test_bonjean.py` | 158 passed |
| `test_transom.py` + `test_wigley.py` + `test_cli.py` + `test_curvature.py` + `test_mesh.py` | 179 passed |
| `test_export` + `test_artifact_export` + `test_gate` + `test_design` + `test_benchmark` + `test_resolution` + `test_mesh_conditioning` + `test_suboff` + `test_suboff_appendages` + `test_voxel` + `test_geometry_field` + `test_adapters_fem` + `test_adapters_ocp` + `test_mirror_contract` | 274 passed, 11 skipped, **2 failed** |

The two failures are
`test_geometry_field.py::test_the_learned_revision_speaks_the_same_protocol` and
`::test_the_learned_field_gradient_is_autograd_and_points_outward`, both
`ModuleNotFoundError: No module named 'safetensors'` — an optional-dependency gap of
this sandbox, present on `main` as well (`test_cli_neural.py`,
`test_neural_pipeline.py`, `test_neural_store.py`, `test_neural_geometry.py`,
`test_neural_format_contract.py`, `test_extract.py` cannot even be collected for the
same reason).  Neither test touches `hull`, `hydrostatics` or the analytic
references.

Total: **727 passed, 11 skipped, 2 environment failures**, no failure attributable to
the merge; the pre-merge branch could not even be imported together with `main`'s
hull module.

## Numerical evidence for the union

`examples/merge_m2925/verify_waterplane_inertia.py` →
`docs/merge_m2925/{results.json,inertia_split.png}` (exits non-zero if a claim
breaks; 14/14 pass).  Case: `L = 10 m`, `B = 2 m`, `T = 0.8 m`, symmetric
`e = 0.7` and asymmetric `(p_fore, p_aft) = (0.5, 0.9)`.

| quantity | analytic | independent quadrature | rel. difference |
|---|---|---|---|
| `I_T`, symmetric | 3.2707254343 m⁴ | 3.2707254343 m⁴ | < 1e-12 |
| `I_L`, symmetric | 76.870781 m⁴ | same | < 1e-10 |
| `I_T`, asymmetric | 3.3329192519 m⁴ | 3.3329192519 m⁴ | < 1e-12 |
| `I_L`, asymmetric (about LCF = 5.258 m) | 77.8807225 m⁴ | same | 1.2e-13 |

Mesh record vs those references (inscribed chords, one sign throughout):

| grid | `I_T` symmetric | `I_L` symmetric | `I_T` asymmetric | `I_L` asymmetric |
|---|---|---|---|---|
| 41x13 | -0.139% | -0.604% | -0.154% | -0.823% |
| 81x25 | -0.034% | -0.198% | -0.038% | -0.297% |
| 161x49 | -0.009% | -0.064% | -0.010% | -0.106% |
| 321x97 | -0.002% | -0.020% | -0.002% | -0.038% |

The exact numbers are in `results.json` (programmatic, nothing transcribed by
hand); the commit message and the mission report quote them from there.

## What is *not* claimed

* No new physics: the merge does not add a hull family, only makes the two
  branches' quantities coexist and be correct on both families.
* The mesh record's inertia is a *discretisation* of the waterplane, so its error
  is the grid's, not the reference's; the reference carries its own quadrature
  evidence (`I_L`) or is a gamma identity (`I_T`).
* Nothing was pushed by this worktree; the merge commit is left for the platform /
  `merge_governor` to publish on `agent/m2925` (fast-forward from `e7e8b18`).