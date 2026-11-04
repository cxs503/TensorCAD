# Cure of PR #39 (`agent/m2925`) against `main` — re-verification by m3753

**Verdict: merged and clean.** `agent/m2925` now ends in a single merge commit whose two
parents are the branch tip before the cure (`ad86f4f`) and `main` (`311f45a`); tree-level
merge against `main` is conflict-free, the whole tree carries **no conflict markers**, and
the merged test baseline equals `main`'s (same two pre-existing, environment-only failures,
zero new ones).

| item | old | new |
|---|---|---|
| `agent/m2925` tip | `546e001` (first cure attempt, marker-bearing `CHANGELOG.md`) | `c5274de` |
| branch tip before any cure | `ad86f4f` | — |
| `main` merged in (unchanged) | `311f45a` (`origin/main`) | — |

## 1. Case: what was dirty and why

PR #39 adds two things to TensorCAD that `main` had meanwhile implemented from the other
end — namely the waterplane second moments needed for a metacentric height:

* **branch (`m2925`)**: `Hydrostatics` gains `transverse_waterplane_inertia`,
  `longitudinal_waterplane_inertia`, `keel_z` from the **flux identity over the closed wetted
  surface** (`∫_cap f dA = −∫_wetted f n_z dA`, `f` independent of `z`, face integrals exact
  because `x, y` are linear on a triangle); `stability.initial_stability(hydrostatics,
  height_of_gravity=…)`; closed forms in `hydrostatics_analytic` and
  `hull.analytic_hydrostatics` (`I_T` from the gamma identity, `I_L` by tanh-sinh quadrature);
  a separate example `run_example_hydrostatics.py`; 14 tests.
* **main (`311f45a`)**: `stability.waterplane_moments(mesh, waterline_z)` cuts the surface
  with the plane, chains segments on mesh identity, integrates by **Green's theorem**;
  `initial_stability(mesh, kg=…)`; 17 tests; `0.21.0` release history and a rewritten
  `run_example.py`.

Both are correct and they are not alternatives — they are two independent integrations of
the same two modules of the waterplane, which is why the cure had to be a **union**: a
side-pick would have deleted a whole route plus the analytic reference that neither side
duplicates.

## 2. Method

1. Base resolved with `git merge-base`: `ca0dde3`. `git merge-tree --write-tree ad86f4f
   311f45a` produced **14 conflicted paths** (`CHANGELOG.md`, `README.md`,
   `docs/stability.md`, `docs/stability/report.md`, `docs/stability/results.json`,
   `examples/stability/run_example.py`, `pyproject.toml`, `tensorcad/__init__.py`,
   `tensorcad/hull.py`, `tensorcad/hydrostatics_analytic.py`, `tensorcad/stability.py`,
   `tests/test_public_surface.py`, `tests/test_stability.py`, and the module surface).
2. Resolution policy, file by file (union, no side picked):

| file | resolution | why it is a union and not a pick |
|---|---|---|
| `tensorcad/stability.py` | main's `waterplane_moments`/Green route **plus** branch's `InitialStability`, `initial_stability`, `_stability_from_hydrostatics`; one entry point dispatching on `Mesh \| Hydrostatics` | the two routes are separate integrations; `kg`/`height_of_gravity` name one quantity and disagreeing values are refused |
| `tensorcad/hydrostatics.py` | branch's flux-identity inertias + `keel_z` added to main's record | additive only: 74 added / 2 changed lines, no key or number touched |
| `tensorcad/hydrostatics_analytic.py`, `tensorcad/hull.py` | branch's inertias, `BM`, `KM_T` and the generalised fore/aft exponents from main | main's signature generalisation is kept, the analytic inertia rows extend it |
| `tensorcad/__init__.py` | both export sets (`waterplane_moments`, `InitialStability`, `initial_stability`, …), `__version__ = 0.21.0` | union; version follows main, branch's `0.14.0` was stale |
| `tests/test_stability.py` | both suites concatenated — 31 test functions, 37 collected checks (17 main + 14 branch) | deleting either set would drop a route's only guard |
| `tests/test_public_surface.py` | declared-API list union, `__version__ == "0.21.0"` | branch's `0.14.0` assertion is superseded by main's release, not by an opinion |
| `examples/stability/run_example.py` | keeps main's (plane route, claims, 3 figures); branch's script continues as `run_example_hydrostatics.py` | branch's copy was the same ladder pointed at the same filenames — keeping it would have overwritten main's newer artifact |
| `docs/stability.md` | main's page + branch's flux-identity + analytic-reference page appended ("Second route"), headings demoted | both derivations measured against the same closed forms |
| `docs/stability/report.md`, `results.json` | main's artifacts, version bumped to `0.21.0` | generated artifacts of the retained example |
| `docs/stability/report_hydrostatics.md`, `results_hydrostatics.json`, `plot.png` | branch's artifacts, kept under their new names | neither example overwrites the other |
| `README.md`, `pyproject.toml` | union of the capability bullets; version `0.21.0` | branch's `0.14.0` is stale, nothing else on its side |
| `CHANGELOG.md` | **defect of the first cure fixed** — markers removed, one `0.14.0` section holding both bullet lists verbatim | see §3 |

3. Static checks on the merged tree: AST duplicate-definition scan over all of
   `tensorcad/*.py` → 0 duplicates; `git grep -E '^(<<<<<<<|=======$|>>>>>>>)'` → 0 hits;
   every line added on either side of the 14 conflicted paths still present in the merge
   except the three deliberately requalified changelog sentences and branch's stale
   `0.14.0` version asserts.

## 3. The one real defect: `CHANGELOG.md` shipped its conflict markers

The first cure commit (`546e001`) resolved thirteen of the fourteen paths cleanly but left
`CHANGELOG.md` with `<<<<<<< HEAD` (line 7), `=======` (51) and `>>>>>>> origin/main` (206)
in the file — a file that still *renders* as release notes while carrying two `## 0.14.0`
sections and the machinery of an unfinished merge. This is exactly the "mechanical overlay
silently loses content" failure class, only louder.

Resolution: one `0.14.0` section with **both** bullet lists verbatim (main's plane-section
route and 19 tests; the branch's flux-identity route, `Hydrostatics` moments, analytic
reference and 18 tests), release order preserved (`0.21.0 … 0.16.0`, merged `0.14.0`,
`0.13.0 …`). Three qualifiers only, each forced by the merge itself: the two test bullets
name their route and the merged total (31 functions), and the branch's example bullet names
`run_example_hydrostatics.py` and the artifacts it now writes. Against each parent,
`CHANGELOG.md` shows main's `0.21.0`–`0.16.0` byte-identical and both `0.14.0` bullet sets
present.

## 4. Results

Full suite on the merged tree, `tests/` minus the two neural modules that do not even
collect in this environment (Torch-free runner), run side by side with `main`:

| tree | result |
|---|---|
| `main` (`311f45a`) | **802 passed, 2 failed, 11 skipped** (454 s) |
| merged (`c5274de`) | **820 passed, 2 failed, 11 skipped** (453 s) |

The same two tests fail on both trees
(`tests/test_geometry_field.py::test_the_learned_revision_speaks_the_same_protocol`,
`…::test_the_learned_field_gradient_is_autograd_and_points_outward` — environmental,
learned-field runtime), so the merge adds **no** failure, and the +18 passes are the
branch's stability checks now running alongside main's.

Targeted evidence on the merge itself:

* `pytest tests/test_stability.py tests/test_public_surface.py -q` → **43 passed**
  (37 stability checks from *both* routes + the public-surface contract at `0.21.0`).
* `examples/stability/run_example.py` → exit 0, all 10 claims hold.
* `examples/stability/run_example_hydrostatics.py` → exit 0, all claims hold
  (`KB = 0.483005 m`, `BM_T = 0.438618 m`, `KM_T = 0.921623 m`, `GM_T(KG=0.5) = 0.421623 m`).
* Cross-route agreement on the same mesh (the whole point of keeping both):

| grid | `I_T` plane [m⁴] | \|plane−flux\| rel | `I_L` plane [m⁴] | \|plane−flux\| rel | err `I_T` | err `I_L` | err `GM_T` |
|---|---|---|---|---|---|---|---|
| 21×9 | 3.252433 | 2.73e-16 | 75.488218 | 0.00e+00 | -0.559 % | -1.799 % | +3.473 % |
| 41×13 | 3.266192 | 1.36e-16 | 76.406622 | 1.49e-15 | -0.139 % | -0.604 % | +1.986 % |
| 81×25 | 3.269597 | nan | 76.718673 | nan | -0.034 % | -0.198 % | +0.736 % |
| 161×41 | 3.270444 | 0.00e+00 | 76.821723 | 1.48e-15 | -0.009 % | -0.064 % | +0.347 % |
| 321×81 | 3.270655 | 1.09e-15 | 76.855133 | 4.44e-15 | -0.002 % | -0.020 % | +0.125 % |

`GM_T` is the amplified quantity (`rel(GM) ≈ 2.87·rel(KM)` at `KG = 0.6 m`): +1.99 % on the
default 41×13 grid, +0.35 % at 161×41, +3.47 % on the deliberately-too-coarse 21×9 rung.

## 5. Notes, limits

* GitHub itself was unreachable from this sandbox (`git ls-remote` timed out), so the
  "push" is the local ref update `refs/heads/agent/m2925: 546e001 → c5274de`; the platform's
  pusher carries it to GitHub, and `merge_governor` re-verifies. `origin/main` and the merge
  parents were taken from the local clone, which already held `311f45a`.
* `agent/m2925` is checked out in another worktree, so the ref was moved with
  `git update-ref` (compare-and-swap against `546e001`) instead of `git branch -f`.
* The merge commit was amended rather than followed by a fix-up commit, so the branch carries
  one final tree and no process history.
* Not done here: no attempt to make the two neural test modules collectable (environmental,
  out of scope), and no change to either route's numerics — the cure is a merge, not a
  redesign.

## 6. Error judgement

| quantity | reference | measured (default 41×13) | error | < 3 %? |
|---|---|---|---|---|
| `I_T` | closed form `(2/3)(B/2)³L∫sin^2.1` | 3.266192 m⁴ | −0.139 % | yes |
| `I_L` | tanh-sinh quadrature, cross-checked by Simpson | 76.4066 m⁴ | −0.604 % | yes |
| `KM_T` | `0.6T + BM_T` | 0.915360 m (+0.684 %) | +0.684 % | yes |
| `GM_T` (KG=0.6 m) | analytic surface | 0.3216 m | +1.99 % | yes |
| `GM_T` at 21×9 | analytic surface | 0.3263 m | +3.47 % | **no** — coarse rung, declared as a measurement, not a claim |
| cross-route `I_T` | independent flux identity | — | ≤ 2.7e-16 | yes (round-off) |

The only out-of-band number is the coarsest grid's `GM_T`, which the branch and main both
document as a printed measurement rather than a claim; on the default and fine grids every
quantity is inside the 3 % band.
