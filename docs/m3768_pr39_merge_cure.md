# PR #39 cure — semantic union of `main` and `agent/m2925` (m3768)

Repo `TensorCAD`. Pull request #39 (`agent/m2925`) was reported `mergeable=false`
against `main`. This document records the conflict inventory, the per-file resolution,
and the numbers that show both sides still work on the merged tree.

## 1. State before the cure

| ref | SHA | note |
|---|---|---|
| `origin/main` | `311f45a` | “Merge #37 (m2873) via merge-governor” |
| `origin/agent/m2925` (PR #39 head, old) | `e7e8b18` | branch rebased onto main `6c568c3` |
| merge base | `6c568c3` | 1 commit unique to the branch, 10 to `main` |

The PR was dirty because both sides touched the same four contracts after `6c568c3`:
the `Hydrostatics` record, the analytic-reference wrapper, the public-surface test and
the version/CHANGELOG header. `main` also added the transom and Wigley hull families,
the asymmetric fore/aft exponents and a rewritten curvature module.

## 2. Conflict inventory and per-file resolution

Three files conflicted (`git merge` with `merge.conflictStyle` default); five more were
auto-merged and inspected by hand for semantic loss.

| file | conflict | resolution (union, nothing dropped) |
|---|---|---|
| `tensorcad/hull.py` | 2 hunks: `analytic_hydrostatics` centre + basis, and the returned record | kept **both** routes: `centre = analytic_longitudinal_centre(length, fore_exponent=fore, aft_exponent=aft)` (main's asymmetric form) **and** the transverse/longitudinal waterplane inertias plus `keel_z` (m2925). The inertias are emitted only when `fore == aft`; otherwise `nan` + a `basis` clause saying why (§4). `basis` carries main's symmetry-aware wording with m2925's quadrature provenance appended. |
| `CHANGELOG.md` | both sides prepended entries under `## 0.21.0 — Added` | concatenated: m2925's `Hydrostatics`-inertia entry first, then main's transom/Wigley entries and the `## 0.20.0` block. No entry reworded or removed. |
| `README.md` | one “what it does **not** provide” paragraph | merged into a single paragraph: main's updated “two families exist … neither with a transom” clause (true after the transom family landed) **plus** m2925's precise stability caveat (“the righting-arm curve, free-surface or damage stability … initial GM is reported, and GZ(φ) is not”). |
| `tensorcad/__init__.py` | auto (both added exports) | verified: exports of `transom`, `wigley`, `stability`, `ANALYTIC_COMPARABLE_QUANTITIES` and the inertia helpers all present, no duplicate `__all__` entry |
| `tests/test_public_surface.py` | auto (both changed the export list) | verified: every symbol the two sides added is asserted once; no assertion lost |
| `pyproject.toml` | auto (both bumped the version) | single `0.21.0` survives; resolution checked against the merged CHANGELOG heading |
| `tests/test_asymmetric_hull.py` | **semantic**, not a textual conflict | the test pinned the old `basis` string exactly; the merged basis is that string plus the inertia provenance, so the pin was updated to the merged sentence (still an exact-equality assertion, not a weakened `in` test) |
| `docs/hydrostatics_verification.md` | auto (branch only) | kept; it now describes a record that the merged code really produces |

## 3. What the union actually computes

The one place where the two sides genuinely overlap is the waterplane second moments.
m2925 derives them from closed forms that assume a fore/aft-**symmetric** waterplane and
take them about midship; main made the fore/aft exponents independent, which moves the
centre of flotation off midship. Reporting m2925's numbers for such a body would be a
silently wrong metacentric radius, so the merged wrapper withholds them (`nan`) exactly
where the assumption fails, and says so in `basis`. Everything else is available on both
kinds of hull: volume, areas, coefficients, centres, `keel_z`.

### 3.1 m2925's route on the merged tree — mesh vs closed form

| grid (nx × nz) | I_T mesh | I_T closed form | rel. error | I_L mesh | I_L closed form | rel. error | keel_z |
|---|---|---|---|---|---|---|---|
| 21 × 9 | 3.252433 m⁴ | 3.270725 m⁴ | -0.5593 % | 75.4882 m⁴ | 76.8708 m⁴ | -1.7986 % | -0.800000 m |
| 41 × 13 | 3.266192 m⁴ | 3.270725 m⁴ | -0.1386 % | 76.4066 m⁴ | 76.8708 m⁴ | -0.6038 % | -0.800000 m |
| 81 × 25 | 3.269597 m⁴ | 3.270725 m⁴ | -0.0345 % | 76.7187 m⁴ | 76.8708 m⁴ | -0.1979 % | -0.800000 m |
| 161 × 41 | 3.270444 m⁴ | 3.270725 m⁴ | -0.0086 % | 76.8217 m⁴ | 76.8708 m⁴ | -0.0638 % | -0.800000 m |

Errors are signed and monotone (inscribed chords under-integrate), ratios 4.0 / 3.0 per
grid doubling in each direction: observed order ≈ 2 for `I_T` and ≈ 1.6 for `I_L`, whose
integers see the waterline closure singularity.

Barge control (closed box, `L·B³/12` and `B·L³/12`):

| quantity | mesh | closed form | rel. error |
|---|---|---|---|
| I_T | 6.666667 | 6.666667 | -1.33e-16 |
| I_L | 166.666667 | 166.666667 | 1.71e-16 |
| keel_z | -0.800000 | -0.800000 | 0 |

### 3.2 main's route on the merged tree — asymmetric taper and the new families

- `analytic_hydrostatics(10, 2, 0.8)` (symmetric): LCB = LCF = 5.000000 m = L/2.
- `analytic_hydrostatics(10, 2, 0.8, fore=0.5, aft=0.9)`: LCB = LCF = 5.257824 m, i.e. +0.2578 m aft of midship, Cb = 0.473274; the two inertias are `nan`
  (`d moments withheld (nan) because their closed forms assume a fore/aft-symmetric waterplane…`), `keel_z` = -0.800 m still reported.
- Wigley at the default parameters: Cb mesh 0.443396 vs analytic 0.444444 (-0.2360 %), displacement 2.7712 m³ vs 2.7778 m³.
- Transom stern (τ = 0.2): Cb mesh 0.513539 vs reference ladder 0.517604 (-0.7852 %).

## 4. Test evidence

Same interpreter, same command (`pytest -q -n 8`, ignoring the two modules that need the
absent optional `safetensors`), three trees:

| tree | result | failing tests |
|---|---|---|
| `main` `311f45a` (baseline) | 2 failed, 802 passed, 11 skipped | `test_geometry_field.py` ×2 (`ModuleNotFoundError: safetensors`) |
| `agent/m2925` `e7e8b18` alone | 4 failed, 647 passed, 11 skipped | the same two + `test_curvature.py` ×2 (branch predates main's curvature rewrite) |
| **merged (this branch)** | **2 failed, 818 passed, 11 skipped** | the same two, both pre-existing environment failures |

So the merge is +16 tests over the `main` baseline with no new failure, and it repairs the
two `test_curvature.py` failures the branch had on its own. Targeted run over the
contract-relevant modules (`test_hydrostatics_inertia`, `test_asymmetric_hull`,
`test_stability`, `test_hull`, `test_public_surface`, `test_transom`, `test_wigley`,
`test_cli`, `test_bonjean`): 294 passed, 0 failed.

Three new tests were added to `tests/test_hydrostatics_inertia.py` to pin the union
itself, so the resolution cannot rot back into a silent loss:

1. symmetric analytic hull reports finite inertias + `tanh-sinh` provenance + `keel_z`;
2. asymmetric analytic hull still reports main's off-midship LCF, but the inertias are
   `nan` and `basis` says `withheld`;
3. mesh record vs symmetric closed form at 81×25 agrees to 5·10⁻³ / 10⁻² relative.

A standalone harness (`tools/verify_merge_m2925.py`) reproduces every number above and
writes `artifacts/m3768/merge_verification.{json,png}`.

## 5. Notes, limits, what was not done

- The local `main` ref in this cluster was stale (`ca0dde3`, 2019-era tip of a different
  arena branch); all work was done against `origin/main` `311f45a` after a fetch.
- A previous cure attempt (`c5274de`, lineage of the pre-rebase `ad86f4f`) exists locally
  but never reached `origin`; this cure is built on the branch head that GitHub actually
  has (`e7e8b18`), which is why the pushed result is a fast-forward, not a force-push.
- `agent/m2925` is checked out in another worktree, so the merge was prepared on
  `cure/m2925-union` (a descendant of `e7e8b18`) and pushed to `agent/m2925` from there.
- Not done: closed forms for the asymmetric waterplane second moments (they would need
  the x-moment re-taken about a displaced LCF). The merged code declares `nan` rather
  than approximating; extending `hydrostatics_analytic` with the split-taper integrals is
  the natural follow-up and is out of scope for a merge cure.
- The two `test_geometry_field.py` failures are a missing optional dependency
  (`safetensors`), identical on `main`, and are not caused by this merge.

## 6. Error judgement against the < 3 % bar

| key quantity | reference | worst rel. error on the merged tree | < 3 % ? |
|---|---|---|---|
| `I_T` of the parametric hull (41 × 13, the default grid) | closed form `B³L/12 · ∫sin³ᵉ` | -0.1386 % | yes |
| `I_L` of the parametric hull (41 × 13) | tanh-sinh quadrature about midship | -0.6038 % | yes |
| `I_T`, `I_L` of a box barge | `L·B³/12`, `B·L³/12` | 1.7e-16 | yes (exact) |
| `keel_z` (all grids) | `−draft` | 0 | yes |
| Wigley `Cb` | `nm/((n+1)(m+1))` | -0.2360 % | yes |
| Transom `Cb` (τ = 0.2) | reference ladder | -0.7852 % | yes |

The basis strings are asserted as exact equality in tests, so the provenance text is part
of the contract, not a comment.
