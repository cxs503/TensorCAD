## PR #39 cured: `main` and `agent/m2925` merged as a semantic union

### ① Case
`TensorCAD` PR #39 (`agent/m2925`) came back `mergeable=false`. Branch head on GitHub was
`e7e8b18` (rebased on main `6c568c3`); `main` had advanced to `311f45a` (transom-stern and
Wigley hull families, independent fore/aft exponents, curvature rewrite). One commit unique
to the branch, ten to main, and both touched the same contracts: the `Hydrostatics` record,
the analytic wrapper `analytic_hydrostatics`, the public-surface test and the version header.
Geometry of the case: the parametric hull family `L=10 m, B=2 m, T=0.8 m`, breadth
`sin(pi x/L)^0.7 * sqrt(z/T)`, water density 1025 kg/m3.

### ② Method
Three-way merge `e7e8b18 + 311f45a` on a scoped worktree, per-file semantic resolution
(no overlay), then the same pytest command run on three trees (main baseline, branch alone,
merged) so "not below baseline" is measured rather than asserted. Resolution:

| file | conflict | resolution |
|---|---|---|
| `tensorcad/hull.py` | 2 hunks | kept main's asymmetric centre **and** m2925's waterplane inertias + `keel_z`; inertias withheld (`nan` + `basis` clause) where the symmetric closed forms do not hold |
| `CHANGELOG.md` | both prepended under `0.21.0 Added` | concatenated, no entry dropped or reworded |
| `README.md` | one limitations paragraph | single paragraph: main's two-family clause + m2925's stability caveat (GZ, damage/free-surface) |
| `tests/test_asymmetric_hull.py` | semantic pin on the old `basis` string | pin updated to the merged sentence (still exact equality) |
| `__init__.py`, `test_public_surface.py`, `pyproject.toml` | auto | inspected: both sides' exports/version survive once |

The only genuine overlap is the waterplane second moment: m2925 derives it about midship for
a symmetric waterplane, main made fore/aft exponents independent so LCF leaves midship.
Union = report the moments where the assumption holds, `nan` where it does not — never the
symmetric number for an asymmetric body.

### ③ Results (values injected from `artifacts/m3768/merge_verification.json`)

| grid | I_T mesh [m4] | I_T closed form | rel. err | I_L mesh [m4] | I_L closed form | rel. err |
|---|---|---|---|---|---|---|
| 21x9 | 3.252433 | 3.270725 | -0.5593% | 75.4882 | 76.8708 | -1.7986% |
| 41x13 | 3.266192 | 3.270725 | -0.1386% | 76.4066 | 76.8708 | -0.6038% |
| 81x25 | 3.269597 | 3.270725 | -0.0345% | 76.7187 | 76.8708 | -0.1979% |
| 161x41 | 3.270444 | 3.270725 | -0.0086% | 76.8217 | 76.8708 | -0.0638% |

Barge control: I_T 6.666667 vs 6.666667 (-1.3e-16), I_L 166.666667 vs
166.666667 (1.7e-16), keel_z -0.800000.

main's route on the merged tree: symmetric LCB = LCF = 5.000000 m (= L/2);
asymmetric (fore 0.5 / aft 0.9) LCB = LCF = 5.257824 m, +0.2578 m aft of
midship, Cb = 0.473274, inertias `nan` with `basis` stating why, `keel_z` = -0.800 m.
Wigley Cb mesh 0.443396 vs analytic 0.444444 (-0.2360%);
transom (tau=0.2) Cb mesh 0.513539 vs reference 0.517604 (-0.7852%).

Tests, same command on three trees: main `311f45a` = 2 failed / 802 passed / 11 skipped;
branch alone `e7e8b18` = 4 failed / 647 passed; **merged = 2 failed / 818 passed / 11 skipped**.
The 2 failures are `ModuleNotFoundError: safetensors` in `test_geometry_field.py`, identical on
main (missing optional dependency). Merging also repairs the two `test_curvature.py` failures
the branch had alone. Targeted contract modules: 294 passed. Three new tests pin the union.

Push: `agent/m2925` `e7e8b18` -> `9663314` (fast-forward, main tip is now an ancestor, so the
PR is mergeable); `merge_governor` can confirm on its next cycle.

### ④ Figure
Left: signed mesh-vs-closed-form error of both waterplane second moments against grid
resolution (four grids, monotone O(h^2)-ish). Right: main's asymmetric route after the
merge — LCB/LCF sit aft of midship at fore=0.5/aft=0.9 while the symmetric case is exactly L/2.
A field contour is not meaningful here: the task changes no flow field, it reconciles two
analytic-reference routes; the comparison curves against the closed forms are the evidence.

### ⑤ Notes / limits
- Local `main` ref in this cluster was stale (`ca0dde3`); all work against `origin/main` after fetch.
- An earlier cure attempt (`c5274de`) existed only locally on the pre-rebase lineage; this cure is
  built on the head GitHub actually has, so the push is a fast-forward, not a force-push.
- `agent/m2925` is checked out in another worktree, so the merge was prepared on
  `cure/m2925-union` and pushed from there; the local ref of `agent/m2925` is therefore still
  `e7e8b18` while `origin/agent/m2925` is `9663314`.
- Not done: closed forms for the asymmetric waterplane moments (needs the x-moment about a
  displaced LCF). The code declares `nan` instead of approximating; that is a follow-up, not a cure.

### ⑥ Error judgement vs the 3 % bar
| quantity | reference | merged-tree error | <3%? |
|---|---|---|---|
| I_T, hull, default 41x13 | closed form B^3 L/12 * beta(3e) | -0.1386% | yes |
| I_L, hull, default 41x13 | tanh-sinh quadrature | -0.6038% | yes |
| I_T, I_L barge | L B^3/12, B L^3/12 | 1.7e-16 | yes (exact) |
| keel_z all grids | -draft = -0.8 | 0 | yes |
| Wigley Cb | nm/((n+1)(m+1)) | -0.2360% | yes |
| transom Cb | reference ladder | -0.7852% | yes |

No quantity is reported without a reference; where a reference is impossible (asymmetric
waterplane moments) the value is withheld as `nan` and the reason is in `basis`.


*(figure: `artifacts/m3768/merge_verification.png`; payload: `artifacts/m3768/report_payload.json`)*
