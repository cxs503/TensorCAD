# Independent re-check of `tensorcad.section_gz` (adversarial review)

**Mission**: m16430 (selfstart verifier seat) · **Date**: 2026-11-12 · **Branch**: agent/m16430
**Subject under review**: `tensorcad/section_gz.py` + `tests/test_section_gz.py` +
`examples/section_gz_exact/` + `docs/section_gz_exact.md`, merged in m16358 (commit `680e5e6`,
merged as `ef5f838`) — roughly four hours before this mission started.
**Disposition**: **approve with notes** (one low-severity defect fixed, one mislabelled claim
corrected in this document, no numerical result of the merged feature changed).

## 1. Why this subject, and what "independent" means here

The verifier seat's job is *verification debt*: take a feature that is already merged and ask
whether its claims survive a route that was not written by the same author.  `section_gz` is the
newest merge in the repository and it ships a **new oracle** — which is exactly the kind of code
that is dangerous precisely because everything else is checked against it.  If the oracle were
wrong, its own tests would agree with it.

The independent route (`tests/section_gz_oracle.py`) shares no algorithm with the module:

| | library (`section_gz`) | independent route (this mission) |
|---|---|---|
| ellipse cut | closed form of the unit disk, `pi/2 + asin d + d sqrt(1-d^2)`, affine map | rotated conic recovered **numerically** from the quadratic form; horizontal chords read off the conic coefficients |
| ellipse integrals | analytic | composite Simpson after `zeta = z_c + H sin(theta)` (1024 panels, error ~1e-11 relative) |
| polygon cut | Sutherland-Hodgman clip | horizontal-chord scan over the polygon edges |
| polygon integrals | Green's theorem on the clipped vertex list | exact per-interval trapezoid + 3-point Simpson (exact for the quadratic integrand) |
| level solve | bisection, 200 halvings | bisection, 120 halvings, separate implementation |
| box reference | the library's own `box_righting_arms_analytic` | the published wall-sided formula `sin(phi)(GM + BM tan^2(phi)/2)` transcribed from Biran / Rawson & Tupper |

The two implementations do share the *conventions* (`(y,z) -> (y c + z s, -y s + z c)`,
`GZ = KN - KG sin phi`, `KG` above the keel, `KN` the submerged centroid's `y`).  Those are
documented in both files and are independently confirmed by the box row below, where the
published wall-sided formula reproduces the arms to `2.2e-16 m` with no shared code at all.

## 2. Checks performed

Run: `python examples/section_gz_independent/run_example.py` →
`examples/section_gz_independent/outputs/{metrics.json,report.md,plot.png}`.
Tests: `python -m pytest tests/test_section_gz_independent.py tests/test_section_gz.py
tests/test_examples_bootstrap.py -q` → **104 passed** (13 of them new).

| # | check | result | verdict |
|---|---|---|---|
| 1 | 2:1 ellipse, closed form vs independent quadrature, 19 heels 0-90 deg | `2.69e-13 m` | PASS (11 orders inside the 3 % rule) |
| 2 | same over both heel directions (±90 deg) | `2.69e-13 m` | PASS |
| 3 | antisymmetry `GZ(-phi) = -GZ(phi)` on the independent curve | `7.8e-16 m` | PASS |
| 4 | box / V section / hexagon rotated 17 deg, polygon oracle vs exact slice integration | `4.4e-16 / 2.2e-16 / 1.7e-16 m` | PASS |
| 5 | box inside the wall-sided limit, vs the published formula (independent transcription) | `2.2e-16 m` | PASS |
| 6 | box above the limit, wall-sided formula genuinely departs | `2.26e-2 m` at 30 deg (min over 30-90 deg) | PASS (not vacuous) |
| 7 | GM: committed artifact `1.1558433455319073` vs independent waterplane moments `1.1558433455315373` | `3.7e-13 m` | PASS |
| 8 | GM: library `dGZ/dphi(0)` (`1.155843345531`) vs independent waterplane moments | `3.7e-13 m` | PASS |
| 9 | GM: independent finite difference at 0.5 deg (`1.1556831`) vs the same | `1.6e-4 m` | PASS — this residual is the *finite-difference truncation*, not a disagreement |
| 10 | KB: library centroid `0.6813855514894` vs independent `0.6813855514891` | `2.7e-13 m` | PASS |
| 11 | GZ_max: committed `0.4436855620569` (5 deg grid, 45 deg) vs independent on the same grid | `2.5e-13 m` | PASS |
| 12 | vanishing angle at KG = 1.6 m: committed and library `31.7924237378636` vs independent `31.7924237379007` | `3.7e-11 deg` | PASS |
| 13 | mesh route (16…256 inscribed chords) vs the **independent** oracle, not the closed form | orders 2.02, 2.05, 2.05, 1.99; error `2.13e-5 m` at 256 chords | PASS (second order as an inscribed chord predicts) |
| 14 | input contract: 9 documented rejections (concave, zero area, non-finite, heel 180, waterline at keel, over-large displacement, non-positive length, ellipse draft at an end, negative semi-axis) | 9/9 fire with the documented message | PASS |
| 15 | clockwise winding accepted and numerically identical | `0.0 m` | PASS |
| 16 | `ellipse_section` really is inscribed, for `sides` = 5, 16, 17, 18, 20, 32 | quadric max `1.0000000000000002` after the fix | PASS after fix (see F1) |

## 3. Findings

### F1 (fixed in this mission, severity: medium-low) — `ellipse_section` was not an inscribed polygon for side counts without a keel vertex

`ellipse_section` samples `sides` points on the ellipse and then subtracts the lowest sample's
height so the keel sits at `z = 0`.  The sampling started at angle `0`, so a vertex sat on the
keel only when `sides` was a multiple of four.  For any other count the shift *translated* the
polygon through the exact section: it kept the right area ordering (still smaller than the
ellipse) but its boundary left the ellipse.

| sides | max quadric on the boundary, **before** | **after** fix |
|---|---|---|
| 5 | — | 1.000000000 |
| 16 | 1.000000 | 1.000000 (vertex set unchanged) |
| 17 | **1.008513** | 1.000000 |
| 18 | **1.030154** | 1.000000 |
| 20 | 1.000000 | 1.000000 |

A quadric value above 1 means the polygon is outside the ellipse.  The docstring promised an
"inscribed polygon, so the chord error has a known sign"; for `sides = 18` the boundary left the
section by 3.0 %, and the keel shift also moved the draft datum, so any convergence study done
with an odd or `2 mod 4` side count would have measured a *systematic offset* on top of the
chord error.  No committed number was affected: the shipped convergence table uses
16/32/64/128/256, all multiples of four.

Fix: sample at `2 pi i / sides - pi/2 + turn`, which guarantees a keel vertex.  For side counts
divisible by four the vertex set is *identical* (the phase is a whole number of sample steps), so
the merged evidence is untouched — re-running `examples/section_gz_exact/run_example.py` after
the fix reproduces every published number, with differences only in the last ulp of the
convergence rows (`5.854042762372225e-3` → `5.854042762372336e-3`, relative `2e-14`, from `sin`
and `cos` being evaluated at different arguments).  Regression test:
`tests/test_section_gz_independent.py::test_ellipse_section_is_inscribed_in_the_exact_ellipse`.

### F2 (documentation, severity: low) — a `max` is published under the name `min`

`examples/section_gz_exact/run_example.py` computes the departure above the wall-sided limit as
`outside_error = max(...)` and publishes it as `"min_error_outside_limit_m"`, and
`docs/section_gz_exact.md` §3 quotes it as "box, **min** error vs wall-sided formula above the
limit | 3.56e-01 m".  The value `0.3565 m` is the departure at `limit + 15 deg = 41.57 deg`; the
smallest departure just above the limit is `0.0456 m` at `31.57 deg` and `0.0226 m` at `30 deg`
(independent route, table 6 above).  The non-vacuity conclusion is unaffected, but the label
overstates the *minimum* by a factor of ~8 and would mislead anyone quoting it as a bound.  This
document quotes the correct quantity; the merged artifact is left as-is (changing a key name in a
committed `results.json` is a separate, mechanical change, and no consumer reads it).

### F3 (advisory, severity: low) — `GZ_max` is a sampled maximum

`docs/section_gz_exact.md` quotes "GZ_max 0.443686 m at 45 deg".  That is `max_arm()` over the
5 deg grid; on a 1 deg grid the independent curve peaks at **0.444443 m at 43 deg** (0.17 %
higher).  The module's `max_arm()` docstring does say "over the sampled angles", but the doc's
table does not; a reader could take `0.443686` as the true maximum.  Not a defect, an
under-statement of scope.

### F4 (nit, severity: informational) — `area_under_curve(limit_deg)` truncates, it does not interpolate

Passing a limit between two sampled heels drops the straddling interval entirely:
`area_under_curve(45.0) == area_under_curve(47.0)` on a 5 deg grid.  Fine for a dynamic-stability
area computed to a sampled heel, but the docstring should say "up to the last sampled heel not
beyond `limit_deg`".  `SectionArmCurve.angle_of` also requires an exact float match
(`angle_of(45)` raises `KeyError` for the int); a tolerant lookup would be friendlier.

## 4. Checks not run (honest scope)

* **No comparison with experimental or model-test data.**  Prismatic sections have no
  experimental counterpart; the review is oracle-against-oracle plus published closed forms, which
  is what this module claims to be.
* **No non-convex, multi-valued, or free-surface cases.**  The module rejects them by design
  (`_validate_section`), and that rejection is verified (table 14), not its numerical behaviour.
* **No performance or memory characterisation**; not part of the claim.
* **The `tensorcad.gz` mesh integrator itself was not re-derived** — it is exercised only as a
  *consumer* (check 13).  Its own verification is out of scope here.
* **The full test suite was not run** in this mission (load average 109 on a 128-core host, and
  the verifier seat's capacity rule is light work): only the three affected suites (104 tests)
  plus the example scripts were executed.

## 5. Claim-by-claim audit of `docs/section_gz_exact.md`

| claim as written | status after this review |
|---|---|
| V section exact vs mesh, max abs error `4.11e-15 m` over 19 heels | reproduced (`examples/section_gz_exact` re-run, identical to the printed digits) |
| V section displacement residual `5.06e-14` | reproduced |
| waterline crosses a section vertex at `26.565 deg` | reproduced (`26.56505117707799`) |
| hexagon rotated 17 deg, max error `3.33e-16 m` | reproduced; independent route agrees to `1.7e-16 m` |
| box vs wall-sided below the limit `1.11e-16 m` | reproduced; independent transcription of the formula agrees to `2.2e-16 m` |
| box departure above the limit `3.56e-01 m` **as a minimum** | the number is a maximum; the minimum is `2.26e-2 m` at 30 deg (F2) |
| ellipse `dGZ/dphi(0)` equals the waterplane GM to `2.2e-16 m/rad` | reproduced; the waterplane GM computed here independently agrees to `3.7e-13 m` |
| ellipse chord convergence orders `2.02, 2.05, 2.05, 1.99` | reproduced against the **independent** oracle: `2.02, 2.05, 2.05, 1.99` |
| `GZ_max 0.443686 m at 45 deg` | true on the 5 deg grid; `0.444443 m at 43 deg` on a 1 deg grid (F3) |
| vanishing angle `31.79 deg` at KG = 1.6 m | reproduced: `31.7924237378636` vs independent `31.7924237379007` |
| wall-sided formula errs by `65.70 %` at 30 deg on the 2:1 ellipse | reproduced (independent wall-sided transcription) |
| "18 tests" in `tests/test_section_gz.py` | confirmed: 18 collected, 18 passed |
| `ellipse_section` is an inscribed polygon | **false for side counts without a keel vertex** (F1, fixed) |

## 6. Limitations of this review

I am the same class of agent that wrote the merged feature, and the independent route was written
in the same session as the review — inspector and author are not institutionally separated.  The
mitigation is that the two routes share no algorithm and that every claim above is a number a
third party can recompute with one command; the residual risk is a *shared* misreading of the
heel convention, which is why the box row is anchored to a formula from the literature rather
than to either implementation.

## 7. Files

| path | what |
|---|---|
| `tests/section_gz_oracle.py` | the independent route: conic slice quadrature, exact polygon slices, wall-sided formula |
| `examples/section_gz_independent/run_example.py` | the run behind every number here |
| `examples/section_gz_independent/outputs/` | `metrics.json`, `report.md`, `plot.png` |
| `tests/test_section_gz_independent.py` | 13 tests pinning the cross-route agreement and the F1 regression |
| `tensorcad/section_gz.py` | the F1 fix (sampling phase in `ellipse_section`) |
