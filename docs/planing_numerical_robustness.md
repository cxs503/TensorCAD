# The planing wetted-topside integral: two numerical defects and the fix

An adversarial re-verification of `tensorcad.planing` (merged as the planing-hull
family, `docs/planing_hull.md`) against an *independent* high-precision reference.
The conclusion is uncomfortable but small: the closed forms for volume, waterplane
area, LCB, LCF and VCB survive the audit to `1e-7` relative; the **wetted-surface
area** did not. It was wrong by up to 0.85 % over most of the family and returned
`inf` on 4 of 420 valid parameter sets, which means the number the module validated
the mesh *against* was broken on those hulls. Both are fixed in `tensorcad/planing.py`
and pinned by `tests/test_planing_topside_robustness.py` (93 cases).

Everything below is reproduced by

```bash
python tools/planing_topside_robustness.py     # 420-case sweep, 3 routes, exits 1 on any failed claim
python -m pytest tests/test_planing_topside_robustness.py -q
```

which writes `artifacts/planing_numerical_robustness/results.json`,
`planing_topside_robustness.png` and `extreme_flare_convergence.json`.

## 1. What was being verified

The wetted surface of the family is three pieces (`analytic_wetted_area`):

| piece | form |
|---|---|
| V-bottom | `Awp * sec(beta)` (closed form: the ruling lines are straight keel-to-chine) |
| wetted topside | `2 int_0^x0 (T - b tan b) sqrt(1 + b'(x)**2) dx` (**quadrature**) |
| transom | `T**2 / tan(beta)` (closed form) |

with `b(x) = (B/2)(x/L)**q` and `x0` the station where the chine meets the waterline.
Only the middle row was suspect, because it is the only one not closed. The shipped
rule was a fixed **4000-interval composite Simpson** in the substituted variable
`x = x0 t**(2/q)`, with the arc element written as `sqrt(1 + slope**2) * dx/dt`.

## 2. The reference the fix is judged against

A double-precision audit of a quadrature rule is worthless if the auditor's own
integral is evaluated in a variable where the integrand is singular, so the reference
uses **two physical variables**, chosen by regime, neither of them the module's
`t**e` map:

* `q < 1` (knife bow, `b' ~ x**(q-1) -> inf`): integrate over the chine half-breadth
  `u = b(x)`, where `sqrt(1 + b'**2) dx = sqrt(1 + (dx/du)**2) du` and
  `dx/du ~ u**(1/q - 1) -> 0` at the stem — smooth.
* `q >= 1` (`b'` bounded): integrate in `x`.

Both are `mpmath.quad` at 40-50 digits (tanh-sinh). Cross-check:
`scipy.integrate.quad` (QUADPACK, adaptive, in `x`) agrees with the mpmath reference
to **3.1e-10** relative over the whole sweep, worst at
`(deadrise 45 deg, q = 0.6, draft 2.0 m)`. The naive "just quad in `x`" reference is
itself 0.25 % wrong on a `q = 0.05` hull (14.6513 against 14.6877 m²) — the same trap
the defect fell into, which is why it is documented here rather than silently used.

## 3. Defect A — overflow to `inf` (4 of 420 cases)

`slope = b'(x) = c q x**(q-1)`; at `q = 0.05` the chine flare is so steep that near
the stem `slope ~ 1e21`, and `sqrt(1 + slope*slope)` squares it to `inf`, so the
topsides contribution — and therefore the whole wetted area — came back `inf`:

| deadrise [deg] | q | draft [m] | legacy wetted area [m²] | reference [m²] | shipped after fix [m²] |
|---|---|---|---|---|---|
| 30 | 0.05 | 0.05 | `inf` | 2.008660254 | 2.008660254 |
| 45 | 0.05 | 0.05 | `inf` | 1.419213562 | 1.419213562 |
| 59 | 0.05 | 0.05 | `inf` | 1.169638476 | 1.169638476 |
| 59 | 0.05 | 0.20 | `inf` | 4.714602438 | 4.714602438 |

Consequences that mattered in practice, not just in theory:

* `PlaningHullModel.metadata()["analytic_reference"]["wetted_surface_area"]` was
  `float('inf')` — strict JSON (`allow_nan=False`) refuses, Python's default encoder
  happily writes `Infinity`, an invalid JSON token for most consumers;
* the mesh-vs-closed-form check compared a real area to `inf`, so its relative error
  was `nan`: the check *failed* but only as a `warning`, and `generate_planing_hull`
  still returned the model. A user reading `validation.ok` saw a failure with no
  number attached to it.

## 4. Defect B — uniform-grid truncation (0.08-0.85 % low)

With `q >= 2` the module's `x = x0 t**(2/q)` substitution puts the singularity of
`dx/dt = x0 (2/q) t**(2/q - 1)` at `t = 0`, where the fixed uniform grid cannot resolve
it; with `q = 0.05` the arc element instead becomes a `t**39` boundary layer against
`t = 1`. Either way the 4000-interval rule under-integrates:

| deadrise [deg] | q | draft [m] | legacy error | shipped error |
|---|---|---|---|---|
| 45 | 4.0 | 2.0 | **-0.8514 %** | -2e-17 |
| 30 | 4.0 | 2.0 | -0.8143 % | +1e-16 |
| 59 | 4.0 | 2.0 | -0.8506 % | -1e-16 |
| 59 | 3.0 | 0.55 | -0.1129 % | -1.3e-16 |
| 12 | 1.5 | 1.0 | -0.00018 % | +1.5e-14 |

126 of the 420 sweep cases have `q >= 2` and all of them carried this bias. It is
inside the module's own 2 % wetted-area band, which is exactly why it survived the
41 original tests: a biased *reference* and a mesh built from the same section
functions agree with each other. The 0.85 % bias ate 42 % of the 2 % budget, so a real
mesh defect of that size would have been reported as 0.15 %.

## 5. The fix (three changes, one file)

1. **Arc element by `hypot`, not by squaring.** `sqrt(1 + b'**2) dx` is
   `hypot(dx/dt, db/dt)`; `hypot` never forms `slope**2`, so it cannot overflow.
2. **Substitution per regime.** `e = 2/q` for `q < 2` (kills the `b' -> inf`
   singularity), `e = 1` for `q >= 2` (`b'` already vanishes at the stem, and `2/q`
   would manufacture a `t**(2/q - 1)` singularity). Both branches write the chine
   half-breadth exactly as `b(x0) t**(e q)`, because `b` *is* a power law — no
   approximation is introduced.
3. **Adaptive 10-point Gauss-Legendre** instead of the fixed Simpson grid, bisection
   driven by the difference between a panel and its two halves, relative tolerance
   `1e-12` (`TOPSIDE_QUADRATURE_RTOL`), depth cap 24. `analytic_wetted_area` gained an
   `rtol` argument; the old `samples=` keyword is still accepted and ignored.

Plus one defensive change outside the integral: a closed form that is non-finite or
zero is now an **`error`-severity** check (it is a broken reference, not sampling
noise) and `generate_planing_hull` raises on post-hydrostatics errors, so a hull can no
longer be handed back with `inf` in its metadata.

## 6. Verification

`tools/planing_topside_robustness.py`, 420 valid parameter sets
(`deadrise 1-59 deg x q 0.05-4.0 x draft 0.05-2.0 m`, `L = 10 m`, `B = 3 m`):

| claim | statement | result |
|---|---|---|
| C1 | the legacy rule returns a non-finite area on some valid hull | **holds**: 4/420, all `q = 0.05`, `deadrise >= 30 deg`, `draft <= 0.2 m` |
| C2 | the shipped rule is finite and within `1e-9` of the reference everywhere | **holds**: worst `2.02e-14` at `(45 deg, q = 1.5, draft 2.0 m)` |
| C3 | where the legacy rule had converged, the fix is a no-op | **holds**: 248 cases, max shift `2.39e-11` relative |
| C5 | where the legacy rule was wrong, the fix lands on the reference | **holds**: 168 cases, worst shipped error `2.02e-14` |
| C4 | the two reference routes agree | **holds**: SciPy vs mpmath `3.06e-10` worst |

The published closed-form values in `docs/planing_hull.md` are **unchanged** to six
decimals (`chine_immersed` 25.670434 m², `chine_dry` 21.636391 m²): those hulls have
`q = 0.6`, where the old rule had already converged, which C3 quantifies.

### 6.1 The severe case, honestly

Fixing the reference makes the *mesh* side visible at `q = 0.05`, where the bow flare
is nearly vertical and uniform stations are a poor sampler
(`extreme_flare_convergence.json`, `L = 10 m`, `B = 3 m`, `T = 0.05 m`,
`deadrise 30 deg`, `q = 0.05`):

| nx x nz | triangles | volume error | waterplane error | wetted error |
|---|---|---|---|---|
| 21 x 9 | 696 | -2.5000 % | -2.5000 % | -1.4508 % |
| 41 x 17 | 2 672 | -1.2500 % | -1.2500 % | -0.8194 % |
| 81 x 33 | 10 464 | -0.6250 % | -0.6250 % | -0.4911 % |
| 161 x 65 | 41 408 | -0.3125 % | -0.3125 % | -0.3070 % |
| 321 x 129 | 164 736 | -0.1563 % | -0.1563 % | -0.1913 % |

First order, and it is the mesh that is short, not the closed form (which is at
`1e-14`). The default grid therefore *warns* — `planing_volume` and
`planing_waterplane` at `-1.25 %` against a 1 % band — on such a hull, and that warning
is the honest reading: a planing bottom drawn with `q < 0.1` needs a finer grid than
the 41 x 17 default before its hydrostatics mean 1 %. `q = 0.6-1.0` (the shapes
`docs/planing_hull.md` tabulates) converge at the usual order and warn on nothing.

## 7. Limitations of this audit

* The ground truth is a high-precision evaluation of the **same paper integral**. If
  the integrand itself is the wrong physical model (e.g. the topside strip should be
  cambered, or the transom should be a loft rather than a fan), no amount of precision
  finds it. This audit settles arithmetic, not modelling.
* The two reference routes share the closed-form `b(x)`, `x0` and `section_area` with
  the module — deliberately, so that the comparison isolates the quadrature — but that
  means an error inside `section_area` would be invisible to all three routes. It was
  checked separately against a hand derivation (V plus wall moment: `-T**2 b + T tau b**2
  - (1/3) tau**2 b**3`), which is what the module's own `Simpson` route
  (`simpson_section_integrals`) also does.
* The sweep is a grid, not a proof. `q = 0.05` and `q = 4.0` are the declared bounds of
  `PlaningHullParameters`; the pathological corner is `q -> 0.05`, `deadrise -> 59 deg`,
  `draft -> 0.05 m` and it is in the sweep.
* Runtime of the adaptive rule is comparable to the fixed 4000-point Simpson it
  replaces (`tools/planing_topside_robustness.py`: 420 cases x 3 routes in ~20 s), so
  the accuracy was not bought with a slow path — but it is no longer a fixed cost:
  a hull with a very fine bow costs more panels than a straight-chine one.

## 8. Where this leaves `docs/planing_hull.md`

Its section 2 table called the topside a quadrature, and its section 4 "anchors that do
not depend on the integrals" were all volume/waterplane/transom anchors — none of them
touched the wetted area, so that page's claims stand. What it did not say is that the
wetted area was the one quantity in the hydrostatics block with a 0.85 %-wide question
mark on it; it now has 2e-14.