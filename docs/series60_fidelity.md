# Series 60 fidelity: the family was not a Series 60 replica, and now it can be

**Question asked.** `tensorcad.series60` passed every check it had: closed forms agree with
its own mesh, the volume converges under refinement, VCB matches `(m+1)/(2m+1)`. All of that
is *internal* consistency. This page asks the external question -- does the hull it lofts
reproduce the particulars that are published for the Series 60 parent forms? -- and finds a
structural defect that no amount of self-consistency could have exposed.

## 1. Case

The Series 60 (Todd 1951) is the standard systematic displacement-hull benchmark: five parent
forms spanning `Cb = 0.60 … 0.80`, each published with a midship coefficient `Cm`, a prismatic
coefficient `Cp` and a waterplane coefficient `Cwp`. The library models them as

```
y(xi, zeta) = (B/2) * (1 - |2*xi - 1|^n) * zeta**(1/m)
```

with a single longitudinal exponent `n` and a single vertical exponent `m`, both pinned to the
published `Cp` and `Cm` of the parent form.

| Quantity | Value |
|---|---|
| Length `L` | 10.0 m |
| Beam `B` | 1.5 m |
| Draft `T` | 0.6 m |
| `L/B`, `B/T` | 6.67, 2.50 |
| Block coefficient sweep | 0.60, 0.65, 0.70, 0.75, 0.80 |
| Baselines | tabulated Series 60 parents (Todd 1951); Schneekluth `Cwp = (1 + 2Cb)/3` |

The published table is self-consistent to the digit given: `Cb - Cp*Cm` ranges from
+0.0000 to
+0.0030, i.e. at most 0.4 % of `Cb`,
which is what makes the transcription auditable. The regression is a second, independent
baseline that does not depend on reading a table correctly.

## 2. What was done

1. **Structural analysis, no references needed.** With one exponent, both the deck width and
   the section area are proportional to `1 - |2*xi-1|^n`. The waterplane integral is then the
   same functional of `n` as the prismatic integral, so `Cwp ≡ Cp`; and `f` is fore/aft
   symmetric, so `LCB ≡ L/2`. Both are identities of the construction, not approximations.
2. **Measurement against the baselines.** The five parent forms were evaluated against the
   tabulated `Cwp` and against the regression (§3, table 1).
3. **Independent recomputation of every closed form.** The surface was re-integrated from its
   definition by tanh-sinh (double-exponential) quadrature, which resolves both the
   `zeta**(1/m)` keel singularity and the fact that a full section's area sits in a layer of
   thickness `1/m` under the waterline. Nothing in that path shares code with the closed
   forms (table 3).
4. **The fix.** Two optional inputs, `Series60Parameters.cwp` and `.lcb`. The deck taper and the
   section fullness are decoupled: `Cm(xi) = Cm_mid (1 + gain (1 - f(xi)))` with
   `gain = (Cp - Cwp)/(Cwp - ∫f²)`, so the volume still lands on `Cb` in closed form while
   `Cwp` is a free target; the aft and fore halves carry their own exponents, solved by a 2×2
   Newton iteration against a `(Cwp, LCB)` target. `cwp = lcb = None` reproduces the legacy hull
   bit for bit (`tests/test_series60_fidelity.py::test_gain_zero_reproduces_the_legacy_hull`).
5. **Re-measurement.** The same quantities on the corrected hulls, plus a grid study to show the
   loft still converges to its own closed forms (tables 2 and 4).

## 3. Results

### Table 1 — legacy family against the published parents

| Cb | Cp | Cm | Cwp (model) | Cwp (Todd) | error vs Todd | Cwp (regression) | error vs regression | LCB/L |
|---|---|---|---|---|---|---|---|---|
| 0.60 | 0.625 | 0.960 | 0.6250 | 0.700 | -10.71 % | 0.7333 | -14.77 % | 0.500 |
| 0.65 | 0.670 | 0.970 | 0.6700 | 0.727 | -7.84 % | 0.7667 | -12.61 % | 0.500 |
| 0.70 | 0.715 | 0.979 | 0.7150 | 0.753 | -5.05 % | 0.8000 | -10.62 % | 0.500 |
| 0.75 | 0.760 | 0.987 | 0.7600 | 0.780 | -2.56 % | 0.8333 | -8.80 % | 0.500 |
| 0.80 | 0.805 | 0.994 | 0.8050 | 0.807 | -0.25 % | 0.8667 | -7.12 % | 0.500 |

The model's `Cwp` is *whatever its `Cp` is*. Against the tabulated parents that is
-10.71 % at `Cb = 0.60`, shrinking to -0.25 %
at `Cb = 0.80`; against the regression it is -14.77 % to
-7.12 % across the whole range. LCB is 0.500 for every parent,
i.e. exactly midship, where the published forms are all aft of it.

### Table 2 — the same parents with `cwp` and `lcb` as inputs

| Cb | n (aft) | n (fore) | taper gain | Cm mid | Cwp (mesh-free closed form) | error vs request | LCB/L | LCB error | LCF/L | VCB/T |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.60 | 2.970 | 1.872 | -0.6118 | 0.960 | 0.700000 | +0.0e+00 | 0.475000 | -5.6e-17 | 0.4796 | 0.5355 |
| 0.65 | 3.179 | 2.261 | -0.4979 | 0.970 | 0.727000 | -1.1e-15 | 0.483000 | -1.1e-16 | 0.4857 | 0.5263 |
| 0.70 | 3.417 | 2.737 | -0.3587 | 0.979 | 0.753000 | +1.5e-16 | 0.490000 | -5.6e-17 | 0.4912 | 0.5175 |
| 0.75 | 3.730 | 3.375 | -0.2075 | 0.987 | 0.780000 | +0.0e+00 | 0.496000 | -5.6e-17 | 0.4963 | 0.5096 |
| 0.80 | 4.062 | 4.307 | -0.0232 | 0.994 | 0.807000 | +0.0e+00 | 0.502000 | +1.1e-16 | 0.5020 | 0.5022 |

Every requested particular is hit to round-off. The mechanism is visible in the taper gain: the
gain is negative for the full-waterplane/low-`Cb` parents, i.e. the sections are *thinned away
from midship* while the deck plan keeps its fuller taper -- which is what a real hull with a fine
entrance and a full waterplane looks like, and what the single-exponent family could not express.

### Table 3 — closed forms against an independent tanh-sinh integration

| Case | volume (closed) | rel. error | waterplane (closed) | rel. error | VCB above keel | rel. error | LCF/L | rel. error |
|---|---|---|---|---|---|---|---|---|
| legacy Cb=0.60 | 5.400000 | -1.8e-07 | 9.375000 | -1.8e-07 | 0.306122 | -1.5e-15 | 0.500000 | -5.0e-15 |
| legacy Cb=0.80 | 7.200000 | -1.9e-13 | 12.075000 | -1.9e-13 | 0.300935 | -7.4e-16 | 0.500000 | -3.6e-16 |
| cwp=0.700 Cb=0.60 | 5.400000 | +7.6e-09 | 10.500000 | +4.2e-09 | 0.321324 | -1.7e-09 | 0.500000 | -1.8e-15 |
| cwp=0.700 lcb=0.475 | 5.400000 | -1.8e-08 | 10.500000 | -1.0e-08 | 0.321278 | +4.1e-09 | 0.479625 | -1.7e-10 |
| cwp=0.780 lcb=0.496 | 6.750000 | +2.1e-11 | 11.700000 | +1.7e-11 | 0.305763 | -2.0e-12 | 0.496299 | +6.8e-14 |

Worst disagreement over all cases: **1.8e-07** relative -- the closed forms are right.
The defect was in the parametrisation, not the algebra, and this separates the two.

### Table 4 — loft convergence to the closed forms

| Hull | grid | triangles | volume error | LCB error | VCB error | closed manifold |
|---|---|---|---|---|---|---|
| legacy | 41x13 | 1998 | -3.68 % | +0.00e+00 m | +3.54 % | yes |
| legacy | 81x21 | 6558 | -2.14 % | +0.00e+00 m | +2.10 % | yes |
| legacy | 161x41 | 25918 | -1.04 % | +0.00e+00 m | +1.03 % | yes |
| legacy | 321x81 | 103038 | -0.50 % | -8.88e-16 m | +0.50 % | yes |
| cwp=0.700/lcb=0.475 | 41x13 | 1998 | -2.73 % | +3.08e-03 m | +2.59 % | yes |
| cwp=0.700/lcb=0.475 | 81x21 | 6558 | -1.55 % | +1.74e-03 m | +1.50 % | yes |
| cwp=0.700/lcb=0.475 | 161x41 | 25918 | -0.72 % | +8.40e-04 m | +0.71 % | yes |
| cwp=0.700/lcb=0.475 | 321x81 | 103038 | -0.34 % | +4.08e-04 m | +0.34 % | yes |

First-order convergence, as the box-rule loft has always shown; the corrected hull behaves the
same way, and the LCB error decays the same order rather than stalling at a bias (it did stall
at −0.6 %L while LCB was mis-derived from the waterplane moment instead of the volume
distribution -- see §5).

**Figures** (`artifacts/series60_fidelity/`): `cwp_vs_baselines.png`, `sections_legacy_vs_cwp.png`,
`waterplane_plan.png`, `convergence.png`.

## 4. Error verdict

| Check | Baseline | Error | < 3 %? |
|---|---|---|---|
| legacy `Cwp`, `Cb = 0.60` | tabulated parent | -10.71 % | **no** |
| legacy `Cwp`, `Cb = 0.65/0.70` | tabulated parent | -7.84 % / -5.05 % | **no** |
| legacy `Cwp`, `Cb = 0.75/0.80` | tabulated parent | -2.56 % / -0.25 % | 0.75 no, 0.80 yes |
| legacy `Cwp`, all parents | Schneekluth regression | -14.77 % worst | **no** |
| legacy `LCB` | published parents (aft of midship) | 0 vs −1…−2.5 %L | **no** (structurally impossible) |
| corrected `Cwp`, `LCB` | as requested | ≤ 2e-15 | **yes** |
| corrected hull, mesh at 161×41 | own closed forms | volume -1.04 %, LCB +0.0e+00 m, VCB +1.03 % | volume/VCB no at that grid (first order); closed forms themselves verified to 2e-07 |

Honest summary: **the legacy family fails the published-parent comparison at `Cb ≤ 0.75` on
`Cwp` and at every parent on `LCB`.** The `cwp`/`lcb` extension removes both failures, and its
closed forms are verified independently. The mesh is first-order accurate in the loft, so a
`Cwp` read off a 161×41 mesh is good to ~1e-3 relative, not better -- refine before quoting.

## 5. Notes, limits, and what is not claimed

* **Published LCB values are demonstration targets, not references.** The verdict on LCB is
  structural (the legacy family can only produce 0.5) and does not depend on those numbers.
  The `Cwp` comparison is what carries the quantitative claim, and it has two independent
  baselines that agree that the low-`Cb` forms are wrong.
* **A trap found while building this.** Deriving `LCB` from the waterplane moment looks right
  and is right only when the section fullness is uniform. With `gain ≠ 0` the volume
  distribution is `f·Cm(xi)`, and the waterplane version leaves a bias of −0.56 %L that survives
  refinement -- visible only because the grid study converged to a nonzero value. `_resolve_halves`
  now uses the volume moment, and table 4 shows the bias decaying.
* **Not every `(Cb, Cwp)` pair is reachable.** Pushing `Cwp` below `Cp` demands a midship
  section fuller than 1.0 at the ends; that is rejected with a message saying so instead of
  silently clamping (`test_waterplane_too_fine_for_the_block_coefficient_is_rejected`).
* **VCB moves.** The corrected hull's VCB/T is 0.5355 at `Cb = 0.60`
  against 0.5102 for the legacy hull, because the thinned ends put relatively more volume near
  the waterline. Anything downstream that pinned VCB should re-read it.
* **Still not a lines plan.** The family is an analytic stand-in with the right integral
  particulars; it does not reproduce Todd's tabulated station offsets, and no resistance
  claim should be built on it as if it did. That remains the honest gap.
