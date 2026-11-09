# Wigley initial stability: closed forms, two mesh routes, and the GM amplification

**Module** `tensorcad/wigley.py` (`analytic_waterplane_inertia`,
`analytic_waterplane_inertia_longitudinal`, `analytic_metacentric_radius_transverse`,
`analytic_metacentric_radius_longitudinal`, `analytic_initial_stability`,
`WigleyStability`). **Evidence** `artifacts/wigley_stability/` — `results.json`,
`report.md`, `fig_convergence.png`, `fig_amplification.png`. **Generator**
`python examples/wigley_stability/run_example.py` (exits non-zero if a claim stops
holding). **Tests** `tests/test_wigley_stability_verification.py` (17 cases).

## 1. Case: which verification debt this closes

`tensorcad.wigley` has always been a verification case as much as a hull form: volume,
waterplane area, block/prismatic/midship coefficients and `VCB` are exact rationals
(`docs/wigley.md`). The **waterplane second moments were not**. That left the
metacentric chain (`BM = I/∇`, `KM = KB + BM`, `GM = KM − KG`) anchored only on the
curved `tensorcad.hull` family, whose `I_T` reference is a Beta integral of a
`sin(πx/L)^0.7` planform (`docs/stability.md`) — a second family, but a different
integrand, and no anchor at all for a hull whose planform is **polynomial** at both ends
and whose proportions are the ones the resistance/diffraction benchmarks actually use.

That matters because the Wigley reference case is the geometry a boundary-element or
LBM study is set up on: a solver user who needs `GM` (or a `KG` limit) for that hull
had to mesh it and trust a convergence study. This page derives the closed forms and
measures both merged mesh routes against them.

| item | value |
|---|---|
| family | `y = (B/2)(1 − \|2x/L − 1\|^n)(1 − (w/T)^m)`, `w` measured down from the WL |
| reference case | ITTC proportions `L = 10 m`, `B = 1 m`, `T = 0.625 m` (`L/B = 10`, `B/T = 1.6`) |
| exponents | `n = 2, 3` on the ladder; `(n, m) ∈ {(2,2), (3,2), (2,3), (3,3), (2.5,2)}` closed-form-only |
| loading condition | `KG = 0.5 m` above the keel for the headline table, ladder `0.30 … 0.527 m` |
| routes compared | plane cut (`tensorcad.stability`), wetted-surface flux (`tensorcad.hydrostatics`) |

## 2. Method: the closed forms, derived here and not read back from the integrator

The waterplane is `|y| ≤ b(x)`, `b = (B/2)(1 − |u|^n)`, `u = 2x/L − 1`.

| quantity | derivation | closed form | at `n = 2` |
|---|---|---|---|
| `I_T` | `∫₀^L (2/3) b³ dx`, `dx = (L/2)du` | `(B³L/12)·K(n)`, `K = 1 − 3/(n+1) + 3/(2n+1) − 1/(3n+1)` | `4LB³/105` |
| `I_L` | `∫ (x − L/2)² 2b dx`, `s = (L/2)u` | `(BL³/4)(1/3 − 1/(n+3))` | `BL³/30` |
| `BM_T` | `I_T/∇`, `∇ = LBT·nm/((n+1)(m+1))` | `B² (n+1)(m+1) K(n) / (12 T n m)` | `3B²/(35T)` |
| `BM_L` | `I_L/∇` | `L² (n+1)(m+1)(1/3 − 1/(n+3)) / (4 T n m)` | `3L²/(40T)` |
| `KB` | existing `VCB` closed form | `T(m+1)/(2(m+2))` | `5T/8` |

Two things a mesh convergence study would never reveal, and which is why the derivation
is written out: the `I_T` integral carries **`3n`** (the half-breadth is cubed, so the
Beta integral is `3n` and not `n`), and `I_L` is about the **centre of flotation**, so the
parallel-axis term is already removed — the same convention `waterplane_moments` reports.

Three independent references back the claims:

1. **Hand rationals** at `n = 2` (`4LB³/105`, `BL³/30`, `3B²/(35T)`, `3L²/(40T)`), written
   out in the test file.
2. **Composite Simpson** on the analytic half-breadth (20000 panels), including
   `n = 2.5`, where Simpson is no longer integrating a polynomial and therefore measures
   real quadrature error rather than round-off.
3. **Two mesh routes that share no code**: the plane cut with Green's theorem on closed
   loops (`tensorcad.stability.waterplane_moments`) and the divergence theorem over the
   **wetted surface**, which needs no waterplane cap at all
   (`tensorcad.hydrostatics::_waterplane_inertia`). Both are exact on a polygonal
   waterplane, so they must agree to round-off; whatever is left is the mesh's chord
   error, not the integrator's.

## 3. Results

### 3.1 Closed forms (`KG = 0.5 m`; `m` moves `KB` and `∇`, not the inertias)

| n | m | I_T [m⁴] | I_L [m⁴] | KB [m] | BM_T [m] | BM_L [m] | KM_T [m] | GM_T [m] |
|---|---|---|---|---|---|---|---|---|
| 2 | 2 | 0.380952381 | 33.333333 | 0.390625 | 0.137142857 | 12.000000 | 0.527767857 | +0.027768 |
| 3 | 2 | 0.482142857 | 41.666667 | 0.390625 | 0.154285714 | 13.333333 | 0.544910714 | +0.044911 |
| 2 | 3 | 0.380952381 | 33.333333 | 0.375000 | 0.121904762 | 10.666667 | 0.496904762 | **−0.003095** |
| 3 | 3 | 0.482142857 | 41.666667 | 0.375000 | 0.137142857 | 11.851852 | 0.512142857 | +0.012143 |
| 2.5 | 2 | 0.437675070 | 37.878788 | 0.390625 | 0.147058824 | 12.727273 | 0.537683824 | +0.037684 |

The `(2,3)` row is deliberately kept: a fuller keel section (`m = 3`) drops `KB` from
`5T/8` to `3T/8`, and at `KG = 0.5 m` the hull is then **unstable upright** (`GM < 0`).
The closed form reports the sign; nothing clamps it.

Worst closed-form-vs-Simpson disagreement over the five exponent pairs: **3.5e-16**
relative (integer exponents, where Simpson is exact to round-off) and within **1e-9** for
`n = 2.5`. The `n = 2` rationals are reproduced to **2.1e-16**.

### 3.2 The ladder: the merged mesh routes converge to the closed forms

`n = 2` (classic Wigley), `KG = 0.5 m`; errors signed against the closed form:

| grid | I_T [m⁴] | err I_T | err I_L | err KB | err BM_T | err KM_T | err GM_T | route gap |
|---|---|---|---|---|---|---|---|---|
| 41×13 | 0.380397173 | −0.1457 % | −0.1041 % | +0.0348 % | +0.0905 % | +0.0493 % | **+0.9361 %** | 4.3e-15 |
| 81×25 | 0.380813514 | −0.0365 % | −0.0260 % | +0.0087 % | +0.0226 % | +0.0123 % | +0.2337 % | 1.0e-15 |
| 161×41 | 0.380917660 | −0.0091 % | −0.0065 % | +0.0031 % | +0.0104 % | +0.0050 % | +0.0954 % | 2.6e-15 |
| 321×81 | 0.380943700 | −0.0023 % | −0.0016 % | +0.0008 % | +0.0026 % | +0.0013 % | **+0.0239 %** | 2.0e-14 |

`n = 3`, same loading:

| grid | err I_T | err I_L | err KB | err BM_T | err KM_T | err GM_T |
|---|---|---|---|---|---|---|
| 41×13 | −0.1458 % | −0.1873 % | +0.0348 % | +0.1113 % | +0.0565 % | +0.6850 % |
| 81×25 | −0.0365 % | −0.0469 % | +0.0087 % | +0.0278 % | +0.0141 % | +0.1710 % |
| 161×41 | −0.0091 % | −0.0117 % | +0.0031 % | +0.0117 % | +0.0056 % | +0.0675 % |
| 321×81 | −0.0023 % | −0.0029 % | +0.0008 % | +0.0029 % | +0.0014 % | +0.0169 % |

Four things are read off this table rather than asserted:

* **`I_T` is under-predicted at every rung** (`−`), which is what a chord approximation
  of a convex taper must do; the sign is a prediction, and it is checked on all eight
  rungs. The two routes agree to **1e-15…1e-14** relative, i.e. round-off: the
  disagreement in the table is the **mesh**, measured identically by both integrators.
* **The order is 2.00** (fitted on `log|err I_T|` vs `log 1/(nx−1)`; quarters per
  doubling). That is the signature of a piecewise-linear waterline and is asserted as a
  band `1.8 … 2.2`, not as a slogan.
* `KB` errs **the other way** (`+`): the lofted volume sits inside the analytic surface,
  so the buoyancy centroid is a little high. `BM_T = I_T/∇` therefore partly cancels —
  `I_T` low and `∇` low — which is why `KM` is more accurate than `I_T` (0.049 % vs
  0.146 % at the default grid). A report that only quoted `I_T` would have overstated the
  error of the quantity that is actually used.
* **`GM` is the worst number on the page** and is printed as such: +0.94 % at the default
  grid, 19× the `KM` error. See §3.3.

### 3.3 The amplification, measured against the loading condition

`GM = KM − KG` is a difference of nearly equal numbers once `KG` approaches `KM`, so the
same `KM` error becomes a large `GM` error. Default 41×13 grid:

| KG [m] | GM ref [m] | GM mesh [m] | rel GM | rel KM | KM/GM |
|---|---|---|---|---|---|
| 0.300 | +0.227768 | +0.228028 | +0.114 % | +0.0493 % | 2.3 |
| 0.450 | +0.077768 | +0.078028 | +0.334 % | +0.0493 % | 6.8 |
| 0.500 | +0.027768 | +0.028028 | +0.936 % | +0.0493 % | 19.0 |
| 0.520 | +0.007768 | +0.008028 | +3.346 % | +0.0493 % | 67.9 |
| 0.527 | +0.000768 | +0.001028 | **+33.854 %** | +0.0493 % | **687.3** |

The identity `rel(GM) = (KM/GM)·rel(KM)` holds exactly (it is algebra, checked to 1e-9),
so the table is a *loading-condition* warning, not a numerical one: the same mesh that is
comfortably inside 3 % at `KG = 0.5 m` is 34 % wrong at `KG = 0.527 m`. Designs that need
`GM` near the limit must refine the grid or use the closed form — which is now available.

Figures: `artifacts/wigley_stability/fig_convergence.png` (relative error of `I_T`,
`I_L` and `GM_T` against the grid step, both exponents, with a slope-2 guide) and
`artifacts/wigley_stability/fig_amplification.png` (left: `GM` falling one-for-one with
`KG`; right: measured `rel(GM)/rel(KM)` against the predicted `KM/GM`, log scale). They
are regenerated by the example and are not duplicated into `docs/`.

## 4. Error judgement against the 3 % standard

| key quantity | case | reference | error | < 3 %? |
|---|---|---|---|---|
| `I_T` | default 41×13, `n = 2` | closed form `4LB³/105` | −0.1457 % | **yes** |
| `I_T` | 321×81, `n = 2` | closed form | −0.0023 % | **yes** |
| `I_L` | default 41×13, `n = 2` | closed form `BL³/30` | −0.1041 % | **yes** |
| `KB` | default 41×13 | `5T/8` | +0.0348 % | **yes** |
| `BM_T` | default 41×13 | `3B²/(35T)` | +0.0905 % | **yes** |
| `KM_T` | default 41×13 | `0.527767857 m` | +0.0493 % | **yes** |
| `GM_T` | default 41×13, `KG = 0.5` | `+0.027767857 m` | +0.9361 % | **yes** |
| `GM_T` | default 41×13, `KG = 0.527` | `+0.000768 m` | +33.854 % | **no — and it is reported, not excused** |
| closed forms | 5 exponent pairs | Simpson | ≤ 3.5e-16 (≤ 1e-9 at `n = 2.5`) | **yes** |

Every quantity the module claims to compute is inside 3 % at the default grid **for a
`KG` that leaves a usable `GM`**; the `KG = 0.527` row is outside it and is printed as a
measurement rather than promoted to a claim, which is the same policy
`docs/stability.md` applies to its coarse rung.

## 5. Notes, limits, and what is deliberately absent

* **Upright, intact, untrimmed only.** No `GZ` curve, no heel, no trim, no free-surface
  correction, no damage. The free-surface correction is a separate model
  (`docs/free_surface.md`, applied as a `KG` raise in `tensorcad.gz`), and damage
  stability is `tensorcad.damage` (lost buoyancy).
* **`KG` is an input, never a default.** `analytic_initial_stability` raises without it,
  for the same reason `initial_stability` does: a default of zero silently reports
  `GM = KM`.
* **The two mesh routes are independent of each other but not of the mesh.** Their
  1e-15 agreement says the integrators are exact on the polygon they are given; it says
  nothing about the polygon. That is what the closed form is for.
* **`m` (vertical exponent) never enters the inertias** — the waterline is at `w = 0`,
  where the section factor is `1` by construction. `m` moves `∇` and `KB` and hence `BM`;
  the test pins that so a future edit cannot quietly couple them.
* **Exponents outside `n, m ≥ 1`** are refused by the closed forms (`ValueError`); the
  loft refuses them too (`WigleyParameters.validate`), but for a different reason (a
  finite waterline corner), which is why both checks exist.
* The `(2,3)` closed-form row has a negative `GM` at `KG = 0.5 m`: this is a property of
  the family, not an error, and no value is clamped.

## 6. Reproducing

```
python examples/wigley_stability/run_example.py            # artifacts/wigley_stability/
python -m pytest tests/test_wigley_stability_verification.py -q
python -m pytest tests/test_wigley.py tests/test_stability.py -q   # nothing else moved
```