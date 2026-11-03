# Hydrostatics verified against a closed-form reference (0.11)

`tensorcad.hydrostatics` integrates the property of a *mesh*, so every number it
reports carries the mesh's discretisation error. Until 0.11 only one of those numbers
had an independent reference in the tree: the displaced volume, against
`tensorcad.hull.analytic_volume`. The coefficients (`Cb`, `Cp`), the centres (`LCB`,
`LCF`, `VCB`), the waterplane area, the maximum section area and the wetted area had
nothing -- `tests/test_hydrostatics.py` only bounded them
(`block_coefficient < 0.5 < prismatic_coefficient`). A one-percent shift in `Cb` would
have passed the suite.

`tensorcad.hydrostatics_analytic` supplies the missing reference, for the surface the
generator samples:

```
b(x, zeta) = (B/2) * sin(pi*x/L)**e * zeta**0.5      zeta = (z + T)/T in [0, 1]
```

Because the breadth taper and the vertical profile are separated, almost every
hydrostatic property of that surface is a closed form -- not a numerical
approximation of one.

## The reference in closed form

| Quantity | Exact value | Form |
|---|---|---|
| displaced volume | `B*L*T*I*2/3` | 7.512698 m³ for L/B/T = 10/2/0.8 |
| displacement mass | volume × ρ (1025 kg/m³) | 7.70052e3 kg |
| waterplane area | `B*L*I` | 14.086309 m² |
| block coefficient `Cb` | `2*I/3` | 0.469544 |
| prismatic coefficient `Cp` | `I` | 0.704315 |
| midship coefficient `Cm` | `2/3` | 0.666667 |
| waterplane coefficient `Cw` | `I` | 0.704315 |
| maximum section area | `2*B*T/3`, at `x = L/2` | 1.066667 m² |
| `LCB`, `LCF` | `L/2` (the form is fore-aft symmetric) | 5.0 m |
| `VCB` | `-0.4*T` (i.e. `0.6*T` above the keel) | -0.320000 m |
| wetted surface area | `2 * ∫∫ sqrt(T²(1 + b_x²) + b_zeta²) dx dzeta` | 22.912730 m² |

with

```
I = int_0^1 sin(pi*x)**e dx = sqrt(pi) * Gamma((e+1)/2) / (pi * Gamma((e+2)/2))
  = 0.7043154582422714        for e = 0.7
```

`I` is derived here from the gamma function, whereas `tensorcad.hull.TAPER_INTEGRAL`
was derived from the volume identity the generator samples. The two agree to
**0.0 ulp** (`taper_integral() - TAPER_INTEGRAL == 0.0`), so the volume identity the
library has shipped since 0.1 is confirmed by an independent derivation -- pinned in
`test_taper_integral_agrees_with_the_volume_derived_constant`.

`VCB = -0.4*T` follows from the z-weighted mean of `zeta**0.5` on `[-T, 0]`:
`(int z*zeta**0.5 dz) / (int zeta**0.5 dz) = (-4T/15)/(2/3) = -0.4*T`. `LCB = LCF = L/2`
follows because `sin(pi*x/L)**e` is symmetric about midship and the section shape does
not depend on `x`.

Only the wetted area needs quadrature, and it needs one that survives two integrable
endpoint singularities: `b_zeta ~ zeta**-0.5` at the keel and
`b_x ~ xi**(e-1) = xi**-0.3` at the stem and the stern. `tanh_sinh_quadrature` is double
exponential, so those are the cases it converges on exponentially fast; the
substitution `zeta = u**2` removes the keel singularity into a smooth integrand first.
The returned value carries its own evidence: the difference to the previous refinement
is `5.0e-13` relative.

## Mesh against the reference

`examples/hydrostatics_verification/run_verification.py` runs a five-level ladder with
both axes refining together, so the nominal cell size `h = L/(nx-1)` is the abscissa
and the fitted slopes are order estimates rather than per-axis rates. Signed errors,
in percent of the reference:

| nx × nz | vertices | `V` | `Awp` | `S` | `Amax` | `VCB` | `Cb` | `Cp` | `LCB` | `LCF` | worst |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 21 × 9 | 341 | -1.8416 | -0.5684 | -0.5487 | -1.2805 | -1.6449 | -1.8416 | -0.5684 | 0.0000 | 0.0000 | 1.8416 |
| 41 × 13 | 1001 | -0.8804 | -0.1749 | -0.3261 | -0.7067 | -0.9391 | -0.8804 | -0.1749 | 0.0000 | 0.0000 | 0.9391 |
| 81 × 21 | 3281 | -0.3867 | -0.0538 | -0.1730 | -0.3330 | -0.4566 | -0.3867 | -0.0538 | 0.0000 | 0.0000 | 0.4566 |
| 161 × 41 | 12961 | -0.1359 | -0.0166 | -0.0742 | -0.1194 | -0.1685 | -0.1359 | -0.0166 | 0.0000 | 0.0000 | 0.1685 |
| 321 × 81 | 51521 | -0.0477 | -0.0051 | -0.0315 | -0.0426 | -0.0613 | -0.0477 | -0.0051 | -0.0000 | 0.0000 | 0.0613 |

Every quantity is under-integrated: the mesh is an inscribed polygonal surface, so it
can only lose volume, area and area-moment. Nothing over-shoots, on any grid.

Observed orders from the last two refinements (log-log slope of the error):

| Quantity | order at 321×81 | note |
|---|---|---|
| `V`, `Cb` | 1.511 | dominated by the `xi**0.7` ends, not by a smooth-surface `O(h²)` |
| `Awp`, `Cp` | 1.700 | the waterline is a one-dimensional integral of the same taper |
| `S` | 1.235 | the wetted surface also inherits the `zeta**0.5` keel edge |
| `Amax` | 1.486 | section-area integral at one station |
| `VCB` | 1.458 | a ratio of two volume-weighted integrals |
| `LCB`, `LCF` | — | exactly `L/2` at every resolution: fixed by symmetry, not by mesh density |

The sub-quadratic order is a property of the *form*, not a defect of the integrator:
`sin(pi*x/L)**0.7` has an unbounded derivative at both ends, so a trapezoidal-like
sum cannot do better than about `h**1.7`. This is the same reason
`estimate_discretisation_error` exists instead of a promise of second-order accuracy.

## Three independent estimates of the wetted area

The wetted area is the only reference that came out of quadrature, so it is
cross-checked against the mesh series itself:

| Estimate | Value [m²] | Method |
|---|---|---|
| quadrature | 22.912730 | nested tanh-sinh, `zeta = u²` substitution |
| Richardson extrapolation of the mesh series | 22.912932 | order 1.212 estimated from the 161×41 → 321×81 steps |
| agreement | **+8.8e-6 relative** | two estimates that share no code |

## Error judgement

| Claim | Value | Criterion | Verdict |
|---|---|---|---|
| worst quantity at the default 41×13 grid | 0.939 % | < 3 % (audit threshold) | **< 3 % — PASS** |
| worst quantity at 161×41 | 0.169 % | < 0.5 % | **PASS** |
| worst quantity at 321×81 | 0.061 % | < 0.1 % | **PASS** |
| volume at the default grid | -0.880 % | matches README/audit (-0.88 %) | **PASS** |
| `Cb`, `Cp` at the default grid | -0.880 %, -0.175 % | < 3 % | **< 3 % — PASS** |
| `VCB` at the default grid | -0.939 % | < 3 % | **< 3 % — PASS** |
| `LCB`, `LCF` | 0.0 (round-off) | exact | **PASS** |
| wetted area reference self-convergence | 5.0e-13 | < 1e-9 | **PASS** |
| quadrature vs Richardson | 8.8e-6 | < 2e-3 | **PASS** |

So: every hydrostatic quantity the library reports is within 3 % of an independent
closed form on the default grid, and within 0.17 % on the 161×41 grid the audit used.
The 3 % criterion is met with a factor of 3 to 50 of margin, and the remaining error is
documented as discretisation, not as a modelling assumption.

## What this does *not* verify

- **The reference describes this one analytic surface.** It is not a hull-form office
  and says nothing about a real ship. `LCB = LCF = L/2` holds because the form is
  fore-aft symmetric: those two rows verify the mesh integrator and would say nothing
  about the trim of an asymmetric hull (a gap logged in `backlog_tensorcad_section.md`
  as "multi-family / transom / asymmetric").
- **Stability, resistance, structure, class rules: not implemented, not verified.**
  `Cb`/`Cp`/`VCB` are inputs to those calculations, not substitutes for them;
  no metacentric height, no added resistance, no longitudinal strength.
- **The wetted area is for the bare lofted surface** -- no appendages, no roughness, no
  transom immersion. The SUBOFF appendage fields are separate geometry and are verified
  in `docs/suboff_appendages.md`.
- **No per-axis convergence rate.** Both grid axes refine together in the ladder, so a
  separate `h_x`/`h_z` study remains open; the orders above are composite.

## Reproduction

```bash
python examples/hydrostatics_verification/run_verification.py \
    --output artifacts/hydrostatics_verification      # ~9 s, exits non-zero on a broken claim
python -m pytest -q tests/test_hydrostatics_analytic.py   # 47 tests
```

The script writes `hydrostatics_verification.json` (ladder, per-quantity errors, fitted
orders, Richardson cross-check, and the five claims it enforces) plus
`convergence.png` and `offset_field_and_error.png`. Those are the only places the
numbers above come from; nothing in this page is typed by hand from a different run.
matplotlib is optional -- without it the JSON is still written and the figures are
skipped.

## Provenance

| Item | Value |
|---|---|
| Implementation | `tensorcad/hydrostatics_analytic.py` (0.11.0), re-exported as `tensorcad.analytic_hydrostatics`, `tensorcad.AnalyticHydrostatics`, `tensorcad.taper_integral`, `tensorcad.tanh_sinh_quadrature` |
| Reference surface | `tensorcad/hull.py` (`sin(pi*x/L)**0.7`, `sqrt(zeta)`), `TAPER_INTEGRAL` |
| Cross-check | `tensorcad.hull.analytic_volume` (volume), Richardson extrapolation (wetted area) |
| Tests | `tests/test_hydrostatics_analytic.py` (47 tests: closed forms, quadrature, mesh comparison, convergence, example claims) |
| Contract | `AnalyticHydrostatics.as_dict()` schema `tensorcad.hydrostatics.analytic/1` |
| Dependency | none beyond the standard library; NumPy and matplotlib are not imported |