# Initial stability measured against closed forms

TensorCAD 0.22.0 · `examples/stability/run_example.py` · hull 10.0 × 2.0 × 0.8 m, KG = 0.6 m above keel

## Claims enforced by this run

| claim | measured | threshold | holds |
|---|---|---|---|
| box: every waterplane and metacentric quantity is exact to round-off | 2.220e-16 | 1e-12 | PASS |
| hemisphere: KM reproduces the sphere radius R (KB = 5R/8, BM = 3R/8) | 3.208e-03 | 2e-02 | PASS |
| default hull grid (41×13): I_T within 3 % of the closed form | 1.386e-03 | 3e-02 | PASS |
| default hull grid (41×13): I_L within 3 % of quadrature | 6.038e-03 | 3e-02 | PASS |
| default hull grid (41×13): GM_T within 3 % of the analytic surface | 1.986e-02 | 3e-02 | PASS |
| fine hull grid (161×41): I_T within 1 % | 8.604e-05 | 1e-02 | PASS |
| fine hull grid (161×41): GM_T within 1 % | 3.470e-03 | 1e-02 | PASS |
| the waterplane error is second order in the grid step (fitted order > 1.5) | 2.006e+00 | 2e+00 | PASS |
| the I_L reference is rule-independent (Simpson vs tanh-sinh within 1e-6) | 7.772e-16 | 1e-06 | PASS |
| the waterplane section is symmetric: |y_c| is round-off (below 1e-9 m) | 6.747e-18 | 1e-09 | PASS |

## Hull ladder: mesh waterplane moments vs the analytic surface

| nx × nz | A [m²] | err A | I_T [m⁴] | err I_T | I_L [m⁴] | err I_L | BM_T [m] | KB [m] | GM_T [m] | err GM |
|---|---|---|---|---|---|---|---|---|---|---|
| 21×9 | 14.0062 | -0.568% | 3.2524 | -0.559% | 75.4882 | -1.799% | 0.4410 | 0.4853 | 0.3263 | +3.473% |
| 41×13 | 14.0617 | -0.175% | 3.2662 | -0.139% | 76.4066 | -0.604% | 0.4386 | 0.4830 | 0.3216 | +1.986% |
| 81×25 | 14.0787 | -0.054% | 3.2696 | -0.034% | 76.7187 | -0.198% | 0.4366 | 0.4811 | 0.3177 | +0.736% |
| 161×41 | 14.0840 | -0.017% | 3.2704 | -0.009% | 76.8217 | -0.064% | 0.4359 | 0.4805 | 0.3165 | +0.347% |
| 321×81 | 14.0856 | -0.005% | 3.2707 | -0.002% | 76.8551 | -0.020% | 0.4356 | 0.4802 | 0.3158 | +0.125% |

Analytic reference: A = 14.086309 m², I_T = 3.270725 m⁴, I_L = 76.870781 m⁴, ∇ = 7.512698 m³, BM_T = 0.435360 m, KM_T = 0.915360 m.

Fitted orders (log–log slope of |error| against the grid step): waterplane_area 1.70, transverse_second_moment 2.01, longitudinal_second_moment 1.62, displacement_volume 1.32.

`I_L` is quadrature, not elementary: the tanh-sinh value agrees with an independent Simpson rule after the `u = 1/2 − s²` substitution to -7.77e-16 relative, and with tanh-sinh one level up.

## Figures

- `sections.png` — the waterplane chords the mesh actually has against the analytic planform.
- `convergence.png` — |error| vs grid step for A, I_T, I_L and the displacement volume.
- `gm_ladder.png` — GM_T against KG, with the GM = 0 loading limit at KG = KM.
