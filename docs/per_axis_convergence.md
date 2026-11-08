# Per-axis convergence of the analytic hull: where the default grid's error lives

**TensorCAD** · tool `tools/per_axis_convergence.py` · artifact
`artifacts/per_axis_convergence/per_axis_convergence.json` · tests
`tests/test_per_axis_convergence.py`

This page closes the last item under "What this does *not* verify" in
`docs/hydrostatics_verification.md`: *"No per-axis convergence rate. Both grid axes
refine together in the ladder, so a separate `h_x`/`h_z` study remains open; the orders
above are composite."* It is now done, and the composite order that page reported
(1.212) turns out to be a mixture of two cleanly separable powers.

Every number below is written by the tool from the artifact; nothing is typed by hand.

## 1. Case

The surface is `tensorcad.hull`'s analytic loft, the same one the hydrostatics page
verifies:

| parameter | value |
|---|---|
| length `L` | 10.0 m |
| beam `B` | 2.0 m |
| draft `T` | 0.8 m |
| fore / aft taper exponent | 0.7 / 0.7 (`sin(pi x/L)**a`) |
| vertical taper | `sqrt(level/(nz-1))` |

The two grid axes do not enter the geometry the same way, so the error they produce
should not be expected to share an order:

* the **longitudinal** direction clips a `sin(a)` taper whose derivative is singular at
  both ends (a chord under-predicts a concave taper);
* the **vertical** direction clips a `sqrt(z/T)` taper, singular at the keel — the
  module's own docstring says this is deliberate ("keeps the sampled surface
  single-valued in z" and "makes the section derivative singular there").

Trapezoidal error on `f ~ t**a` near an endpoint behaves like `h**(a+1)`. That gives
pre-registered predictions **p_x = 1.7** and
**p_z = 1.5**, before any measurement.

## 2. Method

Three one-axis-at-a-time studies, all against `tensorcad.analytic_hydrostatics` (closed
forms), never against another mesh:

1. **axial ladder** — `nx` in `[11, 21, 41, 81, 161, 321]`, `nz` pinned at `81`;
2. **vertical ladder** — `nz` in `[5, 9, 13, 21, 41, 81]`, `nx` pinned at `321`;
3. **cross grids** — seven mid-ladder grids used to test whether the two errors *add*.

Orders are fitted on the **delta against the reference grid** (`E(nx, nz_ref) -
E(nx_ref, nz_ref)`), not on the total error. That distinction is the trap of this
study: once one axis is fine, the total error is floored by the other axis and the
fitted slope collapses (fitting total volume error over the axial ladder gives 1.08
instead of the true 1.8).

The additive model tested on the cross grids is

```
E_pred(nx, nz) = E_ref + [E(nx, nz_ref) - E_ref] + [E(nx_ref, nz) - E_ref]
```

**Toolchain.** Python 3.11, standard library only for the physics; matplotlib only for
the three figures and optional (the JSON is written without it). The tool exits
non-zero if an invariant claim breaks; pre-registered *hypotheses* are reported with a
verdict and may be falsified without failing the run, because a falsified hypothesis is
a result, not a regression.

## 3. Results

### 3.1 Axial ladder (`nz = 81` fixed, reference grid 321x81)

| nx | h_x [m] | volume | waterplane area | wetted area |
|---|---|---|---|---|
| 11 | 1.0000 | -1.8895e-02 | -1.8477e-02 | -2.8484e-03 |
| 21 | 0.5000 | -6.1076e-03 | -5.6840e-03 | -1.1981e-03 |
| 41 | 0.2500 | -2.1745e-03 | -1.7492e-03 | -6.9494e-04 |
| 81 | 0.1250 | -9.6417e-04 | -5.3837e-04 | -4.8736e-04 |
| 161 | 0.0625 | -5.9166e-04 | -1.6570e-04 | -3.7814e-04 |
| 321 | 0.0312 | -4.7701e-04 | -5.1000e-05 | -3.1508e-04 |

### 3.2 Vertical ladder (`nx = 321` fixed)

| nz | h_z [m] | volume | VCB | max section area |
|---|---|---|---|---|
| 5 | 0.2000 | -3.5125e-02 | -4.1318e-02 | -3.5075e-02 |
| 9 | 0.1000 | -1.2855e-02 | -1.6449e-02 | -1.2805e-02 |
| 13 | 0.0667 | -7.1181e-03 | -9.3910e-03 | -7.0675e-03 |
| 21 | 0.0400 | -3.3809e-03 | -4.5663e-03 | -3.3301e-03 |
| 41 | 0.0200 | -1.2445e-03 | -1.6852e-03 | -1.1936e-03 |
| 81 | 0.0100 | -4.7701e-04 | -6.1322e-04 | -4.2603e-04 |

### 3.3 Per-axis orders vs the endpoint-exponent prediction

| quantity | p_x (measured) | p_z (measured) |
|---|---|---|
| displacement_volume | 1.819 | 1.646 |
| waterplane_area | 1.819 | — (x-independent) |
| vertical_centre_buoyancy | — (z-independent) | 1.572 |
| maximum_section_area | — (z-independent) | 1.646 |
| wetted_surface_area | 1.301 | 1.432 |

All QC invariants: axial_order_matches_endpoint_exponent=PASS; vertical_order_matches_keel_exponent=PASS; axis_independent_quantities=PASS; separable_error_model=PASS; vertical_axis_dominates_default_grid=PASS; equal_cell_split_beats_default=PASS; default_grid_within_three_percent=PASS.

An order of `—` is not a failure of the fit: it is the *exact* statement that the
quantity does not integrate that axis at all. The waterplane area is an x-integral
along the waterline, so it is bit-identical for every `nz` (spread
`0.0e+00`); VCB and the
maximum section area are z-integrals of one section, so they are invariant in `nx`
(`3.5e-16` and
`0.0e+00`). Those
three exact zeros are a stronger check of the integrator than any fitted slope.

### 3.4 Where the default grid's error lives

Default grid 41x13, volume error `-8.8043e-03`:

| component | value | share of total |
|---|---|---|
| axial (41 vs 321) | -1.6975e-03 | 19.3 % |
| vertical (13 vs 81) | -6.6411e-03 | 75.4 % |
| reference floor (321x81) | -4.7701e-04 | 5.4 % |

Three quarters of the default-grid volume error is **vertical**. Refining the same
`nx = 41` mesh from `nz = 13` to `nz = 41` cuts the volume error from
`-8.8043e-03` to
`-2.9407e-03` — a factor
`2.99` —
for 3.3x the cells.

### 3.5 The two axis errors add

| grid | cells | measured | predicted | residual | residual/measured |
|---|---|---|---|---|---|
| 21x13 | 240 | -1.2711e-02 | -1.2749e-02 | +3.7411e-05 | 0.294 % |
| 41x13 | 480 | -8.8043e-03 | -8.8156e-03 | +1.1279e-05 | 0.128 % |
| 41x21 | 800 | -5.0735e-03 | -5.0784e-03 | +4.9318e-06 | 0.097 % |
| 81x21 | 1600 | -3.8667e-03 | -3.8681e-03 | +1.4154e-06 | 0.037 % |
| 161x41 | 6400 | -1.3591e-03 | -1.3591e-03 | +8.8035e-08 | 0.006 % |
| 321x41 | 12800 | -1.2445e-03 | -1.2445e-03 | +0.0000e+00 | 0.000 % |
| 161x81 | 12800 | -5.9166e-04 | -5.9166e-04 | +0.0000e+00 | 0.000 % |

Worst residual `0.29 %` of the measured
error, and it shrinks with the grid (it is the `O(h_x^2 h_z^2)` cross term). The reason is
visible in the geometry: the half-breadth is a *product* `b(x) * sqrt(z/T)`, and the
bilinear interpolant of a product is the product of the two linear interpolants, so the
2-D trapezoid error has no leading cross term at all. This is measured here, not assumed.

### 3.6 Equal cell budget: the split matters less than expected

| grid | cells | volume | VCB | waterplane area |
|---|---|---|---|---|
| 41x13 | 480 | -8.8043e-03 | -9.3910e-03 | -1.7492e-03 |
| 31x17 | 480 | -7.4677e-03 | -6.2677e-03 | -2.8527e-03 |
| 25x21 | 480 | -7.4851e-03 | -4.5663e-03 | -4.1689e-03 |
| 21x25 | 480 | -8.2132e-03 | -3.5193e-03 | -5.6840e-03 |
| 17x31 | 480 | -1.0120e-02 | -2.5541e-03 | -8.3070e-03 |
| 13x41 | 480 | -1.4727e-02 | -1.6852e-03 | -1.3550e-02 |
| 11x49 | 480 | -1.9371e-02 | -1.2932e-03 | -1.8477e-02 |

## 4. Notes, limits, falsified hypothesis

* **Falsified hypothesis (kept on the record).** Before looking at the budget grid I
  predicted *"at a fixed cell budget the best split is nz-heavy (nz > nx) and beats
  41x13 by >= 3x"*. It is **false**: the budget optimum is the *balanced* split
  `31x17`, worth only
  `1.18x`,
  and nz-heavy splits are the *worst* of the seven (11x49 is the worst). The reason both
  statements are true at once is that the default grid is not merely z-starved -- it is
  at the x-heavy corner -- so adding cells to either extreme loses to balancing them.
  The correct reading of §3.4 is therefore: *at fixed nx*, refine nz; *at fixed cells*,
  move towards balance.
* **Wetted area is the slowest axis, not the worst error.** Its per-axis orders
  (`1.30` axially,
  `1.43` vertically) are the
  lowest of the six quantities, because the integrand carries the chord's slope and the
  slope of a `sin**0.7` taper is unbounded at the ends. Its default-grid error is
  nevertheless only `0.33 %`,
  smaller than the volume's `0.88 %`.
  The reason is *not* cancellation: the same additive model reproduces the wetted-area
  error to `0.46 %` (computed in the tool under
  `hypotheses`) -- the axis coefficients are simply smaller. Slow order, small constant.
* **One surface only.** These orders belong to this analytic loft with its two chosen
  exponents. A different family (Series 60, planing, Wigley) will have its own; the tool
  is parameter-free, so re-running it on another `HullParameters` is the way to get
  them.
* **Symmetric fore/aft only.** `fore_exponent == aft_exponent` here, so the longitudinal
  singularity is the same at both ends. An asymmetric hull has two different exponents
  and therefore two different endpoint orders; not measured here.

## 5. Error verdict

Against the closed forms at the default 41x13 grid:

| quantity | below 3 % ? | absolute error |
|---|---|---|
| displacement_volume | < 3 % | 0.880 % |
| waterplane_area | < 3 % | 0.175 % |
| vertical_centre_buoyancy | < 3 % | 0.939 % |
| block_coefficient | < 3 % | 0.880 % |
| wetted_surface_area | < 3 % | 0.326 % |

**Every quantity the audit named is < 3 % at the default grid.** The largest is the VCB
at `0.94 %`, then the
displacement-derived block coefficient at
`0.88 %`; the wetted area is
the smallest. The per-axis
predictions are met to within
`0.12`
orders (axial `1.82` vs 1.7;
vertical `1.65` vs 1.5).

**Caveat, stated plainly:** the fitted orders are *asymptotic estimates over a finite
ladder*, not laws. They scatter by ±0.2 between neighbouring refinements (the reduced
ladder in the test file reproduces the same envelope). What is not fuzzy is the
attribution in §3.4 (three quarters vertical) and the exact zeros in §3.3.

## 6. Reproduction

```bash
python3 tools/per_axis_convergence.py --output artifacts/per_axis_convergence --doc docs/per_axis_convergence.md
python3 -m pytest -q tests/test_per_axis_convergence.py
```

Runtime ~45 s for the tool (26 hydrostatics evaluations, largest 321x81), ~7 s for
the test file (reduced ladder). The tool exits non-zero if any QC invariant breaks.
