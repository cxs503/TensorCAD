# Two hull families, two Bonjean curves, two closed forms

`docs/bonjean.md` measured the station cut on a barge (exact) and on the shipped
`hull` family against `analytic_section_area`. Two families arrived afterwards
without a reference to measure them against:

* the **asymmetric** hull (0.17.0) — `fore_exponent` / `aft_exponent` split the
  sine taper per half; its closed forms covered volume, waterplane, LCB/LCF and
  VCB, but never the station-resolved area that damage stability reads;
* the **Wigley** family (0.20.0) — a closed form for hydrostatics, none for the
  area under a given waterline.

This page closes both, and every number below is written by
`python examples/bonjean_families/run_example.py` into
`docs/bonjean_families/results.json` (with `plot.png` and `convergence.png`);
nothing here is typed from a different run. The script exits non-zero if one of
its 8 claims stops holding, and `tests/test_bonjean_families.py`
pins 79 cases of the same content.

## The two references

The hull family's half breadth is `(B/2) sin(pi x/L)^p sqrt(zeta/T)` for
`zeta = z + T` above the keel, so

    A(x, z) = (2/3) B sin(pi x/L)^p zeta^(3/2) / sqrt(T),   p = fore on x >= L/2, aft below

`p` is the only thing the 0.17.0 split changed, and the vertical taper (`sqrt`)
is shared by every member of the family. `bonjean.taper_exponent_at` holds the
split rule in one place: midship takes the fore exponent, which is what
`hull._breadth_exponent` does (`2 i >= nx-1` is exactly `x >= L/2`), and
`sin(pi/2) = 1` makes the choice immaterial there.

The Wigley surface is `y = (B/2)(1-|2x/L-1|^n)(1-(w/T)^m)` with `w` the depth
below the design waterline. Integrating over the immersed height `H = z + T`,
clamped to `[0, T]`:

    A(x, H) = B (1-|2x/L-1|^n) [H - T/(m+1)(1 - (1 - H/T)^(m+1))]

which collapses to the published `analytic_maximum_section_area` at midship and
`H = T`, and to zero below the keel. Both references are checked against a
Gauss-Legendre quadrature of the *surface itself* (`wigley.half_breadth`, and the
generator's own half-breadth law) — an integration path independent of the
algebra — with `zeta = H s^2` substituted so the `sqrt` taper is analytic and the
comparison is exact to rounding rather than converging like `h^1.5`.

## Result 1 — Wigley, meshed against its closed form

`nx = 81, nz = 13`, ITTC reference proportions (`L = 10 m`, `B = 1 m`,
`T = 0.625 m`, `n = m = 2`), 5 stations x 4 waterlines. First twelve samples:

| station [m] | waterline [m] | mesh [m²] | closed form [m²] | rel. error |
|---:|---:|---:|---:|---:|
| 1.00 | -0.5000 | 0.008349 | 0.008400 | -0.608 |
| 1.00 | -0.3125 | 0.046745 | 0.046875 | -0.278 |
| 1.00 | -0.1562 | 0.094727 | 0.094922 | -0.206 |
| 1.00 | -0.0500 | 0.131800 | 0.132038 | -0.181 |
| 2.50 | -0.5000 | 0.017394 | 0.017500 | -0.608 |
| 2.50 | -0.3125 | 0.097385 | 0.097656 | -0.278 |
| 2.50 | -0.1562 | 0.197347 | 0.197754 | -0.206 |
| 2.50 | -0.0500 | 0.274582 | 0.275080 | -0.181 |
| 5.00 | -0.5000 | 0.023192 | 0.023333 | -0.608 |
| 5.00 | -0.3125 | 0.129847 | 0.130208 | -0.278 |
| 5.00 | -0.1562 | 0.263129 | 0.263672 | -0.206 |
| 5.00 | -0.0500 | 0.366110 | 0.366773 | -0.181 |

Worst over the whole grid: **-0.608 %** at
`x = 9 m`, `z = -0.5 m`. The midship
section at the design waterline is the quantity the hydrostatics module already
compared against, so it pins the new expression to the old one:

| nz | mesh [m²] | closed form [m²] | rel. error |
|---:|---:|---:|---:|
| 13 | 0.415943287 | 0.416666667 | -0.174 |
| 25 | 0.416485822 | 0.416666667 | -0.043 |
| 49 | 0.416621455 | 0.416666667 | -0.011 |

## Result 2 — the asymmetric hull, meshed against the split reference

`fore = 1.4`, `aft = 0.7`,
`81 x 41`, 10 stations x 7 waterlines. Worst error per waterline:

| waterline z [m] | worst rel. error |
|---:|---:|
| -0.60 | +0.924 |
| -0.50 | +0.509 |
| -0.40 | +0.333 |
| -0.30 | +0.239 |
| -0.20 | +0.183 |
| -0.10 | +0.145 |
| -0.02 | +0.124 |

The trend is the family's `sqrt` keel taper: the slope is infinite at the keel,
a polygon in `z` cannot follow it, and the error grows as the waterline
approaches it. Worst above `z = -0.4 m`: **+0.333 %**;
at the coarsest waterline count (`81 x 13`) the deepest curve reaches
**+5.308 %**. That is reported
per waterline rather than averaged away, because the average would hide exactly
the region a damage-stability calculation cares about.

## Result 3 — the split is discriminating, not decorative

Agreement with a reference proves nothing unless a *wrong* reference fails.
Reading the same asymmetric mesh against a symmetric reference — the shipped
exponent on both halves — at `x = 0.9 L`:

| reference | worst rel. error at `z = -0.2 m` |
|---|---:|
| symmetric (shipped exponent everywhere) | +72.758 % |
| split (`fore`, `aft` as generated) | +0.183 % |

A factor of ~400 between the two, so the sub-percent agreement above is evidence
that the fore/aft exponents reach the mesh and the closed form the same way.
Two structural checks come with it: swapping the exponents mirrors the measured
curves to `5.54e-08` relative, and
integrating the Bonjean table over 41 stations recovers the asymmetric
closed-form volume to `-0.142 %`
(`6.787274 m³` against `6.796919 m³`).

## Result 4 — both families are second order in the waterline spacing

The section polygon is linear in `z` between vertex rings, so the error should
fall as the square of the spacing; if it does not, the reference or the cut is
wrong.

| Wigley: nz | worst rel. error |
|---:|---:|
| 13 | +0.608 |
| 25 | +0.158 |
| 49 | +0.039 |

| asymmetric hull: nz | worst rel. error at `z = -0.2 m` |
|---:|---:|
| 13 | +1.078 |
| 25 | +0.389 |
| 49 | +0.140 |
| 97 | +0.050 |

Observed ratios: `['3.84', '4.06']` (Wigley) and
`['2.77', '2.79', '2.80']` (hull). Both are ~4 and ~2.8:
second order, with the Wigley constant smaller because its section is smooth at
the keel for `m = 2` while the hull family's `sqrt` is not.

## Claims

| claim | holds | detail |
|---|---|---|
| wigley sections within 1 % of the closed form | yes | 0.608 % at x=9 m, z=-0.5 m |
| wigley midship equals the known section area | yes | -0.174 % at nz = 13 |
| asymmetric sections within 2 % above z = -0.4 m | yes | 0.333 % at nz = 41 |
| the split reference is discriminating | yes | symmetric 72.76 % vs split 0.18 % |
| swapping the exponents mirrors the curves | yes | max |rel| 5.54e-08 over 10 stations |
| wigley error falls like the square of the waterline spacing | yes | ratios ['3.84', '4.06'] |
| asymmetric error falls with refinement | yes | ratios ['2.77', '2.79', '2.80'] |
| integrating the curves recovers the closed-form volume | yes | 0.142 % over 41 stations |

## What this does not cover

* **Above the design waterline** the mesh stops (the loft is capped) while the
  hull family's closed form keeps growing as `zeta^(3/2)`; only the capped range
  is compared, and the Wigley reference clamps to `H <= T` by definition.
* **Keel-region error is characterised, not removed.** `nz = 13` on the `sqrt`
  family is too coarse near the keel for a 3 % claim; the honest statement is
  the per-waterline table above plus the second-order ladder.
* **Plumb, upright, at rest.** No trim, heel, sinkage, flooding or structural
  load is modelled here — this page verifies the geometry of the cut, and a
  Bonjean curve is only the input to those questions.
