# The Wigley hull form

`tensorcad.hull` lofts exactly one family: `b(x, z) = (B/2) sin(pi x/L)^0.7 sqrt(z/T)`.
It is fore-aft symmetric, has no transom, and -- the reason this page exists -- it is not
the surface anybody publishes numbers against. `tensorcad.wigley` adds the form that is:

```
y(x, w) = (B/2) (1 - |2x/L - 1|^n) (1 - (w/T)^m)          w = depth below the design WL
```

with `w = 0` at the waterline and `w = T` at the keel. `n = m = 2` is the classic Wigley
used as the geometry source in wave-resistance, seakeeping and diffraction comparisons;
the exponents are exposed because the same algebra holds for any `n, m >= 1`.

The point of the family is not that it is a ship. It is that its hydrostatics are **exact
rationals**, so a mesh integrator can be tested against algebra instead of against another
mesh. That makes this module simultaneously a hull-form generator and a verification case
for `tensorcad.hydrostatics`.

## Closed forms

With `taper(x) = 1 - |2x/L - 1|^n` and `section(w) = 1 - (w/T)^m`, the two factors of the
surface separate, so every moment is a product of two one-dimensional integrals:

| integral | value | at `n = m = 2` |
|---|---|---|
| `int_0^L taper dx` | `L n/(n+1)` | `2L/3` |
| `int_0^T section dw` | `T m/(m+1)` | `2T/3` |
| `int_0^T w section dw` | `T^2 m/(2(m+2))` | `T^2/4` |

which gives, for `V = 2 int_0^L int_0^T (B/2) taper(x) section(w) dw dx`:

| quantity | closed form | classic Wigley |
|---|---|---|
| displacement volume | `n m/((n+1)(m+1)) L B T` | `4LBT/9` = `0.4444 LBT` |
| waterplane area | `L B n/(n+1)` | `2LB/3` |
| midship section area | `B T m/(m+1)` | `2BT/3` |
| block coefficient `Cb` | `n m/((n+1)(m+1))` | `4/9` |
| prismatic `Cp` | `n/(n+1)` | `2/3` |
| midship `Cm` | `m/(m+1)` | `2/3` |
| waterplane `Cwp` | `n/(n+1)` | `2/3` |
| VCB above keel | `T - T(m+1)/(2(m+2))` | `5T/8` |
| LCB, LCF | `L/2` (fore-aft symmetry, identically) | `L/2` |
| wetted surface | no elementary primitive; composite Gauss-Legendre | `14.879063105 m^2` at `L/B = 10`, `B/T = 1.6` |

`Cb = Cp Cm` holds by construction, which is the cheapest self-check the module has
(`tests/test_wigley.py` asserts it for four exponent pairs).

## Measured convergence

`examples/wigley/run_example.py` regenerates `docs/wigley/{results.json,report.md,plot.png}`
and exits non-zero if any claim stops holding. Reference case `L = 10 m`, `B = 1 m`,
`T = 0.625 m` (`L/B = 10`, `B/T = 1.6` -- the ITTC proportions):

| grid | triangles | V [m^3] | e_V [%] | e_Awp [%] | e_S [%] | max abs err [%] |
|---|---|---|---|---|---|---|
| 21 x 9 | 678 | 2.7600098 | -0.63965 | -0.25000 | -0.04570 | 0.63965 |
| **41 x 13** (default) | 1 998 | **2.7712221** | **-0.23600** | -0.06250 | -0.02009 | 0.23600 |
| 81 x 21 | 6 558 | 2.7756079 | -0.07812 | -0.01562 | -0.00720 | 0.07812 |
| 161 x 41 | 25 918 | 2.7772353 | -0.01953 | -0.00391 | -0.00180 | 0.01953 |

Observed order fitted across the sweep: volume 2.17, waterplane 2.58, wetted area 2.01 --
second order, as a chord-based loft should be. Every signed error is **negative**: the
loft chords cut inside the analytic surface, so the discretisation is systematic. A sign
flip would mean the generator stopped sampling the surface it documents, and the example
fails on it.

For comparison, the `sin(pi x/L)^0.7 sqrt(z/T)` family on the same 41 x 13 grid is
0.88 % low on volume. The difference is physics of the form, not sampling: `sqrt` gives an
infinite section derivative at the keel, the parabolic Wigley a finite keel angle, so the
trapezoidal height integration converges a whole order faster. This is asserted as a test
(`test_the_smooth_form_beats_the_sine_family_on_the_same_grid`).

## Usage

```python
from tensorcad import generate_wigley, wigley_reference_parameters, export_wigley

model = generate_wigley(wigley_reference_parameters(nx=81, nz=21))
model.hydrostatics.block_coefficient          # 0.4441 (exact 4/9)
model.validation.passed                       # topology + dimension + ratio checks
export_wigley(model, "runs/wigley", stl=True) # wigley.obj / _offsets.csv / .json / .stl
```

Shaping the family (`n` fattens or sharpens the shoulders, `m` the sections):

```python
from tensorcad import WigleyParameters, generate_wigley
from tensorcad.wigley import wigley_analytic_hydrostatics

model = generate_wigley(WigleyParameters(length=12, beam=1.5, draft=0.8,
                                        longitudinal_exponent=3.0, vertical_exponent=2.0))
wigley_analytic_hydrostatics(model.parameters).block_coefficient   # 3*2/((4)(3)) = 0.5
```

The waterline offsets a panel or RANS code wants are `model.offsets` as `(x, breadth)`
pairs; `wigley.obj` plus `wigley.json` carry the surface and its provenance. The metadata
document embeds both the integrated hydrostatics and the closed-form reference under
`analytic`, so a downstream run can be checked without this package.

## Limits

* No appendages, no transom, no parallel middle body, no propeller aperture: this is a
  benchmark hull. A real design still needs `tensorcad.suboff_appendages` or a fitted
  offsets table.
* Fore-aft symmetric, so `LCB = LCF = L/2` **identically** -- this case does not exercise
  the longitudinal moment integrator. `tensorcad.hull` and the SUBOFF case do.
* Wetted area is a quadrature reference, not a closed form. It is exact to round-off for
  `n = m = 2`; for non-integer exponents the `|u|^(n-1)` term kinks at midship, so use an
  even panel count (the default is 4 x 24) and re-run `wetted_area_evidence` before quoting
  a number. Measured: 4e-16 at `n = m = 2`, 8e-10 at `n = 1.5, m = 3`.
* Hydrostatics only. Stability, resistance, structural and rule checks are not implemented
  here or anywhere else in the package.
* `nx` must be odd for the realised beam to reach the requested beam, exactly as in
  `tensorcad.hull`; an even `nx` is accepted with a `requested_beam_not_achievable` warning
  rather than being silently short.