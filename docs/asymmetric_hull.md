# Asymmetric fore/aft taper — the longitudinal centres become measurable

`HullParameters` gained two fields in 0.13.0:

| field | default | meaning |
|---|---|---|
| `fore_exponent` | 0.7 | sine-taper exponent on the half forward of midship |
| `aft_exponent` | 0.7 | sine-taper exponent on the half abaft midship |

The half-breadth is still separable, `(B/2) * sin(pi*x/L)^p * sqrt((z+draft)/draft)`,
with `p` chosen per half. The defaults reproduce the shipped symmetric surface bit for
bit — same vertices, same faces, same hydrostatics — so nothing stored or pinned moved.

## Why the parameter exists

Up to 0.12.0 the family was mirror-symmetric, and `LCB = LCF = length/2` followed from
that symmetry. The check was worth having (it validates the frame and the sign of the `x`
moments) but it could not detect an antisymmetric bias, because the geometry had none.
Giving the halves separate exponents makes the centre a **quantity** again:

* it has an independent reference (below),
* it satisfies an exact involution — swapping the exponents reflects the centre about
  midship: `LCB(p_aft, p_fore) = L - LCB(p_fore, p_aft)`,
* it responds monotonically to the asymmetry, so the sign of the response is testable.

The maximum section stays at midship for every exponent (`sin(pi*x)^p` peaks at `x = L/2`),
so the midship section area and therefore `Cp` are untouched by the split.

## Reference

Write `h(p) = int_0^(1/2) sin(pi x)^p dx = taper(p)/2` and
`m(p) = int_0^(1/2) x sin(pi x)^p dx`. Then

| quantity | reference |
|---|---|
| displacement volume | `L*B*T*(h(p_aft) + h(p_fore))*2/3` |
| waterplane area | `L*B*(h(p_aft) + h(p_fore))` |
| maximum section area | `2/3*B*T` (independent of the exponents) |
| VCB | `-2/5*T` (the `sqrt` profile is unchanged) |
| `Cb`, `Cp` | `volume/(L*B*T)`, `volume/(L*Am)` |
| LCB, LCF | `L*(m(p_aft) + h(p_fore) - m(p_fore)) / (h(p_aft) + h(p_fore))` |
| wetted area | Gauss-Legendre quadrature, one element per half |

`h(p) = Gamma((p+1)/2) / (sqrt(pi) Gamma(p/2+1))` is closed form. `m(p)` is not — the
`x`-weighted integral of `sin^p` is not a Beta function, only the unweighted one is — so
it is computed by Gauss-Legendre after `x = t^3/2` (which regularises the tip: the
integrand behaves like `t^(3p+5)`) and independently by an exponential substitution
`x = exp(-v)/2` with a plain Gauss rule on `v in [0, 50]`. The two share neither the
substitution nor the rule and agree to `5e-15` relative over `p in [0.2, 3]`. `basis` in
the returned `Hydrostatics` says which numbers are formula and which are quadrature, so a
caller comparing mesh numbers knows what the reference rests on.

## Measured (`examples/asymmetric_hull/run_example.py`, `docs/asymmetric_hull/`)

13 claims, all passing; the numbers are in `results.json` and `report.md`.

* **Symmetric form untouched**: analytic volume equals `B*L*T*TAPER_INTEGRAL*2/3` exactly,
  the centre is exactly `L/2`, and `taper_integral(0.7)` *is* `TAPER_INTEGRAL`.
* **Centre converges onto the reference**: for `(p_f, p_a) = (0.5, 0.9)` the mesh LCB error
  runs `-7.79 mm` at 41x13, `-2.92` at 81x25, `-1.08` at 161x49, `-0.39` at 321x97 — one
  sign, shrinking, so it is a discretisation bias and not noise. At 161x49 the error is
  `0.011 %` of the length.
* **Involution holds on the mesh to 1e-6 m**: swapping the exponents mirrors the centre,
  while the mesh volume and waterplane area are unchanged to `1e-15` relative. What breaks
  the identity is the discretisation error, not the geometry.
* **`LCB = LCF` becomes two independent mesh integrals.** On the symmetric family the
  equality is one identity; here the buoyancy moment and the waterplane moment are
  different sums that have to land on the same number, and they do to `1e-9`.
* **Monotone response**: a finer forebody (`p_fore` up) moves the centre aft, over
  `p_fore in [0.3, 2.0]` at `p_aft = 0.9`; the reference is continuous through equal
  exponents.

## Limitations

* The family is still synthetic: two exponents on one sine taper, no transom, no bulb, no
  appendages, no parallel middle body. An asymmetric *real* hull is not represented; what
  is represented is a surface whose longitudinal centres are known exactly.
* `m(p)` is quadrature, not a formula. Its accuracy is established by two independent
  schemes and by node-count stability, not by an external source.
* Nothing here changes resistance, stability, seakeeping or class compliance: the model
  still reports geometry and hydrostatic integrals only.
* Fairness is unchanged and still poor at the ends: the taper's slope is singular at the
  stem and the post for every `p < 1` (see `docs/fairness.md`), and an asymmetric split
  moves that singularity but does not remove it.
