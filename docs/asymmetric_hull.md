# The asymmetric hull family: LCB that a formula has to predict

`tensorcad.hull` samples one symmetric surface. Its LCB and LCF are `L/2` *by
construction*, so the moment test in that module can catch a sign or a frame error
but cannot catch an antisymmetric bias, because the geometry has none to catch.
This page is the case that can: the two halves of the length get their own taper
exponent, so where the volume sits becomes a prediction instead of an identity.

Implementation: `tensorcad/asymmetric_hull.py` (0.12.0). Regression: 21 of the 46
tests in `tests/test_asymmetric_hull.py` are pinned to the numbers below.

## The family

```
b(x, z) = (B/2) · sin(pi·x/L)^p(x) · sqrt((z + T)/T)
p(x)   = fore_exponent   for x ≥ L/2
         aft_exponent    for x < L/2
```

`sin(pi·x)` is 1 at midship, so the surface is continuous across the switch
whatever the exponents are; only the waterline slope breaks there, which is what a
bow/stern asymmetry looks like on a real hull. The vertical profile stays the
`sqrt` of the symmetric family, so VCB and the section-based coefficients keep
their closed forms and the *only* new physics is longitudinal.

Exponents are restricted to `[0.4, 3.0]`: the wetted-surface rule and the moment
substitution are calibrated inside that range, and the constructor refuses values
outside it rather than returning a number of unknown accuracy.

## The reference, and what is analytic about it

The surface stays separable, so every volume integral factors:

| quantity | reference | kind |
|---|---|---|
| displaced volume | `L·B·T·(2/3)·(I_f + I_a)` | closed form |
| waterplane area | `L·B·(I_f + I_a)` | closed form |
| midship section area | `2/3·B·T` | closed form |
| VCB | `−2/5·T` | closed form |
| Cb, Cp | `(2/3)(I_f + I_a)`, `I_f + I_a` | closed form |
| LCB, LCF | `L·(M_a + I_f − M_f)/(I_a + I_f)` | quadrature (`M_p`) |
| wetted surface area | per half, summed | quadrature |

with `I_p = ∫_0^{1/2} sin(pi·x)^p dx` and `M_p = ∫_0^{1/2} x·sin(pi·x)^p dx` in
normalised `x`. `I_p` has a Beta-function closed form
`I_p = J(p)/pi`, `J(p) = sqrt(pi)·Γ((p+1)/2)/(2Γ((p+2)/2))`, evaluated with
`math.lgamma` so the core stays standard-library only. The fore half contributes
`∫_{1/2}^{1} x·f dx = I_f − M_f` after mirroring onto `[0, 1/2]`; both halves share
the normalisation `I_a + I_f`.

`LCF == LCB` for this family because the waterplane and the volume are weighted by
the **same** `p(x)`. That is a property of the family, not of the integrator, and
it is asserted as such -- it is not presented as independent agreement.

`M_p` has no elementary antiderivative. It is integrated with Gauss-Legendre after
`x = (1/2)·tau^m`, `m = max(1, 2/p)`: the pointed members of the range (`p < 2`)
get their `x^p` tip turned into `tau^2`, the already-regular members (`p ≥ 2`) stay
on plain Gauss-Legendre. The `2/p` rule is a measured optimum, not a guess -- at
`p = 3` the `2/p` power *loses* three orders (9.7e-12 versus 1.1e-16
node-refinement drift).

## Measured: mesh against the closed form

Case `L = 10 m`, `B = 2 m`, `T = 0.8 m`, `p_fore = 0.7`, `p_aft = 1.6`, seawater
1025 kg/m³ -- a fine bow over a fuller run, which is the shape the symmetric family
cannot express. Signed relative error of the mesh integrals
(`compute_hydrostatics`, divergence theorem on the triangles) against the reference:

| quantity | 11×5 | 21×9 | 41×13 | 81×21 | 161×41 |
|---|---|---|---|---|---|
| displacement volume (= Cb) | −4.559 % | −1.605 % | **−0.806 %** | −0.364 % | **−0.129 %** |
| wetted surface area | −1.111 % | −0.506 % | −0.308 % | −0.167 % | −0.072 % |
| waterplane area (= Cp) | −1.089 % | −0.328 % | −0.100 % | −0.031 % | −0.009 % |
| LCB (= LCF) | −0.824 % | −0.262 % | −0.082 % | −0.026 % | −0.008 % |
| VCB | −4.132 % | −1.645 % | −0.939 % | −0.457 % | −0.169 % |
| midship section area | −3.508 % | −1.281 % | −0.707 % | −0.333 % | −0.119 % |

Reference values at that case: `V = 6.655797633 m³`, `Awp = 12.479620561 m²`,
`S = 22.006225480 m²`, `LCB = LCF = 5.413279060 m`, `VCB = −0.320 m`,
`Cb = 0.415987352`, `Cp = 0.623981028`, `Am = 1.066666667 m²`.

Every quantity is inside 3 % at the 41×13 grid and inside 0.2 % at 161×41, the
error shrinks monotonically with station spacing, and the LCB error keeps its sign
across the sweep -- the sign is a prediction of the tip under-integration, not
noise. VCB is again the largest default-grid error, for the reason recorded in
`docs/hull_analytic_reference.md`: the trapezoidal section rule is high near the
keel, where the `sqrt` profile is steepest.

`LCB` sits `+0.413279 m` (4.13 % of `L`) forward of midship because the fore half
is the fuller one; swapping the exponents gives `4.586720940 m`, and the two sum to
`10.000000000000 m`. The symmetric family would report `5.000000 m` for both and
consider itself verified.

## How the reference itself is checked

Nothing here is checked only against itself:

* `I_p` at `p = 0.7` reproduces the constant shipped in 0.11.0
  (`hull.TAPER_INTEGRAL/2`) to 2.2e-16 relative, and matches a plain Simpson rule
  within 1e-7 across the calibrated range.
* `M_p` is integrated with two different substitution powers (agreement < 1e-7) and
  to quadruple the node count (agreement < 1e-12).
* The whole reference against the mesh path: two independent integrators, tabulated
  above.
* The symmetric limit against the shipped 0.11.0 reference: volume, waterplane,
  VCB, Cb, Cp, Am to 5.6e-16 relative; LCB to 2.2e-16; wetted area to 6.5e-11 (the
  two rules use different substitution exponents). The generated mesh at equal
  exponents is vertex-identical to `hull.build_mesh` (worst coordinate difference
  1.1e-16).
* The wetted-surface quadrature: doubling the node count changes it by 2.3e-12
  relative at the default case, and by < 1e-9 at both exponent extremes.

## Limitations

* **Still one family.** No transom, bulb, flare, knuckle or appendages; `p(x)` is
  the only asymmetry knob, and the switch has to be at midship.
* **LCB is the only new signal.** Cb/Cp/Awp/VCB/A_m are integrals of the same
  separable structure as before; their error behaviour is the symmetric family's.
* **LCF is not independent.** As above: it coincides with LCB by construction.
* **Wetted area is quadrature, not a formula.** Its accuracy is established by two
  quadrature families and node refinement, not by an external measurement.
* **Waterline at `z = 0`.** Trim and heel are out of scope, as everywhere in this
  package.
* **No hydrodynamic claim.** Geometry and integration only: nothing here predicts
  resistance, stability, seakeeping or class compliance, and the exponents are not
  fitted to any measured hull.

## Reproduction

```bash
python examples/asymmetric_hull/run_asymmetric_check.py     # ~15 s, exits 1 if a claim breaks
python -m pytest -q tests/test_asymmetric_hull.py           # 46 tests
```

`examples/asymmetric_hull/run_asymmetric_check.py` writes
`docs/asymmetric_hull/{results.json,report.md,convergence.png}`; the tables above
are transcribed from that `results.json`, and the script fails rather than
reporting a claim it can no longer support.