# The planing hull family: constant deadrise, chine planform power law

Savitsky's planing method (Savitsky 1964, *Hydrodynamic design of planing hulls*,
Marine Technology 1(1)) is written for one shape family: a bottom of **constant
deadrise** `beta`, a **chine** whose breadth grows along the length, and a **transom**.
Everything that method asks for -- deadrise, chine beam, wetted length-to-beam ratio --
is a property of that surface, and the surface itself is what a study of more than one
parent form has to loft. `tensorcad.series60` covers the displacement side of the yard
and `tensorcad.wigley` the benchmark side; a planing bottom was the remaining gap in the
hull-family list. This page documents `tensorcad.planing`, the closed-form hydrostatics
of the family, and the evidence that the loft and the closed forms agree.

Everything on this page is measured by `python examples/planing/run_example.py`, which
writes `artifacts/planing_hull/results.json` and exits non-zero if any claim stops
holding. `tests/test_planing.py` pins the same numbers (41 cases).

## 1. The geometry

Bow at `x = 0`, transom at `x = L`, waterline `z = 0`, keel `z = -T`, centreplane
`y = 0`. The chine half-breadth is a power law,

    b(x) = (B/2) (x/L)**q,        q > 0

`q = 1` is a straight chine (a true prism with constant deadrise); `q < 1` hollows the
bow and fills the stern, which is how a planing bottom is usually drawn. The section is
a V of deadrise `beta` from the keel up to the chine at
`z_chine(x) = -T + b(x) tan(beta)`, with a vertical topside above it, so the submerged
half-breadth is

    y(x, z) = min( (z + T)/tan(beta),  b(x) )

Two regimes therefore live inside one hull. Near the bow the chine is *below* the
waterline, the topside is wetted and the section is V plus a wall. Aft of

    x0 = L (T / ((B/2) tan(beta)))**(1/q)          (clipped to [0, L])

the chine is *above* the waterline and the waterline cuts the V, so the section is the
plain triangle `T**2 / tan(beta)`. A knife bow (`b(0) = 0`) makes the first regime
unavoidable at the stem and the second one likely at a wide transom, so both branches
exist in every hull of the family and both are integrated.

| Case | `L` [m] | `B` [m] | `T` [m] | `beta` [deg] | `q` | `x0` [m] | waterline beam [m] | triangles | chine immersed at transom |
|---|---|---|---|---|---|---|---|---|---|
| `chine_immersed` | 10 | 3.0 | 0.55 | 18 | 0.6 | 10.000 | 3.0000 | 5312 | yes |
| `chine_dry` | 10 | 3.0 | 0.40 | 18 | 0.6 | 7.194 | 2.4621 | 5312 | no |
| `straight_chine` | 10 | 3.0 | 0.60 | 18 | 1.0 | 10.000 | 3.0000 | 5312 | yes |
| `deep_vee` | 10 | 3.0 | 0.60 | 24 | 0.5 | 8.071 | 2.6952 | 5312 | no |

The mesh is the **wetted** hull: the V-bottom, the wetted strip of topside, and the
transom, closed by the planar waterline cap -- the convention
`tensorcad.hydrostatics` states. That is why the reported waterline beam of a
chine-dry hull is `2T/tan(beta) = 2.462 m` and not `B = 3 m`: the dry part of the chine
is deliberately not in the mesh. Two construction details matter for exactness:

* the z-grid of each station is **snapped onto `z_chine`**, so the chine kink is not
  rounded by the sampling;
* the waterline `z = 0` is a vertex line, so no face straddles it and the waterplane
  integral of the divergence theorem is exact.

The transom is triangulated as a fan from an interior point of the last section (the
wetted section is convex, so the fan has no degenerate facet), which is what makes a
planing hull different from the benchmark families that taper both ends to a knife edge.

## 2. The closed forms

With `c = (B/2)/L**q` and `tau = tan(beta)`:

| quantity | expression |
|---|---|
| section area, `x < x0` | `2 T b - tau b**2` |
| section area, `x > x0` | `T**2 / tau` |
| volume | `2T int_0^x0 b dx - tau int_0^x0 b**2 dx + (T**2/tau)(L - x0)` |
| waterplane area | `2 int_0^x0 b dx + 2 (T/tau)(L - x0)` |
| `LCB` | `[2T int_0^x0 x b dx - tau int_0^x0 x b**2 dx + (T**2/(2 tau))(L**2 - x0**2)] / V` |
| section first moment about the waterline | `-T**2 b + T tau b**2 - (1/3) tau**2 b**3` |
| `VCB` | `[int_0^x0 (above) dx - (T**3/(3 tau))(L - x0)] / V` |
| wetted V-bottom | `sec(beta) * Awp` |
| wetted topside | `2 int_0^x0 (T - b tau) sqrt(1 + b'(x)**2) dx` |
| wetted transom | `T**2 / tau` |

Every `int b**k dx` is a power integral, so the volume, the waterplane and both centres
are closed form. Two of them are worth spelling out because they are where a plausible
mistake hides:

* the **section first moment** uses the triangle centroid at *one third of the way up
  from its base*, i.e. `-T + (2/3) b tau`, not `-T + b tau/3`. The wrong choice is a
  23 % error in `VCB` and it survives every sanity band `0 < VCB < T`, which is how it
  was caught here: the first version of the closed form disagreed with the mesh by
  -22.9 % while volume and waterplane area agreed to 0.04 %.
* the **wetted topside is not its vertical projection**. The side follows the chine
  curve, so its surface element carries `sqrt(1 + b'**2)`; with `q < 1` this adds 2.6 %
  of the topside area (`5.037 m**2` against `4.908 m**2` on `chine_immersed`) and it is
  why the projected estimate sat 0.47 % below the mesh until the factor was put in.
  That integral has no elementary closed form once `b'` diverges at the stem, so it is
  evaluated by composite Simpson in the substituted variable `x = x0 t**(2/q)`, which
  removes the endpoint singularity. It is therefore *quadrature*, not closed form, and
  it is labelled as such.

## 3. Results

`chine_immersed` (`T = 0.55 m`, chine wetted to the transom) and `chine_dry`
(`T = 0.40 m`, waterline on the V aft of `x0 = 7.194 m`), both at `81 x 17`:

| quantity | closed form | mesh | mesh - closed form | Simpson quadrature | quadrature - closed form |
|---|---|---|---|---|---|
| `displacement_volume` [m³] | 6.989458 | 6.986746 | -0.0388 % | 6.989457 | -0.0000 % |
| `waterplane_area` [m²] | 18.750000 | 18.745511 | -0.0239 % | -- | -- |
| `wetted_surface_area` [m²] | 25.670434 | 25.662552 | -0.0307 % | -- | -- |
| `LCB` [m] | 5.810983 | 5.812984 | +0.0344 % | 5.810984 | +0.0000 % |
| `VCB` [m] | -0.204943 | -0.204917 | -0.0129 % | -0.204943 | -0.0000 % |
| `LCF` [m] | 6.153846 | 6.155233 | +0.0225 % | -- | -- |

`chine_dry`:

| quantity | closed form | mesh | mesh - closed form | Simpson quadrature | quadrature - closed form |
|---|---|---|---|---|---|
| `displacement_volume` [m³] | 4.199655 | 4.197745 | -0.0455 % | 4.199654 | -0.0000 % |
| `waterplane_area` [m²] | 17.978946 | 17.974094 | -0.0270 % | -- | -- |
| `wetted_surface_area` [m²] | 21.636391 | 21.628910 | -0.0346 % | -- | -- |
| `LCB` [m] | 5.600142 | 5.602486 | +0.0419 % | 5.600142 | +0.0000 % |
| `VCB` [m] | -0.141550 | -0.141526 | -0.0172 % | -0.141550 | -0.0000 % |
| `LCF` [m] | 6.029457 | 6.030851 | +0.0231 % | -- | -- |

The Simpson column is a third route: it integrates the *sections* (the same
`section_area` the loft is built from) at 8192 intervals and never touches the mesh. It
agrees with the paper integrals to `1e-5`, which is what makes the closed form
trustworthy independent of the mesh; the mesh is then a separate discretisation of the
same geometry.

Coefficients, closed form: `chine_immersed` `Cb = 0.4236`, `Cp = 0.7606`, `Cwp = 0.625`,
`VCB/T = 0.627`, transom area `0.9189 m²`; `chine_dry` `Cb = 0.3500`, `Cp = 0.8528`,
`Cwp = 0.5993`, `VCB/T = 0.646`, transom area `0.4924 m²`. `Cwp = 1/(q+1)` exactly when
the chine is immersed everywhere, which is the anchor below.

![verification](artifacts/planing_hull/planing_verification.png)

(a) sections of both regimes with the waterline; (b) mesh error against the closed form
under refinement; (c) error budget per quantity and case.

## 4. Anchors that do not depend on the integrals

| anchor | closed forms | paper formula | agreement |
|---|---|---|---|
| volume, chine immersed everywhere, `q = 0.6` | 6.989458 | `2T c L**(q+1)/(q+1) - tau c**2 L**(2q+1)/(2q+1)` = 6.989458 | exact |
| `Cwp`, chine immersed everywhere | 0.625 | `1/(q+1)` = 0.625 | exact |
| `Cwp`, straight chine `q = 1`, chine immersed | 0.500 | 1/2 | exact |
| transom area, chine dry | 0.492429 | `T**2/tau` = 0.492429 | exact |

## 5. Refinement

`nx x nz` refined by `2n-1` from `21 x 9` to `641 x 257`:

| nx | nz | triangles | volume error | waterplane error | wetted-area error |
|---|---|---|---|---|---|
| 21 | 9 | 696 | -0.33917 % | -0.21148 % | -0.17146 % |
| 41 | 17 | 2672 | -0.11212 % | -0.07136 % | -0.07756 % |
| 81 | 33 | 10464 | -0.03719 % | -0.02394 % | -0.03597 % |
| 161 | 65 | 41408 | -0.01228 % | -0.00800 % | -0.01692 % |
| 321 | 129 | 164736 | -0.00404 % | -0.00266 % | -0.00776 % |
| 641 | 257 | 657152 | -0.00134 % | -0.00088 % | -0.00351 % |

The observed order is **1.60**, not 2, and that is expected rather than a defect: the
trapezoidal rule in `x` sees `b(x) = c x**q` with an infinite slope at the stem, and for
such an integrand the global error is `O(h**(q+1)) = O(h**1.6)`. The straight-chine case
(`q = 1`) is the control: its worst error is 0.0072 % against 0.039-0.074 % for the
`q < 1` hulls at the same grid, i.e. the loss of order is a property of the knife bow and
not of the module. Even at the coarse default the worst disagreement over the four cases
is `7.4e-4`, three orders inside the 1 % gate the module enforces as a validation check.

## 6. What this is not

* **No resistance prediction.** Savitsky's lift, drag and running-trim relations are a
  dynamic problem: they need speed, weight and a trim equilibrium. The geometry and the
  wetted length/beam those relations consume are here; the equilibrium solve is not.
* **No freeboard, sheer, strakes, steps, spray rails, or longitudinal camber.** The
  mesh stops at the waterline, so above-water geometry has no representation.
* **No appendages, no structural or class-rule checks.** Deadrise variation along the
  length is not supported; `beta` is constant, as the Savitsky family assumes.
* The wetted topside area carries a 1-D quadrature (section 2), not a closed form, and
  the transom is triangulated as a fan rather than a loft, so its area is the planar
  polygon area of the last section.

The wetted-area quadrature was re-audited against a high-precision reference and came
out 0.85 % low at `q >= 2` and `inf` on four valid hulls at `q = 0.05`; the rule is now
adaptive Gauss-Legendre with an overflow-free arc element and agrees with the reference
to `2e-14` over 420 cases (`docs/planing_numerical_robustness.md`). The values tabulated
above are unchanged because their `q = 0.6` hulls had already converged.

## 7. Where it plugs in

`tensorcad.planing` returns the same `PlaningHullModel` shape as the other families
(`mesh`, `offsets`-free, `hydrostatics`, `validation`, `metadata()` with the parameters,
the generator description and the analytic reference side by side), so a parameter
sweep over `(L, B, T, beta, q)` is a loop over `generate_planing_hull`. The metadata
records `engineering_status = "geometry and hydrostatics only"` so a downstream consumer
cannot mistake a lofted bottom for a resistance prediction.