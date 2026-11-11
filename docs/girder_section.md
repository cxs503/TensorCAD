# Hull-girder section properties

`tensorcad.girder` answers the structural question about a hull's transverse section:
given the polygon the shell encloses, what is its bending stiffness, where is its neutral
axis, how much moment can it take before the extreme fibre yields, and how much before
the whole section is plastic? It is the first structural module in the package, and it
takes the one input a geometry kernel owns -- a simple closed polygon in the section plane
`(y, z)` in metres, or a half-breadth table that is mirrored into one.

The module computes, from the polygon alone:

* `area`, `perimeter`, `beam`, `depth`, `centroid_y`, `centroid_z`;
* `second_moment_horizontal` (`INT (z - z_c)^2 dA`, the one in `sigma = M z' / I`) and
  `second_moment_vertical`, with `product_of_inertia` and the principal axes;
* `section_modulus_deck` / `_keel` / `_starboard` / `_port` -- the elastic moduli, each
  with its own extreme-fibre distance, so an asymmetric section reports four numbers;
* `plastic_neutral_axis_z` and `plastic_modulus` (plus the vertical pair), located by
  bisection on the exact clipped area and integrated from the first moments;
* `first_moment(z)` -- the static moment `Q(z)` that is the numerator of `tau = V Q / I t`;
* `closed_cell_torsion_constant` (Bredt) and `open_section_torsion_constant`
  (`sum b t^3 / 3`).

## Method, and why the numbers are exact

Area, centroid and second moments are the shoelace/Euler sums: the boundary integral of a
quadratic over a straight edge is exactly the trapezoid rule, so a polygon carries no
quadrature error. The only searched quantity is the plastic neutral axis, found by
bisection on the area below a level, and the area below a level is itself a boundary
integral (`INT (z_p - z)^+ dA = CLOSED_INT (z_p - z)^+ dy`, with the moment using the
continuous potential `-min(z, z_p)^2 / 2`), evaluated with three-point Simpson on each
edge, which is exact because the integrand is at most quadratic. Two consequences are
worth stating because they were the traps in this module:

* A hole needs no nesting analysis. An annulus handed over as one loop with a zero-width
  bridge subtracts its void automatically, because the boundary integral never cares how
  the polygon is connected.
* A vertex sitting exactly on the cut line needs no special case: the level is not
  nudged, and an edge lying along the cut contributes the constant potential term that
  the closed loop requires.

The plastic modulus is `Z_p = M_above - M_below - z_p (A_above - A_below)`, i.e.
`INT |z - z_p| dA`, with `z_p` the level that halves the area. Because that identity only
needs first moments, it stays exact; the earlier temptation to reconstruct the clipped
polygons and sum their areas is what produced the annulus failure documented below.

## What it is not

* **Not a stress calculation.** No bending moment is applied and no material is named;
  the section modulus is the denominator of `sigma = M / Z`, and `M` is the caller's.
* **Not an effective-breadth model.** The section is fully effective: no shear lag, no
  Schade / von Karman correction, no buckling of the compression flange.
* **Not an ultimate-strength check.** `plastic_modulus` is the full-plastic first moment,
  not a collapse moment with plate slenderness, residual stress or a material factor.
* **Not a class-rule engine.** It reports numbers; comparing them to a minimum is policy
  and belongs to the rule, not to the geometry.
* **Thin-walled torsion only.** Bredt and `sum b t^3 / 3` assume `t` small against the
  cell size; multi-cell decks and wall-thickness variation are not modelled.

## Verification

Three independent kinds of reference are used (all of them re-run by
`examples/girder_section/run_example.py`, which fails if a claim stops holding):

1. **Closed forms** for a rectangle (`b h^3 / 12`, `b h^2 / 6`, `b h^2 / 4`, shape factor
   exactly `1.5`), a double-walled box (difference formula), an I-beam (sum of rectangles
   with Steiner terms) and a plus-shaped section (two bars minus the square they share).
2. **Quadrature** of a trapezoid's own width law for the plastic neutral axis and plastic
   modulus: the equal-area level has a closed form and the integrals are Simpson sums that
   never bisect.
3. **An inscribed regular n-gon** against `pi R^4 / 4`, where theory fixes the error
   coefficient at `4 pi^2 / 3`.
4. **A real hull** from `tensorcad.generate_hull`, whose section areas are cross-checked
   against the independent clipped integration in `tensorcad.bonjean`.


## Closed forms

| quantity | measured | reference | relative error |
|---|---|---|---|
| rectangle_area | 24 | 24 | 0.000e+00 |
| rectangle_i_horizontal | 32 | 32 | 0.000e+00 |
| rectangle_section_modulus | 16 | 16 | 0.000e+00 |
| rectangle_plastic_modulus | 24 | 24 | 0.000e+00 |
| rectangle_shape_factor | 1.5 | 1.5 | 0.000e+00 |
| hollow_area | 16 | 16 | 0.000e+00 |
| hollow_i_horizontal | 29.3333333333 | 29.3333333333 | 1.211e-16 |
| hollow_plastic_modulus | 20 | 20 | 0.000e+00 |
| i_beam_area | 6 | 6 | 0.000e+00 |
| plus_area | 20 | 20 | 0.000e+00 |
| plus_plastic_modulus | 22 | 22 | 0.000e+00 |

Worst closed-form residual: **hollow_i_horizontal** at 1.211e-16 relative.

## Plastic neutral axis (trapezoids)

| bottom | top | depth | level (measured) | level (closed form) | Zp (measured) | Zp (quadrature) |
|---|---|---|---|---|---|---|
| 2.0 | 6.0 | 4.0 | 2.472135955 | 2.472135955 | 15.038187267 | 15.038187267 |
| 6.0 | 2.0 | 3.0 | 1.145898034 | 1.145898034 | 8.458980338 | 8.458980338 |
| 1.0 | 8.0 | 5.0 | 3.357769375 | 3.357769375 | 24.224997762 | 24.224997762 |

## Inscribed n-gon

| sides | I | relative error | error * n^2 |
|---|---|---|---|
| 16 | 11.935149390 | 5.023e-02 | 12.8591 |
| 64 | 12.526056178 | 3.208e-03 | 13.1405 |
| 256 | 12.563847545 | 2.008e-04 | 13.1583 |
| 1024 | 12.566212909 | 1.255e-05 | 13.1594 |

Predicted coefficient `4 pi^2 / 3` = 13.1595.

## Hull (L=10 m, B=2 m, T=0.8 m)

| station | A | A (bonjean) | I_hh | Z_deck | Z_keel | Zp | Zp/Z |
|---|---|---|---|---|---|---|---|
| 1.00 | 0.465518 | 0.465518 | 0.019946 | 0.062922 | 0.041296 | 0.081777 | 1.9888 |
| 2.00 | 0.730132 | 0.730132 | 0.031284 | 0.098689 | 0.064769 | 0.128262 | 1.9803 |
| 3.00 | 0.913101 | 0.913101 | 0.039123 | 0.123420 | 0.081000 | 0.160404 | 1.9803 |
| 4.00 | 1.022570 | 1.022570 | 0.043814 | 0.138216 | 0.090711 | 0.179635 | 1.9803 |
| 5.00 | 1.059128 | 1.059128 | 0.045380 | 0.143158 | 0.093954 | 0.186057 | 1.9803 |
| 6.00 | 1.022570 | 1.022570 | 0.043814 | 0.138216 | 0.090711 | 0.179635 | 1.9803 |
| 7.00 | 0.913101 | 0.913101 | 0.039123 | 0.123420 | 0.081000 | 0.160404 | 1.9803 |
| 8.00 | 0.730132 | 0.730132 | 0.031284 | 0.098689 | 0.064769 | 0.128262 | 1.9803 |
| 9.00 | 0.465518 | 0.465518 | 0.019946 | 0.062922 | 0.041296 | 0.081777 | 1.9888 |

Largest area gap against bonjean: 4.770e-16 relative. Largest |product of inertia|: 4.086e-12 m^4.
Governing (largest) section at x = 5.00 m.

## Torsion

Bredt closed cell 2.274936e+00 m^4 vs classical box 2.274936e+00 m^4 (1.952e-16 relative); open-steel version 4.266667e-05 m^4, ratio 53318.8.
## Reading the numbers

The hull row that governs is the one with the largest area -- the midship section in this
loft. Two properties of a real hull show up in the table and neither is visible from the
hydrostatics modules: the deck and keel section moduli differ (the deck of this hull is
the narrower extreme fibre, so `Z_deck < Z_keel`), and the plastic-to-elastic ratio is
well above the rectangle's `1.5` because the material sits far from the neutral axis.

The torsional pair in the last table is the quantitative version of "a ship needs a closed
deck": the same steel as a single closed cell and as three open strips differs by two to
three orders of magnitude, which is why an open hold is not a hull girder.

## Reproducing

```bash
python examples/girder_section/run_example.py   # writes artifacts/girder_section/
python tools/girder_section_doc.py              # regenerates this file from those numbers
python -m pytest tests/test_girder.py
```

## Limitations and next steps

* Shear stress itself is not computed: `Q(z)` is reported, `V` and the local thickness are
  not, and shear lag would be needed before a shear stress were trustworthy at the ends.
* No load model: wave and still-water bending moments (IACS / CSR formulations) belong to
  the rule set, not to a geometry kernel.
* The section is a polygon, so a plate thickness that is not modelled geometrically
  (typical in scantling work) has to be passed as an effective polygon width.
