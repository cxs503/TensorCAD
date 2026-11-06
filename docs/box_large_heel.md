# Beyond the wall-sided window: an exact righting-arm reference for a rectangular barge

`tensorcad.gz` carries a one-line closed form for a box barge,
`GZ(phi) = sin(phi)(GM + BM/2 tan^2 phi)`, and it is exactly right -- inside its window. Past
that window the heeled waterline has reached a bilge or a deck edge, the formula stops
representing the ship, and the module says so and refuses
(`box_righting_arms_analytic(..., depth=...)` raises `GzError`). That refusal is honest and it is
also where large-angle stability actually lives: the peak of `GZ`, the range of stability and the
area under the curve are all *beyond* the limit. Measured on the barge below, the closed form is
**0.0955 m out at 40 degrees** -- 56 % of the
righting arm -- so extrapolating it is not an option.

What was missing out there was not another solver. `tensorcad.gz` integrates the heeled mesh with
Newton on the wetted area, `tensorcad.righting_arm` (PR #40) does the same job face-by-face with
bracket-and-bisect, and the two agree to round-off. Two integrators on the same mesh agreeing is
*not* evidence about the physics: it cannot see a modelling error they share. This page adds the
third answer that can.

## 1. What is added

`tensorcad.box_gz` -- exact righting arms for a rectangular section, valid in all four regimes up
to 89 degrees, in the same style as the rest of the package (standard library only, no mesh):

| symbol | what it answers |
|---|---|
| `box_section(breadth, draft, depth, heel_deg)` | waterline height, immersed area and the buoyancy centre of the heeled rectangle holding `B*T`, by polygon clip and shoelace sums |
| `box_gz(breadth, draft, depth, kg, heel_deg)` | `GZ = KN - KG sin(phi)` with `KN = Y_B`, plus the regime the section is in |
| `box_gz_curve(...)` | the whole curve as `BoxGzPoint`s |
| `box_regime_bounds_deg(breadth, draft, depth)` | the heels of deck-edge immersion, bilge emergence and of losing both walls |
| `box_vanishing_angle_deg(breadth, draft, depth, kg)` | where `GZ` returns to zero, or `inf` if it does not below 89 degrees |

The four regimes, and why they have to be named: the *shape* of the immersed section changes at
each boundary, and a reference is only trustworthy if it is continuous across them.

| regime | waterline cuts | entered when |
|---|---|---|
| `wall_sided` | both side walls | `phi = 0` |
| `deck_immersed` | the deck and the rising-side wall | `tan(phi) -> 2(D-T)/B` |
| `bilge_emerged` | the bottom and the dipping-side wall | `tan(phi) -> 2T/B` (deck permitting) |
| `bilge_and_deck` | the bottom and the deck | both of the above |

## 2. Method: exact in exact arithmetic

The rectangle is written in the body frame (`|y| <= B/2`, `0 <= z <= D`) and rotated into the
world frame with the same convention as `tensorcad.righting_arm`
(`Y = y cos phi + z sin phi`, `Z = -y sin phi + z cos phi`). At a trial waterline the immersed
part is the rectangle clipped against `Z <= c` -- a convex polygon of at most five vertices,
obtained by one Sutherland-Hodgman pass. Its area and first moments are the shoelace sums, which
are *exact* polynomials, not quadrature. Then the only unknown is `c`: holding the immersed area at
`B*T` fixes it, `A(c)` is strictly increasing, and `dA/dc` is the length of the waterline chord, so
Newton converges quadratically to a residual at the last bit. Nothing is discretised, so for this
geometry the "skipped discretisation error" reported elsewhere in the package is zero.

The frame choice is deliberate: clipping the rotated polygon against a *horizontal* line keeps the
arithmetic well conditioned at zero heel, where the body-frame form `z = h + y tan(phi)` has an
infinite slope.

## 3. Results

Four barges, chosen so that both regime orderings appear (`W=wall_sided`, `D=deck_immersed`,
`B=bilge_emerged`, `X=bilge_and_deck`). The four tables below are rewritten in place by
`tools/verify_box_large_heel.py` from its own run (matched on their header rows), so none of
their numbers is typed in; `tests/test_box_gz_exact.py` re-derives them from
`artifacts/m4864/box_large_heel.json` and fails if the page and the evidence part ways.

| geometry | heel of deck-edge immersion | bilge emergence | both | GZ peak | vanishing angle |
|---|---|---|---|---|---|
| doc barge (L=10.0 B=2.0 T=0.8 FB=0.4 KG=0.55) | 21.801 deg | 41.987 deg | 41.987 deg | 0.175849 m at 46 deg | beyond 89 deg |
| wide shallow (L=6.0 B=4.0 T=0.5 FB=0.3 KG=0.9) | 8.531 deg | 14.931 deg | 14.931 deg | 0.442404 m at 18 deg | 55.355 deg |
| square deep (L=3.0 B=3.0 T=1.0 FB=0.5 KG=1.2) | 18.435 deg | 36.870 deg | 36.870 deg | 0.033989 m at 22 deg | 34.364 deg |
| tall barge (L=8.0 B=2.0 T=0.8 FB=1.6 KG=0.6) | 60.945 deg | 38.660 deg | 60.945 deg | 0.592828 m at 85 deg | beyond 89 deg |

Note the `tall barge` row: with enough freeboard the deck never interferes, the section goes
`W -> B -> X`, and there the computed bilge-emergence heel has to equal the classical textbook
`atan(2T/B)`. It does, which is the cleanest independent check of the new solver:

| geometry | wall-sided window ends at | computed bilge emergence | classical atan(2T/B) | difference |
|---|---|---|---|---|
| doc barge | 21.8014 deg | 41.9872 deg | 38.6598 deg | +3.3274 deg |
| wide shallow | 8.5308 deg | 14.9314 deg | 14.0362 deg | +0.8952 deg |
| square deep | 18.4349 deg | 36.8699 deg | 33.6901 deg | +3.1798 deg |
| tall barge | 38.6598 deg | 38.6598 deg | 38.6598 deg | -0.0000 deg |

When a deck immerses first (the other three barges) the bilge stays wet **longer** than
`atan(2T/B)`, because capping the volume above the deck raises the waterline. That difference is
real, not a defect, and it is why the classical limit is only a bound.

### 3.1 The reference against the closed form, and both integrators against the reference

| geometry | wall-sided form vs exact, inside its window | wall-sided form at 40 deg | Newton vs exact, 0-85 deg | bisect vs exact, 0-85 deg |
|---|---|---|---|---|
| doc barge | 8.33e-17 m | 0.0955 m | 3.33e-16 m | 3.73e-13 m |
| wide shallow | 1.11e-16 m | 1.6745 m | 3.44e-15 m | 1.21e-12 m |
| square deep | 2.50e-16 m | 0.2309 m | 4.44e-16 m | 6.89e-13 m |
| tall barge | 1.11e-16 m | 0.0007 m | 1.82e-14 m | 4.64e-13 m |

Inside its window the new solution reproduces the wall-sided closed form to
2.5e-16 m -- round-off, on righting arms of order 0.1 m. That is the
check that licenses using it outside the window. Outside, the two merged integrators follow it over
the whole 0-85 degree sweep: Newton to 1.8e-14 m, bracket-and-bisect to
1.2e-12 m. Until this page, neither had an independent reference past the
wall-sided limit.

| geometry | area residual of the reference | mesh volume residual (bisect route) |
|---|---|---|
| doc barge | 2.78e-16 (relative) | 6.36e-14 (relative) |
| wide shallow | 2.22e-16 (relative) | 9.42e-14 (relative) |
| square deep | 2.96e-16 (relative) | 1.19e-13 (relative) |
| tall barge | 2.78e-16 (relative) | 1.73e-13 (relative) |

![GZ across all four regimes](../artifacts/m4864/box_large_heel_gz.png)

*Left*: the exact curve for the doc barge (black), the wall-sided closed form
(red dashed, valid only in the shaded window), and the two mesh integrators sampled every five
degrees -- they sit on the curve. The dotted verticals are the deck-edge immersion and bilge
emergence heels. *Right*: `|integrator - exact|` per heel on a log scale for all four geometries;
the dotted line is 1e-9 m.

![heeled sections](../artifacts/m4864/box_large_heel_sections.png)

The reference's own geometry: the section at 0, 21.8, 42, 60 and 80 degrees with the immersed part
shaded and the buoyancy centre marked. The regime label under each panel is what
`box_gz(...).regime` returns.

## 4. Error judgement

Key quantities against their references, and whether the target (<3 %) is met:

| quantity | reference | measured | verdict |
|---|---|---|---|
| `GZ` inside the wall-sided window | `sin(phi)(GM + BM/2 tan^2 phi)` | <= 2.5e-16 m on arms of 0.02-0.11 m, i.e. < 1e-13 % | **yes, far below 3 %** |
| `GZ(0..85 deg)`, Newton route | `tensorcad.box_gz` exact | <= 1.8e-14 m on arms up to 0.59 m (~1e-13 %) | **yes** |
| `GZ(0..85 deg)`, bisect route | `tensorcad.box_gz` exact | <= 1.2e-12 m (~1e-11 %) | **yes** |
| bilge-emergence heel, deck-free geometry | classical `atan(2T/B)` | agrees to < 1e-9 deg | **yes** |
| immersed area at every heel | `B*T` | <= 3.0e-16 relative | **yes** |
| wall-sided form at 40 deg (the old documented number) | exact | 0.0955 m on 0.171 m = 56 % | **no, and that is the finding**: the formula is not usable there |

The last row is the reason for the work, not a regression: the closed form's own 21.8 degree
validity limit was already documented, and 9.55 cm reproduces the 9.6 cm on record.

### Negative control

A comparison that cannot fail is not a comparison. Shifting `KG` by 1 mm moves `GZ(40 deg)` by
6.43e-04 m and a 0.5 % displacement shift moves `KN` by
1.04e+00 m -- six orders of magnitude above the
1e-9 m agreement claimed above, so the sweep would have flagged either.

## 5. Note points

* **Double rotation.** The first version of `box_gz` computed `KN = Y_B cos(phi) + Z_B sin(phi)`
  from an *already rotated* centroid, giving `GZ(5 deg) = 0.0577 m` against the closed form's
  `0.0234 m`. The wall-sided test caught it immediately. `KN` is the world abscissa of the buoyancy
  centre; rotating it again is the classic frame error, and it is the reason every check here is
  against a *frame-independent* closed form first.
* **Regime classification must sort the two walls.** `tan(phi)` changes sign at negative heel, so
  the signed body-frame expression swaps the walls: a port heel of 25 degrees was labelled
  `wall_sided` while the starboard heel was `deck_immersed`. The righting arm was correct to
  round-off throughout -- only the label was wrong -- which is exactly the kind of defect that a
  symmetry test (`GZ(-phi) = -GZ(phi)` *and* equal regime) finds and a value test does not.
* **Validity.** The section model ends at 89 degrees (the rectangle degenerates to a line) and has
  no appendages, no flooding, no free surface (that stays in `gz.with_freeboard` and
  `gz.free_surface_correction`). Two of the four barges never vanish below 89 degrees and that is
  reported as `inf`, not as an extrapolated angle.
* **What this does not prove.** The reference is a 2-D section solution; it says nothing about the
  three-dimensional integrators' handling of *non-box* hulls, where the waterline solve is the hard
  part. The `81x25` hull and the sphere keep their own references (`docs/righting_arm.md`).
