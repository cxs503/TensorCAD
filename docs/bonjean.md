# Bonjean curves: a station's immersed area against the waterline

Every other quantity in this package is a whole-hull number. `compute_hydrostatics`
integrates the wetted surface and returns one displaced volume, one waterplane, one
centre; `section_area_at` returns the area of a station the mesh happens to have a
vertex ring on, at the full depth of that ring. Neither can answer the question a
naval architect asks next:

> at station *x*, how much area of the section is under the waterline *z*?

That table — area against waterline, one curve per station — is the **Bonjean
curve**. It is the input to damage stability (the lost buoyancy of one flooded
compartment is the area between two Bonjean curves) and to still-water shear force
and bending moment (the buoyancy load per unit length is `rho g A(x)` at the design
waterline, and the section modulus of the hull girder has to be checked against the
integral of `A(x) * (x - x_section)`).

`tensorcad/bonjean.py` adds it, and this page is the measurement of it: every number
below is written by `python examples/bonjean/run_example.py` into
`docs/bonjean/{results.json,plot.png}`, which exits non-zero if a claim stops
holding. `tests/test_bonjean.py` pins the same behaviour.

## Method

A cut of a closed triangle mesh by the plane `x = station` is a set of closed
polylines: each straddling triangle contributes exactly one segment, and the
segments chain end to end because the mesh is a manifold. Those segments are
matched **by their shared edge key**, so the coordinates on both sides of a shared
edge are the same floats rather than "close enough" floats — chaining cannot drift
apart on a 200 k-triangle mesh.

The plane is placed at `station ± 1e-9 * scale`, the sign chosen to stay inside the
body. Reason: the hull loft puts vertex rings exactly on the station planes, and a
plane through vertices is a degenerate cut (a triangle with two vertices on the
plane has an edge, not a crossing). The nudge converts that degeneracy into an
ordinary crossing whose area differs from the exact station by `O(1e-9)` of the
station spacing — nine orders of magnitude below the discretisation error of any
mesh this package produces. It is stated, not hidden, and it is why the end stations
still report their (zero-area) stem section.

The immersed area is then the area of the cut polygon clipped to `z <= waterline`
(Sutherland–Hodgman, then shoelace), with the first moments taken by the same
pairing so the centroid comes out of the area pass rather than a second one.

## Result 1 — exact on a barge

A wall-sided box has a rectangular section, so nothing is discretised:

| quantity | measured | exact | error |
|---|---|---|---|
| section area, midship, full draft | 1.600000000 m² | `B * T` = 1.6 | 0 (bit-exact) |
| section centroid `z` | −0.400000000 m | `−T/2` | 0 (bit-exact) |
| section perimeter (cut polygon) | 5.600000000 m | `2B + 2T` | 0 (bit-exact) |
| `∫ A(x) dx` over 11 stations | 16.000000000 m³ | `L B T` = 16 | 2.2e-16 relative |

## Result 2 — the parametric hull against its closed form

The generator's half breadth is `(B/2) sin(pi x/L)^0.7 sqrt(zeta/draft)`, so the
section area is analytic and the mesh never saw it:

    A(x, zeta) = (2/3) * B * sin(pi x / L)^0.7 * zeta^(3/2) / sqrt(draft)

(`zeta = z + draft`, measured up from the keel. At midship and the design waterline
this reduces to `(2/3) B T`, the `analytic_maximum_section_area` the hydrostatics
module already compares against — but it now varies with station *and* waterline.)

Worst relative error over 5 stations × 4 waterlines on the default 41×13 grid:
**−5.31 %**, at the shallowest waterline `z = −0.6 m`. That is not a wrong formula,
it is the piecewise-linear skin standing in for `sqrt(z)` where only three of the
thirteen levels are immersed — and it converges:

| nz | section area error at midship, `z = −0.6 m` |
|---:|---:|
| 13 | −5.308 % |
| 21 | −2.539 % |
| 41 | −0.924 % |
| 81 | −0.333 % |

Successive error ratios 2.09 / 2.75 / 2.77 — monotone, second order once the waterline
is resolved by more than three segments. All errors are negative: a linear skin under
a concave profile always under-estimates, never over-estimates, and the test asserts
that sign rather than only the magnitude.

## Result 3 — two independent integrations of one body

The loft is linear in `x` between stations, so the exact volume of the
piecewise-linear body is the trapezoid rule on its exact station areas. Therefore the
section integral must equal the **surface-flux** displacement from
`compute_hydrostatics` (a sum of tetrahedra that never builds a section). It does:

| nx (nz = 41) | Simpson over stations | trapezoid over stations | surface flux | vs analytic |
|---:|---:|---:|---:|---:|
| 11 | 7.469350 m³ (−0.577 %) | 7.365083 m³ | 7.365083 m³ | −1.965 % |
| 21 | 7.493080 m³ (−0.261 %) | 7.461080 m³ | 7.461080 m³ | −0.687 % |
| 41 | 7.500448 m³ (−0.163 %) | 7.490606 m³ | 7.490606 m³ | −0.294 % |
| 81 | 7.502720 m³ (−0.133 %) | 7.499692 m³ | 7.499692 m³ | −0.173 % |

`sections − flux` is **1.1e-16 relative** at every grid: the two integrations agree to
round-off, so the two are measuring the same body and the residual against
`analytic_volume` = 7.512698 m³ is surface sampling, not quadrature. Simpson lands
closer to the analytic volume than the (exact-for-the-mesh) trapezoid because the
station areas are a slightly low *sample* of the analytic area; both numbers are
reported so the choice is visible rather than selected.

## What this buys downstream

* **Damage stability**: the lost buoyancy of a compartment is the area between the
  Bonjean curves at its two bounding stations, from the damage waterline up —
  available now, per station, without re-running a whole-hull hydrostatics.
* **Longitudinal strength**: `A(x)` at the design waterline is the buoyancy curve;
  the still-water shear force and bending moment are its running integrals against
  the weight curve. Before this module the only station-wise number available was
  the area of the stations the mesh happened to have vertices on, at one waterline.
* **The bounding-box shortcut is now measurable** rather than assumed: the section
  table shows which stations are under-resolved.

## Limits, honestly

* One waterline cut at a time; a healed/healed *and* trimmed waterplane (an inclined
  plane) is not supported — the cut plane is `x = const` and the waterplane is
  `z = const`. Heel/trim hydrostatics remains the open geometry item for GZ curves.
* Multi-loop sections (a twin-hull, a flooded compartment with a hole) are summed by
  area, not resolved as separate compartments; `SectionProperties.loops` reports the
  count so a caller can refuse.
* The nudge means the reported station is a nanometre-scale shift of the requested
  one; the section at a *discontinuity* in `x` (a transom, a step) is reported from
  the inside of the body.
* Pre-existing limitation, not new: both the closed form and the mesh share the
  `sqrt(zeta/draft)` profile, so this page says nothing about whether that profile is
  the right ship — only that the mesh represents it to the accuracy tabulated above.

## Error judgement against the standards used in this package

| key quantity | reference | error | < 3 % ? |
|---|---|---:|---|
| barge section area, centroid, `∫A dx` | exact closed form | 0 (bit-exact) | yes |
| hull section area, default 41×13 grid, worst waterline | closed form (2/3) | −5.31 % | **no** (and stated as convergence-limited) |
| hull section area, nz = 41, same waterline | closed form | −0.92 % | yes |
| hull section area, nz = 81 | closed form | −0.33 % | yes |
| displacement from the Bonjean integral, nx = 41, nz = 41 | `analytic_volume` | −0.16 % | yes |
| displacement from the Bonjean integral, nx = 81 | `analytic_volume` | −0.13 % | yes |
| section integral vs surface flux | `compute_hydrostatics` | 1.1e-16 | yes (identity) |

The one "no" is the shallow waterline on the coarsest vertical grid, and it is a
resolution statement, not a formula statement: the same formula is 0.33 % accurate
at `nz = 81`.