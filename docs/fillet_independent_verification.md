# Independent verification of the fillet/chamfer kernel (`tensorcad.fillet`)

**Mission** m16561 (selfstart, verifier seat) · **Module under review**
`tensorcad/fillet.py` (landed by m16453, merge `1a40aa0`, anchor `GAP-CAD-012`) ·
**Code changed** none in the kernel (verification-only mission) ·
**Reproduce** `python examples/fillet_independent_verification/run_example.py`
(53 checks, exit non-zero on regression) · **Tests**
`tests/test_fillet_independent_verification.py` (12 cases) · **Artifacts**
`artifacts/fillet_independent_verification/independent_verification.json`,
`fig1_independent_routes.png`, `fig2_closed_form_agreement.png` ·
**Machine-readable table** `examples/fillet_independent_verification/outputs/{report.md,checks.json}`

---

## 1. Case: what was verified, and why it needed verifying

`tensorcad.fillet` rounds (or chamfers) the corners of a 2-D section and sweeps the result into a
closed solid. It claims *exact* area, perimeter, centroid and second moments, an *exact* swept
volume for a radius that varies linearly along the axis, and a *feasibility boundary*
`r <= L / (tan|phi1|/2 + tan|phi2|/2)` per edge. Those claims are what a hull designer consumes:
the section moments are the waterplane inertia behind `BM = I_T / V`, and the swept volume is a
displacement.

The module's own defence is a cross-check between **two derivations that both live inside the
same file**: the Green's-theorem path integral over lines and arcs, and the per-corner closed form
`dA = -sign(phi) r^2 (tan(|phi|/2) - |phi|/2)`. That catches a coding slip in one of them, but a
mistake shared by both — a wrong sign convention, a wrong integration domain, a wrong tangency
length — survives it. `docs/fillet_chamfer.md` reports the residual as `1e-16`, which is evidence
of *consistency*, not of *correctness*. This mission supplies the missing independent routes.

| property | how the kernel checks it today | route added here |
|---|---|---|
| area, perimeter | Green over the path vs per-corner closed form (same file) | hand-derived rounded-rectangle closed form; signed-corner law on a re-entrant L |
| centroidal second moments | no reference value published | hand-derived closed form **and** brute-force midpoint integration |
| "exact" vs the mesh | sampled area reported | Richardson ratio over `arc_points` (the exact value must be the mesh's limit) |
| swept volume | closed form quoted in the docstring | 5-point Gauss-Legendre of `A(r(x))` + divergence-theorem mesh volume with refinement |
| feasibility | exception when over the limit | circle limit `r = r_max` at machine precision, over-limit refusal, sliver/notch/near-limit cases |
| tangency `g1_error_deg` | analytic identity of the construction | see §5 finding F1 |

Reference geometry is a `4 m x 2 m` rectangle with corner radius `r`, a `2 m x 1 m` L section with
one 90 deg re-entrant corner (base area 3 m², 5 convex + 1 reflex corner), a `2 m x 2 m` square at
its exact limit `r = 1 m` (the circle), and a `3 m x 1 m` section with an interior slot
(two reflex corners facing each other). Sweeps are `L = 3 m`.

## 2. Method: five routes that share no code with the kernel

* **R1 hand closed form.** Rounded rectangle of `w x h` with radius `r`: `A = wh - (4-pi) r^2`,
  `P = 2(w+h) - (8-2pi) r`. For the centroidal second moment a corner bite removes its `r x r`
  square and adds the quarter disk back; in the centroidal coordinate the disk centre sits at
  `d = extent/2 - r` and the quarter disk contributes `d^2 pi r^2/4 + 2 d r^3/3 + pi r^4/16`, so
  `I = other*extent^3/12 - 4 (strip - quarter)` with `strip = r[(extent/2)^3 - (extent/2 -
  r)^3]/3`. The re-entrant section uses the signed-corner law directly.
* **R2 brute force.** Midpoint rule on the rounded body's *indicator function*: it never sees an
  arc, a tangent length or a Green integral — only "is (u,v) inside?". `1600 x 1600` cells for the
  published numbers, `400 x 400` in the test.
* **R3 quadrature.** `A(r(x))` is exactly quadratic for a linear radius profile, so five-point
  Gauss-Legendre integrates it exactly; the mesh volume comes from the divergence theorem on the
  triangle soup.
* **R4 Richardson.** The chorded outline is what a mesh receives, so the "exact" number must be
  the `arc_points -> inf` limit of `sampled_moments()`; a wrong exact value and a correct mesh
  would show up as a non-zero limit or a wrong slope.
* **R5 adversarial geometry.** Circle limit, `r_max + 2 %` refusal, sliver `4 x 0.2` at
  `r = 0.0999`, a slotted section whose two reflex corners face each other, clockwise vs
  counter-clockwise input, and translation of the whole section by up to `1e8` m.

## 3. Results

### 3.1 Rounded rectangle against the hand closed form

| r (m) | area (kernel) | area (closed form) | rel. | I_uu (kernel) | I_uu (closed form) | rel. |
|---:|---:|---:|---:|---:|---:|---:|
| 0.05 | 7.997853981626 | 7.997853981626 | 0.0e+00 | 10.658178007151 | 10.658178007151 | 1.7e-16 |
| 0.40 | 7.862654824574 | 7.862654824574 | 1.1e-16 | 10.164502543893 | 10.164502543893 | 5.2e-16 |
| 0.80 | 7.450619298325 | 7.450619298325 | 3.6e-16 | 8.831924210606 | 8.831924210606 | 0.0e+00 |
| 0.999 | 7.143308609878 | 7.143308609878 | 2.5e-16 | 7.931847653063 | 7.931847653063 | 1.0e-15 |

Perimeter agrees to `<= 3.5e-16` and `I_vv` to `<= 4.2e-16` at every radius
(`2.551995263163` vs `2.551995263163` at `r = 0.4`). Numbers program-injected from
`artifacts/fillet_independent_verification/independent_verification.json`.

### 3.2 Brute force (`r = 0.4`, 1600 x 1600 midpoint cells)

| quantity | brute force | kernel | rel. |
|---|---:|---:|---:|
| area | 7.862700000230 | 7.862654824574 | 5.8e-06 |
| centroid u | 1.999999999940 | 2.0 | 3.0e-11 |
| `int (u-cu)^2 dA` | 10.164647467600 | 10.164502543893 | 1.4e-05 |
| `int (v-cv)^2 dA` | 2.552026511660 | 2.551995263163 | 1.2e-05 |

The brute-force residual is `O(h)` on the curved boundary — it is a property of the grid, not of
the kernel (see Fig. 2's dashed line).

### 3.3 Re-entrant L section and the circle limit

| case | measured | reference | rel. |
|---|---:|---:|---:|
| L area, r=0.1 / 0.3 / 0.5 (m²) | 2.991415926536 / 2.922743338819 / 2.785398163397 | same | 0 / 3.0e-16 / 0 |
| L perimeter, r=0.3 (m) | 7.227433388231 | 7.227433388231 | 2.5e-16 |
| circle-limit area (`r = r_max = 1`), m² | 3.141592653590 | pi | 2.8e-16 |
| circle-limit perimeter (m) | 6.283185307180 | 2 pi | 2.8e-16 |
| `max_fillet_radius(2x2)` (m) | 1.000000000000 | 1 | 2.2e-16 |

As `r` grows the L section **loses** area at its five convex corners and **gains** at the reflex
one; at `r = 0.5` (the edge limit) the result is `2.785398 = 3 - 5(0.1) + 0.1` to round-off, i.e.
the signed-corner law holds through the limit. Requesting `r = 1.02` on the `2 x 2` square raises
`FilletError` naming `the largest feasible size is 1.0` — there is no silent repair.

### 3.4 The "exact" values are the mesh's own limit (R4)

| chords/arc | sampled area - exact (rel.) | error ratio |
|---:|---:|---:|
| 4 | -1.6305e-03 | — |
| 8 | -4.0999e-04 | 3.977 |
| 16 | -1.0265e-04 | 3.994 |
| 32 | -2.5671e-05 | 3.999 |
| 64 | -6.4182e-06 | 4.000 |

Swept mesh volume (`L = 3 m`, `r = 0.4 m`): `1.6305e-03 -> 4.0999e-04 -> 1.0265e-04 -> 2.5671e-05`
for `arc_points = 4, 8, 16, 32`, mean doubling ratio **3.990**; every mesh is watertight,
`V - E + F = 2`, zero zero-area faces. The chordal crease, reported honestly by the module, falls
`22.5 deg -> 2.81 deg` with the same refinement — the exact volume and the mesh volume are two
different numbers and the module never conflates them.

### 3.5 Sweep volume against Gauss-Legendre quadrature

| `r0 -> r1` (m) | exact volume (m³) | Gauss quadrature | rel. | mesh volume rel. |
|---|---:|---:|---:|---:|
| 0.2 -> 0.6 | 23.553628179888 | 23.553628179888 | 0.0e+00 | 1.1e-04 |
| 0.6 -> 0.2 | 23.553628179888 | 23.553628179888 | 1.5e-16 | 1.1e-04 |
| 0.1 -> 0.999 | 23.048969642534 | 23.048969642534 | 0.0e+00 | 2.4e-04 |

Reversing the taper gives the same volume, as it must; `A(r(x))` is exactly quadratic so the
five-point rule is exact and the agreement is at round-off, not at quadrature tolerance.

### 3.6 Adversarial geometry (R5)

| case | `r_max` | sampled vs exact area | orientation preserved | note |
|---|---:|---:|---|---|
| sliver 4 x 0.2, r=0.0999 | 0.100000 | -2.54e-04 | yes | chorded outline stays simple |
| rectangle 4 x 2, r=0.999 | 1.000000 | -2.81e-03 | yes | four arcs nearly meet |
| L section, r=0.5 (limit) | 0.500000 | -1.81e-03 | yes | reflex arc grows material |
| slotted 3 x 1, r=0.29 | 0.300000 | -7.28e-04 | yes | two facing reflex corners |

Clockwise input gives the identical area (`7.862654824574` both ways, rel. 0). All sampled
outlines keep the original orientation and remain simple polygons (checked by segment-segment
intersection counting in the example script).

## 4. Figures

* `artifacts/fillet_independent_verification/fig1_independent_routes.png` — left: chorded section
  area and swept mesh volume vs chords per arc, log-log, with the `O(h²)` slope drawn through;
  right: relative area error under translation, with the `eps * shift / size` line.
* `artifacts/fillet_independent_verification/fig2_closed_form_agreement.png` — rounded-rectangle
  area and `I_uu` against the hand closed form over `r = 0.05 .. 0.999`, plus the brute-force
  `I_uu` level as the horizontal reference.

## 5. Findings

**F1 (non-blocking, documentation).** `SectionShape.g1_error_deg` for a fillet is an *identity of
the construction*: it takes the cosine between the analytic arc-start tangent and the incoming
edge direction, so it is `0.0` by construction and cannot fail. The geometry the caller actually
receives is the chorded outline, which is only G1 to within half the chord angle; the module
reports that number as `sampled_crease_deg` for **solids** (`turn/arc_points`, `22.5 deg` at 4
chords) but a *section* has no equivalent field. Recommendation: state in the docstring that
`g1_error_deg` is analytic, and expose the sampled crease on `SectionShape` too. Not a numerical
defect — the module never claims the chorded outline is G1.

**F2 (documented limit).** The exactness claims are **conditioning-bound**. The integrals are
evaluated in absolute coordinates, so the relative area error grows like `eps * shift / size`:

| translation / body size | 2.5 | 2.5e2 | 2.5e5 | 2.5e7 |
|---|---:|---:|---:|---:|
| relative area error | 0.0e+00 | 2.8e-14 | 2.8e-11 | 7.3e-09 |

For ship sections (`B/2 ~ 20 m`, bilge radius `~ 0.5 m`) the ratio stays near 40 and the error
stays at `1e-15`; a caller who evaluates a millimetre feature at a kilometre-scale offset should
translate the ring to its centroid first. This is now asserted in the test suite so the number
cannot drift silently.

**F3 (methodological, about this mission).** My first hand derivation of the centroidal second
moment was **wrong** — I dropped the sign of the `2 d r^3/3` term, which made the "reference"
disagree with the kernel by `2.76e-02` relative and would have produced a false defect report.
Route R2 (brute force) is what settled it: it reproduced the kernel's `10.16450` to `1.4e-05`
while my formula sat at `9.89144`. Lesson (also for the next verifier): a closed form is only an
anchor once an independent numerical route has confirmed it; when a hand derivation and the code
disagree, the brute force decides, not seniority.

No correctness defect was found in the kernel. Disposition: **approve with notes** (F1, F2).

## 6. Error judgement against the references

| key quantity | reference | measured | error | `< 3 %`? |
|---|---|---:|---:|---|
| rounded-rect area, `r = 0.4` | hand closed form 7.862654824574 | 7.862654824574 | 1.1e-16 | **yes** |
| rounded-rect perimeter, `r = 0.4` | hand closed form 11.313274122872 | 11.313274122872 | 0.0e+00 | **yes** |
| centroidal `I_uu`, `r = 0.4` | hand closed form 10.164502543893 | 10.164502543893 | 5.2e-16 | **yes** |
| centroidal `I_vv`, `r = 0.4` | hand closed form 2.551995263163 | 2.551995263163 | 1.7e-16 | **yes** |
| `I_uu` vs brute force (1600² cells) | 10.164647467600 | 10.164502543893 | 1.4e-05 | **yes** |
| L-section area, `r = 0.3` | signed-corner law 2.922743338819 | 2.922743338819 | 3.0e-16 | **yes** |
| circle-limit area, `r = r_max` | pi | 3.141592653590 | 2.8e-16 | **yes** |
| chordal convergence slope | 2 (order) | 1.9978 (mean ratio 3.9923) | 0.11 % | **yes** |
| swept volume, `r 0.1 -> 0.999` | Gauss quadrature 23.048969642534 | 23.048969642534 | 0.0e+00 | **yes** |
| swept mesh volume at 32 chords | Gauss quadrature | 23.587358955 | 2.6e-05 | **yes** |

Every reference above is either an exact closed form (error at round-off) or, for the two
discretisation statements, a convergence *order* whose measured value is within 0.5 % of the
theoretical one. The brute-force comparison carries an `O(h)` grid error of `1.4e-05`, which is
three orders below the 3 % acceptance line and is attributed to the grid, not to the kernel.

## 7. Limitations of this verification

* Nothing here re-derives the *arithmetic* of the arc integrals in a second symbolic form; R1/R2
  validate their outputs, not their internals. Both would have to be wrong in the same way to
  escape, which is unlikely but not excluded.
* Only 2-D sections and single-loop outlines are covered; the module's stated refusals
  (self-intersecting outlines, non-finite coordinates, zero-area outlines) were exercised at their
  boundary but not fuzzed.
* The sweep is verified as a *volume* and as a *manifold* (`V-E+F=2`, `t=2` per edge, no zero-area
  faces). Surface normals are checked only implicitly, through the sign of the divergence-theorem
  volume; an inverted normal would show as a negative volume and did not.
* `chamfer_solid` with a varying setback was not re-derived independently; only the zero-variation
  and small-variation cases pass through the same quadrature route.
* Runtime: the published example takes ~40 s single-core (the 1600² brute-force grid dominates);
  the test module re-runs it at 400² and takes under a second.

## 8. References

* `tensorcad/fillet.py` — kernel under test; `docs/fillet_chamfer.md` — the author's own report;
  `tests/test_fillet.py` — the author's 66 self-consistency cases.
* `docs/oss_inventory_cad/GAP-ANCHORS.md`, anchor `GAP-CAD-012` — the gap this module closes.
* Closed forms used: rounded-rectangle area/perimeter (elementary), centroidal second moment of a
  rectangle minus corner bites (derived in §2 R1), signed-corner area/perimeter law (module
  docstring, re-derived independently here), Gauss-Legendre 5-point nodes/weights (Abramowitz &
  Stegun table 25.4).