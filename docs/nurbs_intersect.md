# NURBS x NURBS surface intersection and trim domains

`GAP-CAD-A1-NURBS-INTERSECT`: two *parametric* patches meeting each other, and the trim domain that meeting cuts on each of them. `surface_intersection` closed `GAP-CAD-002` for parametric-x-implicit; a B-rep kernel also needs patch-x-patch, because that is what a fillet, a hull/appendage seam and a plate boundary are.

Every number on this page is written by `examples/nurbs_intersect/run_example.py` from `results.json`; none is typed in by hand.

Environment: Python 3.11.6, Linux-6.6.52-1.7.lns23.loongarch64-loongarch64-with-glibc2.38, revision `f43e888`.

## 1. Method

Two patches meet in one dimension, so the problem collapses to a root find as soon as one patch is used as the parameter space. With `Q` the closest point of `S2` to `S1(u, v)` and `n2` the unit normal of `S2` there, the oriented distance

```
g(u, v) = (S1(u, v) - Q) . n2
```

vanishes exactly on the intersection and is a signed *length*, so a residual quoted in units of the patch scale means what it says. Marching squares traces its zero level set on the flat parameter rectangle -- branch and loop structure come from the cell topology, not from a proximity test -- and every crossing is then polished on the full system

```
F(u, v, s, t) = S1(u, v) - S2(s, t) = 0
```

with a minimum-norm Gauss-Newton step on the 3x4 Jacobian `[S1_u, S1_v, -S2_s, -S2_t]`. The march is the seed, the polish is the accuracy: that is why the vertices sit on both surfaces to round-off instead of to grid precision. `refine_levels` inserts parameter-space midpoints and re-projects them, which halves the chord step per level and moves the length down as `O(h**2)` without touching the branch structure.

Trim domains come from the same field: the kept side of a patch is the set where the oriented distance to the other patch is negative. That stays well defined even when the loop runs around a seam, which a polygon containment test does not. The loop is also returned as a polygon in parameter space, with `area_reliable` telling the caller whether its shoelace area means anything.

## 2. Cases

| case | patch A | patch B | analytic reference |
|---|---|---|---|
| sphere x plane | rational biquadratic sphere, R = 2 | plane patch z = 0.5, 6 x 6 m | `r = sqrt(R^2 - h^2)`, `2 pi r` |
| sphere x sphere (equal) | sphere R1 = 1 | sphere patch R2 = 1, centre z = 1.2 | `x = (d^2 + R1^2 - R2^2)/2d`, `r = sqrt(R1^2 - x^2)` |
| sphere x sphere (unequal) | sphere R1 = 1 | sphere patch R2 = 1.2, centre z = 1 | same closed form |
| cylinder x inclined plane | NURBS cylinder R = 1, h = 4, base z = -2 | plane tilted 30 deg | ellipse, `4 b E(m)` |
| cylinder x plane (closure study) | analytic periodic cylinder patch R = 1 | plane z = 0.3 | `2 pi R` |

## 3. Results

### 3.1 Circles: radius and arc length

| pair | r measured | r analytic | \|r - r_a\| / R | length | length analytic | rel. error | closed by |
|---|---|---|---|---|---|---|
| sphere R=2 x plane z=0.5 | 1.936491673104 | 1.936491673104 | 5.55e-16 | 12.167258038 | 12.167336028 | 6.41e-06 | seam |
| sphere R1=1 x sphere R2=1, centre z=1.2 | 0.800000000000 | 0.800000000000 | 4.44e-16 | 5.026516027 | 5.026548246 | 6.41e-06 | seam |
| sphere R1=1 x sphere R2=1.2, centre z=1 | 0.960000000000 | 0.960000000000 | 5.55e-16 | 6.031819232 | 6.031857895 | 6.41e-06 | seam |

Every vertex of the sphere-x-plane curve lies on the sphere to **1.11e-15 m** and on the plane to **6.66e-16 m**; the largest `|S1 - S2|` over that branch is **1.63e-15 m** against a patch scale of 8.500 m.

### 3.2 A tilted NURBS x NURBS pair

A cylinder cut by a plane tilted 30 degrees is an ellipse of semi-axes 1 and `1/cos 30`. Measured perimeter **6.777877122**, reference **6.777917613**, relative error **5.97e-06**. The reference integral itself is cross-checked: AGM `E(m) = 1.467462209339427` against a midpoint quadrature `1.467462209339420`.

### 3.3 Two closure paths

| patch | closes by | length | analytic 2 pi R | rel. error | vertices |
|---|---|---|---|---|---|
| analytic cylinder patch (periodic u) | loop | 6.283145881 | 6.283185307 | 6.27e-06 | 512 |
| NURBS cylinder (clamped seam) | seam | 6.283145033 | 6.283185307 | 6.41e-06 | 513 |

A *periodic* angle parameter closes the loop through grid topology (the seam is the same parameter twice). A *clamped* NURBS seam has two different parameters that coincide in space, so closure is detected on position and only the seam gap may be added as a chord -- refining that last segment would average two parameters a whole domain apart. Getting this wrong loses exactly one segment, which looks like a discretisation error and is not one.

### 3.4 Cross-check against a disjoint implementation

`surface_surface_intersection` marches the zero set of an *implicit* function on its own grid and never forms the 4-D Newton system used here.

| pair | S x S length | implicit-march length | rel. difference |
|---|---|---|---|
| sphere R=2 x plane z=0.5 | 12.167258038 | 12.167258035 | 2.64e-10 |
| sphere R=2 x sphere R=1.2 at z=1 | 5.729732211 | 5.729732210 | 2.76e-10 |
| NURBS cylinder x plane z=0.3 | 6.283145033 | 6.283145032 | 2.49e-10 |

Worst disagreement: **2.76e-10**.

### 3.5 Trim domains

Plane patch with unit parameter axes (so parameter area *is* surface area) trimmed by the sphere: the kept region is a disc of radius 0.458258, area `pi r^2 = 0.659734`.

| quadrature cells | measured area | reference | rel. error |
|---|---|---|---|
| 12 | 0.663194 | 0.659734 | 5.24e-03 |
| 16 | 0.659180 | 0.659734 | 8.41e-04 |
| 24 | 0.661024 | 0.659734 | 1.96e-03 |
| 32 | 0.659424 | 0.659734 | 4.71e-04 |

The masked integral is first order on the region boundary -- the boundary is a curve, not an area -- so the column above oscillates around the truth at the half-percent level instead of marching down monotonically. Splitting the boundary cells (4^3 sub-samples each) removes most of that bias; the loop polygon's shoelace area below is the better-conditioned number for the same region.

Loop polygon: shoelace area 0.659718, perimeter 2.879299, `area_reliable = True`, `seam_spanning = False`. `contains` is True at the disc centre and False at a patch corner.

Sphere (lat/long patch) trimmed by the same plane: the kept region starts at the south pole, so its *parameter* area is exactly `1.160861` and its surface area is `4 pi R^2 - 2 pi R (R - h) = 31.415927`.

| cells (u, v) | parameter area | reference | rel. error | surface area | reference | rel. error |
|---|---|---|---|---|---|---|
| 8, 16 | 1.156250 | 1.160861 | 3.97e-03 | 30.963473 | 31.415927 | 1.44e-02 |
| 16, 32 | 1.171875 | 1.160861 | 9.49e-03 | 31.483797 | 31.415927 | 2.16e-03 |
| 32, 64 | 1.171875 | 1.160861 | 9.49e-03 | 31.467590 | 31.415927 | 1.64e-03 |

Its loop has `area_reliable = False` and `seam_spanning = True`: a latitude line encloses no parameter area, and the module says so rather than printing a number.

### 3.6 OpenCASCADE parity

**WITHHELD** -- `neither the OCP (cadquery-ocp) nor the pythonocc bindings import in this environment; no OCCT numbers were generated and none are reported. The in-repo cross-check against the implicit-field marcher (see docs/nurbs_intersect.md section 3.4) is a different, weaker claim and is not a substitute for it.`

No OpenCASCADE number is reported here because none could be generated. The in-repo cross-check in section 3.4 is a *different and weaker* claim: it shows two independent implementations of the same mathematics agree, not that the result matches an industry kernel. The harness exists and is gated on the optional wheels, so the claim becomes measurable the moment they are installed; until then it is not made.

### 3.7 Convergence

| u cells | length | rel. error | order to next row |
|---|---|---|---|
| 32 | 6.081174137 | 4.10e-04 | 2.00 |
| 64 | 6.083044188 | 1.03e-04 | 2.00 |
| 128 | 6.083512034 | 2.56e-05 | 2.00 |
| 256 | 6.083629018 | 6.41e-06 | - |

| refine_levels | vertices | rel. error |
|---|---|---|
| 0 | 65 | 4.10e-04 |
| 1 | 129 | 1.03e-04 |
| 2 | 257 | 2.56e-05 |
| 3 | 513 | 6.41e-06 |

The chord argument predicts `O(h**2)`; the measured orders in the table are the check on that prediction, not a restatement of it.

## 4. Figures

*`curves_vs_analytic.png`* -- Returned intersection curves (dots) against the analytic circles (dashed) for sphere x plane and sphere x sphere. The vertices lie on both surfaces to round-off, so the two coincide at plot scale.

*`convergence.png`* -- Left: arc-length error against grid refinement and against midpoint refinement, with the O(h^2) chord prediction as a dotted reference. Right: trimmed-region area against the analytic disc and band areas.

*`trim_and_radius_error.png`* -- Left: the intersection curve written in the plane patch's own parameter space -- the trim loop bounding the kept disc. Right: per-vertex deviation from the analytic radius against the 1e-8 R acceptance bound.

## 5. Error verdict

| quantity | measured error | acceptance bound | within bound? | vs the 3% rule |
|---|---|---|---|---|
| circle radius (worst of three pairs) | 5.55e-16 R | 1e-8 R | yes | yes, by 5e+13x |
| arc length (worst of the analytic pairs) | 6.41e-06 relative | 3% (`GAP-CAD-002`) | yes | yes, by 5e+03x |
| cross-check vs implicit marcher | 2.76e-10 relative | 1e-4 (test bound) | yes | yes |
| trim area, plane disc | 4.71e-04 relative | 2% | yes | yes |
| trim area, sphere band | 1.64e-03 relative | 2% | yes | yes |
| OpenCASCADE parity | WITHHELD | 1e-4 length / 1e-4 scale | not claimed in this environment | n/a |

**Is every key quantity within 3%?** Yes: the worst relative error of any number on this page is 1.64e-03, five orders of magnitude inside the 3% rule. The tighter, task-specific bound of the anchor (circle radius < 1e-8 R) is met at 5.55e-16 R.

## 6. Limitations (measured, not asserted)

* **Tangential contact is not a curve.** A plane at z = R grazes the unit sphere in a point: the march returns 0 branches, the round-off cloud of crossings is collapsed and counted (1 collapsed), the field minimum on the grid is 0.0e+00 m against a patch scale of 5.74 m, and `tangential_candidate = True`. Disjoint patches also return no branches, but with the flag `False`, so the two situations are distinguishable.
* **Coincident patches raise.** Two identical spheres: `the two patches coincide over 100% of the sampled domain: the intersection is a surface, not a curve` -- the intersection there is a surface, and a cloud of branches would be a lie.
* **A finite patch truncates the curve.** A narrow plane patch across the same sphere gives 1 open branch (`closed = False`), length 3.738500, cut where the patch ends at x = 0.500 and x = 0.499. 1 crossing failed the residual test and was counted rather than returned.
* **A feature narrower than one cell is invisible.** No method can detect what it never sampled; `samples` is the caller's dial and the grid actually used is echoed in `diagnostics`.
* **Rank-deficient contact is damped, not solved exactly.** Where the two normals are parallel the 3x4 Jacobian loses rank; the step is Tikhonov-damped and the count is reported as `diagnostics['damped_steps']`. The curves around such a point are still correct, but the contact point itself is a point, not a branch.
* **No B-rep stitching.** This returns curves and trim regions; converting them into topological edges and faces of a solid is the next step, not this one.

## 7. Files

| path | what |
|---|---|
| `tensorcad/nurbs_intersect.py` | the oriented-distance field, the march, the 4-D polish, the trim domains |
| `tests/test_nurbs_intersect.py` | 18 tests: analytic anchors per pair, two closure paths, trim area, cross-check, honest failures, contract |
| `examples/nurbs_intersect/run_example.py` | the studies, the figures, this page |
| `tools/xcheck_occt_nurbs_intersect.py` | OpenCASCADE (pythonocc/OCP) parity harness; records `occt_available` and refuses to invent numbers when the wheels are absent |

## 8. Reproduction

```
python examples/nurbs_intersect/run_example.py --output artifacts/nurbs_intersect \
    --doc docs/nurbs_intersect.md
python -m pytest tests/test_nurbs_intersect.py -q
```

