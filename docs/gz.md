# Large-angle stability: the `GZ` curve

`docs/stability.md` ends with "no large-angle stability … `GZ` needs a heeled cut". This is
that cut. `tensorcad.gz` reports the righting arm of a loading condition,

    GZ(phi) = KN(phi) - KG * sin(phi)

with `KN` the transverse distance from the centreline to the vertical through the centre of
buoyancy of the **heeled** hull, and `KG` the loading input, measured above the keel, exactly
as in `initial_stability`. Positive heel means starboard down and positive `GZ` means righting;
both conventions are asserted in `tests/test_gz.py` rather than assumed.

Everything on this page comes from `python examples/gz/run_example.py --output docs/gz`,
which writes `results.json`, `report.md` and `figures/gz-curve.png` and **exits non-zero if
one of its claims stops holding**. The numbers below are read from that JSON.

## What has to be built for it

Two things `stability.py` deliberately does not do:

* **A heeled waterplane.** The model is rotated about a longitudinal axis through the keel
  centreline and cut with the horizontal plane `z = waterline_z` in the rotated frame.
* **A sinkage solve.** Displacement is conserved as the hull goes over, so the waterline is
  moved until the submerged volume equals the upright displaced volume: Newton on
  `V(z) - V_target` with the waterplane area as `dV/dz`, safeguarded by bisection. A column
  of the results tables is that solved waterline, so the conservation is visible rather than
  asserted.

The submerged volume, its centroid and the waterplane area are integrated by clipping the
surface against the half-space, not by sampling:

* triangles wholly below the waterline are kept, wholly above ones dropped, crossing ones
  clipped (Sutherland-Hodgman) and fan-triangulated;
* the **cap is not triangulated and the section is not chained into loops.** Each crossing
  triangle's outward normal gives the section direction `z_hat x n`, which walks the material
  counter-clockwise seen from above, so the oriented segments' signed sum is the area with
  holes already subtracted. The cap's contribution to the volume integrals is then closed
  form: `2*wl*A` to the six-volume sum, `6*wl*x_bar*A`, `6*wl*y_bar*A`, `6*wl^2*A` to the
  first moments. This is exact for a concave section, and it avoids a failure mode that took
  a hand-trace to find: the loop-chaining in `stability.py` keys a crossing point on the
  mesh's own edge identity, and on a facetted cylinder two faces meeting at a shallow angle
  key the *same point* differently -- the chain broke into three loops and the area came out
  at exactly half the truth. The segment orientation cannot break that way.
* the sinkage solve's **initial guess** used to be the bracket midpoint, and the transom
  family merged in alongside this work emits a ring at half draft -- exactly that midpoint.
  A trial waterline on a ring has no area, so the whole curve raised. The start now steps
  off mesh rings by `1e-4` of the span (`gz._start_waterline`); displacement is continuous
  in the waterline height, so the converged height and its conservation are unchanged, and
  `submerged_solid` still refuses an on-edge waterline that is *requested*.
* a plane that runs **along** mesh edges has no section at all. The closure guard counts how
  often each section point appears; an odd count raises `GZError` instead of reporting half an
  area. That is a real case: a facetted circle with a vertex exactly at half draft, cut at
  that draft, puts the whole section on edges.

## The two exact references

| reference | why it is exact | measured |
|---|---|---|
| wall-sided barge (`L=30, B=6, T=2, D=10`, KG=1 m) | while the waterline still cuts both side walls (`tan(phi) < 2T/B = 0.667`, i.e. `phi < 33.7 deg`) the arm is `sin(phi)(GM + BM tan^2(phi)/2)` with `KB=T/2`, `BM=B^2/(12T)` | worst deviation **3.3e-16 m** over 7 angles |
| round bilge (`R=3, L=20, 360 facets`, KG=1.5 m) | a chord through the centre is the only area-bisecting chord of a circle, so B is on the vertical through the centre at every angle: `KN = R sin(phi)` and the floatation is `z = R cos(phi)` | worst `GZ` deviation **9.9e-5 m**, floatation **5.7e-5 m** -- the facetting, not the cut |

The wall-sided case is a polyhedron and the cut is exact, so the agreement is at rounding
level. The round bilge exercises the general path (curved billows, an oblique cut at every
angle) and is asserted to 2e-3, which is where 360 facets put it.

Cross-module consistency: on a hull capped at its design waterline, `gz_curve(...).gm_upright`
equals `initial_stability(...).gm_transverse` to the last bit (`tests/test_gz.py`), and for
the barge `dKN/dphi` at zero heel reproduces `KM = KB + BM` to 2e-3.

## IMO screening

`imo_intact_stability(curve)` screens the six numerical criteria of the 2008 IS Code Part A
2.2 against the measured curve: area 0-30 deg >= 0.055 m·rad, 0-40 deg >= 0.090, 30-40 deg
>= 0.030, `GZ(30) >= 0.20 m`, `GZ_max >= 0.20 m` at not less than 25 deg, and `GM >= 0.15 m`.
The areas are integrated on the piecewise-linear curve (trapezoid on the segments is the
integral, not an approximation of it) and the curve must span 0-40 deg, because an area
extrapolated past the last sampled angle would be a number nobody measured. On the barge the
screen passes with margin (`docs/gz/report.md`); raising `KG` to 2.45 m (GM = 0.05 m) makes it
fail on `gm_upright` and on the area terms, which is the test that the screen is not a
rubber stamp.

This is a **screening**, not a stability booklet: no weather criterion, no free-surface or
grain correction, no damaged/floodable-length case, no class approval.

## What it does not do

* **Constant trim.** The waterline is kept parallel to the still-water line: this is the
  classical cross curve (`KN`), not a free-to-trim solve. A hull that trims as it heels will
  differ, and the difference is second order in the trim angle.
* **Free-surface corrections, grain and heeling-moment allowances, damage stability, class
  rules** -- all still open, as `README.md` says.
* **`generate_hull` and `generate_wigley` cap the model at the design waterline, i.e. they
  have no freeboard.** Heel the model and the deck edge is under water from the first degree,
  so a `GZ` curve on such a mesh is the *no-flooding idealisation* and is only quoted over
  0-20 deg in the example. Large-angle work needs a model with freeboard -- that is why the
  large-angle vehicle here is a barge, and why "a hull family with freeboard" is the next
  item on the list rather than a detail.
* The section is read off the *surface*: a mesh that is not watertight has no volume, and
  `hydrostatics` refuses it first.

## Files

| artefact | what it is |
|---|---|
| `tensorcad/gz.py` | `submerged_solid`, `solve_waterline`, `heel_mesh`, `gz_curve`, `imo_intact_stability` |
| `tests/test_gz.py` | 31 cases: two exact references, IMO pass/fail, sign conventions, error paths |
| `examples/gz/run_example.py` | produces the tables, the figure and the claim ledger above |
| `docs/gz/{results.json,report.md,figures/gz-curve.png}` | the artefacts of that run |