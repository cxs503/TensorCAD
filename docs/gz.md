# Large-angle stability: the righting-arm curve `GZ(phi)`

`tensorcad.stability` answers "how stable is it *at* the upright?" and returns one number,
`GM`, which is the slope of the righting-arm curve there. That is one point of one curve.
What a stability booklet, an intact-stability check and a damage-stability calculation
consume is the curve: the **maximum righting arm** and the heel at which it occurs, the
**range of positive stability**, and the **area under the arm curve** up to 30 and 40
degrees. Until 0.21 that curve was listed as not implemented, in this repository and in the
capability audit (see `docs/stability.md`, "No large-angle stability").

`tensorcad/gz.py` adds it. Every number on this page is written by
`python examples/gz/run_example.py` into `docs/gz/results.json` and
[`docs/gz/tables.md`](gz/tables.md) (generated, not typed), with the figure in
`docs/gz/plot.png`; the script exits non-zero if a claim stops holding.
`tests/test_gz.py` pins the same behaviour (39 tests).

## Method

Heel the body by `phi` about the longitudinal axis, keep the waterline horizontal, and hold
the displacement volume **fixed at its upright value** -- the definition of a stability
calculation:

```
y' =  y cos(phi) + z sin(phi)          positive phi puts starboard (+y) down
z' = -y sin(phi) + z cos(phi)

GZ(phi) = y_B(phi) - z_G sin(phi),     z_G = keel_z + kg
```

Three things are deliberate:

* **The waterline height is solved for, not assumed.** Constant displacement makes it an
  unknown; it is found by a bracketed root solve (Illinois regula falsi, ~15 volume
  evaluations) on "immersed volume below a plane", to a relative displacement residual
  `<= 1e-12`. Every station therefore carries its own convergence evidence
  (`GZCurve.max_volume_error`) instead of being smoothed in from two end states.
* **The immersed volume and centre of buoyancy come from the divergence theorem over the
  clipped surface**, with the apex of the tetrahedra *on the waterplane*:
  `V = (1/3) oint x.n dA`, `M_i = (1/2) oint x_i^2 n_i dA`. Every waterplane cap triangle
  passes through the apex, so it is degenerate: the cap contributes exactly zero volume and
  zero first moment and never has to be built. No cap polygon, no winding fix-up, no
  "close the loop" step -- and the result is exact for the triangulation given, so the only
  error is the hull's own discretisation.
* **A body whose displacement equals its volume is refused.** At exact full submergence the
  waterplane area vanishes and `GZ` is undefined; the module raises with a message naming
  freeboard, rather than returning a smooth meaningless number. The parametric hull of
  `tensorcad.hull` is exactly this case (its lid sits *on* the design waterline).

## Reference 1 -- the wall-sided barge (exact, over a *range* of angles)

A box barge `10 x 2 x 0.8 m` with `0.4 m` freeboard and `KG = 0.6 m` obeys the classical
wall-sided formula exactly while the waterline crosses only its two flat sides:

```
GZ(phi) = sin(phi) * (GM + BM tan^2(phi) / 2),   BM = B^2/(12 T),  KB = T/2
```

The measured arm reproduces it to `1e-13` relative at every station from 1 to 20 degrees --
round-off, not "good agreement" -- which pins the clipping, the constant-displacement solve
and the lever arm simultaneously. The validity limit is not hand-waved: deck-edge immersion
is at `tan(phi) = 2 f / B = 21.80 deg`, and the example runs the ladder past it to show the
measured arm falling *below* the formula (the formula has no standing there).

Beyond the formula's validity the same run reports the curve's engineering content: peak
arm, angle of peak, vanishing angle and the areas to 30 and 40 degrees -- see
[`docs/gz/tables.md`](gz/tables.md). The barge's vanishing angle comes out at exactly 90
degrees, and that is checkable by hand: heeled 90 degrees the section presents a full
rectangle to the waterline, so `y_B` lands on the vertical through `G`.

## Reference 2 -- the circular cylinder (exact at *every* angle)

A circular cylinder floating half-submerged is invariant under the heel rotation, so its
waterplane stays through the axis at every angle, its centre of buoyancy stays `4R/(3 pi)`
below the axis, and its metacentre sits exactly on the axis (`KB + BM = 4R/(3pi) + 4R/(3pi)
= R`). With `G` at `g` above the axis:

```
GZ(phi) = -g sin(phi)      for every phi, not just small ones
```

The measured arms follow it from 5 to 120 degrees, and the only error is the polygon's
approximation of the circle: it falls as `n^-2` (`n = 24 -> 96`), which is what a second-order
geometric error must do. This is the reference that tests the machinery at *large* angles,
where a small-angle expansion cannot hide a bug.

## Reference 3 -- the hull with freeboard (independent route to the same `GM`)

The parametric hull `10 x 2 x 0.8 m` is floated at `z = -0.11 m`, i.e. with `0.11 m` of
freeboard under its flat lid, `KG = 0.7 m`. Two numbers are computed by code paths that
share nothing but the mesh:

* `GM_waterplane = I_T / V - BG` from `tensorcad.stability.waterplane_moments` (a 2-D
  section integral) and `gz.immersed_moments` (a clipped volume integral);
* `GM_curve`, the slope of the arm curve at the first positive station.

They agree to `6e-4` relative on the default grid, and the arm at 10 degrees is grid
converged to a few percent between `41x13` and `81x25`. This is also the worked example of
the loading the freeze-the-upright-waterline convention cannot express: a body that has
freeboard at all.

## Notes, traps and limits

* **`compute_hydrostatics` / `initial_stability` assume the mesh lid *is* the waterline.**
  On a closed mesh they integrate the whole body: the same barge with `0.4 m` freeboard
  reports `V = 24 m3` instead of `16 m3` and `GM` follows it. That is not a defect of this
  module but the convention of that one, and it is why `GZ` carries its own clipping. For a
  body with freeboard, use `gz.immersed_moments` (clips) or the curve slope for `GM`.
* **A waterline through a vertex ring breaks the section code.** `waterplane_moments` raises
  "segments that do not close into a loop" when the plane hits a vertex ring exactly;
  `gz` does not, because Sutherland-Hodgman clipping treats an on-plane vertex as a proper
  member of the clipped polygon. Pick a waterline off the rings (the example uses `-0.11 m`
  for that reason) when calling the section code.
* **No free-surface correction**, no trim, no damage stability, no class rules: `GM` here is
  the solid value and the curve is a heel-only curve.
* **The deck is whatever the mesh closes with.** The packaged hull closes with a flat lid on
  the design waterline, i.e. zero freeboard, and is refused for exactly that reason; a real
  righting-arm curve for that hull family needs a freeboard/deck model first.
* **Trapezoid areas** carry `O(h^2)` error: on a 5 degree ladder the area to 20 degrees is
  3.3e-3 relative away from the closed-form integral. The sampled stations are part of the
  curve, so the spacing is auditable.

## Error verdict against the analytic references

| quantity | reference | measured error | under 3 %? |
|---|---|---|---|
| barge `GZ`, 1..20 deg | wall-sided closed form | `<= 1.4e-13` relative | yes (by 11 orders) |
| barge `GM` from the curve slope | `BM - BG` | `2.9e-4` relative, and it is the predicted `O(phi^2)` term | yes |
| barge area to 20 deg | `integral sin(phi)(GM + BM tan^2/2) dphi` | `3.3e-3` relative (trapezoid) | yes |
| cylinder `GZ`, 5..120 deg | `-g sin(phi)` | `3.2e-3` at 96 segments, second order | yes |
| hull `GM` | waterplane integral (independent code) | `6.2e-4` relative at `41x13` | yes |
| hull `GZ(10 deg)` | grid convergence `41x13 -> 81x25` | a few percent (mesh discretisation, no analytic reference) | yes, but see note |
| hull vanishing angle | no reference | — | declared: no closed form exists for this body |

The one quantity without an analytic reference is the hull's arm near its vanishing angle,
where the arm is the difference of two nearly equal numbers and is grid sensitive; the
example says so instead of quoting it as a result.