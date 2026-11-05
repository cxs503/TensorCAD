# Large-angle righting arms (`GZ`)

`tensorcad.stability` stops at the upright question: `GM = KB + BM - KG`, the slope of
the righting-arm curve at zero heel. `tensorcad.righting` is the curve — hold the
displaced volume fixed, heel the hull, and measure how far the line of buoyancy has
moved from the centre of gravity.

```python
from tensorcad import box_mesh, generate_wigley, righting_arm, righting_arm_curve

curve = righting_arm_curve(generate_wigley().mesh, (0, 10, 20, 30, 40), kg=0.35, waterline_z=-0.1)
curve.max_arm, curve.angle_of_max_arm_deg, curve.area(30.0)   # metres, degrees, m rad
```

## The definition, spelled out

In the body frame (`x` forward, `y` starboard, `z` up) a starboard heel of `phi` puts
world-up at `w = (0, -sin phi, cos phi)`, so the submerged region is the half-space
`w . p <= c`. The horizontal direction perpendicular to `w` (in the heeled frame) is
`e_h = (0, w_z, -w_y)`, and

```
GZ(phi) = (B - G) . e_h = (y_B - y_G) cos(phi) + (z_B - z_G) sin(phi)
```

with the righting moment `= displacement * GZ`. The sign convention makes a starboard
heel `phi > 0` and a positive `GZ` a *righting* arm; at small angles the expression
collapses to `GM sin(phi)`, which is the derivative check the test suite runs.

`kg` is measured **above the keel**, matching `initial_stability`; `waterline_z` sets
the reference displacement (the level waterline of the unheeled body), and
`displacement_volume=` overrides it for a displacement the caller already knows.

## How the integral is done

* **The waterline offset is solved, not assumed.** A `GZ` curve is drawn at constant
  displacement, so for every heel angle the plane has to rise or fall until the
  submerged volume matches. `solve_plane_offset` does that with a Newton iteration
  whose derivative is exact — `dV/dc` is the cut area — safeguarded by bisection on
  `[min(w.p), max(w.p)]`. Six iterations or fewer is typical; the residual volume
  error is reported on every point (`HeeledFloat.volume_residual`) rather than assumed
  away, and it lands at `3e-14` relative in the example.
* **The cut face is never triangulated.** Each triangle is clipped against the plane
  (Sutherland–Hodgman, winding preserved, so the fan keeps the outward orientation),
  and the tetrahedra are taken about a point *on* the waterplane. A tetrahedron with
  four coplanar corners has no volume, so the planar cut contributes exactly zero and
  the submerged part of the hull surface alone is a sufficient closed boundary. That
  removes fan-triangulating a non-convex section, and hole/parity bookkeeping with it.
* **The cut area is a flux, not a polygon.** `submerged_properties(...).cut_area` is
  the magnitude of the flux of that same surface, which equals the waterplane area by
  the divergence theorem — exact for the same reason. It is a genuine check on the
  clipping: at 5° heel the box's in-plane chord is `B / cos(5°)`, 0.38 % longer than
  the breadth, and the reported area matches it to 1e-12.

## Verification, and what it covers

`examples/righting/run_example.py` writes `docs/righting/{results.json,report.md,plot.png}`
and exits non-zero if any of its ten claims breaks. The references are independent of
the code under test:

| reference | what it pins down | measured |
|---|---|---|
| 2D section clip of a box (`analytic_box_righting_arm`) | the whole curve, every angle to 90°, including past deck-edge immersion | `2.2e-16 m` over 14 angles |
| wall-sided closed form `sin(phi)(GM + BM tan²(phi)/2)` | the curve from two upright quantities inside its validity window | `4e-17 m` below `26.565 deg`; `1.1e-2 m` departure just above it |
| `dGZ/dphi(0)` vs `GM = KB + BM - KG` | that the *heeled* integral agrees with the upright waterplane moments (`waterplane_moments`, Green's theorem) | `7e-9` relative on the box, `5e-5 … 2e-4` on the hull family (the residual is the finite-difference step and the faceting) |
| odd symmetry `GZ(-phi) = -GZ(phi)` | sign conventions and the offset solve on a non-prismatic hull | `2e-16 m` |
| mesh refinement on the Wigley family (41x13, 81x25, 161x41) | that the arm has converged, there being no closed form | `0.062484 -> 0.062520 -> 0.062528 m` at 20° (successive differences falling) |

The wall-sided comparison is deliberately reported on *both* sides of
`bilge_emergence_deg = atan(2T/B)`: reproducing the formula where it applies and
departing from it where one side of the section runs dry is what shows the check has
teeth. (`docs/righting/report.md` also carries the column past that limit, where
`tan²(phi)` diverges — that divergence is the formula's, not the code's.)

## What this is not

* **Solid hull, still water, free heel about a longitudinal axis.** No free-surface
  correction (`sum(rho i_x) / displacement` must be subtracted by whoever knows the
  tankage), no floodable compartments or damage cases, no trim, no wind, no lifting
  surfaces.
* **No stability criterion is evaluated.** `RightingArmCurve.area()` integrates the
  curve (`m rad`) so a caller can apply `IS Code 2008` style limits, but the library
  does not know KG, the tankage or the rule; nothing here is a compliance statement.
* **A mesh capped at the design waterline has no freeboard.** The `hull` and `wigley`
  families are closed with a lid exactly at `z = 0`, so floating them at `waterline_z
  = 0` leaves no reserve buoyancy and heel at that displacement is impossible
  (`solve_plane_offset` says so instead of returning a number). Float them lower —
  `waterline_z = -0.1` in the example — and the lid overhangs the waterline, which is
  the geometry of the *solid*, not of a ship with a bulwark.
* **`compute_hydrostatics` is the upright, meshed-at-the-waterline tool.** It
  integrates every triangle of the closed surface, so its displacement is the whole
  enclosed volume regardless of `waterline_z`. That is the right convention for a mesh
  capped at the waterline and the wrong one for a body with freeboard; when in doubt,
  use `upright_displacement(mesh, waterline_z)` here, which cuts the surface.