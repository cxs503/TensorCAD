# Large-angle stability: the righting-arm curve `GZ(φ)`

`tensorcad.stability` answers the small-angle question (`KB`, `BM`, `GM` from the waterplane
moments) and states plainly that `GM` is only the slope of the righting-arm curve at
`φ = 0`. This module supplies the curve itself, by integrating the **heeled displacement**
rather than a waterplane moment.

```
GZ(φ) = KN(φ) − KG·sin φ,        KN(φ) = y_B cos φ + z_B sin φ
```

`(y_B, z_B)` is the centre of buoyancy of the *heeled* displaced volume in the body frame,
measured from the **keel datum**, and `KG` is the centre of gravity above the same keel —
the same convention `tensorcad.stability.initial_stability` takes. `φ > 0` heels the `+y`
side down; `GZ > 0` is a righting moment.

## Method

1. Rotate the mesh about the `x` axis through the keel line.
2. Cut it with the horizontal plane `z = c`; integrate the part below with the divergence
   theorem, apex on the waterline plane so the closing cap contributes exactly zero.
3. Solve `c` such that the displaced volume equals the upright displacement — Newton's
   method with the exact derivative `dV/dc = A(c)`, the waterplane area, which the same
   pass accumulates as the projected area of the wetted triangles. Two to five iterations
   reach a relative volume residual below `1e-13`.

For a polyhedral mesh this is exact: the only error is double precision. The volume
residual is reported per curve point (`volume_error`, `max_volume_error`), because a curve
whose points float at different displacements is not a righting-arm curve.

Cost: one clip pass per Newton iteration per angle. A 41×13 hull mesh (≈2 000 triangles) at
19 heel angles runs in a few seconds on one core; nothing here is parallel or heavy.

## Verification against two closed forms

Both references are exact for the polygonal mesh, so the comparison is of the *integrator*,
not of two approximations.

### Box barge — the wall-sided trapezoid

The submerged section at heel `φ` is a trapezoid with parallel sides
`h1 = T − (B/2)tan φ` and `h2 = T + (B/2)tan φ`, whose mean height is the draft. Integrating
`∫z dA` and `∫y dA` over it:

| quantity | closed form |
|---|---|
| `KB(φ)` | `T/2 + (B²/24T)tan²φ = KB₀ + ½BM₀tan²φ` |
| `y_B(φ)` | `(B²/12T)tan φ = BM₀ tan φ` |
| `KN(φ)` | `sin φ (KB₀ + BM₀ + ½BM₀tan²φ)` |
| `GZ(φ)` | `sin φ (GM₀ + ½BM₀tan²φ)` |

with `KB₀ = T/2`, `BM₀ = B²/(12T)`, `GM₀ = KB₀ + BM₀ − KG`. This is the classical
wall-sided formula, *derived* here rather than quoted.

Validity ends at bilge emergence on the rising side (`tan φ = 2T/B`) or deck-edge immersion
on the dipping side (`tan φ = 2(D − T)/B`); `box_wall_sided_limit_deg` returns the smaller.

Measured (10 × 2 × 2 m barge, draft 0.8 m, KG 0.5 m, limit 38.66°):

| heel [deg] | GZ integrated [m] | wall-sided [m] | abs. error [m] |
|---|---|---|---|
| 5 | 0.027738300 | 0.027738300 | 2.1e-17 |
| 20 | 0.117745731 | 0.117745731 | 5.6e-17 |
| 30 | 0.193055556 | 0.193055556 | 5.6e-17 |
| 35 | 0.240219858 | 0.240219858 | 2.8e-17 |

worst over 0–35°: **6.9e-17 m**. Past the limit the closed form must *stop* holding, and it
does — 4.0e-04 m at 45°, 1.25 m at 70° — which is what makes the agreement above meaningful
rather than vacuous.

### Log — circular section, half submerged

Any line through the centre of a centrally symmetric section splits it into equal areas, so
the waterline stays through the section centre at every heel. The submerged centroid then
sits a distance `d` below that centre along the waterline normal, and rotating it back into
the body frame gives `KN = R sin φ + s`, with `R` the centre's height above the keel and `s`
the centroid's offset along the cut line. For a circle `d` is constant and `s` vanishes:

```
GZ(φ) = (R − KG)·sin φ        for 0 ≤ φ ≤ 90°
```

i.e. a round bilge is stable only with `G` below the section centre — the same statement as
`GM = KB + BM − KG = R − KG`, because `KB = R − 4R/(3π)` and `BM = 4R/(3π)`
(`I_T = L(2R)³/12`, `V = πR²L/2`) sum to exactly `R`.

Measured on a 64-sided polygon (`R = 1 m`, `KG = 0.5 m`) at 33 heel angles from 0 to 90°,
each a symmetry axis of the polygon: worst `|GZ − (R−KG)sin φ| = 1.1e-16 m`.

## Cross-check against the waterplane route

KB and the slope of the curve at zero must reproduce `initial_stability`, which computes
them by Green's theorem on the waterplane loops — a different integral of different
quantities:

| quantity | volume-integral route | waterplane route |
|---|---|---|
| `KB` (41×13 hull) | 0.483005125 m | 0.483005118 m |
| `dGZ/dφ` at 0 (hull) | 0.421622996 m/rad | `GM` = 0.421622996 m (1.8e-07 %) |

On a polygonal hull, `GZ'(0)` has a curvature *kink*: at `φ = 0` the waterline coincides with
the hull's waterline vertex ring, so a central difference converges only linearly in the
step. `examples/gz/run_example.py` removes that first-order term (Richardson, step ratio 10)
before comparing; with it the agreement is at the `1e-9` level.

## Free surface

`free_surface_correction(tank_inertias, displacement_volume)` returns `Σi_x / V`, the loss of
*metacentric height* a slack tank costs, and `rectangular_tank_inertia(l, b) = l·b³/12`
supplies `i_x` for a rectangular tank. A loss of `GM` is a raise of `KG`, so the curve loses
`Σi_x/V · sin(φ)` — not `Σi_x/V` at every ordinate: a flat subtraction puts
`GZ(0) = −Σi_x/V` on an upright symmetric ship and overstates the low-angle areas. A
2.0 × 0.5 m tank in the 7.45 m³ hull example costs 2.80 mm of `GM`; the corrected `A(0-30)`
loses 3.7e-4 m rad where the flat model took 1.5e-3 m rad, 3.9× too much.

`free_surface_model="flat"` reproduces the historical behaviour for an old curve.
`free_surface_model="exact"` (with `slack_tanks=[SlackTank(...)]`) uses the finite-angle shift
of [`tensorcad.free_surface`](free_surface.md) instead of the first-order one; that page
measures the difference — 16.7 % more loss at 30° for a 2.0 × 0.5 × 1.0 m tank half full,
because the free surface widens as `b/cos(φ)`.

## `with_freeboard`, and why a generated hull cannot be heeled

`generate_hull` and `generate_wigley` close the body at the **design waterline**. Heel it and
the waterline immediately leaves the modelled surface on the dipping side, so the displaced
volume silently stops meaning anything. `with_freeboard(mesh, freeboard)` removes the
waterline cap, extrudes the boundary ring upwards with vertical walls and closes it with a
fan; the upright displacement is unchanged (asserted to round-off in the tests) and the body
stays watertight, outward-oriented and Euler-2 (`V − E + F = 2`).

## Limitations

- **Sealed body.** `GZ` is for an intact, sealed hull. Deck-edge immersion
  (`deck_edge_immersion_deg`, `atan(freeboard / B/2)`) does not flood anything, so the curve
  keeps rising past it. Real ships lose righting arm there; that is a flooding criterion, not
  a geometry question, and is not modelled.
- **Symmetry assumed for the sign.** `GZ(−φ) = −GZ(φ)` holds only for a symmetric hull (the
  tests assert it for one and it is the documented convention); the integrator itself handles
  asymmetric waterplanes because `y_B` is measured, not assumed zero.
- **No tanks, no damage.** `KG` is the only loading input. No compartment contents, no
  flooding, no grain heeling, no statutory maxima, no class rules.
- **No weather criterion here.** This module computes no wind heeling arm and no IMO/IS-Code
  check; the areas `area_under_curve` returns are the input to such a check. The check itself
  lives in [`docs/intact_stability.md`](intact_stability.md) (`tensorcad.intact_stability`,
  0.23.0: IS Code 2008 part A 2.2.1-2.2.6 and 2.3), which reads these curves and does not
  touch the hydrostatics.
- **Mesh resolution is the accuracy.** The integrator is exact for the polyhedron it is
  given; the hull's chord-wise waterplane (a 41-station mesh) is what limits `KN`, exactly as
  it limits `KB`/`BM`.

## API

```python
from tensorcad import righting_arm_curve, with_freeboard, heeled_hydrostatics

banded = with_freeboard(hull.mesh, 1.0)               # close the hull above the waterline
curve = righting_arm_curve(banded, kg=0.5, angles_deg=range(0, 91, 5))
curve.max_gz, curve.angle_of_max_gz_deg, curve.range_of_stability_deg
curve.area_under_curve(30.0)                          # dynamical stability, m·rad
curve.to_dict()                                       # JSON-ready, per-point volume residuals
point = heeled_hydrostatics(banded, 30.0)             # one angle: V, KB, KN, A_w, residuals
```

`examples/gz/run_example.py` exercises all of the above, writes `results.json`,
`report.md` and `plot.png`, and exits non-zero if a claim breaks.