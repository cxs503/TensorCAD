# Large-angle stability: the righting-arm curve GZ(φ)

`docs/stability.md` closes with a warning about itself: *"Nor does it give a righting-arm
curve `GZ(phi)` for finite heel"*. That was the honest boundary of `tensorcad/stability`,
and this page is what sits on the other side of it.

    GZ(phi) = (B(phi) - G) . m,   m = (0, cos phi, sin phi),
    B(phi) = centroid of {p : p . n <= d},   n = (0, -sin phi, cos phi),
    V(d) = displacement volume                 (solved for d at every heel)

`tensorcad/gz.py` is the integration; `docs/gz/report.md`, `docs/gz/results.json` and the
three figures in `docs/gz/` are produced by

```bash
python examples/gz/run_example.py --output docs/gz
```

which **exits non-zero if one of its seven claims stops holding**, so the numbers below
cannot drift away from the code. `tests/test_gz.py` pins 45 cases of the same content
(`pytest tests/test_gz.py` → 45 passed in 21 s).

## How it is measured

* **The heeled body is a half-space clip, not a re-loft.** Each triangle is clipped with
  `p . n <= d` into a convex polygon of at most four vertices; the clipped surface plus its
  planar cap is closed, and the cap contributes zero volume and zero first moment when the
  tetrahedral fan's apex is placed *in the cap plane*. So `V` and `M` need no cap, no
  re-meshing and no special case at a deck edge.
* **The displacement is held, not the draft.** For a fixed loading `V` does not change with
  heel, so `d` is the root of `V(d) = V_0`, found by bisection. Bisection rather than Newton
  on purpose: `dV/dd` is the waterplane area, which collapses at a deck edge — exactly where
  a large-angle curve spends its time.
* **The lever direction is pinned by the small-angle limit.** With `y_B ~ BM tan phi` the
  formula above degenerates to `GZ ~ GM sin(phi)`, so the curve's slope at the origin *is*
  the `GM` of `tensorcad/stability`. That is a check, not a definition: `GM` there is
  `KB + BM - KG` built from a displaced volume and a waterplane's second moment, while the
  slope here differentiates a heeled displacement integration.
* **A second, independent route.** `cross_curve_righting_arm` gets the same number from the
  classical cross curves: cut each station into a `(y, z)` outline (`bonjean.section_loops`),
  clip it with the *inclined* waterline, integrate area and centroid along `x`. It shares
  nothing with the 3-D route but the mesh.

## What it measures

**A wall-sided box** (40 × 10 × 8 m, floating at 4 m, KG = 3 m; KB = 2 m, BM = 2.0833 m,
GM = 1.0833 m). The textbook result is `GZ = sin(phi) (GM + BM tan²(phi) / 2)` and a box's
faces are exactly representable, so the measured arm must equal it to round-off:

| heel [deg] | closed form [m] | measured [m] | rel. error |
|---:|---:|---:|---:|
| 5 | +0.095114 | +0.095114 | −2.3e−13 |
| 10 | +0.193743 | +0.193743 | −7.2e−14 |
| 20 | +0.417719 | +0.417719 | −1.1e−13 |
| 30 | +0.715278 | +0.715278 | −1.4e−13 |

Worst relative error over the valid window is **2.3e−15** (the table is read from
`docs/gz/results.json`; the four rows above are a subset). And the window is stated rather
than assumed: the closed form holds only while the waterline still meets the two vertical
sides, i.e. `tan(phi) < T/(B/2)` = **38.66°** for this box. Beyond it the comparison is not
made, because the reference itself stops applying — at 90° the wall-sided formula diverges
while the measured arm stays finite.

**The parametric hull** (12 × 2.4 × 0.9 m, floated at 0.31 m below the deck, KG = 0.55 m) on
three meshes. Two facts, both mesh-converging:

| resolution | faces | volume [m³] | GM [m] | slope of GZ at 0 [m/rad] | rel. error | max GZ [m] |
|---|---:|---:|---:|---:|---:|---:|
| 41×13 | 1 998 | 6.363932 | 0.372439 | 0.372447 | +2.3e−5 | 0.1634 |
| 81×25 | 7 838 | 6.425778 | 0.365372 | 0.365380 | +2.3e−5 | 0.1601 |
| 161×41 | 25 918 | 6.444410 | 0.363230 | 0.363239 | +2.3e−5 | 0.1591 |

* the slope at the origin reproduces `GM = KB + BM − KG` to **2.3e−5 relative** (worst of the
  three), computed from independent integrals;
* the 3-D route and the cross curves agree to **1.2e−3 relative** (worst, heels ≤ 60°), which
  is this page's honest error bar on `GZ` for a faceted hull;
* `GZ(φ) = −GZ(−φ)` to **8.6e−16 m**, and the heeled displacement stays on target to
  **1.1e−15 relative** at every angle sampled.

**A KG ladder.** `GZ_max` falls monotonically as `KG` rises while the loading is still stable
(0.1634 → 0.1134 → 0.0704 → 0.0334 m for KG = 0.55, 0.65, 0.75, 0.85 m), the `GM = 0`
crossing sits at `KM = 0.9154 m`, and past it the arms are negative at 5° (`GZ(5°) = −0.0545 m`
at KG = 1.55 m) with no range of stability to report. `area_under_curve(30)` / `(40)` is the
dynamical stability `∫ GZ dφ` used by statutory criteria — computed, but never read as a
pass/fail.

## Figures

* `docs/gz/box_vs_closed_form.png` — the box curve against the closed form, and the residual
  on a log axis (round-off, four decades below the 3% acceptance line).
* `docs/gz/hull_curves_and_ladder.png` — the hull's `GZ` from both routes, and `GZ_max`
  against `GM` for the KG ladder.
* `docs/gz/convergence.png` — both routes against mesh size, and the slope-vs-`GM` residual.

## Notes, limits and pitfalls found on the way

* **Not a rule check.** No free-surface correction (`Σ ρ i_x / displacement` must be
  subtracted by whoever knows the tankage), no wind heeling arm, no damage/permeability case,
  no `SOLAS`/class criterion. Longitudinal righting (trim) is not integrated.
* **`waterline_z` is a draft, and the displacement follows it.** The displaced volume is the
  part of the hull *below* the waterline, so it is taken from the same clip. This matters:
  `compute_hydrostatics` integrates the *whole* closed surface, so a hull mesh whose top face
  is well above the waterline would otherwise contribute its full volume as displacement. The
  package's own convention (the loft caps the hull at the design waterline) is exactly the
  case where the two agree, which is why `initial_stability` is only self-consistent at that
  waterline: for the hull floated at −0.31 m it reports `GM = 0.2867 m` against the
  metacentric `GM = 0.3632 m` from the clipped volume (`KB` 0.5406 vs 0.3547 m, `BM` 0.2961 vs
  0.5585 m). Not a bug in *its* contract, but a sharp edge worth knowing before comparing the
  two modules.
* **A waterline exactly on a vertex ring is degenerate** for the waterplane-moment cut
  (`waterplane_moments` raises "segments do not close"), while the half-space clip is happy
  with it. The hull's depth rings sit at multiples of 0.075 m, so `z = −0.30` fails and
  `z = −0.31` works. Pick a draft off the rings when comparing against the moment machinery.
* **Faceting is the only error left, and it is first order in the mesh but small.** Doubling
  the hull's stations/levels moves `GZ_max` by −2.0% (41×13 → 81×25) and then by −0.65%
  (81×25 → 161×41); the volume moves +1.0% then +0.3% (6.364 → 6.426 → 6.444 m³). So the
  coarse end is *not* converged: any `GZ` quoted from a 41×13 hull is a statement about the
  mesh, not the ship — 2% is inside the 3% acceptance band, but it is not zero.
* **`GZ_max` is reported at the sample grid's resolution**, and `range_of_stability` is a
  first-order interpolation between two sampled angles — a diagnostic, not a certificate.
* **Not covered:** heel angles beyond 90° (a capsized body needs the deck side solved too),
  and multi-hull cross sections whose station cuts are several disjoint loops nested through
  a tunnel (the station route sums them, the 3-D route does not care — untested).

## Where this sits

`hydrostatics` → `stability` (GM) → **`gz` (large-angle GZ, dynamical stability)**. The
natural next increments are the free-surface correction (needs tankage, an input this package
does not have yet), a wind heeling arm, and the `SOLAS` / `IS Code 2008` weather-criterion
check that consumes these curves; the ice and resistance lines stay in
`ship-performance-platform`, not here.