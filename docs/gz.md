# Large-angle stability: the righting-arm curve

`tensorcad.stability` measures the *slope* of the righting-arm curve at the upright and
says so: `GM` is that slope, and stability at a finite heel is a different integration.
`tensorcad.gz` is that integration. A hull heeled to starboard by `phi` floats on the
plane whose up-normal in the ship's axes is exactly `n = (0, -sin phi, cos phi)`, so the
submerged set is the half-space `p . n <= d`; the displacement does not change with heel,
so `d` is not free — it is the root of `V(d) = displacement volume`, solved by bisection
on `V(d)`, which is monotone for a closed body. With `B` the centroid of the submerged
volume and `G = (0, 0, keel + KG)`,

    GZ(phi) = (B - G) . (0, cos phi, sin phi)

Two independent routes to the same number are shipped and compared:

* **the 3-D route** (`righting_arm`, `righting_arm_curve`) — clip every triangle of the
  closed surface with the inclined plane, integrate volume and first moments as a fan of
  tetrahedra whose apex sits *in* the cutting plane (the cap contributes nothing and is
  never built), and bisect the offset on volume. Exact for the faceted body: the only
  discretisation in the chain is the mesh's own faceting.
* **the cross-curve route** (`cross_curve_righting_arm`, `cross_curve_righting_arm_curve`)
  — the classical construction: cut the hull at each station into `(y, z)` outlines
  (`bonjean.section_loops`), clip each outline with the *inclined* waterline, integrate
  area and centroid along `x`. It shares nothing with the 3-D route except the mesh, so
  agreement between them is evidence about the integration rather than about one
  implementation's habits.

## What it measures

**A wall-sided box** (40 x 10 x 8 m, floating at 4 m, KG = 3 m; KB = 2.0000 m,
BM = 2.0833 m, GM = 1.0833 m). The textbook result is
`GZ = sin(phi) (GM + BM tan^2(phi) / 2)` and a box's faces are exactly representable, so
the measured arm must equal it to round-off. Worst relative error over the valid window
(heels below the wall-sided limit 38.66 deg) is
**2.33e-15**, and the heeled displacement stays on target to
**1.14e-15** relative:

| heel [deg] | closed form [m] | measured [m] | rel. error | volume rel. error |
|---:|---:|---:|---:|---:|
| 5 | +0.095114 | +0.095114 | -2.33e-15 | 1.14e-15 |
| 10 | +0.193743 | +0.193743 | -7.16e-16 | 1.14e-15 |
| 20 | +0.417719 | +0.417719 | -1.06e-15 | 9.95e-16 |
| 30 | +0.715278 | +0.715278 | -1.40e-15 | 9.95e-16 |

Beyond 38.66 deg the waterline leaves the vertical sides and the closed form
stops applying (it diverges at 90 deg while the measured arm stays finite), so no
comparison is claimed there — the full grid is in `docs/gz/report.md`.

**The parametric hull** (12 x 2.4 x 0.9 m, floated at -0.31 m, KG = 0.55 m)
on three meshes. The slope of `GZ` at the origin is checked against
`GM = KB + BM - KG` built from a *different* pair of integrals (`KB` from the displaced
volume, `BM` from the waterplane's second moment):

| resolution | faces | volume [m3] | KB [m] | BM [m] | GM [m] | dGZ/dphi at 0 [m/rad] | rel. error | max GZ [m] |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 41x13 | 1998 | 6.363932 | 0.357971 | 0.564468 | 0.372439 | 0.372447 | 2.28e-05 | 0.1634 |
| 81x25 | 7838 | 6.425778 | 0.355515 | 0.559858 | 0.365372 | 0.365380 | 2.29e-05 | 0.1601 |
| 161x41 | 25918 | 6.444410 | 0.354732 | 0.558499 | 0.363230 | 0.363239 | 2.28e-05 | 0.1591 |

* slope at the origin reproduces `GM` to **2.29e-05** relative (worst of the three);
* the 3-D route and the cross curves agree to **1.22e-03** relative (worst,
  heels <= 60 deg) — this page's honest error bar on `GZ` for a faceted hull;
* `GZ(phi) = -GZ(-phi)` to **8.60e-16 m**, and the heeled displacement stays on target to
  machine precision at every angle sampled.

**A KG ladder.** `GZ_max` falls monotonically as `KG` rises while the loading is still
stable (0.1634, 0.1134, 0.0704, 0.0334 m for KG = 0.55, 0.65, 0.75, 0.85 m),
the `GM = 0` crossing sits at `KM = 0.9154 m`, and past it the arms
are negative at 5 deg with no range of stability to report.
`area_under_curve(30)` / `(40)` is the dynamical stability `int GZ dphi` used by statutory
criteria — computed and reported, never read as a pass/fail:

| KG [m] | GM [m] | GZ_max [m] | at [deg] | GZ(5 deg) [m] | A(0-30) [m rad] | A(0-40) [m rad] | range [deg] |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.55 | +0.3654 | 0.1634 | 30 | +0.0327 | 0.05015 | 0.07838 | — |
| 0.65 | +0.2654 | 0.1134 | 30 | +0.0240 | 0.03677 | 0.05502 | — |
| 0.75 | +0.1654 | 0.0704 | 25 | +0.0152 | 0.02338 | 0.03166 | — |
| 0.85 | +0.0654 | 0.0334 | 20 | +0.0065 | 0.00999 | 0.00830 | — |
| 1.05 | -0.1346 | 0.0000 | 0 | -0.0109 | -0.01679 | -0.03842 | — |
| 1.55 | -0.6346 | 0.0000 | 0 | -0.0545 | -0.08373 | -0.15523 | — |

For the two unstable loadings `GZ_max` is the upright value (0 deg), which is the correct
maximum of the curve and the reason the column pair "GZ_max / at" must be read together:
a maximum at 0 deg with a negative `GZ(5 deg)` is a body that does not float upright.

## Error judgement

| quantity | reference | measured | error | < 3%? |
|---|---|---|---:|---|
| box `GZ`, heels < 38.7 deg | wall-sided closed form | 3-D clip | 2.33e-15 | yes |
| box heeled displacement | upright displacement | 3-D clip | 1.14e-15 | yes |
| hull `dGZ/dphi` at 0 | `KB + BM - KG` | finite difference on the curve | 2.29e-05 | yes |
| hull `GZ`, heels <= 60 deg | cross-curve route | 3-D clip | 1.22e-03 | yes |
| hull `GZ` symmetry | `-GZ(-phi)` | `GZ(phi)` | 8.60e-16 m | yes |

The two routes agree to ~1e-3, not to machine precision, and that gap is *mesh faceting*,
not quadrature: the station cuts are polygons on a faceted surface, the 3-D clip is exact
for the same faceted surface, and their difference is the second-order effect of where the
two discretisations put the waterline on a chord. It is quoted rather than hidden.

## Figures

* `docs/gz/box_vs_closed_form.png` — the box curve against the closed form, and the
  residual on a log axis (round-off, many decades below the 3% acceptance line).
* `docs/gz/hull_curves_and_ladder.png` — the hull's `GZ` from both routes, and `GZ_max`
  against `GM` for the KG ladder.
* `docs/gz/convergence.png` — both routes against mesh size, and the slope-vs-`GM` residual.

## Notes, limits and pitfalls found on the way

* **Not a rule check.** No free-surface correction (`sum(rho i_x) / displacement` must be
  subtracted by whoever knows the tankage), no wind heeling arm, no damage/permeability
  case, no `SOLAS` / `IS Code` criterion. Longitudinal righting (trim) is not integrated.
* **`waterline_z` is a draft, and the displacement follows it.** The displaced volume is
  the part of the hull *below* the waterline, taken from the same clip; `compute_hydrostatics`
  instead integrates the *whole* closed surface. For the package's own convention (the loft
  caps the hull at the design waterline) the two agree; for a hull floated elsewhere they do
  not, which is why `initial_stability` and this module are only directly comparable at the
  design waterline.
* **A hull lofted with its deck at the waterline has no freeboard.** `generate_hull` caps
  the body at the design waterline, so `righting_arm_curve(mesh, kg=..)` with the default
  `waterline_z = 0` floats the hull at its deck edge: heeling submerges nothing and the arms
  come out negative (the hull is a closed box whose top face is the waterline). Float it at a
  real draft (`waterline_z` below the deck) before reading a curve; the example here uses
  `waterline_z = -0.31` on a hull whose deck sits at `z = 0`.
* **A waterline exactly on a vertex ring is degenerate** for the waterplane-moment cut
  (`waterplane_moments` raises "segments do not close") while the half-space clip is happy
  with it. Pick a draft off the rings when comparing against the moment machinery.
* **Faceting is the only error left, and it is first order in the mesh but small.** Doubling
  the hull's stations/levels moves `GZ_max` by about -2.0% (41x13 -> 81x25) and then about
  -0.6% (81x25 -> 161x41). A `GZ` quoted from a 41x13 hull is a statement about the mesh,
  not the ship — inside the 3% band, but not zero.
* **`GZ_max` is reported at the sample grid's resolution**, and `range_of_stability` is a
  first-order interpolation between two sampled angles — a diagnostic, not a certificate.
* **Bisection, not Newton.** `dV/dd` is the waterplane area, which collapses to zero at a
  deck edge — exactly where a large-angle curve spends its time; bisection cannot be trapped
  by that flat spot.
* **Not covered:** heel angles beyond 90 deg (a capsized body needs the deck side solved
  too), and multi-hull cross sections whose station cuts are several disjoint loops nested
  through a tunnel (the station route sums them, the 3-D route does not care — untested).

## API

`DEFAULT_HEEL_ANGLES`, `GZError`, `HeelPoint`, `RightingArmCurve`, `righting_arm`,
`righting_arm_curve`, `cross_curve_righting_arm`, `cross_curve_righting_arm_curve`,
`submerged_properties`, `solve_waterline`, `keel_height`, `station_offsets`,
`wall_sided_righting_arm` — all re-exported from `tensorcad`.

## Where this sits

`hydrostatics` -> `stability` (GM) -> **`gz` (large-angle GZ, dynamical stability)**.
The natural next increments are the free-surface correction (needs tankage, an input this
package does not have yet), a wind heeling arm, and the `SOLAS` / `IS Code 2008`
weather-criterion check that consumes these curves; the ice and resistance lines stay in
`ship-performance-platform`, not here.

Reproduce every number on this page with
`python examples/gz/run_example.py --output docs/gz`; the mechanical claim table is
`docs/gz/report.md`, the raw record `docs/gz/results.json`.
