# Righting arms at finite heel

`initial_stability` answers the upright question with a waterplane moment. `righting_arm_curve`
answers the one a designer actually plots: at each heel angle the body is rotated about its
longitudinal axis and the **waterplane height is solved** so that the immersed part of the closed
solid displaces the same volume as at the design waterline; the buoyancy centre of that clipped
solid gives the arm

    GZ(phi) = KN(phi) - KG sin(phi),

positive when restoring, with `G` on the centreline at `KG` above the keel, `KN` the buoyancy
centre's lever about the keel line. Nothing here assumes a wall-sided body or a small angle.

Evidence: `tests/test_gz_curve.py` (62 tests: every number below is asserted, not just printed),
`examples/gz_curve/run_example.py` (exits non-zero if one of its own claims breaks),
`docs/gz_curve/results.json`, `docs/gz_curve/gz_curves.png`.

## How the integrals are built — and where the code is

The heeled-displacement integration is **one** implementation, in `tensorcad/gz.py`:
Sutherland–Hodgman cut of every triangle against the plane (winding preserved, because the
integrals are signed), fan triangulation, and the divergence-theorem sums

    V        = (1/3) * closed_integral  r . n dA
    int Y dV = (1/2) * closed_integral Y^2 n_Y dA
    int Z dV = (1/2) * closed_integral Z^2 n_Z dA

which are exact for the polyhedron, not a quadrature of it. The waterplane cap needs no mesh:
after shifting `z` so the plane is `z = 0`, the cap's normal is vertical and its `z` is zero, so
all three integrands vanish on it. The plane height is found by Newton on `dV/dc = A(c)` with a
bisection fallback whenever the step leaves the bracket; the residual driven to is `1e-13` of
the displacement, and the *achieved* volume is reported per point as `displacement_volume` —
which is what makes the "displacement is held fixed" invariant testable instead of assumed.

`tensorcad.stability.righting_arm_curve` reuses those primitives rather than growing a second
copy of the numerics. What it adds is the stability-side reading: the arm spelled in `KG`, the
point-level flags a curve is judged by (`deck_immersed`, `immersion_depth`, `waterplane_height`),
the `GM` cross-check against `initial_stability`, and the validation a loading-condition input
needs. A test asserts the two readings produce the *same* arm, buoyancy centre and volume to
`1e-13`, so they cannot drift apart.

## The box: exact section, and the window a closed form lives in

`10 x 2 x (0.5 draft + 1.5 freeboard) m`, `KG = 0.45 m`, displacement 10.000000 m^3.
The reference is an **independent 2-D clipping**: the rectangle's four corners are rotated into
the waterplane frame, cut by a horizontal line, integrated with the shoelace rule, and the line
is bisected until the area is `B*T`. No formula or code is shared with the 3-D path.

| heel [deg] | `GZ` [m] | exact section [m] | `\|err\|` | wall-sided [m] | wall-sided rel. |
|---:|---:|---:|---:|---:|---:|
| 0 | 0.000000000 | 0.000000000 | 0.0e+00 | 0.000000 | — |
| 15 | 0.126976345 | 0.126976345 | 2.8e-17 | 0.126976 | 2.2e-16 |
| 30 | 0.282835427 | 0.282835427 | 0.0e+00 | 0.288889 | 2.1e-02 |
| 45 | 0.388908730 | 0.388908730 | 0.0e+00 | 0.565685 | 3.1e-01 |
| 60 | 0.468478546 | 0.468478546 | 2.2e-16 | 1.270171 | 6.3e-01 |
| 90 | 0.550000000 | 0.550000000 | 0.0e+00 | 1e31 (tan -> inf) | 1.0 |
| 120 | 0.484149399 | 0.484149399 | 1.1e-16 | 1.270171 | 6.2e-01 |
| 180 | 0.000000000 | 0.000000000 | 4.9e-17 | 0.000000 | 6.1e-01 |

* **Worst `|measured - exact section|` = 1.44e-15 m over 0..180 deg**, i.e.
  round-off: the integrator is exact on the polyhedron and every digit of the residual belongs to
  the geometry.
* **The wall-sided formula is a reference only inside its own window.** Two conditions must hold:
  both ends of the waterline stay on the vertical sides, `tan(phi) < 2T/B = 0.5` (`phi <
  26.57 deg`), and the waterline stays between keel and deck. Inside it,
  agreement is **1.2e-15 relative** — the classical
  `GZ = sin(phi)(BM + BM tan^2(phi)/2 - BG)` recovers from an independent route. One grid step
  outside it the formula is **2.1 % high** and by 45 deg it is
  31 % high, because the leeward waterline has crossed the keel edge and the section is no longer
  wall-sided. A test asserts the departure as well as the agreement, so a formula applied
  everywhere would fail.
* **Deck edge immersion is at 65 deg, not at the naive `atan(2F/B) = 56.3 deg`**:
  the body sinks as it heels (fixed displacement, widening immersed section), so the waterline
  reaches the deck later. The flag is checked point by point against the reference's own
  waterplane height.
* `GZ_max = 0.554948 m` at 85 deg, `dynamic_stability_area =
  1.099302 m*rad`. The arm only returns to zero at the inversion, so
  `vanishing_angle_deg` is `None` — reported as "no crossing inside the grid" rather than as 180.

## GM: the same number from two different integrals

| quantity | measured | reference | relative error |
|---|---:|---:|---:|
| box `GM` | 0.466666667 m | `T/2 + B^2/(12T) - KG` = 0.466666667 | 2.2e-16 |
| box `dGZ/dphi(0)`, 5-deg secant | 0.468623 | `GM` | +0.42 % (`O(phi^2)`) |
| hull `GM` (`righting_arm_curve`) | 0.321622996 m | `initial_stability` 0.321622996 m | 0.0e+00 |
| hull `dGZ/dphi(0)`, 1-deg secant | 0.318439 m/rad | `GM` | -0.99 % |
| Wigley `GM` | 0.128027805 m | — | — |

The hull agreement is the interesting one: `BM` comes from the waterplane's **second moment** and
`KB` from the **clipped buoyancy centroid**, while the curve's slope comes from the heeled
displacement integration. They are the same arithmetic only if both paths are right, so the two
numbers must be equal — and they are, to the last digit.

## The hull families: wall-sided topsides

`tensorcad.hull` and `tensorcad.wigley` produce *wetted surfaces* closed at the design waterline.
A body capped there displaces its whole enclosed volume upright, so any heel would submerge it
completely and leave no waterplane at all: the curve refuses that case with an explicit error
instead of returning zeros. `freeboard=<h>` adds vertical, wall-sided topsides over the design
waterplane — the classical assumption, recorded in `basis` and never presented as measured
geometry. `y_B` is the buoyancy centre in the body frame, positive toward the side that dips.

Parametric hull `10 x 2 x 0.8 m`, `KG = 0.6 m`, 0.4 m topsides, displacement
7.446554 m^3:

| heel [deg] | `GZ` [m] | immersion depth [m] | `y_B` [m] | deck immersed |
|---:|---:|---:|---:|---|
| 0 | -0.000000000 | 0.800 | -0.000 | no |
| 10 | 0.051611168 | 0.785 | 0.072 | no |
| 20 | 0.097832346 | 0.775 | 0.138 | no |
| 30 | 0.137128275 | 0.756 | 0.196 | yes |
| 45 | 0.157530088 | 0.762 | 0.249 | yes |
| 60 | 0.149351577 | 0.863 | 0.279 | yes |
| 90 | 0.093730969 | 1.078 | 0.301 | yes |

`GZ_max = 0.157530 m` at 45 deg, `dynamic_stability_area =
0.183279 m*rad`, deck edge immersed from 25 deg onward, and no vanishing
angle inside 0..90 deg.

## What this is not

No free-surface correction, no flooding, no damage case, no class rule, no water on deck. Beyond
the heel at which the deck edge immerses (`deck_immersed` per point) the arm belongs to a closed
watertight body and is **optimistic** — that flag, not the number alone, is the honest range
marker. A twin-hull or bulbous section needs its topsides in the mesh: the wall-sided extension
handles a single waterplane ring and refuses more rather than guessing.

## Reproducing

```bash
python examples/gz_curve/run_example.py --output docs/gz_curve
python -m pytest tests/test_gz_curve.py
```
