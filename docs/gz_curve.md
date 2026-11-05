# Righting arms at finite heel

`initial_stability` answers the upright question with a waterplane moment. `righting_arm_curve`
answers the one a designer actually plots: at each heel angle the body is rotated about its
longitudinal axis and the **waterplane height is solved** so that the immersed part of the closed
solid displaces the same volume as at the design waterline; the buoyancy centre of that clipped
solid gives the arm

    GZ(phi) = Y_G(phi) - Y_B(phi),      Y_G(phi) = -(keel_z + KG) sin(phi),

positive when restoring, with `G` on the centreline at `KG` above the keel. Nothing here assumes a
wall-sided body or a small angle.

Evidence: `examples/gz_curve/run_example.py` (exits non-zero if any claim below breaks),
`tests/test_gz_curve.py` (61 tests), `docs/gz_curve/results.json`, `docs/gz_curve/gz_curves.png`.

## How the integrals are built

The immersed volume and the buoyancy centre are divergence-theorem surface integrals over the
triangles of the *clipped* solid:

    V      = (1/3) * closed_integral  r . n dA
    int Y dV = (1/2) * closed_integral Y^2 n_Y dA
    int Z dV = (1/2) * closed_integral Z^2 n_Z dA

Each triangle is cut by the waterplane with Sutherland–Hodgman (winding preserved, because the
integrals are signed) and fan-triangulated. The per-triangle formulae are exact for polynomials of
degree two, so the result is the exact integral over the polyhedron, not a quadrature of it. The
waterplane cap needs no mesh: after shifting `z` so the plane is `z = 0`, the cap's normal is
vertical and its `z` is zero, so all three integrands vanish on it — the same trick that makes
`hydrostatics.py` cheap, applied to a cut surface instead of a pre-classified one.

The plane height is found by regula falsi with the Illinois anti-stalling correction (~30 exact
evaluations, bracket maintained throughout). Convergence tolerance is `1e-11` of the displacement,
reported per point as `displacement_volume`, which is how the "displacement is held fixed" invariant
is testable.

## The box: exact section, and the window a closed form lives in

`10 x 2 x (0.5 draft + 1.5 freeboard) m`, `KG = 0.45 m`, displacement `10 m^3`. The reference is an
**independent 2-D clipping**: the rectangle's four corners are rotated into the waterplane frame,
cut by a horizontal line, integrated with the shoelace rule, and the line is bisected until the area
is `B*T`. No formula or code is shared with the 3-D path.

| heel [deg] | `GZ` [m] | exact section [m] | `\|err\|` | wall-sided [m] | wall-sided rel. |
|---:|---:|---:|---:|---:|---:|
| 0 | 0.000000000 | 0.000000000 | 2.8e-17 | 0.000000 | — |
| 15 | 0.126976299 | 0.126976299 | 9.7e-17 | 0.126976 | 4.4e-16 |
| 30 | 0.282835427 | 0.282835427 | 1.5e-14 | 0.288889 | 2.1e-02 |
| 45 | 0.388908730 | 0.388908730 | 0.0 | 0.565685 | 3.1e-01 |
| 60 | 0.468478546 | 0.468478546 | 1.5e-14 | 1.270171 | 6.3e-01 |
| 90 | 0.550000000 | 0.550000000 | 0.0 | 1e31 (tan -> inf) | 1.0 |
| 120 | 0.484149399 | 0.484149399 | 1.5e-14 | 1.270171 | 6.2e-01 |
| 180 | 0.000000000 | 0.000000000 | 1.9e-17 | 0.000000 | 5.7e-01 |

* **Worst `|measured - exact section|` = 1.49e-14 m over 0..180 deg**, i.e. round-off: the
  integrator is exact on the polyhedron and every digit of the residual belongs to the geometry.
* **The wall-sided formula is a reference only inside its own window.** Two conditions must hold:
  both ends of the waterline stay on the vertical sides, `tan(phi) < 2T/B = 0.5` (`phi < 26.57
  deg`), and the waterline stays between keel and deck. Inside it, agreement is **4.4e-16 relative**
  — the classical `GZ = sin(phi)(BM + BM tan^2(phi)/2 - BG)` recovers from an independent route.
  One grid step outside it the formula is **2.1 % high** and by 45 deg it is 46 % high, because the
  leeward waterline has crossed the keel edge and the section is no longer wall-sided. A test
  asserts the departure as well as the agreement, so a formula applied everywhere would fail.
* **Deck edge immersion is at 65 deg, not at the naive `atan(2F/B) = 56.3 deg`**: the body sinks as
  it heels (fixed displacement, widening immersed section), so the waterline reaches the deck later.
  The flag is checked point by point against the reference's own waterplane height.
* `GZ_max = 0.554948 m` at 85 deg, `dynamic_stability_area = 1.099302 m*rad`. The arm only returns to
  zero at the inversion, so `vanishing_angle_deg` is `None` — reported as "no crossing inside the
  grid" rather than as 180.

## GM: the same number from two different integrals

| quantity | measured | reference | relative error |
|---|---:|---:|---:|
| box `GM` | 0.466666667 m | `T/2 + B^2/(12T) - KG` = 0.466666667 | 0.0 |
| box `dGZ/dphi(0)`, 5-deg secant | 0.468623 | `GM` | +0.42 % (`O(phi^2)`) |
| hull `GM` (`righting_arm_curve`) | 0.321622996 m | `initial_stability` 0.321622996 m | 0.0 |
| hull `dGZ/dphi(0)`, 1-deg secant | 0.318439 m/rad | `GM` | −0.99 % |
| Wigley `GM` | 0.128028 m | — | — |

The hull agreement is the interesting one: `BM` comes from the waterplane's **second moment** and
`KB` from the **clipped buoyancy centroid**, while the curve's slope comes from the heeled
displacement integration. They are the same arithmetic only if both paths are right, so the two
numbers must be equal — and they are, to the last digit.

## The hull families: wall-sided topsides

`tensorcad.hull` and `tensorcad.wigley` produce *wetted surfaces* closed at the design waterline. A
body capped there displaces its whole enclosed volume upright, so any heel would submerge it
completely and leave no waterplane at all: the curve refuses that case with an explicit error
instead of returning zeros. `freeboard=<h>` adds vertical, wall-sided topsides over the design
waterplane — the classical assumption, recorded in `basis` and never presented as measured geometry.

Parametric hull `10 x 2 x 0.8 m`, `KG = 0.6 m`, 0.4 m topsides, displacement `7.446554 m^3`:

| heel [deg] | `GZ` [m] | immersion depth [m] | `y_B` [m] | deck immersed |
|---:|---:|---:|---:|---|
| 0 | 0.000000000 | 0.800 | 0.000 | no |
| 10 | 0.051611168 | 0.785 | -0.017 | no |
| 20 | 0.097832346 | 0.775 | -0.029 | no |
| 30 | 0.137128275 | 0.756 | -0.037 | yes |
| 45 | 0.157530088 | 0.762 | -0.016 | yes |
| 60 | 0.149351577 | 0.863 | 0.024 | yes |
| 90 | 0.093730969 | 1.078 | 0.106 | yes |

`GZ_max = 0.157530 m` at 45 deg, `dynamic_stability_area = 0.183279 m*rad`, deck edge immersed from
25 deg onward, and no vanishing angle inside 0..90 deg.

## What this is not

No free-surface correction, no flooding, no damage case, no class rule, no deck-immersion water on
deck. Beyond the heel at which the deck edge immerses (`deck_immersed` per point) the arm belongs to
a closed watertight body and is **optimistic** — that flag, not the number alone, is the honest
range marker. A twin-hull or bulbous section needs its topsides in the mesh: the wall-sided
extension handles a single waterplane loop and refuses more rather than guessing.

## Reproducing

```bash
python examples/gz_curve/run_example.py --output docs/gz_curve
python -m pytest tests/test_gz_curve.py
```