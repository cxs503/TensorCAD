# IS Code 2008 intact-stability criteria on a righting-arm curve

`tensorcad.statutory` answers the one question the rest of the stability stack deliberately
does not: **does this GZ curve comply with the general intact-stability criteria?** The two
curve routes (`tensorcad.gz`, `tensorcad.righting_arm`) and the small-angle route
(`tensorcad.stability`) each say in their own docstrings that a statutory check is out of
their scope. This page documents the module that closes that gap, the rule it applies, and
the evidence that the numbers it compares are not the sloppy part of the comparison.

## 1. The rule

IMO Resolution MSC.267(85), *International Code on Intact Stability, 2008*, Part A,
chapter 2, section 2.2 (general criteria). Thresholds are held as data in
`INTACT_STABILITY_CRITERIA`, so a caller can print the rule next to the numbers:

| code | requirement | threshold |
|---|---|---|
| 2.2.1a | area under GZ, 0 to 30 deg | 0.055 m rad |
| 2.2.1b | area under GZ, 0 to 40 deg or to the flooding angle `theta_f` if that is less | 0.090 m rad |
| 2.2.1c | area under GZ, 30 deg to 40 deg (or to `theta_f`) | 0.030 m rad |
| 2.2.2 | GZ at 30 deg | 0.20 m |
| 2.2.3 | heel at which GZ reaches its maximum | 25 deg |
| 2.2.4 | initial metacentric height `GM0` | 0.15 m |

`theta_f` is an input, not an assumption: passing it in is what makes 2.2.1b/2.2.1c *the
rule* rather than the 40-degree special case. Without it the assessment says so in the check
note (`flooding angle not supplied`).

```python
from tensorcad import check_intact_stability, generate_hull, HullParameters
from tensorcad.gz import righting_arm_curve

hull = generate_hull(HullParameters(length=10.0, beam=2.0, draft=0.8))
curve = righting_arm_curve(hull.mesh, kg=0.9, angles_deg=(0, 10, 20, 30, 40), waterline_z=0.8)
assessment = check_intact_stability(curve, flooding_angle_deg=45.0)
print(assessment.verdict(), assessment.to_dict()["checks"])
```

Accepted curve spellings: a callable `f(angle_deg) -> GZ`, an `(angles, levers)` pair, a
`tensorcad.gz.RightingArmCurve`, or a `tensorcad.righting_arm.RightingArmCurve`.

## 2. Method

**Quadrature.** Callable curves go through adaptive Simpson to a requested tolerance with the
`|S2 - S1|/15` error estimate accumulated; the area is returned together with that error.
Tabulated curves are read through monotone cubic Hermite interpolation (Fritsch-Carlson
slopes, which do not overshoot the peak), integrated adaptively, and the *spread against the
plain trapezoidal integral on the same samples* is reported as the uncertainty. On a table
the sampling is the limiting factor, so the module also reports a self-convergence error for
`GM0`: the full table against every second sample.

**Three-valued verdicts.** `passed` is `True`, `False`, or `None`. `None` means the curve
could not answer a rule -- it stops before the angle the rule is about, GZ vanishes first, or
the flooding angle is inside the required range. A margin inside the stated numerical
uncertainty yields `indeterminate`, which also makes `passed` `None`: exact equality at a
threshold passes (the rule reads "at least", and binary round-off must not undo that), but a
value that only *might* be above the threshold is not signed off.

**Peak search.** 2.2.3 is a criterion on an angle, and a 5-degree ladder can miss a peak
sitting between two samples by whole degrees. The maximum is found by golden-section search
(the GZ curve of an intact hull is unimodal up to vanishing) over the range the ship can
actually be heeled through -- truncated at the flooding or vanishing angle -- and a maximum
found at that ceiling is reported with a note saying so.

**GM0.** Supplied from `tensorcad.stability.initial_stability` it is used as given; otherwise
it is read from the curve's slope at upright by Richardson extrapolation
(`GZ/sin(phi) = GM0 + c phi^2 + ...`), on the curve's own small-angle samples for a table.
`gm0_source` always names which of the two it was.

## 3. Results

All numbers below are printed by `tools/intact_stability_validation.py` (artefacts:
`docs/intact_stability/intact_stability_validation.json` and `.png`); nothing is transcribed
by hand.

### 3.1 Area quadrature against the closed forms

References: `GZ = GM sin(phi)` integrates to `GM (1 - cos phi)`; the wall-sided barge
`GZ = sin(phi) (GM0 + BM0/2 tan^2 phi)` integrates to
`GM0 (1 - cos) + BM0/2 (sec + cos - 2)` (B = 6 m, T = 3 m, KG = 2 m, so GM0 = 0.5 m).

| curve | interval [deg] | computed [m rad] | closed form [m rad] | rel. error |
|---|---|---|---|---|
| GZ = 0.5 sin(phi) | 0-30 | 0.066987298108 | 0.066987298108 | 4.14e-16 |
| GZ = 0.5 sin(phi) | 0-40 | 0.116977778441 | 0.116977778441 | 1.19e-16 |
| GZ = 0.5 sin(phi) | 30-40 | 0.049990480333 | 0.049990480333 | 6.94e-16 |
| wall-sided barge | 0-30 | 0.077350269190 | 0.077350269190 | 1.79e-16 |
| wall-sided barge | 0-40 | 0.152703644666 | 0.152703644666 | 1.82e-16 |
| wall-sided barge | 30-40 | 0.075353375477 | 0.075353375477 | 1.84e-16 |

### 3.2 GM0 read from the slope at upright

| curve | computed GM0 [m] | reference [m] | rel. error |
|---|---|---|---|
| GZ = 0.5 sin(phi) | 0.50000000 | 0.50000000 | 0.00e+00 |
| wall-sided barge | 0.49999988 | 0.50000000 | 2.48e-07 |

The wall-sided residual is the `phi^4` term Richardson leaves behind, not quadrature error.

### 3.3 Table coarsening: the error bar is real

Area 0-30 deg against the closed form 0.077350269190 m rad, as the sampling is coarsened:

| step [deg] | samples | area [m rad] | error [m rad] | reported uncertainty [m rad] | covered |
|---|---|---|---|---|---|
| 20 | 4 | 0.077589369 | +2.39e-04 | 2.24e-03 | True |
| 10 | 6 | 0.077356601 | +6.33e-06 | 1.16e-03 | True |
| 5 | 11 | 0.077350634 | +3.65e-07 | 2.92e-04 | True |
| 2.5 | 21 | 0.077350292 | +2.23e-08 | 7.33e-05 | True |
| 1 | 51 | 0.077350270 | +5.68e-10 | 1.17e-05 | True |

The reported uncertainty is conservative (it bounds the trapezoidal error) and it shrinks
with the sampling, so a verdict cannot be decided by the ladder the caller happened to use.

### 3.4 Verdicts, including the mesh route

| case | verdict | statuses |
|---|---|---|
| compliant analytic (GZ = 0.5 sin phi) | `pass` | all six pass |
| wall-sided barge, analytic curve | `pass` | all six pass |
| wall-sided barge, mesh -> `righting_arm_curve` -> check | `pass` | all six pass |
| tender control (GZ = 0.1 sin phi) | `fail` | 2.2.1a/b/c, 2.2.2, 2.2.4 fail |

The mesh route is a closed 12-triangle box (B = 6 m, T = 3 m, D = 8 m) heeled at 0/10/20/30/40
deg by `gz.righting_arm_curve` at constant displacement; its GZ reproduces the closed form to
8.9e-16 m, and the section 2.2 areas then land within 3e-4 relative of the closed-form areas:

| heel [deg] | mesh GZ [m] | closed form [m] | \|error\| [m] |
|---|---|---|---|
| 0 | 0.000000000000000 | 0.000000000000000 | 0.00e+00 |
| 10 | 0.089523554302420 | 0.089523554302420 | 8.33e-17 |
| 20 | 0.193664516574472 | 0.193664516574472 | 2.78e-17 |
| 30 | 0.333333333333333 | 0.333333333333333 | 1.11e-16 |
| 40 | 0.547683387507425 | 0.547683387507424 | 8.88e-16 |

## 4. What this module does not do

* **Section 2.3** (severe wind and rolling / weather criterion) needs the windage area, its
  vertical lever and a gust model -- not geometry this package holds.
* **Sections 2.4-2.5** (timber deck cargo, grain, fishing vessels, containerships, offshore
  supply vessels, ...) amend the 2.2 thresholds; none of them are implemented.
* **Damage stability** is outside for the same reason it is outside `tensorcad.gz`: no
  flooding, no permeability, no compartments.
* **Free-surface correction** is not applied here. Add `sum(i_x)/V` to the curve first
  (`gz.free_surface_correction`) -- a slack tank shifts every ordinate, and the rule is
  applied to the corrected curve, not to the bare one.

Every assessment carries this list in its `scope` field, so a JSON consumer cannot mistake a
section 2.2 pass for a statutory pass.

## 5. Pointers

| artefact | path |
|---|---|
| module | `tensorcad/statutory.py` |
| tests (28) | `tests/test_statutory.py` |
| numbers + figure generator | `tools/intact_stability_validation.py` |
| machine-readable evidence | `docs/intact_stability/intact_stability_validation.json` |
| figure | `docs/intact_stability/intact_stability_validation.png` |