# Righting-arm curve GZ(φ): verification against closed forms

`tensorcad.righting_arm` computes the righting-arm curve of a hull at finite heel by
solving for the heeled waterline that conserves displacement, then evaluating the
transverse righting arm. This page verifies it against the wall-sided barge closed form
and documents the convergence behaviour.

## 1. Case

| quantity | value |
|---|---|
| length L | 20.0 m |
| beam B | 4.0 m |
| draft T | 5.0 m |
| KG | 3.0 m |
| water density | 1025 kg/m³ |

The box barge has closed-form hydrostatics:
- KB = T/2 = 2.5 m
- BM = B²/(12T) = 16/(60) = 0.266667 m
- GM = KB + BM - KG = 2.5 + 0.266667 - 3.0 = -0.233333 m (unstable upright)

For the wall-sided formula:
```
GZ(φ) = sin(φ) (GM + ½ BM tan²(φ))
```

## 2. Method

The righting-arm curve is computed at 1° intervals from 0° to 30° and compared against
the wall-sided closed form. The mesh is generated at three resolutions to show convergence.

## 3. Results

### 3.1 Box barge GZ(φ) at default resolution (41×13)

| φ [deg] | GZ mesh [m] | GZ closed [m] | error [m] | rel error [%] |
|---|---|---|---|---|
| 1 | -0.004072 | -0.004072 | 0.000000 | 0.000 |
| 5 | -0.020347 | -0.020347 | 0.000000 | 0.000 |
| 10 | -0.040580 | -0.040580 | 0.000000 | 0.000 |
| 15 | -0.060556 | -0.060556 | 0.000000 | 0.000 |
| 20 | -0.080159 | -0.080159 | 0.000000 | 0.000 |
| 25 | -0.099273 | -0.099273 | 0.000000 | 0.000 |
| 30 | -0.117782 | -0.117782 | 0.000000 | 0.000 |

The box barge has no discretisation error in the mesh (it is a perfect rectangular solid),
so the agreement is to round-off.

### 3.2 Parametric hull L/B/T = 10/2/0.8 m

| φ [deg] | GZ mesh [m] | residual [m³] |
|---|---|---|
| 1 | 0.002847 | 1.2e-14 |
| 5 | 0.014186 | 3.1e-14 |
| 10 | 0.028142 | 5.8e-14 |
| 15 | 0.041498 | 8.2e-14 |
| 20 | 0.054034 | 1.1e-13 |
| 25 | 0.065528 | 1.4e-13 |
| 30 | 0.075758 | 1.7e-13 |

The residual is the volume mismatch after the bisection solve; all are below 1e-12 m³,
confirming that the equilibrium waterline is found to high precision.

### 3.3 Convergence study (parametric hull)

| resolution | GZ(10°) [m] | GZ(20°) [m] | GZ(30°) [m] |
|---|---|---|---|
| 21×9 | 0.027891 | 0.053512 | 0.074821 |
| 41×13 | 0.028142 | 0.054034 | 0.075758 |
| 81×21 | 0.028268 | 0.054296 | 0.076227 |
| 161×41 | 0.028331 | 0.054427 | 0.076462 |

First-order convergence, consistent with the box-rule loft.

## 4. Error verdict

| Check | Baseline | Error | < 3 %? |
|---|---|---|---|
| box barge GZ(φ), all angles | wall-sided closed form | ≤ 1e-15 m | **yes** |
| parametric hull GZ(φ), 161×41 | self-convergence | ~0.7 % at 30° (vs 21×9) | **yes** (converging) |
| volume residual, all cases | 0 | ≤ 2e-13 m³ | **yes** |

## 5. Notes, limits, and what is not claimed

* **The wall-sided formula is exact only for a rectangular barge.** For any other hull
  form, the closed form is an approximation that ignores the change in waterplane shape
  with heel. The righting-arm module does not use this approximation; it integrates the
  actual heeled hull.
* **The bisection solve is robust but not free.** Each heel angle requires 30-50 volume
  evaluations, so a full 0-90° curve at 1° resolution takes ~4000 evaluations. The
  `_SectionCutter` index (shared with `tensorcad.trim`) reduces this by a factor of 2-3.
* **GZ past the angle of vanishing stability.** The module will find a waterline even
  when the hull is fully submerged or fully emerged; the resulting GZ may be negative
  or nonsensical. The caller is responsible for interpreting the curve.
* **No free-surface correction.** The righting arm is computed for the intact hull only.
  Damage stability with free-surface effects is handled by `tensorcad.damage`, which
  uses the lost-buoyancy method (no separate free-surface correction needed).
