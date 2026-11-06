# Intact-stability criteria: the `GZ` curve, judged

`tensorcad.stability` stops at `GM`, `tensorcad.gz` / `tensorcad.righting_arm` produce the
righting-arm curve, and both say the curve is only as useful as the check run on it. This
module is that check: **`tensorcad.intact_stability`** evaluates the general intact-stability
criteria of the **IMO Intact Stability Code 2008, Part A, §2.2** (successor of resolution
A.749(18) §3.1.2) and returns a per-clause verdict with achieved and required side by side.

```
from tensorcad import assess_intact_stability, assess_righting_arm_curve

report = assess_intact_stability(angles_deg, gz, initial_gm=0.62, flooding_angle_deg=48.0)
report.all_satisfied        # False when any clause fails *or* cannot be answered
report.failed               # ('IS-2008-PartA-2.2.1', ...)
report.undetermined         # clauses the sampled curve is too short to answer
report.to_dict()            # JSON-ready, schema tensorcad.stability.intact_assessment/1
```

## The six clauses and their defaults

| id | requirement | default threshold |
|---|---|---|
| `IS-2008-PartA-2.2.1` | area under `GZ` up to 30° | 0.055 m·rad |
| `IS-2008-PartA-2.2.2` | area under `GZ` up to 40°, or to the flooding angle θf if smaller | 0.090 m·rad |
| `IS-2008-PartA-2.2.3` | area between 30° and 40° (or θf) | 0.030 m·rad |
| `IS-2008-PartA-2.2.4` | `GZ` at 30° | 0.200 m |
| `IS-2008-PartA-2.2.5` | heel at which `GZ` is greatest | 25° |
| `IS-2008-PartA-2.2.6` | initial metacentric height `GM_0` | 0.150 m |

Thresholds live in `IntactThresholds` and are **replaceable**, because an Administration,
a Code revision or a ship type can carry different numbers; the ids in the report say which
clause each number came from. Nothing here paraphrases the regulation or certifies anything.

## Three verdicts, not two

`status` is `satisfied`, `not_satisfied`, or **`undetermined`**. The third is not a soft
failure: it means the sampled curve does not reach the angle the clause is quoted at (or the
caller withheld the flooding angle), so any verdict would be a guess. `achieved` is `None`
exactly when the status is `undetermined`, and `all_satisfied` is `True` only when every
clause is *positively* satisfied — a curve sampled to 20° cannot pass by omission.

Three consequences worth knowing:

- **A short curve is reported, not silently truncated.** `area_between` raises when asked
  for an interval outside the table, and the assessment turns that into `undetermined` with
  the note "widen the curve instead of assuming the arm".
- **A curve that does not start at 0° answers nothing about the areas from 0°.** The same
  rule applies at the other end: `GZ` below the first sample is unknown, not zero.
- **The peak of 2.2.5 is read on the sampled grid.** When the greatest sampled arm *is* the
  last sample, the note says the true peak lies at or beyond it — that is evidence of a
  rising curve, not of a peak.

## How the areas are measured

`int GZ dφ` is the trapezoid rule on the table, with `GZ` linearly interpolated where a
clause boundary (30°, 40°, θf) falls between samples. Angles are converted to radians
*inside* the module: a curve quoted in degrees and integrated in degrees is off by a factor
of 57.3, which is the classic way this number is got wrong. The integral is **signed** — an
interval where `GZ < 0` subtracts — which is what dynamical stability means and what stops
an early loss of stability from being flattered.

## Verification

`examples/intact_stability/run_example.py` (`PYTHONPATH=. python
examples/intact_stability/run_example.py --output docs/intact_stability --plots`) runs six
claims, exits non-zero if any stops holding, and writes `results.json`, `report.md` and
three figures. Measured on the record of 2026-11-06:

1. **The area rule on a curve with an exact integral.** For `GZ = A sin(2φ)`,
   `int_a^b GZ dφ = A/2 (cos 2a - cos 2b)`. Trapezoid on a 0.05° table is within
   **5.2e-8 m·rad**, and the fitted log–log slope of error against step is **2.000**, i.e.
   second order. (First run of this claim *failed* — the closed form in the script was wrong
   by a factor 2 and the check said so. That is the check earning its keep.)
2. **A real box barge against the wall-sided closed form.** For a wall-sided hull
   `GZ = sinφ (GM_0 + BM_0/2 tan²φ)`, which integrates to

   ```
   ∫₀^θ GZ dφ = GM_0 (1 - cos θ) + BM_0/2 (sec θ + cos θ - 2)
   ```

   derived in the module docstring (differentiate it and you get the arm back). A 10 × 1.6 ×
   2.2 m barge floating at 0.8 m draft with `KG = 0.5 m` (so `GM_0 = 0.1667 m`,
   `BM_0 = 0.2667 m`) reproduces it to **1.8e-4 m·rad worst error** over steps of 5, 2.5 and
   1.25°, on intervals 0–30° and 0–40° — against clause thresholds of 0.055 and 0.090 m·rad.
   The barge is deliberately 1.6 m wide and 2.2 m deep so its wall-sided validity limit
   (bilge emergence, `tanφ = 2T/B = 1`, i.e. 45°) sits *above* the Code's 40° boundary.
3. **Negative control.** Raising `KG` by 0.1 m must move the 0–30° area by exactly
   `-0.1 (1 - cos 30°) = -0.013397460 m·rad`, independent of `BM`; measured
   **-0.013397375 m·rad** (7e-8 apart). Any area not tracking the curve would fail this.
4. **The barge's actual verdicts.** `GM_0 = 0.1667 m` passes 2.2.6, and the 0–30° area
   `0.0251 m·rad` fails 2.2.1 by a factor 2.2. A low-freeboard box is not an IS-Code ship —
   reported as numbers, not a shrug.
5. **A hull with a real waterplane, end to end.** A Wigley hull (`L/B = 10`, `t = 0.625 m`)
   is loaded to put `GM` exactly on the 0.15 m floor (`KG = KM - 0.15`), heeled to 90° at
   5° steps; all six clauses are evaluated (`max_volume_error` of the curve is in
   `results.json`). It **fails 2.2.1–2.2.4**: area to 30° is 0.0017 m·rad and `GZ(30°)` is
   0.006 m. That is not a bug in the check: a bare Wigley has a watertight deck at the
   waterline and a very small waterplane inertia, so its arms are tiny. Reading the number
   as "this hull needs freeboard and a `θf`" is the correct reading.
6. **Truncation by a flooding angle.** With `θf = 25°` the 0–40° clause is measured to 25°
   (0.0012 m·rad) and 2.2.3 is reported `undetermined` with the note that the interval is
   empty.

Figures: `figure1_barge_areas.png` (integrated curve over the closed form, with the two
clause areas shaded), `figure2_criteria_ratios.png` (achieved/required per clause, both
hulls, with the 1.0 limit), `figure3_convergence_control.png` (second-order convergence and
the `KG`-shift control).

## What this is not

- **No §2.3 weather criterion**, no parametric rolling or other second-generation
  (Part B) failure modes, no damage stability, no grain/timber-deck rules.
- **No `θf` computation.** The flooding angle belongs to the subdivision, not to the hull
  surface; it is an input. Passing no `θf` assesses the bare hull against the full 40° clauses
  and says so in the report.
- **No free-surface correction of the curve itself.** Build the curve with
  `free_surface_i_x=...` and it carries the flat `Σi_x/V` reduction; the module then judges
  the corrected curve.
- **A watertight deck is not a freeboard.** `gz` integrates the hull it is given; if that
  hull is capped at the waterline, the arms above the deck edge are those of a closed body,
  which is optimistic for a real ship with openings. Feeding `with_freeboard(...)` (or a real
  deck) and a `θf` is the caller's job.
- **`initial_gm` should come from `initial_stability`.** When it is not supplied the module
  reads the chord `GZ(φ)/sinφ` at the first sample ≤ 10° and *says so in the note*: that
  value carries the curvature bias `BM_0/2 tan²φ` and therefore reads high (1 mm at 5° for
  `BM_0 = 0.27 m`; the test pins the bias formula exactly).