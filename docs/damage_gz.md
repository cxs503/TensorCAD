# The damaged ship's righting arm

> Mission m16549 (selfstart generalist seat), TensorCAD 0.24.0.
> Code: `tensorcad/damage.py` (`FloodedArm`, `FloodedArmCurve`, `flooded_righting_arm`,
> `flooded_righting_arm_curve`). Tests: `tests/test_damage_gz.py`, oracle in
> `tests/damage_gz_oracle.py`. Numbers: `examples/damage_gz/outputs/metrics.json`
> (generated, never typed). Figures: `docs/damage_gz/*.png`.

## 1. The gap this closes

`tensorcad.damage` answered one question: **where does the ship float when a compartment
is open to the sea** -- sinkage, trim and list from the lost-buoyancy balance, plus a
damaged `GM`. `tensorcad.gz` / `tensorcad.righting_arm` answered the *intact* ship's
`GZ(phi)`. Between them sat a hole that the module docstring said out loud:

> *Not implemented, and said so: probabilistic subdivision (SOLAS part B-1, index A vs R),
> progressive flooding through internal openings, cross-flooding, the intermediate stages
> of flooding, and **`GZ` curves of the damaged ship**.*

That last one is not a niche. The residual-stability question ("she is holed at the
midship compartment and listing -- what arm is left at 20 degrees?") is asked with a
damaged `GZ` curve, and every damage-stability criterion and every `s`-factor in a
subdivision calculation is a functional of exactly this curve. So the curve had to exist
before anything downstream of it (probabilistic index `A`, residual-stability criteria)
could be written at all.

This mission adds the curve and nothing else: **no rule thresholds, no probabilistic
index.** Those need the `p`-factor tables of the regulation transcribed correctly, which
is a different job with a different kind of risk, and inventing them would be worse than
leaving the gap open.

## 2. Method

The damage module's equilibrium solver at a free heel solves three balances for three
unknowns (`z0`, `tan(trim)`, `tan(heel)`). At a **prescribed** heel there is no third
balance to solve -- the ship is held at that heel by an external moment -- so the system
becomes two equations in two unknowns:

```
V_ext - mu V_c        = Delta/rho      (buoyancy)
(LCB of the net body) = x_G            (longitudinal moment)
```

and the arm is formed from the converged state rather than rotated out of the upright one:

```
GZ(phi) = (y_B - y_G) cos(phi) + (KB - KG) sin(phi)
```

with `y_B = TCB`, `KB` measured above the keel and `y_G`, `KG` the loading inputs. Two
consequences of writing it this way, both of them tests:

* at `phi = 0` the transverse balance is *not* enforced, so `GZ(0) = y_B - y_G` is the
  ship's list, measured. A centred compartment gives exactly zero; a port-side one gives
  a number;
* the small-angle slope of the curve must be the `GM` the *equilibrium* solver reports,
  `dGZ/dphi|0 = KB + I_T/(Delta/rho) - KG = damage_gm(state, KG)`. Two different code
  paths, one number -- 1e-11 relative agreement at `phi = 1e-4` deg.

The two-variable Newton is a separate routine, not the three-variable one with a column
deleted: the residual *count* changes, so the system that gets eliminated changes too.
Warm starting from the previous heel keeps a 13-point curve to a handful of iterations per
point and makes the curve continuous in `heel_deg`.

## 3. Results

Case: watertight box barge `L=10 m`, `B=2 m`, `T=0.8 m`, deck `1.0 m` above the design
waterline, `Delta/rho = 16 m^3`, `KG = 0.5 m`; midships compartment `x = 4..6 m`, full
breadth, keel to the calm waterline, permeability `mu = 0.6`. The reference is the exact
closed-form oracle of §4 -- no grid, no refinement.

| heel [deg] | library GZ [m] | exact oracle GZ [m] | residual [m] | intact GZ [m] | wall-sided GZ [m] | waterline z [m] | KB [m] |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | -0.000000000 | 0.000000000 | -2.26e-18 | 0.000000000 | 0.000000000 | 0.096000 | 0.453760 |
| 5 | 0.032423793 | 0.032423793 | 0.00e+00 | 0.027738300 | 0.027738300 | 0.096000 | 0.455355 |
| 10 | 0.064242983 | 0.064242983 | 1.30e-14 | 0.056113367 | 0.056113367 | 0.094871 | 0.460064 |
| 15 | 0.096323362 | 0.096323362 | 3.15e-14 | 0.085830692 | 0.085830692 | 0.092556 | 0.468060 |
| 20 | 0.130194439 | 0.130194439 | 2.48e-14 | 0.117745731 | 0.117745731 | 0.089804 | 0.479932 |
| 25 | 0.166975204 | 0.166975204 | 1.65e-14 | 0.152973973 | 0.152973973 | 0.086731 | 0.496536 |
| 30 | 0.208106642 | 0.208106642 | 8.55e-15 | 0.193055556 | 0.193055556 | 0.083318 | 0.519173 |
| 35 | 0.255660001 | 0.255660001 | 5.55e-17 | 0.240219858 | 0.240219858 | 0.079493 | 0.549824 |
| 40 | 0.312771645 | 0.312771645 | -1.05e-14 | 0.297092478 | 0.297836736 | 0.075133 | 0.591571 |
| 45 | 0.375917800 | 0.375917800 | -2.66e-15 | 0.353553391 | 0.371231060 | 0.067136 | 0.643464 |
| 50 | 0.421575322 | 0.421575322 | 1.49e-14 | 0.399539808 | 0.469246369 | 0.055433 | 0.689273 |
| 55 | 0.451694972 | 0.451694972 | -4.36e-14 | 0.429380395 | 0.607470652 | 0.041005 | 0.728816 |
| 60 | 0.469985584 | 0.469985584 | 2.88e-14 | 0.447000439 | 0.815507255 | 0.022457 | 0.763795 |

Summary (`metrics.json`):

| quantity | value |
|---|---|
| intact `GM = KB + BM - KG` | 0.316666666667 m |
| damaged `GM` (free heel, three-equation solve) | 0.370426666667 m |
| damaged sinkage at the flooded condition | 0.096 m |
| damaged `GZ` at 20 deg | 0.130194439314 m |
| damaged `GZ` at 60 deg | 0.469985583686 m |
| curve area to peak, `INT GZ dphi` | 0.240059475468 m rad |
| worst |library - oracle| over the ladder | **4.36e-14 m** |

**Why the damaged arm is *larger* than the intact one here**, which looks wrong until the
bookkeeping is unpacked: the flooded compartment is full-breadth and its top is at the calm
waterline, so it removes buoyancy from the *bottom* of the hull and does not touch the
waterplane at all (`I_T` unchanged, `BM` unchanged). Losing bottom buoyancy and sinking
0.096 m to recover it raises `KB` by 0.054 m, so `KM` rises by 0.054 m and `GM` with it.
The case was chosen for that transparency, not because flooding usually helps a ship: the
opposite happens as soon as the flooded volume reaches the waterplane, and
`test_a_side_compartment_lists_the_ship` pins the one-sided case where flooding *does* cost
stability.

## 4. The exact oracle (`tests/damage_gz_oracle.py`)

The library integrates clipped section polygons at Gauss stations, so every number it
produces carries a discretisation argument. For a box hull with a midships box compartment
there is nothing to argue about: the heeled plane is independent of `x`, so the net
buoyancy is a **one-dimensional** integral over `y` of clamped-linear functions, and each
piece contributes

```
INT h dy   = a (yb-ya) + b (yb^2-ya^2)/2
INT y h dy = a (yb^2-ya^2)/2 + b (yb^3-ya^3)/3
INT h^2 dy = a^2 (yb-ya) + a b (yb^2-ya^2) + b^2 (yb^3-ya^3)/3
```

in closed form (`h(y) = a + b y`). Volume, `x`-, `y`- and `z`-moments follow by summing
the pieces, the floating attitude is a bisection on volume (monotone), and the arm comes
out of the centroid. The only assumption is zero trim -- which is the *answer* whenever
hull and compartment are both symmetric about midships, and the library is separately
asserted to return `lcb == x_G` and a vanishing trim there.

Agreement: **4.4e-14 m worst over 0..60 deg**, i.e. the section-clipping route is exact on
this geometry to the round-off of the closed forms. The independence of the panel count
(4/8/16 Gauss panels per interval) is checked separately on a *curved* hull.

## 5. Figures

* `docs/damage_gz/arm_curves.png` -- the damaged arm against the exact oracle (indistinguishable,
  which is the point), the intact arm, and the wall-sided closed form. The vertical line is
  `atan(2T/B) = 38.66 deg`, where the closed form's assumption dies: at 60 deg it is
  **0.346 m** wrong, and the integrator is not.
* `docs/damage_gz/residuals_and_sinkage.png` -- `|library - oracle|` on a log axis
  (round-off, no trend with heel) and the waterline/KB migration that produces the arm.

## 6. Error judgement against the standards

| key quantity | reference | error | <3 %? |
|---|---|---|---|
| damaged `GZ`, 0..60 deg (13 points) | exact closed-form oracle | 4.36e-14 m worst, on arms up to 0.470 m = **9.3e-14 relative** | yes |
| damaged `dGZ/dphi` at upright | `damage_gm(state, KG)` from the independent equilibrium solve | 9.2e-13 m at `phi = 1e-4` deg = 2.5e-12 relative | yes |
| intact `GZ`, 0..35 deg | wall-sided closed form `sin(phi)(GM + BM/2 tan^2 phi)` | 1e-17 m (round-off) | yes |
| intact `GZ`, 40..60 deg | wall-sided closed form | **not a valid reference** past `atan(2T/B)`; the oracle is used instead and still agrees to 2e-14 m | n/a (declared) |
| damaged `GZ` on a curved parametric hull, panels 4 vs 16 | self-convergence | < 5e-4 relative (asserted) | yes |

## 7. Notes, limits, and what is deliberately absent

* **Free-flooding, lost buoyancy, one waterline.** Same model as the rest of
  `tensorcad.damage`: compartments open to the sea, no added-weight cases (an internally
  ruptured tank), no progressive flooding through internal openings, no cross-flooding, no
  intermediate stages.
* **The ladder stops at 60 deg**, and the reason is physical, not numerical: once the
  waterline inside the damaged compartment rises above the damage opening, water drains
  back out and the lost-buoyancy bookkeeping stops describing the case. The integrator
  would happily go further; the *model* would be wrong. `flooded_righting_arm` refuses
  past that with a message that says so.
* **`KG` is required, `lcg`/`tcb` are not.** An arm without a centre of gravity is not an
  arm. `lcg`/`tcb` default to the intact ship's own centres at that displacement -- the
  weight stays where the intact ship carried it.
* **No rule check.** `vanishing_angle_deg`, `positive_range_deg`, `max_arm` and
  `area_to_max_arm` are *measurements of this curve*. Comparing them with a statutory
  threshold is a separate function with a separate loading-condition set, and it is not
  here.
* **Open item, honestly flagged**: the independent cross-check against
  `tensorcad.gz.righting_arm_curve` on a curved parametric hull (intact case, same mesh,
  same displacement) disagrees by **0.36 % at 5 deg heel** (0.0193136 m vs 0.0186471 m),
  far more than either route's own convergence error. The two routes cut the hull
  differently (`section_loops` polygons clipped in `y` versus per-triangle clipping in a
  rotated frame) but they integrate the *same* polygonal mesh, so a 0.36 % gap is
  unexplained and is **not** claimed as agreement. On the box geometry -- where an exact
  answer exists -- this module is exact, so the suspicion falls on the freeboard-band mesh
  handling in one of the two routes. Logged here rather than papered over with a loose
  tolerance; it is the obvious next mission on this line.

## 8. Reproducing

```
python examples/damage_gz/run_example.py            # metrics.json, table.md, figures
python -m pytest tests/test_damage_gz.py tests/test_damage.py -q
```