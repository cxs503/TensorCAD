# Damage stability: the damaged ship's own equilibrium

`tensorcad/damage.py` solves the free-flooding damage case -- a box-shaped space opened
to the sea -- by the lost-buoyancy method, and this page is the measurement of it.
Every number below is written by `python examples/damage/run_example.py` into
`docs/damage/results.json`, and the script exits non-zero if a claim stops holding.

## The case

| quantity | value |
|---|---|
| length L | 10.0 m |
| beam B | 2.0 m |
| design draft T | 0.8 m |
| displacement volume | 16.0 m^3 |
| `KG` (for GM only) | 0.95 m |
| freeboard | 1.0 m |

Compartment geometries: midship (x 4-6 m), aft end (x 0-2 m), port side
(x 4-6 m, y -1-0 m), and the same at `mu = 0.6` and at partial depth.

## How it is solved

The flooded water is inside the hull, so the ship's weight is unchanged and buoyancy is
lost where the water sits:

    V_ext - mu V_c = Delta / rho          (buoyancy)
    centroid(V_ext - mu V_c) = (x_G, y_G) (moment)

Both are evaluated on the *same* waterline plane
`z = z0 + tan(trim)(x - x_ref) + tan(heel) y`, whose three unknowns are found by a damped
Newton solve with a finite-difference Jacobian and a backtracking line search; section
cuts are cached once per station and only re-clipped per iteration. Compartment ends are
panel edges and each panel is integrated with three-point Gauss-Legendre, so the jump in
`A(x)` at a compartment end is never sampled across. `GM = KB + I_T/(Delta/rho) - KG`
uses the *net* waterplane (external minus `mu` times the flooded part), which is the
lost-buoyancy method's answer to the free-surface question.

## Results against the closed forms

A wall-sided barge has no discretisation, so every quantity has an exact value.

| claim | measured | reference | relative error | verdict |
|---|---|---|---|---|
| intact draft [m] | 0.8 | 0.8 | 1.39e-14 % | pass |
| intact KB [m] | 0.4 | 0.4 | 1.39e-14 % | pass |
| intact I_T [m^4] | 6.66666666667 | 6.66666666667 | 1.33e-14 % | pass |
| intact LCB [m] | 5 | 5 | 1.78e-14 % | pass |
| midship damaged draft [m] | 1 | 1 | 0.00e+00 % | pass |
| midship lost buoyancy [m^3] | 4 | 4 | 0.00e+00 % | pass |
| midship damaged I_T [m^4] | 5.33333333333 | 5.33333333333 | 0.00e+00 % | pass |
| midship trim [deg] | -1.08336837024e-14 | 0 | 1.08e-12 % | pass |
| midship GM [m] | -0.116666666667 | -0.116666666667 | 0.00e+00 % | pass |
| mu=0.6 damaged draft [m] | 0.909090909091 | 0.909090909091 | 1.22e-14 % | pass |
| aft compartment z0 [m] | 1.1875 | 1.1875 | 7.76e-12 % | pass |
| aft compartment tan(trim) | -0.1875 | -0.1875 | 4.91e-11 % | pass |
| port compartment z0 [m] | 0.897196261682 | 0.897196261682 | 1.24e-14 % | pass |
| port compartment tan(heel) | -0.14953271028 | -0.14953271028 | 1.86e-14 % | pass |
| stbd heel is the mirror [deg] | 8.50457924088 | 8.50457924088 | 4.18e-14 % | pass |

Worst relative error over the dimensional claims: 7.76e-12 %.

The two damage cases with no closed-form draft (aft end, port side) are checked against
the exact two-equation solution of the same problem on a rectangle: for the aft
compartment `16 z0 + 16 tan = 16` and `96 z0 + 181.3333 tan = 80` give
`z0 = 1.1875 m`, `tan(trim) = -0.1875` (both reproduced to round-off); for the port half
compartment `18 z0 + tan = 16` and `z0 + 6 tan = 0` give the measured
`z0 = 0.897196262 m`, `tan(heel) = -0.149532710`.

## What the figures show

`docs/damage/buoyancy_curve.png` is the damage control picture: the external `A(x)` and
the net curve the ship actually floats on, with the compartment's worth of buoyancy
shaded out. `docs/damage/sweeps.png` sweeps the compartment length: the damaged draft
against the closed form `L T/(L - l)` (exact at every length, not just one), the heel of
a port-side compartment, and `GM`, which goes negative once the damage is long enough --
reported, not clamped.

## Notes, limits, and what is deliberately absent

* **Free-flooding only.** Compartments open to the sea. An internally ruptured tank is
  the added-weight problem with a free surface; passing it here would double-count.
* **Deterministic, not probabilistic.** SOLAS part B-1 (subdivision index `A` against
  `R`), progressive flooding through internal openings, cross-flooding and the
  intermediate stages are not implemented, and no `GZ` curve of the damaged ship is
  produced. `tensorcad.intact_stability` remains the intact-ship rule check.
* **Permeability is an input.** It is never guessed; `mu = 1.0` empty, `0.6` cargo,
  `0.95` machinery, `0` full tank. The state records `lost_buoyancy_volume = mu V_c`.
* **The mesh is the remaining discretisation.** On the curved parametric hull the intact
  solve reproduces `tensorcad.hydrostatics`' own volume exactly and its draft to
  0.5 % -- the difference is the two integrations of the same mesh, not the solve.
* A hull without freeboard cannot be damaged: if the compartment takes away more
  buoyancy than the mesh has above the waterline, the case is refused with
  "the ship would be a submarine" rather than returning a clipped answer.

## Reproducing

```
python -m pytest tests/test_damage.py -q
python examples/damage/run_example.py
```
