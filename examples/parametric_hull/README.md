# First geometry benchmark: parametric hull

A headless, parameterized hull surface that is exported as a **closed, oriented solid** with hydrostatics and a validation record. This directory is the reference for what every later TensorCAD model has to provide: deterministic generation, checks that can fail, machine-readable artifacts, and measured — not decorative — verification.

## Run it

```bash
python -m pip install -e ".[test]"
python examples/parametric_hull/run_example.py --output results/hull-demo
```

The script exits non-zero if any check it makes stops holding. Equivalent CLI runs:

```bash
tensorcad-hull --length 12 --beam 2.4 --draft 0.9 --nx 61 --nz 17 --output results/hull-demo
tensorcad-hull --nx 11 --nz 5 --no-export --json     # report only, no files written
tensorcad-hull --nx 40                                # a coarse grid that admits its own deficit
```

## Artifacts

| File | Contents |
|---|---|
| `hull.obj` | triangle mesh, coordinates written at 17 significant digits so a re-import is bit-identical |
| `offsets.csv` | `x_m,overall_breadth_m` — one waterline breadth row per longitudinal station |
| `metadata.json` | parameters, units, coordinate frame, mesh counts and bounds, validation checks, hydrostatics, discretisation estimate, `engineering_status` |

## Measured output (default parameters, L=10 m, B=2 m, T=0.8 m, nx=41, nz=13)

```
2. closed-solid checks
  [PASS] watertight (every edge used by exactly two faces)
  [PASS] edge-manifold
  [PASS] no zero-area faces
  [PASS] no boundary edges
  [PASS] outward orientation (signed volume > 0): volume = 7.446554 m^3
  [PASS] Euler characteristic V-E+F = 2: V=1001 F=1998 -> 2
3. dimensional fidelity
  [PASS] length / draft / beam realised, waterline plane exactly z=0
4. hydrostatics vs analytic value of the sampled surface
     mesh volume       = 7.446554 m^3
     analytic volume     = 7.512698 m^3      (-0.88 %)
     wetted surface    = 22.838010 m^2
     waterplane area   = 14.061669 m^2
     LCB / LCF / VCB   = 5.0000 / 5.0000 / -0.3170 m
     Cb / Cp           = 0.4654 / 0.7031
5. resolution study (error must shrink as the grid refines)
       nx   nz     volume [m^3]    V err %     Awp [m^2]    Awp err %
       11   5       7.115242     -5.290     13.826032     -1.848
       21   9       7.374345     -1.842     14.006243     -0.568
       41  13       7.446554     -0.880     14.061669     -0.175
       81  21       7.483649     -0.387     14.078726     -0.054
      161  41       7.502488     -0.136     14.083975     -0.017
```

The analytic reference is the exact integral of the generator's own shape functions (`∫₀¹ sin(πx)^0.7 dx = (2/π)·√π·Γ(0.85)/(2Γ(1.35))`, and `∫₀¹ √v dv = 2/3`), so the residual is grid error, not a disagreement between two models. The residual converges as the grid refines, which is the point of printing it.

## Scope and limitations

The breadth distribution is intentionally simple and fore-aft symmetric: no transom, no entry/exit asymmetry, a knuckle at the keel (the `√` profile has an infinite slope there, which is why volume convergence is slower than second order), and no freeboard — the deck is the planar waterline cap.

The checks here prove closure, orientation, manifoldness, index integrity, and dimensional fidelity. They do not prove fairness, hydrostatic accuracy beyond grid error, stability, resistance, or class-rule compliance. STEP/B-rep export and CAD-kernel validity are the next step (issue #4), deliberately kept out of the core package's dependency set.
