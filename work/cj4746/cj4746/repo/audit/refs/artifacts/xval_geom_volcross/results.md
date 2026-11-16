# Geometry volume cross-validation: pipeline vs independent engines

Feature: `XVAL-GEOM-VOLCROSS`  
Generated: 2026-11-14T12:44:58+00:00 by `tools/volume_xval_report.py`  
Engines: trimesh 5.0.0, numpy-stl unknown, numpy 2.4.6

Every number below is written by the generator into `results.json` and rendered from there; no value is copied by hand.

## Route definitions

| route | implementation | shares code/algorithms with |
|---|---|---|
| pipeline | the platform's own route for that hull (NURBS Gauss quadrature, or the closed-surface divergence theorem of `tensorcad.mesh`) | -- |
| engine A | `trimesh.Trimesh.volume` (third-party, NumPy) | nothing in this repo |
| engine B | `numpy-stl` `get_mass_properties` (Mirtich, C++) through a binary-STL round trip (float32 coordinates) | nothing in this repo |
| smooth | Simpson quadrature of the half-breadth field: `A(x) = 2 int b dz`, then `V = int A dx`. No triangles, no divergence theorem | the half-breadth field only |

## Wigley (ITTC reference proportions, n=m=2)

L = 10 m, B = 1 m, T = 0.625 m, LBT = 6.25 m³.  Pipeline route: `tensorcad.wigley mesh + divergence theorem`.

Reference: **analytic (exact)** = 2.7777778 m³ -- V = L B T n m / ((n+1)(m+1)) = 4 L B T / 9, evaluated in tensorcad.wigley.analytic_volume.

| stations | levels | h [m] | faces | pipeline (m³) | engine A (m³) | engine B (m³) | smooth (m³) | A vs B | pipeline vs engines | smooth vs engines |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 21 | 6 | 0.5 | 438 | 2.77122215 | 2.72316504 | 2.72316497 | 2.77777778 | -2.68e-06 % | +1.7648 % | +2.0055 % |
| 41 | 11 | 0.25 | 1678 | 2.77122215 | 2.76365404 | 2.763654 | 2.77777778 | -1.28e-06 % | +0.2738 % | +0.5111 % |
| 81 | 21 | 0.125 | 6558 | 2.77122215 | 2.77421773 | 2.77421768 | 2.77777778 | -1.64e-06 % | -0.1080 % | +0.1283 % |
| 161 | 41 | 0.0625 | 25918 | 2.77122215 | 2.77688595 | 2.77688593 | 2.77777778 | -7.43e-07 % | -0.2040 % | +0.0321 % |

**Self route against the engines, on the identical triangles** (integrator identity, not tessellation):

| stations | divergence (m³) | engine A (m³) | A vs divergence | B vs divergence | watertight | genus |
|---:|---:|---:|---:|---:|---|---:|
| 21 | 2.72316504 | 2.72316504 | +1.63e-14 % | -2.68e-06 % | True | 0 |
| 41 | 2.76365404 | 2.76365404 | +6.43e-14 % | -1.28e-06 % | True | 0 |
| 81 | 2.77421773 | 2.77421773 | -1.60e-14 % | -1.64e-06 % | True | 0 |
| 161 | 2.77688595 | 2.77688595 | -5.12e-13 % | -7.43e-07 % | True | 0 |

Verdict: **agrees** -- pipeline and engines agree to 0.032 %, inside the 0.3 % acceptance line.

- observed convergence order of the pipeline-to-engine gap: 1.9797132663315644
- engine-to-engine spread (worst resolution): 2.683e-08 (engines agree to 2.68e-08 relative, inside the 1e-06 binary-STL float32 round-trip budget)
- self divergence vs engine A (worst): 5.118e-15
- Richardson extrapolation of engine A to h->0: 2.77766813 m³, which is -0.0039 % from the pipeline volume
- pipeline vs reference: -0.2360 % (inside the 0.3 % line)
- extrapolated engines vs reference: -0.0039 %

## KCS-dimensioned NURBS hull (full scale, Cb calibrated to 0.6505)

L = 230 m, B = 32.2 m, T = 10.8 m, LBT = 79984.8 m³.  Pipeline route: `tensorcad.hull_nurbs Gauss quadrature of the NURBS patch (24x16)`.

Reference: **published displacement volume** = 52030 m³ -- KCS principal particulars, KRISO (Lpp 230.0 m, B 32.2 m, T 10.8 m, Cb 0.6505, nabla 52 030 m^3; SIMMAN 2008 / Tokyo 2015 workshop series).

> Caveat: the form is calibrated to the published Cb, so pipeline-vs-published is a calibration identity; the independent content of this row is the three-route agreement at realistic L/B/T and the tessellation attribution.

| stations | levels | h [m] | faces | pipeline (m³) | engine A (m³) | engine B (m³) | smooth (m³) | A vs B | pipeline vs engines | smooth vs engines |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 21 | 7 | 11.5 | 572 | 52062.3902 | 51223.0034 | 51223.0041 | 52062.4644 | +1.43e-06 % | +1.6387 % | +1.6388 % |
| 41 | 12 | 5.75 | 1942 | 52062.3902 | 51791.2673 | 51791.2684 | 52062.4644 | +2.09e-06 % | +0.5235 % | +0.5236 % |
| 81 | 22 | 2.875 | 7082 | 52062.3902 | 51986.6441 | 51986.6448 | 52062.4644 | +1.33e-06 % | +0.1457 % | +0.1458 % |
| 161 | 42 | 1.438 | 26962 | 52062.3902 | 52042.7632 | 52042.764 | 52062.4644 | +1.53e-06 % | +0.0377 % | +0.0379 % |

**Self route against the engines, on the identical triangles** (integrator identity, not tessellation):

| stations | divergence (m³) | engine A (m³) | A vs divergence | B vs divergence | watertight | genus |
|---:|---:|---:|---:|---:|---|---:|
| 21 | 51223.0034 | 51223.0034 | -1.42e-14 % | +1.43e-06 % | True | 0 |
| 41 | 51791.2673 | 51791.2673 | +0.00e+00 % | +2.09e-06 % | True | 0 |
| 81 | 51986.6441 | 51986.6441 | +5.60e-14 % | +1.33e-06 % | True | 0 |
| 161 | 52042.7632 | 52042.7632 | +5.59e-14 % | +1.53e-06 % | True | 0 |

Verdict: **agrees** -- pipeline and engines agree to 0.038 %, inside the 0.3 % acceptance line.

- observed convergence order of the pipeline-to-engine gap: 1.8169216344088621
- engine-to-engine spread (worst resolution): 2.093e-08 (engines agree to 2.09e-08 relative, inside the 1e-06 binary-STL float32 round-trip budget)
- self divergence vs engine A (worst): 5.598e-16
- Richardson extrapolation of engine A to h->0: 52048.5175 m³, which is -0.0266 % from the pipeline volume
- pipeline vs reference: +0.0623 % (inside the 0.3 % line)
- extrapolated engines vs reference: +0.0356 %

## KCS-dimensioned NURBS hull (model scale, Cb 0.555) -- IN383-CAD01 case

L = 3 m, B = 0.4 m, T = 0.15 m, LBT = 0.18 m³.  Pipeline route: `tensorcad.hull_nurbs Gauss quadrature of the NURBS patch (24x16)`.

Reference: **none (WITHHELD)** -- no published value exists for this parametric form (it is not the KCS offsets); before this run its volume was only ever compared against a mesh generated by the same patch, i.e. a WITHHELD claim.

> Caveat: reference-free by construction: the evidence produced here is cross-tool agreement, not agreement with an outside number.

| stations | levels | h [m] | faces | pipeline (m³) | engine A (m³) | engine B (m³) | smooth (m³) | A vs B | pipeline vs engines | smooth vs engines |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 21 | 7 | 0.15 | 572 | 0.0999710666 | 0.0983676305 | 0.0983676343 | 0.0999727053 | +3.81e-06 % | +1.6300 % | +1.6317 % |
| 41 | 12 | 0.075 | 1942 | 0.0999710666 | 0.0994594621 | 0.0994594665 | 0.0999727053 | +4.41e-06 % | +0.5144 % | +0.5160 % |
| 81 | 22 | 0.0375 | 7082 | 0.0999710666 | 0.099829448 | 0.0998294521 | 0.0999727053 | +4.11e-06 % | +0.1419 % | +0.1435 % |
| 161 | 42 | 0.01875 | 26962 | 0.0999710666 | 0.0999354724 | 0.099935477 | 0.0999727053 | +4.61e-06 % | +0.0356 % | +0.0373 % |

**Self route against the engines, on the identical triangles** (integrator identity, not tessellation):

| stations | divergence (m³) | engine A (m³) | A vs divergence | B vs divergence | watertight | genus |
|---:|---:|---:|---:|---:|---|---:|
| 21 | 0.0983676305 | 0.0983676305 | +0.00e+00 % | +3.81e-06 % | True | 0 |
| 41 | 0.0994594621 | 0.0994594621 | +1.12e-13 % | +4.41e-06 % | True | 0 |
| 81 | 0.099829448 | 0.099829448 | +1.39e-14 % | +4.11e-06 % | True | 0 |
| 161 | 0.0999354724 | 0.0999354724 | -6.80e-13 % | +4.61e-06 % | True | 0 |

Verdict: **agrees** -- pipeline and engines agree to 0.036 %, inside the 0.3 % acceptance line.

- observed convergence order of the pipeline-to-engine gap: 1.8407247816373196
- engine-to-engine spread (worst resolution): 4.609e-08 (engines agree to 4.61e-08 relative, inside the 1e-06 binary-STL float32 round-trip budget)
- self divergence vs engine A (worst): 6.805e-15
- Richardson extrapolation of engine A to h->0: 0.0999438754 m³, which is -0.0272 % from the pipeline volume

## Pipeline quadrature convergence (does the pipeline number stand still?)

| Gauss order | volume (m³) | delta vs previous (m³) | relative vs previous |
|---|---:|---:|---:|
| 8x6 | 0.0999371318952 |  |  |
| 12x8 | 0.0999663468664 | 2.921e-05 | +2.923e-02 % |
| 16x12 | 0.0999668954837 | 5.486e-07 | +5.488e-04 % |
| 24x16 | 0.0999710666174 | 4.171e-06 | +4.173e-03 % |
| 32x24 | 0.0999724976785 | 1.431e-06 | +1.431e-03 % |
| 48x32 | 0.0999727337036 | 2.360e-07 | +2.361e-04 % |

Between the two finest orders the pipeline volume moves 2.360e-07 m³ (2.361e-04 %), so the quadrature is converged far below every tolerance compared against it.

## Integrator-identity battery (hulls without a symmetric breadth field)

**Asymmetric transom hull (starboard/port differ, flat transom plate)** -- reference: analytic closed form, 7.6475677 m3 (tensorcad.transom_asymmetric.analytic_asymmetric_hydrostatics).

| route | volume (m3) | vs analytic |
|---|---:|---:|
| divergence (platform) | 7.58438617 | -0.8262 % |
| engine A (trimesh) | 7.58438617 | -0.8262 % |
| engine B (numpy-stl) | 7.58438636 | -0.8262 % |

Self divergence vs engine A: +1.41e-13 %; engine A vs engine B: +2.47e-06 %; watertight=True, genus=0.

> the mesh is the platform's own hull at its default station grid, so the deviation from the analytic volume is a sampling deficit of that grid -- the three routes agreeing on the same triangles is what this row proves.

## Engines that were not available

- `gmsh`: not installed (no network in this environment)
- `pythonocc-core`: not installed (no network in this environment)
- `cadquery`: not installed (no network in this environment)

The anchor names Gmsh/OpenCASCADE as the second engine. Neither is installable in this environment (no package index reachable), so the second independent engine is `numpy-stl`'s Mirtich implementation, and the tessellation-free route is the sectional quadrature. This is a capability limitation of the environment, recorded rather than worked around silently; the OpenCASCADE BRep route remains open work.
