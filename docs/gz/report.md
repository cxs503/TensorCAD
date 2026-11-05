# Large-angle stability: the GZ curve from heeled plane cuts

Produced by `python examples/gz/run_example_gz_curve.py`; every number below is read from
`results.json` written by the same run, and the script exits non-zero if a claim breaks.

## Case

* wall-sided barge: L=30.0 m, B=6.0 m, T=2.0 m, depth 10.0 m, KG=1.0 m
* round bilge: R=3.0 m, L=20.0 m, 360 facets, KG=1.5 m
* hull: `generate_hull(12.0, 2.4, 0.9)`, KG=0.40 m (capped at its waterline: no freeboard)

## Barge: measured against the wall-sided closed form

| heel [deg] | KN [m] | GZ [m] | wall-sided GZ [m] | diff [m] | waterline z [m] |
|---:|---:|---:|---:|---:|---:|
| 0 | -0.000000 | -0.000000 | 0.000000 | -5.26e-17 | 2.000000 |
| 5 | 0.218390 | 0.131234 | 0.131234 | -8.33e-17 | 1.992389 |
| 10 | 0.438170 | 0.264521 | 0.264521 | -1.11e-16 | 1.969616 |
| 15 | 0.660984 | 0.402165 | 0.402165 | -2.22e-16 | 1.931852 |
| 20 | 0.889032 | 0.547012 | 0.547012 | 2.22e-16 | 1.879385 |
| 25 | 1.125467 | 0.702849 | 0.702849 | 0.00e+00 | 1.812616 |
| 30 | 1.375000 | 0.875000 | 0.875000 | -3.33e-16 | 1.732051 |
| 35 | 1.642643 | 1.069066 | 1.071279 | — | 1.637293 |
| 40 | 1.894029 | 1.251241 | 1.303616 | — | 1.509324 |
| 45 | 2.121320 | 1.414214 | 1.590990 | — | 1.342781 |
| 50 | 2.332467 | 1.566423 | 1.965063 | — | 1.139554 |
| 60 | 2.740806 | 1.874781 | 3.247595 | — | 0.625634 |

GM upright = 1.500000 m; displaced volume 360.000 m^3 held constant at every angle.

## Round bilge: measured against circular-segment geometry

| heel [deg] | GZ [m] | analytic (R-KG) sin(phi) [m] | waterline z [m] | R cos(phi) [m] |
|---:|---:|---:|---:|---:|
| 0 | -0.000000 | 0.000000 | 3.000000 | 3.000000 |
| 5 | 0.130724 | 0.130734 | 2.988585 | 2.988584 |
| 10 | 0.260452 | 0.260472 | 2.954425 | 2.954423 |
| 15 | 0.388199 | 0.388229 | 2.897781 | 2.897777 |
| 20 | 0.512991 | 0.513030 | 2.819085 | 2.819078 |
| 25 | 0.633879 | 0.633927 | 2.718934 | 2.718923 |
| 30 | 0.749943 | 0.750000 | 2.598092 | 2.598076 |
| 35 | 0.860299 | 0.860365 | 2.457477 | 2.457456 |
| 40 | 0.964108 | 0.964181 | 2.298160 | 2.298133 |
| 45 | 1.060579 | 1.060660 | 2.121354 | 2.121320 |
| 50 | 1.148979 | 1.149067 | 1.928404 | 1.928363 |
| 60 | 1.298939 | 1.299038 | 1.500057 | 1.500000 |

## IMO IS Code 2008 Part A 2.2 screening

| criterion | requirement | measured | met |
|---|---|---:|:--:|
| `area_0_30` | >= 0.055 m rad | 0.2169 m rad | yes |
| `area_0_40` | >= 0.090 m rad | 0.4030 m rad | yes |
| `area_30_40` | >= 0.030 m rad | 0.1861 m rad | yes |
| `gz_at_30` | >= 0.20 m | 0.8750 m | yes |
| `max_gz` | >= 0.20 m | 1.8748 m | yes |
| `angle_at_max_gz` | >= 25 deg | 60.0000 deg | yes |
| `gm_upright` | >= 0.15 m | 1.5000 m | yes |

Negative control, same barge with `KG = 2.45 m` (`GM = 0.050 m`): `passed = False` on `area_0_30`, `area_0_40`, `gz_at_30`, `gm_upright` -- the screen is not a rubber stamp.

## Claims

| claim | measured | tolerance | met |
|---|---:|---:|:--:|
| wall-sided barge: max |GZ - closed form| over the wall-sided range | 3.331e-16 | 1e-09 | yes |
| round bilge: max |GZ - (R - KG) sin(phi)| at 360 facets | 9.893e-05 | 2e-03 | yes |
| round bilge: max |floatation - R cos(phi)| (chord through the centre) | 5.712e-05 | 2e-03 | yes |
| hull: GM from the curve equal to initial_stability at the same KG | 0.000e+00 | 1e-12 | yes |

## What this does not do

* constant trim: the waterline is kept parallel to the still-water line, so this is the
  classical cross curve, not a free-to-trim solve; a hull that trims as it heels will
  differ, and the difference is second order in the trim angle;
* no free-surface, grain or heeling-moment allowance, no damage stability, no class approval;
* `generate_hull` caps the model at the design waterline, i.e. zero freeboard, so its deck
  edge is under water from the first degree of heel: a GZ curve on such a model is the
  no-flooding idealisation and is only quoted here over 0-20 deg. Large-angle work needs a
  model with freeboard, which is what the barge fixture is;
* a waterline that lies *along* mesh edges has no section at all; the module raises rather
  than reporting half an area (see the closure guard in `tensorcad/gz.py`).
