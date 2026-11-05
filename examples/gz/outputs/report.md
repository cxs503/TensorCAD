# Righting-arm curves: two exact benchmarks and a hull

- box barge 10.0 x 2.0 x 2.0 m, draft 0.8 m, KG 0.5 m:
  worst |GZ - wall-sided| = 6.939e-17 m over 0-38 deg
  (wall-sided limit 38.66 deg); past it the closed form drifts by
  1.252e+00 m at 70 deg
- log R=1.0 m half submerged: worst |GZ - (R-KG)sin(phi)| = 1.110e-16 m over 0-90 deg
- hull: GZ'(0) = 0.421622996 m/rad vs GM = 0.421622996 m (1.77e-07 %)
  max GZ = 0.5589 m at 90 deg, worst displacement residual 8.59e-14
- free surface (2.0 x 0.5 m tank): -2.80 mm on every point

**Verdict: all checks passed**
