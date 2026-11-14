# CASE-CAD-CURVATURE-SPHERE: discrete mean curvature of a sphere vs the analytic `2/R`

**Anchor** `CASE-CAD-CURVATURE-SPHERE` (benchmark case pool, owner 2026-11-12).
**Criterion** relative error below **2 %**, declared before the run.
**Runner** `benchmarks/verified/sphere_mean_curvature/run_benchmark.py` (generated
`result.json`, `gates_verified.json`, `figures/`). Every number in this document is
rendered from `result.json` by `make_report.py`; none is measured by hand.

## 1. Question and reference

Does `tensorcad.curvature.vertex_mean_curvature` -- the cotangent-Laplacian mean
curvature a solver actually consumes -- reproduce the analytic mean curvature of a
sphere of radius `R`?

`kappa_1 = kappa_2 = 1/R` (do Carmo, *Differential Geometry of Curves and Surfaces*,
Prentice-Hall 1976, Sec. 3-3 Example 2; the same statement is already adjudicated in
`docs/fairness.md` and gated at `1e-6` 1/m by `benchmarks/verified/curvature_fairness`
gate G1). The estimator returns the mean-curvature normal of Meyer/Desbrun/Schroeder/Barr
(*VisMath* 2003, Sec. 3.3), ``|K| = 2|H|``, so the quantity it estimates is

```
S = kappa_1 + kappa_2 = 2 / R          (sum convention; S = 2 x vertex_mean_curvature)
```

with `S = 2.000000 1/m` at `R = 1 m`. The value is
evaluated by `reference_sum(R)` in the runner, never copied from a table, and it is
corroborated independently by the integral identity
`int_S S dA = 8 pi R = 25.132741` at `R = 1` -- which is the
`E_bending = 8 pi` closed form already gated by `curvature_fairness` (G2/G3).

Two readings are measured, because "the discrete mean curvature" is ambiguous:

* **integral** -- area-weighted mean of `S`, the reading a bending-energy or wetted-area
  integral consumes, and the one the cotangent operator conserves exactly on a closed mesh;
* **pointwise** -- median of `|S_i|` over vertices, the reading a curvature-driven mesher
  or a fairing loop locally sees.

## 2. Method

1. Analytic reference: `S = 2/R`, plus the scale-free identity `int_S S dA = 8 pi R`.
2. Measurement chain check (G4): `principal_curvatures` on the exact SDF, 40 Fibonacci
   probes plus the six poles, 7 radii -- it must return `2/R` to round-off, otherwise the
   comparison below is measuring the harness.
3. Mesh path, unstructured: `to_mesh(SphereField(R), level, margin=3)` (marching tets),
   levels `[12, 16, 24, 32, 48]`; `S_i =
   2 * vertex_mean_curvature(mesh)`; area weights = one third of the incident triangle area.
4. Mesh path, structured control: a legal UV sphere (one vertex per pole, triangle fan,
   no zero-area face) at `8x16, 16x32, 32x64, 48x96`,
   which separates "the operator is inaccurate" from "the mesher starved the operator".
5. Conditioning probe: `taubin_smooth(10)` + `project_to_field(2)` on the marching-tet
   meshes, to test whether the pointwise error is the mesh or the operator.
6. Radius sweep `0.1 ... 10 m` at level 24 to show the error is a dimensionless property
   of the discretisation.

## 3. Parameters (all dimensionless groups)

| parameter | value | why it is the interesting one |
|---|---|---|
| `S * R` (the dimensionless curvature) | 2 (exact) | the only group in the problem: a sphere has no other length scale |
| normalised mesh size `h/R` (marching tets) | 0.0417 … 0.1667 | resolution in units of the radius |
| normalised mesh size `h/R` (UV) | 0.0654 … 0.3927 | resolution of the structured control |
| radius sweep `R` | 0.1, 0.5, 1, 2, 10 m | scale invariance: `S*R` must stay 2 |
| tolerance | 2 % relative | the case criterion, predeclared |
| probes per radius (field path) | 40 Fibonacci + 6 poles | measurement-chain check only |
| `taubin` iterations / projection iterations | 10 / 2 | the repository's standard conditioning |

## 4. Reference vs measured

| gate | quantity | reference | measured | error | tol | verdict |
|---|---|---|---|---|---|---|
| G4 | field-path kappa_1+kappa_2, sphere R=0.1 m, 40 Fibonacci probes | 20 | 20 | -0.0000 % | 1e-06 rel | ✔ |
| G4 | field-path kappa_1+kappa_2, sphere R=0.25 m, 40 Fibonacci probes | 8 | 7.9999998 | -0.0000 % | 1e-06 rel | ✔ |
| G4 | field-path kappa_1+kappa_2, sphere R=0.5 m, 40 Fibonacci probes | 4 | 3.9999999 | -0.0000 % | 1e-06 rel | ✔ |
| G4 | field-path kappa_1+kappa_2, sphere R=1 m, 40 Fibonacci probes | 2 | 2 | -0.0000 % | 1e-06 rel | ✔ |
| G4 | field-path kappa_1+kappa_2, sphere R=2 m, 40 Fibonacci probes | 1 | 0.99999998 | -0.0000 % | 1e-06 rel | ✔ |
| G4 | field-path kappa_1+kappa_2, sphere R=5 m, 40 Fibonacci probes | 0.4 | 0.39999999 | -0.0000 % | 1e-06 rel | ✔ |
| G4 | field-path kappa_1+kappa_2, sphere R=10 m, 40 Fibonacci probes | 0.2 | 0.2 | -0.0000 % | 1e-06 rel | ✔ |
| G1 | discrete mean curvature, area-weighted integral, 1 m sphere, marching tets level 48 | 2 | 2.0002378 | +0.0119 % | +2.0 % | ✔ |
| G5 | surface integral of kappa_1+kappa_2 over the mesh, level 48 | 25.132741 | 25.124262 | -0.0337 % | +2.0 % | ✔ |
| G2 | pointwise discrete mean curvature, median relative error, structured UV sphere 1 m, 16x32 | 0 | 0.0039387176 | +3.939e-03 (abs) | 0.02 abs | ✔ |
| G3 | observed order of the pointwise estimate on the structured sphere | 2 | 1.9462042 | -5.380e-02 (abs) | 0.5 abs | ✔ |
| G6 | discrete mean curvature, area-weighted integral, sphere R=0.1 m, level 24 | 20 | 20.009798 | +0.0490 % | +2.0 % | ✔ |
| G6 | discrete mean curvature, area-weighted integral, sphere R=0.5 m, level 24 | 4 | 4.0019596 | +0.0490 % | +2.0 % | ✔ |
| G6 | discrete mean curvature, area-weighted integral, sphere R=1 m, level 24 | 2 | 2.0009798 | +0.0490 % | +2.0 % | ✔ |
| G6 | discrete mean curvature, area-weighted integral, sphere R=2 m, level 24 | 1 | 1.0004899 | +0.0490 % | +2.0 % | ✔ |
| G6 | discrete mean curvature, area-weighted integral, sphere R=10 m, level 24 | 0.2 | 0.20009798 | +0.0490 % | +2.0 % | ✔ |
| G7 (diag) | pointwise discrete mean curvature, median relative error, raw marching tets level 32 | 0 | 0.50294436 | +50.29 % | +2.0 % | reported, not gated |

Headline: `S = 2.000238 1/m` against `2/R = 2.000000 1/m`,
**+0.0119 %** relative error, i.e. **inside the 2 % criterion**
and inside 3 %. The pointwise reading on the same mesh is
**+50.12 %** raw and
**+1.94 %** after conditioning;
the structured control is **+0.39 %**
at `16x32`. Those two are the honest limits of the claim: the *integral* passes with two
orders of margin, the *pointwise* value passes only on a well-shaped triangulation.

### Field path (measurement-chain check)

| R [m] | 2/R [1/m] | field-path S [1/m] | relative error |
|---|---|---|---|
| 0.1 | 20 | 20 | -0.000002 % |
| 0.25 | 8 | 7.9999998 | -0.000002 % |
| 0.5 | 4 | 3.9999999 | -0.000002 % |
| 1 | 2 | 2 | -0.000002 % |
| 2 | 1 | 0.99999998 | -0.000002 % |
| 5 | 0.4 | 0.39999999 | -0.000002 % |
| 10 | 0.2 | 0.2 | -0.000002 % |

## 5. Convergence and sensitivity

**Resolution (marching tets, `R = 1 m`).**

| level | h [m] | faces | area err | integral `S` [1/m] | integral err | pointwise median | pointwise p95 | conditioned pointwise median |
|---|---|---|---|---|---|---|---|---|
| 12 | 0.1667 | 3680 | -0.781 % | 2.004167 | +0.2084 % | +48.35 % | +489.8 % | — |
| 16 | 0.1250 | 6656 | -0.429 % | 2.002254 | +0.1127 % | +45.73 % | +390.2 % | +2.08 % |
| 24 | 0.0833 | 15584 | -0.186 % | 2.000980 | +0.0490 % | +48.94 % | +378.2 % | — |
| 32 | 0.0625 | 28032 | -0.104 % | 2.000548 | +0.0274 % | +50.29 % | +441.8 % | +1.98 % |
| 48 | 0.0417 | 63592 | -0.046 % | 2.000238 | +0.0119 % | +50.12 % | +377.3 % | +1.94 % |

The integral error falls by roughly a factor of two per refinement rung
(`+0.2084 %` at level 12,
`+0.0119 %` at level 48): the
marching-tet surface is a facet approximation, so the area (and with it the integral) is
first order in `h`, and the *value* error tracks it. The raw pointwise median does **not**
converge -- it sits at ~50 % at every level, i.e. the sliver population is a property of
the connectivity and refinement does not remove it.

**Mesh family (structured UV sphere, `R = 1 m`).**

| nu x nv | h [m] | vertices | area err | integral err | pointwise median | pointwise p95 | pointwise max |
|---|---|---|---|---|---|---|---|
| 8x16 | 0.3927 | 114 | -3.181 % | -0.00430 % | +1.484 % | +3.498 % | +22.15 % |
| 16x32 | 0.1963 | 482 | -0.801 % | +0.00042 % | +0.394 % | +1.079 % | +24.28 % |
| 32x64 | 0.0982 | 1986 | -0.201 % | +0.00007 % | +0.100 % | +0.283 % | +24.82 % |
| 48x96 | 0.0654 | 4514 | -0.089 % | +0.00002 % | +0.045 % | +0.123 % | +24.92 % |

Observed order of the pointwise median over the first three rungs:
**1.946** against the theoretical 2
(cotangent estimator on a regular triangulation). The structured control is what
proves the operator is second order; the unstructured family above is what proves the
mesher, not the operator, is the limiter.

**Radius sweep (scale invariance, level 24).**

| R [m] | faces | discrete S [1/m] | 2/R [1/m] | integral err | pointwise median |
|---|---|---|---|---|---|
| 0.1 | 15584 | 20.009798 | 20.000000 | +0.0490 % | +48.94 % |
| 0.5 | 15584 | 4.001960 | 4.000000 | +0.0490 % | +48.94 % |
| 1 | 15584 | 2.000980 | 2.000000 | +0.0490 % | +48.94 % |
| 2 | 15584 | 1.000490 | 1.000000 | +0.0490 % | +48.94 % |
| 10 | 15584 | 0.200098 | 0.200000 | +0.0490 % | +48.94 % |

The relative error is *identical* across five decades of radius (spread below `1e-9`):
it is a property of the triangulation in units of `R`, not of the absolute size. The
pointwise median is likewise radius independent at ~49 %.

**Conditioning.** `taubin_smooth(10) + project_to_field(2)` moves the pointwise median
from ~50 % to
+1.94 %
at level 48, and leaves the integral at
+0.0001 %. Conditioning
buys the pointwise value one order of magnitude and still leaves it at the 2 % boundary --
so "smooth the mesh first" is not a way to reach the criterion pointwise.

## 6. Failure modes and artifacts

* **Pointwise on a sliver-ridden mesh fails the criterion and does not converge** (50 %
  median, p95 several hundred percent). The failure is real and is reported as diagnostic
  G7, not gated away: a curvature-driven mesher must not read raw marching-tet vertex
  curvature.
* **The integral reading is nearly free.** It is inside 0.25 % at every level tested and
  0.03 % from level 32, which
  means an integral quantity (bending energy, hydrostatic moment of curvature) can be built
  on the raw mesh; it also means this case cannot detect a connectivity bug through the
  integral alone.
* **The `2 %` criterion is met by the integral, not by the pointwise value.** Any downstream
  claim must say which reading it means; conflating them is a factor of 50.
* **Convention trap.** `vertex_mean_curvature` returns `H = (kappa_1+kappa_2)/2`, so a
  direct comparison against `2/R` is off by exactly 2. The runner multiplies by 2 explicitly
  (`S = 2H`) and documents it; a factor-2 "error" here would have looked like 50 %.
* **Ecosystem friction (not fixed here).** `tests/test_curvature.py::
  test_the_curvature_module_needs_no_third_party_stack` fails on this checkout because
  `tensorcad/__init__.py` imports `mesh_formats`, which imports NumPy, so the
  standard-library-only promise of `tensorcad.curvature` cannot be exercised through the
  package. Pre-existing, unrelated to this case, recorded rather than patched.
* **Artifacts.** `result.json` (all rows, both series, environment, source hashes),
  `gates_verified.json` (verdict per gate), `figures/convergence.png`,
  `figures/readings.png`, `figures/error_cloud.png`.
* **Reproduce.** `python benchmarks/verified/sphere_mean_curvature/run_benchmark.py`
  (41 s, exits non-zero on a broken gate);
  `python -m pytest tests/test_case_cad_curvature_sphere.py -q` pins the same claims in 14 s.
