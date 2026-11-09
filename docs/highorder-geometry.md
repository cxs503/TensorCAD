# High-order geometry: TensorCAD as the ecosystem's p-order truth source

**Campaign node** `p127w2-n01-cad-highorder-geometry` · **branch** `agent/m7706`
· **commit** `8081ada51866` · **generated** 2026-11-09T01:39:32Z
· **environment** Python 3.11.6 on loongarch64, numpy 2.4.6

## 1. Case

A solver that consumes a NURBS patch needs the surface metrics — the area element
`|S_u x S_v|`, the unit normal, the principal curvatures — *together with the polynomial
order they are accurate to*. Today those numbers arrive either from a tessellation
(second-order at best, and usually unstated) or from a CAD-kernel call whose accuracy is a
documented promise rather than a measured number.

This module makes the accuracy claim measurable on the two shapes every hull check starts
from:

* a **sphere of radius 1.0** — an exact rational biquadratic patch, so the truth is
  `4 pi R^2 = 12.566370614359` with no fitting error anywhere in the chain;
* an **ellipsoid hull** with semi-axes `a, b, c = 1.0, 0.6, 0.42`, i.e.
  `L = 2a = 2.00`, `B = 2b = 1.20`, draft `T = c = 0.42`.
  Its wetted surface is half the closed ellipsoid area, obtained independently from
  Legendre's triaxial formula with Carlson symmetric elliptic integrals, and its block
  coefficient is analytic: `Cb = V / (L B T) = pi/6 = 0.523598775598`.

| quantity | analytic value | how the truth was obtained |
|---|---|---|
| sphere area (R = 1.0) | 12.566370614359 | 4 pi R^2 |
| sphere volume | 4.188790204786 | 4/3 pi R^3 |
| ellipsoid closed area | 5.508414766836 | Legendre + Carlson `R_F`/`R_D` (no NURBS code) |
| hull wetted surface S | 2.754207383418 | half of the above |
| hull volume | 0.527787565803 | 2/3 pi a b c |
| hull Cb | 0.523598775598 | V / (L B T) = pi/6 |
| hull Cp, Cwp | 0.333333333333, 0.785398163397 | 1/3 and pi/4 |
| waterplane area | 1.884955592154 | pi a b |

**Resolution scan.** Gauss-Legendre orders 2, 3, 4, 6, 8 (reference 32) per direction per
knot span; curved elements at order 4 on a 2x2 subdivision of every span; STL facet counts
32…128 per direction; voxel resolutions 48…128 per direction. All parameters are constants
in `run_benchmark.py` — no random seeds, no tuning after the fact.

## 2. How it was done

1. **Knot-span Gauss-Legendre quadrature** (`tensorcad/highorder/quadrature.py`). The rule
   is applied per knot span, never across one: a degree-`p` patch *is* a Bezier segment on
   each span, so a rule spanning a knot would integrate a piecewise-smooth integrand as if
   it were smooth. Subdivision refines inside a span without touching the geometry.
2. **Exact second derivatives by Bezier extraction** (`tensorcad/highorder/nurbs_geom.py`).
   `NURBSSurface.derivatives(order=2)` differentiates the control net, which is impossible
   in a direction whose interior knot multiplicity equals the degree — the three-arc circle
   of the rational sphere is exactly that case, and the parent implementation raises rather
   than inventing a number. `restrict_patch` inserts knots until the span of interest is a
   single Bezier segment and differentiates that; knot insertion preserves geometry *and*
   parametrisation, so the extracted segment is the same surface and its second derivatives
   are analytic. That is what makes the finite-difference cross-check independent instead of
   circular.
3. **Closed-form curvature truth.** Sphere: `kappa = 1/R` everywhere. Ellipsoid: the shape
   operator is the Hessian of the implicit function restricted to the tangent plane divided
   by `|grad F|`; with an orthonormal tangent basis the two eigenvalues follow from a
   closed-form quadratic — no iteration and no reuse of the NURBS answer it checks.
4. **An independent area integral** (`tensorcad/highorder/elliptic.py`). Complete integrals
   are implemented twice (AGM series and Carlson `R_F`/`R_D`), incomplete ones once
   (Carlson), then Legendre's triaxial formula. The two complete routes differ by
   6.21e-13 in `E` and
   4.44e-16 in `K`; both are also checked
   against direct Gauss-Legendre quadrature of the defining integrals in the test suite.
5. **Element contract** (`tensorcad/highorder/contract.py`). `(x, w, n, |J|, uv, element_id)`
   as float64 arrays with a JSON header and a SHA-256 over the canonical byte layout,
   written as NPZ, HDF5 and JSON. A consumer recomputes area and volume from the arrays
   alone — that recomputation is the operational meaning of "the contract round-trips".
6. **Negative control.** The same machinery applied to the *linearized* patch (the STL a
   solver would actually get today): area and volume must fall back to second order and must
   miss the 1e-3 volume gate. If they did not, the p-order numbers above would be an
   artefact of the harness rather than a property of the method.

Toolchain: pure-standard-library library code (NumPy/h5py/matplotlib only in the benchmark),
pytest for the gates, and the two commands

```
python benchmarks/verified/cad_highorder_geometry/run_benchmark.py
python benchmarks/verified/cad_highorder_geometry/render_report.py
```

## 3. Results

### 3.1 J1 — quadrature order is exactly 2p − 1, and the integrals converge faster

Degree-of-exactness probe on `[-1, 1]`: the `p`-point rule reproduces `int x^(2p-1) dx`
exactly and misses `int x^(2p) dx`, pinning the order from both sides.

| p | required order 2p−1 | measured exactness degree | error at 2p−1 | error at 2p | fitted convergence order on a curved patch | verdict |
|---|---|---|---|---|---|---|
| 1 | 1 | 1 | 0.00e+00 | 6.67e-01 | 19.04 | PASS |
| 2 | 3 | 3 | 0.00e+00 | 1.78e-01 | 12.03 | PASS |
| 3 | 5 | 5 | 0.00e+00 | 4.57e-02 | 11.20 | PASS |
| 4 | 7 | 7 | 0.00e+00 | 1.16e-02 | 10.47 | PASS |

Relative area error against the order-32 reference (Gauss points 2 | 3 | 4 | 6 | 8 per direction):

| p | target order | 2 | 3 | 4 | 6 | 8 |
|---|---|---|
| 1 | 1 | 1.63e-06 | 1.90e-09 | 2.70e-12 | 6.27e-16 | 9.41e-16 |
| 2 | 3 | 4.89e-03 | 5.13e-05 | 7.88e-06 | 3.12e-08 | 1.70e-10 |
| 3 | 5 | 1.18e-02 | 4.36e-04 | 3.55e-05 | 8.17e-08 | 3.19e-09 |
| 4 | 7 | 1.43e-02 | 9.20e-04 | 8.23e-05 | 2.92e-07 | 1.04e-08 |

The convergence is *spectral* (analytic integrand on a span), so the fitted orders sit far
above the polynomial requirement; the requirement is a floor, and the numbers show how much
headroom there is. The error stalls at ~1e-15, where `observed_orders` deliberately returns
`None` instead of manufacturing an order out of round-off.

### 3.2 J2 — normals and principal curvatures

| case | samples | worst curvature error / kappa_max | worst unit-normal deviation | worst finite-difference error |
|---|---|---|---|---|
| sphere R = 1.5 | 25 | 2.24e-08 | 3.06e-15 | - |
| ellipsoid 1.0 x 0.6 x 0.42 | 49 | 1.88e-14 | 1.31e-15 | 8.79e-06 |

The gate is 1e-3. The ellipsoid curvature route lands at 1.88e-14
— round-off — because both sides are exact: one differentiates the extracted Bezier segment,
the other the implicit function. The sphere is slightly worse
(2.24e-08) because `H ± sqrt(H^2 - K)` loses half the
significant digits when the two principal curvatures coincide, which is the correct
diagnosis rather than a defect. The finite-difference route at
8.79e-06 confirms the second-derivative supply chain independently (central
differences with a 1e-5 parameter step, truncation error `O(h^2)`).

Sign convention, stated because it is a classic trap: with the outward normal `S_u x S_v`,
a sphere *curves away* from its normal, so the patch reports `kappa = -1/R` while the
geometric magnitude is `+1/R`. The comparison uses magnitudes and the tests assert the sign
relation explicitly.

### 3.3 J3 — isoparametric curved element volume

| case | order | elements | volume | relative error vs analytic | area relative error |
|---|---|---|---|---|---|
| sphere_r1 | 2 | 24 | 4.187126302102 | 3.97e-04 | 3.97e-04 |
| sphere_r1 | 3 | 24 | 4.188775466372 | 3.52e-06 | 3.52e-06 |
| sphere_r1 | 4 | 24 | 4.188791067546 | 2.06e-07 | 2.06e-07 |
| sphere_r1 | 6 | 24 | 4.188790205132 | 8.24e-11 | 8.24e-11 |
| ellipsoid_1_0.6_0.42 | 2 | 24 | 1.055155828130 | 3.97e-04 | 2.11e-04 |
| ellipsoid_1_0.6_0.42 | 3 | 24 | 1.055571417526 | 3.52e-06 | 5.38e-06 |
| ellipsoid_1_0.6_0.42 | 4 | 24 | 1.055575349022 | 2.06e-07 | 3.90e-07 |
| ellipsoid_1_0.6_0.42 | 6 | 24 | 1.055575131693 | 8.24e-11 | 4.88e-09 |

At the gate's setting (order 4, 2x2 subdivisions per knot span) the worst volume error is
2.06e-07 against the 1e-6 requirement. One element per span is *not* enough
(7.50e-05); the tests state that explicitly rather than letting the gate be met by an
unreported extra refinement.

### 3.4 J4 — hull `Cb` and wetted surface against the elliptic-integral truth

Geometry gate first (L/B/T/Cb mutual check to 1%), then three independent routes:

| route | wetted surface S | error vs truth | Cb | error vs truth |
|---|---|---|---|---|
| NURBS quadrature (order 16 per span) | 2.754207383455 | 0.000000% | 0.523598775598 | 0.000000% |
| STL facets, 32 per direction (2016 triangles) | 2.741558743868 | 0.4592% | - | - |
| STL facets, 64 per direction (8128 triangles) | 2.751038448338 | 0.1151% | - | - |
| STL facets, 96 per direction (18336 triangles) | 2.752798405450 | 0.0512% | - | - |
| STL facets, 128 per direction (32640 triangles) | 2.753414720725 | 0.0288% | - | - |
| voxel grid, 48^3 | - | - | 0.523908546505 | 0.0592% |
| voxel grid, 64^3 | - | - | 0.523781250475 | 0.0349% |
| voxel grid, 96^3 | - | - | 0.523684198614 | 0.0163% |
| voxel grid, 128^3 | - | - | 0.523649558165 | 0.0097% |
| AGM/Carlson elliptic integral (independent truth) | 2.754207383418 | 0.0000% | 0.523598775598 | 0.0000% |

Worst route error is 0.4592%, against the 0.5% allowance; the NURBS route is
at round-off (it is the same rational surface the truth describes), while the STL and voxel
routes converge to it from their own discretizations — the STL facets from below (inscribed),
the voxel partial-fill estimate from above and below alternately.

**Sensitivity: curvature error → wetted-surface error.** A 1% curvature error does *not*
imply a fixed wetted-surface error; the ratio depends on the metric perturbation that
produced it.

| perturbation | max principal curvature change | wetted surface change | ratio (area % / curvature %) |
|---|---|---|---|
| isotropic +1% | 0.990% | 2.010% | 2.03 |
| isotropic -1% | 1.010% | 1.990% | 1.97 |
| beam +1% | 0.237% | 0.685% | 2.89 |
| beam -1% | 0.229% | 0.683% | 2.98 |
| draft +1% | 0.826% | 0.452% | 0.55 |
| draft -1% | 0.828% | 0.450% | 0.54 |

For isotropic scaling the ratio is 2.0 by construction (`kappa ~ 1/R`, `S ~ R^2`); shape
changes at fixed scale give 0.55 (draft) to 2.9 (beam). The measured band is
0.54–2.98, so a curvature tolerance of 1% only bounds the wetted
surface to roughly 0.6–2.9% — which is why this module gates the area itself (J3/J4) instead
of inferring it from curvature.

### 3.5 J5 — the `(x, w, n, |J|)` contract

| format | case | rows | max array difference | bit identical | content hash |
|---|---|---|---|---|---|
| .npz | sphere_r1 | 384 | 0.0e+00 | True | fd006094f9b6de9f |
| .h5 | sphere_r1 | 384 | 0.0e+00 | True | fd006094f9b6de9f |
| .json | sphere_r1 | 384 | 0.0e+00 | True | fd006094f9b6de9f |
| npz_as_written | sphere_r1 | 384 | - | - | fd006094f9b6de9f |
| .npz | ellipsoid_1_0.6_0.42 | 384 | 0.0e+00 | True | 464498bb4db15ea4 |
| .h5 | ellipsoid_1_0.6_0.42 | 384 | 0.0e+00 | True | 464498bb4db15ea4 |
| .json | ellipsoid_1_0.6_0.42 | 384 | 0.0e+00 | True | 464498bb4db15ea4 |
| npz_as_written | ellipsoid_1_0.6_0.42 | 384 | - | - | 464498bb4db15ea4 |

Every array returns bit-identical and the recomputed area/volume differ by
0.0e+00 — the campaign allows 1e-14, the format delivers 0.
The NPZ writer is deliberately uncompressed so the file bytes (and therefore the file
SHA-256) are reproducible across zlib versions; the recorded file digests are in
`result.json` under `gates/J5/reports`.

### 3.6 J6 — negative control

| facets per direction | triangles | area relative error | volume relative error |
|---|---|---|---|
| 8 | 112 | 6.925e-02 | 1.355e-01 |
| 16 | 480 | 1.783e-02 | 3.547e-02 |
| 32 | 1984 | 4.493e-03 | 8.975e-03 |
| 64 | 8064 | 1.126e-03 | 2.250e-03 |

Fitted order 1.982 (area) and 1.972
(volume): the linearized patch is second order, so it cannot meet `2p-1 = 3` for `p = 2`,
and at 8 facets per direction the volume error is
1.355e-01 — three orders of magnitude above the J3
gate. **The negative control is red, therefore the positive results are informative.**

### 3.7 J7 — CI budget

`pytest tests/test_highorder_geometry.py` took **5.8 s** wall time
(44 passed, 1 warning in 4.65s) against the 120 s budget; the whole benchmark, including all three figures
and the artifact writes, took 19.7 s.

| gate | verdict | criterion |
|---|---|---|
| J1 | PASS | exactness order = 2p-1 (exact there, wrong one degree beyond) and fitted quadrature order >= 2p-1 for p = 1..4 |
| J2 | PASS | L_inf curvature error / kappa_max <= 1e-3 and unit-normal deviation <= 1e-3, against the closed-form ellipsoid/sphere and against a finite-difference route |
| J3 | PASS | order-4 curved elements on 2x2 subdivisions per knot span: volume relative error <= 1e-6 against 4/3 pi R^3 and 4/3 pi a b c |
| J4 | PASS | wetted surface S and block coefficient Cb from NURBS quadrature, STL tessellation and voxel grid all within 0.5% of the independent AGM/Carlson elliptic-integral value |
| J5 | PASS | write-read-recompute is bit identical (recomputed integrals agree to <= 1e-14) and the content hash survives NPZ, HDF5 and JSON |
| J6 | PASS | faceted patch: observed order <= 2 for area and volume, and volume relative error > 1e-3 at 8 and 16 facets per direction (i.e. J1/J3 cannot be met by facets) |
| J7 | PASS | pytest tests/test_highorder_geometry.py wall time <= 120 s and all green |

## 4. Figures

![convergence.png](../benchmarks/verified/cad_highorder_geometry/figures/convergence.png)

*Left: Gauss-Legendre area error against points per direction for degree-1..4 Bezier patches (spectral decay; fitted order in the legend). Right: the same patch linearized -- second order only, and above the 1e-3 volume gate at every affordable facet count.*

![curvature_error.png](../benchmarks/verified/cad_highorder_geometry/figures/curvature_error.png)

*Curvature and unit-normal error over the parameter square of the ellipsoid, both against the closed-form implicit-surface value; several orders of magnitude below the 1e-3 gate away from the poles.*

![cross_validation.png](../benchmarks/verified/cad_highorder_geometry/figures/cross_validation.png)

*Left: STL wetted surface and voxel Cb against the AGM/Carlson elliptic-integral truth (the NURBS route sits at round-off). Right: wetted-area change against curvature change for isotropic and shape perturbations of the same hull -- the ratio is metric dependent, so a curvature guarantee does not replace an area gate.*

## 5. Notes, traps and limitations

1. **Knot insertion is where the bug was.** The first version of `insert_knot` appended the
   tail as `points[span+1:]` instead of `points[span:]`; with existing-knot splits (where the
   multiplicity is already equal to the degree) nothing was inserted, so the defect stayed
   invisible until a test inserted a *new* knot value and compared knot counts. Lesson: test
   the primitive on a value that is not already a knot.
2. **A C0 direction cannot be differentiated twice — and must not be.** The rational sphere's
   three-arc circle has interior multiplicity equal to the degree, so
   `NURBSSurface.derivatives(order=2)` raises. Bezier extraction is the fix; returning zero
   or a finite difference there would silently poison every curvature downstream.
3. **Curvature sign is a convention, not a number.** `kappa` is signed with respect to
   `S_u x S_v`; a convex body with an outward normal reports negative principal curvatures.
   Anything comparing curvature to a textbook value must state which normal it used.
4. **Spectral convergence beats the stated order — and hides failures.** Because the
   integrand is analytic on a span, the error drops by orders of magnitude per refinement; a
   bug that slows convergence to, say, order 4 would still pass `2p-1` for `p <= 2`. The
   negative control (J6) and the exactness probe (J1a) are what keep the gate honest from the
   other side.
5. **`|J|` is smooth, the parametrisation is not.** On the ellipsoid the area element is
   regular but the parametrisation has a seam at `u = 0/1` and two poles; curvature is only
   evaluated at interior parameters, and every quadrature is per span. Pole-adjacent elements
   are legitimate for integration (they carry the pole's share of area) but not for pointwise
   curvature.
6. **The voxel route is a partial-fill estimate, not an SDF integration.** `VoxelGrid.volume`
   uses the plane-cut rule with the implicit function treated as a length; the hard occupancy
   count is reported next to it and is about an order of magnitude coarser. The 0.5% claim
   rests on the partial-fill number, and the raw resolution series is in `result.json`.
7. **Scope.** Axis-aligned tensor-product patches only: no trimmed surfaces, no multi-patch
   stitching, no boolean geometry. A real hull needs multi-patch topology; this module gives
   the per-patch truth source that such a hull would be assembled from, not the assembly.
8. **`1e-6` volume is not a free lunch.** It needs order-4 elements *and* a 2x2 subdivision of
   each span; the test suite asserts that one element per span fails the gate, so the number
   cannot be quoted without the resolution that produced it.

## 6. Error determination against reference values

| quantity | measured | reference | error | within 3%? |
|---|---|---|---|---|
| sphere area, order-8 per span | 12.566370598413 | 12.566370614359 | 1.27e-09 | **yes** |
| sphere volume, order-4 elements, 2x2 | 4.188791067546 | 4.188790204786 | 2.06e-07 | **yes** |
| ellipsoid volume, order-4 elements, 2x2 | 1.055575349022 | 1.055575131606 | 2.06e-07 | **yes** |
| hull wetted surface, NURBS | 2.754207383455 | 2.754207383418 | 1.33e-11 | **yes** |
| hull wetted surface, STL 128/dir | 2.753414720725 | 2.754207383418 | 2.88e-04 | **yes** |
| hull Cb, voxel 128^3 | 0.523649558165 | 0.523598775598 | 9.70e-05 | **yes** |
| hull Cb, NURBS | 0.523598775598 | 0.523598775598 | 1.06e-15 | **yes** |
| curvature, ellipsoid | (49 samples) | closed-form Hessian | 1.88e-14 of kappa_max | **yes** |
| curvature, sphere | (25 samples) | 1/R | 2.24e-08 of kappa_max | **yes** |
| faceted volume, 8 facets/dir (negative control) | 4.756287 | 4.188790204786 | 1.355e-01 | **no — by design, and that is the point** |

Every quantity that the campaign gates is **below 3%**, and every one of them is in fact
below 0.05% except the deliberately linearized control. The one honest caveat: the 3%
question is not answerable for the *curvature-to-area* mapping, because no fixed ratio
exists (section 3.4) — which is precisely why the area is gated directly.

## 7. Reproducing

```
python benchmarks/verified/cad_highorder_geometry/run_benchmark.py   # gates + figures + JSON
python benchmarks/verified/cad_highorder_geometry/render_report.py   # this document
pytest tests/test_highorder_geometry.py -q                           # 41 gates, ~4 s
```

Provenance: module digest `3865fd3412b56456`, benchmark script digest
`1c42e62f0214bcc3`, test file digest `386e48090fb44c47`,
library source dirty at generation: `no`.
