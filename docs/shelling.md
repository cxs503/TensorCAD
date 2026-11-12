# Shelling: hollow a closed solid, and prove the wall is there

`tensorcad.shell` closes `GAP-CAD-011` in `docs/oss_inventory_cad/GAP-ANCHORS.md`, whose
acceptance criterion is explicit: *"a box shelled at thickness `t` matches the analytic
outer-minus-inner volume to 2 %, and a sweep of `t/L` from 0.01 to 0.1 does not fail"*. Before
this module the package could smooth, decimate and project a mesh but not hollow it -- every
operation produced a surface, none produced a **wall**.

```python
from tensorcad.shell import shell, hollow_body
from tensorcad.section_gz import prism_mesh

box = prism_mesh(((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)), 1.0)  # 1 m cube
result = shell(box, 0.1)               # refuses if the wall cannot be achieved
result.shell_volume                    # 0.488 = 1 - 0.8^3, the analytic wall
hollow_body(result).mesh               # outer + reversed inner: watertight export
```

## Method

Every triangle carries a plane `n . x = d` with `n` the outward unit normal, so a solid is
locally the intersection of half-spaces `n . x <= d` and its inward offset is the same solid
with `d` replaced by `d - t`. At a vertex the incident facets hand over their (deduplicated)
planes and the offset vertex is the angle-weighted least-squares point

    minimise  sum_i w_i (n_i . x - (d_i - t))^2 ,     w_i = interior angle of face i at v

solved through the 3x3 normal equations (Cramer, no NumPy: the core stays standard-library
only). For a box corner the three planes are orthogonal and the solution is `d_i - t` along
each with **no** dependence on the weights; on a flat or edge vertex the system is
rank-deficient and a ridge term `1e-14 * trace/3` selects the minimum-norm displacement, which
is exact in the constrained directions because the right-hand side lies in the span of the
normals. The ridge is what the box error of 3.0e-14 is made of.

What the module spends its remaining lines on is **honesty about failure**, because a collided
offset still looks like a mesh:

* the inner signed volume must stay positive (a fully inverted offset does not);
* inner faces must keep the winding of the outer ones (`min dot(n_inner, n_outer)` is reported
  and a non-positive value is refused);
* the achieved wall is *measured*: for every inner vertex the distance to the **whole** outer
  surface is computed through a uniform triangle grid, and a shortfall beyond
  `wall_tolerance` (default 5 %) fails the call. This is the check the plane solve cannot make
  -- two opposite walls of a thin arm, each vertex at `t` from its own faces and closer than
  `t` to the far one.

## Verification

### 1. Box: the acceptance criterion, run literally

| t/L | inner [m^3] | analytic [m^3] | rel. error | wall [m] |
|---|---|---|---|---|
| 0.01 | 0.941191999999972 | 0.941192000000000 | 2.97e-14 | 0.010000000000 |
| 0.02 | 0.884735999999974 | 0.884736000000000 | 2.97e-14 | 0.020000000000 |
| 0.03 | 0.830583999999975 | 0.830584000000000 | 2.99e-14 | 0.030000000000 |
| 0.04 | 0.778687999999977 | 0.778688000000000 | 3.01e-14 | 0.040000000000 |
| 0.05 | 0.728999999999978 | 0.729000000000000 | 3.05e-14 | 0.050000000000 |
| 0.06 | 0.681471999999980 | 0.681472000000000 | 3.00e-14 | 0.060000000000 |
| 0.07 | 0.636055999999981 | 0.636056000000000 | 3.02e-14 | 0.070000000000 |
| 0.08 | 0.592703999999982 | 0.592704000000000 | 2.96e-14 | 0.080000000000 |
| 0.09 | 0.551367999999984 | 0.551368000000000 | 2.98e-14 | 0.090000000000 |
| 0.10 | 0.511999999999985 | 0.512000000000000 | 3.04e-14 | 0.100000000000 |

Worst relative error over the whole sweep: **3.05e-14** against a
2 % threshold -- nine orders inside it, and the wall is the requested
0.100000000000 m at `t/L = 0.10`.

The box alone cannot validate the weighting (three orthogonal planes make every weighting
equal), so the next case uses facets that are **not** mutually orthogonal.

### 2. Pentagon prism vs independent 2-D half-plane clipping

The eroded solid is also the intersection of the offset half-planes. `tests/shell_oracle.py`
computes it the other way round: Sutherland-Hodgman clipping of the section by every edge's
offset half-plane, times the shortened extrusion. It never touches a plane system, a vertex or
a least-squares solve.

| t [m] | inner [m^3] | clipped ref [m^3] | rel. error | corners |
|---|---|---|---|---|
| 0.02 | 1.100325623474289 | 1.100325623474316 | 2.46e-14 | 5 |
| 0.05 | 0.923302866069921 | 0.923302866069944 | 2.49e-14 | 5 |
| 0.10 | 0.671989350204373 | 0.671989350204390 | 2.48e-14 | 5 |
| 0.15 | 0.470598690123310 | 0.470598690123322 | 2.45e-14 | 5 |

Worst relative error: **2.49e-14**. Two routes, one number.

### 3. Sphere: faceting is the only error left

| meridians x parallels | faces | shell [m^3] | analytic [m^3] | rel. error | order |
|---|---|---|---|---|---|
| 12 x 6 | 120 | 1.905196 | 2.044130 | 6.797% |  |
| 24 x 12 | 528 | 2.008892 | 2.044130 | 1.724% | 1.98 |
| 48 x 24 | 2208 | 2.035288 | 2.044130 | 0.433% | 1.99 |

Observed convergence order **1.98 / 1.99**,
i.e. second order, which is the signature of faceting (`O(h^2)`) rather than of a modelling
error, and the achieved wall is the requested `0.2 m` to `1e-12`.

### 4. The wall measurement is falsifiable

Feeding `min_wall_thickness` a *deliberately under-offset* surface must make it report the
smaller distance it was given; that is what makes the number it reports for a correct shell
mean something.

| offset applied [m] | measured wall [m] | rel. error |
|---|---|---|
| 0.020 | 0.020000000000 | 1.02e-14 |
| 0.050 | 0.050000000000 | 9.99e-15 |
| 0.070 | 0.070000000000 | 1.01e-14 |
| 0.100 | 0.100000000000 | 9.99e-15 |

### 5. Failures are first-class

| case | refused | reason |
|---|---|---|
| plate_0.1_shelled_at_0.08 | True | the inward offset of 0.08 m inverts the solid: the inner signed volume is -0.042335999999998736 m^3 (outer 0.09999999999999999 m^3), so the thickness exceeds half the local wall somewhere. Reduce the thickness or split the body |
| taper_shelled_at_0.10 | True | the inward offset folds the surface: min(dot(n_inner, n_outer)) = -1.0 over the faces, so some offset faces are wound backwards |

The plate and the taper are the two genuine geometries the collapse checks exist for: the plate
inverts the inner solid outright, the taper folds faces locally. An open surface, an
inward-oriented surface and a non-positive thickness are refused before any geometry is
attempted.

## The export convention

`hollow_body(result).mesh` is outer faces plus the inner faces **reversed**, so the surface
bounds the material: its divergence-theorem volume is `0.488000000000 m^3` =
material `0.488000000000` = outer 1.0 minus cavity
`0.512000000000`. It is watertight (True) and edge-manifold
(True). A consumer who assumes "both components outward" expects the
*sum* and is wrong by twice the cavity volume; the field names (`material_volume`,
`cavity_volume`) exist to make that impossible to conflate.

## What it does not do

* **No face removal and no rim stitching.** Opening a face and bridging the gap needs B-rep
  topology, recorded as `GAP-CAD-001`; the result here is a closed double-walled body.
* **Curved bodies are offset through their facets**, so the inner surface of a coarse sphere is
  short of the true offset by `O(h^2)` (measured above, and reported as `wall_shortfall`, not
  hidden). The default 5 % wall tolerance is what a mesh has to meet to be accepted; a coarser
  body is refused with both possible causes named.
* **Collision detection is geometric, not topological.** A self-intersection between two
  different patches that neither inverts a face nor shortens a wall is not detected by design;
  `tensorcad.intersect.find_self_intersections` is the separate tool for that, and it is not
  run here because it is `O(n^2)` on the candidate pairs.

## Files

* `tensorcad/shell.py` -- the module (standard library only).
* `tests/test_shell.py`, `tests/shell_oracle.py` -- 48 tests and the independent 2-D oracle.
* `examples/shelling/` -- this report's numbers (`results.json`, `report.md`, `shelling.png`),
  generated programmatically; no number in this document is typed by hand.
