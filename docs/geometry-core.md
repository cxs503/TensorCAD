# Continuous Geometry Core (CGC) v0.1

TensorCAD now has two ways to hold a surface: an **analytic field** (a formula) and a
**learned field** (a `.ngo` object). Before 0.4 they had nothing in common — the
learned one could be trained, edited and extracted, the analytic ones did not exist.
CGC v0.1 is the protocol that makes them interchangeable, plus four exact primitives
to serve as the reference every other representation is judged against.

Everything below is a measured statement. The numbers come from
`python examples/geometry_core/run_benchmark.py`, which fails if any of its own
claims stops being true.

## The protocol

`tensorcad.geometry_field.GeometryField` is the contract. A representation
implements six things and inherits the rest:

| Member | Contract |
|---|---|
| `signed_distance(points) -> list[float]` | the field at an `(n, 3)` batch; the surface is where it is zero |
| `gradient(points) -> list[Point]` | analytic gradient, or a documented approximation (a learned field uses autograd) |
| `bounds() -> BoundingBox \| None` | region containing the surface; `None` means unbounded (a plane) |
| `state() -> dict[str, float]` | the editable controls and their current values |
| `with_state(**values) -> GeometryField` | a **new** field with some controls changed; never mutates |
| `parameters() -> dict` | the full construction set, for the reproducibility record |

Inherited, so two representations cannot drift apart:

| Member | Contract |
|---|---|
| `contains(points)` | inside/outside by the declared `inside_sign` convention |
| `numerical_gradient(points, step=None)` | central differences, for fields with no closed form |
| `metadata()` | kind, quantity, sign, units, frame, bounds, state, parameters, provenance |
| `provenance()` | where the field came from (analytic formula, weights hash) |
| `numpy_callable()` | vectorised evaluator for grid sampling, or `None` (pure-Python still works) |
| `__call__` | alias of `signed_distance` |

Declared per class, not guessed:

| Attribute | Meaning |
|---|---|
| `kind` | short name used in reports (`sphere`, `neural_geometry_object`, …) |
| `quantity` | `signed_distance`, `implicit_function`, or `truncated_signed_distance` |
| `inside_sign` | `"negative"` (default) or `"positive"` |
| `surface_isovalue` | 0.0 by default |
| `units`, `coordinate_frame` | carried, never assumed: mm vs m and object-local vs ship-fixed are real bug sources |

`quantity` is the honesty field. An ellipsoid's `sqrt(Σ(x_i/a_i)²) − 1` has an exact
zero set and an exact sign, but its magnitude is *not* a distance — so it declares
`implicit_function`. A learned TSDF saturates outside its truncation band, so it
declares `truncated_signed_distance`. Code that compares field magnitudes is only
correct for the first of those.

## The four analytic primitives

All four are frozen dataclasses that validate in `__post_init__` and raise
`FieldError` (a `ValueError`) — never `AssertionError`, never a NaN that travels
silently.

| Field | Formula | `quantity` | Measured exactness |
|---|---|---|---|
| `PlaneField` | `n̂·x − offset` | `signed_distance` | `0.00e+00` |
| `SphereField` | `‖x − c‖ − r` | `signed_distance` | `0.00e+00` |
| `CylinderField` | capped cylinder distance (any axis) | `signed_distance` | `0.00e+00` |
| `EllipsoidField` | `√Σ((x−c)ᵢ/aᵢ)² − 1` | `implicit_function` | `0.00e+00` on the surface |

"Exactness" is the benchmark recomputing each field independently — for the sphere and
plane from `math.dist`, for the cylinder from its own closed form (`hypot(max(radial,0),
max(axial,0)) + min(max(radial, axial), 0)`) — and comparing every probe. The
difference is zero in double precision, which is the point: a reference that is 0.1 %
off cannot judge a field that is 0.1 % off.

Gradients are analytic and agree with central differences to `2.7e-11 … 2.0e-10`
(sphere/plane, cylinder, ellipsoid). The capped cylinder needs care at the corner
region: there the nearest feature is the rim circle, and the gradient is the
normalised sum of the radial and axial directions; on the barrel it is radial, on a
cap it is axial. A zero-area corner test pins that down.

`surface_point(field, direction)` finds the surface along a ray by bisecting sign
changes. It exists because **representations must be compared by surface position,
not by field value** — an exact SDF and a truncated learned field have unrelated
magnitudes.

## Edits

`with_state` returns a new field; the original is untouched. The acceptance criterion
is that the edit moves what it names and nothing else:

| Edit | Measured effect |
|---|---|
| sphere `radius` 0.5 → 0.65 | surface moves 0.150000 m in all 24 probe directions; centre unchanged |
| cylinder `height` 1.0 → 1.6 | lid distance moves exactly 0.3 m; barrel distance bit-identical |
| ellipsoid `semi_y` 0.3 → 0.45 | `semi_axes == (0.5, 0.45, 0.2)`; other axes untouched |
| plane `offset` 0.0 → 0.4 | normal identical; value at the origin −0.4 |

Unknown controls raise (`FieldError: sphere has no control(s) ['thickness']`) rather
than being silently ignored — a typo in a design script must not become a
silently-unedited model.

The learned field takes the same edits, with no retraining, because the radius is
field *state* rather than mesh geometry: state 0.45 / 0.55 / 0.60 measured 451.5 /
550.4 / 599.9 mm mean surface radius (+0.33 / +0.08 / −0.01 %).

## Export: field → closed mesh → one acceptance pipeline

`to_mesh(field, resolution)` samples the field, builds a watertight mesh, and runs the
*existing* hull pipeline (`inspect_geometry` + a new `inspect_mesh_quality`) over it.
One pipeline for every representation is the whole point.

Measured at 40 cells per axis:

| Field | Triangles | Zero-area faces | Volume error vs closed form | Slivers (<2°) |
|---|---|---|---|---|
| sphere r = 0.5 | 43,200 | 0 | −0.129 % | 2.33 % |
| cylinder r = 0.3, h = 1.0 | 59,920 | 0 | −0.167 % | 2.32 % |
| ellipsoid 0.5/0.3/0.2 | 42,400 | 0 | −0.134 % | 3.11 % |

All three are watertight, edge-manifold, `V − E + F = 2`, and have no hard
validation errors. The volume bias is the expected marching-tetrahedra deficit and is
reported, not asserted away.

### A defect this work found and fixed

A cylinder of radius 0.3 sampled in a ±0.3 box is exactly 20 cells across, so the
level set ran *through* grid vertices. Marching tetrahedra then emitted **8,108
exactly-zero-area triangles** out of 65,144 (12.5 %), and the sliver share was 22.4 %.
A finite-element code is entitled to reject that mesh. `to_mesh` now offsets its
sampling box by `1/φ²` of one cell, which breaks the coincidence without changing
resolution, units or margin: zero-area faces dropped to **0** and the cylinder's
sliver share to 2.3 %. The regression test states the old number so the fix cannot be
quietly undone.

### Handoff record

`SampledField.handoff()` returns a JSON-ready document describing exactly what a
solver is being given: `schema`, `mesh_sha256` (over float64 vertex bits and indices,
so it identifies the mesh in memory), units, coordinate frame, counts, watertightness,
manifoldness, Euler characteristic, signed volume, surface area, bounds, the full
validation report and every quality check. This is the interface TensorFEM will read;
nothing in it is inferred from a filename.

`inspect_mesh_quality` reports three warnings — smallest interior angle, worst aspect
ratio, sliver share — because a closed solid can still be a poor discretisation. They
are warnings on purpose: the geometry is valid, the mesh is merely not yet solver-grade.

## The learned representation, through the same protocol

`field_from_object(path)` opens a published `.ngo` bundle and returns a
`NeuralGeometryField` — the revision-2 sphere-with-dent committed in this repository.

| Measurement | Result |
|---|---|
| surface vs analytic sphere, 23 unedited directions | mean 3.979 mm, worst 9.330 mm |
| dent axis radius (target 440.000 mm) | 441.730 mm, error +1.730 mm |
| gradient inside the truncation band (autograd) | `|∇|` 0.82–1.19, outward cosine ≥ 0.983 |
| gradient beyond the band | outward cosine 0.721 — the truncation caveat, measured |
| analytic SDF gradient | exactly unit length |
| throughput, 20,000 points, Python loop | 170 k queries/s |
| throughput, vectorised NumPy path | 703 k queries/s (4.1×) |
| throughput, learned field (torch) | 145 k queries/s |
| peak memory | 1.2 MiB vectorised, 4.9 MiB learned |

Two conclusions that matter for anyone using these fields:

1. **`|∇f|` is not a distance estimator.** It drifts from 1 (0.65–1.39 on this bundle)
   because the training objective does not enforce the eikonal condition. Distance
   queries must read the field value; surface queries must read the zero set.
2. **A float32 field must not be differentiated by finite differences at small steps.**
   At `step = 10⁻⁶` the central difference of the network is dominated by float32
   rounding (it disagrees with autograd by up to 0.023); at `step = 10⁻⁴` it agrees to
   `3e-4`. The adapter therefore implements `gradient` with torch autograd and keeps
   `numerical_gradient` only as the protocol's fallback, with the step documented.

### One number that looked wrong

The bundle records `dent_peak_measured_m = 0.052690`, while re-measuring the same
weights through the protocol gives 0.058268 m of radial displacement from the analytic
sphere. Both are right: the recorded peak is measured **against the parent revision**
at the same ray directions (that is the point of the transaction gates — measure the
*edit*, not the fit). It implies a parent axis radius of 494.420 mm, i.e. the fitted
parent sphere sat 5.58 mm inside r = 0.5 on that axis, consistent with the fit's own
recorded boundary error (mean 4.09 mm, max 15.77 mm). The benchmark prints both
numbers side by side so nobody has to guess which baseline a figure used.

## Interfaces left open for what comes next

- **SUBOFF** will be another `GeometryField`: `kind = "suboff_hull"`, analytic or
  offset-table based, with `units` and `coordinate_frame` declared (ship-fixed,
  x aft→fore). It inherits validation, export, handoff and the edit tests for free; the
  datum/units provenance has to travel in `parameters()`.
- **TensorFEM** consumes `to_mesh(field, resolution).handoff()`: units, frame, digest,
  topology, measured volume and quality metrics. A B-rep/STEP adapter becomes a third
  implementation of the same six methods.

## Known limitations

- Marching tetrahedra still produces thin triangles where the surface is nearly
  tangent to a cell edge (2–3 % under 2°); `mesh_min_angle` and `mesh_aspect_ratio`
  report it as warnings. Smoothing/decimation is not implemented.
- `PlaneField` is unbounded and therefore cannot be meshed; that is a property of the
  geometry, not a missing feature.
- The learned field's accuracy (3.98 mm mean, ~0.8 % of the radius) is a property of
  the trained bundle, not of the protocol, and is reported per bundle.
- Only four primitives and one learned family exist. Booleans, sweeps and B-rep are
  not part of v0.1.
