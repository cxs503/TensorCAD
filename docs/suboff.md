# SUBOFF bare hull

The first non-trivial shape on the Continuous Geometry Core: a real submarine hull
form with a published offset table, implemented as an analytic signed-distance field
so it can be validated, edited, meshed and handed to a solver exactly like the sphere
and the cylinder.

Everything here is measured by `python examples/suboff/run_example.py`, which exits
non-zero if any of its claims stops holding. `tests/test_suboff.py` pins the same
numbers, plus the profile against its external source.

## Where the geometry comes from

| Item | Value |
|---|---|
| Source | Groves, Huang & Chang (1989), *Geometric Characteristics of DARPA SUBOFF Models*, DTRC/SHD-1298-01 |
| Transcription | TensorLBM `src/tensorlbm/suboff_cad.py::suboff_radius_profile`, read 2026-10-24 |
| Units in the source | feet: `L = 14.291667 ft`, `R_max = 0.8333333 ft`; `1 ft = 0.3048 m` exactly |
| Metres | `L = 4.356100 m`, `R = 0.254000 m`, `L/D = 8.575` |
| Segments | bow `[0, 3.333333)`, midbody `[3.333333, 10.645833]`, stern taper `(10.645833, 13.979167]`, stern cap `(13.979167, L]` (ft) |
| Frame | `x` from bow tip (`x = 0`) to stern tip (`x = length`); `y`, `z` from the centreline; `units = "m"`, `coordinate_frame = "object_local"` |
| Claimed accuracy | `|s|` equals the distance to the meridian to `2.4e-8 m` (0.09 ppm of R) |

Those five rows are the "units and provenance" item of the acceptance list, and they
travel in `field.parameters()` / `field.provenance()` rather than living in this page.

The profile agrees with TensorLBM's transcription **exactly** at every station tested
(`max |difference| = 0.00e+00` over 6 stations including 0.05, 0.2, 0.5, 0.7449, 0.9,
0.99), and the vectorised sampling path agrees with the scalar documentation path to
`1.1e-16`.

### One real difference between the two implementations

TensorLBM's `SuboffConfig.r_over_l` defaults to `1/(2·8.57) = 0.058343`, i.e.
`R = 0.254148 m`; the offset table itself gives `R/L = 0.8333333/14.291667 = 0.058310`,
i.e. `R = 0.254000 m`. This module uses the table.

The consequence, measured both ways with the *same* quadrature:

| Radius convention | R (m) | Volume (m³) | Wetted area (m²) |
|---|---|---|---|
| offset table (this module) | 0.254000 | 0.699211 | 5.988264 |
| TensorLBM default | 0.254148 | 0.700027 | 5.991830 |

At TensorLBM's radius this implementation reproduces its
`suboff_analytic_wetted_area(4.3561) = 5.991830060618851` to **+0.00000 %** — the
quadrature is the same, so the entire published difference (0.060 % in area,
0.117 % in volume) is the radius default. Against the published experimental wetted
area of 5.998 m², the offset-table value is −0.162 %.

## Field properties

### Continuity at the four segment joints

| Joint | C0 jump | One-sided slope (step 10⁻⁵) |
|---|---|---|
| bow / midbody (ξ = 0.233236) | `0.00e+00` | `+1.3e-07` / `0.0` |
| midbody / stern taper (ξ = 0.744898) | `4.4e-16` | `0.0` / `-4.1e-07` |
| stern taper / stern cap (ξ = 0.978134) | `-1.2e-13` | `-2.8e-03` / `-3.7e-03` |

C0 is exact (the hull closes) and C1 holds (the slopes shrink with the step, so their
limit is zero from both sides). What is *not* continuous is curvature: the taper and
the cap meet tangentially but with different second derivatives. The bow tip is a
genuine sharp point (`ρ ~ ξ^0.476`, so the meridian slope is infinite there) and the
stern tip is an ellipse pole (also infinite slope). Those two tips are where the
meshing cost concentrates — see the sliver column below.

### Distance, gradient, edits

| Property | Measurement |
|---|---|
| `|s|` vs a 400,001-sample brute-force meridian (150 random points) | max error `2.4e-08 m` (0.09 ppm of R) |
| sign vs closed-form containment | 0 disagreements |
| `|∇s|` off the medial axis | `1.000000` |
| analytic gradient vs central differences | `1.7e-08` |
| throughput, 20,000-point batch | 124 k queries/s (8.1 µs per point) |
| `with_state(length=2L)` | midbody barrel distance unchanged bit-for-bit; stern tip 4.356 → 8.712 m; R unchanged |
| `with_state(radius=1.2R)` | a point at 1.1 R flips `+0.0254 m → -0.0254 m`; L unchanged |

The gradient is exact rather than a finite difference: for an axisymmetric solid the
nearest surface point lies in the query's meridian half-plane, so the field is
`±` the distance to the meridian curve, and the gradient is the unit vector along
`(p − q)` outside and `(q − p)` inside, with `q` the nearest meridian point.

## Exports

`to_mesh(field, cells)` accepts three cell counts, which is what a slender body needs.
At `margin=2`:

| Cells | Triangles | Volume error | Wetted-area error | Min angle | Slivers < 2° | Time |
|---|---|---|---|---|---|---|
| `(109, 13, 13)` | 31,304 | −0.567 % | −0.367 % | 0.038° | 0.76 % | 0.7 s |
| `(218, 26, 26)` | 131,700 | −0.138 % | −0.066 % | 0.032° | 4.32 % | 3.6 s |

Both are watertight, edge-manifold, `V − E + F = 2`, with no hard validation errors.
The `(218, 26, 26)` mesh carries 65,852 vertices, a signed volume of 0.698248 m³ and
a surface area of 5.984292 m².

Why anisotropic resolution: an isotropic `64³` grid (262,144 cells — 78 % more than
`(218, 26, 26)`) reaches −0.089 % / −0.080 % but has only 64 stations along a 4.36 m
hull and produces 10.0 % slivers; the anisotropic grid reaches −0.138 % / −0.066 %
with 4.3 % slivers at 147,368 cells. The cube spends its cells in empty water.

### What the mesh does not reach

Both ends are points of zero area, so the mesh stops short of them by about the local
resolution: at `(109, 13, 13)` the bow tip is within 0.07 axial cells and the stern tip
0.068 m short — 1.7 *transverse* cells, because the stern cap closes faster than the
cross-section grid can follow. The volume and wetted-area errors stay under 0.6 %
because there is almost no volume there. This is reported rather than smoothed away:
`handoff["bounds_min_m"]`, `bounds_max_m` and the volume are what the solver sees.

The mesh-quality checks (smallest angle 0.03°, worst aspect ratio ≫ 25, sliver share
4.3 %) are **warnings**: the solid is valid, the discretisation is merely not
solver-grade. No smoothing or decimation exists yet, and the warning says so instead
of hiding it.

## Handoff record

`to_mesh(...).handoff()` produces the record TensorFEM will read:

```json
{"schema": "tensorcad.geometry.handoff/1", "mesh_sha256": "0a9e6493...",
 "units": "m", "coordinate_frame": "object_local",
 "vertices": 65852, "triangles": 131700,
 "watertight": true, "edge_manifold": true, "euler_characteristic": 2,
 "signed_volume_m3": 0.698248, "surface_area_m2": 5.984292,
 "bounds_min_m": [0.0, -0.254, -0.254], "bounds_max_m": [4.3559, 0.254, 0.254],
 "validation": {...}, "quality": {...}}
```

`field.metadata()` carries the geometry provenance (source, transcription, units,
frame, profile description); `field.json` and `handoff.json` are written by
`run_example.py --output DIR` next to the OBJ and STL.

## Not in this scope

- Appendages (sail, rudder/fins) and the AFF-3 / AFF-8 configurations: this is the
  bare hull (AFF-1) only. The Sail is a separate solid and belongs to the boolean
  work that does not exist yet.
- A `Cp`/`Cb`-parameterised hull family: the profile is the published one, editable
  only through `length` and `radius`. Deforming the *lines* (block coefficient,
  entrance angle, transom) is the next geometry step, not this one.
- Resistance or hydrodynamics of any kind. Nothing here predicts flow; it defines
  and verifies the shape.

## Reproduction

```bash
python examples/suboff/run_example.py                 # all checks, ~9 s
python examples/suboff/run_example.py --output results/suboff   # + OBJ/STL/JSON
python -m pytest -q tests/test_suboff.py              # the pinned numbers
```
