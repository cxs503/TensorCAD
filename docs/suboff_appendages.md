# SUBOFF appendages (AFF-8)

The sail and the four stern fins on top of the bare hull from
[docs/suboff.md](suboff.md). Both are published closed-form predicates in the DARPA
report, so neither needs a boolean engine or a CAD kernel: the sail is three
polynomial segments with an elliptical cap, and each fin is a swept NACA section on a
trapezoidal planform. `tensorcad.suboff_appendages` implements those formulas
directly, in the same feet they are published in, and the results are checked point by
point against TensorLBM's torch predicates.

Everything on this page is measured by `python examples/suboff_appendages/run_example.py`,
which exits non-zero if any claim stops holding. `tests/test_suboff_appendages.py`
pins the same numbers.

## Why this is not "the boolean work"

The union of the hull, the sail and the fins *is* three solids overlapping, and the
obvious implementation is a CSG tree. It is also unnecessary: each piece is a
predicate `f(x) < 0`, and the union of predicates is their minimum. The only thing a
boolean engine would add is a B-rep with an explicit intersection curve, which is a
different deliverable (STEP, [issue #4](https://github.com/cxs503/TensorCAD/issues/4)).
What matters for a solver-facing surface is the sign and the zero set, and `min` gives
both exactly. The one trade-off is recorded in the field's `quantity`:
`implicit_function`, not `signed_distance`, because a min-union is a *lower* bound on
the distance in the concave creases where an appendage meets the hull.

## Where the geometry comes from

| Item | Value |
|---|---|
| Source | Groves, Huang & Chang (1989), *Geometric Characteristics of DARPA SUBOFF Models*, DTRC/SHD-1298-01 |
| Transcription | TensorLBM `src/tensorlbm/suboff_cad.py` (`_SAIL_*`, `_FIN_*`, `_NACA_COEFFS`, `suboff_sail_contains_points`, `suboff_fins_contain_points`), read 2026-10-24 |
| Units in the source | feet: hull `L = 14.291667`, `R_max = 0.8333333`; sail `Zmax = 0.109375` half width, `YTMP = 1.507813` body top; fins `H = 13.146284` trailing edge, `r ∈ [0.075, 0.825]` |
| Frame | identical to the hull: `x` bow → stern, `y` transverse, `z` up; `units = "m"`, `coordinate_frame = "object_local"` |
| Scaling | every length is converted with `feet = coordinate · 14.291667 / length`, so the appendages scale with the hull and stay attached |

The sail's planform, in feet of hull length:

| Segment | Interval (ft) | Half width `a(x)` |
|---|---|---|
| entrance | `(3.032986, 3.358507)` | `Zmax·√(2.094759A + 0.2071781B + C)`, `d = 3.0720(x − 3.032986)`, `A = 2d(d−1)⁴`, `B = d²(d−1)³/3`, `C = 1 − (d−1)⁴(4d+1)` |
| middle | `(3.358507, 3.559028]` | `Zmax` |
| exit | `(3.559028, 4.241319]` | `Zmax(2.238361 e f⁴ + 3.106529 e² f³ + 1 − f⁴(4e+1))`, `e = (4.241319 − x)/0.6822917`, `f = e − 1` |

The cross-section is a rectangle of half width `a` up to `z = YTMP` (measured from the
hull axis, so the lower part is buried in the hull) topped by a half-ellipse of
semi-axes `(a, a/2)`. Both pieces are one radial coordinate,

    ρ(x, y, z) = |y|                            z ≤ YTMP
               = √(y² + (2(z − YTMP))²)         z > YTMP

so the sail is `{0 < z, 3.032986 < x < 4.241319, ρ < a(x)}` — exactly the report's
predicate, and no cap-and-body union is needed.

Each fin is a swept NACA 4-digit section: chord `c(r) = −0.466308 r + 0.88859` at span
radius `r`, half thickness `h(s) = 0.2969√s − 0.126s − 0.3516s² + 0.2852s³ − 0.1045s⁴`
at chord position `s = (x − 13.146284)/c(r) + 1`. Because `c` is linear in `r`, the
planform is a convex quadrilateral in `(x, r)`; the four blades are the same predicate
applied twice, once with `(r, t) = (|y|, |z|)` and once with `(|z|, |y|)`.

## Section accuracy against the source

| Quantity | Stations | Max `|error|` vs TensorLBM |
|---|---|---|
| sail half width `a(x)` | 7 (3.05 … 4.20 ft) | `0.00e+00 ft` (0.0000 ppm of `Zmax`) |
| NACA half thickness `h(s)` | 7 (s = 0.02 … 0.98) | `6.9e-18 ft` |

The sail's two segment joints are C0 to `0.00e+00 ft` at the 1e-9 ft probe step, and
both tips close exactly (`a = 0`). The entrance tip is a genuine knife edge — `a` goes
as `√(x − 3.032986)`, so `da/dx` is unbounded there (the finite value at x = 3.10 ft is
`+0.426`). That is a property of the published section, not a defect of this
implementation; it is also where the faint smearing of the reported magnitude below
comes from.

## Containment: the published predicate, point by point

18 points in the DARPA frame were labelled by the **actual torch predicates**
(`suboff_sail_contains_points`, `suboff_fins_contains_points`) and are asserted here:

| Points | Sail labels | Fin labels |
|---|---|---|
| 18 (incl. 0.03 ft above the cap, 3.1e-2 ft forward of a leading edge, the axis column) | 18/18 | 18/18 |

Plus 20,000 random points in `[0.70, 4.45] × [−0.45, 0.45] × [−0.05, 0.60] m` against a
verbatim Python transcription of both predicates, at `scale = 1.0`, `0.85` and `1.2`:
**0 disagreements** (60,000 comparisons per field, 180,000 including the union). The
field's sign agrees with `contains` on every one of them, so a sampler and a predicate
cannot drift apart.

## The magnitude, measured rather than claimed

The sail's cap and the fin's thickness surface have no closed-form nearest point, so
the reported value is the first-order normalised gap `g/|∇g|`: exact on the flat faces,
first-order near them. Measured by walking off the analytic surface along its normal
(1,092 sail normals, 560 fin normals) and comparing `f` with the true offset:

| Field | Offset | Median error | p95 | Max |
|---|---|---|---|---|
| sail | 0.1 mm | 0.00005 mm | 0.00021 mm | 0.00100 mm |
| sail | 0.5 mm | 0.00116 mm | 0.00510 mm | 0.02054 mm |
| sail | 2.0 mm | 0.01575 mm | 0.06718 mm | 0.18125 mm |
| fins | 0.1 mm | 0.00000 mm | 0.00001 mm | 0.00004 mm |
| fins | 0.5 mm | 0.00001 mm | 0.00033 mm | 0.00107 mm |
| fins | 2.0 mm | 0.00018 mm | 0.00553 mm | 0.01839 mm |

`|∇f|` is `p5 … median … max` = `0.9801 … 0.9956 … 1.0000` for the sail and
`0.9987 … 1.0000 … 1.0000` for the fins, i.e. it is bounded above by 1 by
construction and centres on it. On faces that really are planes it is exact, not
first-order: the sail's midbody side is `±a` to `1e-12` relative, and the fin's
trailing edge and tip cylinder likewise.

Where it is weakest: the sail's knife-edge tips. There `da/dx → ∞`, the normalised
gap under-reports, and a sampler sees a slightly rounded tip; the *sign* stays exact,
which is what decides the topology of the extracted surface.

## Volumes: quadrature over the sections, not the sampler

| Quantity | Value |
|---|---|
| sail, total (rectangle + half ellipse integrated over x) | 0.00916043 m³ |
| sail, outside the hull | 0.00423434 m³ (46.2 % of it) |
| fins, total (`2h(s)` over four trapezoids) | 0.00785718 m³ |
| fins, outside the hull | 0.00409551 m³ |
| fin planform area, both sides, four blades | 0.18917 m² |
| bare hull (AFF-1, from `suboff_reference_geometry`) | 0.69921079 m³ |
| **full model = hull + the parts that stick out** | **0.70754064 m³** (+1.191 %) |

That last number is the acceptance reference for the exported mesh: the union adds
1.19 % of volume, and the sampler is judged against the formulas rather than against
its own earlier output.

## Exports

`margin=2`, anisotropic since the model is 8.6:1 with a 0.033 m half-width sail:

| Model | Cells | Triangles | Volume error | Min angle | Slivers < 2° | Euler | Time |
|---|---|---|---|---|---|---|---|
| sail | `(72, 12, 92)` | 121,176 | −0.169 % | 0.109° | 2.31 % | 2 | 2.6 s |
| sail | `(96, 16, 120)` | 212,128 | −0.103 % | 0.179° | 0.40 % | 2 | 4.3 s |
| fins | `(110, 64, 64)` | 191,912 | −0.368 % | 0.005° | 3.21 % | **0** | 4.0 s |
| full | `(218, 26, 38)` | 140,752 | −0.160 % | 0.006° | 2.40 % | 2 | 4.5 s |
| full | `(436, 52, 76)` | 569,056 | −0.059 % | 0.005° | 1.46 % | 2 | 23.6 s |

All are watertight and edge-manifold, with no hard validation errors. The full model's
fine mesh carries 284,530 vertices, a signed volume of 0.707120 m³ (−0.059 % against
0.70754064) and a surface area of 6.350646 m² (the hull alone is 5.988264 m²).

### The fins' Euler characteristic is 0, and that is correct

A bare hull gives `V − E + F = 2`. The fin set alone gives **0**, and the reason is
geometric: the four blades start at `r = 0.075 ft` (not at the axis) and are up to
`0.13 ft` thick, so opposite blades touch each other near the root, while the column
`r < 0.075 ft` at the axis stays empty. The result is one closed shell with a hole
through it — a torus, genus 1, `V − E + F = 2 − 2g = 0`. On the real boat that column
is filled with the hull, which is why the full model is genus 0 again. A test asserts
both values, and `run_example.py` asserts them with that explanation attached rather
than "fixing" the fins into something the report does not describe.

### Where the slivers are

The fin and full-model sliver shares (1.5–3.2 %) live on the NACA leading and trailing
edges: a `√s` nose and a 4th-power tail are knife edges, and marching tetrahedra cuts
them into needles no matter how fine the grid is (3.21 % at both 110 and 160 cells per
axis for the fins). The sail's share, by contrast, halves per resolution step
(2.31 % → 0.40 %) because its only thin feature is the elliptical cap. No smoothing or
decimation exists yet; the numbers are reported as warnings, and removing them is the
next piece of work rather than something this page claims.

## Throughput

| Query | Rate |
|---|---|
| full model, 20,000-point batch | 77,162 queries/s (12.96 µs per point) |

The full model is `min(hull, min(sail, fins))`: the hull pays for its golden-section
nearest-point search and the fins for their four-edge planform distance, so this is
~1.6× the cost of the bare hull alone (124 k queries/s).

## Handoff record

Same schema as every other model (`tensorcad.geometry.handoff/1`), so TensorFEM reads
one shape of record:

```json
{"schema": "tensorcad.geometry.handoff/1", "units": "m", "coordinate_frame": "object_local",
 "vertices": 284530, "triangles": 569056,
 "watertight": true, "edge_manifold": true, "euler_characteristic": 2,
 "signed_volume_m3": 0.707120, "surface_area_m2": 6.350646,
 "bounds_min_m": [0.00006, -0.25396, -0.25400], "bounds_max_m": [4.35472, 0.25396, 0.47625]}
```

`field.metadata()` for `SuboffFullField` nests all three pieces' parameters and their
provenance (`pieces: [hull, appendages]`), because a reviewer asking "which sail
formula" should not have to read this page.

## Not in this scope

- The AFF-3 configuration (hull + sail, no fins) exists as a `SuboffSailField` union
  you can write yourself; it is not a separate class, because it is not a different
  formula.
- The AFF-2 / AFF-4 / AFF-7 variants differ in appendage *layout* (stern planes,
  fairwater planes, flow straighteners), not in this section's formulas. Adding them
  means transcribing their predicates, not extending this abstraction.
- Explicit intersection curves, and therefore STEP: `min` of fields is a valid
  boundary representation only as an implicit surface. Trimming it into a B-rep is
  the B-rep/STEP work.
- Fairing the NACA knife edges into something a finite-element mesh likes: that is the
  smoothing/decimation step, and it belongs on the extracted mesh, not in the published
  field.

## Reproduction

```bash
python examples/suboff_appendages/run_example.py                  # all checks, ~36 s
python examples/suboff_appendages/run_example.py --output results/aff8   # + OBJ/STL/JSON
python -m pytest -q tests/test_suboff_appendages.py                # the pinned numbers
```

The torch comparison that produced the 18 labels and the 0-mismatch counts is not part
of the test suite (it needs TensorLBM); the labels and the verbatim oracle are checked
in instead, so CI verifies this module against the same evidence without the
dependency.
