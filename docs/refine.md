# Refinement to a target size field

`tensorcad.mesh_conditioning` only coarsens — every step of that pipeline (Taubin
smoothing, projection, quadric decimation) keeps a vertex or removes it. A consumer that
needs a **finer** discretisation than the sampler produced had no route up. This page
documents the route: `tensorcad.refine`, which splits triangles until every edge is at or
below the local target size, conformingly, and measures what it did.

The two claims are deliberately separated, because conflating them is how meshing claims
usually go wrong:

| Claim | Where it comes from | Evidence |
|---|---|---|
| Splitting is **geometry-neutral** | a new vertex is the midpoint of an existing edge, so every sub-triangle lies in its parent's plane | volume/area/centroid relative change ≤ 2.2e-12 % on an 8 648-face sphere, 8 648 → 25 744 faces |
| Resolution alone is **not accuracy** | the faceted hull stays equally wrong at four times the face count | projected vs unprojected columns below: only the projected one moves |
| Projection is what **adds accuracy** | bounded Newton step onto the analytic zero level set, new vertices only | sphere volume error −0.655 % → −0.236 %; Suboff −2.350 % → −1.312 % vs published 0.699211 m³ |
| Splitting is **conforming** | one midpoint per marked edge, both incident faces use it (red-green 0/1/2/3-edge patterns) | watertight, edge-manifold, Euler = 2 before and after; zero zero-area faces |
| The target is **local** | an edge is judged at its own midpoint against `h(x, y, z)` | adaptive Suboff: 27.1 % fewer faces than uniform at the same fine target |

Everything here is produced by
`python examples/refine/run_example.py --output results/refine`, which exits non-zero if a
claim stops holding; `tests/test_refine.py` pins the same behaviour on fast meshes.

## API

```python
from tensorcad import refine
mesh, report = refine(mesh, 0.08)                        # uniform target, 8 cm
mesh, report = refine(mesh, lambda x, y, z: 0.08 if x < 2.2 else 0.5)
mesh, report = refine(mesh, 0.08, project=field)         # field: signed_distance + gradient
report.to_dict()      # JSON-ready: before/after metrics + projection residuals
```

`report.converged` is `False` when the pass budget runs out before every edge meets its
target; `report.after.max_edge` then says what was actually achieved, so a truncated run
cannot be mistaken for a finished one.

## Measured: a sphere, against the closed form

Input: `to_mesh(SphereField(radius=1.0), (18, 18, 18))`, 8 648 triangles, volume error
**−0.6553 %** against 4πR³/3.

| target [m] | faces | max edge [m] | plain ΔV | plain ΔA | plain centroid shift [m] | projected volume error | residual [m] |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.15 | 9 544 | 0.1498 | −2.13e-16 | 0.0 | 7.3e-17 | −0.5913 % | 1.11e-16 |
| 0.10 | 13 520 | 0.1000 | +5.76e-15 | 0.0 | — | −0.4281 % | 2.22e-16 |
| 0.07 | 25 744 | 0.0700 | +2.11e-14 | 0.0 | — | −0.2358 % | 2.22e-16 |

Read the "plain" columns as the invariance claim: the mesh has four times the triangles
and the same displacement to fifteen digits. Read the last column as the accuracy claim:
projecting the *new* vertices onto the sphere halves-to-thirds the faceting error, and
the residual after projection is machine epsilon — the Newton step is exact for an SDF
with a unit gradient, which is what makes this a clean test of the operator rather than of
a tolerance.

## Measured: the DARPA Suboff bare hull, against a third-party number

Input: `to_mesh(SuboffHullField(), (28, 6, 8))`, 4 176 triangles, displacement
**0.682780 m³** vs the published AFF-1 volume **0.699211 m³** (Groves et al. 1989) —
**−2.35 %**.

| target [m] | faces | max edge [m] | displacement [m³] | vs published | residual [m] |
|---:|---:|---:|---:|---:|---:|
| input | 4 176 | 0.2366 | 0.682780 | −2.3500 % | — |
| 0.12 | 5 560 | 0.1199 | 0.685967 | −1.8942 % | 9.5e-16 |
| 0.08 | 10 114 | 0.0800 | 0.690035 | −1.3124 % | 9.8e-16 |

The error falls, but not by 4× per level: the residual is dominated by the flat deck and
the waterline edge of the sampled hull, where the facets are *correct* chords and
projection has nothing to fix. That is why the example reports the number instead of an
order of accuracy. The point of the case is that the operator works on a real hull
geometry with an external anchor, not that it is a second-order surface approximant; for
that see `docs/hull_analytic_reference.md` and the grid-refinement rows in `benchmark.py`.

## Measured: adaptivity is worth what it costs

Same Suboff mesh, target `h = 0.08` forward of midship and `h = 0.50` aft:

| variant | faces | vs uniform |
|---|---:|---:|
| uniform `h = 0.08` | 10 100 | — |
| size field (0.08 forward, 0.50 aft) | 7 358 | **−27.1 %** |

Both runs converge. A 27 % saving on a hull where half the surface wants a 6× coarser
element is a modest but real number, and it is the honest one for this shape: the Suboff
aft body is as curved as its fore body, so a size field that halves resolution aft can
only win where the caller actually wants that resolution.

## Where the split patterns come from

Edges are marked by length. A triangle then has 0, 1, 2 or 3 marked edges and is split by
the classical red-green patterns:

| marked edge(s) | new triangles |
|---|---|
| none | 1 (unchanged) |
| one | 2 (bisect toward the opposite vertex) |
| two sharing a vertex | 3 (corner triangle + split quad) |
| three (red) | 4 (corner triangles + inner triangle) |

Every marked edge contributes exactly one midpoint, computed from the *canonical*
(ordered) vertex pair, so both incident faces look up the same coordinate tuple and the
builder welds them into one vertex. That is the conformity argument, and it is why
`V − E + F` survives: red leaves it unchanged, and each green pattern adds
`ΔV = ΔF` vertices more than faces while adding `ΔE = 2ΔF + ΔV − ΔV = …` edges — the
test checks the invariant on the actual mesh rather than trusting the arithmetic.

## Honest limitations

- **No curvature sensor.** The target is the caller's field. A curvature- or
  wavelength-driven `h` has to be supplied by whoever knows the physics; this module does
  not pretend to infer it.
- **No anisotropic splitting.** Long thin triangles are split isotropically, so a
  boundary-layer style aspect ratio cannot be requested; `max_aspect` in the report shows
  what the input had and what the split did to it (it never improves the worst sliver,
  because a sliver's short edges are not marked).
- **No projection of existing vertices.** They are left exactly where the caller put
  them, which is what makes the "original vertices untouched" test possible and keeps the
  operator composable with `mesh_conditioning`, which *does* move them.
- **Cost is linear in the output.** 8 648 → 25 744 faces on the sphere takes ~1 s of pure
  Python; the projection adds one field evaluation and one gradient per new vertex per
  iteration. `max_faces` (default 4 M) refuses a runaway refinement instead of exhausting
  memory.
- **Not a substitute for re-sampling.** Refining a mesh can never recover a feature the
  sampler missed; `resolve_resolution` and `to_mesh` are the tools for that.