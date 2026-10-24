# Mesh conditioning

Marching tets on an analytic field produce a *closed* surface, which is what a solver
needs, and a *bad* triangulation, which is what a solver complains about: needles,
caps, and slivers along every knife edge. The AFF-8 model is the worst case in this
repository, because the NACA fin sections and the sail's leading tip are genuine sharp
edges — the analytic surface really does have a 0° corner there, and the mesher
faithfully reproduces it.

`tensorcad.mesh_conditioning` repairs the mesh without moving the shape. It never edits
geometry: every step either slides vertices along the surface (Taubin smoothing) or puts
them back onto it (projection), and `condition_mesh` refuses to return a mesh whose
volume or area moved beyond a declared tolerance.

Everything on this page is measured by
`python examples/mesh_conditioning/run_example.py`, which exits non-zero if any claim
stops holding. `tests/test_mesh_conditioning.py` pins the same behaviour.

## The pipeline

```
sample ──► Taubin smoothing ──► projection onto the field ──► decimation ──► handoff
           (band-pass)          (bounded Newton, line search)  (quadric error)
```

| Step | What it does | What it costs |
|---|---|---|
| `taubin_smooth(mesh, iterations)` | `lambda = 0.5` followed by `mu = -0.53`, an alternating band-pass | volume drift vs the analytic reference 0.13 % at 30 sweeps; 0.4 s for 65 k triangles |
| `project_to_field(mesh, field, iterations)` | one Newton step per vertex along the field gradient, step bounded by a fraction of that vertex's own mean edge length, 4-step backtracking line search | 1.8 s per iteration at 65 k vertices (0.8 s for the normals plus the field queries) |
| `decimate(mesh, target_faces=…, max_edge_length=…)` | quadric-error edge collapse with the link condition, orientation, duplicate-face and aspect-ratio guards | 5–12 s for 65 k → 20 k |
| `condition_mesh(...)` | the whole thing, with the tolerances enforced | 35 s at 45 k faces |

The order is a decision, not an accident. Smoothing first is what makes the *input*
dense enough for decimation to be safe: a sliver's quadric is nearly degenerate, so
collapsing in the sliver field is where the worst triangles are chosen. The test suite
covers both orders and the example prints both; the numbers are below.

## Measured: the recommended order

Input: `to_mesh(SuboffFullField(), (150, 18, 26), margin=2)` — 65,812 triangles, 1.664 %
of them slivers (thinner than 2°). Reference: the analytic AFF-8 volume
0.70754064 m³ and wetted area 6.350646 m².

| Stage | tris | min angle | max aspect | slivers | ΔV | ΔA | s |
|---|---|---|---|---|---|---|---|
| input | 65,812 | 0.0080° | 7,126 | 1.664 % | −0.323 % | −0.378 % | — |
| Taubin ×30 | 65,812 | 8.3202° | 8.3 | 0.000 % | −0.127 % | −0.729 % | 0.4 |
| + projection ×30 | 65,812 | 3.6590° | 17.6 | 0.000 % | −0.117 % | −0.198 % | 53.6 |
| + decimate → 45,000 | 45,000 | 2.7038° | 25.0 | 0.000 % | −0.118 % | −0.198 % | 5.4 |
| + decimate → 30,000 | 30,000 | 2.5641° | 25.0 | 0.000 % | −0.136 % | −0.196 % | 9.7 |
| + decimate → 20,000 | 20,000 | 2.4021° | 25.0 | 0.000 % | −0.209 % | −0.162 % | 11.6 |

The minimum angle improves by a factor of 1,000 in the smoothing step, then gives some
of it back to the projection, which is correct: the projection is what makes the surface
*the analytic surface* again, and the analytic surface has sharp edges. The worst
displacement over the whole pipeline is 25.6 mm (smoothing) and 18.4 mm (projection) on
a 4.36 m hull, i.e. 0.6 % and 0.4 % of the length, and the worst residual `|f|` after
projection is 0.57 mm.

Max aspect ratio stops at exactly 25.0 from 45 k faces on. That is not a measurement,
it is the gate: `decimate` may not *create* a triangle worse than `max_aspect_ratio`, so
once the gate is the binding constraint the reported number is its value. Reading it as
"the mesher improved to 25" would be wrong.

## Measured: the other order, and what it costs

`condition_mesh(mesh, field=..., target_faces=...)` decimates first, then smooths and
projects. It is the API the tests exercise directly, and at the same face budget it is
measurably worse:

| Stage | tris | min angle | max aspect | slivers | ΔV | ΔA | s |
|---|---|---|---|---|---|---|---|
| `condition_mesh(45,000)` | 45,000 | 2.2075° | 43.8 | 0.000 % | −0.127 % | −0.170 % | 40.3 |
| `condition_mesh(30,000)` | 30,000 | 0.1487° | 514.8 | 0.083 % | −0.299 % | −0.227 % | 34.7 |

At 30 k faces the decimate-first path leaves 28 slivers and a 0.15° minimum angle,
because collapsing the sliver field first is exactly how a decimator picks the worst
triangles to keep. Both results pass their volume/area gates — the gates are about the
*shape*, and both orders preserve the shape — but only the smooth-first order removes
the slivers at that budget. Larger budgets are kinder: at (218, 26, 38) input
(140,752 triangles) the decimate-first path held 0.466° / aspect 211 / 0.0010 % slivers
at 100 k faces and 0.076° / 1,328 / 0.069 % at 70 k.

## What is gated, and what is only reported

Gated (`condition_mesh` raises `MeshConditioningError`, the example exits non-zero):
watertight, edge-manifold, Euler characteristic unchanged, no zero-area faces, volume
drift within 0.5 %, area drift within 1 %, and zero slivers in the recommended order.

Reported, not gated: the minimum angle and the aspect ratio. A decimator that is allowed
to fold edges cannot promise a minimum angle, and a knife edge's neighbourhood has no
triangle with a good angle to offer — the analytic corner is 0°. Publishing a threshold
for those numbers and then loosening it whenever it failed would make the gate
meaningless. The honest statement is the table above: `> 2.4°` at 20 k faces for the
recommended order, `0.15°` for the decimate-first path at 30 k, and the sliver count,
which is a *topological* defect rather than a taste judgement, gated at zero.

The knife edges are also why the example reports, rather than hides, the residual
quality: the fins' trailing edges and the sail's leading tip are thin by construction.

## Long runs

At the repository's recommended export resolution (218, 26, 38 → 140,752 triangles) the
same pipeline gives, without any decimation: Taubin ×30 plus projection ×30 in 120 s →
min 3.868°, aspect 24, no slivers, V = 0.707089 m³ (−0.064 %), A = 6.348993 m², maximum
displacement 20.4 mm; another 60 sweeps (249 s) reaches 4.056°, aspect 27,
V = 0.707093 m³, A = 6.351633 m². Two hours of smoothing is not worth 0.19° on a hull
whose analytic corners are 0°, which is why the example defaults to the smaller grid.

## What this does not do

- It does not change the geometry. Every vertex stays on (or within a micrometre of) the
  analytic surface; features are neither created nor removed.
- It does not recover an analytic corner. A marching-tet surface at a knife edge has
  triangles with a near-zero angle; smoothing moves them, projection puts them back.
- It does not re-mesh to a target size field. `decimate` only coarsens; going the other
  way needs a refinement step, which this module does not implement.
- It does not touch the topology. The link condition, the orientation check and the
  duplicate-face check together guarantee that every collapse preserves the surface as a
  manifold with the same Euler characteristic, and the gates verify that on the output.

## API

```python
from tensorcad.mesh_conditioning import (
    ConditionedMesh, MeshConditioningError, MeshMeasurements, condition_mesh,
    decimate, measure_mesh, project_to_field, taubin_smooth,
)

result = condition_mesh(
    mesh, field=SuboffFullField(), target_faces=45_000,
    smooth_iterations=20, projection_iterations=2,
    volume_tolerance=0.005, area_tolerance=0.01,
)
result.passed          # gates on topology and shape drift
result.summary()       # {"schema": "tensorcad.mesh.conditioning/1", before/after, moves}
result.mesh            # the surface to hand off
```

`measure_mesh(mesh)` is one pass over the mesh returning
`MeshMeasurements(triangles, watertight, edge_manifold, signed_volume, surface_area,
min_angle_degrees, max_aspect_ratio, sliver_share, zero_area_faces)`, using the same
thresholds as `validation.inspect_mesh_quality`; there is no second, laxer definition of
a good triangle in this package.
