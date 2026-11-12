# Multi-patch assembly: stitching a NURBS patch network into one conforming shell

**Module** `tensorcad/patches.py` · **Anchor** `GAP-CAD-017` (proposed: multi-patch
NURBS assembly) · **Tests** `tests/test_patches.py` (29 cases) · **Evidence**
`examples/patch_network/run_example.py` → `examples/patch_network/outputs/{results.json,report.md,figures/*.png}`
(every number below is injected from that JSON by the example; none is typed by hand)
· **Composes with** `tensorcad.edge_alignment` (whose first named limitation this
closes), `tensorcad.nurbs` (`NURBSSurface.tessellate`, `circular_arc`),
`tensorcad.mesh` (topology and integrals).

## 1. Case: the step `edge_alignment` explicitly left out

`tensorcad.edge_alignment` makes two patch edges *compatible*: equal degree, equal
knot vector, therefore identical control points and weights when the two edges
describe the same curve. Its own module docstring names the boundary:

> "Curves only (edges), not whole patches: aligning a patch network means aligning the
> edge curves it shares and using them as the patch data, which is exactly what a
> stitching operation needs; **multi-patch face assembly itself stays out of scope**,
> as it does in `tensorcad.nurbs`."

Real ship geometry is a network of patches: a bulbous bow, a parallel midbody, a
transom, a deck, appended fins. None of them fits in one tensor-product patch, and a
solver needs the assembled shell -- closed, edge-manifold, consistently wound, with a
displaced volume it can trust. This module is that assembly step, and it is written to
keep the two questions separate that a tolerance-based stitcher fuses:

1. **Which nodes are the same node?** Answer: the *declaration*. A `Seam` says that
   side `u1` of patch 3 is identified with side `v0` of patch 7, and in which
   direction. Welding is then an index operation (union-find over node slots), so a
   seam becomes a topological identity by construction.
2. **How far apart did the two patches actually put those nodes?** Answer: the
   *measurement*. `Stitched.seam_report()` reports the largest and mean positional gap
   at paired parameters **before** welding, plus the dihedral angle between the two
   patches' normals at those paired collar points.

Welding on position instead (which is what a single patch's `tessellate` does, keyed on
the exact coordinate tuple) requires the two patches to evaluate their common edge to
the same bits, and independent patches generally do not -- a meridian-parametrised
edge and an arc-parametrised edge reach the same point set at different parameters.
`edge_alignment` already measured that failure mode at 3.9 % of the radius.

## 2. Method

* **Network** = `PatchNetwork(surfaces, seams)`, validated at construction: a side
  claimed by two seams raises (one node cannot be welded to two), a patch index out of
  range raises, a self-seam raises. A pole side may not be claimed at all -- it must be
  left unclaimed so the stitcher can collapse it.
* **Sampling**: one common `(nu, nv)` grid per patch, taken on each patch's own
  parameter domain. A seam that joins a `u`-run side (`v0`/`v1`) to a `v`-run side
  (`u0`/`u1`) needs `nu == nv`; the mismatch is refused with a message naming both
  counts, instead of silently sampling one side more finely.
* **Welding**: union-find over node slots; seams unify paired slots, pole sides unify
  all their slots into one. One vertex per weld class, quads emitted from the unified
  grid with `S_u x S_v` winding (the convention `NURBSSurface.tessellate` uses).
* **`derive_seams`**: for networks whose topology is a consequence of the geometry
  (six box faces, the quadrants of a revolved body) the declarations can be derived
  from coincident side endpoints. It is deliberately conservative: pole sides are
  never candidates, and a side with more than one candidate is left **unmatched**
  rather than guessed, because a wrong weld is worse than a reported open side.

## 3. Results

### 3.1 Box 2 x 3 x 4 m -- six bilinear patches, exact tessellation

A bilinear patch is exact under tessellation, so any error here is a *seam* error, not
a discretisation error. Twelve seams, no unclaimed side, one component:

| quantity | measured | analytic | relative error |
|---|---|---|---|
| volume [m³] | 24 | 24 | 0.0 (to the last bit) |
| area [m²] | 52 | 52 | 0.0 (to the last bit) |
| watertight / edge-manifold / boundary edges | True / True / 0 | -- | -- |
| seam positional gap max [m] | 0.0 | 0 | -- |
| seam dihedral angle | 90.0° at every one of the 12 | 90° (G0 by construction) | -- |

The 90° is not a defect: a box edge is a G0 seam, and the report says so instead of
hiding an honest dihedral behind a "watertight" flag.

### 3.2 Torus R = 2 m, r = 0.5 m -- 16 rational biquadratic patches, 32 seams

Analytic `V = 2 pi^2 R r^2 = 9.869604401 m³`, `A = 4 pi^2 R r = 39.4784176044 m²`.
All 32 seams are tangent seams of one smooth surface, so they must be **G1**:

| samples/patch | triangles | volume [m³] | volume err | area [m²] | area err | observed order |
|---|---|---|---|---|---|---|
| 8 | 2048 | 9.742651548 | -1.2863e-02 | 39.2557065677 | -5.6413e-03 | -- |
| 16 | 8192 | 9.837731861 | -3.2294e-03 | 39.4226065518 | -1.4137e-03 | 1.994 |
| 32 | 32768 | 9.861627839 | -8.0819e-04 | 39.4644564864 | -3.5364e-04 | 1.998 |
| 64 | 131072 | 9.867609733 | -2.0210e-04 | 39.4749268023 | -8.8423e-05 | 2.000 |

measured seam gap max **6.12e-16 m** (round-off), measured seam dihedral max
**1.48e-06 °**, `inset_normals = 0`. Second-order convergence is what a bilinear
triangle discretisation of a smooth surface must show; the observed orders are
1.997 / 2.000 / 2.000 (volume) and 1.997 / 1.999 / 2.000 (area).

**Error judgement** (the explicit question): against the analytic torus the finest
refinement is at -0.0202% volume and -0.0088% area, i.e.
**< 3 %**; the coarsest (8 samples per patch, 2048 triangles) is at
-1.2863% / -0.5641%, also < 3 %. The observed orders are
1.994, 1.998, 2.000 (volume) and
1.997, 1.999, 2.000 (area).

### 3.3 Double cone R = 1 m, half-height h = 1 m -- 8 patches, two apex poles collapsed

Analytic `V = 2 pi R^2 h / 3 = 2.094395102 m³` (h = 1 m), `A = 2 pi R sqrt(R^2 + h^2) = 8.885766 m²`.
The apex sides are claimed by no seam; eight pole sides collapse to two vertices:

| samples/patch | triangles | volume [m³] | volume err | area err | observed order |
|---|---|---|---|---|---|
| 8 | 960 | 2.080881389 | -6.4523e-03 | -4.0309e-03 | -- |
| 16 | 3968 | 2.091010586 | -1.6160e-03 | -1.0099e-03 | 1.997 |
| 32 | 16128 | 2.093548592 | -4.0418e-04 | -2.5260e-04 | 1.999 |
| 64 | 65024 | 2.094183451 | -1.0106e-04 | -6.3160e-05 | 2.000 |

Watertight at every refinement; pole sides collapsed = 2, open sides = 0; measured base-rim dihedral 90.0° (a cone with R = h has
a 45° side slope, so the two cones meet at 90° across the base circle -- again a real
G0 feature of the shape, not an artefact). Finest volume error -0.0101%, coarsest -0.6452%, both **< 3 %**.

### 3.4 Negative controls (all from the same run, `results.json:negative_controls`)

| control | measurement |
|---|---|
| one face moved 1e−3 m | `derive_seams` leaves 8 sides unclaimed (it refuses to guess); with the declarations kept, welding is index-level so the shell is still watertight **and** the report says gap = 9.9999999999989e−04 m, volume error +2.81e−04 |
| one patch removed | 4 open sides, `watertight = False`, 16 boundary edges |
| one seam left unclaimed but bit-identical | 2 open sides reported, mesh still watertight (the position key of `MeshBuilder` welds it) -- so read `open_sides`, not just the watertight flag |
| `reversed` flag omitted on a legitimately reversed seam | measured gap jumps from 0.0 to 1.0 m = the shared edge length, i.e. the direction is declared, never guessed |
| `u1`/`v0` seam with `nu != nv` | refused: `"a u1/v0 seam needs nu == nv (got 4 x 8)"` |

## 4. What this does *not* do

* **No trimming.** Every side is a whole boundary curve. A bulbous-bow outline cut out
  of a hull patch is a trimmed curve, not a patch edge, and is out of scope -- the same
  boundary `tensorcad.nurbs` draws.
* **No smoothness repair.** The output is a mesh: seams are welded, not re-fitted to
  G1. The dihedral measurement tells a caller what smoothness the patches already had.
* **No cap synthesis.** A hole is reported (`open_sides`), not filled; filling a planar
  hole is `tensorcad.csg`/caller territory, a curved one needs trimming.
* **Orientation is inherited** from each patch's `S_u x S_v`; if the patches disagree,
  run `tensorcad.mesh.orient_outward` on the result (the sign of `signed_volume_m3`
  reports it).
