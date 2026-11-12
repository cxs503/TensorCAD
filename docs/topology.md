# B-rep cell topology (`tensorcad.topology`)

**Anchor** `GAP-CAD-001` in `docs/oss_inventory_cad/GAP-ANCHORS.md` (F01, P1: "no
face/edge/vertex topology objects, no topological naming, no Euler operators").
**Tests** `tests/test_topology.py` (30 cases). **Module** `tensorcad/topology.py`.

## 1. What problem this solves

`tensorcad` could answer triangle questions (`tensorcad.mesh`: watertight, edge-manifold,
edge--face counts) and implicit-field questions (`tensorcad.geometry_field`). It could not
answer the question a modelling kernel or a mesher actually asks about *faces*:

* how many faces, edges and vertices does this shell have, and which of them touch which;
* which two faces does this edge separate (so a solver can put a boundary condition on
  one side without naming triangles);
* which faces meet at this vertex (so a feature tree can re-attach a fillet);
* and do those names survive a rebuild at a different sampling density.

`docs/shelling.md` names exactly this as the prerequisite for removing a shell's outer
wall and stitching a rim, and `docs/geometry-core.md` uses triangle topology only.

## 2. Method: labels decide faces, geometry decides nothing

A B-rep face is a **declaration**, not a guess. The builder takes the triangle mesh plus
one label per triangle:

```python
complex_ = from_labeled_mesh(mesh, labels)   # labels: one region name per triangle
complex_.counts()                            # {'vertices': 8, 'edges': 12, 'faces': 6}
complex_.consistency_report()                # Euler characteristic through both routes
```

1. faces = connected components of same-label triangles (union-find over shared edges);
2. boundary arcs = mesh edges whose incident triangles lie in different faces (or one);
3. vertices = nodes where the chain structure changes: degree != 2 in the boundary graph,
   or degree 2 with *different* face pairs on its two arcs. A chain that closes with no
   such node gets one deterministic start vertex, because an edge needs endpoints for
   `V - E + F` to mean anything;
4. edges = the chains between those vertices.

Nothing is merged by proximity and nothing is inferred from dihedral angles: a face that
looks flat but was not declared stays two faces, and a declared face that is not connected
is one face with a disconnected triangle set (the Euler check will show it).

## 3. The check that cannot be gamed

`V - E + F` is a topological invariant. So the face/edge/vertex structure and the triangle
mesh it came from must agree on it, even though their counts are completely different
(a box is 8/12/6 as a B-rep and 386/1152/768 as a 8x8 tessellation). `validate()` also
requires every edge to border one or two faces, every face boundary to be a single closed
chain, and every edge to be traversed in opposite directions by its two faces (winding).

## 4. Measured results (from `tests/test_topology.py`, exact integers)

| geometry | V | E | F | χ | mesh χ | closed | genus |
|---|---:|---:|---:|---:|---:|---|---:|
| box 2x3x4 m, six bilinear patches, 8x8 | 8 | 12 | 6 | 2 | 2 | yes | 0 |
| sphere R=12 m, two hemispheres (24x12) | 1 | 1 | 2 | 2 | 2 | yes | 0 |
| hull L=100 m B=18 m T=6 m, shell+deck, 41x13 | 1 | 1 | 2 | 2 | 2 | yes | 0 |

* Box: 6 faces `x0 x1 y0 y1 z0 z1`, 12 edges `E:<face>|<face>`, 8 vertices
  `V:<f>|<f>|<f>`; every edge borders exactly two faces; every vertex is a three-face,
  three-edge corner; face areas match `bc, bc, ac, ac, ab, ab` to 1e-12 relative.
* Sphere: the equator is one **closed chain** (`loop=True`, both endpoints the same
  vertex) between `north` and `south` -- `1 - 1 + 2 = 2` is exactly why the deterministic
  loop vertex exists.
* Hull: `deck` and `shell` are separated by the closed deck edge (`vg` ring), which runs
  through the stem and the post and therefore has no corner anywhere: 1 vertex, 1 edge,
  2 faces. Coarse (21x9) and fine (81x25) grids give identical counts and names.

Naming stability: the same box at `nu = 4, 8, 16` and with the six patches handed in the
order `(3, 0, 5, 1, 4, 2)` produces the identical `stable_names()` triple. Names come from
labels plus scale-relative rounded coordinates, never from indices.

## 5. Euler operators

| operator | effect | invariant |
|---|---|---|
| `split_face` (mef) | +1 face, +1 edge, V unchanged | χ and validity |
| `split_edge` (mev) | +1 vertex, +1 edge | χ and validity |
| `merge_faces` (kef) | -1 face, -1 edge | **not yet verified** -- see below |

Both verified operators refuse inputs that are not what they claim: a split across two
vertices that are not both on the face raises, a second split that would duplicate a name
is uniquified rather than silently overwritten, and every operator re-checks χ before it
returns a result.

**Known defect (honest FAIL, not covered up).** `merge_faces` refuses the box
split/merge round trip: after joining `z1` with its split half, the merged face's edge
set does not chain into a single boundary walk, so `_assert_preserved` raises
`TopologyError` instead of returning a wrong complex. The case is recorded as a strict
`xfail` in `tests/test_topology.py::TestEulerOperators::test_merge_faces_undoes_the_split`,
so it cannot pass silently and cannot be forgotten. The defect is in the chain-joining
step (which of the two faces' walks is reversed), not in the accounting: χ preservation
is checked before validity.

## 6. Scope and limitations

* **Faces are declared, not discovered.** A feature tree knows its faces; this module
  does not fit them to curvature. That is deliberate (a guess is not reviewable) and it
  means an unlabelled triangle is an error rather than a face of its own.
* **Single manifold shell per call.** Non-manifold edges are refused, because a boundary
  chain through them is not well defined. Multiple shells are separate calls.
* **Derived complexes carry no mesh provenance.** After an Euler operator, the structure
  no longer corresponds to the source tessellation, so `mesh_counts` is dropped and the
  χ-vs-mesh check is replaced by the operator's own before/after check.
* **Not closed by this work:** the other half of `GAP-CAD-001` is *topological naming
  stability across feature regeneration* (a fillet surviving a change of `a` in a box).
  Here names are stable under resampling and reordering, which is what the acceptance
  criterion tests; stability under a parameter change that moves geometry needs a named
  feature history, which this repository does not have yet.
* No example directory: every number in this document comes from the test file above,
  which is the artifact.
