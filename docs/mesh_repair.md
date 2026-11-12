# Mesh repair: hole filling, manifoldisation, winding, and the proof afterwards

`docs/oss_inventory_cad/GAP-ANCHORS.md` records **GAP-CAD-009** as *partial*: the package
could weld, drop degenerate facets and **detect** a defect, but there was no path back to a
mesh a solver can use. The detector said "the shell is open" or "this pair of facets
intersects", and the caller was left holding the diagnosis. `tensorcad.mesh_repair` is the
missing path, and it is written in this repository's register: every step is measured, and a
repair that would trade one defect for another is **refused with the reason attached**
instead of being returned as a success.

Everything on this page is produced by `python examples/mesh_repair/run_example.py`, which
writes `artifacts/mesh_repair/results.json`, `plot.png` and `report.md` and exits non-zero if
any claim stops holding. `tests/test_mesh_repair.py` (20 tests) pins the same behaviour.
Numbers below are read from that artifact, not copied by hand.

## 1. What the module is, and what it is not

| step | what it does | refuses when |
|---|---|---|
| `weld` | merges vertices closer than the weld tolerance (1e-6 m default) | — |
| `drop_degenerate_faces` | deletes zero-area facets | — |
| `drop_duplicate_faces` | keeps one copy of coincident facets (the audit's `DUPLICATE_FACE`) | — |
| `manifoldize` | splits each vertex once **per face sheet**, so a bow-tie edge or vertex becomes two manifold sheets | never — it only splits |
| `orient_consistently` | floods a winding across manifold edges, then turns each closed sheet outward | a non-orientable sheet (a Möbius band has no consistent winding) |
| `fill_holes` | closes each simple boundary loop with an ear-clipped **planar** fill over the loop's own vertices | a loop that is not simple, is collinear, is longer than `max_hole_edges`, or whose fill has zero area |

Two design decisions carry the honesty of the module:

* **The fill invents no vertex.** The patch re-uses the loop's existing vertices, so it
  shares its boundary edges with the shell *exactly*; closure is a topological statement,
  not a tolerance. What that costs is shape: a hole in a curved surface is closed by a chord
  patch, and the distance from that plane is reported per hole
  (`HoleFill.holes[i]["planar_deviation_m"]`) rather than smoothed over.
* **Refusal returns the input.** `repair_mesh` is a question ("is this mesh usable?"), and
  `accepted=False` plus a reason is a more useful answer than a traceback or a quietly
  degraded mesh. The returned mesh is then the *original*, so a caller cannot mistake a
  refusal for a repair.

Explicitly absent, each with the reason it cannot be faked:

* **Self-intersection removal.** That is a remeshing problem (fTetWild-class: local
  re-triangulation plus a projection pass). This module detects, refuses, and says so. It is
  the one clause of GAP-CAD-009 that remains open, and it is reported as open.
* **Re-tessellation / refinement.** It can repair topology and winding, not a facet that is
  simply too coarse.
* **A guarantee that edge-manifoldness implies an intersection-free surface.** It does not:
  two sheets passing within the weld tolerance are topologically separate after
  `manifoldize` and still pierce each other geometrically. That is why
  `measure_mesh(...).self_intersections` is a separate number and part of the acceptance
  test.

## 2. How it was measured

Three analytic fixtures, each with a volume that does not depend on the triangulation:

| fixture | faces | volume | reference |
|---|---:|---:|---|
| box 1×1×1, 2×2 quads per face | 48 | 1.000000000000 m³ | analytic 1 m³ |
| UV sphere r = 0.5 m, 12 rings × 24 segments | 528 | 0.508819 m³ | analytic 4/3πr³ = 0.523599 m³ (**the fixture's own** −2.82 % tessellation error) |
| parametric hull L=10, B=2, T=0.8, 41×13 | 1998 | 7.446554 m³ | analytic 7.512698 m³ (−0.88 %, the surface error the hull audit records) |

Defects are then **injected** into a clean fixture — facets deleted, one facet wound
backwards, the whole shell turned inside out, a facet duplicated, a pole pulled through the
surface — and the repair has to put the measured numbers back.

## 3. Results

| case | measured |
|---|---|
| box, 1 facet removed | 47 → 48 facets, volume 1.000000000000 m³, relative change **0.0** |
| box, one quad (2 facets) removed | 46 → 48 facets, volume 1.000000000000 m³, relative change **0.0** |
| box, two holes on different faces | 47 → 48 facets, volume 1.000000000000 m³, relative change **0.0** |
| box, duplicated facet | 3 non-manifold edges → 0, volume 1.0 m³ |
| box, one backwards facet | 3 winding conflicts → 0, **1** facet flipped (minimal, not the 95 the flood start suggests) |
| box, inside-out shell | volume −1.0 → +1.0 m³, 48 facets flipped |
| sphere, 2 facets removed | closed, no self-intersection, volume 0.508819 m³ = the tessellated value (**+0.000 %**) |
| hull 41×13, 3 facets removed | 1995 → 1998 facets, open edges 5 → 0, volume 7.446529 m³ vs analytic −0.8808 % |

The box rows are the zero-error control: a hole bounded by a planar facet is closed by a
triangle in the same plane, so the volume is restored to round-off, not to 1 %. The sphere
and hull rows are the interesting ones, and they say the thing this repository insists on: a
chord patch closes the surface and gives back the *tessellation's* volume, not the analytic
one. The remaining −2.82 % (sphere) and −0.881 % (hull) were already in the fixture before
the hole was cut; the repair neither hides them nor adds to them.

## 4. Refusals, with the reason the module gives

| injected defect | outcome | reason as returned |
|---|---|---|
| sphere with the north pole pulled through the body (self-intersecting) | `accepted=False`, original mesh returned | "16 self-intersecting face pair(s); removing them is a remeshing problem, not a repair this module performs" |
| 12-facet hole against a 6-edge cap | `accepted=False`, hole left open and recorded in `holes.refused` | "loop has 14 edges, above the 6 cap; a planar fill is not a repair at that size" |
| any repair that moves the volume past the tolerance | `accepted=False` | "volume moved … against the reference (limit …)" — exercised in the test suite at `volume_tolerance=1e-9` |

The refusals are the point. A repair API that returns a mesh for every input is
indistinguishable from a repair API that returns the input for every input; the acceptance
record (`before`/`after` metrics, per-step log, per-hole records, notes on which checks were
skipped) is what separates them.

## 5. Error judgement against GAP-CAD-009

| criterion (GAP-CAD-009) | target | measured | verdict |
|---|---|---|---|
| closed + manifold after repair, on the injected defects | 100 % | 8 / 8 repairable cases closed and manifold | **met for the injected set** |
| volume change after repair | < 1 % | 0.0 % (box), +0.000 % (sphere), 0.88 % (hull, and that is the surface's own error) | **yes** |
| self-intersection set empty after repair | required | empty on every accepted case; a mesh that arrives self-intersecting is refused | **yes for accepted meshes; removal itself is not implemented** |
| the public fTetWild fixture corpus | 100 % repair rate | **not run** — the corpus is not vendored in this repository, and the synthetic injections above are substituted; this is a substitution, not an equivalent | **no** |

That last row is the honest gap: the anchor asks for the fTetWild corpus, this delivery
measures injected defects on analytic fixtures. The substitution is stated here rather than
left for a reader to discover, and it is the first thing a follow-up should close.

## 6. Notes, limits, and what to do next

* **Cost.** The steps are linear-ish in faces except `fill_holes`'s loop simplicity test
  (O(n²) per loop) and `measure_mesh`'s self-intersection grid; at 2 k facets the whole
  example runs in a second. Nothing here is intended for a 10⁶-facet import.
* **The fill is planar.** On a hull, a hole that spans a highly curved region is closed by a
  chord patch whose deviation is worth reading before accepting: it is printed per hole.
* **`orient_consistently` minimises the flips it reports** by inverting any sheet whose flood
  inverted more than half of it. Without that, a single backwards facet costs a whole-sheet
  flip plus a whole-sheet outward flip (95 changes for one wrong facet) — a real defect the
  test suite caught.
* **Next step.** Local re-triangulation for the self-intersecting case, and the fTetWild
  corpus, in that order: the first is the remaining functional gap, the second is what would
  turn the acceptance row above from *no* to *yes*.