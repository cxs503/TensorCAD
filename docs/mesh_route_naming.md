# Naming facet-scale defects on the mesh side

`tensorcad.geomquality` audits a BRep and names fourteen defect codes in patch terms. The
mesh route (`tensorcad.validation`, `tensorcad.geometry_audit`, `tensorcad.mesh`) reads the
triangulation a solver actually consumes. The 2026-11-09 cross-route audit
(`docs/audit_crosscheck.md`) ran both over one 3-hull x 14-injection corpus and found that
six *topology* families survive on the mesh side (open / non-manifold / duplicate / winding
deltas) while four *facet-scale* families -- `SLIVER`, `MICRO_EDGE`, `CURVATURE_JUMP`,
`TANGENT_JUMP` -- all collapsed into the one sentence "the shell is open". That audit
recorded the gap rather than papering over it: *"口径缺失（如实记录，未粉饰）"*.

`tensorcad/geomquality/mesh_defects.py` closes the part of that gap a triangulation can
close. Every number on this page is written by `python examples/mesh_route_naming/run_example.py`
into `artifacts/mesh_route_naming/results.json`, which exits non-zero if a claim fails;
`tests/test_mesh_defects.py` pins the synthetic contracts.

## 1. What the module adds, and what it refuses to add

| code | criterion | measured on the corpus |
|---|---|---|
| `MICRO_EDGE` | an edge shorter than the weld tolerance (1e-6 m -- the tolerance the BRep route calls a micro edge) **and** more than 1 % of the edges | clean 0.0-0.4 % of edges, injected 7.6-7.9 % |
| `SLIVER` | a face at or below `(1e-4 L)**2` -- the BRep route's own sliver expression -- **and** more than 2 % of the faces (the limit `inspect_mesh_quality` uses) | clean 0-1.2 %, injected 20.5-21.2 % |
| `ISOLATED_PATCH` | more than one face component under shared-edge adjacency | clean 1 component, injected 2 |

It **does not** re-implement `FREE_EDGE`, `OPEN_TOP`, `NON_MANIFOLD`, `ORIENTATION` or
`DUPLICATE_FACE`: those belong to `audit_watertightness`, and a second dialect for one fact
is the disease the cross-route audit documented. Open and non-manifold edge counts are
reported as *context* in `metrics`, never as codes.

It **does not** claim `CURVATURE_JUMP` or `TANGENT_JUMP`. Measured at the corpus density the
two are indistinguishable from the facets: both give 1504 faces on the Wigley III (992
clean), 2214 interior edges (2210 and 2214), and 376 vs 408 edges above 5 degrees. The
information that separates them is the second derivative of the surface, and a 16x8
triangulation does not carry it. The honest output for those two, and for the ripple, is
`ISOLATED_PATCH`: "this surface is not one connected shell", with the class *not*
recoverable. The example asserts exactly that (`seam_and_ripple_families_are_not_overclaimed`).

## 2. Why shares, not counts

The BRep route fires on a single micro edge or sliver, which is right for a patch: a CAD
face is a design object. A tessellation is not. All three reference hulls produce faces at
or below `(1e-4 L)**2` at the stem, where the surface closes to zero width (12 of 988 faces
on the Series 60 parent, 0 on the Wigley), and two of them produce a handful of
sub-micrometre edges at the same place. A count-based detector would therefore fire on a
**clean hull** -- the failure mode this repository treats as worse than a miss. The share
limits sit an order of magnitude above the worst clean measurement and an order of
magnitude below the weakest injected one, and the example prints both ends so the margin is
checkable rather than promised.

The declaration the mesh route cannot do without is the same one the cross-route audit
identified: `declared_length`. Without it the sliver criterion is **not evaluated** and the
report says so (`sliver_area_limit is None`) instead of inventing a scale.

## 3. Repair: remove facet damage, or refuse in writing

`repair_facet_defects` performs exactly two operations -- weld vertices within the weld
tolerance, drop facets whose area falls to or below the sliver limit -- and then accepts
only if the result has no new codes, no more open edges, and moves the signed volume and
surface area by less than a declared relative tolerance. On the corpus, of the six
`SLIVER`/`MICRO_EDGE` models:

* every **accepted** repair ends with zero codes (`an_accepted_repair_never_leaves_a_facet_code`);
* every **refused** repair names the number that tripped it (open-edge count, or a volume or
  area change), and returns the original mesh unchanged;
* a defect that welding cannot fix (an island, a curvature jump) is refused with
  `no facet-scale repair available for <codes> -- that is geometry, not a facet`, so an
  accepted repair can never be an empty no-op that reads like a fix.

## 4. Limits of this page

1. The corpus tessellates at 16x8 (`brep.tessellate(16, 8)`), the density the cross-route
   audit used. A finer triangulation resolves a knuckle better and may separate
   `TANGENT_JUMP` from `CURVATURE_JUMP`; that study is not run here.
2. The corpus's facet injections add an **unsewn** patch, so `ISOLATED_PATCH` fires on them.
   The area and length criteria are therefore also pinned on *sewn* synthetic meshes in
   `tests/test_mesh_defects.py` (`sewn_sliver`, `sewn_micro_edge`), where the island shortcut
   is unavailable and the code still fires. Without that test the corpus result would be
   weaker than it looks.
3. `repair_facet_defects` does not sew islands and does not touch curvature: welding vertices
   cannot fix a G1 seam, and silently deleting an island would delete real appendages.
4. No hydrodynamic quantity is computed here, so this page makes no claim about a solver's
   answer -- only about what a gate can name before a solver sees the mesh.
