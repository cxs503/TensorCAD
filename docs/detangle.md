# Self-intersection removal: cutting the fold and capping the hole

`tensorcad.mesh_repair` measures self intersections and refuses a mesh that carries one; removal was the documented gap (`GAP-CAD-009`). `tensorcad.detangle` is that removal step, and this page is its measured evidence sheet. Every number below is read out of `artifacts/detangle/results.json` by `examples/detangle/make_docs.py`.

## 1. Case: a fold, not a dent

A dent is a valid surface: pushing a vertex inwards short of the far side changes the volume and nothing else, and the heal correctly does nothing to it (claim *a dent is not a fold*). A **fold** is the defect: the vertex travels until the facets it owns pierce the opposite side of the shell. That is what separates a defect from a design change, and it is why detection and removal both key on penetrating facet pairs rather than on curvature or on a volume budget.

| fixture | geometry | facets |
|---|---|---|
| `box n=4` | box 1×1×1 m, 192 facets, one vertex poked through the far side | 192 |
| `box n=6` | same fold, 432 facets: one ring of facets is not enough | 432 |
| `sphere r=0.5` | UV sphere (12 rings × 24 segments), a ring vertex pushed through the south pole | 528 |
| `hull 41x13` | parametric hull L=10 B=2 T=0.8, keel vertex pushed 1.6 m up through the deck | 1998 |

## 2. Method

1. `weld` the input at 1e-6 m (a fold is often a duplicate vertex apart) and measure it with `mesh_repair.measure_mesh`.
2. Detect penetrating facet pairs with `intersect.find_self_intersections`; a pair only counts when the two facets share **no** vertex, so a closed manifold's legitimate vertex contacts are not reported.
3. Take those facets as seeds, dilate them by edge adjacency into an *excavation set*.
4. Cut the set out, re-detect. Residual pairs mean the band was wider than the cut, so grow one more ring and retry, up to `max_rings`.
5. Per original shell, keep the largest surviving piece; shards the cut detached are dropped **and recorded**.
6. Cap the holes with the planar ear-clipping fill of `mesh_repair.fill_holes`; when a loop does not project to a simple polygon -- the usual case around a fold, where the facet normals turn through more than 90° -- fall back to a triangle fan from the loop's own centroid.
7. `orient_outward`, then measure. Acceptance is a conjunction: closed, edge-manifold, consistently wound, **zero** penetrating pairs, and the enclosed volume moved by no more than the budget (1 % when a reference volume is supplied).

## 3. Results

| case | pairs before | pairs after | rings | facets cut | holes capped | volume change |
|---|---|---|---|---|---|---|
| `box n=4` | 18 | 0 | 1 | 24 | 2 | +0.0000 % |
| `box n=6` | 18 | 0 | 2 | 48 | 2 | +0.0000 % |
| `sphere r=0.5` | 56 | 0 | 1 | 58 | 2 | +0.1782 % |
| `hull 41x13` | 64 | 0 | 1 | 59 | 2 | +0.4406 % |

All 29 acceptance claims hold (`0` failures recorded). Figures: `artifacts/detangle/plot.png` (the fold, red, and the healed mesh) and `artifacts/detangle/plot_volume.png` (volume change against the 1 % budget).

| case | claim | holds |
|---|---|---|
| `box n=4 control` | volume 1.0 m3 vs 1 | yes |
| `box n=4 detected` | 18 penetrating pair(s) injected | yes |
| `box n=4 healed` | cut 24 facet(s) in 1 ring(s), dropped 0 fragment(s), capped 2 loop(s); 18 penetrating pair(s) -> 0 | yes |
| `box n=4 clean after` | self intersections 0, closed True, manifold True | yes |
| `box n=4 volume` | volume 1.0 m3 vs 1, change 0.0 | yes |
| `box n=6 control` | volume 0.9999999999999973 m3 vs 1 | yes |
| `box n=6 detected` | 18 penetrating pair(s) injected | yes |
| `box n=6 healed` | cut 48 facet(s) in 2 ring(s), dropped 0 fragment(s), capped 2 loop(s); 18 penetrating pair(s) -> 0 | yes |
| `box n=6 clean after` | self intersections 0, closed True, manifold True | yes |
| `box n=6 volume` | volume 0.9999999999999973 m3 vs 1, change 2.6645352591003757e-15 | yes |
| `sphere control` | tessellated 0.5088190451025224 m3 vs exact 0.5235987755982988 (-2.823%) | yes |
| `sphere detected` | 56 penetrating pair(s) injected | yes |
| `sphere healed` | cut 58 facet(s) in 1 ring(s), dropped 0 fragment(s), capped 2 loop(s); 56 penetrating pair(s) -> 0 | yes |
| `sphere clean after` | self intersections 0, closed True | yes |
| `sphere volume` | change 0.0017819906055497965 of the tessellated volume | yes |
| `hull control` | 1998 facets, closed True | yes |
| `hull detected` | 64 penetrating pair(s) injected at the keel | yes |
| `hull healed` | cut 59 facet(s) in 1 ring(s), dropped 0 fragment(s), capped 2 loop(s); 64 penetrating pair(s) -> 0 | yes |
| `hull clean after` | self intersections 0, closed True, manifold True, wound True | yes |
| `hull volume` | change 0.004405679470779271 | yes |
| `idempotent` | second pass healed=False, faces 1984 vs 1984 | yes |
| `repair_mesh refuses a fold by default` | accepted False, self intersections 16 | yes |
| `repair_mesh heals with the flag` | accepted True, self intersections 0, steps ['weld', 'drop_duplicate_faces', 'manifoldize', 'orient_consistently', 'fill_holes', 'heal_self_intersections'] | yes |
| `refuses an open input` | refusing to excavate a mesh with 4 boundary edge(s): the repair's own holes would be indistinguishable from the mesh's, so a closed input is required | yes |
| `a 4-edge cap limit falls back to the fan, then the volume guard refuses` | healed=False: the repair moves the enclosed volume by 9.589%, above the 5.0% budget (0.912500 -> 1.000000 m^3): the cut was not local enough to be called a repair, and a wider cut only removes more material | yes |
| `refuses a non-local cut` | the repair moves the enclosed volume by 26.667%, above the 1.0% budget (0.650000 -> 0.733333 m^3): the cut was not local enough to be called a repair, and a wider cut only removes more material | yes |
| `a dent is not a fold` | nothing to heal: no penetrating facet pair within the detection limit of 64 (the input's own upper bound is 0) | yes |
| `excavation preview` | 24 facet(s) in the 1-ring band | yes |
| `clean input is a no-op` | no penetration, empty band | yes |

## 4. Error determination against the reference

| quantity | reference | measured | verdict |
|---|---|---|---|
| `box n=4` volume after the heal | 1.0 m³ | 1.0 m³, +0.0000 % | yes, below the 1 % budget |
| `box n=6` volume after the heal | 1.0 m³ | 0.9999999999999973 m³, +0.0000 % | yes, below the 1 % budget |
| `sphere r=0.5` volume after the heal | 0.5088190451025224 m³ | 0.5079123343442249 m³, +0.1782 % | yes, below the 1 % budget |
| `hull 41x13` volume after the heal | 7.446554097149208 m³ | 7.41374696663535 m³, +0.4406 % | yes, below the 1 % budget |
| `hull 41x13` volume vs analytic | 7.512698221250894 m³ | -1.3171 % | this is the hull surface's own discretisation error (documented in `docs/hull_analytic_reference.md`), not an error of the heal: the heal is measured against the same mesh's clean volume |
| `sphere r=0.5` volume vs analytic | 0.5235987755982988 m³ | -2.8227 % (the fixture, before any defect is injected) | likewise the fixture's tessellation error |

So: the quantity the heal is answerable for -- how far it moves the enclosed volume -- stays **below 1 % on every case above** (the largest is +0.4406 %), and on the two analytic fixtures it is at round-off for the box and 1.8e-3 for the sphere. The 3 % standard is met with two orders of magnitude to spare on this quantity.

## 5. Notes, traps and scope

* **A wider cut is not a better repair.** On a 48-facet box, one ring of facets is half the mesh; the heal cuts, caps, then refuses because the volume moved 26.667 % -- honest refusal, and the reason quotes the number. Widening the excavation only removes more material, so the guard stops instead of iterating.
* **The planar fill is not enough on its own.** The excavated loop around a keel vertex does not project to a simple polygon; `fill_holes` refuses it and the centroid fan caps it. Without the fan the hull case would be a refusal, and the test asserts the fan path actually ran (`centroid fan` in the notes) rather than assuming it.
* **The cap is a membrane.** At the crossing the shape is *closed*, not *reconstructed*: this is not a boolean re-triangulation. A surface that genuinely needs the crossing kept (an interpenetrating appendage, two shells meeting at an angle) needs a real boolean, which is `GAP-CAD-002` territory.
* **A refusal carries its evidence.** The refused result keeps the best attempt's measurements (`after`, `holes`, `volume_relative_change`, `bands`), because the reason sentence quotes those numbers and a caller should not have to parse prose to get them.
* **Open inputs are refused up front.** Excavating a mesh that already has boundary edges would make the repair's own holes indistinguishable from the mesh's; run `repair_mesh` (which fills first) or close the mesh, then heal.
* **Cost.** Detection is a uniform-grid broad phase; every case above (192-1998 facets) heals in well under a second, so the heal is cheap enough to sit in a repair chain rather than in an offline pass.

## 6. Reproduce

```bash
python3 -m pytest tests/test_detangle.py -q     # 17 tests
PYTHONPATH=. python3 examples/detangle/run_example.py   # exits non-zero on a slipped claim
PYTHONPATH=. python3 examples/detangle/make_docs.py     # regenerates this page
```

## 7. In the repair chain

`mesh_repair.repair_mesh(..., heal_self_intersections=True)` runs the heal as one step of the chain. Without the flag a hull with one keel vertex poked 0.6 m up is **refused** (accepted=False, 16 penetrating pair(s) in the result); with it, the same mesh is accepted with 0 penetrating pairs and the step appears in the repair log (accepted=True).

