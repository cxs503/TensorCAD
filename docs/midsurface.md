# Mid-surface extraction and defeaturing (`GAP-CAD-013`)

`tensorcad.midsurface` turns a **thin-walled closed solid** into the **mid-surface sheet** a
shell structural model needs, and it reports the thickness it measured, the area it covered,
and the mesh it saved. Before this module the package could hollow a solid (`shell`) and repair
a mesh (`mesh_repair`), but nothing produced the sheet itself.

## Why a ship needs it

A double bottom, a bulkhead, a wing tank or a stiffened panel arrives as a solid whose wall is
a few millimetres in a hull a few hundred metres long. A solid-element model of that is
unrunnable (aspect ratios ~1e4); the structural model is the sheet plus one thickness per
patch. The module is the geometric half of that pre-processing step:

```python
from tensorcad.midsurface import extract_midsurface, defeature
from tensorcad.shell import shell, hollow_body

body  = hollow_body(shell(box_mesh, 0.012)).mesh   # or any imported thin-walled solid
sheet = extract_midsurface(body)                   # MidSurface: mesh + panels + thickness map
sheet.area(), sheet.volume(), sheet.face_reduction(), sheet.thickness_map.coverage()
report = sheet.as_dict()                           # JSON-ready
```

## Method (three steps, each with its own number)

1. **Measure the wall.** From every triangle centroid a ray is cast along the inward normal
   through a uniform triangle grid (Möller–Trumbore). The first hit is the opposed wall if
   `-dot(n_src, n_hit) >= opposition` and the distance is within `max_thickness`; the pair is
   kept only if it is **mutual** (the reverse ray returns the original face). Output: a
   `ThicknessMap` — thickness and partner per face, plus a named reason for every failure.
2. **Group faces into panels.** Faces are linked when they share a coincident vertex (positions
   matched within `weld_tolerance`) *and* their planes agree within `angle_tolerance` /
   `offset_tolerance`; panels are the connected components. This is what absorbs a CAD import
   that split one plane into sixteen patches.
3. **Place the sheet.** Each mutual pair is used **once** (its smaller-indexed face), so a
   hollow body yields one sheet and not two coincident ones. The panel plane is offset inward
   by `t/2`; the patch boundary starts as the convex hull of the projected panel vertices and is
   then clipped by the **adjacent panels' mid-planes**. That clip is what a corner needs: a box
   shelled at `t` gives six patches meeting on `x = y = t/2`, i.e. exactly the closed sheet of
   the eroded box. A neighbour that is not a panel (a plate's edge strip, whose opposed face is
   `L` away, not `t`) clips nothing — and a plate's sheet should indeed run out to the edge.

`defeature(midsurface, ...)` re-runs steps 2–3 with relaxed tolerances and an optional area
floor: import seams a hair apart are bridged within `weld_tolerance` (relax `convexity_tolerance`
with it — the bridged sliver is precisely what the convexity check measures) and features below
`min_panel_area` are dropped. Nothing is promised about the saving; it is *measured*:

| quantity | meaning |
|---|---|
| `face_reduction()` | `1 - sheet_faces / source_faces` |
| `removed_area()` | source area that entered no panel |
| `thickness_map.coverage()` | qualified (mutual-pair) area / total area |
| `watertight()` | is the sheet a closed surface (closed body) or a patch (open body) |

## Verification (`tests/test_midsurface.py`, 88 tests with `test_shell.py`)

| check | reference | measured |
|---|---|---|
| box shelled at `t`, `t/L = 0.01…0.1` | area `6 (L-t)^2`, volume `(L-t)^3` | rel. error **< 1e-14** |
| box sheet corners | the six lines `x = y = t/2` | corner set identical to 1e-9 |
| regular 6/12-gon prism | independent 2-D section erosion, extruded (`midsurface_oracle.py`) | rel. error **5e-16** |
| the oracle itself | `n r^2 tan(pi/n)`, `2 n r tan(pi/n)` | rel. error < 1e-12 |
| sheet position | brute-force point-to-triangle distance to the source surface | every vertex at `t/2` within **1e-9** (box, refined box) |
| thickness guard | `thickness_tolerance` shrunk to 1e-20 on a correct wall | fires → the guard is live |
| face-count criterion (plate `8x8`) | sheet/source ≤ 1/10 | **0.0026** (768 → 2 triangles) |
| feature suppression | ≥ 50 % mesh cut | seam bridging **75 %**, small-feature suppression **50 %** |
| plate mid-sheet | area `L*B`, thickness `t`, open sheet | rel. error < 1e-12, `watertight False` |

## Refusals (a wrong sheet still looks like a sheet)

* open or inward-oriented surface → `MidSurfaceError` naming the fix (`weld`, `orient_outward`);
* no opposed thin walls within `max_thickness` → names the qualified fraction;
* wall thickness not constant inside one panel (> `thickness_tolerance`) → names the spread;
* a panel that is not convex once flattened → names the measured hull overshoot;
* `min_coverage` not met → names the fraction.

## Limitations (honest)

* **B-rep rim, face removal, stitching back into a solid** are not implemented (`GAP-CAD-001`).
* **Non-convex panels are refused, not split.** An L-shaped cap (an L-bracket, a bracket in a
  double bottom) needs the outline and an ear-clipping/boolean triangulation; the current hull
  would over-cover it, so the call fails and says so. Splitting such a panel upstream is the
  workaround.
* **A solid cube is indistinguishable from a very thick-walled box** by measurement alone, so
  `max_thickness` is the caller's statement of what a wall is. The default
  (`suggest_max_thickness` = 1.5 × smallest bounding extent) is a guess about the *input*, and it
  accepts a solid cube as "a wall of its own side"; the report always carries the value used.
* **Variable thickness** inside a panel is refused rather than averaged, and a wedge's slanted
  wall does not pair at all (its two faces are not parallel), so wedges must be modelled as
  constant-thickness panels.
* **A refined double-walled body** (the offset surface of a subdivided solid) exposes long
  spurious pairs when `max_thickness` is loose (1.5 × the smallest extent accepts the opposite
  wall of a 1 m box as a "wall"), so pass a `max_thickness` close to the real wall for those.
* **Faceted geometry**: for a faceted tube the sheet is exact for the *faceted* solid and
  `O(h^2)` short of the analytic tube. One measured open case: a regular dodecagon prism with
  `t = 0.05` (inradius 1) reports **1.04 %** low on area, above the 1 % criterion; the hexagon at
  the same `t` and the dodecagon at `t = 0.02` are exact. The corner trim of short facets is the
  suspect and is open work.

## A defect this work exposed in `shell`

`tensorcad.shell` solved the inward offset at each vertex for the **absolute point** with a
ridge term, so a vertex whose incident facet planes span only 1 or 2 dimensions (any vertex on a
flat facet interior or on a straight edge) was pulled toward `(0, 0, 0)` — invisible to every
existing fixture, because a box/prism/sphere corner always has three independent planes.
`_solve_vertex` now solves for the **displacement** from the vertex, which is what the ridge was
meant to mean, and the "inverted solid" refusal now also catches a *collapsed* cavity (at exactly
half the wall the inner volume is 0, and a rounding sign is no basis for calling that a shell).
`tests/test_shell.py` (48 tests) still passes.
