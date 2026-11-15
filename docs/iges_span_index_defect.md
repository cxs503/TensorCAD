# IGES B-spline knot-span index: a silent wrong answer in the reader (m21680)

**Anchor** `GAP-CAD-T1-BASELINE-FLEET` · **Repository** TensorCAD · **Branch** `agent/m21680`
**Site** `tensorcad/iges.py::_find_span` · **Code** `tensorcad/iges.py` (fix), `tools/baseline_fleet/iges_span_defect.py` (measurement)
**Evidence** `runs/m21680_span_defect/evidence.json` · **Figure** `docs/figs/m21680_span_defect.png`
**Tests** `tests/test_iges_span_index.py` (9 passed) · `tests/test_iges_pipeline_gate.py` (updated, see §5)
**Reproduce** `python tools/baseline_fleet/iges_span_defect.py --out runs/m21680_span_defect/evidence.json`

## 1. The defect

`_find_span` picks the knot span that contains a parameter value; `_basis` and therefore
`BSplineSurface.evaluate` / `BSplineCurve.evaluate` index their control points from it. The
pre-fix early-out was

```python
lo, hi = degree, n                      # n = count - 1
if value >= knots[hi]:                  # knots[n], should be knots[n+1]
    return hi - 1                       # n - 1, should be n
```

Two errors in three lines:

* the threshold is the *second-to-last* support knot, not the last (`knots[n+1]`), so the
  early-out fires early;
* it returns `n-1`, so `_basis` starts its index window at `span - degree + i = -1`.

Python does not raise on `control[-1]`: it wraps to the **end** of the knot/control vectors.
The basis functions that result are `0.0` (measured on every affected patch here), the
denominator in `evaluate` is `0.0`, and the function returns `(0.0, 0.0, 0.0)` — a **silent
wrong answer**, not an exception, not a warning.

A direction whose knot vector has a single span (a Bezier patch, `count == degree + 1`, e.g.
`[a]*6 + [b]*6` with `count = 6, degree = 5`) hits the early-out for **every** parameter, so
the whole patch evaluates to the origin. That is exactly the KVLCC2 archive: 5 of its 19
type-128 surfaces are degree-5 Bezier patches in one direction.

## 2. Measured cost on the official archives

Per-surface health, 13x13 samples per surface, same bytes, same process, only the span index
switched (`runs/m21680_span_defect/evidence.json`). "all-zero" = samples returned as the
origin; "outside" = samples outside the surface's own control-net bounding box, which a
correct B-spline evaluation cannot produce (convex-hull property, positive weights).

| ship | file | surfaces | all-zero samples (pre → post) | outside-net samples (pre → post) | surfaces with a defect (pre → post) |
|---|---|---:|---:|---:|---:|
| KCS | `KCS_hull_SVA.igs` (246 164 B) | 55 | 0 → 0 | 699 → **0** | 54 → **0** |
| KVLCC2 | `kvlcc2.igs` (6 270 129 B) | 19 | 773 → **0** | 870 → **0** | 19 → **0** |
| JBC | `jbc.igs` (1 713 554 B) | 45 | 5 343 → **0** | 5 079 → **0** | 41 → **0** |
| ONRT | `onr5613-full_appendage_deckhous_mergedR.igs` (3 374 136 B) | 132 | 17 968 → **2** | 17 867 → **0** | 128 → **2** |

The KCS row is the quietest and the most instructive: no surface evaluated to the origin, yet
**54 of 55** surfaces had samples outside their own control net — every patch's last row and
column of samples was evaluated on the wrong span. The archive "worked" because a boundary
ring of a lofted patch is a small fraction of its area.

The curve path is affected too (same `_basis`), and that is not cosmetic: on KVLCC2 the
type-142 trim loops are composite type-126 curves, and one patch's trim polygon could not be
decided by `tensorcad.fleet.trimfield.trim_polygons` before the fix. Measured:

| quantity (KVLCC2, 19 trimmed patches) | pre-fix | post-fix |
|---|---:|---:|
| trim patches decided | 18 / 19 | **19 / 19** |
| tessellated triangles at grid 40 inside the loops | 9 446 | **50 754** |
| height-field columns (2 m cells) | 485 | **2 218** |

## 3. Effect on the one volume that can be checked today

KCS is the fleet's only ship with an accepted official-geometry volume (∇ = 52 196.1 m³,
+0.319 % vs the published 52 030 m³, built through `iges_field` at dx 1.0 / dz 0.2 in the
m19948 mission). Two independent height-field quadratures of the same 55 patches after the fix
(`tools/baseline_fleet/iges_span_defect.py`, z-buffer of the lowest surface over an (x, y)
raster, waterplane contribution included, half hull doubled):

| quadrature | volume (m³) | vs 52 030 m³ |
|---|---:|---:|
| grid 24 / cell 0.50 m | 52 432.1 | **+0.773 %** |
| grid 32 / cell 0.25 m | 52 366.3 | **+0.646 %** |

Both are inside the 1 % acceptance band and agree with the established `iges_field` number to
0.5 % absolute — the sign and size of the residual are consistent with the `max|y|` binning
bias that the `iges_field` route documents, plus first-order raster error. This quadrature is
a **cross-check, not a certification**: it is resolution-dependent at the 0.13 % level between
the two settings, and the certified KCS number remains the `iges_field` one.

## 4. What is *not* fixed, stated plainly

* **KVLCC2 volume: still not delivered.** The reader now decodes all 19 patches and all 19 trim
  loops, but the trim *semantics* are still wrong for this file: with the loops applied, the
  kept strips cover only part of the hull (2 218 height-field columns against the ~9 700 a
  320 x 58 m hull needs), so the resulting volume is far below 312 622 m³. The strip pattern
  (e.g. patch 2 keeping only x ∈ [-244, -163]) says the accepted loop is a *sub*-loop of the
  true patch boundary, not the boundary itself — an open question about composite (type 102)
  loop traversal, not about this defect. No number is claimed for KVLCC2.
* **JBC / ONRT**: with the reader repaired, JBC's cloud is sound but the file's frame places no
  hull inside `z ∈ [0, T]` (draft window empty), so the `iges_field` route still refuses it;
  ONRT's 132 patches still leave 2 all-zero samples (2 surfaces, 1 sample each), and its
  geometry spans 0.18 m in z — the appendage/deckhouse archive is not a hull-only file.
  `tests/test_iges_pipeline_gate.py::test_jbc_fleet_baseline_status` was updated from "asserts a
  WITHHELD mesh with a > 1 % volume error" to "asserts the named refusal", because the mesh it
  used to assert was built from a broken reader.
* **DTC / Series 60**: unchanged, sources not located (`WITHHELD`).

## 5. Why this was reachable at all

Nothing in the pipeline refused a surface that evaluates to the origin. The route that *did*
refuse (DTC, "sampling collapsed to a point") refuses only when the **whole cloud** collapses;
a file with 5 of 19 dead patches plus 14 live ones produces a plausible-looking cloud and a
plausible-looking wrong hull. `tools/baseline_fleet/iges_span_defect.py::surface_health` is the
check that would have caught it (per-surface all-zero and outside-control-net counts); it is
kept as a tool rather than wired into the gates because wiring it in belongs with the fleet
gate ownership, not with this fix.

## 6. Error judgement

| quantity | reference | measured | error | < 1 %? |
|---|---|---:|---:|---|
| KVLCC2 surfaces evaluating to the origin (13x13 samples) | 0 | 773 pre-fix, **0** post-fix | — | n/a (count) |
| KCS surfaces with samples outside their control net | 0 | 54 pre-fix, **0** post-fix | — | n/a (count) |
| KCS displacement volume (grid 32 / cell 0.25 m) | 52 030 m³ | 52 366.3 m³ | **+0.646 %** | yes |
| KCS displacement volume (grid 24 / cell 0.50 m) | 52 030 m³ | 52 432.1 m³ | **+0.773 %** | yes |
| KVLCC2 displacement volume | 312 622 m³ | not computable (trim semantics open) | — | **no claim** |

The two volume figures are a *quadrature cross-check*; the defect claim itself is the sample
counts, which are exact integers produced by the same code on the same bytes.