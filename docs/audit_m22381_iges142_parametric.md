# GAP-CAD-T1-IGES-PIPE — entity 142 read in parameter space (round m22381)

Anchor: `GAP-CAD-T1-IGES-PIPE` (IGES/STEP → watertight mesh import pipeline, gold
regression = `DTMB5415_Hull_exact`). This round takes the first item of the return-set
written by the m17343 audit:

> 1. entity 142/141 的参数空间路径（当前 3-D 反演在 CATIA 导出上残差 ~1e2 mm，只有 1/8 环可反演）

Test asset: third-party trimmed IGES `KP505_modelScale_250mm.igs` (Kongsberg KP505
propeller, 1 042 550 B, SHA-256 `0bb54ae7618e8d8b1b7c23de8d384e04b0465071e34fe32e1d606ddbb5df462f`).
Not redistributable → the asset lives in the working tree under `benchmarks/anchors/sources/`
(git-ignored) and the integration test skips when it is absent.

## 1. What was wrong (measured on the unmodified reader)

The file census is `{100:8, 102:12, 108:2, 110:12, 120:2, 124:8, 126:53, 128:4, 142:8, 144:8}`:
eight trimmed faces, each with a single 142 record. Four of those records carry **two**
curve references — the 3-D curve and its image in the surface's parameter space. The old
`_curves_of_142` took *every* curve reference as a 3-D candidate and the face builder kept
all of them, so those faces acquired a second, invented boundary loop:

| face | entity | surface | loops before | worst loop residual | verdict |
|---|---|---|---|---|---|
| 0 | 144 @ DE 15 | plane (108) | 1 | 0.0 | meshed |
| 1 | 144 @ DE 45 | revolution (120) | 1 | 6.05e-08 mm | meshed |
| 2 | 144 @ DE 77 | revolution (120) | 1 | 6.21e-09 mm | meshed |
| 3 | 144 @ DE 93 | plane (108) | 1 | 3.55e-15 mm | meshed |
| **4** | 144 @ DE 123 | NURBS (128, 6×6) | **2** | **1.2606e+02 mm** (8.3e5 × gate) | **refused** |

(`_loop_to_uv` subsample probe of the unmodified reader; gate = 1e-6 × extent =
1.518e-04 mm, extent = 151.77.)

The bogus loop is the parameter-space image read as a 3-D curve: its coordinates are
`(u, v, 0) = (-1.945 … 0.998, -0.286 … 0.284, 0)` inside a blade that sits at
`(13.5 …, 113.4 …, 49.9 …)`. It is a perfectly valid closed 3-D loop, which is why no
positional check catches it — only a measurement does. Face 4's parametric image maps to
`(13.512 …, 113.453 …, 49.912 …)`, i.e. onto its own face.

## 2. What was changed (2 files, ~150 lines)

`tensorcad/geomquality/iges_trim.py`

* `TrimmedFace` may now carry, per loop, the **exact parameter-space image** (`uv_loops`)
  and the measured agreement with the 3-D curve (`cross_checks`); `to_dict` reports
  `trim_routes` and `parametric_cross_check`.
* `_curves_of_142` classifies every curve reference **by measurement**: a reference is a
  parameter-space image when its coordinates fall inside the surface domain *and*
  `evaluate(u, v)` reproduces another reference's 3-D polyline within
  `max(2 × chord tolerance, 1e-12)`; everything else stays a 3-D boundary. A reference
  that sits in the domain but does **not** reproduce the 3-D curve is recorded as a
  refusal with the measured disagreement instead of being turned into a boundary loop.
  A record whose only reference is inside the domain cannot be classified without an
  inversion probe (minutes on the 126×54 patch of this file), so it keeps the 3-D route
  and is judged by the mesh-stage gate — stated, not hidden.
* Agreement is measured **point-to-segment**, not point-to-vertex: the 3-D reference is
  itself a chord sampling, and the vertex metric reported 0.265 mm of "disagreement" for
  the same curve (0.043 mm by segment). Both numbers are pinned by unit tests.

`tensorcad/geomquality/trim_mesh.py`

* A loop with a parameter-space image is tessellated from it directly: no inversion, no
  branch ambiguity, no minutes-long walk over the 1 825-point 3-D curve of the same face.
  Per-face statistics gain `trim_routes` and `parametric_cross_check`; the inversion gate
  is unchanged for every loop that has no image (no silent pass).

## 3. Results after the change (KP505, third-party CATIA export)

Reader: 8 faces in 0.6 s, all 142 records resolved; the four surplus references are
recorded as `refused` with the text *"3-D curve … is the same trim as its
parameter-space image … (point-to-polyline agreement 4.285e-02 <= 4.540e-02 model units)"*.

Tessellation of the five-face subset (faces 0–4; the remaining three belong to the
126×54-control-point patches whose inversion is minutes of CPU and was **not run** this
round):

| face | surface | triangles | inversion residual | trim route | parametric cross-check |
|---|---|---|---|---|---|
| 0 | plane | 126 | 0.0 | inverted | — |
| 1 | revolution | — | 2.213e-04 mm (gate 1.146e-04) | inverted | — |
| 2 | revolution | 7 368 | 2.98e-05 mm | inverted | — |
| 3 | plane | 126 | 3.55e-15 mm | inverted | — |
| **4** | NURBS | **11** | **0.0** | **parametric** | **4.285e-02 mm** |

Total 7 631 triangles, 2.5 s. Face 4 went from *refused with a 126.06 mm residual* to
*meshed with a 0.043 mm agreement against an independent 3-D representation of the same
trim*.

## 4. Points to watch (limitations, stated)

1. **Face 1 is refused on a marginal gate**, not on a geometry error: residual
   2.21e-04 mm against a gate of 1.146e-04 mm — i.e. 1.9 × the gate, 1.5e-9 of the model
   extent. This is **pre-existing** (identical before and after the change) and is a
   gate-calibration question, not a trim-read question. Not silently relaxed here.
2. **Faces 5–7 were not tessellated** in this round: their in-domain reference disagrees
   with the 3-D curve of the same 142 record by **36.7 / 31.4 / 3.8 mm**, so it is not
   used as a trim, and inverting the remaining 178/228/122-point 3-D loops on those
   patches costs minutes. What *is* fixed is that the disagreement no longer invents a
   second boundary loop; those faces are now judged on their 3-D loop alone. Why the two
   representations disagree there is an **open question** (candidate explanations: the
   3-D curve is shared with a neighbouring patch; the parameter-space image traces the
   patch boundary rather than the trim). Untested — recorded, not explained away.
3. **No watertight solid yet.** Nine of the ten checks of the m17343 audit stay as they
   were: STEP branch still `NOT_SUPPORTED`, `tensor_benchmark_cases` row still absent,
   `wigley_iii.iges` still untracked, convergence table still missing.
4. The KP505 asset is not redistributable, so the integration test is asset-gated.

## 5. Error judgement (each key quantity vs its baseline)

| quantity | baseline | measured | error | < 3 % ? |
|---|---|---|---|---|
| parametric image vs 3-D curve, face 4 (cross-check) | 3-D curve of the same 142 (independent representation) | 4.285e-02 mm | 2.8e-04 of model extent (151.77) | **yes** |
| the same, measured to vertices instead of segments | — | 2.647e-01 mm | 1.7e-03 of extent | yes (but 6× the true value → metric rejected) |
| face 4 trim route | must not invent a loop | 1 loop, `parametric` | 126.06 mm → 4.3e-02 mm | **yes** (8.3e5× gate → within 2× chord tol) |
| face 4 triangles | > 0 (face must mesh) | 11 | 0 → 11 | **yes** |
| face 1 inversion residual | gate 1.146e-04 mm | 2.213e-04 mm | 93 % of the gate | **no** (pre-existing, unchanged) |
| faces 5–7 mesh | — | not run | — | **not applicable** (stated, not claimed) |
| full-shell watertightness of KP505 | trimesh/manifold | not run | — | **not applicable** |

## 6. Reproduce

```bash
cp <asset> benchmarks/anchors/sources/KP505_modelScale_250mm.igs   # git-ignored
python -m pytest tests/test_iges_142_parametric.py -q              # 7 passed
```

`artifacts/m22381_iges142/kp505_faces0_4.json` holds the per-face numbers above
(program-generated, not transcribed).

## 7. Next round (return set)

1. Explain and close the 36.7 mm disagreement on faces 5–7 (compare the 3-D curve's
   inversion residual on the patch; check whether it is shared with a neighbouring 144).
2. Give the mesh stage a budgeted inversion (cap seeds/iterations per point) so a
   126×54 patch cannot cost minutes; then run the full eight-face file.
3. Register the anchor row in `tensor_benchmark_cases` (metric = trimmed-face success
   rate and cross-check deviation, tolerance explicit) — the L4 gate still cannot pass.
4. Re-calibrate the residual gate for revolved surfaces (face 1 case) with a stated
   argument, not by loosening it until green.