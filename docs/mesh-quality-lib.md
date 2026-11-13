# Mesh-quality library: canonical measures, differential evidence, and the CI gate

*Gap:* `GAP-CAD-D1-MESH-QUALITY-LIB`. *Mission:* m17404 (round 2 — round 1, m17355,
left no artefact and no failure reason, so this page front-loads both the artefact and the
reproduction command.)*

## 1. Problem, and where the reference values come from

Before this module the repository answered *is this triangulation good enough?* in four
dialects, all of them named **aspect ratio**, and none of them the number a mesher or a
solver preprocessor reports for the same mesh. Measured on the meshes of section 4
(`dialect_audit` in the artefact):

| mesh | `gate._facet_aspect` l_max²/A | `validation` l_max²/2A | `parametric_hull` l_max/l_min | canonical `aspect_ratio` | canonical `edge_ratio` |
|---|---|---|---|---|---|
| box_unit | 4 | 2 | n/a | 1.39385 | 1.41421 |
| series60 | 12.3357 | 6.16787 | n/a | 3.8244 | 6.06143 |
| sliver_patch | 800.005 | 400.002 | n/a | 231.23 | 400.001 |
| suboff_coarse | 29652.2 | 14826.1 | n/a | 8559.93 | 9211.34 |
| suboff_taubin30 | 16.8562 | 8.42811 | n/a | 4.93372 | 4.56865 |
| wigley | 10.1932 | 5.09659 | n/a | 3.20676 | 4.9916 |

`gate._facet_aspect` is `2`–`2` times `validation.measure_quality`'s "aspect ratio" across the 6 meshes: the same shape under the same name with a different constant.
That is the gap — not a missing formula, a missing *contract*.

Reference definitions, stated plainly:

* the meshmagick measure set is what `meshmagick.Mesh.print_quality` reports
  (`meshmagick/mesh.py:459-540`, source staged at `/tmp/refs/meshmagick-master`); it
  delegates to `vtk.vtkMeshQuality`, a transcription of the Verdict library (*Verdict*
  *Library Manual, rev A* — the manual meshmagick's own output points at).
* the Gmsh measures quoted here are `gamma` (quality) and the `minSICN` family
  (*Gmsh reference manual*, "Mesh quality"; `gmsh.model.mesh.getElementQualities`).
* **what is not claimed:** this host has no VTK, no Gmsh and no runnable meshmagick
  (`upstream` in the artefact records each import failure verbatim), so the measures are
  *transcribed*, then bounded by closed forms and exact arithmetic — **not** diffed against
  the upstream binaries. `tools/mesh_quality_crosscheck.py` runs them the moment `vtk` or
  `gmsh` imports; the code path exists and is not hand-skipped.

## 2. Method and implementation

Every measure is a formula in the edge lengths `l`, the area `A` and the semiperimeter
`s`, normalised so the equilateral triangle scores 1, so the value is scale free and
comparable across hull sizes.

| measure | definition (source of truth: `METRIC_DEFINITIONS`) | ideal |
|---|---|---|
| `area` | A, the triangle area (model units squared). | n/a |
| `min_angle_deg` | smallest interior angle, degrees. | 60 |
| `max_angle_deg` | largest interior angle, degrees. | 60 |
| `edge_ratio` | l_max/l_min -- Verdict 'Edge Ratio' (tri), ideal 1. | 1 |
| `aspect_ratio` | l_max*(l_a+l_b+l_c)/(4*sqrt(3)*A) -- the edge-length form of Verdict 'Aspect Ratio' (tri), ideal 1. NOT the radius form R/(2r): the two are not equal for a triangle (3-4-5 example in docs/mesh-quality-lib.md: 1.4434 vs 1.2500), so the radius form is reported separately as 1/radius_ratio. | 1 |
| `aspect_frobenius` | (l_a^2+l_b^2+l_c^2)/(4*sqrt(3)*A) -- Verdict 'Aspect Frobenius' (tri), ideal 1. | 1 |
| `gamma` | 4*sqrt(3)*A/(l_a^2+l_b^2+l_c^2) -- Gmsh 'gamma' (quality), ideal 1. | 1 |
| `radius_ratio` | 2*r/R = 8*A^2/(s*l_a*l_b*l_c) with s the semiperimeter -- Verdict 'Radius Ratio' (tri) layout, ideal 1. | 1 |
| `scaled_jacobian` | 2*A/(l_mid*l_max) -- scaled Jacobian of the reference->physical map (a.k.a. min scaled Jacobian), ideal sqrt(3)/2. | 0.866025 |
| `condition_number_inverse` | sigma_min/sigma_max of the 2x2 Jacobian [[v1-v0],[v2-v0]] (Gmsh 'minSICN' family), computed from the two invariants sigma1*sigma2 = 2A and sigma1^2+sigma2^2 = |u|^2+|v|^2; ideal 1/sqrt(3). | 0.57735 |

Code: `tensorcad/geomquality/quality.py` (`triangle_quality`, `mesh_quality`,
`quality_gate`); independent route `tools/quality_oracle.py` (Heron area, law of cosines
from the sides only, `2r/R`, min over the three apices for the scaled Jacobian, LAPACK SVD
for the conditioning, and a `fractions.Fraction` path up to a single square root);
differential driver `tools/mesh_quality_crosscheck.py`; fast CI subset
`tests/test_mesh_quality_lib.py`.

One definition actively corrected during this work: the edge-length form `l_max(l_a+l_b+l_c)/(4√3 A)` is **not** the radius form `R/(2r)` for a triangle (3-4-5: 1.44338 vs 1.25). Both are exposed (`aspect_ratio`, and `1/radius_ratio`) because "aspect ratio" alone is ambiguous, and an ambiguous number cannot be audited.

## 3. Parameters

| tolerance | value |
|---|---|
| independent route | 1e-06 (mission acceptance) |
| exact arithmetic | 1e-12 |
| exact arithmetic, angles | 1e-09 (see §5) |
| closed form | 1e-12 |

| gate measure | sense | limit | why |
|---|---|---|---|
| `aspect_ratio` | max | 25 | the repository's existing decimation guard (25), canonically expressed |
| `edge_ratio` | max | 100 | micro-edge detector; clean hulls sit two orders below it |
| `gamma` | min | 0.05 | 1/20 of the ideal element — below this a facet is a cap |
| `min_angle_deg` | min | 2 | the repository's existing sliver threshold (2°) |

## 4. Comparison table (every number injected from the artefact)

| mesh | tris | worst |library − independent route| | min angle ° | max aspect_ratio | min gamma | s |
|---|---|---|---|---|---|---|
| box_unit | 12 | 2.132e-14 | 45 | 1.39385 | 0.866025 | 0.001 |
| series60 | 1998 | 3.908e-14 | 9.49395 | 3.8244 | 0.281748 | 0.054 |
| sliver_patch | 3 | 2.555e-11 | 0.143239 | 231.23 | 0.0043301 | 0.000 |
| suboff_coarse | 9088 | 2.405e-08 | 0.003865 | 8559.93 | 0.000116839 | 0.271 |
| suboff_taubin30 | 9088 | 3.553e-14 | 10.3079 | 4.93372 | 0.268713 | 0.261 |
| wigley | 1998 | 4.263e-14 | 11.5566 | 3.20676 | 0.339958 | 0.059 |

Exact arithmetic on the rational meshes:

| mesh | worst |library − exact| |
|---|---|
| box_unit | 2.132e-14 |
| sliver_patch | 2.544e-12 |

Gate verdicts, with the measure that decided them:

| mesh | verdict | failing measure (limit, measured, element) |
|---|---|---|
| box_unit | PASS | — |
| series60 | PASS | — |
| sliver_patch | FAIL | `min_angle_deg` min 2, got 0.143239 at face 2; `gamma` min 0.05, got 0.0043301 at face 2; `aspect_ratio` max 25, got 231.23 at face 2; `edge_ratio` max 100, got 400.001 at face 2 |
| suboff_coarse | FAIL | `min_angle_deg` min 2, got 0.003865 at face 4711; `gamma` min 0.05, got 0.000116839 at face 4711; `aspect_ratio` max 25, got 8559.93 at face 4711; `edge_ratio` max 100, got 9211.34 at face 4718 |
| suboff_taubin30 | PASS | — |
| wigley | PASS | — |

The raw marching-tets Suboff shell fails (min angle 0.003865°) — the documented defect the smoother exists for — and 30 Taubin sweeps take it to 10.31°, flipping the gate to PASS. A gate that cannot separate those two meshes is not a gate.

## 5. Sensitivity

Gate verdict vs `min_angle_deg` limit:

| mesh | 0.5 | 1.0 | 10.0 | 2.0 | 5.0 |
|---|---|---|---|---|---|
| box_unit | PASS | PASS | PASS | PASS | PASS |
| series60 | PASS | PASS | FAIL | PASS | PASS |
| sliver_patch | FAIL | FAIL | FAIL | FAIL | FAIL |
| suboff_coarse | FAIL | FAIL | FAIL | FAIL | FAIL |
| suboff_taubin30 | PASS | PASS | PASS | PASS | PASS |
| wigley | PASS | PASS | PASS | PASS | PASS |

Gate verdict vs `aspect_ratio` limit:

| mesh | 10.0 | 100.0 | 25.0 | 50.0 |
|---|---|---|---|---|
| box_unit | PASS | PASS | PASS | PASS |
| series60 | PASS | PASS | PASS | PASS |
| sliver_patch | FAIL | FAIL | FAIL | FAIL |
| suboff_coarse | FAIL | FAIL | FAIL | FAIL |
| suboff_taubin30 | PASS | PASS | PASS | PASS |
| wigley | PASS | PASS | PASS | PASS |

No mesh changes side between 2° and 5° of minimum angle, or between aspect ratio 25 and
50, so the defaults are not sitting on a cliff edge.

## 6. Failure modes, limits, artefacts

* **Upstream binaries absent** — VTK, Gmsh and meshmagick are all absent here, so the
  acceptance is `|library − independent route| < 1e-6`, not `< 1e-6` against gmsh. Stated,
  not papered over; the differential code for both upstream tools is written.
* **Heron loses digits where slivers are.** On the height-1e-4 triangle the float Heron
  route reports a 1e-8 relative area error, which reaches 1e-5 in the aspect ratio; the
  exact route settles it and the library's cross-product route is the accurate one.
* **`acos` conditioning.** Angles are bounded at 1e-9° rather than 1e-12: at 0.01° an
  1e-16 relative error in the cosine becomes ~1e-12 degrees. Stated, not hidden.
* **Quads are not implemented**: `Mesh` stores triangles only, so a quad measure would be
  untested dead code. `minSIGE`, `eta` and warping are likewise named as missing rather
  than assumed.

Reproduce:

```
python tools/mesh_quality_crosscheck.py           # differential + gate, writes artefacts
python -m pytest tests/test_mesh_quality_lib.py   # fast CI subset
python tools/quality_oracle.py                    # the independent route alone
python tools/make_mesh_quality_doc.py             # regenerate this page
```

Artefacts: `artifacts/mesh_quality_lib/crosscheck.json`,
`artifacts/mesh_quality_lib/gate.json`. The driver exits 1 when any tolerance is
violated, so it is the CI door itself.

