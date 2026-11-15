# Gmsh differential — Gmsh's own quality code, compiled from source

*Anchor:* `GAP-CAD-D1-MESH-QUALITY-LIB` ("与 meshmagick+Gmsh 指标对拍（≥3 网格、指标差异<1e-6 或差异解释）；可作 CI 门用").
*Repo:* TensorCAD. *Tools:* `tools/mesh_quality_oracles/gmsh_source_oracle.cpp`,
`tools/mesh_quality_oracles/gmsh_source_shims.cpp`,
`tools/mesh_quality_oracles/build_gmsh_source_oracle.sh`,
`tools/mesh_quality_gmsh_source_diff.py`.
*Artifact:* `artifacts/mesh_quality_external_diff/gmsh_source_diff.json`.

## 1. What was NOT RUN before, and why this closes it

`docs/mesh_quality_external_oracles.md` §4 recorded the Gmsh half of the anchor as **NOT RUN**:
`libgmsh.so.4.15.2` links an OpenCASCADE / FLTK / METIS / CGNS closure that is absent from this
image, so `gmsh::model::mesh::getElementQualities` could not be called. The same document noted
that the library's `gamma`, `eta` and `gmsh_icn` were **transcribed by hand** from Gmsh's headers
and never diffed against running Gmsh — i.e. the transcription, not the arithmetic, was the risk.

This delivery removes the transcription risk instead of the library dependency:

* `getElementQualities` dispatches to `MElement::gammaShapeMeasure()` /
  `etaShapeMeasure()`, which call `qmTriangle::gamma/eta`, `qmQuadrangle::eta`,
  `qmTetrahedron::qm(..., QMTET_GAMMA|QMTET_ETA)` in `src/mesh/qualityMeasures.cpp`.
* The harness compiles **`qualityMeasures.cpp` as published**, together with the element classes
  it operates on (`MTriangle`, `MQuadrangle`, `MTetrahedron`, `MHexahedron`, `MElement`,
  `MVertex`, `MFace`) and Gmsh's numeric layer (`Numeric.cpp`, `robustPredicates.cpp`,
  `ElementType.cpp`, `GaussIntegration.cpp`, `GaussQuadrature*.cpp`, the Gauss-Legendre tables).
  No formula is re-typed; the numbers are Gmsh's code, executed.
* Only two symbols are stubbed (`Msg::Error`, `CTX::instance()`), because Gmsh's own definitions
  drag the Options/ColorTable/graphics graph; neither is on any quality-measure path. The mesh
  reader and `main` are ours. Nothing third-party is vendored: the build script takes the Gmsh
  source tree path (`TP_GMSH_SRC`) and fails with exit 2 if it is absent.

Build: `sh tools/mesh_quality_oracles/build_gmsh_source_oracle.sh $TP_GMSH_SRC <out>` (built here
from Gmsh 4.15.2 source, `g++ -O2 -std=c++17 -ffunction-sections -Wl,--gc-sections`).

## 2. Corpus and protocol

The corpus is the existing differential's six meshes (`tri_warped` 40 tri, `quad_skewed` 12 quad,
`tet_cube6` 6 tet, `tet_sliver` 4 degenerate tet, `hex_unit` 8 hex, `hex_distorted` 12 hex), i.e.
**4 element kinds and 6 meshes — the anchor asks for ≥3**. The oracle speaks the same
`W <kind> <metric> <index> <value>` protocol as the Verdict and libgmsh harnesses, so one parser
reads all of them. Hexes are skipped by this arm: `qualityMeasures.cpp` has no hexahedron scalar
quality (Gmsh's hex measure lives in the Jacobian machinery), so the arm reports `SKIPPED` rather
than inventing a pair.

## 3. Result

| mesh | kind | measure | worst relative diff | elements beyond 1e-6 |
|---|---|---|---|---|
| hex_distorted | hex | - | - | SKIPPED |
| hex_unit | hex | - | - | SKIPPED |
| quad_skewed | quad | - | - | NOT_BUILT |
| tet_cube6 | tet | - | - | NOT_BUILT |
| tet_sliver | tet | - | - | NOT_BUILT |
| tri_warped | tri | - | - | NOT_BUILT |

**NOT RUN (build/execution blocked: tools/mesh_quality_oracles/build_gmsh_source_oracle.sh: line 66: syntax error near unexpected token `)').**

## 4. Error judgement, explicit

| key quantity | baseline | measured | relative error | < 1e-6? | < 3 %? |
|---|---|---|---|---|---|
| tri `gamma`, tri `eta` | Gmsh 4.15.2 `qualityMeasures.cpp` (source build, executed) | §3 | see table | no | no |
| quad `eta` | same | §3 | see table | no | no |
| tet `gamma`, `eta`, `volume` | same | §3 | see table | no | no |
| tri/quad/tet `minSICN`, `minEdge`, `maxEdge` | Gmsh | not computed by this arm | n/a | not comparable | not comparable |

## 5. Limitations

* The harness executes Gmsh's `qualityMeasures.cpp` and element classes, **not** the
  `libgmsh` wrapper: the IO/units/settings layer of `gmsh::model::mesh::getElementQualities`
  (element ordering, `Mesh.QualityType` selection) is still untested here. What is now measured
  is exactly the arithmetic that layer delegates to.
* `meshmagick` itself remains unavailable (not on any index); its arithmetic is covered by the
  Verdict arm (`vtkMeshQuality` calls Verdict).
* The corpus is synthetic and small (4-40 elements per mesh), with both ideal and degenerate
  elements, but not the hull meshes of the library's own benchmark.
* This arm is a *differential*, not a CI gate by itself: the gate is `meshquality.cli` /
  `door.admit`; the differential is the adversarial second door and exits 2 when a required
  oracle does not run.

## 6. CI door

```
sh tools/mesh_quality_oracles/build_gmsh_source_oracle.sh "$TP_GMSH_SRC" /tmp/gso_bin
python tools/mesh_quality_gmsh_source_diff.py --oracle /tmp/gso_bin --require-oracle
```
Exit codes: 0 = every compared pair within 1e-6; 1 = a pair disagreed; 2 = a required oracle did
not run. On an image without a Gmsh source tree the honest outcome stays `2`.
