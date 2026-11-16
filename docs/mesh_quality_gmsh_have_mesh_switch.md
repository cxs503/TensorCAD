# The Gmsh source arm now measures: a config switch, not a dispatch bug

*Anchor:* `GAP-CAD-D1-MESH-QUALITY-LIB` (acceptance: "与 meshmagick+Gmsh 指标对拍（≥3 网格、
指标差异<1e-6 或差异解释）；可作 CI 门用").
*Artifacts under review:* `tools/mesh_quality_oracles/build_gmsh_source_oracle.sh`,
`gmsh_source_oracle.cpp`, `gmsh_source_shims.cpp`, `tools/mesh_quality_gmsh_source_diff.py`,
`artifacts/mesh_quality_external_diff/gmsh_source_diff.json`,
`tests/test_mesh_quality_gmsh_source_differential.py`.
*Prior record:* `docs/audit_m22616_gmsh_source_oracle.md` (no binary), `docs/m22768_gmsh_source_oracle_build.md`
(binary links, tet `volume` agrees, everything else "NOT COMPARABLE YET").

## 1. The previous diagnosis was wrong, and it is recorded as wrong

`docs/m22768_gmsh_source_oracle_build.md` §2 concluded that `gamma`/`eta` came back as exactly
`0.0` because the harness "holds the element by its **concrete** type" and "the no-argument form
falls through to the base default", with the fix being "one line per call site" to call
`qmTriangle::gamma(&el)` directly.

That cannot be the cause. `MTriangle`, `MQuadrangle` and `MTetrahedron` **do** override the
virtual (`geo/MTriangle.h:56-57`, `geo/MQuadrangle.h:60-61`, `geo/MTetrahedron.h:138,141`;
definitions in `geo/MTriangle.cpp:106,115`, `geo/MQuadrangle.cpp:314,323`,
`geo/MTetrahedron.cpp:94,104`), and the no-argument `el.gammaShapeMeasure()` is exactly the call
Gmsh itself makes (`common/gmsh.cpp:2433-2436` in `getElementQualities`). Re-dispatching through
`qmTriangle::gamma(el)` would have run the same function and still returned 0.

The actual cause is one line earlier in the build:

```
$ sed -n 17p src/geo/MTriangle.cpp
#if defined(HAVE_MESH)
$ grep -n HAVE_MESH src/common/GmshConfig.h.in
37:#cmakedefine HAVE_MESH
```

Every shape measure in Gmsh is `#if defined(HAVE_MESH) ... #else return 0.;`:
`MTriangle::gammaShapeMeasure/etaShapeMeasure`, `MQuadrangle::eta/gamma`,
`MTetrahedron::qm`, and `MElement::signedInvCondNumRange` / `signedInvGradErrorRange`
(i.e. `minSICN` and `minSIGE` too). The build script generated its own `GmshConfig.h` by
turning `#cmakedefine X` into `/* undef X */`, so **`HAVE_MESH` was never defined** and every
measure took the `return 0.;` branch. `volume` agreed because `MTetrahedron::getVolume()` is
unconditional — which is precisely why the failure looked like a dispatch problem.

That is the worst failure mode an oracle has: a binary that runs, exits 0, and answers zero for
everything. So the switch is now set explicitly and the reason is written next to it:

```sh
DEFINE_ON="HAVE_MESH"          # every shape measure is guarded by it; undefined => silent zeros
```

Defining it pulls three more translation units into the link (`numeric/CondNumBasis.cpp`,
`mesh/qualityMeasuresJacobian.cpp`, `geo/gmshQuadrics.cpp` for `tetcircumcenter`) and one more
message shim (`Msg::Info`, reached from `_getQualityFunctionSpace`), all named in
`$EXTRA_SOURCES` / `gmsh_source_shims.cpp` with existence checks.

## 2. Measured result — Gmsh's arithmetic vs this library

`sh tools/mesh_quality_oracles/build_gmsh_source_oracle.sh /tmp/gb/gmsh-4.15.2-source /tmp/gso_bin`
→ `built /tmp/gso_bin` (≈90 s, Gmsh 4.15.2 source). Then
`python tools/mesh_quality_gmsh_source_diff.py --oracle /tmp/gso_bin --require-oracle` (exit 0):

| mesh | kind | elements | measure | worst relative difference | verdict |
|---|---|---:|---|---:|---|
| `tri_warped` | tri | 40 | `gamma` | 5.157e-16 | AGREE |
| `tri_warped` | tri | 40 | `eta` | 5.951e-16 | AGREE |
| `quad_skewed` | quad | 12 | `eta` | 1.185e-16 | AGREE |
| `quad_skewed` | quad | 12 | `area` (= Gmsh `getVolume()` on a quad) | 3.331e-16 | AGREE |
| `tet_cube6` | tet | 6 | `gamma` | 0 | AGREE |
| `tet_cube6` | tet | 6 | `eta` | 0 | AGREE |
| `tet_cube6` | tet | 6 | `volume` | 0 | AGREE |
| `tet_sliver` | tet | 4 | `gamma` | 1.853e-16 | AGREE |
| `tet_sliver` | tet | 4 | `eta` | 2.844e-16 | AGREE |
| `tet_sliver` | tet | 4 | `volume` | 1.367e-16 | AGREE |
| `hex_unit` | hex | 8 | — | — | SKIPPED (§3) |
| `hex_distorted` | hex | 12 | — | — | SKIPPED (§3) |

Ten metric groups over four meshes, 134 element-measure pairs, worst relative difference
5.95e-16 against the 1e-6 acceptance threshold — agreement at double-precision round-off, on the
degenerate corpus too (`tet_sliver` element 0: `gamma` = 7.14285669966799e-08 on both sides;
`eta` = 0.00304990741724199 on both sides). Sample values at the worst element are recorded per
metric in the artifact, so nothing here is a transcribed number.

This is the first time this anchor has obtained Gmsh's own `gamma`/`eta` values. It is a
statement about **Gmsh 4.15.2's arithmetic**, not about `meshmagick`; the `meshmagick`/VTK half of
the acceptance stays with the Verdict oracle (`tools/mesh_quality_oracles/verdict_oracle.cpp`,
6 meshes, AGREE ≤ 2.4e-15), because `meshmagick.print_quality()` delegates to
`vtk.vtkMeshQuality`, and neither a VTK wheel nor a `meshmagick` distribution exists for
loongarch64.

## 3. What Gmsh does **not** have, measured rather than assumed

Hexahedra are where the acceptance criterion cannot be met against Gmsh, and the driver now says
so with numbers instead of silence:

* `MHexahedron` overrides neither `gammaShapeMeasure()` nor `etaShapeMeasure()` → Gmsh's own
  `getElementQualities("gamma")` returns **0** for hexes (it is the base `return 0.;`).
* `MElement::getVolume()` returns **0** for a hexahedron (no override; the base implementation
  does not integrate hexes) — on `hex_unit`, whose true volume is 0.125 per cell, Gmsh says 0.
* `minSICNShapeMeasure()` / `minSIGEShapeMeasure()` return **0** / **DBL_MAX** on the same cells.

So the hexagonal arm's external reference remains the Verdict/VTK definitions, and the Gmsh arm
of the gate covers tri/quad/tet only. `minSICN` is likewise excluded from the Gmsh mapping: with
`HAVE_MESH` now defined, Gmsh's `minSICNShapeMeasure()` returns `nan` on the tri/quad corpora and
`0` on the tet corpora (`CondNumBasis` path on elements built outside a Gmsh model), so there is
no Gmsh number to compare — our `min_sicn`/`gmsh_icn` remain backed by the analytic reference
cases. Both exclusions are in `NOT_COMPARABLE` in the driver and in every emitted artifact.

## 4. The gate could never fail (blocking, fixed)

The differential's whole value is the exit code. It was structurally always 0:

* the loop advanced only on `info.get("status") != "COMPARED"`, but `compare()` emits
  `AGREE`/`DIFFER`/`NOT_RETURNED`/`LENGTH_MISMATCH`/`EMPTY` — `"COMPARED"` never occurs, so
  `pairs_compared` stayed **0** even when every pair agreed;
* the disagreement counter read `info.get("beyond_tolerance", 0)` while `compare()` writes
  `n_beyond_tol`, so the count was always 0 → `return 1 if disagreed else 0` could not return 1.

Both are fixed, and `tests/test_mesh_quality_gmsh_source_differential.py` now pins them: a fake
oracle that lies with `0.0` must produce exit **1** with `elements_beyond_tolerance > 0` and
every compared metric `DIFFER`; an absent oracle must be exit **2**, never 0. (Same file pins
`DEFINE_ON="HAVE_MESH"` in the build script and the harness's metric set, so neither defect can
come back silently.)

## 5. Reproduce

```sh
git clone --depth 1 --branch gmsh_4.15.2 https://gitlab.onelab.info/gmsh/gmsh.git /tmp/gb/gmsh-4.15.2-source
sh tools/mesh_quality_oracles/build_gmsh_source_oracle.sh /tmp/gb/gmsh-4.15.2-source /tmp/gso_bin
python tools/mesh_quality_gmsh_source_diff.py --oracle /tmp/gso_bin --require-oracle   # exit 0
python -m pytest tests/test_mesh_quality_gmsh_source_differential.py -q
```

## 6. Limits, stated

* The comparison is against Gmsh **4.15.2 source**, compiled by a harness that links Gmsh's own
  `qualityMeasures.cpp`; the mesh reader and `main` are ours, and two unreachable symbols are
  stubbed (`CTX::instance`, `Msg::*`, plus `tetcircumcenter`, which now answers NaN so an
  accidental future use is visible rather than plausible). Element construction goes through
  Gmsh's real element classes, so `MVertex` bookkeeping is Gmsh's.
* Four meshes, four element kinds, 134 element-measure pairs — enough for the acceptance's "≥3
  meshes", not a statistical claim about arbitrary real-world meshes.
* Hexagonal metrics have no Gmsh counterpart (§3); `minSICN`/`minSIGE` are not usable from a
  Gmsh source build on free-standing elements (§3).
* Pre-existing red test found while running the suite, unrelated to this change and left alone:
  `tests/test_meshquality.py::test_external_availability_reports_versions` expects exactly
  `{"gmsh","meshmagick","trimesh"}` while the availability report now also carries
  `verdict_oracle` (`tests/test_meshquality.py:296`).
