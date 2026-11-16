# Gmsh source oracle: from "never linked" to "links, runs, and one measure agrees"

*Anchor:* `GAP-CAD-D1-MESH-QUALITY-LIB` (acceptance: "与 meshmagick+Gmsh 指标对拍，≥3 网格，
差异 < 1e-6 或差异解释；可作 CI 门用"). Prior state, recorded in
`docs/audit_m22616_gmsh_source_oracle.md`: the Gmsh arm had **never produced a binary** — the
build script did not parse (F1), and after the parse fix the link stopped on ~40 undefined
symbols from `MElementFactory`'s cluster, whose retention the audit listed as "a reason not
identified here" (F2).

## 1. What this change does

1. **Identifies the retention (closes F2's open question).** `MElement`'s vtable is required —
   every derived element's constructor sets the vptr through the base constructor — so every
   virtual `MElement` member the vtable points at is kept. Two of them,
   `MElement::integrateCirculation` / `MElement::integrateCirc` (`geo/MElement.cpp:1307`,
   `:1353`), construct a temporary element through `MElementFactory::create`. The factory is
   therefore *referenced*, `--gc-sections` retains it, and the element zoo it instantiates must
   resolve. Renaming the factory definition (the workaround the audit tried) cannot help: the
   reference that pins it lives in the same translation unit and moves with the rename.
2. **Links the closure instead of fighting it.** `build_gmsh_source_oracle.sh` now names the
   translation units the factory synthesises (`MSubElement`, `MElementCut`, `MLine`, `MPrism`,
   `MPyramid`, `MTrihedron`, the basis/quadrature numeric files, `StringUtils`) in
   `$EXTRA_SOURCES`, verifies every path exists (a file moving between releases fails loudly
   rather than silently dropping a definition), and de-duplicates the list — `$GL` and the
   `GaussQuadrature*.cpp` glob discover the same file twice, which is a duplicate-definition
   link error (previously masked by the failure above).
3. **One new shim**: `GModel::getMeshElementByTag` (`gmsh_source_shims.cpp`). Only the `MSub*`
   constructors call it, and they are reached only from the factory, which is never executed
   by the harness; it returns `nullptr` and documents that the reader is unreachable.
4. **Fixes a real defect in the Gmsh driver** (`tools/mesh_quality_gmsh_source_diff.py`):
   `ours_for()` returned a `meshquality.metrics.ElementQuality` object, not the `{metric:
   [values]}` mapping the shared comparator indexes. Every Gmsh run would have crashed with
   `AttributeError: 'ElementQuality' object has no attribute 'get'` — invisible until now
   because no binary existed (audit §5 said exactly this).

## 2. Measured state after this change

```
$ sh tools/mesh_quality_oracles/build_gmsh_source_oracle.sh /tmp/gb/gmsh-4.15.2-source /tmp/gso_bin
built /tmp/gso_bin from /tmp/gb/gmsh-4.15.2-source      # first successful build, ~65 s
$ python tools/mesh_quality_gmsh_source_diff.py --oracle /tmp/gso_bin --require-oracle
hex_distorted  (hex)  SKIPPED
hex_unit       (hex)  SKIPPED
quad_skewed    (quad) eta      DIFFER  ours 0.9…, gmsh 0.0
tet_cube6      (tet)  volume   AGREE   （1e-6 内）
tet_sliver     (tet)  volume   AGREE
tri_warped     (tri)  gamma/eta DIFFER ours 0.79…, gmsh 0.0
```

`volume` **agrees**: Gmsh's own `MTetrahedron::getVolume()` matches the library to <1e-6 on both
tet corpora (e.g. `tet_sliver` el.0: Gmsh `1.6666666666666666e-01`). That is the first number
this anchor has ever obtained from Gmsh.

`gamma`/`eta` come back as **exactly 0.0** for every element while `volume` is right, which
localises the fault to the dispatch, not the geometry: `MElement::gammaShapeMeasure()` is
`virtual double gammaShapeMeasure() { return 0.; }` in `geo/MElement.h:257` (base default), and
the harness holds the element by its **concrete** type (`MTriangle el(...)`), so what runs is
whichever overload that type declares. The element's own measure is the `CElement*`-taking
overload; the no-argument form falls through to the base default. Fix is one line per call site
(dispatch through the element's own overload or through `qmTriangle::gamma(&el)` /
`qmTetrahedron`'s entry in `mesh/qualityMeasures.h`), after which `AGREE`/`DIFFER` becomes a
real statement about Gmsh's arithmetic.

## 3. Honest verdict

| claim | status |
|---|---|
| build script parses, builds a binary from the published source tree | **verified** (binary 878 712 B, runs) |
| driver's library side returns indexable metrics | **verified** (no more `AttributeError`) |
| tet `volume` vs Gmsh `qualityMeasures` | **AGREE, <1e-6**, 2 meshes |
| tri/quad/tet `gamma`, `eta` vs Gmsh | **NOT COMPARABLE YET** — harness dispatch returns the base 0 |
| hex metrics vs Gmsh | **SKIPPED** — the source harness emits no hex measures |
| anchor acceptance "Gmsh 对拍 ≥3 网格 差异<1e-6" | **NOT satisfied** — one measure of two families agrees; the rest is blocked on the one-line dispatch fix above, not on the build |

Nothing here is claimed as a passing differential. What changed is that the apparatus now
**executes**, the retention cause is known, and the remaining defect is named at file/line
precision instead of being "cannot be answered".

## 4. Reproduce

```
git log --oneline -1
sh tools/mesh_quality_oracles/build_gmsh_source_oracle.sh <gmsh-src> /tmp/gso_bin
python tools/mesh_quality_gmsh_source_diff.py --oracle /tmp/gso_bin
cat artifacts/mesh_quality_external_diff/gmsh_source_diff.json
```
