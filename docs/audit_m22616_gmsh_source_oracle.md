# Audit — the Gmsh *source* oracle: declared, never executed

*Anchor:* `GAP-CAD-D1-MESH-QUALITY-LIB` (mesh-quality metric library; acceptance item
"与 meshmagick+Gmsh 指标对拍，≥3 网格，差异 < 1e-6 或差异解释").
*Reviewed artifacts:* `tools/mesh_quality_oracles/gmsh_source_oracle.cpp`,
`gmsh_source_shims.cpp`, `build_gmsh_source_oracle.sh`,
`tools/mesh_quality_gmsh_source_diff.py`, `docs/mesh_quality_external_oracles.md` §2/§4.
*Independence limitation:* the Gmsh half could not be executed, so nothing here says the
library disagrees with Gmsh — only that the apparatus meant to answer that question never
ran. The Verdict half of the differential (6 meshes, `AGREE` ≤ 2.4e-15) is not re-audited
here.

## 1. Finding F1 (blocking) — the build script does not parse

```
$ sh tools/mesh_quality_oracles/build_gmsh_source_oracle.sh /tmp/gb/gmsh-4.15.2-source
tools/mesh_quality_oracles/build_gmsh_source_oracle.sh: line 66: syntax error near unexpected token `)'
```

Line 66 was an editing fragment left after the `GL=` assignment — a stray argument list
(`"$G/numeric" "$G/mesh" "$G/common" … | tr '\n' ' ')`) with no command and an unbalanced
quote in `tr '\\n' ' '`. Consequence: the script exits at *parse* time, before any
compilation, so no oracle binary has ever been produced by it — the `bin/` directory it
defaults to does not exist. Fixed in this change.

## 2. Finding F2 (blocking) — after the parse fix, the link fails; the failure is in code
that is not on the quality path

With the syntax fixed, `g++` compiles every Gmsh translation unit and the harness, and the
link stage reports ~40 undefined symbols, all of one origin: `MElementFactory::create`
(`src/geo/MElement.cpp`), the constructor dispatch for *every* Gmsh element type. Its call
graph reaches the high-order and sub-element machinery — `MSubPoint`, `MSubLine`,
`MSubTriangle`, `MSubTetrahedron`, `MPrism15`, `MLineN`, `MPolygon`/`MPolygonBorder`
(`isInside`, `getIntegrationPoints`, typeinfo), `BasisFactory`, `JacobianBasis`,
`FuncSpaceData`, `GradientBasis`, `MEdgeN`/`MPolyhedron` — none of which the measures below
touch:

* the harness builds its elements with the concrete constructors it needs
  (`MTriangle`, `MQuadrangle`, `MTetrahedron`), never through the factory;
* the measures called are `gammaShapeMeasure`, `etaShapeMeasure`, `getVolume` from
  `src/mesh/qualityMeasures.cpp` and the element classes.

So the intended `-ffunction-sections`/`--gc-sections` isolation was expected to drop the
factory. Measured: it does not. Isolated, `MElement.cpp` + an empty `main` **does** drop the
factory section; in the full link it is retained. The audit stops at the measurement: the
factory's `.rodata` neighbours (`_ZTV14MPolygonBorder`, `_ZTI14MPolygon`) survive for a
reason not identified here, and every unresolved symbol hangs off that retained cluster.
A workaround tried and *not* kept: compiling `MElement.cpp` with
`-DMElementFactory=MElementFactoryNotOnQualityPath` so the definition is locally renamed
and unreferenced — the same cluster is still retained and still unresolved, which is what
rules out "a stray caller in the harness" (there is none; the harness contains no
`Factory::create` reference at all).

## 3. Finding F3 (honesty defect) — documentation claims compilation that never happened

`docs/mesh_quality_external_oracles.md` §4 states: "Harness written **and compiled**:
`tools/mesh_quality_oracles/gmsh_oracle.cpp` …". The harness cannot have been compiled,
because its build script never parsed (§1). §4's blocker paragraph then correctly explains
why *`libgmsh`* cannot be loaded, so the reader is left with the impression that the
source-compilation route was tried and failed when in fact it was never run. This change
keeps the `NOT RUN` verdict and replaces the "and compiled" claim with the measured
status.

## 4. What is fixed here, and what is still open

Fixed, and verified by the commands quoted above:

1. the parse error in `build_gmsh_source_oracle.sh` (posix `sh -n` clean);
2. the build split into two phases (compile `MElement.cpp` separately, then link) so the
   isolation question above is answerable rather than implicit;
3. `GModel::current(int)` stubbed in `gmsh_source_shims.cpp` — `MVertex`'s constructor
   (genuinely on the executed path) asks the current model for vertex numbering through
   inline accessors (`GModel.h:224`, `:236`); gmsh's own definition drags
   `GModel.cpp`/Mesh/Geo. The stub hands back a zero-initialised model image, and only
   inline bookkeeping is touched;
4. `tests/test_mesh_quality_gmsh_source_build_script.py` — pins the two honest behaviours
   that are testable on this host: the script parses under `sh -n`, and with no source
   tree it exits 2 with a message instead of pretending to build.

Still open, and *not* claimed as done: the Gmsh differential itself. The remaining work is
mechanical and is now located precisely — either link the factory's closure
(`geo/MSubElement.cpp`, `MElementCut.cpp`, `BasisFactory`, `JacobianBasis`, …) or break the
retention, then run

```
python tools/mesh_quality_gmsh_source_diff.py --oracle <bin>
```

which already maps `gamma`, `eta`, tet `volume` and compares at 1e-6. Until then the
anchor's "vs Gmsh" column stays **NOT RUN**, exactly as `docs/mesh_quality_external_oracles.md`
§4/§6 record it, and the acceptance item is **not** satisfied by this change.

## 5. Checks not run

* No oracle executed ⇒ no Gmsh number appears anywhere in this document.
* The Verdict oracle was not rebuilt or re-run (its recorded artifact is unchanged).
* `tools/mesh_quality_gmsh_source_diff.py` was not exercised end-to-end (no binary), so its
  parser/mapping is covered only by its own unit test, not by a measurement.