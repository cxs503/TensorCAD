# Verdict oracle: the meshmagick arm of the mesh-quality acceptance criterion, run

> Mission `m22372`, anchor `GAP-CAD-D1-MESH-QUALITY-LIB` (TensorCAD).
> Numbers in this document are written by `examples/mesh_quality/verdict_oracle_crossval.py`
> into `artifacts/mesh_quality/verdict_crossval.json`; the tables here quote that artifact.

## 1. What the ledger asked for, and what was missing

The acceptance text is *"与 meshmagick+Gmsh 指标对拍（≥3 网格、指标差异 < 1e-6 或差异解释）；
可作 CI 门用"* — a cross-check against meshmagick and Gmsh on at least three meshes, with
differences below `1e-6` **or explained**, plus a usable CI gate.

`docs/mesh_quality.md` §5 recorded both oracles as *not run*, and the reason was
environmental, not methodological:

* **gmsh**: the wheel is published for `x86_64`/`aarch64` only; this host is `loongarch64`
  (`pip download gmsh` → *no matching distribution*). Still not run — see §6.
* **meshmagick**: not on the package index and GPL-3.0; its `print_quality()` delegates to
  `vtk.vtkMeshQuality`, and no VTK wheel exists for `loongarch64`. Mirroring the Verdict
  formulas a third time in-tree would not be a cross-check — it would be the same algebra
  written twice by the same hand, which is exactly what an oracle is supposed to avoid.

The gap was therefore not "the library lacks a metric"; it was **"the library is still its
own judge for the metrics that matter to a solver"**. That is what this work closes.

## 2. What was done

`vtkMeshQuality` is a wrapper; the element metrics live in **Verdict**, VTK's
`ThirdParty/verdict`. So:

1. `tools/verdict_oracle/verdict_oracle.cpp` — a C ABI wrapper over Verdict. It contains
   **no metric algebra**: cell kind → node count, metric name → Verdict function, flat
   `double[n][3]` buffer → the layout Verdict expects. Every number is Verdict's.
2. `tools/verdict_oracle/build.py` — fetches VTK at a pinned release with a sparse,
   blob-filtered clone (only `ThirdParty/verdict` is transferred), instantiates
   `verdict_config.h.in` without CMake, compiles the ten Verdict sources plus the wrapper
   into `build/libverdict_oracle.so`, and writes `provenance.json` (VTK ref + commit,
   compiler, flags, per-file and library digests). The shared library is **not** committed;
   the recipe and the provenance are.
3. `tensorcad/meshquality/verdict_oracle.py` — the Python side: `available()`,
   `provenance()`, `metric_names()`, `metrics()`. Optional by construction, like `gmsh` and
   `trimesh` in `reference.py`: nothing imports it at module load, and `available()` is
   `False` on a host without a compiler. `reference.external_availability()` now reports it.
4. `examples/mesh_quality/verdict_oracle_crossval.py` — the comparison on the meshes the
   repository already generates, with an explicit judgement per pair (below).
5. `tests/test_meshquality_verdict_oracle.py` — 8 tests; skipped (recorded as *not run*)
   when no oracle is built.

Build provenance for the run quoted here: VTK `v9.3.1`
(commit `4f02ae83179c6a1542ac3c335c73dea976578c63`), `g++ (GCC) 12.3.0`,
`libverdict_oracle.so` sha256 prefix — see `artifacts/mesh_quality/verdict_crossval.json`
field `oracle`.

## 3. How a pair is judged (no tuning to agreement)

| mode | meaning | criterion |
|---|---|---|
| `identity` | same definition in both libraries | `|a-b| <= 1e-6 * max(1,|b|)` |
| `factor` | same shape measure, two normalisations | the **ratio** must equal the exact algebraic factor |
| `relative` | same physical quantity, two constructions | reported, informational |

Comparable pairs: triangles and quadrilaterals — `area`, `scaled_jacobian`, `min_angle`,
`max_angle`, `edge_ratio`; tetrahedra — `volume`, `edge_ratio`, `scaled_jacobian` (factor);
hexahedra — `volume`, `scaled_jacobian`, `edge_ratio` (plus `volume` relative on warped
cells). Metrics both libraries compute with *different definitions* are listed in the
artifact under `different_by_definition` and are not compared: `tet.min_angle` (dihedral vs
face angle), `tet.min_sicn` (here `d|det J|^{2/d}/||J||_F^2`, Verdict `||J||_F||J^{-1}||_F/d`),
`hex.min_angle` (Verdict has no hex face-angle metric), `skewness` (Fluent equiangle skew vs
Verdict's per-kind reference angles).

## 4. Results

### 4.1 Meshes and outcomes

Nine meshes (`ideal_tri/quad/tet/hex`, `box_ideal_hex`, `box_stretched_hex`, `hull_hex`,
`hull_tet`, `hull_surface_tri`) — well above the "≥3 meshes" the ledger asks for. Every
gating pair passes:

| mesh.kind | cells | worst gating difference | mode |
|---|---:|---:|---|
| `ideal_tri.tri` | 1 | 4.44e-16 | identity |
| `ideal_quad.quad` | 1 | 0 | identity |
| `ideal_tet.tet` | 1 | 0 (`ratio` = 1/√2 to 0) | factor |
| `ideal_hex.hex` | 1 | 3.33e-16 | identity |
| `box_ideal_hex.hex` | 64 | 8.88e-16 | identity |
| `box_stretched_hex.hex` | 64 | 3.55e-15 | identity |
| `hull_hex.hex` | 3840 | 3.33e-16 (`scaled_jacobian`, `edge_ratio`) | identity |
| `hull_tet.tet` | 23040 | 1.04e-14 (`scaled_jacobian` ratio) | factor |
| `hull_surface_tri.tri` | 1998 | 2.2e-14 | identity |

**Worst gating difference over all meshes and pairs: 1.04e-14** — six orders of magnitude
below the `1e-6` the acceptance text asks for. Exit status of the script is the CI contract:
`0` when every gating pair is within `1e-6`, `1` otherwise.

### 4.2 Two explained differences (the ledger's "或差异解释")

**(a) Tetrahedron scaled Jacobian: exactly `1/√2`.** This tree normalises by the regular
tetrahedron (ideal = 1), Verdict by the orthogonal one (ideal = 1). The observed ratio
`tensorcad / Verdict` over all 21528 tetrahedra of `hull_tet` lies in
`[0.7071067811865401, 0.7071067811865528]` — `1/√2 = 0.7071067811865476` to 14 digits. The
difference is a normalisation, not a disagreement, and the test asserts the factor rather
than the values.

**(b) Verdict refuses 496 of the 3840 `hull_hex` cells.** For those cells Verdict returns
`VERDICT_DBL_MAX` (1e30), its "will not score this" sentinel, because a corner Jacobian is
non-positive. This tree instead reports the negative value, so an inverted/collapsed cell
stays visible instead of being replaced by a sentinel (`docs/mesh_quality.md` §4.5 relies on
that: the gate counts zero-volume cells that `VolumeMesh.quality()` silently drops). On the
3344 cells Verdict does score, `scaled_jacobian` and `edge_ratio` agree to **0** and
**2.2e-16**.

### 4.3 Hexahedron volume on warped cells is a different question, not a discrepancy

The gating `volume` pair is compared on the affine meshes (`ideal_hex`,
`box_ideal_hex`, `box_stretched_hex`, `hull_tet`), where it agrees to `1e-15`. On
`hull_hex`, whose cells have non-planar faces, this tree integrates the six faces
(divergence theorem, the physically meaningful cell volume) while Verdict evaluates the
trilinear Jacobian **at the cell centre**. There the two agree only for affine cells:

| quantity on `hull_hex` | value |
|---|---:|
| median relative difference | 3.9e-2 |
| worst relative difference | 5.37 |
| median absolute difference | 6.0e-6 |
| worst absolute difference | 1.27e-4 |

That is why the pair is informational rather than gating: on a warped hexahedron "the
volume" depends on what you are measuring. Anyone comparing mesh-quality volume numbers
between tools should expect this, and the artifact states it instead of hiding it.

## 5. Reproduce

```
python tools/verdict_oracle/build.py                # sparse VTK v9.3.1 clone + g++ build
python examples/mesh_quality/verdict_oracle_crossval.py
python -m pytest tests/test_meshquality_verdict_oracle.py -q      # 8 passed
```

`build.py --check` reports whether a library is present (exit 1 if not), and
`TENSORCAD_VERDICT_ORACLE=/path/to/libverdict_oracle.so` overrides the location.

## 6. Still not run, and why (honest)

* **gmsh**: still unavailable on this host. The gmsh 4.x/5.x sources *can* be fetched
  (`https://github.com/live-clones/gmsh` mirror) and CMake does configure, but the build is
  a heavy compute job on a shared 128-core node and the mission's compute discipline sends it
  to the scheduler; it was not completed in this window. `reference.gmsh_element_qualities`
  is already written against `getElementQualities('minSICN')` and runs unchanged on a host
  with gmsh, so the arm is a *host* gap, not a code gap. The Gmsh units (`gamma`, `eta`,
  `gmsh_icn`) remain verified only against their analytic reference values and against the
  in-tree second implementation.
* **meshmagick itself**: not installed (GPL-3.0, not on the index, needs VTK). What is now
  verified is the library meshmagick calls — the coverage is arguably stronger than a
  meshmagick run, which would only have re-exported these same Verdict numbers.
* `vtk.vtkMeshQuality`'s *own* wrappers (thresholds, the quality measure → array mapping)
  are not exercised: only Verdict, the metric engine, is.

## 7. Scope limits of this evidence

* Verdict's `hex_*` metrics use only the eight corner Jacobians; the comparison says nothing
  about the interior of a warped cell (see §4.3).
* The tetrahedra compared are the Kuhn split of a structured grid. Degenerate cells
  (zero measure) are excluded and counted — 1512 of 23040 in `hull_tet` — because both
  acos-based angle evaluations and the `1/√2` ratio are ill-conditioned there.
* One platform, one compiler (`g++ 12.3.0`, `loongarch64`). Floating-point agreement at
  `1e-14` is not a promise of bit-identical results under `-ffast-math` or another libm.