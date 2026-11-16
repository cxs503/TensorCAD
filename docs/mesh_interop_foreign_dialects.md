# Foreign mesh dialects: binary STL / MSH 2.2 / VTK XML against a real third-party reader

Anchor: `GAP-CAD-B2-MESH-FORMAT-MATRIX` — *meshio 全格式矩阵进出（STL/OBJ/VTK/SU2/CGNS/Gmsh）*,
acceptance: **往返保真：体积/拓扑/取向不变（6 格式 × 3 几何 = 18 组合全绿表）**.

Status of the acceptance: **met** by `docs/mesh_format_matrix_thirdparty/`
(18 combinations green, later 6 × 4 = 24) and independently verified in
`merge #664 (m22715)`. This document does **not** redo that table. It closes the
two gaps that the verification report itself declares, and nothing else.

## 1. The stated ceiling of the previous evidence

`docs/mesh_format_matrix_thirdparty/verification_report.md` §5:

> **No external CAD/solver tool was executed**: `gmsh`, Salome, ParaView are not
> installed in this worktree, so "another tool can read it" is evidenced by
> spec-conformance of the bytes plus independent parsing, not by running that
> tool. That is the honest ceiling of this verification.

> **Scope**: only the dialects `write_meshio` actually emits are parsed (ASCII
> STL, legacy VTK 5.1 binary, MSH 4.1 binary, ASCII SU2, HDF5 CGNS). Binary STL,
> ASCII MSH 4.1, MSH 2.2, VTK XML (`.vtu`) are **not** verified here and the
> readers raise `IndependentReadError` rather than guessing.

Two consequences, both material for a ship-design exchange:

1. **Ingest was single-dialect.** Every CAD system (CATIA, SolidWorks, Fusion,
   FreeCAD, Rhinoceros) writes **binary** STL, and most FVM decks ship **MSH 2.2
   ASCII**; ParaView's default unstructured container is **`.vtu`**. A hull that a
   real CAD tool exports could not enter TensorCAD at all — the reader refused a
   file the rest of the world considers valid.
2. **No foreign implementation ever ran.** Spec-conformance of our own bytes, read
   by our own parser, cannot catch a shared misreading of the specification.

## 2. What was added

| module | role |
|---|---|
| `tensorcad/mesh_foreign_dialects.py` | readers **and** spec-level writers for `stl_binary`, `gmsh22`, `vtu`; byte-identity facts; the campaign |
| `tensorcad/mesh_interop_thirdparty.py` | trimesh oracle: `trimesh_probe` (third-party reader) / `trimesh_write` (third-party writer); never imports meshio |
| `tensorcad/mesh_interop_independent.py` | dispatch extended with the three dialects (lazy, no import cycle) |
| `tests/test_mesh_interop_foreign_dialects.py` | 59 passed / 9 skipped (skips are "this tool cannot produce that dialect", stated not assumed) |
| `scripts/mesh_interop_foreign_matrix.py` | regenerates `runs/m22862_foreign_dialects/` + the figure |

Dialects parsed, all from the specification (no meshio/h5py/trimesh in the reader):

* **binary STL** — 80-byte header, `uint32` facet count, 50-byte records,
  float32 positions, per-facet normals kept and cross-checked against the winding;
* **Gmsh MSH 2.2 ASCII** — `$Nodes`/`$Elements`, element-type widths, tag
  continuity `1..N`, declared-vs-present counts;
* **VTK XML `.vtu`** — ascii, inline-base64 and **zlib-compressed multi-block**
  inline binary (header `3 + nr_blocks` integers, blocks concatenated).

The third-party tool is `trimesh 5.0.0` (its own STL/OBJ loaders and
its own mass properties), executed in both directions: it reads TensorCAD-written
files, and *its* files are read back by our readers.

## 3. Results

Campaign: `python scripts/mesh_interop_foreign_matrix.py`,
`rows_ok = 24/24`, `consumer_ok = 9/9`,
figures `artifacts/m22862_foreign_dialects/foreign_dialects.png`, JSON + written meshes + per-file sha256 in
`runs/m22862_foreign_dialects/`.

### 3.1 Independent reader on foreign dialects (24 rows)

| geometry | dialect | cells | producer | consumer | cells | V [m^3] | V ref [m^3] | rel. err | tol | winding | count | ok |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| cube | stl_binary | tri | tensorcad | independent reader | 12/12 | 1.000000000000 | 1.000000000000 | 0.00e+00 | 1e-06 | PASS | PASS | PASS |
| cube | stl_binary | tri | trimesh | independent reader | 12/12 | 1.000000000000 | 1.000000000000 | 0.00e+00 | 1e-06 | PASS | PASS | PASS |
| cube | gmsh22 | tri | tensorcad | independent reader | 12/12 | 1.000000000000 | 1.000000000000 | 0.00e+00 | 1e-12 | PASS | PASS | PASS |
| cube | gmsh22 | tri | meshio | independent reader | 12/12 | 1.000000000000 | 1.000000000000 | 0.00e+00 | 1e-12 | PASS | PASS | PASS |
| cube | vtu | tri | tensorcad | independent reader | 12/12 | 1.000000000000 | 1.000000000000 | 0.00e+00 | 1e-12 | PASS | PASS | PASS |
| cube | vtu | tri | meshio | independent reader | 12/12 | 1.000000000000 | 1.000000000000 | 0.00e+00 | 1e-12 | PASS | PASS | PASS |
| cube | gmsh22 | tetra | tensorcad | independent reader | 12/12 | 1.000000000000 | 1.000000000000 | 0.00e+00 | 1e-12 | PASS | PASS | PASS |
| cube | vtu | tetra | tensorcad | independent reader | 12/12 | 1.000000000000 | 1.000000000000 | 0.00e+00 | 1e-12 | PASS | PASS | PASS |
| capsule | stl_binary | tri | tensorcad | independent reader | 1536/1536 | 3.639460513058 | 3.639460644635 | 3.62e-08 | 1e-06 | PASS | PASS | PASS |
| capsule | stl_binary | tri | trimesh | independent reader | 1536/1536 | 3.639460513058 | 3.639460644635 | 3.62e-08 | 1e-06 | PASS | PASS | PASS |
| capsule | gmsh22 | tri | tensorcad | independent reader | 1536/1536 | 3.639460644635 | 3.639460644635 | 0.00e+00 | 1e-12 | PASS | PASS | PASS |
| capsule | gmsh22 | tri | meshio | independent reader | 1536/1536 | 3.639460644635 | 3.639460644635 | 0.00e+00 | 1e-12 | PASS | PASS | PASS |
| capsule | vtu | tri | tensorcad | independent reader | 1536/1536 | 3.639460644635 | 3.639460644635 | 0.00e+00 | 1e-12 | PASS | PASS | PASS |
| capsule | vtu | tri | meshio | independent reader | 1536/1536 | 3.639460644635 | 3.639460644635 | 0.00e+00 | 1e-12 | PASS | PASS | PASS |
| capsule | gmsh22 | tetra | tensorcad | independent reader | 1536/1536 | 3.639460644635 | 3.639460644635 | 0.00e+00 | 1e-12 | PASS | PASS | PASS |
| capsule | vtu | tetra | tensorcad | independent reader | 1536/1536 | 3.639460644635 | 3.639460644635 | 0.00e+00 | 1e-12 | PASS | PASS | PASS |
| wigley_hull | stl_binary | tri | tensorcad | independent reader | 1998/1998 | 2.771222144484 | 2.771222149884 | 1.95e-09 | 1e-06 | PASS | PASS | PASS |
| wigley_hull | stl_binary | tri | trimesh | independent reader | 1998/1998 | 2.771222144484 | 2.771222149884 | 1.95e-09 | 1e-06 | PASS | PASS | PASS |
| wigley_hull | gmsh22 | tri | tensorcad | independent reader | 1998/1998 | 2.771222149884 | 2.771222149884 | 0.00e+00 | 1e-12 | PASS | PASS | PASS |
| wigley_hull | gmsh22 | tri | meshio | independent reader | 1998/1998 | 2.771222149884 | 2.771222149884 | 0.00e+00 | 1e-12 | PASS | PASS | PASS |
| wigley_hull | vtu | tri | tensorcad | independent reader | 1998/1998 | 2.771222149884 | 2.771222149884 | 0.00e+00 | 1e-12 | PASS | PASS | PASS |
| wigley_hull | vtu | tri | meshio | independent reader | 1998/1998 | 2.771222149884 | 2.771222149884 | 0.00e+00 | 1e-12 | PASS | PASS | PASS |
| wigley_hull | gmsh22 | tetra | tensorcad | independent reader | 1996/1996 | 2.771222149884 | 2.771222149884 | 0.00e+00 | 1e-12 | PASS | PASS | PASS |
| wigley_hull | vtu | tetra | tensorcad | independent reader | 1996/1996 | 2.771222149884 | 2.771222149884 | 0.00e+00 | 1e-12 | PASS | PASS | PASS |

`cells` in the third column is the cell kind (`tri` triangles, `tetra` tetrahedra);
the last six rows of each geometry block carry the tetrahedralised mesh, which is
what SU2/CGNS-style volume exchange needs.

### 3.2 Third-party implementation (trimesh) reading TensorCAD files (9 rows)

| geometry | dialect | trimesh volume [m^3] | reference [m^3] | rel. err | watertight | winding consistent | faces | ok |
|---|---|---|---|---|---|---|---|---|
| cube | stl | 1.000000000000 | 1.000000000000 | 0.00e+00 | True | True | 12 | True |
| cube | stl_binary | 1.000000000000 | 1.000000000000 | 0.00e+00 | True | True | 12 | True |
| cube | obj | 1.000000000000 | 1.000000000000 | 0.00e+00 | True | True | 12 | True |
| capsule | stl | 3.639460644635 | 3.639460644635 | 1.22e-16 | True | True | 1536 | True |
| capsule | stl_binary | 3.639460513058 | 3.639460644635 | 3.62e-08 | True | True | 1536 | True |
| capsule | obj | 3.639460644635 | 3.639460644635 | 1.22e-16 | True | True | 1536 | True |
| wigley_hull | stl | 2.771222149884 | 2.771222149884 | 1.60e-16 | True | True | 1998 | True |
| wigley_hull | stl_binary | 2.771222144484 | 2.771222149884 | 1.95e-09 | True | True | 1998 | True |
| wigley_hull | obj | 2.771222149884 | 2.771222149884 | 1.60e-16 | True | True | 1998 | True |

## 4. Tolerances, and why they are per dialect

| dialect | precision | tolerance | justification |
|---|---|---|---|
| `gmsh22`, `vtu` | text/float64, 17 significant digits | `1e-12` | exact representation; measured error is **0.0** |
| `stl_binary` | float32 positions (the format has no float64) | `1e-06` | `~3·eps32`; measured `3.6e-8` (capsule), `2.0e-9` (Wigley) |

The binary-STL bound is **used, not decorative**: a test asserts
`0 < rel_err < 1e-6` on the capsule and the Wigley hull, i.e. float32 really does
cost precision and the tolerance is not hiding it. The cube is exactly
representable and gives `0.0`, which is also asserted.

## 5. Bugs found while building this (both real, both now regression-tested)

1. **Silent short read in the VTU decoder.** The first implementation assumed
   "one base64 chunk per compressed block, each with a size prefix". The actual
   layout (VTK writer and meshio 5.x agree) puts *all* sizes in one header of
   `3 + nr_blocks` integers and concatenates the compressed blocks afterwards.
   The wrong assumption did not raise — it decoded an 8-point mesh to **2
   values**, i.e. a smaller mesh that still "had" a volume. Every length is now
   cross-checked against the header, and
   `test_vtu_rejects_short_compressed_payload` reproduces the failure mode.
2. **`meshio` has no binary STL writer.** `file_format="stl_binary"` is not a
   valid meshio format, so the "our writer vs their writer" comparison for the
   most common CAD dialect required implementing the binary writer from the
   specification *and* using trimesh as the second producer.

A methodological point that the matrix keeps demonstrating: cube facet 0 lies in
the plane `z = 0`, which passes through the origin, so flipping it changes the
divergence-theorem volume by **exactly zero**. Two tests
(`test_volume_alone_is_blind_to_a_facet_through_the_origin`,
`test_harness_detects_flipped_facet_in_vtu`) assert `ΔV == 0` *and* that the
per-cell winding signature still catches the flip. A volume-only gate would pass
a no-longer-closed mesh.

## 6. Residual gaps (NOT_SUPPORTED here, stated rather than guessed)

| gap | status |
|---|---|
| ASCII MSH 4.1, binary MSH 2.2 | not parsed (`IndependentReadError`); binary MSH 4.1 is covered by the 6-format matrix |
| VTU `format="appended"` (raw/appended-raw), LZMA compressor | not parsed, explicit error |
| quad/hex/polyhedral cells | rejected, not silently converted (`meshio_exchange` is tri/tet only) |
| binary STL topology | triangle soup: facet-level topology and winding only, vertex identity is not preserved by the format |
| third-party oracle for VTK/Gmsh/SU2/CGNS | trimesh reads neither; no installed tool reads them, so those four formats still rest on spec-conformance + independent parsing only |
| PLY/OFF/tetgen/med/Nastran | out of scope of the anchor |

## 7. Error verdict

| quantity | reference | worst observed | criterion | verdict |
|---|---|---|---|---|
| volume, text/float64 dialects (16 rows) | reference mesh volume | **0.0** | < 1e-12 | PASS |
| volume, binary STL (8 rows) | reference mesh volume | 3.62e-8 (float32) | < 1e-6 | PASS |
| volume, trimesh reading our files (9 rows) | reference mesh volume | 3.62e-8 (binary STL), 1.6e-16 (float64) | per dialect | PASS |
| cell count | reference count | 0 cells lost | exact | PASS |
| per-cell winding | winding signature | identical multiset | exact | PASS |
| watertightness (trimesh's own test) | — | 9/9 watertight, 9/9 winding-consistent | exact | PASS |
| volume vs analytic (capsule) | 3.665191429 | 0.70 % | < 1 % | PASS (discretisation) |
| volume vs analytic (Wigley L=10, B=1, T=0.625) | 2.777778 | 0.24 % | < 1 % | PASS (discretisation) |

**Is the key quantity within 3 %?** Yes. Exchange error is exactly **0** for every
float64 dialect, and `3.6e-8` for binary STL (float32 ceiling, 5 orders of
magnitude inside the 1 % budget). The only deviations from an analytic value are
0.70 % (capsule) and 0.24 % (Wigley) — faceting of the fixtures, identical in all
dialects and producers, i.e. discretisation, not exchange.

## 8. Reproduction

```
python -m pytest tests/test_mesh_interop_foreign_dialects.py -q   # 59 passed, 9 skipped
python -m pytest tests/test_mesh_format_matrix_thirdparty.py tests/test_meshio_format_matrix.py -q  # 195 passed, 4 skipped
python scripts/mesh_interop_foreign_matrix.py                     # regenerates table, JSON, figure
```
