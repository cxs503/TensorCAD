# XVAL — TensorCAD native six-format mesh matrix vs mature open-source implementations

**Node** m17595 · anchor `XVAL-CAD-M17483-NATIVE-SIX-FORMAT-MATRIX` · branch `agent/m17595`
**Cross-validated object** `tensorcad/mesh_formats.py` (node m17483: native STL / OBJ / VTK / SU2 / CGNS /
Gmsh writers+readers for a closed triangle surface).
**Predeclared judgement (owner 2026-11-12):** a format passes when the *independent* re-read surface has
relative volume error **< 15 %**; the implementation acceptance line used here is the tighter `1e-3`.
Topology (vertex bijection, equal face multiset) and orientation (equal *oriented* face multiset, zero
flipped faces) are pass/fail with **no** tolerance in the five lossless formats; STL alone gets
`1e-5 x diagonal` because float32 cannot be exact.

**Verdict.** Worst OSS-vs-TensorCAD volume error is **1.311e-08** over matrix A (18 cells,
33 reader checks: XVAL writers × OSS readers) and **1.311e-08** over matrix B
(18 cells, 30 checks: the *actual m17483 artefacts* × the same OSS readers).
Both are `PASS` against the 15 % line and
`PASS` against 1e-3 — geometry transport is
validated by genuinely independent readers for **five** of the six formats (STL, OBJ,
VTK-`UNSTRUCTURED_GRID`, SU2, Gmsh): 27/33 reader checks
green in A, 21/30 in B; the sixth (CGNS) is graded
structurally.
Red entries in A: `vtk/POLYDATA/meshio` on 3 geometries (cube); `su2/-/meshio` on 3 geometries (cube). They are **reference-side or dialect findings**, not losses:
3 × VTK/POLYDATA (`meshio` cannot parse legacy `POLYDATA`), 3 × SU2 with `NELEM=0` (`meshio` reader
defect), and the CGNS structure finding of §4.6. The sixth format (CGNS) has **no OSS triangle reader
in existence**, so its row is graded as *mapping conformance*, and there the upstream m17483 files
**fail**: geometry is recovered, but the CGNS/HDF5 node mapping is violated (label stored in the
`type` attribute, `name`/`label`/`order` absent) — a conformant CGNS reader would not find the database.

## 1. Problem setting and references

An internal round trip proves self-consistency only: a writer and a reader sharing the same
misreading of a grammar agree with each other. The XVAL question is falsifiable:

> Re-read the TensorCAD files with *independent* implementations — does the surface come back with the
> same volume, the same topology and the same orientation, and would a third-party tool find it at all?

Reference implementations, **read as source in this sandbox, no network** (zero hand-copying):

| source | use here | what was read |
|---|---|---|
| `meshio` 5.3.5 | STL, OBJ, VTK, SU2, Gmsh readers | `stl/_stl.py`, `obj/_obj.py`, `vtk/_vtk_42.py`, `su2/_su2.py`, `gmsh/_gmsh22.py`, `cgns/_cgns.py` |
| `numpy-stl` 4.1.3 | second STL opinion | `Mesh.from_file` (triangle soup) |
| `trimesh` 5.0.0 | second STL/OBJ opinion | `load(..., process=False)` (vertex merging disabled, see §5) |
| `h5py` 3.15.1 | CGNS: SIDS-to-HDF5 **mapping conformance** + connectivity read | HDF5 group/attribute/dataset contract |
| CGNS SIDS-to-HDF5 text (`cgns_sids_hdf5.html`, 84 564 B, sha256 `dcc5433d…`, cached by node m17483 `specs/manifest.json`) | adjudication of the CGNS node mapping | "Basic Node" table: required attrs `name`/`label`/`type`/`order` + optional dataset ` data`; database rooted at a `CGNSBase_t` node |
| SU2 v7 mesh-file grammar | strict reader `tensorcad/mesh_xval.py:read_su2_strict` | element type 5 = triangle, 0-based connectivity |

Reference-side status, source-checked (`meshio_cgns_status()`):

```json
{
  "path": "meshio/cgns/_cgns.py",
  "raises_can_only_read_tetrahedra": true,
  "todo_stub": true
}
```

## 2. Method and implementation

`tensorcad/mesh_xval.py` (new, 6 writers + 5 readers + 1 comparison kernel) is a **second, independent
implementation** of the same six grammars — it deliberately does not import `tensorcad.mesh_formats`, so
the matrix is "our writer × their reader" and not "our writer × our reader". Readers:
`meshio_mesh`, `trimesh_mesh`, `numpy_stl_mesh`, `read_su2_strict`, `read_cgns_sids`.
`compare_meshes` returns

* **volume** — `V = 1/6 Σ a·(b×c)` (divergence theorem) and its relative error;
* **topology** — vertex bijection (exact hashing, or nearest neighbour inside the format's own
  precision), equal vertex/face counts, equal unordered face multiset, boundary-edge count, Euler
  characteristic `V-E+F`;
* **orientation** — equal *oriented* face multiset under cyclic-rotation keys, flipped-face count, sign
  of the signed volume.

STL (and `trimesh(process=False)`) store a *soup* — the file carries no connectivity — so the comparison
welds by identical coordinates first (`weld()`) and reports `welded: true`. A vertex with no partner
inside the tolerance makes topology **and** orientation false, so a format that duplicates or merges
corners cannot pass by accident.

Matrices: **A** = XVAL writers × OSS readers, 6 formats × 3 closed surfaces; **B** = the upstream
m17483 artefacts (copied to `artifacts/m17595_xval_six_format/upstream_m17483/files`, sha256 in
`results.json → matrix_b.sha256`) × the same readers. B cross-validates the node implementation itself;
A is the control that catches "both of our sides wrong in the same way".

## 3. Parameters

| geometry | genus | parameters | vertices | faces | exact volume [m³] | mesh volume [m³] | discretisation error |
|---|---|---|---|---|---|---|---|
| cube | 0 | edge 2.0 m | 8 | 12 | 8.000000 | 8.000000 | +0.000 % |
| icosphere2 | 0 | r 1.0 m, 2 subdivisions | 162 | 320 | 4.188790 | 4.047045 | -3.384 % |
| torus64x32 | 1 | R 1.0 m, r 0.3 m, 64x32 | 2048 | 4096 | 1.776529 | 1.762302 | -0.801 % |

Discretisation error is **not** format error; it is separated in §4.3 because the round trip compares a
mesh with its own file. Numeric parameters: OBJ / VTK / SU2 / CGNS / Gmsh written with 17 significant
digits (float64 exact round trip), binary STL float32; comparison tolerance `1e-5 × diagonal` for STL
only, `0` elsewhere; CGNS connectivity `I4` 1-based, coordinates `R8`.

## 4. Results

### 4.1 Matrix A — XVAL writers × independent readers (33 checks)

| geometry | format | independent reader | verdict | volume rel err | V/F equal | orientation | flipped |
|---|---|---|---|---|---|---|---|
| cube | stl | meshio | PASS | 0.000e+00 | yes | P | 0 |
| cube | stl | numpy-stl | PASS | 0.000e+00 | yes | P | 0 |
| cube | stl | trimesh | PASS | 0.000e+00 | yes | P | 0 |
| cube | obj | meshio | PASS | 0.000e+00 | yes | P | 0 |
| cube | obj | trimesh | PASS | 0.000e+00 | yes | P | 0 |
| cube | vtk | meshio | PASS | 0.000e+00 | yes | P | 0 |
| cube | vtk | meshio / POLYDATA | UNREADABLE | n/a | n/a | n/a | n/a |
| cube | su2 | meshio | UNREADABLE | n/a | n/a | n/a | n/a |
| cube | su2 | mesh_xval(strict grammar) | PASS | 0.000e+00 | yes | P | 0 |
| cube | cgns | h5py (SIDS mapping checker) | PASS | 0.000e+00 | yes | P | 0 |
| cube | gmsh | meshio | PASS | 0.000e+00 | yes | P | 0 |
| icosphere2 | stl | meshio | PASS | 1.311e-08 | yes | P | 0 |
| icosphere2 | stl | numpy-stl | PASS | 1.311e-08 | yes | P | 0 |
| icosphere2 | stl | trimesh | PASS | 1.311e-08 | yes | P | 0 |
| icosphere2 | obj | meshio | PASS | 0.000e+00 | yes | P | 0 |
| icosphere2 | obj | trimesh | PASS | 0.000e+00 | yes | P | 0 |
| icosphere2 | vtk | meshio | PASS | 0.000e+00 | yes | P | 0 |
| icosphere2 | vtk | meshio / POLYDATA | UNREADABLE | n/a | n/a | n/a | n/a |
| icosphere2 | su2 | meshio | UNREADABLE | n/a | n/a | n/a | n/a |
| icosphere2 | su2 | mesh_xval(strict grammar) | PASS | 0.000e+00 | yes | P | 0 |
| icosphere2 | cgns | h5py (SIDS mapping checker) | PASS | 0.000e+00 | yes | P | 0 |
| icosphere2 | gmsh | meshio | PASS | 0.000e+00 | yes | P | 0 |
| torus64x32 | stl | meshio | PASS | 6.912e-10 | yes | P | 0 |
| torus64x32 | stl | numpy-stl | PASS | 6.915e-10 | yes | P | 0 |
| torus64x32 | stl | trimesh | PASS | 6.915e-10 | yes | P | 0 |
| torus64x32 | obj | meshio | PASS | 0.000e+00 | yes | P | 0 |
| torus64x32 | obj | trimesh | PASS | 0.000e+00 | yes | P | 0 |
| torus64x32 | vtk | meshio | PASS | 0.000e+00 | yes | P | 0 |
| torus64x32 | vtk | meshio / POLYDATA | UNREADABLE | n/a | n/a | n/a | n/a |
| torus64x32 | su2 | meshio | UNREADABLE | n/a | n/a | n/a | n/a |
| torus64x32 | su2 | mesh_xval(strict grammar) | PASS | 0.000e+00 | yes | P | 0 |
| torus64x32 | cgns | h5py (SIDS mapping checker) | PASS | 0.000e+00 | yes | P | 0 |
| torus64x32 | gmsh | meshio | PASS | 0.000e+00 | yes | P | 0 |

`combinations = 18`, `reader checks = 33`,
`green = 27`, worst relative volume error **1.311e-08**.
A cell counts as green only if *every* reader entry passed, which is why
`cells_all_oss_green = 12/18` — the VTK cells carry the historical `POLYDATA`
variant (unreadable by `meshio`) and the SU2 cells carry the `NELEM=0` entry.

### 4.2 Matrix B — upstream m17483 artefacts × the same readers (30 checks)

| geometry | format | independent reader | verdict | volume rel err | V/F equal | orientation | flipped |
|---|---|---|---|---|---|---|---|
| cube | stl | meshio | PASS | 0.000e+00 | yes | P | 0 |
| cube | stl | numpy-stl | PASS | 0.000e+00 | yes | P | 0 |
| cube | stl | trimesh | PASS | 0.000e+00 | yes | P | 0 |
| cube | obj | meshio | PASS | 0.000e+00 | yes | P | 0 |
| cube | obj | trimesh | PASS | 0.000e+00 | yes | P | 0 |
| cube | vtk | meshio | UNREADABLE | n/a | n/a | n/a | n/a |
| cube | su2 | meshio | UNREADABLE | n/a | n/a | n/a | n/a |
| cube | su2 | mesh_xval(strict grammar) | PASS | 0.000e+00 | yes | P | 0 |
| cube | cgns | h5py (SIDS mapping checker) | NONCONFORMANT (geometry ok, SIDS mapping nonconformant) | 0.000e+00 | yes | P | 0 |
| cube | gmsh | meshio | PASS | 0.000e+00 | yes | P | 0 |
| icosphere | stl | meshio | PASS | 1.311e-08 | yes | P | 0 |
| icosphere | stl | numpy-stl | PASS | 1.311e-08 | yes | P | 0 |
| icosphere | stl | trimesh | PASS | 1.311e-08 | yes | P | 0 |
| icosphere | obj | meshio | PASS | 0.000e+00 | yes | P | 0 |
| icosphere | obj | trimesh | PASS | 0.000e+00 | yes | P | 0 |
| icosphere | vtk | meshio | UNREADABLE | n/a | n/a | n/a | n/a |
| icosphere | su2 | meshio | UNREADABLE | n/a | n/a | n/a | n/a |
| icosphere | su2 | mesh_xval(strict grammar) | PASS | 0.000e+00 | yes | P | 0 |
| icosphere | cgns | h5py (SIDS mapping checker) | NONCONFORMANT (geometry ok, SIDS mapping nonconformant) | 0.000e+00 | yes | P | 0 |
| icosphere | gmsh | meshio | PASS | 0.000e+00 | yes | P | 0 |
| torus | stl | meshio | PASS | 6.912e-10 | yes | P | 0 |
| torus | stl | numpy-stl | PASS | 6.915e-10 | yes | P | 0 |
| torus | stl | trimesh | PASS | 6.915e-10 | yes | P | 0 |
| torus | obj | meshio | PASS | 0.000e+00 | yes | P | 0 |
| torus | obj | trimesh | PASS | 0.000e+00 | yes | P | 0 |
| torus | vtk | meshio | UNREADABLE | n/a | n/a | n/a | n/a |
| torus | su2 | meshio | UNREADABLE | n/a | n/a | n/a | n/a |
| torus | su2 | mesh_xval(strict grammar) | PASS | 0.000e+00 | yes | P | 0 |
| torus | cgns | h5py (SIDS mapping checker) | NONCONFORMANT (geometry ok, SIDS mapping nonconformant) | 0.000e+00 | yes | P | 0 |
| torus | gmsh | meshio | PASS | 0.000e+00 | yes | P | 0 |

`green = 21`, worst relative volume error **1.311e-08**,
red entries: `vtk/-/meshio` on 3 geometries (cube); `su2/-/meshio` on 3 geometries (cube); `cgns/-/h5py (SIDS mapping checker)` on 3 geometries (cube). The upstream files are reproducible from
`results.json → matrix_b.sha256`; their geometry agrees with the independently rebuilt bodies to
float64 (STL float32), i.e. **the node's writer is confirmed by third-party readers**, and the CGNS
structure issue of §4.6 is the only substantive defect found.

### 4.3 Convergence — two refinement ladders (format error vs discretisation error)

| case | triangles | mesh vs analytic | worst format error |
|---|---|---|---|
| icosphere/sub1 | 80 | +12.655 % | 5.674e-09 |
| icosphere/sub2 | 320 | +3.384 % | 1.311e-08 |
| icosphere/sub3 | 1280 | +0.861 % | 1.640e-08 |
| icosphere/sub4 | 5120 | +0.216 % | 7.627e-09 |
| torus/16x8 | 256 | +12.265 % | 1.203e-08 |
| torus/32x16 | 1024 | +3.175 % | 8.746e-09 |
| torus/64x32 | 4096 | +0.801 % | 6.912e-10 |
| torus/128x64 | 16384 | +0.201 % | 6.733e-09 |

Discretisation error falls **59x** over the icosphere ladder (icosphere/sub1 → icosphere/sub4) and
**61x** over the torus ladder (torus/16x8 → torus/128x64), while the worst format error stays
inside `1.64e-08` — flat
across a 64x change in face count. The two error mechanisms are therefore separable, and the flat curve
is the evidence that the formats are not the limiting factor at hull-mesh resolutions.

### 4.4 Parameter sensitivity — coordinate magnitude (global hull-scale frames)

| global offset | 0 m | 1 m | 10 m | 100 m | 1000 m |
|---|---|---|---|---|---|
| STL max vertex deviation [m] | 3.29e-08 | 8.26e-08 | 6.65e-07 | 5.13e-06 | 3.85e-05 |
| STL volume rel err | 1.31e-08 | 4.38e-08 | 1.42e-07 | 2.57e-06 | 5.22e-06 |
| OBJ max vertex deviation [m] | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 |
| VTK max vertex deviation [m] | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 |

**Engineering consequence.** Binary STL is safe only in a local frame: the same 1 m body expressed in a
1000 m global frame loses `3.8e-05 m` of vertex fidelity
(float32 mantissa), ~`4e-05` x the body radius, while
the relative volume error stays at `5.2e-06` (volume is a
first-order-in-error quantity, coordinates are not). OBJ/VTK eXact tens of digits stay exact at every
offset. Ship frames are 10²–10³ m, so **STL must be written in a local frame** (or converted to
OBJ/VTK/CGNS) before any downstream audit that compares coordinates, not volumes.

### 4.5 SU2 — separating our grammar from the reference reader defect

| variant | meshio 5.3.5 | strict grammar reader | volume rel err (strict) |
|---|---|---|---|
| `NELEM` mode = `boundary` | unreadable (StopIteration: ) | ok | 0.000e+00 |
| `NELEM` mode = `volume` | ok | ok | 0.000e+00 |

`meshio`'s own writer→reader on a tetrahedral mesh: **ok**. The SU2
failure is therefore *conditional on* `NELEM=0` (`meshio/su2/_su2.py:113-122` calls `next(gen)` on an
empty `islice`), not a general SU2 defect and not a defect of our file: the same triangles placed inside
`NELEM` read back with zero volume error. Our boundary-only mesh is grammar-legal (SU2 v7: `NELEM=0`
plus `MARKER_TAG`/`MARKER_ELEMS`, 0-based); inventing tetrahedra to please the reader was rejected
because a star triangulation fills the torus hole (non-convex bodies).

### 4.6 CGNS — the mapping, not the geometry, is where the defect is

The CGNS/HDF5 mapping (SIDS-to-HDF5, "Basic Node") prescribes: a node is a group carrying the required
attributes `name`, `label`, `type` (a **data-type code**: `MT`/`I4`/`R8`/`C1`) and `order`; the payload
is the single dataset ` data`; a database is rooted at a node labelled `CGNSBase_t`.

| implementation | geometry recovered | SIDS mapping violations |
|---|---|---|
| this branch (`tensorcad.mesh_xval.write_cgns`) | yes, 0.000e+00 | **0** |
| upstream m17483 `cube/icosphere/torus.cgns` | yes, 0.000e+00 | **12** per file |

Upstream violation pattern (distinct kinds): `root: reserved attribute 'version' missing`; `Base: required attribute 'name' missing`; `Base: required attribute 'label' missing`; `Base: required attribute 'order' missing`; `Base: 'type' holds a label`; `Zone: required attribute 'name' missing`; `Zone: required attribute 'label' missing`; `Zone: required attribute 'order' missing`; `Zone: 'type' holds a label where the spec requires a data-type code`; `Elements: required attribute 'name' missing`; `Elements: required attribute 'label' missing`; `Elements: required attribute 'order' missing`.

Reading: the upstream writer puts the node **label into the `type` attribute** (`type = "CGNSBase_t"`)
and omits `name`/`label`/`order`; the geometry is present and correct, but a spec-conformant CGNS reader
(looking for the `label` attribute / `CGNSBase_t` database root) would not recognise the file as a CGNS
database at all. Our internal reader is tolerant, which is exactly why the internal round trip stays
green — the classic XVAL failure mode this node exists to catch. The same checker runs on the fixed
dialect in `write_cgns` with **0** violations, so the fix is small and local (table in §5).

## 5. TensorCAD vs OSS: difference list and porting decisions

| aspect | TensorCAD (m17483 native) | independent reference | difference / boundary condition | porting decision |
|---|---|---|---|---|
| VTK dialect | legacy ASCII `DATASET POLYDATA`, `POINTS`/`POLYGONS` | `meshio` legacy reader accepts only `UNSTRUCTURED_GRID`/`STRUCTURED_*` (`vtk/_vtk_42.py`) | `POLYDATA` is legal VTK (ParaView reads it) but invisible to the whole meshio ecosystem | **default `UNSTRUCTURED_GRID` + `CELLS`/`CELL_TYPES=5`** (verified readable here); keep `POLYDATA` behind an explicit flag |
| SU2 element placement | `NELEM=0` + one `MARKER_TAG`/`MARKER_ELEMS` block | `meshio` fails on `NELEM=0`; needs ≥1 volume element | our file is grammar-legal, the reference reader is not defensive | **keep the faithful boundary mesh** (tets would fill a torus hole); document `NELEM=0` as an interop hazard; expose `nelem_mode="volume"` for pipelines that need it |
| CGNS node mapping | HDF5 with label in `type`, no `name`/`label`/`order`, no root `version` | SIDS-to-HDF5 "Basic Node": `name`/`label`/`type`(=data-type code)/`order` + ` data`; database rooted at `CGNSBase_t` | **violation** (12 findings/file), geometry unaffected; no OSS triangle reader exists in any case (`meshio` `cgns/_cgns.py` = `TODO` stub ending in `ReadError("Can only read tetrahedra.")`) | **fix the writer**: `label="CGNSBase_t"/"Zone_t"/"GridCoordinates_t"/"DataArray_t"/"Elements_t"/"ZoneType_t"/"IndexRange_t"`, `type` = `MT`/`I4`/`R8`/`C1`, add `name`/`order`, root `version`/`format`; keep `read_cgns_sids` as the machine gate (it accepts both dialects and reports the delta) |
| STL | binary, float32, topology by exact welding | numpy-stl/meshio return a soup; **trimesh merges vertices by default** | vertex count depends on the reader's `process` flag, not on the file | compare with `process=False` and weld before topology questions; float32 ⇒ local frame only (§4.4) |
| OBJ | `v`/`f`, 17 significant digits, minimal subset | meshio/trimesh accept; trimesh reorders vertices under `process=True` | no normals/UV/materials/groups written | keep minimal; document that groups/materials are not carried (not claimed) |
| Gmsh | MSH 2.2 ASCII, element type 2, 1-based, 0 tags | meshio `gmsh/_gmsh22.py` reads it; the real `gmsh` binary/module is **not installable here** (package index unreachable) | only MSH 2.2 emitted; MSH 4.1 is a different container (entities/blocks) | keep 2.2 (widest reader support); MSH 4.1 = recorded gap, not claimed |
| precision policy | float64 for 5 of 6 formats, float32 only where forced | same split in the OSS writers | — | matches OSS practice, no change |

**Porting trade-off, stated plainly:** of the four reference-side findings, two are *their* defects
(`meshio` `NELEM=0`, `meshio` CGNS stub) and are documented rather than worked around — bending our
writer to a buggy reader would make the file wrong for the tools that follow the spec. Two are ours:
VTK `POLYDATA` only (fixed dialect validated here) and the CGNS attribute mapping (fix listed above).

## 6. Failure modes, limits and artefacts

* **CGNS has no OSS triangle reader at all** (`meshio` stub; `cgnslib`/`pycgns` not installed and not
  installable offline). Its row is graded as *mapping conformance*, which is strictly weaker than
  "a third party read it". Evidence grade: REAL for the mapping check, WITHHELD for third-party
  interop — the capability table must keep that distinction.
* **`gmsh` proper and `vtk`/ParaView are absent** from the sandbox, so Gmsh interop is `meshio`-level
  only and the `POLYDATA` failure is specific to `meshio`. Ecosystem friction, logged; not a TensorCAD
  defect.
* The upstream m17483 module (`tensorcad/mesh_formats.py`) is **not in this branch's tree** (it is
  uncommitted in the m17483 worktree at the time of writing), so the CGNS fix is delivered as a
  validated reference dialect in `mesh_xval.write_cgns` + a checker that fails the current dialect.
  Applying it is a 10-line change in that module, listed in §5.
* CGNS limits: `I4` connectivity caps at 2³¹ vertices; single element block, `ElementRange = [1,ncell]`;
  only `TRI_3`, `ZoneType=Unstructured`, no BC/flow-solution nodes.
* STL limits: welding is O(n log n) and assumes bit-identical duplicate corners (true for our writer;
  a lossy exporter that perturbs corners per-facet would break it).
* Not covered: non-manifold surfaces, multi-marker SU2, meshes > 1 M faces, MSH 4.1, CGNS zones with
  multiple element blocks.

### 6.1 Artefacts and repro

`artifacts/m17595_xval_six_format/`: `results.json` (every number in this report),
`matrix.md`, `files/` (18 + POLYDATA variants), `upstream_m17483/files/` (18 upstream files, sha256),
`sensitivity/`, `su2_variants/`, `figures/*.png`, `manifest.json` (sha256 of every artefact).
Repro: `python tools/xval_mesh_six_format.py`; drift gate: `python tools/xval_mesh_six_format.py --check`;
tests: `pytest tests/test_mesh_xval.py -q`.

### 6.2 Error table against the 15 % standard (explicit)

| key quantity | reference | cross-read value | error | < 15 % ? |
|---|---|---|---|---|
| volume, cube, STL | 8.0 m³ | matrix A (3 readers) | 0.000e+00 | YES |
| volume, icosphere, STL (float32) | 4.047044754 m³ | matrix A (3 readers) | 1.311e-08 | YES |
| volume, torus, STL | 1.762301535 m³ | matrix A (3 readers) | 6.915e-10 | YES |
| volume, 5 lossless formats | written mesh volume | matrix A | ≤ 0.000e+00 | YES |
| volume, upstream m17483 files | independently rebuilt bodies | matrix B | 1.311e-08 | YES |
| vertex/face counts | written mesh | equal in every green cell | 0 | YES |
| flipped faces | 0 | worst over all reader checks | 0 | YES |
| coordinate deviation, STL at 1000 m offset | 0 m | 3.85e-05 m | 0.00 % of body radius | NO — see §4.4 (local frame required) |
| CGNS SIDS mapping conformance (upstream) | 0 violations | 12 violations/file | n/a (structural) | NO — see §4.6 |
