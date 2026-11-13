# Native six-format mesh matrix (m17907)

Anchor: `XVAL-CAD-M17483-NATIVE-SIX-FORMAT-MATRIX` (acceptance: relative error < 15 %, owner standard).
Module under test: `tensorcad/mesh_formats.py`. Artefacts: `artifacts/m17907_native_six_format/`.
Repro: `python3 tools/native_six_format_matrix.py` (drift gate: `--check`); tests:
`python3 -m pytest tests/test_mesh_formats.py -q` (**38 passed**).

---

## 1. Case: what a "six-format matrix" has to prove

A surface mesh that a hull/solver pipeline hands downstream must survive six different
on-disk grammars without losing geometry or topology, and the *file* (not our reader) must be
readable by a third party. The three bodies are chosen so the reference number is analytic and
the topology is non-trivial:

| body | vertices | triangles | reference volume | why this body |
|---|---|---|---|---|
| `cube(2.0)` | 8 | 12 | **8.0 m³ exactly** | exact float64 oracle; catches any coordinate scaling/transposition |
| `icosphere(1.0, 2)` | 162 | 320 | `4/3·π = 4.18879` (inscribed, O(h²) low) | curvature: float32 damage becomes visible |
| `torus(1.0, 0.3, 64, 32)` | 2048 | 4096 | `2π²Rr² = 1.77653` (O(h²) low) | **genus 1** — Euler χ = 0, so a topology defect cannot hide behind a sphere |

The bodies themselves are pinned by a convergence test: the icosphere volume error falls
3.39e-2 → 8.62e-3 → 2.17e-3 over subdivisions 2→3→4, i.e. ratios **3.93 / 3.98** (second order,
`test_volume_matches_analytic_reference_and_converges_second_order`). If the generator were wrong,
the comparison reference would be wrong too, so it is tested first.

## 2. Method

* **Writer under test** — `tensorcad/mesh_formats.py` (new in this mission): six writers, six own
  readers, a CGNS SIDS-to-HDF5 conformance checker, and one comparison contract
  (`compare_meshes` → volume, vertex bijection, face multiset, *oriented* face multiset, Euler
  characteristic, boundary edges, `pass`).
* **Oracle 1 — own reader**: write → read back → compare. Exercises the grammar, not a library.
* **Oracle 2 — third party**: `meshio 5.3.5`, `numpy-stl` (`import stl`), `trimesh 5.0.0` read the
  written file; a reader that cannot parse it is recorded as **unreadable**, never skipped.
* **Oracle 3 — second implementation**: `tensorcad.mesh_xval` writes the same six grammars with an
  independently written writer set (its docstring states it does not import this module). Writer ×
  writer agreement is a genuine cross-check of both writers and both readers.
* **Gate**: `impl_threshold = 1e-3` relative volume error per cell, topology must pass, CGNS must
  have `sids_violations() == []`. Owner line 0.15 is reported alongside.

Error is measured **against the mesh we wrote** (format fidelity) and **separately against the
analytic value** (discretisation, which no file format can fix). Conflating the two is how
"format" errors get invented or hidden.

## 3. Results

### 3.1 Round trip, own reader (18 cells)

| body | format | written V [m³] | read-back V [m³] | rel err | topo | bytes |
|---|---|---|---|---|---|---|
| cube | stl | 8 | 8 | 0.00e+00 | PASS | 684 |
| cube | obj | 8 | 8 | 0.00e+00 | PASS | 200 |
| cube | vtk | 8 | 8 | 0.00e+00 | PASS | 321 |
| cube | su2 | 8 | 8 | 0.00e+00 | PASS | 245 |
| cube | cgns | 8 | 8 | 0.00e+00 | PASS | 20 040 |
| cube | gmsh | 8 | 8 | 0.00e+00 | PASS | 303 |
| icosphere2 | stl | 4.047005367355 | 4.047005348804 | **4.584e-09** | PASS | 16 084 |
| icosphere2 | obj/vtk/su2/cgns/gmsh | 4.047005367355 | 4.047005367355 | 0.00e+00 | PASS | — |
| torus64x32 | stl | 1.762301534515 | 1.762301535733 | **6.915e-10** | PASS | 204 884 |
| torus64x32 | obj/vtk/su2/cgns/gmsh | 1.762301534515 | 1.762301534515 | 0.00e+00 | PASS | — |

`results.json` holds every cell (18 rows), the byte sizes and the sha256 of each file; the numbers in
this table are read from it, not typed in.

Summary (`results.json.summary`):

```
cells                              18
own_reader_topology_pass           18 / 18
worst_own_reader_rel_err           4.5837e-09   (icosphere2, stl — float32)
own_reader_all_within_impl_threshold  true  (1e-3)
own_reader_all_within_owner_threshold true  (0.15)
third_party_readable_cells         15 / 18      (3 × cgns unreadable: meshio stub)
third_party_worst_rel_err          4.5837e-09   (stl, vs written mesh volume)
cross_writer_agreement             18 / 18      (vs tensorcad.mesh_xval writers)
cgns_native_violations             0
cgns_legacy_violations             1
```

### 3.2 Third-party interop (15 readable / 3 unreadable)

| format | reader | result |
|---|---|---|
| STL | `numpy-stl` | ok, volume error ≤ 4.6e-09 (float32) |
| OBJ | `meshio` | ok, exact |
| VTK (`UNSTRUCTURED_GRID`, the default) | `meshio` | ok, exact |
| SU2 (`nelem_mode="volume"`, the default) | `meshio` | ok, exact |
| Gmsh MSH 2.2 | `meshio` | ok, exact |
| CGNS | `meshio` | **unreadable** — `meshio/cgns/_cgns.py` is a stub that expects `Base/Zone1` and tetrahedra; it ends in `ReadError`/`SystemExit(1)` |

The CGNS cell is honest-FAIL on *interop*, not on *conformance*: no OSS triangle reader for CGNS
exists in this environment (`cgnslib`/`pycgns` not installable offline), so the capability is
graded as **mapping conformance** (REAL) and third-party consumption is **WITHHELD** for CGNS. The
distinction is carried in the tables above rather than blurred into a headline.

### 3.3 CGNS SIDS-to-HDF5 conformance — the defect this mission removes

The previous reference dialect (`tensorcad.mesh_xval.write_cgns`) satisfied *name/label/type/order*
and the reserved root attributes, but had **no `CGNSLibraryVersion_t` node**, so "SIDS conformant"
was an overstatement (board audit m17634, finding F4). The native writer adds the node and the
checker now requires it:

| file | violations |
|---|---|
| `files/sids_probe_native.cgns` (native writer) | `[]` |
| `files/sids_probe_legacy.cgns` (`mesh_xval` dialect) | `["root: CGNSLibraryVersion_t node missing"]` |

The negative control is a test, not a claim:
`test_checker_rejects_a_dialect_without_library_version` builds a file that differs *only* by the
missing node, asserts exactly that one violation, and asserts `read_cgns` refuses to launder it.
`test_sids_checker_flags_a_wrong_attribute_mapping` deletes one `order` attribute and expects a
violation naming it — so the checker is not vacuous.

### 3.4 Writer × writer (18/18 agree)

For all 18 cells the independently written `mesh_xval` file is read back and compared: relative
volume error ≤ 4.584e-09 and topology PASS everywhere. The three CGNS cells agree geometrically
while being structurally non-conformant (the legacy dialect's violation list is recorded in the
same row), which is precisely the separation the audit asked for.

## 4. Design decisions and their measured basis

| decision | alternative | measured consequence of the alternative |
|---|---|---|
| VTK default `UNSTRUCTURED_GRID` + `CELLS`/`CELL_TYPES=5` | legacy `POLYDATA` | `meshio` legacy reader rejects `POLYDATA` (`vtk/_vtk_42.py`) ⇒ the file is invisible to the Python ecosystem. `POLYDATA` kept behind `dataset="POLYDATA"` and covered by a test |
| SU2 default `nelem_mode="volume"` | faithful boundary subset (`NELEM=0`) | `meshio/su2/_su2.py:113` raises `StopIteration` on `NELEM=0`. Boundary mode still emitted, legal, and read by our own reader (`test_su2_boundary_mode_is_legal_and_parsed_here`) |
| CGNS with `CGNSLibraryVersion_t` | version-less dialect | 1 conformance violation per file (measured above) |
| STL float32 + weld before topology | compare raw facet soup | raw soup has 3·n_faces vertices; topology questions are meaningless until welded (`test_stl_is_a_facet_soup_and_needs_welding`) |
| Gmsh MSH 2.2 only | MSH 4.1 | 4.1 is a different container (entities/blocks); refused with an explicit error instead of a half-parse |

## 5. Limits, honest

* **CGNS third-party consumption is untested** (no `cgnslib`/`pycgns` offline). Grade: mapping
  conformance REAL, external consumption WITHHELD.
* **CGNS subset**: one unstructured zone, one `TRI_3` block, `I4` connectivity (vertex count < 2³¹),
  no BC/flow-solution nodes, no multi-zone.
* **STL local frame only**: float32 ⇒ at a 1000 m coordinate offset the representation error is
  ~4e-5 m — the same 15 %-table row that m17595 marked NO; not re-measured here (out of scope for
  this anchor, and unchanged by this module).
* **`_bijection(tol>0)` is O(n²)** in the worst case; fine for the ≤ 4096-triangle bodies here and
  for typical hull panels, not for million-face meshes.
* **OBJ minimal subset**: no normals/UV/materials/groups are written (not claimed).
* **Not covered**: non-manifold surfaces, multi-marker SU2, MSH 4.1, meshes > 1 M faces.
* The `tensorcad/__init__.py` re-export makes the module discoverable
  (`tensorcad.mesh_formats`) without shadowing the top-level `TriMesh`-free namespace.

## 6. Error judgement against the standards (explicit)

| key quantity | reference | got | error | < 3 % ? | < 15 % ? |
|---|---|---|---|---|---|
| volume, cube, all six formats | 8.0 m³ (analytic) | 8.0 | 0.00e+00 | **yes** | yes |
| volume, icosphere2, STL | 4.047005367355 (written mesh) | 4.047005348804 | 4.584e-09 (4.6e-7 %) | **yes** | yes |
| volume, icosphere2, 5 float64 formats | 4.047005367355 | 4.047005367355 | 0.00e+00 | **yes** | yes |
| volume, torus, STL | 1.762301534515 | 1.762301535733 | 6.915e-10 | **yes** | yes |
| volume, torus, 5 float64 formats | 1.762301534515 | 1.762301534515 | 0.00e+00 | **yes** | yes |
| third party (numpy-stl/meshio) vs written mesh | written V | written V | ≤ 4.58e-09 | **yes** | yes |
| cross-writer (mesh_xval) vs written mesh | written V | written V | ≤ 4.58e-09 | **yes** | yes |
| CGNS SIDS conformance, native | 0 violations | 0 | — | **yes** | yes |
| CGNS SIDS conformance, legacy dialect | 0 violations | 1 violation | structural | **no** (documented, fixed by this module) | n/a |
| CGNS third-party read | readable | unreadable (meshio stub) | — | **unknown / not run** (no cgnslib) | n/a |
| icosphere volume vs analytic sphere | 4.18879 | 4.04701 (sub2) | 3.38 % | no by construction | yes — this is O(h²) discretisation of the body, not a format error |

## 7. Artefacts

```
artifacts/m17907_native_six_format/
  results.json         every number above (18 rows + cross-writer + sids + env + sha256 index)
  matrix.md            generated tables
  files/               18 native files + 18 mesh_xval files + 2 SIDS probes
  figures/fig_native_six_format.png
                       (A) log-scale round-trip error per format/body with the 1e-3 gate and the
                           15 % owner line; (B) torus body, χ = 0
```