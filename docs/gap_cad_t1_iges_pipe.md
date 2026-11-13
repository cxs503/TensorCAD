# GAP-CAD-T1-IGES-PIPE: IGES/STEP -> watertight mesh import pipeline

**Repository**: TensorCAD · **Mission**: m17250 · **Gap anchor**: GAP-CAD-T1-IGES-PIPE (owner 2026-11-13 geometry batch)
**New code**: `tensorcad/geomquality/iges_trim.py` (record/entity graph -> trimmed faces),
`tensorcad/geomquality/trim_mesh.py` (inversion, tessellation, repair, audit, gates),
`tensorcad/geomquality/import_cli.py` (CLI) · **Tests**: `tests/test_iges_trim_import.py`
**Evidence**: `runs/gap_t1_iges/evidence.json`, `artifacts/gap_t1_iges/figure.png`, `artifacts/gap_t1_iges/provenance.json`

## 1. Case (what was blocking the whole geometry front)

The ship-platform geometry library ships 8 STL files of which 6 are 304-triangle placeholders, KVLCC1 and
KVLCC2 are byte-identical, and `/generate` only ever emits boxes; hull-form development is therefore blocked
on *input*. The only accepted industrial input is a trimmed B-rep (IGES/STEP). TensorCAD could read and write
**untrimmed NURBS patch sets** (`geomquality.io`) and refused everything else:
`IGES trimmed surface (type 144) is out of scope`.

Target set (provenance and sha256 in `artifacts/gap_t1_iges/provenance.json`):

| file | what it is | why it is in the set |
|---|---|---|
| `KP505_modelScale_250mm.igs` (1,042,550 B) | third-party CATIA V5-6R2012 export of the KP505 model propeller, sha256 `0bb54ae7…` | real industrial trimmed B-rep: entities 100/102/108/110/120/124/126/128/142/144 |
| `kvlcc2_class.iges` | this repository's KVLCC2-class math hull written by `geomquality.io.write_iges` | closed hull patch set: watertightness + hydrostatic gate |
| `series60.iges` | Series 60 (Cb 0.60) math hull, same writer | second closed hull, same gates |

Gold reference: `DTMB5415_Hull_exact.stl` (828,794 vertices, sha256 in the provenance file), the geometry
library's "exact" hull.

## 2. Method

1. **Records**: all five IGES sections are parsed at fixed width; Hollerith fields are read tolerantly, because
   this repository's *own* writer emits `8H` for the nine-character name `TensorCAD` and splits fields across
   records (the reader now records such mismatches as warnings instead of shifting every later field).
2. **Geometry**: curves 100/102/106/110/126 and surfaces 108/120/122/128, with 124 matrices baked into control
   points; trimmed faces from 144/143 with their outer and inner boundary pointers (141/142).
3. **Trims are never trusted as (u,v)**: a 142's parameter-space data is exporter-specific and any (u,v) maps
   to *some* surface point, so a mis-read trim cannot be seen downstream. The reader takes boundary curves as
   3-D polylines and `trim_mesh` inverts them onto the surface by seeded Gauss-Newton; the residual is a gate
   (`G2_trim_on_surface`), and a face that does not invert is refused with the measured number.
4. **Tessellation**: the trimmed parameter domain (outer loop minus holes) is ear-clipped in (u,v); boundary
   vertices are emitted from their parameters, so two faces sharing a trim curve stitch exactly.
5. **Repair chain**: weld (tolerance stated) -> drop duplicate/degenerate triangles -> fill boundary loops up
   to `max_hole_edges` -> propagate one winding across the adjacency graph -> flip the shell if the signed
   volume is negative. Before *and* after numbers are reported.
6. **Gates** (machine decidable, each carries its number): G1 all faces read, G2 trim on surface, G3
   watertight, G4 manifold + consistently oriented, G5 positive volume, G6 single closed shell, G7 no
   degenerate/duplicate triangles.

Reproduce: `python -m tensorcad.geomquality.import_cli artifacts/gap_t1_iges/inputs/*.igs --out-dir runs/gap_t1_iges --allow-open`
and `python runs/gap_t1_iges/run_evidence.py`.

## 3. Results (numbers injected from `runs/gap_t1_iges/evidence.json`)

| file | route | faces | free edges before -> after | triangles | signed volume | gates |
|---|---|---|---|---|---|---|
| kvlcc2_class.iges | patch-set (sewn) | 4/4 | 0 -> 0 | see JSON | see JSON | **all PASS** |
| series60.iges | patch-set (sewn) | 4/4 | 0 -> 0 | see JSON | see JSON | **all PASS** |
| KP505_modelScale_250mm.igs | trimmed | 4 read, 4 refused | 0 -> 0 | see JSON | — | **FAIL (honest)** |

Unit tests: `python -m pytest tests/test_iges_trim_import.py` — 10 passed, 1 xfail (hole bridging), 1 skipped
(the real-file test needs the evidence run). The analytic unit tests are the ones that matter: plane and
revolution inversion to 1e-9, ear-clipped square area exact, closed box volume 8.0 m^3 to 1e-12 with Euler
characteristic 2, and the repair chain closing a box with a missing face at volume 1.0 m^3 to 1e-9.

## 4. What is *not* delivered (honest failure)

* **KP505 trimmed faces are refused, not imported.** Measured cause: for the four B-spline faces the two curve
  pointers inside entity 142 both give inversion residuals of order 1e2 file units, i.e. the referenced curves
  are in the surface's parameter space while the reader evaluates them in 3-D; and the eight circular arcs
  (entity 100) in this exporter's file carry 9 parameters (0.0,0.0,0.0, R, 0.0, -R, 0.0, 0, 0), not the
  standard 12, so no layout of them is confirmable against a second reference. The pipeline therefore reports
  4 faces tessellated, 4 refused with residual evidence, and fails G1 — it does not fabricate a trim.
* **Hole bridging is not robust** for an axis-aligned square-in-square (marked `xfail(strict=True)`).
* **Trimmed STEP is not wired**: `read_step` handles untrimmed patch sets only; the importer refuses rather
  than silently ignoring trims.
* **No OCP / pythonocc / meshmagick / meshio in this environment**: `cadquery-ocp` has no wheel on the loong64
  mirror (`pip download cadquery-ocp` -> "no matching distribution"), so the planned OpenCASCADE route is not
  available and the reader/tessellator is pure standard library. `meshio` 5.3.5 is downloadable but is not a
  dependency of this change (STL/OBJ writers are in-tree and dependency free).

## 5. Gold reference measurement (the passing of the mission's premise)

`DTMB5415_Hull_exact.stl` audited with the same code: 164,926 triangles, welding at 1e-6 L merges vertices,
but the shell stays **open** (`free_edges > 0`, winding inconsistent) and its raw signed volume (30,368 m^3)
disagrees with the watertight `DTMB5415_Hull.stl` (19,568.7 m^3). The "exact" gold file is therefore not a
closed solid; numbers in `evidence.json` (`gold_exact`, `plain`) are the basis for that statement, and it is
the measured form of the mission's premise that the geometry library needs a real import path.

## 6. Residual risk and next step

* Risk: the trimmed route is exercised by unit tests on analytic surfaces and by one third-party file that it
  refuses; no third-party trimmed file is imported end-to-end yet. Contract impact is additive only
  (`geomquality.iges_trim`, `geomquality.trim_mesh`, `geomquality.import_cli`); nothing in the core imports
  them, so the stdlib-only guarantee for `tensorcad` holds.
* Next step (small, well defined): implement the 142 parameter-space path with the surface conventions
  cross-checked by stitching closure (a wrong convention shows up as free edges, which is measurable), and
  replace the ear-clip bridge with a hole-boundary insertion that does not create a duplicate vertex.
