# TensorCAD meshio format matrix — extended to a real ship hull

**Mission**: m18268 — anchor `GAP-CAD-B2-MESH-FORMAT-MATRIX`
**Date**: 2026-11-13 · **Repo**: TensorCAD @ `agent/m18268`
**Run evidence**: `artifacts/m18268_meshio_matrix_hull/matrix18.json` (this document is
generated from it by `scripts/gen_m18268_doc.py`; no number below is hand-copied)

---

## 1. Case introduction

### Physical background
Hull geometry has to leave TensorCAD and reach CFD/structural solvers. meshio is
the unified reader/writer used for that, and the format matrix asks a single
question per (format, geometry) pair: *does the geometry survive the trip out and
back with its volume, topology and orientation intact?*

The matrix already existed for 18 combinations (6 formats × 3 geometries), but all
three geometries were **convex toy bodies** (cube, unit sphere, capsule), while the
module docstring claimed ship hulls work "in practice". A hull is the case that
matters and the case that was not measured, so it is now the fourth geometry.

### Geometry added: Wigley benchmark hull

| Item | Value |
|---|---|
| Source | `tensorcad.wigley.generate_wigley()` (the package's own validated generator, not a hand-rolled fixture) |
| Hull form | Wigley parabolic: half breadth `(B/2)(1-&#124;2x/L-1&#124;^n)(1-(w/T)^m)`, n = m = 2 |
| Main dimensions | L = 10.0 m, B = 1.0 m, T = 0.625 m |
| Discretisation | nx = 41 stations × nz = 13 levels, closed by a deck cap |
| Mesh | 1001 vertices / 1998 triangles, watertight = True, edge-manifold = True |
| Analytic volume | `4LBT/9` = 2.777778 m³ (`tensorcad.wigley.analytic_volume`) |
| Discrete volume | 2.771222 m³ → **0.24 %** low (grid under-resolves the ends) |

### Formats

| Format | Topology it can carry | Verification status |
|---|---|---|
| STL | triangles | REAL (round trip) |
| OBJ | triangles | REAL (round trip) |
| VTK | triangles | REAL (round trip) |
| Gmsh | triangles | REAL (round trip) |
| SU2 | tets only | REAL (round trip, after tetrahedralisation) |
| CGNS | tets only | REAL (round trip) · third-party CGNS reader conformance WITHHELD (no `cgnslib` on host) |

## 2. What was done

1. **Measured the missing case first, honestly.** Feeding the hull to the existing
   `tetrahedralise_surface` raised

   ```
   TetrahedralisationError: no interior point sees the whole surface: best candidate
   (5.0, 0.0, -0.3125) still gives a worst-case tetrahedron volume 0.000e+00
   ```

   Diagnosis (not a guess — the offending facets were printed): where the Wigley
   section ring collapses to zero breadth at the stem and the post, 2 of the 1998
   facets lie **exactly in the centreplane y = 0**. Every off-plane apex *inverts*
   those cells (3 inverted tets at 5 % beam offset) while an on-plane apex sees them
   **edge-on** (zero-volume cells). The old rule demanded strictly positive cells, so
   the hull was refused and reported as "not star-shaped" — the docstring claim about
   hulls was simply false.

2. **Replaced the strict-positivity heuristic with a certificate.** A candidate apex
   is accepted when (a) **no** tetrahedron is inverted and (b) the tets it spans
   reproduce the surface volume to 1e-9 relative. Edge-on facets span zero-volume
   cells; those cells are dropped and counted in the new
   `TetrahedralMesh.edge_on_facets`. Dropping zero-measure cells is exact — the
   volume identity in (b) is what proves it, and the certificate is what a
   non-star body (torus, new regression test) cannot satisfy.
   Among accepted candidates the one leaving the fewest edge-on facets wins, so the
   existing cube/sphere/capsule results are unchanged (`edge_on_facets = 0`).

3. **Extended the matrix to 24 combinations** (6 formats × 4 geometries) with the
   hull in every format, and re-ran the evidence generator into mission-owned
   directories (`--out` / `--fig` options added, so the m17802 evidence is untouched).

4. **Pinned the residual limitation** as a characterisation test rather than prose:
   the delivered hull volume mesh drops 2 skin facets, leaving 6 unpaired rim edges
   (slit area 1.3021e-02 m² of the hull's 2.15386e+01 m² skin, measured). Volume conservation is exact;
   a solver-ready watertight volume mesh is **WITHHELD** for the hull.

## 3. Results

### 3.1 Round-trip matrix — 24 combinations, 24 PASS

| Format | Geometry | Topology sent | V original [m³] | V round-trip [m³] | rel. err | Topology | Orientation | Status |
|---|---|---|---|---|---|---|---|---|
| stl | cube | triangle (8v/12c) | 1.000000000 | 1.000000000 | 0.0e+00 | ok | ok | PASS |
| stl | sphere | triangle (1106v/2208c) | 4.158971436 | 4.158971436 | 0.0e+00 | ok | ok | PASS |
| stl | capsule | triangle (770v/1536c) | 3.639460645 | 3.639460645 | 0.0e+00 | ok | ok | PASS |
| stl | wigley_hull | triangle (1001v/1998c) | 2.771222150 | 2.771222150 | 0.0e+00 | ok | ok | PASS |
| obj | cube | triangle (8v/12c) | 1.000000000 | 1.000000000 | 0.0e+00 | ok | ok | PASS |
| obj | sphere | triangle (1106v/2208c) | 4.158971436 | 4.158971436 | 0.0e+00 | ok | ok | PASS |
| obj | capsule | triangle (770v/1536c) | 3.639460645 | 3.639460645 | 0.0e+00 | ok | ok | PASS |
| obj | wigley_hull | triangle (1001v/1998c) | 2.771222150 | 2.771222150 | 0.0e+00 | ok | ok | PASS |
| vtk | cube | triangle (8v/12c) | 1.000000000 | 1.000000000 | 0.0e+00 | ok | ok | PASS |
| vtk | sphere | triangle (1106v/2208c) | 4.158971436 | 4.158971436 | 0.0e+00 | ok | ok | PASS |
| vtk | capsule | triangle (770v/1536c) | 3.639460645 | 3.639460645 | 0.0e+00 | ok | ok | PASS |
| vtk | wigley_hull | triangle (1001v/1998c) | 2.771222150 | 2.771222150 | 0.0e+00 | ok | ok | PASS |
| gmsh | cube | triangle (8v/12c) | 1.000000000 | 1.000000000 | 0.0e+00 | ok | ok | PASS |
| gmsh | sphere | triangle (1106v/2208c) | 4.158971436 | 4.158971436 | 0.0e+00 | ok | ok | PASS |
| gmsh | capsule | triangle (770v/1536c) | 3.639460645 | 3.639460645 | 0.0e+00 | ok | ok | PASS |
| gmsh | wigley_hull | triangle (1001v/1998c) | 2.771222150 | 2.771222150 | 0.0e+00 | ok | ok | PASS |
| su2 | cube | tetra (9v/12c) | 1.000000000 | 1.000000000 | 0.0e+00 | ok | ok | PASS |
| su2 | sphere | tetra (1107v/2208c) | 4.158971436 | 4.158971436 | 0.0e+00 | ok | ok | PASS |
| su2 | capsule | tetra (771v/1536c) | 3.639460645 | 3.639460645 | 0.0e+00 | ok | ok | PASS |
| su2 | wigley_hull | tetra (1002v/1996c) | 2.771222150 | 2.771222150 | 0.0e+00 | ok | ok | PASS |
| cgns | cube | tetra (9v/12c) | 1.000000000 | 1.000000000 | 0.0e+00 | ok | ok | PASS |
| cgns | sphere | tetra (1107v/2208c) | 4.158971436 | 4.158971436 | 0.0e+00 | ok | ok | PASS |
| cgns | capsule | tetra (771v/1536c) | 3.639460645 | 3.639460645 | 0.0e+00 | ok | ok | PASS |
| cgns | wigley_hull | tetra (1002v/1996c) | 2.771222150 | 2.771222150 | 0.0e+00 | ok | ok | PASS |

Measured with meshio 5.3.5, numpy 2.4.6, Python 3.11.6 on
`Linux-6.6.52-1.7.lns23.loongarch64-loongarch64-with-glibc2.38`. Per-file sha256 and the exported files themselves are in
`artifacts/m18268_meshio_matrix_hull/matrix18.json`, plus the exported mesh
files themselves under `artifacts/m18268_meshio_matrix_hull/files/`.

### 3.2 Geometry gate (a round-trip claim is only meaningful on a valid solid)

| Geometry | Vertices | Faces | Watertight | V discrete [m³] | V analytic [m³] | rel. err | Tets | Edge-on facets dropped | ΣV_tet − V_surface |
|---|---|---|---|---|---|---|---|---|---|
| cube | 8 | 12 | True | 1.000000 | 1.000000 | 0.00 % | 12 | 0 | 0.00e+00 |
| sphere | 1106 | 2208 | True | 4.158971 | 4.188790 | 0.71 % | 2208 | 0 | -4.26e-14 |
| capsule | 770 | 1536 | True | 3.639461 | 3.665191 | 0.70 % | 1536 | 0 | -4.00e-15 |
| wigley_hull | 1001 | 1998 | True | 2.771222 | 2.777778 | 0.24 % | 1996 | 2 | 1.11e-14 |

The hull's tetrahedralisation now closes to `1.11e-14` m³ of the
surface integral (floating-point noise), i.e. the 1996 emitted tets tile the body exactly.

### 3.3 Figures

| Figure | Path | Content |
|---|---|---|
| Volume relative error, all combinations | `artifacts/m18268_meshio_matrix_hull/volume_relative_error.png` | per-combination relative volume error, log scale, with the 1e-12 lossless line |
| Status matrix | `artifacts/m18268_meshio_matrix_hull/status_matrix.png` | format × geometry grid with PASS/FAIL, topology and cell counts |
| Geometry gate | `artifacts/m18268_meshio_matrix_hull/geometry_gate.png` | discrete vs analytic volume per fixture, hull included |

## 4. Notes / limitations

* **The claim that was wrong.** "Star tessellation covers ship hulls whose centroid
  sees the whole shell" was untested when it was written. Against `tensorcad.wigley`
  it fails for a structural reason (centre-plane tip facets), and the fix is a
  certificate, not a looser tolerance: ~5 % of beam off-centre-plane inverts cells and
  is still refused (`test_hull_off_centreplane_apex_is_refused`).
* **Open boundary slit.** The hull volume mesh boundary is the hull skin minus the 2
  edge-on facets (6 unpaired rim edges). Volume-exact, not solver-watertight. A real
  volume mesher (gmsh-style) remains the right tool for a CFD-grade hull volume mesh;
  what this module now delivers for a hull is a volume-exact, round-trip-stable tet
  mesh.
* **CGNS scope.** Round trip through meshio is REAL; conformance for a third-party
  CGNS reader stays WITHHELD (no `cgnslib` here, and meshio writes a single uniform
  tetra block without the full CGNS/HDF5 metadata).
* **STL is lossless in this environment, by accident of the writer.** All STL rows
  measure rel. error 0 because meshio 5.3.5 writes **ASCII** STL; a binary
  (float32) STL writer would lose ~1e-7 relative volume. The 1e-3 STL tolerance in the
  tests is the right guard, but the earlier "STL = float32 detour" figure caption was
  not what the measurement showed and has been corrected.
* Prior docs numbers for a "Wigley" fixture (V ≈ 0.0344 m³ vs a claimed 0.075 m³) came
  from a hand-rolled, non-manifold shelf and are superseded by this run; the fixture
  was replaced earlier by the capsule and is now joined by the package generator's hull.

## 5. Error judgement against references

| Quantity | Reference | Measured | Error | < 3 % ? |
|---|---|---|---|---|
| Hull displaced volume (discrete surface) | `4LBT/9` = 2.777778 m³ (closed form) | 2.771222 m³ | 0.24 % | **yes** |
| Hull volume after tetrahedralisation | surface integral 2.771222 m³ | 2.771222 m³ | 4.0e-15 (rel.) | **yes** (exact to machine precision) |
| Hull volume after SU2 / CGNS / STL / OBJ / VTK / Gmsh round trip | surface integral, per format | identical to 9 printed digits | 0.0e+00 measured | **yes** |
| Sphere / capsule volume vs analytic (fixture validity) | 4π/3, πr²L+4/3πr³ | see §3.2 | 0.71 % / 0.70 % | **yes** |

No quantity in this work is above the 3 % acceptance line. The one thing that is
*not* claimed: solver-grade watertightness of the hull volume mesh (WITHHELD, see §4).

## 6. Reproduce

```bash
python scripts/meshio_matrix18.py --out artifacts/m18268_meshio_matrix_hull/files \
    --fig artifacts/m18268_meshio_matrix_hull
cp artifacts/m18268_meshio_matrix_hull/files/matrix18.json artifacts/m18268_meshio_matrix_hull/
python -m pytest tests/test_meshio_format_matrix.py -q     # 141 passed, 4 skipped
```
