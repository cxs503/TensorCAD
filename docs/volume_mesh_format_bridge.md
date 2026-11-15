# Volume-mesh → tetrahedra → SU2/CGNS: closing the format-matrix workaround

**Mission**: m19950 — anchor `GAP-CAD-B2-MESH-FORMAT-MATRIX`
**Module**: `tensorcad.meshio_exchange.volume_mesh_to_tetrahedral`
**Date**: 2026-11-14

## 1. The gap this closes

`GAP-CAD-B2-MESH-FORMAT-MATRIX` was delivered by m17350/m17802: 6 formats × 3 geometries,
round-trip volume/topology/orientation preserved. Its report stated a residual limit, in
its own words:

> *Star tets are not a general tetrahedraliser. […] The package's own `hull.build_mesh`
> output is refused […] For such bodies use `tensorcad.volume_mesh` (body-fitted hex → tet)
> and feed the tetrahedra instead.*

That remedy was **unbacked**. `tensorcad.volume_mesh` does build the body-fitted grid and
does Kuhn-split each hexahedron into six tetrahedra, and `meshio_exchange.TetrahedralMesh`
is the only object SU2/CGNS accept — but nothing converted one into the other. The sentence
had exactly the defect it was written to replace: a promise with no code and no test behind
it. For a real hull there was therefore **no path at all** into the volume formats:

| body | star tessellation (`tetrahedralise_surface`) | body grid (`volume_mesh`) |
|---|---|---|
| cube / sphere / capsule | exact | n/a |
| package hull (`hull.build_mesh`) | **refused** (rings collapse at stem/post) | built — but not exportable |

## 2. What was added

`volume_mesh_to_tetrahedral(volume_mesh, *, units, coordinate_frame, tolerance)` —
~80 lines, no new dependency (`volume_mesh` is pure standard library, meshio stays optional
and is imported only by the write/read path).

What it does, and what it refuses to do:

* **Classifies every tetrahedron on its signed volume.** Strictly positive → kept.
  Negative → `MeshioFormatError` naming the cell and its signed volume. The adapter never
  reorders a cell to make it positive: a negative cell is a builder defect, and silently
  fixing it is how a solver ends up with an inverted mesh that converges to nonsense.
* **Drops zero-volume cells and counts them** in the new
  `TetrahedralMesh.degenerate_cells` field. Every structured hull grid contains them —
  the keel line and both end stations collapse to zero breadth, so whole cell columns have
  no measure. Dropping is exact (zero measure contributes zero), and the count is carried
  so it is visible rather than silent.
* **Asserts conservation instead of hoping for it.** The kept tetrahedra must reproduce
  `volume_mesh.volume()` within a scale-aware tolerance (default `1e-12 · V`); otherwise it
  raises. A silent loss of measure fails loudly.

## 3. Measured results

`python examples/meshio_volume_bridge/run_example.py` →
`runs/m19950_volume_bridge/bridge_matrix.json` (every number below is read from that file).

| body | hexes | tets kept | zero-volume dropped | V tets [m³] | V analytic [m³] | split vs grid | vs analytic | SU2 | CGNS |
|---|---|---|---|---|---|---|---|---|---|
| box 2×1×0.6 | 72 | 432 | 0 | 1.200000000000 | 1.200000000000 | 1.1e-15 | 1.5e-15 | 0.0 | 0.0 |
| Wigley body L10×B1×T0.625 | 1536 | 8280 | 936 | 2.762123390480 | 2.762123390480 | 7.7e-15 | 8.0e-15 | 0.0 | 0.0 |

Three separate anchors, deliberately not conflated:

1. **Split conservation** — tetrahedra vs the hex grid's own measure. Machine precision
   (1e-15), as it must be: the six Kuhn paths tile a hexahedron exactly.
2. **Trapezoid integral** — `2·∫b dx dz`, exact for the piecewise-linear cells the body is
   built from. Machine precision. This is an independent code path (quadrature, not cells).
3. **Continuum closed form** — `(4/9)·B·L·T = 2.777778 m³` for
   `b = (B/2)(1−(2x/L)²)(1−(z/T)²)`. Measured **0.564 %** low. That is the *geometry
   discretisation* of the body (piecewise-linear breadth vs a parabola, 25 stations ×
   9 levels), **not** an error of the bridge; it is reported separately for exactly that
   reason. Under 3 %, and it shrinks with refinement.

Format round trip through meshio preserves volume **exactly** (0.0 relative, both SU2 and
CGNS), with cell and vertex counts unchanged and all tetrahedra still positive.

## 4. Reproduce

```bash
python examples/meshio_volume_bridge/run_example.py      # table + JSON + figure, exits non-zero if any claim fails
python -m pytest tests/test_volume_mesh_format_bridge.py -q   # 14 passed
```

## 5. Limitations, stated not hidden

1. **CGNS metadata is still WITHHELD**, unchanged from the anchor's report: meshio's writer
   emits a single uniform tetra block without `CGNSLibraryVersion`/`ElementType`/`ZoneType`.
   Round-trip through meshio is REAL; conformance for a third-party reader is unverified —
   no `cgnslib` on this host. Do not claim "valid CGNS file" without that check.
2. **SU2/CGNS remain volume-only.** Surface triangle meshes are still refused, which is a
   property of the formats, not a missing feature.
3. **VTK is not reachable from here.** `vtk` is classified as a *surface* format in this
   module, so the topology gate refuses a tetrahedral mesh for it. meshio's VTK writer would
   happily take tets; changing that classification is a separate, larger decision (it would
   change the meaning of the existing 18-combination matrix) and was left alone.
4. **The grid must be split at build time** (`split_tets=True`, the default). A grid built
   with `split_tets=False` is refused with that remedy in the message rather than re-split
   here, so there is exactly one place that decides how a hexahedron is cut.

## 6. Second pass (m22029): the file a solver can actually load

### 6.1 The residual gap, measured

Section 2 delivered *cells*. It did not deliver a *boundary*, and a volume format without one
is not a solver file. Measured on the box of §3 (`test_unmarked_su2_has_no_markers`):

```
NDIME= 3
NPOIN= 140
NELEM= 432
NMARK= 0          <-- no marker at all
```

The cause is not carelessness, it is a missing binding: meshio's SU2 writer emits
`MARKER_ELEMENTS` only for triangle blocks that carry integer cell data, and
`TetrahedralMesh` is tetrahedra and nothing else. So `write_meshio(..., file_format="su2")`
writes a file that a solver config cannot bind `MARKER_HEATFLUX= (hull, 0.0)` to.

### 6.2 What was added

`tensorcad.meshio_exchange`:

* `volume_mesh_patches(volume_mesh, *, tetrahedral=None, tolerance=None) -> VolumePatches`
  — **welds coincident points**, takes the outer skin of the surviving tetrahedra, cancels
  the zero-thickness slit, and names the remaining triangles from the grid's own station and
  level planes (`symmetry`, `bottom`, `top`, `fore`, `aft`, and `hull` as the catch-all for
  the body surface). Refuses to return a skin that is not watertight or whose divergence-
  theorem volume differs from the tetrahedra's.
* `write_su2_volume_mesh(path, volume_mesh, *, patches=None) -> Path` — writes the SU2 ASCII
  mesh directly, with **string** `MARKER_TAG`s. It is written here rather than through meshio
  because meshio cannot keep string tags: it renumbers them 1, 2, 3… and warns as it goes.
* `read_su2_markers(path) -> dict[str, tuple[Tri, ...]]` — reads the marker sections back
  with the tags as names, closing the loop on a file meshio would flatten.

Why welding is not cosmetic: `volume_mesh` keeps a distinct node per
`(station, level, breadth)` index even where the breadth collapses. Dropping the zero-volume
cells then exposes the *interior* faces of the surviving cells. On the Wigley body of §3,
**1050 of the 2184 triangles that look like boundary are that
slit** — coincident pairs of a zero-thickness sheet. Exporting them would give a solver a
doubly covered, non-manifold wall. Welding merges 328 coincident points, the pairs cancel, and what is
left is a closed surface.

### 6.3 Measured results

`python examples/meshio_volume_bridge/run_example.py` →
`runs/m19950_volume_bridge/bridge_matrix.json` (numbers below read from that file; the two
`.su2` files it writes are in the same directory).

| body | boundary tri (split) | skin tri | slit tri cancelled | points welded | NPOIN | NELEM | NMARK | hull area [m²] | anchor [m²] | rel | watertight |
|---|---|---|---|---|---|---|---|---|---|---|---|
| box 2x1x0.6 | 216 | 216 | 0 | 0 | 140 | 432 | 5 | 2.400000 | 2.400000 | 0.00% | yes |
| Wigley body L10xB1xT0.625 | 2184 | 1134 | 1050 | 328 | 1696 | 8280 | 2 | 14.839762 | 14.879052 | 0.26% | yes |

Patch table:

| body | patch | triangles | area [m²] |
|---|---|---|---|
| box 2x1x0.6 | `bottom` | 48 | 2.000000000 |
| box 2x1x0.6 | `top` | 48 | 2.000000000 |
| box 2x1x0.6 | `fore` | 36 | 0.600000000 |
| box 2x1x0.6 | `aft` | 36 | 0.600000000 |
| box 2x1x0.6 | `hull` | 48 | 2.400000000 |
| Wigley body L10xB1xT0.625 | `bottom` | 368 | 6.655092593 |
| Wigley body L10xB1xT0.625 | `hull` | 766 | 14.839761791 |

Four independent checks, kept separate:

1. **Closure.** Every edge of the skin is used exactly once in each direction, and the
   divergence-theorem volume of the skin reproduces the tetrahedra's volume to
   8.7e-15 (Wigley) / 1.9e-16 (box). A slit face would fail both.
2. **Box areas are exact.** A box has closed-form sides: `bottom = top = L·B = 2.0`,
   `fore = aft = B·H = 0.6`, `hull = 2·L·H = 2.4`, total surface `7.6 m²`. Measured to
   1.9e-16 relative — a zero-error control for the patch split and for the winding.
3. **Wigley area against a quadrature of the breadth function.**
   `2·∬sqrt(1 + b_x² + b_z²) dx dz = 14.879052 m²` (midpoint, 400×200, converged to 1e-6);
   the mesh's `hull` patch measures `14.839762 m²`, **0.26% low**. That is the piecewise-linear
   body discretisation (25 stations × 9 levels), not a defect of the export: refining the grid
   25×9 → 49×17 → 97×25 closes the gap 0.264 % → 0.066 % → 0.023 %.
4. **Independent reader.** meshio reads the file we wrote and recovers 8280 tetrahedra whose
   volume equals the patch object's to 1e-12, plus 1134 boundary triangles equal, as
   geometry, to ours. Our own reader recovers the tags as *names*
   (`bottom, hull` for the Wigley body, `bottom, top, fore, aft, hull` for the box).

### 6.4 What is still withheld

* **Solver-side acceptance is WITHHELD.** No `SU2_CFD` binary exists on this host (checked:
  `which SU2_CFD` empty), so "SU2 runs this file" is not claimed. What is claimed and tested
  is that the file follows the SU2 mesh-file layout (element types 10 and 5, 0-based node
  indices, `NMARK`/`MARKER_TAG`/`MARKER_ELEMS` sections) and that two independent parsers
  read it back consistently.
* **String tags are a deliberate deviation from this package's meshio path.** meshio
  renumbers them; SU2's specification uses names. Code that must round-trip through meshio
  should expect `su2:tag` integers, not the names — use `read_su2_markers` when the names
  matter.
* **Patches are grid planes, not physics.** `hull` is "the skin that is not on the
  centreplane or on a station/level cap"; a real case would merge or split patches
  per component (rudder, skeg, deck) and that is per-case configuration, not geometry.
* The CGNS metadata limitation of §5.1 and the VTK classification of §5.3 are unchanged.

### 6.5 Reproduce

```bash
python -m pytest tests/test_volume_mesh_su2_markers.py -q   # the marker contract
python examples/meshio_volume_bridge/run_example.py         # table + JSON + 4-panel figure
```
