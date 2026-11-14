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