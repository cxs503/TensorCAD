# GAP-CAD-B1 — STEP AP203/AP214 read pipeline (`tensorcad.step_brep`)

Status: **delivered**, evidence REAL for the supported entity table, WITHHELD for
third-party-kernel cross-validation (see Limitations). Last verified: 2026-11-13,
mission m17349, repository TensorCAD, branch `agent/m17349`.

## 1. What the gap was

`tensorcad.exchange.read_step_facets` reads only `FACETED_BREP`: a triangle soup that
somebody else already tessellated. A real CAD handover is
`MANIFOLD_SOLID_BREP` / `ADVANCED_FACE` over *analytic* surfaces, and that was refused
outright (`_BREP_ENTITIES` guard). The warehouse therefore could not open the file it is
most likely to be sent, and no volume could be measured from it.

## 2. What was built

`tensorcad/step_brep.py` (pure Python, no OCP, no numpy):

* ISO-10303-21 parser: comments, strings, complex instances, typed parameters
  (`LENGTH_MEASURE`), enumeration flags — `Enum` is a `str` subclass, so every boolean
  goes through `_flag()`; `bool(Enum('.F.'))` is `True` and that bug inverted two faces
  and produced a mesh that looked plausible and had the wrong volume.
* Geometry: `AXIS2_PLACEMENT_3D` frames; `PLANE`, `CYLINDRICAL_SURFACE`,
  `CONICAL_SURFACE`, `SPHERICAL_SURFACE`, `TOROIDAL_SURFACE` with analytic
  `evaluate`/`inverse`; `LINE`, `CIRCLE`, `ELLIPSE`, `B_SPLINE_CURVE_WITH_KNOTS`
  (de Boor, rational or not) for edge sampling.
* Trim handling: `EDGE_LOOP -> ORIENTED_EDGE -> EDGE_CURVE`; the traversal order is
  `ORIENTED_EDGE.orientation XOR EDGE_CURVE.same_sense`, which is what makes a seam
  edge pushed twice come back the way it went.
* Tessellation: parameter-space triangulation (hole bridging + ear clipping, with an
  area-coverage assertion), a strip grid for full-revolution trims, and longest-edge
  bisection until the angular span is below `angular_tolerance`.
* Measurement: divergence-theorem volume, watertightness, Euler characteristic, and an
  honest `volume_m3 = None` whenever the mesh is partial or open.
* `exchange.read_step()` now dispatches: `FACETED_BREP` to the old reader, everything
  else to the analytic reader; `supported_routes()["step_brep_read"]` was updated from
  `NOT_SUPPORTED` to the real level *for the named table*.

## 3. Coverage census (machine-readable: `artifacts/gap-cad-b1-step-reader/coverage.json`)

| file | schema | instances | faces | solids | supported surfaces | supported curves | refused types |
|---|---|---|---|---|---|---|---|
| `bushing.step` | AUTOMOTIVE_DESIGN { 1 0 1030 | 132 | 4 | 1 | PLANEx2, CYLINDRICAL_SURFACEx2 | LINEx2, CIRCLEx8 | - |
| `l_bracket.step` | AUTOMOTIVE_DESIGN { 1 0 1030 | 322 | 8 | 1 | PLANEx8 | LINEx36 | - |
| `plate_with_hole.step` | AUTOMOTIVE_DESIGN { 1 0 1030 | 277 | 7 | 1 | PLANEx6, CYLINDRICAL_SURFACEx1 | LINEx25, CIRCLEx4 | - |
| `stepped_shaft.step` | AUTOMOTIVE_DESIGN { 1 0 1030 | 143 | 5 | 1 | PLANEx2, CYLINDRICAL_SURFACEx2, CONICAL_SURFACEx1 | LINEx3, CIRCLEx8 | - |
| `2010_30p30n_thik_te_18inches.stp` | CONFIG_CONTROL_DESIGN | 2283 | 0 | 0 | - | B_SPLINE_CURVE_WITH_KNOTSx11 | - |

Supported surfaces: PLANE, CYLINDRICAL_SURFACE, CONICAL_SURFACE, SPHERICAL_SURFACE,
TOROIDAL_SURFACE. Supported curves: LINE, CIRCLE, ELLIPSE, B_SPLINE_CURVE_WITH_KNOTS.

Refused **by name** (never guessed): `B_SPLINE_SURFACE[_WITH_KNOTS]`,
`RATIONAL_B_SPLINE_SURFACE`, `SURFACE_OF_REVOLUTION`, `SURFACE_OF_LINEAR_EXTRUSION`,
`RECTANGULAR_TRIMMED_SURFACE`, `CURVE_BOUNDED_SURFACE`, `OFFSET_SURFACE`, `POLYLINE`,
`TRIMMED_CURVE`, `COMPOSITE_CURVE`, `SURFACE_CURVE`, `SEAM_CURVE`. A face whose single
closed boundary loop on a doubly curved surface is ambiguous (two caps fit the same
loop) is also refused: the volume is then `None`, not a plausible number.

## 4. Volume acceptance (criterion: mesh volume vs CAD nominal < 1 %)

The three solid files under `artifacts/gap-cad-b1-step-reader/samples/` are AP214
`MANIFOLD_SOLID_BREP` files written by `tools/step_samples.py`, whose geometry is
defined in closed form there, so the nominal is a formula and not a copied number
(`nominal_volumes.json`). Provenance of each file is in `provenance.json`.

| file | nominal (m^3) | measured (m^3) | error | triangles | Euler | criterion |
|---|---|---|---|---|---|---|
| `plate_with_hole.step` | 0.000454867 | 0.000454926 | +0.0129 % | 228 | 0 | PASS |
| `stepped_shaft.step` | 0.000336936 | 0.000336147 | -0.2341 % | 420 | 2 | PASS |
| `bushing.step` | 0.000183783 | 0.000183353 | -0.2341 % | 424 | 0 | PASS |
| `l_bracket.step` | 1.8e-05 | 1.8e-05 | +0.0000 % | 20 | 2 | PASS |

All four are inside the 1 % band; the worst (0.234 %) is the inscribed-polygon bias of
the circular trims, which is predicted by `2*pi^2/(3*n^2)` for `n = 2*pi/0.1186 = 53`
samples per revolution — the measured error matches the prediction, it is not luck.

## 5. Convergence (why the residual error is trustworthy)

| angular step (rad) | plate | shaft | bushing |
|---|---|---|---|
| 0.48 | +0.1836 % | -3.3234 % | -3.3234 % |
| 0.24 | +0.0497 % | -0.9001 % | -0.9001 % |
| 0.12 | +0.0129 % | -0.2341 % | -0.2341 % |
| 0.06 | +0.0033 % | -0.0597 % | -0.0597 % |
| 0.03 | +0.0008 % | -0.0149 % | -0.0149 % |

Error ratio per halving: plate 3.999, shaft
3.999, bushing 3.999 against the
second-order value 4.0. Second-order convergence in a *geometric* quantity is the
evidence that the trim evaluator is right, not merely self-consistent.

## 6. The genuinely third-party file

`artifacts/gap-cad-b1-step-reader/real/2010_30p30n_thik_te_18inches.stp` — the BANC
30P-30N high-lift geometry, written by Siemens NX 7.5 (ST-Developer v11), schema
`CONFIG_CONTROL_DESIGN`, 2283 instances. Provenance is
recorded in `real/provenance.json` (source repo, upstream path, sha256).

* It parses: `['B_SPLINE_CURVE_WITH_KNOTS']`, `LENGTH_UNIT` resolved to
  **39.3701 file units per metre** (`CONVERSION_BASED_UNIT('INCH', ...)`),
  i.e. the inch conversion is honoured — the bounding box comes out at 0.556 m,
  consistent with an 18-inch stowed chord plus slat and flap.
* It is a `GEOMETRICALLY_BOUNDED_WIREFRAME_SHAPE_REPRESENTATION`: 0
  solids, 0 faces. Asking it for a mesh is refused with
  "no MANIFOLD_SOLID_BREP ... a wireframe or sheet-model file has no solid to measure",
  which is the honest answer, not a silent empty mesh.

## 7. Limitations and honest gaps

1. **No third-party CAD kernel on this machine** (OCP/FreeCAD have no loong64 wheels;
   `pip download cadquery-ocp` finds nothing). So the volumes are validated against
   closed-form nominals plus second-order convergence, not against OCCT's own
   `GProp` volume. This is exactly what `supported_routes()["step_file_external_validation"]`
   still reports as `WITHHELD`.
2. **The three solid samples are authored fixtures**, not downloaded CAD parts: no
   non-wireframe third-party STEP file exists anywhere on this box (the only `.stp`
   files are copies of the 30P30N wireframe). Their provenance says so openly.
3. **B-spline surfaces are refused.** They need a Newton inversion of a tensor-product
   patch; refusing is better than auditing the tessellator.
4. Sphere/torus faces that close on a single loop (pole caps) are refused as ambiguous.
5. `angular_step`/`angular_tolerance` are knobs with a measured effect, not magic
   numbers; default 0.12/0.10 rad gives the 0.23 % bias above.

## 8. Bugs this work fixed on the way (each found by a machine check)

* `bool(Enum('.F.')) == True` inverted `same_sense`/`orientation`; caught by the
  watertight + Euler check, not by inspection.
* Two independently-sampled faces of one rim did not weld; caught by
  `boundary_edges > 0`, fixed with a neighbour-cell weld.
* Ear clipping with a duplicated bridge vertex emitted **overlapping** triangles: the
  plate's area came out 1.79x and the bushing's 2.46x too large. Found by comparing the
  triangle-area sum with the parameter-region area; that comparison is now an assertion
  in `_triangulate_face`, so the failure mode cannot come back silently.
* An inch file was scaled by `25.4 * 1000` instead of `1000 / 25.4`; caught by the
  bounding-box check on the real NX file (0.00086 m instead of 0.457 m).

## 9. Reproduce

```bash
python tools/step_samples.py --out artifacts/gap-cad-b1-step-reader/samples   # writes the 4 STEP files + nominal volumes
python tools/step_reader_audit.py                                            # census + volumes + convergence + figures, exit 0 == pass
python -m pytest tests/test_step_brep_reader.py -q                            # 17 tests
python -m tensorcad.step_brep census artifacts/gap-cad-b1-step-reader/real/2010_30p30n_thik_te_18inches.stp
python -m tensorcad.step_brep mesh artifacts/gap-cad-b1-step-reader/samples/plate_with_hole.step --obj /tmp/plate.obj
```

## 10. Ecosystem friction (recorded, not hidden)

* The repo had no STEP fixture at all, and downloading one is blocked by the sandbox
  network policy; the only real CAD file available was inside another repository's
  (`TensorFVM`) test data, so the cross-repository read had to be done by copying the
  file with its provenance. A shared `ecosystem-cad-samples` location would have saved
  that.
* The full TensorCAD suite did not finish inside the 120 s shell window (2 % in 100 s),
  so the evidence here is the targeted suites (`tests/test_step_brep_reader.py`,
  `tests/test_exchange_verification.py`) plus the audit script.
