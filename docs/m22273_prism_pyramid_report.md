# Prism/Pyramid Mesh Quality — Re-verification Report (m22273)

## 1. Case Introduction

Ship CFD meshes use **prism (wedge6) layers** on hull surfaces for boundary-layer resolution
and **pyramids (5-node)** at hex/tet transitions. The `tensorcad.meshquality` library
previously only measured triangles, quadrilaterals, tetrahedra, and hexahedra. PR#581
(m19964) added full metric support for prisms and pyramids; it was merged via #460 but
the audit artifacts were not committed. This mission re-opens the deliverable, verifies
correctness end-to-end, and commits the artifacts.

### Geometry & Parameters

| Element | Nodes | Ideal geometry | Key closed-form values |
|---|---|---|---|
| Prism (wedge6) | 6 | Equilateral-triangle base (side=1), unit height | Volume = √3/4 ≈ 0.4330, min_sicn = 1.0, scaled_jacobian = sin60° ≈ 0.8660 |
| Pyramid (5-node) | 5 | Unit-square base, apex at height 1 | Volume = 1/3 ≈ 0.3333, min_sicn = 1.0, scaled_jacobian = √(2/3) ≈ 0.8165 |

Test meshes: ideal prism, ideal pyramid, sheared prism (non-affine), wedge prism
(boundary-layer collapse), warped pyramid (non-planar base), prism plate layer
(3-layer BL, 24 cells), pyramid cube fan (6 pyramids from cube centre), collapsed
prism layer, inverted pyramid fan.

## 2. Method

1. **Metrics implementation** (`tensorcad/meshquality/metrics.py`): `prism_metrics()`
   and `pyramid_metrics()` compute volume (conforming tet split), edge ratio, min/max
   angle, equiangle skewness, scaled Jacobian, and min SICN. The SICN normaliser is
   the ICN of the ideal element at each sample point, so an ideal prism/pyramid scores
   1.0. The pyramid apex (singular by construction) is excluded from sampling.
2. **Independent second route** (`tensorcad/meshquality/reference.py`):
   `reference_prisms()` and `reference_pyramids()` use finite-difference Jacobians and
   Cayley-Menger volumes — deliberately different arithmetic from the library path.
3. **Gate policy** (`tensorcad/meshquality/report.py`): `DEFAULT_LIMITS` extended with
   prism (min_sicn ≥ 0.1, scaled_jacobian ≥ 0.1, min_angle ≥ 10°) and pyramid
   (min_sicn ≥ 0.1, scaled_jacobian ≥ 0.0, min_angle ≥ 10°). Pyramid skewness is NOT
   gated because the ideal square pyramid has skewness ≈ 0.197, not 0.
4. **Audit tool** (`tools/prism_pyramid_audit.py`): Runs 5 test meshes, compares
   library vs second route, checks ideal-element closed forms, and writes
   `artifacts/mesh_quality_prism_pyramid/`.
5. **Test suite** (`tests/test_meshquality_prism_pyramid.py`): 21 tests covering
   closed-form values, rigid-motion invariance, inversion detection, flat-element
   handling, independent-route agreement, measure()/gate()/CLI integration, and
   exchange-format round-trips.

## 3. Results

### Ideal-element verification (library vs closed form)

| Element | Metric | Closed form | Library | |Δ| |
|---|---|---|---|---|
| Ideal prism | volume | 0.4330127018922193 | 0.4330127018922193 | 0 |
| Ideal prism | min_sicn | 1.0 | 1.0 | 0 |
| Ideal prism | scaled_jacobian | 0.8660254037844386 | 0.8660254037844387 | 1.1e-16 |
| Ideal prism | min_angle | 60.0 | 60.0 | 7.1e-15 |
| Ideal prism | skewness | 0.0 | 1.2e-16 | 1.2e-16 |
| Ideal prism | edge_ratio | 1.0 | 1.0 | 2.2e-16 |
| Ideal pyramid | volume | 0.3333333333333333 | 0.3333333333333333 | 0 |
| Ideal pyramid | min_sicn | 1.0 | 1.0 | 0 |
| Ideal pyramid | scaled_jacobian | 0.816496580927726 | 0.8164965809277261 | 1.1e-16 |
| Ideal pyramid | min_angle | 48.18968510422141 | 48.1896851042214 | 1.4e-14 |
| Ideal pyramid | skewness | 0.19683858159630982 | 0.19683858159631004 | 2.2e-16 |
| Ideal pyramid | edge_ratio | 1.224744871391589 | 1.224744871391589 | 0 |

**ideal_element_failures: []** (empty — all match to machine precision)

### Library vs independent second route (5 meshes)

| Mesh | Kind | n_cells | max|lib−ref| | Gate | Failing checks |
|---|---|---|---|---|---|
| prism_plate_layer | prism | 24 | 7.851e-12 | PASS | — |
| pyramid_cube_fan | pyramid | 6 | 2.415e-12 | PASS | — |
| sheared_prism_stack | prism | 3 | 1.937e-11 | PASS | — |
| collapsed_prism_layer | prism | 24 | 6.217e-10 | FAIL | prism.min_sicn, prism.scaled_jacobian, prism.min_angle |
| inverted_pyramid_fan | pyramid | 6 | 1.554e-12 | FAIL | pyramid.min_sicn, pyramid.scaled_jacobian |

Worst per-metric difference over all 5 meshes:

| Metric | max |library − second route| |
|---|---:|
| volume | 6.217e-10 |
| min_sicn | 1.258e-11 |
| scaled_jacobian | 1.937e-11 |
| min_angle | 1.066e-14 |
| skewness | 0.0 |
| edge_ratio, min_edge, max_edge, max_angle | 0.0 |

### Test suite

```
pytest tests/test_meshquality_prism_pyramid.py tests/test_meshquality.py -q
47 passed, 1 warning in 1.63s
```

All 21 prism/pyramid-specific tests and all 26 general meshquality tests pass.

## 4. Figures

See `artifacts/mesh_quality_prism_pyramid/quality.png` — the audit tool's summary plot
showing per-mesh metric distributions and gate verdicts.

## 5. Notes & Limitations

- **No parity against gmsh/meshmagick binaries.** Neither imports on this host.
  Definitions were transcribed and cross-validated by a second route; the
  parsed-upstream comparison remains open.
- **Pyramid min_sicn/scaled_jacobian sample 5 points (base corners + base centre),
  not Gauss points or the apex.** A metric sampling the singular apex would be
  constant zero.
- **Volume is the conforming tet split**, exact for planar faces. For a warped base,
  the 0-2 diagonal split is reported.
- **Non-orthogonality is not extended to prisms/pyramids.**
  `topology.non_orthogonality` still accepts only hex/tet.
- **FD-route precision is 6e-10 on degenerate cells** (collapsed prism), not 1e-15.
  Acceptance tolerance is 1e-6.
- **F1 — Pyramid winding matters.** The base must be wound counter-clockwise seen
  from the apex (reverse of a hex face). Wrong winding gives scaled_jacobian = -1.0
  with all other metrics identical.
- **F2 — Real boundary layers sit close to the min_sicn limit.** The 3-layer prism
  plate (aspect ratio up to 28) scores min_sicn = 0.1493, passing ≥ 0.1 with 49%
  margin but no margin if the layer is one refinement thinner. This is the intended
  reading; edge_ratio is NOT gated for prisms.

## 6. Error Assessment

| Quantity | vs Reference | Error | <3%? |
|---|---|---|---|
| Ideal prism volume | Closed form √3/4 | 0 (exact) | YES |
| Ideal prism min_sicn | Closed form 1.0 | 0 (exact) | YES |
| Ideal prism scaled_jacobian | Closed form sin60° | 1.1e-16 | YES |
| Ideal pyramid volume | Closed form 1/3 | 0 (exact) | YES |
| Ideal pyramid min_sicn | Closed form 1.0 | 0 (exact) | YES |
| Ideal pyramid scaled_jacobian | Closed form √(2/3) | 1.1e-16 | YES |
| Library vs 2nd route (all metrics, all meshes) | Independent FD implementation | max 6.2e-10 | YES |

**All key quantities agree with their references to far better than 3% — the worst
residual is 6.2e-10 on a degenerate cell, i.e. 1e-8 of the 1% band.** The quantities
are algebraically identical, not approximately equal.
