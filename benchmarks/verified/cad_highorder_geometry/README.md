# Verified case: `cad_highorder_geometry`

Campaign node `p127w2-n01-cad-highorder-geometry`. TensorCAD as the ecosystem's p-order
geometry truth source: quadrature, normals, principal curvatures, wetted surface and `Cb` on
an exact rational sphere and an exact ellipsoid, each reported **with the polynomial order it
is accurate to** and cross-checked against an implementation that shares no code with it.

The report is `docs/highorder-geometry.md`; it is *generated* from `result.json`, so no number
in it is transcribed by hand.

## Contents

| file | what it is |
|---|---|
| `run_benchmark.py` | recomputes gates J1–J7, writes `result.json`, `gates_verified.json`, `figures/`, `artifacts/` |
| `render_report.py` | renders `docs/highorder-geometry.md` from `result.json` |
| `result.json` | every measurement, the environment, and the git/module/script/test hashes |
| `gates_verified.json` | one line per gate: criterion, verdict, pointer into `result.json` |
| `figures/*.png` | convergence (with the faceted negative control), curvature error maps, cross-validation and sensitivity |
| `artifacts/contract_*.npz`, `.h5`, `.json` | a real `(x, w, n, \|J\|)` element contract as consumed by a solver, in all three formats |

The three formats carry the same arrays and the same content hash: `result.json` records each
file's size and the hash a solver should see, so the committed bytes are checkable rather
than decorative. `matplotlib` and `h5py` are optional -- on a host without them the benchmark
records a skip for that output and still writes complete gate evidence (the CI/702 cluster
node has neither).

## Reproduce

```bash
python benchmarks/verified/cad_highorder_geometry/run_benchmark.py   # gates + figures + JSON
python benchmarks/verified/cad_highorder_geometry/render_report.py   # docs/highorder-geometry.md
pytest tests/test_highorder_geometry.py -q                           # 44 gates, ~4 s
```

## Measured on two independent hosts

Both runs, from the same working tree, on loongarch64 nodes (control plane 99 and cluster
node n315):

| gate | quantity | value |
|---|---|---|
| J1 | fitted quadrature order, p = 1…4 | 19.04, 12.03, 11.20, 10.47 (requirement 2p−1 = 1, 3, 5, 7) |
| J2 | curvature error / κ_max, ellipsoid | 1.884e-14 (gate 1e-3) |
| J2 | unit-normal deviation, ellipsoid | 1.310e-15 (gate 1e-3) |
| J3 | order-4 element volume error, 2×2 subdivisions | 2.06e-07 (gate 1e-6) |
| J4 | NURBS wetted surface vs elliptic truth | 1.33e-11 |
| J4 | STL wetted surface, 128 facets/direction | 2.88e-04 |
| J4 | voxel Cb, 128³ | 9.70e-05 |
| J5 | contract round-trip | bit identical, worst array difference 0.0 |
| J6 | negative control order (area, volume) | 1.982, 1.972 (gate ≤ 2) |
| J7 | `pytest tests/test_highorder_geometry.py` | 3.1 s on n315, 3.5 s on the control plane (gate 120 s) |

## Honest limitations

* Axis-aligned tensor-product patches only: no trimmed surfaces, no multi-patch stitching.
* The 1e-6 volume figure needs order-4 elements **and** a 2×2 subdivision per knot span; one
  element per span misses it (`7.5e-05`). The test suite asserts that, so the number cannot be
  quoted without its resolution.
* A curvature accuracy of 1% does **not** bound the wetted area to a fixed percentage: the
  measured ratio spans 0.55–2.9 depending on which perturbation produced the curvature error,
  which is why the area is gated directly instead of inferred (`result.json#/gates/J4/sensitivity`).
* The voxel `Cb` uses the partial-fill volume rule; the hard occupancy count is reported
  alongside and is about an order of magnitude coarser.