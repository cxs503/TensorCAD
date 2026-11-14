# Verified case: `sphere_mean_curvature`

Anchor `CASE-CAD-CURVATURE-SPHERE` — **discrete mean curvature of a sphere against the
analytic `2/R`**, criterion **2 %**.

| file | what it is |
|---|---|
| `run_benchmark.py` | recomputes gates G1–G6 plus diagnostic G7, writes `result.json`, `gates_verified.json`, `figures/*.png` |
| `result.json` | every reference/measured pair, both convergence series, the radius sweep, environment and source hashes |
| `gates_verified.json` | one line per gate: criterion, reference, measured, error, verdict |
| `figures/` | `convergence.png`, `readings.png`, `error_cloud.png` (per-vertex error field on the sphere) |

```
python benchmarks/verified/sphere_mean_curvature/run_benchmark.py           # ~45 s, exits non-zero on a broken gate
python benchmarks/verified/sphere_mean_curvature/run_benchmark.py --quick    # ~15 s smoke run
```

## Reference

`kappa_1 = kappa_2 = 1/R`, so the quantity the cotangent Laplacian estimates,
`S = kappa_1 + kappa_2`, is `2/R` (sum convention: `|K| = 2|H|`, see
`tensorcad.curvature.vertex_mean_curvature`). The value is **evaluated** by
`reference_sum(R)`, never copied: do Carmo *Differential Geometry of Curves and
Surfaces* (1976) Sec. 3-3 Example 2 for the closed form; Meyer/Desbrun/Schroeder/Barr
(VisMath 2003) Sec. 3.3 for the discrete convention. Independent corroboration:
`int_S S dA = 8 pi R`, i.e. `8 pi` at `R = 1`, which is the `E_bending` closed form
already gated by `benchmarks/verified/curvature_fairness` (G2/G3).

## Two readings, both reported

* **integral** = area-weighted mean of the discrete `S`. The cotangent operator is a
  discrete divergence, so this is the reading that is insensitive to triangle shape.
* **pointwise** = median `|S_i|` over vertices. Second order on a structured
  triangulation, non-convergent on a sliver-ridden marching-tet mesh.

## Gates

| gate | quantity | reference | tolerance |
|---|---|---|---|
| G1 | discrete `S`, area-weighted integral, marching tets level 48 (`1 m`) | `2/R` | 2 % |
| G2 | pointwise `S` median, structured UV sphere 16x32 | `2/R` | 2 % |
| G3 | observed order of the pointwise estimate on the structured family | 2 | ±0.5 absolute |
| G4 | field-path `kappa_1+kappa_2` (7 radii, finite-difference Hessian) | `2/R` | `1e-6` relative |
| G5 | `int_S S dA` over the mesh, level 48 | `8 pi R` | 2 % |
| G6 | area-weighted integral, 5 radii `0.1 … 10 m` at level 24 | `2/R` | 2 % each |
| G7 | **diagnostic**: raw pointwise median, marching tets level 32 | `2/R` | 2 % (reported, not gated) |

## Result (committed run)

`S = 2.00024 1/m` vs `2/R = 2.0` → **+0.012 %**, inside the 2 % criterion. Pointwise on
the same mesh: **50.1 %** raw, **1.94 %** after `taubin_smooth(10) + project_to_field(2)`.
Structured UV sphere pointwise: 1.48 % → 0.39 % → 0.10 % → 0.045 % over
`8x16 … 48x96`, observed order **1.95**.

Details, figures and the failure modes: `docs/sphere_mean_curvature.md`.