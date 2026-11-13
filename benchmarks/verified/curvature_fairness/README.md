# Verified case: `curvature_fairness`

Anchor `GAP-CAD-A3-CURVATURE-FAIRNESS` — curvature field (principal / Gaussian / mean)
and smoothing energy, measured against references that are not mesh estimates.

| file | what it is |
|---|---|
| `run_benchmark.py` | recomputes gates G1–G6 + diagnostic G7, writes `result.json`, `gates_verified.json`, `figures/curvature_fairness.png` |
| `result.json` | every reference/measured pair, the convergence series, the invariance probes, the three-case cross-check, environment and source hashes |
| `gates_verified.json` | one line per gate: criterion, verdict, failing rows |

## What is gated, and against what

| gate | quantity | reference | tolerance |
|---|---|---|---|
| G1 | sphere principal curvature, 7 radii (0.1 … 10 m) | closed form `k = 1/R` | `1e-6` 1/m absolute (the A3 clause) |
| G2 | sphere `E_bending` / `E_willmore` / `E_gauss`, field path, level 24 | `8π` / `4π` / `4π` | 0.1 % |
| G3 | the same three on a conditioned mesh (level 32) | `8π` / `4π` / `4π` | 6 % |
| G4 | Wigley smooth `E_bending`, `E_willmore` at 121×37 | `analytic_flank_energy` (Gauss-Legendre over the analytic flank) | 15 % |
| G5 | Wigley Richardson-extrapolated `E_bending`; the analytic reference's own order-doubling change | analytic flank; 0 | 0.5 %; `1e-12` |
| G6 | invariance of `E_bending` under scaling / cell reordering / 1.6 km translation | 0 | `1e-10` absolute |
| G7 | *diagnostic only*: three-case cross-check, median pointwise deviation from the exact field, plus the Wigley flank integral | exact analytic field / analytic flank half | 15 % reported, **not gated** |

G7 is deliberately not a gate. The acceptance clause asked for a `meshmagick`
comparison on three cases; `meshmagick` cannot be obtained in this environment
(`pip download` finds no version on the internal index, `git+https` is blocked), so
the harness runs the three cases against `trimesh`'s normal-cycle measures
(Cohen-Steiner & Morvan) and the closed forms instead, and records the substitution
in `result.json` → `cross_check.provenance`. Medians:

| case | H (TensorCAD mesh vs exact) | K (TensorCAD mesh vs exact) |
|---|---|---|
| sphere, level 24 | 1.58 % | 0.048 % |
| ellipsoid, level 24 | 12.83 % | **22.43 %** |
| Wigley flank patch (integral) | — | `E_bending` −12.87 %, area −13.61 % vs the analytic half |

The ellipsoid Gaussian-curvature row **breaks the 15 % standard** and is reported as
a failing diagnostic rather than absorbed: the cotangent-Laplacian `H` and the
angle-defect `K` stop being a consistent pair on a coarse marching-tet ellipsoid,
where `H² − K` can even go negative (`negative_discriminant_vertices` in
`result.json`). The field path (G2) and the conditioned sphere mesh path (G3) show
the estimators themselves are right; the mesh path on a coarse anisotropic body is
where the error lives.

## Reproduce

```bash
python benchmarks/verified/curvature_fairness/run_benchmark.py    # ~50 s, gates + figure + JSON
pytest tests/test_fairness_entry.py tests/test_fairness_energy.py -q
```

One command for the whole evidence bundle (benchmark + CLI records + three-case
cross-check + tests + sha256 manifest):

```bash
bash runs/m17465/repro_cmd.sh
```