# Benchmark: parametric hull geometry and hydrostatics (ref vs got)

Node `p316w1-n03-cad-parametric-hull` (repository `TensorCAD`).
This directory is the committed run record: it is regenerated, never edited by hand.

```
python benchmarks/verified/parametric_hull/run_benchmark.py
```

The run (a) executes the verification campaign
[`examples/parametric_hull_dag/run_verification.py`](../../../examples/parametric_hull_dag/run_verification.py),
which measures every acceptance criterion of the node, and (b) drives the command-line
entry point end to end (`python -m tensorcad.parametric_hull_cli wigley ... --gmsh ...`),
then re-reads the exported gmsh 2.2 file and reports the volume it carries. Outputs:

| file | what it is |
|---|---|
| `result.json` | the ref-vs-got table below, the CLI handoff record, gate verdicts |
| `verification.json` | the full criterion record (convergence study, quality scan, Cb sweep) |
| `wigley_L10.msh` | the exported gmsh 2.2 surface (`sha256 b31f72e19a59f331...`) |
| `wigley_L10.json` | the CLI record for that export |
| [`docs/parametric_hull/figures/`](../../../docs/parametric_hull/figures/) | the six report figures |

Reference hull: L = 10.0 m, B = 1.0 m,
T = 0.625 m (Wigley III, Cb = 4/9). Mesh 97 x 16
(mesh_density_for(L, T, n_v=16) -> dx = 2.5 dz).

## ref vs got

| case | quantity | reference | measured | error | threshold | verdict |
|---|---|---|---|---|---|---|
| `a_wigley_ordinates` | max relative error of hull ordinates y(x,z) at 10 stations x 5 waterlines | `y = (B/2)(1-(2x-L)^2/L^2)(1-(z/T)^2), Wigley III closed form` | `3.42798e-05` | `3.4280e-05` | `0.001` | PASS |
| `b_wigley_displacement` | displaced volume [m^3] | `V = p q L B T/((p+1)(q+1)) = 2.777778` | `2.77367` | `-1.4793e-03` | `0.005` | PASS |
| `c_wigley_waterplane` | waterplane area [m^2] | `Aw = p L B/(p+1) = 6.666667` | `6.66595` | `-1.0793e-04` | `0.01` | PASS |
| `d_suboff_offsets` | max relative radius error over 20 SUBOFF stations | `DARPA SUBOFF four-segment radius table (Groves, Huang & Chang 1989)` | `0.000931037` | `9.3104e-04` | `0.005` | PASS |
| `e_min_scaled_jacobian` | minimum scaled Jacobian over the hull triangles | `criterion floor 0.30 (equilateral = 0.866)` | `0.365326` | `6.5326e-02` | `0.3` | PASS |
| `e_max_aspect_ratio` | maximum triangle aspect ratio | `criterion ceiling 5.0` | `2.73728` | `2.2627e+00` | `5` | PASS |
| `f_parametric_cb_sweep` | max relative volume error over the Cb sweep (Cb 0.35 -> 0.55) | `closed-form volume at each Cb` | `0.00243237` | `2.4324e-03` | `0.005` | PASS |
| `cli_gmsh_volume` | displaced volume of the exported gmsh file [m^3] | `V = 4/9 L B T = 2.777778 for L=10, B=1, T=0.625` | `2.77367` | `-1.4793e-03` | `0.005` | PASS |

`all_passed = True`; observed order of accuracy
p = 2.07 (Richardson-extrapolated displacement error
0.0e+00).

## What the reference is, and what it is not

* Wigley III closed-form hydrostatics (`V = p q L B T/((p+1)(q+1))`, `Aw = p L B/(p+1)`):
  textbook values, pinned independently in `tests/test_parametric_hull_dag.py`.
* DARPA SUBOFF: the published four-segment radius table of Groves, Huang & Chang (1989),
  transcribed in `tensorcad/suboff.py` and pinned by `tests/test_suboff.py`.
* **No experimental data is used, and none applies**: geometry, displacement and
  waterplane area of these bodies are closed-form quantities, not measured ones. The
  external anchor is the published SUBOFF offset table; everything else is exact
  mathematics. A hydrodynamic quantity (resistance, sinkage) would need towing-tank data
  and is out of this module's scope.
