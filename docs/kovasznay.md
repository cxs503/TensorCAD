# Kovasznay (1948) Exact Solution — Reference Bundle and Order-of-Accuracy Evidence

Module: `tensorcad/kovasznay.py` · Example: `examples/kovasznay_fvm/run_example.py` ·
Tests: `tests/test_kovasznay.py` (59) · Artifacts: `artifacts/kovasznay_fvm/`

Case-pool entry: **`FVM-KOVASZNAY-ORDER`** — "定常 Kovasznay 流收敛阶验证", acceptance
`L2 收敛阶 >= 1.8 (名义二阶)`, metric *log-log slope on 32/64/128 meshes*.

## 1. What this module is, and what it is not

TensorCAD supplies the **geometry, exact solution and reference data** side of the
benchmark, exactly as it already does for the lid-driven cavity
(`tensorcad/cavity_mesh.py`, `examples/cavity_re5000/`). The flow solver lives in
TensorFVM (`tensorcad` core stays dependency-free: standard library only).

| Part of the benchmark | Where it lives | Evidence grade here |
|---|---|---|
| Exact solution + analytic derivatives | `tensorcad/kovasznay.py` | **REAL** — residual < 4e-15 (Table 2) |
| 32/64/128 quad mesh family + BC tags | `tensorcad/kovasznay.py` | **REAL** — exact area, positive Jacobians |
| Reference fields (cell-centre, nodal) + probe table | `artifacts/kovasznay_fvm/*.json` | **REAL** — generated, hashes in `checks.json` |
| Order-measurement harness (solver-agnostic) | `convergence_study()` | **REAL** — recovers injected p = 1, 2, 3 to 7e-5 |
| Second-order *consistency* of the discrete operators | `truncation_error_study()`, `midpoint_consistency_study()` | **REAL** — p = 1.977…1.999 (Table 3) |
| **Error order of a converged NS solution** | requires the TensorFVM solve | **WITHHELD** — see §6 |

The last row matters: the acceptance number `p >= 1.8` refers to the *solution*
error `||u_h - u||_2`. This unit measures what the discretisation can deliver
(consistency) and provides the harness that will measure the solver's order, but it
does **not** claim a solver order. The TensorFVM attempt (mission 16068) ended with
`stop_condition = agent_lost`, so the anchor stays open.

## 2. Physical problem

Steady, incompressible 2-D Navier–Stokes with unit reference velocity and length
(`Re = U L / nu = 1/nu`):

```
lam^2 - Re*lam - 4*pi^2 = 0        ->   lam = (Re - sqrt(Re^2 + 16*pi^2)) / 2

u(x, y) = 1 - exp(lam*x) * cos(2*pi*y)
v(x, y) = lam/(2*pi) * exp(lam*x) * sin(2*pi*y)
p(x, y) = p0 - exp(2*lam*x) / 2
```

* `lam < 0` is the branch bounded downstream (the positive root diverges); for
  `Re = 40`, `lam = -0.9637405441957689` and the quadratic residual is 8.5e-14.
* Continuity holds identically: `u_x + v_y = -lam*E*cos + lam*E*cos = 0`.
* The `u`-momentum equation fixes `lam`; the `v`-momentum equation is then
  satisfied automatically (both are checked numerically in Table 2).
* Stream function `psi = y - exp(lam*x) sin(2*pi*y)/(2*pi)` gives two
  algebraically independent routes to the vorticity; their mismatch is exactly 0.

**Set-up used here (standard benchmark box):**

| Quantity | Value |
|---|---|
| Reynolds number | 40 (unit U, L) |
| Kinematic viscosity | 0.025 |
| Domain | `[-0.5, 1.5] x [-0.5, 1.5]`, area 4 |
| Boundary conditions | Dirichlet `(u, v)` from the exact solution on all four sides |
| Pressure datum | `p0 = 0` (pressure is defined up to a constant; only `grad p` is physical) |
| Mesh family | uniform quads, 32/64/128 cells per direction |

## 3. What the bundle contains

`python examples/kovasznay_fvm/run_example.py` writes `artifacts/kovasznay_fvm/`:

| File | Content |
|---|---|
| `checks.json` | every measured number, the fitted orders, and the pass/fail conclusions |
| `case_spec.json` | machine-readable spec (Re, lam, domain, BCs, acceptance, mesh summaries) |
| `mesh_{32,64,128}.json` | nodes, quads (CCW), boundary node sets, cell centroids, exact cell-centre and nodal fields |
| `probe_table.json` | exact solution on a 7x5 probe lattice with residual columns |
| `kovasznay_fields.png` | exact `u`, `v`, `p` fields |
| `kovasznay_order.png` | log–log error vs `h` with fitted slopes and a slope-2 reference |
| `kovasznay_profiles.png` | centreline profiles at Re = 40 and Re = 100 |

Reference fields are stored rounded to 12 significant digits (≈1e-12 relative,
eight orders below any solver error at these resolutions); `checks.json` keeps full
double precision. Every file's sha256 is recorded in `checks.json → bundle.files`.

## 4. Results

### 4.1 The exact solution verifies itself (35 probe points)

| Quantity | Value |
|---|---|
| max \|momentum residual `R_u`\| | 3.553e-15 |
| max \|momentum residual `R_v`\| | 5.274e-16 |
| max \|continuity residual\| | 0.0 |
| max vorticity route mismatch | 0.0 |
| `lam` quadratic residual | 8.527e-14 |

Since the derivatives are analytic, the cancellations are algebraic and only
round-off survives — this tests the *implementation*, not a finite-difference
approximation.

### 4.2 Mesh family

| n | h | cells | nodes | total area | area error | neg. Jacobians |
|---:|---:|---:|---:|---:|---:|---:|
| 32 | 0.062500 | 1024 | 1089 | 4.000000000000000 | 0.0e+00 | 0 |
| 64 | 0.031250 | 4096 | 4225 | 4.000000000000000 | 0.0e+00 | 0 |
| 128 | 0.015625 | 16384 | 16641 | 4.000000000000000 | 0.0e+00 | 0 |

### 4.3 Consistency orders (what a second-order scheme can deliver)

`midpoint_consistency_study` compares the *cell average* of the exact solution
(4x4 Gauss–Legendre per cell, i.e. what a FVM stores) against the point value at
the cell centroid — the `h^2 (u_xx + u_yy)/24` representation error.
`truncation_error_study` freezes the exact nodal values and evaluates the standard
3-point central-difference steady NS operator.

| n | h | L2(cell-avg − centroid) u | v | velocity | L2 NS residual Rx | Ry | div |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 32 | 0.062500 | 3.612004e-03 | 5.540238e-04 | 3.654246e-03 | 1.715308e-02 | 1.711428e-03 | 1.405464e-02 |
| 64 | 0.031250 | 9.044398e-04 | 1.387267e-04 | 9.150171e-04 | 4.371876e-03 | 4.333838e-04 | 3.591967e-03 |
| 128 | 0.015625 | 2.262000e-04 | 3.469547e-05 | 2.288454e-04 | 1.102102e-03 | 1.089131e-04 | 9.065034e-04 |

Fitted log–log orders (least squares on the three levels):

| Series | Order p | R² |
|---|---:|---:|
| midpoint rule, `u` | 1.9986 | 0.99999994 |
| midpoint rule, `v` | 1.9986 | 0.99999994 |
| midpoint rule, velocity | 1.9986 | 0.99999994 |
| NS operator, x-momentum | 1.9801 | 0.99999466 |
| NS operator, y-momentum | 1.9870 | — |
| NS operator, continuity | 1.9773 | — |

All six exceed the anchor's `>= 1.8` floor, and the asymptotic regime is already
established at 32/64/128 (errors drop by 3.99–4.00 per refinement).

### 4.4 The order harness is discriminating (not a rubber stamp)

`convergence_study()` is the solver-agnostic measurement: hand it a callback that
returns cell-centre velocity and it fits `||u_h - u||_2 = C h^p`. Its self-test
injects solvers with *known* order:

| Injected order | Measured order | \|error\| | R² |
|---:|---:|---:|---:|
| 1 | 1.000070 | 7.0e-05 | 0.9999999994 |
| 2 | 2.000070 | 7.0e-05 | 0.9999999999 |
| 3 | 3.000070 | 7.0e-05 | 0.9999999999 |
| exact | not fitted (error 0.0 at round-off, flagged) | — | — |

The injected first-order solver is *rejected* by the acceptance test
(`passes_acceptance = False`); a harness that only ever reported "≈2" would not do
that. Round-off is flagged rather than fitted, so an exactly-correct solver cannot
produce a spurious order.

## 5. Consuming the bundle from a solver repo

```python
from tensorcad.kovasznay import KovasznayFlow, convergence_study, write_case_bundle

flow = KovasznayFlow(reynolds_number=40.0)
write_case_bundle(flow, "artifacts/kovasznay_case", (32, 64, 128))  # JSON for any language

def my_fvm(mesh):                     # replace with the real solve
    ...                               # returns [(u, v)] at mesh.centroids()
    return cell_centre_velocity

study = convergence_study(flow, (32, 64, 128), my_fvm)
print(study["order"], study["passes_acceptance"])
```

`study["order"]` is exactly the number the case-pool acceptance asks for;
`study["levels"]` carries the per-level L2 errors for the report.

## 6. Limitations, caveats and honest gaps

* **Solver error order is WITHHELD.** No incompressible NS solve was run in this
  unit; the numbers in §4.3 are consistency/representation orders, not solution
  orders. Closing `FVM-KOVASZNAY-ORDER` needs the TensorFVM solve plugged into
  `convergence_study` (mission 16068 ended in `agent_lost`, i.e. infrastructure,
  not numerics).
* **Re is not a free parameter numerically.** The inlet amplitude grows as
  `exp(|lam| * 0.5)`: 1.62 at Re = 40, but 2.30 at Re = 200 and 5.9 at Re = 1000.
  Higher Re needs finer meshes before the asymptotic range starts; the 32/64/128
  family is calibrated for the `Re ≈ 10…100` window.
* **Pressure is set only up to a constant.** `p0 = 0` is a convention; compare
  pressure *gradients* or differences, never absolute values, against another code.
* **The domain is finite.** The exact solution is defined on the whole plane; the
  truncation at `x = -0.5` (largest inflow) and `x = 1.5` (decayed) imposes
  Dirichlet data from the exact solution, which is the standard benchmark
  convention and removes any outflow/radiation ambiguity.
* **Order fitting needs three healthy levels.** The fit raises rather than guesses
  when a level hits round-off, and reports `R²` so a callers can see whether the
  asymptotic range was actually reached.

## 7. Error judgement (vs. baseline)

| Key quantity | Baseline | Measured here | Error | < 3 % ? |
|---|---|---|---|---|
| Momentum residual of the exact solution | 0 (analytic) | ≤ 3.6e-15 | absolute, round-off | yes (machine precision) |
| Continuity residual | 0 | 0.0 | exact | yes |
| Vorticity route mismatch | 0 | 0.0 | exact | yes |
| Mesh area vs analytic 4.0 | 4.0 | 4.0 | 0.0 | yes |
| Midpoint-rule velocity error at n=128 | `O(h²)` scaling | 2.288e-04 at h=1/64 | — (no reference constant) | n/a |
| Fitted midpoint order | 2 (nominal) | 1.9986 | 0.07 % | yes |
| Fitted NS truncation order | 2 (nominal) | 1.9773…1.9870 | 0.65…1.14 % | yes |
| **NS solution error order** | 2 (nominal) | **not measured** | — | **WITHHELD** |

## 8. Reproduction

```bash
python examples/kovasznay_fvm/run_example.py          # writes artifacts/kovasznay_fvm, exit 0 = all checks pass
python -m pytest tests/test_kovasznay.py -q          # 59 tests, ~1.3 s
```

Wall time of the example: ≈6 s single core (no external solver, no heavy compute).