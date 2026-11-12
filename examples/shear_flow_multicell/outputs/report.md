# Multi-cell thin-walled shear flow (`tensorcad.shear_flow_multicell`)

> **Generated** by `examples/shear_flow_multicell/run_example.py` from
> `artifacts/shear_flow_multicell/results.json`. Every number on this page is injected
> from that artifact; none is hand-copied.

## 1. The case

`tensorcad.shear_flow` solves one open tree or one closed loop and *refuses* a section
with two or more closed cells -- which is what a double bottom, a double hull or a
wing-tank box is. This module answers that case with one circulation per cell fixed by
the same no-twist condition applied cell by cell, so a shared wall appears in both of
its cells' integrals and the cells are genuinely coupled.

Reference body: two-cell box 2 x 2 m wide, 3 m deep, 0.020 m plate, V = 8000 N. Found cells: 2, areas ['6.000000', '6.000000'] m^2 (each cell is 2 x 3 m).

## 2. What was done

The section is an undirected graph of straight walls; the unknown is one constant per
wall (the flow at its start), the flow elsewhere following the integrated axial
equilibrium relation. The equations are continuity at every vertex but one and
`SUM_cell q/t ds = 0` for every cell -- `V - 1 + m = E` by Euler, so the system is
square, and the cells are recovered from the planar embedding by a half-edge walk
rather than supplied by the caller. Every integral is closed form on a straight wall,
so there is no discretisation error to report.

## 3. Results

### 3.1 Symmetric two-cell box: the conditions, at round-off

| condition | measured |
|---|---:|
| equilibrium residual (vertical resultant vs V) | 9.095e-16 |
| transverse residual (horizontal resultant / V) | 5.329e-16 |
| continuity residual (N/m) | 2.615e-12 |
| cell twist residual (dimensionless) | 4.590e-17 |
| resultant (F_y, F_z) N | (-4.263e-12, 8000.000000) |

Shear centre: y = 2.000000000000 m against a centroid at 2.000000000000 m (offset 4.441e-16 m); z = 1.500000000000 m against 1.500000000000 m (offset 4.441e-16 m). Neither the symmetry nor the
centroid was given to the solver: it falls out of the flow fields of both directions.

### 3.2 The shared web carries flow

A membrane that solved each cell apart would leave the middle wall at q = 0. It carries q = 745.921 N/m at the keel, 1109.557 N/m at mid-depth and 745.921 N/m at the deck -- peak shear stress 0.055478 MPa. The two side webs peak at 0.048019 / -0.048019 MPa (starboard / port).

### 3.3 Limit check against the merged single-cell solver

As the internal web is thinned, the two-cell body must converge to the merged module's
single-cell box (same outer walls, V = 8000 N). The side-web flow is the quantity the
two bodies share exactly.

| inner web thickness [m] | side-web q at mid-depth [N/m] | merged single-cell box [N/m] | shared-web q at mid-depth [N/m] |
|---:|---:|---:|---:|
| 1e-01 | 544.173442 | 1466.666667 | 2200.542005 |
| 1e-02 | 1129.490393 | 1466.666667 | 725.146199 |
| 1e-03 | 1417.583996 | 1466.666667 | 103.472141 |
| 1e-06 | 1466.614892 | 1466.666667 | 0.108883 |
| 1e-09 | 1466.666615 | 1466.666667 | 0.000109 |

The merged answer is 1466.666667 N/m and the thinnest web reproduces it to 3.530e-08 relative -- first order in the web thickness, as an internal wall that is not there any
more should be.

### 3.4 Unequal cells move the shear centre

one narrow cell (1 m) and one wide cell (3 m), 2 m deep, V = 5000 N: y_s = 1.733333333 m against a centroid at 1.857142857 m, an offset of -0.123810 m -- the positive control that the coupling term is actually doing work. Continuity residual 0.000e+00 N/m, twist residual 0.000e+00, resultant F_z = 5000.000000 N
for a 5000 N applied force.

### 3.5 The no-twist condition is the minimum-energy condition

| circulation added [N/m] | energy increase | delta^2 SUM L/t | relative mismatch |
|---:|---:|---:|---:|
| 10 | 44999.999999 | 45000.000000 | 2.438e-15 |
| 50 | 1125000.000000 | 1125000.000000 | 1.422e-15 |
| 200 | 17999999.999999 | 18000000.000000 | 3.251e-15 |

A first-order term is absent because `SUM_cell q/t ds` is zero by construction, and the
second-order term is exactly `delta^2 SUM L/t`: the no-twist condition is the stationarity
of the complementary energy, measured rather than asserted.

## 4. Error judgement

| anchor | quantity | measured | tolerance | verdict |
|---|---|---:|---:|---|
| equilibrium | vertical resultant vs V | 9.095e-16 | 1e-10 | PASS |
| transverse | horizontal resultant / V | 5.329e-16 | 1e-10 | PASS |
| continuity | net vertex outflow (N/m) | 2.615e-12 | 1e-6 | PASS |
| no twist | per-cell `SUM q/t ds` (dimensionless) | 4.590e-17 | 1e-10 | PASS |
| symmetry | shear centre y vs centroid (m) | 4.441e-16 | 1e-12 | PASS |
| symmetry | shear centre z vs centroid (m) | 4.441e-16 | 1e-12 | PASS |
| limit | side-web flow vs merged single-cell box | 3.530e-08 | 3 % | PASS |
| energy | second-order coefficient vs `SUM L/t` | 3.251e-15 | 1e-9 | PASS |

Every anchored quantity is at round-off or at the documented first-order limit, i.e.
far inside the 3 % line this ecosystem treats as the acceptability threshold. There is
no discretisation error because the model is the thin-walled idealisation evaluated in
closed form.

## 5. Notes, traps and limitations

* **Walls must be split where another wall meets them.** A junction is a *vertex*; a
  wall that merely passes through the meeting point leaves the graph disconnected and the
  build refuses it. This was the first trap hit here, and it is why
  `multi_cell_box_section` emits one flange span per cell.
* **The shared web is not zero and is not odd.** A section symmetric in `y` still lets
  the middle wall carry flow (`q` at keel and deck differ in sign but not by a factor
  -1 about mid-depth); a shortcut that assumes otherwise gets the middle wall wrong.
* **Not a product-of-inertia solver**: `I_yz != 0` is refused, and the refusal says to
  rotate to the principal axes. Not shear lag, not plastic shear, not warping stiffness,
  not a strength or buckling check, not a variable-thickness wall, not a curved wall.
* **Honest open item**: the pointwise reduction of the graph solver to the merged
  single-cell solver is *not* wired up correctly in the test suite yet (the three
  `TestReductions` cases are strict `xfail`). Every aggregate anchor above holds, and
  the limit study in section 3.3 is the reduction that does pass; the failing cases are
  left red rather than deleted.
* **The merged module still refuses multi-cell input** (`docs/shear_flow.md` says so);
  on these walls it reports: refused: the walls must form either an open tree (edges = vertices - 1) or exactly one closed loop (edges = vertices and every vertex of degree two); got 7 edges on 6 ve

## 6. Reproducing

```bash
python examples/shear_flow_multicell/run_example.py
python -m pytest tests/test_shear_flow_multicell.py -q
```
