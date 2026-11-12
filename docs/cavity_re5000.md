# Lid-Driven Cavity Re=5000 Benchmark

## Overview

The lid-driven cavity is the canonical benchmark for incompressible Navier–Stokes solvers. This document describes the TensorCAD implementation, which provides:

1. **Structured mesh generation** for the unit square domain
2. **Reference data** from Ghia et al. (1982) for validation
3. **Boundary condition metadata** for solver consumption

## Physical Problem

### Geometry
- **Domain**: Unit square [0, 1] × [0, 1] (non-dimensional)
- **Lid**: Top wall (y = 1) moves with velocity (u, v) = (1, 0)
- **Walls**: Bottom, left, and right walls are stationary (no-slip)

### Parameters
- **Reynolds number**: Re = 5000
- **Characteristic velocity**: U = 1 (lid velocity)
- **Characteristic length**: L = 1 (domain size)
- **Kinematic viscosity**: ν = U·L/Re = 2 × 10⁻⁴

### Flow Features at Re=5000
- Steady, laminar flow
- Primary vortex centered near (0.5, 0.5)
- Secondary vortices in bottom-left and bottom-right corners
- Thin boundary layers along all walls
- Well-documented centreline velocity profiles

## Reference Data

The reference data is from:

> Ghia, U., Ghia, K. N., & Shin, C. T. (1982). "High-Re Solutions for Incompressible Flow Using the Navier-Stokes Equations and a Multigrid Method." *Journal of Computational Physics*, 48, 387–411.

### Centreline Velocities

Two profiles are provided:

1. **u-velocity along x = 0.5** (vertical centreline)
   - u(y) from y = 0 (bottom wall, u = 0) to y = 1 (lid, u = 1)
   - Captures the primary vortex structure

2. **v-velocity along y = 0.5** (horizontal centreline)
   - v(x) from x = 0 (left wall, v = 0) to x = 1 (right wall, v = 0)
   - Captures the horizontal flow structure

The data is digitised from Table 1 of Ghia et al. (1982) and stored in:
- `tensorcad/cavity_mesh.py`: `GHIA_RE5000_U_CENTERLINE`, `GHIA_RE5000_V_CENTERLINE`
- `results/cavity_re5000/ghia_re5000_reference.json`

## Mesh Generation

### Structured Quadrilateral Mesh

The mesh is a conforming structured grid on [0, 1]²:

```python
from tensorcad.cavity_mesh import generate_cavity_mesh, CavityParameters

# Uniform mesh
mesh = generate_cavity_mesh(nx=64, ny=64)

# Graded mesh (finer near walls)
mesh = generate_cavity_mesh(nx=128, ny=128, wall_grading=3.0)
```

### Mesh Properties

| Property | Description |
|----------|-------------|
| `nx`, `ny` | Number of cells in x and y directions |
| `n_nodes` | Total nodes: (nx+1) × (ny+1) |
| `n_elements` | Total quads: nx × ny |
| `nodes` | List of (x, y) coordinates |
| `elements` | Quad connectivity (counter-clockwise) |
| `boundary` | Boundary condition metadata |

### Wall Grading

The `wall_grading` parameter controls cell size distribution:

- `wall_grading = 1.0`: Uniform spacing
- `wall_grading > 1.0`: Smaller cells near walls, larger in the centre
- Typical values: 2.0–5.0 for boundary-layer resolution

The grading is symmetric: cells are refined towards all four walls.

## Usage

### Generate Mesh and Reference Data

```bash
python examples/cavity_re5000/run_example.py --nx 128 --output results/cavity_re5000
```

This produces:
- `cavity_re5000_mesh.json`: Mesh data (nodes, elements, BCs)
- `ghia_re5000_reference.json`: Reference velocity profiles
- `cavity_mesh_summary.txt`: Human-readable summary
- `cavity_mesh_quality.png`: Visualisation

### Access Reference Data in Python

```python
from tensorcad.cavity_mesh import ghia_u_at_y, ghia_v_at_x

# Interpolate u-velocity at x=0.5 for a given y
u = ghia_u_at_y(0.5)  # u ≈ 0.01944 at centre

# Interpolate v-velocity at y=0.5 for a given x
v = ghia_v_at_x(0.5)  # v ≈ 0.00126 at centre
```

### Export Mesh for Solvers

```python
from tensorcad.cavity_mesh import generate_cavity_mesh, write_mesh_json

mesh = generate_cavity_mesh(nx=128, ny=128, wall_grading=2.0)
write_mesh_json(mesh, "cavity_mesh.json")
```

The JSON format includes:
- Node coordinates
- Element connectivity
- Reynolds number and viscosity
- Boundary condition assignments

## Validation

### Mesh Quality Checks

The example script verifies:
1. **Node count**: (nx+1) × (ny+1)
2. **Element count**: nx × ny
3. **Connectivity bounds**: All indices in [0, n_nodes)
4. **Domain bounds**: All nodes in [0, L]²
5. **Reference data integrity**: Endpoints, monotonicity

### Solver Validation

To validate a solver against the Ghia reference:

1. Run the solver on the generated mesh
2. Extract centreline velocities:
   - u along x = 0.5
   - v along y = 0.5
3. Compare with `ghia_u_at_y()` and `ghia_v_at_x()`
4. Acceptance criterion: <8% error (per feature ledger)

Example validation script:

```python
import json
from tensorcad.cavity_mesh import ghia_u_at_y

# Load solver results
with open("solver_results.json") as f:
    results = json.load(f)

# Compare u-profile
max_error = 0.0
for y_solver, u_solver in results["u_centerline"]:
    u_ref = ghia_u_at_y(y_solver)
    error = abs(u_solver - u_ref) / max(abs(u_ref), 1e-8)
    max_error = max(max_error, error)

print(f"Max error vs Ghia: {max_error:.2%}")
assert max_error < 0.08, "Validation failed"
```

## Integration with Tensor Ecosystem

### TensorLBM
The structured mesh can be converted to a lattice grid for LBM solvers. The uniform spacing simplifies the mapping.

### TensorFVM
The quad elements are directly usable as finite volumes. The boundary metadata specifies which faces have which conditions.

### TensorFEM
The mesh can be converted to triangular elements (each quad → 2 triangles) for FEM solvers.

## Limitations

1. **2D only**: This module provides 2D meshes. For 3D cavity flow, extend to a box domain.
2. **Steady Re=5000**: The reference data is for steady flow. At higher Re (>10000), the flow becomes unsteady.
3. **No turbulence**: Re=5000 is laminar. For turbulent cavity flow (Re>10000), different reference data is needed.
4. **Digitised data**: The Ghia reference values are digitised from the published table, not re-computed. Small discrepancies (<0.1%) may exist vs. the original.

## References

1. Ghia, U., Ghia, K. N., & Shin, C. T. (1982). "High-Re Solutions for Incompressible Flow Using the Navier-Stokes Equations and a Multigrid Method." *JCP*, 48, 387–411.
2. Botella, O., & Peyret, R. (1998). "Benchmark Spectral Results on the Lid-Driven Cavity Flow." *Computers & Fluids*, 27(4), 421–433.
3. Bruneau, C.-H., & Saad, M. (2006). "The 2D Lid-Driven Cavity Problem Revisited." *Computers & Fluids*, 35(3), 326–348.

## Feature Ledger

- **Slug**: `CASE-CAVITY-RE5000`
- **Title**: 算例池·顶盖方腔 Re=5000
- **Acceptance**: <8% error vs Ghia et al. 1982 reference
- **Status**: Implemented (mesh + reference data)
- **Next**: Solver validation (TensorLBM/TensorFVM)
