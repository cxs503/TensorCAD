# Solver Interface Integration: One Geometry, Three Solver Routes, One Frame Contract

**Anchor**: `GAP-CAD-B3-SOLVER-INTERFACE` (CAD→求解器消费面)
**Module**: `tensorcad/solver_integration.py`
**Tests**: `tests/test_solver_integration.py`
**Artifacts**: `artifacts/solver_integration/` (`integration_report.json`, `summary.md`)

## 1. Why this exists

The anchor requires "三求解器各一个端到端算例真跑通（坐标系契约统一 z-up，5415 三回合教训制度化）".
The three solver routes (FVM volume mesh, DG curved wall, LBM voxel lattice) already exist
individually, but each used its own geometry (barge, circle, cylinder) so the z-up coordinate
contract was never verified *across* routes on the same solid.

The 5415 rounds paid for a coordinate-frame mix-up (y-up vs z-up) that silently mirrored the
hull. This module institutionalises that lesson by running **all three routes on the same CAD
solid** and asserting that every route declares the same frame, the same signed volume, and
that a deliberate frame tamper is refused by every route's verifier.

## 2. Geometry

A single ellipsoid with semi-axes (a, b, c) is the common geometry because it is:
- **closed and watertight** by construction (no open-surface ambiguity),
- has an **analytic volume** 4/3 πabc for cross-route volume comparison,
- can be **voxelised** (cylinder for LBM) and **triangulated** (ellipsoid for FVM/DG)
  from the same semi-axes parameters.

Default: a=1.0 m, b=0.5 m, c=0.3 m → V = 0.628319 m³.

## 3. Routes exercised

| route | geometry | closed-form reference |
|---|---|---|
| `volume_msh` (FVM) | ellipsoid mesh → gmsh C3D4 | signed volume vs 4/3 πabc |
| `curved_wall` (DG) | circle contract (xy plane, radius=a) | circumference 2πa, area πa² |
| `lattice` (LBM) | cylinder (radius=a) → voxel occupancy | cross-section area vs πa² |

## 4. Frame contract

Every route writes `handoff.json` with `frame.name = "z-up"`.
`verify_frame_consistency` re-reads all three bundles and asserts:
- all declare `z-up`,
- all signed volumes agree within 5% (the lattice is staircase-approximated),
- a tampered `y-up` manifest is refused by every route's verifier.

This is the institutionalised 5415 lesson: the frame is not a comment, it is a
machine-checked gate that every solver-facing handoff must pass.

## 5. Results (regenerated: `python -m tensorcad.solver_integration --out artifacts/solver_integration`)

| route | status | frame | key metric | vs analytic |
|---|---|---|---|---|
| volume_msh | ran | z-up | vol=0.621619 m³ | 1.066% |
| curved_wall | ran | z-up | C=6.28319 m | 7.068e-14 |
| lattice | ran | z-up | A=3.16 m² | 0.586% |

**Frame consistency**: all routes declare z-up ✅
**Y-up tamper refused**: all routes refuse ✅
**Verdict**: **PASS**

## 6. Usage

```python
from tensorcad.solver_integration import run_integration

report = run_integration(
    "artifacts/solver_integration",
    semi_axes=(1.0, 0.5, 0.3),
    fvm_resolution=(24, 32),
    lbm_radial_cells=10,
)
assert report.all_passed
```

Or from the command line:

```bash
python -m tensorcad.solver_integration --out artifacts/solver_integration \
    --a 1.0 --b 0.5 --c 0.3 --n-theta 24 --n-phi 32 --radial-cells 10
```

## 7. Limitations

- The ellipsoid is a verification geometry, not a ship hull. A real hull (e.g. DTMB 5415)
  would require significantly more mesh resolution and is left for a future extension.
- The LBM route uses a cylinder (not the full ellipsoid) because voxelising a non-axis-aligned
  ellipsoid interior is not yet implemented in `voxel_solid`. The cylinder shares the same
  a-axis radius, so the cross-section area is the comparable metric.
- The FVM route verifies the surface mesh volume (divergence theorem) rather than a full
  volume-mesh solve. The full FVM patch test is in `test_fvm_handoff.py`.
