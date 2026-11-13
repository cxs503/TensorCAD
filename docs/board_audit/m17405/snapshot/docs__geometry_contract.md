# TensorCAD Geometry Contract

**Version:** 1.0
**Date:** 2026-11-13
**Anchor:** GAP-CAD-D3-GEOM-CONTRACT

## Purpose

This document defines the four conventions that every TensorCAD mesh must honour.
A mesh that passes all four checks is *contract-compliant*; downstream consumers
(solvers, viewers, export pipelines) can rely on the declared units, frame, scale
and orientation without guessing.

## The Four Conventions

### 1. Coordinate System

**Rule:** The coordinate system must be right-handed, Cartesian.

**Rationale:** All numerical methods in the TensorCAD ecosystem (BEM, FEM, FVM,
LBM) assume a right-handed basis. A left-handed system would flip the sign of
cross products, surface normals, and volume integrals, silently corrupting every
downstream calculation.

**Check:** The `coordinate_frame` string must contain the token `right-handed`
(case-insensitive). If the frame is `unspecified`, a warning is emitted but the
mesh is not rejected (the default assumption is right-handed).

**Violation:** `left-handed` in the frame string → **error**.

### 2. Units

**Rule:** The mesh must declare its length unit as one of the recognised SI tokens.

**Allowed tokens:** `m`, `mm`, `cm`, `dm`, `km`

**Rationale:** STL files carry no unit metadata. The unit must be declared by the
caller when reading the file, and it must be a standard SI token so that scale
conversions are unambiguous.

**Check:** `Mesh.units` must be in the allowed set.

**Violation:** Any other token (e.g., `inch`, `ft`, `unknown`) → **error**.

### 3. Scale

**Rule:** The bounding-box diagonal, expressed in metres, must lie within a
plausible range for a ship-scale model.

**Range:** 0.05 m ≤ diagonal ≤ 500 m

**Rationale:** A model whose diagonal falls outside this band is almost certainly
mislabelled (e.g., a millimetre model declared as metres, or vice versa). The
lower bound of 5 cm excludes toy-scale models; the upper bound of 500 m exceeds
the length of any real ship.

**Check:** Compute the bounding-box diagonal from vertex coordinates, convert to
metres using the declared unit, and verify it lies within the range.

**Violation:** Diagonal outside the range → **error**.

### 4. Orientation

**Rule:** For a closed (watertight) surface, the signed volume must be positive
(outward-facing normals). For an open surface, the majority of face normals must
point into the same hemisphere (consistent winding).

**Rationale:** The divergence-theorem volume integral `V = (1/6) Σ (a · (b × c))`
is positive if and only if the faces are wound counter-clockwise when viewed from
outside. A negative volume means the normals point inward, which would cause
solvers to compute negative pressures, inverted forces, and other catastrophes.

**Check:**
- Closed surface: `signed_volume > 0`
- Open surface: ≥ 50% of face normals agree with the majority direction

**Violation:** Negative volume on a closed mesh → **error**.
Inconsistent winding on an open mesh → **warning**.

## Validator API

### Library

```python
from tensorcad.geom_contract import validate_contract, validate_stl_file

# Validate a Mesh object directly
report = validate_contract(mesh, source="my_model.stl")
if not report.passed:
    for v in report.violations:
        print(f"[{v.severity}] {v.convention}: {v.message}")

# Validate an STL file (units and frame must be declared)
report = validate_stl_file("hull.stl", units="m",
    coordinate_frame="x aft->fore, y port->starboard, z up; right-handed")
```

### CLI

```bash
# Human-readable output
python -m tensorcad.geom_contract hull.stl --units m --frame "x aft->fore; right-handed"

# JSON output (for pipelines)
python -m tensorcad.geom_contract hull.stl --units m --json

# Multiple files
python -m tensorcad.geom_contract *.stl --units m
```

## Compliance Table

| Asset | Status | Units | Diagonal (m) | Volume (m³) | Watertight | Violations |
|---|---|---|---|---|---|---|
| wigley.stl | PASS | m | 10.07 | 2.77 | yes | (warning: frame unspecified) |
| series60_cb60.stl | PASS | m | 20.28 | 52.01 | yes | — |
| series60_cb75.stl | PASS | m | 20.28 | 64.75 | yes | — |
| box_barge.stl | PASS | m | 31.45 | 1200.0 | yes | — |
| sphere_r1.stl | PASS | m | 3.46 | 4.17 | yes | — |
| inverted_wigley.stl | **FAIL** | m | 10.07 | -2.77 | yes | orientation: negative volume |
| wrong_units.stl | **FAIL** | mm | 0.01 | 2.77 | yes | scale: diagonal below minimum |
| left_handed.stl | **FAIL** | m | 10.07 | 2.77 | yes | coordinate_system: left-handed |

**5/8 assets passed.** The 3 failures are deliberate red-team injections proving
the validator catches real defects.

## Red-Team Evidence

The three failing assets demonstrate that the validator detects:

1. **Inverted normals** (`inverted_wigley.stl`): Signed volume = -2.77 m³.
   The mesh is geometrically identical to the valid Wigley but with all face
   windings reversed. The orientation check catches this immediately.

2. **Wrong units** (`wrong_units.stl`): The Wigley geometry (diagonal ≈ 10 m)
   is declared as millimetres, giving a diagonal of 0.01 m in "metres". The
   scale check flags this as implausible for a ship model.

3. **Left-handed frame** (`left_handed.stl`): The coordinate frame declaration
   explicitly states "left-handed". The coordinate system check rejects this
   because all TensorCAD solvers assume right-handed coordinates.

## Provenance

- **Contract module:** `tensorcad/geom_contract.py`
- **Test assets:** `artifacts/geom_contract_assets/`
- **Compliance report:** `artifacts/geom_contract_assets/compliance_report.json`
- **Asset generator:** `tools/generate_contract_assets.py`
- **Compliance runner:** `tools/run_contract_compliance.py`
