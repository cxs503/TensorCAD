# Placeholder STL Replacement Plan — GAP-CAD-T1-BASELINE-FLEET

**Repository** TensorCAD · **Anchor** GAP-CAD-T1-BASELINE-FLEET
**Date** 2026-11-13 · **Author** agent/m18255

## 1. Problem Statement

The ship-platform geometry library shipped 8 STL files of which 6 were 304-triangle
placeholder boxes, and two of the "different" ships (KVLCC1/KVLCC2) were byte-identical.
No build-time check detected this because "a file exists" was the only predicate evaluated.

This plan defines the **replacement pipeline** that turns official workshop geometry
(IGES / offset tables) into the watertight, manifold, volume-verified solids the fleet
gate requires, and the **intake gate** that prevents placeholders from ever reaching a
solver again.

## 2. Replacement Pipeline

### 2.1 Source Acquisition

| Ship | Official Source | Format | URL |
|---|---|---|---|
| KCS | NMRI / SIMMAN 2014 | IGES (B-spline patches) | https://simman2014.dk/ |
| KVLCC2 | SIMMAN 2008 / NMRI | IGES (B-spline patches) | https://simman2014.dk/ |
| ONRT | DTMB 5613 / Gothenburg 2010 | IGES / offset table | https://simman2014.dk/ |
| JBC | Tokyo 2015 Workshop | IGES / offset table | https://tokyo2015.k.u-tokyo.ac.jp/ |
| DTC | HSVA / el Moctar et al. 2012 | IGES / offset table | https://www.uni-due.de/ship/dtc.php |
| Series 60 | DTMB Report 1710 (Todd 1963) | Offset table (published) | https://ntrl.ntis.gov/ |

### 2.2 Geometry Routes (implemented in `tensorcad.fleet.build`)

Five routes, one contract: each returns a `Mesh` or raises `RouteRefusal` with a
measured reason.

| Route | Input | How it guarantees watertightness | When to use |
|---|---|---|---|
| `series60_form` | Cb, L, B, T | Structured waterline slices from analytic section law | Series 60 hulls, Cb ∈ [0.55, 0.85] |
| `parametric_sections` | Cb, Cm, LCB, L, B, T | Beta section-area law solved to published particulars | Any displacement hull with published Cm/LCB |
| `offset_table` | stations × waterline half-breadths | Loft by construction: shared keel/deck vertices, fan closure | Official offset tables from workshops |
| `imported_stl` | binary/ASCII STL file | Weld duplicate vertices, orient outward, verify | Vendor STL files, legacy meshes |
| `iges_field` | IGES file (B-spline patches) | Half-breadth field → structured waterline slices | Official IGES from NMRI/SIMMAN/Gothenburg |

The **IGES field route** is the primary path for KCS, KVLCC2, ONRT, JBC and DTC:
1. Read IGES B-spline surfaces (`tensorcad.iges.read_iges`)
2. Sample each patch on an adaptive grid (`sample_patches`)
3. Collapse to half-breadth field `b(x,z) = max|y|` (`field_from_points`)
4. Build watertight solid from structured waterline rings (`build_solid`)

This guarantees watertightness **by construction** because topology is built from the
field grid, not from patch seams. The cost is explicit: `max` binning biases the surface
outwards by O(cell²/radius), bounded by the reported `field_cover` metric.

### 2.3 Intake Gate (implemented in `tensorcad.fleet.gates`)

Nine machine-decidable gates, evaluated in order. A mesh must pass all to be admitted.

| Gate | Name | What it checks | Tolerance |
|---|---|---|---|
| G0 | provenance | source URL, license, sha256, retrieval stated | complete |
| G1 | not_placeholder | triangle count > 304, non-zero volume, not a box | min 1000 triangles |
| G2 | no_duplicate_geometry | sha256 unique across fleet | exact match |
| G3 | watertight | zero boundary edges | 0 |
| G4 | manifold_oriented | no non-manifold edges, consistent orientation, no zero-area/duplicate faces | 0 |
| G5 | single_closed_shell | Euler characteristic = 2, positive signed volume | exact |
| G6 | particulars_consistent | Cb·L·B·T agrees with declared displacement | 3 % |
| G7 | volume_vs_reference | mesh volume vs reference volume | 1 % |
| G8 | principal_dimensions | mesh L/B/T span declared values | 0.1 % |

### 2.4 Negative Controls (implemented in `tensorcad.fleet.controls`)

Four controls that **must be rejected**, proving the gate is not a rubber stamp:

| Control | Defect it reproduces | Expected verdict |
|---|---|---|
| `placeholder_304_triangles` | 304-triangle box (the original defect) | Rejected by G1 |
| `byte_identical_pair` | Two ships sharing one STL file (KVLCC1=KVLCC2) | Rejected by G2 |
| `clipped_cb` | Cb raised into form-family window but volume unchanged | Rejected by G7 |
| `incoherent_particulars` | Declared displacement contradicts Cb·L·B·T | Rejected by G6 |

All four are rejected in the current test suite (36 tests, all passing).

## 3. Replacement Procedure

### Step 1: Acquire official geometry
Download IGES or offset table from the official source. Record:
- Source URL
- Publisher
- License
- SHA-256 of the downloaded file
- Retrieval timestamp

### Step 2: Update the manifest
Add or update the ship record in `tensorcad/fleet/data/fleet_manifest.json`:
```json
{
  "name": "KCS",
  "geometry_route": "iges_field",
  "geometry_path": "data/iges/KCS_hull.iges",
  "geometry_sha256": "<sha256 of IGES file>",
  "provenance": {
    "source": "NMRI T2015 workshop",
    "url": "https://www.t2015.nmri.go.jp/",
    "license": "public benchmark data",
    "retrieval": "downloaded",
    "verified": true
  }
}
```

### Step 3: Run the fleet gate
```bash
python -m tensorcad.fleet.cli --out-dir runs/fleet_verification --nx 481 --nz 97
```

### Step 4: Verify the verdict
Check `runs/fleet_verification/evidence.json`:
- `ships[].passed` must be `true`
- `ships[].topology.watertight` must be `true`
- `ships[].volume_relative_error` must be `< 0.01`
- `negative_controls[].rejected` must all be `true`

### Step 5: Replace the placeholder STL
Once a ship passes all gates, replace the placeholder STL in the ship-platform geometry
library with the verified mesh:
```bash
python -c "
from tensorcad.fleet import load_manifest, build_ship_geometry, write_binary_stl
manifest = load_manifest()
ship = manifest.by_name('KCS')
built = build_ship_geometry(ship, nx=481, nz=97)
write_binary_stl(built.mesh, 'ship_performance_platform/data/KCS_Hull.stl')
"
```

## 4. Current Status

| Ship | Route | V_ref [m³] | V_mesh [m³] | Error | Watertight | Status |
|---|---|---|---|---|---|---|
| KCS | series60_form | 52030.0 | 51850.7 | 0.34 % | ✅ | Admitted (surrogate) |
| KVLCC2 | parametric_sections | 312622.0 | 310946.1 | 0.54 % | ✅ | Admitted (surrogate) |
| ONRT | parametric_sections | 8518.1 | 8493.5 | 0.29 % | ✅ | Admitted (surrogate) |
| JBC | parametric_sections | 178378.2 | 177417.2 | 0.54 % | ✅ | Admitted (surrogate) |
| DTC | series60_form | 155834.4 | 155407.2 | 0.27 % | ✅ | Admitted (surrogate) |
| Series60 | series60_form | 4266.6 | 4248.9 | 0.42 % | ✅ | Admitted (surrogate) |

All six ships pass all gates at production resolution (481×97). The worst volume error
is 0.54 % (JBC), well within the 1 % mission tolerance.

**What is WITHHELD**: Official IGES/offset files could not be downloaded in this
environment (no outbound network). The current meshes are surrogates built from published
particulars (Cb, Cm, LCB). The gate is proven on surrogates; admitting a real file only
needs `geometry_route: "iges_field"` or `"offset_table"` plus its sha256 in the manifest.

## 5. What This Plan Prevents

| Defect | How it is prevented |
|---|---|
| 304-triangle placeholder | G1: triangle count > 304, non-zero volume, not a box |
| Byte-identical files for different ships | G2: sha256 collision detection across fleet |
| Unclosed surface reaching a solver | G3: zero boundary edges required |
| Non-manifold mesh | G4: no non-manifold edges, consistent orientation |
| Wrong ship geometry | G7: volume vs reference < 1 % |
| Incoherent particulars | G6: Cb·L·B·T must agree with declared displacement |
| Untraceable geometry | G0: provenance (source, URL, license, sha256) required |

## 6. Reproduction

```bash
# Run the full fleet gate with negative controls
python -m tensorcad.fleet.cli --out-dir runs/fleet_verification --nx 481 --nz 97

# Run the test suite
python -m pytest tests/test_fleet_baseline.py -q

# Verify a single ship
python -c "
from tensorcad.fleet import load_manifest, build_ship_geometry, measure
manifest = load_manifest()
ship = manifest.by_name('Series60')
built = build_ship_geometry(ship, nx=481, nz=97)
facts = measure(built.mesh)
print(f'Volume: {facts.signed_volume:.2f} m³')
print(f'Watertight: {facts.boundary_edge_count == 0}')
print(f'Euler: {facts.euler_characteristic}')
"
```

## 7. Future Work

1. **Official geometry ingestion**: When outbound network is available, download IGES
   files from NMRI/SIMMAN/Gothenburg and run them through the `iges_field` route.
2. **Appendage handling**: The current fleet is bare hull only. Appendages (rudder,
   propeller, bilge keels) should be declared as separate geometry files with their
   own provenance.
3. **Free-surface mesh**: The current mesh is capped at the waterline. A full
   free-surface mesh would extend above the waterline and include the air domain.
4. **Class-society compliance**: Add gates for class-society rules (e.g., minimum
   bow height, freeboard, intact stability criteria).
