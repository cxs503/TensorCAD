# Standard Hull Baseline Fleet — Verification Report

**Anchor**: `GAP-CAD-T1-BASELINE-FLEET`  
**Mission**: m21168  
**Date**: 2026-11-15  
**TensorCAD version**: 0.4.0+ (agent/m21168 branch)  
**Run directory**: `runs/m21168_fleet/`

---

## 1. Case Introduction

### 1.1 Physical Background

The baseline fleet comprises six standard hulls that form the backbone of computational ship hydrodynamics benchmarking:

| Ship | Class | Workshop / Source | Lpp [m] | B [m] | T [m] | Cb |
|------|-------|-------------------|--------:|------:|------:|---:|
| KCS | KRISO Container Ship | Gothenburg 2010 / SIMMAN 2014 | 230.0 | 32.2 | 10.8 | 0.651 |
| KVLCC2 | MOERI VLCC Tanker | SIMMAN 2008 / 2014 | 320.0 | 58.0 | 20.8 | 0.8098 |
| ONRT | ONR Tumblehome (DTMB 5613) | Gothenburg 2010 | 154.0 | 18.78 | 5.494 | 0.535 |
| JBC | Japan Bulk Carrier (SR266) | Tokyo 2015 | 280.0 | 45.0 | 16.5 | 0.858 |
| DTC | Duisburg Test Case Container Ship | El Moctar et al. 2012 | 355.0 | 51.0 | 14.5 | 0.661 |
| Series 60 | DTMB methodical series parent | Todd 1963, DTMB Report 1710 | 100.0 | 13.33 | 5.33 | 0.60 |

These hulls span the full range of commercial and naval displacement forms: fast container ships (KCS, DTC), full-form tankers (KVLCC2, JBC), a naval tumblehome (ONRT), and a methodical-series parent (Series 60). Their published particulars are the reference against which all CFD and potential-flow solvers are validated.

### 1.2 Geometry Routes

Each ship is built through one of two parametric geometry routes:

- **`series60_form`** — Todd's Series 60 section law, validated for Cb ∈ [0.55, 0.85]. Used for KCS, DTC, Series 60.
- **`parametric_sections`** — Generalized Beta section-area law solved to the record's Cb, Cm, and LCB. Used for KVLCC2, ONRT, JBC (Cb outside the Series 60 window).

Both routes produce watertight meshes by construction: waterline slices are closed rings, consecutive rings are joined by quads, and end stations are fanned shut. No repair or welding step is needed.

### 1.3 Acceptance Criteria

| Criterion | Threshold | Gate |
|-----------|-----------|------|
| Watertight (boundary edges = 0) | 0 | G3 |
| Edge-manifold (non-manifold edges = 0) | 0 | G4 |
| Consistent outward orientation | 0 inconsistent edges | G4 |
| Single closed shell (Euler χ = 2) | 2 | G5 |
| Volume vs reference | < 1% | G7 |
| Provenance complete (source + URL + license) | — | G0 |
| Not a placeholder (> 304 triangles) | > 1000 | G1 |
| No duplicate geometry (unique SHA-256) | — | G2 |

---

## 2. Method

### 2.1 Build Pipeline

```
ShipRecord (L, B, T, Cb, Cm, LCB)
    │
    ├─ series60_form ──→ Series60Parameters → section law → waterline slices → closed mesh
    │
    └─ parametric_sections ──→ Beta shape solve (a, b, k) → section-area distribution
                                → waterline slices → closed mesh
    │
    ▼
MeshBuilder → orient_outward → Mesh
    │
    ▼
Topology gates (G1–G5) + Volume gate (G7) + Provenance gate (G0)
    │
    ▼
Verdict: ADMITTED / REFUSED
```

### 2.2 Production Resolution

All ships are built at **nx = 481 stations, nz = 97 waterline levels** (production resolution). Convergence is verified on three grids: coarse (121×25), medium (241×49), and fine (481×97), using Richardson extrapolation to estimate the exact volume and observed convergence order.

### 2.3 Independent Volume Cross-Check

The mesh volume (divergence theorem over triangles) is independently verified against the hydrostatics solver (`tensorcad.hydrostatics.compute_hydrostatics`), which integrates the submerged volume by a different numerical method. Agreement better than 10⁻⁶ relative is required.

### 2.4 Artifact Round-Trip

Each admitted mesh is written to binary STL, then re-read and re-gated. The re-ingested mesh must reproduce the original topology and volume to float32 precision (relative error < 10⁻⁶).

---

## 3. Results

### 3.1 Fleet Summary

| Metric | Value |
|--------|-------|
| Ships registered | 6 |
| Ships admitted | 6 |
| Ships refused/failed | 0 |
| Negative controls all rejected | ✅ (5/5) |
| Cross-route agreement (KCS) | ✅ (< 1%) |
| Volume integration agreement | ✅ (all < 10⁻⁶) |
| Artifacts written | 6 STL files |
| Artifacts re-ingest OK | ✅ (6/6) |
| Source anchors complete | ✅ (6/6) |

### 3.2 Per-Ship Volume Agreement

| Ship | Route | V_ref [m³] | V_mesh [m³] | Error [%] | Richardson Error [%] | Order p |
|------|-------|-----------:|-----------:|----------:|---------------------:|--------:|
| KCS | series60_form | 52,030.0 | 51,850.7 | **0.3445** | 0.0958 | 1.066 |
| KVLCC2 | parametric_sections | 312,622.0 | 310,946.1 | **0.5361** | −0.0049 | 1.058 |
| ONRT | parametric_sections | 8,518.1 | 8,493.5 | **0.2887** | −0.2036 | 1.382 |
| JBC | parametric_sections | 178,369.9 | 177,417.2 | **0.5341** | 0.0009 | 1.063 |
| DTC | series60_form | 173,467.0 | 173,051.7 | **0.2394** | 0.2061 | 1.065 |
| Series 60 | series60_form | 4,266.6 | 4,248.9 | **0.4161** | −0.0017 | 1.077 |

**All six ships are within the 1% mission tolerance.** The worst error is KVLCC2 at 0.536%, well below the threshold.

### 3.3 Richardson Extrapolation

The Richardson-extrapolated volumes (estimates of the exact volume at zero grid spacing) are all within **0.21%** of the reference, confirming that the discretisation error is well-characterised and vanishing at the expected first-order rate.

### 3.4 Topology Verification

All six ships pass all topology gates:

| Ship | Triangles | Vertices | Boundary Edges | Non-manifold | Euler χ | Signed Volume [m³] |
|------|----------:|---------:|---------------:|-------------:|--------:|-------------------:|
| KCS | 185,278 | 92,641 | 0 | 0 | 2 | 51,850.7 |
| KVLCC2 | 185,278 | 92,641 | 0 | 0 | 2 | 310,946.1 |
| ONRT | 185,278 | 92,641 | 0 | 0 | 2 | 8,493.5 |
| JBC | 185,278 | 92,641 | 0 | 0 | 2 | 177,417.2 |
| DTC | 185,278 | 92,641 | 0 | 0 | 2 | 173,051.7 |
| Series 60 | 185,278 | 92,641 | 0 | 0 | 2 | 4,248.9 |

### 3.5 Negative Controls

| Control | Expected | Result |
|---------|----------|--------|
| 304-triangle placeholder box | Rejected | ✅ Rejected (placeholder_signature) |
| Byte-identical hull pair | Rejected | ✅ Rejected (G2_no_duplicate_geometry) |
| Clipped block coefficient | Rejected | ✅ Rejected (G7_volume_vs_reference) |
| Wrong volume record | Rejected | ✅ Rejected (G6_particulars_consistent) |
| Malformed manifest | Rejected | ✅ Rejected (ManifestError) |

---

## 4. Figures

### 4.1 Volume Agreement, Richardson Error, and Convergence Order

![Fleet Verification](../runs/m21168_fleet/fleet_verification.png)

*Left*: Volume agreement at production resolution. All ships below 1% (red dashed line). Green bars = literature reference, blue bars = closure reference. *Center*: Richardson-extrapolated error (estimate of exact-volume error). All below 0.21%. *Right*: Observed convergence order. All close to p = 1.0 (expected for linear chord interpolation).

### 4.2 Principal Particulars

![Fleet Particulars](../runs/m21168_fleet/fleet_particulars.png)

*Principal particulars of the six standard hulls: Lpp, beam, draft, and block coefficient (×100).*

### 4.3 Convergence Curves

![Fleet Convergence](../runs/m21168_fleet/fleet_convergence.png)

*Volume convergence on three grids for each ship. All show monotonic convergence toward the reference volume.*

---

## 5. Discussion

### 5.1 Geometry Route Selection

- **KCS, DTC, Series 60** use `series60_form` because their Cb values (0.651, 0.661, 0.60) fall within the validated window [0.55, 0.85].
- **KVLCC2, ONRT, JBC** use `parametric_sections` because KVLCC2 (Cb = 0.8098) is near the upper edge, ONRT (Cb = 0.535) is below the window, and JBC (Cb = 0.858) is above it. The Beta section-area law handles the full range Cb < Cm ≤ 1.

### 5.2 Surrogate vs Official Geometry

All six ships currently use **parametric surrogates** rather than official IGES geometry. The reasons:

| Ship | Official IGES Status | Issue |
|------|---------------------|-------|
| KCS | ✅ VERIFIED (0.30% error) | IGES field route works |
| KVLCC2 | ❌ WITHHELD | "fewer than two usable waterline slices" — IGES patches don't evaluate below DWL |
| JBC | ❌ WITHHELD | −18.28% volume error; 29/45 patches degenerate |
| ONRT | ❌ WITHHELD | −98.91% error; model-scale bytes or units mismatch |
| DTC | ❌ Not fetched | No official IGES in custody |
| Series 60 | ❌ Not digitised | Offset table from DTMB Report 1710 not digitised |

The parametric surrogates carry the published Cb, Cm, and LCB and reproduce the reference volume to < 1%, but they lack hull-specific features (bulbous bows, transom sterns, local fairing). This is an honest limitation recorded in each ship's provenance (`reference_basis: literature_value_recalled` or `closure`).

### 5.3 Convergence Behaviour

All ships show first-order convergence (p ≈ 1.0–1.4), consistent with linear chord interpolation between waterline stations. The Richardson-extrapolated volumes are all within 0.21% of the reference, confirming that the discretisation error is well-characterised.

### 5.4 Limitations

1. **No official IGES geometry** for 5 of 6 ships (KVLCC2, JBC, ONRT, DTC, Series 60). The parametric surrogates are validated against published displacement but do not capture hull-specific features.
2. **All provenance records have `verified = false`** because this environment has no outbound network to re-fetch workshop pages.
3. **Series 60 reference is closure-based** (Cb×L×B×T) rather than a published displacement, because the parent form's displacement is not independently quoted in the literature.

---

## 6. Error判定

| Ship | V_mesh vs V_ref | < 1%? | Richardson Error | < 0.3%? |
|------|----------------|-------|-----------------|---------|
| KCS | 0.3445% | ✅ Yes | 0.0958% | ✅ Yes |
| KVLCC2 | 0.5361% | ✅ Yes | −0.0049% | ✅ Yes |
| ONRT | 0.2887% | ✅ Yes | −0.2036% | ✅ Yes |
| JBC | 0.5341% | ✅ Yes | 0.0009% | ✅ Yes |
| DTC | 0.2394% | ✅ Yes | 0.2061% | ✅ Yes |
| Series 60 | 0.4161% | ✅ Yes | −0.0017% | ✅ Yes |

**All six ships pass the 1% volume tolerance and the 0.3% Richardson tolerance.**

---

## 7. Artefacts

| File | Description |
|------|-------------|
| `runs/m21168_fleet/evidence.json` | Full run evidence (gates, topology, convergence, hydrostatics) |
| `runs/m21168_fleet/ref_vs_got.csv` | Machine-readable volume comparison table |
| `runs/m21168_fleet/manifests/*.json` | Per-ship provenance manifests |
| `runs/m21168_fleet/*.stl` | Binary STL geometry files (6 × ~9 MB) |
| `runs/m21168_fleet/fleet_verification.png` | Volume agreement + Richardson error + convergence order |
| `runs/m21168_fleet/fleet_particulars.png` | Principal particulars bar chart |
| `runs/m21168_fleet/fleet_convergence.png` | Per-ship convergence curves |

---

## 8. Reproduction

```bash
# Run the full fleet gate at production resolution
python -m tensorcad.fleet.cli --out-dir runs/m21168_fleet --nx 481 --nz 97

# Run the test suite
python -m pytest tests/test_fleet_baseline.py -v
```

---

*Report generated by TensorCAD fleet run (m21168). All numbers are programmatic — no manual transcription.*
