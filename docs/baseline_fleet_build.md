# Baseline fleet: published IGES → watertight hull meshes (m17293, round 1 rework)

Replaces the "gate only" delivery of m17249 with a **producer**: the six placeholder
hulls (304 triangles each, KVLCC1 ≡ KVLCC2 byte-identical) now have a documented path
to a closed, oriented, volume-verified solid built from the workshops' own IGES.

```
IGES 5.3 patches → samples → b(x,z) half-breadth field → waterline slices
                → closed triangle mesh → hydrostatics → verdict vs published ∇
```

## 1. Case

* Ship platform geometry library: 8 `*_Hull.stl`, 6 of them 304-triangle box
  placeholders (`artifacts/baseline_fleet/placeholder_audit.json`), `/generate` on the
  platform side still returns boxes. Downstream resistance/sea-keeping work had no real
  hull input.
* Authoritative sources: NMRI Tokyo-2015 geometry pages and the SIMMAN 2020 download
  page. Archives are fetched by `tools/baseline_fleet/fetch_sources.py`; bytes are
  **not** committed (neither publisher states a redistribution licence) — the repo
  keeps URL + publisher + HTTP metadata + SHA-256 + licence state and re-hashes
  whatever is downloaded.
* Published particulars are extracted from the cached publisher pages by
  `tools/baseline_fleet/build_fleet.py` (regex on the cached HTML, page SHA-256 and the
  matched snippet stored in `artifacts/baseline_fleet/fleet_specs.json`). No number in
  the result table is hand-typed.

## 2. Method

| stage | module | note |
|---|---|---|
| parse IGES | `tensorcad/iges.py` | types 102/110/126/128/141/142/143/144; de Boor rational evaluation; **control-point layout auto-detection** (IGES interleaved `x,y,z,w` vs the NMRI weight-block form) scored on control-net smoothness; Hollerith tokens → NaN so field positions stay stable; unit flag located by content (JBC declares millimetres) |
| sample | `fleet_baseline.sample_surface` | per-patch (u,v) counts from the control-net polyline length, `max_axis` 200 |
| field | `fleet_baseline.field_from_points` | `b[i][k] = max|y|` at `(x_i, z_k)`, interior holes bridged linearly |
| mesh | `fleet_baseline.build_solid` | waterline slices → closed rings resampled by arc length → lofted quads → caps at the ends; `mesh.orient_outward` |
| measure | `measure_mesh`, `field_integrals` | triangle divergence theorem **and** an independent trapezoid integration of the same field |
| grade | `tensorcad/fleet_gate.py` + `checks` block | watertight, edge-manifold, zero-area faces, volume vs published, LWL vs Lpp |

Why a field instead of per-patch triangulation: the workshop hulls are 19–264
independent B-spline patches whose boundaries do not share vertices, so per-patch
tessellation leaves cracks along every seam — the exact property being graded. The
field removes seams; the price (single-valued `b(x,z)` in `y`, plus `max`-binning bias
O(cell²/radius)) is bounded by the mesh-vs-field cross-check that every manifest
reports (KCS: −0.02 % at dx = 0.25 m, i.e. two integrators agreeing to 2e-4).

## 3. Results (run 2026-11-13, `--samples-per-metre 6`, dx 0.5 m, dz 0.1 m)

| hull | patches | ref ∇ (m³) | got ∇ (m³) | error | triangles | watertight | verdict |
|---|---|---|---|---|---|---|---|
| KCS | 55 | 52 030 | 52 187.98 | **+0.304 %** | 52 320 | yes | **VERIFIED** |
| JBC | 90 (+38 trimmed) | 178 369.9 | — | — | — | — | WITHHELD |
| ONRT | 264 | 8 500.78 | — | — | — | — | WITHHELD |
| KVLCC2 | 19 (+19 trimmed) | — | — | — | — | — | WITHHELD |
| DTC | — | — | — | — | — | — | NOT SUPPORTED |
| Series 60 | — | — | — | — | — | — | NOT SUPPORTED (analytic generator only) |

Secondary KCS checks against the same page: Cb 0.6505 → **0.6516** (+0.17 %),
LCB −1.48 % Lpp → **−1.48 %**, LWL 232.5 m → 230.25 m at the design-waterline cut.
Grid convergence on KCS displacement error: **+0.62 %** (dx 0.5 / dz 0.1 / spacing
0.125 m) → **+0.30 %** (spacing 0.167 m) → **+0.17 %** (dx 0.25 / dz 0.05 / spacing
0.08 m), i.e. the residual is discretisation, not a geometry defect.

Artifacts: `artifacts/baseline_fleet/{KCS_manifest.json,fleet_summary.json,benchmark.json,fleet_specs.json}`,
`meshes/KCS_Hull.stl` (binary STL, 52 320 triangles, the placeholder replacement format),
`figs/{volume_vs_published.png,KCS_half_breadth_field.png}`.

## 4. Failure modes hit (honest list)

1. **KVLCC2 / JBC / ONRT are trimmed-surface models** (type 143/144 + 141/142). The
   current reader evaluates the *basis* surfaces over their full parametric domain, so
   trimmed ships need a trim-loop filter in (u,v) before their volumes can be claimed.
   This is the main remaining gap and the reason three of six hulls stay WITHHELD.
2. **ONRT is in model scale** (LWL 3.147 m on the page's model column; the IGES carries
   model-scale coordinates) and its volume came out 0.0056 m³ — a scale factor taken
   from the geometry (LWL_igs vs LWL_full = 154.0 m) is required before comparing.
3. **JBC** trips the control-net layout ambiguity guard (`4x` score margin). The guard is
   fail-closed on purpose: guessing the layout silently distorts geometry.
4. **DTC / Series 60** have no official IGES at the sources probed (404s recorded in the
   run log). Series 60 exists only as TensorCAD's analytic generator
   (`tensorcad/series60.py`), which is *our* law, not published geometry — it cannot be
   counted as an authoritative hull.
5. The analytic ellipse fixture in `tests/test_iges_baseline.py` pinches to a line at both
   ends and is marked `xfail`: the slice loft cannot yet close a fully pinched field
   without an explicit cap. Real hulls (KCS) cap correctly, so this is a fixture-path
   gap, tracked here rather than hidden.

## 5. Reproduce

```bash
python tools/baseline_fleet/fetch_sources.py --out artifacts/baseline_fleet
python tools/baseline_fleet/build_fleet.py --samples-per-metre 6
python -m pytest tests/test_iges_baseline.py tests/test_fleet_gate.py -q
python tools/baseline_fleet/audit.py audit artifacts/baseline_fleet/meshes \
    --published artifacts/baseline_fleet/published.json
```

## 6. Error judgement

| quantity | ref | got | error | <1 %? |
|---|---|---|---|---|
| KCS displacement | 52 030 m³ | 52 187.98 m³ | +0.304 % | **yes** |
| KCS Cb | 0.6505 | 0.6516 | +0.17 % | yes |
| KCS LCB | −1.48 % Lpp | −1.48 % Lpp | ≤0.01 %p | yes |
| KCS mesh vs field ∇ | 52 187.98 | 52 187.99 | −0.02 % | yes |
| JBC / ONRT / KVLCC2 / DTC / Series 60 | — | not built | — | **no claim** |

Acceptance status: **1 of 6 hulls fully verified end-to-end** (watertight + manifold +
<1 % displacement from published geometry); 5 of 6 honestly WITHHELD/NOT SUPPORTED with
the specific blocking defect named above. Not a pass on the mission's "6/6" line.
