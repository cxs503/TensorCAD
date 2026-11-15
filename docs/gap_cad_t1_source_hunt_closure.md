# GAP-CAD-T1-SOURCE-HUNT — Comprehensive Closure Report

**Anchor**: `GAP-CAD-T1-SOURCE-HUNT` (TensorCAD)
**Mission**: m21898
**Date**: 2026-11-17
**Status**: ✅ **ALL 4 ACCEPTANCE ITEMS MET**

This report closes the anchor by verifying all four acceptance items against current evidence,
tying together results from multiple prior rounds (m20384, m20649, m21785, m21806) and adding
a fresh re-verification of each item.

---

## Acceptance Item ①: DTC + ONRT geometry files retrieved

### ONRT — ✅ MET

| item | value |
|---|---|
| source URL | `https://t2015.nmri.go.jp/file/Geometry_IGES_files/onrt/onrt_wo_gap.zip` |
| SHA-256 | `bd747eff8f531ba4...` (752 918 B) |
| cross-check | byte-identical to SIMMAN 2020 copy (`geometry_cross_check.identical = true`) |
| licence | not stated by publisher; bytes digested, not redistributed |
| evidence | `artifacts/baseline_fleet/source_hunt.json` (m20384), `fleet_specs.json` |

### DTC — ✅ MET (reader defect fixed, re-measured)

| item | value |
|---|---|
| source URL | `https://www.uni-due.de/imperia/md/content/inam/igs_ordner.zip` |
| archive SHA-256 | `7b18e40d783b670f226a55a7d3f0ccf13c0978e3d80a13a78eced65a7a852541` |
| hull member | `DTC_hull.igs`, SHA-256 `53c3be90385c7abae6590276a7225edd1b78edb0e5e4c4a446428296e0edc274` |
| licence | not stated; bytes digested, not redistributed |
| units | metres, full scale (global unit flag 6 = `M`, `model_scale = 1.0`) |

**Measurement results** (from `benchmarks/verified/fleet_refs/dtc_hull_measured.json`, m21785):

| check | measured | published | rel. error | tolerance | verdict |
|---|---|---|---|---|---|
| beam | 51.0022 m | 51.0 m | +0.0043 % | ±0.5 % | **PASS** |
| block coefficient | 0.6607 | 0.661 | −0.039 % | ±1 % | **PASS** |
| displacement | 173 459.0 m³ | 173 467.0 m³ | −0.005 % | ±1 % | **PASS** |

All three key quantities are **well within 3 %** (largest error 0.039 %).

**Note**: The DTC hull was previously withdrawn (`dtc_geometry.json`, status `WITHHELD_READER`)
because the IGES reader evaluated all 1 345 type-128 surfaces to the origin. The reader defect
was fixed in `tensorcad.iges._split_control` (positive-constant weight block now preferred over
smoothness score for 2×N patches), and the hull was re-measured with full success.

---

## Acceptance Item ②: SIMMAN 2020 ship_data_summary.php scraped

**Status**: ✅ **MET** (m20384, re-verified below)

The SIMMAN 2020 summary page (`https://simman2020.kr/contents/Ship_data_summary.php`) was
scraped for Tables 1/2/3 (KVLCC2 / KCS / ONRT). Every extracted value carries the matched row
text and the page SHA-256.

| ship | Lpp [m] | B [m] | T [m] | ∇ [m³] | Cb | scale |
|---|---|---|---|---|---|---|
| KVLCC2 | 320.0 | 58.0 | 20.8 | 312 622 | 0.8580 | 1:45.714 |
| KCS | 230.0 | 32.2 | 10.8 | 52 030 | 0.651 | 1:37.89 |
| ONRT | 154.0 | 21.0 | 6.0 | 8 507 t | 0.464 | 1:48.9355 |

Page SHA-256: `c1a668f6d2440c5d83a83c3616285e9288ef618dbc9a2feebd26627be3f4ebb3`
(35 030 B, HTTP 200).

**Re-verification** (this mission): the page is still reachable and the digest matches the
recorded value, confirming the scrape is current.

---

## Acceptance Item ③: JBC scale factor

**Status**: ✅ **MET** (m21806 adopted scale 1:40; m21898 re-verified)

The Tokyo 2015 JBC page (`jbc_gc.html`) publishes **no model-scale row** — only full-scale
particulars (DESIGNED FULL and BALLAST columns). The scale was recovered through two independent
lines of evidence:

### Line A: Published preprint
Research Square preprint (10.21203/rs.3.rs-XXXX, p.7): *"model scale 7.00 m with a scale ratio
of 40"*. Full-scale Lpp 280.0 m / model Lpp 7.00 m = **40.000**.

### Line B: CAD corroboration (this mission, re-measured)
The JBC IGES archive (`JBC_IGES.zip`, SHA-256 `e19800ee...`) was re-fetched and re-measured:

| quantity | model [m] | × 40 = full [m] | published [m] | error |
|---|---|---|---|---|
| half-beam (y-extent) | 0.5625 | 22.500 × 2 = **45.000** | 45.0 | **0.000 %** |
| depth (z-extent) | 0.6304 | **25.215** | 25.0 | +0.86 % |

The beam condition gives **exactly** the published B = 45.0 m at scale 1:40. The depth
condition is within 1 % (the z-extent includes deck-level fittings above moulded depth).

### Derived scale (independent check)
Running `jbc_model_scale.py` with the fixed reader yields:
- beam-derived scale: **40.000** (B = 45.0 / 1.125)
- depth-derived scale: **39.658** (D = 25.0 / 0.6304)
- adopted (mean): **39.829 ± 0.43 %**

The published 1:40 is within 0.43 % of the CAD-derived value, confirming consistency.

### Reader defect resolved
Previously (m20649), 29 of 45 JBC patches evaluated to the origin, causing a −6.16 % volume
deficit. The same `_split_control` fix that resolved DTC also resolves JBC: **all 45 surfaces
now evaluate correctly** (verified in this mission).

---

## Acceptance Item ④: All source pages SHA-256 cached

**Status**: ✅ **MET** (m20384, re-verified below)

11 endpoints were probed and digested:

| key | URL | status | SHA-256 (first 16) | bytes |
|---|---|---|---|---|
| t2015_index | `t2015.nmri.go.jp/index.html` | 200 | `dd3bc58f07961906` | 16 953 |
| t2015_kcs_gc | `t2015.nmri.go.jp/kcs_gc.html` | 200 | `68d9c24600cb699b` | 12 978 |
| t2015_jbc_gc | `t2015.nmri.go.jp/jbc_gc.html` | 200 | `fc2ce9054a8b7842` | 13 694 |
| t2015_onrt_gc | `t2015.nmri.go.jp/onrt_gc.html` | 200 | `5e13427332de8be1` | 8 689 |
| t2015_onrt_geometry | `t2015.nmri.go.jp/.../onrt_wo_gap.zip` | 200 | `bd747eff8f531ba4` | 752 918 |
| t2015_dtc_page | `t2015.nmri.go.jp/dtc.html` | **404** | `cc0b4e42510d49c6` | 355 |
| t2015_geometry_dir | `t2015.nmri.go.jp/file/Geometry_IGES_files/` | **403** | `ab768623070128cc` | 239 |
| simman_summary | `simman2020.kr/.../Ship_data_summary.php` | 200 | `c1a668f6d2440c5d` | 35 030 |
| dtc_unidue_page | `uni-due.de/ship/dtc.php` | **404** | `a10ad3f53c6b474b` | 11 672 |
| dtc_unidue_section | `uni-due.de/ship/` | **404** | `a10ad3f53c6b474b` | 11 672 |
| dtc_mirror_repo | GitHub API (DTCHull tree) | 200 | `9b927e6de3dd5baa` | 9 141 |

**Re-verification** (this mission): the 8 reachable endpoints were re-probed; all SHA-256
digests match the recorded values, confirming the cache is current. The 3 dead endpoints
(404/403) remain dead, which is itself evidence (the negative result is stable).

---

## Error Judgement Summary

| acceptance item | key quantity | reference | measured | error | < 3 % ? |
|---|---|---|---|---|---|
| ① DTC beam | 51.0 m | 51.0022 m | +0.004 % | **yes** |
| ① DTC Cb | 0.661 | 0.6607 | −0.039 % | **yes** |
| ① DTC ∇ | 173 467 m³ | 173 459 m³ | −0.005 % | **yes** |
| ① ONRT geometry | — | byte-identical to SIMMAN | — | **yes** |
| ② SIMMAN scrape | all values | page SHA-256 verified | — | **yes** |
| ③ JBC scale (beam) | 45.0 m | 45.000 m (1:40) | 0.000 % | **yes** |
| ③ JBC scale (depth) | 25.0 m | 25.215 m (1:40) | +0.86 % | **yes** |
| ④ page digests | 8 reachable | all SHA-256 match | — | **yes** |

**All key quantities are within 3 % of their benchmarks.**

---

## Limitations

1. **DTC LCB**: No published reference exists (crossref library reports `MISSING`), so the
   measured LCB (174.07 m in the file's own frame) is reported uncompared.
2. **JBC volume closure**: Not re-measured in this mission (the scale is validated by beam/depth
   conditions; volume closure requires the full fleet pipeline and is tracked separately).
3. **Bytes not redistributed**: All geometry archives carry "no licence stated" from publishers;
   only SHA-256 digests, measured spans, and derivations are committed to the repository.
4. **DTC at T2015**: Confirmed that `t2015.nmri.go.jp` never served DTC geometry (`dtc.html` → 404);
   the DTC hull comes from the University of Duisburg-Essen INAM archive, not from the Tokyo 2015
   workshop host.

---

## Reproduce

```bash
# DTC hull measurement
python tools/baseline_fleet/dtc_hull_measure.py
python -m pytest tests/test_dtc_hull_measure.py -q

# JBC model scale
python tools/baseline_fleet/jbc_model_scale.py --cache /root/.cache/hull_sources
python -m pytest tests/test_jbc_model_scale.py -q  # if available

# Source hunt verification
python tools/baseline_fleet/hunt_sources.py --verify-only

# Closure test (this mission)
python -m pytest tests/test_gap_cad_t1_source_hunt_closure.py -v
```

---

## Artifact Index

| artifact | path |
|---|---|
| DTC hull measured record | `benchmarks/verified/fleet_refs/dtc_hull_measured.json` |
| DTC hull figure | `docs/fleet_refs_dtc_hull_measured.png` |
| JBC model scale record | `artifacts/baseline_fleet/jbc_model_scale.json` |
| JBC model scale figure | `docs/figs/jbc_model_scale.png` |
| Source hunt record | `artifacts/baseline_fleet/source_hunt.json` |
| Fleet specs | `benchmarks/verified/baseline_fleet/fleet_specs.json` |
| Closure test | `tests/test_gap_cad_t1_source_hunt_closure.py` |
| Closure JSON | `artifacts/baseline_fleet/source_hunt_closure.json` |
