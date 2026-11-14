# GAP-CAD-T1-SOURCE-HUNT — geometry source hunt (T2015 / SIMMAN 2020 / DTC)

Round: m20384 (TensorCAD). Tool: `tools/baseline_fleet/hunt_sources.py`.
Record: `artifacts/baseline_fleet/source_hunt.json` (regenerate with
`python tools/baseline_fleet/hunt_sources.py --out artifacts/baseline_fleet --fetch-geometry`;
CI runs `--verify-only`, tests `tests/test_source_hunt.py`).

## 1. What was asked and what happened

| acceptance item | outcome | evidence |
|---|---|---|
| ① DTC + ONRT geometry retrieved, digest/URL/licence in a fleet manifest | **ONRT: MET. DTC: NOT MET (WITHHELD)** | ONRT archive from the T2015 host, SHA-256 `bd747ef…95 685`, 752 918 B, byte-identical to the SIMMAN copy on record (`geometry_cross_check.identical = true`). DTC: every publisher host probed is dead (table below) |
| ② SIMMAN 2020 `Ship_data_summary.php` main particulars scraped | **MET** | Tables 1/2/3 parsed for KVLCC2 / KCS / ONRT; every value carries the matched row text and the page SHA-256 `c1a668f…0c5d` |
| ③ JBC scale factor | **NOT MET (WITHHELD)** | `t2015_jbc_gc.html` has *no* model-scale row: its two numeric columns are "design" and "ballast", both full scale (∇ 178 369.9 / 89 185.9 m³). The offset must be measured from the IGES control-cloud span, which this round did not run |
| ④ every source page digested | **MET** | 11 endpoints recorded, 8 reachable with full SHA-256, 3 failures with digested error bodies |

## 2. HTTP outcome of every probe (2026-11 round, `source_hunt.json`)

| key | status | bytes | SHA-256 (first 16) |
|---|---|---|---|
| t2015_index | 200 | 16 953 | dd3bc58f07961906 |
| t2015_kcs_gc | 200 | 12 978 | 68d9c24600cb699b |
| t2015_jbc_gc | 200 | 13 694 | fc2ce9054a8b7842 |
| t2015_onrt_gc | 200 | 8 689 | 5e13427332de8be1 |
| t2015_onrt_geometry (zip) | 200 | 752 918 | bd747eff8f531ba4 |
| t2015_dtc_page (`/dtc.html`) | **404** | 355 | cc0b4e42510d49c6 |
| t2015_geometry_dir | **403** | 239 | ab768623070128cc |
| simman_summary | 200 | 35 030 | c1a668f6d2440c5d |
| dtc_unidue_page (`uni-due.de/ship/dtc.php`) | **404** | 11 672 | a10ad3f53c6b474b |
| dtc_unidue_section (`uni-due.de/ship/`) | **404** | 11 672 | a10ad3f53c6b474b |
| dtc_mirror_repo (GitHub tree API) | 200 | 9 141 | 9b927e6de3dd5baa |
| web.archive.org | no response (3 attempts, 30–90 s) | — | not reachable from this network |
| Zenodo API `"Duisburg Test Case"` | 200, 2 hits, none geometry | — | rejected |
| GitHub repo search `Duisburg Test Case` | 200, 0 items | — | rejected |

## 3. Findings

1. **The anchor's DTC lead is falsified for the T2015 host.** The workshop site serves
   KCS / JBC / KVLCC2 / ONRT / DTMB 5415 pages; `dtc.html` is 404 and the geometry
   directory is not listable (403). DTC was never a T2015 benchmark hull.
2. **The URL the DTC paper names is gone.** `uni-due.de/ship/dtc.php` and the enclosing
   section both return the university's 404 page (identical digest → a generic error
   page, so the institute section itself was removed, not just the file).
3. **ONRT is closed properly, and by two independent publishers.** The T2015-hosted
   archive is *byte-identical* to the SIMMAN 2020 copy already in the fleet registry,
   so the earlier "model-scale coordinates, no scale factor obtained" WITHHELD can be
   retired: the SIMMAN summary publishes the scale row (1 : 48.9355) and the T2015
   page publishes model 3.147 m / full 154.0 m LWL, matching it (154.0/3.147 = 48.94).
4. **Model scales for KCS / KVLCC2 / ONRT are now machine-read from one page**:
   1:37.89, 1:45.714, 1:48.9355; each is cross-checked against the Lpp pair in the same
   table in `tests/test_source_hunt.py::test_model_scale_factors`.
5. **DTC residual**: the only publicly reachable DTC hull surface found is a 21.2 MB
   `constant/triSurface/DTC-scaled.stl` inside a third-party GitHub copy of the OpenFOAM
   `DTCHull` tutorial. It carries no publisher digest and no publisher licence, so it is
   recorded as a lead and **not adopted** (`open_questions[0].residual`).

## 4. Error judgement (vs published references)

| quantity | reference | read | error | <3 % ? |
|---|---|---|---|---|
| KCS Lpp | 230.0 m (SIMMAN Tab. 2) | 230.0 m | 0 % | yes (identity by construction) |
| KVLCC2 ∇ | 312 622 m³ (Tab. 1) | 312 622 m³ | 0 % | yes |
| ONRT ∇ | 8507 t (T2015 `onrt_gc`) | 8507 t | 0 % | yes |
| ONRT ∇, SIMMAN vs T2015 | 8507 t | 8507 t | 0 % | yes — the earlier 2.43 % gap came from comparing a mass-derived column with the Cb·L·B·T closure, not from the publishers disagreeing |
| ONRT scale, SIMMAN vs LWL pair | 48.9355 | 154.0 / 3.147 = 48.9355 | 9×10⁻⁶ | yes |

No *mesh* was rebuilt in this round, so the "mesh volume vs published ∇ < 1 %" row of
`benchmarks/verified/baseline_fleet/README.md` stays **not run** here — this round
touched provenance only.

## 5. Limitations

* DTC geometry: **WITHHELD**. No publisher-hosted copy reachable; archive.org unreachable
  from this network, so a historical snapshot could not be tested.
* JBC scale: **WITHHELD**. Needs the IGES control-cloud span; the diagnostics file
  (`benchmarks/verified/baseline_fleet/geometry_diagnostics.json`) records per-hull spans
  but its paths point at another work tree (`/root/tp-agent/work/m20163/...`), so the
  derivation is deferred rather than guessed.
* No bytes redistributed: pages and archives live in `artifacts/baseline_fleet/raw/`
  (git-ignored); git carries the digests, the matched snippets and the fetch script.
