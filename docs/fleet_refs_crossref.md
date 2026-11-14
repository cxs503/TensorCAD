# Fleet reference library — six benchmark hulls, ≥2 independent sources each

Anchor **REPRO-GEOM-FLEET-REFS** (TensorCAD). Machine artefact: `benchmarks/verified/fleet_refs/crossref.json`; this document is rendered from it by `tools/baseline_fleet/crossref_report.py` (no number is transcribed by hand).

## 1. Problem setting and where the references come from

The fleet benchmark compares a reconstructed displacement against **one** published number per hull (KCS verified at +0.304 % against the NMRI page, everything else withheld). A single number cannot separate a reconstruction error from a source typo. This library therefore holds every principal particular that a publisher states for the six hulls — ∇, Cb, LCB (%Lpp), Lwl, and the Lpp/Bwl/T the coefficients are defined on — reads them out of the fetched bytes, and asks per quantity whether **two independent publishers** state it and agree.

Sources (digest pinned; the bytes stay out of the repository, see the licence column):

| source | publisher | family | status | sha256 (16) |
|---|---|---|---|---|
| `t2015_kcs` | — | — | read | 68d9c24600cb699b |
| `t2015_jbc` | — | — | read | fc2ce9054a8b7842 |
| `t2015_onrt` | — | — | read | 5e13427332de8be1 |
| `simman2020_summary` | — | — | read | c1a668f6d2440c5d |
| `simman2014_kvlcc2` | — | — | read | bfc9101e3a4afa4f |
| `simman2014_kcs` | — | — | read | 3cbc3a22daf3d494 |
| `ittc_benchmarks` | — | — | read | da101e9a352108e4 |

## 2. Method

`tensorcad/fleet/crossref.py` reads a table row by row (`official.table_rows`), matches the quantity label **in the first cell only** (`LABEL_PATTERNS`) and takes the numeric cell of the publisher's full-scale column. Two guards keep the reading honest: a row whose label is a derived quantity (`izz/Lpp`) is skipped, and a row whose configured full-scale column is empty is kept for the audit trail but flagged `fallback_column` and never allowed to vote. Independence is judged at the *lineage* level, not the URL: the SIMMAN 2014 and SIMMAN 2020 editions print the identical KVLCC2 table including the identical `Bwl 580` defect, so they count as one source. Each ship/quantity is then classified `CROSS_VERIFIED` (≥2 lineages, relative spread ≤ tolerance), `CONFLICT` (≥2 lineages, spread above tolerance), `SINGLE_SOURCE` or `MISSING`.

## 3. Parameter table (as published, full scale)

| hull | quantity | published value(s) by lineage | verdict |
|---|---|---|---|
| KCS | ∇ (m³) | NMRI-T2015 52,030.0, SIMMAN-SERIES 52,030.0 | CROSS_VERIFIED |
| KCS | Cb | NMRI-T2015 0.6505, SIMMAN-SERIES 0.651 | CROSS_VERIFIED |
| KCS | LCB (% Lpp, fwd+) | NMRI-T2015 -1.48, SIMMAN-SERIES -1.48 | CROSS_VERIFIED |
| KCS | Lwl (m) | NMRI-T2015 232.5, SIMMAN-SERIES 232.5 | CROSS_VERIFIED |
| KCS | Lpp (m) | NMRI-T2015 230, SIMMAN-SERIES 230 | CROSS_VERIFIED |
| KCS | Bwl (m) | NMRI-T2015 32.2, SIMMAN-SERIES 32.2 | CROSS_VERIFIED |
| KCS | T (m) | NMRI-T2015 10.8, SIMMAN-SERIES 10.8 | CROSS_VERIFIED |
| KVLCC2 | ∇ (m³) | SIMMAN-SERIES 312,622.0 | SINGLE_SOURCE |
| KVLCC2 | Cb | SIMMAN-SERIES 0.8098 | SINGLE_SOURCE |
| KVLCC2 | LCB (% Lpp, fwd+) | SIMMAN-SERIES 3.48 | SINGLE_SOURCE |
| KVLCC2 | Lwl (m) | SIMMAN-SERIES 325.5 | SINGLE_SOURCE |
| KVLCC2 | Lpp (m) | SIMMAN-SERIES 320 | SINGLE_SOURCE |
| KVLCC2 | Bwl (m) | SIMMAN-SERIES 580 | SINGLE_SOURCE |
| KVLCC2 | T (m) | SIMMAN-SERIES 20.8 | SINGLE_SOURCE |
| JBC | ∇ (m³) | NMRI-T2015 178,369.9 | SINGLE_SOURCE |
| JBC | Cb | NMRI-T2015 0.858 | SINGLE_SOURCE |
| JBC | Lwl (m) | NMRI-T2015 285 | SINGLE_SOURCE |
| JBC | Lpp (m) | NMRI-T2015 280 | SINGLE_SOURCE |
| JBC | Bwl (m) | NMRI-T2015 45 | SINGLE_SOURCE |
| JBC | T (m) | NMRI-T2015 16.5 | SINGLE_SOURCE |
| ONRT | ∇ (m³) | NMRI-T2015 8,299.5, SIMMAN-SERIES 8,507.0 | CONFLICT |
| ONRT | Cb | NMRI-T2015 0.535, SIMMAN-SERIES 0.535 | CROSS_VERIFIED |
| ONRT | LCB (% Lpp, fwd+) | NMRI-T2015 1.625, SIMMAN-SERIES -1.68 | CONFLICT |
| ONRT | Lwl (m) | NMRI-T2015 154, SIMMAN-SERIES 154 | CROSS_VERIFIED |
| ONRT | Bwl (m) | NMRI-T2015 18.78, SIMMAN-SERIES 18.78 | CROSS_VERIFIED |
| ONRT | T (m) | NMRI-T2015 5.494, SIMMAN-SERIES 5.494 | CROSS_VERIFIED |

## 4. Cross-check result — does each hull have ≥2 independent sources?

| hull | ∇ | Cb | LCB | Lwl |
|---|---|---|---|---|
| KCS | CROSS_VERIFIED (2) | CROSS_VERIFIED (2) | CROSS_VERIFIED (2) | CROSS_VERIFIED (2) |
| KVLCC2 | SINGLE_SOURCE (1) | SINGLE_SOURCE (1) | SINGLE_SOURCE (1) | SINGLE_SOURCE (1) |
| JBC | SINGLE_SOURCE (1) | SINGLE_SOURCE (1) | MISSING (0) | SINGLE_SOURCE (1) |
| ONRT | CONFLICT (2) | CROSS_VERIFIED (2) | CONFLICT (2) | CROSS_VERIFIED (2) |
| DTC | MISSING (0) | MISSING (0) | MISSING (0) | MISSING (0) |
| Series60 | MISSING (0) | MISSING (0) | MISSING (0) | MISSING (0) |

Coverage: **11/42** cells cross-verified, 13 single-source, 2 conflict, 16 missing. Hulls with all four key quantities (∇, Cb, LCB, Lwl) backed by ≥2 lineages: **KCS**. Hulls still below two sources: **DTC, JBC, KVLCC2, Series60**.

### 4.1 Reference vs reconstruction

| hull | ref ∇ (m³) | got ∇ (m³) | error | ref verdict | benchmark |
|---|---|---|---|---|---|
| KCS | 52,030.0 | 52,188.0 | +0.30 % | CROSS_VERIFIED | VERIFIED |
| KVLCC2 | 312,622.0 | — | — | SINGLE_SOURCE | WITHHELD |
| JBC | 178,369.9 | 145,769.0 | -18.28 % | SINGLE_SOURCE | WITHHELD |
| ONRT | 8,299.5 | 92.9489 | — | CONFLICT | WITHHELD (got artefact not comparable) |

### 4.2 Publisher self-consistency (Cb closure)

| hull | lineage | Cb·Lpp·Bwl·T (m³) | published ∇ (m³) | rel. error | verdict |
|---|---|---|---|---|---|
| KCS | NMRI-T2015 | 52,030.1 | 52,030.0 | +0.000 | ok |
| KCS | SIMMAN-SERIES | 52,070.1 | 52,030.0 | +0.001 | ok |
| KVLCC2 | SIMMAN-SERIES | 3,126,216.7 | 312,622.0 | +9.000 | **defect** — implied Bwl 58.00 m |
| JBC | NMRI-T2015 | 178,378.2 | 178,369.9 | +0.000 | ok |

## 5. Sensitivity / what the classification is sensitive to

Two knobs decide the verdicts, and both are exercised in `tests/test_fleet_refs_crossref.py`:

* **tolerance** — 0.5 % for ∇/Cb/Lwl/Lpp/Bwl/T and 2 % for LCB (a small percent-of-Lpp number); at 0.5 % the KCS Cb pair (0.6505 vs 0.651, 0.077 %) is `CROSS_VERIFIED`, and lowering the tolerance to 0.05 % would flip it to `CONFLICT`: the two publishers round differently, which is corroboration, not disagreement.
* **lineage grouping** — treating SIMMAN 2014 and SIMMAN 2020 as two sources would report KVLCC2 as cross-verified on all seven quantities; the closure test (section 4.2) shows the two editions share the `Bwl 580` digit, so they are one source. This is the single most consequential choice in the library.

## 6. Failure modes and artefacts

Honest gaps, with the reason each document could not be used:

| wanted | url | outcome |
|---|---|---|

* **DTC** and **Series 60**: no numeric source reachable at all in this environment (0 cells). DTC's originator page is retired and its IGES bundle (`uni-due.de/imperia/md/content/inam/igs_ordner.zip`, sha256 recorded in the registry) carries geometry but no particulars; Series 60's primary table (Todd) is behind a 403.
* **JBC**: one lineage only (NMRI Tokyo 2015). Its ∇ 178 369.9 m³ is *not* yet external in the sense this anchor requires; the reconstruction currently reads 145 769 m³ (−18.3 %) with the control-net layout guard fired, so the hull is `WITHHELD` in the benchmark for a reason that must not be blamed on the reference.
* **ONRT**: ∇ and LCB are `CONFLICT` — the NMRI page prints displacement as a *mass* (8 507 t → 8 299.5 m³ at 1.025 t/m³) while SIMMAN prints 8 507 m³, and LCB is printed on different references (the NMRI cell is a metre offset that the parser reads as 1.625, the SIMMAN cell is −1.68 %Lpp). Both cells are kept, unflagged as agreement, because picking either one would be a flattering choice.

## Reproduce

```bash
python tools/baseline_fleet/crossref_refs.py fetch --out artifacts/fleet_refs   # network; digests checked
python tools/baseline_fleet/crossref_refs.py build --raw artifacts/fleet_refs/raw
python tools/baseline_fleet/crossref_report.py \
    --library benchmarks/verified/fleet_refs/crossref.json \
    --doc docs/fleet_refs_crossref.md --figs docs/figs/fleet_refs
python -m pytest tests/test_fleet_refs_crossref.py -q
```

