# Fleet reference blockers - anchor REPRO-GEOM-FLEET-REFS (round 9)

Generated 2026-11-15T22:54:25+00:00 from `benchmarks/verified/fleet_refs/crossref.json` by `tools/baseline_fleet/refs_blockers.py`. Machine-readable form: `benchmarks/verified/fleet_refs/blockers.json`.

## What changed this round

This round adds no reference numeral. It converts the remaining acceptance failure from a verdict into a per-cell, evidence-backed work item: for every quantity that is not cross-verified it names *which kind* of gap it is and *which external document* would close it, using only what the pinned library already records.

## Coverage

* cells checked: **42** (6 ships x 7 quantities)
* cross-verified: **21**
* still blocked: **21**
* blockers by kind: `CONVERTED_ONLY` 1, `DISPUTED` 1, `NO_DOCUMENT` 9, `SINGLE_LINEAGE` 10
* ships still failing the anchor's >=2-source criterion on the four named quantities: **DTC, JBC, KVLCC2, ONRT, Series60**

## Per-ship status on the four named quantities (displacement, Cb, LCB, Lwl)

| ship | passes >=2 independent sources | quantities short |
|---|---|---|
| DTC | no | lcb_percent, lwl |
| JBC | no | displacement_volume, cb, lcb_percent, lwl |
| KCS | yes | - |
| KVLCC2 | no | lcb_percent, lwl |
| ONRT | no | displacement_volume, lcb_percent |
| Series60 | no | displacement_volume, cb, lcb_percent, lwl |

## Blockers

| ship | quantity | verdict | blocker | families | missing | voting sources | converted reading |
|---|---|---|---|---|---|---|---|
| DTC | lcb_percent | MISSING | `CONVERTED_ONLY` | 0 | 2 | - | kvlcc2_ugent (converted, not a vote) |
| DTC | lwl | MISSING | `NO_DOCUMENT` | 0 | 2 | - | - |
| JBC | displacement_volume | SINGLE_SOURCE | `SINGLE_LINEAGE` | 1 | 1 | jbc_cfdletters, jbc_memorial, jbc_scse, t2015_jbc | - |
| JBC | cb | SINGLE_SOURCE | `SINGLE_LINEAGE` | 1 | 1 | jbc_cfdletters, jbc_scse, t2015_jbc | - |
| JBC | lcb_percent | SINGLE_SOURCE | `SINGLE_LINEAGE` | 1 | 1 | jbc_scse | - |
| JBC | lwl | SINGLE_SOURCE | `SINGLE_LINEAGE` | 1 | 1 | jbc_cfdletters, jbc_memorial, jbc_scse, t2015_jbc | - |
| JBC | length_pp | SINGLE_SOURCE | `SINGLE_LINEAGE` | 1 | 1 | jbc_cfdletters, jbc_memorial, jbc_scse, t2015_jbc | - |
| JBC | beam | SINGLE_SOURCE | `SINGLE_LINEAGE` | 1 | 1 | jbc_cfdletters, jbc_memorial, jbc_scse, t2015_jbc | - |
| JBC | draft | SINGLE_SOURCE | `SINGLE_LINEAGE` | 1 | 1 | jbc_cfdletters, jbc_memorial, jbc_scse, t2015_jbc | - |
| KVLCC2 | lcb_percent | SINGLE_SOURCE | `SINGLE_LINEAGE` | 1 | 1 | kvlcc2_osaka, simman2014_kvlcc2, simman2020_summary | kvlcc2_ugent (converted, not a vote) |
| KVLCC2 | lwl | SINGLE_SOURCE | `SINGLE_LINEAGE` | 1 | 1 | kvlcc2_osaka, simman2014_kvlcc2, simman2020_summary | - |
| ONRT | displacement_volume | CONFLICT | `DISPUTED` | 2 | 0 | simman2020_summary, t2015_onrt | - |
| ONRT | lcb_percent | SINGLE_SOURCE | `SINGLE_LINEAGE` | 1 | 1 | simman2020_summary | t2015_onrt (converted, not a vote) |
| ONRT | length_pp | MISSING | `NO_DOCUMENT` | 0 | 2 | - | - |
| Series60 | displacement_volume | MISSING | `NO_DOCUMENT` | 0 | 2 | - | - |
| Series60 | cb | MISSING | `NO_DOCUMENT` | 0 | 2 | - | - |
| Series60 | lcb_percent | MISSING | `NO_DOCUMENT` | 0 | 2 | - | - |
| Series60 | lwl | MISSING | `NO_DOCUMENT` | 0 | 2 | - | - |
| Series60 | length_pp | MISSING | `NO_DOCUMENT` | 0 | 2 | - | - |
| Series60 | beam | MISSING | `NO_DOCUMENT` | 0 | 2 | - | - |
| Series60 | draft | MISSING | `NO_DOCUMENT` | 0 | 2 | - | - |

## Targets and why the pinned set is not enough

### DTC / lcb_percent (`CONVERTED_ONLY`)

- why the pinned set is not enough: both pinned DTC papers were re-read this round: dtc_strath Table 1 and dtc_ije 'Geometrical Characteristics' print Lpp / Bwl / T / volume / CB / WSA only - neither prints LCB or Lwl; dtc_ugent (JSR-2017) prints LCG for the DTC design condition, i.e. a converted reading, not a stated LCB.
- would close it: El Moctar, Shigunov & Zorn (2012), 'Duisburg Test Case: Post-Panamax Container Ship for Benchmarking', Ship Technology Research 59(3) - the DTC primary description; states the design-condition LCB/LCG.
- would close it: A thesis that tabulates DTC including LCB and Lwl (e.g. the Ghent record biblio.ugent.be/publication/8516555 noted by the previous round).
- evidence read: dtc_strath: Table 1 'Main particulars of DTC hull' = lambda, Lpp 355.0, Bwl 51.0, T 14.5, volume 173467.0, CB 0.661, WSA 22032.0 (no LCB, no Lwl)
- evidence read: dtc_ije: 'Geometrical Characteristics' table = Lpp 355, Bwl 51, T 14.5, volume 173468, CB 0.661, WSA 22032 (no LCB, no Lwl)
- reachable from this environment: no

### DTC / lwl (`NO_DOCUMENT`)

- why the pinned set is not enough: no pinned DTC document prints Lwl; the Lwl would have to be measured off the geometry file, which is a different kind of statement from a published cell.
- would close it: DTC primary description (El Moctar et al. 2012) - waterline length at design draft.
- would close it: DTC hull description published with the geometry package (INAM / uni-due).
- evidence read: same two Table 1 reads as DTC/lcb_percent
- reachable from this environment: no

### JBC / displacement_volume (`SINGLE_LINEAGE`)

- why the pinned set is not enough: only one lineage states the quantity in the pinned set
- would close it: a second independent publisher stating this quantity
- reachable from this environment: no

### JBC / lcb_percent (`SINGLE_LINEAGE`)

- why the pinned set is not enough: jbc_cfdletters / jbc_scse / jbc_memorial reproduce the NMRI table (lineage NMRI-T2015), so their agreement is transcription, not independence; a document from a different upstream (e.g. the JBC CFD workshop's own outcome report) is needed.
- would close it: NMRI Tokyo 2015 JBC geometry page (or the NMRI report behind it) - the page as pinned does not print LCB.
- would close it: A JBC paper that prints LCB from its own source: the three JBC PDFs pinned this round all carry lineage NMRI-T2015, so they cannot break the single lineage.
- evidence read: crossref documents: jbc_cfdletters, jbc_scse, jbc_memorial all lineage=NMRI-T2015
- evidence read: t2015_jbc quantities_read has no lcb_percent
- reachable from this environment: no

### KVLCC2 / lcb_percent (`SINGLE_LINEAGE`)

- why the pinned set is not enough: all printed values are SIMMAN series copies (3.48); the JSR-2017 paper prints LCG 171.1 m, which converts to 3.4688 %Lpp - agreement at 0.32 %, but a converted numeral.
- would close it: A second explicitly printed LCB (%Lpp, fwd+) for KVLCC2 from a lineage other than the SIMMAN series.
- evidence read: simman2014_kvlcc2 / simman2020_summary / kvlcc2_osaka: 3.48 (one lineage)
- evidence read: kvlcc2_ugent (JSR-2017): LCG 171.1 m -> 3.46875 %Lpp (unit_basis=derived)
- reachable from this environment: no

### KVLCC2 / lwl (`SINGLE_LINEAGE`)

- why the pinned set is not enough: every pinned Lwl reading is 325.5 m from the SIMMAN series (2014 page, 2020 summary, Osaka reprint); they are one upstream table.
- would close it: A second printed Lwl (m) for KVLCC2 from outside the SIMMAN series.
- evidence read: crossref: KVLCC2/lwl entries all family SIMMAN-SERIES, value 325.5
- reachable from this environment: no

### ONRT / lcb_percent (`SINGLE_LINEAGE`)

- why the pinned set is not enough: one publisher (SIMMAN 2020) prints LCB directly; the NMRI reading is a converted distance aft of FP, and a converted numeral is corroboration, never a vote.
- would close it: A second *printed* LCB for ONRT (the T2015 page states LCB as a distance aft of FP, which the library converts to -1.7045 %Lpp).
- evidence read: t2015_onrt cell: 'LCB (m) aft of FP 1.625' -> -1.704545 %Lpp (unit_basis=derived)
- evidence read: simman2020_summary: LCB (%) fwd+ = -1.68 (votes); residual 1.4 % relative, inside the 2 % LCB tolerance, so this is corroboration pending a printed second value
- reachable from this environment: no

### Series60 / displacement_volume (`NO_DOCUMENT`)

- why the pinned set is not enough: the library holds no Series 60 document at all; three Series-60 PDFs present in this environment's scratch cache (digests recorded in the dossier) were opened this round and print no parent-hull particulars - they are CFD papers about a Cb=0.6 model, not the series description.
- would close it: Todd, F.H. (1963), 'Series 60 Methodical Experiments with Models of Single-Screw Merchant Ships', David Taylor Model Basin Report 1712 - the parent-form tables.
- would close it: ITTC 7.5-02-06-06 benchmark-data appendix (the pinned ITTC document is the procedure index and yields no numeral).
- evidence read: crossref documents: no source with ship Series60 except the ITTC index (0 cells)
- evidence read: scratch-cache PDFs read this round: no 'LCB' / 'particulars' table for the parent form
- reachable from this environment: no

### Series60 / cb (`NO_DOCUMENT`)

- why the pinned set is not enough: no pinned Series 60 document states any numeral.
- would close it: Todd (1963) DTMB Report 1712 parent-form table (Cb = 0.60 parent).
- evidence read: same as Series60/displacement_volume
- reachable from this environment: no

### Series60 / lcb_percent (`NO_DOCUMENT`)

- why the pinned set is not enough: no pinned Series 60 document states any numeral.
- would close it: Todd (1963) DTMB Report 1712 - LCB (%Lpp) of the parent form.
- would close it: Validation studies of the Series 60 parent hull that tabulate LCB (e.g. towing-tank papers) as a *second lineage*.
- evidence read: same as Series60/displacement_volume
- reachable from this environment: no

### Series60 / lwl (`NO_DOCUMENT`)

- why the pinned set is not enough: no pinned Series 60 document states any numeral.
- would close it: Todd (1963) DTMB Report 1712 - waterline length of the parent form at design draught.
- evidence read: same as Series60/displacement_volume
- reachable from this environment: no

## Scope and limits

* Nothing here upgrades a verdict. `CONVERTED_ONLY` cells stay `SINGLE_SOURCE` in the library: a converted numeral is arithmetic the library performed, so promoting it would launder a transcription slip into 'two sources agree'.
* The dossier cannot tell whether an unpinned document exists; it only states that the *pinned* set does not print the quantity. That is a statement about evidence held, not about the world.
* Error judgement: no accuracy claim is made this round, because no reference value was added or changed. The one numeric check performed is a reconciliation - blocked cells must equal `acceptance.coverage.cells - acceptance.coverage.cross_verified`, and the generator refuses to emit a dossier when that identity fails (a drift guard, not a measurement).

## Next actions

1. Obtain the DTC primary description (El Moctar et al. 2012) for DTC LCB/Lwl - the only blocker pair that a single document closes.
2. Obtain Todd (1963) DTMB Report 1712 for Series 60: it is the only route to four of the ten remaining blockers, and Series 60 currently has no document at all.
3. For KVLCC2/ONRT LCB, hunt a *printed* percentage from a non-SIMMAN lineage; the converted readings already agree to 0.32 % and 1.4 % respectively, so a printed value would settle them in one fetch.

## Reproduction

```
python tools/baseline_fleet/refs_blockers.py \
    --library benchmarks/verified/fleet_refs/crossref.json \
    --acceptance benchmarks/verified/fleet_refs/acceptance.json \
    --out benchmarks/verified/fleet_refs/blockers.json \
    --doc docs/fleet_refs_round9_blockers.md \
    --figure docs/figs/fleet_refs_round9_blockers.png
python -m pytest tests/test_fleet_refs_blockers.py -q
```
