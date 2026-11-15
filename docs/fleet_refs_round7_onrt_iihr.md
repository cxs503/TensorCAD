# Round 7 — an independent ONRT witness and the basis of the one disputed cell

Anchor `REPRO-GEOM-FLEET-REFS`, criterion ① (two independent lineages per quantity).

## Document pinned this round

* University of Iowa (IIHR) thesis 2016, Bottiglieri - Table 2.1 main particulars
* URL `https://iro.uiowa.edu/view/delivery/01IOWA_INST/12730586600002771/13730793940002771`
* SHA-256 `6fd4b42f0f0b873b96114fddeac77dce962a73d7ddc11ea55470170dcf91961d` — matches pin: **True**
* independence family `IIHR-2016`, lineage `IIHR-2016` (the IIHR basin measured this model; the table is not a reprint of a workshop page)

## Numerals read (verbatim quote kept per cell)

| quantity | value | raw | basis |
|---|---:|---|---|
| lwl | 154.0 | `154` | volume |
| beam | 18.78 | `18.78` | volume |
| draft | 5.494 | `5.494` | volume |
| displacement_volume | 8299.5122 | `8507` | mass |
| cb | 0.535 | `0.535` | volume |

## The disputed cell, decided by arithmetic (not by averaging)

| family | value | basis as stored |
|---|---:|---|
| SIMMAN-SERIES | 8507.0 | volume |
| NMRI-T2015 | 8299.5 | mass |

Closure of the newly read rows: Cb·Lwl·B·T = 8500.8 m³
→ residual **0.0731 %** if 8507 is a volume, **2.425 %** if it is a mass.
Verdict: `UNIT_BASIS_CONFLICT` — the IIHR table's numeral 8507 closes as a volume (residual 0.0731 % vs 2.425 % if read as mass); read as a mass it converts to 8299.5 m3, i.e. the NMRI-T2015 value. The disagreement is therefore a basis question, not two measurements that differ: no published value is changed here, and the cell must keep both readings until a document states the basis explicitly.

## LCB does not vote

`LCB LCB (m) aft of FP 1.625 NA` — quoted model value 1.625 m, which would be 79.625 m aft of FP at 1/49; the full-scale column reads `NA`, so `ONRT/lcb_percent` remains `SINGLE_DOCUMENT`. A future round needs a full-scale LCB in print (SIMMAN 2020 ONRT detail page or a DTMB 5613 model report).

## What still blocks criterion ①

| ship | quantity | status after this round | next document wanted |
|---|---|---|---|
| DTC | lcb_percent, lwl | MISSING | Moctar/Shigunov/Zorn 2012 primary table, or a thesis reprint |
| JBC | lcb_percent | SINGLE_DOCUMENT | NMRI 2015 JBC report (LCB is not on the geometry page) |
| ONRT | lcb_percent | SINGLE_DOCUMENT | full-scale LCB in print |
| Series60 | ∇, Cb, LCB, Lwl | MISSING | Todd's parent-form table reprinted by an independent body |

Nothing above is a guess: each row is the state of `gap_ledger.json` after the round, and a candidate is a hunch until its SHA-256 is pinned and its quote read.

## Reproduce

```
python tools/baseline_fleet/refs_round7_iihr.py --raw artifacts/baseline_fleet/raw
python -m pytest tests/test_fleet_refs_round7_iihr.py
```

Figure: `docs/figs/fleet_refs_round7_onrt_iihr.png` (published ONRT displacement by family, with the closure-implied volume drawn as a band, not as a vote).
