# Round 4 — the numerals the reference library was dropping

Anchor **REPRO-GEOM-FLEET-REFS** (TensorCAD). Machine artefact: `benchmarks/verified/fleet_refs/derived_cells.json` (regenerate: `python tools/baseline_fleet/refs_round4_derived.py`). Every number below is rendered from that file.

Documents re-read from the pinned bytes (no network): 6 (t2015_kcs, t2015_jbc, t2015_onrt, simman2020_summary, simman2014_kvlcc2, simman2014_kcs).

## 1. Cells that do not vote, and what they say

| ship | quantity | basis | derived value | voters | Δ vs voter | within tolerance |
|---|---|---|---|---|---|---|
| ONRT | LCB (% Lpp, fwd+) | `derived_length_aft_fp` | -1.6365 | SIMMAN-SERIES=-1.68 | +0.0435 | False |

## 2. Derivations in full

* ONRT LCB (% Lpp, fwd+): LCB = (L/2 - x_aft_FP)/L * 100 with L = Lwl = 154 m, x_aft_FP = 1.625 m * 48.9355 = 79.52 m
  * note: 1.625 m aft of FP at model scale -> 79.52 m full scale (48.9355) -> -1.6365 % of Lwl = 154 m, midship at L/2

## 3. Why these cells do not vote

A published numeral votes only when the publisher prints it in a column whose scale is known. A model-scale numeral cannot vote for a full-scale quantity (that is the `fallback_column` rule, enforced by `tensorcad.fleet.crossref.votes`). The derived readings above use the *document's own* model/full ratio, so they are arithmetic on one publisher's numbers and one length assumption -- evidence about the document, not a second opinion about the ship. None of them is counted in criterion ①.

