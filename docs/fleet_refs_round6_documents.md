# Fleet reference library — the document axis (round 6)

Anchor **REPRO-GEOM-FLEET-REFS** (TensorCAD). Machine artefacts:
`benchmarks/verified/fleet_refs/crossref.json` (new `documents` block),
`benchmarks/verified/fleet_refs/gap_ledger.json`.
Regenerate: `python tools/baseline_fleet/refs_round6_documents.py`.
Every digit below is computed from those files.

## 1. Why a second axis

Criterion ① is judged on *lineages*: two documents that print one upstream table are
one source, because their agreement is a statement about copying, not about the hull.
That rule is kept. But it hid the documents: `JBC` was reported `SINGLE_SOURCE` while
three publishers state its particulars. The two axes answer different questions —

* *documents*: is the numeral printed by more than one body?
* *lineages*: is there an origin independent of the first one?

## 2. Per-cell verdict (`verdict` = documents × lineages, `lineage_verdict` = library)

| ship | ∇ (m³) | Cb | LCB (% Lpp, fwd+) | Lwl (m) | Lpp (m) | Bwl (m) | T (m) |
|---|---|---|---|---|---|---|---|
| DTC | INDEPENDENT (3d/2l) | INDEPENDENT (3d/2l) | MISSING (0d/0l) | MISSING (0d/0l) | INDEPENDENT (3d/2l) | INDEPENDENT (3d/2l) | INDEPENDENT (3d/2l) |
| JBC | LINEAGE_LIMITED (4d/1l) | LINEAGE_LIMITED (3d/1l) | SINGLE_DOCUMENT (1d/1l) | LINEAGE_LIMITED (4d/1l) | LINEAGE_LIMITED (4d/1l) | LINEAGE_LIMITED (4d/1l) | LINEAGE_LIMITED (4d/1l) |
| KCS | INDEPENDENT (3d/2l) | INDEPENDENT (3d/2l) | INDEPENDENT (3d/2l) | INDEPENDENT (3d/2l) | INDEPENDENT (3d/2l) | INDEPENDENT (3d/2l) | INDEPENDENT (3d/2l) |
| KVLCC2 | INDEPENDENT (5d/2l) | INDEPENDENT (5d/2l) | LINEAGE_LIMITED (3d/1l) | LINEAGE_LIMITED (3d/1l) | INDEPENDENT (6d/3l) | INDEPENDENT (5d/2l) | INDEPENDENT (6d/3l) |
| ONRT | CONFLICT (2d/2l) | INDEPENDENT (2d/2l) | SINGLE_DOCUMENT (1d/1l) | INDEPENDENT (2d/2l) | MISSING (0d/0l) | INDEPENDENT (2d/2l) | INDEPENDENT (2d/2l) |
| Series60 | MISSING (0d/0l) | MISSING (0d/0l) | MISSING (0d/0l) | MISSING (0d/0l) | MISSING (0d/0l) | MISSING (0d/0l) | MISSING (0d/0l) |

`Nd/Ml` = N voting documents from M lineages. `INDEPENDENT` needs ≥2 lineages;
`LINEAGE_LIMITED` means ≥2 documents but one origin — corroboration of transcription,
never a second measurement.

## 3. The anchor's four named quantities

| ship | documents on ∇/Cb/LCB/Lwl | independent (≥2 lineages) | lineage-limited | single document | missing | documents in library |
|---|---|---|---|---|---|---|
| DTC | 2/4 | 2/4 | 0 | 0 | 2 | 3 |
| JBC | 3/4 | 0/4 | 3 | 1 | 0 | 4 |
| KCS | 4/4 | 4/4 | 0 | 0 | 0 | 3 |
| KVLCC2 | 4/4 | 2/4 | 2 | 0 | 0 | 6 |
| ONRT | 3/4 | 2/4 | 0 | 1 | 0 | 2 |
| Series60 | 0/4 | 0/4 | 0 | 0 | 4 | 0 |

Ships with ≥2 documents on all four named quantities: KCS, KVLCC2
.
Ships independent (≥2 lineages) on all four: KCS
.
The first list is the owner's literal wording (`≥2 独立来源`, journal papers and
workshop pages both named as sources); the second is the strict reading the library
already enforced. Neither is presented as the other.

## 4. What this round did *not* achieve

8 named cells still carry fewer than two documents. The
search actually run is recorded per cell in `gap_ledger.json`; no candidate was
fetched or pinned, so no verdict moved. In particular Series 60 still has no pinned
document, and DTC has no document at all for LCB or Lwl.

| ship | quantity | documents | verdict | search run this mission | candidates |
|---|---|---|---|---|---|
| DTC | LCB (% Lpp, fwd+) | 0 | MISSING | Duisburg Test Case DTC container ship LCB waterline length | 2 |
| DTC | Lwl (m) | 0 | MISSING | Duisburg Test Case DTC container ship LCB waterline length | 2 |
| JBC | LCB (% Lpp, fwd+) | 1 | SINGLE_DOCUMENT | — | 0 |
| ONRT | LCB (% Lpp, fwd+) | 1 | SINGLE_DOCUMENT | ONR tumblehome DTMB 5613 main particulars LCB | 0 |
| Series60 | ∇ (m³) | 0 | MISSING | Series 60 Cb 0.6 hull main particulars table resistance | 4 |
| Series60 | Cb | 0 | MISSING | Series 60 Cb 0.6 hull main particulars table resistance | 4 |
| Series60 | LCB (% Lpp, fwd+) | 0 | MISSING | Series 60 Cb 0.6 hull main particulars table resistance | 4 |
| Series60 | Lwl (m) | 0 | MISSING | Series 60 Cb 0.6 hull main particulars table resistance | 4 |

## 5. How this was verified

```
python tools/baseline_fleet/refs_round6_documents.py
python -m pytest tests/test_fleet_refs_document_axis.py tests/test_fleet_refs_crossref.py \
                 tests/test_fleet_refs_adjudication.py tests/test_fleet_refs_freshness.py -q
python tools/baseline_fleet/refs_acceptance.py check
```

![documents vs lineages](docs/figs/fleet_refs_round6_documents.png)
