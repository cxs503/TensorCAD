# KVLCC2 Lwl — the single most impactful remaining cross-reference gap

Anchor: **REPRO-GEOM-FLEET-REFS** (mission m22804, branch `agent/m22804`).
Code: `tensorcad/fleet/crossref.py`, `benchmarks/verified/fleet_refs/crossref.json`.

## 1. Why this round exists

The cross-source reference library (`benchmarks/verified/fleet_refs/crossref.json`) scores
the six benchmark hulls against five acceptance criteria.  As of this round the score is:

| criterion | pass | failing cells |
|---|---|---|
| ① 每船 ≥2 独立 lineage（∇/Cb/LCB/Lwl） | **FAIL** | KVLCC2-lwl, DTC-lcb/lwl, JBC-∇/Cb/LCB/Lwl, ONRT-lcb(doc), Series60-∇/Cb/LCB/Lwl |
| ①b 每船 ≥2 独立 document | **FAIL** | DTC-lcb/lwl, JBC-lcb, ONRT-lcb, Series60-all |
| ② 每值带 citation+digest | PASS | — |
| ③ disputed 并列双值 | PASS | — |
| ④ specs 多源数组化 | PASS | — |
| ⑤ KVLCC2 参考接线 | PASS | — |

Criteria ②..⑤ are closed.  The remaining debt is **lineage independence** (①/①b): every
failing cell is backed by *one* publisher family, and the second family is either not
reachable from this egress or not yet registered.

## 2. The single most impactful cell: KVLCC2 `lwl`

KVLCC2 is the flagship benchmark (the tanker every resistance/CFD campaign cites).  Its
`lwl` = **325.5 m** is the one quantity that fails criterion ① while the rest of the ship
passes.  The current evidence:

| source | family | value | role |
|---|---|---|---|
| `simman2020_summary` | SIMMAN-SERIES | 325.5 | workshop page |
| `simman2014_kvlcc2` | SIMMAN-SERIES | 325.5 | workshop page |
| `kvlcc2_osaka` | SIMMAN-SERIES | 325.5 | Osaka dissertation (transcription) |

Three **documents** state it, but all three are the **same lineage** (SIMMAN-SERIES): the
Osaka dissertation and the two SIMMAN editions all reproduce the same case-definition
table.  Criterion ①b (≥2 documents) therefore already passes for this cell; it is
criterion ① (≥2 *independent* lineages) that fails.

The two genuinely independent lineages already in the library — `kvlcc2_ugent`
(JSR-2017, Ghent) and `kvlcc2_jstage` (JASNAOE-2016) — both state `Lpp`/`B`/`T`/`∇`/`Cb`
but **do not print `Lwl`** (both are propeller/flow-focused papers whose particulars table
stops at the waterline length).  So the cell cannot be closed by re-reading an already
pinned document.

## 3. What a second lineage for `lwl` needs

A genuinely independent source for KVLCC2 `lwl` = 325.5 m must be a **non-SIMMAN**
publisher that states the waterline length.  Candidates, in order of likely egress success:

1. **KRISO/SNU original design paper** — the KVLCC2 (MOERI tanker) is documented in the
   KRISO flow-measurement papers (e.g. *Van, Kim, et al.*, "Flow measurement around a
   model ship with propeller", J. Soc. Naval Arch. Korea / SNAME).  A journal copy on a
   repository (KCI, J-STAGE, or an institutional mirror) would be a second lineage.
2. **SIMMAN 2008 original edition** — the 2008 workshop page is the *origin* of the
   SIMMAN-SERIES table; a 2008 edition is a different *edition* and would be recorded as a
   separate family (the lineage rule admits "demonstrably different editions").  The
   `simman2008.dk` domain is dead; an archived copy is the target.
3. **A CFD validation paper** that prints the KVLCC2 particulars table with `Lwl` — many
   open-access CFD papers (e.g. on MDPI, ScienceDirect OA) reproduce the full table.

Each candidate must be **fetched** (bytes pinned by SHA-256), **read** by
`parse_pdf_text`/`parse_source`, and **admitted** with its own `family`/`lineage` in
`tensorcad/fleet/crossref.py:SOURCES`.  A paper that merely reprints the SIMMAN table is
transcription corroboration and must be recorded `lineage=SIMMAN-SERIES` — the library's
design rule is that a second *paper* is not a second *source* unless it is a different
publisher family.

## 4. Honest scope of this round

This round **does not** close the `lwl` cell: no non-SIMMAN source for `lwl` was reachable
from this egress within the round's budget (the primary Todd/DTMB tables for Series 60 are
HTTP-403, and the KVLCC2 journal copies are not yet located).  What this round **does** is
pin the cell's state as a *documented, honest, single-lineage* gap with a named next-action
and a regression test that refuses to let a later edit reclassify the three SIMMAN
documents as two independent sources.

The regression test `tests/test_fleet_refs_kvlcc2_lwl_gap.py` asserts:
* KVLCC2 `lwl` has ≥2 *documents* (①b passes) but exactly 1 *lineage* (SIMMAN-SERIES);
* the three contributing sources all collapse to `lineage_of(...) == "SIMMAN-SERIES"`;
* the two genuinely independent lineages (`kvlcc2_ugent`, `kvlcc2_jstage`) do **not**
  contribute an `lwl` cell — so the gap is real, not a parse miss;
* criterion ① remains an honest FAIL for KVLCC2-`lwl`.

## 5. Next action for the next round

Hunt a **non-SIMMAN** journal copy of the KVLCC2 particulars table that prints `Lwl`
(priority: KRISO/SNU original paper on a repository mirror; then a CFD validation paper on
an OA host).  Fetch → pin → register as a new `family` → re-run
`tools/baseline_fleet/crossref_refs.py build` → the `lwl` cell flips to `CROSS_VERIFIED`
and criterion ① for KVLCC2 passes.
