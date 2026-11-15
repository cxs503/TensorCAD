# What the six-hull reference library actually determines — transcription fidelity and determinability

Anchor **REPRO-GEOM-FLEET-REFS** (TensorCAD). Machine artefacts:
`benchmarks/verified/fleet_refs/fidelity.json` (regenerate:
`python tools/baseline_fleet/refs_fidelity.py --library benchmarks/verified/fleet_refs/crossref.json --raw artifacts/fleet_refs/raw`);
code `tensorcad/fleet/fidelity.py`; tests `tests/test_fleet_refs_fidelity.py`;
generated tables `docs/fleet_refs_fidelity_tables.md`; figure `docs/figs/fleet_refs_fidelity.png`.
Every number quoted below is rendered from that JSON — nothing is transcribed by hand.

## 1. The question this round answers

The library built in the previous rounds answers *"do two independent publishers state this
quantity, and do they agree?"* It reports 42 cells over six hulls, and the acceptance record is
honest about the outcome: **one hull (KCS) has all four core quantities cross-verified; five do
not.** Two rounds of open-access source hunting (13 recorded probes in
`tools/baseline_fleet/crossref_refs.py::PROBES`) already showed *why*: a benchmark case
definition is published once by its originator workshop, and every later document reprints that
table. Copying is not independence.

That left an unmeasured claim in the middle of the report: *how good is a "single-lineage" cell,
really?* Three documents repeating one table either agree digit for digit — in which case the
transcription channel is clean and the residual risk is "the originator was wrong", not "the
literature disagrees" — or they do not, in which case there are typos in the corpus that a
downstream reconstruction would be blamed for. This round measures that, and then classifies
every cell by what the retrievable corpus *determines*, so the gap list becomes work items
instead of a pass/fail.

## 2. Method

Pure functions of the library JSON plus the fetched bytes (`tensorcad/fleet/fidelity.py`):

1. **Transcription fidelity.** Per cell, over the entries allowed to vote (a
   `fallback_column` reading has no known scale and never votes): the distinct values, the
   spread `(max−min)/min`, a single-linkage modal cluster at the quantity's tolerance, and every
   entry outside that cluster listed with its printed `raw` text. One finding per
   `(source, value, raw)` — a document printing the same typo twice is one defect.
2. **Determinability.** Each cell is classified
   `NO_DETERMINATION` / `SINGLE_LINEAGE_SINGLE_DOCUMENT` / `SINGLE_LINEAGE_TRANSCRIPTION` /
   `MULTI_LINEAGE_AGREE` / `MULTI_LINEAGE_DISPUTED`. The first three still block the anchor's
   "≥ 2 independent sources".
3. **Reader-gap scan.** For each fetched document, count occurrences of each quantity's own
   tokens and compare with what the parser extracted from that document. A mention with no
   extraction is a *reader* gap (fixable here) or a prose mention (not a particular); a quantity
   no document mentions is a *source* gap (needs the network). 11 documents scanned.

## 3. Results

Tables are generated into `docs/fleet_refs_fidelity_tables.md`; the numbers that matter:

* 42 cells, **30 stated by ≥ 2 documents** — so the corpus is not thin, it is *vertically*
  thin: many copies, few lineages.
* **2 cells carry a transcription outlier**, and both are the defects the adjudication round
  found by hand: KVLCC2 `Bwl` printed as `580` in both SIMMAN editions (900 % off the modal
  `58.0`, a power-of-ten slip) and ONRT ∇ printed as `8507 m³` where NMRI prints a *mass*
  (2.44 %, the CONFLICT cell). The measurement rediscovers them independently — that is the
  test that the measurement is not a tautology.
* Excluding those two, the widest disagreement in the whole corpus is **0.077 %** (KCS Cb over
  three documents), i.e. **one to two orders below the 0.5 % tolerance**. Transcription, where
  documents exist at all, is faithful.
* Determinability: `MULTI_LINEAGE_AGREE`/`DISPUTED` for 6 of 7 quantities on KCS, KVLCC2 and
  DTC and 5 of 7 on ONRT; JBC is `SINGLE_LINEAGE_TRANSCRIPTION` on 6 cells and
  `SINGLE_LINEAGE_SINGLE_DOCUMENT` on LCB; **Series 60 is `NO_DETERMINATION` on all 7** — the
  library has no Series 60 document at all.
* Reader gaps found: 4 (e.g. `jbc_memorial` mentions LCB three times but the parser read no LCB
  from it; `ittc_benchmarks` mentions Cb four times in prose). These are *not* new sources — they
  are cases where a document may state more than the library currently records, so they are the
  cheapest next increment.

## 4. What this says about the anchor's acceptance criterion ①

Criterion ① ("each hull: displacement/Cb/LCB/LWL from ≥ 2 independent sources") is **not met**
and, on this egress, cannot be met for five hulls. The evidence is now per cell rather than
per hull:

* **DTC `LCB`, `LWL`** — no retrievable document prints them. The three DTC documents the
  library holds (`dtc_ije`, `dtc_strath`, `dtc_ugent`) were scanned for the tokens `LCB`/`Lwl`:
  **zero occurrences in all three**. This is a source gap, not a parser gap.
* **Series 60** — the primary tables (Todd, DTMB Report 1258) return HTTP 403 from the
  Biodiversity Heritage Library mirror, DTIC does not resolve, ScienceDirect blocks; the
  JASNAOE full texts reachable from this host (`s60_1996`, `s60_1975`, `s60_2002`) state no
  principal-particulars table with LCB.
* **JBC** — all four documents descend from the Tokyo 2015 workshop table (the three later ones
  declare it), so the lineage count is 1 by construction, not by failure of the search.
* **KVLCC2 `LCB`, `LWL`** and **ONRT `LCB`** — one lineage each (SIMMAN); the NMRI originator
  page for KVLCC2 (`t2015.nmri.go.jp/kvlcc2.html`) is HTTP 404 today, recorded as a probe.

Honest consequence, stated rather than papered over: **of the corpus' 42 cells only KCS' core
four are backed by two independent lineages**; what is now demonstrated is that the *transcription* channel is clean
to ≤ 0.077 % and that the corpus contains exactly two typographic defects, both already
adjudicated. A single-lineage cell in this library therefore means "the originator's table was
copied faithfully", which is a weaker claim than "two publishers agree" and is now labelled as
such on every cell.

Next-round suggestions (each is a single, checkable piece of work): (i) close the 4 reader gaps
above by extending the per-document `quotes`, which is a parser change with no network
dependency; (ii) for DTC/Series 60, treat an *independent determination from the originator's
released geometry* (pinned SHA in the asset library) as a documented second channel and label it
`GEOM-INTEGRATION` — never as literature; (iii) if the anchor's wording is to stand unchanged,
record in the ledger that criterion ① is unachievable for 5/6 hulls and cite this document.

## 5. Caveats

* The scan depends on which documents are present in `artifacts/fleet_refs/raw` (git-ignored,
  re-fetched on demand). 11 were present for this run; a document that could not be re-fetched
  is reported as missing rather than as clean.
* "Mentions without extraction" is deliberately over-inclusive: prose mentions count. Each of the
  4 rows needs a look at the text before it is called a defect; none is claimed as one here.
* The spread statistic is computed on full-scale values only; a model-scale column is excluded
  by the same rule that denies it a vote.