# Is the missing reference data a retrieval gap or a reader gap? — reservoir sweep

Anchor **REPRO-GEOM-FLEET-REFS** (TensorCAD). Machine artefact
`benchmarks/verified/fleet_refs/reservoir_sweep.json` (regenerate:
`python tools/baseline_fleet/refs_reservoir_sweep.py --out benchmarks/verified/fleet_refs/reservoir_sweep.json
--markdown docs/fleet_refs_reservoir_tables.md --reservoir artifacts/fleet_refs/raw --reservoir /tmp/hunt
--reservoir /tmp/refs --reservoir /tmp/oa --reservoir /tmp/xref/raw`); code `tensorcad/fleet/reservoir.py`;
tests `tests/test_fleet_refs_reservoir.py`; generated tables `docs/fleet_refs_reservoir_tables.md`;
figure `docs/figs/fleet_refs_reservoir.png`. Every number below is rendered from that JSON.

## 1. The question this round answers

The fidelity round closed with a work list: five hulls still fail the anchor's criterion ①
("each hull: displacement / Cb / LCB / LWL from ≥ 2 independent sources"), and the report could
not say, for any single missing cell, whether the next worker should **fetch a new document** or
**parse a document that is already on disk**. Those are very different jobs, and the owner's
repetition rule (2026-11-09) makes the second kind the more expensive mistake: a new hunt over
ground a previous hunt already covered is work thrown away.

The reservoir answers it offline. Earlier rounds fetched far more bytes than they admitted — hunts
for a second JBC source, the DTC originator paper, Series 60 tables, KVLCC2 editions — and stopped
when the *admitted* set was still young. Those bytes are still in the scratch caches. This round
scans every one of them against exactly the cells the library is missing.

## 2. Method

Pure functions of (document bytes, library JSON) — `tensorcad/fleet/reservoir.py`:

1. **Roster.** Every PDF/HTML/text file under the given roots, SHA-256 pinned, empty files dropped.
   A document is marked `in_library` when its digest is one the library already cites (collected by
   walking the library JSON), so an admitted source is never mistaken for a new lead.
2. **Degeneration filters.** A `.txt` sidecar of a PDF is the *same* document: documents are
   collapsed by a whitespace-normalised text digest, and non-printable blobs (one 19 MB `.bin` in
   the cache) are recorded as *not text* instead of being scanned into gibberish. 21 candidate files
   were excluded this way — the exclusions are listed in the artefact, not silently dropped.
3. **Blocking cells.** Read from the library itself: `MISSING` cells plus every cell whose
   determinability class (from `fidelity.json`) is still blocking. 20 cells over five hulls.
4. **Mention scan.** For each reservoir document, the hulls it names and every *numeral-bearing*
   mention of a blocking quantity, quoted with its line. Three precision rules, each of which was
   added because the first run's output was unreadable: reference-list lines are skipped
   (`Proceedings|Vol.|pp.|doi|…`); numerals outside a wide plausibility window are rejected
   (`LCB 2021.8` is not an LCB); and LCB is also matched over a **two-line** window, because PDF
   extraction of a particulars table drops the numeral onto the next line — that one rule is the
   reason this round's only real lead was found at all.
5. **Verdict per cell.** `CANDIDATE_FIELD` (some reservoir document appears to state it — a parser
   job, no network) / `SOURCE_GAP_CONFIRMED` (the hull is in the reservoir, the quantity is nowhere
   stated) / `NOT_IN_RESERVOIR`.

**Nothing is admitted.** A candidate is a digest plus a quoted line; it never becomes a library
value, never changes a verdict, and cannot: criterion ② of the anchor needs a citation *and* a
digest, and a file whose publisher URL was not recorded when it was fetched has only half of that.

## 3. Results

See the generated tables (`docs/fleet_refs_reservoir_tables.md`) and the figure. The counts:

* **38 documents scanned** (digest-deduped, 8 of them already library sources, 30 outside it);
  21 files excluded as duplicates or non-text.
* **20 blocking cells**: 14 classify as `CANDIDATE_FIELD`, 6 as `SOURCE_GAP_CONFIRMED`,
  0 as `NOT_IN_RESERVOIR` — every hull the library owes is named by *some* document in the cache,
  so no cell is blocked by the reservoir being silent about the ship.
* **The one lead that matters.** JBC `LCB` — the only cell in the library carried by a *single*
  document (`SINGLE_LINEAGE_SINGLE_DOCUMENT`) — is stated by `d8991e.pdf`
  (sha256 `61c19a141a…`), "Numerical Simulation with Overset Approach for a Bulk Carrier (JBC) in
  Shallow Water", Table 1 particulars of the JBC model: `Lpp 4.8276 | Lwl 4.9138 | B 0.7759 |
  d 0.2845 | ∇ 0.9155 | Cb 0.8580 | LCB (%Lpp), fwd + 2.5475`. The digits are the Tokyo-2015 table's
  own (Cb 0.858, Cm 0.9981, ∇ scaled by λ³), i.e. a *transcription* witness, not a new lineage —
  admitting it would move the cell from `SINGLE_LINEAGE_SINGLE_DOCUMENT` to
  `SINGLE_LINEAGE_TRANSCRIPTION` (still blocking, but corroborated) and add nothing to the
  "≥ 2 independent sources" count.
* **DTC `LCB`, `LWL`** — `SOURCE_GAP_CONFIRMED`: three DTC documents are in the reservoir and none
  states either quantity. The previous round's claim ("no retrievable document prints them") is now
  measured over the whole fetched corpus, not just the admitted subset.
* **KVLCC2 `LCB`/`LWL`, ONRT `LCB`, the seven Series 60 cells** — `CANDIDATE_FIELD` for most of
  them, but on hand review (below) none of those leads is a genuine statement: the scan is
  document-level, so a JBC paper that cites KVLCC2 in its references contributes "hits" to KVLCC2
  cells.

### 3.1 Hand review of the 14 candidate cells (this is a review, not a measurement)

| cell | lead | review outcome |
|---|---|---|
| JBC `lcb_percent` | `d8991e.pdf`, quoted above | **genuine** — the numeral is in the JBC model table itself |
| JBC `cb`,`draft`,`displacement_volume`,`length_pp`,`lwl` | `jbc_maasch.pdf`, `jbc_researchsquare.pdf`, `jbc_semarak.pdf` | genuine *rows* but redundant: these documents are already represented for those cells or descend from the same Tokyo-2015 table; no verdict moves |
| KVLCC2 `lcb_percent`, ONRT `lcb_percent` | same JBC/ONRT papers | **rejected** — the numeral belongs to another hull's table (cross-hull attribution, limitation below) |
| Series 60 (5 cells) | `s60_*` papers and bibliography survivors | **rejected** — the Series 60 papers *set* `Lpp = 150 m` for their own computations or cite "Series 60 CB=0.60 model" in a reference list; none prints a principal-particulars table with ∇ and LCB |

So the honest bottom line: **the reservoir contains exactly one genuinely new field, and it is a
transcription witness for the one cell that was single-document. It does not close criterion ① for
any hull.** Series 60 remains `NO_DETERMINATION` on all seven quantities, which is a *retrieval*
gap: the primary tables (Todd, DTMB Report 1258) are behind HTTP 403 on this egress and the
reachable JASNAOE full texts do not print the table.

## 4. Why the lead cannot be admitted this round

`d8991e.pdf` sits in a scratch cache from an earlier hunt whose fetch log recorded neither the URL
nor the retrieval time. Admitting it would put a value in the library whose only provenance is a
local file — exactly the failure mode criterion ② exists to prevent. The artefact therefore keeps
it as `CANDIDATE_FIELD` with the digest and the quoted line, so the *first* thing a future worker
with working egress does is re-fetch this one document (title and digest both pinned) instead of
re-running a hunt.

## 5. Caveats and limitations

* **Document-level, not table-level.** A numeral is attributed to every hull the document names.
  That is the source of the KVLCC2/ONRT false leads above, and it is why the artefact labels leads
  "not admitted values" instead of "values".
* **Precision was not measured.** I reviewed the 14 candidate cells by hand; a precision figure
  would need a labelled sample. Marked *not_run*, not "high precision".
* **The reservoir is whatever the scratch caches held at scan time.** Caches are not part of the
  repository and can be garbage-collected; the *digests and quotes* in the artefact survive, the
  bytes do not. A cell classified `SOURCE_GAP_CONFIRMED` therefore means "absent from the 38
  documents listed here", not "absent from the literature".
* **Roots are explicit.** The scan only reads directories passed with `--reservoir`; it never
  wanders the filesystem, and it writes nothing outside `benchmarks/verified/` and `docs/`.
* The `.txt` sidecars in the caches were produced by an earlier round's extractor; where a PDF and
  its sidecar were both present the PDF was kept (first in sort order) and the sidecar recorded as a
  duplicate of it, so a document never votes twice.

## 6. Error judgement (anchor's numerical frame does not apply — stated explicitly)

This round moves **no numerical quantity**, so "is the key quantity within 3 % of a baseline?" has
no subject: there is no reconstruction, no reference value and no benchmark case here, only a
document census. The two checkable statements it does make are exact, not approximate:

| claim | how it is checked | result |
|---|---|---|
| every scanned document's identity is pinned | 38/38 documents carry a 64-hex SHA-256 in the artefact; test `test_committed_artefact_keeps_its_policy_and_digests` enforces it | 0 exceptions |
| the JBC `LCB` numeral quoted is the one in the document | quote printed verbatim from the extracted text (`fwd + 2.5475`), digest `61c19a141a…` | matches |
| no candidate changes a library verdict | the artefact's `policy` field and the absence of any write path from `reservoir.py` into `crossref.py`; test asserts the policy string | holds |
| DTC `LCB`/`LWL` truly absent | 3 DTC documents × 2 quantities, mention counts in the generated table | 0 statements |

Whether the *anchor* is closer to done: no. This round removes re-work (the next hunt now knows
that 37 of the 38 cached documents do not contain the missing cells) and converts "source gap" from
an assertion into a census, but criterion ① is exactly as unmet as it was in the fidelity round.