# Fleet reference library: re-binding every digit to its own bytes

Anchor: `REPRO-GEOM-FLEET-REFS` (six benchmark hulls, ≥2 independent sources per quantity).
Tool: `tools/baseline_fleet/digest_rebind.py`; module: `tensorcad/fleet/rebind.py`.
Artifact: `benchmarks/verified/fleet_refs/rebind_audit.json`.
Tests: `tests/test_fleet_refs_rebind.py` (10 cases, invariants not smoke).

## 1. why this audit exists

The library (`benchmarks/verified/fleet_refs/crossref.json`) prints, for every published
particular, a source key, a URL and the SHA-256 of the page or PDF the digit was read from.
That is the anchor's acceptance clause ② ("every value carries a citation plus a cached
digest"). Two half-checks existed and they never met:

| check | what it verified | what it could not see |
|---|---|---|
| `pin_audit` → `digest_audit.json` | 7 round-one pages: pinned digest == fetched bytes | whether the numeral the library carries is in those bytes |
| `oa_ingest.json` pins | 10 open-access PDFs digested at ingest | nothing ever re-read them: a cell whose *arithmetic* changed under a later code revision stayed published |

A library can therefore be all-green on digests and still print a number its own document does
not support. `rebind` closes the chain **bytes → sha256 → text → numeral → library value** for
every cell, document by document, offline.

Reading rules are the *current* ones (`parse_source` for HTML tables, `parse_pdf_text` for the
PDF quotes) — deliberately, because the question is "does this digit follow from these bytes
today", not "did it once".

## 2. method

1. `SOURCES` (17 documents: 7 round-one pages/reports + 10 open-access PDFs) each get their
   cached bytes from `/tmp/xref/raw` (HTML/PDF) or `/tmp/oa` (OA PDFs). A document with no
   cached copy is `NOT_CACHED` — **not** a pass.
2. Digest check with a *stability mask* (see §3.1).
3. Cells are re-derived from the bytes with the current rules and compared cell by cell.
4. Cell verdicts:
   * `MATCH` — the digit comes back out of its document;
   * `TRANSCRIPTION_DRIFT` — the library misquotes a printed numeral (never auto-edited: a
     repair would dress a new reading as an old one);
   * `STALE_DERIVED` — the cell is arithmetic on two numerals of the same document and the
     arithmetic changed (`unit_basis` starts with `derived`); safe to recompute from the same
     bytes, `repair` mode does exactly that and nothing else;
   * `NOT_IN_BYTES` — the library carries a cell the document no longer yields, or the document
     states a quantity the library lacks (a *shrunk* library becoming visible).
5. `repair` rewrites stale derived cells, re-derives the cached `closure`/`coverage`/
   `adjudication` blocks (`refresh_derived`) and drops the derived-vs-voter arithmetic the
   repair invalidated.

## 3. results

```
documents 17/17 read   cells rebound 105/107
documents_by_status  DIGEST_MATCH 17
cell statuses        MATCH 105, NOT_IN_BYTES 2
repairs              ONRT lcb_percent  -1.704545 -> -1.6365  (derived, ratio 49 -> 48.9355)
```

### 3.1 Finding A — two pinned digests can never be stable (page carries a request counter)

The SIMMAN-2014 KVLCC2 and KCS pages differ between fetches by exactly one line:

```
- ML:768-PU:145-PLA:19
+ ML:768-PU:146-PLA:19
```

66 359 bytes on both sides, same length, one counter byte-range. Pinning the *raw* bytes of
such a page guarantees a red freshness check on every future re-fetch, which is how a reader
learns to ignore drift. `rebind.stable_digest` masks publisher request counters declared in
`VOLATILE_PATTERNS` and reports which masks fired, so "the bytes moved" and "the document
changed" stop being the same sentence:

* `DIGEST_MATCH` — bytes identical to the pin;
* `DIGEST_VOLATILE_ONLY` — raw digest moved, masked digest identical (counter only);
* `DIGEST_DRIFT` — the content moved; the document must be re-read deliberately.

The mask is anchored on the publisher's own widget (`ML:<n>-PU:<n>-PLA:<n>`), never on a
generic number pattern, so it cannot eat a table numeral. Two tests pin both directions: the
counter-only case reports `VOLATILE_ONLY`, changing `3.48` to `3.51` still reports `DRIFT`.

### 3.2 Finding B — the published ONRT LCB was stale arithmetic (repaired)

`crossref.json` carried, for ONRT / LCB (% Lpp) / NMRI-T2015, the derived value **−1.704545**
("LCB 1.625 × 49 m … (74.375 − 77.000)/154 = −1.7045"). The document's own length ratio is
`Lwl 154 m / 3.147 m = 48.9355`, so the current reading rules give **−1.6365**, and
`derived_cells.json` — written in a later round — already said −1.6365. Two published
artifacts of the same anchor therefore disagreed by 0.068 %Lpp (4.1 % of the derived value)
and nothing detected it, because the cached blocks are recomputed from `report` while the cell
itself lives inside `report`.

After `repair` the library, its `derived` block and the audit agree: −1.6365, one basis, one
derivation string. The test `test_the_committed_library_is_internally_consistent` fails if a
later round reintroduces a derived ONRT LCB that its own document does not reproduce.

**Consequence for the anchor (honest, not smoothed over).** With the correct scale the two
independent routes to ONRT LCB disagree by 2.59 % relative (−1.68 from the SIMMAN 2020 summary
table vs −1.6365 from NMRI-2015's "1.625 m aft of FP" plus the pinned scale), i.e. above both
the 2 % tolerance used for that quantity and the 0.5 % threshold in acceptance clause ③. The
derived route is a *conversion*, not an independent publication of LCB in %Lpp, so it does not
vote — but the disagreement is material and is now recorded in the artifact instead of
disappearing into a stale cell. Recommendation: keep ONRT LCB `SINGLE_SOURCE` and add
`discordant_derivation: -1.6365` rather than promoting the conversion to a second source.

### 3.3 Finding C — two derived cells are no longer reproducible (open)

`kvlcc2_ugent` (JSR 2017, Table 1) still yields the five Table-1 columns, but the current rules
no longer produce an `lcb_percent` from it, while `crossref.json` carries derived cells
`DTC −0.969296` and `KVLCC2 3.46875` sourced to that document. They are audit-trail cells (the
document prints LCB, not LCB in %Lpp), so they do not vote — but a published digit that its
own document does not reproduce under the current reading rules is exactly what an audit is for.
They are reported, not deleted: deleting them would silently shrink the evidence. Next round
should either restore the derivation with a recorded rule or move them to a `superseded` block.

## 4. limits

* The audit verifies **that the digit follows from the bytes**; it cannot verify that the
  publisher is right, that the two sources are truly independent (that is the lineage axis, in
  `tensorcad/fleet/crossref.py`), or that a reprinted table was not copied from the other source.
* The reference bytes are the ones this environment still holds (`/tmp/xref/raw`, `/tmp/oa`).
  They are not redistributed (no publisher states a licence) and they are not in the repository,
  so a reviewer without the cache sees `NOT_CACHED` — which is the honest answer, never a pass.
* `--fetch` distinguishes cached-vs-live and is **off** by default so the audit is reproducible
  offline. In this session a live re-fetch was performed by hand for five documents
  (`simman2020_summary`, `t2015_kcs`, `t2015_onrt`, `simman2014_kvlcc2`, `simman2014_kcs`):
  the first three still hash to their pins, the two SIMMAN-2014 pages differ only by the
  request counter (§3.1). A full 17-document `--fetch` sweep exceeded this environment's
  120 s shell budget, so the committed artifact is the offline one — the live axis is a
  spot-check, not a claim.
* DTC/JBC `lcb_percent` and Series 60 (all quantities) remain **open** gaps of the anchor: they
  are missing *documents*, and no arithmetic can invent one.