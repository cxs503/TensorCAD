# Round 8 — arithmetic corroboration of ONRT's LCB and the residual-gap probe register

Anchor `REPRO-GEOM-FLEET-REFS`. Machine-readable copy:
`benchmarks/verified/fleet_refs/round8_lcb_corroboration.json`.

## 1. Why an arithmetic corroboration exists

The library counts a quantity `CROSS_VERIFIED` only when two independent **publishers**
state it. The IIHR thesis (`onrt_iihr`, sha256
`6fd4b42f0f0b873b…`) measures ONRT's LCB in metres at 1/49 scale and prints
`NA` in the full-scale column:

> `LCB LCB (m) aft of FP 1.625 NA`

Round 7 correctly refused to vote that cell. But refusing to *use* it throws away an
independent measurement, so this round converts it with the document's own stated scale and
the library's own Lpp — and keeps the result outside the census (`votes = false`).

| quantity | value | source |
|---|---|---|
| measured LCB (model, aft of FP) | 1.625 m | onrt_iihr |
| model scale | 1/49 | same document |
| **full-scale LCB (aft of FP), converted** | **79.625 m** | this tool's arithmetic |
| published LCB (SIMMAN-2020, fwd +) | -1.68 % Lpp | simman2020_summary |
| length used for the conversion | 154.0 m | simman2020_summary (`lwl`) |
| published LCB expressed aft of FP | 79.587 m | arithmetic |
| residual | +0.038 m (0.047 %) | |
| verdict | **CORROBORATES** (tolerance 0.5 %) | votes: False |

So the published percentage now has a second, *independent measurement* behind it that
agrees to 0.047 % — but the library verdict stays
`SINGLE_SOURCE`: a conversion is not a publisher, and a transcription slip in
the scale or in Lpp would otherwise be laundered into false agreement.

## 2. Where the other cells were looked for (negative results, recorded)

| document | targets | outcome |
|---|---|---|
| kvlcc2_ugent | KVLCC2/lcb_percent, KVLCC2/lwl, DTC/lcb_percent, DTC/lwl | no row states LCB or Lwl; the table's rows are Lpp, B, CB, TFP/TAP, volume displacement, LCG, VCG, GMT, rxx, ryy, rzz |
| jbc_memorial | JBC/lcb_percent | seven rows (LOA, B, LPP, LWL, draught, displacement, WSA); LCB is absent, so the single JBC LCB document is not corroborated |
| dtc_strath | DTC/lcb_percent, DTC/lwl | the text layer contains no LCB / longitudinal-centre-of-buoyancy / waterline-length row (zero matches); the document states Lpp, Bwl, T, volume, CB only |

The three documents are already pinned in the library, so each row is reproducible from the
digest it was read under; the Strathclyde row is a *zero-match* text search, not an OCR read.

Search provenance (never evidence, per the library's own rule):

| query | candidates | why not closed |
|---|---|---|
| Duisburg Test Case container ship (OpenAlex, open-access filter) | UGent 8516555 'open model test data' (item page, no PDF link answered); Strathclyde 67912 (already pinned as dtc_strath) | no candidate answers with a full text that prints both LCB and Lwl for the DTC |
| Japan Bulk Carrier JBC main particulars LCB block coefficient | Memorial 2025 thesis (already pinned), Sciendo 10.2478/pomr-2023-0002 (L/B ratio study, no JBC table hit) | the only JBC table reachable restates the NMRI rows and has no LCB |
| Series 60 Todd parent hull principal particulars LCB | Todd 1963 'Methodical experiments with models of single-screw merchant ships', Biodiversity Heritage Library item 102062 | the BHL item-PDF endpoint answered 403 to this egress and the catalogue PDF is a page scan without a text layer; Series 60 therefore still has zero documents and stays an open gap |

## 3. What this round does not claim

* No cell's verdict changes. `ONRT/lcb_percent` remains `SINGLE_SOURCE`; the four other
  open cells remain open.
* No numeral is transcribed by hand: every input is read from `crossref.json` or from the
  committed round-7 artifact.
* `Series60` still has zero documents. The 1963 Todd report is the right primary source and
  its BHL item-PDF endpoint refused this egress (403); that is recorded above so the next
  round starts from there rather than from scratch.
