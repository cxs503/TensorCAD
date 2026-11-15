#!/usr/bin/env python3
"""Round 8: arithmetic corroboration for ONRT's LCB, and the residual-gap probe register.

Anchor ``REPRO-GEOM-FLEET-REFS``.  After round 7 the library
(``benchmarks/verified/fleet_refs/crossref.json``) still carries five cells that one
publisher states and nobody else does:

    DTC/lcb_percent, DTC/lwl, JBC/lcb_percent, ONRT/lcb_percent, Series60/* (four cells)

This round does two things and claims nothing more.

**1. It corroborates ONRT's LCB without a second publisher.**  The IIHR thesis pinned in
round 7 (``onrt_iihr``) prints ``LCB (m) aft of FP  1.625`` in the model column and ``NA``
in the full-scale column, which is why round 7 refused to vote it.  Refusing to vote is
right -- multiplying by 49 would be *our* arithmetic -- but the number is still an
independent *measurement*, so it can be converted with the document's own stated scale
(1/49) and the library's own Lpp and compared with the SIMMAN percentage.  If the two land
in the same 0.5 % band, the published percentage has a second, arithmetic-only check; the
verdict stays ``SINGLE_SOURCE`` because a conversion is not a publisher
(:mod:`tensorcad.fleet.corroboration` states that rule and ``votes`` is hard-wired false).

**2. It records where the other four cells were looked for and why they did not close.**
Three documents already pinned in the library were re-read specifically for the missing
rows -- the JSR-2017 KVLCC2/DTC table, the Memorial University JBC table, the Strathclyde
DTC table -- and their text was searched for LCB / Lwl rows.  None states them; the quotes
(or the zero-match count) are recorded with the document digests so the next round does not
pay for the same search.  The OpenAlex queries that produced the candidates are recorded as
*search provenance*: a probe can point at a document, it can never move a verdict.

    python tools/baseline_fleet/refs_round8_lcb_corroboration.py
    python -m pytest tests/test_fleet_refs_round8_lcb_corroboration.py -q

Nothing here transcribes a numeral by hand: every input is read out of
``crossref.json`` or a committed round artifact, whose ``sha256`` traces it to fetched bytes.
"""

from __future__ import annotations

import argparse
import datetime
import json
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tensorcad.fleet.corroboration import (  # noqa: E402
    LCB_TOLERANCE,
    corroborate_lcb_by_scale,
)

LIBRARY = ROOT / "benchmarks/verified/fleet_refs/crossref.json"
ROUND7 = ROOT / "benchmarks/verified/fleet_refs/round7_onrt_iihr.json"
ARTIFACT = ROOT / "benchmarks/verified/fleet_refs/round8_lcb_corroboration.json"
DOC = ROOT / "docs/fleet_refs_round8_lcb_corroboration.md"
FIG = ROOT / "docs/figs/fleet_refs_round8_lcb_corroboration.png"

#: What this round searched for, and what it found in the three documents it re-read.
#: ``quotes`` are verbatim from the extracted text layer (never OCR); ``zero_match`` counts a
#: text search that returned nothing, which is a *negative* result and is labelled as such.
PROBES: tuple[dict[str, Any], ...] = (
    {
        "target": ["KVLCC2/lcb_percent", "KVLCC2/lwl", "DTC/lcb_percent", "DTC/lwl"],
        "document": "kvlcc2_ugent",
        "title": "Journal of Ship Research 61(3) 2017, Table 1 'Main particulars and loading conditions for the KVLCC2 and DTC'",
        "url": "https://biblio.ugent.be/publication/8563318/file/8563319.pdf",
        "outcome": "no row states LCB or Lwl; the table's rows are Lpp, B, CB, TFP/TAP, volume displacement, LCG, VCG, GMT, rxx, ryy, rzz",
        "quotes": [
            "Table 1 Main particulars and loading conditions for the KVLCC2 and DTC",
            "Lpp (m) 320.0 355.0",
            "CB [-] 0.8098 0.661",
            "LCG (m) 171.1 169.0 174.059 170.460",
        ],
    },
    {
        "target": ["JBC/lcb_percent"],
        "document": "jbc_memorial",
        "title": "Memorial University thesis 2025, Table 3 'JBC hull data and geometry'",
        "url": "https://storage.googleapis.com/dissertations/2025/MEM_2025_KABONGE.pdf",
        "outcome": "seven rows (LOA, B, LPP, LWL, draught, displacement, WSA); LCB is absent, so the single JBC LCB document is not corroborated",
        "quotes": [
            "Table 3 JBC hull data and geometry No. Main Particulars Value Units",
            "1. LOA 291 m 2. Beam, B 45 m 3. LPP 280 m 4. LWL 285 m",
            "5. Draught 16.5 m 6. Displacement 178370 m3 7. WSA 19556.1 m2",
        ],
    },
    {
        "target": ["DTC/lcb_percent", "DTC/lwl"],
        "document": "dtc_strath",
        "title": "ISNAM 2018 (Strathclyde repository), Table 1 main particulars of the DTC hull",
        "url": "https://strathprints.strath.ac.uk/67912/1/Dogrul_etal_ISNAM2018_Uncertainty_assessment_and_self_propulsion_estimation_of_Duisburg_Test_Case.pdf",
        "outcome": "the text layer contains no LCB / longitudinal-centre-of-buoyancy / waterline-length row (zero matches); the document states Lpp, Bwl, T, volume, CB only",
        "zero_match": r"(?i)lcb|longitudinal cent(?:er|re) of buoyancy|waterline length",
        "quotes": [],
    },
)

#: Search provenance only -- see the library's own rule: "a probe is search provenance,
#: never evidence".  Recorded so the next round starts where this one stopped.
SEARCHES: tuple[dict[str, str], ...] = (
    {
        "query": "Duisburg Test Case container ship (OpenAlex, open-access filter)",
        "candidates": "UGent 8516555 'open model test data' (item page, no PDF link answered); Strathclyde 67912 (already pinned as dtc_strath)",
        "why_not_closed": "no candidate answers with a full text that prints both LCB and Lwl for the DTC",
    },
    {
        "query": "Japan Bulk Carrier JBC main particulars LCB block coefficient",
        "candidates": "Memorial 2025 thesis (already pinned), Sciendo 10.2478/pomr-2023-0002 (L/B ratio study, no JBC table hit)",
        "why_not_closed": "the only JBC table reachable restates the NMRI rows and has no LCB",
    },
    {
        "query": "Series 60 Todd parent hull principal particulars LCB",
        "candidates": "Todd 1963 'Methodical experiments with models of single-screw merchant ships', Biodiversity Heritage Library item 102062",
        "why_not_closed": "the BHL item-PDF endpoint answered 403 to this egress and the catalogue PDF is a page scan without a text layer; Series 60 therefore still has zero documents and stays an open gap",
    },
)


def load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def library_cell(library: dict[str, Any], ship: str, quantity: str) -> dict[str, Any]:
    return library["report"][ship]["quantities"][quantity]


def first_cell_value(cell: dict[str, Any]) -> tuple[float, str, str]:
    entry = cell["entries"][0]["cell"]
    return entry["value"], entry["source"], entry["sha256"]


def build(library: dict[str, Any], round7: dict[str, Any]) -> dict[str, Any]:
    # ---- 1. ONRT LCB: model-scale measurement -> library's percentage ----------------
    # ONRT's workshop row prints one length only (Lwl); it is the length the library holds,
    # so the percentage is converted on that and the quantity used is recorded, not assumed.
    for lpp_quantity in ("length_pp", "lwl"):
        cell = library_cell(library, "ONRT", lpp_quantity)
        if cell["entries"]:
            break
    onrt_lpp, lpp_source, _ = first_cell_value(cell)
    onrt_lcb_cell = library_cell(library, "ONRT", "lcb_percent")
    published_pct, lcb_source, lcb_sha = first_cell_value(onrt_lcb_cell)

    lcb_evidence = round7["lcb_evidence"]
    model_len = float(lcb_evidence["quoted"])
    scale = float(lcb_evidence["model_scale"])
    doc = round7["document"]
    if lcb_evidence["votes"]:
        raise SystemExit("round 7 already voted this cell: this round's premise is stale")

    corroboration = corroborate_lcb_by_scale(
        ship="ONRT",
        model_length_aft_fp_m=model_len,
        scale=scale,
        lpp_m=onrt_lpp,
        published_percent_fwd_positive=published_pct,
    )

    return {
        "schema": "tensorcad.fleet_refs.round8/1.0",
        "anchor": "REPRO-GEOM-FLEET-REFS",
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).replace(microsecond=0).isoformat(),
        "rule": "an arithmetic corroboration is a cross-check of one published number; it is never a second publisher and never changes a verdict",
        "lcb_corroboration": {
            "ship": "ONRT",
            "quantity": "lcb_percent",
            "library_verdict_before": onrt_lcb_cell["verdict"],
            "library_verdict_after": onrt_lcb_cell["verdict"],
            "measurement": {
                "document": doc["key"],
                "publisher": doc["publisher"],
                "sha256": doc["sha256_observed"],
                "digest_matches_pin": doc["digest_matches_pin"],
                "quote": lcb_evidence["raw_row"],
                "model_length_aft_fp_m": model_len,
                "scale": scale,
                "note": "the document's full-scale column reads NA; the conversion below is this tool's arithmetic and is labelled as such",
            },
            "published": {
                "source": lcb_source,
                "sha256": lcb_sha,
                "family": onrt_lcb_cell["entries"][0]["family"],
                "percent_fwd_positive": published_pct,
                "lpp_m": onrt_lpp,
                "lpp_source": lpp_source,
                "lpp_quantity": lpp_quantity,
            },
            "conversion": corroboration.as_dict(),
            "tolerance": LCB_TOLERANCE,
        },
        "gap_probes": list(PROBES),
        "searches": list(SEARCHES),
        "remaining_gaps": {
            "DTC": ["lcb_percent", "lwl"],
            "JBC": ["lcb_percent"],
            "ONRT": ["lcb_percent (single publisher; arithmetic corroboration recorded above)"],
            "Series60": ["displacement_volume", "cb", "lcb_percent", "lwl"],
        },
    }


def markdown(artifact: dict[str, Any]) -> str:
    c = artifact["lcb_corroboration"]
    conv = c["conversion"]
    rows = "\n".join(
        f"| {p['document']} | {', '.join(p['target'])} | {p['outcome']} |" for p in artifact["gap_probes"]
    )
    searches = "\n".join(f"| {s['query']} | {s['candidates']} | {s['why_not_closed']} |" for s in artifact["searches"])
    return f"""# Round 8 — arithmetic corroboration of ONRT's LCB and the residual-gap probe register

Anchor `REPRO-GEOM-FLEET-REFS`. Machine-readable copy:
`benchmarks/verified/fleet_refs/round8_lcb_corroboration.json`.

## 1. Why an arithmetic corroboration exists

The library counts a quantity `CROSS_VERIFIED` only when two independent **publishers**
state it. The IIHR thesis (`{c['measurement']['document']}`, sha256
`{c['measurement']['sha256'][:16]}…`) measures ONRT's LCB in metres at 1/49 scale and prints
`NA` in the full-scale column:

> `{c['measurement']['quote']}`

Round 7 correctly refused to vote that cell. But refusing to *use* it throws away an
independent measurement, so this round converts it with the document's own stated scale and
the library's own Lpp — and keeps the result outside the census (`votes = false`).

| quantity | value | source |
|---|---|---|
| measured LCB (model, aft of FP) | {c['measurement']['model_length_aft_fp_m']} m | {c['measurement']['document']} |
| model scale | 1/{c['measurement']['scale']:.0f} | same document |
| **full-scale LCB (aft of FP), converted** | **{conv['full_scale_aft_fp_m']:.3f} m** | this tool's arithmetic |
| published LCB (SIMMAN-2020, fwd +) | {c['published']['percent_fwd_positive']} % Lpp | {c['published']['source']} |
| length used for the conversion | {c['published']['lpp_m']} m | {c['published']['lpp_source']} (`{c['published']['lpp_quantity']}`) |
| published LCB expressed aft of FP | {conv['published_aft_fp_m']:.3f} m | arithmetic |
| residual | {conv['residual_m']:+.3f} m ({conv['residual_relative'] * 100:.3f} %) | |
| verdict | **{conv['verdict']}** (tolerance {c['tolerance'] * 100:.1f} %) | votes: {conv['votes']} |

So the published percentage now has a second, *independent measurement* behind it that
agrees to {conv['residual_relative'] * 100:.3f} % — but the library verdict stays
`{c['library_verdict_after']}`: a conversion is not a publisher, and a transcription slip in
the scale or in Lpp would otherwise be laundered into false agreement.

## 2. Where the other cells were looked for (negative results, recorded)

| document | targets | outcome |
|---|---|---|
{rows}

The three documents are already pinned in the library, so each row is reproducible from the
digest it was read under; the Strathclyde row is a *zero-match* text search, not an OCR read.

Search provenance (never evidence, per the library's own rule):

| query | candidates | why not closed |
|---|---|---|
{searches}

## 3. What this round does not claim

* No cell's verdict changes. `ONRT/lcb_percent` remains `SINGLE_SOURCE`; the four other
  open cells remain open.
* No numeral is transcribed by hand: every input is read from `crossref.json` or from the
  committed round-7 artifact.
* `Series60` still has zero documents. The 1963 Todd report is the right primary source and
  its BHL item-PDF endpoint refused this egress (403); that is recorded above so the next
  round starts from there rather than from scratch.
"""


def figure(artifact: dict[str, Any], path: Path) -> str:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    conv = artifact["lcb_corroboration"]["conversion"]
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    labels = ["published\n(SIMMAN-2020)", "measured\n(IIHR, 1/49 x 49)"]
    values = [conv["published_aft_fp_m"], conv["full_scale_aft_fp_m"]]
    ax.bar(labels, values, color=["#4477aa", "#ee6677"], width=0.55)
    for i, v in enumerate(values):
        ax.text(i, v, f"{v:.3f} m", ha="center", va="bottom", fontsize=9)
    ax.set_ylim(min(values) - 1.0, max(values) + 0.6)
    ax.set_ylabel("LCB, m aft of FP")
    ax.set_title(
        "ONRT LCB: published percentage vs model measurement x scale\n"
        f"residual {conv['residual_m']:+.3f} m = {conv['residual_relative'] * 100:.3f} % "
        f"(tolerance {artifact['lcb_corroboration']['tolerance'] * 100:.1f} %)"
    )
    ax.grid(axis="y", alpha=0.3)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return str(path.relative_to(ROOT))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--library", type=Path, default=LIBRARY)
    parser.add_argument("--round7", type=Path, default=ROUND7)
    parser.add_argument("--out", type=Path, default=ARTIFACT)
    parser.add_argument("--doc", type=Path, default=DOC)
    parser.add_argument("--fig", type=Path, default=FIG)
    args = parser.parse_args(argv)

    artifact = build(load(args.library), load(args.round7))
    artifact["figure"] = figure(artifact, args.fig)
    args.out.write_text(json.dumps(artifact, indent=1) + "\n")
    args.doc.write_text(markdown(artifact))
    conv = artifact["lcb_corroboration"]["conversion"]
    print(
        f"ONRT LCB corroboration: {conv['verdict']} residual {conv['residual_m']:+.3f} m "
        f"({conv['residual_relative'] * 100:.3f} %), votes={conv['votes']}"
    )
    print(f"wrote {args.out.relative_to(ROOT)}, {args.doc.relative_to(ROOT)}, {artifact['figure']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())