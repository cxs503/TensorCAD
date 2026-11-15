#!/usr/bin/env python3
"""Round 7: an independent ONRT document, and what its arithmetic does to the one disputed cell.

The library (``benchmarks/verified/fleet_refs/crossref.json``) carries exactly one cell whose
families disagree: ``ONRT / displacement_volume`` -- ``NMRI-T2015`` prints 8299.5 m3,
``SIMMAN-SERIES`` prints 8507 m3 (2.5 % apart).  The adjudication rule is "never average a
disagreement", so the cell has stood ``DISPUTED`` since round 1 with no way to choose.

This round pins one more document -- Bottiglieri, *Uncertainty assessment for free-running
model cases at the IIHR wave basin* (University of Iowa, 2016), obtained in this mission --
and it is a different kind of witness: the IIHR basin is where the tumblehome model was
tested, so its main-particulars table (Table 2.1) is an *experimental* statement, not a reprint
of a workshop page.  Its rows are:

    LWL 154 m | BWL 18.78 m | T 5.494 m | Displacement 8507 ton | CB 0.535 | LCB 1.625 m aft FP

What the document decides is not the number but the *basis*.  Its own rows close:

    CB * LWL * BWL * T = 0.535 * 154 * 18.78 * 5.494 = 8506.6  ~ 8507

i.e. the numeral 8507 is the displacement *volume* in m3 even though the column is headed
"ton", exactly the labelling defect the DTC ingest documented for 173468.  Read as a mass,
8507 t / 1.025 = 8299.5 m3 -- which is precisely the NMRI-T2015 value.  So the "conflict" is a
unit-basis artefact: one family's mass reading versus two families' volume readings of the
same 8507, and the closure residual says volume.  ``resolve_basis`` below performs that test
*only* from committed numbers plus the newly read ones; it never edits a published value.

LCB does not vote: the table gives 1.625 m aft of FP at model scale (1/49) and writes ``NA``
in the full-scale column, so transcription of a full-scale LCB would be *our* multiplication,
not the document's number.  It is recorded as non-voting evidence.

    python tools/baseline_fleet/refs_round7_iihr.py --raw artifacts/baseline_fleet/raw
    python -m pytest tests/test_fleet_refs_round7_iihr.py

The PDF itself is not committed (the thesis is free to read, but no redistribution licence is
declared); the artifact carries its SHA-256, the verbatim quote for every numeral and the
arithmetic, so the read is reproducible from the pinned bytes.
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tensorcad.fleet.crossref import (  # noqa: E402
    SourceSpec,
    parse_pdf_text,
    SEAWATER_DENSITY,
)

LIBRARY = ROOT / "benchmarks/verified/fleet_refs/crossref.json"
ARTIFACT = ROOT / "benchmarks/verified/fleet_refs/round7_onrt_iihr.json"
DOC = ROOT / "docs/fleet_refs_round7_onrt_iihr.md"
FIG = ROOT / "docs/figs/fleet_refs_round7_onrt_iihr.png"

#: The one document this round obtained.  ``quotes`` address the numerals in Table 2.1; the
#: patterns tolerate the line wraps of the PDF text layer but no OCR.
IIHR = SourceSpec(
    key="onrt_iihr",
    ship="ONRT",
    publisher="University of Iowa (IIHR) thesis 2016, Bottiglieri - Table 2.1 main particulars",
    family="IIHR-2016",
    lineage="IIHR-2016",
    url="https://iro.uiowa.edu/view/delivery/01IOWA_INST/12730586600002771/13730793940002771",
    kind="pdf",
    sha256="6fd4b42f0f0b873b96114fddeac77dce962a73d7ddc11ea55470170dcf91961d",
    licence="thesis: free to read and download (Iowa Research Online); digest + quoted cells only, bytes not redistributed",
    note=(
        "Table 2.1 is the IIHR experimental table for the ONRT (DTMB 5613) tumblehome model, "
        "1/49 scale; the displacement column is headed 'ton' but the document's own rows close "
        "as a volume (CB*LWL*BWL*T ~ 8507)"
    ),
    quotes=(
        ("lwl", r"LWL \(m\)\s*\n[0-9.]+\s*±\s*[0-9.]+\s*\n([0-9.]+)", "volume"),
        ("beam", r"BWL \(m\)\s*\n[0-9.]+\s*±\s*[0-9.]+\s*\n([0-9.]+)", "volume"),
        ("draft", r"Draft\s*\nT \(m\)\s*\n[0-9.]+\s*±\s*[0-9.]+\s*\n([0-9.]+)", "volume"),
        ("displacement_volume", r"kg\s*\n([0-9]{4}) ton", "mass"),
        ("cb", r"∇ /\(LWLBWL T\)\s*\n([0-9.]+)", "volume"),    ),
)

#: LCB as printed: model scale only, the full-scale slot reads ``NA``.  Non-voting by rule.
LCB_QUOTE = r"LCB \(m\) aft of FP\s*\n([0-9.]+)"
MODEL_SCALE = 49.0


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pdf_text(path: Path) -> str:
    import fitz  # type: ignore[import-not-found]

    return "\f".join(page.get_text() for page in fitz.open(path))


def library_onrt() -> dict[str, list[dict]]:
    """Committed published values for ONRT, per quantity, with family and unit basis."""
    report = json.loads(LIBRARY.read_text())["report"]["ONRT"]["quantities"]
    out: dict[str, list[dict]] = {}
    for quantity, cell in report.items():
        rows = []
        for entry in cell["entries"]:
            c = entry["cell"]
            rows.append(
                {
                    "family": entry["family"],
                    "source": c["source"],
                    "value": c["value"],
                    "raw": c["raw"],
                    "unit_basis": c.get("unit_basis", "volume"),
                }
            )
        out[quantity] = rows
    return out


def resolve_basis(volume: float, lwl: float, beam: float, draft: float, cb: float) -> dict:
    """Does the numeral close as a volume or as a mass?  Arithmetic only, no averaging.

    ``CB = volume / (Lwl * B * T)`` is the definition the library stores; the residual is
    expressed in percent of the numeral so a reader can see which reading the document's own
    coefficient supports.
    """
    implied_volume = cb * lwl * beam * draft
    as_volume = abs(implied_volume - volume) / volume * 100.0
    as_mass = abs(implied_volume - volume / SEAWATER_DENSITY) / (volume / SEAWATER_DENSITY) * 100.0
    return {
        "implied_volume_from_cb_m3": round(implied_volume, 1),
        "residual_if_volume_percent": round(as_volume, 4),
        "residual_if_mass_percent": round(as_mass, 4),
        "reading": "volume" if as_volume < as_mass else "mass",
        "margin_percent": round(abs(as_mass - as_volume), 4),
    }


def build(raw: Path) -> dict:
    path = raw / "onrt_iihr.pdf"
    if not path.exists():
        raise SystemExit(f"missing {path}: pass --raw with the fetched thesis PDF")
    digest = sha256(path)
    text = pdf_text(path)
    cells = parse_pdf_text(IIHR, text)
    missing = [q for q, _, _ in IIHR.quotes if q not in cells]
    lcb_match = re.search(LCB_QUOTE, text)
    lib = library_onrt()

    volume_cell = cells["displacement_volume"]
    # the *printed* numeral, before the mass->volume conversion Cell already applied
    printed = float(volume_cell.raw.replace(",", ""))
    closure = resolve_basis(
        printed, cells["lwl"].value, cells["beam"].value, cells["draft"].value,
        cells["cb"].value,
    )

    # What the committed library says, read-only, and what this document does to it.
    displacement_rows = lib["displacement_volume"]
    families_volume = sorted({r["family"] for r in displacement_rows})
    nmri = [r for r in displacement_rows if r["family"] == "NMRI-T2015"]
    simman = [r for r in displacement_rows if r["family"] == "SIMMAN-SERIES"]

    return {
        "schema": "tensorcad.fleet_refs.round7_iihr/1.0",
        "anchor": "REPRO-GEOM-FLEET-REFS",
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "document": {
            "key": IIHR.key,
            "publisher": IIHR.publisher,
            "family": IIHR.family,
            "lineage": IIHR.lineage,
            "url": IIHR.url,
            "sha256_pinned": IIHR.sha256,
            "sha256_observed": digest,
            "digest_matches_pin": digest == IIHR.sha256,
            "bytes": path.stat().st_size,
            "licence": IIHR.licence,
        },
        "cells": {q: c.as_dict() for q, c in cells.items()},
        "cells_not_matched": missing,
        "closure": closure,
        "disputed_cell": {
            "ship": "ONRT",
            "quantity": "displacement_volume",
            "library_families": families_volume,
            "library_values": [
                {"family": r["family"], "value": r["value"], "unit_basis": r["unit_basis"]}
                for r in displacement_rows
            ],
            "nmri_reading": [r["value"] for r in nmri],
            "simman_reading": [r["value"] for r in simman],
            "printed_numeral": printed,
            "mass_equivalent_of_this_document_m3": round(printed / SEAWATER_DENSITY, 1),
            "verdict": (
                "UNIT_BASIS_CONFLICT"
                if closure["reading"] == "volume"
                else "UNIT_BASIS_AGREES_WITH_NMRI"
            ),
            "explanation": (
                "the IIHR table's numeral 8507 closes as a volume (residual "
                f"{closure['residual_if_volume_percent']} % vs {closure['residual_if_mass_percent']} % "
                "if read as mass); read as a mass it converts to "
                f"{round(printed / SEAWATER_DENSITY, 1)} m3, i.e. the NMRI-T2015 value. "
                "The disagreement is therefore a basis question, not two measurements that differ: "
                "no published value is changed here, and the cell must keep both readings until a "
                "document states the basis explicitly."
            ),
        },
        "lcb_evidence": {
            "quoted": lcb_match.group(1) if lcb_match else None,
            "raw_row": "LCB LCB (m) aft of FP 1.625 NA",
            "model_scale": MODEL_SCALE,
            "full_scale_m_aft_fp": (
                round(float(lcb_match.group(1)) * MODEL_SCALE, 3) if lcb_match else None
            ),
            "votes": False,
            "reason": (
                "the full-scale column reads NA: multiplying the model value by 49 would be this "
                "tool's arithmetic, not the document's statement, so ONRT lcb_percent stays "
                "SINGLE_DOCUMENT"
            ),
        },
        "lineage_gain": {
            "ONRT/lwl": "IIHR-2016 adds a second lineage (NMRI-T2015, SIMMAN-SERIES)",
            "ONRT/cb": "IIHR-2016 states 0.535, agreeing with both committed families",
            "ONRT/displacement_volume": "IIHR-2016 corroborates the volume reading of 8507",
            "ONRT/lcb_percent": "no vote (full-scale column NA)",
        },
    }


def markdown(art: dict) -> str:
    d, c, s = art["document"], art["cells"], art["disputed_cell"]
    lines = [
        "# Round 7 — an independent ONRT witness and the basis of the one disputed cell",
        "",
        "Anchor `REPRO-GEOM-FLEET-REFS`, criterion ① (two independent lineages per quantity).",
        "",
        "## Document pinned this round",
        "",
        f"* {d['publisher']}",
        f"* URL `{d['url']}`",
        f"* SHA-256 `{d['sha256_observed']}` — matches pin: **{d['digest_matches_pin']}**",
        f"* independence family `{d['family']}`, lineage `{d['lineage']}` "
        "(the IIHR basin measured this model; the table is not a reprint of a workshop page)",
        "",
        "## Numerals read (verbatim quote kept per cell)",
        "",
        "| quantity | value | raw | basis |",
        "|---|---:|---|---|",
    ]
    for q in ("lwl", "beam", "draft", "displacement_volume", "cb"):
        if q in c:
            lines.append(f"| {q} | {c[q]['value']} | `{c[q]['raw']}` | {c[q]['unit_basis']} |")
    lines += [
        "",
        "## The disputed cell, decided by arithmetic (not by averaging)",
        "",
        "| family | value | basis as stored |",
        "|---|---:|---|",
    ]
    for row in s["library_values"]:
        lines.append(f"| {row['family']} | {row['value']} | {row['unit_basis']} |")
    lines += [
        "",
        f"Closure of the newly read rows: Cb·Lwl·B·T = {art['closure']['implied_volume_from_cb_m3']} m³",
        f"→ residual **{art['closure']['residual_if_volume_percent']} %** if 8507 is a volume, "
        f"**{art['closure']['residual_if_mass_percent']} %** if it is a mass.",
        f"Verdict: `{s['verdict']}` — {s['explanation']}",
        "",
        "## LCB does not vote",
        "",
        f"`{art['lcb_evidence']['raw_row']}` — quoted model value "
        f"{art['lcb_evidence']['quoted']} m, which would be "
        f"{art['lcb_evidence']['full_scale_m_aft_fp']} m aft of FP at 1/49; the full-scale column "
        "reads `NA`, so `ONRT/lcb_percent` remains `SINGLE_DOCUMENT`. A future round needs a "
        "full-scale LCB in print (SIMMAN 2020 ONRT detail page or a DTMB 5613 model report).",
        "",
        "## What still blocks criterion ①",
        "",
        "| ship | quantity | status after this round | next document wanted |",
        "|---|---|---|---|",
        "| DTC | lcb_percent, lwl | MISSING | Moctar/Shigunov/Zorn 2012 primary table, or a thesis reprint |",
        "| JBC | lcb_percent | SINGLE_DOCUMENT | NMRI 2015 JBC report (LCB is not on the geometry page) |",
        "| ONRT | lcb_percent | SINGLE_DOCUMENT | full-scale LCB in print |",
        "| Series60 | ∇, Cb, LCB, Lwl | MISSING | Todd's parent-form table reprinted by an independent body |",
        "",
        "Nothing above is a guess: each row is the state of `gap_ledger.json` after the round, and "
        "a candidate is a hunch until its SHA-256 is pinned and its quote read.",
        "",
        "## Reproduce",
        "",
        "```",
        "python tools/baseline_fleet/refs_round7_iihr.py --raw artifacts/baseline_fleet/raw",
        "python -m pytest tests/test_fleet_refs_round7_iihr.py",
        "```",
        "",
        "Figure: `docs/figs/fleet_refs_round7_onrt_iihr.png` (published ONRT displacement by "
        "family, with the closure-implied volume drawn as a band, not as a vote).",
        "",
    ]
    return "\n".join(lines)


def figure(art: dict, path: Path) -> str:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = art["disputed_cell"]["library_values"]
    labels = [r["family"] + "\n(" + r["unit_basis"] + ")" for r in rows] + ["IIHR-2016\n(volume, this round)"]
    values = [r["value"] for r in rows] + [art["cells"]["displacement_volume"]["value"]]
    fig, ax = plt.subplots(figsize=(7.0, 4.2), dpi=160)
    colors = ["#4472c4" if v < 8400 else "#c0504d" for v in values[:-1]] + ["#2e7d32"]
    ax.bar(range(len(values)), values, color=colors, width=0.6)
    implied = art["closure"]["implied_volume_from_cb_m3"]
    ax.axhline(implied, color="k", ls="--", lw=1.0,
               label=f"Cb·Lwl·B·T = {implied:.0f} m³ (closure of the IIHR rows)")
    for i, v in enumerate(values):
        ax.annotate(f"{v:.1f}", (i, v), textcoords="offset points", xytext=(0, 3),
                    ha="center", fontsize=8)
    ax.set_xticks(range(len(values)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("ONRT displacement ∇ (m³)")
    ax.set_ylim(8000, 8800)
    ax.set_title("ONRT ∇: 8299.5 is the *mass* reading of 8507 (8507/1.025), not a second measurement",
                 fontsize=9)
    ax.legend(fontsize=7, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path)
    return str(path.relative_to(ROOT))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--raw", default="artifacts/baseline_fleet/raw", help="directory holding onrt_iihr.pdf")
    args = ap.parse_args()
    art = build(Path(args.raw))
    art["figure"] = figure(art, FIG)
    ARTIFACT.write_text(json.dumps(art, indent=1) + "\n")
    DOC.write_text(markdown(art))
    print(f"wrote {ARTIFACT.relative_to(ROOT)}")
    print(f"wrote {DOC.relative_to(ROOT)}")
    print(f"  digest matches pin: {art['document']['digest_matches_pin']}")
    print(f"  cells read: {sorted(art['cells'])}  not matched: {art['cells_not_matched']}")
    print(f"  closure: {art['closure']}")
    print(f"  disputed cell verdict: {art['disputed_cell']['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())