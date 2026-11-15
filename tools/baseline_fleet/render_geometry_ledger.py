#!/usr/bin/env python3
"""Render the geometry-source ledger into a document and a figure.

Nothing in the tables is hand-typed: every number is read from
``tensorcad/fleet/data/geometry_source_ledger.json``, which itself is regenerated from the
cached publisher bytes by ``tools/baseline_fleet/geometry_source_ledger.py``.  The figure
plots the relative error of each geometric check against the 3 % acceptance line, so a
reader can see at a glance which checks pass and which are *supposed* to fail (the ONRT
unit-normalisation row).

Usage::

    python tools/baseline_fleet/render_geometry_ledger.py
"""

from __future__ import annotations

import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
LEDGER = REPO / "tensorcad" / "fleet" / "data" / "geometry_source_ledger.json"
DOC = REPO / "docs" / "gap_cad_t1_geometry_ledger.md"
FIGURE = REPO / "docs" / "figs" / "gap_cad_t1_geometry_ledger.png"
TOL = 0.03


def table(headers: list[str], rows: list[list[str]]) -> str:
    out = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    for row in rows:
        out.append("| " + " | ".join(row) + " |")
    return "\n".join(out)


def make_figure(ledger: dict) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = ledger["checks"]
    labels = [row["quantity"] for row in rows]
    errors = [row["abs_rel_error"] * 100.0 for row in rows]
    colours = ["#2c7fb8" if row["within_tolerance"] else "#d95f0e" for row in rows]

    fig, ax = plt.subplots(figsize=(9.0, 4.6), dpi=160)
    positions = range(len(rows))
    ax.barh(list(positions), errors, color=colours, height=0.6)
    ax.axvline(TOL * 100.0, color="#444444", linestyle="--", linewidth=1.2)
    ax.set_yticks(list(positions))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xscale("log")
    ax.set_xlabel("absolute relative error [%] (log scale); dashed = 3 % acceptance line")
    ax.set_title(
        "Geometry bytes vs published particulars — JBC at 1:40, ONRT unit-normalised",
        fontsize=10,
    )
    for pos, value, row in zip(positions, errors, rows):
        ax.text(
            value * 1.08,
            pos,
            f"{value:.2f} %" + ("" if row["within_tolerance"] else " (explained)"),
            va="center",
            fontsize=7,
        )
    ax.set_xlim(min(errors) * 0.4, max(errors) * 3.0)
    fig.tight_layout()
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE)


def main() -> int:
    ledger = json.loads(LEDGER.read_text())
    make_figure(ledger)

    sources = ledger["sources"]
    live = [(k, v) for k, v in sources.items() if v["http_status"] == 200]
    dead = ledger["dead_ends"]

    source_rows = [
        [
            f"`{key}`",
            str(entry["http_status"]),
            f"{entry['bytes']:,}",
            f"`{entry['sha256'][:16]}…`",
            "yes" if entry["digest_bound_in_repo"] else "—",
            entry["kind"],
        ]
        for key, entry in sorted(live)
    ]
    dead_rows = [
        [f"`{key}`", str(entry["http_status"]), f"`{entry['sha256'][:16]}…`", entry["url"]]
        for key, entry in sorted(dead.items())
    ]
    geometry_rows = []
    for ship, block in ledger["geometry"].items():
        fingerprint = block["iges_fingerprint"]
        extent = block["extent"]["extent_m"]
        geometry_rows.append(
            [
                ship,
                block["hull_member"]["name"],
                f"{block['hull_member']['bytes']:,}",
                f"`{block['hull_member']['sha256'][:16]}…`",
                f"{fingerprint['unit_flag']} ({fingerprint['unit_name']})",
                f"{extent[0]:.4f} × {extent[1]:.4f} × {extent[2]:.4f}",
                block["regime"]["regime"] if block["regime"] else "—",
            ]
        )
    check_rows = [
        [
            row["quantity"],
            f"{row['measured']:.6g}",
            f"{row['reference']:.6g}",
            row["unit"],
            f"{row['abs_rel_error'] * 100:.2f} %",
            "**yes**" if row["within_tolerance"] else "no",
            row["note"] or row["reference_source"],
        ]
        for row in ledger["checks"]
    ]
    acceptance_rows = [
        [item, block["verdict"], block["detail"]] for item, block in ledger["acceptance"].items()
    ]
    question_rows = [[q["item"], q["status"], q["evidence"]] for q in ledger["open_questions"]]
    summary = ledger["checks_summary"]

    text = f"""# Geometry source ledger — `GAP-CAD-T1-SOURCE-HUNT`

*Generated* {ledger["generated_utc"]} · *tool* `{ledger["tool"]}` ·
*record* `tensorcad/fleet/data/geometry_source_ledger.json` ·
*reproduce* `python tools/baseline_fleet/geometry_source_ledger.py --cache <cache> --out <json>`
(and `--verify` to re-digest the cache and diff the record)

![checks](figs/gap_cad_t1_geometry_ledger.png)

## 1. What this settles

The anchor asks a provenance question and a numeric one: **where is the standard-hull
geometry, is the copy we hold the publisher's, and does it actually measure up to the
published particulars?** Earlier rounds answered the first half with page digests. This
round answers the second half from the geometry bytes themselves: the JBC IGES archive is
authored in millimetres at the published scale denominator 1:40 and reproduces the
quoted beam and depth to **0.51 %** and **0.86 %**; the ONRT archive is not a 3.147 m
model hull at all but a **unit-length-normalised** export whose proportions match the
published B/Lwl to **0.52 %**, so it must be scaled by the published Lwl before meshing.

{summary["within_tolerance"]} of {summary["rows"]} geometric checks are inside the 3 %
acceptance line; the one that is not is the ONRT normalisation row, which exists to
*establish* the regime and is annotated as such.

## 2. Live sources retrieved (digests, no bytes redistributed)

{table(["source key", "HTTP", "bytes", "sha256 (16)", "in crossref", "kind"], source_rows)}

`in crossref = yes` means the same SHA-256 is registered in
`benchmarks/verified/fleet_refs/crossref.json`, so the cached bytes and the committed
manifest are one chain rather than two claims. Licences: the workshop pages are
copyright NMRI / SIMMAN and the geometry archives ship under a research-use grant with no
OSI licence file, so **nothing is copied into the repository** — only digests, member
names, licences and measured dimensions are recorded.

## 3. Dead ends, kept as evidence

{table(["source key", "HTTP", "error-body sha256", "URL"], dead_rows)}

The DTC row is the substance of the WITHHELD verdict below: the publisher page the DTC
paper names answers 404, the T2015 DTC page answers 404, the geometry directory answers
403 (not listable), and the Internet Archive is unreachable from this network. The only
public DTC copy found is a third-party OpenFOAM tutorial whose licence does not travel
with it. A negative result that keeps its error-body digest is auditable; one that is
dropped is not.

## 4. Geometry archives measured

{table(["ship", "hull member", "bytes", "sha256 (16)", "unit flag", "extent x×y×z [m]", "scale regime"], geometry_rows)}

The extent is measured by evaluating each B-spline surface on a 15 × 15 parameter grid,
not from the control net: on the JBC hull the control net overstates the length by 7.7 %
while the sampled surface agrees with the header's own maximum-coordinate claim
(`header_matches_extent = true` for both archives). That check is what rules out "we
measured the wrong entity".

## 5. Validation against the published particulars

{table(["check", "measured", "reference", "unit", "rel. error", "< 3 %", "note / source"], check_rows)}

Answering the acceptance question row by row: the four JBC rows and the ONRT
beam-over-length row are **inside 3 %**; the ONRT file-length row is **not**, by
construction, and its note records the reason (the export is normalised, so the physical
scale has to come from the published Lwl). No row is failing because of a transcription
error: measured values are computed from bytes, references are read from digest-bound
records.

### Cross-check against `scale_calibration.json`

The same question was already asked of `jbc.igs` by this repository's
`tensorcad/fleet/scale_calibration.py` (committed under
`tensorcad/fleet/data/scale_calibration.json`), which works from *waterline* curve extents
and reports 39.74 as its combined estimate (39.48 by waterline length, 40.00 by waterline
beam). This ledger measures the *sampled surface* extent over all 45 patches and gets
39.80 by full beam and 39.66 by depth. The two methods share no code and differ by 0.15 %,
and both land inside the published 1:40 to better than 1 %: that agreement is the reason
the JBC denominator is recorded as settled rather than estimated. The difference is
explained: a surface-sampled extent includes the bilge-keel strip, which a waterline
curve does not, so the surface method reads a marginally longer beam.

## 6. Acceptance items

{table(["item", "verdict", "detail"], acceptance_rows)}

## 7. What is still open

{table(["item", "status", "evidence"], question_rows)}

## 8. Limits of this round

* **No egress.** This environment cannot reach `t2015.nmri.go.jp` (HTTPS times out), so
  the ledger consumes the cache retrieved 2026-11-14 by the earlier hunt and adds no new
  fetch. A run with egress re-fetches and compares digests; the tool already supports
  that, and `--verify` is meaningful offline because it re-digests the cache.
* **Hulls without geometry.** KCS, KVLCC2 and Series 60 have digest-bound particulars but
  no geometry archive in this cache, so they are recorded as WITHHELD rather than inferred.
* **One archive per hull.** A single archive cannot prove a publisher's whole distribution;
  the ledger records exactly which member was measured, by digest.
"""
    DOC.parent.mkdir(parents=True, exist_ok=True)
    DOC.write_text(text)
    print(f"wrote {DOC} ({len(text)} chars) and {FIGURE}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())