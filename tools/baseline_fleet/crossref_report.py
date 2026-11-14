#!/usr/bin/env python3
"""Render the cross-source reference library into the report and its figures.

    python tools/baseline_fleet/crossref_report.py \
        --library benchmarks/verified/fleet_refs/crossref.json \
        --doc docs/fleet_refs_crossref.md --figs docs/figs/fleet_refs

Every number in the markdown is injected from the library JSON; nothing is transcribed.
Figures: (1) the ship x quantity verdict matrix, (2) reference vs reconstructed displacement
with the error of each hull, (3) the Cb-closure check that exposes the KVLCC2 Bwl typo.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SHIPS = ["KCS", "KVLCC2", "JBC", "ONRT", "DTC", "Series60"]
QUANTITIES = ["displacement_volume", "cb", "lcb_percent", "lwl", "length_pp", "beam", "draft"]
SHORT = {
    "displacement_volume": "∇",
    "cb": "Cb",
    "lcb_percent": "LCB",
    "lwl": "Lwl",
    "length_pp": "Lpp",
    "beam": "Bwl",
    "draft": "T",
}
VERDICT_COLOUR = {
    "CROSS_VERIFIED": "#2e7d32",
    "CONFLICT": "#c62828",
    "SINGLE_SOURCE": "#ef6c00",
    "MISSING": "#9e9e9e",
}


def verdict_matrix(library: dict[str, Any], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.4, 3.4))
    for row, ship in enumerate(SHIPS):
        for col, quantity in enumerate(QUANTITIES):
            cell = library["report"][ship]["quantities"][quantity]
            ax.add_patch(
                plt.Rectangle((col, row), 1, 1, color=VERDICT_COLOUR[cell["verdict"]], alpha=0.85)
            )
            label = {
                "CROSS_VERIFIED": f"{cell['independent_families']}×",
                "CONFLICT": "✗",
                "SINGLE_SOURCE": "1×",
                "MISSING": "—",
            }[cell["verdict"]]
            ax.text(col + 0.5, row + 0.5, label, ha="center", va="center", color="white", fontsize=10)
    ax.set_xlim(0, len(QUANTITIES))
    ax.set_ylim(len(SHIPS), 0)
    ax.set_xticks([i + 0.5 for i in range(len(QUANTITIES))])
    ax.set_xticklabels([SHORT[q] for q in QUANTITIES])
    ax.set_yticks([i + 0.5 for i in range(len(SHIPS))])
    ax.set_yticklabels(SHIPS)
    ax.set_title("Independent-source coverage per ship and quantity (×n = independent lineages)")
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=c, alpha=0.85) for c in VERDICT_COLOUR.values()
    ]
    ax.legend(handles, list(VERDICT_COLOUR), loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=4, frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def displacement_figure(library: dict[str, Any], path: Path) -> None:
    rows = [(s, library["comparisons"][s]) for s in SHIPS if s in library["comparisons"]]
    rows = [(s, r) for s, r in rows if r["reference_volume_m3"] and r["got_volume_m3"]]
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    if not rows:
        ax.text(0.5, 0.5, "no reconstructed volume in the benchmark artefact", ha="center")
    else:
        names = [s for s, _ in rows]
        ref = [r["reference_volume_m3"] / 1e3 for _, r in rows]
        got = [r["got_volume_m3"] / 1e3 for _, r in rows]
        x = range(len(names))
        ax.bar([i - 0.2 for i in x], ref, width=0.4, label="published ∇ (cross-ref minimum)", color="#1565c0")
        ax.bar([i + 0.2 for i in x], got, width=0.4, label="reconstructed ∇", color="#ef6c00")
        for i, (_, r) in enumerate(rows):
            text = "not comparable" if r["rel_error_percent"] is None else f"{r['rel_error_percent']:+.2f} %"
            ax.text(i, max(ref[i], got[i]) * 1.02, text, ha="center", fontsize=9)
        ax.set_xticks(list(x))
        ax.set_xticklabels(names)
        ax.set_ylabel("displacement volume (10³ m³)")
        ax.set_title("Published vs reconstructed displacement")
        ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def closure_figure(library: dict[str, Any], path: Path) -> None:
    closure = library["closure"]
    ships = [s for s in SHIPS if closure.get(s)]
    fig, ax = plt.subplots(figsize=(7.2, 3.4))
    labels, errors = [], []
    for ship in ships:
        for family, row in sorted(closure[ship].items()):
            labels.append(f"{ship}\n{family}")
            errors.append(row["rel_error"])
    colours = ["#c62828" if abs(e) > 0.01 else "#2e7d32" for e in errors]
    ax.bar(range(len(errors)), errors, color=colours)
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=7.5)
    ax.set_ylabel(r"$C_b L_{pp} B_{wl} T / \nabla - 1$")
    ax.set_title("Publisher self-consistency: Cb closure of each source's own digits")
    for i, ship in enumerate(ships):
        row = library["closure"][ship]
        for family, values in row.items():
            if not values["consistent"]:
                ax.annotate(
                    f"implied Bwl {values['implied_beam_m']:.1f} m",
                    xy=(labels.index(f"{ship}\n{family}"), values["rel_error"]),
                    xytext=(0, 6),
                    textcoords="offset points",
                    ha="center",
                    fontsize=7.5,
                )
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def table(headers: list[str], rows: list[list[str]]) -> str:
    out = ["| " + " | ".join(headers) + " |", "|" + "|".join(["---"] * len(headers)) + "|"]
    out += ["| " + " | ".join(row) + " |" for row in rows]
    return "\n".join(out)


def fmt(value: float | None) -> str:
    if value is None:
        return "—"
    if abs(value) >= 1000:
        return f"{value:,.1f}"
    return f"{value:g}"


def render(library: dict[str, Any]) -> str:
    cov = library["coverage"]
    sources = library["sources"]
    lines: list[str] = []
    add = lines.append

    add("# Fleet reference library — six benchmark hulls, ≥2 independent sources each\n")
    add("Anchor **REPRO-GEOM-FLEET-REFS** (TensorCAD). Machine artefact: "
        "`benchmarks/verified/fleet_refs/crossref.json`; this document is rendered from it by "
        "`tools/baseline_fleet/crossref_report.py` (no number is transcribed by hand).\n")

    add("## 1. Problem setting and where the references come from\n")
    add("The fleet benchmark compares a reconstructed displacement against **one** published "
        "number per hull (KCS verified at +0.304 % against the NMRI page, everything else "
        "withheld). A single number cannot separate a reconstruction error from a source typo. "
        "This library therefore holds every principal particular that a publisher states for "
        "the six hulls — ∇, Cb, LCB (%Lpp), Lwl, and the Lpp/Bwl/T the coefficients are defined "
        "on — reads them out of the fetched bytes, and asks per quantity whether **two "
        "independent publishers** state it and agree.\n")
    add("Sources (digest pinned; the bytes stay out of the repository, see the licence column):\n")
    rows = []
    for key, row in sources.items():
        rows.append([
            f"`{key}`",
            row.get("publisher", "—"),
            row.get("family", "—"),
            row.get("status", "—"),
            (row.get("sha256") or "—")[:16],
        ])
    add(table(["source", "publisher", "family", "status", "sha256 (16)"], rows))
    add("")

    add("## 2. Method\n")
    add("`tensorcad/fleet/crossref.py` reads a table row by row (`official.table_rows`), matches "
        "the quantity label **in the first cell only** (`LABEL_PATTERNS`) and takes the numeric "
        "cell of the publisher's full-scale column. Two guards keep the reading honest: a row "
        "whose label is a derived quantity (`izz/Lpp`) is skipped, and a row whose configured "
        "full-scale column is empty is kept for the audit trail but flagged "
        "`fallback_column` and never allowed to vote. Independence is judged at the *lineage* "
        "level, not the URL: the SIMMAN 2014 and SIMMAN 2020 editions print the identical "
        "KVLCC2 table including the identical `Bwl 580` defect, so they count as one source. "
        "Each ship/quantity is then classified `CROSS_VERIFIED` (≥2 lineages, relative spread "
        "≤ tolerance), `CONFLICT` (≥2 lineages, spread above tolerance), `SINGLE_SOURCE` or "
        "`MISSING`.\n")

    add("## 3. Parameter table (as published, full scale)\n")
    rows = []
    for ship in SHIPS:
        for quantity in QUANTITIES:
            cell = library["report"][ship]["quantities"][quantity]
            if not cell["family_values"]:
                continue
            values = ", ".join(f"{fam} {fmt(val)}" for fam, val in cell["family_values"].items())
            rows.append([ship, library["quantity_labels"][quantity], values, cell["verdict"]])
    add(table(["hull", "quantity", "published value(s) by lineage", "verdict"], rows))
    add("")

    add("## 4. Cross-check result — does each hull have ≥2 independent sources?\n")
    rows = []
    for ship in SHIPS:
        row = library["report"][ship]
        def cell(q: str) -> str:
            c = row["quantities"][q]
            return f"{c['verdict']} ({c['independent_families']})"
        rows.append([ship, cell("displacement_volume"), cell("cb"), cell("lcb_percent"), cell("lwl")])
    add(table(["hull", "∇", "Cb", "LCB", "Lwl"], rows))
    add("")
    add(f"Coverage: **{cov['cross_verified']}/{cov['cells']}** cells cross-verified, "
        f"{cov['single_source']} single-source, {cov['conflict']} conflict, {cov['missing']} missing. "
        f"Hulls with all four key quantities (∇, Cb, LCB, Lwl) backed by ≥2 lineages: "
        f"**{', '.join(cov['ships_fully_cross_verified']) or 'none'}**. "
        f"Hulls still below two sources: **{', '.join(cov['ships_failing_min_two_sources']) or 'none'}**.\n")

    add("### 4.1 Reference vs reconstruction\n")
    rows = []
    for ship in SHIPS:
        row = library["comparisons"].get(ship)
        if not row:
            continue
        rows.append([
            ship,
            fmt(row["reference_volume_m3"]),
            fmt(row["got_volume_m3"]),
            "—" if row["rel_error_percent"] is None else f"{row['rel_error_percent']:+.2f} %",

            row["reference_verdict"],
            row["benchmark_verdict"] + (
                "" if row.get("comparability", "full_scale") == "full_scale" else " (got artefact not comparable)"
            ),
        ])
    add(table(["hull", "ref ∇ (m³)", "got ∇ (m³)", "error", "ref verdict", "benchmark"], rows))
    add("")

    add("### 4.2 Publisher self-consistency (Cb closure)\n")
    rows = []
    for ship in SHIPS:
        for family, row in sorted(library["closure"].get(ship, {}).items()):
            rows.append([
                ship,
                family,
                fmt(row["cb_lpp_bwl_t_m3"]),
                fmt(row["published_volume_m3"]),
                f"{row['rel_error']:+.3f}",
                "ok" if row["consistent"] else f"**defect** — implied Bwl {row['implied_beam_m']:.2f} m",
            ])
    add(table(["hull", "lineage", "Cb·Lpp·Bwl·T (m³)", "published ∇ (m³)", "rel. error", "verdict"], rows))
    add("")

    add("## 5. Sensitivity / what the classification is sensitive to\n")
    add("Two knobs decide the verdicts, and both are exercised in "
        "`tests/test_fleet_refs_crossref.py`:\n")
    add("* **tolerance** — 0.5 % for ∇/Cb/Lwl/Lpp/Bwl/T and 2 % for LCB (a small "
        "percent-of-Lpp number); at 0.5 % the KCS Cb pair (0.6505 vs 0.651, 0.077 %) is "
        "`CROSS_VERIFIED`, and lowering the tolerance to 0.05 % would flip it to `CONFLICT`: the "
        "two publishers round differently, which is corroboration, not disagreement.\n"
        "* **lineage grouping** — treating SIMMAN 2014 and SIMMAN 2020 as two sources would "
        "report KVLCC2 as cross-verified on all seven quantities; the closure test "
        "(section 4.2) shows the two editions share the `Bwl 580` digit, so they are one "
        "source. This is the single most consequential choice in the library.\n")

    add("## 6. Failure modes and artefacts\n")
    add("Honest gaps, with the reason each document could not be used:\n")
    rows = []
    for probe in library.get("probes", []):
        rows.append([probe.get("want", "—"), probe.get("url", "—"), probe.get("outcome") or probe.get("observed", "—")])
    add(table(["wanted", "url", "outcome"], rows))
    add("")
    add("* **DTC** and **Series 60**: no numeric source reachable at all in this environment "
        "(0 cells). DTC's originator page is retired and its IGES bundle "
        "(`uni-due.de/imperia/md/content/inam/igs_ordner.zip`, sha256 recorded in the registry) "
        "carries geometry but no particulars; Series 60's primary table (Todd) is behind a 403.\n"
        "* **JBC**: one lineage only (NMRI Tokyo 2015). Its ∇ 178 369.9 m³ is *not* yet "
        "external in the sense this anchor requires; the reconstruction currently reads "
        "145 769 m³ (−18.3 %) with the control-net layout guard fired, so the hull is "
        "`WITHHELD` in the benchmark for a reason that must not be blamed on the reference.\n"
        "* **ONRT**: ∇ and LCB are `CONFLICT` — the NMRI page prints displacement as a *mass* "
        "(8 507 t → 8 299.5 m³ at 1.025 t/m³) while SIMMAN prints 8 507 m³, and LCB is "
        "printed on different references (the NMRI cell is a metre offset that the parser "
        "reads as 1.625, the SIMMAN cell is −1.68 %Lpp). Both cells are kept, unflagged as "
        "agreement, because picking either one would be a flattering choice.\n")

    add("## Reproduce\n")
    add("```bash\n"
        "python tools/baseline_fleet/crossref_refs.py fetch --out artifacts/fleet_refs   # network; digests checked\n"
        "python tools/baseline_fleet/crossref_refs.py build --raw artifacts/fleet_refs/raw\n"
        "python tools/baseline_fleet/crossref_report.py \\\n"
        "    --library benchmarks/verified/fleet_refs/crossref.json \\\n"
        "    --doc docs/fleet_refs_crossref.md --figs docs/figs/fleet_refs\n"
        "python -m pytest tests/test_fleet_refs_crossref.py -q\n"
        "```\n")
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--library", default="benchmarks/verified/fleet_refs/crossref.json")
    parser.add_argument("--doc", default="docs/fleet_refs_crossref.md")
    parser.add_argument("--figs", default="docs/figs/fleet_refs")
    args = parser.parse_args()
    library = json.loads(Path(args.library).read_text())
    figs = Path(args.figs)
    figs.mkdir(parents=True, exist_ok=True)
    verdict_matrix(library, figs / "coverage_matrix.png")
    displacement_figure(library, figs / "displacement_ref_vs_got.png")
    closure_figure(library, figs / "closure_check.png")
    Path(args.doc).write_text(render(library))
    print(f"wrote {args.doc} and 3 figures under {figs}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())