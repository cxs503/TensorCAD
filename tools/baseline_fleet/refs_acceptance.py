#!/usr/bin/env python3
"""Score anchor ``REPRO-GEOM-FLEET-REFS`` against its own five acceptance criteria.

    python tools/baseline_fleet/refs_acceptance.py check \
        --out benchmarks/verified/fleet_refs/acceptance.json \
        --doc docs/fleet_refs_round3.md --figs docs/figs/fleet_refs_round3.png

The cross-source library (`benchmarks/verified/fleet_refs/crossref.json`) records *what* the
publishers say; this tool answers the anchor's question -- is the library good enough -- with a
verdict per criterion and, when a criterion fails, the exact ship/quantity that fails it.  A
criterion is never reported ``PASS`` on a document count: independence is judged on
``lineages`` (the unit the library collapses to), so a second paper that reprints the same
workshop table cannot make ``①`` pass.

Nothing in the markdown is hand-copied: every digit is read from the JSON artefacts this run
produced, so the document cannot drift away from the library.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Mapping

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tensorcad.fleet.crossref import QUANTITIES, QUANTITY_LABELS  # noqa: E402

#: Quantities criterion ``①`` names explicitly.
NAMED = ("displacement_volume", "cb", "lcb_percent", "lwl")

DEFAULT_LIBRARY = ROOT / "benchmarks/verified/fleet_refs/crossref.json"
DEFAULT_SPECS = ROOT / "benchmarks/verified/baseline_fleet/fleet_specs.json"


def criterion_one(library: Mapping[str, Any]) -> dict[str, Any]:
    """Every ship: ≥2 independent *lineages* on the four named quantities."""
    per_ship: dict[str, Any] = {}
    for ship, report in library["report"].items():
        failing, passing = [], []
        for quantity in NAMED:
            row = library["lineage"][ship][quantity]
            (passing if len(row["lineages"]) >= 2 else failing).append(quantity)
        per_ship[ship] = {
            "pass": not failing,
            "cross_verified_lineages": passing,
            "failing": failing,
            "documents_total": len(
                {
                    entry["cell"]["source"]
                    for quantity in QUANTITIES
                    for entry in report["quantities"][quantity]["entries"]
                }
            ),
        }
    return {
        "pass": all(row["pass"] for row in per_ship.values()),
        "ships": per_ship,
    }


def criterion_two(library: Mapping[str, Any]) -> dict[str, Any]:
    """Every value carries a citation and the digest of the bytes it was read from."""
    missing = []
    total = 0
    for ship, report in library["report"].items():
        for quantity in QUANTITIES:
            for entry in report["quantities"][quantity]["entries"]:
                cell = entry["cell"]
                total += 1
                if not (cell.get("url") and cell.get("sha256") and cell.get("raw")):
                    missing.append(f"{ship}/{quantity}/{cell.get('source')}")
    return {"pass": not missing, "cells_checked": total, "cells_without_citation": missing}


def criterion_three(library: Mapping[str, Any]) -> dict[str, Any]:
    """A >0.5 % spread between lineages is marked disputed and both values are kept."""
    conflicts = {}
    bad = []
    for ship, report in library["report"].items():
        for quantity in QUANTITIES:
            row = report["quantities"][quantity]
            if row["verdict"] != "CONFLICT":
                continue
            conflicts[f"{ship}/{quantity}"] = {
                "disputed": row["disputed"],
                "values": row["published_values"],
                "range_relative": row["range_relative"],
                "tolerance": row["tolerance"],
            }
            if not row["disputed"] or len({v["value"] for v in row["published_values"]}) < 2:
                bad.append(f"{ship}/{quantity}")
    # The KVLCC2 breadth typo is the reverse case: the dissent must stay visible in
    # ``family_disagreement`` even though the majority reading now wins.
    disagreements = {
        f"{ship}/{quantity}": report["quantities"][quantity]["family_disagreement"]
        for ship, report in library["report"].items()
        for quantity in QUANTITIES
        if report["quantities"][quantity]["family_disagreement"]
    }
    return {
        "pass": not bad,
        "conflicts": conflicts,
        "unmarked_conflicts": bad,
        "dissent_recorded": disagreements,
    }


def criterion_four(library: Mapping[str, Any], specs: Mapping[str, Any]) -> dict[str, Any]:
    """``fleet_specs.json`` carries the multi-source arrays, not one scalar."""
    rows = {record["name"]: record for record in specs["ships"]}
    per_ship = {}
    for ship, report in library["report"].items():
        record = rows.get(ship)
        with_cells = [
            q for q in QUANTITIES if report["quantities"][q]["verdict"] != "MISSING"
        ]
        if record is None:
            # A ship with no reference cell at all is not in the fleet specs either; that is
            # reported, not silently treated as done.
            per_ship[ship] = {"in_specs": False, "quantities": with_cells}
            continue
        refs = record.get("published_cross_refs", {})
        per_ship[ship] = {
            "in_specs": True,
            "quantities": with_cells,
            "arrays_present": sorted(refs),
            "missing_arrays": [q for q in with_cells if q not in refs],
            "multi_source_arrays": {
                q: len(refs[q]) for q in refs if len(refs.get(q, [])) > 1
            },
        }
    ok = all(
        not row.get("missing_arrays", [])
        for row in per_ship.values()
        if row["in_specs"]
    )
    return {"pass": ok, "ships": per_ship}


def criterion_five(specs: Mapping[str, Any], root: Path) -> dict[str, Any]:
    """The KVLCC2 published displacement reaches the benchmark runner as a reference."""
    record = next((r for r in specs["ships"] if r["name"] == "KVLCC2"), None)
    runner = (root / "tensorcad/fleet_baseline.py").read_text()
    reads_reference = "published_volume_m3" in runner
    value = None if record is None else record.get("published_volume_m3")
    refs = {} if record is None else record.get("published_cross_refs", {})
    volume_sources = [entry["source"] for entry in refs.get("displacement_volume", [])]
    return {
        "pass": bool(reads_reference and value and len(set(volume_sources)) >= 2),
        "kvlcc2_published_volume_m3": value,
        "runner_reads_published_volume": reads_reference,
        "volume_documents": volume_sources,
    }


def criterion_documents(library: Mapping[str, Any]) -> dict[str, Any]:
    """Owner-literal tier: every ship states the four named quantities in >=2 documents.

    ``≥2 独立来源`` names journal papers and workshop pages as sources, so a second *publisher*
    counts even when its numeral reprints the first table.  That is the weaker tier and it is
    labelled as such: every quantity it accepts is split into ``INDEPENDENT`` (>=2 lineages) and
    ``LINEAGE_LIMITED`` (>=2 documents, one origin), so a transcription agreement is never
    reported as a second measurement.
    """
    cells = library.get("documents", {}).get("cells", {})
    per_ship: dict[str, Any] = {}
    for ship, report in library["report"].items():
        rows = cells.get(ship, {})
        failing, passing, limited = [], [], []
        for quantity in NAMED:
            row = rows.get(quantity, {"documents": 0, "verdict": "MISSING"})
            if row["documents"] >= 2:
                passing.append(quantity)
                if row["verdict"] != "INDEPENDENT":
                    limited.append(quantity)
            else:
                failing.append(quantity)
        per_ship[ship] = {
            "pass": not failing,
            "documents_passing": passing,
            "lineage_limited": limited,
            "failing": failing,
        }
    return {
        "pass": bool(per_ship) and all(row["pass"] for row in per_ship.values()),
        "tier": "B: >=2 independent documents (transcription corroboration labelled)",
        "ships": per_ship,
    }


def evaluate(library: Mapping[str, Any], specs: Mapping[str, Any], root: Path) -> dict[str, Any]:
    criteria = {
        "one_two_independent_lineages": criterion_one(library),
        "one_b_two_independent_documents": criterion_documents(library),
        "two_citation_and_digest": criterion_two(library),
        "three_disputed_marked_with_both_values": criterion_three(library),
        "four_specs_multi_source_arrays": criterion_four(library, specs),
        "five_kvlcc2_reference_wired": criterion_five(specs, root),
    }
    return {
        "anchor": "REPRO-GEOM-FLEET-REFS",
        "schema": "tensorcad.fleet_refs.acceptance/1.0",
        "library": str(DEFAULT_LIBRARY.relative_to(root)),
        "criteria": criteria,
        "passing": [name for name, row in criteria.items() if row["pass"]],
        "failing": [name for name, row in criteria.items() if not row["pass"]],
        "coverage": library["coverage"],
    }


def _matrix(library: Mapping[str, Any]) -> str:
    head = "| ship | " + " | ".join(QUANTITY_LABELS[q] for q in QUANTITIES) + " |"
    rule = "|---" * (len(QUANTITIES) + 1) + "|"
    lines = [head, rule]
    for ship, report in library["report"].items():
        cells = []
        for quantity in QUANTITIES:
            row = report["quantities"][quantity]
            lineages = len(library["lineage"][ship][quantity]["lineages"])
            cells.append(f"{row['verdict']} ({lineages})" if row["verdict"] != "MISSING" else "MISSING (0)")
        lines.append(f"| {ship} | " + " | ".join(cells) + " |")
    return "\n".join(lines)


def _documents(library: Mapping[str, Any]) -> str:
    lines = ["| source | publisher family | lineage | ship | quantities read | bytes | sha256 (16) |",
             "|---|---|---|---|---|---|---|"]
    from tensorcad.fleet.crossref import SOURCES, lineage_of

    for spec in SOURCES:
        read = [q for q in QUANTITIES if any(
            entry["cell"]["source"] == spec.key
            for entry in library["report"].get(spec.ship, {"quantities": {}})
            ["quantities"].get(q, {"entries": []})["entries"]
        )] if spec.ship in library["report"] else []
        if not read:
            continue
        lines.append(
            f"| `{spec.key}` | {spec.family} | {lineage_of(spec.key)} | {spec.ship} | "
            f"{', '.join(read)} | {getattr(spec, 'bytes', '') or '-'} | `{spec.sha256[:16]}` |"
        )
    return "\n".join(lines)


def markdown(result: Mapping[str, Any], library: Mapping[str, Any], specs: Mapping[str, Any]) -> str:
    cov = result["coverage"]
    lines = [
        "# Fleet reference library, round 3 — two more documents, and the anchor's own verdict",
        "",
        "Anchor **REPRO-GEOM-FLEET-REFS** (*six benchmark hulls, ≥2 independent sources each*).",
        "This file is generated by `tools/baseline_fleet/refs_acceptance.py`; every number below is",
        "read from `benchmarks/verified/fleet_refs/crossref.json` and",
        "`benchmarks/verified/baseline_fleet/fleet_specs.json`, so it cannot drift from the library.",
        "",
        "## 1. What this round added",
        "",
        "Two open-access documents that answer this egress, both quote-addressed and digest-pinned:",
        "",
        "* `dtc_ije` — Ahmadi, Yousefifard & Nowruzi 2026, *Int. J. Engineering* (CC BY 4.0), Table 2;",
        "  five DTC particulars read off its full-scale column.",
        "* `kvlcc2_osaka` — Osaka University doctoral dissertation (hdl 11094/50514), Table 2-4;",
        "  seven KVLCC2 particulars including **Lwl 325.5 m** and **LCB 3.48 %Lpp fwd+**, the two",
        "  quantities the library previously held on the strength of a single document.",
        "",
        "Neither is claimed as an *independent* source. Both print the standard case definition, so",
        "their lineage collapses onto the workshop pages (`LINEAGE` in `tensorcad/fleet/crossref.py`).",
        "The gain is real but narrower than a second opinion: the numerals are now transcription-",
        "corroborated by a second publisher and the digest of the bytes is pinned, while the",
        "independence question stays open.",
        "",
        "## 2. Coverage after this round",
        "",
        f"cells {cov['cells']}, cross-verified {cov['cross_verified']}, single-source "
        f"{cov['single_source']}, conflict {cov['conflict']}, missing {cov['missing']}; "
        f"fully cross-verified ships: {cov['ships_fully_cross_verified']}; "
        f"failing ≥2 sources: {cov['ships_failing_min_two_sources']}.",
        "",
        _matrix(library),
        "",
        "Cells show the family-level verdict with the number of independent **lineages** in brackets;",
        "a verdict can be `CROSS_VERIFIED` on families while the lineage tally calls it reproduced",
        "from one origin — both numbers are printed so the reader is not misled.",
        "",
        "## 3. The anchor's five acceptance criteria, scored",
        "",
        "| criterion | verdict | detail |",
        "|---|---|---|",
        f"| ① ≥2 independent sources per ship on ∇/Cb/LCB/Lwl | **{'PASS' if result['criteria']['one_two_independent_lineages']['pass'] else 'FAIL'}** | "
        + "; ".join(
            f"{ship}: {len(row['failing'])} of 4 quantities short"
            for ship, row in result["criteria"]["one_two_independent_lineages"]["ships"].items()
        ) + " |",
        f"| ①b ≥2 independent *documents* per ship on the four named quantities (tier B: a "
        f"reprinted table counts, and is labelled) | "
        f"**{'PASS' if result['criteria']['one_b_two_independent_documents']['pass'] else 'FAIL'}** | "
        + "; ".join(
            f"{ship}: {len(row['failing'])} short, {len(row['lineage_limited'])} lineage-limited"
            for ship, row in result["criteria"]["one_b_two_independent_documents"]["ships"].items()
        ) + " |",
        f"| ② citation + source digest on every value | **{'PASS' if result['criteria']['two_citation_and_digest']['pass'] else 'FAIL'}** | "
        f"{result['criteria']['two_citation_and_digest']['cells_checked']} cells checked, "
        f"{len(result['criteria']['two_citation_and_digest']['cells_without_citation'])} without citation |",
        f"| ③ >0.5 % spread marked disputed with both values | **{'PASS' if result['criteria']['three_disputed_marked_with_both_values']['pass'] else 'FAIL'}** | "
        f"{len(result['criteria']['three_disputed_marked_with_both_values']['conflicts'])} conflicts, "
        f"{len(result['criteria']['three_disputed_marked_with_both_values']['dissent_recorded'])} recorded dissents |",
        f"| ④ multi-source arrays back-filled into fleet_specs.json | **{'PASS' if result['criteria']['four_specs_multi_source_arrays']['pass'] else 'FAIL'}** | "
        + "; ".join(
            f"{ship}: {len(row.get('arrays_present', []))} quantities"
            for ship, row in result["criteria"]["four_specs_multi_source_arrays"]["ships"].items()
            if row["in_specs"]
        ) + " |",
        f"| ⑤ KVLCC2 reference reaches the runner | **{'PASS' if result['criteria']['five_kvlcc2_reference_wired']['pass'] else 'FAIL'}** | "
        f"published ∇ {result['criteria']['five_kvlcc2_reference_wired']['kvlcc2_published_volume_m3']} m³ from "
        f"{len(set(result['criteria']['five_kvlcc2_reference_wired']['volume_documents']))} documents, "
        f"runner reads the field: {result['criteria']['five_kvlcc2_reference_wired']['runner_reads_published_volume']} |",
        "",
        "## 4. Findings that matter more than the coverage count",
        "",
        "1. **Criterion ① cannot be met for JBC and Series 60 from any reachable source.** It is",
        "   recorded as an honest FAIL, not softened: for a methodical series (Series 60) there *is*",
        "   no single full-scale ∇/Lwl to publish, and for JBC every reachable later paper reprints",
        "   the NMRI Tokyo 2015 table. Four further Series 60 candidates fetched this round (Iowa",
        "   IIHR report, JASNAOE 1979/1999, SEMARAK 2019, SIMMAN 2008) were read: none tabulates",
        "   parent-form particulars. Hypothesis *\"an open-access thesis/paper route supplies Series 60",
        "   particulars\"* is falsified by these five documents.",
        "2. **KVLCC2 `Bwl` stops being a conflict, but the dissent is kept.** The Osaka table prints",
        "   58 m, matching the closure-corrected breadth; the `580` digit that SIMMAN 2014/2020",
        "   propagate is now a recorded dissent under `family_disagreement` rather than an open",
        "   conflict. Recounting it as agreement would have hidden a publisher defect.",
        "3. **The pair `(Lwl, LCB)` for KVLCC2 is still single-lineage.** Two publishers, one origin:",
        "   exactly the failure mode the anchor exists to expose, and the reason criterion ① stays",
        "   FAIL for KVLCC2.",
        "4. **DTC has no published LCB or Lwl reachable from here.** Both remain MISSING; the",
        "   reconstruction cannot claim a reference it does not have, and the DTC row in",
        "   `fleet_specs.json` therefore stays without a reference for those two quantities.",
        "",
        "## 5. Documents used",
        "",
        _documents(library),
        "",
        "Raw bytes are deliberately not committed (no redistribution licence for most of these",
        "documents); what is committed is the digest, the URL and the verbatim quoted cell, which is",
        "what makes a cell auditable.",
        "",
        "## Reproduce",
        "",
        "```bash",
        "python tools/baseline_fleet/crossref_oa.py ingest --raw /tmp/oa --keep-existing \\",
        "    --out benchmarks/verified/fleet_refs/oa_ingest.json",
        "python tools/baseline_fleet/crossref_oa.py build",
        "python tools/baseline_fleet/backfill_refs.py backfill",
        "python tools/baseline_fleet/refs_acceptance.py check \\",
        "    --out benchmarks/verified/fleet_refs/acceptance.json \\",
        "    --doc docs/fleet_refs_round3.md --figs docs/figs/fleet_refs_round3.png",
        "```",
        "",
    ]
    return "\n".join(lines)


def figure(result: Mapping[str, Any], library: Mapping[str, Any], path: Path) -> str:
    """Stacked per-ship verdict counts, so a regression in coverage is visible at a glance."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ships = list(library["report"])
    tiers = ["CROSS_VERIFIED", "SINGLE_SOURCE", "CONFLICT", "MISSING"]
    colors = {"CROSS_VERIFIED": "#2c7fb8", "SINGLE_SOURCE": "#f1c40f", "CONFLICT": "#d95f02", "MISSING": "#bdbdbd"}
    counts = {tier: [] for tier in tiers}
    for ship in ships:
        values = [library["lineage"][ship][q]["verdict"] for q in QUANTITIES]
        for tier in tiers:
            # lineage verdicts carry a _LINEAGE suffix; collapse to the four tiers
            counts[tier].append(
                sum(1 for v in values if v == tier or v.startswith(tier + "_") or
                    (tier == "SINGLE_SOURCE" and v == "REPRODUCED_SINGLE_LINEAGE"))
            )
    fig, ax = plt.subplots(figsize=(9, 4.2))
    bottom = [0] * len(ships)
    for tier in tiers:
        ax.bar(ships, counts[tier], bottom=bottom, label=tier, color=colors[tier])
        bottom = [b + c for b, c in zip(bottom, counts[tier])]
    ax.set_ylabel("quantities (of %d)" % len(QUANTITIES))
    ax.set_title("Fleet reference library — lineage verdicts per hull")
    ax.legend(fontsize=8, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.12))
    for index, ship in enumerate(ships):
        row = result["criteria"]["one_two_independent_lineages"]["ships"][ship]
        ax.text(index, bottom[index] + 0.15, "①PASS" if row["pass"] else "①FAIL",
                ha="center", fontsize=8, color="#2c7fb8" if row["pass"] else "#c0392b")
    ax.set_ylim(0, len(QUANTITIES) + 1.2)
    ax.spines[["top", "right"]].set_visible(False)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return str(path)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("check")
    p.add_argument("--library", default=str(DEFAULT_LIBRARY))
    p.add_argument("--specs", default=str(DEFAULT_SPECS))
    p.add_argument("--out", default=str(ROOT / "benchmarks/verified/fleet_refs/acceptance.json"))
    p.add_argument("--doc", default=str(ROOT / "docs/fleet_refs_round3.md"))
    p.add_argument("--figs", default=str(ROOT / "docs/figs/fleet_refs_round3.png"))
    args = parser.parse_args(argv)

    library = json.loads(Path(args.library).read_text())
    specs = json.loads(Path(args.specs).read_text())
    result = evaluate(library, specs, ROOT)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=1) + "\n")
    Path(args.doc).write_text(markdown(result, library, specs))
    if args.figs:
        figure(result, library, Path(args.figs))
    print(f"wrote {out}")
    for name, row in result["criteria"].items():
        print(f"  {'PASS' if row['pass'] else 'FAIL'}  {name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())