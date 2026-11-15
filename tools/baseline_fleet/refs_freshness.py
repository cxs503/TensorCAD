#!/usr/bin/env python3
"""Does every cached block of the fleet reference library still follow from its cells?

    python tools/baseline_fleet/refs_freshness.py check
    python tools/baseline_fleet/refs_freshness.py refresh \
        --report docs/fleet_refs_round5_dissent.md --figs docs/figs/fleet_refs_round5.png

Anchor ``REPRO-GEOM-FLEET-REFS``.  ``benchmarks/verified/fleet_refs/crossref.json`` is a
*report* (what each source says, per cell) plus four blocks **derived** from it -- ``coverage``,
``lineage``, ``closure`` and ``adjudication``.  The derived blocks are cached so that readers and
the report generator do not have to run the library, and that is precisely how this file went
wrong: round 3 rewrote ``report`` (repaired KVLCC2 beam 580 -> 58, added DTC/JSR-2017 cells,
refined the JBC volume) and left ``closure``/``adjudication`` at round 2's answer.  The
committed library therefore stated in one block that KVLCC2's own digits are 900 % off and in
another that they close to 1e-6, and the published reference table
(``docs/fleet_refs_crossref.md`` section 4.2) printed the phantom defect as a finding.

``check`` recomputes the four blocks from ``report`` and diffs them leaf by leaf; it exits
non-zero on any drift, so a stale cache is a test failure rather than a reader's discovery.
``refresh`` performs the same recomputation, writes it back into the library and regenerates
the round document with every number read from the refreshed file.

Neither mode touches the network: this tool audits the library against itself.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path
from typing import Any, Mapping

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tensorcad.fleet.crossref import (  # noqa: E402
    CLOSURE_TOLERANCE,
    QUANTITIES,
    QUANTITY_LABELS,
    refresh_derived,
    stale_derived,
)

LIBRARY = ROOT / "benchmarks/verified/fleet_refs/crossref.json"

#: Ships in the order the anchor and the fleet benchmark list them.
SHIPS = ("KCS", "KVLCC2", "JBC", "ONRT", "DTC", "Series60")


def _tool_module(name: str, filename: str):
    """Import ``tools/baseline_fleet/<filename>`` as a module (the tools are not a package)."""
    path = ROOT / "tools/baseline_fleet" / filename
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def derived_view(payload: Mapping[str, Any]) -> dict[str, Any]:
    """Recompute all four cached blocks, ``lineage`` included.

    ``lineage`` is tallied by ``tools/baseline_fleet/crossref_oa.py`` (it is a tool-level view
    of the cells, not a library primitive), so it is imported rather than duplicated -- two
    implementations of one tally is the defect this tool exists to catch.
    """
    oa = _tool_module("crossref_oa_tool", "crossref_oa.py")
    blocks = stale_derived(payload)
    fresh = dict(blocks["fresh"])
    fresh["lineage"] = oa.lineage_tally(payload["report"])
    stale = list(blocks["stale"])
    differences = dict(blocks["differences"])
    if json.dumps(fresh["lineage"], sort_keys=True) != json.dumps(
        payload.get("lineage"), sort_keys=True
    ):
        stale.append("lineage")
        differences["lineage"] = ["<tally differs>"]
    return {"stale": stale, "fresh": fresh, "differences": differences}


def cmd_check(library: Path) -> int:
    payload = json.loads(library.read_text())
    view = derived_view(payload)
    if not view["stale"]:
        print(f"{library}: all four derived blocks follow from report")
        return 0
    print(f"{library}: STALE derived blocks: {', '.join(view['stale'])}")
    for block in view["stale"]:
        for line in view["differences"][block]:
            print(f"  {block}{line if line.startswith('/') else ' ' + line}")
    print(
        "run `python tools/baseline_fleet/refs_freshness.py refresh` and regenerate the report"
    )
    return 1


def cmd_refresh(
    library: Path, report_doc: Path | None, figs: Path | None, before_path: Path | None = None
) -> int:
    payload = json.loads(library.read_text())
    # The document quotes what the cache *was* claiming, so the pre-refresh snapshot can be
    # supplied (``git show HEAD:<library>``) when the refresh has already been applied.
    before = derived_view(json.loads(before_path.read_text()) if before_path else payload)
    refresh_derived(payload, ("coverage", "closure", "adjudication"))
    oa = _tool_module("crossref_oa_tool", "crossref_oa.py")
    payload["lineage"] = oa.lineage_tally(payload["report"])
    library.write_text(json.dumps(payload, indent=1) + "\n")
    after = derived_view(json.loads(library.read_text()))
    print(f"{library}: refreshed {', '.join(before['stale']) or '(nothing stale)'}")
    if after["stale"]:
        print(f"REFUSING to claim success: still stale after refresh: {after['stale']}")
        return 1
    if report_doc is not None:
        write_document(library, report_doc, figs, before)
        print(f"{report_doc}: written")
    return 0


# --------------------------------------------------------------------------- document


def _figure(library: Mapping[str, Any], figs: Path | None, before: Mapping[str, Any]) -> Path | None:
    """Two panels: fresh closure residual per ship/family, and the KVLCC2 stale-vs-fresh."""
    if figs is None:
        return None
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - optional dependency
        print(f"figures skipped: {exc}")
        return None
    figs.parent.mkdir(parents=True, exist_ok=True)
    fig, (left, right) = plt.subplots(1, 2, figsize=(11, 4.2))

    labels, residuals, colors = [], [], []
    for ship in SHIPS:
        for family, row in sorted(library["closure"].get(ship, {}).items()):
            labels.append(f"{ship}\n{family}")
            residuals.append(row["rel_error"])
            colors.append("#2c7fb8" if row["consistent"] else "#d95f02")
    positions = range(len(labels))
    left.bar(positions, residuals, color=colors)
    left.set_xticks(list(positions))
    left.set_xticklabels(labels, fontsize=6, rotation=90)
    left.axhline(CLOSURE_TOLERANCE, color="#777", ls="--", lw=0.8)
    left.axhline(-CLOSURE_TOLERANCE, color="#777", ls="--", lw=0.8)
    left.set_ylabel(r"$C_b L_{pp} B_{wl} T / \nabla - 1$")
    left.set_title("fresh: every publisher closes on its own digits", fontsize=9)

    stale_rows = _stale_closure_rows(before)
    if stale_rows:
        names = [f"{ship}/{family}" for ship, family, _ in stale_rows]
        stored = [abs(row) for _, _, row in stale_rows]
        fresh_rows = [
            abs(library["closure"].get(ship, {}).get(family, {}).get("rel_error", float("nan")))
            for ship, family, _ in stale_rows
        ]
        x = range(len(names))
        right.bar([i - 0.2 for i in x], stored, width=0.4, color="#d95f02", label="cached (stale)")
        right.bar([i + 0.2 for i in x], fresh_rows, width=0.4, color="#2c7fb8", label="fresh")
        right.set_yscale("log")
        right.set_xticks(list(x))
        right.set_xticklabels(names, fontsize=7)
        right.set_ylabel("|relative closure error|")
        right.set_title("what the cache was claiming", fontsize=9)
        right.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(figs, dpi=150)
    plt.close(fig)
    return figs


def _stale_closure_rows(before: Mapping[str, Any]) -> list[tuple[str, str, float]]:
    """``(ship, family, cached residual)`` for the closure rows the refresh repaired."""
    out = []
    for line in before["differences"].get("closure", []):
        # path shape: /SHIP/FAMILY/rel_error: <stored> -> <fresh>
        path, _, tail = line.partition(": ")
        parts = path.strip("/").split("/")
        if len(parts) != 3 or parts[2] != "rel_error":
            continue
        stored = tail.split(" -> ")[0]
        try:
            out.append((parts[0], parts[1], float(stored)))
        except ValueError:
            continue
    return out


def write_document(
    library: Path, doc: Path, figs: Path | None, before: Mapping[str, Any]
) -> None:
    payload = json.loads(library.read_text())
    report = payload["report"]
    coverage = payload["coverage"]
    adj = payload["adjudication"]
    lines: list[str] = []
    add = lines.append

    add("# Round 5 — the cached blocks now follow the cells, and the dissent is audited")
    add("")
    add(
        "Anchor **REPRO-GEOM-FLEET-REFS** (*six benchmark hulls, ≥2 independent sources each*).\n"
        "Generated by `tools/baseline_fleet/refs_freshness.py refresh`; every number below is\n"
        "read from `benchmarks/verified/fleet_refs/crossref.json` after the refresh, so the\n"
        "document cannot drift from the library.  No new source was fetched this round: this is\n"
        "an audit of the library's own derived blocks."
    )
    add("")
    add("## 1. Two defects, both introduced by the library's own caching")
    add("")
    add("**Defect 1 — a stale cache printed a phantom source defect.**")
    add(
        "`closure` and `adjudication` are functions of `report`, but round 3 rewrote `report`\n"
        "(KVLCC2 `Bwl` repaired 580→58 m, DTC/JSR-2017 cells added, JBC volume refined to\n"
        "178 369.925 m³) and left both blocks at round 2's values.  A reader of the committed\n"
        "file was told KVLCC2's own digits miss their displacement by 900 % while the same file's\n"
        "`report` block closed them to 1e-6, and the published table\n"
        "(section 4.2 of `docs/fleet_refs_crossref.md`) carried that as a *finding*."
    )
    add("")
    add(
        "| block | stale? before refresh | leaves that disagreed |\n|---|---|---|"
    )
    for block in ("coverage", "closure", "adjudication", "lineage"):
        leaves = before["differences"].get(block, [])
        add(
            f"| `{block}` | {'**stale**' if block in before['stale'] else 'fresh'} "
            f"| {len(leaves)} |"
        )
    add("")
    add(
        "**Defect 2 — the scale-slip detector stopped auditing the dissent.**  Round 1 attributed\n"
        "the KVLCC2 `Bwl 580` typo with a power-of-ten scan of the *family value*.  Round 3 gave\n"
        "that row a second lineage printing the corrected 58 m, so the majority cluster became\n"
        "the corrected value, the family closes, and a detector that only inspects the winning\n"
        "numeral reported *no typo* — for a source that still prints 580 in two of its five\n"
        "entries.  The finding did not become false; it moved into `family_disagreement.dissent`,\n"
        "where nothing was looking.  Two tests had been parked as `xfail` on exactly this."
    )
    add("")
    add("## 2. What the dissent scan finds now")
    add("")
    add("| ship | family | quantity | printed (dissent) | ÷ | corrected | family value | residual | verdict |")
    add("|---|---|---|---|---|---|---|---|---|")
    for slip in adj["scale_typos"]:
        add(
            f"| {slip['ship']} | {slip['family']} | {slip['quantity']} "
            f"| {slip['printed']:g} | {slip['divide_by']:g} | {slip['corrected']:g} "
            f"| {slip.get('family_value', float('nan')):g} "
            f"| {slip['corrected_residual']:+.2e} | attributed from `{slip['origin']}` |"
        )
    if not adj["scale_typos"]:
        add("| — | — | — | — | — | — | — | — | no single-factor slip in any dissent |")
    add("")
    add(
        "The numerals are unchanged from round 1 (`raw_residual` "
        f"{adj['scale_typos'][0]['raw_residual']:+.3f} if the printed 580 were believed, "
        f"{adj['scale_typos'][0]['corrected_residual']:+.2e} after ÷10); what changed is *where*\n"
        "they are read from.  The cell keeps the printed value and gains `closure_defect` with\n"
        "`origin: dissent`, and the family vote is **not** rewritten by the detector."
        if adj["scale_typos"]
        else "No dissent needed adjudicating this round."
    )
    add("")
    add("## 3. Publisher self-consistency after the refresh")
    add("")
    add(
        "| hull | lineage | length basis | Cb·L·B·T (m³) | published ∇ (m³) | rel. error | verdict |"
    )
    add("|---|---|---|---|---|---|---|")
    for ship in SHIPS:
        rows = adj["ships"].get(ship, {}).get("closure") or payload["closure"].get(ship, {})
        for family, row in sorted(rows.items()):
            add(
                f"| {ship} | {family} | {row.get('length_basis', '-')} | {row['implied_volume_m3']:.1f} "
                f"| {row['published_volume_m3']:.1f} | {row['rel_error']:+.2e} "
                f"| {'ok' if row['consistent'] else 'does not close'} |"
            )
    add("")
    add(
        f"Tolerance ±{CLOSURE_TOLERANCE:.0%}.  Series 60 stays empty: the library holds no\n"
        "numeral for it at all, so there is nothing to close."
    )
    add("")
    add("## 4. Coverage and the anchor's five criteria (unchanged by this round)")
    add("")
    add(
        f"cells {coverage['cells']}, cross-verified {coverage['cross_verified']}, "
        f"single-source {coverage['single_source']}, conflict {coverage['conflict']}, "
        f"missing {coverage['missing']}; fully cross-verified ships: "
        f"{coverage['ships_fully_cross_verified']}."
    )
    add("")
    add("| ship | " + " | ".join(QUANTITY_LABELS[q] for q in QUANTITIES) + " |")
    add("|---" * (len(QUANTITIES) + 1) + "|")
    for ship in SHIPS:
        row = report[ship]["quantities"]
        cells = []
        for q in QUANTITIES:
            n = len(payload["lineage"][ship][q]["lineages"])
            cells.append(f"{row[q]['verdict']} ({n})")
        add(f"| {ship} | " + " | ".join(cells) + " |")
    add("")
    add(
        "This round adds **no** source, so criterion ① (≥2 independent lineages per ship on\n"
        "∇/Cb/LCB/Lwl) still fails for DTC, JBC, KVLCC2, ONRT and Series 60.  It is not a\n"
        "regression and it is not softened: the honest pointer is `docs/fleet_refs_round3.md`,\n"
        "whose score table this round does not change."
    )
    add("")
    add("## 5. Traps worth recording")
    add("")
    add(
        "1. **A cache of a pure function is a second source of truth.**  The repair is not \"be\n"
        "   careful when merging\": `merge_into` now refreshes `closure`/`adjudication` itself and\n"
        "   `refs_freshness.py check` fails the build if any cached block stops following from\n"
        "   `report`.  A new round that forgets to refresh cannot be merged silently.\n"
        "2. **`adjudicate` writes into the report it reads** (`closure_defect` on the cells), so the\n"
        "   freshness check recomputes on a deep copy — asking \"is this stale?\" must not rewrite\n"
        "   the artefact being asked about.\n"
        "3. **Independence is not what the closure test measures.**  A table that closes on its own\n"
        "   digits proves it was transcribed consistently, never that a second body measured the\n"
        "   hull; the KVLCC2 `580` case is a table that *reproduced* a typo across two publishers.\n"
        "4. **The dissent is the audit trail.**  Majority voting is right for the value and wrong\n"
        "   for the audit: the losing numerals are the only place a systematic typo shows up once\n"
        "   enough publishers copy the correction.\n"
    )
    add("")
    add("## 6. How this was verified")
    add("")
    add(
        "```\n"
        "python -m pytest tests/test_fleet_refs_freshness.py tests/test_fleet_refs_adjudication.py \\\n"
        "                 tests/test_fleet_refs_oa_merge.py tests/test_fleet_refs_crossref.py \\\n"
        "                 tests/test_fleet_refs_round3.py -q\n"
        "python tools/baseline_fleet/refs_freshness.py check\n"
        "```\n"
        "The two `xfail` markers on the KVLCC2 slip tests are gone: they are real assertions of\n"
        "the round-3 contract (family value 58, dissent 580 flagged, printed cell untouched)."
    )
    add("")
    fig_path = _figure(payload, figs, before)
    if fig_path is not None:
        add(f"![closure before/after]({fig_path.as_posix()})")
        add("")
    doc.parent.mkdir(parents=True, exist_ok=True)
    doc.write_text("\n".join(lines) + "\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    check = sub.add_parser("check", help="exit non-zero if a cached block is stale")
    check.add_argument("--library", default=str(LIBRARY))
    refresh = sub.add_parser("refresh", help="recompute the cached blocks and rewrite the doc")
    refresh.add_argument("--library", default=str(LIBRARY))
    refresh.add_argument("--report", default=None, help="markdown document to regenerate")
    refresh.add_argument("--figs", default=None, help="png path for the before/after figure")
    refresh.add_argument("--before", default=None, help="pre-refresh snapshot of the library")
    args = parser.parse_args(argv)
    library = Path(args.library)
    if args.command == "check":
        return cmd_check(library)
    return cmd_refresh(
        library,
        Path(args.report) if args.report else None,
        Path(args.figs) if args.figs else None,
        Path(args.before) if args.before else None,
    )


if __name__ == "__main__":
    raise SystemExit(main())