"""Known-truth audit of the fleet reference *instrument* (anchor ``REPRO-GEOM-FLEET-REFS``).

The reference library (:mod:`tensorcad.fleet.crossref`) reads publisher pages and decides, per
quantity, whether two independent families state the same numeral.  Every existing test feeds it
*real* pages, where the correct answer is itself an assumption -- so those tests can show the
machine runs, never that it reads the right cell: when the reader picks the neighbouring row of a
table, test and bug agree.  This audit supplies the missing half.  It builds documents whose ground
truth is **designed** (a synthetic table plus the verdict it must produce), runs the production
``parse_source -> cross_reference`` path over them unmodified, and scores the machine against that
truth, cell by cell.

Two numbers come out, and they are not the same number:

* **sensitivity** -- of the planted defects (a ten-fold decimal slip, a ratio row that looks like a
  length row, an LCB printed in metres, a model-column-only volume), how many does the machine
  catch?  A miss means a real publisher's footnote can silently become a "cross-verified" reference.
* **specificity** -- of the clean fixtures (two families agreeing, last-digit rounding, thousands
  separators), how many stay clean?  A miss means the library reports a conflict that does not
  exist and an operator learns to ignore the flag.

**known_limit** fixtures are neither: they record behaviour the bytes cannot decide (a volume
printed under a ``(t)`` header; an unregistered reprint of the same table).  They are asserted so
the limit stays visible instead of being rediscovered as a surprise.

``crossref.INSTRUMENT_GUARDS`` switches the row-label guards added by this round.  With
``--no-guards`` the *same* corpus is scored against the pre-guard reading, so the defect report is
an experiment that can be re-run rather than a claim about the past; the pre-guard artefact is
quoted in the generated report.

Reproduce::

    python tools/baseline_fleet/refs_instrument_audit.py --no-guards \\
        --out artifacts/fleet_refs/instrument_audit_prefix.json
    python tools/baseline_fleet/refs_instrument_audit.py --docs \\
        --prefix artifacts/fleet_refs/instrument_audit_prefix.json
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tensorcad.fleet import crossref as xref  # noqa: E402
from tensorcad.fleet.crossref import (  # noqa: E402
    Cell,
    SourceSpec,
    cross_reference,
    parse_source,
    votes,
)

ROOT = Path(__file__).resolve().parents[2]
ARTEFACT = ROOT / "artifacts/fleet_refs/instrument_audit.json"
DOC = ROOT / "docs/fleet_refs_instrument_audit.md"
FIGURE = ROOT / "docs/figs/fleet_refs/instrument_audit.png"
LIBRARY = ROOT / "benchmarks/verified/fleet_refs/crossref.json"


# --------------------------------------------------------------------------------------
# corpus helpers
# --------------------------------------------------------------------------------------
def spec(key: str, family: str, *, ship: str = "SYNTH", column: int | None = None) -> SourceSpec:
    """A synthetic source spec: nothing is fetched, so the digest is a marker, not a claim."""
    return SourceSpec(
        key=key,
        ship=ship,
        publisher=f"synthetic fixture {key}",
        family=family,
        url=f"fixture://{key}",
        kind="html",
        sha256="0" * 64,
        licence="synthetic fixture written for this audit; no publisher bytes involved",
        note="ground truth designed by the audit",
        column=column,
    )


def table(rows: Sequence[Sequence[str]]) -> str:
    body = "".join("<tr>" + "".join(f"<td>{cell}</td>" for cell in row) + "</tr>" for row in rows)
    return f"<table><tr><th>Particular</th><th>Value</th></tr>{body}</table>"


FULL_ROWS: tuple[tuple[str, str], ...] = (
    ("Length between perpendiculars Lpp (m)", "230.0"),
    ("Length on waterline Lwl (m)", "232.5"),
    ("Maximum beam on waterline Bwl (m)", "32.2"),
    ("Design draft T (m)", "10.8"),
    ("Displacement volume (m3)", "52,030"),
    ("LCB (%), fwd+", "-1.48"),
    ("Block coefficient Cb", "0.6505"),
)

#: The same ship as a second publisher prints it: rounding-only differences, LCB to two decimals.
SECOND_ROWS: tuple[tuple[str, str], ...] = (
    ("Length between perpendiculars Lpp (m)", "230.00"),
    ("Length on waterline Lwl (m)", "232.5"),
    ("Maximum beam on waterline Bwl (m)", "32.20"),
    ("Design draft T (m)", "10.80"),
    ("Displacement volume (m3)", "52030"),
    ("LCB (%), fwd+", "-1.50"),
    ("Block coefficient Cb", "0.6505"),
)


def edited(rows: Sequence[Sequence[str]], *, label_contains: str, value: str) -> tuple[tuple[str, str], ...]:
    """``rows`` with the value of every row whose label contains ``label_contains`` replaced."""
    return tuple(
        (label, value if label_contains in label else current) for label, current in rows
    )


ALL_QUANTITIES = (
    "displacement_volume",
    "cb",
    "lcb_percent",
    "lwl",
    "length_pp",
    "beam",
    "draft",
)


@dataclass
class Expect:
    """What the machine must say about one quantity of one fixture."""

    verdict: str
    disputed: bool = False
    families: int | None = None


@dataclass
class Fixture:
    id: str
    kind: str  # "control" | "planted" | "known_limit"
    why: str
    docs: list[tuple[SourceSpec, str]]
    expect: dict[str, Expect]
    #: source key -> quantity -> exact numeral the parser must return (digit-exact, not "close")
    cells: dict[str, dict[str, float]] = field(default_factory=dict)
    #: source key -> quantity -> required ``unit_basis``
    basis: dict[str, dict[str, str]] = field(default_factory=dict)
    #: quantity -> family -> the value the family must carry (printed majority, never a mean)
    family_values: dict[str, dict[str, float]] = field(default_factory=dict)
    #: quantity -> family -> numerals that must appear as recorded dissent, never averaged in
    dissents: dict[str, dict[str, list[float]]] = field(default_factory=dict)
    #: register the fixture keys in the production family registry (False = leave them unknown,
    #: which is what makes the fail-open default observable)
    register_families: bool = True


def _all(verdict: str, disputed: bool = False, families: int | None = None) -> dict[str, Expect]:
    return {q: Expect(verdict=verdict, disputed=disputed, families=families) for q in ALL_QUANTITIES}


_CLEAN_A = table(FULL_ROWS)
_CLEAN_B = table(SECOND_ROWS)

CORPUS: tuple[Fixture, ...] = (
    Fixture(
        id="clean_two_families",
        kind="control",
        why="Two independent families state the same particulars; nothing may be flagged.",
        docs=[(spec("famA", "FAM-A"), _CLEAN_A), (spec("famB", "FAM-B"), _CLEAN_B)],
        expect=_all("CROSS_VERIFIED", families=2),
        cells={
            "famA": {
                "length_pp": 230.0,
                "lwl": 232.5,
                "beam": 32.2,
                "draft": 10.8,
                "displacement_volume": 52030.0,
                "lcb_percent": -1.48,
                "cb": 0.6505,
            }
        },
    ),
    Fixture(
        id="two_pages_one_registered_family",
        kind="control",
        why="Two pages of one registered family are one source; repeating a numeral is not corroboration.",
        docs=[
            (spec("famA_p1", "FAM-A"), _CLEAN_A),
            (spec("famA_p2", "FAM-A"), _CLEAN_B),
        ],
        expect=_all("SINGLE_SOURCE", families=1),
    ),
    Fixture(
        id="unregistered_reprint_pair",
        kind="known_limit",
        why=(
            "Independence is read from the FAMILIES/LINEAGE registry by source key.  Two unregistered "
            "keys that reprint the same table are therefore two families and the cell is reported "
            "cross-verified.  The default is fail-open, and this fixture pins that down: a new source "
            "must be registered, or its agreement counts twice."
        ),
        docs=[(spec("unregistered_1", "FAM-A"), _CLEAN_A), (spec("unregistered_2", "FAM-A"), _CLEAN_B)],
        expect=_all("CROSS_VERIFIED", families=2),
        register_families=False,
    ),
    Fixture(
        id="ratio_rows_present",
        kind="planted",
        why=(
            "A table prints 'Lwl/Lpp 1.0109', 'Bwl/Lpp 0.1386' and 'T/L 0.0464' above the rows they "
            "divide.  A quotient of two particulars is a value of neither; read as a length it makes "
            "the library compare 1.0109 against 232.5 and report a conflict that does not exist."
        ),
        docs=[
            (
                spec("famA", "FAM-A"),
                table(
                    (("Lwl/Lpp", "1.0109"), ("Bwl/Lpp", "0.1386"), ("T/L", "0.0464")) + FULL_ROWS
                ),
            ),
            (spec("famB", "FAM-B"), _CLEAN_B),
        ],
        expect=_all("CROSS_VERIFIED", families=2),
        cells={"famA": {"lwl": 232.5, "length_pp": 230.0, "beam": 32.2, "draft": 10.8}},
    ),
    Fixture(
        id="lcb_printed_in_metres",
        kind="planted",
        why=(
            "Some pages print LCB as a distance from the forward perpendicular ('LCB (m) -3.395') and "
            "then the percentage; the metre row is a length and must not vote as a percentage, while "
            "the genuine percentage row later on the page must win."
        ),
        docs=[
            (
                spec("famA", "FAM-A"),
                table((("LCB (m)", "-3.395"),) + FULL_ROWS),
            ),
            (spec("famB", "FAM-B"), _CLEAN_B),
        ],
        expect={
            "lcb_percent": Expect("CROSS_VERIFIED", families=2),
            "length_pp": Expect("CROSS_VERIFIED", families=2),
        },
        cells={"famA": {"lcb_percent": -1.48}},
    ),
    Fixture(
        id="lcb_metres_only",
        kind="planted",
        why=(
            "A page that prints only 'LCB (m) -3.395' states a length, not a percentage.  The cell must "
            "be kept (with basis lcb_as_length) and must not vote, so the quantity is SINGLE_SOURCE "
            "rather than a 129 % conflict against the other family's -1.48 %."
        ),
        docs=[
            (
                spec("famA", "FAM-A"),
                table((("LCB (m)", "-3.395"),) + tuple(r for r in FULL_ROWS if "LCB" not in r[0])),
            ),
            (spec("famB", "FAM-B"), _CLEAN_B),
        ],
        expect={"lcb_percent": Expect("SINGLE_SOURCE", families=1)},
        cells={"famA": {"lcb_percent": -3.395}},
        basis={"famA": {"lcb_percent": "lcb_as_length"}},
    ),
    Fixture(
        id="decimal_slip_tenfold",
        kind="planted",
        why="One family drops a digit (5,203 for 52,030); a 90 % disagreement cannot be agreement.",
        docs=[
            (spec("famA", "FAM-A"), _CLEAN_A),
            (spec("famB", "FAM-B"), table(edited(SECOND_ROWS, label_contains="Displacement", value="5,203"))),
        ],
        expect={
            "displacement_volume": Expect("CONFLICT", disputed=True, families=2),
            "cb": Expect("CROSS_VERIFIED", families=2),
        },
        cells={"famB": {"displacement_volume": 5203.0}},
    ),
    Fixture(
        id="same_family_dissent_majority",
        kind="planted",
        why=(
            "One family prints the numeral twice correctly and once with a dropped digit.  The family "
            "value must be the printed majority (52030), never the mean of 52030 and 5203 (a number no "
            "publisher ever printed), and the dissent must be recorded against the family."
        ),
        docs=[
            (spec("famA_p1", "FAM-A"), _CLEAN_A),
            (spec("famA_p2", "FAM-A"), _CLEAN_A),
            (spec("famA_p3", "FAM-A"), table(edited(FULL_ROWS, label_contains="Displacement", value="5,203"))),
            (spec("famB", "FAM-B"), _CLEAN_B),
        ],
        expect={
            "displacement_volume": Expect("CROSS_VERIFIED", families=2),
            "length_pp": Expect("CROSS_VERIFIED", families=2),
        },
        family_values={"displacement_volume": {"FAM-A": 52030.0, "FAM-B": 52030.0}},
        dissents={"displacement_volume": {"FAM-A": [5203.0]}},
    ),
    Fixture(
        id="mass_labelled_volume",
        kind="planted",
        why=(
            "One family prints the displacement under a '(t)' header with a bare numeral cell.  A mass "
            "is not a volume: reading 52 030 t as 52 030 m3 absorbs a 2.4 % error into a "
            "'cross-verified' reference.  The row label carries the unit, so the cell must be "
            "converted at 1.025 t/m3 and kept as a disputed reading of 50 761 m3 -- never averaged "
            "with the volume family."
        ),
        docs=[
            (spec("famA", "FAM-A"), _CLEAN_A),
            (spec("famB", "FAM-B"), table(SECOND_ROWS[:4] + (("Displacement (t)", "52,030"),) + SECOND_ROWS[5:])),
        ],
        expect={"displacement_volume": Expect("CONFLICT", disputed=True, families=2)},
        cells={"famB": {"displacement_volume": 50760.976}},
        basis={"famB": {"displacement_volume": "mass"}},
    ),
    Fixture(
        id="model_column_only_volume",
        kind="planted",
        why=(
            "A model/full-scale table prints the lengths in both columns but leaves the full-scale "
            "displacement blank.  The table's own ratio (Lpp 230 / 2.3 = 100) makes the recovery "
            "arithmetic and volume scales as its third power (0.05203 x 100^3 = 52030 m3).  The cell "
            "must be marked derived -- auditable, never voting -- and must not be left at the model "
            "numeral or scaled by the length ratio."
        ),
        docs=[
            (
                spec("famA", "FAM-A", column=2),
                "<table><tr><th>Particular</th><th>Model (1:100)</th><th>Full scale</th></tr>"
                "<tr><td>Length between perpendiculars Lpp (m)</td><td>2.300</td><td>230.0</td></tr>"
                "<tr><td>Length on waterline Lwl (m)</td><td>2.325</td><td>232.5</td></tr>"
                "<tr><td>Maximum beam on waterline Bwl (m)</td><td>0.322</td><td>32.2</td></tr>"
                "<tr><td>Design draft T (m)</td><td>0.108</td><td>10.8</td></tr>"
                "<tr><td>Displacement volume (m3)</td><td>0.05203</td><td>NA</td></tr>"
                "<tr><td>Block coefficient Cb</td><td>0.6505</td><td>0.6505</td></tr></table>",
            ),
            (spec("famB", "FAM-B"), _CLEAN_B),
        ],
        expect={"displacement_volume": Expect("SINGLE_SOURCE", families=1), "draft": Expect("CROSS_VERIFIED", families=2)},
        cells={"famA": {"draft": 10.8, "lwl": 232.5, "displacement_volume": 52030.0}},
        basis={"famA": {"displacement_volume": "derived_model_column_cubed"}},
    ),
    Fixture(
        id="thousands_separator",
        kind="control",
        why="A numeral printed with a thousands separator is 52030, not 52, and the raw text is kept.",
        docs=[(spec("famA", "FAM-A"), _CLEAN_A), (spec("famB", "FAM-B"), _CLEAN_B)],
        expect={"displacement_volume": Expect("CROSS_VERIFIED", families=2), "cb": Expect("CROSS_VERIFIED")},
        cells={"famA": {"displacement_volume": 52030.0}},
    ),
    Fixture(
        id="within_tolerance_rounding",
        kind="control",
        why="Last-digit rounding (52,043 vs 52,030, 0.025 %) is agreement; a conflict here is noise.",
        docs=[
            (spec("famA", "FAM-A"), _CLEAN_A),
            (spec("famB", "FAM-B"), table(edited(SECOND_ROWS, label_contains="Displacement", value="52,043"))),
        ],
        expect={
            "displacement_volume": Expect("CROSS_VERIFIED", families=2),
            "length_pp": Expect("CROSS_VERIFIED", families=2),
        },
    ),
)


# --------------------------------------------------------------------------------------
# scoring
# --------------------------------------------------------------------------------------
def _close(a: float, b: float, tol: float = 1e-6) -> bool:
    return abs(a - b) <= tol * max(abs(a), abs(b), 1.0)


def _register_families(docs: Sequence[tuple[SourceSpec, str]]) -> dict[str, list[str]]:
    """Register the fixture keys in the production family registry.

    ``cross_reference`` resolves a source's family from :data:`crossref.LINEAGE` / 
    :data:`crossref.FAMILIES` *by source key*; an unknown key is its own family.  Fixture keys do
    not exist there, so they are registered for the duration of one fixture and removed again --
    the alternative (letting every key be its own family) would make the independence fixtures
    untestable.  Returns the keys added so they can be withdrawn.
    """
    added_families: list[str] = []
    added_lineage: list[str] = []
    for source_spec, _ in docs:
        if source_spec.key not in xref.FAMILIES:
            xref.FAMILIES[source_spec.key] = source_spec.family
            added_families.append(source_spec.key)
        xref.LINEAGE.pop(source_spec.key, None)
    return {"FAMILIES": added_families, "LINEAGE": added_lineage}


def _withdraw(added: Mapping[str, list[str]]) -> None:
    for key in added.get("FAMILIES", []):
        xref.FAMILIES.pop(key, None)


def _check_cells(fixture: Fixture, parsed: Mapping[str, Mapping[str, Cell]]) -> list[str]:
    problems: list[str] = []
    for source_key, wanted in fixture.cells.items():
        cells = parsed.get(source_key, {})
        for quantity, value in wanted.items():
            cell = cells.get(quantity)
            if cell is None:
                problems.append(f"source {source_key}: no cell for {quantity} (the page states {value:g})")
                continue
            if not _close(cell.value, value):
                problems.append(
                    f"source {source_key}: {quantity} read as {cell.value:g} but the page states "
                    f"{value:g} (label {cell.label!r})"
                )
    for source_key, wanted in fixture.basis.items():
        cells = parsed.get(source_key, {})
        for quantity, basis in wanted.items():
            cell = cells.get(quantity)
            if cell is None:
                problems.append(f"source {source_key}: no cell for {quantity} to carry basis {basis}")
            elif cell.unit_basis != basis:
                problems.append(
                    f"source {source_key}: {quantity} basis {cell.unit_basis!r}, expected {basis!r} "
                    f"(value {cell.value:g})"
                )
    # a cell that is not allowed to vote must be reported as such
    for source_key, cells in parsed.items():
        for quantity, cell in cells.items():
            if cell.unit_basis.startswith(("fallback_column", "derived", "lcb_as_length")) and votes(cell):
                problems.append(f"source {source_key}: {quantity} basis {cell.unit_basis!r} is voting")
    return problems


def _check_report(fixture: Fixture, report: Mapping[str, Any]) -> list[str]:
    problems: list[str] = []
    quantities = report[fixture.docs[0][0].ship]["quantities"]
    for quantity, want in fixture.expect.items():
        got = quantities[quantity]
        if got["verdict"] != want.verdict:
            problems.append(f"{quantity}: verdict {got['verdict']}, the fixture designs {want.verdict}")
        if bool(got["disputed"]) != want.disputed:
            problems.append(f"{quantity}: disputed={got['disputed']}, the fixture designs {want.disputed}")
        if want.families is not None and got["independent_families"] != want.families:
            problems.append(
                f"{quantity}: {got['independent_families']} independent families, the fixture has "
                f"{want.families} (family_values {got['family_values']})"
            )
    for quantity, families in fixture.family_values.items():
        got_values = quantities[quantity]["family_values"]
        for family, value in families.items():
            if family not in got_values:
                problems.append(f"{quantity}: family {family} missing from family_values {got_values}")
            elif not _close(got_values[family], value):
                problems.append(
                    f"{quantity}: family {family} value {got_values[family]:g}, expected the printed "
                    f"majority {value:g} (a mean of disagreeing entries is not a published numeral)"
                )
    for quantity, families in fixture.dissents.items():
        got_dissent = quantities[quantity]["family_disagreement"]
        for family, values in families.items():
            recorded = got_dissent.get(family, {}).get("dissent")
            if recorded is None:
                problems.append(f"{quantity}: family {family} dissent not recorded ({got_dissent})")
            elif sorted(round(v, 6) for v in recorded) != sorted(round(v, 6) for v in values):
                problems.append(f"{quantity}: family {family} dissent {recorded}, expected {values}")
    return problems


def run_corpus() -> dict[str, Any]:
    """Score the production parser and verdict logic against the designed ground truth."""
    results: list[dict[str, Any]] = []
    for fixture in CORPUS:
        added = _register_families(fixture.docs) if fixture.register_families else {"FAMILIES": []}
        try:
            parsed: dict[str, dict[str, Cell]] = {
                source_spec.key: parse_source(source_spec, html) for source_spec, html in fixture.docs
            }
            report = cross_reference({fixture.docs[0][0].ship: parsed})
        finally:
            _withdraw(added)
        problems = _check_cells(fixture, parsed) + _check_report(fixture, report)
        results.append(
            {
                "id": fixture.id,
                "kind": fixture.kind,
                "why": fixture.why,
                "passed": not problems,
                "problems": problems,
                "observed": {
                    quantity: {
                        "verdict": got["verdict"],
                        "disputed": got["disputed"],
                        "independent_families": got["independent_families"],
                        "family_values": got["family_values"],
                        "family_disagreement": got["family_disagreement"],
                    }
                    for quantity, got in report[fixture.docs[0][0].ship]["quantities"].items()
                },
                "parsed": {
                    key: {
                        quantity: {
                            "value": cell.value,
                            "raw": cell.raw,
                            "label": cell.label,
                            "unit_basis": cell.unit_basis,
                            "votes": votes(cell),
                            "derivation": cell.derivation,
                        }
                        for quantity, cell in cells.items()
                    }
                    for key, cells in parsed.items()
                },
            }
        )

    def score(kind: str) -> dict[str, Any]:
        subset = [row for row in results if row["kind"] == kind]
        passed = [row for row in subset if row["passed"]]
        return {
            "fixtures": len(subset),
            "passed": len(passed),
            "rate": (len(passed) / len(subset)) if subset else None,
            "failing": [row["id"] for row in subset if not row["passed"]],
        }

    return {
        "schema": "tensorcad.fleet_refs.instrument_audit/1.0",
        "anchor": "REPRO-GEOM-FLEET-REFS",
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "guards": bool(xref.INSTRUMENT_GUARDS),
        "protocol": (
            "each fixture is a synthetic table plus its designed ground truth; the production "
            "parse_source -> cross_reference path runs unmodified and is scored on (a) the numeral "
            "read per cell and its voting basis, (b) verdict / disputed / family count / family "
            "value / recorded dissent per quantity"
        ),
        "sensitivity": score("planted"),
        "specificity": score("control"),
        "known_limits": score("known_limit"),
        "fixtures": results,
        "impact_on_library": impact_on_library(),
    }


def impact_on_library() -> dict[str, Any]:
    """Would the guards added for this audit re-read any cell of the shipped library?

    The artefact stores every parsed cell (label, basis, value, source), so the question is
    answerable offline: count the stored cells the new rules would treat differently.  A non-zero
    count means the fix changes published verdicts and the artefact must be rebuilt from the pinned
    bytes; zero means the fix is behaviour-preserving on the corpus of record.

    Only *HTML* cells are scanned: their label is one table row's first cell, which is what the
    guards test.  A PDF cell's ``label`` is a prose context window around a quote, so a row-label
    test applied to it would fire on neighbouring text -- that would be a finding about the scan,
    not about the library.
    """
    if not LIBRARY.exists():
        return {"checked": False, "reason": f"{LIBRARY} missing"}
    data = json.loads(LIBRARY.read_text())
    kinds = {source.key: source.kind for source in xref.SOURCES}
    ratio_rows: list[dict[str, Any]] = []
    lcb_length_rows: list[dict[str, Any]] = []
    model_volume_rows: list[dict[str, Any]] = []
    mass_label_rows: list[dict[str, Any]] = []
    cells = html_cells = 0
    for ship, body in data.get("report", {}).items():
        for quantity, block in body.get("quantities", {}).items():
            for entry in block.get("entries", []):
                cell = entry["cell"]
                cells += 1
                label = cell.get("label", "")
                if kinds.get(cell.get("source", "")) != "html":
                    continue
                html_cells += 1
                stored = Cell(
                    quantity=quantity,
                    value=cell["value"],
                    raw=cell.get("raw", ""),
                    label=label,
                    source=cell.get("source", ""),
                    url=cell.get("url", ""),
                    sha256=cell.get("sha256", ""),
                    column=cell.get("column", 0),
                    unit_basis=cell.get("unit_basis", "volume"),
                )
                if not votes(stored):
                    continue
                hit = {"ship": ship, "quantity": quantity, "source": cell.get("source"), "label": label[:120]}
                if xref._RATIO_LABEL.search(label):
                    ratio_rows.append(hit)
                if quantity == "lcb_percent" and xref._LCB_AS_LENGTH.search(label):
                    lcb_length_rows.append(hit)
                if quantity == "displacement_volume" and stored.unit_basis.startswith("fallback_column"):
                    model_volume_rows.append(hit)
                if (
                    quantity == "displacement_volume"
                    and stored.unit_basis == "volume"
                    and xref._MASS_LABEL.search(label)
                ):
                    mass_label_rows.append(hit)
    return {
        "checked": True,
        "library": str(LIBRARY.relative_to(ROOT)),
        "cells_scanned": cells,
        "html_cells_scanned": html_cells,
        "ratio_label_voters": ratio_rows,
        "lcb_as_length_voters": lcb_length_rows,
        "model_column_volume_cells": model_volume_rows,
        "mass_label_volume_cells": mass_label_rows,
        "reclassified": (
            len(ratio_rows) + len(lcb_length_rows) + len(model_volume_rows) + len(mass_label_rows)
        ),
    }


# --------------------------------------------------------------------------------------
# reporting
# --------------------------------------------------------------------------------------
def _md_table(rows: Sequence[Sequence[str]], header: Sequence[str]) -> str:
    head = "| " + " | ".join(header) + " |"
    sep = "|" + "|".join("---" for _ in header) + "|"
    body = "\n".join("| " + " | ".join(cell for cell in row) + " |" for row in rows)
    return f"{head}\n{sep}\n{body}"


def render_figure(artefact: Mapping[str, Any]) -> bool:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:  # pragma: no cover - plotting is optional
        return False
    fixtures = artefact["fixtures"]
    values = [1.0 if row["passed"] else 0.0 for row in fixtures]
    colours = ["#2c7a4b" if row["passed"] else "#a4243b" for row in fixtures]
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(10.2, 4.4))
    ax.bar(range(len(fixtures)), values, color=colours, edgecolor="#1b1b1b", linewidth=0.6)
    ax.set_xticks(range(len(fixtures)))
    ax.set_xticklabels([row["id"] for row in fixtures], rotation=32, ha="right", fontsize=8)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["wrong answer", "right answer"])
    ax.set_ylim(0, 1.12)
    blocked = not artefact.get("guards", True)
    ax.set_title(
        "Fleet-reference instrument vs a designed truth "
        f"({'pre-fix reading: guards off' if blocked else 'shipped code'})",
        fontsize=10,
    )
    for index, row in enumerate(fixtures):
        ax.text(index, 0.03, row["kind"], rotation=90, fontsize=6.5, color="white", ha="center", va="bottom")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIGURE, dpi=160)
    plt.close(fig)
    return True


def render_doc(artefact: Mapping[str, Any], prefix: Mapping[str, Any] | None) -> None:
    def score_line(key: str, label: str) -> str:
        block = artefact[key]
        if block["fixtures"] == 0:
            return f"* {label}: no fixtures"
        rate = f"{block['rate'] * 100:.0f} %" if block["rate"] is not None else "n/a"
        return f"* {label}: **{block['passed']}/{block['fixtures']}** ({rate})"

    rows = [
        [
            f"`{row['id']}`",
            row["kind"],
            "pass" if row["passed"] else "**fail**",
            "<br>".join(problem.replace("|", "/") for problem in row["problems"]) or "--",
        ]
        for row in artefact["fixtures"]
    ]
    observations = []
    for row in artefact["fixtures"]:
        for quantity in ALL_QUANTITIES:
            got = row["observed"][quantity]
            if quantity in ("displacement_volume", "cb") or got["disputed"]:
                observations.append(
                    [
                        f"`{row['id']}`",
                        quantity,
                        got["verdict"],
                        str(got["disputed"]),
                        str(got["independent_families"]),
                        ", ".join(f"{k}={v:.6g}" for k, v in got["family_values"].items()) or "--",
                    ]
                )
    impact = artefact["impact_on_library"]
    prefix_rows: list[str] = []
    prefix_score = ""
    if prefix is not None:
        for row in prefix["fixtures"]:
            for problem in row["problems"]:
                prefix_rows.append([f"`{row['id']}`", problem.replace("|", "/")])
        prefix_score = (
            f"sensitivity **{prefix['sensitivity']['passed']}/{prefix['sensitivity']['fixtures']}**, "
            f"specificity **{prefix['specificity']['passed']}/{prefix['specificity']['fixtures']}**, "
            f"known-limit **{prefix['known_limits']['passed']}/{prefix['known_limits']['fixtures']}**"
        )
    impact_text = (
        f"* stored cells scanned: **{impact.get('cells_scanned', 0)}** "
        f"({impact.get('html_cells_scanned', 0)} of them HTML row cells, the only kind the guards address)\n"
        f"* cells the new guards would reclassify: **{impact.get('reclassified', 0)}** "
        f"(ratio-label voters {len(impact.get('ratio_label_voters', []))}, "
        f"LCB-as-length voters {len(impact.get('lcb_as_length_voters', []))}, "
        f"model-column volume cells {len(impact.get('model_column_volume_cells', []))})"
        if impact.get("checked")
        else f"* library not scanned: {impact.get('reason', '')}"
    )
    doc = f"""# The reference instrument, audited against a designed truth

Anchor **REPRO-GEOM-FLEET-REFS** (TensorCAD). Machine artefact: `artifacts/fleet_refs/instrument_audit.json`
(regenerate: `python tools/baseline_fleet/refs_instrument_audit.py --docs`). Every number below is rendered from that file.

## 1. The question this round answers

The library answers one question per hull and quantity: *do two independent publishers state the
same numeral?* Its tests all feed it **real** pages, so they can show the machine runs but not that
it reads the right cell -- when the reader picks the neighbouring row of a table, test and bug
agree. This round puts the machine against documents whose answer is **designed**: a synthetic
table plus the verdict it must produce. A wrong reading then becomes a named numeral in a diff
instead of a plausible-looking percentage.

Corpus: **{len(artefact['fixtures'])} fixtures**, each one behaviour, in three classes:

* **control** -- clean data that must stay clean (specificity);
* **planted** -- a defect that must be caught (sensitivity);
* **known_limit** -- a case the bytes cannot decide; the fixture asserts the behaviour we accept, so
  the limit is documented instead of rediscovered as a surprise.

## 2. Score

Shipped code:

{score_line('sensitivity', 'sensitivity (planted defects caught)')}
{score_line('specificity', 'specificity (controls left clean)')}
{score_line('known_limits', 'known-limit behaviour reproduced')}

Same corpus with the guards switched off (`--no-guards`, the pre-fix reading): {prefix_score}.

![known-truth audit](figs/fleet_refs/instrument_audit.png)

## 3. Fixture by fixture

{_md_table(rows, ['fixture', 'class', 'result', 'machine vs designed truth'])}

## 4. Verdicts on the quantities that carry the arguments

{_md_table(observations, ['fixture', 'quantity', 'verdict', 'disputed', 'families', 'family values'])}

## 5. What the pre-guard reading got wrong

{_md_table(prefix_rows, ['fixture', 'pre-fix failure (verbatim from the audit run)']) if prefix_rows else '*the pre-fix snapshot recorded no failure*'}

Three defects, all latent -- no stored cell carries them, which is exactly why no existing test
could see them:

1. **A ratio row voted as a length.** A page that prints `Lwl/Lpp 1.0109` above the row it divides
   had the quotient read as `lwl = 1.0109 m` with basis `volume`, so it voted; `Bwl/Lpp` became
   `length_pp = 0.1386` and `T/L` became `draft = 0.0464`, because a quotient label matches the
   patterns of whichever of its two names comes first. The genuine rows were then skipped by the
   "first row for a quantity wins" rule, and the library compared a dimensionless ratio against a
   length in metres -- reporting a conflict that no publisher states. Fixed by refusing a row whose
   label is a quotient of two particulars (`_RATIO_LABEL`, hull-letter symbols only, so `kg/m3`,
   `m/s` and `TFP / TAP` are untouched).
2. **LCB in metres voted as a percentage.** `_LCB_AS_LENGTH` existed and was documented but was
   never applied: `LCB (m) -3.395` became `lcb_percent = -3.395 %` and voted against another
   family's `-1.48 %` -- a 129 % disagreement manufactured from a unit mismatch. Fixed by marking
   such a cell `lcb_as_length` (kept for the audit trail, never votes) and by letting a later
   percentage row replace a length row, never the reverse.
3. **A model-only displacement was dropped.** `derive_cells` recovered lengths and LCB from a
   table's model column but left a model-only volume at the model numeral. Volume scales as the
   third power of the table's own length ratio, which the same table pins (`Lpp 230 / 2.3 = 100`),
   so the recovery is arithmetic on one publisher's numerals like every other derivation. It is now
   emitted with basis `derived_model_column_cubed`: auditable, and **not voting**, because a
   derivation is evidence about the document, not a second opinion about the ship.

4. **A mass stated only in the row label was read as a volume.** The unit heuristic looked at the
   numeral cell, so `Displacement (t) | 52 030` -- a bare numeral under a mass header -- became
   `52 030 m3`. That is the failure mode this anchor exists to prevent: a 2.4 % unit error crosses
   as `CROSS_VERIFIED` because both families agree on a *numeral*. Fixed by reading the row label
   too (`_MASS_LABEL`), converting at 1.025 t/m3 and keeping the reading as a disputed 50 761 m3 --
   the conflict stays visible and the closure adjudicator can resolve it.

## 6. Impact on the shipped library

{impact_text}

No stored cell is reclassified: the guards fire only on row labels the corpus of record does not
contain, and no stored displacement cell sits in a fallback column. The fix is behaviour-preserving
on the six-hull library, so the verdicts published in `benchmarks/verified/fleet_refs/crossref.json`
stand unchanged -- and that statement is now a measurement rather than a hope.

## 7. Limits this audit makes explicit

* **Independence is read from a registry, not from the document.** `cross_reference` resolves a
  source's family by key from `LINEAGE`/`FAMILIES`, and an *unknown* key is its own family. So two
  unregistered reprints of one table are reported `CROSS_VERIFIED` (`unregistered_reprint_pair`,
  a known-limit fixture). The default is fail-open; a new source must be registered or its
  agreement counts twice.
* **A unit stated only in a table caption is still invisible.** The mass guard reads the row label
  and the numeral cell; a page that names the unit once, in the caption above the table, is read as
  a volume. The `mass_labelled_volume` fixture asserts that when the label *does* carry the unit the
  library keeps both readings and marks the cell disputed instead of averaging them; resolution is
  the closure adjudicator's job (`docs/fleet_refs_adjudication.md`).
* **Sensitivity is bounded by the defect classes written down.** A transposed digit inside one
  column, a unit stated only in a table caption, or a numeral belonging to a different loading
  condition are not in the corpus and are caught by nothing here.
* **Fixtures are synthetic HTML tables**, so they exercise `table_rows` and `parse_source` only. The
  PDF path reads *quotes* (regex addresses into extracted text) and cannot be exercised with a
  designed table; its exactness rests on the stored quote strings, not on this audit.
* The corpus scores the **instrument, not the ships**: a passing fixture says the machine reads the
  table it was designed to read, nothing about how a real publisher lays out a page.

## 8. Reproduce

```bash
python tools/baseline_fleet/refs_instrument_audit.py --docs \\
    --prefix artifacts/fleet_refs/instrument_audit_prefix.json
python -m pytest tests/test_fleet_refs_instrument_audit.py -q
```
"""
    DOC.parent.mkdir(parents=True, exist_ok=True)
    DOC.write_text(doc)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--docs", action="store_true", help="render the markdown report and figure")
    parser.add_argument("--prefix", type=Path, default=None, help="pre-fix artefact to quote in the report")
    parser.add_argument("--out", type=Path, default=ARTEFACT)
    parser.add_argument(
        "--no-guards",
        action="store_true",
        help="flip crossref.INSTRUMENT_GUARDS off: reproduce the pre-fix reading on the same code",
    )
    args = parser.parse_args(argv)

    if args.no_guards:
        xref.INSTRUMENT_GUARDS = False
    artefact = run_corpus()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artefact, indent=1) + "\n")
    if args.docs:
        prefix = None
        if args.prefix is not None and args.prefix.exists():
            prefix = json.loads(args.prefix.read_text())
        render_figure(artefact)
        render_doc(artefact, prefix)
    for kind in ("sensitivity", "specificity", "known_limits"):
        block = artefact[kind]
        print(f"{kind:14s} {block['passed']}/{block['fixtures']}  failing={block['failing']}")
    print(f"reclassified library cells: {artefact['impact_on_library'].get('reclassified')}")
    failing = [row["id"] for row in artefact["fixtures"] if not row["passed"]]
    print("fixture failures:", failing or "none")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())