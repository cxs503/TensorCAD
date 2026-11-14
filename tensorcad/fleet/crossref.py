"""Cross-source reference library for the six benchmark hulls.

Anchor ``REPRO-GEOM-FLEET-REFS``.  The fleet benchmark
(:mod:`tensorcad.fleet.official`, ``benchmarks/verified/baseline_fleet``) compares a
reconstructed displacement against **one** published number per hull.  A single number
cannot tell a reconstruction error from a source typo, so this module keeps a *library* of
published principal particulars -- displacement volume ∇, block coefficient Cb, LCB (% Lpp,
forward positive) and length on the waterline Lwl, plus the Lpp/Bwl/T the coefficients are
defined on -- and asks, per quantity, whether at least two **independent publishers** state
it and whether they agree.

Design rules that the tests enforce
-----------------------------------
* Nothing is transcribed by hand.  :func:`parse_source` reads a *fetched byte string* whose
  SHA-256 is pinned in :data:`SOURCES`; every cell carries the raw cell text, the row it
  came from and the source digest, so a value can always be traced back to the bytes.
* Independence is a property of the *publisher family*, not of the URL.  Two pages served
  by the same workshop committee are one family; the same workshop series in two years is
  still one family unless the pages are demonstrably different editions (recorded in
  ``family``).  A quantity backed by one family is ``SINGLE_SOURCE`` even when three URLs
  repeat it -- that is the failure mode this anchor exists to expose.
* A disagreement between families is reported as ``CONFLICT`` with per-family values; it is
  never averaged away.
"""

from __future__ import annotations

import re
from dataclasses import asdict, dataclass, field
from typing import Any, Iterable, Mapping, Sequence

from .official import table_rows

#: Quantities the library tracks.  ``displacement_volume`` is always a volume in m^3;
#: a publisher that states a mass is converted at :data:`SEAWATER_DENSITY` and the raw
#: mass is kept in ``Cell.raw`` together with ``Cell.unit_basis == "mass"``.
QUANTITIES: tuple[str, ...] = (
    "displacement_volume",
    "cb",
    "lcb_percent",
    "lwl",
    "length_pp",
    "beam",
    "draft",
)

#: Human labels used in the report, in report order.
QUANTITY_LABELS: dict[str, str] = {
    "displacement_volume": "∇ (m³)",
    "cb": "Cb",
    "lcb_percent": "LCB (% Lpp, fwd+)",
    "lwl": "Lwl (m)",
    "length_pp": "Lpp (m)",
    "beam": "Bwl (m)",
    "draft": "T (m)",
}

#: Seawater density used only when a publisher states a *mass* (t) instead of a volume.
SEAWATER_DENSITY = 1.025

#: Agreement tolerance per quantity.  Published reproductions of the same principal
#: particular differ in the last printed digit; 0.5 % covers rounding, not a different ship.
TOLERANCE: dict[str, float] = {
    "displacement_volume": 0.005,
    "cb": 0.005,
    "lcb_percent": 0.02,  # percent-of-Lpp is a small number; 2 % of it is 0.03 %Lpp
    "lwl": 0.005,
    "length_pp": 0.005,
    "beam": 0.005,
    "draft": 0.005,
}

_NUMBER = re.compile(r"[-+]?\d+(?:\.\d+)?")

#: Row-label patterns.  Matched case-insensitively against the *whole* flattened row so a
#: page that prints "Maximum beam of waterline B WL (m)" and one that prints "Bwl (m)" are
#: both understood.  Order matters: the first pattern that matches wins, so LCB is tested
#: before Cb and Lwl before Lpp.
LABEL_PATTERNS: tuple[tuple[str, str], ...] = (
    ("lcb_percent", r"\blcb\b|longitudinal cent(?:er|re) of buoyancy"),
    ("lwl", r"\bl\s?wl\b|waterline length|length (?:of|on) (?:the )?water\s?line|length at (?:the )?waterline"),
    ("length_pp", r"\bl\s?pp\b|length between perpendiculars"),
    ("beam", r"\bb\s?wl\b|maximum beam|breadth"),
    ("draft", r"\bdraft\b|\bdraught\b|^\s*t\s*(?:\(|\b)"),
    ("displacement_volume", r"displacement"),
    ("cb", r"\bc\s?b\b|block coefficient"),
)


@dataclass(frozen=True)
class SourceSpec:
    """One fetched publisher document and how to read it."""

    key: str
    ship: str
    publisher: str
    #: Independence family.  Two sources from the same family never count as two sources.
    family: str
    url: str
    kind: str  # "html" | "pdf"
    #: SHA-256 of the fetched bytes (pinned after the first fetch, verified on every re-fetch).
    sha256: str
    licence: str
    note: str = ""
    #: Column holding the *full-scale* value; ``None`` = first numeric cell after the label.
    column: int | None = None
    #: For multi-ship pages: the sub-table heading that starts this ship's block.
    section: str | None = None


#: The document set.  Digests are filled in by ``tools/baseline_fleet/fetch_crossref.py --pin``.
SOURCES: tuple[SourceSpec, ...] = (
    SourceSpec(
        key="t2015_kcs",
        ship="KCS",
        publisher="NMRI Tokyo 2015 workshop (G2015 geometry pages)",
        family="NMRI-T2015",
        url="https://www.t2015.nmri.go.jp/kcs_gc.html",
        kind="html",
        sha256="68d9c24600cb699b16394fc5c32762ebc875c856d14ef8460e1d0ebf4e7c78a4",
        licence="page copyright NMRI; digests + extracted cells only, page bytes not redistributed",
        column=2,
    ),
    SourceSpec(
        key="t2015_jbc",
        ship="JBC",
        publisher="NMRI Tokyo 2015 workshop (G2015 geometry pages)",
        family="NMRI-T2015",
        url="https://www.t2015.nmri.go.jp/jbc_gc.html",
        kind="html",
        sha256="fc2ce9054a8b7842162bac94365853fa076da2527e8c7e455832f9710f58f875",
        licence="page copyright NMRI; digests + extracted cells only",
        column=2,
    ),
    SourceSpec(
        key="t2015_onrt",
        ship="ONRT",
        publisher="NMRI Tokyo 2015 workshop (G2015 geometry pages)",
        family="NMRI-T2015",
        url="https://www.t2015.nmri.go.jp/onrt_gc.html",
        kind="html",
        sha256="5e13427332de8be17b9dc2245dc5830300baf2e2340efcd8a80eee2269c43a11",
        licence="page copyright NMRI; digests + extracted cells only",
        column=3,
    ),
    SourceSpec(
        key="simman2020_summary",
        ship="KCS|KVLCC2|ONRT",
        publisher="SIMMAN 2020 workshop (KRISO/SNU) — Table 1/2/3 main particulars",
        family="SIMMAN-2020",
        url="https://simman2020.kr/contents/Ship_data_summary.php",
        kind="html",
        sha256="c1a668f6d2440c5d83a83c3616285e9288ef618dbc9a2feebd26627be3f4ebb3",
        licence="page copyright SIMMAN2020; digests + extracted cells only",
    ),
    SourceSpec(
        key="simman2014_kvlcc2",
        ship="KVLCC2",
        publisher="SIMMAN 2014 workshop (FORCE Technology) — KVLCC2 geometry page",
        family="SIMMAN-2014",
        url="https://simman2014.dk/ship-data/moeri-kvlcc2-tanker/geometry-and-conditions-moeri-kvlcc2-tanker/",
        kind="html",
        sha256="bfc9101e3a4afa4fbe02750de3f0ec830ccc948074265ad4e63184f5b67d310b",
        licence="page copyright FORCE Technology; digests + extracted cells only",
    ),
    SourceSpec(
        key="simman2014_kcs",
        ship="KCS",
        publisher="SIMMAN 2014 workshop (FORCE Technology) — MOERI container ship geometry page",
        family="SIMMAN-2014",
        url="https://simman2014.dk/ship-data/moeri-container-ship/geometry-and-conditions-moeri-container-ship/",
        kind="html",
        sha256="3cbc3a22daf3d494ffb26ce81a1bc951cd5dca6cb04c27dd2fd478f9a5dc9bf1",
        licence="page copyright FORCE Technology; digests + extracted cells only",
    ),
    SourceSpec(
        key="ittc_benchmarks",
        ship="KCS|KVLCC2|JBC|ONRT|Series60|DTC",
        publisher="ITTC Benchmark Repository — 'List of Benchmarks' (ITTC, Zürich)",
        family="ITTC",
        url="https://ittc.info/media/11250/list-of-benchmarks-2.pdf",
        kind="pdf",
        sha256="da101e9a352108e4cc73eb5f871fbe27d585f8e7bd4fb21af4261edd85fbaca1",
        licence="ITTC publication; bibliographic endorsement of the hull set (no numerals read)",
    ),
)


@dataclass
class Cell:
    """One quantity as stated by one source."""

    quantity: str
    value: float
    raw: str
    label: str
    source: str
    url: str
    sha256: str
    column: int
    unit_basis: str = "volume"
    note: str = ""

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


def split_ships(spec: SourceSpec, payload: bytes | str) -> list[tuple[str, str]]:
    """Split a multi-ship page into one markup block per ship.

    The workshop summary page prints one numbered table per hull ("Table 2 Main Particulars
    of KCS"); the block a ship owns runs from its own heading to the next heading.  A
    single-ship page is returned whole.
    """
    markup = payload.decode("utf-8", "replace") if isinstance(payload, bytes) else payload
    ships = spec.ship.split("|")
    if len(ships) == 1:
        return [(ships[0], markup)]
    hits = [
        (m.start(), m.group(1).upper())
        for m in re.finditer(r"Table\s*\d+\s*Main Particulars of\s+([A-Za-z0-9]+)", markup, flags=re.I)
    ]
    out: list[tuple[str, str]] = []
    for index, (start, name) in enumerate(hits):
        end = hits[index + 1][0] if index + 1 < len(hits) else len(markup)
        out.append((name, markup[start:end]))
    return out


def _split_section(markup: str, spec: SourceSpec) -> str:
    """Restrict a page to the block that belongs to ``spec.ship``."""
    if not spec.section:
        return markup
    start = markup.find(spec.section)
    if start < 0:
        return ""
    nxt = re.search(r"<table", markup[start + len(spec.section):], flags=re.I)
    end = start + len(spec.section) + (nxt.start() if nxt else len(markup) - start)
    return markup[start:end]


def _pick(cells: Sequence[str], column: int | None) -> tuple[float, str, int, str] | None:
    """The numeric cell that carries the full-scale value, with its index and basis."""
    order: Iterable[int]
    if column is not None:
        order = [column] + [i for i in range(1, len(cells)) if i != column]
    else:
        order = range(1, len(cells))
    for index in order:
        if index >= len(cells):
            continue
        text = cells[index]
        match = _NUMBER.search(text.replace(",", ""))
        if match is None:
            continue
        basis = "volume"
        if re.search(r"(?i)\btonns?\b|\bton\b|\bt\s*$|\(t\)", text):
            basis = "mass"
        return float(match.group(0)), text, index, basis
    return None


def parse_source(spec: SourceSpec, payload: bytes | str) -> dict[str, Cell]:
    """Read the tracked quantities of ``spec.ship`` out of a fetched document."""
    if spec.kind != "html":
        return {}
    markup = payload.decode("utf-8", "replace") if isinstance(payload, bytes) else payload
    block = _split_section(markup, spec)
    out: dict[str, Cell] = {}
    for row in table_rows(block):
        if not row:
            continue
        # The *label* lives in the first cell.  Reading it from the flattened row instead
        # lets derivative rows win ("izz/Lpp 0.25" is not Lpp), which is how a page can
        # silently hand back the wrong quantity.
        label = row[0]
        if re.search(r"(?i)izz|ixx|hub|\bic\b", label):
            continue
        quantity = None
        for name, pattern in LABEL_PATTERNS:
            if re.search(pattern, label, flags=re.I):
                quantity = name
                break
        if quantity is None or quantity in out:
            continue
        picked = _pick(row, spec.column)
        if picked is None:
            continue
        value, raw, index, basis = picked
        if spec.column is not None and index != spec.column:
            # The configured full-scale column was empty, so a neighbour was read.  Its scale
            # is unknown: keep it for the audit trail, never let it vote.
            basis = "fallback_column"
        out[quantity] = Cell(
            quantity=quantity,
            value=value,
            raw=raw,
            label=row[0][:120],
            source=spec.key,
            url=spec.url,
            sha256=spec.sha256,
            column=index,
            unit_basis=basis,
        )
    if "displacement_volume" in out:
        cell = out["displacement_volume"]
        if cell.unit_basis == "mass":
            cell.value = round(cell.value / SEAWATER_DENSITY, 1)
            cell.note = f"mass {cell.raw} converted at {SEAWATER_DENSITY} t/m3"
    return out


def spread(values: Sequence[float]) -> tuple[float, float]:
    """Absolute range and range relative to the mean (guard for a zero mean)."""
    lo, hi = min(values), max(values)
    mean = sum(values) / len(values)
    rel = abs(hi - lo) / abs(mean) if mean else float("inf")
    return hi - lo, rel


def cross_reference(extractions: Mapping[str, Mapping[str, Mapping[str, Cell]]]) -> dict[str, Any]:
    """Per ship, per quantity: the independent families that state it and whether they agree.

    ``extractions`` maps ``ship -> source_key -> quantity -> Cell``.
    """
    report: dict[str, Any] = {}
    for ship in sorted(extractions):
        per_quantity: dict[str, Any] = {}
        for quantity in QUANTITIES:
            entries = []
            for source_key, cells in sorted(extractions[ship].items()):
                cell = cells.get(quantity)
                if cell is None:
                    continue
                family = LINEAGE.get(source_key, FAMILIES.get(source_key, source_key))
                entries.append({"family": family, "cell": cell.as_dict()})
                if cell.unit_basis.startswith("fallback_column"):
                    continue
            families: dict[str, list[float]] = {}
            for entry in entries:
                # A cell read from a fallback column has an unknown scale: it stays in
                # ``entries`` for the audit trail but never votes.  The original loop put the
                # guard at the *end* of the collection loop, where it excluded nothing, so a
                # neighbouring-column reading could outvote (or silently corroborate) the real
                # full-scale value -- two tests had been failing on exactly that.
                if entry["cell"]["unit_basis"].startswith("fallback_column"):
                    continue
                families.setdefault(entry["family"], []).append(entry["cell"]["value"])
            family_values = {fam: sum(vals) / len(vals) for fam, vals in sorted(families.items())}
            if len(family_values) >= 2:
                delta, rel = spread(list(family_values.values()))
                tol = TOLERANCE[quantity]
                verdict = "CROSS_VERIFIED" if rel <= tol else "CONFLICT"
            else:
                delta, rel = 0.0, 0.0
                verdict = "SINGLE_SOURCE" if family_values else "MISSING"
            ordered = [
                {"family": fam, "value": value}
                for fam, value in sorted(family_values.items())
            ]
            per_quantity[quantity] = {
                "verdict": verdict,
                # A disagreement above the tolerance is *disputed*, not averaged: both
                # family values stay in ``ordered`` so a consumer must choose explicitly.
                "disputed": verdict == "CONFLICT",
                "published_values": ordered,
                "independent_families": len(family_values),
                "family_values": family_values,
                "range": delta,
                "range_relative": rel,
                "tolerance": TOLERANCE[quantity],
                "entries": entries,
            }
        # Same-family repeats are counted separately: they are corroboration, not independence.
        repeated = sum(
            1
            for quantity in QUANTITIES
            for src in [per_quantity[quantity]]
            if len(src["entries"]) > src["independent_families"] > 0
        )
        report[ship] = {
            "quantities": per_quantity,
            "cross_verified": sum(1 for q in per_quantity.values() if q["verdict"] == "CROSS_VERIFIED"),
            "single_source": sum(1 for q in per_quantity.values() if q["verdict"] == "SINGLE_SOURCE"),
            "conflict": sum(1 for q in per_quantity.values() if q["verdict"] == "CONFLICT"),
            "missing": sum(1 for q in per_quantity.values() if q["verdict"] == "MISSING"),
            "corroborated_by_repeat_family": repeated,
            "disputed_quantities": [
                q for q in QUANTITIES if per_quantity[q]["verdict"] == "CONFLICT"
            ],
        }
    return report


#: ``source key -> independence family``, derived from :data:`SOURCES`.
FAMILIES: dict[str, str] = {spec.key: spec.family for spec in SOURCES}

#: Independence *lineage*: the unit that may count as one independent source.  It is
#: stricter than :data:`FAMILIES`, because a later edition of the same workshop reproduces
#: the same upstream table: SIMMAN 2014 and SIMMAN 2020 state KVLCC2 as ``Bwl 580`` where the
#: ship is 58.0 m wide -- the identical defect in both editions is the evidence that one
#: table was copied from the other, so their agreement is *not* corroboration.  KCS is
#: backed by two lineages (the NMRI Tokyo 2015 page and the SIMMAN series), KVLCC2 by one.
LINEAGE: dict[str, str] = {
    "t2015_kcs": "NMRI-T2015",
    "t2015_jbc": "NMRI-T2015",
    "t2015_onrt": "NMRI-T2015",
    "simman2020_summary": "SIMMAN-SERIES",
    "simman2014_kvlcc2": "SIMMAN-SERIES",
    "simman2014_kcs": "SIMMAN-SERIES",
    "ittc_benchmarks": "ITTC",
}


def closure_check(report: Mapping[str, Any], ships: Sequence[str] | None = None) -> dict[str, Any]:
    """Does each source's own ∇ match its own ``Cb·Lpp·Bwl·T``?

    A publisher table is internally consistent or it is not; the check needs no external
    truth.  When it fails, the *implied* breadth is returned, which is what exposes the
    KVLCC2 ``Bwl 580`` defect (the implied breadth is 58.0 m, i.e. the printed 580 is a
    factor-of-ten typo, and it travels with the table into every later edition).
    """
    out: dict[str, Any] = {}
    for ship in ships or sorted(report):
        cell = report[ship]["quantities"]
        volume = cell["displacement_volume"]["family_values"]
        cb = cell["cb"]["family_values"]
        lpp = cell["length_pp"]["family_values"]
        beam = cell["beam"]["family_values"]
        draft = cell["draft"]["family_values"]
        rows: dict[str, Any] = {}
        for family in sorted(volume):
            if not all(family in d for d in (cb, lpp, beam, draft)):
                continue
            product = cb[family] * lpp[family] * beam[family] * draft[family]
            rel = (product - volume[family]) / volume[family] if volume[family] else float("inf")
            implied = volume[family] / (cb[family] * lpp[family] * draft[family]) if cb[family] * lpp[family] * draft[family] else float("nan")
            rows[family] = {
                "cb_lpp_bwl_t_m3": product,
                "published_volume_m3": volume[family],
                "rel_error": rel,
                "consistent": abs(rel) <= 0.01,
                "implied_beam_m": implied,
            }
        out[ship] = rows
    return out


def coverage_summary(report: Mapping[str, Any]) -> dict[str, Any]:
    """Fleet-level tallies used by the gate and by the report."""
    ships = sorted(report)
    cells = [report[s]["quantities"][q] for s in ships for q in QUANTITIES]
    return {
        "ships": ships,
        "cells": len(cells),
        "cross_verified": sum(1 for c in cells if c["verdict"] == "CROSS_VERIFIED"),
        "single_source": sum(1 for c in cells if c["verdict"] == "SINGLE_SOURCE"),
        "conflict": sum(1 for c in cells if c["verdict"] == "CONFLICT"),
        "missing": sum(1 for c in cells if c["verdict"] == "MISSING"),
        "ships_fully_cross_verified": [
            s
            for s in ships
            if all(report[s]["quantities"][q]["verdict"] == "CROSS_VERIFIED" for q in QUANTITIES)
        ],
        "ships_failing_min_two_sources": [
            s
            for s in ships
            if any(
                report[s]["quantities"][q]["independent_families"] < 2
                for q in ("displacement_volume", "cb", "lcb_percent", "lwl")
            )
        ],
    }

#: A pinned digest is a SHA-256 hex string, nothing else.  ``pin_audit`` exists because a
#: digest that was typed rather than fetched still *looks* like a fingerprint: two of the
#: seven source digests in the first version of this registry shared their leading eight
#: hex digits with the real bytes and then diverged, which is what a hand-written
#: continuation looks like.  The audit recomputes every digest from the bytes and says so.
DIGEST_PATTERN = re.compile(r"^[0-9a-f]{64}$")


def pin_audit(fetched: Mapping[str, bytes | str]) -> dict[str, dict[str, Any]]:
    """Per source: is the pinned digest the digest of the bytes that were actually read?

    ``fetched`` maps source key -> payload.  A key that is not in ``fetched`` is reported as
    ``not_fetched``; it never counts as verified.  ``hand_written`` marks a pinned digest
    that is not even 64 lowercase hex digits, and ``drift`` marks a fetch whose digest no
    longer matches the pin -- a source change that must be re-read deliberately.
    """
    import hashlib

    out: dict[str, dict[str, Any]] = {}
    for spec in SOURCES:
        payload = fetched.get(spec.key)
        row: dict[str, Any] = {"url": spec.url, "family": spec.family}
        if payload is None:
            row["verdict"] = "not_fetched"
            row["pinned"] = spec.sha256
            out[spec.key] = row
            continue
        data = payload.encode("utf-8") if isinstance(payload, str) else payload
        digest = hashlib.sha256(data).hexdigest()
        row.update(pinned=spec.sha256, fetched=digest, bytes=len(data))
        if not spec.sha256:
            row["verdict"] = "unpinned"
        elif not DIGEST_PATTERN.match(spec.sha256):
            row["verdict"] = "hand_written"
        elif spec.sha256 == digest:
            row["verdict"] = "match"
        else:
            row["verdict"] = "drift"
        out[spec.key] = row
    return out


@dataclass(frozen=True)
class Reference:
    """One published quantity, resolved to a single value *and* to its provenance."""

    ship: str
    quantity: str
    value: float
    disputed: bool
    families: tuple[str, ...]
    #: Every family value, so a disputed reference can be quoted both ways.
    values: tuple[tuple[str, float], ...]
    #: ``(source_key, sha256)`` for every entry that backs the value.
    provenance: tuple[tuple[str, str], ...]

    def as_dict(self) -> dict[str, Any]:
        return {
            "ship": self.ship,
            "quantity": self.quantity,
            "value": self.value,
            "disputed": self.disputed,
            "families": list(self.families),
            "values": [{"family": fam, "value": value} for fam, value in self.values],
            "provenance": [{"source": key, "sha256": sha} for key, sha in self.provenance],
        }


def resolve_reference(
    report: Mapping[str, Any],
    ship: str,
    quantity: str = "displacement_volume",
    *,
    preferred_family: str | None = None,
) -> Reference | None:
    """Resolve one published quantity to a value plus provenance, never to an average.

    * agreed sources -> the agreed value;
    * one family -> that family's value, ``disputed=False`` (honest single source);
    * several families that disagree -> ``disputed=True`` with *all* values kept.  The
      returned ``value`` is the preferred family's, else the family with the most entries
      (ties broken alphabetically) -- a deterministic rule, not a mean, and the caller can
      see from ``disputed`` that the number is one of two.
    * nothing published -> ``None`` (the caller must fall back to a closure, visibly).
    """
    quantity_row = report[ship]["quantities"][quantity]
    families = quantity_row["family_values"]
    if not families:
        return None
    entries = quantity_row["entries"]
    counts: dict[str, int] = {}
    provenance: list[tuple[str, str]] = []
    for entry in entries:
        fam = entry["family"]
        counts[fam] = counts.get(fam, 0) + 1
        provenance.append((entry["cell"]["source"], entry["cell"]["sha256"]))
    if preferred_family is not None:
        if preferred_family not in families:
            raise KeyError(f"{preferred_family!r} states no {quantity} for {ship}")
        chosen = preferred_family
    elif len(families) == 1:
        chosen = next(iter(families))
    else:
        chosen = sorted(families, key=lambda fam: (-counts[fam], fam))[0]
    return Reference(
        ship=ship,
        quantity=quantity,
        value=families[chosen],
        disputed=quantity_row["verdict"] == "CONFLICT",
        families=tuple(sorted(families)),
        values=tuple(sorted(families.items())),
        provenance=tuple(sorted(set(provenance))),
    )
