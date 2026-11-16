"""Re-bind every published digit of the fleet reference library to the bytes it was read from.

The library (``benchmarks/verified/fleet_refs/crossref.json``) claims, for every published
particular, a *source*, a *URL* and the SHA-256 of the page or PDF the digit was read from.
Until this module existed that claim rested on two half-checks that never met in the middle:

* ``pin_audit`` recomputed digests of the round-one pages -- but it says nothing about the
  numerals the library carries; a digest can match while the cell was typed by hand.
* ``oa_ingest`` pinned the digests of the open-access PDFs at ingest time -- and nothing ever
  re-read those bytes, so a cell whose arithmetic changed under a code revision stayed in the
  published library unnoticed.

:func:`audit_cells` closes the chain end to end: bytes -> digest -> text -> numeral -> library
value, document by document.  A cell that does not come back out of its own document is
reported, never quietly accepted; and the count of cells that *did* come back is what the
anchor's acceptance clause ("every value carries a citation plus a cached source-page digest")
may honestly claim.

Two classes of stored cell are distinguished, because they deserve different treatment:

``transcribed``
    the numeral is printed in the document.  A disagreement means the library misquotes its
    own witness.  It may never be repaired automatically -- the repair would be a new reading
    dressed up as an old one -- so it is reported as ``TRANSCRIPTION_DRIFT`` for a human to
    re-read.
``derived`` (``unit_basis`` starts with ``derived``)
    the numeral is arithmetic on two printed numerals of the same document (e.g. a model-scale
    distance converted with the table's own length ratio).  A disagreement means the code's
    arithmetic changed, not that evidence moved, so it is reported as ``STALE_DERIVED`` and is
    safe to recompute from the same bytes.

Digest stability is treated with the same care.  A page that embeds a request counter has no
stable byte digest, and pinning one guarantees a red freshness check forever -- which trains a
reader to ignore drift.  :func:`stable_digest` masks counted nonces that were observed in the
pinned pages and reports which masks it applied, so "the bytes changed" and "the content
changed" stop being the same sentence.
"""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass, field
from typing import Any, Iterable, Mapping, Sequence

SCHEMA = "tensorcad.fleet_refs.rebind/1.0"

#: Byte patterns that change on every request and carry no document content.  Each mask is
#: anchored on the publisher's own widget so it can never eat a numeral of a table: the
#: SIMMAN-2014 pages emit ``ML:768-PU:<request counter>-PLA:19`` on two consecutive fetches
#: that are otherwise byte-identical (see ``docs/fleet_refs_rebind.md``).
VOLATILE_PATTERNS: tuple[tuple[str, str], ...] = (
    ("simman_request_counter", r"ML:\d+-PU:\d+-PLA:\d+"),
)

#: Cell verdicts.
MATCH = "MATCH"
STALE_DERIVED = "STALE_DERIVED"
TRANSCRIPTION_DRIFT = "TRANSCRIPTION_DRIFT"
NOT_IN_BYTES = "NOT_IN_BYTES"
UNCHECKED = "UNCHECKED"

#: Document verdicts.
DIGEST_MATCH = "DIGEST_MATCH"
DIGEST_VOLATILE_ONLY = "DIGEST_VOLATILE_ONLY"
DIGEST_DRIFT = "DIGEST_DRIFT"
NOT_CACHED = "NOT_CACHED"
NO_PINNED_DIGEST = "NO_PINNED_DIGEST"


def _mask(payload: bytes) -> tuple[bytes, list[str]]:
    """Replace counted nonces; return the masked bytes and the names of the masks applied."""
    applied: list[str] = []
    masked = payload
    for name, pattern in VOLATILE_PATTERNS:
        new = re.sub(pattern.encode("ascii"), b"<nonce>", masked)
        if new != masked:
            applied.append(name)
        masked = new
    return masked, applied


def stable_digest(payload: bytes | str) -> tuple[str, list[str]]:
    """SHA-256 of ``payload`` with request counters masked, plus the masks that fired."""
    data = payload.encode("utf-8") if isinstance(payload, str) else payload
    masked, applied = _mask(data)
    return hashlib.sha256(masked).hexdigest(), applied


def digest_verdict(
    pinned: str | None,
    cached: bytes | str | None,
    live: bytes | str | None = None,
) -> dict[str, Any]:
    """Do the bytes still answer to the pinned digest, and did the *content* move?

    ``cached`` is the copy the library was read from; ``live`` is a fresh fetch, when the
    caller has one.  The two are kept apart because a page can change bytes on every request
    without changing a single digit: comparing *masked* digests (see :func:`stable_digest`)
    separates "the bytes moved" from "the document changed", and only the second is drift.

    ``NOT_CACHED`` and ``NO_PINNED_DIGEST`` are not passes: an unread document has not been
    checked, and saying otherwise is how an unverified library looks verified.
    """
    data = None
    if cached is not None:
        data = cached.encode("utf-8") if isinstance(cached, str) else cached
    live_data = None
    if live is not None:
        live_data = live.encode("utf-8") if isinstance(live, str) else live
    row: dict[str, Any] = {"pinned": pinned, "status": NOT_CACHED}
    if data is None and live_data is None:
        return row
    if data is None:
        row["status"] = NOT_CACHED
        row["live_only"] = True
    actual = hashlib.sha256(data).hexdigest() if data is not None else None
    row.update(
        cached=actual,
        cached_stable=stable_digest(data)[0] if data is not None else None,
        cached_bytes=len(data) if data is not None else None,
        masks=stable_digest(data)[1] if data is not None else [],
    )
    if live_data is not None:
        row.update(
            live=hashlib.sha256(live_data).hexdigest(),
            live_stable=stable_digest(live_data)[0],
            live_bytes=len(live_data),
            live_masks=stable_digest(live_data)[1],
        )
    if not pinned:
        row["status"] = NO_PINNED_DIGEST
        return row
    if actual is not None and actual == pinned:
        row["status"] = DIGEST_MATCH
    elif actual is not None:
        row["status"] = DIGEST_DRIFT
    if live_data is not None:
        if row.get("live") == pinned:
            row["live_status"] = DIGEST_MATCH
        elif data is not None and row.get("live_stable") == row.get("cached_stable"):
            row["live_status"] = DIGEST_VOLATILE_ONLY
            if row["status"] == DIGEST_DRIFT:
                row["status"] = DIGEST_VOLATILE_ONLY
        else:
            row["live_status"] = DIGEST_DRIFT
    return row


@dataclass(frozen=True)
class CellVerdict:
    """One library digit tested against the document it cites."""

    source: str
    ship: str
    quantity: str
    stored: float | None
    fresh: float | None
    status: str
    basis: str = ""
    detail: str = ""
    fresh_cell: Mapping[str, Any] | None = field(default=None, compare=False, repr=False)

    @property
    def ok(self) -> bool:
        return self.status == MATCH

    def as_dict(self) -> dict[str, Any]:
        out: dict[str, Any] = {
            "source": self.source,
            "ship": self.ship,
            "quantity": self.quantity,
            "stored": self.stored,
            "fresh": self.fresh,
            "status": self.status,
            "basis": self.basis,
        }
        if self.detail:
            out["detail"] = self.detail
        return out


def _stored_value(entry: Any) -> tuple[float | None, str]:
    """``(value, unit_basis)`` out of a library cell, tolerant of a bare number."""
    if isinstance(entry, Mapping):
        value = entry.get("value")
        basis = str(entry.get("unit_basis") or "")
        return (float(value) if isinstance(value, (int, float)) else None), basis
    if isinstance(entry, (int, float)):
        return float(entry), ""
    return None, ""


def _is_derived(basis: str) -> bool:
    return basis.startswith("derived")


def compare_cells(
    source: str,
    ship: str,
    stored: Mapping[str, Any],
    fresh: Mapping[str, Any],
    tol: float = 1e-9,
) -> list[CellVerdict]:
    """Every stored cell re-read out of the document, plus every numeral the document adds.

    ``stored`` maps quantity -> library cell (mapping or bare number); ``fresh`` maps
    quantity -> a :class:`tensorcad.fleet.crossref.Cell` (or its ``as_dict()`` form).  A
    quantity the document states but the library never carried is reported too -- that is how
    a *shrunken* library (a document read while a parser was weaker) becomes visible.
    """
    out: list[CellVerdict] = []
    for quantity in sorted(set(stored) | set(fresh)):
        stored_value, basis = _stored_value(stored.get(quantity))
        fresh_cell = fresh.get(quantity)
        fresh_value = None
        if fresh_cell is not None:
            fresh_value = (
                fresh_cell.get("value")
                if isinstance(fresh_cell, Mapping)
                else getattr(fresh_cell, "value", None)
            )
            fresh_value = float(fresh_value) if isinstance(fresh_value, (int, float)) else None
        if quantity not in stored:
            out.append(
                CellVerdict(
                    source=source,
                    ship=ship,
                    quantity=quantity,
                    stored=None,
                    fresh=fresh_value,
                    status=NOT_IN_BYTES,
                    basis="",
                    detail="the document states this quantity; the library does not carry it",
                )
            )
            continue
        if fresh_cell is None:
            out.append(
                CellVerdict(
                    source=source,
                    ship=ship,
                    quantity=quantity,
                    stored=stored_value,
                    fresh=None,
                    status=NOT_IN_BYTES,
                    basis=basis,
                    detail="the library carries this quantity; the document no longer yields it",
                )
            )
            continue
        if (
            stored_value is not None
            and fresh_value is not None
            and abs(stored_value - fresh_value) <= tol
        ):
            out.append(
                CellVerdict(
                    source=source,
                    ship=ship,
                    quantity=quantity,
                    stored=stored_value,
                    fresh=fresh_value,
                    status=MATCH,
                    basis=basis,
                    fresh_cell=fresh_cell if isinstance(fresh_cell, Mapping) else fresh_cell.as_dict(),
                )
            )
            continue
        derived = _is_derived(basis)
        out.append(
            CellVerdict(
                source=source,
                ship=ship,
                quantity=quantity,
                stored=stored_value,
                fresh=fresh_value,
                status=STALE_DERIVED if derived else TRANSCRIPTION_DRIFT,
                basis=basis,
                detail=(
                    f"library {stored_value} vs document {fresh_value}"
                    + (" (derived arithmetic: the code changed, not the evidence)" if derived else "")
                ),
                fresh_cell=fresh_cell if isinstance(fresh_cell, Mapping) else fresh_cell.as_dict(),
            )
        )
    return out


def repairable(verdicts: Iterable[CellVerdict]) -> list[CellVerdict]:
    """Only arithmetic may be recomputed; a transcription is never edited by a tool."""
    return [v for v in verdicts if v.status == STALE_DERIVED and v.fresh_cell is not None]


def apply_repair(
    payload: dict[str, Any],
    verdicts: Sequence[CellVerdict],
    *,
    ship_blocks: Mapping[str, str] = {"report": "quantities"},
) -> list[dict[str, Any]]:
    """Rewrite stale *derived* cells in a library payload from the re-read document.

    Transcribed cells and unmatched quantities are left untouched on purpose: a library that
    silently rewrites its own quotations is worse than a stale one.
    """
    done: list[dict[str, Any]] = []
    for verdict in repairable(verdicts):
        block_path = ship_blocks["report"]
        ship = payload.get("report", {}).get(verdict.ship, {}).get(block_path, {})
        cell_block = ship.get(verdict.quantity)
        if not isinstance(cell_block, dict):
            continue
        rows = cell_block.get("entries") or []
        changed = False
        for row in rows:
            cell = row.get("cell") or {}
            if cell.get("source") != verdict.source:
                continue
            if not _is_derived(str(cell.get("unit_basis") or "")):
                continue
            if abs(float(cell.get("value") or 0.0) - float(verdict.stored or 0.0)) > 1e-12:
                continue  # a different reading of the same source: not ours to touch
            before = cell.get("value")
            fresh = dict(verdict.fresh_cell or {})
            cell.update({k: fresh[k] for k in ("value", "raw", "note", "derivation") if k in fresh})
            done.append(
                {
                    "ship": verdict.ship,
                    "quantity": verdict.quantity,
                    "source": verdict.source,
                    "before": before,
                    "after": cell.get("value"),
                    "basis": cell.get("unit_basis"),
                }
            )
            changed = True
        if changed and "derived" in cell_block:
            cell_block["derived"] = {
                "cells": [
                    {
                        "family": row["cell"].get("family", row.get("family", "")),
                        "value": row["cell"]["value"],
                        "unit_basis": row["cell"].get("unit_basis", ""),
                        "derivation": row["cell"].get("derivation", ""),
                        "note": row["cell"].get("note", ""),
                    }
                    for row in rows
                    if row.get("cell", {}).get("unit_basis", "").startswith("derived")
                ]
            }
            # The derived-vs-voter arithmetic is cached beside the cell; drop the stale copy so
            # a later reader recomputes it instead of trusting a number this repair invalidated.
            cell_block["derived"].pop("vs_voters", None)
    return done


def summarize(documents: Sequence[Mapping[str, Any]], cells: Sequence[CellVerdict]) -> dict[str, Any]:
    """Counts and an explicit all-or-nothing verdict; unread documents are never a pass."""
    by_status: dict[str, int] = {}
    for verdict in cells:
        by_status[verdict.status] = by_status.get(verdict.status, 0) + 1
    docs_by_status: dict[str, int] = {}
    for row in documents:
        status = str(row.get("status"))
        docs_by_status[status] = docs_by_status.get(status, 0) + 1
    checked_sources = {
        row.get("key") for row in documents if row.get("status") != NOT_CACHED
    }
    unread = sorted(
        {verdict.source for verdict in cells if verdict.status == UNCHECKED}
        | {str(row.get("key")) for row in documents if row.get("status") == NOT_CACHED}
    )
    blocked = [verdict.as_dict() for verdict in cells if verdict.status != MATCH]
    return {
        "documents": len(documents),
        "documents_checked": len(checked_sources),
        "documents_by_status": docs_by_status,
        "cells": len(cells),
        "cells_by_status": by_status,
        "cells_rebound": by_status.get(MATCH, 0),
        "non_rebound": blocked,
        "unchecked_sources": unread,
        "all_cells_rebound": not blocked and not unread,
    }