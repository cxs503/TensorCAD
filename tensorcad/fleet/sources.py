"""Per-ship source anchors: the citation behind every manifest digit, and its check state.

The fleet manifest carries numbers; this registry carries *where they came from*, in a
form a reader can re-derive them from: citation, URL, licence, and the locator (the
table or section the digits would be read from). The two files are joined by ship name
and the join is gated: a ship with no anchor cannot pass ``G0_provenance``, because a
number with no named source is precisely the defect the fleet anchor exists to remove.

The registry is also honest about the difference between "published" and "checked":
every entry records ``retrieval`` (was the source re-read in this environment?) and
``verified`` (was it confirmed against the source?), and both are false here because
this environment has no egress. The arithmetic check that *is* possible offline lives
in :func:`tensorcad.fleet.gates.gate_reference_consistency`: a cited displacement has
to agree with the record's own ``Cb*L*B*T``.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping

REGISTRY_PATH = Path(__file__).with_name("data") / "source_registry.json"
SCHEMA = "tensorcad.fleet.source_registry.v1"


class RegistryError(ValueError):
    """The source registry is missing, malformed, or does not cover the manifest."""


@dataclass(frozen=True)
class SourceAnchor:
    """Where one ship's published particulars come from, and how far that was checked."""

    ship: str
    anchor_id: str
    citation: str
    url: str
    license: str
    locator: str
    anchors: tuple[str, ...]
    retrieval: str
    verified: bool
    note: str = ""

    @property
    def complete(self) -> bool:
        return all(
            str(value).strip() for value in (self.citation, self.url, self.license, self.locator, self.anchor_id)
        )

    def as_dict(self) -> dict[str, Any]:
        return {
            "ship": self.ship,
            "anchor_id": self.anchor_id,
            "citation": self.citation,
            "url": self.url,
            "license": self.license,
            "locator": self.locator,
            "anchors": list(self.anchors),
            "retrieval": self.retrieval,
            "verified": self.verified,
            "note": self.note,
        }


@dataclass(frozen=True)
class SourceRegistry:
    sources: tuple[SourceAnchor, ...]
    sha256: str | None = None
    path: Path | None = None

    def for_ship(self, ship: str) -> SourceAnchor | None:
        for source in self.sources:
            if source.ship == ship:
                return source
        return None

    def missing_for(self, ships: Iterable[str]) -> list[str]:
        return [name for name in ships if self.for_ship(name) is None]


def load_registry(path: str | Path | None = None) -> SourceRegistry:
    import hashlib

    target = Path(path) if path is not None else REGISTRY_PATH
    payload = target.read_bytes()
    try:
        document = json.loads(payload)
    except json.JSONDecodeError as exc:  # pragma: no cover - malformed file
        raise RegistryError(f"{target}: not valid JSON ({exc})") from exc
    if document.get("schema") != SCHEMA:
        raise RegistryError(f"{target}: schema {document.get('schema')!r} is not {SCHEMA!r}")
    anchors: list[SourceAnchor] = []
    for entry in document.get("sources", []):
        anchor = SourceAnchor(
            ship=str(entry["ship"]),
            anchor_id=str(entry["anchor_id"]),
            citation=str(entry["citation"]),
            url=str(entry["url"]),
            license=str(entry["license"]),
            locator=str(entry["locator"]),
            anchors=tuple(entry.get("anchors", ())),
            retrieval=str(entry.get("retrieval", "unknown")),
            verified=bool(entry.get("verified", False)),
            note=str(entry.get("note", "")),
        )
        if not anchor.complete:
            raise RegistryError(f"{target}: anchor for {anchor.ship!r} is incomplete")
        anchors.append(anchor)
    return SourceRegistry(
        sources=tuple(anchors),
        sha256=hashlib.sha256(payload).hexdigest(),
        path=target,
    )