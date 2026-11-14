"""Pointer hygiene for :mod:`tensorcad.fleet_baseline` spec files.

A spec file is the contract between the geometry bytes and the benchmark. Its pointers
used to be *paths into a mission work tree* (``/root/tp-agent/work/m18488/...``) and
*paths relative to that work tree* (``artifacts/baseline_fleet/raw/pages/*.html``).
Both die with the mission, and they die silently: the benchmark cannot be re-run, but
nothing in the repository says so.

Three states are distinguished, and only three:

``resolved``
    the pointer resolves to bytes whose SHA-256 equals the SHA-256 the spec records;

``dead``
    the pointer resolves to nothing, and the recorded SHA-256 still identifies the
    byte — so the pointer can be repaired *by content*, not by guessing a path;

``withheld``
    the pointer resolves to nothing and no durable byte exists anywhere (the publisher
    page cache). Repair keeps the loss visible: the path becomes the publisher URL
    (a stable citation the spec did not carry before) and the entry is marked
    ``pointer_status: "cache-lost"`` with the original path kept as ``legacy_path``.

A pointer that is already an ``http(s)`` URL is ``cited``: it names a durable source
rather than a local byte, so there is nothing left to repair and nothing to verify.

Repair never invents a path: geometry is re-pointed only when
``iges_sha256`` finds the byte in the store, and the rewritten pointer is an
``asset://`` URI (machine-independent) rather than an absolute path.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping

from .registry import AssetRegistry
from .store import AssetError, AssetStore, asset_id_of, asset_uri, sha256_file

GEOMETRY = "geometry"
PAGE = "page"
ARCHIVE = "archive"
PUBLISHED_SOURCE = "published_source"


@dataclass(frozen=True)
class Pointer:
    """One pointer inside a spec file, its state, and what it should become."""

    where: str
    kind: str
    value: str
    status: str
    sha256: str = ""
    asset_id: str = ""
    replacement: str = ""
    detail: str = ""

    @property
    def ok(self) -> bool:
        return self.status == "resolved"

    @property
    def settled(self) -> bool:
        """True when the pointer needs no repair (resolves, or cites a durable URL)."""
        return self.status in {"resolved", "cited"}

    def as_dict(self) -> dict[str, Any]:
        return {
            "where": self.where,
            "kind": self.kind,
            "value": self.value,
            "status": self.status,
            "sha256": self.sha256,
            "asset_id": self.asset_id,
            "replacement": self.replacement,
            "detail": self.detail,
        }


@dataclass
class PointerReport:
    path: str
    pointers: list[Pointer] = field(default_factory=list)

    def by_status(self, status: str) -> list[Pointer]:
        return [p for p in self.pointers if p.status == status]

    @property
    def dead(self) -> list[Pointer]:
        return self.by_status("dead")

    @property
    def withheld(self) -> list[Pointer]:
        return self.by_status("withheld")

    @property
    def resolved(self) -> list[Pointer]:
        return [p for p in self.pointers if p.ok]

    @property
    def cited(self) -> list[Pointer]:
        return self.by_status("cited")

    def fraction(self, predicate) -> float:
        if not self.pointers:
            return 0.0
        return round(sum(1 for p in self.pointers if predicate(p)) / len(self.pointers), 6)

    def as_dict(self) -> dict[str, Any]:
        return {
            "path": self.path,
            "pointers": len(self.pointers),
            "resolved": len(self.resolved),
            "cited": len(self.cited),
            "dead": len(self.dead),
            "withheld": len(self.withheld),
            "resolvable_fraction": self.fraction(lambda p: p.ok),
            "detail": [p.as_dict() for p in self.pointers],
        }


def load_spec_payload(path: str | Path) -> dict[str, Any]:
    return json.loads(Path(path).read_text())


def iter_pointers(payload: Mapping[str, Any]) -> Iterable[tuple[str, str, str, str]]:
    """Yield ``(where, kind, value, sha256)`` for every pointer the spec carries."""
    ships = payload.get("ships", []) if isinstance(payload, Mapping) else payload
    for index, ship in enumerate(ships or []):
        name = str(ship.get("name", f"#{index}"))
        if ship.get("iges_path"):
            yield (f"ships[{index}]({name}).iges_path", GEOMETRY, str(ship["iges_path"]), str(ship.get("iges_sha256", "")))
        if ship.get("published_source"):
            yield (
                f"ships[{index}]({name}).published_source",
                PUBLISHED_SOURCE,
                str(ship["published_source"]),
                str((payload.get("evidence", {}).get(name, {}) or {}).get("page", {}).get("sha256", "")),
            )
    evidence = payload.get("evidence", {}) if isinstance(payload, Mapping) else {}
    for name, blob in (evidence or {}).items():
        block = blob or {}
        iges = block.get("iges") or {}
        if iges.get("path"):
            yield (f"evidence.{name}.iges.path", GEOMETRY, str(iges["path"]), str(iges.get("sha256", "")))
        archive = block.get("archive") or {}
        if archive.get("path"):
            yield (f"evidence.{name}.archive.path", ARCHIVE, str(archive["path"]), str(archive.get("sha256", "")))
        page = block.get("page") or {}
        if page.get("path"):
            yield (f"evidence.{name}.page.path", PAGE, str(page["path"]), str(page.get("sha256", "")))
        for field_name, value in (block.get("published_fields") or {}).items():
            target = (value or {}).get("page")
            if target:
                yield (
                    f"evidence.{name}.published_fields.{field_name}.page",
                    PAGE,
                    str(target),
                    str((value or {}).get("page_sha256", "")),
                )


def _candidate_roots(spec_path: Path, repo_root: Path | None) -> list[Path]:
    roots = [Path(spec_path).resolve().parent]
    if repo_root is not None:
        roots.append(Path(repo_root).resolve())
    roots.append(Path.cwd().resolve())
    return roots


def _locate(value: str, roots: list[Path], store: AssetStore) -> Path | None:
    asset_id = asset_id_of(value)
    if asset_id:
        path = store.path_of(asset_id)
        return path if path.is_file() else None
    raw = Path(value)
    if raw.is_absolute():
        return raw if raw.is_file() else None
    for root in roots:
        candidate = (root / value).resolve()
        if candidate.is_file():
            return candidate
    return None


def scan_spec_pointers(
    spec_path: str | Path,
    *,
    store: AssetStore | None = None,
    repo_root: str | Path | None = None,
    registry: AssetRegistry | None = None,
) -> PointerReport:
    """Classify every pointer in ``spec_path`` without modifying it."""
    store = store or AssetStore((registry.store_root if registry else "") or None)
    payload = load_spec_payload(spec_path)
    roots = _candidate_roots(Path(spec_path), None if repo_root is None else Path(repo_root))
    report = PointerReport(path=str(spec_path))
    for where, kind, value, sha in iter_pointers(payload):
        if str(value).startswith(("http://", "https://")):
            report.pointers.append(
                Pointer(
                    where=where,
                    kind=kind,
                    value=value,
                    status="cited",
                    sha256=sha,
                    detail="durable publisher URL (no local byte to verify)",
                )
            )
            continue
        located = _locate(value, roots, store)
        asset_id = ""
        if located is not None and sha:
            actual = sha256_file(located)
            if actual != sha:
                report.pointers.append(
                    Pointer(
                        where=where,
                        kind=kind,
                        value=value,
                        status="dead",
                        sha256=sha,
                        detail=f"resolves to {located} but sha256 {actual} != recorded {sha}",
                    )
                )
                continue
        if located is not None:
            if asset_id_of(value):
                asset_id = asset_id_of(value) or ""
            report.pointers.append(
                Pointer(
                    where=where,
                    kind=kind,
                    value=value,
                    status="resolved",
                    sha256=sha,
                    asset_id=asset_id,
                    detail=str(located),
                )
            )
            continue
        replacement = ""
        found = ""
        if kind == GEOMETRY and sha:
            if registry is not None:
                found = registry.find_by_sha(sha, store=store) or ""
            else:
                found = store.find_by_sha(sha) or ""
            if found:
                replacement = asset_uri(found)
        status = "dead" if replacement else "withheld"
        detail = (
            f"absent; byte recovered by sha256 as {replacement}"
            if replacement
            else "absent; no durable byte and no SHA-256 to recover it by"
        )
        report.pointers.append(
            Pointer(
                where=where,
                kind=kind,
                value=value,
                status=status,
                sha256=sha,
                asset_id=found,
                replacement=replacement,
                detail=detail,
            )
        )
    return report


def _hull_of_where(where: str) -> str:
    """Ship name behind a ``where`` label (``ships[0](KCS).x`` / ``evidence.KCS.x``)."""
    if where.startswith("ships[") and "(" in where and ")" in where:
        return where.split("(", 1)[1].split(")", 1)[0]
    if where.startswith("evidence."):
        return where.split(".")[1]
    return ""


def _set_pointer(payload: dict[str, Any], where: str, replacement: str, *, status: str, legacy: str | None) -> None:
    """Rewrite one pointer addressed by the ``where`` strings :func:`iter_pointers` emits."""
    if where.startswith("ships["):
        index = int(where[len("ships[") : where.index("]")])
        ship = payload["ships"][index]
        key = where.split(").", 1)[1]
        if status == "cache-lost":
            ship["pointer_status"] = {**(ship.get("pointer_status") or {}), key: "cache-lost"}
            ship[key + "_legacy"] = legacy
        ship[key] = replacement
        return
    parts = where.split(".")
    if parts[0] != "evidence":
        raise KeyError(f"unknown pointer {where!r}")
    hull, tail = parts[1], parts[2:]
    node: Any = payload["evidence"][hull]
    for key in tail[:-1]:
        node = node[key]
    last = tail[-1]
    if status == "cache-lost":
        container = node if isinstance(node, dict) else None
        if container is not None:
            statuses = container.setdefault("pointer_status", {})
            if isinstance(statuses, dict):
                statuses[last] = "cache-lost"
            container["legacy_path"] = legacy
    node[last] = replacement


def repair_spec_file(
    spec_path: str | Path,
    *,
    out_path: str | Path | None = None,
    store: AssetStore | None = None,
    registry: AssetRegistry | None = None,
    repo_root: str | Path | None = None,
    url_map: Mapping[str, str] | None = None,
    archive_url_map: Mapping[str, str] | None = None,
    dry_run: bool = False,
) -> dict[str, Any]:
    """Re-point dead pointers and return a change log.

    * geometry with a recorded SHA-256 present in the store -> ``asset://`` URI;
    * page/archive pointers with a known publisher URL -> the URL, marked
      ``cache-lost`` with the original path kept in ``legacy_path``;
    * anything else is *left alone* and reported as ``withheld``.
    """
    spec_path = Path(spec_path)
    payload = load_spec_payload(spec_path)
    report = scan_spec_pointers(spec_path, store=store, repo_root=repo_root, registry=registry)
    url_map = dict(url_map or {})
    archive_url_map = dict(archive_url_map or {})
    evidence = payload.get("evidence", {}) or {}
    changes: list[dict[str, Any]] = []
    unresolved: list[dict[str, Any]] = []
    for pointer in report.pointers:
        if pointer.settled:
            continue
        replacement = pointer.replacement
        status = "asset-uri"
        if not replacement and pointer.kind in {PAGE, PUBLISHED_SOURCE}:
            hull = _hull_of_where(pointer.where)
            page_key = str(((evidence.get(hull, {}) or {}).get("page", {}) or {}).get("key", ""))
            replacement = url_map.get(page_key, "")
            status = "cache-lost"
        if not replacement and pointer.kind == ARCHIVE:
            basename = Path(pointer.value).name
            replacement = archive_url_map.get(basename, "")
            status = "cache-lost"
        if not replacement:
            unresolved.append(pointer.as_dict())
            continue
        _set_pointer(payload, pointer.where, replacement, status=status, legacy=pointer.value)
        changes.append(
            {
                "where": pointer.where,
                "kind": pointer.kind,
                "from": pointer.value,
                "to": replacement,
                "status": status,
                "sha256": pointer.sha256,
            }
        )
    result = {
        "schema": "tensorcad.assetlib.pointer_repair.v1",
        "spec": str(spec_path),
        "changes": changes,
        "unresolved": unresolved,
        "before": report.as_dict(),
    }
    if not dry_run:
        target = Path(out_path) if out_path is not None else spec_path
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload, indent=2) + "\n")
        after = scan_spec_pointers(target, store=store, repo_root=repo_root, registry=registry)
        result["after"] = after.as_dict()
        result["output"] = str(target)
    return result