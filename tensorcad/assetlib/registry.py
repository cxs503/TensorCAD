"""The geometry asset registry: a git-versioned index over the durable store.

Why a registry at all
---------------------
Pointers rot, and a rotted pointer looks exactly like a working one until somebody
opens it. A registry replaces *paths* with *claims that can be checked*: for every
asset it records the SHA-256 and the byte count, so "the file is gone" and "the file
came back different" become different, detectable states rather than a traceback in a
mission report. It also joins three things that used to live apart:

* the store's per-hull ``MANIFEST.json`` (what the importer claimed),
* ``benchmarks/verified/baseline_fleet/fleet_specs.json`` (what the benchmark cites),
* the benchmark verdict (what was actually proven with those bytes).

The registry carries **no bytes**. The workshop archives are redistributed by nobody
(see ``docs/asset_library.md``); what the repo gets is the fingerprint, which is
enough to detect loss and to re-import from a legitimate source.

Growth metrics
--------------
:meth:`AssetRegistry.growth` differentiates two registry snapshots. The metrics are
deliberately boring and hard to game:

``assets`` / ``bytes``
    how many verifiable files the library holds, and their size;
``verified_fraction``
    share of assets whose bytes matched the recorded SHA-256 *at build time*;
``spec_bindings`` / ``coverage``
    how many fleet-spec ships resolve to a verified asset — the number that says
    whether downstream physics can still be re-run from durable bytes;
``resolvable_pointer_fraction``
    share of geometry pointers in the specs that resolve (no dead pointers).
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping

from .store import AssetError, AssetStore, asset_id_of, asset_uri, hull_of, sha256_file

REGISTRY_SCHEMA = "tensorcad.assetlib.registry.v1"
KIND_GEOMETRY = "geometry"


class RegistryError(ValueError):
    """The registry is malformed, inconsistent with the store, or incomplete."""


def _now(clock: str | None) -> str:
    if clock:
        return clock
    import datetime as _dt

    return _dt.datetime.now(_dt.timezone.utc).replace(microsecond=0).isoformat()


@dataclass(frozen=True)
class AssetEntry:
    """One verifiable file in the durable store."""

    asset_id: str
    kind: str
    sha256: str
    bytes: int
    hull: str = ""
    status: str = "unknown"
    source: str = ""
    duplicate_of: str = ""
    note: str = ""

    @property
    def uri(self) -> str:
        return asset_uri(self.asset_id)

    @property
    def verified(self) -> bool:
        return self.status == "verified"

    def as_dict(self) -> dict[str, Any]:
        return {
            "asset_id": self.asset_id,
            "uri": self.uri,
            "kind": self.kind,
            "hull": self.hull,
            "sha256": self.sha256,
            "bytes": self.bytes,
            "status": self.status,
            "verified": self.verified,
            "source": self.source,
            "duplicate_of": self.duplicate_of,
            "note": self.note,
        }


@dataclass(frozen=True)
class ShipBinding:
    """A fleet-spec ship and the asset its geometry SHA-256 points at."""

    ship: str
    asset_id: str
    iges_sha256: str
    sha_match: bool
    benchmark_verdict: str = ""
    benchmark_rel_error_percent: float | None = None

    @property
    def uri(self) -> str:
        return asset_uri(self.asset_id)

    def as_dict(self) -> dict[str, Any]:
        return {
            "ship": self.ship,
            "asset_id": self.asset_id,
            "uri": self.uri,
            "iges_sha256": self.iges_sha256,
            "sha_match": self.sha_match,
            "benchmark_verdict": self.benchmark_verdict,
            "benchmark_rel_error_percent": self.benchmark_rel_error_percent,
        }


@dataclass
class AssetRegistry:
    """The index over the store, plus the operations consumers need."""

    assets: dict[str, AssetEntry] = field(default_factory=dict)
    bindings: dict[str, ShipBinding] = field(default_factory=dict)
    hulls: dict[str, dict[str, Any]] = field(default_factory=dict)
    snapshot: dict[str, Any] = field(default_factory=dict)
    store_root: str = ""
    schema: str = REGISTRY_SCHEMA

    # ------------------------------------------------------------------ lookups
    def by_kind(self, kind: str) -> list[AssetEntry]:
        return [entry for entry in self.assets.values() if entry.kind == kind]

    def by_hull(self, hull: str) -> list[AssetEntry]:
        return [entry for entry in self.assets.values() if entry.hull == hull]

    def entry(self, asset_id: str) -> AssetEntry:
        key = asset_id_of(asset_id) or asset_id.strip("/")
        try:
            return self.assets[key]
        except KeyError:
            raise RegistryError(f"asset {asset_id!r} is not registered") from None

    def ship(self, ship: str) -> ShipBinding:
        try:
            return self.bindings[ship]
        except KeyError:
            raise RegistryError(f"ship {ship!r} has no geometry binding in the registry") from None

    # -------------------------------------------------------------- resolution
    def resolve(self, name: str, *, store: AssetStore | None = None, require_verified: bool = False) -> str:
        """Resolve an asset id/URI or a ship name to a verified absolute path.

        Byte integrity is not optional: a consumer holding a name gets either bytes
        whose SHA-256 equals the registry record, or an explicit refusal (``AssetError``
        with ``reason`` in ``{missing, sha_mismatch, size_mismatch, unverified}``).

        ``require_verified`` is a *separate* gate. Byte integrity and physics
        acceptance are different claims: KVLCC2's bytes are hash-clean while its
        displacement is still unfixed (``status: pending`` in the store manifest), and a
        reader that needs the geometry must still get it. Gates that publish a number
        pass ``require_verified=True`` and refuse anything not yet proven.
        """
        store = store or AssetStore(self.store_root or None)
        if asset_id_of(name) or "/" in name:
            asset_id = asset_id_of(name) or name.strip("/")
        else:
            asset_id = self.ship(name).asset_id
        entry = self.entry(asset_id)
        canonical = self.entry(entry.duplicate_of) if entry.duplicate_of else entry
        path = store.verify_asset(asset_id, entry.sha256, entry.bytes)
        asset_id, entry = canonical.asset_id, canonical
        if require_verified and not entry.verified:
            # Bytes agree with the registry; the manifest itself never declared them
            # verified. Refuse to launder "present" into "proven".
            raise AssetError(
                f"{asset_id} is present and hash-clean but not verified ({entry.note or 'no note'})",
                reason="unverified",
            )
        return path

    def find_by_sha(self, sha256: str, *, store: AssetStore | None = None) -> str | None:
        for entry in self.assets.values():
            if entry.sha256 == sha256:
                return entry.asset_id
        store = store or AssetStore(self.store_root or None)
        return store.find_by_sha(sha256)

    # ------------------------------------------------------------ verification
    def verify(self, *, store: AssetStore | None = None) -> dict[str, Any]:
        """Re-hash every registered byte; report per-asset state, never raise per asset."""
        store = store or AssetStore(self.store_root or None)
        results: dict[str, str] = {}
        for asset_id, entry in sorted(self.assets.items()):
            try:
                store.verify_asset(asset_id, entry.sha256, entry.bytes)
                results[asset_id] = "ok" if entry.verified else "present_unverified"
            except AssetError as exc:
                results[asset_id] = exc.reason
        counts: dict[str, int] = {}
        for state in results.values():
            counts[state] = counts.get(state, 0) + 1
        return {
            "schema": "tensorcad.assetlib.verify.v1",
            "store_root": str(store.root),
            "assets": len(results),
            "states": dict(sorted(counts.items())),
            "problems": {aid: st for aid, st in results.items() if st not in {"ok", "present_unverified"}},
            "detail": results,
            "ok": not any(st in {"missing", "sha_mismatch", "size_mismatch"} for st in results.values()),
        }

    # ------------------------------------------------------------------ growth
    def growth(self, previous: Mapping[str, Any] | None = None) -> dict[str, Any]:
        """Differentiate this snapshot against ``previous`` (defaults to ``snapshot``)."""
        current = self.snapshot_dict()
        before = dict(previous or {})
        prev_assets = int(before.get("assets", 0) or 0)
        prev_bytes = int(before.get("bytes", 0) or 0)
        prev_ids = set(before.get("asset_ids", []) or [])
        cur_ids = set(self.assets)
        added = sorted(cur_ids - prev_ids) if prev_ids else sorted(cur_ids)
        removed = sorted(prev_ids - cur_ids)
        cur_bytes = int(current["bytes"])
        cov_prev = before.get("coverage", None)
        return {
            "schema": "tensorcad.assetlib.growth.v1",
            "current": current,
            "previous": before or None,
            "delta": {
                "assets": int(current["assets"]) - prev_assets,
                "bytes": cur_bytes - prev_bytes,
                "verified": int(current["verified"]) - int(before.get("verified", 0) or 0),
                "spec_bindings": int(current["spec_bindings"]) - int(before.get("spec_bindings", 0) or 0),
                "coverage": None
                if cov_prev is None
                else round(float(current["coverage"]) - float(cov_prev), 6),
            },
            "added": added,
            "removed": removed,
            "basis": "two registry snapshots; bytes = sum of registered file sizes",
        }

    def snapshot_dict(self) -> dict[str, Any]:
        """A compact, comparable summary of the library right now."""
        total_bytes = sum(entry.bytes for entry in self.assets.values())
        by_kind: dict[str, dict[str, int]] = {}
        for entry in self.assets.values():
            slot = by_kind.setdefault(entry.kind, {"assets": 0, "bytes": 0, "verified": 0})
            slot["assets"] += 1
            slot["bytes"] += entry.bytes
            slot["verified"] += int(entry.verified)
        by_hull: dict[str, dict[str, int]] = {}
        for entry in self.assets.values():
            if not entry.hull:
                continue
            slot = by_hull.setdefault(entry.hull, {"assets": 0, "bytes": 0})
            slot["assets"] += 1
            slot["bytes"] += entry.bytes
        verified = sum(int(entry.verified) for entry in self.assets.values())
        bound = sum(int(binding.sha_match) for binding in self.bindings.values())
        return {
            "assets": len(self.assets),
            "bytes": total_bytes,
            "verified": verified,
            "verified_fraction": round(verified / len(self.assets), 6) if self.assets else 0.0,
            "by_kind": {k: by_kind[k] for k in sorted(by_kind)},
            "by_hull": {k: by_hull[k] for k in sorted(by_hull)},
            "spec_bindings": bound,
            "spec_ships": len(self.bindings),
            "coverage": round(bound / len(self.bindings), 6) if self.bindings else 0.0,
            "asset_ids": sorted(self.assets),
        }

    # --------------------------------------------------------------------- I/O
    def as_dict(self) -> dict[str, Any]:
        return {
            "schema": self.schema,
            "store_root": self.store_root or "$TENSORCAD_ASSET_ROOT",
            "license_note": (
                "Fingerprints and pointers only: the geometry bytes are publisher material "
                "without a redistribution grant and are NOT stored in this repository."
            ),
            "assets": [self.assets[key].as_dict() for key in sorted(self.assets)],
            "bindings": [self.bindings[key].as_dict() for key in sorted(self.bindings)],
            "hulls": {key: self.hulls[key] for key in sorted(self.hulls)},
            "snapshot": self.snapshot,
        }

    def to_json(self, *, indent: int = 1) -> str:
        return json.dumps(self.as_dict(), indent=indent, sort_keys=False) + "\n"


def _load_json(path: str | Path) -> Any:
    return json.loads(Path(path).read_text())


def build_registry(
    *,
    store: AssetStore | None = None,
    specs_path: str | Path | None = None,
    benchmark_path: str | Path | None = None,
    clock: str | None = None,
) -> AssetRegistry:
    """Read the store, join it with the fleet specs and the benchmark, verify hashes.

    Nothing is assumed: the per-hull manifests are re-checked against the bytes, and a
    spec ship is bound only when its recorded ``iges_sha256`` equals the bytes on disk.
    A file present in the store but absent from every spec is still registered (it is
    library growth), it simply carries no ship binding.
    """
    store = store or AssetStore()
    assets: dict[str, AssetEntry] = {}
    hulls: dict[str, dict[str, Any]] = {}
    problems: list[str] = []

    specs_ships: list[dict[str, Any]] = []
    if specs_path is not None and Path(specs_path).is_file():
        payload = _load_json(specs_path)
        specs_ships = list(payload.get("ships", []) if isinstance(payload, dict) else payload)

    verdicts: dict[str, dict[str, Any]] = {}
    if benchmark_path is not None and Path(benchmark_path).is_file():
        payload = _load_json(benchmark_path)
        verdicts = dict(payload.get("hulls", {}) or {})

    spec_sha = {str(entry.get("name", "")): str(entry.get("iges_sha256", "")) for entry in specs_ships}
    spec_hull_owner = {sha: ship for ship, sha in spec_sha.items() if sha}

    for hull_dir in store.hull_dirs():
        hull = hull_dir.name
        try:
            check = store.verify_manifest(hull_dir)
        except AssetError as exc:
            problems.append(f"{hull}: {exc}")
            continue
        manifest_path = Path(check.manifest_path)
        manifest = _load_json(manifest_path)
        declared_status = str(manifest.get("status", "unknown"))
        note = str(manifest.get("verify", ""))
        source = str(manifest.get("source", ""))
        # A hull is "verified" only when its own manifest says so AND its bytes match.
        hull_status = "verified" if declared_status == "verified" and check.ok else declared_status
        if not check.ok:
            hull_status = "broken"
        first_seen: dict[str, str] = {}
        for rel, digest in sorted(check.actual.items()):
            asset_id = f"geometry/hulls/{hull}/{rel}"
            if digest in first_seen:
                dup, status, entry_note = first_seen[digest], "duplicate", f"same bytes as {first_seen[digest]}"
            else:
                first_seen[digest] = asset_id
                dup, status, entry_note = "", hull_status, note
            assets[asset_id] = AssetEntry(
                asset_id=asset_id,
                kind=KIND_GEOMETRY,
                sha256=digest,
                bytes=(hull_dir / rel).stat().st_size,
                hull=hull,
                status=status,
                source=source,
                duplicate_of=dup,
                note=entry_note,
            )
        for rel in check.missing:
            problems.append(f"{hull}: declared {rel} is missing from the store")
        for rel in check.mismatched:
            problems.append(f"{hull}: {rel} sha256 differs from the manifest claim")
        hulls[hull] = {
            "hull": hull,
            "status": hull_status,
            "manifest": str(manifest_path),
            "manifest_declared_status": declared_status,
            "declared_files": len(check.declared),
            "files_on_disk": len(check.actual),
            "verified": hull_status == "verified",
            "license": str(manifest.get("license", "")),
            "provenance": str(manifest.get("provenance", "")),
            "source": source,
            "verification_note": note,
            "captured": str(manifest.get("captured", "")),
            "ship_bound": hull in spec_sha,
        }

    bindings: dict[str, ShipBinding] = {}
    for entry in specs_ships:
        ship = str(entry.get("name", ""))
        sha = str(entry.get("iges_sha256", ""))
        if not ship:
            continue
        match_id = ""
        if sha:
            for asset_id, asset in sorted(assets.items(), key=lambda kv: (bool(kv[1].duplicate_of), kv[0])):
                if asset.sha256 == sha and asset.hull == hull_of(asset_id):
                    match_id = asset_id
                    break
            if not match_id:
                match_id = store.find_by_sha(sha) or ""
        verdict = verdicts.get(ship, {}) or {}
        bindings[ship] = ShipBinding(
            ship=ship,
            asset_id=match_id,
            iges_sha256=sha,
            sha_match=bool(match_id),
            benchmark_verdict=str(verdict.get("verdict", "")),
            benchmark_rel_error_percent=verdict.get("rel_error_percent"),
        )
        if not match_id:
            problems.append(f"spec ship {ship}: iges_sha256 {sha or '<none>'} matches no asset in the store")

    registry = AssetRegistry(
        assets=assets,
        bindings=bindings,
        hulls=hulls,
        store_root=str(store.root),
    )
    snapshot = registry.snapshot_dict()
    snapshot["generated_at"] = _now(clock)
    snapshot["problems"] = sorted(problems)
    registry.snapshot = snapshot
    return registry


def load_registry(path: str | Path) -> AssetRegistry:
    payload = _load_json(path)
    if payload.get("schema") != REGISTRY_SCHEMA:
        raise RegistryError(f"{path}: schema {payload.get('schema')!r} != {REGISTRY_SCHEMA!r}")
    assets = {}
    for row in payload.get("assets", []):
        entry = AssetEntry(
            asset_id=str(row["asset_id"]),
            kind=str(row.get("kind", KIND_GEOMETRY)),
            sha256=str(row["sha256"]),
            bytes=int(row["bytes"]),
            hull=str(row.get("hull", "")),
            status=str(row.get("status", "unknown")),
            source=str(row.get("source", "")),
            duplicate_of=str(row.get("duplicate_of", "")),
            note=str(row.get("note", "")),
        )
        assets[entry.asset_id] = entry
    bindings = {}
    for row in payload.get("bindings", []):
        binding = ShipBinding(
            ship=str(row["ship"]),
            asset_id=str(row.get("asset_id", "")),
            iges_sha256=str(row.get("iges_sha256", "")),
            sha_match=bool(row.get("sha_match")),
            benchmark_verdict=str(row.get("benchmark_verdict", "")),
            benchmark_rel_error_percent=row.get("benchmark_rel_error_percent"),
        )
        bindings[binding.ship] = binding
    return AssetRegistry(
        assets=assets,
        bindings=bindings,
        hulls=dict(payload.get("hulls", {})),
        snapshot=dict(payload.get("snapshot", {})),
        store_root=str(payload.get("store_root", "")),
        schema=str(payload.get("schema", REGISTRY_SCHEMA)),
    )


def write_registry(registry: AssetRegistry, path: str | Path, *, history: str | Path | None = None) -> Path:
    """Write the registry, and optionally append the snapshot to a JSONL history."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(registry.to_json())
    if history is not None:
        hist = Path(history)
        hist.parent.mkdir(parents=True, exist_ok=True)
        line = json.dumps(
            {
                key: value
                for key, value in registry.snapshot.items()
                if key != "asset_ids"  # keep the history small; ids are in the registry
            },
            sort_keys=True,
        )
        existing = hist.read_text().splitlines() if hist.is_file() else []
        if line not in existing:
            with open(hist, "a", encoding="utf-8") as handle:
                handle.write(line + "\n")
    return target


def latest_history(path: str | Path) -> dict[str, Any] | None:
    hist = Path(path)
    if not hist.is_file():
        return None
    lines = [line for line in hist.read_text().splitlines() if line.strip()]
    if not lines:
        return None
    return json.loads(lines[-1])