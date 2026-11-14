"""Access to the durable geometry store.

The store is a plain directory tree outside git::

    $TENSORCAD_ASSET_ROOT/geometry/hulls/<HULL>/<file...>
    $TENSORCAD_ASSET_ROOT/geometry/hulls/<HULL>/MANIFEST.json

``MANIFEST.json`` is written by whoever imported the hull and carries the SHA-256 of
every file plus the provenance and the verification note. It is treated here as the
*authoritative claim*; :meth:`AssetStore.verify_manifest` checks that claim against
the bytes, because a manifest that nobody re-checks is just a comment.

Asset ids are store-relative POSIX paths (``geometry/hulls/KCS/KCS_hull_SVA.igs``)
and asset *URIs* are the same string behind the ``asset://`` scheme, which is what a
committed JSON file is allowed to contain: it names the byte without claiming a
machine-local absolute path.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping

SCHEME = "asset://"
DEFAULT_ROOT = "/root/tp-agent/assets"
ENV_VAR = "TENSORCAD_ASSET_ROOT"


class AssetError(LookupError):
    """An asset is missing, malformed, or does not match its recorded SHA-256."""

    def __init__(self, message: str, *, reason: str = "unavailable") -> None:
        super().__init__(message)
        self.reason = reason


def default_store_root() -> Path:
    """Where the platform keeps the durable store (overridable for tests/CI)."""
    return Path(os.environ.get(ENV_VAR, DEFAULT_ROOT))


def sha256_file(path: str | Path, chunk: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            block = handle.read(chunk)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def asset_uri(asset_id: str) -> str:
    return SCHEME + asset_id.strip("/")


def asset_id_of(pointer: str) -> str | None:
    """Return the asset id when ``pointer`` is an ``asset://`` URI, else ``None``."""
    text = str(pointer).strip()
    if text.startswith(SCHEME):
        return text[len(SCHEME) :].strip("/")
    return None


def hull_of(asset_id: str) -> str:
    """``geometry/hulls/KCS/x.igs`` -> ``KCS`` (empty string if not a hull asset)."""
    parts = asset_id.strip("/").split("/")
    if len(parts) >= 3 and parts[0] == "geometry" and parts[1] == "hulls":
        return parts[2]
    return ""


@dataclass(frozen=True)
class ManifestCheck:
    """Result of re-checking one store manifest against the bytes on disk."""

    hull: str
    manifest_path: str
    declared: Mapping[str, str]
    actual: Mapping[str, str]
    missing: tuple[str, ...]
    mismatched: tuple[str, ...]
    undeclared: tuple[str, ...]

    @property
    def ok(self) -> bool:
        return not (self.missing or self.mismatched)

    def as_dict(self) -> dict[str, object]:
        return {
            "hull": self.hull,
            "manifest": self.manifest_path,
            "declared_files": len(self.declared),
            "missing": list(self.missing),
            "mismatched": list(self.mismatched),
            "undeclared": list(self.undeclared),
            "ok": self.ok,
        }


class AssetStore:
    """A directory tree of durable asset bytes plus the per-hull manifests."""

    def __init__(self, root: str | Path | None = None) -> None:
        self.root = Path(root) if root is not None else default_store_root()

    # ------------------------------------------------------------------ layout
    def exists(self) -> bool:
        return self.root.is_dir()

    def hulls_dir(self) -> Path:
        return self.root / "geometry" / "hulls"

    def hull_dirs(self) -> list[Path]:
        base = self.hulls_dir()
        if not base.is_dir():
            return []
        return sorted(p for p in base.iterdir() if p.is_dir())

    def path_of(self, asset_id: str) -> Path:
        rel = asset_id.strip("/")
        if rel.startswith("/") or ".." in Path(rel).parts:
            raise AssetError(f"asset id {asset_id!r} escapes the store", reason="bad_id")
        return self.root / rel

    def relpath(self, path: str | Path) -> str:
        """Store-relative POSIX path (the asset id) for an absolute path."""
        return Path(path).resolve().relative_to(self.root.resolve()).as_posix()

    # ------------------------------------------------------------- verification
    def verify_manifest(self, hull_dir: Path | None = None, *, hull: str | None = None) -> ManifestCheck:
        target = hull_dir
        if target is None:
            if hull is None:
                raise AssetError("verify_manifest needs hull_dir or hull", reason="bad_call")
            target = self.hulls_dir() / hull
        manifest_path = Path(target) / "MANIFEST.json"
        if not manifest_path.is_file():
            raise AssetError(f"{manifest_path} is missing", reason="no_manifest")
        payload = json.loads(manifest_path.read_text())
        declared = {str(k): str(v) for k, v in dict(payload.get("files", {})).items()}
        actual: dict[str, str] = {}
        for path in sorted(Path(target).rglob("*")):
            if path.is_file() and path.name != "MANIFEST.json":
                actual[path.relative_to(Path(target)).as_posix()] = sha256_file(path)
        missing = tuple(sorted(set(declared) - set(actual)))
        mismatched = tuple(sorted(name for name in set(declared) & set(actual) if declared[name] != actual[name]))
        undeclared = tuple(sorted(set(actual) - set(declared)))
        return ManifestCheck(
            hull=Path(target).name,
            manifest_path=str(manifest_path),
            declared=declared,
            actual=actual,
            missing=missing,
            mismatched=mismatched,
            undeclared=undeclared,
        )

    def verify_all(self) -> list[ManifestCheck]:
        checks = []
        for hull_dir in self.hull_dirs():
            try:
                checks.append(self.verify_manifest(hull_dir))
            except AssetError:
                continue
        return checks

    def verify_asset(self, asset_id: str, expected_sha256: str, expected_bytes: int | None = None) -> str:
        """Return the resolved path when bytes match, else raise :class:`AssetError`."""
        path = self.path_of(asset_id)
        if not path.is_file():
            raise AssetError(f"{asset_id} is missing from the store", reason="missing")
        if expected_bytes is not None and path.stat().st_size != int(expected_bytes):
            raise AssetError(
                f"{asset_id}: {path.stat().st_size} bytes on disk != {expected_bytes} recorded",
                reason="size_mismatch",
            )
        actual = sha256_file(path)
        if actual != expected_sha256:
            raise AssetError(
                f"{asset_id}: sha256 {actual} != {expected_sha256} recorded",
                reason="sha_mismatch",
            )
        return str(path)

    def find_by_sha(self, sha256: str, *, prefer: Iterable[str] = ()) -> str | None:
        """Locate an asset id by content hash (survives renames, not deletions)."""
        if not sha256:
            return None
        preferred = [pref for pref in prefer if pref]
        ordered: list[Path] = []
        for pref in preferred:
            base = self.path_of(pref)
            if base.is_file():
                ordered.append(base)
            elif base.is_dir():
                ordered.extend(sorted(p for p in base.rglob("*") if p.is_file()))
        ordered.extend(
            sorted(p for p in self.root.rglob("*") if p.is_file() and p.name != "MANIFEST.json")
        )
        seen: set[Path] = set()
        for path in ordered:
            if path in seen:
                continue
            seen.add(path)
            if sha256_file(path) == sha256:
                return self.relpath(path)
        return None

    # ------------------------------------------------------------------ writing
    def import_asset(
        self,
        source: str | Path,
        asset_id: str,
        *,
        expected_sha256: str | None = None,
        overwrite: bool = False,
    ) -> dict[str, object]:
        """Copy one file into the store atomically, iff its SHA-256 is what we expect.

        Returns a record ``{asset_id, path, sha256, bytes, action}``. ``action`` is
        ``copied``, ``identical`` (bytes already there) or ``skipped``; a hash
        mismatch raises and writes nothing, because importing the wrong bytes under a
        trusted name is worse than not importing at all.
        """
        src = Path(source)
        if not src.is_file():
            raise AssetError(f"source {src} is not a file", reason="missing_source")
        digest = sha256_file(src)
        if expected_sha256 is not None and digest != expected_sha256:
            raise AssetError(
                f"{src}: sha256 {digest} != expected {expected_sha256}",
                reason="sha_mismatch",
            )
        dest = self.path_of(asset_id)
        record = {
            "asset_id": asset_id.strip("/"),
            "path": str(dest),
            "sha256": digest,
            "bytes": src.stat().st_size,
        }
        if dest.is_file():
            if sha256_file(dest) == digest:
                return {**record, "action": "identical"}
            if not overwrite:
                raise AssetError(
                    f"{asset_id} already exists with different bytes (pass overwrite=True)",
                    reason="conflict",
                )
        dest.parent.mkdir(parents=True, exist_ok=True)
        handle, tmp_name = tempfile.mkstemp(dir=str(dest.parent), prefix=".import-")
        os.close(handle)
        try:
            shutil.copyfile(src, tmp_name)
            if sha256_file(tmp_name) != digest:
                raise AssetError(f"{src}: copy does not reproduce sha256 {digest}", reason="copy_corrupt")
            os.replace(tmp_name, dest)
        finally:
            if os.path.exists(tmp_name):
                os.unlink(tmp_name)
        return {**record, "action": "copied"}