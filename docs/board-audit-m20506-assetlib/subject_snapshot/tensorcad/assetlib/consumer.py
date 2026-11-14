"""What a downstream module calls to get geometry bytes.

This is deliberately tiny. The contract is:

``geometry_path(ship)``      verified path for a fleet ship, from the registry;
``resolve_pointer(value, sha)``  verified path for a spec pointer, whatever form it
                             takes (``asset://`` URI, absolute path, dead work-tree path);
``open_geometry(ship)``      a binary handle on those bytes.

Resolution order for a spec pointer: ``asset://`` (registry/store) -> existing path
with a matching SHA-256 -> *content lookup by SHA-256 in the store* -> the pointer
itself, unchanged. The last step is on purpose: a spec that was written before the
store existed must keep working, and a genuinely absent file must still fail where it
always did (in the geometry reader, with its own message), not be papered over here.
"""

from __future__ import annotations

from pathlib import Path
from typing import IO, Any

from .registry import AssetRegistry, load_registry
from .store import AssetError, AssetStore, asset_id_of, sha256_file

REPO = Path(__file__).resolve().parents[2]
REGISTRY_PATH = REPO / "assets" / "geometry_registry.json"


def _default_registry(path: str | Path | None = None) -> AssetRegistry | None:
    target = Path(path) if path is not None else REGISTRY_PATH
    if not target.is_file():
        return None
    try:
        return load_registry(target)
    except Exception:
        return None


def resolve_pointer(
    value: str,
    sha256: str = "",
    *,
    store: AssetStore | None = None,
    registry: AssetRegistry | None = None,
) -> str:
    """Turn a spec pointer into a path whose SHA-256 equals ``sha256`` (when given)."""
    registry = registry if registry is not None else _default_registry()
    store = store or AssetStore((registry.store_root if registry else "") or None)
    asset_id = asset_id_of(value)
    if asset_id:
        if registry is not None:
            return registry.resolve(asset_id, store=store)
        if sha256:
            return store.verify_asset(asset_id, sha256)
        return str(store.path_of(asset_id))
    path = Path(value)
    if path.is_file():
        if sha256 and sha256_file(path) != sha256:
            raise AssetError(
                f"{value}: sha256 {sha256_file(path)} != spec {sha256}", reason="sha_mismatch"
            )
        return str(path)
    if sha256:
        found = (registry.find_by_sha(sha256, store=store) if registry else None) or store.find_by_sha(sha256)
        if found:
            return store.verify_asset(found, sha256)
    return value


def geometry_path(
    ship: str,
    *,
    registry: AssetRegistry | None = None,
    store: AssetStore | None = None,
    require_verified: bool = False,
) -> str:
    """Path to ``ship``'s registered geometry whose bytes match the registry SHA-256.

    ``require_verified=True`` additionally demands that the store manifest calls the
    hull verified (physics accepted), not merely hash-clean; see
    :meth:`tensorcad.assetlib.registry.AssetRegistry.resolve`.
    """
    registry = registry if registry is not None else _default_registry()
    if registry is None:
        raise AssetError("no asset registry available", reason="no_registry")
    return registry.resolve(ship, store=store, require_verified=require_verified)


def open_geometry(ship: str, *, registry: AssetRegistry | None = None, mode: str = "rb") -> IO[Any]:
    return open(geometry_path(ship, registry=registry), mode)