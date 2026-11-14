"""Asset library: durable geometry store, registry, growth metrics, consumer surface.

Why this package exists
-----------------------
Mission evidence used to live in *work trees*: ``fleet_specs.json`` pointed at
``/root/tp-agent/work/m18488/artifacts/baseline_fleet/raw/KCS/KCS_hull_SVA.igs``.
Work trees are disposable, so every one of those pointers becomes a dead pointer the
moment the mission closes, and a benchmark that cites a dead pointer cannot be
re-run. The platform now keeps the bytes in a durable store
(``$TENSORCAD_ASSET_ROOT/geometry/hulls/<hull>/``, default
``/root/tp-agent/assets``). This package is the layer that makes that store
*addressable and checkable* from git:

``store``
    Read-side and write-side access to the store: resolve the root, hash bytes,
    read the per-hull ``MANIFEST.json`` the store already carries, import a file
    with an expected SHA-256 and an atomic rename.
``registry``
    The versioned index (``assets/geometry_registry.json``) that joins store bytes,
    their SHA-256, the fleet-spec bindings and the benchmark verdicts; it can
    verify itself against the bytes, resolve an asset or a ship name, and compute
    the growth metrics that show whether the library is actually growing.
``specs``
    Pointer hygiene for ``fleet_specs.json``: find pointers that do not resolve,
    re-point them at ``asset://`` URIs *by SHA-256* (never by guesswork), and report
    what could not be recovered instead of pretending.
``cli``
    ``python -m tensorcad.assetlib`` → build / verify / resolve / growth / repair.

The bytes are deliberately **not** in this repo: the workshop archives carry
publisher terms without a redistribution grant (see every per-hull ``MANIFEST.json``
and ``docs/asset_library.md``). What is in the repo is the registry: SHA-256, sizes,
provenance and pointers, which is what makes a lost byte detectable instead of
silent.
"""

from __future__ import annotations

from .registry import (
    AssetEntry,
    AssetRegistry,
    RegistryError,
    REGISTRY_SCHEMA,
    build_registry,
    load_registry,
    write_registry,
)
from .store import AssetError, AssetStore, default_store_root, sha256_file
from .specs import PointerReport, repair_spec_file, scan_spec_pointers

__all__ = [
    "AssetEntry",
    "AssetError",
    "AssetRegistry",
    "AssetStore",
    "PointerReport",
    "REGISTRY_SCHEMA",
    "RegistryError",
    "build_registry",
    "default_store_root",
    "load_registry",
    "repair_spec_file",
    "scan_spec_pointers",
    "sha256_file",
    "write_registry",
]