"""Tests for the geometry asset library: store, registry, pointer repair, consumer.

The invariant under test is the one the mission cares about: *a name either resolves to
bytes with the recorded SHA-256 or it fails loudly with a reason*. Every test that
matters therefore tampers with bytes or deletes them and asserts the failure mode,
rather than asserting that a function returned without raising.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from tensorcad.assetlib.consumer import resolve_pointer  # noqa: E402
from tensorcad.assetlib.registry import (  # noqa: E402
    AssetRegistry,
    build_registry,
    load_registry,
    write_registry,
)
from tensorcad.assetlib.specs import repair_spec_file, scan_spec_pointers  # noqa: E402
from tensorcad.assetlib.store import AssetError, AssetStore, sha256_bytes, sha256_file  # noqa: E402

STORE_SPECS = REPO / "benchmarks" / "verified" / "baseline_fleet" / "fleet_specs.json"
STORE_BENCH = REPO / "benchmarks" / "verified" / "baseline_fleet" / "result.json"
REGISTRY = REPO / "assets" / "geometry_registry.json"
REAL_STORE = AssetStore()


def _hull(root: Path, hull: str, files: dict[str, bytes], *, status: str = "verified") -> dict[str, str]:
    """Create ``root/geometry/hulls/<hull>/`` with a manifest, return the sha map."""
    base = root / "geometry" / "hulls" / hull
    base.mkdir(parents=True, exist_ok=True)
    shas: dict[str, str] = {}
    for name, payload in files.items():
        target = base / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(payload)
        shas[name] = sha256_bytes(payload)
    (base / "MANIFEST.json").write_text(
        json.dumps(
            {
                "hull": hull,
                "files": shas,
                "status": status,
                "verify": f"test fixture, status={status}",
                "source": "unit-test fixture",
                "license": "none",
                "captured": "2026-11-14",
            },
            indent=1,
        )
    )
    return shas


def _specs(path: Path, ships: list[dict], evidence: dict) -> Path:
    path.write_text(json.dumps({"ships": ships, "evidence": evidence}, indent=2))
    return path


@pytest.fixture()
def store(tmp_path: Path) -> AssetStore:
    root = tmp_path / "assets"
    _hull(root, "ALPHA", {"alpha.igs": b"ALPHA-IGES", "alpha/mirror.igs": b"ALPHA-IGES", "memo.txt": b"memo"})
    _hull(root, "BETA", {"beta.igs": b"BETA-IGES"}, status="pending")
    return AssetStore(root)


def test_registry_records_every_byte_and_marks_duplicates(store: AssetStore) -> None:
    registry = build_registry(store=store)
    assert set(registry.assets) == {
        "geometry/hulls/ALPHA/alpha.igs",
        "geometry/hulls/ALPHA/alpha/mirror.igs",
        "geometry/hulls/ALPHA/memo.txt",
        "geometry/hulls/BETA/beta.igs",
    }
    mirror = registry.entry("geometry/hulls/ALPHA/alpha/mirror.igs")
    assert mirror.status == "duplicate"
    assert mirror.duplicate_of == "geometry/hulls/ALPHA/alpha.igs"
    # the canonical copy is the one the registry trusts
    assert registry.entry("geometry/hulls/ALPHA/alpha.igs").verified
    snapshot = registry.snapshot
    assert snapshot["assets"] == 4
    assert snapshot["bytes"] == sum(len(p) for p in (b"ALPHA-IGES", b"ALPHA-IGES", b"memo", b"BETA-IGES"))
    assert snapshot["by_kind"]["geometry"]["assets"] == 4
    assert snapshot["by_hull"]["ALPHA"]["assets"] == 3


def test_resolve_verifies_bytes_and_honours_the_verified_gate(store: AssetStore) -> None:
    registry = build_registry(store=store)
    assert Path(registry.resolve("geometry/hulls/ALPHA/alpha.igs")) == store.root / "geometry/hulls/ALPHA/alpha.igs"
    assert resolve_pointer("asset://geometry/hulls/ALPHA/alpha.igs", registry=registry).endswith("alpha.igs")
    # BETA's bytes are hash-clean but its manifest never claimed physics verification:
    # reading is allowed, publishing a number from it is not.
    registry.resolve("geometry/hulls/BETA/beta.igs")
    with pytest.raises(AssetError) as excinfo:
        registry.resolve("geometry/hulls/BETA/beta.igs", require_verified=True)
    assert excinfo.value.reason == "unverified"


def test_tampered_byte_is_detected_with_the_reason(store: AssetStore) -> None:
    registry = build_registry(store=store)
    target = store.root / "geometry" / "hulls" / "ALPHA" / "alpha.igs"
    original = target.read_bytes()
    # same length, one byte different: only a hash can catch this
    target.write_bytes(original[:-1] + b"X")
    report = registry.verify(store=store)
    assert report["ok"] is False
    assert report["problems"]["geometry/hulls/ALPHA/alpha.igs"] == "sha_mismatch"
    with pytest.raises(AssetError) as excinfo:
        registry.resolve("geometry/hulls/ALPHA/alpha.igs", store=store)
    assert excinfo.value.reason == "sha_mismatch"
    # a different length is caught one step earlier, and named differently
    target.write_bytes(original + b"!")
    assert registry.verify(store=store)["problems"]["geometry/hulls/ALPHA/alpha.igs"] == "size_mismatch"


def test_deleted_byte_is_detected_and_named(store: AssetStore) -> None:
    registry = build_registry(store=store)
    (store.root / "geometry" / "hulls" / "ALPHA" / "memo.txt").unlink()
    report = registry.verify(store=store)
    assert report["problems"]["geometry/hulls/ALPHA/memo.txt"] == "missing"
    # the store manifest check flags the same loss from the other direction
    check = store.verify_manifest(hull="ALPHA")
    assert not check.ok
    assert check.missing == ("memo.txt",)


def test_growth_is_computed_from_two_snapshots(store: AssetStore) -> None:
    registry = build_registry(store=store)
    previous = dict(registry.snapshot)
    previous.update({"assets": 2, "bytes": 10, "verified": 1, "spec_bindings": 0, "coverage": 0.0,
                     "asset_ids": ["geometry/hulls/ALPHA/alpha.igs", "geometry/hulls/OLD/old.igs"]})
    growth = registry.growth(previous)
    assert growth["delta"]["assets"] == 2
    assert growth["delta"]["bytes"] == registry.snapshot["bytes"] - 10
    assert growth["delta"]["coverage"] == 0.0  # no spec ships bound in this fixture
    assert "geometry/hulls/BETA/beta.igs" in growth["added"]
    assert growth["removed"] == ["geometry/hulls/OLD/old.igs"]


def test_scan_classifies_and_repair_repoints_by_content(store: AssetStore, tmp_path: Path) -> None:
    shas = {
        name: sha256_file(store.root / "geometry" / "hulls" / "ALPHA" / name)
        for name in ("alpha.igs", "memo.txt")
    }
    specs = _specs(
        tmp_path / "fleet_specs.json",
        ships=[
            {
                "name": "ALPHA",
                "draft": 10.0,
                "length_pp": 100.0,
                "beam": 20.0,
                "iges_sha256": shas["alpha.igs"],
                "iges_path": "/root/tp-agent/work/m00000/artifacts/raw/ALPHA/alpha.igs",
                "published_source": "artifacts/raw/pages/alpha.html",
            }
        ],
        evidence={
            "ALPHA": {
                "iges": {"path": "artifacts/raw/ALPHA/alpha.igs", "sha256": shas["alpha.igs"]},
                "page": {"key": "alpha_page", "path": "artifacts/raw/pages/alpha.html", "sha256": "ab" * 32},
                "archive": {"path": "artifacts/raw/ALPHA.zip", "sha256": "cd" * 32},
                "published_fields": {"length_pp": {"value": 100.0, "page": "artifacts/raw/pages/alpha.html",
                                                   "page_sha256": "ab" * 32}},
            }
        },
    )
    registry = build_registry(store=store)
    before = scan_spec_pointers(specs, store=store, registry=registry, repo_root=tmp_path)
    assert {p.where: p.status for p in before.pointers} == {
        "ships[0](ALPHA).iges_path": "dead",
        "ships[0](ALPHA).published_source": "withheld",
        "evidence.ALPHA.iges.path": "dead",
        "evidence.ALPHA.page.path": "withheld",
        "evidence.ALPHA.archive.path": "withheld",
        "evidence.ALPHA.published_fields.length_pp.page": "withheld",
    }
    result = repair_spec_file(
        specs,
        store=store,
        registry=registry,
        repo_root=tmp_path,
        url_map={"alpha_page": "https://example.org/alpha.html"},
        archive_url_map={"ALPHA.zip": "https://example.org/ALPHA.zip"},
    )
    assert result["unresolved"] == []
    repaired = json.loads(specs.read_text())
    assert repaired["ships"][0]["iges_path"] == "asset://geometry/hulls/ALPHA/alpha.igs"
    assert repaired["ships"][0]["published_source"] == "https://example.org/alpha.html"
    assert repaired["ships"][0]["published_source_legacy"] == "artifacts/raw/pages/alpha.html"
    assert repaired["ships"][0]["pointer_status"]["published_source"] == "cache-lost"
    assert repaired["evidence"]["ALPHA"]["archive"]["path"] == "https://example.org/ALPHA.zip"
    after = scan_spec_pointers(specs, store=store, registry=registry, repo_root=tmp_path)
    assert after.dead == [] and after.withheld == []
    assert len(after.resolved) == 2 and len(after.cited) == 4


def test_repair_refuses_to_invent_a_path_when_no_sha_is_recorded(store: AssetStore, tmp_path: Path) -> None:
    specs = _specs(
        tmp_path / "specs.json",
        ships=[{"name": "GHOST", "draft": 1.0, "length_pp": 1.0, "beam": 1.0,
                "iges_sha256": "", "iges_path": "/root/tp-agent/work/dead/ghost.igs"}],
        evidence={},
    )
    result = repair_spec_file(specs, store=store, registry=build_registry(store=store), repo_root=tmp_path)
    assert len(result["unresolved"]) == 1
    assert json.loads(specs.read_text())["ships"][0]["iges_path"] == "/root/tp-agent/work/dead/ghost.igs"


def test_registry_roundtrip(store: AssetStore, tmp_path: Path) -> None:
    registry = build_registry(store=store)
    target = write_registry(registry, tmp_path / "registry.json", history=tmp_path / "history.jsonl")
    loaded = load_registry(target)
    assert loaded.snapshot_dict()["assets"] == registry.snapshot_dict()["assets"]
    assert loaded.entry("geometry/hulls/ALPHA/alpha.igs").sha256 == registry.entry(
        "geometry/hulls/ALPHA/alpha.igs"
    ).sha256
    lines = (tmp_path / "history.jsonl").read_text().splitlines()
    assert len(lines) == 1 and json.loads(lines[0])["assets"] == 4


def test_import_asset_refuses_the_wrong_bytes(store: AssetStore, tmp_path: Path) -> None:
    source = tmp_path / "new.igs"
    source.write_bytes(b"NEW-BYTES")
    with pytest.raises(AssetError) as excinfo:
        store.import_asset(source, "geometry/hulls/GAMMA/new.igs", expected_sha256="00" * 32)
    assert excinfo.value.reason == "sha_mismatch"
    assert not (store.root / "geometry" / "hulls" / "GAMMA" / "new.igs").exists()
    record = store.import_asset(source, "geometry/hulls/GAMMA/new.igs", expected_sha256=sha256_bytes(b"NEW-BYTES"))
    assert record["action"] == "copied"
    assert (store.root / "geometry" / "hulls" / "GAMMA" / "new.igs").read_bytes() == b"NEW-BYTES"
    assert store.import_asset(source, "geometry/hulls/GAMMA/new.igs")["action"] == "identical"


def test_store_rejects_ids_that_escape_the_root(store: AssetStore) -> None:
    with pytest.raises(AssetError):
        store.path_of("../etc/passwd")


# ---------------------------------------------------------------- committed state
needs_store = pytest.mark.skipif(not REAL_STORE.exists(), reason="durable asset store not mounted")


@needs_store
def test_committed_specs_have_no_dead_pointers() -> None:
    registry = load_registry(REGISTRY)
    report = scan_spec_pointers(STORE_SPECS, store=REAL_STORE, registry=registry, repo_root=REPO)
    assert report.dead == [], [p.as_dict() for p in report.dead]
    assert report.withheld == [], [p.as_dict() for p in report.withheld]
    geometry = [p for p in report.pointers if p.kind == "geometry"]
    assert geometry and all(p.status == "resolved" for p in geometry)


@needs_store
def test_committed_registry_matches_the_store_bytes() -> None:
    registry = load_registry(REGISTRY)
    report = registry.verify(store=REAL_STORE)
    assert report["ok"], report["problems"]
    # a registry that has drifted from the store is a lie in git
    rebuilt = build_registry(store=REAL_STORE, specs_path=STORE_SPECS, benchmark_path=STORE_BENCH)
    assert rebuilt.snapshot_dict()["asset_ids"] == registry.snapshot_dict()["asset_ids"]
    assert rebuilt.snapshot_dict()["bytes"] == registry.snapshot_dict()["bytes"]


@needs_store
def test_every_committed_geometry_pointer_carries_the_spec_sha() -> None:
    payload = json.loads(STORE_SPECS.read_text())
    registry = load_registry(REGISTRY)
    for ship in payload["ships"]:
        pointer = ship["iges_path"]
        assert pointer.startswith("asset://"), pointer
        resolved = resolve_pointer(pointer, ship["iges_sha256"], registry=registry)
        assert sha256_file(resolved) == ship["iges_sha256"], ship["name"]


@needs_store
def test_fleet_baseline_reads_geometry_through_the_registry() -> None:
    from tensorcad.fleet_baseline import specs_from_json

    specs = specs_from_json(STORE_SPECS)
    assert specs, "no specs parsed"
    for spec in specs:
        assert Path(spec.iges_path).is_file(), spec.iges_path
        assert sha256_file(spec.iges_path) == spec.iges_sha256, spec.name
        assert Path(spec.iges_path).is_relative_to(REAL_STORE.root), spec.iges_path