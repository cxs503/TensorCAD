"""Hull-primary resolution: ask the library for a hull by name, not by file name.

The gap these tests close is a *silent* one. ``assets/geometry_registry.json`` held six
hulls, but only the four ships a fleet spec binds (KCS/KVLCC2/JBC/ONRT) were reachable by
name; ``consumer.geometry_path("DTC")`` raised ``RegistryError: no geometry binding`` even
though the DTC hull block was in the registry and marked verified. The library therefore
looked complete (``coverage: 1.0``) while two hulls it stored could not be used, because
``coverage`` counted spec ships only.

Two rules make those hulls addressable, and the same two rules refuse to guess:

``unique_geometry``   the hull owns exactly one non-duplicate geometry asset;
``hull_form_name``    exactly one asset's stem is ``<hull>`` or ``<hull>_hull``.

Everything else is a ``RegistryError`` that names the candidates. ``DTC_bilgekeel.igs``
is not a hull, and a resolver that hands it out because it sorts first would be worse than
one that fails.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from tensorcad.assetlib.consumer import geometry_path  # noqa: E402
from tensorcad.assetlib.registry import (  # noqa: E402
    AssetRegistry,
    RegistryError,
    build_registry,
    load_registry,
)
from tensorcad.assetlib.store import AssetError, AssetStore, sha256_bytes  # noqa: E402

REGISTRY = REPO / "assets" / "geometry_registry.json"
REAL_STORE = AssetStore()
needs_store = pytest.mark.skipif(not REAL_STORE.exists(), reason="durable asset store not mounted")


def _hull(root: Path, hull: str, files: dict[str, bytes], *, status: str = "verified") -> dict[str, str]:
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


@pytest.fixture()
def store(tmp_path: Path) -> AssetStore:
    root = tmp_path / "assets"
    # DELTA mirrors DTC: a hull byte plus appendages whose stems also start with the hull
    # name, so a prefix match would be ambiguous and only the exact-stem rule is honest.
    _hull(
        root,
        "DELTA",
        {
            "delta_hull.igs": b"DELTA-HULL",
            "delta_bilgekeel.igs": b"DELTA-BILGEKEEL",
            "delta_skeg.igs": b"DELTA-SKEG",
        },
    )
    # GAMMA mirrors KP505: one geometry file, no name relationship to the hull directory.
    _hull(root, "GAMMA", {"gamma_modelScale_250mm.igs": b"GAMMA-IGES"}, status="pending")
    # OMEGA has two hull-sized geometries and no naming rule: it must stay unresolvable.
    _hull(root, "OMEGA", {"omega_a.igs": b"OMEGA-A", "omega_b.igs": b"OMEGA-B"})
    _hull(root, "ALPHA", {"alpha.igs": b"ALPHA-IGES"})
    return AssetStore(root)


def _spec(tmp_path: Path, store: AssetStore) -> Path:
    """A fleet spec that binds ALPHA to a *different* byte than the naming rule would."""
    other = store.root / "geometry" / "hulls" / "DELTA" / "delta_skeg.igs"
    payload = {
        "ships": [
            {
                "name": "ALPHA",
                "iges_path": "asset://geometry/hulls/DELTA/delta_skeg.igs",
                "iges_sha256": sha256_bytes(other.read_bytes()),
            }
        ]
    }
    path = tmp_path / "fleet_specs.json"
    path.write_text(json.dumps(payload, indent=1))
    return path


def test_hull_form_name_rule_picks_the_hull_and_not_an_appendage(store: AssetStore) -> None:
    registry = build_registry(store=store)
    entry, rule = registry.hull_primary("DELTA")
    assert rule == "hull_form_name"
    assert entry.asset_id == "geometry/hulls/DELTA/delta_hull.igs"
    assert registry.resolve("DELTA").endswith("delta_hull.igs")


def test_unique_geometry_rule_covers_hulls_with_one_file(store: AssetStore) -> None:
    registry = build_registry(store=store)
    entry, rule = registry.hull_primary("GAMMA")
    assert rule == "unique_geometry"
    assert entry.asset_id == "geometry/hulls/GAMMA/gamma_modelScale_250mm.igs"
    # hash-clean but never declared verified: reading is fine, publishing is not
    registry.resolve("GAMMA")
    with pytest.raises(AssetError) as excinfo:
        registry.resolve("GAMMA", require_verified=True)
    assert excinfo.value.reason == "unverified"


def test_ambiguity_is_refused_with_the_candidates_named(store: AssetStore) -> None:
    registry = build_registry(store=store)
    with pytest.raises(RegistryError) as excinfo:
        registry.resolve("OMEGA")
    message = str(excinfo.value)
    assert "omega_a.igs" in message and "omega_b.igs" in message
    assert "no unambiguous hull-form geometry" in message


def test_unknown_hull_is_still_refused_by_name(store: AssetStore) -> None:
    registry = build_registry(store=store)
    with pytest.raises(RegistryError):
        registry.resolve("NOT_A_HULL")


def test_a_spec_binding_outranks_every_naming_rule(tmp_path: Path, store: AssetStore) -> None:
    registry = build_registry(store=store, specs_path=_spec(tmp_path, store))
    result = registry.resolve_explained("ALPHA")
    assert result["via"] == "spec_binding"
    assert result["asset_id"] == "geometry/hulls/DELTA/delta_skeg.igs"


def test_resolve_explained_reports_the_rule_and_the_fingerprint(store: AssetStore) -> None:
    registry = build_registry(store=store)
    result = registry.resolve_explained("DELTA")
    assert result["via"] == "hull_primary" and result["rule"] == "hull_form_name"
    assert result["sha256"] == registry.entry(result["asset_id"]).sha256
    assert Path(result["path"]).read_bytes() == b"DELTA-HULL"
    # an explicit asset id needs no inference at all
    assert registry.resolve_explained("asset://geometry/hulls/DELTA/delta_hull.igs")["via"] == "asset_id"


def test_coverage_counts_only_spec_ships_so_hull_coverage_is_reported(store: AssetStore, tmp_path: Path) -> None:
    """The honest metric: 4 of 5 hulls resolvable while ``coverage`` still says 1.0."""
    registry = build_registry(store=store, specs_path=_spec(tmp_path, store))
    snapshot = registry.snapshot_dict()
    assert snapshot["coverage"] == 1.0  # the one spec ship binds
    assert snapshot["hulls"] == 4
    assert snapshot["hulls_resolvable"] == 3  # OMEGA is ambiguous
    assert snapshot["hull_coverage"] == 0.75


def test_growth_carries_the_hull_coverage_delta(tmp_path: Path, store: AssetStore) -> None:
    registry = build_registry(store=store)
    current = registry.snapshot_dict()
    previous = dict(current, hull_coverage=0.5)
    growth = registry.growth(previous)
    assert growth["delta"]["hull_coverage"] == pytest.approx(current["hull_coverage"] - 0.5)


# --------------------------------------------------------------------------- #
# the committed registry against the durable store
# --------------------------------------------------------------------------- #
@needs_store
def test_committed_registry_reaches_every_registered_hull() -> None:
    registry = load_registry(REGISTRY)
    snapshot = registry.snapshot
    assert snapshot["hull_coverage"] == 1.0, snapshot
    for hull in registry.hulls:
        path = registry.resolve(hull)
        assert Path(path).is_file(), hull
        assert Path(path).is_relative_to(REAL_STORE.root), hull


@needs_store
def test_dtc_resolves_to_the_hull_and_not_an_appendage() -> None:
    """DTC has no fleet-spec binding; the store holds its hull byte plus two appendages."""
    path = geometry_path("DTC")
    assert Path(path).name == "DTC_hull.igs"
    registry = load_registry(REGISTRY)
    assert registry.resolve("DTC", require_verified=True) == path  # manifest says verified
    assert Path(Path(path).parent / "DTC_bilgekeel.igs").is_file()  # the trap that was avoided


@needs_store
def test_dtc_resolution_returns_the_measured_byte() -> None:
    """The name points at the *measured* byte, not merely at a plausible file.

    ``benchmarks/verified/fleet_refs/dtc_hull_measured.json`` is the record that measured
    DTC's displacement; if the resolver returned a different byte, the measurement would
    silently stop describing what a caller gets.
    """
    record = json.loads((REPO / "benchmarks/verified/fleet_refs/dtc_hull_measured.json").read_text())
    measured_sha = record["source"]["hull_sha256"]
    result = load_registry(REGISTRY).resolve_explained("DTC")
    assert result["sha256"] == measured_sha
    assert result["asset_id"] == "geometry/hulls/DTC/DTC_hull.igs"
    displacement = record["measured"]["checks"]["displacement_volume"]
    assert displacement["within_tolerance"]
    assert abs(displacement["rel_error_percent"]) < 3.0  # the platform's hard line


@needs_store
def test_kp505_resolves_by_unique_geometry() -> None:
    result = load_registry(REGISTRY).resolve_explained("KP505")
    assert result["rule"] == "unique_geometry"
    assert result["asset_id"] == "geometry/hulls/KP505/KP505_modelScale_250mm.igs"
    assert result["verified"] is False  # not in the fleet baseline yet, and it says so


@needs_store
def test_cli_resolve_shows_how_the_name_was_tied_to_the_bytes() -> None:
    proc = subprocess.run(
        [sys.executable, "-m", "tensorcad.assetlib", "resolve", "DTC"],
        cwd=REPO,
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, proc.stderr
    payload = json.loads(proc.stdout)
    assert payload["resolved"] is True
    assert payload["via"] == "hull_primary" and payload["rule"] == "hull_form_name"
    assert payload["asset_id"] == "geometry/hulls/DTC/DTC_hull.igs"