"""Transaction and metadata tests for the neural geometry object store (no torch)."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from tensorcad.neural_geometry import (
    COMMITTED,
    ROLLED_BACK,
    Lineage,
    NeuralGeometryFormatError,
    ObjectStore,
    PayloadRef,
    TransactionError,
    TransactionRecord,
    build_manifest,
    dumps_canonical,
    sha256_bytes,
    write_json_atomic,
)

PAYLOAD = b"weights-v1"


def make_manifest(
    *,
    revision: int,
    payload: bytes = PAYLOAD,
    parent: int | None = None,
    operation: str = "fit",
) -> dict[str, object]:
    return build_manifest(
        object_id="urn:tensorcad:test:sphere",
        revision=revision,
        weights=PayloadRef("weights.safetensors", sha256_bytes(payload), len(payload)),
        architecture_config={"coordinate_dim": 3, "output_dim": 1, "state_dim": 1,
                             "hidden_width": 8, "hidden_layers": 2},
        field_spec={"quantity": "truncated_signed_distance", "inside_sign": "negative",
                    "surface_isovalue": 0.0, "length_unit": "m", "truncation_distance": 0.15},
        coordinates={"dimension": 3, "frame": "object_local",
                     "bounds": {"min": [-0.6, -0.6, -0.6], "max": [0.6, 0.6, 0.6]}},
        state={"schema_version": "0.2.0", "field_state": {"radius": 0.5}, "controls": {}},
        training={"schema_version": "0.2.0", "steps": 1},
        validation={"schema_version": "0.2.0", "status": "passed", "revision": revision},
        lineage=Lineage(parent_revision=parent, root_revision=1, operation=operation),
        provenance={"producer": "test"},
    )


def bundle_files(payload: bytes = PAYLOAD) -> dict[str, bytes | str]:
    return {
        "weights.safetensors": payload,
        "state.json": dumps_canonical({"schema_version": "0.2.0", "field_state": {"radius": 0.5}}),
        "training.json": dumps_canonical({"schema_version": "0.2.0", "steps": 1}),
        "validation.json": dumps_canonical({"schema_version": "0.2.0", "status": "passed"}),
    }


def record(revision: int, parent: int | None, status: str = COMMITTED, **extra: object) -> TransactionRecord:
    return TransactionRecord(
        revision=revision,
        parent_revision=parent,
        status=status,
        operation="fit" if parent is None else "dent_edit",
        weights_sha256=sha256_bytes(PAYLOAD) if status == COMMITTED else None,
        recorded_utc="2026-10-24T00:00:00+00:00",
        **extra,  # type: ignore[arg-type]
    )


def publish_first(store: ObjectStore) -> None:
    store.publish(
        files=bundle_files(),
        manifest=make_manifest(revision=1),
        expected_parent=None,
        record=record(1, None),
    )


def test_commit_publishes_revision_and_writes_history(tmp_path: Path) -> None:
    store = ObjectStore(tmp_path / "sphere.ngo")
    assert store.active_revision() is None
    publish_first(store)
    assert store.active_revision() == 1
    assert store.exists
    entries = store.history()
    assert [entry["status"] for entry in entries] == [COMMITTED]
    assert entries[0]["revision"] == 1
    # no staging or backup directories left behind
    leftovers = [path.name for path in tmp_path.iterdir() if path.name.startswith(".")]
    assert leftovers == []


def test_second_commit_keeps_previous_history(tmp_path: Path) -> None:
    store = ObjectStore(tmp_path / "sphere.ngo")
    publish_first(store)
    payload = b"weights-v2"
    store.publish(
        files=bundle_files(payload),
        manifest=make_manifest(revision=2, payload=payload, parent=1, operation="dent_edit"),
        expected_parent=1,
        record=record(2, 1),
    )
    assert store.active_revision() == 2
    assert [entry["revision"] for entry in store.history()] == [1, 2]
    assert store.manifest().lineage.parent_revision == 1


def test_stale_parent_is_rejected_and_nothing_changes(tmp_path: Path) -> None:
    store = ObjectStore(tmp_path / "sphere.ngo")
    publish_first(store)
    before = store.active_weights_sha256()
    payload = b"weights-v2"
    staging = store.stage(bundle_files(payload), make_manifest(revision=2, payload=payload, parent=1))
    with pytest.raises(TransactionError, match="stale candidate"):
        store.commit(staging, expected_parent=0, record=record(2, 1))
    assert not staging.exists(), "the rejected staging directory must be removed"
    assert store.active_revision() == 1
    assert store.active_weights_sha256() == before


def test_rollback_keeps_active_revision_byte_identical(tmp_path: Path) -> None:
    store = ObjectStore(tmp_path / "sphere.ngo")
    publish_first(store)
    before = (store.root / "weights.safetensors").read_bytes()
    store.rollback(
        None,
        record=record(2, 1, status=ROLLED_BACK, metrics={"dent_peak_error_m": float("inf")},
                      failed_checks=("dent_peak_error",), reason="gates not met"),
    )
    assert (store.root / "weights.safetensors").read_bytes() == before
    assert store.active_revision() == 1
    entries = store.history()
    assert entries[-1]["status"] == ROLLED_BACK
    assert entries[-1]["failed_checks"] == ["dent_peak_error"]
    assert entries[-1]["metrics"]["dent_peak_error_m"] is None, "non-finite metrics must become JSON null"


def test_commit_refuses_non_committed_record(tmp_path: Path) -> None:
    store = ObjectStore(tmp_path / "sphere.ngo")
    staging = store.stage(bundle_files(), make_manifest(revision=1))
    with pytest.raises(TransactionError, match="refusing to commit"):
        store.commit(staging, expected_parent=None, record=record(1, None, status=ROLLED_BACK))
    assert not store.exists


def test_rollback_requires_an_active_revision(tmp_path: Path) -> None:
    store = ObjectStore(tmp_path / "missing.ngo")
    with pytest.raises(TransactionError, match="no active revision"):
        store.rollback(None, record=record(1, None, status=ROLLED_BACK))


def test_revision_mismatch_between_manifest_and_record_is_rejected(tmp_path: Path) -> None:
    store = ObjectStore(tmp_path / "sphere.ngo")
    staging = store.stage(bundle_files(), make_manifest(revision=3))
    with pytest.raises(TransactionError, match="does not match record"):
        store.commit(staging, expected_parent=None, record=record(1, None))


def test_stage_rejects_payload_that_does_not_match_the_manifest(tmp_path: Path) -> None:
    store = ObjectStore(tmp_path / "sphere.ngo")
    manifest = make_manifest(revision=1)
    files = bundle_files(b"tampered")
    with pytest.raises(NeuralGeometryFormatError, match="SHA-256 mismatch"):
        store.stage(files, manifest)
    leftovers = [path.name for path in tmp_path.iterdir() if path.name.startswith(".")]
    assert leftovers == []


def test_publish_surfaces_commit_failures_without_leaving_state(tmp_path: Path) -> None:
    store = ObjectStore(tmp_path / "sphere.ngo")
    publish_first(store)
    with pytest.raises(TransactionError):
        store.publish(
            files=bundle_files(b"weights-v2"),
            manifest=make_manifest(revision=2, payload=b"weights-v2", parent=1),
            expected_parent=7,
            record=record(2, 1),
        )
    assert store.active_revision() == 1
    leftovers = [path.name for path in tmp_path.iterdir() if path.name.startswith(".")]
    assert leftovers == []


def test_json_writes_are_deterministic_and_reject_non_finite(tmp_path: Path) -> None:
    first = tmp_path / "a.json"
    second = tmp_path / "b.json"
    write_json_atomic(first, {"b": 1, "a": float("nan")})
    write_json_atomic(second, {"a": float("nan"), "b": 1})
    assert first.read_text() == second.read_text()
    assert json.loads(first.read_text())["a"] is None


def test_derived_records_must_stay_inside_the_bundle() -> None:
    from tensorcad.neural_geometry import derived_record

    with pytest.raises(NeuralGeometryFormatError, match="unsafe bundle path"):
        derived_record(kind="mesh", fmt="obj", path="../escape.obj", digest="0" * 64,
                       size_bytes=1, metrics={})


def test_manifest_summary_and_payload_reference(tmp_path: Path) -> None:
    store = ObjectStore(tmp_path / "sphere.ngo")
    publish_first(store)
    manifest = store.manifest()
    summary = manifest.summary()
    assert summary["revision"] == 1
    assert summary["weights_bytes"] == len(PAYLOAD)
    assert summary["parent_revision"] is None
    assert summary["derived_artifacts"] == []
    assert PayloadRef("w", "0" * 64).to_dict() == {"path": "w", "sha256": "0" * 64}


def test_store_exposes_state_validation_and_file_listing(tmp_path: Path) -> None:
    store = ObjectStore(tmp_path / "sphere.ngo")
    publish_first(store)
    assert store.state()["field_state"]["radius"] == 0.5
    assert store.validation()["status"] == "passed"
    files = {path.name for path in store.iter_files()}
    assert {"manifest.json", "weights.safetensors", "state.json", "history.json"} <= files
    assert ObjectStore(tmp_path / "absent.ngo").active_weights_sha256() is None


def test_read_manifest_can_report_a_tampered_payload(tmp_path: Path) -> None:
    from tensorcad.neural_geometry import read_manifest

    store = ObjectStore(tmp_path / "sphere.ngo")
    publish_first(store)
    (store.root / "weights.safetensors").write_bytes(b"tampered")
    with pytest.raises(NeuralGeometryFormatError, match="SHA-256 mismatch"):
        read_manifest(store.root)
    manifest = read_manifest(store.root, verify_payload=False)
    assert manifest.revision == 1


def test_derived_artifacts_survive_a_round_trip(tmp_path: Path) -> None:
    from tensorcad.neural_geometry import derived_record

    store = ObjectStore(tmp_path / "sphere.ngo")
    record_entry = derived_record(
        kind="surface_mesh", fmt="obj", path="derived/surface.obj",
        digest="1" * 64, size_bytes=128, metrics={"watertight": True},
    )
    manifest = make_manifest(revision=1)
    manifest["derived"] = [record_entry]
    store.publish(
        files={**bundle_files(), "derived/surface.obj": b"o mesh"},
        manifest=manifest,
        expected_parent=None,
        record=record(1, None),
    )
    loaded = store.manifest()
    assert loaded.derived[0]["format"] == "obj"
    assert loaded.derived[0]["metrics"]["watertight"] is True
    assert (store.root / "derived/surface.obj").read_bytes() == b"o mesh"


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (lambda d: d.update(format_version="0.2"), "semantic major.minor.patch"),
        (lambda d: d["model"]["weights"].update(bytes=-1), "non-negative integer"),
        (lambda d: d.update(derived="surface.obj"), "list of artifact records"),
        (lambda d: d["lineage"].update(history_file="../history.json"), "unsafe bundle path"),
        (lambda d: d.update(state_file="/etc/passwd"), "unsafe bundle path"),
    ],
)
def test_manifest_rejects_malformed_v2_fields(mutate: object, message: str) -> None:
    from tensorcad.neural_geometry import parse_manifest

    data = make_manifest(revision=1)
    mutate(data)  # type: ignore[operator]
    with pytest.raises(NeuralGeometryFormatError, match=message):
        parse_manifest(data)


def test_manifest_parses_a_minimal_v01_document_for_backwards_compatibility() -> None:
    from tensorcad.neural_geometry import parse_manifest

    legacy = make_manifest(revision=1)
    for key in ("lineage", "derived", "state", "training", "validation", "provenance"):
        legacy.pop(key, None)
    parsed = parse_manifest(legacy)
    assert parsed.revision == 1
    assert parsed.lineage.history_file == "history.json"
    assert parsed.derived == ()
