"""Tests for safe neural-geometry metadata loading."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from tensorcad.neural_geometry import NeuralGeometryFormatError, load_bundle, parse_manifest


def manifest_for(weights: bytes) -> dict[str, object]:
    return {
        "format": "tensorcad.neural-geometry",
        "format_version": "0.1.0",
        "object_id": "urn:tensorcad:test:sphere",
        "revision": 1,
        "model": {
            "family": "conditional_implicit_field",
            "architecture_id": "tensorcad.mlp_fourier.v1",
            "architecture_config": {
                "coordinate_dim": 3,
                "state_dim": 1,
                "hidden_width": 16,
                "hidden_layers": 2,
                "coordinate_encoding": "fourier",
                "output_dim": 1,
            },
            "weights": {"path": "weights.safetensors", "sha256": hashlib.sha256(weights).hexdigest()},
        },
        "field": {
            "quantity": "truncated_signed_distance",
            "inside_sign": "negative",
            "surface_isovalue": 0.0,
            "length_unit": "m",
            "truncation_distance": 0.1,
        },
        "coordinates": {"dimension": 3, "frame": "object_local"},
        "state_file": "state.json",
        "training_file": "training.json",
        "validation_file": "validation.json",
    }


def write_bundle(path: Path, *, weights: bytes = b"test weight bytes") -> None:
    path.mkdir()
    manifest = manifest_for(weights)
    (path / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    (path / "weights.safetensors").write_bytes(weights)
    (path / "state.json").write_text('{"schema_version":"0.1.0","controls":{}}', encoding="utf-8")
    (path / "training.json").write_text('{"schema_version":"0.1.0","steps":0}', encoding="utf-8")
    (path / "validation.json").write_text('{"schema_version":"0.1.0","status":"not_run"}', encoding="utf-8")


def test_parse_manifest_accepts_valid_core_contract() -> None:
    parsed = parse_manifest(manifest_for(b"weights"))
    assert parsed.object_id == "urn:tensorcad:test:sphere"
    assert parsed.revision == 1
    assert parsed.weights.sha256 == hashlib.sha256(b"weights").hexdigest()


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (lambda d: d.update(format="unknown"), "unsupported format"),
        (lambda d: d.update(format_version="1.0.0"), "unsupported NGOM major"),
        (lambda d: d.update(revision=0), "positive integer"),
        (lambda d: d["model"].update(architecture_id="python.module:run"), "not allow-listed"),
        (lambda d: d["model"]["weights"].update(path="../outside.safetensors"), "unsafe bundle path"),
        (lambda d: d["model"]["weights"].update(sha256="bad"), "SHA-256"),
        (lambda d: d["field"].update(inside_sign="positive"), "inside_sign"),
        (lambda d: d["field"].update(truncation_distance=float("nan")), "non-finite"),
    ],
)
def test_rejects_malformed_manifest(mutate: object, message: str) -> None:
    data = manifest_for(b"weights")
    mutate(data)  # type: ignore[operator]
    with pytest.raises(NeuralGeometryFormatError, match=message):
        parse_manifest(data)


def test_bundle_loader_verifies_hash_and_metadata(tmp_path: Path) -> None:
    bundle = tmp_path / "sphere.ngo"
    write_bundle(bundle)
    loaded = load_bundle(bundle)
    assert loaded.object_id == "urn:tensorcad:test:sphere"
    assert loaded.weights.path == "weights.safetensors"


def test_bundle_loader_rejects_weight_hash_mismatch(tmp_path: Path) -> None:
    bundle = tmp_path / "sphere.ngo"
    write_bundle(bundle)
    (bundle / "weights.safetensors").write_bytes(b"tampered")
    with pytest.raises(NeuralGeometryFormatError, match="SHA-256 mismatch"):
        load_bundle(bundle)


def test_bundle_loader_rejects_missing_metadata(tmp_path: Path) -> None:
    bundle = tmp_path / "sphere.ngo"
    write_bundle(bundle)
    (bundle / "state.json").unlink()
    with pytest.raises(NeuralGeometryFormatError, match="cannot read valid JSON"):
        load_bundle(bundle)


def test_bundle_loader_rejects_non_object_json(tmp_path: Path) -> None:
    bundle = tmp_path / "sphere.ngo"
    write_bundle(bundle)
    (bundle / "state.json").write_text("[]", encoding="utf-8")
    with pytest.raises(NeuralGeometryFormatError, match="must contain a JSON object"):
        load_bundle(bundle)
