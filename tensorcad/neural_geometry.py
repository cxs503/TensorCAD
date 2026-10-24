"""Safe metadata and bundle loading for experimental neural geometry objects.

This module intentionally uses only the Python standard library. It validates
metadata and payload integrity; it does not deserialize executable model objects
or claim that a neural field has been trained.
"""
from __future__ import annotations

import hashlib
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

_FORMAT = "tensorcad.neural-geometry"
_VERSION_RE = re.compile(r"^(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)$")
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_ALLOWED_ARCHITECTURES = frozenset({"tensorcad.mlp_fourier.v1"})


class NeuralGeometryFormatError(ValueError):
    """Raised when an NGOM bundle is malformed or unsafe."""


@dataclass(frozen=True)
class PayloadRef:
    path: str
    sha256: str


@dataclass(frozen=True)
class NeuralGeometryManifest:
    format_version: str
    object_id: str
    revision: int
    architecture_id: str
    architecture_config: dict[str, Any]
    weights: PayloadRef
    field: dict[str, Any]
    coordinates: dict[str, Any]
    state_file: str
    training_file: str
    validation_file: str
    raw: dict[str, Any]


def _reject_constant(value: str) -> None:
    raise NeuralGeometryFormatError(f"non-standard JSON numeric constant is forbidden: {value}")


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"), parse_constant=_reject_constant)
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise NeuralGeometryFormatError(f"cannot read valid JSON from {path.name}: {exc}") from exc
    if not isinstance(value, dict):
        raise NeuralGeometryFormatError(f"{path.name} must contain a JSON object")
    _ensure_finite(value, path.name)
    return value


def _ensure_finite(value: Any, where: str) -> None:
    if isinstance(value, float) and not math.isfinite(value):
        raise NeuralGeometryFormatError(f"non-finite number in {where}")
    if isinstance(value, dict):
        for key, child in value.items():
            _ensure_finite(child, f"{where}.{key}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _ensure_finite(child, f"{where}[{index}]")


def _safe_relative_path(raw: Any, field_name: str) -> str:
    if not isinstance(raw, str) or not raw.strip():
        raise NeuralGeometryFormatError(f"{field_name} must be a non-empty relative path")
    candidate = Path(raw)
    if candidate.is_absolute() or ".." in candidate.parts or candidate.as_posix() != raw:
        raise NeuralGeometryFormatError(f"unsafe bundle path in {field_name}: {raw!r}")
    if raw in {".", ""}:
        raise NeuralGeometryFormatError(f"invalid bundle path in {field_name}")
    return raw


def _required_mapping(parent: dict[str, Any], key: str) -> dict[str, Any]:
    value = parent.get(key)
    if not isinstance(value, dict):
        raise NeuralGeometryFormatError(f"manifest field {key!r} must be an object")
    return value


def parse_manifest(data: dict[str, Any]) -> NeuralGeometryManifest:
    """Validate a manifest's core contract without opening any payload."""
    _ensure_finite(data, "manifest")
    if data.get("format") != _FORMAT:
        raise NeuralGeometryFormatError(f"unsupported format: {data.get('format')!r}")
    version = data.get("format_version")
    if not isinstance(version, str) or not (match := _VERSION_RE.fullmatch(version)):
        raise NeuralGeometryFormatError("format_version must be semantic major.minor.patch")
    if int(match.group(1)) != 0:
        raise NeuralGeometryFormatError(f"unsupported NGOM major version: {version}")
    object_id = data.get("object_id")
    if not isinstance(object_id, str) or not object_id.strip():
        raise NeuralGeometryFormatError("object_id must be a non-empty string")
    revision = data.get("revision")
    if not isinstance(revision, int) or isinstance(revision, bool) or revision < 1:
        raise NeuralGeometryFormatError("revision must be a positive integer")

    model = _required_mapping(data, "model")
    if model.get("family") != "conditional_implicit_field":
        raise NeuralGeometryFormatError("model.family must be 'conditional_implicit_field'")
    architecture_id = model.get("architecture_id")
    if architecture_id not in _ALLOWED_ARCHITECTURES:
        raise NeuralGeometryFormatError(f"architecture is not allow-listed: {architecture_id!r}")
    config = _required_mapping(model, "architecture_config")
    if config.get("coordinate_dim") != 3 or config.get("output_dim") != 1:
        raise NeuralGeometryFormatError("architecture must declare coordinate_dim=3 and output_dim=1")
    for key in ("hidden_width", "hidden_layers"):
        value = config.get(key)
        if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
            raise NeuralGeometryFormatError(f"architecture_config.{key} must be a positive integer")
    state_dim = config.get("state_dim")
    if not isinstance(state_dim, int) or isinstance(state_dim, bool) or state_dim < 0:
        raise NeuralGeometryFormatError("architecture_config.state_dim must be a non-negative integer")

    weights_data = _required_mapping(model, "weights")
    weights_path = _safe_relative_path(weights_data.get("path"), "model.weights.path")
    digest = weights_data.get("sha256")
    if not isinstance(digest, str) or not _SHA256_RE.fullmatch(digest):
        raise NeuralGeometryFormatError("model.weights.sha256 must be a lowercase SHA-256 hex digest")

    field = _required_mapping(data, "field")
    if field.get("quantity") != "truncated_signed_distance":
        raise NeuralGeometryFormatError("field.quantity must be truncated_signed_distance")
    if field.get("inside_sign") != "negative":
        raise NeuralGeometryFormatError("field.inside_sign must be negative")
    if field.get("surface_isovalue") != 0.0:
        raise NeuralGeometryFormatError("field.surface_isovalue must be 0")
    if not isinstance(field.get("length_unit"), str) or not field["length_unit"]:
        raise NeuralGeometryFormatError("field.length_unit is required")
    truncation = field.get("truncation_distance")
    if not isinstance(truncation, (int, float)) or isinstance(truncation, bool) or truncation <= 0:
        raise NeuralGeometryFormatError("field.truncation_distance must be positive")

    coordinates = _required_mapping(data, "coordinates")
    if coordinates.get("dimension") != 3:
        raise NeuralGeometryFormatError("coordinates.dimension must be 3")
    if not isinstance(coordinates.get("frame"), str) or not coordinates["frame"]:
        raise NeuralGeometryFormatError("coordinates.frame is required")
    paths = {
        key: _safe_relative_path(data.get(key), key)
        for key in ("state_file", "training_file", "validation_file")
    }
    return NeuralGeometryManifest(
        format_version=version,
        object_id=object_id,
        revision=revision,
        architecture_id=architecture_id,
        architecture_config=config,
        weights=PayloadRef(weights_path, digest),
        field=field,
        coordinates=coordinates,
        state_file=paths["state_file"],
        training_file=paths["training_file"],
        validation_file=paths["validation_file"],
        raw=data,
    )


def load_bundle(bundle_path: str | Path) -> NeuralGeometryManifest:
    """Load and verify an NGOM directory bundle without deserializing model code."""
    root = Path(bundle_path).resolve()
    if not root.is_dir():
        raise NeuralGeometryFormatError(f"bundle directory does not exist: {root}")
    manifest = parse_manifest(_read_json(root / "manifest.json"))
    payload = (root / manifest.weights.path).resolve()
    try:
        payload.relative_to(root)
    except ValueError as exc:
        raise NeuralGeometryFormatError("weights path escapes bundle root") from exc
    if not payload.is_file():
        raise NeuralGeometryFormatError(f"weights payload is missing: {manifest.weights.path}")
    actual_hash = hashlib.sha256(payload.read_bytes()).hexdigest()
    if actual_hash != manifest.weights.sha256:
        raise NeuralGeometryFormatError("weights SHA-256 mismatch")

    # Require all declared metadata files to exist and be valid JSON. Payload
    # formats remain architecture-specific and are never executed by this loader.
    for relative in (manifest.state_file, manifest.training_file, manifest.validation_file):
        file_path = (root / relative).resolve()
        try:
            file_path.relative_to(root)
        except ValueError as exc:
            raise NeuralGeometryFormatError(f"metadata path escapes bundle root: {relative}") from exc
        _read_json(file_path)
    return manifest
