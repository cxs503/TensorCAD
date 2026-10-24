"""Safe metadata, versioned bundle I/O and edit transactions for neural geometry.

A *neural geometry object* (NGO) is a directory bundle whose canonical content is
a learned continuous field plus the metadata needed to interpret it:

    part.ngo/
      manifest.json      # entry point: architecture, field semantics, lineage
      weights.safetensors
      state.json         # named editable controls and their current values
      training.json      # how this revision was produced
      validation.json    # measured acceptance metrics and tolerances
      history.json       # every commit/rollback attempt, oldest first
      derived/           # optional meshes/STL: outputs, never canonical truth

This module is deliberately standard-library only. It validates metadata, hashes
payloads and publishes revisions atomically; it never deserializes executable
model objects. The PyTorch field runtime lives in :mod:`tensorcad.neural`, which
is an optional extra — importing this module must not pull it in.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import re
import shutil
import tempfile
import uuid
from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from dataclasses import field as dataclass_field
from pathlib import Path
from typing import Any

_FORMAT = "tensorcad.neural-geometry"
#: Format version written by this release. Readers accept any ``0.x.y``.
CURRENT_FORMAT_VERSION = "0.2.0"
_FORMAT_VERSION_RE = re.compile(r"^(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)$")
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_ALLOWED_ARCHITECTURES = frozenset({"tensorcad.mlp_fourier.v1"})

MANIFEST_FILE = "manifest.json"
HISTORY_FILE = "history.json"
STATE_FILE = "state.json"
TRAINING_FILE = "training.json"
VALIDATION_FILE = "validation.json"
WEIGHTS_FILE = "weights.safetensors"
DERIVED_DIR = "derived"

COMMITTED = "committed"
ROLLED_BACK = "rolled_back"


class NeuralGeometryFormatError(ValueError):
    """Raised when an NGO bundle is malformed or unsafe."""


class TransactionError(RuntimeError):
    """Raised when a revision cannot be published safely."""


@dataclass(frozen=True)
class PayloadRef:
    path: str
    sha256: str
    bytes: int | None = None

    def to_dict(self) -> dict[str, Any]:
        data: dict[str, Any] = {"path": self.path, "sha256": self.sha256}
        if self.bytes is not None:
            data["bytes"] = self.bytes
        return data


@dataclass(frozen=True)
class Lineage:
    """Where this revision came from and what operation produced it."""

    parent_revision: int | None = None
    root_revision: int | None = None
    operation: str = "fit"
    history_file: str = HISTORY_FILE

    def to_dict(self) -> dict[str, Any]:
        return {
            "parent_revision": self.parent_revision,
            "root_revision": self.root_revision,
            "operation": self.operation,
            "history_file": self.history_file,
        }


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
    lineage: Lineage = dataclass_field(default_factory=Lineage)
    derived: tuple[dict[str, Any], ...] = ()
    validation_status: str = "unknown"
    provenance: dict[str, Any] = dataclass_field(default_factory=dict)

    @property
    def units(self) -> str:
        value = self.field.get("length_unit", "m")
        return str(value)

    @property
    def coordinate_frame(self) -> str:
        return str(self.coordinates.get("frame", "object_local"))

    def summary(self) -> dict[str, Any]:
        """Compact, human-facing description used by ``tensorcad-neural info``."""
        return {
            "object_id": self.object_id,
            "revision": self.revision,
            "format_version": self.format_version,
            "architecture": self.architecture_id,
            "parent_revision": self.lineage.parent_revision,
            "operation": self.lineage.operation,
            "validation_status": self.validation_status,
            "units": self.units,
            "coordinate_frame": self.coordinate_frame,
            "weights_sha256": self.weights.sha256,
            "weights_bytes": self.weights.bytes,
            "truncation_distance": self.field.get("truncation_distance"),
            "derived_artifacts": [artifact.get("kind") for artifact in self.derived],
            "provenance": self.provenance,
        }


# --------------------------------------------------------------------------- #
# hashing and JSON helpers
# --------------------------------------------------------------------------- #


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: str | Path) -> str:
    """Hash a file in chunks so a large payload never has to fit in memory."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_safe(value: Any) -> Any:
    """Replace non-finite floats with ``None`` so the output is standard JSON.

    A failed edit can legitimately produce ``inf``/``nan`` (no valid ray, an
    empty waterplane, ...). Such a run must still be able to *report* why it
    failed, so the values are surfaced as ``null`` rather than crashing the
    serializer or emitting a non-standard ``NaN`` token.
    """
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return value


def dumps_canonical(data: Mapping[str, Any]) -> str:
    """Serialise deterministically: sorted keys, no NaN, trailing newline."""
    return json.dumps(_json_safe(data), indent=2, sort_keys=True, allow_nan=False) + "\n"


def _fsync_dir(path: Path) -> None:
    """Make a directory entry durable (rename visibility after a crash)."""
    try:
        handle = os.open(path, os.O_RDONLY)
    except OSError:  # pragma: no cover - platforms without directory fds
        return
    try:
        os.fsync(handle)
    finally:
        os.close(handle)


def atomic_write_text(path: str | Path, text: str) -> Path:
    """Write text via a sibling temporary file plus ``os.replace``."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(f".{target.name}.{uuid.uuid4().hex}.tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        stream.write(text)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, target)
    _fsync_dir(target.parent)
    return target


def write_json_atomic(path: str | Path, data: Mapping[str, Any]) -> Path:
    """Write a JSON object deterministically and atomically."""
    return atomic_write_text(path, dumps_canonical(data))


def write_json_array_atomic(path: str | Path, data: list[Any]) -> Path:
    """Write a JSON array deterministically and atomically (used for history)."""
    return atomic_write_text(path, json.dumps(_json_safe(data), indent=2, allow_nan=False) + "\n")


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


# --------------------------------------------------------------------------- #
# manifest parsing
# --------------------------------------------------------------------------- #


def parse_manifest(data: dict[str, Any]) -> NeuralGeometryManifest:
    """Validate a manifest's core contract without opening any payload.

    Unknown keys are preserved in :attr:`NeuralGeometryManifest.raw` rather than
    rejected: forward-compatible readers must be able to open objects written by
    a newer minor version without silently dropping metadata they do not model.
    """
    _ensure_finite(data, "manifest")
    if data.get("format") != _FORMAT:
        raise NeuralGeometryFormatError(f"unsupported format: {data.get('format')!r}")
    version = data.get("format_version")
    if not isinstance(version, str) or not (match := _FORMAT_VERSION_RE.fullmatch(version)):
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
    weights_bytes = weights_data.get("bytes")
    if weights_bytes is not None and (
        not isinstance(weights_bytes, int) or isinstance(weights_bytes, bool) or weights_bytes < 0
    ):
        raise NeuralGeometryFormatError("model.weights.bytes must be a non-negative integer")

    field_spec = _required_mapping(data, "field")
    if field_spec.get("quantity") != "truncated_signed_distance":
        raise NeuralGeometryFormatError("field.quantity must be truncated_signed_distance")
    if field_spec.get("inside_sign") != "negative":
        raise NeuralGeometryFormatError("field.inside_sign must be negative")
    if field_spec.get("surface_isovalue") != 0.0:
        raise NeuralGeometryFormatError("field.surface_isovalue must be 0")
    if not isinstance(field_spec.get("length_unit"), str) or not field_spec["length_unit"]:
        raise NeuralGeometryFormatError("field.length_unit is required")
    truncation = field_spec.get("truncation_distance")
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

    lineage_data = data.get("lineage")
    lineage = Lineage(
        parent_revision=(
            lineage_data.get("parent_revision") if isinstance(lineage_data, dict) else None
        ),
        root_revision=(lineage_data.get("root_revision") if isinstance(lineage_data, dict) else None),
        operation=(lineage_data.get("operation", "fit") if isinstance(lineage_data, dict) else "fit"),
        history_file=(
            _safe_relative_path(lineage_data.get("history_file", HISTORY_FILE), "lineage.history_file")
            if isinstance(lineage_data, dict)
            else HISTORY_FILE
        ),
    )
    derived_data = data.get("derived", [])
    if not isinstance(derived_data, list):
        raise NeuralGeometryFormatError("derived must be a list of artifact records")
    derived = tuple(record for record in derived_data if isinstance(record, dict))

    validation_data = data.get("validation")
    validation_status = "unknown"
    if isinstance(validation_data, dict):
        validation_status = str(validation_data.get("status", "unknown"))
    elif VALIDATION_FILE in paths:
        validation_status = "recorded"

    return NeuralGeometryManifest(
        format_version=version,
        object_id=object_id,
        revision=revision,
        architecture_id=architecture_id,
        architecture_config=config,
        weights=PayloadRef(weights_path, digest, weights_bytes),
        field=field_spec,
        coordinates=coordinates,
        state_file=paths["state_file"],
        training_file=paths["training_file"],
        validation_file=paths["validation_file"],
        raw=data,
        lineage=lineage,
        derived=derived,
        validation_status=validation_status,
        provenance=dict(data.get("provenance", {})) if isinstance(data.get("provenance"), dict) else {},
    )


def _verify_payload(root: Path, manifest: NeuralGeometryManifest) -> None:
    payload = (root / manifest.weights.path).resolve()
    try:
        payload.relative_to(root)
    except ValueError as exc:
        raise NeuralGeometryFormatError("weights path escapes bundle root") from exc
    if not payload.is_file():
        raise NeuralGeometryFormatError(f"weights payload is missing: {manifest.weights.path}")
    actual_hash = sha256_file(payload)
    if actual_hash != manifest.weights.sha256:
        raise NeuralGeometryFormatError("weights SHA-256 mismatch")
    if manifest.weights.bytes is not None and payload.stat().st_size != manifest.weights.bytes:
        raise NeuralGeometryFormatError("weights size does not match the manifest")


def read_manifest(bundle_path: str | Path, *, verify_payload: bool = True) -> NeuralGeometryManifest:
    """Parse a bundle's manifest.

    ``verify_payload=False`` is for tools that want to *report* an integrity
    failure (``tensorcad-neural verify``) rather than turn it into an exception
    before they can describe it.
    """
    root = Path(bundle_path).resolve()
    if not root.is_dir():
        raise NeuralGeometryFormatError(f"bundle directory does not exist: {root}")
    manifest = parse_manifest(_read_json(root / MANIFEST_FILE))
    if verify_payload:
        _verify_payload(root, manifest)
    return manifest


def load_bundle(bundle_path: str | Path) -> NeuralGeometryManifest:
    """Load and verify an NGO directory bundle without deserializing model code."""
    root = Path(bundle_path).resolve()
    manifest = read_manifest(root, verify_payload=False)
    _verify_payload(root, manifest)

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


# --------------------------------------------------------------------------- #
# manifest construction
# --------------------------------------------------------------------------- #


def build_manifest(
    *,
    object_id: str,
    revision: int,
    weights: PayloadRef,
    architecture_config: dict[str, Any],
    field_spec: Mapping[str, Any],
    coordinates: Mapping[str, Any],
    state: Mapping[str, Any],
    training: Mapping[str, Any],
    validation: Mapping[str, Any],
    lineage: Lineage,
    provenance: Mapping[str, Any],
    derived: list[Mapping[str, Any]] | None = None,
) -> dict[str, Any]:
    """Assemble a v0.2 manifest.

    Top-level ``state_file``/``training_file``/``validation_file`` are kept for
    readers of the 0.1 layout, while the richer nested records carry the actual
    values, so both an old and a new reader agree about what the object is.
    """
    return {
        "format": _FORMAT,
        "format_version": CURRENT_FORMAT_VERSION,
        "object_id": object_id,
        "revision": revision,
        "lineage": lineage.to_dict(),
        "model": {
            "family": "conditional_implicit_field",
            "architecture_id": "tensorcad.mlp_fourier.v1",
            "architecture_config": dict(architecture_config),
            "weights": weights.to_dict(),
        },
        "field": dict(field_spec),
        "coordinates": dict(coordinates),
        "state": dict(state),
        "training": dict(training),
        "validation": dict(validation),
        "provenance": dict(provenance),
        "derived": [dict(record) for record in (derived or [])],
        "state_file": STATE_FILE,
        "training_file": TRAINING_FILE,
        "validation_file": VALIDATION_FILE,
    }


def derived_record(
    *, kind: str, fmt: str, path: str, digest: str, size_bytes: int, metrics: Mapping[str, Any]
) -> dict[str, Any]:
    """Describe one non-canonical artifact (mesh, STL, image) derived from the field."""
    return {
        "kind": kind,
        "format": fmt,
        "path": _safe_relative_path(path, "derived.path"),
        "sha256": digest,
        "bytes": size_bytes,
        "metrics": dict(metrics),
    }


# --------------------------------------------------------------------------- #
# versioned object store: stage -> validate -> commit / rollback
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class TransactionRecord:
    revision: int
    parent_revision: int | None
    status: str
    operation: str
    weights_sha256: str | None
    metrics: dict[str, Any] = dataclass_field(default_factory=dict)
    failed_checks: tuple[str, ...] = ()
    reason: str | None = None
    recorded_utc: str = ""
    command: list[str] | None = None
    elapsed_seconds: float | None = None

    def to_dict(self) -> dict[str, Any]:
        data: dict[str, Any] = {
            "revision": self.revision,
            "parent_revision": self.parent_revision,
            "status": self.status,
            "operation": self.operation,
            "weights_sha256": self.weights_sha256,
            "metrics": _json_safe(self.metrics),
            "recorded_utc": self.recorded_utc,
        }
        if self.failed_checks:
            data["failed_checks"] = list(self.failed_checks)
        if self.reason:
            data["reason"] = self.reason
        if self.command:
            data["command"] = list(self.command)
        if self.elapsed_seconds is not None:
            data["elapsed_seconds"] = self.elapsed_seconds
        return data


def _utc_now() -> str:
    from datetime import datetime, timezone

    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


class ObjectStore:
    """A directory bundle with append-only history and atomic revision swaps.

    Invariants this class is responsible for:

    * an unfinished edit never becomes visible — candidates are written to a
      sibling staging directory and only ``commit`` renames them into place;
    * a failed validation leaves the active revision byte-identical;
    * every attempt (commit *and* rollback) is appended to ``history.json``;
    * ``expected_parent`` makes a stale writer fail instead of clobbering a
      revision that someone else published in the meantime.
    """

    def __init__(self, root: str | Path) -> None:
        self.root = Path(root)

    # -- reading ---------------------------------------------------------- #

    @property
    def exists(self) -> bool:
        return (self.root / MANIFEST_FILE).is_file()

    def manifest(self) -> NeuralGeometryManifest:
        return load_bundle(self.root)

    def active_revision(self) -> int | None:
        """Revision currently on disk, or ``None`` for an object that is not published yet."""
        if not self.exists:
            return None
        return parse_manifest(_read_json(self.root / MANIFEST_FILE)).revision

    def active_weights_sha256(self) -> str | None:
        if not self.exists:
            return None
        weights = self.root / WEIGHTS_FILE
        return sha256_file(weights) if weights.is_file() else None

    def history(self) -> list[dict[str, Any]]:
        path = self.root / HISTORY_FILE
        if not path.is_file():
            return []
        entries = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(entries, list):
            raise NeuralGeometryFormatError("history.json must contain a JSON array")
        return [entry for entry in entries if isinstance(entry, dict)]

    def state(self) -> dict[str, Any]:
        data = _read_json(self.root / STATE_FILE)
        return data

    def validation(self) -> dict[str, Any]:
        return _read_json(self.root / VALIDATION_FILE)

    # -- staging ---------------------------------------------------------- #

    def stage(self, files: Mapping[str, bytes | str], manifest: Mapping[str, Any]) -> Path:
        """Write a candidate bundle to a hidden staging directory.

        Nothing here is visible to readers: the staging directory sits next to
        the bundle root under a dot-prefixed name and is renamed into place by
        :meth:`commit` only.
        """
        parsed = parse_manifest(dict(manifest))
        self.root.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(
            tempfile.mkdtemp(prefix=f".{self.root.name}.staging-", dir=self.root.parent)
        )
        try:
            for relative, payload in files.items():
                target = staging / _safe_relative_path(relative, "staged file")
                target.parent.mkdir(parents=True, exist_ok=True)
                if isinstance(payload, str):
                    target.write_text(payload, encoding="utf-8")
                else:
                    target.write_bytes(payload)
            write_json_atomic(staging / MANIFEST_FILE, dict(manifest))
            _verify_payload(staging, parsed)
        except Exception:
            shutil.rmtree(staging, ignore_errors=True)
            raise
        return staging

    # -- publishing ------------------------------------------------------- #

    def commit(
        self,
        staging: Path,
        *,
        expected_parent: int | None,
        record: TransactionRecord,
    ) -> NeuralGeometryManifest:
        """Atomically make ``staging`` the active revision.

        ``expected_parent`` is optimistic concurrency control: pass the revision
        the candidate was trained from, and the commit fails if the on-disk
        revision moved on in the meantime.
        """
        current = self.active_revision()
        if current != expected_parent:
            shutil.rmtree(staging, ignore_errors=True)
            raise TransactionError(
                f"stale candidate: expected parent revision {expected_parent}, on-disk revision is {current}"
            )
        if record.status != COMMITTED:
            shutil.rmtree(staging, ignore_errors=True)
            raise TransactionError(f"refusing to commit a {record.status!r} transaction")

        previous_history = self.history()
        entries = [*previous_history, record.to_dict()]
        write_json_array_atomic(staging / HISTORY_FILE, entries)
        staged_manifest = parse_manifest(_read_json(staging / MANIFEST_FILE))
        if staged_manifest.revision != record.revision:
            shutil.rmtree(staging, ignore_errors=True)
            raise TransactionError(
                f"staged manifest revision {staged_manifest.revision} does not match record {record.revision}"
            )

        backup: Path | None = None
        if self.root.exists():
            backup = self.root.with_name(f".{self.root.name}.previous-{uuid.uuid4().hex}")
            os.replace(self.root, backup)
        try:
            os.replace(staging, self.root)
            _fsync_dir(self.root.parent)
        except Exception:  # pragma: no cover - filesystem failure path
            if backup is not None and not self.root.exists():
                os.replace(backup, self.root)
            raise
        if backup is not None:
            shutil.rmtree(backup, ignore_errors=True)
        return self.manifest()

    def rollback(
        self,
        staging: Path | None,
        *,
        record: TransactionRecord,
    ) -> TransactionRecord:
        """Discard a candidate and append the failed attempt to history.

        The active revision is never touched: only ``history.json`` is rewritten,
        and that rewrite is itself atomic.
        """
        if staging is not None:
            shutil.rmtree(staging, ignore_errors=True)
        if record.status != ROLLED_BACK:
            raise TransactionError(f"refusing to record a {record.status!r} transaction as rolled back")
        if not self.exists:
            raise TransactionError("cannot roll back: no active revision exists")
        entries = [*self.history(), record.to_dict()]
        write_json_array_atomic(self.root / HISTORY_FILE, entries)
        return record

    def publish(
        self,
        *,
        files: Mapping[str, bytes | str],
        manifest: Mapping[str, Any],
        expected_parent: int | None,
        record: TransactionRecord,
    ) -> NeuralGeometryManifest:
        """``stage`` + ``commit`` in one call, with cleanup on any failure."""
        staging = self.stage(files, manifest)
        try:
            return self.commit(staging, expected_parent=expected_parent, record=record)
        except Exception:
            shutil.rmtree(staging, ignore_errors=True)
            raise

    def derived_dir(self) -> Path:
        path = self.root / DERIVED_DIR
        path.mkdir(parents=True, exist_ok=True)
        return path

    def iter_files(self) -> Iterator[Path]:
        for path in sorted(self.root.rglob("*")):
            if path.is_file():
                yield path
