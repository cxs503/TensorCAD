"""Object-level operations: fit, edit, validate and extract.

This is the layer an application calls. It deliberately keeps the three concerns
separate:

* :mod:`tensorcad.neural.field` knows about networks and measurements;
* :mod:`tensorcad.neural_geometry` knows about bundles, hashing and transactions;
* this module turns one into a committed revision of the other, and refuses to
  publish anything that has not passed its measured gates.
"""
from __future__ import annotations

import math
import shutil
import sys
import tempfile
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, ClassVar

from ..export import write_obj, write_stl
from ..geometry import BoundingBox, Mesh
from ..geometry_field import (
    TRUNCATED_SIGNED_DISTANCE,
    FieldError,
    GeometryField,
    Point,
    _points,
)
from ..neural_geometry import (
    COMMITTED,
    HISTORY_FILE,
    ROLLED_BACK,
    STATE_FILE,
    TRAINING_FILE,
    VALIDATION_FILE,
    WEIGHTS_FILE,
    Lineage,
    NeuralGeometryFormatError,
    ObjectStore,
    PayloadRef,
    TransactionRecord,
    build_manifest,
    derived_record,
    sha256_bytes,
)
from .extract import ExtractionReport, extract_surface, summarise
from .field import (
    ARCHITECTURE_ID,
    DEFAULT_TRUNCATION,
    ConditionalTSDF,
    FieldConfig,
    field_callable,
    fine_tune_dent,
    gradient_callable,
    load_weights,
    measure_dent_edit,
    measure_sphere,
    normalise_direction,
    save_weights,
    train_base_field,
)

SCHEMA_VERSION = "0.2.0"


class NeuralObjectError(RuntimeError):
    """Raised when an object cannot be fitted, edited or extracted."""


@dataclass(frozen=True)
class Tolerances:
    """Acceptance gates. These are declared per object and recorded with the result."""

    field_mae: float = 0.008
    boundary_error_mean: float = 0.010
    ray_coverage: float = 0.95
    dent_peak_error: float = 0.025
    dent_region_error_mean: float = 0.020
    protected_mean_displacement: float = 0.012
    protected_p99_displacement: float = 0.020
    protected_max_displacement: float = 0.030

    def to_dict(self) -> dict[str, float]:
        return {
            "field_mae_m": self.field_mae,
            "boundary_error_mean_m": self.boundary_error_mean,
            "ray_coverage": self.ray_coverage,
            "dent_peak_error_m": self.dent_peak_error,
            "dent_region_error_mean_m": self.dent_region_error_mean,
            "protected_mean_displacement_m": self.protected_mean_displacement,
            "protected_p99_displacement_m": self.protected_p99_displacement,
            "protected_max_displacement_m": self.protected_max_displacement,
        }


@dataclass(frozen=True)
class Check:
    name: str
    passed: bool
    measured: float
    limit: float
    comparison: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "passed": self.passed,
            "measured": self.measured,
            "limit": self.limit,
            "comparison": self.comparison,
        }


def _check(name: str, measured: float, limit: float, *, maximum: bool = True) -> Check:
    finite = isinstance(measured, float) and math.isfinite(measured)
    passed = finite and (measured <= limit if maximum else measured >= limit)
    return Check(
        name=name,
        passed=bool(passed),
        measured=float(measured) if isinstance(measured, float) else float("nan"),
        limit=float(limit),
        comparison="<=" if maximum else ">=",
    )


def _failed(checks: Sequence[Check]) -> tuple[str, ...]:
    return tuple(check.name for check in checks if not check.passed)


@dataclass(frozen=True)
class FitResult:
    store: Path
    revision: int
    weights_sha256: str
    metrics: dict[str, Any]
    checks: tuple[Check, ...]
    elapsed_seconds: float
    final_loss: float

    @property
    def passed(self) -> bool:
        return not _failed(self.checks)


@dataclass(frozen=True)
class EditResult:
    status: str
    parent_revision: int | None
    active_revision: int | None
    weights_sha256: str | None
    metrics: dict[str, Any]
    checks: tuple[Check, ...]
    elapsed_seconds: float
    reason: str | None = None
    derived: tuple[dict[str, Any], ...] = field(default=())

    @property
    def failed_checks(self) -> tuple[str, ...]:
        return _failed(self.checks)

    @property
    def committed(self) -> bool:
        return self.status == COMMITTED


def _utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def _command_line(command: Sequence[str] | None) -> list[str]:
    return list(command) if command is not None else list(sys.argv)


def _weights_bytes(model: ConditionalTSDF) -> tuple[bytes, str]:
    """Serialise weights to a temporary file so the digest and blob agree exactly."""
    directory = Path(tempfile.mkdtemp(prefix="tensorcad-weights-"))
    try:
        path = directory / WEIGHTS_FILE
        digest, _size = save_weights(model, path)
        return path.read_bytes(), digest
    finally:
        shutil.rmtree(directory, ignore_errors=True)


def field_state_value(manifest_state: dict[str, Any]) -> float:
    """Read the numeric state the field is conditioned on from ``state.json``."""
    try:
        return float(manifest_state["field_state"]["radius"])
    except (KeyError, TypeError, ValueError) as exc:
        raise NeuralGeometryFormatError(
            "state.field_state.radius is required to evaluate a field"
        ) from exc


def _state_document(radius: float, *, controls_extra: dict[str, Any] | None = None) -> dict[str, Any]:
    controls: dict[str, Any] = {
        "radius": {
            "value": radius,
            "unit": "m",
            "role": "field_state",
            "edit_mode": "conditioned",
            "description": "radius the field is conditioned on; continuous, not a mesh parameter",
        }
    }
    if controls_extra:
        controls.update(controls_extra)
    return {
        "schema_version": SCHEMA_VERSION,
        "field_state": {"radius": radius},
        "controls": controls,
    }


def _coordinates(bounds_min: Sequence[float], bounds_max: Sequence[float]) -> dict[str, Any]:
    return {
        "dimension": 3,
        "frame": "object_local",
        "units": "m",
        "bounds": {
            "min": [float(value) for value in bounds_min],
            "max": [float(value) for value in bounds_max],
        },
    }


def _field_spec(truncation: float) -> dict[str, Any]:
    return {
        "quantity": "truncated_signed_distance",
        "inside_sign": "negative",
        "surface_isovalue": 0.0,
        "length_unit": "m",
        "truncation_distance": truncation,
        "continuity": "continuous",
        "resolution": "implicit; bounded only by the query point density",
    }


def _sphere_bounds(
    radius: float, margin: float = 0.16
) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    half = radius * (1.0 + margin)
    return (-half, -half, -half), (half, half, half)


# --------------------------------------------------------------------------- #
# fit
# --------------------------------------------------------------------------- #


def fit_sphere(
    output: str | Path,
    *,
    radius: float = 0.5,
    steps: int = 1800,
    seed: int = 503,
    tolerance: Tolerances | None = None,
    config: FieldConfig | None = None,
    truncation: float = DEFAULT_TRUNCATION,
    command: Sequence[str] | None = None,
    rays: int = 2048,
    log: Any = None,
) -> FitResult:
    """Fit a conditional sphere family and publish it as revision 1."""
    gates = tolerance or Tolerances()
    store = ObjectStore(output)
    if store.exists:
        raise NeuralObjectError(
            f"{store.root} already holds revision {store.active_revision()}; "
            "fit writes revision 1 into a new object, use edit to continue it"
        )
    if radius <= 0:
        raise ValueError("radius must be positive")

    outcome = train_base_field(
        steps=steps, seed=seed, truncation=truncation, config=config, log_every=50, log=log
    )
    metrics = measure_sphere(
        outcome.model, radius, seed=seed + 1, truncation=truncation, rays=rays
    )
    measured = {
        **metrics.to_dict(),
        "final_batch_mse": outcome.final_loss,
        "training_seconds": outcome.elapsed_seconds,
    }
    checks = (
        _check("field_mae", metrics.field_mae, gates.field_mae),
        _check("boundary_error_mean", metrics.boundary_error_mean, gates.boundary_error_mean),
        _check("ray_coverage", metrics.ray_coverage, gates.ray_coverage, maximum=False),
    )
    if _failed(checks):
        raise NeuralObjectError(
            "sphere fit did not meet its gates: " + ", ".join(_failed(checks))
        )

    payload, digest = _weights_bytes(outcome.model)
    bounds_min, bounds_max = _sphere_bounds(radius)
    revision = 1
    training = {
        "schema_version": SCHEMA_VERSION,
        "method": "supervised_tsdf",
        "optimizer": "Adam",
        "learning_rate": 1e-3,
        "steps": steps,
        "seed": seed,
        "batch": 2048,
        "final_batch_mse": outcome.final_loss,
        "elapsed_seconds": outcome.elapsed_seconds,
        "equations": "s(x) = clip(|x| - r, -tau, tau); r sampled from [0.35, 0.65]",
    }
    validation = {
        "schema_version": SCHEMA_VERSION,
        "revision": revision,
        "weights_sha256": digest,
        "status": "passed",
        "metrics": measured,
        "tolerances": gates.to_dict(),
        "checks": [check.to_dict() for check in checks],
        "failed_checks": [],
    }
    state = _state_document(radius)
    manifest = build_manifest(
        object_id=f"urn:tensorcad:geometry:sphere-r{radius:g}",
        revision=revision,
        weights=PayloadRef(WEIGHTS_FILE, digest, len(payload)),
        architecture_config=(config or FieldConfig()).to_dict(),
        field_spec=_field_spec(truncation),
        coordinates=_coordinates(bounds_min, bounds_max),
        state=state,
        training=training,
        validation=validation,
        lineage=Lineage(parent_revision=None, root_revision=revision, operation="fit"),
        provenance={
            "producer": "TensorCAD",
            "experiment": "neural-geometry-v0.2",
            "created_utc": _utc_now(),
            "command": _command_line(command),
            "architecture_id": ARCHITECTURE_ID,
        },
    )
    record = TransactionRecord(
        revision=revision,
        parent_revision=None,
        status=COMMITTED,
        operation="fit",
        weights_sha256=digest,
        metrics=measured,
        recorded_utc=_utc_now(),
        command=_command_line(command),
        elapsed_seconds=outcome.elapsed_seconds,
    )
    files: dict[str, bytes | str] = {
        WEIGHTS_FILE: payload,
        STATE_FILE: _json_text(state),
        TRAINING_FILE: _json_text(training),
        VALIDATION_FILE: _json_text(validation),
    }
    published = store.publish(
        files=files, manifest=manifest, expected_parent=None, record=record
    )
    return FitResult(
        store=store.root,
        revision=published.revision,
        weights_sha256=digest,
        metrics=measured,
        checks=checks,
        elapsed_seconds=outcome.elapsed_seconds,
        final_loss=outcome.final_loss,
    )


def _json_text(data: dict[str, Any]) -> str:
    from ..neural_geometry import dumps_canonical

    return dumps_canonical(data)


# --------------------------------------------------------------------------- #
# edit
# --------------------------------------------------------------------------- #


def apply_dent_edit(
    target: str | Path,
    *,
    depth: float,
    sigma: float,
    direction: Sequence[float] = (1.0, 0.0, 0.0),
    steps: int = 1000,
    seed: int = 503,
    tolerance: Tolerances | None = None,
    command: Sequence[str] | None = None,
    extract_resolution: int | None = None,
    rays: int = 4096,
) -> EditResult:
    """Apply a local dent as a transaction: train, measure, commit or roll back.

    On rejection the active revision is left byte-identical and the attempt is
    appended to ``history.json`` with the measured metrics and the names of the
    failed gates — a failed edit is evidence, not silence.
    """
    gates = tolerance or Tolerances()
    store = ObjectStore(target)
    parent_manifest = store.manifest()
    parent_revision = parent_manifest.revision
    parent_path = store.root / parent_manifest.weights.path
    parent_hash_before = store.active_weights_sha256()
    config = FieldConfig.from_dict(parent_manifest.architecture_config)
    truncation = float(parent_manifest.field["truncation_distance"])
    parent = load_weights(parent_path, config)
    radius = field_state_value(parent_manifest.raw.get("state", {}))
    normalise_direction(direction)

    outcome = fine_tune_dent(
        parent,
        radius=radius,
        depth=depth,
        sigma=sigma,
        steps=steps,
        seed=seed,
        truncation=truncation,
    )
    metrics = measure_dent_edit(
        parent,
        outcome.model,
        radius=radius,
        depth=depth,
        sigma=sigma,
        seed=seed + 99,
        rays=rays,
    )
    checks = (
        _check("ray_coverage", metrics["ray_coverage"], gates.ray_coverage, maximum=False),
        _check("dent_peak_error", metrics["dent_peak_error_m"], gates.dent_peak_error),
        _check(
            "dent_region_error_mean", metrics["dent_region_error_mean_m"], gates.dent_region_error_mean
        ),
        _check(
            "protected_mean_displacement",
            metrics["protected_mean_displacement_m"],
            gates.protected_mean_displacement,
        ),
        _check(
            "protected_p99_displacement",
            metrics["protected_p99_abs_displacement_m"],
            gates.protected_p99_displacement,
        ),
        _check(
            "protected_max_displacement",
            metrics["protected_max_abs_displacement_m"],
            gates.protected_max_displacement,
        ),
    )
    failed = _failed(checks)
    next_revision = parent_revision + 1
    recorded = {
        **metrics,
        "training_seconds": outcome.elapsed_seconds,
        "final_batch_mse": outcome.final_loss,
        "dent_depth_m": depth,
        "dent_sigma_radians": sigma,
    }

    if failed:
        record = TransactionRecord(
            revision=next_revision,
            parent_revision=parent_revision,
            status=ROLLED_BACK,
            operation="dent_edit",
            weights_sha256=None,
            metrics=recorded,
            failed_checks=failed,
            reason="validation gates not met: " + ", ".join(failed),
            recorded_utc=_utc_now(),
            command=_command_line(command),
            elapsed_seconds=outcome.elapsed_seconds,
        )
        store.rollback(None, record=record)
        assert store.active_weights_sha256() == parent_hash_before
        return EditResult(
            status=ROLLED_BACK,
            parent_revision=parent_revision,
            active_revision=store.active_revision(),
            weights_sha256=parent_hash_before,
            metrics=recorded,
            checks=checks,
            elapsed_seconds=outcome.elapsed_seconds,
            reason=record.reason,
        )

    payload, digest = _weights_bytes(outcome.model)
    validation = {
        "schema_version": SCHEMA_VERSION,
        "revision": next_revision,
        "weights_sha256": digest,
        "status": "passed",
        "metrics": recorded,
        "tolerances": gates.to_dict(),
        "checks": [check.to_dict() for check in checks],
        "failed_checks": [],
        "comparison_basis": (
            "protected displacement is measured against the parent revision at the same ray "
            "directions, so the field's own approximation error is not counted as edit drift"
        ),
    }
    state = _state_document(
        radius,
        controls_extra={
            "dent": {
                "value": depth,
                "unit": "m",
                "role": "online_edit",
                "edit_mode": "online_fit",
                "sigma_radians": sigma,
                "direction": [float(component) for component in normalise_direction(direction)],
            }
        },
    )
    training = {
        "schema_version": SCHEMA_VERSION,
        "method": "online_fit_with_replay",
        "optimizer": "Adam",
        "learning_rate": 2e-4,
        "steps": steps,
        "seed": seed,
        "batch": 2048,
        "replay_weight": 3.0,
        "elapsed_seconds": outcome.elapsed_seconds,
        "final_batch_mse": outcome.final_loss,
        "parent_revision": parent_revision,
        "parent_weights_sha256": parent_hash_before,
        "equation": "surface target r(a) = R - depth*exp(-a^2/(2*sigma^2))",
    }

    # Derived artifacts are produced in a scratch directory and shipped inside the
    # staged bundle: writing them into the live object before the commit would
    # make a rejected transaction leave files behind and would let a reader see a
    # mesh no manifest describes.
    carried_records, derived_files = _carry_forward_derived(
        store, parent_manifest, parent_revision
    )
    new_records: list[dict[str, Any]] = []
    if extract_resolution:
        scratch = Path(tempfile.mkdtemp(prefix="tensorcad-derived-"))
        try:
            _mesh, _report, records, payloads = _extract_and_write(
                manifest_state=state,
                weights_bytes=payload,
                config=config,
                truncation=truncation,
                bounds_min=parent_manifest.raw["coordinates"]["bounds"]["min"],
                bounds_max=parent_manifest.raw["coordinates"]["bounds"]["max"],
                resolution=extract_resolution,
                directory=scratch,
                formats=("obj", "stl"),
            )
        finally:
            shutil.rmtree(scratch, ignore_errors=True)
        new_records = list(records)
        derived_files.update(payloads)

    # A regenerated artifact replaces the stale record for the same path.
    replaced = {record["path"] for record in new_records}
    derived_records = [record for record in carried_records if record["path"] not in replaced]
    derived_records.extend(new_records)
    derived = tuple(derived_records)

    manifest = build_manifest(
        object_id=parent_manifest.object_id,
        revision=next_revision,
        weights=PayloadRef(WEIGHTS_FILE, digest, len(payload)),
        architecture_config=config.to_dict(),
        field_spec=_field_spec(truncation),
        coordinates=parent_manifest.raw["coordinates"],
        state=state,
        training=training,
        validation=validation,
        lineage=Lineage(
            parent_revision=parent_revision,
            root_revision=parent_manifest.lineage.root_revision or parent_revision,
            operation="dent_edit",
        ),
        provenance={
            "producer": "TensorCAD",
            "experiment": "neural-geometry-v0.2",
            "created_utc": _utc_now(),
            "command": _command_line(command),
            "architecture_id": ARCHITECTURE_ID,
        },
        derived=list(derived),
    )
    files: dict[str, bytes | str] = {
        WEIGHTS_FILE: payload,
        STATE_FILE: _json_text(state),
        TRAINING_FILE: _json_text(training),
        VALIDATION_FILE: _json_text(validation),
        **derived_files,
    }

    record = TransactionRecord(
        revision=next_revision,
        parent_revision=parent_revision,
        status=COMMITTED,
        operation="dent_edit",
        weights_sha256=digest,
        metrics=recorded,
        recorded_utc=_utc_now(),
        command=_command_line(command),
        elapsed_seconds=outcome.elapsed_seconds,
    )
    published = store.publish(
        files=files, manifest=manifest, expected_parent=parent_revision, record=record
    )
    return EditResult(
        status=COMMITTED,
        parent_revision=parent_revision,
        active_revision=published.revision,
        weights_sha256=digest,
        metrics=recorded,
        checks=checks,
        elapsed_seconds=outcome.elapsed_seconds,
        derived=derived,
    )


# --------------------------------------------------------------------------- #
# extraction
# --------------------------------------------------------------------------- #


def build_field(manifest: Any, weights_path: str | Path) -> Any:
    """Load the weights of a manifest and return a numpy field callable."""
    config = FieldConfig.from_dict(manifest.architecture_config)
    model = load_weights(weights_path, config)
    radius = field_state_value(manifest.raw.get("state", {}))
    return field_callable(model, radius)


def _extract_and_write(
    *,
    manifest_state: dict[str, Any],
    weights_bytes: bytes,
    config: FieldConfig,
    truncation: float,
    bounds_min: Sequence[float],
    bounds_max: Sequence[float],
    resolution: int,
    directory: Path,
    formats: Sequence[str],
) -> tuple[Mesh, ExtractionReport, tuple[dict[str, Any], ...], dict[str, bytes]]:
    """Extract a mesh, write it under ``directory`` and return the file bytes.

    The bytes are what a transaction stages: callers never point the extraction at
    the live bundle, because a candidate that is later rejected must not leave
    anything behind.
    """
    directory.mkdir(parents=True, exist_ok=True)
    scratch = Path(tempfile.mkdtemp(prefix="tensorcad-extract-"))
    try:
        weights_path = scratch / WEIGHTS_FILE
        weights_path.write_bytes(weights_bytes)
        model = load_weights(weights_path, config)
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
    radius = field_state_value(manifest_state)
    field = field_callable(model, radius)
    mesh = extract_surface(
        field,
        bounds_min,
        bounds_max,
        resolution,
        units="m",
        coordinate_frame="object_local",
    )
    report = summarise(mesh, resolution, expected_units="m")
    metrics = report.to_dict()
    records: list[dict[str, Any]] = []
    payloads: dict[str, bytes] = {}
    for fmt in formats:
        if fmt == "obj":
            path = write_obj(mesh, directory / "surface.obj")
        elif fmt == "stl":
            path = write_stl(mesh, directory / "surface.stl")
        else:  # pragma: no cover - guarded by the CLI choices
            raise ValueError(f"unsupported derived format {fmt!r}")
        data = path.read_bytes()
        relative = f"derived/{path.name}"
        payloads[relative] = data
        records.append(
            derived_record(
                kind="surface_mesh",
                fmt=fmt,
                path=relative,
                digest=sha256_bytes(data),
                size_bytes=len(data),
                metrics=metrics,
            )
        )
    return mesh, report, tuple(records), payloads


def _carry_forward_derived(
    store: ObjectStore, manifest: Any, parent_revision: int
) -> tuple[list[dict[str, Any]], dict[str, bytes]]:
    """Keep the parent revision's derived artifacts instead of deleting them.

    A commit swaps the whole bundle directory, so anything the new manifest does
    not reference disappears. Carrying the records forward -- marked with the
    revision they were derived from -- keeps an export discoverable and keeps it
    honest about being stale rather than silently dropping it.
    """
    records: list[dict[str, Any]] = []
    payloads: dict[str, bytes] = {}
    for artifact in manifest.derived:
        relative = artifact.get("path")
        if not isinstance(relative, str):
            continue
        source = store.root / relative
        if not source.is_file():
            continue
        payloads[relative] = source.read_bytes()
        records.append({**artifact, "source_revision": parent_revision})
    return records, payloads


def extract_object(
    target: str | Path,
    *,
    resolution: int = 96,
    formats: Sequence[str] = ("obj", "stl"),
    output: str | Path | None = None,
) -> tuple[Mesh, ExtractionReport, tuple[dict[str, Any], ...]]:
    """Derive a mesh from the object's active revision.

    Without ``output`` the artifacts are published *into the object* under
    ``derived/`` as part of a new revision, so ``manifest.derived`` describes them
    and a later commit cannot orphan them. With ``output`` they are written to a
    directory of the caller's choosing and are not registered in the bundle --
    that is an export, not a version.

    Either way the field stays the source of truth: the records carry the hash of
    the exact bytes and the metrics of the mesh that was extracted.
    """
    store = ObjectStore(target)
    manifest = store.manifest()
    weights_bytes = (store.root / manifest.weights.path).read_bytes()
    config = FieldConfig.from_dict(manifest.architecture_config)
    scratch = Path(tempfile.mkdtemp(prefix="tensorcad-derived-"))
    try:
        mesh, report, records, payloads = _extract_and_write(
            manifest_state=manifest.raw.get("state", {}),
            weights_bytes=weights_bytes,
            config=config,
            truncation=float(manifest.field["truncation_distance"]),
            bounds_min=manifest.raw["coordinates"]["bounds"]["min"],
            bounds_max=manifest.raw["coordinates"]["bounds"]["max"],
            resolution=resolution,
            directory=scratch,
            formats=formats,
        )
    finally:
        shutil.rmtree(scratch, ignore_errors=True)

    if output is not None:
        directory = Path(output)
        directory.mkdir(parents=True, exist_ok=True)
        for relative, data in payloads.items():
            (directory / Path(relative).name).write_bytes(data)
    else:
        store.publish_derived(files=payloads, derived=records, operation="extract")
    return mesh, report, records


class NeuralGeometryField(GeometryField):
    """A published ``.ngo`` revision exposed through the shared field protocol.

    The point is interchangeability: benchmark code, export code and validation
    code should not care whether a surface came from ``|x|-r`` or from a network.
    Two differences are stated rather than hidden:

    * ``quantity`` is ``"truncated_signed_distance"`` -- outside the truncation
      band the magnitude is meaningless, so a comparison against an exact SDF has
      to compare *surfaces*, not values;
    * ``gradient`` is a central difference of the network (torch autograd would be
      exact but is not available on the NumPy-only extraction path).

    ``with_state`` edits the conditioning value of the frozen weights, which is the
    cheap parameter edit the format was designed around -- no retraining happens.
    """

    kind: ClassVar[str] = "neural_geometry_object"
    quantity: ClassVar[str] = TRUNCATED_SIGNED_DISTANCE

    def __init__(
        self, manifest: Any, weights_path: str | Path, *, state_value: float | None = None
    ) -> None:
        config = FieldConfig.from_dict(manifest.architecture_config)
        self.manifest = manifest
        self._weights_path = Path(weights_path)
        self._model = load_weights(self._weights_path, config)
        self._state_value = (
            field_state_value(manifest.raw.get("state", {})) if state_value is None else state_value
        )
        self._evaluate = field_callable(self._model, self._state_value)
        self._gradient: Callable[[Any], Any] | None = gradient_callable(self._model, self._state_value)
        bounds = manifest.raw["coordinates"]["bounds"]
        self._bounds = BoundingBox(tuple(bounds["min"]), tuple(bounds["max"]))
        self.units = str(manifest.field.get("length_unit", "m"))
        self.coordinate_frame = str(manifest.coordinates.get("frame", "object_local"))
        self.truncation_distance = float(manifest.field.get("truncation_distance", 0.0))

    def signed_distance(self, points: Iterable[Sequence[float]]) -> list[float]:
        import numpy as np

        batch = [(float(x), float(y), float(z)) for x, y, z in points]
        if not batch:
            return []
        values = np.asarray(self._evaluate(np.asarray(batch, dtype=np.float64)), dtype=np.float64)
        return [float(value) for value in values]

    def gradient(self, points: Iterable[Sequence[float]]) -> list[Point]:
        """Exact gradient of the network, via autograd -- not a finite difference.

        A float32 network makes central differences unreliable below ~1e-4 of the
        object size, so the protocol's numerical fallback is kept only as a
        fallback.
        """
        import numpy as np

        batch = _points(points)
        if not batch:
            return []
        if self._gradient is not None:
            values = np.atleast_2d(self._gradient(np.asarray(batch, dtype=np.float64)))
            return [(float(row[0]), float(row[1]), float(row[2])) for row in values]
        return self.numerical_gradient(batch, step=1e-4 * max(self._bounds.sizes))

    def bounds(self) -> BoundingBox:
        return self._bounds

    def state(self) -> dict[str, float]:
        return {"radius": self._state_value}

    def with_state(self, **values: float) -> NeuralGeometryField:
        unknown = set(values) - {"radius"}
        if unknown:
            raise FieldError(f"{self.kind} has no control(s) {sorted(unknown)}")
        return NeuralGeometryField(
            self.manifest, self._weights_path, state_value=float(values.get("radius", self._state_value))
        )

    def parameters(self) -> dict[str, Any]:
        return {
            "architecture_id": self.manifest.architecture_id,
            "weights_sha256": self.manifest.weights.sha256,
            "revision": self.manifest.revision,
            "truncation_distance": self.truncation_distance,
        }

    def provenance(self) -> dict[str, Any]:
        return {
            "producer": "TensorCAD",
            "representation": "neural_field",
            "object_id": self.manifest.object_id,
            "revision": self.manifest.revision,
            "weights_sha256": self.manifest.weights.sha256,
        }

    def numpy_callable(self) -> Any:
        return self._evaluate


def field_from_object(target: str | Path) -> NeuralGeometryField:
    """Open a bundle and present its active revision as a :class:`GeometryField`."""
    store = ObjectStore(target)
    manifest = store.manifest()
    return NeuralGeometryField(manifest, store.root / manifest.weights.path)


def history_entries(target: str | Path) -> list[dict[str, Any]]:
    return ObjectStore(target).history()


__all__ = [
    "HISTORY_FILE",
    "Check",
    "EditResult",
    "FitResult",
    "NeuralGeometryField",
    "NeuralObjectError",
    "Tolerances",
    "apply_dent_edit",
    "build_field",
    "extract_object",
    "field_from_object",
    "field_state_value",
    "fit_sphere",
    "history_entries",
]
