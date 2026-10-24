# TensorCAD Neural Geometry Object Metadata (NGOM) v0.1

Status: proposed schema for the first executable prototype. This is a versioned object format, not yet a finalized public standard.

## Design principles

1. The canonical object is a continuous neural field plus its interpretation metadata.
2. Weights alone are insufficient: architecture, tensor names/shapes/dtypes, field semantics, coordinate frame, units, state, and training/validation provenance must be recoverable.
3. One field-query contract applies to every geometry. Shape labels are descriptive metadata, never runtime backend selectors.
4. Metadata is declarative. Loading an object must not execute code from the object bundle.
5. Updates are versioned transactions. An unvalidated candidate must never overwrite the last valid object.
6. Meshes and STEP/IGES exports are derived artifacts; they are not the canonical geometry state.

## Suggested package layout

    part.ngo/
      manifest.json
      weights.safetensors
      state.json
      training.json
      validation.json
      artifacts/
        preview.glb
        solver_mesh.vtu

The .ngo directory/package name is provisional. Start with a directory bundle; a single-file container can come later.

## manifest.json

The manifest is the entry point. It describes how to interpret every payload without embedding executable code.

    {
      "format": "tensorcad.neural-geometry",
      "format_version": "0.1.0",
      "object_id": "urn:tensorcad:geometry:example-sphere",
      "revision": 1,
      "model": {
        "family": "conditional_implicit_field",
        "architecture_id": "tensorcad.mlp_fourier.v1",
        "architecture_config": {
          "coordinate_dim": 3,
          "state_dim": 1,
          "hidden_width": 64,
          "hidden_layers": 4,
          "coordinate_encoding": "fourier",
          "output_dim": 1
        },
        "weights": {
          "path": "weights.safetensors",
          "sha256": "..."
        }
      },
      "field": {
        "quantity": "truncated_signed_distance",
        "inside_sign": "negative",
        "surface_isovalue": 0.0,
        "length_unit": "m",
        "truncation_distance": 0.1
      },
      "coordinates": {
        "dimension": 3,
        "frame": "object_local",
        "axis_convention": "right_handed_xyz",
        "normalization": {
          "kind": "affine",
          "center": [0.0, 0.0, 0.0],
          "scale": [1.0, 1.0, 1.0]
        }
      },
      "state_file": "state.json",
      "training_file": "training.json",
      "validation_file": "validation.json",
      "created_utc": "2026-10-24T00:00:00Z",
      "producer": {
        "name": "TensorCAD",
        "version": "0.1.0"
      }
    }

The hash is a placeholder and must be generated from the actual weight bytes. The architecture identifier resolves only through an installed, allow-listed implementation registry. The object file must not provide arbitrary Python code, dynamic imports, or executable model definitions.

## state.json

State records editable controls and their semantic constraints, independently of model weights.

    {
      "schema_version": "0.1.0",
      "controls": {
        "radius": {
          "value": 0.5,
          "unit": "m",
          "role": "global_parameter",
          "minimum": 0.05,
          "maximum": 5.0,
          "edit_mode": "conditioned"
        }
      },
      "pose": {
        "translation": [0.0, 0.0, 0.0],
        "rotation_xyzw": [0.0, 0.0, 0.0, 1.0]
      },
      "bounds": {
        "minimum": [-0.6, -0.6, -0.6],
        "maximum": [0.6, 0.6, 0.6]
      }
    }

edit_mode is descriptive metadata for the edit engine, not permission to switch geometry representations. Initial values may be conditioned (change state with frozen weights) and online_fit (update the same model weights against constraints). An operation can combine them: update state, then train and validate.

For local controls, use stable IDs and local coordinate frames. Avoid baking UI labels or screen coordinates into canonical state.

## training.json

Record enough to reproduce or audit the most recent training/edit transaction.

    {
      "schema_version": "0.1.0",
      "method": "supervised_tsdf",
      "optimizer": "Adam",
      "seed": 503,
      "steps": 1200,
      "learning_rate": 0.001,
      "losses": {
        "field_regression": 0.0,
        "surface": 0.0,
        "preservation": 0.0,
        "smoothness": 0.0
      },
      "source_constraints": [
        {
          "kind": "analytic_reference_samples",
          "description": "sphere reference used to generate labels"
        }
      ],
      "parent_revision": 0,
      "compute": {
        "device": "cpu",
        "elapsed_seconds": 0.0
      }
    }

These numbers are illustrative only; production output must contain measured values. If exact continuation is required, optimizer state and replay data belong in separately identified payloads with hashes. Do not put large training arrays directly in JSON.

## validation.json

Validation is explicit, machine-readable, and tied to the exact weight hash and revision.

    {
      "schema_version": "0.1.0",
      "revision": 1,
      "weights_sha256": "...",
      "status": "passed",
      "metrics": {
        "boundary_error_mean": 0.0,
        "boundary_error_p95": 0.0,
        "boundary_error_max": 0.0,
        "normal_error_degrees_mean": 0.0,
        "protected_region_displacement_max": 0.0,
        "inside_outside_accuracy": 1.0,
        "mesh_components": 1,
        "degenerate_triangle_count": 0
      },
      "tolerances": {
        "boundary_error_max": 0.001,
        "protected_region_displacement_max": 0.001
      },
      "checks": [
        {"name": "finite_field_values", "passed": true},
        {"name": "surface_extraction", "passed": true},
        {"name": "edit_preservation", "passed": true}
      ]
    }

Zeros and tolerance values above are examples, not established engineering acceptance thresholds. Tolerances must be selected for geometry scale and use case.

## Query contract

Every implementation must support the same logical interface:

- evaluate(points, object_revision) -> field_values
- gradient(points, object_revision) -> spatial_gradients
- extract_surface(object_revision, extraction_config) -> derived_mesh
- apply_edit(object_revision, constraint_spec) -> candidate_revision
- validate(candidate_revision, validation_profile) -> validation_report

The interface must not dispatch to an analytic sphere evaluator, a SUBOFF formula evaluator, or a separate mesh-based geometry backend. Explicit formulas can exist in data-generation and validation code only.

## Safe update protocol

1. Read and verify the current manifest, payload hashes, schema versions, and architecture allow-list.
2. Create a candidate revision without modifying the active revision.
3. Apply edit constraints and train the same field model.
4. Write candidate weights and metadata to new immutable payloads.
5. Validate finite values, field semantics, geometric errors, protected regions, and requested engineering checks.
6. If all required checks pass, atomically mark the candidate as active. Otherwise preserve the previous active revision and record the failure report.
7. Never treat a missing validation report as success.

## Schema validation rules

- Reject unknown major format versions; permit compatible minor versions only under explicit rules.
- Reject NaN/infinity in numeric metadata.
- Require declared units, coordinate frame, field sign convention, isovalue, architecture ID, and tensor payload hash.
- Verify all weight tensors match the allow-listed architecture's names, shapes, and dtypes.
- Do not load arbitrary serialized Python objects.
- Require a revision and parent revision for edits.
- Bind validation reports to exact weight hashes.
- Keep artifact paths relative to the bundle root and reject path traversal.
- Do not trust status=passed without re-running required validation in the consuming environment.

## First prototype deliverables

1. A typed Python model for manifest/state/training/validation records.
2. JSON Schema or equivalent validation for each metadata file.
3. A safe bundle loader that validates hashes and rejects unknown architecture IDs.
4. A sphere object generated, saved, reloaded, queried, and compared against reference values.
5. A revisioned online-edit transaction with validation and rollback.
6. Tests for malformed metadata, hash mismatch, unknown architecture, coordinate/unit mistakes, and failed edits.

The format is intentionally small. Evolve it from measured needs in the sphere-and-dent experiment, not speculative fields for every future CAD feature.
