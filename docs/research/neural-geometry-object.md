# Neural Geometry Object: Runtime Learning Mechanism

This document refines the unified-field experiment into an executable system concept. It does not claim that the mechanism is implemented yet.

## Core proposition

A TensorCAD geometry object is a learned continuous field, persisted primarily as model weights plus a versioned geometry state and model specification. Meshes and CAD exchange files are derived artifacts.

For a point x and object state c, the network predicts a scalar field:

    s = F_theta(x, c)

The zero set is the surface. Use a truncated signed-distance target initially, with negative values inside and positive values outside. Do not label predictions as exact signed distances unless measurements justify it.

## A geometry object is more than a weight file

Persist a small, versioned object bundle:

- model.safetensors (or another safe tensor-only format): network weights;
- geometry.json: object ID, units, coordinate frame, model architecture/version, named controls and current values, bounds, TSDF convention/truncation, training state and provenance;
- training_state (optional): optimizer state and replay samples if continued training must resume exactly;
- validation.json: measured boundary, normal, edit, and mesh-extraction metrics;
- derived exports (optional): mesh/STEP/IGES as compatibility outputs, never canonical truth.

Never serialize executable Python or load arbitrary pickled objects from untrusted geometry files.

## One model family

Use one shared conditional neural field architecture for all objects in the experiment. Shape identity and editable controls are inputs, not a switch to a shape-specific evaluator.

The first proof-of-concept may use a compact coordinate MLP with positional/Fourier features and a small state encoder. The exact width/depth is an experimental parameter, not a design dogma. Keep the field output scalar and use the same evaluator for all shapes.

To make the object compact and editable, distinguish:

- shared architecture and optional pretrained initialization;
- per-object weights or compact per-object learned parameters;
- named controls/state;
- observations/constraints used to update the object.

The first prototype should choose the simplest option that can prove frozen-weight parameter edits. Do not assume one global network will automatically generalize to every geometry class.

## Initialization and editing are training transactions

All geometry-changing operations enter through one training transaction:

1. Read current object weights and state; create a checkpoint.
2. Convert the requested operation into observations/constraints in the common field target: positive/negative samples, near-surface distance samples, boundary points, fixed regions, or user-provided points/curves.
3. Build a batch containing new constraints plus replay constraints sampled from the current object. Replay is important: fitting only the new indentation can destroy the rest of the shape (catastrophic forgetting).
4. Optimize the same field model for a bounded number of steps with a bounded time/compute budget.
5. Validate both the edited region and protected regions against explicit acceptance tolerances.
6. Commit a new object version only if validation passes; otherwise roll back to the prior weights and report failure.

This transaction is the core mechanism. The user-facing tool merely produces a constraint specification.

## Sphere experiment

For radius r, use the analytic sphere only as a reference generator:

    d(x) = norm(x - center) - r

Sample points around the surface and in a narrow volume band, train the neural field, and measure field/boundary error. Test on radii excluded from training.

There are two distinct tests:

A. **State-conditioned edit:** train one conditional field over a radius interval, freeze all weights, change r to an unseen value, and query the same network. This proves the network learned a useful control mapping.

B. **Online edit:** start with a fitted sphere and change its radius or add a local indentation by adding new constraints and fine-tuning the same object's weights. Validate that the intended region changes while protected regions remain stable.

A successful A does not prove B; both must be measured. The user’s proposed “train it live” mechanism is primarily B.

## Local indentation experiment

Represent “make a dent here” as a constraint operation, not as a separate geometry backend:

- center and frame of the edit region;
- target depth and radius/footprint;
- boundary samples or a target field in the affected region;
- protected-region constraints outside the edit region;
- smoothness/gradient regularization;
- tolerance and compute budget.

Training updates the field weights. The edit succeeds only if the dent reaches the target depth, the surrounding shape remains within tolerance, and the field does not introduce unintended components or unstable normals.

The edit tool does not directly manipulate network layers or construct a separate analytic dent primitive as the canonical shape. It specifies desired geometric constraints.

## Training losses

The first experiment should keep losses interpretable and log each term separately:

- field regression to TSDF targets;
- surface consistency for known boundary samples;
- spatial-gradient regularization near the surface, used cautiously;
- protected-region preservation/replay;
- optional smoothness regularization;
- optional edit-locality penalty.

Use explicit analytic references to generate labels and compute test errors, not as runtime fallbacks. PINN-style differential constraints may be explored when a meaningful governing equation or geometric constraint exists; a PINN is not automatically needed to learn a sphere.

## Safety and correctness gates

A training run must not overwrite the last valid geometry. Require:

- deterministic seeds and recorded model/config versions;
- frozen validation sets and explicit error tolerances;
- pre-edit checkpoint and rollback;
- no silent acceptance when validation fails;
- unit and coordinate-frame checks;
- mesh extraction checks for disconnected components, degenerate triangles, and normal consistency;
- separate reporting of field error and extracted-mesh error.

For engineering use, downstream solvers should receive a mesh extracted from a validated field. The network file itself is the canonical editable object, but it is not automatically a substitute for a watertight, solver-ready discretization.

## First acceptance criteria

1. Train a compact neural field for a sphere using generated reference data.
2. Save and reload the object bundle and reproduce query values within tolerance.
3. Change the radius to an unseen value with frozen weights and report boundary error.
4. Fine-tune the same object to produce a local indentation; report target dent depth, protected-region displacement, normal error, training time, and inference time.
5. Demonstrate rollback when a deliberately invalid edit fails validation.
6. Run the same field API and test harness for slab, cylinder, ellipsoid, then SUBOFF.

Do not declare the research successful based only on a visually plausible mesh. The decisive result is a compact, persistent neural geometry object that can be edited through constraints, validated quantitatively, and safely versioned.
