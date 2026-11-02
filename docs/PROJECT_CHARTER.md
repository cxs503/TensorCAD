# TensorCAD — computational geometry charter

**Canonical ecosystem vision:** [Tensor Ecosystem: Vision for the AI Era](https://github.com/cxs503/TensorPlatform/blob/main/docs/architecture/TENSOR_ECOSYSTEM_AI_ERA_VISION.md)

## Mission
Provide agent-callable computational geometry representations, parameterized editing, validation, and solver-ready artifacts, including exploration of continuous/neural geometry objects.

## Owns
- Geometry schemas, analytic/continuous fields, editable state and geometry-specific transactions.
- Mesh generation/conditioning, topology and shape validation, geometric metrics and export contracts.
- Geometry-specific reproducibility, provenance and solver handoff manifests.

## Does not own
- A general-purpose agent platform or a full traditional CAD desktop application by default.
- Downstream fluid/structural solver execution or claims of engineering performance based only on geometry checks.
- Treating neural-field output as accurate without comparison against declared references and thresholds.

## Integration
Expose stable geometry queries, explicit units and coordinate conventions, artifact hashes, validation reports and solver handoffs. TensorFEM/TensorLBM and other consumers own their solver interpretation; TensorPhysics owns scientific workflow validation; TensorPlatform owns shared general infrastructure.

## Quality gates
Use independent analytic/geometric references, topology/manifold checks, resolution studies, round-trip export checks, shape-drift measurements and explicit limitations. Geometry validity does not establish hydrodynamic or structural suitability.

## Agent rules
Read this charter and README; state ownership and affected interfaces; avoid private cross-project coupling; add reproducible tests and measured benchmarks; report failures and untested paths; separate geometry verification from downstream engineering validation.

## Maturity
This charter states boundaries; feature maturity and engineering readiness must be assessed from current benchmark evidence.
