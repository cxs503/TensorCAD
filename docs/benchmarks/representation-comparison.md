# Geometry Representation Benchmark

**Status:** benchmark contract / implementation plan. This document defines the comparison; it does not claim the full cross-representation benchmark has already run.

## Purpose

TensorCAD already has measured geometry-core tests and a neural geometry-object edit experiment. The next step is to combine them under a common benchmark contract, then add optional CAD-kernel and solver handoff paths. Avoid duplicating measurements or comparing incompatible quantities.

The benchmark answers: **which representation should an agent choose for a given engineering edit, and what evidence supports that choice?**

## Existing measured entry points

Run these independently first and preserve their current acceptance gates:

- `python examples/geometry_core/run_benchmark.py` — analytic primitives and the learned field through the common GeometryField protocol.
- `python examples/neural_geometry/run_object_pipeline.py --output results/object` — neural object fit, local edit, commit/rollback, artifact hashes and derived mesh.
- `python examples/step_export/run_example.py` — Gmsh/JSON FEM handoff; optional OCP sewing and STEP round-trip checks.
- `python examples/parametric_hull/run_example.py` — parametric hull mesh validity and hydrostatics/resolution study.

The unified benchmark should consume structured outputs, not parse console text, where practical.

## Representations to compare

| ID | Representation | Role / limitation |
|---|---|---|
| `analytic_field` | Exact analytic field | Ground truth and exact parameter edits where a closed form exists |
| `neural_field` | Versioned `.ngo` neural implicit object | Learned shape and local edits; quantify approximation and locality errors |
| `mesh` | Explicit triangle mesh | Mesh editing/conditioning baseline; no hidden claim of design intent |
| `ocp_step` | OCP/Open CASCADE export | Optional kernel path; label faceted solids as faceted, never as analytic feature B-rep |
| `parametric_hull` | Existing parameterized hull model | Engineering case with measured hydrostatics and solver-facing outputs |

A capability unavailable for a representation must be reported as `not_supported`, not silently approximated or omitted.

## Benchmark cases

### Case A — localized edit with analytic ground truth

Start with the existing sphere + dent experiment. Use a fixed seed and probe set. Measure the target dent and protected surface separately. Preserve existing gates and add the unified result record.

### Case B — geometry export and round trip

Use the existing SUBOFF AFF-8 surface. Validate topology, component/patch identity, volume/area drift, and export/import integrity. Distinguish mesh-to-faceted-STEP from a true editable feature model.

### Case C — solver handoff

After A and B pass, hand one geometry to a real Tensor solver (prefer TensorFEM for the first deterministic case). Record solver version/configuration, boundary conditions, mesh, convergence status, physical observables, and comparison against a documented reference. A screenshot or successful export is not a solver result.

## Common metrics

Every metric must name its units, reference, sampling method, and tolerance.

- **Shape:** symmetric surface-distance mean/P95/max; bounding-box and key-dimension error; relative volume and area error where meaningful.
- **Edit semantics:** target-region edit error; protected-region displacement mean/P95/max; locality ratio; parameter recovery; committed/rejected/rolled-back outcome.
- **Validity:** watertightness, edge/vertex manifoldness, orientation, Euler characteristic, degenerate/sliver faces, and kernel validity where available.
- **Cost:** wall time by stage, peak memory when measurable, training steps/time, query throughput, and artifact bytes.
- **Handoff/physics:** round-trip displacement, patch/region preservation, solver convergence, and selected physical quantities versus reference.
- **Provenance:** input/config hash, geometry and weight hashes, seed, Python/package/kernel/solver versions, device, and command line.

Do not compare raw field values between an exact signed-distance field and a truncated signed-distance field. Compare extracted surface positions and record each field's `quantity`.

## Result format

The unified runner should emit:

1. `results.json` — stable machine-readable schema with run metadata, per-representation metrics, gate outcomes, and artifact hashes.
2. `report.md` — generated summary table, limitations, and measured recommendation.
3. Plots generated from `results.json`, e.g. surface error versus runtime and protected-region displacement.
4. Optional exported artifacts in a run-specific output directory; avoid committing large generated meshes/weights.

Each metric should carry `value`, `unit`, `reference`, `tolerance` when defined, and `status`. Missing metrics use `not_supported` or `not_run`; never encode missing data as zero.

## External benchmark to learn from

[CADGenBench](https://github.com/huggingface/cadgenbench) provides precedents for validity-gated scoring, surface similarity, interface matching, and topology metrics. TensorCAD should borrow applicable metric definitions, while adding local-edit preservation, transactional rollback, provenance, solver handoff, and engineering observables.

## Gates and CI

- CI runs a small deterministic CPU smoke case with fixed probes and numerical gates.
- Neural training, OCP/STEP, and solver runs may be optional jobs, but their commands and dependencies must be documented.
- A gate failure fails the relevant case and remains visible in the report; aggregate scores must not hide validity or protected-region failures.
- Outputs are generated by the code under test, never hand-entered or mocked.
- Keep expensive performance numbers informational unless hardware/runtime are controlled; correctness gates should be portable.

## Implementation sequence

1. Re-run the existing entry points above and record failures/dependency availability.
2. Define a versioned JSON schema and adapt existing benchmark scripts without weakening their gates.
3. Implement Case A first: analytic reference + neural field + mesh baseline.
4. Add optional OCP/STEP round-trip measurements.
5. Add one real TensorFEM solver case only after geometry/handoff gates pass.
6. Publish measured results and a representation-selection recommendation.

Related: [Issue #19](https://github.com/cxs503/TensorCAD/issues/19), [Issue #6](https://github.com/cxs503/TensorCAD/issues/6), [Issue #4](https://github.com/cxs503/TensorCAD/issues/4).
