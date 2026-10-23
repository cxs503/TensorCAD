# TensorCAD

**TensorCAD is the agent-callable geometry and engineering-model foundation for the Tensor scientific-computing ecosystem.**

The goal is not to reproduce every feature of a traditional desktop CAD system. TensorCAD should let a research agent define geometry from explicit parameters, validate it, export solver-ready artifacts, and reproduce the exact model used by a simulation.

## Position in the Tensor ecosystem

- **TensorCAD** — reusable geometry, parameterization, validation, and pre-processing interfaces.
- **TensorShipDesign** — planned downstream ship-design application built on TensorCAD, starting with parametric hull geometry and hydrostatics.
- **TensorLBM / TensorFEM / TensorDEM / TensorBEM** — simulation engines that consume suitable geometry, meshes, and boundary/configuration data.
- **TensorPlatform** — orchestrates agent tasks, model generation, simulation, evaluation, versioning, and iterative design.

## Design principles

1. **Agent-first, headless-first:** core operations must be callable from Python and automation, without requiring interactive GUI use.
2. **Reproducibility:** record parameters, units, coordinate frame, backend/version, validation results, and generated artifacts.
3. **Geometry before optimization:** establish deterministic geometry and validity checks before adding differentiable optimization or learned surrogates.
4. **Reuse mature kernels:** evaluate Open CASCADE/OCP or other mature geometry libraries rather than rewriting a full B-rep kernel.
5. **Clear boundaries:** keep generic geometry services separate from ship-specific design logic and solver-specific adapters.
6. **Real verification:** examples and benchmarks must run actual code and report measured outputs; no mock results.

## Initial milestone

Deliver a minimal Python package and a headless, parameterized geometry example with tests, documented installation, and export/validation where supported. The first application scenario will be ship-hull geometry, with naval-architecture validation added incrementally.

## Proposed package boundaries

- `tensorcad.parameters` — parameter definitions, units, validation.
- `tensorcad.geometry` — geometry backend abstraction and model results.
- `tensorcad.validation` — geometric and engineering checks.
- `tensorcad.export` — artifact export and metadata.
- `tensorcad.adapters` — optional solver/meshing integrations.
- `examples/` and `tests/` — runnable examples and regression tests.

These are initial boundaries, not a commitment to prematurely freeze every module.

## Status

New repository; implementation and backend compatibility are to be verified. The initial benchmark must not claim engineering validity solely from visual appearance.
