# Agent instructions

Before making changes, read the [Tensor ecosystem vision](https://github.com/cxs503/TensorPlatform/blob/main/docs/architecture/TENSOR_ECOSYSTEM_AI_ERA_VISION.md) and this repository's [project charter](docs/PROJECT_CHARTER.md). These documents define ownership boundaries; the source code and tests define what is actually implemented.

## Required development behavior

- State the goal, why it belongs in this repository, and which public or cross-repository contracts it affects.
- Preserve the repository boundary: TensorCAD owns computational geometry representation, editable geometry state, mesh/topology checks and solver handoff artifacts. Read docs/PROJECT_CHARTER.md and the ecosystem vision. Do not implement a general agent platform or claim hydrodynamic/structural suitability from geometry checks. Use independent references, round-trip tests, shape drift and explicit export/engineering limitations.
- Inspect relevant code, tests, recent changes, and existing benchmark records before editing.
- Prefer stable, versioned public interfaces; avoid private cross-repository imports and duplicated sources of truth.
- Include runnable tests/examples and quantitative evidence appropriate to the change. Report exact commands, results, failures, untested paths, assumptions, and limitations.
- Distinguish implementation, unit tests, integration tests, benchmarks, domain validation, and production deployment evidence. Never claim more than the evidence establishes.
- Keep changes scoped, reviewable, and reversible. Do not optimize for PR count, code volume, or workflow completion alone.
