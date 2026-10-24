# TensorCAD

**TensorCAD is the agent-callable geometry and engineering-model foundation for the Tensor scientific-computing ecosystem.**

The goal is not to reproduce every feature of a traditional desktop CAD system. TensorCAD should let a research agent define geometry from explicit parameters, validate it, export solver-ready artifacts, and reproduce the exact model used by a simulation.

[![tests](https://github.com/cxs503/TensorCAD/actions/workflows/tests.yml/badge.svg)](https://github.com/cxs503/TensorCAD/actions/workflows/tests.yml)

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

Principle 6 is enforced in CI: `examples/parametric_hull/run_example.py` prints measured numbers and fails the build if any of its own claims stop being true.

## What is implemented

`tensorcad.hull` generates a symmetric displacement-hull surface from principal dimensions and hands back a **closed, oriented solid** — not an open shell:

| Property of the exported mesh | Status |
|---|---|
| Watertight (every edge used by exactly two faces) | checked, and generation raises if it fails |
| Edge- and vertex-manifold | checked (the zero-breadth stem/post is reported, not hidden) |
| Consistently outward-oriented (signed volume > 0) | checked |
| Euler characteristic V−E+F = 2 | checked |
| No zero-area faces | checked |
| Realised L/B/T agrees with the request | checked, with the grid-sampling deficit named explicitly |
| Discretisation error measured against a refined grid | reported in `metadata.json` |

Hydrostatics are integrated from the mesh itself (divergence theorem), so the numbers describe the artifact a solver will receive rather than a formula it never got:

```python
from tensorcad import HullParameters, generate_hull

model = generate_hull(HullParameters(length=10.0, beam=2.0, draft=0.8, nx=41, nz=13))
model.hydrostatics.displacement_volume      # 7.446554 m^3  (analytic form: 7.512698)
model.hydrostatics.block_coefficient        # 0.4654
model.hydrostatics.prismatic_coefficient    # 0.7031
model.hydrostatics.longitudinal_centre_buoyancy
model.validation.warnings                   # [] at nx=41; non-empty when the grid hides something
```

Measured at the default resolution, the sampled surface sits 0.88 % below the analytic volume of its own generating functions; at `nx=161, nz=41` the gap is 0.14 %. That is reported per export rather than asserted in prose.

## Neural geometry objects (0.3)

A hull mesh is one way to store geometry. A `.ngo` object is another: the geometry
is a learned continuous field plus the editable state that conditions it, and every
mesh — OBJ, STL, whatever a solver wants — is a derived export.

```python
from tensorcad.neural.objects import fit_sphere, apply_dent_edit, extract_object

fit_sphere("results/sphere.ngo", radius=0.5, steps=1200, seed=503)          # revision 1
edit = apply_dent_edit("results/sphere.ngo", depth=0.06, sigma=0.22)        # transaction
edit.status                    # 'committed' -> revision 2, or 'rolled_back' and nothing changes
mesh, report, artifacts = extract_object("results/sphere.ngo", resolution=96)
report.watertight, report.euler_characteristic                              # True, 2
```

```bash
tensorcad-neural fit     --output results/sphere.ngo --radius 0.5 --steps 1200
tensorcad-neural edit    --object results/sphere.ngo --depth 0.06 --sigma 0.22 --resolution 96
tensorcad-neural extract --object results/sphere.ngo --resolution 96 --format obj,stl
tensorcad-neural info    --object results/sphere.ngo        # manifest + every attempt
tensorcad-neural verify  --object results/sphere.ngo        # re-open, re-hash, re-query
```

An object is a directory bundle: `manifest.json`, `weights.safetensors`, `state.json`,
`training.json`, `validation.json`, `history.json`, and optional `derived/` meshes.
Edits are transactions: a candidate is trained in a hidden staging directory, measured
against declared gates, and only then renamed into place; a rejection leaves the active
revision byte-identical and records the failed gates and their measurements in
`history.json`. Measured on the reference run: dent peak 0.0527 m against a 0.06 m
target, dent-region error 0.0028 m, protected-region displacement mean 0.0012 m /
P95 0.0030 m / max 0.0097 m, extracted mesh watertight with V−E+F = 2.

The metadata layer (`tensorcad.neural_geometry`) is standard-library only and never
deserializes executable code. The field runtime (`tensorcad.neural`) needs the
`neural` extra; importing `tensorcad` never pulls in torch or numpy. Format details,
field semantics and the transaction model are in
[docs/neural-geometry-format.md](docs/neural-geometry-format.md).

## Package layout

The boundaries below are the ones the code actually has (they were proposals in 0.1):

- `tensorcad.parameters` — parameter specs, units, ranges, dimensionless ratio checks.
- `tensorcad.geometry` — `Mesh`, `BoundingBox`, `MeshBuilder` (exact-position vertex sharing).
- `tensorcad.mesh` — generic algorithms: weld, orient outward, boundary/manifold/watertight analysis, surface integrals.
- `tensorcad.validation` — named geometric checks, errors versus warnings, `ValidationReport`.
- `tensorcad.hydrostatics` — mesh-integrated displacement, areas, centres, coefficients.
- `tensorcad.export` — OBJ/CSV/JSON writers and an OBJ reader (round-trip exact).
- `tensorcad.hull` — the one ship-specific model.
- `tensorcad.adapters` — optional solver/kernel bridges; the core must never import from here.
- `tensorcad.neural_geometry` — neural geometry object format: metadata, hashing, atomic
  revision publishing, edit transactions (standard library only).
- `tensorcad.neural` — optional field runtime: training, local edits, marching-tetrahedra
  extraction to a closed mesh (`neural` extra, consumes torch/numpy).
- `tensorcad.cli_neural` — `tensorcad-neural` entry point for fit/edit/extract/info/verify.
- `examples/` and `tests/` — runnable benchmark and 129 tests.

The core package depends on nothing outside the standard library. That is a tested guarantee, not a preference: a B-rep kernel stays an optional backend (see issue #4), and `tensorcad.adapters.ocp` is the planned seam.

## Install and run

```bash
python -m pip install -e ".[test,dev]"
pytest --cov
python examples/parametric_hull/run_example.py --output results/hull-demo
tensorcad-hull --length 12 --beam 2.4 --draft 0.9 --nx 61 --nz 17 --output results/hull-demo
tensorcad-hull --json            # full reproducibility record on stdout

# neural geometry objects (needs the torch-backed extra)
python -m pip install -e ".[neural]"
python examples/neural_geometry/run_object_pipeline.py --output results/object
```

## Status and limitations

0.3.0 adds the neural geometry object layer described above. 0.2.0 delivered the hull milestone:
a headless, parameterized geometry model with tests, CI, export, and validation that fails loudly.

What it does **not** provide: any shape family beyond a sphere for the neural objects (the
format and transaction layer are shape-agnostic, the fitted families are not), a field
accurate enough to replace a CAD surface (radial error is about 3 % of the radius), B-rep/STEP geometry, a real hull form (the breadth distribution is `sin(πx/L)^0.7 · √(z/draft)`, fore-aft symmetric, with a knuckle at the keel and no transom), fairness or surface-fairing checks, lightship weight, stability, resistance, or any classification-rule compliance. `metadata.json` carries that as `engineering_status`, and no number in this repository is reported without a check that produced it.
