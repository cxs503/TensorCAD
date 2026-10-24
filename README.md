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

## Continuous geometry core (0.4)

Every representation now answers the same questions. `tensorcad.geometry_field`
defines one protocol — batch signed distance, gradient, bounds, editable state,
parameters, sign convention, units, frame — and ships four exact primitives as the
reference the others are measured against:

```python
from tensorcad import SphereField, CylinderField, to_mesh

sphere = SphereField(radius=0.5)
sphere.signed_distance([(0.9, 0.0, 0.0)])        # 0.4  (exactly |x|-r, not an approximation)
wider = sphere.with_state(radius=0.65)           # a new field; the original is untouched
wider.center == sphere.center                    # True -- the edit moved the radius only

sampled = to_mesh(sphere, 40)                    # field -> closed mesh -> validation
sampled.report.watertight, sampled.report.euler_characteristic     # True, 2
sampled.report.signed_volume / (4 / 3 * 3.14159 * 0.125) - 1       # -0.0013  (-0.13 %)
sampled.handoff()["mesh_sha256"]                 # what a solver is actually being given
```

The primitives (plane, sphere, capped cylinder, ellipsoid) are exact to `0.00e+00` in
double precision, so a learned field can be compared against them without wondering
which side of the comparison is wrong. An ellipsoid declares
`quantity="implicit_function"` and a learned field `"truncated_signed_distance"`:
the magnitude of neither is a distance, and the protocol says so instead of letting a
caller assume otherwise.

A published `.ngo` object is a `GeometryField` too:

```python
from tensorcad.neural.objects import field_from_object

learned = field_from_object("docs/experiments/2026-10-24-neural-geometry-object/bundle")
learned.quantity                                 # 'truncated_signed_distance'
learned.with_state(radius=0.55)                  # conditioned edit, no retraining
```

Measured through this protocol, the reference bundle's surface sits 3.98 mm (mean) /
9.33 mm (worst) from the analytic sphere away from its dent, and a conditioned radius
of 0.45/0.55/0.60 m is reproduced within 0.33/0.08/0.01 %. Details, including the
zero-area-triangle defect this work found in the sampler, are in
[docs/geometry-core.md](docs/geometry-core.md).

## SUBOFF bare hull (0.5)

The first real hull form on the geometry core. `SuboffHullField` implements the DARPA
SUBOFF bare hull (AFF-1) from its published offset table — four segments in feet,
carried as a true signed distance in metres:

```python
from tensorcad import SuboffHullField, to_mesh

hull = SuboffHullField()                          # L = 4.3561 m, D = 0.5080 m, L/D = 8.575
hull.length, hull.max_radius                      # 4.35610, 0.25400   (offset table, not a rounded L/D)
hull.signed_distance([(2.0, 0.30, 0.0)])          # +0.0460  (outside the barrel)
hull.with_state(radius=1.2 * hull.max_radius)     # a fatter hull, same length -- no drift

sampled = to_mesh(hull, (218, 26, 26))            # cells per axis: a slender body is not a cube
sampled.report.watertight, sampled.report.euler_characteristic      # True, 2
sampled.report.signed_volume / 0.699211 - 1       # -0.00138  (-0.138 % vs the analytic volume)
sampled.report.surface_area / 5.988264 - 1        # -0.00066  (-0.066 % vs the analytic wetted area)
```

Measured, with the details in [docs/suboff.md](docs/suboff.md):

| Check | Result |
|---|---|
| `|s|` vs a brute-force meridian | `2.4e-08 m` (0.09 ppm of R) |
| joints | C0 exact, C1 continuous at all three; curvature jumps at the stern cap |
| gradient | unit length, matches central differences to `1.7e-08` |
| profile vs TensorLBM's transcription | `0.00e+00` at every station |
| wetted area at TensorLBM's radius | reproduces its 5.991830 m² to `+0.00000 %` |
| offset-table wetted area vs the experiment | −0.162 % of 5.998 m² |
| export `(218,26,26)` | watertight, genus 0, volume −0.138 %, area −0.066 % |

The tips are points of zero area, so the mesh stops short of them by about a cell and
4.3 % of the triangles there are thinner than 2°; that is reported as a warning rather
than smoothed away. `python examples/suboff/run_example.py` re-measures every row and
fails the build if one stops holding.

## SUBOFF appendages (0.6)

The sail and the four stern fins (AFF-8) on top of that hull — from the published
predicates, so no boolean engine and no CAD kernel is involved:

```python
from tensorcad import SuboffFullField, SuboffSailField, SuboffFinsField, to_mesh

full = SuboffFullField()                       # hull + sail + cruciform fins, L = 4.3561 m
full.contains([(1.10, 0.0, 0.35)])             # [True]  -- inside the sail, above the hull
full.signed_distance([(1.10, 0.0, 0.35)])      # -0.0333 m (a lower bound: min-union)
full.with_state(sail_scale=1.5, fin_scale=1.2) # bigger appendages, still on the hull

sail = SuboffSailField()                       # the sail alone: 3 polynomials + elliptical cap
sail.height                                    # 0.47625 m above the axis (0.22225 above the deck)
SuboffFinsField().parameters()["chord_root_m"] # 0.26018 m at the root, 0.15358 m at the tip

sampled = to_mesh(full, (436, 52, 76))         # anisotropic: 8.6:1 model, 0.033 m sail
sampled.report.signed_volume / 0.70754064 - 1  # -0.00059  (-0.059 % vs quadrature)
```

Measured, with the details in
[docs/suboff_appendages.md](docs/suboff_appendages.md):

| Check | Result |
|---|---|
| sail half width / NACA half thickness vs TensorLBM's transcription | `0.00e+00 ft` / `6.9e-18 ft` |
| containment vs TensorLBM's torch predicates | 18/18 hand-picked labels; 20,000 random points × 3 scales with **0 mismatches** (checked externally at 400,000 per scale, also 0) |
| normalised gap, 0.5 mm off the analytic surface | sail p50 0.0012 mm, fins p50 0.00001 mm |
| `\|∇f\|` | sail 0.980 … 0.996 … 1.000, fins 0.999 … 1.000 … 1.000 |
| volume from quadrature over the sections | sail 0.00916043, fins 0.00785718, full 0.70754064 m³ (+1.19 % over the hull) |
| export `(436,52,76)` | watertight, genus 0, volume −0.059 %, 569,056 triangles, 6.350646 m² |
| fins alone | watertight genus **1** (the blades touch near the root; the hull fills the axis) |
| throughput | 77 k queries/s on a 20,000-point batch |

## Mesh conditioning (0.7)

Marching tets gives a closed surface and a bad triangulation. `tensorcad.mesh_conditioning`
fixes the triangulation without moving the shape: Taubin smoothing (a band-pass, so
unlike a plain Laplacian it does not shrink the body), then a bounded Newton projection
back onto the analytic field, then quadric-error decimation with the link, orientation,
duplicate-face and aspect-ratio guards.

```python
from tensorcad import SuboffFullField, condition_mesh, measure_mesh, to_mesh

mesh = to_mesh(SuboffFullField(), (150, 18, 26), margin=2).mesh
measure_mesh(mesh).sliver_share                # 0.0166 -- 1.66 % of triangles are slivers
result = condition_mesh(mesh, field=SuboffFullField(), target_faces=45_000)
result.passed                                  # topology and shape-drift gates
measure_mesh(result.mesh).sliver_share         # 0.0
result.after.signed_volume / 0.70754064 - 1    # -0.00127  (0.13 % below the analytic volume)
```

Measured at `(150, 18, 26)`, with the details in
[docs/mesh_conditioning.md](docs/mesh_conditioning.md):

| Stage | tris | min angle | max aspect | slivers | ΔV | ΔA | s |
|---|---|---|---|---|---|---|---|
| input | 65,812 | 0.0080° | 7,126 | 1.664 % | −0.323 % | −0.378 % | — |
| Taubin ×30 | 65,812 | 8.3202° | 8.3 | 0.000 % | −0.127 % | −0.729 % | 0.4 |
| + projection ×30 | 65,812 | 3.6590° | 17.6 | 0.000 % | −0.117 % | −0.198 % | 53.6 |
| + decimate → 45,000 | 45,000 | 2.7038° | 25.0 | 0.000 % | −0.118 % | −0.198 % | 5.4 |
| + decimate → 20,000 | 20,000 | 2.4021° | 25.0 | 0.000 % | −0.209 % | −0.162 % | 11.6 |

Volume/area drift is measured against the analytic AFF-8 references, not against the
input mesh. The minimum angle is reported, not gated: a knife edge's neighbourhood has no
good triangle to offer, and the example prints the decimate-first path's worse numbers
(0.15° at 30 k faces) rather than hiding them. The sliver share *is* gated, at zero,
because a sliver is a degenerate fold rather than a taste judgement.

## STEP and the TensorFEM handoff (0.7)

Two exits, for two consumers, and neither is imported by the core:

```python
from tensorcad.adapters.fem import fem_surface, write_gmsh, write_fem_json
from tensorcad.adapters.ocp import have_ocp, write_step

surface = fem_surface(sampled, components={"hull": SuboffHullField(),
                                           "sail": SuboffSailField(),
                                           "fins": SuboffFinsField()})
write_gmsh(surface, "aff8_surface.msh")     # gmsh 2.2, one physical tag per patch
write_fem_json(surface, "aff8_surface.json") # same tables + handoff manifest

have_ocp()                                   # False without the optional extra
record = write_step(mesh, "aff8.step", round_trip=True)
```

Measured by `examples/step_export/run_example.py`, details in
[docs/step_export.md](docs/step_export.md):

| Check | Result |
|---|---|
| gmsh 2.2 round trip (7,386 nodes) | worst node displacement **0.0 m**, worst per-patch area difference 2.2e−14 m² |
| JSON round trip | digest identical to the source surface |
| patches on the AFF-8 surface | hull 96.75 %, sail 1.72 %, fins 1.53 % of 6.105 m²; worst centroid-to-component distance 30–32 mm on a 60 mm cell |
| OCP sewing (14,764 triangles) | 12.8 s, 0 free edges, 0 multiple edges, all shells closed, `BRepCheck_Analyzer` valid |
| kernel volume vs mesh signed volume | 0.694708267 m³ vs 0.694708267 m³ (**+5.3e−15 relative**) |
| STEP round trip | re-imported volume −3.1e−13 relative; 37.14 MB, one planar face per triangle |
| disconnected islands | written as a compound with every island counted (`shells`/`solids`), never dropped |

The STEP output is a **faceted solid**, not a B-rep model: 14,764 planar faces, no
fillets, no analytic surfaces, no feature tree. That is what a solver-facing exchange file
needs and what a designer does not want; a feature model is [issue #4](https://github.com/cxs503/TensorCAD/issues/4).
The optional extra is `pip install "tensorcad[ocp]"`; without it `write_step` raises
`OcpUnavailable` with that command in the message, and the FEM handoff still runs,
because it needs nothing but the standard library.

## Design loop (0.8)

Everything above answers "what is this surface?". This answers the question an agent
actually receives — *"give me a hull that displaces 12.06 m³ with a B/T of 2.67 and a
wetted area under 33.3 m², and prove it"*:

```python
from tensorcad import DesignProblem, DesignTargets, HullParameters, write_solver_manifest

solution = DesignProblem(
    base=HullParameters(length=10.0, beam=2.0, draft=0.8),
    targets=DesignTargets(displacement_m3=12.063418, beam_to_draft=2.4 / 0.9,
                          wetted_area_m2_max=33.3),
    variables=("length", "beam", "draft"),
).solve()

solution.status            # 'converged'
solution.solved            # converged AND every geometric gate green
solution.failed_checks     # names, never a message to parse
solution.convergence_order # 1.940 -- measured, above the finite-difference floor
write_solver_manifest(solution, "results/solver_manifest.json")
```

Damped Newton on `log` of the free dimensions, finite-difference Jacobian, backtracking
line search, bounds from the hull model's own specs, every iteration kept and measured.
Standard library only, like the hull model it drives.

| Check | Measurement |
|---|---|
| recover a known hull (L 12, B 2.4, T 0.9) from its measured displacement + ratios | all three dimensions to **1e-7 relative** in 6 evaluations, 1.1 s |
| convergence order, from the residual triple (9.1e-3, 4.1e-5, 8.3e-10) | **1.940** (Newton is 2; the estimate stops above the 1e-10 finite-difference floor) |
| Jacobian sanity | B/T has sensitivity **+1.000000 / −1.000000**, displacement **+1.000000** in beam and draft |
| replay | same targets → same mesh digest `814f4770b0189c58` |
| re-open with a 15 % heavier brief | converges from the previous answer in 5 evaluations, 13.872930 m³ |
| grid sensitivity | (41,13) → (81,25): displacement **+0.577 %**, area +0.18 % |
| a design that *converges and is refused* | 200 m³ at L = 10 m: ‖r‖ = 1.13e-12, and the verdict is `solved=False` naming `mesh_min_angle`, `mesh_aspect_ratio`, `mesh_sliver_share` |

The manifest carries the target, every iteration, the full `handoff/1` record, the
reference quantities a solver needs and `"solver": {"status": "ready", "executed_here":
false}`. TensorCAD measures and hands off geometry; running the flow or the structure is
TensorFEM's and TensorLBM's job, and no number here claims otherwise. Details:
[docs/design_loop.md](docs/design_loop.md).

## Representation benchmark (0.9)

Analytic formula, explicit triangles, a learned field and a kernel export, all through one
harness on the same rays, the same closed form and the same gates:

```bash
python examples/benchmarks/run_benchmark.py --quick --output results/bench
```

Case A is a sphere with a Gaussian dent of known depth; Case B is a design target, a converged
hull, and the gmsh artifact a solver reads. A row whose optional extra is missing reports
`not_supported` with the command that would measure it — never a zero.

| Case A, full sampling (512 rays, one edit, one reference) | analytic | mesh edit | neural field |
|---|---|---|---|
| cone-mean depth error | 0 (reference) | **2.3e-6 m** | 1.9e-3 m |
| dent surface error, mean | 0 | 4.5e-5 m | 2.1e-3 m |
| leak outside the edit (vs its own parent) | 0 | **0 m** | 1.3e-3 m (imported) |
| surface error over every ray (mean / P95) | 0 | 1.0e-4 / 1.8e-4 m | 3.4e-3 / 9.2e-3 m |
| wall clock for 512 rays (informational, see the doc) | 0.05 s | 3.1 s | 0.47 s |

The measured verdict is that a plain vertex move beats the trained field by 13–840× on depth
and surface accuracy for this task, with an exactly-zero leak outside its support; the field
keeps a 222 KB topology-free artifact and a continuous surface, which this benchmark does not
score, and answers rays faster only because the harness's own mesh ray caster has no spatial
index. Both rows face the **same 1 %-of-depth threshold**, so the neural row fails the mesh
row's gate (exit code 1) instead of being given a looser one. The 40×80 grid that happens to
pass is rejected in favour of 80×160 because its chord-sag floor (1.9e-4 m) must sit well
below the gate. Details, the resolution ladder and the threats to validity:
[docs/representation_benchmark.md](docs/representation_benchmark.md); the committed runs are in
`docs/benchmarks/{quick,full,neural}/`.

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
- `tensorcad.geometry_field` — the continuous-geometry protocol, four exact primitives,
  field sampling to a validated mesh, and the solver handoff record.
- `tensorcad.suboff` — the DARPA SUBOFF bare hull as an analytic field, with its
  offset-table provenance, analytic reference geometry and anisotropic export.
- `tensorcad.suboff_appendages` — the AFF-8 sail and cruciform stern fins from their
  published predicates, plus the hull-sail-fin union field.
- `tensorcad.design` — the design loop: targets in, damped Newton over the free hull
  dimensions, every iteration measured and verified, and a solver manifest out.
- `tensorcad.mesh` — generic algorithms: weld, orient outward, boundary/manifold/watertight analysis, surface integrals.
- `tensorcad.mesh_conditioning` — Taubin smoothing, field projection, quadric decimation,
  and the `condition_mesh` gates that refuse to hand on a surface whose shape drifted.
- `tensorcad.validation` — named geometric checks (including mesh-quality slivers and
  aspect ratios), errors versus warnings, `ValidationReport`.
- `tensorcad.hydrostatics` — mesh-integrated displacement, areas, centres, coefficients.
- `tensorcad.export` — OBJ/CSV/JSON writers and an OBJ reader (round-trip exact).
- `tensorcad.hull` — the one ship-specific model.
- `tensorcad.adapters` — solver/kernel bridges the core must never import: `adapters.fem`
  (node/element tables, gmsh 2.2 + JSON, no extra) and `adapters.ocp` (sew → solid → STEP,
  `ocp` extra).
- `tensorcad.neural_geometry` — neural geometry object format: metadata, hashing, atomic
  revision publishing, edit transactions (standard library only).
- `tensorcad.neural` — optional field runtime: training, local edits, marching-tetrahedra
  extraction to a closed mesh (`neural` extra, consumes torch/numpy).
- `tensorcad.cli_neural` — `tensorcad-neural` entry point for fit/edit/extract/info/verify.
- `tensorcad.benchmark` — the representation benchmark: metric schema with explicit missing
  states, ray caster, exact references, and `results.json` / `report.md` generation.
- `examples/` and `tests/` — runnable benchmarks and 394 tests, 15 of which need the
  `neural` extra; the OCP bridge adds 10 more when its wheels are installed.

The core package depends on nothing outside the standard library. That is a tested guarantee, not a preference: a B-rep kernel stays an optional backend (see issue #4), and `tensorcad.adapters.ocp` is the planned seam.

## Install and run

```bash
python -m pip install -e ".[test,dev]"
pytest --cov
python examples/parametric_hull/run_example.py --output results/hull-demo
python examples/geometry_core/run_benchmark.py --resolution 40
python examples/suboff/run_example.py --output results/suboff
python examples/design_loop/run_example.py --output results/design
python examples/mesh_conditioning/run_example.py --output results/conditioning
python examples/step_export/run_example.py --output results/step   # STEP part needs .[ocp]
tensorcad-hull --length 12 --beam 2.4 --draft 0.9 --nx 61 --nz 17 --output results/hull-demo
tensorcad-hull --json            # full reproducibility record on stdout

# neural geometry objects (needs the torch-backed extra)
python -m pip install -e ".[neural]"
python examples/neural_geometry/run_object_pipeline.py --output results/object
```

## Failure behaviour

The point of a headless geometry library is that a caller can tell "bad input" from
"bad output" from "I could not finish", so each phase has a defined answer:

- `generate_hull` raises on invalid parameters or a mesh that fails a check, and reports
  every problem at once;
- an under-sampled grid is a **warning** that names the measured error against the
  analytic volume (`nx=3, nz=3` is −35.7 %, the default 41×13 is 0.88 %);
- `tensorcad-hull` wraps every phase: `--json` prints `{"error": ..., "stage": ...}` on
  stdout, humans get one line on stderr, and the exit code is 2 — no tracebacks, and no
  half-written artifact directory (the metadata document is built before anything is
  written);
- `tensorcad-neural verify` separates integrity (manifest, hashes, history, derived
  artifacts — standard library only) from the live field query, which reports `skipped`
  when the `neural` extra is not installed;
- a rejected edit leaves the active revision byte-identical and appends the failed gates
  and their measurements to `history.json`.

## Status and limitations

0.9.0 adds the representation benchmark: one harness that puts an analytic formula, a mesh
edit, a learned field and an optional kernel export on the same rays, the same closed form
and the same thresholds, and publishes `results.json` + `report.md` for every run. 0.8.0 added the design loop: a displacement and shape-ratio target, a damped Newton search
over the free hull dimensions, a verdict that refuses a converged design whose triangles
fail validation, and a solver manifest with the geometry handoff attached. 0.7.0 adds mesh conditioning (smoothing, projection and decimation with measured shape
drift) and the two solver exits: a standard-library FEM handoff (gmsh 2.2 + JSON with named
boundary patches) and an optional OpenCASCADE STEP bridge behind the `ocp` extra. 0.6.0 adds the SUBOFF appendages (sail and cruciform fins) on top of the bare hull,
from their published closed-form predicates rather than a boolean engine. 0.5.0 added the SUBOFF bare hull on top of the geometry core. 0.4.0 added the continuous geometry core: one protocol for analytic and learned fields,
four exact primitives, mesh-quality checks, a solver handoff record, and a benchmark that
enforces every claim above. 0.3.1 fixed the defects the 0.3.0 review found (a legal model
could crash the convergence estimate, a failed commit could leave derived meshes in the
live bundle, and extraction could be deleted by the next edit). 0.3.0 added the neural
geometry object layer. 0.2.0 delivered the hull milestone: a headless, parameterized geometry model with
tests, CI, export, and validation that fails loudly.

What it does **not** provide: the AFF-2/3/4/7 appendage *layouts* (only AFF-8's sail and
fins are transcribed), a `Cb`-parameterised hull family (SUBOFF is the published hull,
editable in length and radius only), any shape family beyond a sphere for the neural objects (the
format and transaction layer are shape-agnostic, the fitted families are not), a field
accurate enough to replace a CAD surface (mean radial error 0.8 % of the radius, worst
1.9 %), booleans or sweeps, a hydrodynamic objective (the design loop matches displacement and shape
ratios, not resistance or stability), a B-rep/STEP *feature model* (the STEP export is a faceted solid),
mesh refinement or remeshing to a target size field (conditioning only coarsens), a real hull form (the breadth distribution is `sin(πx/L)^0.7 · √(z/draft)`, fore-aft symmetric, with a knuckle at the keel and no transom), fairness or surface-fairing checks, lightship weight, stability, resistance, or any classification-rule compliance. `metadata.json` carries that as `engineering_status`, and no number in this repository is reported without a check that produced it.
