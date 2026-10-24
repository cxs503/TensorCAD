# Changelog

All notable changes to TensorCAD are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
semantic versioning.

## 0.7.0 — 2026-10-24

Mesh conditioning and the two solver exits. The 0.6.0 release closed with "mesh
smoothing/decimation: 1.5–3.2 % of triangles are thinner than 2°; repairing them is the
next piece of work"; this release repairs them, measures the shape drift that repair
causes, and hands the result to TensorFEM (and to any tool that reads gmsh or STEP)
without adding a mandatory dependency.

### Added
- **`tensorcad.mesh_conditioning`** — `taubin_smooth` (λ 0.5 then μ −0.53, a band-pass so
  the body does not shrink), `project_to_field` (one Newton step per vertex along the
  field gradient, with a per-vertex step limit and a 4-step backtracking line search),
  `decimate` (quadric-error collapse under the link, orientation, duplicate-face and
  aspect-ratio guards, with every rejected collapse counted), `measure_mesh` and
  `condition_mesh`, which refuses to return a surface whose watertightness, manifolding,
  Euler characteristic, volume (0.5 %) or area (1 %) moved beyond its declared tolerance.
- **Measured on AFF-8 at `(150, 18, 26)`** (65,812 triangles, 1.664 % slivers):
  Taubin ×30 → min angle 0.0080° → 8.3202°, no slivers, volume −0.127 % of the analytic
  0.70754064 m³ in 0.4 s; + projection ×30 (53.6 s) → 3.6590°, aspect 17.6, volume
  −0.117 %, area −0.198 %, worst displacement 18.4 mm, worst residual `|f|` 0.57 mm;
  + decimation to 45,000 / 30,000 / 20,000 faces in 5.4 / 9.7 / 11.6 s → min angle
  2.70 / 2.56 / 2.40°, zero slivers, volume drift −0.12 / −0.14 / −0.21 %.
- **The decimate-first order is reported, not recommended.** `condition_mesh` at the same
  budgets reaches 2.21° (45 k) and 0.149° with 28 slivers (30 k); at 45 k faces it costs
  40 s against 5 s smoothing-first. Both pass the shape gates; the example prints both
  tables.
- **`tensorcad.adapters.fem`** — solver-neutral handoff needing no extra: node/element
  tables, boundary patches assigned by nearest analytic component, gmsh 2.2 ASCII with one
  physical tag per patch, and a JSON manifest. Round trip: worst node displacement **0.0 m**,
  worst per-patch area difference 2.2e−14 m², digest identical to the source surface.
- **`tensorcad.adapters.ocp`** — optional STEP bridge (`pip install "tensorcad[ocp]"`):
  sew → shell → solid → STEP, with the free-edge count, shell closure, `BRepCheck_Analyzer`
  verdict, kernel volume, written size and re-imported volume all reported. Measured on
  14,764 AFF-8 triangles: 12.8 s sewing, 4.4 s writing, 0 free edges, kernel volume
  0.694708267 m³ against a mesh signed volume of 0.694708267 m³ (**+5.3e−15 relative**),
  re-import −3.1e−13 relative, 37.14 MB. Detached islands are kept and counted
  (`shells`/`solids`) rather than dropped.

### Notes
- The STEP file is a **faceted solid**: one planar face per triangle (14,764 faces). That
  is what a solver-facing exchange file needs; a B-rep feature model is issue #4.
  `write_step` raises `OcpUnavailable` with the pip command when the extra is absent — it
  never writes a mesh file with a `.step` extension.
- 357 tests in the default lane (336 without the `neural` extra), plus 10 OCP tests that
  run when the wheels are installed; 93.2 % line coverage overall, with
  `tensorcad/adapters/ocp.py` at 32 % in the default lane because its tests skip there.
  Strict `mypy` clean on 22 files (OCP is an untyped import behind a documented
  override), `ruff` clean.

### Still missing
- Mesh *refinement* or remeshing to a target size field: conditioning only coarsens.
- A B-rep/STEP feature model, booleans, fillets, sweeps.

## 0.6.0 — 2026-10-24

The DARPA SUBOFF appendages (AFF-8): the sail and the four cruciform stern fins, on top
of the 0.5.0 bare hull. Both are published closed-form predicates, so they are
implemented as formulas — no boolean engine, no CAD kernel, no CSG tree. The three
solids are unioned as `min` of their fields, which is exact for the sign and the zero
set. `examples/suboff_appendages/run_example.py` enforces every claim below;
`docs/suboff_appendages.md` records the numbers.

### Added
- **`tensorcad.suboff_appendages.SuboffSailField`** — the published three-segment sail:
  entrance polynomial, parallel middle, exit polynomial, plus the elliptical cap
  (`Zmax = 0.109375 ft` half width, body top `1.507813 ft`, axial
  `[3.032986, 4.241319] ft` of the 14.291667 ft hull). The section is one radial
  coordinate `ρ = |y|` below the body top and `√(y² + (2(z − YTMP))²)` above it, so the
  cap and the body need no union. `height` (the cap apex, 0.47625 m) is a property.
- **`tensorcad.suboff_appendages.SuboffFinsField`** — the swept NACA 4-digit cruciform:
  chord `c(r) = −0.466308 r + 0.88859` ft, thickness
  `(0.2969, 0.126, 0.3516, 0.2852, 0.1045)`, span `r ∈ [0.075, 0.825] ft`, trailing edge
  at `x = 13.146284 ft`. The planform is a convex quadrilateral in `(x, r)`, whose
  distance is computed exactly (four edge distances), not as a gap.
- **`SuboffAppendagesField` and `SuboffFullField`** — the sail ∪ fins, and the whole
  AFF-8 model hull ∪ sail ∪ fins. `sail_scale` / `fin_scale` reproduce TensorLBM's
  similarity axes (sail about its footprint centre, centreplane and deck plane; fins
  about their trailing edge and root radius), so a scaled appendage stays attached.
- **Verified against the source twice.** The section values match TensorLBM's
  transcription exactly (`0.00e+00 ft` on the sail, `6.9e-18 ft` on the NACA section);
  containment matches the actual torch predicates on 18 hand-picked labels and on
  20,000 random points at each of three scales with **0 mismatches**, against a
  verbatim in-repo transcription of those predicates that CI can run.
- **An independent volume reference.** Quadrature over the published sections gives
  sail 0.00916043 m³ (46.2 % of it outside the hull), fins 0.00785718 m³, and the full
  model 0.70754064 m³ (+1.191 % over the bare hull). The exported
  `(436, 52, 76)` mesh lands at 0.707120 m³ (−0.059 %) with 569,056 triangles,
  6.350646 m² of surface, watertight, edge-manifold, genus 0.
- **The magnitude is measured, not claimed.** The cap and the thickness surface have no
  closed-form nearest point, so the fields declare `quantity = "implicit_function"` and
  report the first-order normalised gap. Along the analytic normals: sail 0.0012 mm
  median error at a 0.5 mm offset, fins 0.00001 mm; `|∇f|` is 0.996 (sail) and 1.000
  (fins) at the median, and exact on the flat faces.

### Notes
- The fins alone are a **genus-1** closed shell (`V − E + F = 0`) and that is the
  published geometry, not a meshing bug: the blades start at `r = 0.075 ft` and touch
  near the root, leaving the axis column empty. The hull fills it, so the full model is
  genus 0. Both are asserted.
- Throughput: 77 k queries/s on a 20,000-point batch for the full model (~1.6× the
  bare hull's cost).
- 314 tests (295 without the `neural` extra), 95.9 % line coverage; `suboff_appendages`
  99 %. Strict `mypy` clean on 19 files, `ruff` clean.

### Still missing
- Mesh smoothing/decimation: 1.5–3.2 % of the triangles are thinner than 2°, on the
  NACA knife edges. Reported as warnings; repairing them is the next piece of work.
- The AFF-2/3/4/7 appendage layouts, explicit intersection curves, and therefore STEP.

## 0.5.0 — 2026-10-24

The first real hull form on the geometry core: the DARPA SUBOFF bare hull (AFF-1) as an
analytic signed-distance field, with its offset-table provenance, an analytic reference,
anisotropic meshing, and a measured export. `examples/suboff/run_example.py` enforces
every claim below; `docs/suboff.md` records the numbers.

### Added
- **`tensorcad.suboff.SuboffHullField`** — the four-segment DARPA offset table (bow
  polynomial, parallel midbody, sextic stern taper, elliptical stern cap) in its
  published feet, carried as a true signed distance in metres. `L = 4.356100 m`,
  `R = 0.254000 m`, `L/D = 8.575`, bow at `x = 0`, stern at `x = length`.
  `|s|` matches a brute-force 400,001-sample meridian to `2.4e-08 m` (0.09 ppm of R);
  the gradient is the exact unit normal (`1.000000`, agreeing with central differences
  to `1.7e-08`).
- **Provenance in the record, not in a document.** `parameters()` / `provenance()`
  carry the source (Groves, Huang & Chang 1989, DTRC/SHD-1298-01), the transcription
  (TensorLBM `suboff_cad.suboff_radius_profile`), the unit convention, the axial origin
  and the profile description. The profile agrees with that transcription to
  `0.00e+00` at every station tested.
- **`suboff_reference_geometry()`** — analytic volume, wetted area, midship area and
  prismatic coefficient by quadrature, reported at two resolutions with their
  difference (volume 0.000 ppm, area 0.005 ppm between 10⁵ and 10⁶ samples).
  `V = 0.699211 m³`, `A = 5.988264 m²`, `Cp = 0.79194`.
- **Anisotropic export.** `to_mesh` and `extract_surface` now accept three cell counts
  per axis, and `sample_grid` marches a non-cubic grid. A 8.575:1 hull at
  `(218, 26, 26)` (147,368 cells) reaches −0.138 % volume / −0.066 % wetted-area
  error with 4.3 % slivers, where an isotropic `64³` (262,144 cells) reaches
  −0.089 % / −0.080 % with 10.0 % slivers and 64 stations along the hull.

### Fixed
- **`test_public_surface` poisoned the interpreter.** It deleted `numpy` from
  `sys.modules` and reloaded `tensorcad` to prove the core needs no third-party stack,
  which left torch's numpy bridge broken for every later test in the session — a
  real order-dependent failure, observed when the file was run before
  `test_geometry_field.py`. The check now runs in a fresh interpreter, like the
  neighbouring opt-in tests.
- **A scale-dependent "degenerate face" criterion.** 0.4.0 counted faces with an area
  below 10⁻¹² m², which is a different geometric claim on a 0.3 m cylinder and a
  4.4 m hull. Tests and benchmarks now use area below 10⁻⁶ of the longest edge
  squared. The corrected 0.4.0 record: the cylinder at 40 cells emitted 8,108 faces
  below 10⁻¹² m² (12.4 %), 2,936 of them with an aspect ratio above 10⁶ — not, as
  0.4.0's notes said, "exactly zero area".

### Notes
- 281 tests (262 without the `neural` extra), 95 % line coverage.
- TensorLBM's `r_over_l` default is `1/(2·8.57) = 0.058343` (`R = 0.254148 m`); the
  offset table gives `0.058310` (`R = 0.254000 m`). This module uses the table, which
  makes its wetted area 0.060 % smaller than TensorLBM's published value. At
  TensorLBM's radius the two implementations agree to `+0.00000 %`, so the whole
  difference is the radius convention.
- Mesh quality is reported, not repaired: at `(218, 26, 26)` 4.3 % of triangles are
  thinner than 2°, concentrated at the two tips where the meridian slope is infinite.
  There is no smoothing or decimation yet.

## 0.4.0 — 2026-10-24

The continuous geometry core (CGC v0.1): analytic and learned geometry now answer the
same questions, so they can be interchanged, compared and handed to a solver without
each consumer learning a second API. Everything here is measured by
`examples/geometry_core/run_benchmark.py`, which fails the build if a claim stops
being true; see `docs/geometry-core.md` for the protocol and the numbers.

### Added
- **One protocol for every representation.** `tensorcad.geometry_field.GeometryField`
  defines batch `signed_distance`, `gradient`, `bounds`, editable `state`/`with_state`,
  `parameters`, `provenance`, `metadata`, `contains`, `numerical_gradient` and an
  optional `numpy_callable` fast path. Declared per field, not assumed: `kind`,
  `quantity`, `inside_sign`, `surface_isovalue`, `units`, `coordinate_frame`.
- **Four exact primitives** — `PlaneField`, `SphereField`, `CylinderField` (capped,
  arbitrary axis), `EllipsoidField`. Measured against independently recomputed closed
  forms, the difference is `0.00e+00` in double precision. Gradients are analytic and
  match central differences to `2.7e-11 … 2.0e-10`. `EllipsoidField` declares
  `quantity="implicit_function"`: its zero set is exact, its magnitude is not a distance.
- **Edits that provably do not drift.** `with_state` returns a new field and leaves the
  original untouched. Measured: sphere radius 0.5 → 0.65 moves the surface 0.150000 m in
  all 24 probe directions with the centre unchanged; a cylinder height edit leaves the
  barrel value bit-identical; an ellipsoid axis edit leaves the other two axes alone;
  unknown controls raise `FieldError` instead of being ignored.
- **`to_mesh(field, resolution)`** — field → closed mesh → the existing hull acceptance
  pipeline. At 40 cells per axis: sphere −0.129 %, capped cylinder −0.167 %, ellipsoid
  −0.134 % volume error; all watertight, edge-manifold, `V − E + F = 2`, no hard errors.
- **`inspect_mesh_quality`** — smallest interior angle, worst aspect ratio and sliver
  share as *warnings*: a closed solid can still be a poor discretisation.
- **`SampledField.handoff()` and `mesh_digest()`** — a self-describing record (schema,
  mesh SHA-256 over float64 bits, units, frame, counts, topology, volume, area, bounds,
  validation, quality) for the TensorFEM handoff.
- **`surface_point(field, direction)`** — ray-bisected surface position, so
  representations are compared by where their surfaces are rather than by field values
  that are not comparable.
- **`NeuralGeometryField` / `field_from_object`** — a published `.ngo` revision behind
  the same protocol, with `quantity="truncated_signed_distance"`, autograd gradients and
  conditioned edits. Measured on the committed bundle: surface 3.98 mm mean / 9.33 mm
  worst from the analytic sphere away from the dent; dent axis 441.73 mm against a
  440.00 mm target; conditioned radius 0.45/0.55/0.60 reproduced within 0.33/0.08/0.01 %.

### Fixed
- **Squashed triangles from grid-vertex coincidence.** A radius-0.3 cylinder in a ±0.3
  box is exactly 20 cells across, so the level set passed through grid vertices and
  marching tetrahedra emitted 8,108 faces with an area below 10⁻¹² m² (12.4 % of
  65,144), 2,936 of them with an aspect ratio above 10⁶, and a 22.4 % sliver share.
  `to_mesh` now offsets the sampling box by an irrational fraction of a cell per axis:
  every count → 0, sliver share → 2.3 %, volume error unchanged at −0.17 %.
- **Finite differences are the wrong derivative for a float32 field.** At `step=1e-6`
  a central difference of the network is dominated by rounding (0.023 disagreement with
  autograd); at `1e-4` it agrees to `3e-4`. `NeuralGeometryField.gradient` now uses
  torch autograd, and `numerical_gradient`'s docstring states the step requirement.

### Notes
- The learned field's `|∇f|` ranges 0.65–1.39, i.e. it is not a metric: distance queries
  must read the field value and surface queries the zero set. Beyond the truncation band
  the outward cosine drops from 0.983 to 0.721 — the truncation caveat, measured.
- The bundle's recorded `dent_peak_measured_m = 0.052690` is parent-relative; re-measuring
  against the analytic sphere gives 0.058268 m. Both are correct and both are printed
  side by side; the difference is the parent's own 5.58 mm axis error.
- 262 tests (245 without the `neural` extra), 95 % line / 91 % branch coverage.

## 0.3.1 — 2026-10-24

The review of 0.3.0 found four defects that could crash a legal model, corrupt an
active bundle or delete an export. They are fixed here, with regression tests that
fail on the old behaviour; nothing was re-thresholded to make a gate pass.

### Fixed
- **Convergence estimate crashed on a legal model.** `estimate_discretisation_error`
  guarded the requested grid (`nx*nz`), then built the *refined* parameters, which
  could exceed the `4001 x 513` limits (`nx=2002, nz=3` raised
  `nx must be <= 4001 stations, got 4003`). It now reports `skipped` with the
  refined grid in the reason. The refined parameter set also keeps
  `water_density` instead of silently falling back to 1025.
- **Half-written export on failure.** `export_hull` wrote the OBJ and CSV before
  the metadata, so a failure during reporting left a directory with a mesh and no
  record of the parameters that produced it. The metadata document is now built
  first and all three files are written afterwards.
- **`tensorcad-hull` had a half-open error boundary.** Only `generate_hull` was
  inside a `try`; export and metadata failures escaped as tracebacks. Every phase
  now reports `stage` + `error` (JSON on stdout under `--json`, human text on
  stderr) and exits 2. The CLI also computes the metadata document once instead of
  generating the refined grid twice.
- **`tensorcad-neural verify` required torch.** The `build_field` import sat
  outside the error boundary, so a torch-free environment got a `ModuleNotFoundError`
  traceback even though the manifest, hash and history checks need no neural
  extra. Those checks now run everywhere, `field_query` reports `skipped` when the
  extra is absent, and `verify` gained derived-artifact integrity checks
  (presence, SHA-256, size) as warnings.
- **Derived meshes could be written into the live bundle before publishing.**
  `apply_dent_edit` extracted straight into `bundle/derived/` and then committed,
  so a rejected or failed commit left files behind that no manifest described.
  Extraction now happens in a scratch directory and travels inside the staged
  bundle.
- **Extraction was invisible and destructible.** `extract_object` wrote
  `derived/surface.obj|stl` without registering them, and the next commit -- which
  swaps the whole bundle directory -- deleted them silently. Extraction into the
  object now publishes its own revision (`lineage.operation = "extract"`, weights
  untouched) and registers the records; an edit carries the parent's registered
  artifacts forward, marked with `source_revision`. `extract --output DIR` remains
  an unregistered export, which is now stated in the docstring.
- **`parse_manifest` did not validate `derived[]`.** Readers accepted an escaping
  `../../../../etc/passwd` path, a non-hex digest or a negative size, so `info`
  would print them. Records are now validated exactly like the writer builds them.
  `ObjectStore.history()` also rejects `NaN`/`Infinity` instead of accepting the
  tokens this project refuses to write.
- **`max_section_area` was quadratic on non-grid meshes.** It re-scanned every
  vertex per station, which is invisible on a lofted hull and ruinous on an
  extracted implicit surface: 14.1 s of a 14.1 s `compute_hydrostatics` at 64^3.
  One bucketing pass makes it `O(V + stations)`; the same call now takes 0.34 s and
  returns the identical numbers (a test pins it against the old definition).

### Added
- Coarse-grid warning: `generate_hull` compares the mesh-integrated volume against
  the analytic volume of the generating functions and warns above 2 % (`nx=3, nz=3`
  is -35.7 % and used to report `warnings == []`). `analytic_volume` is public, and
  the benchmark example now imports it instead of duplicating the formula.
- `tensorcad/py.typed`, so the strict typing this package is developed under is
  actually visible to PEP 561 consumers.
- Regression tests for every item above, including a parametrised overflow matrix
  and a "derived files survive an edit" test.

### Changed
- Test suite: 194 tests without the neural extra, 213 with it (was 175 / 189);
  96 % branch coverage.
- Corpus bundles: `verify` reports recorded-but-absent derived meshes as warnings.
  The committed experiment bundle intentionally omits the 20 MB exports, so it now
  shows two warnings -- which is the honest description of what is in the repo.

## 0.3.0 — 2026-10-24

Geometry can now be stored as an object rather than a mesh: a learned continuous
field plus the metadata, versioning and transactions needed to trust it. STL/OBJ
become derived exports instead of the source of truth.

### Added
- `tensorcad.neural_geometry` (standard library only) — the NGO v0.2 bundle format:
  manifest with lineage, field semantics, state, training, validation, provenance
  and derived-artifact records; SHA-256 verification; path-escape rejection;
  `history.json`; `ObjectStore` with stage → validate → commit/rollback, atomic
  directory swaps, optimistic concurrency (`expected_parent`), and JSON-safe
  history records for failed attempts.
- `tensorcad.neural` (optional `neural` extra) — `ConditionalTSDF` field runtime,
  base training, replay-protected local edits, measurement of dent peak / dent-region
  error / protected-region mean, P95, P99, max, and a marching-tetrahedra
  `extract_surface` that emits an already-outward-oriented, watertight mesh.
- `tensorcad.cli_neural` / `tensorcad-neural` — `fit`, `edit`, `extract`, `info`,
  `verify`; every command prints measured numbers, supports `--json`, and exits
  non-zero when a gate fails (a rolled-back edit exits 1, not 0).
- `write_stl` (binary STL) for derived meshes, alongside the existing OBJ writer.
- `examples/neural_geometry/run_object_pipeline.py` — the end-to-end benchmark:
  fit → edit → extract → verify → deliberate rejection, with assertions.
- `schemas/neural-geometry-manifest-v0.2.schema.json` and
  `docs/neural-geometry-format.md` describing what the format stores, what it
  replaces, and what it does not.
- Coverage gate (`fail_under = 90`, currently 96 %) and a torch-free `neural`
  CLI help check in the default test job.

### Changed
- Manifest v0.2 keeps the 0.1 top-level file keys, so older bundles still load;
  unknown fields are preserved rather than dropped.
- Protected-region gates now include P99 and max, not only the mean, and failed
  edits report `failed_checks` instead of a single generic reason.
- The dent metrics are split into `dent_peak_*` (on the axis) and
  `dent_region_error_mean_m` (inside the cone), so a peak target is never compared
  against a region average.
- The dent cone is always sampled with a deterministic stratified grid, so the
  dent metrics do not silently become infinite when few random rays land in it.
- Torch↔NumPy conversion falls back to list conversion when the installed wheels
  have mismatched ABIs (torch 2.2 with NumPy 2.x) rather than failing a run.
- Test suite: 189 tests (from 142), 96 % branch coverage on the measured modules.

### Known limits
- The only fitted shape family is the sphere (with radial state and local dents);
  ellipsoids, cylinders and hulls are not wired into the object layer yet.
- Field accuracy is about 3 % of the radius at 0.5 m; extraction at 96³ takes
  several seconds on CPU.

## 0.2.0 — 2026-10-24

The exported surface becomes a solid a solver can actually use, and the checks
that matter are now enforced rather than documented.

### Fixed
- **Inverted winding.** The port half of the hull was lofted with the same face
  order as the starboard half, so its normals pointed inward and the mesh's
  signed volume cancelled to exactly zero. Solved structurally: one closed
  section ring per station, plus a generic `orient_outward` pass, instead of two
  mirrored sheets.
- **Tolerance welding artefacts.** Merging the two mirrored sheets by a floating
  point tolerance produced coincident duplicates along the keel and one
  non-manifold edge. Vertices are now shared by exact position at build time.
- **Open, unclosed hull.** The surface ended at the waterline with no deck face;
  it is now capped, closed, edge-manifold, and has Euler characteristic V-E+F = 2.
- **Non-exact waterline plane.** `z = -draft + draft * (nz-1)/(nz-1)` left the
  waterline vertices about 1e-16 off the plane, breaking coplanarity assumptions;
  the mapping is now `z = -draft * (1 - fraction)`, which is exact.
- **Public API.** `export_hull` was not importable from `tensorcad` although it was
  the documented entry point; it is re-exported with the rest of the API.
- **OBJ precision silently altered the model.** Files were written with 9
  significant digits; export now uses 17 so `read_obj(write_obj(mesh))` returns the
  identical vertex table (a test locks this in).

### Added
- `tensorcad.geometry` / `tensorcad.mesh` / `tensorcad.parameters` /
  `tensorcad.validation` / `tensorcad.export` / `tensorcad.hydrostatics` /
  `tensorcad.adapters` — the module boundaries the README promised in 0.1 are now
  real, with generic mesh algorithms kept out of the ship-specific model.
- Validity checks that a geometry foundation has to provide: watertightness, edge
  and vertex manifoldness, outward orientation, zero-area faces, Euler
  characteristic, volume-bounded-by-bbox, and realised-versus-requested L/B/T.
  `generate_hull` raises instead of returning a model that fails them.
- Hydrostatics integrated from the mesh with the divergence theorem: displaced
  volume, displacement mass, wetted surface, waterplane area, LCB, LCF, VCB, Cb,
  Cp, and maximum section area — all reported in `metadata.json`.
- Measured discretisation evidence: `estimate_discretisation_error` regenerates on
  a refined grid and records the relative change (default grid: 0.88 % under the
  analytic volume, 0.14 % at nx=161/nz=41).
- Dimensionless ratio checks (`L/B`, `B/T`, `T/L`) reported as warnings, and hard
  bounds on `nx`/`nz` so an accidental `1e8` fails fast instead of allocating.
- Explicit `nx` parity reporting: with an even station count no station lands on
  midship, so the realised beam falls short; that is now a named warning instead of
  a silent ~1 % dimensional error.
- CLI: `--density`, `--json`, `--no-export`, `--no-convergence`, non-zero exit
  codes with the full validation list on stderr.
- `examples/parametric_hull/run_example.py`: a runnable benchmark that verifies
  closure, orientation, dimensional fidelity, hydrostatics against the analytic
  value, and a resolution study — and exits non-zero if a claim fails.
- `LICENSE` (Apache-2.0 text was declared in `pyproject.toml` but absent),
  `CHANGELOG.md`, `.gitignore`, ruff + mypy (strict) + coverage (90 % gate, currently 99 %)
  configuration, and a Python 3.13 CI matrix entry with a byte-identical
  reproducibility check.
- 129 tests, up from 8.

### Changed
- `HullModel` holds a `Mesh` (`vertices`/`faces` remain as accessors), plus
  `validation` and `hydrostatics`; `metadata()` no longer re-runs validation on every
  call.
- Face indices in the API are 0-based; the OBJ writer applies the 1-based offset.
- `validate_hull(model)` keeps its signature and returns the list of error strings.
- `metadata.json` schema version is `2.0`.

## 0.1.0 — 2026-10-23

- Initial package: deterministic parametric sampled hull surface, OBJ/CSV/JSON
  export, CLI entry point, structural validation, and CI.
