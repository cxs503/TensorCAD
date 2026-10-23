# Changelog

All notable changes to TensorCAD are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
semantic versioning.

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
