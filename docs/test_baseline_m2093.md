# Test baseline — 2026-11-02 (m2093, TensorCAD 0.9.0 + gate/export work)

Environment: Python 3.11.6, numpy 2.4.6, torch 2.7.1, **no `safetensors`**, no
`cadquery-ocp`. Commands were run from a checkout (`python3 -m pytest -q
--continue-on-collection-errors`); the collection flag is required because the two
neural-extra modules cannot even be imported here (see "Red" below).

| run | tree | result | wall |
|---|---|---|---|
| baseline | clean `git archive HEAD` snapshot (no local changes) | **3 failed, 380 passed, 11 skipped, 2 errors** | 304 s |
| after this branch | worktree with the gate/export work | **3 failed, 411 passed, 11 skipped, 2 errors** | 309 s |

The failure/error/skip sets are identical, so the 31 new tests all pass and no
pre-existing result changed. `tests/test_gate.py` (12), `tests/test_resolution.py`
(6) and `tests/test_artifact_export.py` (13) are the additions.

## Red, honestly

| test | cause | introduced by this branch? |
|---|---|---|
| `test_geometry_field.py::test_the_learned_revision_speaks_the_same_protocol` | `ModuleNotFoundError: safetensors` (neural extra not installed) | no — fails identically on clean HEAD |
| `test_geometry_field.py::test_the_learned_field_gradient_is_autograd_and_points_outward` | same missing extra | no |
| `test_suboff_appendages.py::test_the_scalar_path_walks_the_sail_entrance_segment` | scalar and numpy code paths differ in the last bit (`0.00017595222258749225` vs `…227`); the test uses exact `==` | no — fails identically on clean HEAD |
| `ERROR tests/test_cli_neural.py`, `ERROR tests/test_neural_pipeline.py` | collection-time import of `safetensors` | no |

The first two and the two collection errors are an environment gap: installing the
`neural` extra (`torch`, `safetensors`, `numpy`) is enough. The suboff failure is a
genuine (pre-existing) numerical-strictness defect in that test, not touched here —
fixing it means deciding which of the two equivalent formulas is the reference, and
that is outside this task's scope. It is left red and reported.

## Slow (single-run wall clock, measured before the change)

| module | wall | note |
|---|---:|---|
| `test_benchmark.py` | 69 s | full benchmark sampling |
| `test_suboff_appendages.py` | 54 s | analytic field evaluation |
| `test_design.py` | 52 s | damped Newton loop, many hull generations |
| `test_suboff.py` | 46 s | offset tables |
| `test_geometry_field.py` | 23 s | Marching-tetrahedra extraction |
| whole suite | ~305 s | dominated by the four rows above |

The new suites add ~11 s in total: `test_artifact_export.py` pays for the
resolution search (three hull generations up to 81×33) and the STL round trip;
`test_resolution.py` regenerates grids up to 81×21.

## Skips

11 skips: 10 are the OCP bridge (wheels absent), 1 is the benchmark's
`ocp_step` row. Nothing is skipped silently to make the run green.

## Run it

```bash
python3 -m pytest -q --continue-on-collection-errors        # whole suite
python3 -m pytest -q tests/test_gate.py tests/test_resolution.py tests/test_artifact_export.py
python3 examples/hull_stl_export/run_example.py --output artifacts/hull_stl_export
```