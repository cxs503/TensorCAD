# PR #84 CI red — root cause and fix (mission m4091)

## What was red

`tests.yml` matrix (`test (3.10) … test (3.13)`) failed on PR #84 (branch `agent/m3858`,
head `3dd1df8`). Only `test (3.13)` reported `failure`; the other three reported
`cancelled`, which is the GitHub Actions matrix default (`fail-fast: true`): the first job
to fail cancels its siblings. So one job's failure explains all four red/grey cells.

The failing step is the **first** step of the job, `Lint` (`ruff check .`), which runs
before `mypy` and before `pytest`:

```
$ git checkout 3dd1df8 && ruff check .
Found 26 errors.
[*] 9 fixable with the --fix option (15 hidden fixes can be enabled with --unsafe-fixes)
```

`mypy` and the test suite were therefore never reached in that job — the other three
matrix legs were killed while still installing, hence `cancelled`.

## The 26 lint findings (all introduced by the #84 content, none in unrelated files)

| rule | count | where | fix |
|---|---|---|---|
| `B905` `zip()` without `strict=` | 8 | `tensorcad/righting_arm.py`, `examples/righting_arm/run_example.py`, `artifacts/m3858/make_figure.py` | `strict=True` where the sequences are equal-length by construction |
| `RUF007` prefer `itertools.pairwise` | 2 | same | replaced the `zip(xs, xs[1:])` idiom with `itertools.pairwise` |
| `B007` unused loop control variable | 3 | `_solve_waterline`, example, figure script | `for _` / drop the dead `for … : pass` loop |
| `F841` unused local | 2 | `righting_arm._triangle_terms` (`mean_x`, `mean_y`) | deleted (the moments actually used are `mean_z`, `mean2_*`) |
| `F401` unused import | 1 | `righting_arm.py` (`Face`) | removed the symbol from the `.geometry` import |
| `RUF059` unpacked but unused | 1 | `RightingArmCurve.to_dict` (`vanishing`) | `_, crossed = self._vanishing()` |
| `RUF022` / `RUF023` unsorted `__all__` / `__slots__` | 2 | `tensorcad/__init__.py`, `righting_arm.py` | sorted |
| `F541` f-string without placeholders | 1 | example | plain string |
| `C416` needless comprehension | 1 | example | `list(range(0, 91, 5))` |
| `W292` no newline at EOF | 5 | `righting_arm.py`, `tests/test_righting_arm.py`, example, 2 artifact scripts | trailing newline |

Ten of the fixes are mechanical (`ruff check --fix`: `W292`, `F401`, `F541`, `RUF022`,
`RUF023`); the other sixteen were applied by hand because they change control flow or
argue about intent.

## Behavioural safety of the non-mechanical edits

* `zip(..., strict=True)` in `_HeeledSection.integrals` — `self.faces` and `self.terms`
  are built in the same comprehension from the same face list, so their lengths are
  equal by construction; `strict=True` turns a hypothetical desynchronisation into a
  loud `ValueError` instead of a silent truncation.
* `_solve_waterline`: `for iterations in range(1, max_iterations + 1)` →
  `while iterations < max_iterations: iterations += 1`. Same iteration count, same
  bracket updates, same returned `iterations` (including the degenerate
  `max_iterations <= 0`, which returns the initialised `0`). Verified by the whole
  `tests/test_righting_arm.py` module.
* `pairwise(xs)` ≡ `zip(xs, xs[1:])` for sequences (both stop after `len-1` pairs), used in
  `area_under_curve`, `righting_arm_curve`'s monotonicity check and the example's
  refinement-monotonicity claim.
* `examples/righting_arm/run_example.py`: a `for row in results["sphere"]["convergence"]:
  pass` dead loop was deleted; it had no side effect.
* `artifacts/m3858/make_figure.py`: `enumerate(bars.items())` → `enumerate(bars.values())`
  (the dict key was unused; bar positions come from the enumeration index, unchanged).

No public API, signature, schema, number or docstring claim was touched.

## Verification (local, Python 3.11.6 — the runner image's 3.13 was not available offline)

| check | before (3dd1df8) | after |
|---|---|---|
| `ruff check .` | **26 errors** | **All checks passed!** |
| `mypy` (strict, 33 files) | clean | **clean** |
| `pytest tests/test_righting_arm.py` | n/a | **41 passed** |
| `pytest tests/test_stability.py tests/test_wigley.py` | n/a | **49 passed** |
| `pytest tests/test_transom.py tests/test_curvature.py tests/test_bonjean.py` | n/a | **158 passed** |
| `pytest tests/test_voxel.py tests/test_gate.py tests/test_analytic_hydrostatics.py tests/test_hydrostatics_analytic.py` | n/a | **100 passed** |
| `pytest tests/test_hydrostatics_inertia.py tests/test_asymmetric_hull.py tests/test_mirror_contract.py` | n/a | **71 passed** |
| `pytest tests/test_artifact_export.py tests/test_cli.py tests/test_hull.py` | n/a | **66 passed** |
| `pytest tests/test_neural_format_contract.py tests/test_suboff_appendages.py` | n/a | **67 passed** |
| `pytest tests/test_mesh.py tests/test_export.py tests/test_design.py tests/test_geometry_field.py tests/test_parameters.py` | n/a | **126 passed** |
| `pytest tests/test_mesh_conditioning.py tests/test_neural_geometry.py` | n/a | **35 passed** |
| `pytest tests/test_resolution.py tests/test_suboff.py tests/test_adapters_fem.py tests/test_adapters_ocp.py` | n/a | **46 passed, 10 skipped** |
| `pytest tests/test_benchmark.py` (33 tests, split into three selections) | n/a | **all passed** |

## Limitations / what this round did not do

* The 3.13 leg was not reproduced locally (no 3.13 interpreter offline); the failure is
  version-independent by construction (`ruff` is a static source check pinned by
  `pyproject.toml`, `target-version = "py310"`), so a 3.11 reproduction of the same step
  is conclusive for the root cause.
* The full suite was executed in selections because the sandbox caps a single command at
  ~120 s; the `quick_run`-fixture benchmark tests are the slowest single group
  (~110 s on their own) and were run as their own selection.
* `pytest --cov` (`fail_under = 90`) was not re-measured end to end for the same time
  reason. It never ran in the red job — `Lint` failed first — and none of the edits
  removes an executed line, so the coverage arithmetic is unchanged by this fix.