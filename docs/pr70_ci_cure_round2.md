# PR #70 round-2 CI cure (agent/m3767)

## What was red
`test (3.12)=failure`, `test (3.10/3.11/3.13)=cancelled` (fail-fast) on PR #70
(`agent/m3767`, head `cb8035a`, wet-surface-uncertainty evidence).

## Root cause found in the branch itself
`tests/test_wetted_uncertainty.py` mutated the **global import path at collection
time**:

```python
sys.path.insert(0, str(REPO / "tools"))
from wetted_uncertainty import (...)
```

`sys.path.insert` in a test module is executed during collection, so `tools/`
stays at `sys.path[0]` for *every* subsequently executed test module of the
session.  Two failure modes follow, and they are order dependent:

* any later module that resolves a name which also exists in `tools/` gets the
  `tools/` copy (silent shadowing);
* conversely the module works only while `tools/` happens to be on the path, so a
  test that is fine per-file can fail in the full-suite run, which is exactly the
  CI configuration (`pytest --cov --cov-report=term-missing`).

Reproduced locally by removing the global insert and re-running the module: the
doc test then raised `ModuleNotFoundError: No module named
'wetted_uncertainty_doc'` (the import had always depended on the leaked path).

## Fix (minimal, one file)
Both `tools/` helpers are now loaded **by file location** with
`importlib.util.spec_from_file_location`, registered in `sys.modules` under their
own names and never put on `sys.path`:

```python
def _load_tools_module(name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, TOOLS / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module
```

No other file is touched.  The suite therefore no longer depends on collection
order and cannot shadow stdlib/project modules from `tools/`.

## Branch hygiene
`agent/m3767` was cut from `fa9c95c` and was 2 commits behind `main` (`#39`,
`#45`).  This round merges the current `main` (`67e1d5c`) into the PR branch so
CI evaluates a current tree; the merge is conflict free and leaves the previous
PR deliverables intact.

## Verification (local, Python 3.11.6, numpy 2.4.6)
| check | result |
| --- | --- |
| `ruff check .` | All checks passed |
| `mypy` | Success: no issues found in 33 source files |
| `pytest -q tests/test_wetted_uncertainty.py` | 11 passed (was: 1 failed after the path fix, proving the module depended on the leak) |
| `pytest -q` per file over the merged tree | all files pass (adapters, analytic_hydrostatics, hull, transom, voxel, suboff, mesh_conditioning, gz, wigley, stability, validation, benchmark subset …) |
| `python examples/...` | unchanged by this PR |

## Limitations (honest)
The GitHub Actions log for the failing job was not retrievable from this sandbox
(outbound API calls are policy-blocked), so the diagnosis is derived from the
code and from local reproduction of the ordering hazard rather than from the
runner's log.  The slowest module `tests/test_neural_pipeline.py` needs 106 s of
CPU per test on the current (heavily loaded) control node and could not be run to
completion end-to-end here; it was verified test-by-test.
