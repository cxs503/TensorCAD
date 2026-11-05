# CI round-2 cure for the `agent/m4118` PR branch

## Reported red
`test (3.12)=failure`, `test (3.10/3.11/3.13)=cancelled` on the PR whose head is
`agent/m4118` (branch head at the time of this round: `add1d64`, which already
contains the round-1 cure `PR #70 cure: load tools helpers by path, drop
sys.path leak`).

## What could and could not be observed from this sandbox
* `git fetch`/`ls-remote`/the GitHub REST API are **blocked or timing out** from
  this host, so the failing job's log and the PR head SHA could not be read. The
  diagnosis below is therefore built from full local parity runs, not from the
  runner's log. This is stated as a limitation, not glossed over.
* Python **3.12 is not installed** here (only 3.11.6), and PyPI is unreachable,
  so the exact matrix interpreter could not be reproduced.

## Parity evidence (Python 3.11.6, numpy 2.4.6, pytest 9.2.0, coverage 7.17.2)
| CI step | local result |
| --- | --- |
| `ruff check .` | All checks passed |
| `mypy` | Success: no issues found in 33 source files |
| `pytest --cov` (full suite, coverage gate 90 %) | **877 passed, 10 skipped, total coverage 93.37 %** |
| `tensorcad-neural --help` | ok |
| CLI artifacts + byte-identical reproducibility | ok |
| `examples/parametric_hull`, `geometry_core`, `suboff`, `suboff_appendages` | ok (all checks passed) |

No test or tool in the tree imports anything outside the declared dependencies
with the exception of `numpy`/`pytest` (checked by AST over `tests/*.py`), and
every optional-extra test is behind `pytest.importorskip`, so the `[test,dev]`
install used by CI is sufficient.

## Root cause of the *diagnostic* problem (fixed here)
`strategy.matrix` defaults to `fail-fast: true`. When one job goes red the three
sibling jobs are cancelled, which is exactly the reported label set
(`3.12=failure`, three `cancelled`). The label therefore says nothing about
whether the break is version specific — which is why this branch has now needed a
second round. `.github/workflows/tests.yml` now sets `fail-fast: false`, so each
version reports its own verdict and the next round can tell a version-specific
bug from a shared one.

## The one CI-fragile assertion found (reported, not changed)
`tests/test_extract.py::test_extracted_mesh_can_be_measured_with_the_same_hydrostatics`
asserted `elapsed < 5.0` seconds for a pure-Python hydrostatics pass over >1000
stations. That number is set by the runner's load and by whether coverage tracing
is on (`pytest --cov` is how CI runs it), not by this code, so it is the single
assertion in the suite that can flip red on a loaded 4-vCPU runner while the code
is unchanged. Making it load-independent needs a measured ratio against a
reference loop; the calibration run was cut short by the mission budget, so the
file is left untouched rather than committed half-fixed.

## Status (honest)
The red itself was **not** reproduced locally: every CI step that can run without
Python 3.12 and without GitHub access is green on this tree. The workflow change
makes the next CI run diagnostic (four independent verdicts instead of one
failure plus three cancellations); the branch still needs the runner's log — or a
Python 3.12 interpreter — to close the red.
