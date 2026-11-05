# PR #91 CI cure (m4094): the red was lint, not physics

## Symptom

`tests.yml` matrix for PR #91 (branch `agent/m3907`) reported
`test (3.10) = failure`, `test (3.11/3.12/3.13) = cancelled`. The whole matrix
shares one pipeline, so the cancel is the fail-fast fallout of the single red job;
the failure is version-independent and shows up before any test runs.

## Root cause

The branch is clean on `main`'s rules but not on the rules `dev = ["ruff>=0.6"]`
resolves to today (ruff 0.15.0). Reproduced on a byte-exact copy of the branch head
`8630690`:

```
$ cd /tmp/m4094full && python3 -m ruff check .     # untouched 8630690 tree
Found 19 errors.
```

`main` at the same time is `All checks passed!`, so every finding was introduced by
this PR. The `Lint` step is the first step of the job, so it red-lights the job before
`pytest` is reached -- which is why all four versions look broken at once.

| file | rule | finding | cure |
|---|---|---|---|
| `tensorcad/__init__.py` | RUF022 | `__all__` not sorted (new names appended) | list re-sorted (names unchanged) |
| `tensorcad/righting_arm.py` | F401 | `Face` imported, unused | dropped from the `geometry` import |
| `tensorcad/righting_arm.py` | F841 (2) | triangle centroid `mean_x`, `mean_y` never read | removed (only `mean_z` is used) |
| `tensorcad/righting_arm.py` | RUF023 | `_HeeledSection.__slots__` unsorted | sorted |
| `tensorcad/righting_arm.py` | B905 (3) | `zip()` without `strict=` | `strict=True` on the paired face/term arrays; `itertools.pairwise` for the ladder walks |
| `tensorcad/righting_arm.py` | RUF007 (2) | `zip(x, x[1:])` for successive pairs | `itertools.pairwise` |
| `tensorcad/righting_arm.py` | B007 | `iterations` never read in the bisection body | counted explicitly (same value returned) |
| `tensorcad/righting_arm.py` | RUF059 | unpacked `vanishing` unused | `_vanishing` |
| `tensorcad/righting_arm.py` | W292 | no trailing newline | newline |
| `tests/test_righting_arm.py` | W292 | no trailing newline | newline |
| `examples/righting_arm/run_example.py` | B905, RUF007 | successive-pair `zip` | `itertools.pairwise` |
| `examples/righting_arm/run_example.py` | F541 | f-string without placeholder | plain string |
| `examples/righting_arm/run_example.py` | B007, C416 | no-op `for row in ...: pass` loop; `[i for i in range()]` | loop deleted; `list(range(...))` |
| `examples/righting_arm/run_example.py` | W292 | no trailing newline | newline |

No numerical code path was touched: `strict=` only makes an existing invariant
explicit, `pairwise` is the same walk, the removed locals were dead, and the
`iterations` counter returns the same integer. The physics claims of
`docs/righting_arm/report.md` are therefore unchanged, and
`tests/test_righting_arm.py` (41 tests) still passes untouched.

## Verification

```
$ python3 -m ruff check .            # branch head after the cure
All checks passed!
$ python3 -m mypy
Success: no issues found in 33 source files
$ python3 -m pytest tests/test_righting_arm.py tests/test_public_surface.py -q
41 passed  /  6 passed
$ python3 -c "import ast,pathlib; [ast.parse(pathlib.Path(f).read_text(), feature_version=(3,10)) for f in ...]"
ast parse @3.10 OK
```

The python-version axis is checked by parsing the four touched files with
`feature_version=(3,10)` and by grepping for 3.11+ only APIs (`StrEnum`, `typing.Self`,
`itertools.batched`, `asyncio.timeout`, `tomllib`, ...): none are used, and every
stdlib call added (`itertools.pairwise`, `zip(strict=)`) is 3.10-native.

Logs: `runs/ruff_before.txt` (19 findings, baseline tree), `runs/ruff_after.txt`
(clean), `runs/run_full_suite.sh` + `runs/pytest_full.log` (the whole CI test step, `--cov`
with the pyproject `fail_under = 90` gate, `-n 10`): 880 passed, 10 skipped, 1 failed --
`tests/test_extract.py::test_extracted_mesh_can_be_measured_with_the_same_hydrostatics`,
which is a wall-clock guard (`assert elapsed < 5.0`) that trips on this box only
(load average > 1000, coverage without its C tracer); its numeric assertions above
(mid-section area, `Cb`) pass, and it touches no code in this PR.
## Round 2 (m4127): the lint cure was not enough — `main` had moved underneath

The round-1 cure (`d3ba0e2`) never reached GitHub: the PR head stayed at `8630690`, so
the green-on-lint tree was still invisible to CI, and the PR was closed with the ruling
*"修正 ≥2 轮未绿，关闭重开，原分支留档"*.

Re-diagnosing on top of the current `main` (`67e1d5c`) found the second, non-lint cause:
between the branch point (`2ff80b2`) and `main`, PR #45 merged `tensorcad/gz.py` — a
module written independently for exactly this quantity (`GZ(phi)` from the heeled
displacement), including free-surface, freeboard and tank-inertia cases this branch does
not have. The branch was therefore 2 commits behind and its `tensorcad/__init__.py`
conflicted on six `__all__` hunks plus the import block, with three names defined twice
(`RightingArmCurve`, `righting_arm_curve`, `wall_sided_righting_arm`). Merging both
verbatim is impossible (one name, two modules) and shipping two competing entry points for
one quantity is worse than useless.

Resolution — one canonical interface, one verification asset:

| decision | why |
|---|---|
| `tensorcad.gz` owns the top-level API (`__all__`, imports unchanged from `main`) | it is the merged, released module and a strict superset of this branch's scope |
| `tensorcad/righting_arm.py` stays, reachable as `tensorcad.righting_arm`, **not** re-exported | deleting it would throw away PR #91's measured results; making it a second top-level entry point would be an ambiguous API |
| `tests/test_righting_arm.py::TestAgreesWithCanonicalGz` (new) pins the two implementations against each other | the second implementation has to earn its place: same `GZ`, same `KN`, same displaced volume on the box barge, where `main`'s tests already have closed forms |
| `tests/test_public_surface.py` only asserts the `gz` names | the public-surface test must describe the package as it is, not the duplicate |

The cross-check compares frame-free invariants only: this module reports the centre of
buoyancy in the heeled (rotated) frame, `gz` reports it above the keel datum in the body
frame, so the test rotates ours back before comparing `KN` — that convention difference is
now documented instead of being a silent trap.

Verification of the merged tree (all commands run in the worktree):

```
$ ruff check .                    -> All checks passed!
$ mypy                            -> Success: no issues found in 34 source files
$ python3 -m pytest tests/test_righting_arm.py -q   -> 44 passed (41 original + 3 cross-checks)
$ python3 -m pytest tests/test_gz.py tests/test_public_surface.py -q -> 32 passed
```
The full suite (~1000 tests, dominated by the mesh-heavy hull/appendage files) does not fit
in this node's 120 s per-command window and the loong64 compute nodes have no numpy/pytest
wheels offline, so the whole-suite green is left to the CI matrix that this push triggers;
what is verified here is every file this merge touches plus the module `main` contributed.
