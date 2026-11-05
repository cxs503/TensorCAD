# PR #77 cure (mission m4088): righting-arm line made ruff-clean

Branch `agent/m4088`, commit `8277e85`, rebuilt from PR #77 head
`65d7a94a6510865fafb473d96908dec842acd2e5` (branch `agent/m3817`, base `main`
`67e1d5c`) on top of current `main`.

## 1. What CI said

`test (3.13)` = failure, `3.10/3.11/3.12` = cancelled (fail-fast). The failing
step was the **first** step, `Lint` (`ruff check .`), exit code 1; `Type check`,
`Tests with coverage gate` and every example step were reported `skipped`
(`.github/workflows/tests.yml`). So the red was a lint red, and everything behind
it was unverified.

Reproduced locally with ruff 0.15.0 on the PR tree: **19 findings** in four files
that the PR adds or rewrites:

| file | findings |
|---|---|
| `tensorcad/righting_arm.py` | F401 unused `Face` import, 2×F841 unused locals (`mean_x`, `mean_y`), 3×B905 `zip()` without `strict=`, 2×RUF007 `zip(x, x[1:])`, B007 loop variable `iterations` unused in body, RUF059 unpacked `vanishing` unused, W292 missing final newline |
| `examples/righting_arm/run_example.py` | B905, RUF007, B007 (`row`), C416, F541, W292 |
| `tensorcad/__init__.py` | RUF022 `__all__` not sorted |
| `tests/test_righting_arm.py` | W292 |

## 2. How the branch was rebuilt (ancestry proof)

`git fetch` of the branch fails in this sandbox, so the PR head tree was taken
from the public codeload tarball of `65d7a94` and unpacked over the worktree
(which sits on `main` = `67e1d5c`, the PR's base). The six files that `main` has
and the PR deletes (`tensorcad/gz.py`, `examples/gz/*`, `tests/test_gz.py`,
`docs/gz.md`) were removed, then

```
git add -A && git write-tree   ->  5b320f25dad5f74a778a83a3a9de9845a93a5f1c
GitHub commit 65d7a94 tree      ->  5b320f25dad5f74a778a83a3a9de9845a93a5f1c
```

i.e. the rebuild is **byte-identical** to the PR head tree, and the cure commit
adds exactly the four lint fixes on top of it. `git diff --stat` of the single
cure commit: 19 files, +2531 / -3018 — the PR content itself, not unrelated work.

## 3. The fixes (behaviour preserving)

* `tensorcad/righting_arm.py`
  * drop the unused `Face` import and the unused `mean_x` / `mean_y` triangle
    centroids (the terms actually returned use `mean_z`, `mean2_*`);
  * sort `_HeeledSection.__slots__` (order is not semantic for `__slots__`) and
    sort `__all__` in `__init__.py`;
  * `zip(self.faces, self.terms, strict=True)` — the two tuples are built
    together, so the assertion can never fire, but a silent truncation can no
    longer hide;
  * `itertools.pairwise` in place of `zip(points, points[1:])` for the
    trapezoidal stability integral and for the "heels strictly increasing" check;
  * the displacement bisection now counts with `while iterations < max_iterations`
    instead of `for iterations in range(1, max_iterations + 1)` whose target was
    only read after the loop — same iteration count, same return tuple;
  * `crossed = self._vanishing()[1]` instead of unpacking an unused name.
* `examples/righting_arm/run_example.py`: `pairwise`, `list(range(...))`, removal
  of a dead `for row in results["sphere"]["convergence"]: pass` loop (it shadowed
  `row` and did nothing), `"GZ(30) = "` without a placeholder-less f-string,
  final newline.
* `tests/test_righting_arm.py`: final newline.

## 4. Verification run locally (Python 3.11.6, ruff 0.15.0, mypy 1.18.2)

| check | result |
|---|---|
| `ruff check .` (before / after) | 19 errors / **All checks passed!** |
| `mypy` (strict, `packages = ["tensorcad"]`) | **Success: no issues found in 33 source files** |
| `pytest tests/test_righting_arm.py` | **41 passed** (18.5 s) |
| `pytest tests/test_public_surface.py tests/test_stability.py tests/test_validation.py` | **39 passed** |
| `pytest tests/test_{hull,hydrostatics,hydrostatics_analytic,hydrostatics_inertia,mesh,mirror_contract}.py` | **156 passed** |
| `pytest tests/test_[a-g]*.py -n 32` | no failure in the 91 % that completed inside the 110 s slice |

## 5. Limits (honest)

The whitebox host was shared and saturated (load average ≈ 1300 on 128 vCPU),
so the **full** suite could not be completed inside the tool's 110–120 s slices:
two full runs reached 88–96 % with **one** failure marker in the 25–32 % region
that I could not localise within budget, while every file group in that region
that I did re-run individually passed (156 + 39 tests). That residual marker is
either a load/timing artefact of the 16–32 worker runs or a test in
`test_[h-w]` after `test_mirror_contract.py`; it is left to CI, which runs the
same command on a dedicated runner. The coverage gate (`fail_under = 90`) is
likewise not reproduced here — coverage was measured on the same aborted runs.

## 6. Numeric error vs baseline

This cure changes no numerics: no expression, tolerance or summation order was
touched (the bisection loop keeps the identical iteration sequence, `pairwise`
is the same pair sequence as `zip(x, x[1:])`, `strict=True` never triggers). The
PR's own analytic baselines therefore keep their values exactly — sphere
`GZ` vs `GM·sin φ` (gate: worst |error| < 1e-3 m) and the wall-sided barge GZ
vs its closed form — and their assertions are part of the 41 tests that pass.
Deviation of the cure commit vs PR head `65d7a94`: **0 (identical computation,
identical asserted thresholds); 0 %**.