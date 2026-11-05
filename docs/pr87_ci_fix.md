# PR #87 (agent/m3872) — CI red: what failed and what was changed

The `tests` workflow (`.github/workflows/tests.yml`, matrix 3.10–3.13) failed in
`test (3.11)`; the other three jobs were cancelled by the matrix's default
`fail-fast`. The job failed at its **Lint** step, before any test ran, and the
**Type check** step was red as well. Both were reproduced locally in a worktree
cut from `agent/m3872`.

## Reproduced failures

| Step | Command | Before | After |
|---|---|---|---|
| Lint | `ruff check .` | `Found 24 errors` | `All checks passed!` |
| Type check | `mypy` | `Found 4 errors in 1 file (checked 33 source files)` | `Success: no issues found in 33 source files` |
| Tests | `pytest tests/test_trim.py` | `33 passed` | `33 passed` (+ 7 guard tests) |

### Lint (24 → 0)

`examples/trim/run_example.py`, `tensorcad/trim.py`, `tests/test_trim.py`:

* `RUF100` (6) — `# noqa: E402` directives that the per-file ignore already covers.
* `F401` (2) — `bonjean.section_properties` and `stability.StabilityError` imported
  but never used.
* `UP035`/`UP037` (3) — `Sequence` from `typing`, quoted annotations.
* `RUF022`/`RUF023` (2) — unsorted `__all__` and `__slots__`.
* `W292` (2) — missing newline at end of file.
* `B905` (3) — `zip()` without `strict=`: the station/value pair and the
  area/section pair are equal length by construction, so `strict=True` states it.
* `SIM108` (2) — the same regula-falsi step was written twice as an `if`/`else`
  block. Both are now one helper, `_false_position(a, fa, b, fb)`, with the flat
  residual fallback and its `# pragma: no cover` in a single place.
* `RUF059` — an unused `height_b` unpack replaced by `fb = _residual(b)[0]`.
* `E501` (2 + 1) — three lines over 110 characters, wrapped.

### Type check (4 → 0), `mypy --strict`

* `unused` was annotated `set[tuple[int, int]]` but stores **edges** between
  crossing keys; the annotation is now `set[tuple[tuple[int, int], tuple[int, int]]]`.
* `adjacency.get(current, ())` → `.get(current, [])` (the value type is `list`).
* The local `key` in the chain walk was reused for two different shapes; the edge
  variables are now `edge`.
* **A real annotation bug**: `trimmed_hydrostatics(..., preserved=...)` was typed
  `Sequence[float] | None` while the body calls `float(preserved)`. `preserved` is
  a displaced volume in m³, so it is now `float | None`. Runtime behaviour is
  unchanged (every caller passed a scalar); `mypy` was telling the truth.

## Coverage gate

`pyproject.toml` gates `pytest --cov` at `fail_under = 90`. Measured on
`agent/m3872` with the trim suite alone (`coverage run --branch -m pytest
tests/test_trim.py`), `tensorcad/trim.py` — the 1008-line module this PR adds —
stood at **89 %** (447 statements / 136 branches). The missing lines were almost
all refusal guards (`TrimError` for a non-finite waterline, a non-finite trim
angle, degenerate box arguments) plus the trapezoid fallback and the `to_dict`
surfaces, i.e. paths the numbers-only suite never enters.

`tests/test_trim_guards.py` closes those: `_round` on non-finite/non-float input,
`_integrate` on a one-station and on a non-uniform sample, an indexed station whose
bucket is empty, the four `BoxParticulars`/`box_trim_tangent` argument refusals,
and `hydrostatic_particulars` / `trimmed_hydrostatics` refusing a non-finite
waterline or trim angle. Seven tests, no change to `tensorcad/`.

## Also in this branch

`CHANGELOG.md` had been truncated mid-entry: the PR's copy stopped in the middle of
the 0.14.0 section and dropped **715 lines** of release history (0.14.0 onward and
every earlier release). The new 0.23.0 entry is kept and the history from 0.14.0
down is restored verbatim from the branch point, so the diff of this file is now
purely additive.

## Not verified here

The full-project `pytest --cov` percentage was not re-measured locally: the suite
runs several minutes and the measurement host was saturated by concurrent jobs, so
each attempt was killed before coverage wrote its data file. The authoritative
total is CI's. Everything above that *was* measured is quoted from a local command
in this worktree.
