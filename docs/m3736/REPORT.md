# Cure of TensorCAD PR #37 (`agent/m2873`) — semantic union into `main`

Mission m3736 · repo TensorCAD · branch `agent/m3736` · target PR branch `agent/m2873`.

## 1. State on entry and the shape of the conflict

| item | value |
|---|---|
| `origin/main` (ours) | `8389c7e` "Merge #33 (m2575) via merge-governor", version `0.20.0` |
| `origin/agent/m2873` (theirs, PR #37) | `206698a` "m2873: rebuild on current main 6c568c31", version `0.19.0` |
| merge base | `6c568c3` — the branch was rebuilt on an *older* main, hence `mergeable=false` |
| files main touched since base | 27 (Wigley family, asymmetric hull docs, CLI, analytic hydrostatics, tests) |
| files branch touched since base | 15 (transom-stern family, docs, example, tests) |
| genuine three-way conflicts | 4: `CHANGELOG.md`, `pyproject.toml`, `tensorcad/__init__.py`, `tests/test_public_surface.py` |

Both sides only **added**; nothing was deleted on either side. The conflict is therefore
purely a "same three lines, different content" collision on the version string, the
`__all__` list, the expected-public-API set and the changelog head — i.e. exactly the
files where a mechanical overlay would silently drop one family.

## 2. Per-file resolution (semantic union, nothing dropped)

| file | main side | branch side | resolution |
|---|---|---|---|
| `pyproject.toml` | `version = "0.20.0"` | `version = "0.19.0"` | **`0.21.0`** — the branch's feature lands *after* main's 0.20.0, so the merged release is one minor step above both |
| `tensorcad/__init__.py` (imports) | `from .wigley import (...)` | `from .transom import (...)` | git auto-merged both import blocks; verified both present and alphabetically ordered (`suboff → transom → validation → voxel → wigley`) |
| `tensorcad/__init__.py` (`__all__`) | `"generate_wigley"` | `"generate_transom_hull"` | **both kept**, sorted: `generate_transom_hull`, then `generate_wigley`; result is the exact set union, 120 names |
| `tensorcad/__init__.py` (`__version__`) | `0.20.0` | `0.19.0` | `0.21.0` (same rule as pyproject) |
| `tests/test_public_surface.py` (`EXPECTED_API`) | 10 Wigley names | 13 transom names | **union of 23**, each block labelled with its family and release |
| `tests/test_public_surface.py` (version assert) | `0.20.0` | `0.19.0` | `0.21.0` |
| `tests/test_public_surface.py` (module list) | adds `"wigley"` | adds `"transom"` | both, in one tuple |
| `CHANGELOG.md` | 0.20.0 (Wigley) + 0.17.0 sections | 0.19.0 (transom) section | new `## 0.21.0` heading carrying the branch's transom entry **verbatim** (renumbered, with an explicit provenance note), followed by main's `0.20.0` and `0.17.0` sections unchanged |

No file was resolved by taking one side wholesale. `tools/report_m2873.py`,
`docs/m2873/outputs/transom_stern/*`, `docs/transom_stern.md`, `examples/transom_stern/`,
`tensorcad/transom.py`, `tests/test_transom.py` came from the branch untouched;
`docs/wigley/*`, `tensorcad/wigley.py`, `tests/test_wigley.py`, `docs/asymmetric_hull/*`,
`tensorcad/hull.py`, `tensorcad/cli.py`, `tests/test_cli.py` came from main untouched.

## 3. Completeness audit (machine, not by eye)

`docs/m3736/verify_union.py` walks `git diff --name-status base..side` for **both** sides and
compares each such path's blob against the merged working tree, flagging every addition,
deletion or content difference that does not match the owning side.

```
counts: ours_only=27  theirs_only=11  both=4
silent_loss: []          ok: true
```

So: 27 main-only paths byte-identical to main, 11 branch-only paths byte-identical to the
branch, 4 conflict files hand-resolved (reported above). Machine-readable:
`docs/m3736/union_audit.json`.

Surface equality (import the merged package and diff against both revisions' `__all__`):

```
version: 0.21.0
__all__ size: 120   expected union: 120
missing: []   extra: []   duplicates: 0
EXPECTED_API == union: True    __all__ == EXPECTED_API: True
```

## 4. Test evidence

Full suite, same interpreter, same flags, `--ignore` of two collection errors that are
pre-existing environment gaps (`safetensors` not installed; they abort collection on main
too, so they are excluded from *both* runs):

| run | tree | passed | failed | skipped | log |
|---|---|---|---|---|---|
| baseline | `origin/main` @ `8389c7e` (exported with `git archive`) | **706** | 2 | 11 | `runs/baseline_main_pytest.log` |
| merged | this branch (semantic union) | **802** | 2 | 11 | `runs/merged_pytest.log` |

* Same two failures in both runs, both pre-existing and unrelated:
  `tests/test_geometry_field.py::test_the_learned_revision_speaks_the_same_protocol` and
  `::test_the_learned_field_gradient_is_autograd_and_points_outward` (torch-backed learned
  field). Verdict: **not below the main baseline** (identical failure set, +96 net new
  passing tests, which are the branch's transom suite plus the widened public-surface test).
* Targeted: `pytest tests/test_transom.py tests/test_public_surface.py tests/test_wigley.py`
  → **132 passed** (22 s).
* Example re-run in the merged tree: `PYTHONPATH=. python3 examples/transom_stern/run_example.py`
  → exit 0, `runs/transom_example.log`, artefacts reproduced byte-for-byte in content.

## 5. Numbers the merged tree reproduces (from the example, no hand-typing)

Incomplete taper integral `∫_τ^1 sin(πu)^p du` — `sinc^p` power series vs two independent
evaluations, 26 points: worst relative disagreement **1.11e-15** (Γ closed form at τ=0,
tanh-sinh quadrature at τ>0).

Mesh vs analytic hydrostatics, `L=10 m, B=2 m, T=0.8 m, p=0.7`:

| τ | 41×13 V err | 321×81 V err | A_wp err (finest) | LCB err |
|---|---|---|---|---|
| 0.00 | −0.880 % | −0.048 % | −0.005 % | ±1.8e-16 |
| 0.05 | −0.813 % | −0.045 % | −0.003 % | −0.002 % |
| 0.20 | −0.785 % | −0.045 % | −0.002 % | −0.002 % |
| 0.35 | −0.770 % | −0.044 % | −0.002 % | −0.002 % |

Design space (`A_t/A_m → Cb`): p=0.7 gives Cb 0.4701→0.4825, p=3.0 gives 0.3097→0.3443.

## 6. Error judgement vs the acceptance criteria

| key quantity | reference | error | <3 %? |
|---|---|---|---|
| mesh volume, finest grid (all τ) | analytic closed form | 0.044 % | **yes** |
| mesh volume, coarsest 41×13 | analytic closed form | 0.77–0.88 % | **yes** |
| waterplane area, finest | analytic closed form | ≤0.005 % | **yes** |
| LCB, finest | analytic closed form | ≤0.002 % | **yes** |
| taper integral series | Γ form / tanh-sinh | 1.1e-15 (≈5 ulp) | **yes** (by 13 orders) |
| version / API surface | union of both sides | 0 missing, 0 extra | **yes** |
| silent content loss | `verify_union.py` | 0 paths | **yes** |
| test suite | main baseline (2 known failures) | same 2 failures, +96 passes | **yes** (not below baseline) |

Honest caveats: (a) the mesh→analytic error is an *inscribed-chord* discretisation error and
the log-log slope against `nx·nz` is only −0.73 over the tested ladder (≈1.5 in h), i.e. the
table is not yet in the asymptotic range for the coarsest grid — the finest grid is what the
0.044 % refers to; (b) the two `test_geometry_field.py` failures exist on `main` before this
merge and are untouched here; (c) no experimental or literature benchmark is involved — this
is pure geometry/hydrostatics against closed forms, which is the strongest available
reference for this family; (d) the merged tree's `.ruff_cache` and the example's
root-level `outputs/` directory are run artefacts and are not part of the deliverable.

## 7. Reproduction

```bash
git checkout agent/m3736                 # merge of origin/main + origin/agent/m2873
python3 docs/m3736/verify_union.py       # completeness audit (exit 1 on silent loss)
python3 -m pytest -q --ignore=tests/test_cli_neural.py --ignore=tests/test_neural_pipeline.py
PYTHONPATH=. python3 examples/transom_stern/run_example.py
PYTHONPATH=. python3 docs/m3736/make_figure.py     # docs/m3736/merge_evidence.png
```
