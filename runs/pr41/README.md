# PR #41 cure evidence (m3213, 管到底 round 1)

`pr41_union_check.py` is the gate: 30 checks, exit 0 only if the merge into
`origin/main` kept both sides and the PR delta sits under the governor's 200 KB
artifact cap. Run it from the repo root; it writes `union_check.{tsv,json}`.

| area | what it proves |
|---|---|
| A1/A2/A2b | no main path disappears; every branch-owned path survives byte-identical, except the three figures that were re-encoded (PSNR >= 40 dB, measured 45.0/45.7/45.8 dB) |
| B | main's lines remain a sub-multiset of the merged file for every file the merge touched (only documented allowances: 1 re-wrapped line in `tests/test_public_surface.py`) |
| C | branch inventory: `heel.py`, `test_heel.py`, `examples/gz_curve/run_example.py`, `docs/gz_curve*`, 117 `__all__` names, 21 module-tuple members, exactly one CHANGELOG bullet, `cli.py` = main's blob |
| D | artifact-class bytes in the PR delta: 50.7 KB < 204.8 KB |

Method note: this clone's history is grafted (`--is-shallow-repository` -> true), so
`git merge-base 292d794 origin/main` cannot resolve and `runs/cure_pr41_*.py` from the
earlier round (base-dependent) cannot be replayed here. The base-free whole-tree
criterion above is authoritative for this round; the older scripts are kept unmodified
as the round-1 record.

Test parity (`head_pytest.log` = merged tree, `main_pytest.log` = origin/main@2ff80b2,
both `pytest -q tests` minus the two collection-error files that import `safetensors`):

| tree | passed | skipped | failed |
|---|---|---|---|
| origin/main | 818 | 11 | 2 (`test_geometry_field`: `safetensors` not installed) |
| merged | 848 | 11 | 2 (same two, same cause) |

+30 tests, all from `tests/test_heel.py`, all green; the failing pair is an environment
gap present identically on main, so the merged tree is at or above main's baseline.
