# Golden geometry gate: mutation red-team

- schema: `tensorcad.golden.redteam/1`
- gate mode: quick (coarse level, as the git hook runs it)
- null control (unmodified copy): PASS
- mutants killed: **7/7** (kill rate 100%), 96.53 s wall clock

| mutant | kind | patch site | expected finding | observed findings | verdict |
|---|---|---|---|---|---|
| `reference-wigley-waterplane` | reference | `tensorcad/wigley.py` | `reference_not_reproducible` | `reference_not_reproducible` | KILLED |
| `reference-series60-waterplane` | reference | `tensorcad/series60.py` | `reference_not_reproducible` | `reference_mismatch`, `reference_not_reproducible` | KILLED |
| `numeric-wigley-grid-silently-coarsened` | numeric | `tensorcad/golden.py` | `geometry_changed` | `geometry_changed` | KILLED |
| `identity-orientation-winding` | identity | `tensorcad/neural/extract.py` | `geometry_changed` | `geometry_changed` | KILLED |
| `numeric-sphere-radius` | numeric | `tensorcad/golden.py` | `drift`, `geometry_changed` | `drift`, `geometry_changed`, `reference_not_reproducible` | KILLED |
| `numeric-5415-mirror-axis` | numeric | `tensorcad/golden.py` | `invariant_violated` | `invariant_violated` | KILLED |
| `baseline-reference-hand-edited` | reference | `docs/golden/tensorcad_geometry_golden.json` | `reference_not_reproducible` | `reference_not_reproducible` | KILLED |

## What each mutant proves

- **reference-wigley-waterplane** (reference): the closed form of the Wigley waterplane area is perturbed by 1e-3; the stored reference no longer equals what its declared source computes, so a hand-copied criterion cannot survive a source change
- **reference-series60-waterplane** (reference): same proof on the second analytic family (Series 60 exponents)
- **numeric-wigley-grid-silently-coarsened** (numeric): the fine grid is silently coarsened without re-blessing: in full mode the fine measurements move (drift) and in quick mode, where only the coarse level is re-measured, the fine mesh fingerprint is what fails
- **identity-orientation-winding** (identity): marching tetrahedra stops repairing inverted winding: the surface area takes the norm, the sphere/ellipsoid volume the absolute value and the extents are unchanged, so every measured number is bit-identical and only the geometry fingerprint can fail
- **numeric-sphere-radius** (numeric): the case definition itself changes (sphere radius 0.05 % larger): the reference moves with it, so only the drift lock and the fingerprint fire
- **numeric-5415-mirror-axis** (numeric): the case declares the mirror axis as x although the shipped shell mirrors about y: the measured extents do not change at all, so only the invariant identity can fail -- a check that is neither drift nor fingerprint
- **baseline-reference-hand-edited** (reference): someone edits the blessed Wigley displacement by +0.4 % by hand; the gate recomputes the reference from the source and refuses the pass

## How to reproduce

```
python -m tensorcad.golden --check            # the gate itself, on the committed tree
python tools/golden/mutation_redteam.py \
    --json docs/golden/redteam_quick.json --markdown docs/golden/redteam.md
```

The machine-readable companion is `docs/golden/redteam_quick.json` (schema
`tensorcad.golden.redteam/1`). Mutants are applied to a temporary copy of `tensorcad/`,
the golden baseline and the DTMB 5415 fixture; the worktree is never modified. The null
control run is included so that a harness that fails for an unrelated reason cannot be
mistaken for a kill.

