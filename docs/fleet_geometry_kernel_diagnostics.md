# Fleet geometry kernel: silent control-net misreads made visible (GAP-CAD-T1-BASELINE-FLEET)

Subject: `tensorcad/iges.py` (type-128 reader), `tools/baseline_fleet/diagnose_sources.py`,
`tests/test_iges_layout_tiebreak.py`. Machine-readable evidence:
`benchmarks/verified/baseline_fleet/geometry_diagnostics.json` (produced by the tool
from the fetched bytes, no hand-entered numbers).

Prior state (mission 18488): 1 of 6 hulls `VERIFIED` — KCS, +0.304 % on 52 030 m³;
KVLCC2 `WITHHELD` ("field produced fewer than two usable waterline slices"); JBC
−18.28 %; ONRT −98.91 %; DTC / Series 60 `NOT_SUPPORTED`. Three root causes had been
written down for those three failures. This mission attacks the *diagnosis* first,
because two of the three turned out to be wrong or incomplete.

## 1. Falsified hypothesis: ONRT needs IGES entity 124 (transformation matrices)

Mission 18488 §4.1 blamed ONRT on "the reader documents no transformation-matrix
processing", i.e. per-surface local frames. Measured on the fetched bytes
(`geometry_diagnostics.json → hulls.ONRT.transformation_matrices`):

| file | type-124 records |
|---|---|
| KCS_hull_SVA.igs | 0 |
| kvlcc2.igs | 0 |
| jbc.igs | 0 |
| onr5613-full_appendage_deckhous_mergedR.igs | **0** |

There is no transformation matrix in the ONRT file to ignore, so implementing entity
124 could not have moved that hull by a single millimetre. The file's entity types are
`126 ×820, 128 ×132, 141 ×98, 143 ×80, 314 ×5, 406 ×135`. Implementing 124 remains a
reasonable robustness item; it is *not* the ONRT root cause.

The measurement that does characterise the ONRT file: its control-point cloud spans
**1.076 m × 1.061 m × 1.065 m** — i.e. the geometry is non-dimensionalised by the hull
length, not shipped at the published model scale of 1 : 48.9355. Two independent
agreements with the published particulars at a scale of 154 m (= L_PP):
half-beam 0.0610 × 154 = 9.39 m vs published B_wl/2 = 9.39 m, and keel 0.0357 × 154 =
5.49 m vs published T = 5.494 m. So ONRT carries the `"published model scale"` policy
(`geometry_scale = 48.9355`) while its bytes are `length = 1`, which is exactly the
0.276 factor mission 18488 measured (42.5 m vs 154 m).

## 2. Defect fixed: the layout heuristic is undefined on a degenerate control net

A type-128 record stores its control points either as IGES 5.3 `x, y, z, w`
quadruples or as a block of `n` weights followed by `3n` coordinates (NMRI /
Tokyo-2015 writer). The reader chose by "smallest sum of squared second differences".
For a net with no interior control point (2×2, 2×n) *every* second difference is
identically zero, so every candidate scores 0.0 and the choice fell to tuple sort
order — `"interleaved" < "weight_block"` — with `layout_ambiguous = False`.

Reproduced on the real sources (before the fix): 15 of 132 ONRT patches and 2 of 45 JBC
patches were read as `interleaved`, producing control points inside the unit cube
(`max|y| = 1.000 m` on an 18.78 m-beam hull whose half-beam is 9.39 m). The reader
already had a quarantine hook for "ambiguous" patches; the defect was that this case
never reported itself as ambiguous.

Fix (`tensorcad/iges.py`):

1. `_uniform_weights()` — a candidate whose weights are a strictly positive *constant*
   block is preferred outright. This is a fact about the writers, not a heuristic: all
   four authoritative files carry `weights ≡ 1.0` in every patch (KCS 55/55, KVLCC2
   19/19, JBC 45/45, ONRT 132/132 patches).
2. The `rational` argument (PROP3) is now used. It was accepted and ignored, and the
   plain-triples layout was offered for rational records whenever the tail happened to
   fit, so a rational patch could be parsed as truncated coordinates.
3. When the weight test is silent *and* every candidate is perfectly smooth — the
   degenerate-net case — the patch is returned with `layout_ambiguous = True` instead of
   a silent choice. Fail closed, do not guess.

Effect, measured by the same tool (patches whose layout is now *reported* as undecidable):

| hull | patches | ambiguous before | ambiguous after | quarantined by `partition_surfaces` |
|---|---|---|---|---|
| KCS | 55 | 0 | 0 | 0 |
| KVLCC2 | 19 | 0 | 1 | — (build fails earlier, see §3) |
| JBC | 45 | 1 | 4 | 1 → **4** (41 usable) |
| ONRT | 132 | 0 | 26 | 0 → 26 |

That is the honest shape of the increment: the volume gate did **not** move, because
for these patches the truth is not recoverable from the record alone — 2×2 patches
carry no interior control point, so the byte stream genuinely cannot distinguish the
two layouts. What changed is that the pipeline now says so, and
`partition_surfaces` drops those patches instead of meshing a possibly mis-read net.
KCS keeps 55/55 usable and its `VERIFIED` verdict is untouched (its patches are all
non-degenerate; 41 + 3 pre-existing tests around the fleet stay green).

## 3. KVLCC2: two independent bytes-level faults, one still ambiguous

`geometry_diagnostics.json → hulls.KVLCC2`:

* global field 15 declares **`MM`** (`unit_flag = 2`), so the reader applies
  0.001 m/unit → the whole hull is 0.334 m long;
* the control cloud spans 333.6 / 40.0 / 29.2 in *file* units. With 320 m L_PP, 58 m
  beam and 30 m depth, the only self-consistent reading is *metres*: axis0 = 333.6 m
  (L_PP 320 plus fore/aft overhang, +4.3 %), and one of the two transverse/vertical
  axes spans 29.18 which matches B/2 = 29.0 to **+0.6 %**. The unit flag is therefore
  wrong (a millimetre flag on metres of data), which by itself is enough to make the
  draft-based waterline cut fail: with the 0.001 factor the "hull" is 28 mm deep, so the
  field has no usable waterline slice — the mission-18488 error message, now with a
  named cause.
* Which of axis1 / axis2 is vertical is **not** decided here: depth 30.0 and half-beam
  29.0 are within 3 % of each other, so the extents cannot separate them. Stated as
  open, not guessed. The decisive next measurement is the published KVLCC2 section
  offsets at one station (a few typed numbers from the same SIMMAN page), not more
  geometry heuristics.

## 4. Verification

```
python -m pytest tests/test_iges_layout_tiebreak.py tests/test_iges_units_and_patch_quarantine.py \
      tests/test_iges_baseline.py tests/test_iges_gold_regression.py tests/test_iges_trim_import.py \
      tests/test_fleet_gate.py -q
# 44 passed, 3 skipped, 2 xfailed
python tools/baseline_fleet/diagnose_sources.py --json benchmarks/verified/baseline_fleet/geometry_diagnostics.json --summary
```

`tests/test_iges_units_and_patch_quarantine.py` was updated from `(44, 1)` to
`(41, 4)` usable/quarantined for JBC — the three extra quarantined patches are the
2×2 nets that used to be read as `interleaved` and used silently. That assertion is
the regression lock for this fix.

## 5. Error judgement against the acceptance line

| criterion | target | result | met |
|---|---|---|---|
| six hulls watertight + manifold | 6/6 | 1/6 (KCS, unchanged) | **no** |
| mesh volume vs published displacement | < 1 % | KCS +0.304 % (unchanged); JBC/ONRT/KVLCC2 still `WITHHELD` | **1/6** |
| per-hull manifest | 6/6 | 4/6 sources, licence not stated by either publisher | **partial** |
| placeholder STL replacement | plan + gate | unchanged (gate + plan) | met |
| **no silent control-net misread** | 0 patches | 0 silent; 31 patches now flagged undecidable and quarantined | **met** |

No volume claim is made for JBC/ONRT/KVLCC2 in this mission; a rebuild with the
quarantined patches removed is the next step and is *not* reported as a pass here.

## 6. Notes / limits

* Publisher bytes are not redistributed (neither NMRI T2015 nor SIMMAN 2020 states a
  licence), so the tool reads `artifacts/baseline_fleet/raw/` (`gitignore`d) and the
  committed JSON carries only counts, layouts, extents and SHA-256.
* The ONRT scale finding above is evidence, not a fix: changing
  `GEOMETRY_SCALE_POLICY["ONRT"]` from `("published", …)` to `("inferred", None)` must be
  graded by a rebuild (and by hull-only patch selection, since that file also carries
  deckhouse and rudder surfaces) — not by a policy edit alone.
* The weight test is a writer fact, not a law: a genuinely rational patch with
  varying weights falls back to the smoothness score, and by construction reports
  ambiguity when the net is degenerate. Both branches are covered by
  `tests/test_iges_layout_tiebreak.py`.