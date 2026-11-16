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

> **Superseded in part — see §7.** Those "ambiguous after" counts were measured with the
> reader as it stood at m20163. m20487 added the rejection of non-positive weights, which
> decides these records outright: on the current reader all four hulls report **0**
> ambiguous patches and **0** quarantines. The table above is kept as the historical
> measurement, not as a description of `main`.

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
| **no silent control-net misread** | 0 patches | 0 silent; 31 patches flagged undecidable and quarantined at m20163, 0 flagged on the current reader because every patch has a unique admissible reading (§7) | **met** |

No volume claim is made for JBC/ONRT/KVLCC2 in this mission; a rebuild with the
quarantined patches removed is the next step and is *not* reported as a pass here.

## 6. Notes / limits

* Publisher bytes are not redistributed (neither NMRI T2015 nor SIMMAN 2020 states a
  licence), so the tool reads `artifacts/baseline_fleet/raw/` (`gitignore`d) and the
  committed JSON carries only counts, layouts, extents and SHA-256. The SHA-256 was
  `null` in the original record (the reader never computed one, so the field could not
  be filled); §7 repairs that and every hull now carries the digest of its bytes.
* The ONRT scale finding above is evidence, not a fix: changing
  `GEOMETRY_SCALE_POLICY["ONRT"]` from `("published", …)` to `("inferred", None)` must be
  graded by a rebuild (and by hull-only patch selection, since that file also carries
  deckhouse and rudder surfaces) — not by a policy edit alone.
* The weight test is a writer fact, not a law: a genuinely rational patch with
  varying weights falls back to the smoothness score, and by construction reports
  ambiguity when the net is degenerate. Both branches are covered by
  `tests/test_iges_layout_tiebreak.py`.

## 7. Correction (this revision): the record disagreed with the reader, and the digests were empty

Re-running the tool on the same four publisher files against the current reader produced
different numbers from the committed record. Two independent defects, both of which made
the evidence unusable for a reviewer, plus one overclaim in the field names.

### 7.1 The committed record is stale w.r.t. m20487

`_split_control` gained a rule after m20163: a candidate reading whose weight slots hold a
non-positive value is **inadmissible**, because no writer in these archives emits a zero or
negative weight — such a value means the tail was split at the wrong offset. That rule
decides the records m20163 could only flag:

* JBC DE 386 (3x3 keel plate): the `interleaved` reading puts coordinates in the weight
  slots and yields weights `m20163` could only call ambiguous; measured now, the minimum of
  that candidate's weights is negative (`-1295.938`), so only `weight_block` is admissible.
  It decodes to a **3.406 m x 0.518 m** plate, which is what the keel piece is.
* ONRT/JBC 2x2 nets: same mechanism. Reading a weight block `1, 1, 1, 1` as `x, y, z, w`
  invents the control point (1, 1, 1) — the signature of the original defect — while the
  competing `weight_block` reading of a genuine quadruple tail needs a non-positive weight,
  so exactly one reading survives.

Measured on the publisher bytes, current reader, tool output in
`benchmarks/verified/baseline_fleet/geometry_diagnostics.json` (regenerated by
`python tools/baseline_fleet/diagnose_sources.py`):

| hull | units | m/unit | patches | layouts | why the layout holds | ambiguous | usable | quarantined | SHA-256 (first 16) |
|---|---|---|---|---|---|---|---|---|---|
| KCS | M | 1 | 55 | weight_block 55 | polynomial records, smoothness separates (no weight block in the tail) | 0 | 55 | 0 | `6f5556b2b5e58e14` |
| KVLCC2 | MM | 1e-3 | 19 | weight_block 19 | as KCS | 0 | 19 | 0 | `ed06c287c7a858a1` |
| JBC | MM | 1e-3 | 45 | weight_block 45 | 10 by unique admissible reading, 35 by smoothness | 0 | 45 | 0 | `2de08d290e98cc58` |
| ONRT | M | 1 | 132 | weight_block 132 | 32 by unique admissible reading, 100 by smoothness | 0 | 132 | 0 | `ff061ab8584e5b0a` |

The consequences are corrections, not improvements:

* `tests/test_iges_units_and_patch_quarantine.py` locked `(41 usable, 4 quarantined)` for
  JBC; that assertion **failed** against the current reader (measured `(45, 0)`) and is now
  the measured value, with the plate extents above as the proof that the accepted reading is
  the right one. The test only runs when the cached publisher bytes are present, which is
  why `main` stayed green while the lock was broken — the offline half of the same lock
  lives in `tests/test_iges_diagnostics_tool.py`.
* The control-point cloud extents in the old record were partly produced by misread patches
  and must not be quoted as hull measurements. JBC axis spans were `7.506 / 4.255 / 4.255`
  (two misread patches) and are now `7.506 / 1.125 / 0.630`; ONRT was
  `1.076 / 1.061 / 1.065` and is now `1.011 / 0.150 / 0.179`. In particular the ONRT scale
  corroboration in §1 above used `0.0610` as an axis extremum of that polluted cloud; it is
  not reproducible from the current decode (the extreme half-beam now reads `0.0752` on two
  mirrored 2x2 patches, `0.0638` on the hull proper). The conclusion of §1 — the file is
  non-dimensionalised by hull length, not shipped at the published model scale — rests on
  the length axis and on the entity census, both of which are unchanged; the half-beam
  agreement is withdrawn as evidence.

### 7.2 The SHA-256 field could never be non-null

`diagnose()` read `parsed.__dict__.get("sha256")` but `IgesFile` had no such attribute, so
all four digests in the committed JSON were `null` and nothing tied the counts to the bytes.
Fixed in the reader (`tensorcad.iges.read_iges` now hashes the bytes it parses, so every
consumer — the manifest, the fleet report, this tool — can cite the file it measured) and
the tool now reports that digest. `tests/test_iges_diagnostics_tool.py` fails if a hull in
the committed record has an empty digest, and fails if the digest is not cited in this
document.

### 7.3 The remaining overclaim: `decided_by_weight_test`

The field was computed as "the weights are constant", which is not the same as "the weight
test decided". It read `55/55` for KCS, whose patches are *polynomial* records: they have no
weight block, the weight prior never applied, and the layout came from the smoothness score.
Renamed to `weights_uniform` (what it measures). The tool now also reports, per patch,
`admissible_layouts` and `layout_evidence` — one of `unique admissible reading`,
`weight prior (constant positive weights)`, `smoothness`, or
`smoothness tie, flagged ambiguous` — so a reviewer can see *why* a patch is not ambiguous
instead of trusting a flag. That is the table column "why the layout holds" above, and the
counts are asserted in `tests/test_iges_diagnostics_tool.py`.

### 7.4 What this correction does not do

No volume or watertightness claim is re-graded here for JBC / ONRT / KVLCC2: not being
ambiguous is not the same as being verified, and the hull-level acceptance numbers of §5 are
unchanged (KCS `VERIFIED` at +0.304 %; the other three still `WITHHELD`). The two ONRT
patches whose half-beam exceeds the published `B_wl/2` by 23 % are recorded as an open
question for the hull-level rebuild, not as a reader defect: their alternative reading is
the one that produces the (1, 1, 1) artifact point.

Verified with (offline, no publisher bytes needed):
`python -m pytest tests/test_iges_diagnostics_tool.py tests/test_iges_layout_tiebreak.py -q`
and, with the cached bytes staged under `artifacts/baseline_fleet/raw/`,
`python -m pytest tests/test_iges_units_and_patch_quarantine.py -q` (9 passed).
