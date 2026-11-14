# Baseline fleet — two IGES reader defects found on the authoritative bytes

Case: `GAP-CAD-T1-BASELINE-FLEET` (standard-hull authority geometry set).
Scope of this note: **the reader**, not the fleet verdict. No acceptance item of the
anchor is claimed complete below.

## 1. What was inspected

Cached publisher bytes from the m18488 acquisition run (`/tmp/fleet_fetch/iges/`,
not redistributed: neither publisher states a licence; digests live in
`benchmarks/verified/baseline_fleet/fleet_specs.json`):

| file | units | entities |
|---|---|---|
| `KCS_hull_SVA.igs` | m | 55 × type 128 |
| `jbc.igs` | mm | 45 × 128, 19 × 141, 19 × 143 |
| `kvlcc2.igs` | mm | 19 × 128, 19 × 142, 19 × 144 |
| `onr5613-full_appendage_deckhous_mergedR.igs` (ONRT) | m | 132 × 128, 98 × 141, 80 × 143, 820 × 126 |

## 2. Defect A — the control-tail layout tie-break was decided by string order

`_split_control` chooses between three readings of a type-128 control-point tail
(interleaved `x,y,z,w`, a weight block followed by coordinates, and plain triples).
The smoothness score that was supposed to decide is **identically 0.0 on a 2×2 or
2×n net**, because such a net has no interior control point. The previous code then
kept the first candidate in insertion order, i.e. `"interleaved"`.

Measured consequence on ONRT DE 27 (`128,1,1,1,1,0,0,1,0,0 …`, P records 14–23):
the weight block `1,1,1,1` was read as a control point (1.0, 1.0, 1.0) and the true
coordinates became weights, one of them **negative** (−0.052).

| file | patches with a non-positive weight before | flag `layout_ambiguous` before → after |
|---|---|---|
| ONRT | 15 / 132 | 26 → **0** |
| JBC | 2 / 45 | 4 → **0** |
| KVLCC2 | 0 / 19 | 1 → **0** |
| KCS | 0 / 55 | 0 → 0 |

Fix (two rules, both in `tensorcad/iges.py`):
1. a candidate carrying a non-positive weight is rejected — no rational patch in
   these archives has one, so a zero/negative member proves the tail was split at
   the wrong offset;
2. when several readings claim constant weights, the one whose constant block is
   **present in the record** (`weight_block`) beats the one that merely assumes an
   all-ones vector (`polynomial`).

Model-space bounding box after the fix (model units): ONRT `1.0112 × 0.1503 ×
0.1794` (before: `1.0756 × 1.0610 × 1.0655`, the unit-cube artifact).

## 3. Defect B — DE pointers use first-record numbering in the ONRT file

The reader keys every directory entry by its **second** 80-column record
(`second[73:80]`, so keys are even), while the ONRT/Rhino export writes pointers to
the **first** record (`141,1,3,27,…` names DE 27, whose key here is 28). Every
cross-reference in that file therefore missed, and `trimmed_surfaces()` aborted with
`boundary DE 1 missing`. Added `IgesFile.resolve()`: exact key first, then
`pointer + 1` — but only for an odd pointer, so a lookup can never walk onto an
unrelated entity. `trimmed_surfaces()` now uses it for the basis surface, the
boundary and the boundary curves.

## 4. Falsified hypothesis (negative result, recorded deliberately)

m18488 recommended, as step 1, "implement IGES entity 124 (transformation matrix)
in `tensorcad/iges.py` … rebuild ONRT". **The ONRT file contains no type-124 entity
at all** — its entity histogram is `{126: 820, 128: 132, 141: 98, 143: 80, 314: 5,
406: 135}`. Entity 124 cannot be the reason ONRT's cloud spans 1.006 m, so that
recommendation is withdrawn rather than implemented.

Likewise, the layout defect above does **not** explain the displacements. Rebuilding
through `tensorcad.fleet_baseline.build_ship` before and after the fix, same grid
(`dx = 1.0 m`, `dz = 0.4 m`, `samples_per_metre = 1.5`):

| hull | ∇ field [m³] before | after | vs published | verdict |
|---|---|---|---|---|
| KCS | 52 081.082 | 52 081.082 | +0.120 % | VERIFIED (no regression) |
| JBC | 150 276.616 | 150 276.616 | −16.59 % | WITHHELD |
| ONRT | 104.090 | 104.090 | −97.74 % | WITHHELD |

Byte-identical to 15 digits, because the patches the fix re-reads were already
quarantined by `partition_surfaces` and so never reached the field. The JBC gap is
therefore **not** a mis-read control net: its half-breadth max is 21.09 m against a
published B/2 of 22.5 m (−6.3 %), i.e. breadth is missing in the geometry that *is*
read, exactly as m18488 suspected.

## 5. What remains (concrete, measured)

* ONRT type-143/141 records do not follow the conforming layout: `143, <form>=1, SO,
  count, B1, B2` and `141, <form>=1, N, (type, pointer)×N`. That reading decodes
  consistently for **76 of 80** wrappers (basis pointer → a type-128 DE, boundaries
  → type-141 DEs; `count` is 1 or 2, and the reader's `BoundedSurface` supports only
  one boundary). Implementing it is the next step that can actually move ONRT.
* ONRT is still `WITHHELD` for a second, independent reason: `z` spans −0.0495 …
  +0.0589 m and neither frame convention (keel at 0 / waterline at 0) yields two
  waterline slices, so the builder raises `field produced fewer than two usable
  waterline slices`.
* JBC/HVLCC-class breadth: the field max|y| must be compared against published B/2
  *before* any volume verdict is quoted (a breadth gate, not a volume gate).

## 6. Reproduce

```bash
python -m pytest tests/test_iges_layout_tiebreak.py -q     # 4 passed
python -m pytest tests/test_iges_baseline.py tests/test_iges_units_and_patch_quarantine.py -q
```

The ONRT/JBC inspection snippets need the cached bytes (`fetch_sources.py`); the new
tests need no network and no data, they encode the two records verbatim.
