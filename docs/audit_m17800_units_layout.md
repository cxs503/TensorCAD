# m17800 · GAP-CAD-T1-BASELINE-FLEET: length-scale audit and two reader defects

Anchor: `GAP-CAD-T1-BASELINE-FLEET` (standard-hull baseline fleet, 6 hulls, published
geometry -> watertight mesh -> volume within 1 % of the published displacement).
This node does **not** claim the anchor closed; it reports two defects that were hiding
behind the anchor's verdicts, fixes them, and puts the remaining distance to acceptance
in numbers.

## 1. What was wrong (defects, both in `tensorcad/iges.py`)

**D1 — units flag read from a content scan (silent factor 1000).** The flag was located
by "first real number followed by a small integer" in the global section. Global fields
7-11 are the precision magnitudes (`32,38,6,308,15`), which have exactly that shape, so
the scan returned **6 (metres) for the millimetre JBC file**, whose real flag is 2. Every
JBC length was then scaled by 1000. Fixed by Hollerith-tokenising the global section
(`,,`-style and `1H,,1H;`-style headers now land the flag on the same token), reading
field 13/14/15 positionally, and cross-checking the numeric flag against the units *name*
(disagreement is a recorded warning, the name wins). The model-space scale (field 13) is
now read and applied as well.

**D2 — one ambiguous patch aborted the whole file.** `_split_control` raised when two
control-point layouts scored within 4x in control-net roughness. For JBC that was a
single 3x3 patch (DE 386, a flat ~3.4 m x 0.5 m piece) and it made all 45 patches
unreadable. The reader now returns the smoother net with `layout_ambiguous = True`, and
`fleet_baseline.partition_surfaces` quarantines such a patch (census carried into the
report; the field's interior-hole fill covers it from its neighbours) instead of guessing
or failing.

## 2. New gate: length scale before volume

A model-scale file compared against a full-scale displacement produces a "99.99 %"
volume error that says nothing about the mesh. `build_ship` therefore checks
`mesh x-extent / published L_PP` first (`0.85 … 1.20` accepted) and withholds with the
measured ratio, the implied scale factor and the volume that scaling would give. This is
the check that turns the corpus below into a diagnosis.

## 3. Results (machine-generated from `artifacts/baseline_fleet/fleet_run.json`)

| hull | units/flag (source) | patches used/total | quarantined | mesh ∇ [m³] | err vs published | C_B | extent/L_PP | verdict |
|---|---|---|---|---|---|---|---|---|
| KCS | m/6 (positional) | 55/55 | 0 | 52214.08 | +0.354 % | 0.6528 | 1.0391 | VERIFIED |
| JBC | mm/2 (positional) | 44/45 | 1 | 4.44 | -99.998 % | 0.0000 | 0.0268 | WITHHELD |
| KVLCC2 | mm/2 (positional) | 19/19 | 0 | n/a | n/a | n/a | n/a | WITHHELD |
| ONRT | m/6 (positional) | 132/132 | 0 | 0.01 | -100.000 % | 0.0000 | 0.0065 | WITHHELD |

Reasons recorded by the code (verbatim, `fleet_run.json`):

* **KCS** — (none)
* **JBC** — volume 4 m3 is -100.00% off the published 178370 m3 (tolerance 1%); geometry extent 7.499 m vs published L_PP 280.0 m (ratio 0.0268): model-scale bytes or a units error, so a volume comparison against the full-scale displacement is not meaningful; scaling by 37.339^3 would give 231004.2 m3
* **KVLCC2** — mesh construction failed: ValueError: field produced fewer than two usable waterline slices; declared units mm (flag 2, positional); model-space scale 1.0; raw coordinate extent 333.617 x 29.795 x 28.122 file units vs published L_PP 320.0 m: the control points are not in the ship frame
* **ONRT** — volume 0 m3 is -100.00% off the published 8501 m3 (tolerance 1%); geometry extent 1.000 m vs published L_PP 154.0 m (ratio 0.0065): model-scale bytes or a units error, so a volume comparison against the full-scale displacement is not meaningful; scaling by 154.000^3 would give 20423.7 m3

![scale and volume](figs/m17800_scale_and_volume.png)

## 4. What this means for the anchor's acceptance

* **KCS** stays VERIFIED (`+0.354 %` vs the
  publisher's 52 030 m³, watertight, manifold) — unaffected by either defect.
* **JBC** is now read correctly (millimetres, 44/45 patches, 1 quarantined). Its mesh
  x-extent is `0.0268` of the published 280 m, i.e. the
  distributed IGES is the **model** (1:37.34), so the meshed 4.44 m³ has to be compared
  with the model-scale displacement (178 369.9 / 37.34³ ≈ 3.43 m³): the scaled mesh is
  ≈ +29 % high, so the patch set is *not* the design-condition hull below the waterline
  (it carries more than the hull body). Honest WITHHELD, with the number.
* **ONRT** — the cached file is a 1.00 m model fragment (`extent/L_PP = 0.0065`),
  so no comparison is meaningful; unchanged WITHHELD, now for a stated reason instead of
  a −99.99 % volume.
* **KVLCC2** — meshing fails; the measured record shows the control points of the CATIA
  file are not in the ship frame (extent 333.617 x
  29.795 x
  28.122 file units for a 320 m ship):
  a type-124 transformation / local frame has to be resolved before KVLCC2 can be graded.
  WITHHELD with the measured extent (this is the next node, not a claim).

So: **1 of 6 hulls VERIFIED, 3 WITHHELD for stated, measured reasons**, and the two
reader defects that made the other five uninformative are fixed and regression-tested.

## 5. Verification

```
python -m pytest tests/test_iges_units_and_patch_quarantine.py tests/test_iges_baseline.py \
                 tests/test_fleet_gate.py -q      # 20 passed, 1 xfailed (pre-existing)
python runs/m17800/run_fleet.py                   # writes artifacts/baseline_fleet/fleet_run.json
```

The regression tests bind the fix to the publisher's bytes: the same global fields that
the old scan read as `6` (metres) are asserted to be `2`/`MM`, the real DE-386 tail is
asserted to resolve with `ambiguous = True` instead of raising, and KCS is pinned as
still-metres so the fix cannot silently move the one verified hull.

## 6. Limits / next steps

* No outbound network in this environment: literature particulars entered offline are
  flagged `retrieval=offline_recollection, verified=false` in `fleet/data/fleet_manifest.json`;
  digits are not laundered.
* KVLCC2 needs type-124 transformation support (next node).
* JBC/ONRT need the full-scale (not model) distribution, or an explicit acceptance at
  model scale with the publisher's model displacement.
