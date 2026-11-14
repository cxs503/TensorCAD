# GAP-CAD-T1-BASELINE-FLEET — why two of the six hulls are WITHHELD: the reader drops type-126 surfaces

Mission **m20862** (TensorCAD) · anchor `GAP-CAD-T1-BASELINE-FLEET` · date 2026-11-14
Code: `tensorcad/iges.py` (`IgesFile.surface_census`, `ENTITY_BSPLINE_SURFACE_POLYNOMIAL`),
`tensorcad/fleet_baseline.py` (`build_ship` refusal attribution)
Tool: `tools/baseline_fleet/diagnose_ship.py` · Tests: `tests/test_iges_surface_census.py`

---

## 1. What the ledger said, and what the benchmark actually reported

Acceptance for the anchor: six hulls watertight + manifold, volume within 1 % of the
publisher's displacement, per-ship provenance manifest. `benchmarks/verified/baseline_fleet/result.json`
carries the last recorded verdicts:

| hull | volume got / ref (m³) | rel. error | verdict | recorded reason |
|---|---|---|---|---|
| KCS | 52 187.98 / 52 030 | **+0.30 %** | VERIFIED | — |
| JBC | 145 769.05 / 178 370 | −18.28 % | WITHHELD | volume outside 1 % |
| ONRT | 92.95 / 8 500.78 | −98.91 % | WITHHELD | extent 42.50 m vs L_PP 154 m |
| KVLCC2 | — | — | WITHHELD | `ValueError: field produced fewer than two usable waterline slices` |

Read alone, those three reasons describe *shapes*: one hull too slim, one too small,
one with no waterlines. They do not. They are one reader defect seen three times.

## 2. Measurement: the geometry is in the file, the reader never looks at it

`BSplineSurface` in `tensorcad/iges.py` decodes **one** entity type — 128, the *rational*
B-spline surface. IGES also defines **126**, the *non-rational* B-spline surface, and that
is what CATIA (and the NMRI-2015 archive's CATIA exports) writes in bulk. Census taken
here with `IgesFile.surface_census()` on the durable asset store:

| hull | file | units | type 128 (decoded) | **type 126 (dropped)** | other |
|---|---|---|---|---|---|
| KCS | `KCS_hull_SVA.igs` | m (flag 6) | 55 | **0** | — |
| JBC | `jbc.igs` | mm (flag 2) | 45 | **172** | 19×141 boundary, 19×143 bounded surface, 1×314, 1×406 |
| KVLCC2 | `kvlcc2.igs` | mm (flag 2) | 19 | **304** | 19×144 trimmed surface, 38×102 composite curve, 19×142 curve-on-surface |

The correlation is exact and needs no fitting: **the one hull with zero type-126 entities
is the one hull that verifies** (KCS, +0.30 %); the two hulls carrying the geometry
predominantly in type 126 are exactly the two refused on geometry.

KVLCC2 in detail — the `.igs` global section declares `2H MM` and its data is a CATIA
`kvlcc2_CATIA.igs` export; its 19 decodable type-128 entities sample into a set of
surfaces whose bounding box is `x ∈ [−323.68, 0], y ∈ [−0.42, 29.12], z ∈ [−28.01, 0.00]`
(metres). `x` and `y` reproduce the published L_PP = 320 m and B/2 = 29 m, so the scale
and the frame guard are right; the surface is simply not the hull. Feeding that residue
to the half-breadth field produces a grid with one usable level, and `build_solid` refuses
with the message the benchmark recorded.

## 3. What this mission changed

Two minimal, test-covered changes — the diagnosis, not the decoder:

1. **`IgesFile.surface_census()`** (`tensorcad/iges.py`) returns
   `{"readable_type128": n, "unreadable_type126": m, "by_type": {...}}`. A caller that
   integrates a volume can now see that it is integrating residuum instead of inferring
   it from a downstream symptom. `ENTITY_BSPLINE_SURFACE_POLYNOMIAL = 126` names the gap
   where the number was previously absent.
2. **Refusal attribution** (`build_ship`, `tensorcad/fleet_baseline.py`). When the field
   or the loft fails on a file that carries undecoded type-126 entities, the raised error
   now says so and carries both counts:

   > `KVLCC2: field produced fewer than two usable waterline slices -- the source carries
   > 304 type-126 (non-rational B-spline) surfaces against 19 type-128; this reader decodes
   > type 128 only, so the hull arrived partial and the refusal is a reader gap, not a
   > shape verdict`

   The same attribution is attached to the "no surfaces at all" branch, and the census
   travels in every hull report under `source.surface_census`. A file with nothing
   dropped keeps its plain message — the excuse is not invented when it does not apply
   (`test_refusal_stays_plain_when_no_surfaces_were_dropped`).

Tooling: `tools/baseline_fleet/diagnose_ship.py` walks the chain stage by stage (units →
scale → frame guard → samples below draft → field levels → usable slices) and prints the
stage that broke, so the next round does not have to guess which of six stages refused a
hull.

## 4. Honest status against the anchor's acceptance line

**Not met, unchanged by this mission.** 1 of 6 hulls is watertight + within 1 %
(KCS), 0 of 6 have a complete per-ship manifest signed off this round. The blocker is now
named and bounded instead of being a mystery shape:

| hull | blocker | next step (decoder side) |
|---|---|---|
| KCS | none | — |
| JBC | 172 type-126 entities dropped | decode 126, then re-run; scale already calibrated to 39.7393 (spread 1.30 %) |
| KVLCC2 | 304 type-126 entities dropped | decode 126 (CATIA parameter order, see §5); frame guard, units and L_PP already agree |
| ONRT | 3.146 m model bytes at 1 : 48.9355 → 42.50 m built extent vs L_PP 154 m | re-check the archive member used (`onrt_wo_gap/…mergedR.igs` is a merged appendage/deckhouse file, not the bare hull) before touching scale |
| DTC | no geometry in the store | acquire Gothenburg-2022 member |
| Series60 | own module (`tensorcad/series60.py`), not in this fleet spec | wire into the fleet manifest |

## 5. Why the decoder is not in this mission

Type 126 and type 128 differ in parameter *order*, not just in the rational flag. In 128
the knot vectors precede the control points; in 126 the control points precede the knots,
and the point stride is 3 (polynomial) or 4 (rational) while the existing layout scoring in
`_split_control` was calibrated on 128 tails. Getting that wrong does not fail loudly — it
returns smooth-but-wrong surfaces, which is precisely the failure mode this mission exists
to make visible (`layout_ambiguous` / `partition_surfaces` quarantine, `mesh_vs_field_volume_percent`).
A blind decoder written without a converged volume check on two hulls would trade a
*known* WITHHELD for an *unknown* VERIFIED. The work is one decoder plus a per-hull volume
re-run, and it is scoped as such in the follow-up proposal.

## 6. Reproduction

```bash
python -m pytest tests/test_iges_surface_census.py -q          # 4 passed, offline
python tools/baseline_fleet/diagnose_ship.py --all --specs \
       benchmarks/verified/baseline_fleet/fleet_specs.json     # per-stage refusal report
```

Census for any file, without the fleet chain:

```python
from tensorcad.iges import read_iges
read_iges("/root/tp-agent/assets/geometry/hulls/KVLCC2/kvlcc2.igs").surface_census()
# {'readable_type128': 19, 'unreadable_type126': 304, 'by_type': {...}}
```

Publisher bytes stay out of the repository (they are research material without a
redistribution grant): only URLs, SHA-256 digests and these counts are committed.