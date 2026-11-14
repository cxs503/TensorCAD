# GAP-CAD-T1-SOURCE-HUNT ③ — JBC model scale: from derivation to adoption

Round m21074 (TensorCAD). Anchor item ③ of `GAP-CAD-T1-SOURCE-HUNT` ("JBC needs its
scale factor"). This round closes the **plumbing** half of the item: m21021 derived the
scale, but the record was never committed and the fleet manifest still built JBC with
the fallback fit, so the derivation changed nothing downstream.

* derivation tool/record: `tools/baseline_fleet/jbc_model_scale.py` →
  `artifacts/baseline_fleet/jbc_model_scale.json` (now committed);
* adoption tool: `tools/baseline_fleet/adopt_derived_scale.py` (manifest ⇄ record sync);
* policy: `tools/baseline_fleet/build_fleet.py` `GEOMETRY_SCALE_POLICY["JBC"] = ("derived","jbc")`;
* impact measurement: `tools/baseline_fleet/jbc_scale_impact.py` →
  `artifacts/baseline_fleet/jbc_scale_impact.json`;
* tests: `tests/test_jbc_scale_adoption.py` (7 new checks).

## 1. Case and why it matters

JBC (Japan Bulk Carrier) is the Tokyo-2015 benchmark hull. The publisher page
`jbc_gc.html` (SHA-256 `fc2ce9054a8b7842…`,
13694 B) publishes **full-scale particulars only** — no
model-scale row — so the scale cannot be scraped. `fleet_specs.json` therefore carried
`geometry_scale: null`, and `tensorcad.fleet_baseline.build_ship` responds to `None` by
fitting the point cloud's **x-extent** to the published `L_PP`.

That fit is wrong for this hull. The x-extent is the deck-level extreme (bow flare plus
stern overhang): the record measures it at
7.5428 m model space against an implied model `L_PP` of
7.0478 m — a ratio of
1.0702. Fitting x to `L_PP`
shrinks the ship by that same ratio. **No displacement gate can see this**: the same
factor rescales the volume, so the hull stayed self-consistent while being ~5 % small.

## 2. Method

1. Re-ran the derivation tool to regenerate the record from the publisher CAD
   (network live): 45 B-spline patches evaluated on a
   41×41 grid
   = 75645 points. Reproduced m21021's number to the digit
   (1:39.7285), so the published document is reproducible.
2. Committed the record. Previously its 12 record tests **skipped** in CI, which is why the
   numbers were unverifiable; they now execute.
3. Added the `derived` policy so a manifest regeneration reads the record instead of
   silently reverting to the `None` fit.
4. Synced `fleet_specs.json` from the record with `adopt_derived_scale.py` (no value typed
   by hand; `--check` is asserted by the tests).
5. Measured the impact by building the *same* IGES bytes twice through the production
   route, differing only in `geometry_scale` (coarse grid 0.5×0.1 m, 240 ring points).

## 3. Results

Scale derived from two independent published lengths:

| reference | published | model span | scale |
|---|---|---|---|
| B (max beam of waterline) | 45 m | 1.1307 m | 1:39.7987 |
| D (moulded depth) | 25 m | 0.6304 m | 1:39.6583 |
| **adopted mean** | — | — | **1:39.7285** ±0.18% |

Production route, same bytes, two policies (draft cut at 16.5 m):

| quantity | published | x-extent fit (old) | derived (new) |
|---|---|---|---|
| L_wl | 285 m | 270.50 m (-5.09 %) | **287.00 m (+0.70 %)** |
| B_wl | 45 m | 41.97 m (-6.73 %) | **44.92 m (-0.18 %)** |
| ∇ | 178369.9 m³ | 145769 m³ (-18.28 %) | 165741 m³ (-7.08 %) |
| C_b | 0.8580 | 0.7011 (-18.28 %) | 0.7972 (-7.08 %) |
| scale used | — | 1:37.1216 | 1:39.7285 |
| verdict | — | WITHHELD | WITHHELD |

Figure: `docs/figs/jbc_scale_adoption.png` (left: both policies against the 3 % band;
right: the two scale references, the adopted mean, and the x-extent fit that was used before).

## 4. Error judgement — is it < 3 %?

| quantity | reference | value | error | < 3 % ? |
|---|---|---|---|---|
| L_wl, derived scale | 285 m (publisher) | 287.00 m | +0.70 % | **yes** |
| B_wl, derived scale | 45 m | 44.92 m | -0.18 % | **yes** |
| L_wl, old x-extent fit | 285 m | 270.50 m | -5.09 % | **no** |
| B_wl, old x-extent fit | 45 m | 41.97 m | -6.73 % | **no** |
| ∇, derived scale | 178 369.9 m³ | 165741 m³ | -7.08 % | **no — honest FAIL, unchanged** |
| C_b, derived scale | 0.8580 | 0.7972 | -7.08 % | **no — honest FAIL** |
| beam/depth scale disagreement | — | 1:39.7987 vs 1:39.6583 | 0.18% | **yes** |

The volume and C_b deficits are **not** caused by the scale and are **not** claimed fixed.
They are the known unclosed-surface defect of the published archive (29 of 45 type-128
patches evaluate to a point; see the m20649 note), which under-predicts displacement by
about 7 %. Adopting the derived scale even makes that number slightly *worse* in
appearance (-18.28 % → -7.08 %) because the old
fit also shrank ∇ toward the published value — two errors partially cancelling. That is
exactly why the length set is the honest discriminator here.

## 5. Notes, traps, scope

* **Lengths, not volumes, close this item.** A volume gate cannot detect a scale bias that
  rescales the hull; only independently measured lengths can.
* `mesh.extent_m` must not be used as `L_PP`/`L_wl`: it is clipped at the design
  waterline (it reads 287.0 m here, while the undeformed x-extent at deck is
  299.7 m). Use `derived.length_waterline_m`.
* The scale is **derived**, not workshop-stated: no publisher document gives 1:39.73 for
  JBC. Half-spread 0.18 % between two independent lengths is the honest uncertainty; the
  depth reference is biased slightly low by deck-level fittings above the moulded depth.
* Coarse grid (38 s/build) is sufficient for a 5 % bias, two orders of magnitude above
  discretisation noise; no convergence study is claimed for the volumes.
* Publisher bytes are still not redistributed: only digests, spans and the derivation are
  committed. The build test therefore resolves the IGES from the asset store *or* a
  git-ignored checkout copy, and skips neither — it asserts the committed records only.
* Still open on the anchor: item ① DTC geometry (WITHHELD, all publisher hosts dead) —
  unchanged by this round.
