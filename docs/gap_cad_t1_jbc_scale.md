# GAP-CAD-T1-SOURCE-HUNT — JBC model scale, measured from the published IGES

Round: m20649 (TensorCAD) · Anchor: `GAP-CAD-T1-SOURCE-HUNT` (item ③ "JBC 补 scale 因子")
Tool: `python -m tensorcad.fleet.scale_calibration` · Record: `tensorcad/fleet/data/scale_calibration.json`
Tests: `tests/test_scale_calibration.py` (5 passed) · Registry: `tensorcad/fleet/data/source_registry.json`

## 1. What was open and what this round did

The previous round (m20384) closed items ①(ONRT) ②(SIMMAN scrape) ④(page digests) and left
③ **JBC scale** WITHHELD with the note *"needs the IGES control-cloud span, which this round
did not run"*, because the geometry download directory (`/file/Geometry_IGES_files/`) answers
**403** for directory listing.  That conclusion was too pessimistic: the *file* URLs under it
are served normally.  This round fetched

| item | value |
|---|---|
| archive | `https://t2015.nmri.go.jp/file/Geometry_IGES_files/jbc/JBC_IGES.zip` |
| bytes / sha256 | 2 508 534 B / `e19800ee954f4960ef0ccb5bd764b45774b3200cb5c783fdaf17f6b1661f3f9c` |
| member used | `jbc.igs` (1 713 554 B, sha256 `2de08d29…d61e9ab`) |
| provenance check | that sha256 is **byte-identical** to the digest already pinned in `benchmarks/verified/baseline_fleet/fleet_specs.json` for JBC — an independent fetch reproducing the recorded digest |
| units | IGES global section: unit flag 2 (millimetres), model space scale 1.0 |

Bytes are **not** redistributed: the artifact lives in the git-ignored
`artifacts/baseline_fleet/raw/`; git carries URL + digests + licence note only.

## 2. Method (why not the published table, why not the x-extent shortcut)

The T2015 JBC page (`jbc_gc.html`, sha256 `…`) publishes **no model-scale row**: its two
numeric columns are *DESIGNED FULL* and *BALLAST*, both full scale.  The fleet pipeline
therefore used `geometry_scale = None`, i.e. *normalise the point cloud x-extent to L_PP*.
That is biased: the x-extent of a hull surface is **length overall**, so for JBC
7.5428 m → `s = 280/7.5428 = 37.12`, 7 % low.

Instead the scale is read off quantities that mean the same thing in the file and in the
table — the **design waterline**.  With the keel at z = 0 and a uniform scaling by `s`, the
design waterline sits at model immersion `z = T/s`, so

* `s·L_wl(T/s) = L_WL`  ⇔  `L_wl(z)/z = L_WL/T = 17.2727`
* `s·B_wl(T/s) = B_WL`  ⇔  `B_wl(z)/z = B_WL/T = 2.7273`

Each is one equation in `z`, solved by scan + bisection on the file's own surface.  A hull
yields several level crossings (skeg edge, bilge, deck edge), so the roots are **ranked by the
measured waterline length/beam** and the widest one is taken — the design waterline.  (The
first version took the lowest bracket and returned `s = 741.9` from a skeg edge at z = 22 mm;
the synthetic-hull test now pins that behaviour.)

## 3. Results

| quantity | value |
|---|---|
| scale from beam condition | **39.99999** (model B = 1.1250 m vs B_WL = 45.0 m) |
| scale from length condition | 39.47861 (root z = 417.95 mm) |
| adopted (median) | 39.7393 |
| spread between the two independent conditions | 1.30 % |
| LWL at adopted scale | 286.62 m vs published 285.0 → **+0.57 %** |
| B_WL at adopted scale | 44.71 m vs published 45.0 → **−0.65 %** |
| patches in `jbc.igs` | 45 |
| patches that evaluate to a single point | **29 (64.4 %)** |
| displacement at adopted scale | 167 386 m³ vs published 178 370 → **−6.16 %** |
| Cb from geometry | 0.8051 vs published 0.8580 → **−6.2 %** |
| Cm from geometry | 0.9417 vs published 0.9981 → **−5.6 %** |
| scale that *would* close the volume | 40.59 |

The volume closure is **not** within the 3 % gate, so the record says
`gate = "WITHHELD"` and the verdict string in the JSON states that the scale is validated by
the waterline conditions only, never by displacement.

## 4. The finding that matters more than the scale

`jbc.igs` carries 45 type-128 patches; **29 of them evaluate to (0,0,0) for every parameter**,
not because they are empty but because the reader's evaluation returns the origin for them
(all 29 have a real control net: e.g. patch 12 has 3×3 control points spanning
x ∈ [−2.219, +2.110] m, weights 1.0, knots `[0,0,0,1,1,1]`, and still evaluates to zero).
The −6.16 % volume deficit and the −5.6 % midship-section deficit are consistent with meshing
only 16 of 45 patches, i.e. the fleet's JBC mesh is built from 36 % of the hull surface.

This is recorded as an **open defect**, not fixed here: the fix belongs in
`tensorcad/iges.py` (evaluation of the affected patch family) and must be proven by a
volume/Cm closure, not by "the code no longer returns zero".  The calibration module
deliberately *counts and reports* degenerate patches so the next round cannot lose them
silently.

## 5. Error judgement

| quantity | reference | read | error | < 3 % ? |
|---|---|---|---|---|
| model scale (beam condition) | 45.0 m / 1.1250 m | 40.000 | 0.000 % | yes |
| model scale (length condition) | 285.0 m / L_wl(root) | 39.479 | 1.30 % vs beam | yes |
| LWL at adopted scale | 285.0 m | 286.62 m | +0.57 % | yes |
| B_WL at adopted scale | 45.0 m | 44.71 m | −0.65 % | yes |
| displacement volume | 178 369.9 m³ | 167 386 m³ | **−6.16 %** | **NO** |
| midship section coefficient | 0.9981 | 0.9417 | **−5.6 %** | **NO** |

The two failing rows are the honest result of this round: they falsify the assumption that
the published JBC IGES, as *currently evaluated* by this repository, is a complete
description of the published hull.

## 6. Limitations

* The length condition (39.479) and the beam condition (40.000) differ by 1.30 %; LWL is a
  definition-sensitive quantity (5 m of waterline overhang beyond L_PP), so the beam
  condition — which pins the file's 1.1250 m model beam against B_WL = 45.0 m exactly — is
  quoted as the verdict and the length result is kept as the independent check.
* Volume closure uses a 803 847-point cloud and a (10 mm × 5 mm) model-frame field grid; the
  grid resolution study (dx ∈ {5, 8, 10} mm) moved the closure by ≤ 1.0 percentage points, so
  resolution is not the explanation for a 6 % deficit.
* No DTC geometry was obtained: the T2015 index lists only KCS / JBC / KVLCC2 / ONRT /
  DTMB5415 (`dtc.html` → 404), which re-confirms m20384's falsification.
* Bytes not redistributed (licence: published workshop data, research use, no OSI file).