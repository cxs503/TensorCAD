## Independent cross-check of the m20649 result (round m21021)

This file is the **second, independent derivation** of the same quantity. Round m20649 solved
it through the design waterline (`tensorcad.fleet.scale_calibration`, adopted 1:39.7393,
two-condition spread 1.30 %); this round measures the surface breadth and depth directly and
divides the published B / D by them, without any waterline search:

| derivation | method | value | spread between its two references |
|---|---|---|---|
| m20649 (adopted earlier) | waterline length + waterline beam conditions | 1:39.7393 | 1.30 % |
| m21021 (this record) | surface breadth B, surface depth D | 1:39.7285 | 0.35 % |

**Agreement: 0.027 %** — two different constructions of the same
model-space measurement land on the same scale, which is the strongest statement available
short of a publisher-stated number. The two records are kept separate on purpose: they disagree
about *which* published length is the safe reference (`scale_calibration` uses the waterline,
`jbc_model_scale` uses the maximum breadth / depth), so a reader can see the sensitivity.

Neither derivation closes the displacement: m20649 reports −6.16 % (∇ 167 386 m³ vs 178 370 m³)
and traces it to 29 of the 45 type-128 patches in `jbc.igs` evaluating to a single point, i.e. the
hull surface is not a closed shell in this archive. Nothing in this file changes that; the scale
here is likewise validated by lengths only, never by volume.

---

# GAP-CAD-T1-SOURCE-HUNT (3) — JBC model scale derived from the publisher's own CAD

Round m21021 (TensorCAD). Anchor item ③ of `GAP-CAD-T1-SOURCE-HUNT`
("JBC needs its scale factor"). Tool: `tools/baseline_fleet/jbc_model_scale.py`;
record: `artifacts/baseline_fleet/jbc_model_scale.json` (regenerate with
`python tools/baseline_fleet/jbc_model_scale.py`, re-check offline with `--verify-only`);
tests: `tests/test_jbc_model_scale.py`. Every number below is injected programmatically
from the record — none is typed by hand.

## 1. What the publisher states, and what it does not

The Tokyo-2015 page `jbc_gc.html` (SHA-256 `fc2ce9054a8b7842…`,
13694 B) publishes **two full-scale columns only** —
"DESIGNED FULL" and "BALLAST CONDITION". There is **no model-scale row**, which is why
the previous round recorded the JBC scale as WITHHELD: it cannot be scraped.

| published quantity | value | row text the value was read from |
|---|---|---|
| Lpp | 280.000 m | `Length between perpendiculars / 280.000` |
| Lwl | 285.000 m | `Length of waterline / 285.000` |
| B (waterline, max) | 45.000 m | `Maximum beam of waterline / 45.000` |
| D (moulded depth) | 25.000 m | `Depth / 25.000` |
| T (design draft) | 16.500 m | `Draft(Full) / 16.500` |
| ∇ | 178369.9 m³ | `Displacement volume of a hull / 178369.9` |
| Cb | 0.8580 | `Block coefficient (\(C_B\)) / 0.8580` |

Cb closure check (arithmetic on the publisher's own numbers, no geometry involved):
∇/(Lpp·B·T) = 0.85796 versus published Cb = 0.8580 →
**0.0047 %**. The particular table is internally consistent.

The same page links the geometry archive `JBC_IGES.zip`
(SHA-256 `e19800ee954f4960…`, 2508534 B, 7 members)
and the propeller archive `MP687_particulars.zip`
(SHA-256 `feb298a22caef7ef…`, 248313 B).
Neither archive bytes nor page bytes are redistributed: the publisher states no licence,
so git carries digests, matched row text, measured spans and the derivation only.

## 2. Method

1. Download the three sources, digest them (SHA-256) and cache the bytes under
   `artifacts/baseline_fleet/raw/` (git-ignored).
2. Parse the IGES hull (`jbc.igs`, 45 B-spline patches, header units `MM`,
   `metres_per_unit = 0.001`) with `tensorcad.iges.read_iges_surfaces`,
   and **evaluate** every patch on a 41×41 parameter grid
   (75645 points) to get the model-space extent. Control points are *not* used for the
   extent: the control polygon is the convex hull of the surface and over-reports the length.
3. Fix the scale from two **independent** published full-scale lengths, both of which exist in
   the same model space:
   `s_beam = B_full / (model breadth)` and `s_depth = D_full / (model depth)`.
4. Report the mean, the half-spread between the two references, and the residual of the
   round-number hypotheses a reader might hold (1:40, 1:39.7, 1:39.4).
5. Cross-checks that do **not** feed the derivation: keel plane at z = 0, port/starboard
   symmetry, IGES header units, and the model propeller diameter from `MP687_particulars.txt`
   (a second, independent length in the same model space).

The hull's **x extent is deliberately excluded** as a scale reference. It is the deck-level
extreme (bow flare + stern overhang): 7.5428 m, which against Lpp would give a naive
1:37.122 — about 6 % low. That bias is the main trap in this
derivation, so it is reported as a rejected candidate rather than dropped.

## 3. Results

Model-space extent of the hull surface (metres, after applying `metres_per_unit`):
x ∈ [-3.8562, 3.6866],
y ∈ [-0.5653, 0.5653],
z ∈ [0.0000, 0.6304] →
spans **x 7.5428 m / y 1.1307 m / z 0.6304 m**.

| reference | full-scale value | model span | derived scale |
|---|---|---|---|
| beam B | 45.000 m | 1.1307 m | **1:39.7987** |
| depth D | 25.000 m | 0.6304 m | **1:39.6583** |
| mean (adopted) | — | — | **1:39.7285 ± 0.18 %** |

Implied model hull at the adopted scale: Lpp 7.0478 m, Lwl 7.1737 m,
∇ 2.8446 m³ (rounded: a 7.0 m model with 2845 litres of displacement).

## 4. Figure

`docs/figs/jbc_model_scale.png` — left: side projection of the sampled model extent
(x 7.543 m × z 0.630 m); right: the two scale references against the
1:40 / 1:39.7 / 1:39.4 hypotheses.

## 5. Error judgement (explicit "<3 % ?" answer)

| quantity | reference | derived/read | error | < 3 % ? |
|---|---|---|---|---|
| scale from beam | B_full = 45.000 m | model B = 1.1307 m → 1:39.7987 | internal, by construction 0 % | yes (identity) |
| scale from depth | D = 25.000 m | model D = 0.6304 m → 1:39.6583 | internal, by construction 0 % | yes (identity) |
| disagreement between the two references | — | — | 0.35 % | **yes** (0.35 % < 3 %) |
| hypothesis 1:40 vs derived | 1:40 | 1:39.7285 | beam 0.51 %, depth 0.86 % | yes |
| hypothesis 1:39.7 vs derived | 1:39.7 | 1:39.7285 | beam 0.25 %, depth 0.11 % | yes — the best-fitting round value |
| hypothesis 1:39.4 vs derived | 1:39.4 | 1:39.7285 | beam 1.00 %, depth 0.65 % | yes, but 1:39.7 fits better |
| Cb closure of the published table | Cb = 0.8580 | ∇/(Lpp·B·T) = 0.85796 | 0.0047 % | **yes** |
| model propeller diameter (independent length) | MP687 model D = 203.00 mm | 0.203 m, same model space as the hull | consistency only (no full-scale D in the fetched sources) | yes (unit consistency), no full-scale cross-check available |

**Answer: yes — every comparison that has a reference is below 3 %.** The adopted value is
**1:39.73 ± 0.18 %**; it is a *derived* quantity, not a
workshop-stated one, and the record says so.

## 6. Notes, limits, and what is still open

* **Derived, not published.** No document on the T2015 host states the JBC model scale. Anything
  downstream that needs the scale must cite this derivation (tool + record digest), not a workshop.
* **Depth span includes deck fittings.** The z extent ends at the deck-level geometry, which is
  slightly above the moulded depth (hatch coaming / bulwark), which biases `s_depth` low by
  ~0.1–0.3 %. That is why `s_beam` (a moulded dimension) is the larger of the two and why the two
  are averaged rather than one being picked.
* **Sampling deficit.** Extremes are found on a 41×41 grid per patch; the
  curvature of the bilge radius at the widest section bounds the deficit well below the quoted
  half-spread. Stability was measured, not assumed: the sampled spans of x/y/z in metres are
  [7.54197, 1.13069, 0.63039] at 11×11, [7.54243, 1.13069, 0.63039] at 21×21, [7.54277, 1.13069, 0.63039] at 41×41 and [7.54277, 1.13069, 0.63039] at 81×81,
  i.e. the breadth settles to 0.000 % by 41×41.
* **Licence.** "Not stated" on every page. Therefore no bytes in git: only digests, the matched
  row text, the measured spans and the code.
* **DTC (anchor item ①) unchanged this round.** The same T2015 index
  (`dd3bc58f07961906…`) serves exactly five hull pages — `kcs`, `kvlcc2`, `jbc`, `onrt`,
  `dtmb5415` — and no DTC page, so the anchor's "DTC lives on the T2015 host" lead stays
  **falsified**; `web.archive.org` still does not answer from this network (bounded 120 s
  attempt, no response). DTC geometry therefore remains **WITHHELD**, as recorded by the
  earlier round.
* **Why the tool fails loudly**: `--verify-only` rejects a record whose two references disagree
  by more than 2 %, whose digests are malformed, whose IGES header is not MM, or whose keel plane
  is off z = 0. The same checks caught a real transcription bug in this round: the first
  generated record read ∇ as `3.0` m³ because the label row carries digits inside its unit
  ("(m 3 )"); the parser now only accepts a *following* row that starts with a number, and
  `tests/test_jbc_model_scale.py::test_particulars_parser_ignores_unit_digits` pins it.
