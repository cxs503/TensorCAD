# LCB / LCG reference frames — reading metre-denominated particulars into the fleet library

Anchor **REPRO-GEOM-FLEET-REFS** (TensorCAD). Machine artefact:
`benchmarks/verified/fleet_refs/lcb_frames_out.json`; registry:
`benchmarks/verified/fleet_refs/lcb_frames.json`; code `tensorcad/fleet/lcb_frames.py`;
tests `tests/test_fleet_lcb_frames.py`. Every number below is rendered from the artefact.

## 1. The gap this closes

Four of the six hulls are blocked on `lcb_percent` and/or `lwl`. Two of those blocks are not
source-availability problems: the pinned documents *do* state the number, in metres, from a
datum the library does not read. JSR-2017 prints `LCG (m) 171.1` for KVLCC2 and `174.059` for
DTC under a header that says `all values refer to the origin located at [AP/CL/BL]`; the NMRI
ONRT page prints `LCB (m) aft of FP 1.625` in the model column and `NA` in the full-scale one.
Reading those cells means converting frames, and converting frames honestly means classifying
what the conversion costs.

## 2. Rule set

Canonical frame: **percent_Lpp_fwd_from_lpp_mid**.

| provenance class | votes | why |
|---|---|---|
| `explicit_frame` | yes | full-scale numeral, datum and reference length stated: the conversion only re-addresses the publisher's own number |
| `scale_derived` | no | model-scale numeral; the full-scale figure is our multiplication, so it corroborates |
| `identity_derived` | no | LCG substituted for LCB; valid only where the document states even keel (TFP == TAP) |
| `computed_not_published` | no | the source's own computation, not a published particular: a comparison, not a reference |
| `unstated` | refused | datum or reference length missing — the converter raises rather than guesses |

## 3. Readings

| id | hull | published | frame | canonical | provenance | votes | digest | vs published |
|---|---|---|---|---|---|---|---|---|
| `kvlcc2_jsr2017_lcg` | KVLCC2 | 171.1 m (LCG) | LCG 171.1 | +3.469 | `identity_derived` | no | verified | -0.011 %Lpp (0.32 %) |
| `dtc_jsr2017_lcg` | DTC | 174.059 m (LCG) | LCG 174.059 | -0.969 | `identity_derived` | no | verified | — |
| `onrt_t2015_lcb` | ONRT | 1.625 m (LCB) | LCB 1.625 x49 | -1.705 | `scale_derived` | no | verified | -0.025 %Lpp (1.46 %) |
| `kcs_arxiv_lcb_lwlmidship` | KCS | 2.3219 m (LCB) | LCB 2.3219 | -0.533 | `computed_not_published` | no | verified | +0.947 %Lpp (63.96 %) |

### 3.1 Quoted source text

* **KVLCC2** — Journal of Ship Research 61(3) 2017 (Ghent University), Table 1 (p. 4, sha256 `172f78d49f36be74…`)
  * frame: Table 1 Main particulars and loading conditions for the KVLCC2 and DTC (all values refer to the origin located at [AP/CL/BL])
  * quote: `Table 1 Main particulars and loading conditions for the KVLCC2 and DTC (all values refer to the origin located at [AP/CL/BL]) KVLCC2 DTC Lpp (m) 320.0 355.0 B (m) 58.0 51.0 CB [-] 0.8098 0.661 ... LCG (m) 171.1 169.0 174.059 170.460`
* **DTC** — Journal of Ship Research 61(3) 2017 (Ghent University), Table 1 (p. 4, sha256 `172f78d49f36be74…`)
  * frame: same Table 1 header: all values refer to the origin located at [AP/CL/BL]
  * quote: `LCG (m) 171.1 169.0 174.059 170.460  (KVLCC2 | heavy ballast | design | light ballast)`
* **ONRT** — NMRI Tokyo 2015 workshop, ONRT geometry page (reproducing the IIHR/DTMB model table) (p. 1, sha256 `5e13427332de8be1…`)
  * frame: LCB (m) aft of FP | model 1.625 | full scale NA; the table's columns are the 1/49 model and full scale (LWL 3.147 m model, 154 m full scale)
  * quote: `LCB LCB (m) aft of FP 1.625 NA`
* **KCS** — Keramati et al., 'Localized Generative Design of Ship Hulls with Diffusion Models and High-Fidelity Assessment' (arXiv HTML full text), Table 5 (p. 16, sha256 `f7af76c786f9a7f6…`)
  * frame: The longitudinal center of buoyancy is measured identically on both hulls as a percentage of waterline length aft of the waterline mid-point; the published KCS value of +1.48 %Lpp is referred to the forward perpendicular, which the bulb displaces.
  * quote: `LCB [%L aft] +1.00 ... (Lwl [m] 232.19, Lpp [m] 230.00); 1.00 % of 232.19 m = 2.3219 m aft of the waterline mid-point`

## 4. Worked case: the KCS LCB that two sources print in two frames

KCS is the only hull with an LCB both publishers state and the reconstruction reproduces, so it
is the place to test whether a frame explanation is *sufficient*. The arXiv paper prints its own
computed KCS LCB, -0.533 %Lpp after conversion, and explains the gap to the
published value by the bulb displacing the waterline datum.

| quantity | value |
|---|---|
| Lpp (m) | 230.00 |
| Lwl (m) | 232.50 |
| bulb overhang Lwl − Lpp (m) | 2.50 |
| frame term (half the overhang, as %Lpp) | 0.543 |
| published LCB (%Lpp) | -1.48 |
| paper's computed LCB, canonicalised (%Lpp) | -0.533 |
| difference (%Lpp) | +0.947 |
| explained by the frame (%Lpp) | 0.543 |
| unexplained (%Lpp) | 0.403 |

**Verdict:** the frame term is a minority of the difference: the paper's own computed LCB sits 0.95 %Lpp away from the published value, of which the bulb displaces at most 0.54 %Lpp; the remainder is a geometry / reconstruction difference, not bookmarking.

This is the honest reading of the anchor's third criterion applied to frames: a datum difference
must not be allowed to absorb a disagreement that is really geometric, and the size of the
frame term is computable, so it is stated rather than argued.

## 5. Effect on the library

* No verdict changes. The two readings that carry new information (`kvlcc2_jsr2017_lcg`,
  `dtc_jsr2017_lcg`) are `identity_derived`, i.e. LCG standing in for LCB under a stated even
  keel, which is arithmetic of ours; `onrt_t2015_lcb` is `scale_derived` (model 1/49 column,
  full-scale cell `NA`); the KCS reading is the source's own computation.
* What *does* change when `--merge` is used: those numerals enter the library as non-voting
  cells (`unit_basis='derived'`), and the library's `derived` block records the distance between
  each converted reading and the voting value — ONRT `-1.705` vs voter `-1.68` (delta
  `-0.025 %Lpp`), KVLCC2 `+3.469` vs voter `+3.48` (delta `-0.011 %Lpp`), DTC `-0.969` with no
  voter at all. The corroboration is therefore visible in the machine artefact instead of
  living in a report paragraph.
* What does *not* change: the acceptance verdicts. `tools/baseline_fleet/refs_acceptance.py`
  counts only voting cells towards the document tier, so ONRT stays failing on
  `lcb_percent` and DTC on `lcb_percent` + `lwl`. That is the intended outcome — a converted
  numeral is evidence *about* the published one, not a second publisher, and the acceptance
  record must not be able to reach a pass by arithmetic.
* KVLCC2's converted LCG lands within 0.05 %Lpp of the published percentage (see table): the
  two publishers' numerals are the same number to four digits, which is strong evidence that
  neither the SIMMAN table nor the JSR table was mistyped — but it is *one* measurement, so the
  verdict stays `SINGLE_SOURCE` until a second lineage states it.

## 6. Limitations

* `identity_derived` is only valid at even keel. A trimmed condition (TFP ≠ TAP) makes LCB and
  LCG differ by the trimming moment; the converter raises rather than converting.
* `scale_derived` needs the scale to come from the document. The ONRT scale (1/49) is stated by
  the IIHR thesis and is implied by the two columns of the NMRI page (LWL 3.147 m model, 154 m
  full scale); both are recorded in the registry, and no floor is quietly assumed.
* The `Lwl` cells remain unfilled for KVLCC2 and DTC: neither JSR-2017 nor any other pinned
  document states them, and a waterline length is not derivable from `Lpp` without the hull.
  That gap stays open and is still the reason KVLCC2/DTC cannot pass criterion 1.
* Digest verification depends on the document bytes still being cached on this host; rows whose
  bytes are gone report `bytes_absent`, which is unknown, not verified.

## Reproduce

```bash
python tools/baseline_fleet/lcb_frames_build.py --raw /tmp/h7 --raw /tmp --doc --fig
python tools/baseline_fleet/lcb_frames_build.py --merge   # non-voting cells into the library
python -m pytest tests/test_fleet_lcb_frames.py -q
```
