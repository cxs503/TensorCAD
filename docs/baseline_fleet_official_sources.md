# Baseline fleet on official sources — published digits re-read, records re-verified

Anchor: **GAP-CAD-T1-BASELINE-FLEET** (mission m20045, branch `agent/m20045`).
Code: `tensorcad/fleet/official.py`, `tools/baseline_fleet/official_fleet.py`,
`tests/test_fleet_official.py`.

## 1. Why this round exists (what was still owed)

Every earlier round of this anchor ended with the same honest caveat: *"no egress, so the
publisher pages and the official IGES archives could not be reached; the particulars are
`retrieval=offline_recollection, verified=false`, and five of six hulls are built by
surrogate laws."* That caveat is now **falsified for the digits**: the environment has
outbound TCP/HTTPS from Python (only the `curl` binary is policy-blocked), so both the
workshop pages and the geometry archives are reachable, and all five archives re-fetch with
the digests that were recorded blind in an earlier round.

| ship | archive | bytes | sha256 (first 16) | matches recorded |
|---|---|---|---|---|
| KCS | NMRI T2015 `KCS_hull_SVA.zip` | 74 086 | `b9e334853857b819` | yes |
| KVLCC2 | SIMMAN 2020 `KVLCC2-Hull.zip` | 1 345 336 | `dc032647f00a68d3` | yes |
| JBC | NMRI T2015 `JBC_IGES.zip` | 2 508 534 | `e19800ee954f4960` | yes |
| ONRT | SIMMAN 2020 `onrt_wo_gap.zip` | 752 918 | `bd747eff8f531ba4` | yes |

Page digests are pinned in `tensorcad/fleet/official.py:PAGES` (KCS `68d9c246…`, JBC
`fc2ce905…`, ONRT `5e134273…`, KVLCC2 summary `c1a668f6…`).

## 2. How the digits are read

`parse_particulars(ship, html)` does not transcribe: it flattens the publisher's `<tr>/<td>`
table and reads a quantity out of the **cell of the column that publisher uses for full
scale** (`COLUMN_INDEX = {KCS: 2, JBC: 2, ONRT: 3, KVLCC2: 1}`; ONRT prints its model column
first and its full-scale LCB as `NA`, which the parser reports as `reduced_scale` rather than
silently adopting). Each value keeps its raw cell text, its row, the page URL and the page
sha256, so any reader can re-run the parse on the same bytes.

Two independent checks then run on the extraction:

* **page closure** `Cb·Lpp·Bwl·T` must reproduce the page's own published ∇ (`<0.5 %`);
* **record fidelity**: the packaged fleet record is compared field-by-field at 0.2 %.

## 3. Results (all numbers injected from `published_particulars.json`, no hand copy)

### 3.1 Published particulars vs the packaged records

| ship | quantity | record (offline recall) | publisher page | rel. diff | verdict |
|---|---|---|---|---|---|
| KCS | Lpp, Bwl, T | 230.0 / 32.2 / 10.8 | 230.0 / 32.2 / 10.8 | 0.00 % | agree |
| KCS | Cb | 0.651 | 0.6505 | 0.077 % | agree |
| KCS | **Cm** | 0.971 | **0.9849** | **1.41 %** | **contradicted** |
| KCS | ∇ | 52 030 | 52 030 | 0.00 % | agree |
| KVLCC2 | Lpp, T | 320 / 20.8 | 320 / 20.8 | 0.00 % | agree |
| KVLCC2 | Bwl | 58.0 | 580 (page) | 90 % | page defect (§3.2) |
| KVLCC2 | **Cm** | 0.9956 | 0.998 | 0.24 % | contradicted |
| KVLCC2 | ∇ | 312 622 | 312 622 | 0.00 % | agree |
| JBC | Lpp, Bwl, T | 280 / 45 / 16.5 | 280 / 45 / 16.5 | 0.00 % | agree |
| JBC | **Cm** | 0.992 | 0.9981 | 0.61 % | contradicted |
| JBC | ∇ | 178 378.2 (closure) | **178 369.9** | 0.005 % | now external |
| ONRT | Lwl, Bwl, T | 154 / 18.78 / 5.494 | 154 / 18.78 / 5.494 | 0.00 % | agree |
| ONRT | Cb | 0.535 | 0.535 | 0.00 % | agree |
| ONRT | ∇ | 8 518.1 (recall) | 8 507 t → 8 299.5 m³ | 2.57 % | see §3.3 |

Longitudinal/beam/draft/Cb/∇ of KCS, KVLCC2 and JBC reproduce to ≤0.08 % — the recalled digits
were sound. The **midship coefficients were not**: all three recalled `Cm` values are wrong by
0.24–1.41 %, which matters because `parametric_sections` used `Cm` to shape the section area
law of those surrogates. This is exactly the class of defect the anchor was asking for, and
it was invisible without re-reading the page.

### 3.2 A page-level defect the closure test exposes

KVLCC2's SIMMAN summary row prints `Bwl (m) = 580`. Around that digit the closure gives
**3 126 217 m³** against a published **312 622 m³** (900 % off); with the decimal point
restored (`58.0`) it gives **312 621.7 m³** (1e-6 relative). `resolve_digits()` therefore
admits `beam = 58.0` and records the repair verbatim in the ship's `corrections` list —
page digit, admitted digit, and both closures. No other digit is ever altered, and a page
whose closure agrees is left untouched (KCS/JBC/ONRT have empty correction lists).

### 3.3 ONRT publishes two currencies that disagree by 2.4 %

ONRT publishes a mass (8507 ton) and a Cb defined on `LWL·BWL·T` (0.535). Turning the mass
into a volume at 1.025 t/m³ gives **8 299.5 m³**; the Cb closure gives **8 500.8 m³**. The
tool returns **both** candidate references and reports the mesh error against each, rather
than choosing the flattering one — a closed loop that "validates" a mesh against whichever
reference it happens to fit is not a validation.

### 3.4 Verdict on the shipped fleet records

`fleet_manifest_official.json` (written by the tool, round-trips through `load_manifest`)
now carries, for the four hulls with a published hull IGES:

* `geometry_route = "iges_field"`, `geometry_path` = the archive member, `geometry_sha256`;
* `provenance.retrieval = "fetched"`, `verified = true`, URL and publisher licence quoted
  (both publishers state **no** licence — the bytes are digested, not redistributed);
* `displacement_basis = "literature_value_verified"` (so `reference_is_external()` is true,
  where it was false for every ship before);
* the corrected digits from §3.1–3.2, and the correction list in `notes`.

## 4. What is **not** done (honest)

* **The meshes are not built in this round.** Job `cj2450` wrote the official manifest and
  then ended without `official_build.json` (the IGES→field loft of four hulls at 8 samples/m
  did not complete inside the slot); the coarse re-submission `cj2476` (2 samples/m, dx = 1.0 m,
  dz = 0.25 m, 120 ring points) was still **queued** when this round closed (queue depth 25–41
  jobs). The geometry stage therefore remains **WITHHELD** for this mission, and the anchor's
  "6 hulls watertight + volume < 1 %" acceptance is **not** met yet. Nothing here claims it is.
* **DTC and Series 60 stay surrogate.** `uni-due.de/ship/dtc.php` returns HTTP 404 (the CAD
  used to live there) and no DTC hull IGES is in the workshop archives; Series 60 is a
  published *family*, not a single official file. Both are reported `WITHHELD` with the reason
  rather than built from an invented file.
* **No redistribution.** IGES bytes and page bytes live in `.cache/` (git-ignored); the
  repository carries digests and the fetch script only.

## 5. Error judgement

| quantity | worst case | criterion | verdict |
|---|---|---|---|
| L/B/T/Cb/∇ of KCS, KVLCC2, JBC vs publisher page | 0.077 % (KCS Cb) | <1 % | **pass** |
| Cm of KCS/KVLCC2/JBC vs page | 1.41 % | <0.2 % | **fail** — recalled digits corrected in the official manifest |
| page closure vs published ∇ | 2.16e-6 (KCS), 4.65e-5 (JBC), 1e-6 (KVLCC2 after repair) | <0.5 % | **pass** |
| ONRT mass-derived vs closure volume | 2.43 % | <1 % | **fail by the publisher**, reported as two references |
| mesh volume vs published ∇ | — | <1 % | **not run** (jobs cj2450/cj2476, see §4) |

3 % criterion: the digit-level comparison is at 0.077 % worst case and the single failure
(`Cm`, 1.41 %) is reported as a defect and corrected, so no admitted digit is above 3 %.
The mesh-level criterion is **not evaluated** — no mesh was produced in this round.

## 6. Reproduction

```
python tools/baseline_fleet/fetch_sources.py --out artifacts/baseline_fleet
python tools/baseline_fleet/official_fleet.py extract \
    --pages artifacts/baseline_fleet/raw/pages --out artifacts/baseline_fleet_official
python tools/baseline_fleet/official_fleet.py build \
    --src artifacts/baseline_fleet --out artifacts/baseline_fleet_official \
    --samples-per-metre 2.0 --dx 1.0 --dz 0.25 --ring-points 120
python -m pytest tests/test_fleet_official.py -q      # 20 passed
```

The test suite runs the parser on markup shaped like the publishers' tables (the publisher
pages themselves are not redistributed) and re-checks the pinned page digests and digits on
the real bytes when a fetch has materialised them; without a fetch it **skips**, never fakes.