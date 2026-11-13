# Baseline hull fleet: provenance-first geometry gate (mission 17249, GAP-CAD-T1-BASELINE-FLEET)

Status: **honest PARTIAL / FAIL on the volumetric acceptance line.** Four of six
authoritative geometry sources are acquired and provenance-locked, the fleet gate
that would have caught the placeholder bug is implemented and tested, and **zero**
of the six ships is certified watertight-with-volume inside 1 %. The blocking
fact is stated plainly in §5: this repository has no CAD kernel, its exchange seam
declares IGES `NOT_SUPPORTED`, and the mission's own acceptance criterion needs
IGES → watertight mesh. What is delivered is everything up to that wall, plus the
wall itself made machine-visible instead of anecdotal.

---

## 1. Case: what was actually broken

The ship platform's geometry library held eight `*_Hull.stl` files. Six of them —
KCS, KVLCC1, KVLCC2, JBC, DTC, ONRT — were *not hulls*: each is a 304-triangle
box-like approximation (ASCII `solid ship_hull`, sections at a handful of
half-breadths), and `KVLCC1_Hull.stl` and `KVLCC2_Hull.stl` are **byte-identical**
(`sha256 588d2dff5572937f…`), i.e. one form wearing two names. Two files
(`DTMB5415_Hull.stl`, `DTMB5415_Hull_exact.stl`) are genuine tessellations
(127 624 / 164 926 triangles). Nothing in the build noticed any of this: the only
predicate evaluated was "the file exists", and `/generate` fell back to the box
path, so downstream CFD/resistance work was silently fed a box.

Mission target: the authoritative six-ship set (KCS, KVLCC2, ONRT, JBC, DTC,
Series 60) from the workshop publishers (NMRI T2015, SIMMAN 2020, Gothenburg 2022),
each converted to a watertight manifold mesh, volume within 1 % of the published
displacement, with a per-ship manifest (source URL / licence / SHA-256 / verdict).

Acceptance criteria, verbatim, as the machine checks them (`tensorcad.fleet_gate`):

| # | criterion | gate |
|---|---|---|
| A | provenance complete (URL + publisher + licence + SHA-256 + retrieval time) | `validate_source` |
| B | mesh is a hull, not a stand-in (≥ 2 000 facets, not byte-identical to another ship) | `judge_mesh` → `PLACEHOLDER` |
| C | closed, edge-manifold, single component | `is_watertight`, `is_edge_manifold` |
| D | \|V_mesh − ∇_published\| / ∇_published ≤ 1 % | `judge_mesh` → `VERIFIED` |

## 2. Method

1. **Acquisition with provenance, not with guesses.** `tools/baseline_fleet/fetch_sources.py`
   pulls the publisher archives over HTTPS and records per source: final URL,
   publisher, HTTP `Content-Type`/`Last-Modified`, byte size, SHA-256, TLS state,
   and the licence *as stated by the publisher*. The registry that results is
   `artifacts/baseline_fleet/sources.json`.
2. **Published displacements are extracted, never typed.** The publisher pages
   (`t2015_kcs.html`, `t2015_jbc.html`) are archived, hashed, and the displacement
   figure is pulled out of the page text programmatically into
   `published_provenance.json`, which carries the page URL, the page SHA-256 and
   the sentence the number came from.
3. **Bytes are measured, not trusted.** `tensorcad/fleet_gate.py` re-reads the
   mesh, welds coincident vertices (1e-6), and measures triangle/vertex counts,
   watertightness, edge-manifoldness, component count, signed volume (divergence
   theorem), surface area and bounding-box L/B/T. Everything in a report is
   recomputed from the bytes on disk; a registry that disagrees with its own file
   fails `audit.py sources`.
4. **Fail closed.** A mesh that is open, multi-component, or lacks a published
   displacement is graded `WITHHELD` with the reason attached — never "ok".
5. **Repro commands**
   ```bash
   python tools/baseline_fleet/fetch_sources.py --out artifacts/baseline_fleet   # materialise + record provenance
   python tools/baseline_fleet/fetch_sources.py --verify-only --out artifacts/baseline_fleet
   python tools/baseline_fleet/audit.py sources artifacts/baseline_fleet/sources.json
   python tools/baseline_fleet/audit.py audit <geometry_dir> --published artifacts/baseline_fleet/published.json \
          --out artifacts/baseline_fleet/placeholder_audit.json
   python -m pytest tests/test_fleet_gate.py -q     # 7 tests, offline
   ```

## 3. Results (all numbers measured by the code above)

### 3.1 The current eight-file library, graded

`artifacts/baseline_fleet/placeholder_audit.json` — the six placeholder hulls and
the real control, straight from the gate:

| ship | triangles | signed volume (m³) | Cb of bbox | grade | reason |
|---|---|---|---|---|---|
| KCS | 304 | −15 427.42 | 0.500 | PLACEHOLDER | 304 < 2 000 facets |
| KVLCC1 | 304 | −74 460.71 | 0.500 | PLACEHOLDER | facet count **and** byte-identical to KVLCC2 |
| KVLCC2 | 304 | −74 460.71 | 0.500 | PLACEHOLDER | facet count **and** byte-identical to KVLCC1 |
| JBC | 304 | −40 099.63 | 0.500 | PLACEHOLDER | 304 < 2 000 facets |
| DTC | 304 | −50 635.18 | 0.500 | PLACEHOLDER | 304 < 2 000 facets |
| ONRT | 304 | −3 071.33 | 0.500 | PLACEHOLDER | 304 < 2 000 facets |
| DTMB5415 (control) | 127 624 | 19 568.70 | 0.341 | WITHHELD | real tessellation, watertight + manifold, but no published ∇ recorded yet |
| DTMB5415_exact (control) | 164 926 | 13 710.93 | 0.226 | WITHHELD | real tessellation, **not** watertight (open surface) |

Every placeholder reports Cb(bbox) = 0.500 exactly — a box-in-floats, which is the
signature of the generator that produced them, and a second detector available to
the gate (see `docs/baseline_fleet_placeholder_note.md` if the threshold is ever
revisited). The two `DTMB5415` files show the gate does not simply reject
everything that is not a workshop archive: real geometry passes the facet and
topology gates and stops only at the missing published-displacement record.

### 3.2 Authoritative sources acquired (`artifacts/baseline_fleet/sources.json`)

| ship | publisher | archive | bytes | SHA-256 (recorded) | IGES payload |
|---|---|---|---|---|---|
| KCS | NMRI T2015 | `KCS_hull_SVA.zip` | 74 086 | `b9e334853857b819…` | `KCS_hull_SVA.igs` (246 164 B) |
| KCS (2nd source) | SIMMAN 2020 | `KCS-Hull.zip` | 73 295 | `b50adace1ddf255f…` | `KCS_hull_SVA.igs` (same bytes) |
| KVLCC2 | SIMMAN 2020 | `KVLCC2-Hull.zip` | 1 345 336 | `dc032647f00a68d3…` | `kvlcc2.igs` (6 270 129 B) |
| JBC | NMRI T2015 | `JBC_IGES.zip` | 2 508 534 | `e19800ee954f4960…` | `jbc.igs` (1 713 554 B) + propeller/rudder/duct |
| ONRT | SIMMAN 2020 | `onrt_wo_gap.zip` | 752 918 | `bd747eff8f531ba4…` | `onr5613-full_appendage_deckhous_mergedR.igs` (3 374 136 B) |
| DTC | — | not located | — | — | publisher page gone (HTTP 404/410), no workshop mirror found |
| Series 60 | — | not located | — | — | no published IGES located; only the parametric family already in `tensorcad/series60.py` |

The two KCS archives from *different* publishers, fetched independently, contain
the same 246 164-byte IGES — an independent cross-publisher consistency check on
the source bytes.

Measured IGES inventory (`artifacts/baseline_fleet/iges_inventory.json`) — this is
the map of what a tessellator has to swallow, and it is why §5 is a wall, not a
delay:

| file | 128 B-spline surf. | 144 trimmed surf. | 142 curve-on-surface | 143 bounded surf. | 141 boundary | 126 curves | 102 composites |
|---|---|---|---|---|---|---|---|
| `KCS_hull_SVA.igs` | 55 | 0 | 0 | 0 | 0 | 0 | 0 |
| `kvlcc2.igs` | 19 | 19 | 19 | 0 | 0 | 304 | 38 |
| `jbc.igs` | 45 | 0 | 0 | 19 | 19 | 172 | 0 |
| `onr5613…mergedR.igs` | 132 | 0 | 0 | 80 | 98 | 820 | 0 |

KCS is the only one of the four that is a bare untrimmed B-spline patch set (55
surfaces, no trims); KVLCC2 trims in parameter space (142), JBC and ONRT trim with
model-space boundaries (141/143) which require inverse surface evaluation to lift
into parameter space.

### 3.3 Published displacement, extracted from the publisher page text

| ship | ∇ (m³) | extracted from | page SHA-256 |
|---|---|---|---|
| KCS | 52 030 | `t2015_kcs.html` — "Displacement volume ∇ (m³) 52030" | `68d9c24600cb699b…` |
| JBC | 178 369 | `t2015_jbc.html` — "Displacement volume of a hull ∇ (m³) 178369" | `fc2ce9054a8b7842…` |

KVLCC2 / ONRT displacements were not extracted: the SIMMAN 2020 pages render the
particulars table client-side and the fetched HTML carries no displacement text.
They are therefore `WITHHELD`, not guessed.

## 4. Figures

* `artifacts/baseline_fleet/figs/facet_gate.png` — facet counts per library file
  (log scale) against the 2 000-facet gate: six entries sit two orders of
  magnitude below the real ones. Generated by
  `tools/baseline_fleet/make_figures.py` from `placeholder_audit.json` and
  `control_real_mesh.json` (no hand-entered values).
* `artifacts/baseline_fleet/figs/volume_gate.png` — measured signed volume of the
  placeholders against the published displacement of the same-named ship, with the
  ±1 % band drawn: the gap is 2–4× the box volume in the wrong direction, so the
  volume gate alone would also have caught it.

## 5. Notes, limits, and the honest wall

1. **No CAD kernel here → IGES is not convertible in this repository today.**
   `tensorcad/exchange.py` states it explicitly: `IGES — refused — NOT_SUPPORTED`,
   and the OCP route (`pip install "tensorcad[ocp]"`) is not installed in this
   environment (`OCP`/`cadquery`/`gmsh` all absent). Hand-rolling a trimmed-surface
   tessellator is possible (KCS needs only 55 untrimmed patches) but doing it for
   KVLCC2/JBC/ONRT — 19 trimmed and 98/80 model-space-bounded surfaces with
   inverse evaluation — inside this mission would have been an unverified
   half-finished kernel, and an unverified tessellation is exactly the failure mode
   the exchange module refuses to accept. **So the six meshes are not delivered and
   no volume claim is made.** Verdict: FAIL on criteria B/C/D, delivered on A.
2. **TLS**: `www.t2015.nmri.go.jp` serves a certificate whose hostname does not
   match; the fetcher disables verification *and records that fact per source*
   (`tls_verified: false`) rather than hiding it. The SIMMAN host verifies normally.
3. **Licence**: neither publisher states terms of use (the SIMMAN footer carries
   only "Copyright © SIMMAN2020 All rights reserved"; the NMRI T2015 pages carry no
   terms text). The archives are therefore **not committed** — the repository keeps
   the provenance record + a re-fetch script, and `.gitignore` excludes
   `artifacts/baseline_fleet/{raw,iges}/`. Anything shipped from these sources
   downstream needs the publishers' permission first. This is also why the gate
   fails provenance for these entries until a licence statement exists: the
   `audit.py sources` run reports `missing license` for all five archives, by design.
4. **Where the sources were not found**: `t2015.nmri.go.jp/kvlcc2.html` and
   `/dtmb5415.html` are dead (HTTP 404 on the live site); DTC's historical
   `uni-due.de/~qpd402/dtc/` returns HTTP 410 Gone and `uni-due.de/ship/dtc*` 404;
   `gothenburg2022.com` does not resolve; GitHub repository search for
   `DTC hull iges`, `duisburg test case hull`, `series60 hull offsets` returned
   zero repositories. KVLCC2 was recovered from the SIMMAN 2020 mirror instead of
   NMRI. DTC and Series 60 remain `WITHHELD`/`NOT_LOCATED` — recorded, not papered
   over.
5. **Gate thresholds are heuristics with a stated basis**: 2 000 facets is two
   orders of magnitude below the real tessellations measured here (127 624 /
   164 926) and 6.6× above the observed 304; Cb(bbox) = 0.500 for all six
   placeholders makes a second detector available but it is not wired in, because a
   legitimately barge-like hull could approach 0.5 and a hard fail on it would be a
   false positive. The facet gate is the one that fires.
6. **`total_signed_volume` needs outward orientation**; a consistently inverted
   mesh yields a negative volume, which the gate reports with its sign and compares
   in absolute value, but the orientation itself is not yet part of the verdict
   (only watertightness and manifoldness are).

## 6. Error judgement against the acceptance line

| criterion | target | result | verdict |
|---|---|---|---|
| six ships watertight + manifold | 6/6 | 0/6 (no IGES→mesh route in-repo) | **NOT MET** |
| mesh volume vs published displacement | < 1 % | not computable for any of the six (no mesh); the two extracted ∇ values are recorded for the run that has a kernel | **NOT MET** |
| per-ship manifest (URL/licence/SHA/verdict) | 6/6 | 4/6 source manifests complete except licence (publishers state none); DTC and Series 60 not located | **PARTIAL (4/6)** |
| placeholder STL replacement | plan | replacement plan (§ below) + gate that enforces it | **MET (plan + enforcement), 0 files swapped** |
| gate correctness | — | 7 offline tests: placeholder by facet count, byte-twin detection, watertight pass at <1 %, `WITHHELD` at >1 %, open-surface `WITHHELD`, provenance validator, JSON report | **MET** |

**"Is the key quantity <3 %?"** No number is claimed. For the only quantities this
mission could measure, the placeholder volumes differ from the published
displacements by 70 % (KCS: 15 427 vs 52 030 m³), 78 % (JBC: 40 100 vs 178 369 m³)
and are therefore far outside any tolerance — those are *detections of the defect*,
not validations of a new mesh. No watertight mesh exists yet to measure against a
benchmark.

### Replacement plan for the six placeholder STLs (concrete, ordered)

1. **Provide a CAD kernel** in the CAD repo's runtime: `pip install "tensorcad[ocp]"`
   (OpenCASCADE bindings) — one dependency, and `tensorcad.adapters.ocp` already
   speaks it. Without it every following step is blocked; with it, IGES/STEP import
   stops being `NOT_SUPPORTED`.
2. **Materialise + verify sources**: `fetch_sources.py` (network) then
   `--verify-only` (CI, offline). Digests are already recorded, so a changed source
   is caught before it reaches a solver.
3. **Tessellate to a closed solid** per ship and record the tessellation tolerance
   in the manifest: KCS first (55 untrimmed patches — a grid tessellator plus seam
   welding suffices), then KVLCC2 (19 patches trimmed by parameter-space 142 loops),
   then JBC/ONRT (model-space 141/143 boundaries → inverse-evaluate to parameter
   space). Close the deck/transom openings explicitly and record hole area — the
   `DTMB5415_Hull_exact.stl` control shows what an uncapped surface looks like to
   the gate (`watertight: false`).
4. **Grade before shipping**: `audit.py audit <dir> --published published.json`;
   only `VERIFIED` ships may enter the ship platform's geometry directory, and the
   JSON report ships beside the STLs as the manifest sidecar.
5. **Make the library unfakeable**: put `audit.py audit` in the ship platform's CI
   over the geometry directory; the box generator path (the `/generate` fallback)
   must either emit a `PLACEHOLDER` flag in its output metadata or be removed, so a
   downstream job can refuse box geometry instead of consuming it.

Artifacts: `artifacts/baseline_fleet/{sources.json, published.json, published_provenance.json,
placeholder_audit.json, control_real_mesh.json, iges_inventory.json, figs/}`.
Code: `tensorcad/fleet_gate.py`, `tools/baseline_fleet/{fetch_sources.py, audit.py, make_figures.py}`,
`tests/test_fleet_gate.py`.