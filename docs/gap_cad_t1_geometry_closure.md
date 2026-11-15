# GAP-CAD-T1-SOURCE-HUNT — geometry closure of the fleet (measured, not asserted)

Round: m21253 (TensorCAD) · Anchor: `GAP-CAD-T1-SOURCE-HUNT`
Tool: `python tools/baseline_fleet/geometry_closure.py --all` ·
Record: `tensorcad/fleet/data/geometry_closure.json` ·
Tests: `tests/test_geometry_closure.py` ·
Probes: `python tools/baseline_fleet/geometry_closure.py --verify-only` (offline, what CI runs)

## 1. What was open

The fleet record (`benchmarks/verified/baseline_fleet/fleet_specs.json`) carries, per hull, a
`geometry_scale` and a set of published particulars, but the row *"mesh volume vs published
displacement < 1 %"* in `benchmarks/verified/baseline_fleet/README.md` has been **not run** since
the fleet was built (m20384 says so explicitly: *"No mesh was rebuilt in this round, so the mesh
volume row stays not run"*).  The anchor's item ① also leaves DTC/DTMB5415 open.  This round runs the
row for every hull whose publisher geometry is still reachable, and records the hulls whose geometry
cannot be measured at all as evidence rather than as silence.

## 2. Method (no hand-typed digits)

1. Fetch the T2015 particulars page and the geometry archive, digest both (bytes stay in the
   git-ignored cache; git gets URL + byte count + SHA-256).
2. Read the published table off the *digested* page: the full-scale value of a length/beam/draft row
   is the larger of the row's first two numbers (the KCS rows carry footnote fragments after them,
   `86`, `77`, `0`, `18`, and the first attempt of this tool took `draft = 18 m` from exactly such a
   fragment — the rule now pins the pair, and the implied `full/model` ratio per row is stored as a
   cross-check: [41, 81] and [61, 121] runs both reproduce the workshop's scale
   of ~31.6 for KCS).  The displacement row's bare `3` (the `m 3` exponent) is excluded by a unit rule.
3. Solve two independent waterline conditions on the IGES patch grids
   (`tensorcad.fleet.scale_calibration`): `L_wl(z)/z = L_WL/T` and `B_wl(z)/z = B_WL/T`; the two
   scales must agree (they do for KCS/JBC, they do **not** for ONRT — see §4).
4. Cut the design waterline, integrate the sampled patch cloud below it (trapezoidal, bins derived
   from the hull's own extent), and grade L_WL / B_WL / ∇ / C_b against the published values at 3 %.
5. Every hull is measured at two grid resolutions `[41, 81]` and `[61, 121]`;
   the finer one is reported and both are kept, because the waterline cut is grid-limited (KCS L_WL
   moves from 224.0 m to 226.7 m between the two, i.e. 1.2 %).

## 3. Results

### 3.1 KCS — PASS (file is full scale, digest matches the fleet record)

| quantity | published | measured | error | within 3 % |
|---|---|---|---|---|
| length_wl | 232.5000 | 226.6924 | -2.50 % | yes |
| beam_wl | 32.2000 | 33.1092 | +2.82 % | yes |
| volume | 52030.0000 | 52735.8004 | +1.36 % | yes |
| cb | 0.6505 | 0.6593 | +1.36 % | yes |

File→full scale recovered from the geometry alone: **1.0273** (`1.0055` at
`[41, 81]`), against the fleet record's `geometry_scale = 1.0` — a 2.73 %
disagreement, inside the 3 % band.  The IGES digest of the fetched member is
`6f5556b2b5e58e14…`, byte-identical to `fleet_specs.json`; the page digest is
`68d9c24600cb699b…`.

### 3.2 JBC — FAIL on displacement, and it is *not* reproducible across sampling

| quantity | published | measured | error | within 3 % |
|---|---|---|---|---|
| length_wl | 285.0000 | 286.6181 | +0.57 % | yes |
| beam_wl | 45.0000 | 44.7067 | -0.65 % | yes |
| volume | 178369.9000 | 130214.7070 | -27.00 % | **NO** |
| cb | 0.8580 | 0.6263 | -27.00 % | **NO** |

The waterline rows close; the volume does not, and it does not converge — three sampling settings
give −68.5 %, −27.0 % and −6.2 % (`volume_sensitivity` in the record).  Cause: **29 of 45 patches
(64.4 %) evaluate to a single
point**, so the integrated volume measures the cloud grid, not the hull.  Consequence for the
record: the committed JBC closure of −6.16 % (m20649, `scale_calibration.json`) is one grid's value
and must not be read as a validated displacement; the JBC *scale* (39.739 vs recorded 39.728,
mismatch factor 0.999728) **is** reproduced
independently, which is what item ③ of the anchor asked for.

### 3.3 ONRT — NOT_MEASURABLE (and the recorded scale is not supported by the file)

* digest of the fetched member `ff061ab8584e5b0a…` = `fleet_specs.json` ✓, and the
  archive `bd747eff8f531ba4…` is the cross-host digest the earlier round pinned;
* the file *does* contain 132 type-128 surfaces (`surface_census`), but
  **106 of 132
  (80.3 %) evaluate to a point**, so the two waterline
  conditions disagree by 94.6 % (length gives
  330.4, beam gives 169.8) and no scale may be quoted;
* therefore the recorded `geometry_scale = 48.9355` (the workshop *model* scale) is **unverified
  against the file** — `mismatch_factor = 0.1957`.
  This is declared in the record, not silently corrected: the waterline scale of this file cannot be
  established until the reader can evaluate its surfaces (see §5).

### 3.4 What the host actually serves (anchor item ①, six-hull reading)

| endpoint | status |
|---|---|
| `t2015_index` | 200 |
| `t2015_dtmb5415_page` | 404 |
| `t2015_kvlcc2_page` | 404 |
| `t2015_kvlcc2_gc_page` | 404 |
| `t2015_dtc_page` | 404 |
| `t2015_kvlcc2_iges_guess` | 404 |
| `t2015_dtmb5415_iges_guess` | 404 |
| `t2015_dtc_iges_guess` | 404 |

The workshop index page *links* `dtmb5415.html` and `kvlcc2.html`; both are 404, and no
`kvlcc2/`, `dtmb5415/` or `dtc/` geometry path exists under `/file/Geometry_IGES_files/`.  So the
geometry this host serves is **KCS, JBC, ONRT only** — the six-hull set the anchor assumes is not
retrievable from `t2015.nmri.go.jp`, and the TLS certificate of that host is issued without the
`www.` prefix the anchor quotes (every URL here is the bare name).

## 4. Error judgement (explicit)

| hull | quantity | error | < 3 % ? |
|---|---|---|---|
| KCS | L_WL / B_WL / ∇ / C_b | −2.50 / +2.82 / +1.36 / +1.36 % | yes, all four |
| KCS | file→full scale | +2.73 % vs recorded 1.0 | yes |
| JBC | L_WL / B_WL | +0.57 / −0.65 % | yes |
| JBC | ∇ / C_b | −27.0 % (and −68.5 %, −6.2 % on other grids) | **no** |
| JBC | scale 39.739 vs recorded 39.728 | +0.03 % | yes |
| ONRT | scale | not measurable (conditions disagree 94.6 %) | no basis |

So: the *golden* hull closes inside 3 % on all four quantities; the JBC hull cannot be validated on
displacement by geometry alone; the ONRT file cannot be measured by the current reader.  The
platform's precision standard is met for KCS, honestly failed for JBC-∇ and ONRT.

## 5. Notes, limits, follow-ups

* **Reader defect (open).** 106/132
  ONRT and 29/45
  JBC patches evaluate to a single point while `read_iges` reports the entities as well-formed
  surfaces; on ONRT the parameter list of a `K=M=1` patch carries 41 numbers where the plain type-128
  layout predicts ≤ 26, so the reader and the file disagree about the layout.  This is the single
  blocking item behind the two FAILs and is worth its own round.
* Sampling: the displacement integral is sampling-dominated on incomplete surfaces (JBC §3.2), so a
  "volume closure" claim must always state its cloud spacing — the record now does.
* The 3 % standard is applied to *geometry-measured* quantities; the published columns themselves
  are the references (KCS ∇ 52030 m³, Cb 0.6505, ONRT ∇ 8507, JBC ∇ 178369.9 — each with its page
  digest in the record).
* No bytes redistributed: the archives and pages live under `artifacts/baseline_fleet/raw/closure/`
  (git-ignored); git carries digests, the licence statement and this tool.
