# DTC official hull geometry: measured, not withdrawn (anchor GAP-CAD-T1-SOURCE-HUNT)

Machine record: `benchmarks/verified/fleet_refs/dtc_hull_measured.json`
Tool: `tools/baseline_fleet/dtc_hull_measure.py` (stdlib + this repo only, no hand-entered number)
Guard: `tests/test_dtc_hull_measure.py` (5 tests, ~2 s)
Figure: `docs/fleet_refs_dtc_hull_measured.png`

## 1. What was owed, and what was wrong with the standing answer

The anchor's acceptance item ① asks for the **DTC and ONRT geometry files retrieved, with sha256 +
source URL + licence in a fleet-spec-style manifest**. ONRT was already in the store; DTC was not.
The standing DTC statement was a withdraw with a *measurement* behind it
(`benchmarks/verified/fleet_refs/dtc_geometry.json`): "IGES sampling collapsed to a point: 448 590
samples from 1 345 B-spline surfaces span 0.000 m".

That measurement is **stale, and this round falsifies it**. Re-reading the *same cached bytes*
today yields a real point cloud that spans 372.8 m x
51.00 m x 31.75 m. The
DTC file was never broken: the type-128 control-net reader was. Between the two rounds the layout
choice in `tensorcad.iges._split_control` was replaced (a positive-constant weight block is now the
primary signal, instead of the sum of squared second differences, which is identically zero on the
2x4 nets these patches carry and let the "interleaved" reading win). DTC's 1 345 surfaces are
`count_u = 2` patches: exactly the case that fix addressed. **The fix landed; nobody re-measured the
withdrawn hull.** This round does.

## 2. Source and licence (nothing redistributed)

| item | value |
|---|---|
| archive | `igs_ordner.zip` (`dtc-2012-01-10/`), cached as `/root/.cache/hull_sources/dtc_igs_ordner.zip` |
| archive sha256 | `7b18e40d783b670f226a55a7d3f0ccf13c0978e3d80a13a78eced65a7a852541` |
| archive bytes | 2343841 |
| URL | https://www.uni-due.de/imperia/md/content/inam/igs_ordner.zip (workshop page https://www.uni-due.de/ship/dtc.php) |
| licence | publisher states no licence; the workshop publishes the IGES set for benchmarking. Bytes are digested here, never redistributed. |
| hull member | `DTC_hull.igs`, sha256 `53c3be90385c7abae6590276a7225edd1b78edb0e5e4c4a446428296e0edc274` |

Members imported into the platform asset store as
`asset://geometry/hulls/DTC/{DTC_hull.igs, DTC_rudder_skeg.igs, DTC_bilgekeel.igs}` with a
per-hull `MANIFEST.json` (source URL, licence, digests, measured status) and re-hashed by
`python -m tensorcad.assetlib build-registry`: the registry entry for DTC is `status: verified`,
`verified: true`, `problems: []`. **Bytes stay out of the repository** — only digests and pointers.

### 2.1 Two file-level defects worth pinning

1. **Backslash member paths.** The archive's central directory spells members
   `dtc-2012-01-10\DTC_hull.igs` while the local file headers use `dtc-2012-01-10/DTC_hull.igs`, so
   `zipfile.ZipFile.read()` raises `BadZipFile: File name in directory ... and header ... differ`.
   `extract_members()` normalises `ZipInfo.orig_filename` in memory (the header spelling is the one
   the offsets point at) and leaves the bytes untouched. Anyone re-deriving this record must do the
   same, or they will conclude the archive is corrupt.
2. **Units are metres, and the hull is full scale.** Global unit flag 6 = `M`, `model_scale = 1.0`,
   x-extent 372.8 m against `L_PP` = 355.0 m.
   DTC is *not* a 1/59.407 model-scale file like the ONRT archive: any earlier "model scale 1/59"
   assumption applied to this file would have produced a 59x-too-small hull and a 200 000x volume
   error. No normalisation was applied; the volume gate runs on the publisher's own frame.

## 3. Method

`tensorcad.fleet_baseline` chain, no new hydrostatics code:

1. `read_iges` -> 1 345 type-128 surfaces (`partition_surfaces` quarantined
   0 of them);
2. `sample_patches` at two surface densities (1.0 m and 0.5 m);
3. `field_from_points(points, T = 14.5 m, dx = 0.5 m, dz = 0.1 m)` -> half-breadth field
   `b(x, z)`, with interior-hole bridging reported (58088 nodes at
   the fine density);
4. `field_integrals` -> displacement, LCB, waterplane area, max half-breadth (trapezoid integration,
   independent of the mesh route);
5. comparison against the published particulars **injected** from
   `tensorcad/fleet/data/fleet_manifest.json` (the same quadruple the fleet build trusts) and cross
   reference verdicts from `dtc_geometry.json`.

## 4. Results

| spacing [m] | samples | B measured [m] | &#8711; [m&#179;] | C_B | A_WP [m&#178;] | LCB x [m] | &#8711; error | bridged nodes |
|---|---|---|---|---|---|---|---|---|
| 1.0 | 105359 | 51.0005 | 173272.7 | 0.6600 | 15258 | 174.08 | -0.1120 % | 84312 |
| 0.5 | 315616 | 51.0022 | 173459.0 | 0.6607 | 15305 | 174.07 | -0.0046 % | 58088 |

Published (El Moctar, Shigunov & Zorn 2012, as carried by the fleet manifest):
`L_PP 355.0 m, B 51.0 m, T 14.5 m, C_B 0.661, &#8711; 173467.0 m&#179;`
(the two-source library mean is 173 467.5 m&#179;, 2.9e-6 apart from the manifest digit).

Checks, as written by the tool (tolerance in the record):

| check | measured | published | rel. error | tolerance | verdict |
|---|---|---|---|---|---|
| beam | 51.0022 m | 51.0 m | +0.0043 % | &plusmn;0.5 % | PASS |
| block coefficient | 0.6607 | 0.661 | -0.0394 % | &plusmn;1 % | PASS |
| displacement | 173459.0 m&#179; | 173467.0 m&#179; | -0.0046 % | &plusmn;1 % | PASS |

**Is every key quantity within 3 % of its benchmark? Yes** — the largest is C_B at 0.039 %, and the
displacement is 0.005 % from the published value; the beam is an *independent* zero-parameter check
(the file's own geometry, nothing fitted) at 0.004 %.

Figure `docs/fleet_refs_dtc_hull_measured.png`: (left) section-area curve to T = 14.5 m,
maximum 730 m&#178; near midship; (middle) waterline half-breadth against published B/2 = 25.50 m;
(right) displacement versus sampling density, with the 1 % band (both densities inside it, spread
0.107 % — the record's convergence guard).

## 5. Pitfalls / limits

* **The withdraw was stale, not wrong when written.** Any reader who quotes
  `dtc_geometry.json`'s `route.outcome = refused` should read this note instead; that file is left
  byte-identical on purpose as the historical record of the reader defect.
* The patches are **untrimmed** (the census shows 0 type-144 entities), so the sampling grid includes
  the whole parametric domain: the x-extent 372.8 m is the
  *patch sheet*, not a published L_OA, and is not compared with anything. Beam is compared because
  the maximum |y| is genuinely attained on the hull side.
* **LCB has no published reference** (the crossref library reports `MISSING` for the DTC LCB cell),
  so 174.07 m is reported **in the file's own frame, uncompared** — it is not a
  "% of L_PP from AP" value and must not be quoted as one.
* `field_from_points` bins by `max|y|`, which biases the hull outwards by O(cell²/radius); the -0.005 %
  volume says that bias is negligible at dx = 0.5 m here, and the coarse/fine spread bounds it.
* The appendage files (`DTC_rudder_skeg.igs`, `DTC_bilgekeel.igs`) are imported and digested but
  **not measured**: the hull-only displacement is what the publisher's 173 467 m&#179; refers to.

## 6. Reproduce

```
python tools/baseline_fleet/dtc_hull_measure.py            # rewrites the record, prints the table
python -m pytest tests/test_dtc_hull_measure.py -q          # 5 tests: record + recompute from cache
python -m tensorcad.assetlib verify                         # re-hashes the store entry
```
