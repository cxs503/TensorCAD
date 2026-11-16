# GAP-CAD-T1-IGES-PIPE — patch-set IGES reader: parameter-record keying + type 128 field order

**Mission**: m22830 · **Branch**: agent/m22830 · **Anchor**: GAP-CAD-T1-IGES-PIPE
**Commit**: see `git log -1` on `agent/m22830`
**Files**: `tensorcad/geomquality/io.py`, `tests/test_iges_parameter_record_keying.py`

## 1. What was wrong (found on publisher specimens, not on self-written files)

The untrimmed patch-set route (`geomquality.io.read_iges`, used for every publisher patch
set: KCS 55 patches, DTC 1345 patches) had two independent defects. Both were invisible
on files this package writes for itself, because writer and reader shared the same
private conventions.

**(a) P-section records were keyed on the wrong column.**
`_iges_parameters` grouped parameter lines by columns 65–72 and only then decoded them.
Columns 65–72 of a P line are *per-line entity metadata*, not the record identity: IGES
5.3 §2.2.4 puts a line's own sequence number in columns 74–80 and terminates a parameter
record with `;`, and the D section's parameter-data pointer names the sequence number of
the record's **first line**. A record longer than one line was therefore cut wherever two
consecutive lines disagreed in columns 65–72; the entity's parameters came out as a
fragment, and whatever field happened to land in the "weights" slot was reported as a
non-positive rational weight.

**(b) Type 128 was decoded with this package's private field order.**
`read_iges` read `deg_u, deg_v, count_u, count_v, PROP1, 0, 0, 0`, then control points,
then knots, then a unit-weight block — exactly what `geomquality.io.write_iges` emits.
IGES 5.3 §4.24 puts `K1, K2, M1, M2, PROP1..4` first, then the two knot vectors, then the
control points (interleaved `x,y,z,w`, or a weight block before the coordinates). Every
publisher file uses the standard order, so their records were decoded as if they were the
private dialect.

Consequence in the ledger's own words, from the previous round: KVLCC2 "fails with
`NURBSError: weights[1] = -121.489` — a data quality issue in the source file, not a code
defect". That conclusion was wrong. −121.49 is a KVLCC2 *x coordinate in metres*, read as
a weight because of (a). The publisher file is conformant: all 19 type-128 records decode
with strictly positive weights once the record is keyed correctly.

## 2. The fix (minimal, one implementation of the standard record)

* `_iges_parameters` now sorts the P lines, keys each record on its first line's sequence
  number (columns 74–80), and ends a record at `;`. A numeric value in columns 65–72 is
  kept as an **alias key**, so files written by this package's own `write_iges` still
  decode without a format migration.
* The type-128 field decode is delegated to `tensorcad.iges.BSplineSurface.from_parameters`
  — the decoder the fleet route already scores against KCS/KVLCC2/JBC/ONRT — whenever the
  record starts with entity type 128 (the standard layout). That decoder also settles the
  interleaved-vs-weight-block control-net layout by scoring net smoothness, which is what
  the publisher files need. The private dialect stays on the old branch, selected by the
  record's first value, so round trips of this package's exports are unchanged.

No behaviour change for conformant *self-written* files: `tests/test_geomquality_io.py`,
`tests/test_iges_baseline.py`, `tests/test_iges_diagnostics_tool.py` —
25 passed, 2 skipped, 1 xfailed.

## 3. Results (numbers read back from the runs; no hand-typed field)

### 3.1 Decode of publisher patch sets (patch-set route, `read_iges`)

| specimen | D-section census | faces decoded | time | control-net weights | bbox x [m] | bbox y [m] | bbox z [m] |
|---|---|---:|---:|---|---|---|---|
| KCS (`KCS_hull_SVA.igs`, SVA Potsdam) | 128×110 | 55 | 0.35 s | all 1.0 (polynomial) | −6.00 … 237.84 | −0.01 … 16.15 | −0.02 … 23.50 |
| DTC (`DTC_hull.igs`, uni-due) | 128×2690, 143×48, 141×48, 314/406/126/102 | 1345 | 9.99 s | all 1.0 | −6.74 … 366.07 | −25.52 … 25.52 | −0.01 … 31.75 |

### 3.2 End-to-end import of KCS (`import_cad`, same entry the CLI uses)

| quantity | measured | published / reference | error | <3 % ? |
|---|---:|---:|---:|---|
| faces in file / refused | 55 / 0 | 55 / 0 | — | ✅ |
| triangles after repair | 4314 | — | — | — |
| beam (2 × max \|y\|) | 32.22 m | B = 32.2 m | **+0.06 %** | ✅ |
| length (bbox x-extent) | 243.84 m | Lpp = 230 m | containment only (+6.0 % bbox, *not* scored: Lpp is a datum, KCS has bow/stern overhang) | ✅ (one-sided) |
| G1 all faces read | True | True | — | ✅ |
| G2 trim on surface | True | True | — | ✅ |
| G3 watertight / G6 single shell | False | unreachable for a half-breadth open patch set (no centreplane, no deck) | — | n/a, reported |

| DTC hull | measured | published | error | <3 % ? |
|---|---:|---:|---:|---|
| beam (2 × max \|y\|) | 51.04 m | B = 51.0 m | **+0.08 %** | ✅ |
| length (bbox) | 372.81 m | Lpp = 355 m | containment only (DTC has a large overhang and appendages) | ✅ (one-sided) |

### 3.3 Route status of the three named specimens

| specimen | content | route | verdict this round |
|---|---|---|---|
| KCS | 55 untrimmed 128 patches | patch-set | **decodes + imports, gates G1/G2/G5 pass, 0 refused faces** |
| KVLCC2 | 19 × (144 trimmed + 142 curve-on-surface), 304 × 126 | trimmed | trimmed reader parses all 19 faces in **1.83 s**; the patch-set route *correctly* refuses (`type 144 out of scope`). Remaining cost is tessellation, not parsing |
| ONRT | 264 × 128, 196 × 141, 160 × 143, 1640 × 126 (appendage+deckhouse merged) | patch-set | **still fails**: `NURBSError: parameter 2.6194 is outside the domain [-35.0718, 2.6194]` — a bounded-surface (143) curve parameter outside its own basis-surface domain. Untouched this round |
| DTC | 2690 × 128 + 48 × 143/141 | patch-set | decodes (new specimen, independent of the three above) |

## 4. Honest limits of this round

* **Not fixed**: ONRT's bounded-surface parameter/domain inconsistency (see §3.3); the
  KVLCC2 trimmed tessellation cost. The anchor's "3 real IGES all pass the gate" criterion
  is therefore **not** closed by this round — KCS and DTC are green, KVLCC2 parses but the
  route choice/tessellation budget is still open, ONRT still refuses.
* **Volume is not measured as displacement** for KCS/DTC here: both specimens are
  open, half-breadth patch sets, and the signed volume of an open shell is not a
  displacement. The scored quantities are the ones a publisher dimension can be compared
  with: beam and length containment.
* The `length` comparison is reported one-sided (bbox ≥ Lpp) because Lpp is a hydrostatic
  datum defined on the design waterline, not a bounding-box length; scoring
  `|bbox − Lpp| / Lpp` would manufacture a 6 % "error" on a correctly scaled model.

## 5. Reproduction

```bash
python -m pytest tests/test_iges_parameter_record_keying.py -q      # 3 passed
python -m pytest tests/test_geomquality_io.py tests/test_iges_baseline.py \
                tests/test_iges_diagnostics_tool.py -q              # 25 passed, 2 skipped, 1 xfailed
python -c "from tensorcad.geomquality.io import read_iges; \
print(len(read_iges('/root/tp-agent/assets/geometry/hulls/KCS/KCS_hull_SVA.igs').faces))"  # 55
```

Baseline check: with `tensorcad/geomquality/io.py` stashed to the pre-fix version, all
three new tests fail (`3 failed`), including the real-specimen one.