# Board audit m20506 — GAP-GEOM-ASSET-LIB (subject m20465, TensorCAD)

**Verdict: L3 (system integration), gate expected to rule on the evidence below.**
Subject: mission 20465, branch `agent/m20465` @ `9de9850`, status `succeeded`, platform report id 18636.
Audit contract: three-part kit (run artifacts / benchmark ref-vs-got / executable repro), benchmark judge =
`tensor_benchmark_cases`, then L0–L5 grading where the level is the functional substance, not the page count.

## ① Subject and benchmark provenance

The anchor asks for a *durable geometry byte store + registry + growth metric + consumer surface* in TensorCAD.
The benchmark judge has **no case registered for this anchor**:

```
SELECT count(*) FROM tensor_benchmark_cases WHERE anchor_slug='GAP-GEOM-ASSET-LIB';  -- 0
```

The identical physical metric (KCS displacement vs the publisher's 52030 m³, watertight + edge-manifold)
*is* registered under a sibling anchor and has already passed:

```
id=52  GAP-CAD-T1-BASELINE-FLEET / CASE-CAD-BASELINE-FLEET-KCS-DISPLACEMENT  status=passed
tolerance: <1% relative error on displacement AND watertight AND edge-manifold; <3% ecosystem line
```

So the numeric leg of m20465 can still be judged — against case 52's tolerance — but the anchor itself is
unregistered, and the mission's own `result.server_quality_audit` is
`{gate_valid:false, key_present:false, gate_missing_warning:true}`: it registered no gate key.

## ② How the audit was done

1. Read the mission result, the branch diff (`9de9850`), the 7 run artifacts and `docs/asset_library.md`.
2. Re-ran the subject's own kit in its worktree: `pytest tests/test_assetlib.py` (14), `pytest tests/test_assetlib.py tests/test_fleet_baseline.py` (43),
   `python -m tensorcad.assetlib verify`, `scan-specs`, `growth`.
3. Re-ran the physics repro **independently** (own output path, own invocation) and compared numbers bit-for-bit.
4. Checked integration: is the consumer reached from a main pipeline path, or only from its own tooling?
5. Checked reproducibility: are the bytes that produce the numbers versioned in the repo?
6. Generated every number in this report and in `verdict_m20506.json` programmatically from the artifacts
   (`docs/board-audit-m20506-assetlib/gen_audit.py`); artifact snapshots are copied into `artifacts/` here so the audit
   verdict does not depend on another worktree staying alive.
7. Because the subject branch is not on `main` yet, `subject_snapshot/` holds **byte-identical copies** of the subject's
   implementation (`tensorcad/assetlib/*.py`, `fleet_baseline.py`, `tests/test_assetlib.py`, the repaired `fleet_specs.json`,
   the registry), each with its SHA-256 in `verdict_m20506.json:evidence.snapshot_sha256`; the gate resolves `impl_files`
   against these while it cannot resolve them against the ecosystem checkout.
8. Self-checked the verdict with the gate's own judge function (`maturity_gate.judge(3, evidence, 'TensorCAD', 20506)`
   → `(3, [])`), recorded in `artifacts/audit_gate_selfcheck.txt`, so the claimed level is the level the machine re-check
   should reproduce — not a self-declaration.

## ③ Results

### 3.1 Pointer repair (criterion ①: dead pointers cleared)

| kind | pointers | before dead/withheld/resolved | after dead/withheld/resolved/cited |
|---|---|---|---|
| geometry | 8 | 8 / 0 / 0 | 0 / 0 / 8 / 0 |
| page | 27 | 0 / 27 / 0 | 0 / 0 / 0 / 27 |
| published_source | 4 | 0 / 4 / 0 | 0 / 0 / 0 / 4 |
| archive | 4 | 0 / 4 / 0 | 0 / 0 / 0 / 4 |

Totals: 43 pointers → 43; dead 8→0,
withheld 35→0, resolved 0→8,
cited 35. Auditor reran `scan-specs` and got the same tri-state.
**Caveat:** the 35 "cited" pointers are publisher URLs (`durable publisher URL (no local byte to verify)`) — the pointer
is now honest, the byte is still absent. Geometry pointers only: 8/8 resolved; all-pointer fraction is
8/43 = 0.186. The growth headline "0.00→1.00" is
correct *scoped to geometry pointers*, and must not be read as "43/43 bytes present" (finding D6).

### 3.2 KCS rebuilt from store bytes (criterion ③: the bytes are usable)

Source: `artifacts/assetlib/physics_from_store.json`; reference = published particulars carried by the spec
(52030 m³, Cb 0.6505, LCB −1.48 %Lpp); acceptance = 1.0% (subject's declared line, stricter than the platform's 15%).

| dx [m] | ∇ got [m³] | field volume [m³] | error vs 52030 | Cb got | LCB got [%Lpp] | triangles | verdict |
|---|---|---|---|---|---|---|---|
| 0.125 | 52075.28 | 52110.71 | +0.087% | 0.65106 | -1.493 | 52320 | VERIFIED |
| 0.25 | 52134.65 | 52145.44 | +0.201% | 0.65181 | -1.502 | 52320 | VERIFIED |
| 0.5 | 52214.08 | 52224.23 | +0.354% | 0.65280 | -1.506 | 52320 | VERIFIED |

Auditor rerun (`--ships KCS --dx 0.5`, own output file):
got 52214.0807 m³, error **+0.353797%**, sha256 match
`True`, difference from the subject's archived number
+0.000000 m³ → **reproduced**.

### 3.3 Fleet scope

| ship | sha match | verdict | best |volume error| | reason |
|---|---|---|---|---|
| KCS | True | VERIFIED | -0.087% | - |
| KVLCC2 | True | NO_VERDICT | n/a | ValueError: field produced fewer than two usable waterline slices |
| JBC | True | WITHHELD | -99.998% | - |
| ONRT | True | WITHHELD | -100.000% | - |

### 3.4 Growth metrics

| metric | baseline | current | delta |
|---|---|---|---|
| store-managed bytes | 0 | 20165115 | 20165115 |
| store assets | 0 | 13 | 13 |
| resolvable geometry pointer fraction | 0.0 | 1.0 | 1.0 |
| dead pointers | 8 | 0 | -8 |
| withheld pointers | 35 | 0 | -35 |
| spec ship coverage | 0.0 | 1.0 | 1.0 |

Note the registry's own `verified` flag covers 2
of 13 assets (hash-ok but `present_unverified`);
`spec_ship_coverage` = ships with a binding, **not** ships whose physics verified (that is 1/4), see D6/D2.

## ④ Figures

See `fig_m20506_assetlib_audit.png`:
(a) pointer tri-state per kind, before vs after; (b) KCS volume error vs grid spacing against the 1% and 3% lines,
with the auditor's independent rerun marked; (c) fleet scope — one green VERIFIED, three red WITHHELD with error magnitudes.
The subject's own `artifacts/assetlib/assetlib_growth.png` is copied here unmodified as secondary evidence.

## ⑤ Notes / limits of this audit

- The audit ran on the platform host, where the geometry store `/root/tp-agent/assets` exists.
  **A clean clone cannot reproduce any of these numbers**: the store root is recorded in the registry, but the
  20.17 MB of publisher bytes are deliberately not redistributed (license note in the registry) and are not in git
  (`git ls-files` carries 20 KB of registry JSON and no hull IGES). This is a design choice, honest and documented —
  but it is the reason L4 reproducibility fails (D3).
- The subject's branch was not on `origin` at audit close (0 remote refs, not in `main`); all pointers above resolve in
  `/root/tp-agent/work/m20465` @ `9de9850`. If the gate resolves evidence against `main` only, this bundle will look
  empty — the gate must resolve against the mission branch or after merge (D5).
- Not run by this audit: `build-registry` in write mode against the live store (it would rewrite a shared platform
  registry and change a snapshot history that is evidence of other missions); the tamper-detection paths were exercised
  through the test suite only, not by mutating live bytes.
- The subject's `present_unverified` status for 11 assets was accepted on its own semantics; the audit verified the
  hash-match claim (`verify` → 13/13 hashes agree) but did not re-derive why the registry grades only 2 as `verified`.

## ⑥ Error judgement (explicit)

| quantity | reference | got | error | <3%? |
|---|---|---|---|---|
| KCS displacement (dx=0.5) | 52030 m³ | 52214.08 m³ | +0.354% | **YES** |
| KCS displacement (dx=0.25) | 52030 m³ | 52134.65 m³ | +0.201% | YES |
| KCS displacement (dx=0.125) | 52030 m³ | 52075.28 m³ | +0.087% | YES |
| KCS Cb (dx=0.5) | 0.6505 | 0.65280 | +0.354% | YES |
| KCS LCB (dx=0.5) | −1.48 %Lpp | -1.506 | -0.0256 pp | YES |
| KVLCC2 / JBC / ONRT | 312622 / 178369.9 / 8500.78 m³ | no valid got value | not computable | **NO — WITHHELD**, no <3% claim possible |

Error decreases monotonically with dx (0.125m:+0.087%, 0.25m:+0.201%, 0.5m:+0.354%) — first-order-looking
convergence toward the published value, so the +0.354% at the coarsest grid is discretisation, not a wrong scale.

## ⑦ Level decision

**L3 — system integration.**
Substance behind the level: 1470 LOC of real code (store / registry / specs / consumer / CLI) whose data path is
bytes → SHA-256 → registry → spec pointer → hull rebuild; a runnable CLI the auditor executed; 14 tests (43 with
`fleet_baseline`, no regression); the consumer is reached from the main fleet build path
(`fleet_baseline.specs_from_json` → `tensorcad/fleet/build.py`); one real external anchor met at 0.354% and reproduced
independently.

**Why not L4:** (D1) no benchmark case registered for the anchor and no mission gate key, so the sole judge has nothing
to confirm; (D2) fleet scope 1/4, with three hulls honestly WITHHELD and one (KVLCC2) failing inside the inherited
field reconstruction; (D3) the bytes that carry the numbers are outside version control, so the case is not
reproducible from the repository — reproducibility is the L4 leg this mission cannot satisfy on its own.

**Why not L2:** the work is not a standalone module that runs alone; it changes an existing consumer
(`fleet_baseline`) and is exercised through the fleet build path, i.e. it has integration callers with real behaviour.

Findings D1–D7 are listed in `verdict_m20506.json`. None of them is a fraud finding: the subject disclosed every
failure (m18488 archive drift included) rather than hiding it, which is what the honest-FAIL rule asks for.
