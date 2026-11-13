# Audit m17325 — GAP-CAD-T1-BASELINE-FLEET (subject m17293, TensorCAD)

**Verdict: maturity L2** (independent functional implementation, independently
reproduced) · disposition **request changes** · gate state left for the machine
(`tensor_maturity_verdicts` id written by this audit).

- subject mission: m17293 ("GAP-CAD-T1 返修1：官方 IGES→水密船体网格生产者")
- subject worktree: `/root/tp-agent/work/m17293` @ `404be46`
  (branch `agent/m17293`, `origin/main` = `9ab0897`)
- audit worktree: `/root/tp-agent/work/m17325` (branch `agent/m17325`)
- anchor goal: a 6-hull authority geometry set (KCS/KVLCC2/ONRT/…)
- independence: no tensorcad.* import; volume/topology/geometry recomputed from shipped STL bytes with stdlib only; reference re-extracted from cached publisher HTML, not from the subject's fleet_specs.json

## 1. What was audited (case definition)

The subject claims an IGES (B-spline patch) → half-breadth field → watertight
mesh → displacement pipeline, validated on the KCS hull from the NMRI
Tokyo-2015 / SIMMAN-2020 workshop archive, with five other hulls honestly
WITHDRAWN. The published reference is the workshop page cached in the subject
tree; the auditor re-extracted the reference numbers from the cached bytes and
compared them with the subject's own recorded snippets.

| quantity | published (re-extracted by auditor) | page sha256 match |
|---|---|---|
| L_PP [m] | 230.0 | True |
| B_WL [m] | 32.2 | True |
| T [m] | 10.8 | True |
| displacement ∇ [m³] | 52030.0 | True |
| C_B (−) | 0.6505 | True |
| LCB [%L_PP] | -1.48 | True |

Reference provenance is sound: `sha256(page) = 68d9c24600cb699b16394fc5c32762ebc875c856d14ef8460e1d0ebf4e7c78a4` equals the
sha256 recorded in `fleet_specs.json` (True), and the auditor's independent
re-extraction returns exactly the subject's values (True). Reference
values are machine-extracted from the publisher's page, not hand-typed.

## 2. How (auditor method, no subject imports)

`runs/m17325/audit_independent.py` (stdlib only) recomputes everything from artifact
bytes: divergence-theorem volume over the shipped STL triangles, an undirected-edge
census (every edge must be used exactly twice for a closed 2-manifold), a mid-ship
section area by plane clipping + breadth integration, a design-waterline section, and a
regex re-extraction of the publisher page. Test suite re-run in the subject tree.

```
python3 runs/m17325/audit_independent.py --subject /root/tp-agent/work/m17293 \
        --out runs/m17325/audit_evidence.json
cd /root/tp-agent/work/m17293 && python3 -m pytest tests/test_iges_baseline.py tests/test_fleet_gate.py -q
```

## 3. Results (machine-interpolated)

### 3.1 Primary metric — displacement volume, KCS

| metric | subject claim | auditor recomputation | agreement |
|---|---|---|---|
| ∇ mesh [m³] | 52187.9792 | 52187.9806 | Δ = 1.44e-03 m³ |
| rel. error vs ∇_pub (%) | 0.303631 | 0.303634 | 2.77e-06 pp |
| watertight | True | True | agree |
| edge-manifold | True | True | agree |
| triangles | 52320 | 52320 | agree |

Topology census: 78480 unique edges, all used exactly twice
(78480), 0 boundary edges,
0 non-manifold edges.

### 3.2 Secondary published checks (geometry plausibility)

| quantity | published | auditor from STL | Δ |
|---|---|---|---|
| C_B = ∇/(L_PP·B_WL·T) | 0.6505 | 0.65247 | +0.303 % |
| C_M (mid-ship section, x=L/2) | 0.985 | 0.98441 | -0.059 % |
| design-waterline span at z=T [m] | L_PP 230.0 | 228.52 (x ∈ [-2.23, 226.29]) | -1.48 m |
| overall LOA from STL [m] | 232.5 | 239.00 | +2.80 % |

The integral quantities (∇, C_B) agree to 0.303 %, the mid-ship
section to 0.06 %, so the mesh is the KCS hull and not a
displacement-matched blob. The **local fore-body shape is the weak point**: the mesh's
x-extent runs to 239.0 m (+2.80 % of LOA) while the
design-waterline section stops at x = 226.29 m, i.e.
3.71 m short of FP; 328 vertices sit beyond x = 232.6 m
(up to x = 236.77 m at z = 7.7–7.9 m from the file), which is a field-interpolation bulge at
the stem, not hull geometry.

### 3.3 Honest-FAIL section preserved (not overwritten)

| hull | subject verdict | auditor check | disposition |
|---|---|---|---|
| KCS | VERIFIED | volume/topology/geometry reproduced independently | upheld |
| JBC | WITHHELD | IGES control-net layout guard fired ("candidates score within 4x") | honest, no reference published in cached page |
| ONRT | WITHHELD | STL bbox = 1.000 × 0.1218 × 0.100 m, ∇ = 0.005595 m³ vs published 8500.78 m³ | honest: model-scale bytes, no scale factor obtained |
| KVLCC2 | (no ref) | cached page re-scanned by auditor: no "volume"/"∇"/0.8098 token | claim "no published table reachable" is true, not laziness |
| DTC, Series 60 | NOT_SUPPORTED | — | stated as such in the subject README |

Anchor census: 1/6 hulls VERIFIED
(KCS); verdict census in `benchmark.json` = {'VERIFIED': 1, 'WITHHELD': 2}.
The board's framing that "6 hulls succeeded" is **not supported by the artifacts**; the
subject's own report is honest about this, the board summary is not.

### 3.4 Test evidence re-run by the auditor

`python3 -m pytest tests/test_iges_baseline.py tests/test_fleet_gate.py -q` →
**11 passed, 1 xfailed** (own log `runs/m17325/audit_run.log`). The tests assert
invariants, not mere execution: watertightness, ∇ < 1 % of the published value,
C_B ≈ 0.6505 ± 0.01, LCB ≈ −1.48 ± 0.15 % L_PP, closed-ring area vs πab, and the two
independent volume integrators agreeing to 2 %. The xfail is a disclosed pinch-to-a-line
limitation of the analytic fixture, not a hidden failure.

## 4. Findings (severity-ranked)

| id | severity | finding | evidence |
|---|---|---|---|
| F1 | **blocking** | None of the deliverables is in git: `tensorcad/iges.py`, `fleet_baseline.py`, `fleet_gate.py`, `tools/baseline_fleet/*`, the two test files and `artifacts/baseline_fleet/*` are untracked in the subject worktree and absent from `HEAD` and `origin/main` (`any_file_in_origin_main = False`, `any_file_tracked_in_head = False`). The capability exists only as loose files in one worktree. | `audit_evidence.json → E_landing.files` |
| F2 | material | No public API / no integration: `tensorcad/__init__.py` and `tensorcad/cli.py` mention neither `iges` nor `fleet_baseline` → the producer cannot be called through any stable interface (L3 not reachable). | `E_landing.public_api_wiring` |
| F3 | material | Two documented secondary numbers contradict the subject's own machine artifact: docs claim C_B 0.6516 (+0.17 %) and LWL 230.25 m, while `KCS_manifest.json` records `derived.cb = 0.6524737` and `derived.length_waterline_m = 229.0`; the auditor's independent value is C_B = 0.65247 (+0.303 %). The doc lines are hand-typed and violate the zero-hand-copy rule. | `docs/baseline_fleet_build.md:57,102`; `benchmarks/verified/baseline_fleet/README.md` |
| F4 | material | Fore-body fidelity: LOA overshoot +2.80 % and waterline section 3.71 m short of FP. Fine for displacement, not for stem/appendage-sensitive work. | §3.2 |
| F5 | info (+1) | Honest-FAIL discipline is respected: 5 of 6 hulls are withheld/not-supported with named blockers, and the NaN/zero-volume case is not dressed up. | §3.3 |

## 5. Maturity call

| rung | criterion | met? | evidence |
|---|---|---|---|
| L0 | document only | yes | docs + report exist |
| L1 | math/unit verification | yes | 11 passed / 1 xfail, closed-form ring area, two-integrator agreement |
| L2 | real functional implementation (input → compute → output) | **yes** | IGES reader (672 lines) + producer (747 lines) + gate (384 lines); auditor rebuilt ∇ from the shipped STL to 52187.98 m³, watertight |
| L3 | system integration (solver/platform) | no | F1/F2: untracked, no export, no CLI, no caller |
| L4 | benchmark validation at anchor scope (complete case + error analysis + reproducible) | no | 1/6 hulls verified; no `tensor_benchmark_cases` row existed before this audit (registered by it, id in DB); error analysis contaminated by F3 |
| L5 | research/engineering result | no | — |

**Level = 2.** The single KCS case would support an L4-grade *pipeline* claim if the code
were landed and integrated, but the anchor's scope is a 6-hull fleet and the deliverables
are not in the repository. Level is function substance, not report length.

## 6. Error judgment (explicit)

| quantity | reference | measured | error | < 3 %? |
|---|---|---|---|---|
| KCS displacement ∇ | 52030.0 m³ (publisher) | 52187.98 m³ (auditor) | +0.3036 % | **yes** |
| KCS C_B | 0.6505 (publisher) | 0.65247 (auditor) | +0.3034 % | **yes** |
| KCS C_M mid-ship | 0.985 (publisher) | 0.98441 (auditor) | -0.0594 % | **yes** |
| KCS LOA | 232.5 m (publisher) | 239.00 m (auditor) | +2.80 % | **yes (integral), but the stem profile is locally wrong — see F4** |
| ONRT ∇ | 8500.78 m³ (publisher) | 0.005595 m³ | -99.9999 % | no — correctly reported as WITHHELD, not a silent pass |

## 7. Scope, limits, checks not run

- Not run: the subject's `fetch_sources.py` (network fetch) — the cached bytes were used, and
  their sha256 was verified instead. Not re-run: `build_fleet.py` end-to-end (≈36 s/hull; the
  shipped STL + manifest were audited instead, which is the stronger check for the claim).
- Not verified: that the workshop archive licences redistribution (the subject states bytes are
  not redistributed and ships only derived numbers/meshes).
- The auditor did not fix the subject tree (audit, not rework): F1–F4 are handed back.
- Independence limits: the same machine and the same IGES bytes were used; a share-bug in the
  IGES reader could only shift both the subject and the auditor if the auditor had imported it
  — it did not (artifact-level recomputation only).
