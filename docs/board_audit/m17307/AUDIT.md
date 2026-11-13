# m17307 board audit — m17251 `GAP-CAD-T1-OFFSET-LOFT` (TensorCAD)

**Disposition: approve at L2 (独立功能实现).**  L3 (system integration) and L4 (registered
benchmark) are **refused with named gaps**, not because the work is weak but because those
rungs require things this delivery does not have (a downstream consumer; a benchmark case in
`tensor_benchmark_cases`).  The engineering itself is real, reproducible and unusually
well-instrumented — it is an honest L2 with L4-grade rigour in its own self-checks.

Audited subject: mission **m17251**, branch `agent/m17251`, commit **5bc0db7**
(`offset_loft: offset table -> NURBS loft with Simpson baseline, provenance digests and
station/waterline sensitivity`), 6 files committed + 10 on-disk artifacts.  Snapshot with
per-file sha256: `runs/m17307/subject_m17251/MANIFEST.json`.

## What was checked (independently, not from the report)

| check | method | result |
|---|---|---|
| code is real, not a shell | read `tensorcad/offset_loft.py` (1282 lines) + `offset_forms.py` (284 lines): offset table/provenance bundle → `parametric_hull.loft_nurbs_surface` → Green-theorem section integrals over the patch's own basis → volume/LCB/VCB/A_wp/wetted area, plus a `NURBSSurface.derivatives` slow cross-check and a tessellated mesh route | pass |
| tests | `pytest tests/test_offset_loft.py -q` → **18 passed** (17.5 s); with `tests/test_loft.py` → **49 passed** (27.2 s), auditor rerun | pass |
| reproducibility | reran `examples/offset_loft/run_example.py` in the subject worktree: `results.json` sha256 `c1aeb70e…` **before and after → byte identical**; worktree clean afterwards | pass (`runs/m17307/logs/`) |
| baseline numbers not fabricated | re-implemented composite Simpson from the CSV alone (uniform 1/3 + 3/8 tail; unequal-spacing quadratic rule), **no import of the audited module** | reproduces the reported `table_simpson` to **≤9.3e-16** relative (all 3 hulls) |
| independent quadrature | not-a-knot cubic spline in z per station, then in x (different rule) | Wigley: **3.2e-16** vs analytic `4BLT/9`, **0.0** vs analytic `2LB/3`; Series 60: 5.1e-5 vs the reported NURBS volume; KVLCC2: 2.4e-3 (the hull's own 41-station discretisation error) |
| geometry sanity (Cb/L/B) | Cb from the landed tables | Wigley 0.444444 (`4/9` ✓), Series 60 0.600702 (parent 0.600), KVLCC2 0.810160 vs published 0.8098 → 0.04 % |
| stated acceptance | loft volume **and** waterplane vs the table's own Simpson integral, 0.5 % | worst **0.157 %** (volume, KVLCC2) / **0.237 %** (waterplane, KVLCC2); loft-mesh route worst 0.128 % → **passes with margin**, and both are reproduced from `results.json` |
| station-count sensitivity | read `results.json.station_sensitivity` | KVLCC2 11 / 15 / 21 / 41 / 81 stations → 1.604 % / 0.498 % / 0.152 % / 0.123 % / 0.089 % vs the 81-station table: the 0.5 % band needs ≥15 stations, and that is stated in the delivered report |
| integrity claim ("tampered ordinate rejected") | mutated one interior ordinate (+5 %) and deleted a provenance sidecar on copies of two landed bundles, ran the audited reader | 4/4 probes rejected with `digest mismatch` / `provenance sidecar missing` — the claim is true, not just test-covered (`runs/m17307/tamper_probe.json`) |
| registry判据 | `select count(*) from tensor_benchmark_cases where anchor_slug='GAP-CAD-T1-OFFSET-LOFT'` | **0 rows**; `tensorbench_results` for mission 17251 = **0 rows**; `/tensor/bench/l4-verdict/17251` → `CLAIMED_NOT_PROVEN`. There is no registered ref/tolerance to score against. |
| integration | `git ls-tree HEAD tensorcad/offset_loft.py` in `/root/tensor-ecosystem/TensorCAD`; grep for consumers | empty (branch unmerged); only re-export in `tensorcad/__init__.py`, the companion `offset_forms.py` of the same commit, and the example script. No solver/CLI/API/platform caller. |

## Severity-ranked findings

1. **[material, scope — not fraud] The three "船型值表" are the library's own parametric
   reconstructions, not issuing-body tables.** Every provenance sidecar carries
   `official: false`; the KVLCC2 sidecar states verbatim *"NOT the official KRISO offset
   table … Cb/LCB agreement is a representation consistency check"*. Consequence: the
   acceptance measures **interpolation + quadrature consistency** (plus exact analytic
   identities for Wigley), not fidelity to real ship geometry. The stated motivation of the
   gap — ship-platform geometry inputs are placeholder STLs / KVLCC1≡KVLCC2 — is therefore
   **not** solved: a pipeline now exists, real geometry data does not. Downstream users must
   not read "三型船落地" as "three real ships' offsets ingested".
2. **[material, level-capping] No benchmark registration and no external anchor.**
   `tensor_benchmark_cases` has no row for the anchor → the platform's L4 rule is
   unsatisfiable and `l4-verdict` reports `CLAIMED_NOT_PROVEN`. The available anchors are
   analytic (Wigley `Cb=4/9`, `Cwp=2/3`, `VCB=5T/8`, reproduced to 3e-16) and published
   *particulars* (KVLCC2 L/B/T/Cb/LCB, used for calibration only); there is no measured
   table or published hydrostatic value to compare a `got` against.
3. **[material, level-capping] No system integration.** Re-export in the package `__init__`
   is API surface, not a consumer: no hydrostatics/resistance path, CLI or platform endpoint
   calls the loft, and the branch is not merged, so nothing downstream can call it today.
4. **[minor] Platform report 16253 is the auto-generated `report.md`**, i.e. tables only: it
   has no explicit six-element structure and does not state the "是否<3%" answer (the numbers
   do answer it: worst 0.237 % vs the 0.5 % acceptance, so <3 % holds for every quantity on
   the table-Simpson reference; the clustered-waterline p=1 cases reach 1.77 % vs closed form
   and are disclosed in the same table). Reporting discipline gap, not a numeric defect.
5. **[minor, disclosure only] The acceptance reference is the coarsest route in the study.**
   With 17 waterlines clustered as `ζ^(1/10)` and first-waterline exponent p=1, the table's
   own Simpson integral is less accurate than the loft: KVLCC2 p=1 shows dV 0.230 % vs table
   but 1.153 % vs closed form. The delivered report prints both columns, so the reader cannot
   be misled — but a user comparing only to "table Simpson" would over-trust p=1 tables.
   Recommend the ingestion path reject first-waterline clustering that coarse.

## Checks not run / limitations of this audit

- No independent *geometry* truth for Series 60 / KVLCC2 (none exists in this sandbox);
  Cb checks bound the table, they do not certify the hull shape.
- Mesh route (`mesh_integrals` → `parametric_hull.compute_hydrostatics`) was verified only
  through byte-identical reproduction of `results.json`, not re-implemented from scratch.
- Performance/memory at production table sizes (e.g. 200 × 60) was not measured; the module
  is pure Python.
- `git ls-tree` was run against `/root/tensor-ecosystem/TensorCAD` HEAD at audit time; if the
  branch is merged later, finding 3 changes (integration still absent, merge provenance
  changes).

## Why L2 and not L3/L4

L2 requires real code with input→compute→output: met (`offset_loft.py`, entry CLI, 12 run
artifacts, 18 tests). L3 requires a caller in the solver/platform: absent, so claiming it
would be inflation. L4 requires a benchmark case with a registered ref/tolerance: 0 rows, and
the in-repo "ref" is the same table, so it cannot substitute. Declaring L2 is not a demotion
of the engineering — the study's internal verification (three independent integral routes,
station/waterline sensitivity, digest-checked provenance, tamper probes) is L4-grade in form;
what is missing is *external* evidence.

**What would earn L4:** register a case `CASE-CAD-OFFSET-LOFT` in `tensor_benchmark_cases`
with a published anchor (e.g. KCS or KVLCC2 official hydrostatics: ∇, LCB, A_wp) and a
tolerance, ingest at least one issuing-body offset table with its licence/sha, publish
ref-vs-got for all three quantities, and record the rows through `/tensor/bench/results`.
**What would earn L3:** a consumer — e.g. the ship-platform hydrostatics path or a
`tensorcad` CLI — calling `loft_offset_table` on a table in the pipeline.

## Reproduce this audit

```
python3 runs/m17307/independent_check.py --subject /root/tp-agent/work/m17251
python3 runs/m17307/tamper_probe.py
python3 docs/board_audit/m17307/build_evidence.py     # evidence.json + audit_figures.png
python3 docs/board_audit/m17307/insert_verdict.py     # row into tensor_maturity_verdicts
```

Verdict row: `('TensorCAD', 'm17251 GAP-CAD-T1-OFFSET-LOFT', 17307, 2, evidence, 'board-audit')`.