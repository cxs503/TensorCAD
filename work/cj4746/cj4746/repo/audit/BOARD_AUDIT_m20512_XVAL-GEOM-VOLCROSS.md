# Board audit — m20512 `XVAL-GEOM-VOLCROSS` (repo TensorCAD)

Audit mission: **m20544** · auditor: `board-audit` · subject worktree:
`/root/tp-agent/work/m20512` · branch `agent/m20512` · HEAD `0d829a6`
("TensorCAD: cross-tool volume cross-check reaches the product surface",
19 files / +4480 −2, worktree clean, artifacts tracked).
Anchor = *几何体积跨工具互证*: the platform's own NURBS/mesh volume integral vs
independently implemented engines.

Referee mandate: `tensor_benchmark_cases` is the sole numeric referee
(`SELECT case_name, metric, tolerance … WHERE anchor_slug='XVAL-GEOM-VOLCROSS'`).

---

## 0. Verdict

| item | result |
|---|---|
| level granted | **L3 — 系统集成** (real code + entry + run artifacts + a real caller in `tensorcad.fleet_gate`) |
| L4 | **refused**, one item short: the referee table holds **0 rows** for this anchor (132 cases / 93 anchors total). Substance for L4 (complete case + error analysis + reproducible) is *verified present* — the refusal is registration, not mathematics. |
| L2 / L1 | granted (real code, 40 tests pass, entry exit codes 0/1/2 observed) |
| L5 | refused (no research-grade result beyond the above) |
| disposition | **approve with notes**, L3 confirmed; re-gate to L4 after one `INSERT` (SQL in §5) |

Everything below is either read from the subject's committed artifacts or produced
by `audit/m20544_audit.py` (this worktree); raw numbers in
`audit/m20544_evidence_raw.json`.

---

## 1. Evidence three-piece set — does it exist?

| claimed | found | how verified |
|---|---|---|
| run artifacts | `artifacts/xval_geom_volcross/{results.json,results.md,figures/*.png}` (3 PNG, 902-line JSON, 146-line MD) + `benchmarks/verified/volume_cross_engine/{result.json,gates_verified.json,figures/…png}` | sha256 inventory, all 18 paths present (`inventory` block) |
| benchmark ref vs got | in-repo gates G1–G5 recomputed **from `result.json` by this audit** and equal to `gates_verified.json` to the last bit: G1 2.682978e-08 ≤ 1e-06; G2 3.210584e-04 ≤ 3e-03; G3 1.598721e-16 ≤ 1e-09; G4 2.020500631935433 ∈ [1.5, 2.6]; G5 [] | `gates` block, `matches: true` ×4 |
| **platform referee** | **0 rows** for the anchor → no `ref_source`, no `tolerance` registered by the board | `referee` block (psycopg, direct SQL) |
| repro_cmd | `python benchmarks/verified/volume_cross_engine/run_benchmark.py` → **exit 0**, gates byte-identical to the shipped `result.json`, live rows identical, `40 passed` in 10.07 s for the two test modules | `repro` block, run in disposable copy `/tmp/m20544_repro` |
| entry | `python -m tensorcad.cli volume-xval wigley` → 0; `… wigley --stations 61 --levels 13` → **1** (honest FAIL, `reference_error −0.329 %`); `… mesh /tmp/does_not_exist.stl` → **2** | observed return codes 0/1/2 |

## 2. Independent recomputation (own arithmetic, no `tensorcad` import)

| quantity | my arithmetic | artifact | agreement |
|---|---|---|---|
| Wigley analytic V = 4LBT/9 | 2.7777777777777777 m³ | 2.7777777777777777 | exact |
| Wigley finest rung vs analytic (161×41) | −3.210584e-04 (−0.0321 %) | −0.0321 % | Δ = 0.0 |
| per-rung deviation, all 4 rungs | recomputed | recomputed | max abs Δ = **0.0** |
| per-rung engine spread | own (A−B)/A | `engine_agreement` | max Δ = 5.37e-08 (float32 STL budget) |
| measured deficit order on the artifact ladder | **1.9797** | live ladder order 2.0205 (different level rule `n//4+1` vs `n//4+2`) | both inside the pre-declared window 1.5–2.6 |
| Richardson limit from my own fit | 2.7777922 m³ → +5.20e-06 vs analytic | ships +… same construction | consistent |
| asymmetric-transom battery: engine A vs analytic | −0.8262 % (7.58438617 vs 7.647567686), watertight, genus 0 | −0.8262 % | exact |
| KCS full scale: Cb from the reported geometry volume / LBT | **0.650904** (published 0.6505, +6.20e-04) | report prints +0.0623 % vs 52030 m³ | value disclosed, label imprecise → **F1** |
| KCS modelscale | Cb from volume/LBT = 0.5554; `reference_kind = none (WITHHELD)` | same | honest |

**Geometry gate (Cb/L/B) passes**: LBT 230·32.2·10.8 = 79 984.8 m³ = artifact value;
KCS 0.6509 vs 0.6505 (+0.062 %), modelscale 0.5554 vs 0.555. No sign of a shape error.

## 3. Findings, severity-ranked

**F1 — medium (claim accuracy, geometry anchor).** The mission labels the full-scale
hull "Cb calibrated to 0.6505", but the hull that is actually reported has
Cb = 0.650904. Cause, measured here: `calibrate_form` (tensorcad/hull_nurbs.py:548)
converges at quadrature order **(10,6)** — volume 52 030.088 m³, 4.7e-07 from the
52 030.1124 m³ target — while the reported pipeline number uses **(24,16)** =
52 062.390 m³; the order gap is **6.2e-04 relative**, i.e. 600× the calibration
tolerance of 1e-6. `hull_nurbs.py:557-559` states "the tests verify that the order
dependence is far below the calibration tolerance" — false for the KCS-scale hull
(this audit measured (10,6)=52 030.088, (16,12)=52 063.197, (24,16)=52 062.390,
(32,24)=52 062.433, (48,32)=52 062.405) and contradicted by the repo's own
`tests/test_hull_nurbs.py:271-280`, whose docstring honestly says 7e-4 and whose
assertion allows 2e-3. Impact: the row called a "calibration identity" is not an
identity to 1e-6; the +0.0623 % vs the published 52 030 m³ is the calibration-order
residual. No number is hidden (the report prints +0.0623 %), but the label and the
docstring overstate. Fix: calibrate at the reporting order (or report the achieved
Cb = 0.65090 in the row name).

**F2 — blocking for L4 (process, not maths).** `tensor_benchmark_cases` has **0 rows**
for `XVAL-GEOM-VOLCROSS`, so the referee cannot confirm any got-vs-ref pair. The
in-repo gates are self-authored (tolerance 0.3 % chosen by the subject), which is why
a referee row is required before L4. One `INSERT` fixes it — §5.

**F3 — low (artifact metadata).** `benchmarks/verified/volume_cross_engine/result.json`
records `"numpy_stl": "absent"` although numpy-stl **4.1.3** executed engine B in that
very run (`stl` exposes no `__version__`, `run_benchmark.py:_version` swallows it).
A machine reading the artifact concludes the second engine was missing, contradicting
report §2's toolchain list. Fix: `importlib.metadata.version("numpy-stl")`.

**F4 — low (gate strength).** G5 "committed rows hold" re-reads the committed
`results.json` and checks spreads/band on the *stored* numbers; it does not rebuild
hulls, so G5 is self-consistency, not reproduction. Reproduction is real but comes
from a command outside the declared `repro_cmd`: this audit ran
`tools/volume_xval_report.py --out …` on a disposable copy and got **`results.json`
identical (schema + every numeric leaf)** and `results.md` identical except the
generation stamp. Recommend folding that command into the benchmark run.

**F5 — low (report vs artifact).** Posted report §3 renders the Wigley volumes of the
multi-hull table at one decimal — `2.7` for 2.7231650437231827 (0.85 % rendering
error) and `2.8` for 2.763654038838464 (1.32 %), 2.774218 (0.93 %), 2.776886 (0.83 %)
— so the "程序化注入零手抄" claim holds for §3's first table and the artifacts but
not for that table. KCS rows are fine (52042.8, 7.1e-07).

**F6 — info.** `docs/volume_xval.md` title promises "独立实体引擎" (B-rep solid); the
second engine is a triangle-mesh engine (numpy-stl/Mirtich). The body discloses this
and the gmsh/OpenCASCADE route as `WITHHELD`, so the title is the only overstatement.

## 4. Checks performed / not run

Performed: sha256 inventory of 18 paths; referee SQL (132 cases / 93 anchors / 0 for
this anchor); independent analytic + deviation + order + Richardson recomputation;
KCS Cb/LBT geometry gate; quadrature-order sweep (10,6)→(48,32); gate recomputation
from `result.json`; DB report (`tensor_reports.id=18660`, 1 figure, 6 sections, all 8
spot-checked tokens present) vs artifacts; full `repro_cmd` replay in `/tmp/m20544_repro`
(exit 0, gates identical); `tools/volume_xval_report.py` regeneration (identical);
`pytest tests/test_volume_xval.py tests/test_volume_xval_entry.py` (40 passed);
CLI exit codes 0/1/2; `fit_order` on a zero/constant ladder (returns `None`, no crash).

Not run / unknown (marked, not passed): no B-rep (gmsh / OpenCASCADE / cadquery)
cross-check — not installable here, so "third independent engine" remains open;
trimesh's and numpy-stl's internals were not re-derived, only cross-compared;
the KCS published 52 030 m³ source page was not re-fetched by this audit (taken as
quoted by m20511-era sources); the JBC −18.3 % incident premise was not re-verified.

## 5. Making it L4 (for the board / repo owner, not applied by the auditor)

```sql
INSERT INTO tensor_benchmark_cases (anchor_slug, case_name, ref_source, metric, tolerance, repo)
VALUES ('XVAL-GEOM-VOLCROSS','CASE-CAD-VOLCROSS-WIGLEY-ANALYTIC',
        'V = 4LBT/9 analytic (tensorcad.wigley.analytic_volume); manifest ref 2.7777777777777777 vs got 2.7768859489811852 at 161x41',
        'displacement_volume_m3 rel error at finest rung', '3e-3','TensorCAD'),
       ('XVAL-GEOM-VOLCROSS','CASE-CAD-VOLCROSS-ENGINE-AGREEMENT',
        'two independent engines (trimesh 5.0.0, numpy-stl 4.1.3) on identical triangles; got 2.68e-08',
        'engine-to-engine relative spread', '1e-6','TensorCAD');
```

Then the same evidence bundle (`benchmark` key with ref/got/error, `repro_cmd`) is
re-gateable without touching the code.