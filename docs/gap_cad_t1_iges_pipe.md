# GAP-CAD-T1-IGES-PIPE — IGES/STEP → watertight mesh import pipeline (rework round 1)

**Repo** TensorCAD · **Mission** m17426 (rework of m17333 ← m17250) · **Anchor** GAP-CAD-T1-IGES-PIPE
**Declared maturity: L2** (impl files + entry + run artifacts). **L3/L4 NOT claimed**: the
integration entry and the gold benchmark harness exist, but the gold regression is **not
numerically closed** in this round (honest FAIL, see §6).

## 1. Case

Ship-platform geometry input is blocked: of its 8 STL files six are 304-triangle placeholders,
KVLCC1 ≡ KVLCC2 byte-identical, `/generate` only emits boxes. The only industrial input is a
trimmed B-rep (IGES/STEP). TensorCAD's `geomquality.io` reads **untrimmed** NURBS patch sets
and explicitly refuses entity 144; the previous round added `geomquality/iges_trim.py` +
`trim_mesh.py` + `import_cli.py`, but that round was judged L1 because (a) 0/4 faces of the one
real third-party file tessellated, (b) no run artifacts committed, (c) no benchmark row.

## 2. Method (what changed in this round)

| change | file | evidence |
|---|---|---|
| **Defect fix — IGES entity 100 field layout.** The reader assumed 12 fields (plane point, centre, start, end); IGES 5.3 §5.11 defines `ZT, X1, Y1, X2, Y2, X3, Y3` — one planar Z, then centre/start/end in that plane, normal +Z rotated by the entity transform. | `tensorcad/geomquality/iges_trim.py` | KP505 CATIA file: trimmed faces read **4/8 → 8/8**, refused-curve list empty (was 8 entries, all "type 100 has 9 parameters, field 10 is missing") |
| Tessellation-density knob plumbed end-to-end (`u_samples`, `v_samples`) so the discretisation error of the patch-set route is a measurable, closable term instead of a fixed bias. | `tensorcad/geomquality/trim_mesh.py:import_cad`, `import_patch_set` | convergence sweep runs at 12x8 … 96x48 (`runs/gap_t1_iges/evidence.json`) |
| **Integration entry** `tensorcad.geometry.load_cad_geometry(path)` — documented public API in the mesh module, delegates to the pipeline, returns `(mesh, report)`. | `tensorcad/geometry.py` | benchmark harness imports it; caller chain `geometry.load_cad_geometry → trim_mesh.import_cad → gates` |
| **Gold loft harness**: `DTMB5415_Hull_exact.stl` → transverse sections → separable B-spline interpolation → hull patch + ruled deck + 2 ruled end caps → entity-128 IGES, plus an independent Gauss-Legendre quadrature (`V = 1/3 ∮ S·(S_u×S_v)`) as the reference integral. | `benchmarks/verified/iges_import/gold_loft.py` | `artifacts/gap_t1_iges/inputs/DTMB5415_gold_loft.iges` (132,354 B, sha256 in provenance) |
| Benchmark driver, one-command repro. | `benchmarks/verified/iges_import/run_benchmark.py` | `runs/gap_t1_iges/evidence.json` |

## 3. Results (all numbers read back from `runs/gap_t1_iges/evidence.json`)

Repro:
```bash
PYTHONPATH=$PWD python benchmarks/verified/iges_import/run_benchmark.py --out-dir runs/gap_t1_iges
```

Inputs written and hashed in `artifacts/gap_t1_iges/provenance.json`:
`series60_cb060.iges` (149,121 B), `series60_cb065.iges` (194,724 B), `kvlcc2_class.iges`
(131,301 B), `DTMB5415_gold_loft.iges` (132,354 B); gold mesh
`DTMB5415_Hull_exact.stl` 8,246,384 B, 164,926 triangles, sha256 recorded.

**KP505_modelScale_250mm.igs** (real CATIA V5-6R2012 export, 1,042,550 B,
sha256 0bb54ae7618e8d8b1b7c23de8d384e04b0465071e34fe32e1d606ddbb5df462f), trimmed route:

| metric | before this round | after |
|---|---|---|
| trimmed faces parsed | 4 / 8 | **8 / 8** |
| faces tessellated | 0 / 8 | 2 / 8 |
| triangles out | 0 | 252 |
| free edges | — | 256 (not watertight) |
| gates | G1 F, G2 T, G3 F, G4 T, G5 T, G6 F, G7 T | unchanged: G1/G3/G6 fail |

So the parser defect is fixed and verified; the *whole-file* watertight result for the
propeller is still **not** achieved — the remaining 6 faces are refused with measured
inversion residuals (the trim curves of those faces do not invert onto the reader's
`RevolvedSurface`/`NURBSSurface` within the file's own tolerance). This is a genuine,
measured refusal, not a silent skip.

## 4. Figures

None. Justification: the only converged quantitative run this round is KP505's face
accounting (a 2-row table); the gold regression numbers are invalid (§6), so a plot of them
would present unvalidated data as a result. Plot scripts exist (`run_benchmark.py` writes
per-file `*.report.json`) and a figure is the first task of the next round.

## 5. Notes / traps hit

* **The `DTMB5415_Hull_exact.stl` "gold" is not a validated DTMB5415.** Measured extents
  151.90 × 20.539 × 19.45 m, bow top at z = 16.17 m, 164,930 triangles (trimesh-merged),
  `is_watertight = False`. Published DTMB 5415 particulars are L_pp 142.0 m, B 19.06 m,
  T 6.15 m. The file is therefore used **only as the declared reference surface** (shape
  comparison), never as a volume/`Cb` gold — its divergence volume over the triangle soup is
  13,710.93 m³ while `trimesh.volume` reports 30,368.14 m³, i.e. the volume of an open shell is
  not a number worth regression-testing.
* Reusing an existing trap: the earlier audit already found that a green G1–G7 set says
  nothing about geometric accuracy (G5 only tests `volume > 0`). The quadrature reference was
  added for exactly that reason.
* `numpy 2.x` removed `ndarray.ptp`; `ndarray.reshape(...).view('<f4')` is still the fastest
  binary-STL read (164,926 triangles in 0.03 s).

## 6. Error judgement (explicit)

| key quantity | reference | measured | error | < 3 % ? |
|---|---|---|---|---|
| KP505 faces parsed | 8 faces in file | 8 | 0 | **yes** |
| KP505 faces tessellated | 8 needed | 2 | −75 % | **no (function not achieved)** |
| series60 volume vs quadrature | quadrature of the written patches | harness row errored (`KeyError: component_count`) then swept: 12x8 → +0.0011 %, 24x12 → +2.02 %, 48x24 → +3.18 %, 96x48 → +3.49 % | not monotone ⇒ **the reference or the sweep is wrong** | **unknown / not_run** |
| DTMB5415 gold volume | same quadrature | sweep −129.47 % … −126.51 % | sign/scale defect in the harness quadrature (`same_sense` guessed as all-True) | **no — harness defect, numbers withdrawn** |
| DTMB5415 section deviation vs gold | gold STL sections | not computed (row errored) | — | **not_run** |

Honest reading: this round **fixes a real parser defect with a reproducible before/after
(4/8 → 8/8 faces) and commits the missing artifacts, entry point and integration caller**, but
the gold benchmark (L4) does **not** pass: my quadrature reference and the sampling sweep do
not agree, so no volume/`Cb` error may be claimed. A next round must (1) fix the harness
`component_count` key and the `same_sense` signs of the loft (mirrored/cap faces), (2) show the
sweep converging monotonically to the analytic/qerature value, (3) extend the trimmed reader to
the remaining 6 KP505 faces (142 parameter-space curves).

## 7. Files

| path | role |
|---|---|
| `tensorcad/geomquality/iges_trim.py` | IGES parser (entity 100 layout fixed) |
| `tensorcad/geomquality/trim_mesh.py` | inversion, tessellation, repair, audit, gates, `import_cad` (density knobs) |
| `tensorcad/geomquality/import_cli.py` | CLI entry point |
| `tensorcad/geometry.py` | `load_cad_geometry()` — public integration entry |
| `benchmarks/verified/iges_import/gold_loft.py` | gold STL → NURBS loft → IGES + quadrature reference |
| `benchmarks/verified/iges_import/run_benchmark.py` | one-command benchmark / evidence driver |
| `artifacts/gap_t1_iges/inputs/*.iges`, `provenance.json` | inputs + provenance (path, sha256, bytes, source) |
| `runs/gap_t1_iges/evidence.json`, `*.report.json` | run artifacts |
