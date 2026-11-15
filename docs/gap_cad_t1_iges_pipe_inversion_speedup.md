# GAP-CAD-T1-IGES-PIPE — the import hot spot was rebuilding derivative control nets per evaluation

**Repo** TensorCAD · **Anchor** GAP-CAD-T1-IGES-PIPE · **Branch** `agent/m22464`
**Maturity: L3 maintained.** This round adds **no new claim** about geometry quality; it removes an
algorithmic defect in the surface evaluator that made the large publisher files (KVLCC2, ONRT,
the very files acceptance criterion #2 asks for) unreachable inside the 120 s test budget. The
KVLCC2/ONRT runs are still **not** through the gates, and the reason is diagnosed and measured
below — reported as an open item, not as a pass.

## 1. Case

The anchor's acceptance is 【管道 CLI+库函数；3 个真实 IGES 全过门；金标回归测试绿；依赖锁定】.
`docs/gap_cad_t1_iges_pipe_final_verification.md` closed criteria 1, 3, 4 and left criterion 2 as
**PARTIAL**: `KCS_hull_SVA.igs` imports (all gates, 2 s), while `kvlcc2.igs` (6.27 MB, 19 trimmed
faces) and `onr5613-…mergedR.igs` (3.37 MB) exceed the 120 s per-file budget with no result at all.
The earlier report called this "known performance limitation, O(n²) inversion"; that attribution is
**wrong** and is corrected here — the cost was not in the search strategy but in the evaluator.

## 2. Method — where the time actually goes

`cProfile` on `import_cad(kvlcc2.igs)` interrupted after 95 s (`work/prof_import.py`):

| function | cumtime | note |
|---|---|---|
| `trim_mesh.invert_point` | 81.8 s | 58 loop points → 171 multi-start lattices → 4 327 Gauss–Newton runs |
| `nurbs.NURBSSurface.derivatives` | **84.8 s** | 1 024 794 calls to `nurbs.derivative_net` |
| `trim_mesh._gauss_newton` itself | 0.8 s | the solver is cheap; the *surface* is expensive |

`NURBSSurface.derivatives` rebuilt the homogeneous control net and, for every requested
order, the differentiated control nets **from the raw control net on each call** — O(net) work per
evaluation where O(degree²) is enough. A trim loop inverting a 126×54 patch therefore paid the full
net-differentiation cost thousands of times per point.

## 3. Fix (two minimal, semantics-preserving edits)

1. `tensorcad/nurbs.py` — the homogeneous net and the u/v derivative nets are functions of the
   (immutable, `frozen=True`) surface only, so they are memoised once per surface in
   `_net_cache` (`_homogeneous_net`, `_levels_u`, `_levels_v`). Same arithmetic, built once.
2. `tensorcad/geomquality/trim_mesh.py::_gauss_newton` — for surfaces that expose `partials`, one
   `derivatives(u, v, 1)` call yields the point **and** both tangents; the old code called
   `evaluate` and then `partials`, paying for the value twice per iteration.

No tolerance, gate, seed lattice, iteration count or acceptance criterion was touched.

## 4. Results (programmatic: `python tools/bench_iges_inversion.py` → `artifacts/gap_t1_iges_inversion/perf.json`)

The script measures HEAD **and** the patched tree in one run (patch stashed to a backup, files
checked out from HEAD, measured, restored), so the table is generated, never transcribed.

| quantity (KVLCC2, first 2 faces, 4 loops, 78 loop points) | HEAD | patched | ratio |
|---|---|---|---|
| inversion wall time | 35.9 s | 10.5 s | **×3.4** |
| cost per loop point | 460 µs | 134 µs | **×3.4** |
| worst inversion residual | 1.366e+02 | 1.366e+02 | identical |
| KCS full import | 3.9 s | 2.1 s | ×1.8 |
| KCS triangles / vertices | 4 499 / 2 514 | 4 499 / 2 514 | identical |
| KCS gate verdicts | G1,G2,G5 pass | G1,G2,G5 pass | identical |

The residual and mesh identity rows are the point: the patch is a **cost** change, and the
measured geometry does not move by one vertex.

Regression: `pytest tests/test_gap_cad_t1_iges_pipe.py -k "not KVLCC2 and not ONRT"` → **10 passed**
(CLI, library entry, 4 DTMB5415 gold-regression tests, KCS, dependency-lock checks);
`pytest tests/test_nurbs_surface_degree_elevation.py tests/test_iges_142_parametric.py` → **24 passed**.

## 5. Why KVLCC2/ONRT still do not pass (open, measured)

The 3.4× is not enough, and the reason is a **reader defect, not speed**. Evidence
(`python tools/diag_iges_142.py work/assets/kvlcc2.igs 35`):

* 142 record **#35** reads `('0', '5', '33', '19', '1', '0', '0')`: surface `5` (type 128), then two
  composite curves `33`, `19` (type 102), then the preference flag `1`.
* `_curves_of_142` resolves **every integer token** as an entity pointer, so the flag `1` is taken as
  a reference to **entity 1** — a type-126 curve belonging to *another* face (`z = −28`,
  best distance to this surface 28.4 m). Face 0 therefore gets a phantom 2-point loop; measured
  loop count is 3 where the file's single 142 record supports 2.
* The genuine parameter-space image (`33`, UV box `(−173.1, 122.8) × (−19.4, 11.9)` inside the domain
  `(−176.1, 125.7) × (−19.7, 12.2)`) is cross-checked against the **outer** 3-D boundary `19`, which is a
  *different* loop (deviation 13.2 m ≫ limit 0.158 m). It is therefore demoted to a 3-D boundary,
  pushed through inversion, and reaches residual **136.6 m** ≫ tolerance 3.3e-4 m → the face is
  refused *after* minutes of fruitless multi-start search.

Recommended follow-up (not done here; it changes trim semantics and needs its own before/after):
(a) parse 142 by field position (pointer fields only) instead of scanning tokens, and (b) keep an
in-domain reference as an **exact parametric loop** with the agreement value reported, instead of
using disagreement with an unrelated reference as a veto.

## 6. Notes / traps

* The original report blamed "O(n²) boundary points × 49 seeds". The profile contradicts that: the
  search strategy costs 0.8 s of 95 s; the evaluator cost 85 s. Attaching a cause to a benchmark
  without a profile is how this survived two rounds.
* Memoising on a `frozen=True` dataclass must go through `object.__setattr__` in `__post_init__`;
  the cache is not a dataclass field, so equality/hash of surfaces is unchanged.
* The cost of a *failed* inversion is unchanged in structure: `invert_point` still runs the coarse
  and fine lattices when the walk misses. It is only cheaper per evaluation now. Removing the
  futile work needs the reader fix in §5.
* Reproduce: `python tools/bench_iges_inversion.py` (≈90 s, needs `work/assets/*.igs`).

## 7. Error judgement

| key quantity | reference | measured | error | < 3 % ? |
|---|---|---|---|---|
| KCS mesh triangles | HEAD mesh (4 499) | 4 499 | 0 | **yes** (identity, not a tolerance) |
| KCS worst inversion residual | HEAD (unchanged value) | identical | 0 | **yes** |
| KVLCC2 import wall time in test budget | ≤120 s | still >58 s (not finished) | — | **no** — honest FAIL, cause in §5 |
| DTMB5415 gold regression (volume/BB/gates) | gold STL | 4/4 tests pass | <3 % / <1 % as before | **yes** (unchanged, re-run) |

Acceptance criterion #2 (3 real IGES through all gates) therefore **remains open**; this round
removes a real defect on the path to it and documents the next one with measurements.