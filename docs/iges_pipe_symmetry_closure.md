# Symmetry closure for half-hull imports (GAP-CAD-T1-IGES-PIPE, increment)

**Repo** TensorCAD · **stage** `tensorcad/geomquality/symmetry.py` · **date** 2026-11-16
**Evidence** every number below is copied from a tool's JSON in `artifacts/iges_seam/`
(`tools/iges_seam_sweep.py`, `tools/iges_gap_diag.py`, `tools/iges_boundary_census.py`),
zero hand-computed values.

---

## 1. Why this stage, and what it changes

The import pipeline's `G3_watertight` / `G6_single_shell` gates can never pass on a file
that is the **starboard half** of a hull: the surface legitimately ends on the centre
plane, and no local repair may close that boundary — capping it would build a wall through
the ship. Publisher ship IGES files are frequently exported exactly this way, so the
missing capability is not "better stitching" but **completing the half into the full hull**.

`symmetry_stage` does that: it detects a symmetry plane that carries the open boundary,
mirrors the mesh across it, reverses the winding of the mirrored copy (reflection reverses
handedness), snaps the seam vertices *onto* the plane and welds them so the seam is shared
rather than duplicated. Detection is fail-closed on two counts:

* **one-sidedness** — no vertex may lie more than `band` beyond the candidate plane, so a
  solid straddling its own bounding-box face is never doubled;
* **dominance** — the winning plane must carry ≥ `min_coverage` (default 0.40) of the
  open-boundary *length* and ≥ 3× the runner-up plane's share.

Measurement is the default (`symmetry_closure=False` reports the detected plane and
changes nothing), because mirroring changes the model's extent.

## 2. Analytic verification (independent of the stage's own bookkeeping)

`tests/test_iges_symmetry_closure.py`, 7 tests, green:

| test | assertion | measured |
|---|---|---|
| half cylinder → full | watertight, 1 component, 1 closed component, oriented | `free_edge_count = 0` |
| half cylinder volume | vs closed form `πr²h = 2π` | within **2e-3** (96 segments; the residual is the chordal area error `1 − sin(2π/96)/(2π/96) ≈ 7e-4`, not the mirror's) |
| half sphere volume | vs `4π/3` | within **3e-3** |
| exact doubling | mirrored volume = 2× the half's signed volume | `rel = 1e-9` |
| closed mesh | a mesh with no open boundary is never mirrored | refused (`detected = None`) |
| straddling mesh | a solid that crosses the plane | `y` plane refused |
| measurement-only | mesh untouched | vertices/faces identical |

The doubling test is the important one: it fails if the stage is doing anything other than
an exact reflection (a re-loft, a tolerance fudge, an extra cap would all break it).

## 3. Measured behaviour on a real publisher file (`work/assets/kcs.igs`, 55 trimmed faces)

### 3.1 Seam budget ladder (no symmetry, no capping)

| seam tolerance | gates passed | free edges | components |
|---|---|---|---|
| none (measure) | 3/7 | 803 | 20 |
| 0.02 | 3/7 | 374 | 14 |
| 0.2 | 3/7 | 300 | 58 |

Closing the seam buys free edges but not a solid, and a wide budget costs manifoldness
(58 components at 0.2 against 14 at 0.02): the seam stage is not the missing piece.

### 3.2 What the surviving open boundary actually is (seam = 0.02)

`tools/iges_gap_diag.py` classifies all 374 surviving free edges by the repair that would
close each one:

| class | count |
|---|---|
| closable by a weld at 1e-6 (coincident twin) | **0** |
| closable by a T-junction split | 6 |
| genuine gaps (`> 1e-6` from the nearest free edge) | **368** |

Gap quantiles (distance from an edge midpoint to the nearest other free edge):
`q00 0.25`, `q10 0.52`, `q50 2.20`, `q90 6.95`, `q100 11.41` model units.

So the residual boundary is **not** a tolerance problem: the halves of the model are metres
apart. Nothing in the repair chain may legitimately bridge that.

### 3.3 Where those edges lie, and why the symmetry stage refused

`tools/iges_boundary_census.py` (band = 0.5 % of extent = 1.22):

| region | edges | share |
|---|---|---|
| centre plane `y = 0` | 219 | 58.6 % |
| `elsewhere` (mid-surface) | 117 | 31.3 % |
| deck / transom / keel | 38 | 10.2 % |

That reads like a half hull — which is why the symmetry closure was implemented. But the
*length* carried by a tight band does not support it. Coverage of the open-boundary length
by each bounding-box plane (`artifacts`, diagnostic run at seam = 0.02):

| plane | band 0.244 (1e-3 extent) | band 1.0 | band 2.5 |
|---|---|---|---|
| `y = 0` (centre plane) | **0.348** | 0.467 | 0.473 |
| `z = lo` (keel) | **0.356** | 0.364 | 0.379 |
| `y = hi`, `z = hi`, `x = lo`, `x = hi` | ≤ 0.30 | ≤ 0.32 | ≤ 0.40 |

At the bandwidth that a real cut would satisfy (0.244 on a 244-unit model) the centre plane
carries 34.8 % and the keel plane carries 35.6 %: **the two are indistinguishable**, and the
coarse-band prominence of `y = 0` is an artefact of the metres-wide defect gaps sitting near
the centre plane, not of a clean cut. The stage therefore refuses on this file
(`"no plane carries a dominant share of the open boundary; mesh kept"`), which is the
designed fail-closed outcome: mirroring on this evidence could have fabricated a hull.

**Honest verdict on the anchor's `3 个真实 IGES 全过门`:** not met, and this increment does
not meet it. The measurement says why: `kcs.igs` is an open patch set whose free boundary is
(a) 31 % mid-surface defect gaps of metres, and (b) cut-plane boundary that is not
geometrically planar to better than ~1 % of the model extent. Gate closure needs either a
conformant export or the OCP healing route — capping or mirroring this boundary would
produce a solid whose volume is an artefact.

## 4. What ships

| file | role |
|---|---|
| `tensorcad/geomquality/symmetry.py` | the stage: detect / mirror+reverse / snap+weld, fail-closed |
| `tensorcad/geomquality/trim_mesh.py` | `symmetry_closure`, `symmetry_band`, `symmetry_min_coverage` on both routes; `report["symmetry"]` |
| `tensorcad/iges_pipe_cli.py` | `--seam-tolerance`, `--symmetry-closure`, `--symmetry-min-coverage` (the CLI previously could not reach the seam stage at all) |
| `tools/iges_gap_diag.py` | free-edge classification: weld / T-junction / genuine gap |
| `tools/iges_seam_sweep.py` | seam-budget ladder, one JSON per rung |
| `tools/iges_boundary_census.py` | where the open boundary lies, with the band it used |
| `tests/test_iges_symmetry_closure.py` | 7 analytic tests (§2) |
| `artifacts/iges_seam/*.json` | the measurements tabulated above |

## 5. Limitations / next

* The detector assumes the cut boundary is planar to within `band`. Files whose centre-plane
  boundary is ragged after tessellation (KCS) are refused — by design, but it means the
  stage only helps clean half-hull exports.
* Mirroring assumes the half is the *starboard* side and the plane is a bounding-box face;
  an arbitrarily placed plane is not supported.
* `G7_no_degenerate` (3 degenerate/duplicate triangles) and `G4_manifold` (27) are
  independent of this stage and unchanged.
* Next step with the highest leverage: make the tessellated cut boundary conform to the
  plane (snap the boundary polyline before mirroring), then this stage would apply to KCS.