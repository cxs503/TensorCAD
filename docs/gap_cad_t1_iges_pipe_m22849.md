# GAP-CAD-T1-IGES-PIPE — the seam stage was hiding its own damage (m22849)

**Anchor** GAP-CAD-T1-IGES-PIPE · **Mission** m22849 · **Repo** TensorCAD · **Branch** agent/m22849
**Scope of this round** the repair chain (`seam_stage`) of the trimmed import pipeline: what it
actually does to a publisher patch set, and a guard so a stitch can never be counted as a heal
while it quietly makes the mesh non-manifold.

## 1. Why this round exists

Rounds m17333/m17426/m20797/m21127/m22022/m22381/m22464/m22645/m22790 each closed a real defect
(142 parametric route, singular edges, interior sampling, latin-1 decode, gold reference
quadrature). The anchor is still `open`, and the nearest unmet acceptance clause is "3 个真实 IGES
全过门". Re-running the three publisher files this round gives the current baseline:

| File | bytes | sha256 (head) | route | faces | triangles | free edges (before → after repair) | gates |
|---|---|---|---|---|---|---|---|
| KCS_hull_SVA.igs | 246,164 | `6f5556b2…` | trimmed | 55/55 | 4,314 | 1,697 → 795 | G1 G2 G5 pass; **G3 G4 G6 G7 fail** |
| KP505_modelScale_250mm.igs | 1,042,550 | `0bb54ae7…` | trimmed | 8/8 | 77,429 | 678 → 80 | G1 G2 G5 pass; **G3 G4 G6 G7 fail** |
| DTC_hull.igs | 2,309,858 | `53c3be90…` | trimmed | — | — | — | `not_run`: >120 s per command budget |

(commands: `python -m tensorcad.geomquality.import_cli <file> --out-dir /tmp/{kcs,kp505,dtc}_now`)

So both files now *read* completely (G1/G2 pass — the 142/128 work of the earlier rounds holds),
and both fail exactly on closure. The accepted next step is the seam stage.

## 2. What the seam stage really does (measured, not asserted)

`seam_stage(vertices, faces, tolerance=t)` on the KCS tessellation, audit before/after
(`tests/test_seam_stage_guard.py`, `_final.py` in the mission scratch):

| tolerance (m) | duplicate triangles dropped | free edges before → after | non-manifold edges | orientation conflicts | components |
|---|---|---|---|---|---|
| 0.005 | 5 | 1,697 → 1,456 | 0 → **6** | 0 → **131** | 55 → 25 |
| 0.020 | 22 | 1,697 → 1,335 | 0 → **14** | 0 → **216** | 55 → 28 |
| 0.050 | 42 | 1,697 → 1,275 | 0 → **37** | 0 → **252** | 55 → 35 |

The stitch closes 15–25 % of the open boundary and, in the same pass, introduces 6–37
non-manifold chords, 131–252 orientation conflicts and 5–42 duplicate triangles. A downstream
`audit_mesh` then reported `manifold: False` on a mesh whose only visible evidence used to be
"free edges went down". That is the defect this round fixes: **the repair chain was scoring itself
on the one number it improves.**

Two further measurements show the same thing from the other side. Signed volume of an open shell
is not a measurement: the KCS half hull integrates to 17,671 m³ before stitching and to 34,768 m³
at tolerance 0.05 — a 2× jump produced entirely by closure bookkeeping. Mirroring the half hull
about y = 0 and welding gives 69,537 m³, i.e. +33.5 % against the published displacement
L·B·T·Cb = 230 × 32.2 × 10.8 × 0.651 = 52,083 m³. `tests/test_seam_stage_guard.py` therefore
never asserts a volume for these files.

## 3. The change

`tensorcad/geomquality/trim_mesh.py`

* `_drop_duplicate_triangles(faces)`: a T-junction split re-emits the parent and its children; when
  the same chord is hit twice the split repeats an existing triangle. Repeats are never geometry,
  so they are dropped and counted (`stitch["duplicate_triangles_dropped"]`).
* `seam_stage(..., strict=False)` now records `topology_before` / `topology_after`
  (free edges, non-manifold, orientation conflicts, components) and `manifold_degraded`.
  `strict=True` refuses a stitch that would increase the non-manifold edge count and returns the
  input mesh with `rejected` / `note` set.
  The default keeps the stitch (it still closes most of the boundary, and the existing
  `test_real_kcs_seam_stage_reduces_free_edges` pins that capability) — but the damage is now
  in the report instead of only in a later audit.

No geometry heuristic, tolerance default or gate threshold was changed; behaviour of a
non-degrading stitch is byte-identical.

## 4. Tests

| Test | What it pins | Result |
|---|---|---|
| `test_duplicate_triangles_are_dropped_under_any_rotation` | dedup is rotation-invariant | PASS |
| `test_t_junction_is_stitched_and_reported` | a genuine T-junction is still stitched; record carries before/after topology | PASS |
| `test_measured_stitch_never_returns_duplicate_triangles` | audit of the returned mesh: 0 duplicates | PASS |
| `test_publisher_patch_set_stitch_is_refused_when_it_degrades_manifoldness` | on KCS: default keeps the stitch **and** reports `manifold_degraded`; `strict=True` refuses it | PASS |

`python -m pytest tests/test_seam_stage_guard.py tests/test_iges_seam_stitch.py
tests/test_iges_trim_import.py -q` → **24 passed, 4 skipped, 1 xfailed**.
`tests/test_iges_gold_regression.py` (gold reference = `DTMB5415_Hull_exact.stl`) → **9 passed**.

## 5. Error judgement (vs the 3 % standard, item by item)

| Quantity | Base | Measured | Error | <3 %? |
|---|---|---|---|---|
| Gold regression volume (4-patch DTMB5415 loft, patch-set route) | gold STL | watertight, 1 component, bbox within 1 % | test threshold | **yes** (unchanged by this round) |
| KCS half-hull signed volume | 26,041 m³ (published Cb, half) | 17,671 m³ (open shell) | −32.1 % | **no — and not a valid measurement: the shell has 795 free edges, so the divergence integral has no cap** |
| KP505 closed volume | — | 8.3e-05 m³ shell | — | **not applicable**, 80 free edges remain at the blade root |
| Seam stitch topology | manifold input | +6…+37 non-manifold edges | — | **no** (this is the defect fixed) |

## 6. Limitations and what is left

* Closure of the publisher patch sets is **not** achieved. KCS needs its 55 patches sewn at the
  surface level (a shared discretisation of the common edge curves), not snapped after
  tessellation; KP505 needs the blade-root/hub junction closed (its 80 free edges are confined to
  x ∈ [−27.4, 21.7] mm, y ∈ [15.2, 23.4] mm, z ∈ [−17.2, 12.3] mm and none lies on a symmetry
  plane).
* `DTC_hull.igs` is not measured: >120 s in one command. It needs a queued run, not the control
  plane.
* The greedy post-tessellation stitch is the wrong tool for publisher files and this round only
  makes its cost visible. The next unit of work is patch-level sewing (build the patch adjacency
  graph from boundary polylines, resample each shared curve once, tessellate to the shared
  discretisation), which removes the T-junctions by construction.