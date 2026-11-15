# Fleet intake: measured unit gate + trim-loop sampling (anchor GAP-CAD-T1-IGES-PIPE)

**Repo** TensorCAD · **Anchor** GAP-CAD-T1-IGES-PIPE · **Files**
`tensorcad/fleet_baseline.py`, `tests/test_fleet_trim_intake.py`,
`outputs/fleet_trim_intake/{evidence.json,figure.png}`

## 1. Case

The baseline fleet is the platform's geometry source: four publisher hull archives
(KCS, KVLCC2, JBC, ONRT) that must arrive as *watertight, volume-verified* solids.
Before this round three of the four returned `WITHHELD`. Two distinct defects sat
behind the KVLCC2 refusal, and both are intake defects rather than geometry problems:

| # | defect | evidence |
|---|---|---|
| A | the IGES **global section contradicts the coordinates**. `kvlcc2.igs` carries unit flag 2 (millimetres) over metre coordinates, so the declared metre factor is 1e-3 and a 320 m ship becomes a 0.334 m one — still watertight, still a mesh, of the wrong object. | declared reading: implied length 0.334 m vs published `L_PP` 320 m (ratio 0.001); metre reading: 333.619 m (ratio 1.043) |
| B | **every** basis surface is bound to a trimmed patch (19/19), so the untrimmed domains are not the ship. KVLCC2's untrimmed grid reaches z = −28.07 m for a 20.8 m draft and up to 6× outside the trimmed region — and the half-breadth field takes `max |y|` per `(x, z)` cell, i.e. the stray sheet *sets* the breadth instead of averaging out. | domain of patch 8 is 0…15.29 × 0…8.72 while its trim loop covers 0…15.29 × 0…8.72 only after the parameter-space pointer is chosen correctly (see §5) |

## 2. Method

| change | where | why it is measurable |
|---|---|---|
| `measured_unit_factor()` — both unit readings (declared global section, one-file-unit-is-metre) are scored against the publisher's own `L_PP` inside a wide decade band `UNIT_BAND = (0.5, 2.0)`; a plausible declared reading always wins, a contradicted one is overridden **and recorded**, and when neither is plausible the declared reading is kept with verdict `no_plausible_unit` so the existing scale guard refuses instead of guessing | `tensorcad/fleet_baseline.py` | gate output carries every candidate's implied length, ratio and verdict; the KCS/JBC/ONRT readings are unchanged by construction (ratio 1.06 / 1.07 / 0.007) |
| `sample_trimmed_patches()` + `trim_loops_uv()` — when the census says the patch set is fully trimmed, sample inside the trim loops (even-odd, so an inner loop is a hole) instead of over the untrimmed domain | same file | per-patch `kept`/`grid`/`domain_area_fraction` rows travel in the report; `patches_used` = 19/19 for KVLCC2 |
| frame guard corrected: the waterline-at-origin shift is `-z_min` (keel → z = 0). The previous `-z_max` moved the waterline plane by ~0.05 m and left the hull below the cut — the branch had never fired on a real file because defect A kept every point inside one z-cell | same file | `source.frame_shift_z` = 28.0 m for KVLCC2, 0.0 for KCS |
| one-step `dz` refinement on the trimmed route, reported as `dz_m` / `volume_m3` / `relative_change_percent` | same file | KVLCC2: +6.49 % at dz = 0.1 m vs −0.13 % at dz = 0.05 m; the change (−6.22 %) is in the evidence, not hidden |

Route selection is a data fact, not a heuristic: `basis_surfaces_untrimmed == 0` →
trim loops are the authoritative sheet; otherwise (KCS, 0/55 trimmed) the untrimmed
cloud is used exactly as before, and the trim route is only a fallback.

## 3. Results (read back from `outputs/fleet_trim_intake/evidence.json`)

| hull | route | unit verdict | volume [m³] | published [m³] | error | Cb | watertight | edge-manifold | triangles |
|---|---|---|---|---|---|---|---|---|---|
| KCS | untrimmed_patch_sheet (unchanged) | declared_consistent | 52 211 | 52 030 | **+0.35 %** | 0.6528 | ✓ | ✓ | 136 800 |
| KVLCC2 | trimmed_patch_loops (new) | declared_contradicted_by_geometry | 312 204 | 312 622 | **−0.13 %** | 0.8087 (pub. 0.8098) | ✓ | ✓ | 200 160 |

Before this round KVLCC2 was `WITHHELD: field produced fewer than two usable waterline
slices`; the mesh-vs-field cross-check (`max|y|` binning bias versus an independent
trapezoid integral of the same field) is reported per hull in the evidence file.

## 4. Error judgement

| quantity | reference | computed | error | < 3 % ? |
|---|---|---|---|---|
| KVLCC2 displacement | 312 622 m³ (NMRI T2015) | 312 204 m³ | −0.13 % | **yes** |
| KVLCC2 Cb | 0.8098 | 0.8087 | −0.14 % | **yes** |
| KCS displacement | 52 030 m³ | 52 211 m³ | +0.35 % | **yes** (regression check: identical route, unchanged) |
| KVLCC2 `max|y|` half-breadth | B/2 = 29.0 m | 29.06 m | +0.2 % | yes |

JBC (−18.28 %) and ONRT (−98.91 %) remain `WITHHELD` with their existing reasons; the
unit gate explicitly reports `no_plausible_unit` for ONRT as it arrives rather than
rescaling it, so those two are **not** claimed here. The analytic B-rep STEP/IGES
route via OpenCASCADE stays `NOT_SUPPORTED`: `pythonocc-core`/`OCP` is not installable
in this environment (no PyPI access, `cadquery` present without its `OCP` backend), and
no tessellation is guessed in its place.

## 5. Notes, and one trap worth recording

* **A degenerate trim curve scores a perfect inside-fraction.** The parameter-space
  pointer of a type-142 record is chosen by measurement ("a uv loop lies inside the
  knot domain"), but a curve that samples to one repeated point sitting in a domain
  corner satisfies that test with margin 1.0. `trim_loops_uv` therefore also requires
  the polygon's spread to be non-degenerate; `test_degenerate_candidate_cannot_win`
  pins this, and the 3-D candidate wins instead.
* Untrimmed domains are **not** a conservative fallback: binning `max |y|` makes
  out-of-loop sheet an over-estimate, which is why "it still produced a watertight
  mesh" was not evidence of correctness (+3.66 % for KVLCC2 on the untrimmed route).
* `dz` matters once a set is keel-to-deck deep (28 m) rather than draft-deep: the
  default 0.1 m cell is not converged, hence the reported refinement.
* KCS keeps its historical route deliberately: its 55 surfaces are all untrimmed, and
  changing a verified hull's path to "prove" a fix is how a regression gets shipped.

## 6. Reproduction

```bash
PYTHONPATH=$PWD python -m pytest tests/test_fleet_trim_intake.py -q     # 10 passed
PYTHONPATH=$PWD python -m pytest tests/test_iges_pipeline_gate.py tests/test_real_iges_gate.py \
    tests/test_iges_gold_regression.py -q                              # 21 passed, 1 skipped
```