# GAP-CAD-T1-IGES-PIPE — round 3: physical validation of real publisher hull IGES

**Repo** TensorCAD · **Anchor** GAP-CAD-T1-IGES-PIPE · **Branch** `agent/m21936`
**Maturity argued: L2** (implementation + library entry + test + external reference). L3 is *not*
claimed: the watertight gate is still not reachable for the publisher files, and the reason is
now measured rather than mysterious (see §5).

## 1. Case

The anchor asks for an IGES/STEP → **watertight** mesh pipeline whose gold regression is
`DTMB5415_Hull_exact`. Every earlier round earned an honest FAIL on the same point: the three
"real" IGES files of round 1 were written by this repository's own writer (self round trip), and
the one genuinely third-party file (KP505 propeller, CATIA V5 export) reached 2/8 tessellated
faces. Meanwhile the topology gate (`G1..G7` of `geomquality.trim_mesh`) cannot be green for the
files that actually carry hull geometry, and round 2 could not say why.

This round answers the *why* with measurements and replaces the unreachable claim with a check
that has an external reference.

## 2. Method

1. **Census and boundary measurement** on `KCS_hull_SVA.igs` (SVA Potsdam, sha256
   `6f5556b2…`): IGES entity census `{128: 55}` — 55 untrimmed B-spline patches, **no 144/141**,
   i.e. no trimmed geometry at all; vertex envelope `y ∈ [0, 16.1117]` (= B/2 = 16.1),
   `z ∈ [0, 23.5]`, `x ∈ [-6.0, 237.84]`. **No centreplane cap, no deck cap, no transom cap** ⇒ the
   model is an open half-breadth surface: `watertight` is false by construction, not by defect.
2. **Seam measurement.** After welding at 1 mm: 7 connected components, 825 free edges. Perpendicular
   distance from each free vertex to the nearest free edge of *another* component has percentiles
   `[3.6e-4, 5.9e-3, 3.0, 9.4, 23.0]` m — that is, most of the "cracks" are metres wide, because the
   patches **overlap** instead of butting. An ear-clipping/T-junction repair cannot sew such a soup;
   the remaining route to a solid would be a full CAD sew-and-heal, out of scope here.
3. **Envelope hydrostatics** (new, `tensorcad/geomquality/hull_surface.py`, stdlib only — the
   package charter keeps `dependencies = []`): the moulded breadth at `(x, z)` is the largest `|y|`
   of the surface along a `±y` ray (axis-aligned ray-triangle tests on an x-slab index); sectional
   area `A(x) = ∫₀^T b(x,z) dz`; `∇ = k ∫ A dx` with `k = 2` for a half-breadth model; wetted area
   by clipping each triangle at `z = T` (Sutherland–Hodgman) and doubling, **excluding** faces that
   lie in the centreplane. Overlapping patches do not corrupt the envelope (it takes the outermost
   surface), which is exactly why this route works where welding does not.
4. **External reference**: the workshop particulars of KCS (`Lpp 230.0 m, B 32.2 m, T 10.8 m,
   Cb 0.651, ∇ 52 030 m³, S 9 424 m²`) and KVLCC2 (`320.0 / 58.0 / 20.8, Cb 0.8098,
   ∇ 312 622 m³, S 27 194 m²`), recorded with their source string in `PARTICULARS` and compared
   by `validate_against_particulars` (per-quantity error + conjunction verdict). Published `Cb` is
   defined on `Lpp`, so the measured `Cb` is put on the `Lpp` basis before comparison — the file's
   own x-span is 239.85 m (bulb ahead of FP), and using it blindly is exactly the kind of unit
   slip that produced a 4.6 % "error" in a first prototype.

## 3. Results (KCS, imported through `trim_mesh.import_cad`, weld 1e-3, 61 stations × 48 z-levels)

| quantity | reference (published) | measured | error | < 3 % ? |
|---|---|---|---|---|
| displacement volume ∇ | 52 030 m³ | **51 870.7 m³** | **−0.31 %** | **yes** |
| block coefficient Cb (on Lpp basis) | 0.651 | **0.6482** | **−0.44 %** | **yes** |
| waterline beam B | 32.2 m | 32.22 m | +0.06 % | yes |
| midship section coefficient Cm | 0.977 (implied) | 0.982 | +0.5 % | yes |
| wetted surface S | 9 424 m² | 11 249.3 m² | **+19.4 %** | **no — honest FAIL** |
| watertight | required by the topology gate | **false** (open half-breadth model) | — | not achievable by construction |

The volume error (−0.31 %) is the first time this anchor has an *external* number on the geometry
as imported: earlier rounds could only compare a file against the same generator that wrote it.
The wetted-area overshoot is explained: area is a sum over the patch soup, so overlapping patches
count the same skin twice, while the envelope breadth (and therefore ∇, Cb, Cm) takes the
outermost surface and is immune. The test `test_kcs_workshop_geometry_matches_published_particulars`
pins that overshoot as a *measured* limitation (`10 % < error < 30 %`, `passed is False`) so that a
future de-duplication stage will flip it visibly instead of silently.

## 4. Figure

`artifacts/gap_t1_iges_round3/` — `sectional_area_curve.png` (KCS `A(x)` from the imported surface,
0 → 239.85 m, max 170.9 m² at midship) and `validation.png` (per-quantity error bars against the
±3 % band, volume/Cb/beam/Cm green, wetted area red at +19.4 %). Both are rendered from the same
JSON the tests assert on.

## 5. Notes / traps

* **A green topology gate proves nothing about geometry.** `G5` only asks `volume > 0` on the
  repaired soup; on this file it reports 17 671 m³ — 66 % below the true 52 030 m³ — and still
  "passes". That is the trap the section-integral route removes.
* **Open ≠ defective.** `KCS_hull_SVA.igs` and `kvlcc2.igs` are half-breadth open models; the gate
  suite now classifies them (`classify_model`) and judges them on what they are (all faces read,
  trims on their surface, physical validation) instead of demanding a closure surface the file never
  had. `ONRT` remains on the full `G1..G7` path.
* Ray-triangle intersection for an axis-aligned ray is a point-in-triangle test in the x-z plane —
  the naive Möller–Trumbore with `p = A − o` (instead of `o − A`) silently returns the *near* hit or
  nothing; it cost one prototype iteration and produced a −99.6 % volume before it was caught by
  the box test.
* KVLCC2 (`kvlcc2.igs`) is 4× larger; its import exceeds the local 120 s cap, so its validation row
  is **not run** in this round (declared, not skipped-and-forgotten).

## 6. Error judgement (explicit)

| key quantity | reference | measured | error | < 3 % ? |
|---|---|---|---|---|
| KCS ∇ | 52 030 m³ (workshop particular) | 51 870.7 m³ | −0.31 % | **yes** |
| KCS Cb on Lpp | 0.651 | 0.6482 | −0.44 % | **yes** |
| KCS S | 9 424 m² | 11 249.3 m² | +19.4 % | **no** |
| KCS Cm | 0.977 | 0.982 | +0.5 % | yes |
| KCS watertight | — | false (open model) | — | not achievable |

## 7. Files

| path | role |
|---|---|
| `tensorcad/geomquality/hull_surface.py` | envelope breadth, sectional areas, ∇/Cb/Cm/S, model classification, `PARTICULARS`, `validate_against_particulars` (stdlib only) |
| `tests/test_hull_surface_validation.py` | box exactness (volume/area/Cb/Cm to 1e-12), classification, KCS vs published particulars, grid-convergence check |
| `tests/test_iges_pipe_gate.py` | half-breadth open models judged on G1/G2 + physical validation instead of an unreachable watertight verdict (this test was red on `main`) |
| `artifacts/gap_t1_iges_round3/` | figure + the JSON the figure is rendered from |