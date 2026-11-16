# GAP-CAD-T1-IGES-PIPE — round 5: type-142 read by field position (the ``IT`` flag is not a pointer)

**Repo** TensorCAD · **Anchor** GAP-CAD-T1-IGES-PIPE · **Branch** `agent/m22574`
**Maturity: L3 for the trim-reading quantity.** The defect round 4 identified and left open
(«(a) parse 142 by field position … (b) keep an in-domain reference as an exact parametric
loop») is fixed, with a before/after diff on three archives and a synthetic-layout test.
The anchor's *watertight* verdict for publisher files stays unreachable and is **not**
claimed; acceptance criterion #2 is still open (§6).

## 1. Case

The trimming reader (`tensorcad/geomquality/iges_trim.py`) turns the trim records of a
publisher IGES file into boundary loops. It handled entity 142 (curve on a parametric
surface) by scanning **every integer token** of the record and resolving each one that
happened to name an existing curve. IGES 5.3 §4.25 defines the record as

```
142, CRV, SPTR, BPTR, CPTR, IT
```

with `IT` a **0/1 flag** ("the composite `CPTR` is (0) / the individual curves `CRV` and
`BPTR` are (1) the preferred representation"). A flag of `1` is therefore a valid pointer
value, and on the SIMMAN KVLCC2 archive it resolved to directory entry 1 — a type-126
curve belonging to a *different face*.

Measured field map (both writers, no guessing):

| archive | record | raw parameter data |
|---|---|---|
| KVLCC2 | DE 35 | `142,0,5,33,19,1,0,0` → CRV=0, SPTR=5, BPTR=33, CPTR=19, **IT=1** |
| KP505 | DE 13 | `142,0,1,0,11,2,0,0` → CRV=0, SPTR=1, BPTR=0, CPTR=11, **IT=2** |

Both pad the record with trailing zeros that are not pointers either.

## 2. Method

1. **Root cause, part (a).** `_fields_of_142()` reads the record positionally. The
   alignment (whether the writer repeated the leading `142`, which `parse_iges` strips)
   is decided by `SPTR`, the only field whose role is unmistakable: it must name a
   surface (108/120/122/128). Fields after `IT` are outside the entity definition and
   are never resolved.
2. **Root cause, part (b).** `_curves_of_142()` now honours the spec meaning of `BPTR`:
   when it exists and its samples fall inside the surface's knot domain it **is** the
   parameter-space image and is consumed as an exact parametric loop (no inversion). The
   record's other representation (`CRV`, or the 3-D children of the `CPTR` composite) is
   a **cross-check**: the measured deviation is reported on the loop, and a disagreement
   is recorded as a note instead of vetoing the only usable trim. Records without a
   usable `BPTR` (KP505 writes `BPTR = 0`) keep the previous measurement-based
   classification unchanged.
3. **Measurement, not assertion.** `tools/measure_142_route.py` writes a JSON
   fingerprint (faces, loops per face, trim routes, cross-checks, tessellation refusals)
   for every archive present on the machine. The two fingerprints below were produced by
   `git show HEAD:tensorcad/geomquality/iges_trim.py` (before) and the patched file
   (after), same command, same assets.

## 3. Results

| archive | faces | loops per face | trim route | read wall time |
|---|---|---|---|---|
| KVLCC2 — before | 19 | `1,3,3,3,3,3,1,1,3,1,3,3,3,3,3,3,3,3,3` | 3-D (inverted) | 2.3 s |
| KVLCC2 — after | 19 | `1 × 19` | parametric | 3.0 s |
| KP505 — before | 8 | `1,1,1,1,1,2,2,2` | inverted | 0.37 s |
| KP505 — after | 8 | `1 × 8` | parametric (faces 4-7) | 0.53 s |
| KCS — before / after | 55 | `1 × 55` (unchanged) | inverted | 0.2 s |

Phantom loop on KVLCC2: **15 of 19** faces carried a third loop that resolved to entity 1
(distance 28.4 m from the trimmed surface). With it gone, the veto against the unrelated
3-D composite also disappears — measured disagreements on the same records: DE 35
`1.326e+01`, DE 69 `4.775e-04`, DE 107 `1.408e+01`, DE 141 `1.309e+01` model units. DE 69
is the case where the two representations *are* the same curve; it keeps the old refusal
text and the measured agreement, so the information survives.

Consequence at the mesh stage (KP505 propeller, the fixture that was refused in the
original defect report):

| | before | after |
|---|---|---|
| tessellated triangles | 6 519 | **69 510** |
| refused faces | 1, 5, 6, 7 | **1** |

Full pipeline on KVLCC2 (`import_cad`, weld 1e-3, seam 1e-3): **811.7 s**, 19 faces,
4 813 vertices, 42 579 triangles, gates `G1_all_faces_read` **True** and
`G2_trim_on_surface` **True** (both were unreachable while the faces were refused),
`G5_positive_volume` True, `G7_no_degenerate` True; `G3_watertight`, `G4_manifold`,
`G6_single_shell` False with 249 free edges and 19 components before repair — the model
is an open half-breadth surface by construction, as recorded in round 3.

## 4. Notes / traps

* A positional rule alone is not enough: the phantom loop was removed *and* the veto had
  to stop overriding the parameter-space curve. Fixing only (a) still left each KVLCC2
  face with two 3-D loops and refused them.
* The read is not the cost. 2.9 s of the 811.7 s import is the reader; the rest is
  inversion/tessellation of 19 trimmed patches at the default density. The wall-clock
  budget problem of the acceptance criterion is a meshing-cost problem, not a reader one.
* `BPTR` and the `CPTR` composite genuinely differ on most KVLCC2 patches (13 m). Any
  future "repair" that averages or stitches the two representations would be inventing
  geometry.
* `_fields_of_142` decides the alignment from `SPTR` rather than from token 0, because
  `parse_iges` silently strips a repeated entity-type token; a record that keeps it is one
  field longer and the naive positional read would shift every field by one.

## 5. Error judgement

| key quantity | reference | before | after | < 3 % ? |
|---|---|---|---|---|
| loops per trimmed KVLCC2 face | 1 (one trim record = one loop) | 3 on 15/19 faces | **1 on 19/19** | **yes** (identity, not a tolerance) |
| phantom reference to entity 1 | 0 | 15 records | **0** | **yes** (identity) |
| KP505 refused faces | 0 | 4 | **1** | **no** — improved, not closed (§6) |
| KCS faces / loops (regression) | 55 / 1 each | 55 / 1 each | 55 / 1 each | **yes** (unchanged) |
| KVLCC2 ∇ via `fleet_baseline` route | 312 622 m³ (SIMMAN 2020) | 312 568.0 m³ (−0.017 %) | unchanged (that reader is not touched) | **yes** |

The trim-route volume is **not** used for a 3 % judgement: the imported mesh is open
(249 free edges, 19 components), so its enclosed volume is not a displacement.

## 6. Honest status of acceptance criterion #2

「3 个真实 IGES 全过门」 remains **open**:

* KVLCC2: G1/G2 now pass, G3/G4/G6 cannot pass for a half-breadth publisher surface
  (declared since round 3), and the mesh stage needs 811.7 s — far outside a test budget.
* ONRT: the same reader fix applies to its 142 records, but its mesh stage is larger
  still; not measured in this round (**not_run**, not passed).
* KP505 fixture: 1 of 8 faces is still refused (a genuine off-surface 3-D loop; its
  refusal message is the one the round-2 test pins).

## 7. Files

`tensorcad/geomquality/iges_trim.py` (`_fields_of_142`, `_curves_of_142`) ·
`tests/test_iges_142_fields.py` (5 tests) · `tools/measure_142_route.py` ·
`tools/patch_142_fields.py` · `runs/round5_142_fields/{before,after}.json`

Verification run in this session:

```
python -m pytest tests/test_iges_142_fields.py -q                     -> 5 passed
python -m pytest tests/test_iges_142_parametric.py tests/test_iges_trim_binding.py \
    tests/test_iges_trim_import.py tests/test_trim.py tests/test_trim_guards.py \
    tests/test_gap_cad_t1_iges_pipe.py -q -k "not KVLCC2 and not ONRT" \
  -> 77 passed, 3 skipped, 2 deselected, 1 xfailed, 1 failed
```

The one failure (`test_trim.py::test_the_module_needs_no_third_party_stack`) is
**pre-existing** in this tree: it fails because `numpy` is importable in the environment,
with or without this change (verified on the unmodified file before patching).