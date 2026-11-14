# GAP-CAD-T1-BASELINE-FLEET — type 126 is a curve, and the reader was indexing it wrong

Mission **m20952** (TensorCAD) · anchor `GAP-CAD-T1-BASELINE-FLEET` · date 2026-11-15
Code: `tensorcad/iges.py` (`BSplineCurve.from_parameters`, `CurveOnSurface.from_parameters`,
`IgesFile.surface_census`, `IgesFile.trim_pointer_audit`)
Tool: `tools/baseline_fleet/trim_audit.py` · Data: `artifacts/baseline_fleet/t126_trim_audit.json`
Tests: `tests/test_iges_surface_census.py`, `tests/test_iges_trim_binding.py` (30 passed)

---

## 1. The claim this round had to test

Mission m20862 recorded the blocker for three of the six hulls as an entity-type gap:

> "`BSplineSurface` decodes **one** entity type — 128, the *rational* B-spline surface.
> IGES also defines **126**, the *non-rational* B-spline surface* … the two hulls carrying
> the geometry predominantly in type 126 are exactly the two refused on geometry."

and it added `ENTITY_BSPLINE_SURFACE_POLYNOMIAL = 126` plus a census key
`unreadable_type126`, which `fleet_baseline.build_ship` then used to word its refusals.

**That claim is false, on the entity identity and on the reader.** IGES 5.3 defines no
non-rational surface entity: type 126 is the *rational B-spline curve* (§4.24), type 128
the rational B-spline surface (§4.26). Trim boundaries are curves, so a hull archive is
expected to carry hundreds of 126 records next to its surfaces. Measured census on the
publisher archives (`surfaces_type128` / `curves_type126` / basis surfaces bound to a
trimmed patch):

| hull | file | 128 surfaces | 126 curves | trimmed patches (143/144) | basis surfaces trimmed | untrimmed |
|---|---|---|---|---|---|---|
| KCS | `KCS_hull_SVA.igs` | 55 | 0 | 0 | 0 | 55 |
| JBC | `jbc.igs` | 45 | 172 | 19 | 0 | 45 |
| KVLCC2 | `kvlcc2.igs` | 19 | 304 | 19 | 19 | 0 |
| ONRT | `onr5613-…mergedR.igs` | 132 | 820 | 80 | 0 | 132 |

KCS — the one hull that verifies — carries **zero** 126 records, which is the only fact the
old story fit; JBC carries 172 of them and was refused, but so does ONRT with 820 while
KVLCC2 carries 304 and *does* bind every one of its 19 basis surfaces to a trimmed patch.
The correlation m20862 reported was between 126 count and refusal, and it does not hold
once the entities are counted for what they are.

## 2. The defect the falsification exposed

While checking the curves, the decoder itself failed the check that matters: **a knot
vector is non-decreasing, and coordinates live in the model box.** On JBC DE 8
(K = 18, M = 2, PROP3 = 0) the reader produced

| reading | knot vector (first entries) | first control point |
|---|---|---|
| reader, before | `106.995, −3365.936, −0.0, 107.52, …` | `0.462, 0.462, 0.538` |
| IGES 5.3 field order | `0.3887, 0.3887, 0.3887, 0.4615, …` → `1.0` | `−3326.8, 0.0, 84.5` |

The decoder read the parameter list **tail-first** (`control … knots`) and skipped three
fields into the knot vector, so millimetres landed in the knot slots and the knot vector
landed in the control points. Type 126 was never "undecoded"; it was decoded into garbage,
which is worse, because the garbage has the shape of geometry. The fix reads the spec
order — knots `K+M+2`, then weights `K+1` (rational only), then `K+1` xyz triples — and
after it **all 1 296 type-126 curves in the three archives decode with non-decreasing knot
vectors and 0 records failing** (57 ONRT curves sample outside the model box, which is the
expected signature of parameter-space loops).

## 3. The trim pointer, decided by measurement rather than by field order

Type 142 (`curve on a parametric surface`) carries `CRV, SPTR, BPTR, CPTR, IT`. The reader
assigned `curve_3d`←BPTR and `curve_uv`←CPTR (one field to the right), and a test pinned
that assignment (`curve.curve_uv == 19`) without any evidence that 19 was a parameter-space
loop — a test that asserted the implementation, not the contract.

`IgesFile.trim_pointer_audit` decides it from the data: a parameter-space loop must lie
inside the basis surface's knot domain, a model-space curve cannot. Over the 19 KVLCC2
patches (only archive with 142 records):

| candidate | patches whose sampled points lie inside `u × v` on ≥ half the samples |
|---|---|
| `BPTR` (field 4, index 3) | **18 / 19** (1 undecided within a 10-point margin) |
| `CRV` (field 2, index 1) | 0 / 19 |

So the loop is BPTR, and `CurveOnSurface` now maps the record to the spec
(`curve_3d`←CRV, `curve_uv`←BPTR, `preferred`←CPTR, `preferred_flag`←IT).

## 4. What is still open (unchanged acceptance, no verdict claimed)

* **JBC and ONRT**: 0 of 45 / 0 of 132 basis surfaces bind to a trimmed patch. Their
  wrappers are type 143 with a type-141 boundary, and `resolve()` finds no type-128 basis
  surface behind those pointers under either directory-numbering convention. That is a
  *pointer-binding* gap, measured, and it is the next concrete thing to fix — not an
  entity type that does not exist.
* **KVLCC2**: 19/19 patches are trimmed and their loops are now read in (u, v), but the
  untrimmed integral of the same 19 basis surfaces is **114 817 m³ against the publisher's
  312 622 m³ (−63 %)**, i.e. the patch set as a whole does not close the hull; trimming
  alone will not rescue it. Also measured: this CATIA export declares `2H MM` (unit flag 2)
  while its coordinates are metres (x-span 328.1, half-breadth 29.1, depth 29.2 against
  published L_PP 320 / B 58 / D 30) — a 1000× contradiction between the file's own unit
  declaration and its geometry.
* No hull's gate verdict changes in this round. The anchor stays **open**.

## 5. Repro

```bash
python tools/baseline_fleet/trim_audit.py jbc.igs kvlcc2.igs KCS_hull_SVA.igs onrt.igs \
       --out artifacts/baseline_fleet/t126_trim_audit.json
python -m pytest tests/test_iges_surface_census.py tests/test_iges_trim_binding.py -q
```

Publisher bytes are not redistributed here (neither publisher states a licence); the
digests are pinned in `benchmarks/verified/baseline_fleet/fleet_specs.json`.
