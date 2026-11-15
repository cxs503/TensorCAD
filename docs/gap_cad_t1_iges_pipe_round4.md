# GAP-CAD-T1-IGES-PIPE — round 4: the wetted-area overshoot was duplicate skin, and it is now removed

**Repo** TensorCAD · **Anchor** GAP-CAD-T1-IGES-PIPE · **Branch** `agent/m22083`
**Maturity: L3 for the wetted-area quantity** (implementation + library entry + regression
test + external reference within 3 %, on two independent knobs). The *watertight* verdict of
the anchor stays unreachable for publisher files and is not claimed (round 3 §5: the models
are open half-breadth surfaces by construction).

## 1. Case

Round 3 left the anchor with three quantities inside the ±3 % band against the published KCS
particulars (∇ −0.31 %, Cb −0.44 %, Cm +0.5 %) and one honest FAIL: the wetted surface came out
at 11 249 m³·m⁻²… 11 249 m², **+19.4 %** over the workshop value 9 424 m². It pinned that
overshoot in a test and blamed "overlapping patches count the same skin twice" — an
*explanation*, never a *measurement*. This round measures the duplication, removes it, and
re-derives the number.

Geometry: `KCS_hull_SVA.igs` (SVA Potsdam), 55 untrimmed B-spline patches → 2 510 vertices,
4 324 triangles at default chord tolerance, half-breadth open model (`classify_model` →
`symmetry_factor = 2`), design draft T = 10.8 m.

## 2. Method

1. **Census of the soup** on the imported triangles: 4 324 triangles, only 5 exactly identical
   vertex triples — so the duplication is *geometric*, not literal. Each triangle was given a
   ray from its centroid along its outward normal (outward = away from the bounding-box
   centroid) and tested for an occluder; hidden area was 1 630 m², i.e. **14.5 % of the soup**,
   and 96 % of it lies within 0.5 m of the sheet that hides it. That is the quantitative
   statement round 3 was missing: the patches do not merely touch, they lie on top of each other.
2. **Why the first two estimators failed.** Summing triangle areas over the soup is +19.4 %
   (double skin). The *envelope graph* estimator — area of `b(x,z)` via `√(1+b_x²+b_z²)` — is
   −21.8 %: a flat of bottom and a transom are edge-on in the `(x, z)` projection and vanish
   from a graph representation. An occlusion filter that casts only along `+y` gave −21.5 % for
   the same reason (it discarded 1 556 m² of *legitimate* bottom skin). Both are recorded here
   because each of them looks plausible until it is measured — the estimator must be
   orientation-agnostic.
3. **The stage that shipped** (`tensorcad/geomquality/hull_surface.py`, stdlib only — the
   package charter keeps `dependencies = []`): `wetted_area(vertices, faces, waterline=…,
   dedup=True)` casts the centroid ray along the triangle's own outward normal through a
   uniform 3-D cell index (reach = 4 cells, stride = ½ cell), with two additions that the naive
   version needs to be correct:
   * **coincident sheets** never cross the ray, so a triangle whose centroid lies *on* another
     triangle (plane distance ≤ tol and inside it) is treated as hidden too;
   * **mutual annihilation** — two exactly coincident sheets each lie on the other, so with the
     rule above *both* would vanish. A coincident claim is therefore only honoured from a
     **lower-numbered** triangle, so each coincident cluster keeps exactly one copy. Without
     this tie-break an earlier prototype lost both copies (found by the duplication-invariance
     test, not by the hull file).
   `hull_characteristics(..., dedup=True)` uses it and reports `s_wetted_naive` (raw sum) and
   `s_wetted_discarded` next to `s_wetted`, so the raw defect stays on record.

## 3. Results (KCS, same import as round 3: `import_cad`, weld 1e-3, 61 stations × 48 levels)

| quantity | reference (published) | measured before | measured now | error now | < 3 % ? |
|---|---|---|---|---|---|
| ∇ | 52 030 m³ | 51 870.7 | 51 870.7 | −0.31 % | **yes** |
| Cb (Lpp basis) | 0.651 | 0.6482 | 0.6482 | −0.44 % | **yes** |
| **wetted S** | **9 424 m²** | 11 249.3 (**+19.4 %**, FAIL) | **9 619.2** | **+2.07 %** | **yes** |
| Cm | 0.977 | 0.9825 | 0.9825 | +0.5 % | yes |
| watertight | topology gate | false (open model) | false | — | not achievable, declared |

`validate_against_particulars` conjunction: **passed** (was `False` in round 3).
De-duplication moved only the wetted area: ∇/Cb/Cm are section integrals of the envelope
breadth and were never sensitive to the soup.

## 4. Figure

`artifacts/gap_t1_iges_dedup/kcs_dedup_validation.png` (rendered by `artifacts/make_figure.py`
from the same JSON the tests assert on): left — per-quantity error bars after de-duplication
with the ±3 % band, wetted area now green at +2.07 % (the naive +19.4 % is in the title);
right — the two knobs. Index-cell size is insensitive (span/16 → span/128 gives 1.89 % → 2.14 %),
which is the behaviour an index must show; the coincidence tolerance has a **plateau from
1e-5·span to 1e-4·span (1.20 %…2.87 %) and collapses beyond 3e-4·span** (−8.4 %, then −43 %),
because a tolerance comparable to real geometric separation starts hiding genuine skin.

## 5. Notes / traps

* **The +2.07 % is not tuned per hull, but it is tolerance-sensitive.** The default is
  1e-4·span = 24 mm on a 244 m model; at 1e-5·span it reads +2.87 %, at 3e-4·span −8.4 %. The
  agreement inside the plateau is a property of the geometry, but any single number should be
  quoted with the plateau, not alone.
* Triangle-level decisions mean a triangle straddling the edge of an overlap is kept or dropped
  whole; on this file that is worth a few m², and it is why the result is *not* claimed to
  converge below ~1 %.
* The outward-normal heuristic (away from the bounding-box centroid) is exact for a star-shaped
  hull skin and **is not a solid modeller**: a deeply concave pocket can be hidden in the wrong
  direction. Declared in the docstring of `wetted_area`.
* The old test `test_kcs_workshop_geometry_matches_published_particulars` pinned the FAIL
  (`10 % < error < 30 %`, `passed is False`). It was flipped to assert the de-duplicated value
  **and** that the naive sum is still > 1.15 × reference, so the defect remains measured.

## 6. Tests

| test | invariant |
|---|---|
| `test_box_wetted_area_is_exact_and_unchanged_by_dedup` | box hull: exact area to 1e-12, nothing discarded |
| `test_duplicated_skin_does_not_double_the_area` | **duplication invariance**: an exact copy of the side adds 0 m², and the raw sum shows the +2·L·T defect |
| `test_sheet_hidden_under_the_outer_skin_is_dropped` | a patch 1 cm inside the skin is dropped entirely, area unchanged |
| `test_full_hull_model_is_not_doubled` | full hull: symmetry factor 1, no doubling |
| `test_kcs_dedup_flips_the_pinned_wetted_area_defect` | KCS wetted area within 3 %, ∇ unchanged, raw sum still > 1.15× |
| `test_kcs_workshop_geometry_matches_published_particulars` (flipped) | conjunction `passed is True`, raw defect still on record |

`python -m pytest tests/test_hull_surface_dedup.py tests/test_hull_surface_validation.py -q`
→ **11 passed** (17.9 s).

## 7. Not done / declared

* **KVLCC2** (second workshop hull, S 27 194 m², ∇ 312 622 m³) with the same default knobs:
  its import exceeds the local 110 s cap, so it was submitted to the cluster queue as **cj3898**
  and was still queued when this round closed. The second-hull cross-check is therefore
  **not run**, not passed.
* **Watertight import** for publisher files: unchanged, unreachable (round 3 §5), not claimed.