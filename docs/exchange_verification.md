# Exchange seam verification: units, welding and the refusal guards

**Mission**: m16084 (verifier seat, selfstart) · **Repository**: TensorCAD · **Module**: `tensorcad/exchange.py` · **Evidence**: `artifacts/exchange_verification/figure.png` and `artifacts/exchange_verification/results.json` (every number below is injected from that file by `examples/exchange_verification/run_example.py`).

## 1. Case: what the seam is for, and what it must never do

`tensorcad.exchange` is the place a shipyard's file enters the quality audit: a faceted STEP (`FACETED_BREP` over `POLY_LOOP`), a neutral `.patchset.json`, or a refusal. The audit judges a triangle mesh, so two silent failures are worse than a refusal, because nothing downstream can notice them: **a wrong scale** and **a wrong triangulation**. The module's own docstring already promised both would not happen -- "the scale travels with the data instead of being assumed" and "it never guesses a tessellation" -- and neither promise had a test. This page is the verification that closed that debt, plus the four defects it found.

Fixture geometry is a unit cube (8 vertices, 12 triangles, 6 m² of surface, 1 m³) plus four hand-built STEP files that this repository's writer cannot produce: a per-corner-vertex cube, a face with a hole, a non-convex U-shaped loop, and a truncated cube. The hand-built files are emitted by a fixture writer inside the example, so the code under test is not grading its own output.

## 2. Method

1. **Independent structural parse.** The written file is re-read with a separate grammar (entity id, keyword, reference closure, framing markers) that shares no code with the reader; counts and reference closure are asserted against the mesh.
2. **Declaration round trip.** The writer is run at three scales (1, 1e3, 1e6) and the reader is given *no* scale, so the declared `LENGTH_UNIT` is the only route by which the coordinates can come back right.
3. **Naive-route measurement.** The pre-fix behaviour is reproduced on the current tree by a deliberately unguarded reader inside the example (hard-coded /1000, no framing check, concatenate-and-fan over every loop). Each fixture is read both ways; the error of the naive route is what the guards removed.
4. **Negative controls.** Analytic B-rep, holes, non-convex loops, truncation, dangling references, zero/negative/unrepresentable scales, empty mesh, IGES and unknown extensions all have to raise, and a convex quad has to keep working.

## 3. Results

### 3.1 Silent defects, measured on the naive route and closed by the guards

| defect | quantity | naive | truth | relative error | current outcome |
|---|---|---:|---:|---:|---|
| declared unit ignored on read (metre file divided by 1000) | max |x| of a 1 m cube, metres | 0.001 | 1 | -99.90 % | declared SI_UNIT($,.METRE.) read as metres; caller override still possible |
| per-corner VERTEX_POINTs not welded | closed cube reported watertight (1 = yes) | 0 | 1 | -100.00 % | welded: 8 vertices, watertight=True |
| face with a hole concatenated and fanned | planar area, m^2 | 6.25 | 3 | +108.33 % | refused |
| non-convex loop fanned from its first vertex | planar area, m^2 | 11 | 7 | +57.14 % | refused |
| truncated file accepted (no framing or reference closure check) | file accepted without its closing statements (1 = accepted, must be 0) | 1 | 0 | binary flag | refused |

### 3.2 Round trip at each declared scale (reader given no scale)

| writer `length_scale` | reader `length_scale` | source of that scale | vertex deviation (m) | triangles |
|---:|---:|---|---:|---:|
| 1 | 1 | declared | 0 | 12 |
| 1000 | 1000 | declared | 0 | 12 |
| 1000000 | 1000000 | declared | 0 | 12 |

### 3.3 Welding a per-corner writer's solid

| quantity | naive read | current read | truth |
|---|---:|---:|---:|
| vertices | 36 | 8 | 8 |
| removed duplicate rows | 0 | 28 | 28 |
| watertight | no | yes | yes |

Welding merges only *exactly* equal coordinates, so no position changes; the vertex numbering is remapped, which is why the test compares triangles as coordinate sets rather than index tuples. Passing `weld_duplicates=False` restores entity-level fidelity and the open-mesh reading, for a caller that wants to audit the file rather than the model.

### 3.4 Convex loops keep working

A 2 m x 1 m convex quad reads as 2 triangles with fan area 2 m² against a true area of 2 m² (+0.00e+00 %). The convexity guard is not an over-refusal: a fan over a convex loop is exact.

## 4. What was changed

| # | change | why it is not a workaround |
|---|---|---|
| 1 | writer emits the `SI_UNIT` prefix that matches `length_scale` (`$`/`.DECI.`/`.CENTI.`/`.MILLI.`/`.MICRO.`/`.NANO.`), and refuses a scale with no SI prefix | the file now describes its own numbers; a caller using 25.4 (inch) gets an error instead of a file that lies to every conformant reader |
| 2 | reader derives the divisor from the declared `LENGTH_UNIT`, records `length_scale_source`, and keeps an explicit `length_scale` override | the docstring's "scale travels with the data" claim is now true, and the assumed-millimetre fallback is labelled as an assumption |
| 3 | `weld_duplicates` (default on) merges exactly equal `VERTEX_POINT`s and reports `welded_vertex_rows` / `collapsed_triangles` | a legal per-corner file no longer reads as an open solid, so the audit judges the model instead of the writer's indexing |
| 4 | a face with more than one bound, or a non-convex `POLY_LOOP`, is refused with the face id | a fan over those invents area (see 3.1); refusing is the module's stated policy, and `skipped_faces`/`collapsed_triangles` keep the rest countable |
| 5 | truncation and dangling references are refused; `ADVANCED_FACE` dead code removed | a file missing its closing statements used to read as a smaller valid body |

## 5. Points to watch

1. **Still no third-party validation.** No CAD kernel is installed here, so the writer is validated by an internal round trip and an independent parse of the text, not by a foreign reader. `supported_routes()` keeps that at `WITHHELD`, correctly.
2. **`CONVERSION_BASED_UNIT` is not converted.** A file declaring inches has no SI prefix this reader honours; the millimetre fallback applies and `length_scale_source` says so. Convert or override the scale explicitly.
3. **Welding is exact equality only.** A writer that stores the same corner with different round-off (e.g. `0.1` vs `0.10000000000000001`) will not weld; that is a deliberate choice, since a tolerance would be a guess about the model.
4. **Convexity is tested in the loop's own plane** (Newell normal). A loop whose vertices are not coplanar has no single normal, and the guard then decides on the projected winding; such a face is not a planar facet and belongs to a B-rep, not here.
5. **The reader still fans convex n-gons.** A convex pentagon becomes three triangles; the area is exact, the diagonal choice is not checked for aspect ratio.

## 6. Error judgement

| quantity | benchmark | measured | verdict |
|---|---|---|---|
| round-trip vertex deviation | exact recovery (< 1e-15 m) | 0 m | **PASS (< 3 %)** |
| convex-quad fan area | exact (2.0 m²) | 2 m² | **PASS (< 3 %)** |
| max |x| of a 1 m cube, metres (declared unit ignored on read (metre file divided by 1000)) | naive route error | -99.90 % | **closed** (naive route now exact or refused) |
| closed cube reported watertight (1 = yes) (per-corner VERTEX_POINTs not welded) | naive route error | -100.00 % | **closed** (naive route now exact or refused) |
| planar area, m^2 (face with a hole concatenated and fanned) | naive route error | +108.33 % | **closed** (naive route now exact or refused) |
| planar area, m^2 (non-convex loop fanned from its first vertex) | naive route error | +57.14 % | **closed** (naive route now exact or refused) |
| file accepted without its closing statements (1 = accepted, must be 0) (truncated file accepted (no framing or reference closure check)) | must be refused | accepted=1 on the naive route, refused now | **PASS (refusal)** |

The largest naive-route error is 108.3 % on "face with a hole concatenated and fanned", far outside the 3 % acceptability line; every one of those routes is now either exact or an explicit refusal, which is the point. The claims that remain *unchecked* are the foreign-reader interoperability (no CAD kernel) and the `CONVERSION_BASED_UNIT` conversion, both marked above rather than claimed.

## 7. Reproducing

```bash
python3 examples/exchange_verification/run_example.py   # writes artifacts/ + this page
python3 -m pytest tests/test_exchange_verification.py -q  # 26 cases
```

Figure: left panel is the relative error the guards removed on three quantitative routes (log scale, the dashed line is the 3 % acceptability line); right panel is a closed cube written with one `VERTEX_POINT` per face corner -- 36 vertices and not watertight before welding, 8 vertices and closed after.
