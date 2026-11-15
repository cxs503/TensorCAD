# GAP-CAD-T1-IGES-PIPE — a both-sides hull was measured as half a ship

**Repo** TensorCAD · **Anchor** GAP-CAD-T1-IGES-PIPE · mission **m22344**
**Verdict: defect fixed, closed-form regression added, external beam check inside 3 %; the
displacement comparison on the gold file is grid-limited and is reported as such, not as a pass.**

## 1. Case

`tensorcad.geomquality.hull_surface` is the physical-validation stage of the anchor: it turns
an imported (open, publisher) hull surface into ∇, Cb/Cm, Lwl, beam and wetted area, and
`validate_against_particulars` judges those against published particulars at 3 %. Round 4
(`docs/gap_cad_t1_iges_pipe_round4.md`) validated the *wetted-area* half of it on KCS.

Its section breadth was `envelope_breadth = max |y|` at every station, multiplied afterwards by
`classify_model(...)["symmetry_factor"]`. For the **half-breadth** publisher files that is
right (one side, doubled). For a model that carries **both sides** the factor is 1 and
`max |y|` is exactly **half** the breadth, so:

* A(x) came out half the true section area,
* therefore ∇ half, and beam half,
* while `Cb = ∇ / (Lwl · beam · T)` came out **correct**, because numerator and denominator
  shrink by the same two.

That is why three rounds of validation missed it: the anchor's gold reference
`DTMB5415_Hull_exact.stl` is a both-sides model, every workshop hull used for the external
checks (KCS, KVLCC2) is half-breadth, and the only full-model unit test asserted the *wetted
area* (`test_full_hull_model_is_not_doubled`), which never goes through `envelope_breadth`.

## 2. Method

1. **Closed-form reproducer** (no asset, no mesh, 0.8 s): a box L×B×T = 10×4×2 m carried with
   both sides (y ∈ [−2, 2]). Exact ∇ = 80 m³, exact beam 4 m, Cb = Cm = 1 by construction.
2. **Fix**: `envelope_breadth(..., full=)` and `sectional_area(..., full=)` gain a two-sided
   reading (`y_max − y_min`, initialised to ±∞, `0.0` when the ray misses everything);
   `hull_characteristics` infers the flag from `classify_model` and exposes
   `full_breadth: bool | None` so the pre-fix behaviour is **one keyword away** (a measurement,
   not a claim). The half-breadth path is bit-for-bit unchanged.
3. **Gold file**: `DTMB5415_Hull_exact.stl` (164 926 triangles, both sides, `symmetry_factor`
   = 1) measured with the shipped code path at `stations=8, levels=8`, before and after; the
   independent reference is the **shoelace area of the section polylines** produced by
   `benchmarks/verified/iges_import/gold_loft.section_polyline` (the anchor's own gold tooling),
   closed along the waterline and integrated over the same stations — it shares no code with
   `envelope_breadth`.
4. `scripts/dtmb5415_gold_full_model_check.py` writes the same table at `stations=21,
   levels=24`; runs `runs/dtmb5415_full_model/evidence_local.json` (8×8) is the committed one.

## 3. Results (programmatic, `runs/dtmb5415_full_model/evidence_local.json`)

Closed-form box (L×B×T = 10×4×2):

| quantity | exact | before (`full_breadth=False`) | after | after error |
|---|---|---|---|---|
| ∇ [m³] | 80.0 | **40.0** | **80.0** | 0.0 % (rel 1e-12) |
| beam [m] | 4.0 | **2.0** | **4.0** | 0.0 % |
| Cb [–] | 1.0 | 1.0 | 1.0 | 0.0 % (the defect is invisible here) |

`DTMB5415_Hull_exact.stl`, waterline 6.15 m, 8 stations × 8 levels:

| quantity | before | after | independent reference | after vs reference | published |
|---|---|---|---|---|---|
| ∇ [m³] | 4 084.75 | **8 169.51** | 9 013.72 (shoelace sections) | **−9.37 %** | — |
| beam [m] | 9.538 | **19.076** | 19.076 (= y-span of the mesh) | 0.0 % | **19.06** → **+0.09 %** |
| Cb [–] | 0.5239 | 0.5239 | — | — | — |
| ratio after/before | — | **∇ 2.0000 / beam 2.0000** | — | — | — |

## 4. Notes / traps

* **The 2.0000 ratio is the whole finding.** It is measured on the real gold file, not on a
  synthetic box, and it lands on exactly the factor the algebra predicts (1.99999998 for ∇,
  2.0 for beam), while Cb does not move at all.
* **The residual −9.37 % is grid, not shape.** At 8 stations × 8 levels the z-trapezoid and the
  x-midpoint rule are both coarse, and the *reference* closes every section along the waterline,
  which adds the waterline segment the hull surface itself does not have below a deck. The
  cluster run at 21 × 24 is **cj4057** and was still queued when this mission closed: the
  converged displacement of the gold file is **not established here** and is not claimed.
  The beam is grid-insensitive (it is an extremum) and that is why the one external number this
  round does claim is the beam: 19.076 m vs 19.06 m published = **+0.09 % < 3 %**.
* **My own patch had a defect for ten minutes**: the two-sided span was initialised as
  `best_lo = 0.0`, so a model living entirely at `y > 0` reported `y_max` — the offset-box test
  below caught it, and it is why the fix initialises to ±∞. Recorded because the test, not the
  reasoning, is what found it.
* **`classify_model` cannot resolve an offset both-sides model** (`y_min ≥ 0` ⇒ "half-breadth").
  That is a convention, not a measurement error: no ship model lies entirely at `y > 0`. The
  escape hatch is `full_breadth=True`, which a caller with an unusual frame must combine with
  the factor themselves — documented in the test.
* Values measured before this round on the half-breadth workshop files (KCS +2.07 %, KVLCC2
  −0.017 % ∇) are unaffected: on those files the flag is a no-op, proved by
  `test_half_breadth_path_is_unchanged`.

## 5. Files

`tensorcad/geomquality/hull_surface.py` (fix) ·
`tests/test_hull_characteristics_full_model.py` (7 tests, 0.8 s) ·
`scripts/dtmb5415_gold_full_model_check.py` · `runs/dtmb5415_full_model/{evidence_local.json,figure.png}` ·
`scripts/kvlcc2_wetted_area_check.py` (round-4 §7 debt, submitted as cluster **cj4055**, still
queued at close → **not run, not passed**).

## 6. Tests

`python -m pytest tests/test_hull_characteristics_full_model.py tests/test_hull_surface_validation.py -q`
→ **12 passed, 13.62 s**.

| test | invariant |
|---|---|
| `test_envelope_breadth_returns_the_two_sided_span_for_full_models` | one side = B/2, `full=True` = B |
| `test_full_box_displacement_beam_and_coefficients_are_exact` | ∇ = L·B·T, B, Lwl, Cb = Cm = 1 to 1e-12 |
| `test_forcing_half_breadth_readings_reproduces_the_halving` | ∇, B halved; Cb, Cm unchanged |
| `test_half_breadth_path_is_unchanged` | sym 2, ∇ exact, flag is a no-op |
| `test_a_model_living_entirely_on_positive_y_is_declared_half_breadth` | convention + escape hatch |
| `test_validate_against_particulars_passes_a_full_model_now` | volume error 0.0 %, `passed is True` |

Not run locally (heavy, submitted to the queue): `tests/test_hull_surface_dedup.py` and
`tests/test_iges_gold_regression.py` (the KCS import alone exceeds this shell's 120 s cap).

## 7. Error judgement

| key quantity | baseline | measured | error | < 3 % ? |
|---|---|---|---|---|
| beam of the gold hull | 19.06 m (published DTMB5415) | 19.076 m | **+0.09 %** | **yes** |
| ∇ of the gold hull | shoelace section integral 9 013.72 m³ | 8 169.51 m³ | −9.37 % | **no** — grid-limited (8×8), converged value not established |
| ∇ of a closed-form box | 80.0 m³ | 80.0 m³ | 0.0 % | **yes** |
