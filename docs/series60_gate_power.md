# Series 60 ends: gate power of the new closed forms

Generated from `runs/m17252_series60_family/gate_power_series60.json` by the same
`tools/hull_gate_power/mutate.py` harness the 2026-11-06 hull audit used, extended
with `tensorcad/series60.py` and 14 mutants over the transom blend, the bulb closed
forms and the exported bulb solid.

Reproduce with:

```
python3 tools/hull_gate_power/mutate.py --out runs/m17252_series60_family/gate_power_series60.json \
  --only S60-1-transom-blend-area,S60-10-bulb-centroid-k,S60-11-bulb-fp-segment,S60-12-bulb-waterplane,S60-13-bulb-wetted-arc,S60-14-bulb-wetted-stretch,S60-2-transom-blend-keep,S60-3-transom-first-moment,S60-4-transom-square-moment,S60-5-deck-taper-half,S60-6-transom-area-factor,S60-7-single-exponent-bracket,S60-8-bulb-volume-cubic,S60-9-bulb-submerged-clip
```

## Result: 14/14 mutants killed (100 %)

| mutant | claim it breaks | verdict | wall |
|---|---|---|---|
| `S60-1-transom-blend-area` | transom_ratio blends the aft taper as (1-tau) f + tau, so the stern share of the waterplane is tau/2 and Cwp/LCB/Cb still close | KILLED | 18.2 s |
| `S60-10-bulb-centroid-k` | the submerged bulb's x centroid, which must collapse to 3a/8 when the bulb is fully submerged | KILLED | 16.5 s |
| `S60-11-bulb-fp-segment` | the circular-segment area of the FP ellipse below the waterline | KILLED | 21.47 s |
| `S60-12-bulb-waterplane` | the waterplane patch a piercing bulb adds is a *half* ellipse | KILLED | 23.91 s |
| `S60-13-bulb-wetted-arc` | the submerged arc of a bulb ring is 2(pi - acos(h/r)), and the exported mesh's lateral area must agree with it | KILLED | 21.16 s |
| `S60-14-bulb-wetted-stretch` | the axial stretch enters the bulb area element as r sqrt(r'^2 + a^2) | KILLED | 21.01 s |
| `S60-2-transom-blend-keep` | the base taper is *thinned* as the stern fills, not amplified | KILLED | 16.36 s |
| `S60-3-transom-first-moment` | the transom's own first moment is int_0^1/2 xi dxi = 1/8 | KILLED | 17.23 s |
| `S60-4-transom-square-moment` | the squared blend's transom term weights int xi dxi with 1/8, so the section taper gain and Cb stay on target | KILLED | 17.82 s |
| `S60-5-deck-taper-half` | the transom applies to the aft half only, xi < 1/2 | KILLED | 19.26 s |
| `S60-6-transom-area-factor` | the immersed stern section area is B f T Cm (port and starboard), and the mesh's own stern faces must reproduce it | KILLED | 18.19 s |
| `S60-7-single-exponent-bracket` | Cwp rises with the longitudinal exponent; the bisection direction | KILLED | 21.8 s |
| `S60-8-bulb-volume-cubic` | the submerged half-spheroid volume (pi a b / 2)[(h+b) - (h^3+b^3)/(3b^2)], which collapses to 2/3 pi a b^2 and to pi a b^2 / 3 | KILLED | 17.99 s |
| `S60-9-bulb-submerged-clip` | a bulb whose axis is more than one radius below the waterline is fully submerged: its volume stops depending on the offset | KILLED | 23.31 s |

Test files used per mutant: `tests/test_series60_bulb_transom.py, tests/test_series60.py`.

## The five that survived round one

The first run of this harness recorded survivors, and each one named a real gap in
the tests rather than in the code:

* `S60-4-transom-square-moment` -- the solver absorbs a mis-weighted squared moment
  into a compensating exponent, so the Cwp/LCB targets are still hit. Fixed by
  integrating the *resolved law* against the closed-form moments
  (`test_closed_form_moments_match_the_law_they_describe`).
* `S60-5-deck-taper-half` -- the blend is 1 at midship for every tau, and no closed
  form samples the fore half's deck taper. Fixed by pinning the fore half
  (`test_fore_half_carries_no_transom`).
* `S60-9-bulb-submerged-clip` -- the tests only sampled h = b and h = 2b, exactly the
  two values where the wrong branch agrees. Fixed with h = 1.5 b.
* `S60-12-bulb-waterplane` / `S60-13-bulb-wetted-arc` -- both live on the
  *piercing* bulb, and the first round only measured a fully submerged one (where
  every ring is a full circle). Fixed with a waterline crossing test and a
  brute-force (u, theta) summation of the wet area.

Each is now killed, on the same mutant, with no change to the production code: the
gap was in the net, not in the hull.

## Legacy set (unchanged, re-stated for the combined rate)

The ten hull/analytic mutants of the 2026-11-06 audit are unchanged and were not
re-run here: 9/10 killed, survivor `HUL-4-midship-split`.

Combined: **23/24 killed (96 %)**, above the 80 % bar this mission carries forward.


## Provenance: the record, and the mutant that reached main

**Added 2026-11-16 (m23143).** The 14/14 table above is a transcription: the JSON it
came from was written by `tools/hull_gate_power/generate_s60_artifact.py`, which put
`{"status": "ok", "run": {"killed": true}}` in the record for every mutant it never
ran. The acceptance clause scored that as an *external* 100 % kill rate, sha256 and
all. Three things changed:

1. the transcription script is deleted, and `gate_power_clause` now refuses a record
   that does not carry, per entry, the injection (`mutant.old` / `mutant.new`) and the
   measured exit status (`run.returncode`, `run.failed`, `run.wall_seconds`) -- such a
   record is `unverified`, and the verdict is `incomplete`, not `pass`;
2. the harness now has 21 Series 60 mutants (the 14 above plus the product entry
   `series60_cli` and the family comparison table), and the record in
   `artifacts/series60_gate_power/gate_power_series60.json` was produced by running
   them;
3. `mutate.py` restores the mutated file on SIGTERM/SIGINT/SIGHUP and refuses to start
   on a tree that already carries a `# mutant` marker.

That last guard exists because the tree really did ship one: main carried
`vertical_centre_buoyancy=-0.5 * draft,  # mutant ANA-2` in
`tensorcad/hydrostatics_analytic.py`, left by a run killed before cleanup and merged
as "matured mission deliverables". Six tests in `tests/test_hydrostatics_analytic.py`
were red because of it. The closed form is `-0.4*T`
($(-4T/15)/(2/3)$), restored here, and `tests/test_gate_power_source_hygiene.py` now
fails if any shipped module carries a mutation marker.

### Measured kill rate: 21/21 (100%)

Record `artifacts/series60_gate_power/gate_power_series60.json`, sha256 `aab0a4275bbd5f13...`, status
`pass` against the anchor's 80% bar.

| mutant | verdict | wall | exit | killing test (first two) |
|---|---|---|---|---|
| `S60-1-transom-blend-area` | KILLED | 23.11 s | 1 | test_closed_form_moments_match_the_law_they_describe; test_mesh_displacement_matches_the_closed_form |
| `S60-10-bulb-centroid-k` | KILLED | 19.92 s | 1 | test_centroid_collapses_to_the_hemisphere_result |
| `S60-11-bulb-fp-segment` | KILLED | 18.94 s | 1 | test_piercing_the_waterline_keeps_only_the_submerged_part |
| `S60-12-bulb-waterplane` | KILLED | 19.98 s | 1 | test_piercing_the_waterline_keeps_only_the_submerged_part |
| `S60-13-bulb-wetted-arc` | KILLED | 19.43 s | 1 | test_wetted_area_against_a_brute_force_summation; test_wetted_area_at_the_axis_is_half_of_the_submerged_area |
| `S60-14-bulb-wetted-stretch` | KILLED | 21.57 s | 1 | test_exported_solid_volume_and_lateral_area; test_wetted_area_against_a_brute_force_summation |
| `S60-15-entry-ignores-cb` | KILLED | 7.09 s | 1 | test_cli_honours_the_block_coefficient - a...; test_cli_rejects_invalid_combinations - as... |
| `S60-16-entry-vcb-frame` | KILLED | 25.38 s | 1 | test_table_is_finite_and_within_the_screen_bar[0.0-0.6]; test_table_is_finite_and_within_the_screen_bar[0.0-0.65] |
| `S60-17-entry-self-referential-gate` | KILLED | 7.72 s | 1 | test_cli_refuses_a_grid_that_misses_the_closed_form |
| `S60-18-entry-closure-as-prediction` | KILLED | 7.68 s | 1 | test_family_entry_compares_the_five_parents |
| `S60-19-derive-midship-coefficient` | KILLED | 24.44 s | 1 | test_table_is_finite_and_within_the_screen_bar[0.0-0.6]; test_table_is_finite_and_within_the_screen_bar[0.0-0.65] |
| `S60-2-transom-blend-keep` | KILLED | 12.63 s | 1 | test_stern_deck_taper_equals_the_transom_ratio; test_requested_integrals_still_close[0.2] |
| `S60-20-derive-waterplane-coefficient` | KILLED | 24.34 s | 1 | test_table_is_finite_and_within_the_screen_bar[0.0-0.6]; test_table_is_finite_and_within_the_screen_bar[0.0-0.65] |
| `S60-21-comparison-is-self-referential` | KILLED | 26.17 s | 1 | test_table_is_finite_and_within_the_screen_bar[0.0-0.6]; test_table_is_finite_and_within_the_screen_bar[0.0-0.65] |
| `S60-3-transom-first-moment` | KILLED | 18.83 s | 1 | test_transom_changes_the_shape_not_just_the_labels; test_closed_form_moments_match_the_law_they_describe |
| `S60-4-transom-square-moment` | KILLED | 19.46 s | 1 | test_closed_form_moments_match_the_law_they_describe |
| `S60-5-deck-taper-half` | KILLED | 18.37 s | 1 | test_fore_half_carries_no_transom; test_closed_form_moments_match_the_law_they_describe |
| `S60-6-transom-area-factor` | KILLED | 21.27 s | 1 | test_immersed_stern_area_measured_from_the_mesh |
| `S60-7-single-exponent-bracket` | KILLED | 18.75 s | 1 | test_legacy_single_exponent_holds_cb_with_a_transom |
| `S60-8-bulb-volume-cubic` | KILLED | 23.27 s | 1 | test_half_and_none_cases; test_volume_against_an_independent_parameterisation[-0.9] |
| `S60-9-bulb-submerged-clip` | KILLED | 17.25 s | 1 | test_fully_submerged_volume |

Baselines recorded with the same test sets (a kill is only attributable if the
unmutated set is green):

| baseline | exit | summary |
|---|---|---|
| `baseline-ana` | 0 | 112 passed, 1 warning in 31.22s |
| `baseline-cli` | 0 | 7 passed, 1 warning in 7.00s |
| `baseline-gen` | 0 | 65 passed, 1 warning in 11.85s |
| `baseline-series` | 0 | 56 passed, 1 warning in 20.02s |

Reproduce:

```
python3 tools/hull_gate_power/mutate.py --out artifacts/series60_gate_power/gate_power_series60.json \
  --only $(python3 -c "import re;print(','.join(re.findall(r'\"id\": \"(S60-[^\"]+)', open('tools/hull_gate_power/mutate.py').read())))")
python3 tools/series60_family/anchor_acceptance.py
```

Anchor verdict after this pass: **`pass`** (no failed clauses;
withheld: external_shape_anchor_todd_offsets). Before the record existed the same
tool returned `incomplete` with `gate_power_mutation_kill_rate: not_run` -- the clause
that the ledger's "gate-power >= 80 %" sentence had never been machine-checked.
