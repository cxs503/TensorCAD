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
