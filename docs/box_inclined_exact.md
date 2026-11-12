# An exact attitude reference for a barge: combined heel **and** trim

`tensorcad.box_inclined` (`tensorcad/box_inclined.py`, `tests/test_box_inclined.py`,
`examples/box_inclined/run_example.py` → `examples/box_inclined/outputs/`).

## 1. The gap this closes

Every hydrostatic reference in this package is one-dimensional in attitude:

| route | attitude it can hold | why it cannot hold the other |
|---|---|---|
| `tensorcad.gz.box_righting_arms_analytic` | heel, inside the wall-sided window | single formula, refuses past the bilge/deck edge |
| `tensorcad.box_gz` (`box_gz`, `box_gz_area`) | heel, any angle | a **2-D** section is clipped and integrated; a trim would tilt the section plane and the section is not the section |
| `tensorcad.trim` | upright trim | TPC/MCT from waterplane moments of an upright hull |
| `tensorcad.hydrostatics.compute_hydrostatics` | mesh cut by a **horizontal** plane | a trimmed floating body has no horizontal waterline in body coordinates |

A ship in a turn, and any damage case, floats at a heel *and* a trim at the same time. There
was no exact answer for that attitude, so the mesh route (rotate the mesh, cut it) could only
be checked against other mesh routes — the objection `docs/box_large_heel.md` makes: two
integrators on the same mesh agreeing is not evidence about the physics.

## 2. Method (no discretisation, no case analysis)

The submerged body of a box is `box ∩ {n · p ≤ c}`: a convex polyhedron. Clipping the six
box faces by the waterline plane with Sutherland–Hodgman and closing the cut with the
waterplane polygon gives that polyhedron **exactly**.

* **Volume and buoyancy centre** — divergence theorem over the triangulated clipped faces
  (signed tetrahedra about the origin), so `∇`, `LCB`, `TCB`, `VCB` are exact to round-off.
* **Waterplane** — capsule polygon measured *in the plane* (Green's theorem), giving area,
  centroid and the in-plane second moments; `BM_T = i_transverse/∇`, `BM_L = i_longitudinal/∇`.
* **Righting arm** — `GZ = t · (B − G)` where `t` is the horizontal transverse direction of
  *that* waterplane (body `x` projected into the plane, crossed with the plane normal). Both
  attitudes enter one expression; there is no "heel formula" and "trim formula" to disagree.
* **Attitude** — up-normal `n = normalize(−sin τ, −sin φ cos τ, cos φ cos τ)`; `φ > 0`
  starboard down, `τ > 0` bow down; submerged side `n · p ≤ c`. The displacement is held and
  the waterline offset is solved by bisection on the strictly monotone `volume(c)`.

Because the cut is a general plane, **bilge emergence, deck-edge immersion and a cut running
through a single corner are all the same code path**. That is what makes the reference usable
outside the window where the closed form refuses to answer.

## 3. Measured results

Barge `L × B × D = 10 × 2 × 1.2 m`, displacement `∇ = 16 m³` (draft `T = 0.8 m`), `KG = 0.55 m`.
All numbers below are produced by `examples/box_inclined/run_example.py`; none is typed in.

| check | reference | measured worst relative error | < 3 %? |
|---|---|---|---|
| upright `VCB` | `T/2` | `3.6e-15` | ✅ |
| upright `i_transverse`, `i_longitudinal` | `L B³/12`, `B L³/12` | `0.0` (exact) | ✅ |
| upright `BM_T` | `B²/(12T)` | `3.7e-15` | ✅ |
| trim `LCB`, trims 1–4° (cut stays on the side walls) | `L² tan τ / (12 T)` | `2.3e-15` | ✅ |
| displacement over 12 combined attitudes (heel 0–30°, trim 0–4°) | `16 m³` | `5.7e-15` | ✅ |
| `GZ`, heel 1–15° | wall-sided closed form | `1e-12` relative (test bound) | ✅ |

Two things the numbers also say, and they matter more than the round-off:

* **The closed form is not the same answer past its window.** At 40° the exact arm is
  `0.170156 m` while the wall-sided formula (evaluated outside its validity) is more than
  `0.02 m` away — asserted in `TestHeel::test_the_wall_sided_formula_is_not_vacuous_here`, so
  the agreement inside the window is not a tautology.
* **The trim closed form has a window too.** `L² tan τ/(12T)` assumes the cut stays on the
  side walls; at `τ = 6°` the bow waterline reaches `z = 1.225 m > D = 1.2 m`, the deck is
  immersed, and the exact `LCB` departs from the closed form. The example tabulates both
  regimes and flags `deck_immersed`, rather than quietly extrapolating.
* `GZ` peaks at `0.175598 m` at `45°` on this barge, and the arm is still `0.091252 m` at
  80° — a range of stability that the one-line formula cannot describe at all.

Figures: `examples/box_inclined/outputs/fig1_gz_curves.png` (exact vs wall-sided `GZ`, and the
log-scale relative deviation that shows where the formula stops applying),
`fig2_trim_and_combined.png` (`LCB` vs trim against the closed form; `TCB` drift with heel at
4° trim).

## 4. How it was verified

* **29 tests** in `tests/test_box_inclined.py`, all green (`1.2 s`, no NumPy).
* **Independent routes, and what they are worth.** The tests use two references, neither of
  which is the module's own algorithm:
  1. the wall-sided closed form (transcribed from the literature, upright window only);
  2. a **2-D section oracle written in the test file**: the rectangle is clipped by the
     body-frame waterline and its area/centroid come from a bisection on the offset. This is
     one dimension below the module and shares no code, but it *does* share the half-space
     idea, so what it demonstrates is "two implementations of the same cut agree at 2 × 10⁻⁹ m
     on `TCB`/`VCB`", not "two independent methods agree". Called that way in the test
     docstring rather than dressed up as more.
  3. the analytic trim moment `L² tan τ/(12T)` is a genuinely independent third anchor for
     the trim-only branch (derived from the odd/even split of `∫ x (a + x tanτ) dx`).
* **Structural checks**: the clipped polyhedron closes (`polyhedron_volume_m3 ==
  volume_m3`), every waterplane vertex stays inside the box, mirrored heels mirror `TCB`,
  `dGZ/dφ|₀ = GM` to `1e-6`, the curve rises then falls with an interior peak, and impossible
  requests (volume above the box, negative volume, zero length, NaN attitude, waterline below
  the keel) are refused with `BoxInclinedError`.

## 5. Honest limits

1. **A box.** A shaped section is out of scope by construction; this is the attitude anchor
   the mesh routes get compared *to*, not a replacement for them.
2. **No mesh-route comparison was run in this mission.** The natural next step is
   `compute_hydrostatics` on a box mesh rotated to `(φ, τ)` — the exact route is now available
   for it. Stating it as not-run rather than implying agreement.
3. **No free-surface correction, no sinkage/trim solve, no windage, no criterion.** The
   attitude and the displacement are inputs; `GZ` is a lever, not a stability verdict (that is
   `tensorcad.statutory` / `tensorcad.intact_stability`).
4. **Degenerate attitudes are refused, not approximated**: a waterline that leaves no
   submerged volume raises, and a displacement above the box volume raises.
5. **Round-off floor.** Volume closure is `5.7e-15` relative on this geometry; below that the
   route has nothing more to say, and a hull-scale comparison should not be pushed past ~1e-13.