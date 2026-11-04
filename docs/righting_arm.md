# Righting arms: `GZ(phi)` measured against a closed form and an independent section solve

`docs/stability.md` closes with what it does not do: *"`GM` is the slope of `GZ(phi)` at
`phi = 0`. `GZ` needs a heeled displacement, which is a different integration."* This page
is that second integration, and it is the last open item of the stability backlog.

    GZ(phi) = y_B(phi) - y_G(phi)

with the body heeled, the waterline **solved** so the displaced volume stays at the loading
displacement, and `B` the centroid of the heeled volume. `G` is fixed in the body at `KG`
above the keel. Heel is positive to starboard.

Everything below is produced by `python examples/righting_arm/run_example.py --output
docs/righting_arm`, which writes `results.json`, `report.md` and two figures, and **exits
non-zero if one of its fourteen claims stops holding**. The numbers here are read from that
JSON; `tests/test_righting_arm.py` pins 39 cases of the same content.

## Why the waterline is an unknown

A body that keeps its displacement while it heels cannot keep its waterline: the low side
goes deeper, the high side shallower, and only a change of immersion preserves the
displaced volume. So `GZ` is not a property of a mesh rotated about a fixed plane. The
module bisects the waterline height until `|V - V_loading| <= 1e-13 V_loading` and carries
the residual on every arm, so a curve whose solve failed cannot be mistaken for a converged
one. The lever itself does not care where the plane ends up — a vertical shift moves neither
`y_B` nor `y_G` — which is why the solve only has to make the volume right.

## No cap, and therefore no degeneracy

The immersed moments are fluxes of **z-directed fields that all carry a factor `(z - w)`**:

| quantity | field `F` | `div F` | flux through the waterplane |
|---|---|---|---|
| `volume` | `(0, 0, z - w)` | `1` | `0` |
| `∫x dV` | `(0, 0, x(z - w))` | `x` | `0` |
| `∫y dV` | `(0, 0, y(z - w))` | `y` | `0` |
| `∫z dV` | `(0, 0, (z - w)^2/2)` | `z - w` | `0` |

The cap's normal is `z` and `z = w` on it, so **the cap contributes exactly zero and never
has to be triangulated**. That is worth two things. First, no cut-segment chaining: the
waterplane integrator of `stability.py` chains segments into loops because it *needs* the
waterplane area, but a cap triangulated from chained loops introduces its own failure modes.
Second, the degenerate configurations become free: a face lying exactly in the clipping
plane, or a plane running exactly along a mesh edge, contributes zero to all four integrals
above instead of doubling the cap. (Empirically that matters — asking the waterplane
integrator for a section at `z = -0.8` on a hull whose vertex rings sit at `-0.8` raises
"160 segments that do not close"; the same plane through the immersed-volume integrator is
unremarkable.)

The trap to avoid is putting `x` on the x-component: `F = (x(z - w), 0, 0)` has divergence
`z - w`, not `x`, and its flux through the cap does **not** vanish — the horizontal moments
then come out as `-w · volume`, which is exactly zero at `w = 0` and silently wrong at every
heel that needs a sinkage. That bug was written, caught by the barge's analytic centroid, and
is the reason the table above is in this page.

Each clipped polygon is integrated exactly: a triangle clipped to a half-space is a convex
polygon of at most four vertices, fan-triangulated, and each triangle is integrated with the
three-point edge-midpoint rule, which is exact for the degree-two integrands above. The only
discretisation left is the one already baked into the mesh.

## The barge: where the arm has a closed form

`L = 10`, `B = 4`, `T = 2`, `D = 6` m, `KG = 0.9` m, `V = 80 m³`, `KB = 1`, `BM = 2/3`,
`GM = 0.7667` m. For a body whose waterline still spans the full beam,

    GZ = sin(phi) · (GM + BM tan(phi)² / 2)

is exact (the wall-sided formula), and `wall_sided_validity` reports where that stops: bilge
emergence at `atan(2T/B) = 45°` for this body, deck immersion at `atan(2(D-T)/B) = 63.4°`, so
bilge emergence binds at `45°`.

| heel [deg] | GZ mesh [m] | wall-sided [m] | 2D section [m] | mesh − wall | mesh − 2D |
|---|---|---|---|---|---|
| 2 | 0.026770467 | 0.026770467 | 0.026770467 | `3.4e-16` | `0` |
| 10 | 0.134929913 | 0.134929913 | 0.134929913 | `-8.1e-16` | `0` |
| 30 | 0.438888889 | 0.438888889 | 0.438888889 | `-8.8e-15` | `0` |
| 44 | 0.748507086 | 0.748507086 | 0.748507086 | `7.8e-16` | `0` |
| **46** | 0.807738096 | 0.808615392 | 0.807738096 | `-8.8e-04` | `0` |
| 60 | 1.233691384 | 1.529978213 | 1.233691384 | `-2.96e-01` | `0` |
| 90 | 2.100000000 | — (no claim) | 2.100000000 | — | `0` |

Worst deviation inside the regime: `8.8e-15 m`. At `46°` the closed form is already
`8.8e-4 m` off and at `60°` it is `0.296 m` off (`19.4 %` of the arm), which is asserted as
a claim: a closed form has a *regime*, and the boundary is not cosmetic.

Past that boundary the reference is an **independent two-dimensional computation**: the
section of a prism is a rectangle, so the heeled problem reduces exactly to clipping a
polygon and taking its shoelace centroid, with its own bisection on the section area. It
reaches the same lever through a different dimension and a different integral (Green's
theorem in the plane instead of a divergence theorem in space) and agrees at every heel from
`0°` to `90°` to `2.9e-14 m`.

The rest of the barge's claims:

| quantity | measured | reference | residual |
|---|---|---|---|
| displacement held at every heel | `max rel` | `V = 80 m³` | `1.0e-13` |
| `GZ(0)` and oddness (`GZ(-φ) = -GZ(φ)`) | — | `0` | `2.2e-16 m` |
| `KG + 0.4 m` shift (`ΔGZ = -0.4 sin φ`) | — | identity | `1.1e-16 m` |
| resting heel, `KG = 1.9` (`GM = -0.2333 < 0`) | `39.91603°` | `atan√(-2GM/BM) = 39.91788°` | `-1.8e-3°` |
| area under `GZ` to `45°` | `0.26499310 m rad` | `GM(1-cosφ) + BM/2(cosφ + secφ - 2) = 0.26499158` | `1.5e-6` |

The area residual is *entirely* quadrature: the trapezoid rule on the closed-form arm over
the same 0.25° grid differs from the exact integral by `1.513e-6`, i.e. the measured curve
and the closed form integrate identically. The resting-heel residual is grid-limited too
(the zero crossing is interpolated between 0.5° samples), not an error in the arm.

## A hull with topsides: no closed form, two other checks

The shipped families cap the surface at the design waterline, which leaves **no freeboard**
(see the trap below), so the example generates the deeper body — `draft = 1.6 m` — and floats
it at `z = -0.9 m`: 0.7 m immersed, 0.9 m of freeboard at the centreline. `V = 4.237853 m³`,
`KB = 0.428123`, `BM = 0.222483`, `GM = 0.350607 m` at `KG = 0.30 m`. The topsides are the
analytic continuation of the same `√(z/draft)` section law, not a designed deck — an
idealisation, stated here rather than buried.

**The slope at zero heel.** `GZ/sin(φ)` → `GM` with an `O(φ²)` defect, removed by
Richardson from 2° and 4°: `0.3506131 m` against `KB + I_T/V - KG = 0.3506066 m`, a relative
residual of `1.9e-5`. That `GM` is assembled from `tensorcad.stability`'s waterplane
integrator (Green's theorem on the plane section) and the module's own immersed volume, so
the two halves of the check share no code with the heeled-displacement integration. The same
`KG`-shift identity holds on the hull to `1.7e-16 m`.

**Convergence.** The arm at 20° of heel: `0.125390727` (41×13), `0.122578961` (81×25),
`0.121679056` (161×41) — observed order `2.04`, i.e. the same second order the hydrostatics
shows, but the *constant* is larger: the shipped default grid is **+3.05 %** high on the arm
at 20° where it is `-0.88 %` on the volume. A righting arm is a first moment of the immersed
volume divided by the volume, so it inherits both errors and adds a lever; the arm is the
quantity to refine for, not the displacement.

## The trap: a waterline-capped hull has no large-angle curve

`generate_hull` and `generate_wigley` cap the surface *at the design waterline*. A surface
capped there has no freeboard, so its displacement **is** its whole volume, and no heeled
waterline can hold that: displacing the entire body means submerging it, and the arm of a
fully submerged body is identically zero (its `B` is fixed in the body, so `y_B = y_G`). The
shipped hull therefore has no large-angle curve at all, and asking for one raises:

```
the requested displacement 7.446554097149238 m^3 is the whole volume of the body
(7.446554097149233 m^3), so the body has no freeboard: a closed surface capped at its
waterline cannot heel at constant displacement ... Generate a deeper body and float it at a
lower waterline_z, e.g. HullParameters(draft=2 * t) with waterline_z=-t
```

A row of zeros would have been the wrong answer: it looks like a result and carries none.

**Second, smaller trap**: `immersed_volume` is the *general* routine — it takes the waterline
as an argument. `compute_hydrostatics` sums signed tetrahedra of the whole closed surface and
is only the submerged volume when the mesh is already cut at its waterline
(`docs/hydrostatics_verification.md`), so it cannot be pointed at a heeled body by passing a
different `waterline_z`. On a surface that *is* cut there, the two agree to `1e-12` — which
is a test.

## What this is not

No free-surface correction (a slack tank stays the caller's problem, as in `stability.py`),
no trim coupling (the body is heeled about a fixed longitudinal axis; a real heel changes
trim too), no damage or flooding case, and no statutory verdict: `IS Code 2008` needs wind
heeling arms, roll amplitudes and tank free surfaces, none of which is here.
`range_of_stability` reports where the computed `GZ` first returns to zero inside the
computed grid and refuses to guess past it.

## Reproducing

```bash
python examples/righting_arm/run_example.py --output docs/righting_arm
python -m pytest tests/test_righting_arm.py      # 39 tests
```

| artefact | content |
|---|---|
| `docs/righting_arm/results.json` | every measured number, the claim table, and both curves |
| `docs/righting_arm/report.md` | the same tables rendered |
| `docs/righting_arm/gz_curves.png` | barge against the closed form and the 2D section solve; hull at two loadings against `GM sin φ` |
| `docs/righting_arm/deviation.png` | the two differences at `1e-18` resolution, and the arm's second-order convergence |