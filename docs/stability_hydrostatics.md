# Waterplane inertia and initial stability

> **Merge note (PR #39, ``agent/m2925``).** This page and
> `examples/stability/run_hydrostatics_example.py` describe the
> *hydrostatics-first* route to KM and GM: the two waterplane second moments are
> carried on `Hydrostatics` itself (flux identity, so no waterline cap is needed)
> and `initial_stability_from_hydrostatics` turns that record into KB/BM/KM/GM
> without touching a mesh. `docs/stability.md` describes the *mesh-first* route
> that landed on main first (plane cut + Green's theorem, exact for the chords the
> mesh has, and the only one of the two that handles a heeled or multi-loop
> section). Both are exported; `tests/test_stability.py` checks that they agree on
> the same hull, and the closed forms below are the common reference for both.

The hydrostatics module has reported displacement, areas and centres since 0.4.0. It
had no waterplane **second** moments, so the one number a naval architect asks for
first — the metacentric height — could not be produced, and every page in this
repository said so (`docs/hydrostatics_verification.md`, README "what this is not").
This release closes that gap with the two integrals, and turns them into KB, BM, KM
and GM in `tensorcad.stability`.

It closes it with a reference attached, not just a formula: the same two integrals
are available in closed form on the analytic hull surface, so every mesh number has
something independent to be measured against.

## What is computed

| Quantity | Symbol | Where |
|---|---|---|
| transverse waterplane inertia, `∫ y² dA` about the centreline | `I_T` | `Hydrostatics.transverse_waterplane_inertia` |
| longitudinal waterplane inertia, `∫ (x-LCF)² dA` about the flotation axis | `I_L` | `Hydrostatics.longitudinal_waterplane_inertia` |
| height of the centre of buoyancy over the keel | `KB` | `HydrostaticsInitialStability.height_of_centre_of_buoyancy` |
| transverse / longitudinal metacentric radius, `I/∇` | `BM_T`, `BM_L` | `…_metacentric_radius` |
| height of the metacentre over the keel, `KB + BM` | `KM_T`, `KM_L` | `…_metacentre` |
| transverse / longitudinal metacentric height, `KM - KG` | `GM_T`, `GM_L` | `…_metacentric_height` |

```python
from tensorcad import generate_hull, initial_stability_from_hydrostatics

model = generate_hull()                      # L=10, B=2, T=0.8, 41x13
stability = initial_stability_from_hydrostatics(model.hydrostatics, height_of_gravity=0.5)
stability.transverse_metacentric_radius      # 0.438618 m
stability.height_of_transverse_metacentre    # 0.921623 m
stability.transverse_metacentric_height      # 0.421623 m
stability.is_stable_in_roll                  # True
```

`KG` is an **input**. TensorCAD has no weight schedule, so any KG it derived would
be a guess wearing a result's clothes; the caller supplies the loading, the code
supplies the geometry.

## How the inertias are integrated off the mesh

No waterline cap is needed, and no polygon has to be reconstructed. For any closed
surface the flux of `F = (0, 0, f)` is zero whenever `f` does not depend on `z`, and
every wetted face lies below the waterline, so

```
∫_{waterplane} f dA  =  - ∫_{wetted} f n_z dA
```

With `f = y²` this is `I_T`; with `f = x²` it is the moment about the origin, which
the parallel-axis theorem converts into `I_L` about the centre of flotation:

```
I_L = ∫ x² dA - A_w · LCF²
```

This is the same identity the module already used for the waterplane *area* (with
`f = 1`), so the waterplane quantities stay consistent with each other by
construction, including on a surface whose cap is missing.

The face integrals are exact: `x` and `y` are linear over a triangle, so

```
∫_face f² dA = A/6 · (f₁² + f₂² + f₃² + f₁f₂ + f₂f₃ + f₃f₁)
```

and `face_normal` is twice the area, which makes the accumulated `n_z · ∫f²` a
single multiply. The only error left in the result is the discretisation of the
*surface*, which is what the reference measures.

Two honest limits:

* Faces that straddle the waterline are not sliced (the same restriction the rest of
  the module has), and
* an **open** surface returns `nan` for both inertias rather than a finite number
  that would only reflect where the shell happens to stop — the flux identity needs
  a closed surface to hold at all. `basis` says `open surface (volume not
  trustworthy)` in that case, and `initial_stability_from_hydrostatics` raises.

The sign of the accumulated flux is folded out instead of assumed, so a mesh built
with the opposite winding still reports non-negative moments.

## The reference on the analytic surface

The generating surface has half-breadth `b(x) = (B/2)·sin(πx/L)^e` with `e = 0.7`
(see `docs/hull_analytic_reference.md`), and the waterplane is the band
`|y| ≤ b(x)`. Both integrals collapse:

| Quantity | Closed form | Method |
|---|---|---|
| `I_T` | `(2/3)·(B/2)³·L·∫₀¹ sin(πξ)^{3e} dξ` | gamma identity — exact |
| `I_L` | `B·L³·∫₀¹ sin(πξ)^e (ξ-1/2)² dξ` | tanh-sinh quadrature, self-reported |
| `KB` | `0.6·T` (`VCB = -0.4T`, keel at `-T`) | exact |
| `BM_T`, `BM_L` | `I_T/∇`, `I_L/∇` with `∇ = B·L·T·I·2/3` | exact ratio |
| `KM_T` | `0.6·T + BM_T` | exact |

The longitudinal integrand keeps the `ξ^{-0.3}` endpoint of the taper, so it goes
through `tanh_sinh_quadrature` and carries `value_at_previous_level` with it; the
difference between levels 6 and 7 is `0.0` in float64, i.e. the reference is exact
to the resolution of the mesh values compared against it. The transverse one did not
need quadrature at all: `∫₀¹ sin(πξ)^{2.1} dξ` is the same gamma identity at a
different exponent, and `tests/test_stability.py` checks it against a direct
integration of `∫(2/3)b³ dx` to `1e-12` relative.

## Verification

`python examples/stability/run_hydrostatics_example.py` runs the resolution ladder
`21×9 … 321×81`, compares six quantities against the closed forms, writes
`docs/stability_hydrostatics/{results.json,report.md,plot.png}` and exits non-zero if a claim
fails. Measured on the default hull (L=10, B=2, T=0.8):

| Quantity | 41×13 (default) | 161×41 | 321×81 |
|---|---|---|---|
| `I_T` | −0.139 % | −0.009 % | −0.002 % |
| `I_L` | −0.604 % | −0.064 % | −0.020 % |
| `BM_T` | +0.748 % | +0.127 % | +0.046 % |
| `KM_T` | +0.684 % | +0.120 % | +0.043 % |
| `GM_T` (KG = 0.5 m) | +1.508 % | +0.263 % | +0.095 % |
| `∇` | −0.880 % | −0.136 % | — |

Every quantity is inside 1.5 % at the default grid and inside 0.3 % at 161×41, so
the 3 % acceptance band holds with room to spare. The table is regenerated by the
example; `docs/stability_hydrostatics/report.md` is the same table as emitted by the run.

`GM_T` is the weakest of the six, and for a structural reason worth stating: `BM_T`
is `I_T/∇`, so it inherits the *volume* error with the opposite sign. The mesh
under-integrates the volume by 0.88 % and the transverse inertia by only 0.14 %, so
`BM_T` comes out about 0.75 % high. A user who wants GM to 1 % needs a grid where
the volume error is also near 1 % — 81×21 in this ladder.

Two more checks live in `tests/test_stability.py`, and they are the ones that would
catch a *systematically* wrong identity rather than a discretisation error:

* a **rectangular barge** built vertex by vertex has `I_T = L·B³/12`,
  `I_L = B·L³/12`, `KB = T/2` and `GM = BM_T` when `KG = KB`. All four come out to
  `1e-12` relative on three different tessellations, including a single-element box
  — the flux identity is exact on a flat waterplane, so any leftover would be a bug
  in the bookkeeping, not in the physics.
* a **circular waterplane** (vertical cylinder) approaches `πR⁴/4` from below as the
  segment count doubles, with the error ratios falling by 4 each time (the `O(1/n²)`
  of an inscribed `n`-gon) and no face of it being horizontal-but-for-the-bottom.

## What this is not

Even keel, upright, small angles, one displacement, no free-surface correction, no
damage case, no `KG` model, no righting-arm curve, no IMO or class-rule criterion.
`GM_T > 0` is necessary for roll stability at small angles and says nothing about
the large-angle range: a hull with a large GM can still have a poor range because
the deck edge immerses early. The righting-arm curve needs `GZ(φ)`, which needs the
buoyancy centroid at heel, which needs waterline slicing of the mesh — the next
piece of work, not this one.