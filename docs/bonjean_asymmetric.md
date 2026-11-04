# Bonjean closed forms for the asymmetric hull family

The asymmetric generator (0.17.0) gave the parametric hull a taper exponent per half
(`fore_exponent` / `aft_exponent`) and `analytic_volume`, `analytic_waterplane_area`,
`analytic_longitudinal_centre` and `analytic_wetted_surface_area` followed. The
**Bonjean** module did not: `bonjean.analytic_section_area` still hard-wired
`BREADTH_EXPONENT`, so the one function that answers *how much area sits under this
waterline at this station* could not be checked against a hull whose halves differ —
the half of the family the split introduced had no reference value at all. This page
closes that verification debt and adds the second closed form the table was missing.

Every number below is written by `python examples/bonjean_asymmetric/run_example.py`
into `docs/bonjean_asymmetric/{results.json,plot.png}` (the example exits non-zero if
a claim stops holding); `tests/test_bonjean_asymmetric.py` pins the same behaviour
with 49 tests.

## Case

The parametric hull family, on the ITTC-like proportions used throughout the
verification pages:

| parameter | value |
|---|---|
| length `L` | 10.0 m |
| beam `B` | 2.0 m |
| draft `T` | 0.8 m |
| fore exponent | 0.5 (full forebody) |
| aft exponent | 0.9 (fine afterbody) |
| grids | 41×13, 81×21, 161×41 (station rings × waterline levels) |

Half breadth is `(B/2) sin(pi x/L)^p sqrt(zeta/T)` with `p = fore` forward of midship
and `aft` abaft it; the *area* closed form is therefore per half as well:

    A(x, zeta) = (2/3) B sin(pi x/L)^p zeta^(3/2) / sqrt(T),
    zeta = z + T measured up from the keel.

The second closed form is new here. Because the vertical profile is `sqrt(zeta/T)` at
every station, the area is `c zeta^(3/2)` and its first moment about the keel is
`(2/5) c zeta^(5/2)`, so the section centroid is

    z_c = -T + (3/5) zeta

— independent of the station *and* of the taper exponent. At `zeta = T` that is the
`(3/5) T` of `analytic_vertical_centre_buoyancy`, reached per station instead of per
hull. A centroid that drifts with `x` or with `p` is therefore a cutting error, not a
formula error, which is what makes it a sharp test.

## Method

1. `bonjean.analytic_section_area` gains keyword-only `fore_exponent`/`aft_exponent`;
   the exponent is chosen by `2 * station >= L` — the generator's own rule — and both
   default to `BREADTH_EXPONENT`, so every previously published digit is reproduced by
   the same expression (a test asserts *bit* equality, not 1e-15 closeness).
2. `bonjean.analytic_section_centroid_z(draft, waterline_z)` — the `3/5` closed form.
3. `BonjeanTable.integrated_first_moment` and `BonjeanTable.longitudinal_centre` — the
   `int x A dx` companion of `integrated_volume`, so the longitudinal centre of the
   *Bonjean areas* becomes measurable without a single surface integral. That is an
   independent measurement of the same number `compute_hydrostatics` and
   `analytic_longitudinal_centre` report.
4. Mesh measurement uses the existing exact-cut machinery (`section_properties`:
   plane cut, edge-key chained loops, Sutherland–Hodgman clip, shoelace), unchanged.
5. Still-water comparison at `zeta = T` (design waterline), Simpson along 41 stations
   for the Bonjean reductions.

## Results

### Section area against the split closed form

Worst signed relative error over 5 stations (0.1L, 0.3L, 0.5L, 0.7L, 0.9L) × 4
waterlines (-0.6, -0.4, -0.1, 0.0 m):

| grid | worst area error | worst centroid error |
|---|---|---|
| 41×13 | −5.308 % | +0.939 % |
| 81×21 | −2.539 % | +0.457 % |
| 161×41 | −0.924 % | +0.169 % |

Both ladders are monotone, and both worst cases are the shallowest waterline — the
piecewise-linear skin standing in for `sqrt(zeta)` where only three of thirteen levels
are immersed. All area errors are **negative**: inscribed chords remove area, never add
it. That is the same −5.31 % / −0.924 % the symmetric family reports (`docs/bonjean.md`),
which is the expected result: the split changes the *longitudinal* taper, and the area
error is a *vertical* discretisation effect, so it cannot tell the halves apart.

### The centroid closed form

| grid | measured `z_c` at midship, design waterline | closed form `-T + 3/5 T` | relative |
|---|---|---|---|
| 41×13 | −0.316995 m | −0.32 m | +0.939 % |
| 81×21 | −0.318539 m | −0.32 m | +0.457 % |
| 161×41 | −0.319461 m | −0.32 m | +0.169 % |

The independent part of this test is the *spread*: at 41×13 the centroid measured at
`x = 2, 3.5, 5, 6.5, 8 m` varies by **4.2e-12 m** (3.8e-13 m at 161×41) — the cut is
exact, and the value is one number for the whole hull, as the closed form says it must
be despite the two halves having different exponents. All of the 0.94 % is therefore a
vertical-skin effect, none of it a cutting or averaging artefact. A barge control
(`centroid_z = -0.625 m` at a -0.25 m waterline on a box from z = -1) is exact to
1e-15, confirming the sign convention the hull numbers are compared in.

### Displacement and longitudinal centre from the Bonjean table

Analytic references at the split exponents: volume 7.572387 m³, LCB 5.257824 m.

| grid | Bonjean volume | vs closed form | surface-integral volume | vs closed form | Bonjean LCB | vs closed form |
|---|---|---|---|---|---|---|
| 41×13 | 7.512605 m³ | −0.790 % | 7.500729 m³ | −0.946 % | 5.254159 m | −0.0697 % |
| 81×21 | 7.540882 m³ | −0.416 % | 7.540966 m³ | −0.415 % | 5.254159 m | −0.0697 % |
| 161×41 | 7.557047 m³ | −0.203 % | 7.561214 m³ | −0.148 % | 5.254159 m | −0.0697 % |

Two things in that table are worth more than the percentages.

*The LCB does not move with the grid* — to 1e-9 m across three grids, while the volume
error falls by a factor of four. The Bonjean LCB is a ratio of two integrals of the same
areas, so the vertical discretisation — a near-uniform multiplicative bias on every
`A(x)` — cancels. The residual is therefore **station quadrature, not the mesh**, and
that can be proved directly: the same Simpson rule applied to the *closed-form* areas at
the same 41 stations gives LCB = 5.254159343 m, i.e. within **8.3e-9 m** of the mesh
number. Refining only the stations of the closed form walks the residual down:

| stations on the closed form | LCB | vs analytic |
|---|---|---|
| 41 | 5.254159 m | −0.0697 % |
| 81 | 5.256520 m | −0.0248 % |
| 161 | 5.257360 m | −0.0088 % |

The ratio 2.81 / 2.81 between successive errors is what a `x·A(x)` integrand does when
`A ~ sin(pi x/L)^0.5` leaves a `sqrt`-type endpoint singularity: Simpson degrades from
second order to about 1.5. So the -0.0697 % is fully accounted for, and the chain
Bonjean → volume → LCB carries **no additional model error** beyond the station rule.

*The centre is forward of midship* — 5.254159 m against 5.0 m — because the afterbody
(0.9) is thinner than the forebody (0.5). The measurement could not have produced that
number if the exponent split had silently collapsed to the symmetric hull, and a direct
check compares the two halves' areas at mirrored stations.

### Exact identities

* Mirror: `LCB(p, q) + LCB(q, p) = L` because the two moments telescope into
  `h(p) + h(q)` over the same denominator. Residual 1.78e-15 m at (0.5, 0.9) and
  exactly 0 at (0.4, 1.2) and (1.0, 0.6).
* Midship is a fixed point of the split: `sin(pi/2) = 1` for every exponent, so
  `A(L/2, zeta) = 0.692820323027551 m²` exactly at `zeta = 0.6 m`, asserted bit-exactly.
* The station-wise closed form integrates to `analytic_volume(..., fore, aft)` with a
  4001-station Simpson to 2e-6 relative — the per-station and whole-hull closed forms are
  one model, so the new function cannot drift from the volume it is compared against.

## Error verdict

| quantity | reference | error | < 3 % ? |
|---|---|---|---|
| section area, 161×41 | split closed form | −0.924 % | yes |
| section area, 41×13 | split closed form | −5.308 % (worst, `z = -0.6 m`) | **no** — vertical skin, converges |
| section centroid, 41×13 | `-T + 3/5 zeta` | +0.939 % | yes |
| section centroid, 161×41 | `-T + 3/5 zeta` | +0.169 % | yes |
| displacement from Bonjean, 41×13 | `analytic_volume` split | −0.790 % | yes |
| displacement from Bonjean, 161×41 | `analytic_volume` split | −0.203 % | yes |
| LCB from Bonjean, all grids | `analytic_longitudinal_centre` split | −0.0697 % | yes |
| LCB, mesh vs closed-form 41-station quadrature | same rule on closed form | 8.3e-9 m abs | yes |

The single > 3 % entry is the shallow-waterline section area on the coarsest grid, which
is a known and documented property of the piecewise-linear skin (`docs/bonjean.md`
reports the identical −5.31 % / −0.924 % for the symmetric family) and is reported rather
than hidden; every quantity that enters a design decision is under 1 %.

## Caveats

* The exponent selector is a *convention* (`2 * station >= L` is forward), matching
  `_breadth_exponent`; a hull whose exponents agree with each other is unaffected, but a
  caller who builds a hull with `fore < aft` inherits a centre **aft** of midship — the
  mirror identity says the two directions are exact counterparts, and both are tested.
* `analytic_section_centroid_z` is the *watertight-hull* centroid: it assumes the full
  `sqrt` profile is immersed, so it is a reference for a section, not for a partially
  flooded compartment after damage.
* The Bonjean reduction uses the design waterline at the top of the section; the tables
  themselves interpolate linearly between tabulated waterlines, so a query between levels
  carries an `O(h²)` interpolation error that this page does not measure (the closed form
  does, via the vertical ladder above).