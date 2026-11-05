# Righting arms at finite heel: `GZ(φ)`

`tensorcad.stability` ends its own docstring with the sentence this page answers:

> *Nor does it give a righting-arm curve `GZ(φ)` for finite heel: `GM` is the slope of that
> curve at `φ = 0`, and stability at large angles is a different integration (a heeled
> displacement integration, not a waterplane moment).*

`tensorcad.heel` is that integration. Everything on this page is produced by
`python examples/gz_curve/run_example.py --output docs/gz_curve`, which writes
`results.json`, `report.md` and three figures and **exits non-zero if one of its twelve
claims stops holding**; `tests/test_heel.py` pins 30 cases of the same content.

## What is computed

A body of displacement ∇₀ is heeled by `φ` about the longitudinal axis **at constant
displacement**. The water surface is a plane whose normal in the body frame is

```
u(φ) = (0, −sin φ, cos φ)          # positive φ = starboard side down
```

so the submerged side is `u·r ≤ c`, and `c` is *solved* for each angle so that the
submerged volume equals ∇₀. That solve is the physical content: cutting at the upright
waterline and rotating the answer is a different (wrong) question. The centre of buoyancy
`B` comes out of the same integral, and

```
GZ(φ) = y_B cos φ + (z_B − z_G) sin φ,     z_G = z_keel + KG
```

with `B`, `G` in the body frame and `y` to starboard. `KG` is a required keyword, for the
same reason as in `initial_stability`: it is a loading condition, not a shape property.

**Why it is exact, not merely converged.** Each triangle is clipped against the
half-space (Sutherland–Hodgman on the three signed vertex heights, no geometric tolerance
anywhere), and each clipped polygon is split into tetrahedra with a single apex **on the
water surface**, `o = c·u`. The tetrahedra of the waterplane lid are then degenerate and
contribute exactly zero, so the signed tetrahedron sum over the clipped hull triangles
alone is the exact submerged volume and first moment:

```
V = (1/6) Σ_faces det[a−o, b−o, c−o],      ∫ r dV = Σ_tet V_tet (o+a+b+c)/4
```

This is an identity for a closed surface, so a heeled waterline that leaves the deck, a
non-convex section, or an asymmetric body needs no special case. The only remaining error
is the one-dimensional root-find on `c` (bracketed regula falsi, Illinois damping, residual
`≈ 1e-12 ∇₀`).

## The barge: the integration is exact

`L = 20 m`, `B = 6 m`, `T = 2 m`, depth `4 m` (freeboard 2 m), `KG = 2 m`. The deck edge
immerses at `atan(2f/B) = 33.69°`, and below that angle the wall-sided closed form
`GZ = sin φ (GM + ½ BM tan²φ)` is exact:

| φ [deg] | GZ mesh [m] | wall-sided [m] | relative error | solved offset c [m] |
|---|---|---|---|---|
| 5 | +0.044078206 | +0.044078206 | `+6.3e-16` | `−4.7e-16` |
| 10 | +0.090873287 | +0.090873287 | `+3.1e-16` | `0` |
| 20 | +0.204991739 | +0.204991739 | `+8.1e-16` | `0` |
| 30 | +0.375000000 | +0.375000000 | `−1.0e-15` | `+9.5e-16` |

Worst relative deviation below the deck edge: **`1.04e-15`**; the largest excursion of the
solved waterline from `z = 0` is **`9.47e-16 m`**, which is the wall-sided property (`V`
linear in `c`) *coming out of the solve* rather than being assumed. The elliptical prism
in the tests repeats the same check with a waterplane whose `I_T` is not `L B³/12`
(`π a b³/4`), agreeing to `1e-11 m` — so `BM` really comes from the section, not from a
rectangle shortcut.

## Where the reference stops: the shallow barge

Same hull at depth `3 m` (freeboard 1 m), `KG = 2.2 m` → `GM = +0.300 m`, deck edge at
`18.43°`. Past the deck edge the unbounded-wall formula over-predicts — by **2.424 m at
60°** — and the measured curve turns over:

| quantity | value |
|---|---|
| `GZ_max` | `0.1470 m` at `25°` |
| vanishing angle (`GZ = 0` again) | `45.95°` |

So the wall-sided formula is quoted only where it is a reference, and the change of regime
is *measured*, not assumed. `vanishing_angle_deg` is a linear interpolation between two
samples and is `None` when the curve has not crossed by the last sample — no invented
capsizing angle.

## The sphere: an exact large-angle reference with no deck

For a sphere, heel changes nothing: `B` stays on the vertical through the centre, so
`GZ = sin φ (z_C − z_G)` **at any angle**. With `R = 1 m`, draft `0.6 R`, `KG = 0.3 m`:

| φ [deg] | GZ mesh [m] | `sin φ (z_C − z_G)` [m] | error [mm] | offset c [m] |
|---|---|---|---|---|
| 10 | +0.121273414 | +0.121553724 | `+0.280` | `−0.0061` |
| 40 | +0.449100805 | +0.449951327 | `+0.851` | `−0.0934` |
| 60 | +0.605457528 | +0.606217783 | `+0.760` | `−0.1996` |
| 80 | +0.689065618 | +0.689365427 | `+0.300` | `−0.3300` |

Worst `0.851 mm`. That residual is **the faceting of the mesh, not the integration**, and
the ladder removes the doubt — the error falls by `3.81×` then `4.02×` per halving of the
facet size `h = π/rings`:

| rings | error at 60° [m] | factor per doubling |
|---|---|---|
| 12 | `2.896e-03` | — |
| 24 | `7.603e-04` | `3.81×` |
| 48 | `1.892e-04` | `4.02×` |

Second order, as a faceted sphere should be. Note the offset `c` moving 0.32 m between 0°
and 80° *while ∇ is fixed*: the solve is doing real work, not reusing the upright plane.

## Consistency with `tensorcad.stability`

The two modules compute different objects, so agreeing is evidence:

| check | value |
|---|---|
| 1° chord slope vs `GM` from the waterplane moments | `0.500203` vs `0.500000 m/rad`, `+4.06e-04` relative (the wall-sided `tan²` term a 1° chord already carries) |
| `GM` from the waterplane moments vs `KB + BM − KG` | `0.500000` vs `0.500000 m` (exact) |
| odd symmetry `|GZ(φ) + GZ(−φ)|` | `2.22e-16 m` |
| constant displacement, all sampled angles | `9.79e-13` relative |

## One trap worth stating: `GM` is not `initial_stability`'s `GM` for a freeboard body

`initial_stability` divides `I_T` by *the whole closed volume* and takes `KB` from that
body's centroid — right for a mesh capped at the waterline, wrong for a body with
freeboard. On the barge (240 m³ displaced inside a 480 m³ closed body) it reports
`GM = 0.75 m` where the truth is `0.50 m`, from two errors that only cancel by accident.
`heel.gz_curve` therefore builds `GM` from its own *submerged* upright body and
`waterplane_moments`; `tests/test_heel.py` pins both numbers so the difference cannot
silently disappear.

A related edge case, found by the sphere case and **not** fixed here: `waterplane_moments`
double-counts a section whose plane passes exactly through an interior vertex ring (a UV
sphere with a ring on the equator gives `A = 6.265257 m²` instead of `π` and `loops = 2`,
because the ring is traced once from above and once from below). `GZ` is unaffected — it
never uses the plane section — but a `GM` read from such a waterline is `2×` off. Water
lines that pass through a ring of vertices are therefore avoided in this example, and the
defect is recorded here for `tensorcad.stability` rather than patched under a GZ change.

## What this does not do

- **No free-surface correction**: tanks are a loading property (`Σρ i_x / ∇`), not a shape
  one; subtract it at the caller who knows the tankage.
- **No damage stability, no wind or launching heeling moment, no class rule.**
- `G` is on the centreline; an off-centre weight needs a `y_G` input.
- **The closed body *is* the hull.** A mesh capped at the design waterline has zero
  freeboard, so no heeled equilibrium exists and `HeelError` is raised instead of a number.
  Computing a ship's `GZ` needs the body closed at the deck or bulwark — the next piece of
  geometry work, not a tolerance to loosen.
- Angles are sampled: `vanishing_angle_deg` interpolates between two samples and is `None`
  when there is no crossing in range.