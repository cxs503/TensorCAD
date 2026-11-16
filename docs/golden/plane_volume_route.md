# The waterplane-closure route (open shells): validated integrator, unresolved 5415 gap

Anchor: `GAP-CAD-D2-GOLDEN-REGRESSION`. The golden library's one remaining `WITHHELD` is
the DTMB 5415 displacement: the shipped shell is an open sheet (one boundary loop of 722
edges), so `total_signed_volume` on it measures nothing, and the case publishes the
station-area displacement instead (`tensorcad.dtmb5415.DTMB5415Space.hydrostatics`).

This note adds the missing *independent* route — a surface integral that needs no mesh
repair — proves it against closed-form geometry, and reports honestly that it does **not**
yet close the 5415 gap. The golden baseline was deliberately **not** re-blessed.

## 1. Method — closure supplied analytically, not by repair

For a shell open at the top, the column of water below a plane `z = T` is bounded by:

* the part of the shell with `z <= T` (lateral + bottom boundary),
* the waterplane patch at `z = T` (planar, outward normal `+e_z`),
* every other opening of the sheet (transom plane, longitudinal stem seam), which is
  **vertical**, so it contributes nothing to either integral below.

With `g = ∫ n_z dA` and `h = ∫ z n_z dA` over the clipped shell, the divergence theorem on
that closed surface (`div(z e_z) = 1`, `div(e_z) = 0`) gives

    A_waterplane = -g            V = h - T·g

with the sign of the shell's winding resolved from `g` (reported as `orientation`).
Triangles that the plane crosses are clipped exactly (Sutherland-Hodgman, fan
triangulation), so the result is exact for the given tessellation: no voxelisation, no
volume estimate, no watertight-repair precondition. Implementation:
`tensorcad.mesh.volume_below_plane` → `PlaneVolume`.

## 2. Validation against closed forms — `tests/test_mesh_plane_volume.py`

| case | quantity | computed | closed form | verdict |
| --- | --- | --- | --- | --- |
| box 2×3×4, deck removed, plane at `z=4` | volume | 24.0 | 24.0 | exact (rel 1e-15) |
| same | waterplane area | 6.0 | 6.0 | exact |
| same | submerged area | 46.0 | 6+40 | exact |
| same, plane at `z=2` (8 triangles clipped) | volume / area | 12.0 / 6.0 | 12.0 / 6.0 | exact |
| same, winding inverted | volume / area | identical, `orientation=-1` | — | sign resolved, not assumed |
| closed box at its top plane | volume / area | 24.0 / 0.0 | 24.0 / 0.0 | closure adds nothing |
| lower hemisphere `r=1` (96×48), plane `z=0` | volume | 2.09066 | 2/3π = 2.09440 | −0.18 % (inscribed polyhedron, converges: −0.71 % at 48×24) |
| same | waterplane area | 3.13935 | π = 3.14159 | −0.07 % (−0.29 % at 48×24) |
| shell entirely above the plane | volume / area | 0 / 0 | 0 / 0 | fail-closed, no NaN |

`pytest tests/test_mesh_plane_volume.py -q` → **7 passed** (0.63 s).

## 3. Result on the shipped 5415 sheet — and why it does not lift the WITHHELD

Same mesh the golden case fingerprints (`25802` vertices, `50880` faces), compared against
the station-area route it is meant to cross-check:

| waterplane `z` [m] | `V` closure [m³] | `V` station [m³] | Δ | `A_wp` closure [m²] | `A_wp` station [m²] | Δ |
| --- | --- | --- | --- | --- | --- | --- |
| −2.000 | 23.11 | 23.41 | −1.29 % | 43.09 | 41.52 | +3.78 % |
| 6.157 (design WL) | 8024.43 | 8452.76 | **−5.07 %** | 2010.80 | 2100.33 | **−4.26 %** |

At the waterline the closure route reports `submerged_surface_area = 2986.57 m²` with
`1116` triangles clipped and `orientation = −1` (the shipped winding points *into* the
enclosed column — resolved, not a defect).

A 5 % disagreement between the two routes is far outside any band the golden suite could
assert, so **the `displacement-withheld-open-shell` limitation stands unchanged** and the
baseline was not re-blessed. Recorded as a negative result rather than a criterion.

### Candidate causes (not yet discriminated)

1. **Non-vertical caps.** The closure assumes the sheet's non-waterplane openings are
   vertical planes. Measured off the boundary loop, the stern opening is a *raked* transom
   bottom curve (`x` varies over ≈0.4 m while `z` runs 5.7→10) and the stem seam dips to
   `z ≈ 1.6`; the part of each below the waterline is not a vertical plane. A crude bound
   (few m² × O(0.1) normal tilt) is a few m³ — an order of magnitude too small to explain
   428 m³, so this is a contributor at most.
2. **The tessellated sheet is not the model the station route integrates.** The mesh is a
   bicubic spline of the 2250 published points on a `49 × 133` per-patch grid; the station
   route re-samples the same spline on a refined station grid and clips each section
   exactly. Measured over a z-ladder the deviation is not a constant factor: it is −1.3 % at
   `z = −2 m`, +4.0 % at `z = −1 m`, then a nearly constant ≈ −5 % from `z = 2 m` upwards
   (−5.4 % at 2 m, −5.1 % at the waterline). A sign change below the bilge plus a flat 5 %
   offset above it is not the signature of a single global scale error; it points at the
   whole tessellated surface above the bilge (curvature, deck edge, the two patch seams).
3. **A hole below the waterline.** A missing triangle range (for instance the forefoot
   between the aft/fore patches, `x = 80.13…81.79`) would also show up as an
   under-integral. This is testable and not yet excluded.

### The discriminating test to run next

Per-station section area from the mesh (clip the sheet at `x = const`, integrate the
enclosed polygon at `z <= T`) against `DTMB5415Space._stations` at the same `x`, on a
sweep of stations. That localises the 5 % in `x` and separates cause 2 (smooth, spread over
the hull) from cause 3 (a step at one `x` range). If it turns out to be a hole in the
shipped sheet, the finding belongs to the dataset case, and the closure route then becomes
the independent anchor that makes the 5415 displacement `REAL` instead of `WITHHELD`.

## 4. What the new capability is good for now

Any hull whose surface is a sheet rather than a solid gets a displacement and a waterplane
area with no repair step: box barges, Wigley/Series-60 foils before capping, imported IGES
hull surfaces that are open at the deck. It is also the cheap way to check a *closed* mesh
against itself (`A_waterplane → 0` there, which the test above pins).