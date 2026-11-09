# Boolean operations on fields: union, intersection, difference

Mission m7565 (`TensorCAD`). Code: `tensorcad/csg.py`. Tests: `tests/test_csg.py`
(25 cases). Example + evidence: `examples/csg_boolean/run_example.py`,
`docs/csg_boolean/results.json`, `docs/csg_boolean/plot.png`.

## 1. Case

Every real model is more than one surface. A submarine is a hull plus a sail plus four
fins; a tank is an opening cut out of a compartment; a bulb is the part of a body above
the design waterline. Until 0.24.0 the package could say that only by hand-writing a
class: `SuboffFullField` *is* a hard-coded `min` of hull, sail and fins. "No boolean or
sweep operations" has stood on the capability ledger since the 0.4.0 audit
(`backlog_tensorcad_section.md`), and it is the first item TensorShipDesign's geometry
base needs after a hull loft.

What a Boolean has to deliver, and what this node measured:

| requirement | how it is checked here |
|---|---|
| the composed zero set is exactly the Boolean of the operand solids | containment at 2 000 random points against the operand booleans, exact (`==`), no tolerance |
| the solid's volume is the Boolean's volume | six closed forms from elementary solid geometry, voxel integration at four resolutions |
| it composes with what already exists | `union(hull, appendages)` reproduces `SuboffFullField` **bit for bit**, same SDF grid digest |
| the exported mesh is still a closed solid | marching tetrahedra + the standard `inspect_geometry` gate, and the topology change a through-cut causes is reported, not hidden |

## 2. What was done

Three operators, in the standard `min`/`max` form, returning `GeometryField`s:

```
union(a, b, ...)         min(a, b, ...)            points inside any operand
intersection(a, b, ...)  max(a, b, ...)            points inside every operand
difference(a, b, ...)    max(a, -b, ...)           a with the rest removed
```

* **Sign is exact, magnitude is not.** `min(a,b) < 0` is `a < 0 or b < 0` for any real
  numbers, so the zero set, `contains`, mesh extraction and voxel occupancy are exact
  statements about the composed solid. `|min(a,b)|` is the true distance only where the
  nearer surface belongs to the operand that attains the minimum; inside the other body
  it reports the distance to *that* body's surface. An exact CSG distance needs a
  nearest-point query, so the field declares `quantity = "implicit_function"` — the same
  honesty `EllipsoidField` and `SuboffFullField` already apply. The error is in the safe
  direction: `|value| >= true distance` everywhere, which is what keeps a distance-band
  test, a voxel occupancy array and marching tetrahedra conservative.
* **Gradient** is the active operand's gradient (negated for the subtracted operand).
  Exact away from the seam; on the seam one-sided derivatives are all a field can offer,
  and that is stated rather than smoothed over.
* **Bounds** are per operator: the box around both for a union (`None` if either operand
  is unbounded, e.g. a half-space), the box overlap for an intersection — and an
  intersection of non-overlapping boxes raises `FieldError` instead of sampling an empty
  grid that would silently contain no surface — and the left operand's box for a
  difference.
* **State** is the operands' state, addressed through a `left__` / `right__` prefix, so
  `field.with_state(left__radius=0.6)` edits one leaf and returns a new frozen
  composition. Unprefixed names are an error: ambiguous by design.
* **NumPy path** composes `minimum`/`maximum` over the operands' vectorised forms, which
  select values rather than recompute them, so it is bit-identical to the Python path
  (asserted, not assumed).
* N-ary forms fold left, so `union(a, b, c) == union(union(a, b), c)`.

Nothing in the module knows what a hull is: it is arithmetic on the field protocol, and it
therefore accepts a learned `.ngo` field exactly as it accepts a sphere.

## 3. Results

### 3.1 Volumes against closed forms

Measured with `to_voxel_grid(...).volume()` (cell partial fill, a second-order rule) — an
integration route independent of the field formula being tested. Printed by
`examples/csg_boolean/run_example.py`, stored in `docs/csg_boolean/results.json`.

| composed solid | closed form | 48 cells | 64 cells | 96 cells | 128 cells |
|---|---|---|---|---|---|
| cap: ball R=1 minus half-space z≥0.5 | `πh²(3R−h)/3` = 0.654498 m³ | +0.1264 % | +0.0531 % | +0.0227 % | **+0.0032 %** |
| zone: ball R=1, \|z\|≤0.4 | `2π(R²t−t³/3)` = 2.379233 m³ | +0.1021 % | +0.0407 % | +0.0208 % | **+0.0134 %** |
| bicylinder: two cylinders r=0.6, perpendicular | `16r³/3` = 1.152000 m³ | +0.1703 % | +0.0076 % | +0.0296 % | **−0.0068 %** |
| lens: two balls R=0.75, centres 1.0 apart | `π(4R+d)(2R−d)²/12` = 0.261799 m³ | +0.0398 % | +0.0295 % | +0.0185 % | **+0.0106 %** |
| shell: ball R=1 minus concentric R=0.6 | `4π/3(Ro³−Ri³)` = 3.284012 m³ | +0.0601 % | +0.0095 % | +0.0102 % | **−0.0001 %** |
| ring: ball R=0.6 minus coaxial cylinder r=0.3 | `4π/3(R²−r²)^{3/2}` = 0.587671 m³ | +0.0280 % | +0.0396 % | −0.0056 % | **+0.0170 %** |

Worst error at 128 cells: **0.017 %**. The error is not strictly monotone in the
resolution case by case (voxel phase noise: 0.396 % at 64 then −0.056 % at 96 for the
ring), so the claim the example machine-checks is the trend — every case is closer at 128
cells than at 48 — plus the absolute bound of 1 %.

The `ring` reference is worth a remark because getting it wrong is easy and the code
would have taken the blame. A ball with a coaxial through-hole is **not** "sphere minus a
cylinder of length `2√(R²−r²)`": the hole is a cylinder capped by two spherical caps,
and the naive reading is 7.7 % low. Derived properly, the solid left behind is
`4π/3(R²−r²)^{3/2}`, exactly the ball the hole would carve out of. The first draft of
`tests/test_csg.py` used the naive form and reported a 3.8 % "error" against correct
code; the table above is the corrected version.

### 3.2 Mesh export and the topology a cut can change

`difference(ball R=0.6, coaxial cylinder r=0.3)` extracted by marching tetrahedra:
watertight, edge-manifold, vertex-manifold, no degenerate facets, outward orientation,
signed volume 0.5877 m³ vs the closed form 0.587671 m³ (−0.133 % at 56³, −0.025 % at 128³).

A through-hole makes the solid a **torus**, so `inspect_geometry`'s genus-0 check returns
`V−E+F = 0` and fails — correctly. That is a real semantic every downstream consumer has
to know: a Boolean cut can change the topology of the body, and the validation gate says
so instead of reporting genus 0 to keep a green tick. The same operator with a blind
pocket (cylinder that does not pierce) stays genus 0 and passes every check; both cases
are pinned in `tests/test_csg.py::TestExport`.

### 3.3 The composition the package already hard-coded

`union(SuboffFullField().hull(), SuboffFullField().appendages())` against the existing
`SuboffFullField`, at 400 random points spanning the 4.356 m model and on a 64³ voxel
grid: **identical field values, identical containment, identical volume, identical SDF
SHA-256 digest**. `min` selects an operand's value rather than recomputing it, so
equality is exact rather than "within 1 ULP". The general operator therefore did not
quietly change the SUBOFF geometry the rest of the ecosystem benchmarks against.

## 4. Figures

`docs/csg_boolean/plot.png`: (left) centreplane field of the composed SUBOFF with the
zero level set drawn — the surface a solver would receive; (middle) the field of the
drilled ball with the hole visible and its genus-1 topology; (right) measured volume
error against the six closed forms plus the extracted-mesh error, versus resolution.

## 5. Notes, limits

* **No exact distance inside a composition.** Stated in §2 and in the module docstring;
  measured (`tests/test_csg.py::test_composed_value_never_understates_the_distance`): the
  inequality against a brute-force nearest-surface search holds everywhere, and equality
  holds exactly for points where one operand dominates (0 violations over 200 probes).
* **No smooth blend.** A `smooth_min` with a fillet radius would help
  hull-appendage junctions, but it changes the solid (the union is no longer an exact
  set), so it is a different claim and is not made here.
* **No sweep/extrude.** The ledger item "no boolean *or sweep* operations" is half
  closed: solids can now be combined and cut, but a profile cannot yet be swept into a
  body.
* **Intersection bounds are a box overlap, not the true extent**, and a difference
  reports the left operand's box. Both are supersets — correct for sampling, wasteful for
  a slender intersection of two large bodies (a 2.4:1 waste in the bicylinder case above).
* **The genus check is not "CSG-friendly".** A cut that changes topology trips
  `euler_characteristic_genus_zero`; the check is right, but a consumer that treats any
  failing `inspect_geometry` row as "bad geometry" will misread a legitimate through-hole.
  Left as a documented behaviour rather than special-cased, because the topology of the
  body is the caller's business, not the Boolean's.

## 6. Error determination

Key quantities and their error against their reference, each item explicit:

| quantity | reference | measured | error | < 3 %? |
|---|---|---|---|---|
| cap volume | closed form, 0.654498 m³ | 0.654519 m³ (128³) | +0.0032 % | **yes** |
| zone volume | closed form, 2.379233 m³ | 2.379551 m³ | +0.0134 % | **yes** |
| bicylinder volume | closed form, 1.152000 m³ | 1.151921 m³ | −0.0068 % | **yes** |
| lens volume | closed form, 0.261799 m³ | 0.261827 m³ | +0.0106 % | **yes** |
| shell volume | closed form, 3.284012 m³ | 3.284009 m³ | −0.0001 % | **yes** |
| ring volume (voxel) | closed form, 0.587671 m³ | 0.587771 m³ | +0.0170 % | **yes** |
| ring volume (mesh, 128³) | same | 0.587526 m³ | −0.0247 % | **yes** |
| SUBOFF composed vs `SuboffFullField` | same package, 64³ grid | identical | 0 (bitwise) | **yes** |

Every case is below 0.02 % at 128 cells; the 3 % band is not
binding here. No external literature value is compared against in this node — the
references are closed forms derived in the example file, which is why they count as an
independent route (integration of cells) rather than as a self-consistency check
(field formula re-evaluated).