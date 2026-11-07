# PC6 hull form and the de Long -> Pevek ice channel as one geometry contract

**Mission** IN12-P1 / PC6 geometry service (TensorCAD). **Artifact** `artifacts/geom_pc6.json`
(schema `tensorcad.pc6.geom/1.0`, digest `cf5accb44183b2a392d0c701735a1458...`).
**Status** all acceptance criteria met: watertight closed manifold, beam sweep 24 -> 28 m
regenerating cleanly, 14 pytest cases passing (0 failed), both
negative controls caught.

## 1. Case: what the geometry has to be good for

The ice legs of the campaign need one geometry that four different solvers can share: a DEM
hull-contact surface, a FEM ice belt, and an LBM/FVM flow domain. Getting four slightly
different hulls is the classic way to make four slightly different answers. This node
therefore publishes **one** contract holding both representations of the same body plus its
parameter table, and re-verifies it every time it is loaded.

The hull is a simplified PC6 (Polar Code category 6) form at **L = 200 m, B = 24 m, T = 11 m**:
a forward-raked stem, flared topsides and a wide flat transom, i.e. the waterline body of a
full-form ice-class bulk carrier. The body is closed at the design waterline, so the "deck" of
the closed solid is the waterplane; that is the body a resistance or seakeeping solver meshes,
and it is stated rather than implied because the contract does **not** describe the topsides
above the DWL.

The channel is the de Long (Bennett Island) -> Pevek leg of the East Siberian Sea route: eight
waypoints on a great-circle polyline, a 250 m wide corridor around it and the closed prism of
that corridor down to the hull's design draft.

### Parameter table (from the contract)

| parameter | value | unit | description |
|---|---|---|---|
| aft_exponent | 1.1 | 1 | run fullness exponent |
| beam | 24 | m | moulded breadth at the design waterline B |
| bilge_height | 0.08 | 1 | draft fraction over which the bilge turns |
| bilge_ratio | 0.95 | 1 | half-breadth at the keel / DWL half-breadth |
| draft | 11 | m | keel to design waterline T |
| flare | 0.15 | 1 | cubic topside flare coefficient |
| fore_exponent | 0.22 | 1 | entrance fullness exponent |
| length | 200 | m | length between perpendiculars L |
| longitudinal_stations | 97 | count | requested station count (breaks inserted on top) |
| midship_fraction | 0.62 | 1 | position of the maximum section as a fraction of L |
| stem_rake | 8 | m | horizontal setback of the stem at the keel (forward rake) |
| transom_ratio | 0.88 | 1 | transom half-breadth / maximum half-breadth |
| vertical_levels | 19 | count | requested level count (breaks inserted on top) |
| water_density | 1,025 | kg/m^3 | density turning displacement volume into mass |

### Principal quantities of the design point

| quantity | sampled mesh | independent reference | relative error |
|---|---|---|---|
| displacement volume  ∇ | 41,544.876 | 41,629.884 | -0.204 % |
| mass displacement  Δ = ρ∇ | 42,583.5 | — | at rho = 1025 kg/m^3 |
| waterplane area  A_wp | 4,276.397 | 4,284.019 | -0.178 % |
| block coefficient  C_b | 0.78683 | 0.78844 | -0.204 % |
| maximum section area (C_m x B x T) | 239.368 | 237.868 | +0.630 % |
| prismatic coefficient  C_p | 0.86781 | 0.87506 | derived |
| LCB from AP | 94.407 | 94.586 | -0.189 % |
| LCF from AP | 96.276 | 96.458 | -0.189 % |
| VCB from DWL | -5.334 | -5.336 | -0.041 % |
| wetted surface  S | 8,143.4 | no closed form | integration of the sampled shell |

`C_b = 0.7868`, `C_m = 0.9010`, `C_p = 0.8751`,
`C_w = 0.8925`, LCB at `94.59 m` from the AP
(`47.3 % L`). The set is coherent for a full-form ice-class
bulker (C_w > C_b > C_p C_m holds by construction, VCB near half draft) but it is **not** a
model-tested design: no resistance, stability or ice-load validation is claimed here.

![hull form](pc6_hull_form.png)

*Figure 1 -- the analytic form the mesh is lofted from: profile with the raked stem and flat
transom, body plan at four stations, waterlines from keel to DWL, and the section-area curve
whose mean over the length is C_b.*

## 2. Method

**Parametric surface.** Stations run on `s in [0, 1]`, heights on `u in [0, 1]`:

```
z(u)    = -T (1 - u)
x(s, u) = s * x_stem(z),     x_stem(z) = L + stem_rake * z / T      (forward-raked stem)
y(s, u) = +- 0.5 B * g(s) * sec(u)
sec(u)  = base(u) (1 + flare u^3) / (1 + flare),  base = f_b + (1 - f_b) sqrt(u / u_b)
g(s)    = t0 + (1 - t0) (s/s_m)^e_aft   (s <= s_m);   (1 - (s-s_m)/(1-s_m))^e_fore  (s > s_m)
```

`g(1) = 0` with `sec(u) > 0`, so the starboard and port sheets meet exactly on the raked stem
curve. The mesh therefore closes **by construction**: the two sheets, the flat of bottom, the
waterplane deck and the flat transom share their boundary vertices, with no welding and no
tolerance parameter deciding the topology. Station and level lists insert `s_m` and `u_b`
exactly, so the sampled solid reproduces L, B and T to machine precision rather than
undershooting them by interpolation.

**Implicit field.** The same body is published as `f(x, y, z) = -min(-z, z+T, x, x_stem(z)-x,
b(s,u)-|y|)`: negative inside, zero on the shell, exact sign, magnitude a lower bound on the
distance to the nearest constraint (not a signed-distance function -- using it as one would
bias a level-set solver).

**Independent references.** Three things are never checked against themselves: (i) the
displacement, waterplane area, LCB/LCF/VCB and C_m come from closed-form power integrals of the
analytic surface, a code path the mesh integrator does not share; (ii) the implicit field is
compared with the sign of the *solid angle* (winding number) subtended by the STL at
deterministic random points; (iii) the channel prism volume is compared with plan length x
width x depth, which for mitred offsets must hold exactly.

**Self-intersection.** A watertight mesh can still fold through itself, so a geometric test
runs on every generated hull: edge-versus-face piercing (the standard three-edge formulation,
which is complete for two convex triangles because any intersection segment begins or ends on
a boundary edge of one of them) plus a coplanar containment branch, pruned with a uniform hash
grid over triangle bounding boxes. Faces sharing a vertex are skipped -- a manifold mesh meets
there by construction.

**Channel.** Waypoints are projected onto a local azimuthal-equidistant plane centred on the
first waypoint, so along-track lengths are true great-circle distances (a plate carree at
70-77 N would be off by a factor 1/cos(74 deg) ~ 3.6 in longitude). Walls are mitred
perpendicular offsets, mitre length capped at twice the half-width; the prism is built from the
four corner columns per station, so it is closed for any simple plan shape.

**Toolchain.** Pure standard library for the geometry (`tensorcad.pc6`, `tensorcad.icechannel`,
`tensorcad.intersect`, `tensorcad.pc6_contract`), pytest for the acceptance suite, matplotlib
only for the figures. One command regenerates everything: `python tools/pc6_contract.py`.

## 3. Results

### 3.1 Topology of the published solid

| check | value | note |
|---|---|---|
| vertices / faces / edges | 3900 / 7796 / 11694 | closed quad loft, 2 sheets + 3 caps |
| watertight (every edge in exactly 2 faces) | True | criterion |
| edge manifold | True | criterion |
| boundary edges | 0 | criterion |
| Euler characteristic | 2 | genus-0 shell = 2 |
| connected shells | 1 | one solid |
| zero-area faces | 0 | no degenerate triangles |
| vertex-link defects | 0 | no pinched vertices |
| signed volume (divergence) | 41,544.876 m^3 | positive = outward orientation |
| surface area | 12,419.8 m^2 |  |

### 3.2 Discretisation convergence (independent of any tuning)

| grid (stations x levels) | triangles | volume error | waterplane error | LCB error |
|---|---|---|---|---|
| 41 x 11 | 1988 | -0.545 % | -0.516 % | -0.546 % |
| 57 x 13 | 3216 | -0.382 % | -0.343 % | -0.364 % |
| 73 x 15 | 4700 | -0.288 % | -0.252 % | -0.268 % |
| 97 x 19 | 7796 | -0.204 % | -0.178 % | -0.189 % |
| 129 x 25 | 13464 | -0.144 % | -0.125 % | -0.134 % |
| 161 x 31 | 20668 | -0.110 % | -0.096 % | -0.102 % |

The shipped design point is 97 x 19 (`7796` triangles): its displacement error is
`-0.204 %`, 2.4x inside the 0.5 % gate. The error is a one-sided
inscribed-polyhedron bias that falls monotonically with refinement, which is the signature of
discretisation and not of a modelling error.

### 3.3 Beam sweep 24 -> 28 m (automatic regeneration)

| B [m] | triangles | ∇ mesh [m^3] | ∇ analytic [m^3] | volume error | C_b | A_wp [m^2] | valid |
|---|---|---|---|---|---|---|---|
| 24.0 | 7796 | 41,544.9 | 41,629.9 | -0.204 % | 0.78683 | 4,276.4 | yes |
| 26.0 | 7796 | 45,006.9 | 45,099.0 | -0.204 % | 0.78683 | 4,632.8 | yes |
| 28.0 | 7796 | 48,469.0 | 48,568.2 | -0.204 % | 0.78683 | 4,989.1 | yes |

Each beam goes through the identical validated path (`beam_sweep`), so a beam that broke
watertightness or the analytic displacement check would raise instead of being recorded as a
pass. The sweep meshes are regenerated on demand rather than committed: the numbers and flags
above are the evidence, and one call reproduces any of them.

### 3.4 Two representations, one solid

The implicit field and the STL agree at **24/24** deterministic sample points
(13 of them inside the mesh by winding number); `passed = True`.
The contract also carries a coarse 33 x 9 x 7 sample grid so a solver can
initialise without reimplementing the formula; the analytic definition travels with it.

### 3.5 Negative controls

| input | layer | outcome |
|---|---|---|
| draft = 0 m | parameter validation | rejected: invalid PC6 parameters: draft must be a finite positive number, got 0. |
| draft = -4 m | parameter validation | rejected: invalid PC6 parameters: draft must be a finite positive number, got -4 |
| stem_rake = 1.25 L | geometric self-intersection detector | 64 intersecting face pairs |
| stem_rake = 1.25 L | parameter validation | rejected: invalid PC6 parameters: stem_rake=250.0 m exceeds 35 % of length=200.0 |
| clean design point | geometric self-intersection detector | 0 intersecting face pairs |

The second row of the geometric case is the important one: the folded loft is rejected by the
*geometry*, not merely by the parameter domain check, which is what makes the detector evidence
rather than decoration.

### 3.6 Ice channel

| quantity | value |
|---|---|
| route | de Long (Bennett Island) -> Pevek, 8 waypoints |
| great-circle length | 1061.8 km (573 nm) |
| projected plan length | 1062.0 km (error +0.023 %) |
| corridor width | 250 m (measured min 250.00 m, max 254.60 m at mitred corners) |
| corridor depth | 11 m below the ice surface |
| prism solid | 60 faces, watertight = True, chi = 2 |
| prism volume | 2.9206e+09 m^3 vs plan estimate 2.9206e+09 m^3 |
| plan shape simple | True (0 self-intersecting face pairs) |

![ice channel](pc6_ice_channel.png)

*Figure 3 -- the de Long -> Pevek polyline with the corridor walls (orange, dashed) drawn back
in degrees from the local plane.*

![validation](pc6_hull_validation.png)

*Figure 2 -- left: displacement and waterplane error against the closed-form integrals versus
grid refinement, with the 0.5 % gate; centre: the 24/26/28 m beams against their analytic
volumes; right: the implicit field at midship, blue where f < 0, with the analytic section
superimposed (the boundary of the blue region is the section curve).*

## 4. Notes, pitfalls and limits

* **Real pitfall found in this node.** The triangle-intersection test first called two parallel
  wall faces 250 m apart *coplanar*, because it compared plane offsets computed with raw
  face normals: two back-to-back triangles give antiparallel normals and offsets of opposite
  sign, so an offset difference of 250 m shrank to ~0.03 m and the coplanar branch fired. The
  fix canonicalises the normal (largest-magnitude component positive) before comparing offsets.
  Before the fix the channel reported 3 false self-intersections; after it, zero -- and the
  folded-loft negative control still reports 64, so the fix removed a false positive without
  weakening the true-positive path.
* **The body stops at the design waterline.** Above-DWL flare, bulwarks and the real deck are
  outside this contract. A seakeeping leg that needs the above-water body must extend the
  generator, not reinterpret this mesh.
* **Single-valued sections only.** `y` is a function of `(x, z)`, so there is no bulbous bow, no
  reentrant stern, no skeg, no bilge keel, no appendage and no ice knife. The form is a
  hydrodynamic envelope, not a structural model.
* **C0 at the maximum section.** `g` switches from the run branch to the entrance branch at
  `s_m` with a slope discontinuity: fair in position, not in curvature. Displacement and
  waterplane area do not feel it; a resistance calculation on the analytic surface would.
* **C_m is 0.63 % off its closed form.** That residual is the trapezoidal section integrator of
  `tensorcad.hydrostatics` acting on a section with a square-root bilge, not a mesh error -- it
  shrinks with the vertical level count and is reported rather than smoothed away.
* **Corridor mitering.** At a sharp waypoint the miter cap clips the outer corner at twice the
  half-width; the measured wall-to-wall distance therefore reaches
  254.6 m on straight runs at corners while the perpendicular
  width stays 250 m. The route itself is campaign-scale (0.1 deg rounded
  positions): it is not a surveyed track and must not be used for navigation.
* **Applicability.** Right for geometry hand-off, voxel classification, mesh seeding and
  parameter studies. Not right for anything that needs the real ship: no model tests, no class
  rule check, no ice-load calibration.

## 5. Error verdict, quantity by quantity

The campaign rule is that a discrepancy above 3 % against a reference is unacceptable.

| quantity | reference | relative error | how the reference was made | < 3 % ? |
|---|---|---|---|---|
| displacement volume ∇ | 41,629.884 m^3 | -0.204 % | closed-form integral of the same analytic surface | YES |
| waterplane area | 4,284.019 m^2 | -0.178 % | B L ∫g ds | YES |
| LCB | 94.586 m | -0.189 % | first moment of the analytic integral | YES |
| VCB | -5.336 m | -0.041 % | vertical moment of the analytic integral | YES |
| C_b | 0.78844 | -0.204 % | ∇ / (L B T) | YES |
| C_m (max section) | 237.868 m^2 | +0.630 % | B T ∫sec du; mesh side uses the trapezoidal section integrator of tensorcad.hydrostatics | YES |
| channel prism volume | 2,920,607,332 m^3 | +0.000 % | plan length x width x depth | YES |
| channel route length | 1,061.791 km | +0.023 % | sum of great-circle segments vs the projected plane | YES |

All eight quantities are below 3 %; the two that matter for a solver hand-off (displacement and
the implicit-field/mesh agreement) are at -0.204 % and 0 disagreements out of
24 samples. No quantity in this node is claimed against a model test or a
published resistance/stability reference, and none is claimed at all -- there is no such
reference for a synthetic PC6 form, and inventing one would be worse than saying so.

## 6. Reproduction and audit chain

```
python tools/pc6_contract.py     # regenerate artifacts + figures + acceptance block
python -m pytest tests/test_pc6.py -q
```

| item | value |
|---|---|
| contract | `artifacts/geom_pc6.json` (schema `tensorcad.pc6.geom/1.0`) |
| content digest | `cf5accb44183b2a392d0c701735a145834bf6de2bf30fda4fd7085af57969d52` |
| hull STL | `artifacts/pc6_hull.stl`, 389884 bytes, sha256 `a459732a2c5f401e26346fb7c39812da...` |
| channel STL | `artifacts/pc6_channel.stl`, 3084 bytes |
| acceptance suite | 14 passed, 0 failed (`14 passed, 1 warning in 10.17s`) |
| reload verification | `load_contract` rebuilt the hull from the STL and matched the recorded displacement to +0.000000 % |

Consumer entry point:

```python
from tensorcad.pc6_contract import load_contract

geometry = load_contract("artifacts/geom_pc6.json")   # re-checks topology and displacement
geometry.hull                                        # closed manifold Mesh, metres
geometry.field(x, y, z)                              # implicit margin, negative inside
geometry.inside(x, y, z), geometry.channel           # membership and the corridor prism
```
