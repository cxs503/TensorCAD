# Righting arms: GZ(phi) at constant displacement

`tensorcad.stability` ends where this module begins. `GM` is the *slope* of the
righting-arm curve at zero heel, and its docstring was explicit that a curve at finite
heel is a different integration: the displaced volume no longer sits under a horizontal
waterline through the upright hull, and the buoyancy centre must be found where the
**heeled** hull displaces the same volume. `tensorcad.righting_arm` is that integration,
and this page is what it measures.

Everything below is written by `python examples/righting_arm/run_example.py
--output docs/righting_arm`, which emits `results.json`, `report.md` and three figures and
**exits non-zero if one of its ten claims stops holding**. The tables on this page are the
script's own output; no number is transcribed.

## How the immersed part is measured

For the region `z <= h` of a closed surface, the divergence theorem needs no cap:

| quantity | surface integral | cap at `z = h` |
|---|---|---|
| `V` | `oint (z - h) n_z dA` | contributes `0` |
| `int y dV` | `oint (1/2) y^2 n_y dA` | `n_y = 0` |
| `int z dV` | `oint (1/2)(z - h)^2 n_z dA + h V` | `z - h = 0` |

Every integrand is a polynomial of degree <= 2 over a planar triangle, so each
contribution is **exact for the polygonal hull**; what is left is the mesh's own
discretisation. Two consequences worth naming:

* no loop has to be chained, so a waterline lying *exactly* on a vertex ring is not a
  special case (it is one for a section-integration route -- see the caveat below);
* a face lying *in* the waterplane contributes exactly zero, which is why a hull capped at
  its design waterline needs no case of its own.

The equilibrium is solved, not assumed: `V(h)` is strictly increasing in `h` (raising the
plane only adds volume), so bracket-and-bisect cannot fail while the displacement is
inside `(0, V_total)`, and the residual is reported rather than trusted.

The transverse arm is `GZ = y_B - y_G`, both in the world frame of the heeled body, with
`y_G = (z_keel + KG) sin(phi)` and `phi > 0` heeled to starboard. The heel is a rotation
about the longitudinal axis **through the upright waterline**, so a mesh whose `z` origin
is not its waterline behaves like an origin-centred one.

## The barge: exact, and exact only up to a stated angle

For a wall-sided barge, `GZ = sin(phi) (GM + (1/2) BM tan^2(phi))` is **not** an
approximation: it follows from `y_B = BM tan(phi)` and `z_B = KB + (1/2) BM tan^2(phi)` at
constant displacement, and it holds while the waterline stays on the vertical sides.
`wall_sided_validity_deg` returns where that stops (bilge uncovered, or deck edge
immersed -- the smaller of the two). On `L = 10, B = 2, T = 0.8 m` with `0.4 m` freeboard
and `KG = 0.55 m` the mesh reproduces the closed form to **1e-12 m** up to the deck-edge
angle, and past it the two part company by centimetres while the mesh still agrees with a
2-D polygon clip of the same section -- so the deviation is the *formula's* limit, not the
integration's error (both are tests, not text).

## The sphere: `GZ = GM sin(phi)` at every heel

A floating sphere's buoyancy centre stays on the vertical through the sphere's centre, so
`GZ = GM sin(phi)` exactly for all heel angles; the residual is pure meshing error. The
ladder in `docs/righting_arm/results.json` shows it falling monotonically with the
tessellation, and a sphere with `KG` at its centre reports `GZ = 0` to well under a
millimetre at 60 deg. This is the curved-surface benchmark the barge cannot give.

## The hull: no closed form, so the checks are different

The parametric hull (`sin^0.7 sqrt`) has no closed-form `GZ`. Generated with depth `1 m`
and floated at `z = -0.47 m`, it carries `0.47 m` of freeboard, and the three checks are:

1. **small-angle slope vs the waterplane moments** -- `dGZ/dphi` at 1e-2 deg against
   `KB + I_T/nabla - KG` computed by `tensorcad.stability` from the *same* hull by a
   different integration (a plane section, not a volume integral): agreement to 1e-6
   relative, i.e. the two independent routes to `GM` agree;
2. **refinement** -- `GZ(30 deg)` on `21x9 ... 161x41` converges monotonically;
3. **constant displacement** -- the relative volume residual is below 1e-9 at every
   sampled heel, so the curve really is one loading condition.

`docs/righting_arm/report.md` lists the curve itself with the solved waterline and the
buoyancy centre at each heel, the peak arm, and the dynamical stability
`int GZ dphi`. The example also *computes* the deck-edge immersion angle from the solved
waterlines instead of assuming it: past that angle the arm is the arm of a **closed
solid**, and a real ship would be taking water on deck, so it is an upper bound.

## Caveats, stated rather than smoothed away

* **Freeboard is the caller's job.** A mesh capped *on* its operating waterline has no
  volume above the water and the solved waterline has to climb over that cap; the arm then
  describes a capped solid, not a ship. `heeled_displacement` reports the waterline, which
  is what makes this visible -- a test asserts the failure mode rather than hiding it.
  `generate_hull` produces exactly such a capped hull, so the example generates depth
  `1 m` and floats it at `-0.47 m` to give it freeboard.
* **Closed solid, no downflooding.** Above the deck-edge angle the integral keeps buoyancy
  where a real ship would have water on deck. No damage, no flooding, no freeing ports.
* **No free-surface correction, no trim coupling.** The heel is prescribed and the
  waterline *height* is solved; a free-trim solve is a second root and is not here. A slack
  tank's `sum(rho i_x) / displacement` must still be subtracted by the caller.
* **`KG` is an input.** It is a loading condition, never inferred, exactly as in
  `initial_stability`.
* **A pre-existing trap in the waterplane route.** `stability.waterplane_moments` raises
  "does not close into a loop" when the waterline falls exactly on a vertex ring of a
  `generate_hull` mesh (e.g. `nz = 13`, `z = -0.5`). The volume integration here is
  unaffected; the example dodges it by floating at `-0.47 m` and that is a reported
  finding, not a fix.
* **`|heel| <= 90 deg`.** Past that the transverse heel convention stops describing the
  problem and the call is refused rather than answered.

## Public surface

`righting_arm`, `righting_arm_curve`, `heeled_displacement`, `RightingArm`,
`RightingArmCurve`, `HeeledDisplacement`, `RightingArmError`, `wall_sided_righting_arm`,
`wall_sided_validity_deg`, `DEFAULT_HEELS_DEG`. Standard library only.

At the package surface three of those names collide with :mod:`tensorcad.gz`, which
answers the same question with different numerics (see below), so `main`'s names keep
their meaning and this module's are reachable under a route-explicit alias:

| this module | package alias | `tensorcad.gz` keeps |
|---|---|---|
| `righting_arm_curve` | `tensorcad.bisection_righting_arm_curve` | `righting_arm_curve` (Newton waterline) |
| `RightingArmCurve` | `tensorcad.BisectionRightingArmCurve` | `RightingArmCurve` |
| `wall_sided_righting_arm` | `tensorcad.barge_wall_sided_righting_arm` | `wall_sided_righting_arm` (`gm`, `bm` arguments) |

Everything else -- `righting_arm`, `heeled_displacement`, `HeeledDisplacement`,
`RightingArm`, `RightingArmError`, `wall_sided_validity_deg`, `DEFAULT_HEELS_DEG` -- is
exported unchanged. `tests/test_righting_arm.py` asserts both directions: the aliases point
at this module, and the plain names still point at `tensorcad.gz`.

## The second route, and why keeping two is the cure

`tensorcad.gz` (0.21.0) computes the same `GZ(phi)` by rotating the mesh, cutting it at the
waterline and driving that waterline with Newton's method on the projected wetted area;
this module clips the heeled section face by face and bisects the monotone `V(h)`. Both
integrals are exact for a planar triangle, so on a given mesh the only difference is the
waterline *solver* -- which is exactly what makes agreement between them worth measuring:

| mesh | `max |GZ_bisection - GZ_Newton|` |
|---|---|
| barge `10 x 2 x 0.8 m` | 1.2e-13 m |
| hull `41 x 13` at `z = -0.47 m` | 1.6e-13 m |
| hull `81 x 25` at `z = -0.47 m` | 6.9e-14 m |
| sphere `R = 1 m`, resolution 32 | 2.8e-16 m |

`tools/verify_merge_m4038.py` writes that table, the arms themselves, and a **negative
control**: the same barge compared against a `gz` curve whose displacement target was
moved by 10 %, which shifts `GZ(20 deg)` by 3.2 mm. A cross-check that cannot fail is not a
check, so the mismatch is part of the record. The merge itself is documented in
`docs/m4038_pr40_merge_cure.md`; the numbers live in
`artifacts/m4038/merge_verification.json` and its figure.
