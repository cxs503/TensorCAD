# Large-angle stability: GZ from the submerged volume

`tensorcad.stability` answers the small-angle question. It returns `KB`, `BM`, `KM`,
`GM` from the upright waterplane, which is enough to say whether a loading condition
is stable and nothing more. A `GM` cannot produce a `GZ` curve: the deck-edge
immersion and the loss of waterplane that decide whether a hull survives 30-60 degrees
of heel are second-order in the small-angle expansion, i.e. absent from it.

`tensorcad.large_angle` (new in 0.21.0) computes the righting arm directly:

```
GZ(phi) = KN(phi) - KG sin(phi)
```

with `KN` the transverse lever of the centre of buoyancy about the keel, at
**constant displacement** and level trim. For every heel angle the immersed region is
the hull clipped by a **tilted** waterplane whose sinkage is solved so the displaced
volume is unchanged — exactly the geometry `tensorcad.hydrostatics` declines to handle
("no face straddles the waterline").

## Why the answer is exact rather than grid-convergent

A triangle mesh cut by a plane of constant `z` gives a cross-section whose boundary
vertices move **linearly** in `z`. Its area and first moments are therefore
**quadratic in `z`** between two consecutive vertex heights. Three samples per band, at
`1/4, 1/2, 3/4` of the band (chosen so no sample ever lands on a vertex height, where
the section is degenerate), determine that quadratic exactly, and the band integral is
then evaluated in closed form — including a *partial* band, which is what a bisection on
the waterline needs.

Two practical consequences:

* the sinkage solve is a bisection on a monotone function, not a re-mesh, so a heel
  angle costs no more than a waterline query;
* the numbers are floating-point exact, so the residual against a closed-form reference
  measures the *reference*, not the quadrature. That is what makes the box case below
  worth running: there is nothing else left in the error budget.

Cross-section boundaries are traversed in the direction `z_hat x n` for the face's
outward normal `n`, which is the way a positive-area loop is traversed. Shoelace sums
then come out positive without ever linking segments into loops, so a section that
splits into several closed curves (or develops a hole) needs no special case.

## Conventions

* `+x` forward, `+y` to starboard, `+z` up, upright waterline at `z = 0`. `heel_deg > 0`
  heels **to starboard**: the body is rotated about `x` by `-heel`, so the water stays
  level.
* `kg` is the height of the centre of gravity **above the keel**, and it is a required
  keyword. A default would make a silent `KG = 0` look like a valid loading condition.
* Fixed trim. A trim excursion changes `GZ` by an amount this module does not model.
* `displacement_volume` (hold displacement, let sinkage follow) **or** `waterline_z`
  (hold the waterplane, let displacement follow, as a fixed-carriage model test does).
  With neither, the upright displacement at `z = 0` is the target.

## The freeboard problem, and `extend_freeboard`

Every generator in this package caps the surface **at the design waterline**. That is
right for a wetted-surface calculation and fatal here: such a body displaces its entire
volume upright, so no tilted waterplane can hold that displacement and `GZ` does not
exist. `heeled_hydrostatics` refuses the request with an explanation instead of
returning a number:

```
LargeAngleError: the mesh is capped at its own top surface: it displaces its whole
volume (7.446554097149208 m^3), so no tilted waterplane can hold that displacement and
a heeled righting arm is undefined. This is a hull *without freeboard* -- wrap the mesh
in extend_freeboard() ...
```

`extend_freeboard(mesh, height)` deletes the flat deck faces, extrudes the deck's
boundary loop vertically, and re-caps with the deleted deck's own triangulation lifted.
The result is watertight and **waterline neutral**: no vertex below the deck moves, so
the immersed volume and centre of buoyancy it reports reproduce
`compute_hydrostatics` on the capped mesh to 2.4e-16 relative (tested). `GM` therefore
does not change either.

One boundary to keep in mind: `compute_hydrostatics` and `initial_stability` integrate a
surface *already clipped at the waterline* and treat a closed mesh as the whole wetted
solid, so on a freeboarded mesh they report the entire solid rather than the
displacement. Feed those two the capped mesh and this module the extended one.

## Verification

Run `python examples/large_angle_stability/run_example.py --output docs/large_angle_stability`
(PYTHONPATH=.). Four references, none of which is another implementation of the same
integral:

| reference | what it grades | measured |
|---|---|---|
| wall-sided closed form `sin(phi)(GM + BM tan^2(phi)/2)` for a box barge | the section integral, on a surface with no curvature | worst rel. error **2.88e-15** over 0-30 deg |
| closed form `int_0^Phi GZ dphi = GM(1-cos Phi) + (BM/2)(cos Phi + sec Phi - 2)` | the trapezoid area under a curve | 6.04e-04 relative at 2.5 deg steps, `O(h^2)` on refinement |
| `GM` from `stability.initial_stability` on the capped mesh | the slope `dGZ/dphi(0)` measured on the freeboarded mesh | **7.74e-07** relative |
| rotational symmetry of a circular-section barge | the integrator on a curved surface: `GZ = (R - KG) sin(phi)` at every angle | 5.08e-04 at 32 facets, 4.03e-06 at 128 (facet error, shrinks ~126x) |

### Box barge (`L = 6, B = 2, T = 0.8`, freeboard 0.6 m, `KG = 0.5 m`)

| heel [deg] | GZ measured [m] | wall-sided [m] | rel. error | waterline [m] |
|---:|---:|---:|---:|---:|
| 5 | 0.027738300466 | 0.027738300466 | 3.25e-15 | -0.0000 |
| 10 | 0.056113366873 | 0.056113366873 | 1.36e-15 | -0.0000 |
| 20 | 0.117745730766 | 0.117745730766 | 4.71e-16 | -0.0000 |
| 30 | 0.193055555556 | 0.193055555556 | 2.30e-15 | +0.0000 |

(For a box of constant displacement, the mean immersed depth never changes, so the
waterline offset stays at zero — a free consistency check on the sinkage solve. The
measured arms are exact at 30 deg, just inside the deck-edge limit
`atan(2*0.6/2) = 30.96 deg`; beyond it the reference over-predicts, which
`tests/test_large_angle.py` asserts as a **negative control**: a shallow-freeboard box
at 25 deg returns 0.041 m against a formula value of 0.153 m.)

### Shipped hull family with 0.6 m freeboard (`L = 10, B = 2, T = 0.8`, `KG = 0.5 m`)

`Cb = 0.4654`, `Cp = 0.7031`, `GM = 0.421623 m`, `KM = 0.921623 m`, displacement
7.446554 m^3 (the capped mesh, unchanged by the extension).

| heel [deg] | GZ [m] | KN [m] | waterline [m] | displacement [m^3] |
|---:|---:|---:|---:|---:|
| 0 | 0.000000 | 0.000000 | +0.0000 | 7.446554097 |
| 10 | +0.068976 | 0.155800 | -0.0033 | 7.446554097 |
| 20 | +0.132034 | 0.303044 | -0.0130 | 7.446554097 |
| 30 | +0.193606 | 0.443606 | -0.0293 | 7.446554097 |
| 40 | +0.253147 | 0.574541 | -0.0499 | 7.446554097 |

Peak arm 0.2531 m at 40 deg (the sector ends there, so the peak is not bracketed: the
code says so rather than implying the curve turns over). Displacement is held to
9.78e-15 relative and `GZ` is odd in heel to 2.78e-17 m. The waterline **sinks** as heel
grows — the section lost on the emerging side has to be recovered by pushing the hull
down — which is the physical content of the constant-displacement solve.

A second variant with only 0.1 m of freeboard (`KG = 0.6 m`) peaks at 0.0487 m near
15 deg and vanishes at **48.71 deg** (linear interpolation between the 5-deg
tabulation), which is the behaviour a `GM`-only method cannot represent at all.

## What this is not

No free-surface correction (tank contents are solid here), no wind or launching heels,
no damage or flooded-compartment cases, no IMO/IS-Code compliance statement. `GZ` is a
geometric function of one hull, one displacement and one `KG`; whether a ship satisfies
a criterion is a regulatory question this package does not answer.