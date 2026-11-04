# Curvature and surface fairness

Every other check in this package is topological or metric: closed, outward, right
volume, no slivers. None of them says whether the surface is **fair** -- whether its
curvature is finite, continuous, and free of accidental knuckles. `tensorcad/hull.py`
says so in its own module docstring ("nothing here checks fairness"), `docs/
mesh_conditioning.md` ends up reporting slivers along sharp edges without ever
measuring the curvature that produced them, and the backlog lists the fairness check
as the last open P3 item. `tensorcad.curvature` is that item.

Everything on this page is measured by `python examples/fairness/run_example.py`,
which writes `docs/fairness/{results.json,report.md,plot.png}` and exits non-zero if
any claim stops holding. `tests/test_curvature.py` pins the same numbers.

## Two estimates, on purpose

The curvature of the zero set of a differentiable field is standard differential
geometry: with `n = grad f / |grad f|` and `P = I - n n^T`, the shape operator is

    M = P (Hess f) P / |grad f|

and its two non-zero eigenvalues in the tangent plane are the principal curvatures.
`P` kills the normal direction, so the formula is correct for a true signed distance
*and* for an implicit function that is not one -- which is what makes it usable on the
ellipsoid and on the learned `.ngo` fields, not just on the exact primitives. The
Hessian is a central difference of `field.gradient`, and the seven gradients per probe
(the point plus three axis pairs) go out in a single batched call: fields evaluate
point lists vectorised, and measuring that was worth it (a per-probe call is roughly
an order of magnitude slower on the AFF-8 min-union).

`vertex_mean_curvature` gives the mesh estimate from the cotangent Laplacian. A
consumer of a mesh can recompute that; nobody can recompute the field once it is gone.
The two disagree, and the disagreement is the discretisation error:

| estimate | how accurate | why it exists |
|---|---|---|
| field (`principal_curvatures`) | exact to round-off on the primitives (`1e-15` relative) | the reference; what a designer means by "the shape" |
| mesh, pointwise (`vertex_mean_curvature`) | 49 % median on a raw marching-tet sphere, 1.4 % after conditioning | what a solver actually consumes |
| mesh, area-weighted mean | `0.05 %` on the raw mesh, no conditioning needed | the integral `(1/A) int H dA`, which the cotangent operator conserves exactly |

That last row is the useful surprise. The cotangent operator is a *discrete
divergence*: on a closed mesh the area-weighted sum telescopes and reproduces the
exact integral of mean curvature, so a sliver field cannot move it. Pointwise values
are a different story -- slivers carry almost no area and almost all of the pointwise
error, which is exactly why the raw mesh's median is tens of percent off while its
integral is right.

## What the audit reports

`assess_fairness(mesh, field=...)` returns a `FairnessReport` with three groups:

* **dihedral angles** over interior edges (mean, p50, p95, max, min) and the count and
  length fraction of edges above a crease threshold (15 deg by default, the usual
  ship-form convention: a knuckle or transom edge is far above it, an under-resolved
  smooth patch is far below);
* **mesh curvature** from `vertex_mean_curvature`, summarised as median, p95 and the
  area-weighted mean;
* **curvature jumps** across interior edges -- only with a field, because the jump is
  then a property of the geometry rather than of the triangulation: the mean curvature
  is evaluated at the two face centroids either side of an edge and differenced.

Nothing here is a gate. A hard threshold on creases would reject the DARPA AFF-8
model, whose NACA fins genuinely have 0-degree trailing edges.

## Measured

**Analytic field, exact.** Sphere `1/R, 1/R`, cylinder `1/R, 0`, plane `0, 0`,
ellipsoid at all three poles (`a/b^2, a/c^2` etc.): worst relative error `1.4e-15`,
the plane exactly zero. The estimator's own reported step error stays at `~1e-7 1/m`,
seven orders below the curvatures it measures.

**AFF-8 appendages (80x12x16, 22,060 triangles).** p50 dihedral `0.049 deg`, p95
`11.9 deg`, max `163.5 deg`, min `0 deg`, 1,171 crease edges over 15 deg = `3.04 %`
of the interior edge length. Ten Taubin sweeps do not remove them -- they should not --
but they cut the number that survive at the same threshold, which is the sliver
population.

**The hull family is not G1 at its ends.** `b(x) = B/2 * sin(pi x/L)^0.7` gives
`db/dx ~ x^-0.3`; the fit over `x/L = 1e-2 .. 1e-7` returns `-0.3000`. The waterline
has an infinite slope at the stem and the stern, so the surface is C0 there, its
curvature diverges, and any curvature-driven mesher or fairing loop will fight it.
The keel line has the same problem for the same reason: the profile is
`sqrt(z/T)`, so `dy/dz -> infinity` at `z = 0`. Both are properties of the published
family, not of a mesher, and the honest fix is a finite bilge radius and a fuller end
in the parameterisation -- not more grid.

**A crease threshold on a coarse mesh measures the mesh.** On the hull, p95 dihedral
falls `26.80 -> 3.59 deg` and the crease length fraction `7.55 % -> 1.44 %` between
`41x13` and `161x41`, while the crease *count* rises `232 -> 628` (more edges survive
because refinement puts more of them in the end regions). Judge fairness on the p95
and the length fraction, or on the field; a raw count is not resolution-independent.

## What this does not do

No curvature-driven remeshing, no fairing/optimisation loop (the audit says where, it
does not fix), no principal *directions* exported (the eigenvalues are computed from a
2x2 restriction, so the eigenvectors exist but are not part of the API yet), and no
umbilic/singularity classification. The finite-difference Hessian needs a field whose
gradient is defined at and around the probe point; on a mesh with no field behind it
the report degrades to the dihedral and cotangent parts, which is the intended
fallback and is labelled as such in the JSON.

## Measured tolerances, not aspirations (m3313)

Two tests in this suite used to demand more than the estimators can give, and the
numbers behind their replacement are measured by
`runs/2026-11-04-curvature-adjudication/measure.py` (JSON + figure alongside).

**The unrefined probe is second order, and that is the whole of its error.**
`|kappa_max - 1/R| / (0.5 h^2) = 1.000000` (worst deviation 4.8e-5) over fourteen
steps from `8e-3` to `9.8e-7 m`, and the observed order over the halving ladder is
2.00 at every rung. At the default step (`1e-4` of the bbox diagonal, `3.464e-4 m`
on the unit sphere) the error is `6.0e-8`, which is exactly `0.5 h^2` -- **not** a
float32 bandwidth floor: a float32 field (gradients and point rounded to binary32)
sits on a flat `2.5e-8` floor below `h = 2e-4` and cannot improve, while the double
precision path keeps falling to `4.8e-13` at `h = 9.8e-7`. The test therefore pins
`0.5 * step**2` for `refine=False` and `1e-12` only for the Richardson-extrapolated
`refine=True` path, where the `O(h^2)` term is the thing being removed.

**A sphere fixture must be a legal mesh.** A UV sphere built from `nu + 1` rings has
`nv` coincident vertices at each pole, hence `2 nv` zero-area faces; `assess_fairness`
refuses such a mesh (a face with no normal has no dihedral), and the cotangent
estimator reports `H = 0` at the collapsed vertices -- a silent "this is flat". The
suite's fixture is now a pole-fan sphere (single pole vertex, no zero-area face, all
validator checks pass, area-weighted mean `1.000004` at `16x32`) and the collapsed
mesh is kept as the negative control that pins both halves: 32 zero-area faces, four
failed validator checks, `CurvatureError` from the fairness report.

**The dihedral threshold is resolution-dependent by construction.** On the pole-fan
sphere p95 is `22.80 -> 11.29 -> 5.63 deg` at `8x16 -> 16x32 -> 32x64`, halving with
the grid (ratios 2.020 and 2.005), so the `16x32` assertion is `< 12 deg` *and*
`max < DEFAULT_CREASE_DEGREES` with a separate refinement test pinning the halving.
Note the consequence: at `8x16` the fan's own dihedral p95 (22.8 deg) is inside the
15 deg crease band, i.e. a coarse enough sphere is audited as creased -- the
threshold says nothing about a mesh this under-resolved, which is why the refinement
trend is asserted rather than a single number.
