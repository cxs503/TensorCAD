# Curvature field and smoothing energy (GAP-CAD-A3-CURVATURE-FAIRNESS)

`docs/fairness.md` answers "is the curvature finite and continuous here?" per probe
point. This page covers the rest of A3: a **curvature field** on a mesh (mean,
Gaussian and principal curvature at every vertex, with the inconsistency of the two
discrete operators counted rather than hidden) and the **smoothing (fairness)
energy** -- how much curvature a shape carries in total, whether that number is
stable, and how much of it is a designed crease rather than a fair surface.

Everything on this page is produced by

    python examples/fairness_energy/run_example.py      # writes docs/fairness_energy/{results.json,report.md,plot.png}
    python tools/crosscheck_curvature_thirdparty.py     # writes artifacts/curvature_crosscheck.json
    python -m pytest tests/test_fairness_energy.py      # 20 pins

The runner exits non-zero if any claim stops holding, so the table in
`docs/fairness_energy/report.md` is regenerated, never transcribed.

## The three densities

| name | integrand | SI | anchor |
|---|---|---|---|
| `bending` | `k1^2 + k2^2 = 4H^2 - 2K` | dimensionless | sphere: `8 pi`, any `R` |
| `willmore` | `H^2` | dimensionless | sphere: `4 pi` (Willmore equality case) |
| `gaussian` | `K` | dimensionless | closed genus-0: `4 pi` charge (Gauss-Bonnet) |

They are dimensionless because `dA` carries `1/m^2`, so a 10 m model and a 300 m ship
that are geometrically similar must return the same three numbers. That invariance is
tested, not assumed: doubling `L, B, T` changes `E_bending` by `0.0` (exactly) while
the wetted area grows by `4.0`.

API:

    from tensorcad.fairness_energy import fairness_energy_report, principal_curvature_field
    report = fairness_energy_report(mesh, field=field, crease_angle_degrees=30.0)
    report.bending.smooth_energy, report.bending.crease_energy, report.identity_gap
    field_values = principal_curvature_field(mesh)   # per-vertex H, K, k1, k2, flags

`field=` switches the density evaluation to the exact field at each face centroid;
without it the densities come from the cotangent mean curvature and the angle-defect
Gaussian curvature. `crease_angle_degrees=` splits the surface into smooth cells and
cells within one ring of a crease or an open boundary, and reports the two sums
separately (`smooth_energy`, `crease_energy`, `crease_area_fraction`).

## Why the split matters (the Wigley result)

The lofted Wigley solid is closed by a flat deck cap and by both flanks meeting along
the keel and stem lines, so it has four real creases and four conical corner points
(each carrying a 180 degree angle defect; the discrete defect sum is `4 pi` exactly).
A pointwise mesh estimate of the curvature *at* a crease grows like `1/h`, so the raw
total energy of that hull would diverge under refinement and look like an unstable
estimator. It is not: the smooth part converges and the crease part diverges, and this
capability reports both. Measured (`results.json`):

* smooth `E_bending` climbs with resolution and, extrapolated to `h = 0` from the three
  finest grids, lands **0.06 %** from the analytic reference;
* the crease bucket grows monotonically (`465` -> `1130` over 21x7 -> 121x37), which is
  the signature of designed knuckles and is the number a ship designer wants.

`crease_area_fraction` falls with refinement (~`0.36` at 121x37), so the excluded band
is `O(h)` wide and cannot be hiding a converged contribution.

## The analytic Wigley reference (why it exists)

No published Wigley bending energy exists to compare against, so it is derived and
integrated here: `tensorcad.wigley.analytic_flank_energy` evaluates the standard
graph formulas

    H = ((1 + hw^2) hxx - 2 hx hw hxw + (1 + hx^2) hww) / (2 J^3)
    K = (hxx hww - hxw^2) / J^4,                 J = sqrt(1 + hx^2 + hw^2)

on the analytic flank `h(x, w) = (B/2)(1 - |2x/L - 1|^n)(1 - (w/T)^m)` with a
tensor-product Gauss-Legendre rule (4 panels x 64 nodes by default). Three independent
checks, all in the runner:

1. **quadrature p-extension** (order vs 2*order): `~1e-15` relative;
2. **an independent quadrature of the same surface** (`analytic_wetted_surface_area`,
   composite panels, hand-written Jacobian): wetted area agrees to `4e-16` relative --
   the area is the shared factor of all three energies, so this validates the volume
   element, not just the code path;
3. **the mesh path** converges to it (above).

`E_gauss` of the smooth flank integrates to `~1e-17` -- the crest and trough regions
of the Wigley cancel exactly to quadrature precision -- and the whole `4 pi` charge
sits at the four conical corners. That is consistent with the mesh: the corner
vertices hold 180 degree defects each.

## Numerical stability

`math.fsum` for every accumulation (exactly rounded, so cell order cannot change the
result -- tested on the accumulator and end to end), three separate buckets, no
clipping. Measured on the Wigley solid: scale invariance `0.0`, cell reordering
`2e-16`, translation by 1.6 km `2e-13` (the hard one: cotangent weights are coordinate
differences, so a large offset cancels badly without an exact sum). `condition`
reports `sum|integrand dA| / |sum integrand dA|` -- `1` for the positive densities and
as large as `sum|K| dA / 4 pi` for the signed charge, which is the honest size of the
cancellation the sum survived. Non-positive cells are counted, never clipped:
`4H^2 - 2K < 0` is a discretisation artefact and hiding it would make a mesh look
smoother than it is.

## Limitations, measured

* **Pointwise mesh curvature stays noisy.** Agreement with the exact field
  (`artifacts/curvature_crosscheck.json`, 40 probes per case, ball/normal-cycle
  conversion): sphere median `H` 1.6 %, `K` 0.05 %; ellipsoid median `H` 12.8 %,
  `K` 22.4 %; Wigley flank patch 40 % / 100 %. The energies are integrals and are
  accurate to <0.1 %; single-vertex values are not, and 15 % cannot be claimed for
  them on stretched or open meshes.
* **The `(H, K) -> (k1, k2)` inversion is unusable pointwise.** On a conditioned
  marching-tet sphere, `H^2 - K < 0` at half the vertices because the cotangent
  Laplacian and the angle defect are different discretisations. `principal_curvature_field`
  reports those vertices as `nan` and counts them; the energy path therefore never
  routes through the inversion (`bending := 4H^2 - 2K` directly).
* **Raw marching-tet meshes break the mesh energy** (measured `+35 %`, `+140 %`, `+43 %`
  at levels 16/24/32, non-monotone): slivers carry almost no area but all of the
  pointwise error. Condition first (`taubin_smooth` + `project_to_field`), or use the
  field.
* **The analytic reference covers the classic Wigley (`n = m = 2`) exactly**; for other
  exponents it is quadrature-converged only while `1 < m < 2` keeps `(w/T)^(m-2)`
  integrable at the waterline.

## Provenance and the meshmagick clause (honest status)

The A3 acceptance clause asks for a **meshmagick** comparison on three cases.
`meshmagick` could not be obtained in this environment:

* `import meshmagick` -> `ModuleNotFoundError`;
* `pip download meshmagick --no-deps` -> no versions on the internal index
  (`lpypi.loongnix.cn`), ties to 2026-11-13;
* egress is blocked by the sandbox network policy, so neither a clone nor a
  `git+https` wheel is reachable.

The clause is therefore **WITHHELD**, not satisfied. What was executed instead, and
what `tools/crosscheck_curvature_thirdparty.py` records in
`artifacts/curvature_crosscheck.json`:

* three cases -- sphere, ellipsoid, Wigley flank patch -- against **the exact field**
  (validated to 2e-7 on the sphere) and against an independent **normal-cycle**
  operator (Cohen-Steiner & Morvan, the measure family `trimesh.curvature` implements);
* the harness prefers `meshmagick`, then `trimesh`, then the in-repo normal-cycle
  implementation, and it *refuses to guess* a meshmagick API it has never seen: if the
  package becomes importable the run stops with an instruction instead of pretending.

Result of the substituted comparison: the sphere and the ellipsoid integrals/charges
agree (`K` 0.05 % median on the sphere, total Gaussian measure `= 4 pi` exactly), the
pointwise numbers do not reach 15 % on stretched or open meshes (see limitations).
The honest reading is: **integral-level curvature diagnostics are REAL; pointwise mesh
curvature is a noisy estimate that this task did not improve to a 15 % guarantee.**
