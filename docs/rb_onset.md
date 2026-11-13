# Rayleigh-Benard onset — reference value, criterion and solver bundle

Module: `tensorcad/rb_onset.py` · Example: `examples/rb_onset/run_example.py` ·
Tests: `tests/test_rb_onset.py` (23) · Artifacts: `artifacts/rb_onset/`

Case-pool entry: **`LBM-RB-ONSET-RA1708`** — "Rayleigh-Benard 对流起始临界 Ra_c",
acceptance `临界 Rayleigh 数判定值 vs 1707.76 (Ra 二分扫描+静止态线性稳定性判据)`,
tolerance `< 2%`, reference Chandrasekhar (1961) ch. II (rigid-rigid plates).

## 0. What this module is, and what it is not

TensorCAD supplies the **geometry, boundary-condition and reference-value** side of
the benchmark, exactly as it already does for the lid-driven cavity
(`tensorcad/cavity_mesh.py`, `examples/cavity_re5000/`) and for Kovasznay flow
(`tensorcad/kovasznay.py`). The time-marching solver lives in TensorLBM / TensorFVM
/ TensorDG; `tensorcad` core stays dependency-free (standard library only, verified
by the same constraint the rest of the package obeys).

| Part of the benchmark | Where it lives | Evidence grade here |
|---|---|---|
| Critical Rayleigh number and wavenumber (reference value) | `tensorcad/rb_onset.py` | **REAL** — spectral collocation + Ra bisection, `1.3e-7` from Chandrasekhar (Table 1) |
| Exact stress-free benchmark `(π²+a²)³/a²` | `tensorcad/rb_onset.py` | **REAL** — `3.0e-9` worst relative error (Table 2) |
| Cell geometry, periodic/no-slip/isothermal BC tags, mesh family, monitor contract | `rb_cell()` | **REAL** — consumed by a solver, `artifacts/rb_onset/rb_onset_bundle.json` |
| Critical eigenfunction (mode shape) | `onset_mode()` | **REAL** — dimensionless residual `7e-14` |
| Ra-bisection driver for a solver's growth-rate probe | `bisect_onset()` | **REAL** — self-tested against a planted onset (`1e-7` recovery) |
| The LBM/FVM run of the case itself | TensorLBM / TensorFVM | **not in this repo** — see §6 |

## 1. Physics and the criterion

A Newtonian layer between two infinite horizontal plates, heated from below, is
non-dimensionalised with the plate separation `d`, the temperature drop `ΔT`,
velocity `κ/d` and time `d²/κ`. The conduction state `T = 1 - z` (bottom plate hot at
`z = 0`) is an exact solution of the Boussinesq equations for every `Ra`. It loses
stability by *exchange of stability* (stationary rolls) — the criterion the case
entry names, "静止态线性稳定性判据".

Linearising about the conduction state and taking normal modes `∝ exp(i a x)`:

```
(D² - a²)² Ψ - Ra Θ = 0                       (1a)   Ψ = stream function amplitude
(D² - a²)   Θ + a² Ψ = 0                      (1b)   Θ = temperature amplitude
```

with `D = d/dz`. Eliminating `Θ` recovers the classical sixth-order form
`(D² - a²)³ Ψ = -Ra a² Ψ` (Chandrasekhar 1961, ch. II). Plate conditions:

| plates | kinematic | thermal |
|---|---|---|
| rigid (no-slip) | `Ψ = DΨ = 0` | `Θ = 0` |
| free (stress-free) | `Ψ = D²Ψ = 0` | `Θ = 0` |

`Ra` enters (1a) linearly, so `det M(a, Ra)` is a polynomial in `Ra`. The conduction
state is linearly stable iff the only stationary perturbation is trivial, i.e. iff
`det M ≠ 0`; **the first sign change of `det M` while scanning `Ra` upward is the
marginal state of the lowest mode**, and `Ra_c = min_a Ra(a)`. That is literally the
case-pool criterion: an Ra bisection scan gated by the stationary stability
condition.

The stationary onset is independent of the Prandtl number (no `Pr` appears in (1)) —
so the reference applies to any LBM relaxation time, which is what makes it usable as
an LBM benchmark.

## 2. Method

* Chebyshev-Lobatto collocation in `z ∈ [0,1]` on the coupled system (1a)-(1b); the
  four `Ψ` conditions and the two `Θ` conditions replace the operator rows adjacent
  to each wall ("boundary bordering"), leaving `2n+2` equations for `2n+2` unknowns.
* `det M` is accumulated in log space (LU with partial pivoting); the bordered matrix
  carries entries up to `n⁸`, so a plain determinant overflows double precision.
* `Ra` scan: multiplicative bracket (`×1.25`) from `Ra = 1`, then bisection to a
  relative tolerance (default `1e-6`).
* `min_a`: wavenumber scan → parabola through the three lowest samples (two-pass
  vertex refine) and a golden-section refine; the marginal curve is flat in `a`, so
  `a_c` is quoted to ~`1e-3` relative and guarded against the bisection noise floor.
* Independent cross-check: `σ_min(M)/σ_max(M)` minimised over `Ra` by SVD (no
  determinant, no sign logic) — same discretisation, different criterion.

## 3. Results (programmatically generated, `artifacts/rb_onset/checks.json`)

### 3.1 Case-pool verdict (Chebyshev order 28)

| quantity | measured | reference | relative error | acceptance (<2%) |
|---|---|---|---|---|
| `Ra_c`, rigid-rigid | **1707.761784** | 1707.762 (Chandrasekhar 1961) | **−1.26e-07** | **PASS** |
| `a_c`, rigid-rigid | **3.116750** | 3.11632 | +1.4e-04 | PASS |
| `Ra_c`, stress-free | **657.511433** | 657.511364 (= 27π⁴/4) | +1.04e-07 | PASS |
| `a_c`, stress-free | **2.221953** | 2.221441 (= π/√2) | +2.3e-04 | PASS |

### 3.2 Exact stress-free benchmark of the whole scan

`Ra(a) = (π² + a²)³/a²` is exact for free plates. Measured at `a = 1.5, π/√2, 3.0`:
worst relative error **2.97e-09** — the discretisation, the Ra bisection and the
`min_a` step are all validated against a closed form *before* the rigid-rigid claim is
made.

### 3.3 Chebyshev-order convergence (rigid-rigid, fixed `a = 3.11632`)

| order `n` | `Ra_c` | relative error vs 1707.762 |
|---:|---|---|
| 8 | 1707.563093 | −1.16e-04 |
| 10 | 1707.756643 | −3.14e-06 |
| 12 | 1707.761694 | −1.79e-07 |
| 16 | 1707.761777 | −1.30e-07 |
| 20 … 40 | 1707.761777 | −1.30e-07 |

Spectral convergence down to an **error floor of 1.3e-7 that is the six-digit rounding
of the published reference** (|1707.761777 − 1707.762| = 2.2e-4 absolute), not
discretisation error: the converged value is stable in eight digits from order 16 to
order 40 (spread < 1e-9). Error at order 8 is already 1.2e-4, i.e. 900× inside the 3%
band.

### 3.4 Independent criteria

| route | result | note |
|---|---|---|
| determinant sign flip | `det M(0.9 Ra_c)` and `det M(1.1 Ra_c)` have opposite signs | the destabilisation really is a marginal crossing |
| SVD minimisation of `σ_min/σ_max` over `Ra` | `Ra = 1707.761776` (rel −1.31e-07), `σ_min/σ_max = 5.2e-19` | no determinant, no sign tracking — independent criterion |

### 3.5 Critical mode (compare field by field, not only the eigenvalue)

`onset_mode()` returns `Ψ(z)`, `Θ(z)` on the Chebyshev nodes with dimensionless
residual `‖Mx‖/(‖M‖‖x‖) = 7e-14`; `Ψ` and `Θ` vanish on both plates, the roll core
sits at mid-depth. The 2-D reconstruction `ψ = Ψ(z)cos(a_c x)`,
`θ = Θ(z)cos(a_c x)` is what `artifacts/rb_onset/rb_critical_mode.png` plots.

## 4. Solver-side harness

`bisect_onset(probe, rayleigh_stable=…, rayleigh_unstable=…)` bisects `Ra` with a
solver-supplied growth rate `σ(Ra)` (positive = unstable). The bracket is validated
first: a probe that does not change sign raises `ValueError` instead of returning a
fabricated number. `growth_rate(times, amplitudes)` is the matching log-linear probe
estimator. Self-test with a planted onset at `Ra = 1234.5`: recovered to `1e-7`
relative in 22 iterations — the driver discriminates, so a solver reporting
`Ra_c ≈ 1707.8` is reporting its own discretisation error, not a harness artefact.

## 5. Notes, limits, and traps met on the way

* **Determinant overflow.** `det` of the bordered matrix (entries `~n⁸`) overflows
  double precision for `n ≥ 20`; hence `determinant_sign_log`. A naive `numpy.det`
  returns `inf`/`nan` and the sign scan silently degrades.
* **The marginal curve is flat in `a`** (quadratic with curvature ~2.4e2 near the
  minimum, bisection noise ~2e-3 in `Ra`): `a_c` is therefore only good to ~1e-3
  relative, while `Ra_c` is good to `1e-7`. Quoting `a_c` to five digits would be
  dishonest.
* **A generalized-eigenvalue route via `scipy.linalg.eig(M0, -dM)` is unusable here**:
  the pencil's `B` block is rank deficient (BC rows are empty in `B`) and the Chebyshev
  scaling drives the conditioning through the floor; it returned `1e10`-scale spurious
  eigenvalues. The SVD route replaces it and is stable.
* The Chebyshev derivative chain `z = (1-ξ)/2` needs `d/dz = -2 d/dξ`; `D¹` enters the
  stability operator only squared or as `DΨ = 0`, so a sign slip there is invisible in
  `Ra_c` — it was caught by the polynomial-exactness test, not by the physics checks.
* Scope: two infinite parallel plates, Oberbeck-Boussinesq, stationary onset only
  (no overstability / `Pr`-dependent oscillatory modes, no finite-aspect or sidewall
  effects, no temperature-dependent properties).

## 6. Hand-off

`artifacts/rb_onset/rb_onset_bundle.json` (sha256-protected) contains the case spec,
the measured onset, the cell geometry/BC tags/mesh family/probe points and the mode
profile. A solver run should: (i) build the cell from `case["cell"]` with
`width = 2π/a_c ≈ 2.0165 d`, (ii) start from the conduction state plus a small
`sin`-roll perturbation, (iii) fit `σ` over the exponential window with
`growth_rate()`, (iv) bisect `Ra` with `bisect_onset()`, and (v) check `Nu = 1` below
and `Nu > 1` above the onset. Executing that LBM run is TensorLBM's side of the
case-pool entry, not TensorCAD's.

## References

* Chandrasekhar, S. (1961). *Hydrodynamic and Hydromagnetic Stability*, Oxford
  University Press, ch. II (Benard problem), table III: rigid-rigid
  `Ra_c = 1707.762`, `a_c = 3.116`. doi:10.1016/0368-3281(63)90025-1 (secondary).
* Drazin, P. G. & Reid, W. H. (1981). *Hydrodynamic Stability*, Cambridge University
  Press, ch. II.
* de Vahl Davis, G. (1983). *Int. J. Numer. Methods Fluids* 3, 249 — `Nu` benchmarks
  for the supercritical companion entries `CASE-RAYLEIGH-BENARD-RA1E4/1E5`.
* Berrut, J.-P. & Trefethen, L. N. (2004). *SIAM Review* 46, 501 — barycentric
  interpolation.
* Trefethen, L. N. (2000). *Spectral Methods in MATLAB*, SIAM — Chebyshev operators.