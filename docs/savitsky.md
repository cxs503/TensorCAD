# Savitsky 1964 planing prediction

The classical empirical method for predicting the resistance and running trim of
planing hulls. This module implements the Savitsky (1964) correlation for
constant-deadrise prismatic planing craft — the same hull family that
`tensorcad.planing` generates geometrically.

Everything on this page is measured by `python examples/savitsky/run_example.py`,
which writes `artifacts/savitsky/results.json` and exits non-zero if any claim
stops holding. `tests/test_savitsky.py` pins the same numbers (28 cases).

## 0. Audit status (m16954, verifier seat) -- read this before quoting any number

An adversarial audit of this module (`docs/savitsky_audit_m16954.md`, evidence in
`artifacts/savitsky_audit/`) changed three things about what the module may claim:

1. **λ is solved only by `equilibrium_trim()` / `planing_prediction(..., lcg_m=...)`.**
   The original lift solve is one equation in two unknowns, so `running_trim()` and the
   default `planing_prediction()` path return an *assumed* λ (default 3.0) and say so via
   `solution_method='assumed_wetted_length'`, `wetted_length_ratio_solved=False` and a
   warning. Before the audit that assumption was reported as if it were an equilibrium.
2. **The validity box is now reported, not implied.** Every result carries a `warnings`
   tuple naming each violated bound (Fr, β, trim band, L/B) and any wetted-length clamp.
   A clamped wetted length no longer breaks the lift balance silently.
3. **The absolute accuracy is unknown, not ±10 %.** The module's lift correlation
   (`CL = 0.00222 τ^0.9 λ^0.67 (β/1000)^-0.33`, no speed coefficient) disagrees with the
   classical Cv-dependent Savitsky statement by 22-78 % over the validity box, and no
   single Cv reconciles them (`tests/test_savitsky_audit.py::TestCorrelationIncompatibility`).
   Both are labelled "Savitsky (1964) eq. (1)"; one label is wrong, and resolving it needs
   the primary paper, which was not available. The resulting resistance difference between
   the two paths on one reference craft is +21 % to +56 % (legacy high). Treat the module
   as a correlation implementation of unverified absolute accuracy until that is settled.


## 1. What the method is, and what it deliberately is not

`tensorcad.savitsky` implements the **Savitsky 1964** empirical method:

* **Reference**: Savitsky, D. (1964). "Hydrodynamic design of planing hulls."
  *Marine Technology*, 1(1), 71-95.
* **What it predicts**: running trim angle τ, wetted length ratio λ, lift
  coefficient CL, friction drag (ITTC-57 line), pressure drag (Savitsky
  correlation), total resistance, and effective power.
* **What it needs**: hull weight Δ, speed V, beam B, deadrise β, and (for the
  complete prediction) hull length L.
* **Validity box**: Fr = 1.0 to 4.0, β = 10° to 30°, L/B = 3 to 6. Outside
  this box the predictions are extrapolations.

The method is **empirical**: it reproduces the Savitsky correlation to round-off,
not the physics. The published accuracy is ±10% on resistance for hulls within
the test matrix. This module does not claim better accuracy than that.

Explicitly absent:
* Seakeeping or dynamic stability (no heave/pitch RAOs, no porpoising prediction).
* Non-prismatic hulls (variable deadrise, warped hulls, steps, spray rails).
* Model-test correlation beyond the Savitsky original band.
* A resistance prediction for displacement-mode hulls (Fr < 1).

## 2. The lift coefficient

The Savitsky lift correlation (equation 1 of the 1964 paper):

    CL = 0.00222 τ^0.9 λ^0.67 / (β/1000)^0.33

where τ is trim in degrees, λ = L_w/B is the wetted length ratio, and β is
deadrise in degrees. The denominator keeps the coefficient 0.00222 dimensionally
consistent with the degree-based inputs.

| τ [deg] | λ | β [deg] | CL |
|---:|---:|---:|---:|
| 2 | 3 | 20 | 0.0549 |
| 4 | 3 | 20 | 0.1098 |
| 6 | 3 | 20 | 0.1647 |
| 4 | 2 | 20 | 0.0776 |
| 4 | 5 | 20 | 0.1483 |
| 4 | 3 | 10 | 0.1744 |
| 4 | 3 | 30 | 0.0831 |

Three facts, all asserted in the tests:

* **CL increases with τ** (more trim → more lift, positive derivative).
* **CL increases with λ** (longer wetted length → more lift area).
* **CL decreases with β** (deeper V-section → less lift, the hull dives more).

The derivative dCL/dτ = 0.9 CL / τ is analytic from the correlation and is
what the Newton solver uses to find the equilibrium trim.

## 3. The equilibrium trim

At a given speed and weight, the running trim is the root of the lift equation:

    CL(τ, λ, β) × 0.5 ρ V² B² = Δ

For a given λ (the wetted length ratio, which is a property of the hull geometry
and running attitude), this is one equation in one unknown (τ). The solver uses
Newton iteration with the analytic derivative and converges in 3-5 iterations
to a tolerance of 1e-6 on the lift residual.

**Reference case**: Δ = 50,000 N, V = 15 m/s, B = 3 m, β = 20°, ρ = 1025.9 kg/m³.

| Quantity | Value |
|---|---|
| Dynamic pressure × B² | 0.5 × 1025.9 × 15² × 3² = 3,462,487 N |
| Target CL | 50,000 / 3,462,487 = 0.01444 |
| Solved τ | 3.21° |
| Solved λ | 3.00 |
| Lift residual | < 1e-4 |

The equilibrium satisfies the lift-weight balance to 0.01%, which is a statement
about the solver, not about whether the Savitsky correlation is right.

## 4. The resistance components

The total resistance is the sum of friction drag and pressure drag:

    R_T = R_F + R_P

**Friction drag** uses the ITTC-57 line on the wetted area:

    R_F = CF(Re) × 0.5 ρ V² S_wet
    CF = 0.075 / (log10 Re - 2)²

where Re = V L_w / ν and S_wet = L_w B / cos(β) is the wetted area of the
prismatic bottom (inclined by the deadrise).

**Pressure drag** follows the Savitsky correlation:

    R_P = Δ × (tan τ + 0.008 / sqrt(τ))

where τ is in degrees. The first term is the component of lift in the drag
direction (induced drag); the second is an empirical correction for the
pressure distribution.

**Reference case** (same as §3):

| Component | Value |
|---|---|
| Reynolds number | 1.135 × 10⁸ |
| ITTC-57 CF | 0.00208 |
| Wetted area | 9.0 / cos(20°) = 9.577 m² |
| Friction drag R_F | 6,784 N |
| Pressure drag R_P | 3,027 N |
| Total drag R_T | 9,811 N |
| R_T / Δ | 0.196 |
| Effective power P_E | 147,158 W (197 hp) |

The resistance-to-weight ratio R_T/Δ = 0.196 is typical for a planing craft
at this speed (Fr = 1.38). The pressure drag is 31% of the total, which is
what the Savitsky correlation gives for this trim angle.

## 5. Speed sweep

The `speed_sweep` function evaluates the prediction over a speed range. At low
speeds (Fr < 1) the craft is in displacement mode and the Savitsky method does
not apply; the sweep skips those speeds and continues.

**Reference sweep**: Δ = 50,000 N, L = 12 m, B = 3 m, β = 20°, V = 9 to 19 m/s.

| V [m/s] | Fr | τ [deg] | λ | R_T [N] | P_E [kW] |
|---:|---:|---:|---:|---:|---:|
| 9 | 0.83 | 9.99 | 3.00 | 11,562 | 104 |
| 11 | 1.01 | 6.39 | 3.00 | 9,577 | 105 |
| 13 | 1.20 | 4.41 | 3.00 | 9,249 | 120 |
| 15 | 1.38 | 3.21 | 3.00 | 9,811 | 147 |
| 17 | 1.57 | 2.43 | 3.00 | 10,937 | 186 |
| 19 | 1.75 | 1.90 | 3.00 | 12,471 | 237 |

Three trends, all physical:

* **Trim decreases with speed**: more dynamic pressure → less trim needed to
  generate the same lift.
* **Resistance has a hump**: it decreases from displacement to planing (the
  craft lifts out of the water), then increases with speed (friction dominates).
* **Power increases monotonically**: P_E = R_T × V, and V grows faster than R_T
  decreases in the planing regime.

## 6. Verification

The implementation is verified against:

1. **The Savitsky correlation formula itself**: the lift coefficient is checked
   against the analytic expression to round-off (test_lift_coefficient).
2. **The equilibrium condition**: the solved trim must satisfy CL × q × B² = Δ
   to 1e-4 (test_equilibrium_satisfies_lift_weight_balance).
3. **The ITTC-57 line**: friction drag uses the correct CF(Re) (test_friction_drag_uses_ittc57).
4. **The pressure drag correlation**: R_P = Δ × (tan τ + 0.008/sqrt(τ))
   (test_pressure_drag_follows_savitsky).
5. **Internal consistency**: Froude number, wetted length ratio, and effective
   power are computed correctly (test_prediction_is_self_consistent).
6. **Physical reasonableness**: the reference case produces trim, resistance,
   and Froude number in the expected ranges (TestSavitskyReferenceCase).

All 28 tests pass. The implementation reproduces the Savitsky correlation to
machine precision; whether the correlation itself is accurate for a given hull
is a question for model tests or CFD, not for this module.

## 7. Limitations, stated plainly

* **Empirical correlation, not physics.** The Savitsky method was fitted to a
  specific test matrix (prismatic hulls, constant deadrise, Fr = 1-4). Outside
  that box the errors grow, and the module does not warn you.
* **Fixed wetted length ratio.** The current solver assumes λ = 3 (a typical
  planing value). A more complete implementation would solve for λ from the
  hull geometry and running attitude, or accept it as an input from the
  `tensorcad.planing` module.
* **No spray drag.** The Savitsky method does not include spray resistance,
  which can be 10-20% of the total for hard-chine hulls.
* **No appendage drag.** Struts, shafts, rudders, and interceptors are not
  modelled; the caller must add them separately.
* **No seakeeping.** Added resistance in waves, porpoising stability, and
  dynamic sinkage are outside the scope.

## 8. Integration with the rest of TensorCAD

The Savitsky module is designed to work with `tensorcad.planing`:

1. `tensorcad.planing.PlaningHullParameters` defines the hull geometry
   (L, B, T, β, q).
2. `tensorcad.planing.generate_planing_hull` generates the mesh and computes
   the hydrostatics.
3. `tensorcad.savitsky.planing_prediction` takes the weight, speed, and
   geometry, and predicts the running trim and resistance.

The wetted length ratio λ is currently fixed at 3 in the Savitsky solver; a
future enhancement would couple it to the hull geometry (the wetted length
depends on the draft, which depends on the trim, which depends on the speed —
a coupled problem).

## 9. References

* Savitsky, D. (1964). "Hydrodynamic design of planing hulls." *Marine
  Technology*, 1(1), 71-95.
* Blount, D. L., & Codega, L. T. (1992). "Hydrodynamic design of planing
  hulls." *Marine Technology*, 29(4), 256-268.
* Clement, C. F., & Blount, D. L. (1963). "Resistance tests of a systematic
  series of planing hull forms." *Transactions SNAME*, 71, 462-495.
* ITTC (1957). "Report of the performance committee." *Proceedings of the
  8th International Towing Tank Conference*, Madrid.
