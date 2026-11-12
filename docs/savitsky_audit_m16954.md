# Audit: `tensorcad.savitsky` — does the module solve what it claims? (m16954)

**Mission**: m16954 (verifier seat, selfstart)
**Target**: `tensorcad/savitsky.py` as merged in #359 (mission m16869), `docs/savitsky.md`,
`examples/savitsky/run_example.py`, `tests/test_savitsky.py`
**Method**: independent re-derivation of the execution path (not a re-run of the author's
tests), then probes of the lift closure, the wetted-length solve, the hull-length cap, the
validity box, and the two lift correlations.
**Evidence bundle**: `examples/savitsky_audit/run_example.py` →
`artifacts/savitsky_audit/results.json`, `report.md`, four figures. Every claim below is
asserted by that script (it exits non-zero if a claim stops holding) and pinned by
`tests/test_savitsky_audit.py` (42 tests).
**Independence limitation**: the primary reference (Savitsky 1964, *Marine Technology*
1(1) 71–95) is not available offline in this environment. Correlation coefficients could
therefore be compared against each other and against physical invariants, but not against
the paper. Findings that depend on the paper are marked **WITHHELD**, not "confirmed error".

---

## 1. Why audit this module

The module was merged with 28 green tests and a passing example. Its tests assert
*self-consistency* (the code reproduces its own formula, lift closes on weight, trends
have the right sign) — no test asserted that the *reported wetted length ratio* is a
result of the solve rather than an input. A verification seat exists to catch exactly
that: green tests that execute the code without asserting the claimed capability.

## 2. Findings

Severity: **B**locking (claim false / result unusable), **M**ajor (result misleading),
**m**oderate (diagnostics/hygiene), **W** = evidence withheld.

### F1 — BLOCKING — the "equilibrium wetted length ratio" was never solved

`running_trim()` initialised `lam = 3.0`, solved only the lift equation for τ, and
returned `(tau, lam)` with that untouched initial guess. `planing_prediction()` then
reported `wetted_length_ratio`, `wetted_length_m`, the wetted area, the friction drag and
therefore the total resistance and effective power as functions of a constant.

Evidence (`results.json → wetted_length_legacy`): six physically unrelated cases
(Δ = 8–40 kN, V = 8–18 m/s, B = 2.0–3.0 m, β = 10–30°) return **one identical λ = 3.000**
while τ varies from 0.4° to 14.9°.

Consequence: every wetted-geometry-dependent number in the module was insensitive to the
craft's weight, speed, beam and LCG. `docs/savitsky.md` nonetheless listed "wetted length
ratio λ" among the predictions.

Root cause: one equation (`CL(τ,λ,β) 0.5ρV²B² = Δ`) in two unknowns. The docstring claimed
"Savitsky provides a second equation from the pressure distribution … which gives a
relationship between τ and λ", but no second equation was implemented and no LCG was ever
an input to the module.

### F2 — MAJOR — the trim solve need not have iterated, and its failure mode was misreported

The lift correlation is invertible in closed form:
`τ = (CL_target (β/1000)^0.33 / (0.00222 λ^0.67))^(1/0.9)`.
The Newton loop's real job was to discover, 50 iterations later, that a clamp at 15° had
been hit; it then raised `Newton iteration failed to converge`, which reads as a numerical
problem and is actually "this craft needs 47° of trim" (evidence:
`closed_form_trim.out_of_window_message`). Fixing this changes nothing physical — the
closed form reproduces the lift balance to 0.0e+00 relative error — but it removes a
misleading diagnostic and the wasted iterations.

Related (m): the same defect made many legitimate-looking heavy/narrow cases fail with a
message that gave the caller no way to distinguish "not planing" from "solver broken".

### F3 — MAJOR — the hull-length cap silently broke lift equilibrium

`planing_prediction()` clamped `λ` to `L/B` when `Lw > L` and then recomputed the drag
components at the clamped λ, so the returned `lift_coefficient` no longer balanced the
weight. No warning, no flag. Repro (pre-fix): `planing_prediction(15000, 12, 5.0, 2.5, 20)`
→ λ reported 2.0 with `CL·qB²` ≠ Δ.

After the fix the clamp is announced in `warnings` and the trim is re-solved at the capped
λ, so lift equilibrium holds again (measured residual 2.4e-16 relative) and the caller
learns that the craft is in the hump / fully wetted.

### F4 — MAJOR — the documented validity box was never evaluated

The module docstring, `docs/savitsky.md` and the example all state a validity box
(Fr = 1–4, β = 10–30°, L/B = 3–6) and that "outside this box the predictions are
extrapolations". Nothing in the code computed Fr, or compared β or L/B to the box, or
reported extrapolation. A result at Fr = 0.74 (hump/displacement) was returned
indistinguishable from one at Fr = 1.8. Now every prediction carries a `warnings` tuple
enumerating each violated bound, and a case fully inside the box returns an empty tuple
(tested both ways).

### F5 — MODERATE — a shortened speed sweep was indistinguishable from a complete one

`speed_sweep()` swallowed every `PlaningPredictionError` with an unused `as e` binding.
A caller asking for Fr = 0.5–3.0 received a shorter list with no indication that half the
range had failed. `failures=[...]` now records `(speed, reason)` and
`len(results) + len(failures) == len(speeds)` is asserted.

### F6 — WITHHELD — two incompatible lift correlations are both labelled "Savitsky (1964) eq. (1)"

The module's `lift_coefficient` is
`CL = 0.00222 τ^0.9 λ^0.67 (β/1000)^-0.33`, documented as "equation 1 of the 1964 paper".
It contains no speed coefficient `Cv`. The classical statement of Savitsky's lift is
`CL_β = τ^1.1 (0.0120 λ^0.5 + 0.0055 λ^2.5 / Cv²)` with the deadrise correction
`CL_0 = CL_β − 0.0065 β CL_β^0.6`, which does depend on `Cv`.

Measured (`results.json → lift_models`, figure 1): over τ = 2–10°, λ = 2–4, β = 20° the
module's correlation is **22 % to 78 % lower** than the classical one at Cv = 2, and the
best-fitting single `Cv` still leaves an **88 % worst-case disagreement** over the whole
grid. Two formulas that disagree by tens of percent cannot both be eq. (1) of the same
paper; at least one citation is wrong, and no single-`Cv` reconciliation exists.

I did **not** rewrite the module's formula: I cannot check the primary paper here, and
inventing a coefficient would be worse than flagging one. The incompatibility is pinned by
`tests/test_savitsky_audit.py::TestCorrelationIncompatibility` so that any future
unification forces this document to be revisited.

Same status applies to the pressure-drag term `R_P/Δ = tan τ + 0.008/√τ`
(`resistance_components`): the `0.008/√τ` term is not something this audit can source. Note
only that it is **not small** at low trim — at τ = 1° it is 46 % of the `tan τ` term, at
τ = 8° it is 2 % — so any error in it is concentrated exactly where planing craft run.

### F7 — the practical consequence (what the owner should know)

Comparing the two paths on one craft (Δ = 15 kN, B = 2.5 m, L = 12 m, β = 20°, LCG = 3.0 m;
`report.md`):

| V [m/s] | Fr | solved τ | legacy τ | solved λ | legacy λ | solved R_T | legacy R_T |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 8 | 0.74 | 3.06° | 5.11° | 2.12 | 3.00 | 1998 N | 2904 N |
| 12 | 1.11 | 2.50° | 2.07° | 1.80 | 3.00 | 2807 N | 3825 N |
| 20 | 1.84 | 1.52° | 0.67° | 1.67 | 3.00 | 5488 N | 8571 N |

The legacy path over-predicts resistance by **+21 % to +56 %** and the two paths disagree
on trim by up to **2.3×**. Both are far outside the ±10 % band the module claims for the
method, so the choice of path — not the correlation's published accuracy — dominates the
current error budget. Until F6 is resolved with the paper, `tensorcad.savitsky` should not
be quoted as a <10 % resistance predictor.

## 3. What was changed (minimal, additive-first)

`tensorcad/savitsky.py`:

* **New, and the only path that solves λ**: `equilibrium_trim()` — the classical
  two-equation closure (lift = weight, centre of pressure under the LCG) with
  `classic_lift_coefficient()`, `center_of_pressure_ratio()`, `speed_coefficient()` and a
  `TrimSolution` result carrying both residuals. Bisection, no tuned step sizes; residuals
  4.7e-15 (lift) and 1.2e-10 m (moment) at every speed tested.
* `planing_prediction(..., lcg_m=…)` routes to that closure and reports
  `solution_method='classic_equilibrium'`, `wetted_length_ratio_solved=True`. Without
  `lcg_m` the legacy one-equation behaviour is preserved **verbatim** (default λ = 3.0) but
  is now labelled `'assumed_wetted_length'` and carries a warning. No existing caller
  breaks: the 28 pre-existing tests and `examples/savitsky/run_example.py` still pass
  unchanged.
* `running_trim()`: λ becomes an explicit input (`wetted_length_ratio`, default 3.0) that
  is echoed back, with an honest docstring; the misleading iteration/clamp path is replaced
  by the exact closed form `trim_for_wetted_length()`.
* `planing_prediction()`: the hull-length cap is announced, lift equilibrium is restored
  after clamping, and `validity_warnings()` is folded into `warnings`.
* `speed_sweep(..., lcg_m=…, failures=[…])`: failures are recordable.
* `PlanPrediction` gains `wetted_length_ratio_solved`, `solution_method`, `warnings`
  (defaulted, serialised by `as_dict()`).

`docs/savitsky.md` records the audit status of each claim; `CHANGELOG.md` has the entry.

## 4. Checks run / not run

Run:

* `python -m pytest tests/test_savitsky_audit.py tests/test_savitsky.py -q` → 70 passed.
* `python examples/savitsky_audit/run_example.py` → 10/10 claims PASS, 4 figures.
* `python examples/savitsky/run_example.py` → all pre-existing claims still pass.
* Full suite `python -m pytest tests/ -q` (recorded in the mission report).

**Not run** (and why):

* Comparison against Savitsky (1964) itself, or against model-test data (Sottorf, Kapryan &
  Boyd, Fridsma): no offline copy of any of them in this environment. This is what keeps
  F6 at WITHHELD and the classical coefficients provisional.
* Comparison against `ship-performance-platform`'s planing resistance path: not reachable
  from this work tree; a cross-code check with an independently written Savitsky
  implementation would be the cheapest way to settle F6 and is recommended to the owner.
* Any validation of the pressure-drag `0.008/√τ` term (F6b) — same reason.

## 5. Findings deliberately not "fixed"

* The module's `lift_coefficient` formula is untouched (F6): changing a published
  correlation on the strength of a memory of the paper would launder an unverified number
  into the library.
* `resistance_components`' wetted area `Lw·B/cos β` is untouched: it is a defensible
  first-order idealisation for a constant-deadrise prismatic hull, and this audit has no
  better reference to replace it with. It is now simply one more reason not to quote the
  module as ±10 %.
* The default path remains the legacy one-equation path. Switching the default to a
  WITHHELD correlation would be a bigger change than the evidence supports.

## 6. Disposition

**Request changes, delivered partially → approve with notes.**
The blocking defect (F1) and the three major defects (F2, F3, F4) plus the sweep hygiene
(F5) are fixed and pinned. F6 is escalated as an evidence gap: it needs the primary paper
or a second independent implementation, and until then the module must be described as a
correlation implementation of unknown absolute accuracy, not as a ±10 % predictor.

Recommended next actions, in order:

1. Obtain Savitsky (1964) eqs. (1)–(3) and settle which lift/CP form is correct; if the
   classical form is right, `lift_coefficient` needs replacing — a two-function change
   (`classic_lift_coefficient`, `center_of_pressure_ratio` are already in place) plus a
   re-pin of `TestCorrelationIncompatibility` and `docs/savitsky.md`'s table.
2. Cross-check `equilibrium_trim` against `ship-performance-platform`'s planing module on
   3–5 of the table's cases; agreement to <5 % would upgrade F6's evidence level.
3. Add a `planing` example that drives `lcg_m` end-to-end (geometry → trim → resistance)
   so the closed loop is demonstrated on a hull rather than on a parameter set.

## 7. Reproduction

```bash
python -m pytest tests/test_savitsky_audit.py -q            # 42 audit tests
python examples/savitsky_audit/run_example.py               # claims + artifacts
python -m pytest tests/test_savitsky.py -q                  # 28 pre-existing tests
```