"""Audit tests for :mod:`tensorcad.savitsky` (m16954, verifier seat).

Each test pins one finding of ``docs/savitsky_audit_m16954.md`` so that the defect
cannot come back silently:

* the trim solve used to report λ = 3.0 as an "equilibrium" wetted length ratio;
* the wetted-length cap used to break lift equilibrium without saying so;
* the validity box was documented but never reported;
* the failure mode of the trim solve was reported as "convergence", not as a window
  violation;
* a speed sweep used to drop failing speeds without a trace;
* the module's lift correlation and the classical Savitsky lift correlation differ by
  tens of percent, which cannot both be Savitsky (1964) eq. (1).
"""
from __future__ import annotations

import math

import pytest

from tensorcad.savitsky import (
    DEFAULT_WETTED_LENGTH_RATIO,
    PlaningPredictionError,
    center_of_pressure_ratio,
    classic_lift_coefficient,
    equilibrium_trim,
    lift_coefficient,
    planing_prediction,
    running_trim,
    speed_coefficient,
    speed_sweep,
    trim_for_wetted_length,
    validity_warnings,
)

RHO = 1025.9
G = 9.80665


def dynamic_pressure(weight_n: float, speed_ms: float, beam_m: float) -> float:
    """q B² for the lift equation."""
    return 0.5 * RHO * speed_ms ** 2 * beam_m ** 2


class TestWettedLengthWasNotSolved:
    """Finding F1 (blocking): λ was the initial guess, not a solution."""

    @pytest.mark.parametrize(
        "weight,speed,beam,deadrise",
        [
            (8000.0, 8.0, 2.5, 20.0),
            (12000.0, 10.0, 2.0, 10.0),
            (15000.0, 12.0, 2.5, 20.0),
            (20000.0, 15.0, 3.0, 25.0),
            (26000.0, 9.0, 2.2, 30.0),
            (40000.0, 18.0, 3.0, 15.0),
        ],
    )
    def test_running_trim_echoes_the_assumed_wetted_length_ratio(
        self, weight, speed, beam, deadrise
    ):
        """Six physically different cases, one identical λ: it is an input, not a result."""
        _, lam = running_trim(weight, speed, beam, deadrise)
        assert lam == DEFAULT_WETTED_LENGTH_RATIO

    def test_running_trim_honours_an_explicit_wetted_length_ratio(self):
        for lam_in in (1.5, 2.0, 5.0):
            _, lam = running_trim(15000.0, 12.0, 2.5, 20.0, wetted_length_ratio=lam_in)
            assert lam == lam_in

    def test_prediction_labels_the_assumed_path(self):
        pred = planing_prediction(15000.0, 12.0, 8.0, 2.5, 20.0)
        assert pred.solution_method == "assumed_wetted_length"
        assert pred.wetted_length_ratio_solved is False
        assert any("assumed" in w for w in pred.warnings)

    def test_prediction_labels_the_solved_path(self):
        pred = planing_prediction(15000.0, 12.0, 12.0, 2.5, 20.0, lcg_m=4.5)
        assert pred.solution_method == "classic_equilibrium"
        assert pred.wetted_length_ratio_solved is True


class TestClosedFormTrim:
    """The lift equation is closed-form solvable; the old Newton loop hid that."""

    @pytest.mark.parametrize("lam", [1.0, 2.0, 3.0, 4.5])
    def test_trim_satisfies_the_lift_equation_exactly(self, lam):
        weight, speed, beam, beta = 15000.0, 12.0, 2.5, 20.0
        tau = trim_for_wetted_length(weight, speed, beam, beta, lam)
        lift = lift_coefficient(tau, lam, beta) * dynamic_pressure(weight, speed, beam)
        assert abs(lift - weight) / weight < 1e-12

    def test_required_trim_outside_the_window_is_reported_with_its_value(self):
        """Heavy craft + narrow beam: the message must carry the trim asked for."""
        with pytest.raises(PlaningPredictionError) as excinfo:
            trim_for_wetted_length(40000.0, 6.0, 2.0, 20.0, 3.0)
        message = str(excinfo.value)
        assert "required trim" in message
        assert "outside the implemented window" in message
        assert "failed to converge" not in message


class TestWettedLengthClampIsHonest:
    """Finding F3: the cap used to break lift equilibrium silently."""

    def test_clamped_wetted_length_keeps_the_lift_balance(self):
        # L/B = 2 < 3, so λ = 3 (Lw = 7.5 m) cannot fit in a 5 m hull.
        pred = planing_prediction(15000.0, 12.0, 5.0, 2.5, 20.0)
        assert pred.wetted_length_ratio == pytest.approx(5.0 / 2.5)
        assert pred.wetted_length_m == pytest.approx(5.0)
        lift = pred.lift_coefficient * dynamic_pressure(15000.0, 12.0, 2.5)
        assert abs(lift - pred.weight_n) / pred.weight_n < 1e-9
        assert any("clamped" in w for w in pred.warnings)

    def test_clamp_case_is_not_advertised_as_a_solved_equilibrium(self):
        pred = planing_prediction(15000.0, 12.0, 5.0, 2.5, 20.0)
        assert pred.wetted_length_ratio_solved is False


class TestEquilibriumSolver:
    """The new capability: solve (τ, λ) from lift *and* the centre-of-pressure moment."""

    def test_both_residuals_close(self):
        solution = equilibrium_trim(
            15000.0, 12.0, 2.5, 20.0, 4.5, max_wetted_length_m=12.0
        )
        assert solution.lift_residual < 1e-9
        assert solution.moment_residual_m < 1e-8
        assert solution.method == "classic_savitsky_lift_and_moment"

    def test_centre_of_pressure_sits_under_the_lcg(self):
        for lcg in (3.5, 4.0, 4.5, 5.0):
            solution = equilibrium_trim(
                15000.0, 12.0, 2.5, 20.0, lcg, max_wetted_length_m=12.0
            )
            cp = (
                center_of_pressure_ratio(
                    solution.wetted_length_ratio, solution.speed_coefficient
                )
                * solution.wetted_length_ratio
                * 2.5
            )
            assert cp == pytest.approx(lcg, abs=1e-7)

    def test_wetted_length_shrinks_with_speed(self):
        lams = [
            equilibrium_trim(
                15000.0, v, 2.5, 20.0, 4.5, max_wetted_length_m=12.0
            ).wetted_length_ratio
            for v in (8.0, 10.0, 12.0, 15.0, 18.0, 20.0)
        ]
        assert all(later < earlier for earlier, later in zip(lams, lams[1:])), lams

    def test_moving_the_lcg_aft_raises_the_trim(self):
        tau_aft = equilibrium_trim(15000.0, 12.0, 2.5, 20.0, 3.5, max_wetted_length_m=12.0).trim_deg
        tau_fwd = equilibrium_trim(15000.0, 12.0, 2.5, 20.0, 5.5, max_wetted_length_m=12.0).trim_deg
        assert tau_aft > tau_fwd

    def test_heavier_craft_trims_more(self):
        tau_light = equilibrium_trim(10000.0, 12.0, 2.5, 20.0, 4.5, max_wetted_length_m=12.0).trim_deg
        tau_heavy = equilibrium_trim(25000.0, 12.0, 2.5, 20.0, 4.5, max_wetted_length_m=12.0).trim_deg
        assert tau_heavy > tau_light

    def test_unreachable_lcg_is_an_honest_error_not_a_clamp(self):
        with pytest.raises(PlaningPredictionError) as excinfo:
            equilibrium_trim(15000.0, 12.0, 2.5, 20.0, 6.0, max_wetted_length_m=12.0)
        assert "centre of pressure can only be placed" in str(excinfo.value)

    @pytest.mark.parametrize(
        "kwargs",
        [
            {"weight_n": 0.0},
            {"speed_ms": -1.0},
            {"beam_m": 0.0},
            {"lcg_m": 0.0},
        ],
    )
    def test_invalid_inputs_are_rejected(self, kwargs):
        base = dict(
            weight_n=15000.0, speed_ms=12.0, beam_m=2.5, deadrise_deg=20.0, lcg_m=4.5
        )
        base.update(kwargs)
        with pytest.raises(PlaningPredictionError):
            equilibrium_trim(**base)


class TestValidityBox:
    """Finding F4: the documented validity box is now reported, not implied."""

    def test_low_froude_number_is_flagged(self):
        # Fr = 0.83: solvable, but outside the planing box the module claims for itself.
        pred = planing_prediction(8000.0, 9.0, 12.0, 3.0, 20.0)
        assert any("outside the planing box" in w for w in pred.warnings)

    def test_deadrise_outside_the_fitted_band_is_flagged(self):
        pred = planing_prediction(15000.0, 12.0, 12.0, 2.5, 5.0)
        assert any("deadrise" in w for w in pred.warnings)

    def test_length_beam_ratio_outside_the_matrix_is_flagged(self):
        pred = planing_prediction(15000.0, 12.0, 5.0, 2.5, 20.0)
        assert any("L/B=" in w for w in pred.warnings)

    def test_a_case_inside_the_box_raises_no_warning(self):
        # Fr = 1.11, beta = 20 deg, L/B = 4.8, trim = 3.63 deg: inside every band.
        pred = planing_prediction(25000.0, 12.0, 12.0, 2.5, 20.0, lcg_m=3.0)
        assert pred.warnings == ()

    def test_validity_helper_returns_one_message_per_violation(self):
        assert validity_warnings(2.0, 20.0, 12.0, 2.5) == ()
        assert len(validity_warnings(0.5, 5.0, 5.0, 2.5)) == 3


class TestSpeedSweepIsAuditable:
    """Finding F5: a silently shortened sweep is not evidence."""

    def test_failures_can_be_recorded(self):
        speeds = (2.0, 4.0, 6.0, 8.0, 12.0, 16.0)
        failures: list[tuple[float, str]] = []
        results = speed_sweep(
            60000.0, 12.0, 2.5, 20.0, speeds, failures=failures
        )
        assert len(results) + len(failures) == len(speeds)
        assert failures, "an overloaded craft must fail at low speed"
        assert all(isinstance(msg, str) and msg for _, msg in failures)

    def test_sweep_is_still_backwards_compatible(self):
        results = speed_sweep(15000.0, 12.0, 2.5, 20.0, (8.0, 12.0, 16.0))
        assert all(hasattr(r, "total_drag_n") for r in results)


class TestClassicCorrelation:
    """Properties of the classical Cv-dependent correlation used by the closure."""

    def test_speed_coefficient_matches_definition(self):
        assert speed_coefficient(12.0, 2.5) == pytest.approx(12.0 / math.sqrt(G * 2.5))

    def test_lift_falls_with_speed_coefficient(self):
        cv_low = classic_lift_coefficient(4.0, 3.0, 2.0, 20.0)
        cv_high = classic_lift_coefficient(4.0, 3.0, 4.0, 20.0)
        assert cv_low > cv_high

    def test_lift_rises_with_trim_and_wetted_length(self):
        assert classic_lift_coefficient(6.0, 3.0, 2.0, 20.0) > classic_lift_coefficient(
            4.0, 3.0, 2.0, 20.0
        )
        assert classic_lift_coefficient(4.0, 4.0, 2.0, 20.0) > classic_lift_coefficient(
            4.0, 3.0, 2.0, 20.0
        )

    def test_deadrise_correction_reduces_lift(self):
        flat = classic_lift_coefficient(4.0, 3.0, 2.0, 0.0)
        veed = classic_lift_coefficient(4.0, 3.0, 2.0, 20.0)
        assert flat > veed > 0.0

    def test_centre_of_pressure_ratio_is_a_fraction_of_the_wetted_length(self):
        ratio = center_of_pressure_ratio(3.0, 2.0)
        assert 0.0 < ratio < 1.0
        # The *fraction* falls with wetted length while the absolute centre of pressure
        # moves forward, which is the behaviour the moment solve relies on.
        assert center_of_pressure_ratio(4.0, 2.0) < ratio
        assert center_of_pressure_ratio(4.0, 2.0) * 4.0 > ratio * 3.0


class TestCorrelationIncompatibility:
    """Finding F6 (withheld): two different lift correlations are both labelled
    Savitsky (1964) eq. (1). At least one label is wrong; the difference is material."""

    def test_the_two_correlations_differ_by_a_large_margin(self):
        module = lift_coefficient(4.0, 3.0, 20.0)
        classic = classic_lift_coefficient(4.0, 3.0, 2.0, 20.0)
        relative_difference = (module - classic) / classic
        # Pinned: ~-60%. If someone unifies the two formulas this test must fail and
        # the audit document has to be revisited.
        assert -0.65 < relative_difference < -0.53, relative_difference

    def test_no_single_speed_coefficient_reconciles_them(self):
        grid = [
            (tau, lam, beta)
            for tau in (2.0, 4.0, 6.0, 8.0)
            for lam in (1.0, 2.0, 3.0, 4.0)
            for beta in (10.0, 20.0, 30.0)
        ]
        best_worst_case = min(
            max(
                abs(lift_coefficient(t, lam, b) - classic_lift_coefficient(t, lam, cv, b))
                / classic_lift_coefficient(t, lam, cv, b)
                for t, lam, b in grid
            )
            for cv in [0.5 + 0.1 * i for i in range(46)]
        )
        # Even the best-fitting Cv leaves >50% worst-case disagreement.
        assert best_worst_case > 0.5, best_worst_case


class TestModuleSurface:
    """The audit must not have widened the public surface by accident."""

    def test_forward_and_inverse_trim_agree(self):
        weight, speed, beam, beta, lam = 15000.0, 12.0, 2.5, 20.0, 3.0
        tau, lam_out = running_trim(weight, speed, beam, beta, wetted_length_ratio=lam)
        assert lam_out == lam
        assert tau == pytest.approx(
            trim_for_wetted_length(weight, speed, beam, beta, lam)
        )

    def test_prediction_dict_carries_the_provenance(self):
        pred = planing_prediction(15000.0, 12.0, 12.0, 2.5, 20.0, lcg_m=4.5)
        payload = pred.as_dict()
        assert payload["solution_method"] == "classic_equilibrium"
        assert payload["wetted_length_ratio_solved"] is True
        assert isinstance(payload["warnings"], list)