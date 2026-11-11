"""Resistance lines, the Prohaska fit, and model-to-ship scaling.

The assertions here are of three kinds, and the difference matters:

* **Definition identities** -- ``CF(1e7) == 0.003`` exactly, ``R_T = 0.5 rho S V^2 CT``,
  ``P_E = R_T V``. A failure is a coding error, not a modelling question.
* **Implementation checks against an independent route** -- the Schoenherr root's own
  residual, the monotonicity of a friction line in ``Re``, the recovery of an injected
  ``(1 + k)`` from synthetic data generated *from* it. These can fail because the code
  is wrong even when the physics is fine.
* **Measured band statements** -- how far the 1957 line sits from the Schoenherr line
  as a function of ``Re``. These pin *measurements*, so a change in agreement is a
  visible event rather than a silent drift; they are stated with the band they were
  measured in, and the low-Reynolds divergence is asserted too, because pretending the
  lines agree everywhere would be a false claim.
"""
from __future__ import annotations

import json
import math

import pytest

from tensorcad.mesh import surface_area, triangle_areas
from tensorcad.resistance import (
    FRICTION_LINES,
    MINIMUM_REYNOLDS,
    ModelShipPrediction,
    ProhaskaFit,
    Resistance,
    ResistanceError,
    bare_hull_resistance,
    friction_coefficient,
    froude_number,
    ittc57_friction_coefficient,
    mesh_bare_hull_resistance,
    prandtl_schlichting_friction_coefficient,
    prohaska_fit,
    reynolds_number,
    schoenherr_friction_coefficient,
    ship_from_model_test,
)
from tensorcad.wigley import (
    analytic_wetted_surface_area,
    generate_wigley,
    reference_parameters,
)

REYNOLDS_GRID = (1.0e4, 1.0e5, 1.0e6, 1.0e7, 1.0e8, 1.0e9)


def _deck_area(mesh) -> float:
    """Triangles lying in the design waterline plane ``z = 0``."""
    areas = triangle_areas(mesh)
    return math.fsum(
        area
        for face, area in zip(mesh.faces, areas)
        if max(abs(mesh.vertices[index][2]) for index in face) <= 1e-12
    )


class TestReynoldsAndFroude:
    def test_reynolds_number_is_the_definition(self):
        assert reynolds_number(4.36, 2.0, 1.0e-6) == pytest.approx(4.36 * 2.0 / 1.0e-6, rel=1e-15)

    def test_froude_number_is_the_definition(self):
        assert froude_number(10.0, 0.3 * math.sqrt(9.80665 * 10.0)) == pytest.approx(0.3, rel=1e-13)

    def test_scaling_law_between_the_two(self):
        # Re and Fn are not independent: at fixed fluid, Re = Fn sqrt(g L) L / nu.
        length, velocity, nu = 6.0, 1.8, 1.2e-6
        fn = froude_number(length, velocity)
        assert reynolds_number(length, velocity, nu) == pytest.approx(
            fn * math.sqrt(9.80665 * length) * length / nu, rel=1e-13
        )

    @pytest.mark.parametrize("bad", [0.0, -1.0, math.inf, math.nan])
    def test_non_physical_inputs_are_rejected(self, bad):
        with pytest.raises(ResistanceError):
            bare_hull_resistance(length=bad, velocity=1.0, wetted_surface=1.0)
        with pytest.raises(ResistanceError):
            bare_hull_resistance(length=1.0, velocity=bad, wetted_surface=1.0)
        with pytest.raises(ResistanceError):
            bare_hull_resistance(length=1.0, velocity=1.0, wetted_surface=bad)


class TestFrictionLines:
    def test_ittc57_is_the_closed_form(self):
        for re in REYNOLDS_GRID:
            assert ittc57_friction_coefficient(re) == pytest.approx(
                0.075 / (math.log10(re) - 2.0) ** 2, rel=1e-15
            )

    @pytest.mark.parametrize(
        ("reynolds", "expected"),
        [(1.0e5, 0.008333333333333333), (1.0e6, 0.0046875), (1.0e7, 0.003),
         (1.0e8, 0.0020833333333333333), (1.0e9, 0.001530612244897959)],
    )
    def test_ittc57_values_at_decade_reynolds(self, reynolds, expected):
        # These are the formula's own values to 16 digits; pinning them catches a
        # coefficient or a log base changing under the numerical tests.
        assert ittc57_friction_coefficient(reynolds) == pytest.approx(expected, rel=1e-15)

    def test_schoenherr_root_satisfies_its_implicit_equation(self):
        for re in REYNOLDS_GRID[1:]:
            cf = schoenherr_friction_coefficient(re)
            residual = 0.242 / math.sqrt(cf) - math.log10(re * cf)
            assert abs(residual) < 1e-12, (re, cf, residual)

    def test_schoenherr_bisection_tolerance_is_honoured(self):
        loose = schoenherr_friction_coefficient(1.0e8, tolerance=1e-6)
        tight = schoenherr_friction_coefficient(1.0e8, tolerance=1e-14)
        assert abs(loose - tight) <= 1e-6 * tight
        assert loose != tight

    @pytest.mark.parametrize("line", sorted(FRICTION_LINES))
    def test_every_line_decreases_with_reynolds(self, line):
        values = [friction_coefficient(re, line) for re in REYNOLDS_GRID]
        assert all(later < earlier for earlier, later in zip(values, values[1:]))

    def test_unknown_line_is_rejected_by_name(self):
        with pytest.raises(ResistanceError, match="unknown friction line"):
            friction_coefficient(1.0e7, "ittc-1957")  # a plausible typo must not pass

    def test_below_the_turbulent_domain_is_rejected(self):
        assert MINIMUM_REYNOLDS == 1.0e4
        with pytest.raises(ResistanceError, match="below the turbulent-line domain"):
            ittc57_friction_coefficient(1.0e3)
        # The 1957 formula itself would happily return a number here; it is the
        # domain guard, not the arithmetic, that refuses.
        assert 0.075 / (math.log10(1.0e3) - 2.0) ** 2 == pytest.approx(0.075, rel=1e-3)


class TestLineAgreement:
    """Measured distance between independent flat-plate lines, as a function of Re."""

    def test_ship_range_agreement_with_schoenherr(self):
        # Re >= 1e7 is the range the 1957 line was fitted for; the claim is 2.3 %.
        worst = max(
            abs(ittc57_friction_coefficient(re) - schoenherr_friction_coefficient(re))
            / schoenherr_friction_coefficient(re)
            for re in (1.0e7, 1.0e8, 1.0e9)
        )
        assert worst < 0.023, worst

    def test_low_reynolds_divergence_is_real_and_asserted(self):
        # Not a bug to be fixed: the lines genuinely part below ~1e6 (16 % at 1e5),
        # and a future change that made them agree there would be a different model.
        deviation = abs(ittc57_friction_coefficient(1.0e5) - schoenherr_friction_coefficient(1.0e5)) / (
            schoenherr_friction_coefficient(1.0e5)
        )
        assert 0.15 < deviation < 0.17, deviation

    def test_prandtl_schlichting_bounds_the_1957_line_in_the_ship_range(self):
        for re in (1.0e7, 1.0e8, 1.0e9):
            ittc = ittc57_friction_coefficient(re)
            ps = prandtl_schlichting_friction_coefficient(re)
            assert abs(ittc - ps) / ps < 0.026, (re, ittc, ps)


class TestProhaskaFit:
    MODEL_LENGTH = 4.36
    NOMINAL_K = (1.15, 0.25)  # (form factor, wave coefficient) injected

    def _synthetic(self, form_factor: float, wave: float, *, length: float | None = None):
        length = self.MODEL_LENGTH if length is None else length
        froudes = (0.12, 0.15, 0.18, 0.21, 0.24, 0.27)
        coefficients = []
        for fn in froudes:
            re = reynolds_number(length, fn * math.sqrt(9.80665 * length), 1.1892e-6)
            coefficients.append(form_factor * ittc57_friction_coefficient(re) + wave * fn**4)
        return froudes, tuple(coefficients)

    def test_exact_recovery_of_the_injected_form_factor(self):
        form_factor, wave = self.NOMINAL_K
        froudes, totals = self._synthetic(form_factor, wave)
        fit = prohaska_fit(froudes, totals, length=self.MODEL_LENGTH)
        assert isinstance(fit, ProhaskaFit)
        assert fit.form_factor == pytest.approx(form_factor, rel=1e-12)
        assert fit.wave_coefficient == pytest.approx(wave, rel=1e-12)
        assert fit.max_relative_residual < 1e-14

    def test_a_wrong_form_factor_is_detected_not_absorbed(self):
        # If the fit were merely "close enough", a 20 % injection would still pass.
        froudes, totals = self._synthetic(1.05, 0.10)
        fit = prohaska_fit(froudes, totals, length=self.MODEL_LENGTH)
        assert fit.form_factor == pytest.approx(1.05, rel=1e-12)
        assert abs(fit.form_factor - 1.15) > 0.09

    def test_explicit_friction_coefficients_give_the_same_answers(self):
        froudes, totals = self._synthetic(*self.NOMINAL_K)
        cf = tuple(
            ittc57_friction_coefficient(
                reynolds_number(
                    self.MODEL_LENGTH, fn * math.sqrt(9.80665 * self.MODEL_LENGTH), 1.1892e-6
                )
            )
            for fn in froudes
        )
        with_length = prohaska_fit(froudes, totals, length=self.MODEL_LENGTH)
        explicit = prohaska_fit(froudes, totals, friction_coefficients=cf)
        assert explicit.form_factor == pytest.approx(with_length.form_factor, rel=1e-14)

    def test_scatter_is_reported_rather_than_hidden(self):
        froudes, totals = self._synthetic(*self.NOMINAL_K)
        noisy = tuple(ct * (1.0 + 0.01 * (index % 2)) for index, ct in enumerate(totals))
        fit = prohaska_fit(froudes, noisy, length=self.MODEL_LENGTH)
        assert fit.max_relative_residual > 1e-3
        assert fit.intercept_standard_error > 0.0
        assert fit.slope_standard_error > 0.0

    def test_fitted_line_reproduces_the_data_it_was_fitted_to(self):
        froudes, totals = self._synthetic(*self.NOMINAL_K)
        fit = prohaska_fit(froudes, totals, length=self.MODEL_LENGTH)
        for fn, ct, cf in zip(froudes, totals, fit.friction_coefficients):
            re = reynolds_number(self.MODEL_LENGTH, fn * math.sqrt(9.80665 * self.MODEL_LENGTH))
            assert fit.total_coefficient(fn, re) == pytest.approx(ct, rel=1e-12)

    def test_two_points_are_refused(self):
        froudes, totals = self._synthetic(*self.NOMINAL_K)
        with pytest.raises(ResistanceError, match="at least 3 points"):
            prohaska_fit(froudes[:2], totals[:2], length=self.MODEL_LENGTH)

    def test_length_and_friction_coefficients_cannot_both_be_given(self):
        froudes, totals = self._synthetic(*self.NOMINAL_K)
        cf = tuple(ittc57_friction_coefficient(1.0e7) for _ in froudes)
        with pytest.raises(ResistanceError, match="not both"):
            prohaska_fit(froudes, totals, length=self.MODEL_LENGTH, friction_coefficients=cf)

    def test_neither_route_given_is_refused(self):
        froudes, totals = self._synthetic(*self.NOMINAL_K)
        with pytest.raises(ResistanceError, match="needs friction_coefficients or a model length"):
            prohaska_fit(froudes, totals)

    def test_mismatched_lengths_are_refused(self):
        froudes, totals = self._synthetic(*self.NOMINAL_K)
        with pytest.raises(ResistanceError, match="froude_numbers has"):
            prohaska_fit(froudes, totals[:-1], length=self.MODEL_LENGTH)

    def test_negative_coefficients_are_refused(self):
        froudes, totals = self._synthetic(*self.NOMINAL_K)
        with pytest.raises(ResistanceError, match="must be positive"):
            prohaska_fit(froudes, (-1.0,) + totals[1:], length=self.MODEL_LENGTH)


class TestResistanceComponents:
    def _case(self, **overrides):
        arguments = dict(
            length=10.0,
            velocity=1.5,
            wetted_surface=14.879,
            form_factor=1.15,
            residual_coefficient=4.0e-4,
            correlation_allowance=0.0,
        )
        arguments.update(overrides)
        return bare_hull_resistance(**arguments)

    def test_coefficients_sum_to_the_total(self):
        body = self._case()
        assert body.total_coefficient == pytest.approx(
            body.viscous_coefficient + body.residual_coefficient + body.correlation_allowance,
            rel=1e-15,
        )

    def test_every_resistance_is_its_coefficient_times_dynamic_pressure_area(self):
        body = self._case()
        q = 0.5 * body.density * body.wetted_surface * body.velocity**2
        assert body.dynamic_pressure_area == pytest.approx(q, rel=1e-15)
        assert body.friction_resistance == pytest.approx(q * body.friction_coefficient, rel=1e-15)
        assert body.residual_resistance == pytest.approx(q * body.residual_coefficient, rel=1e-15)
        assert body.total_resistance == pytest.approx(q * body.total_coefficient, rel=1e-15)

    def test_form_resistance_is_the_form_factor_increment(self):
        body = self._case()
        assert body.form_resistance == pytest.approx(
            (body.form_factor - 1.0) * body.friction_resistance, rel=1e-15
        )
        assert body.viscous_resistance == pytest.approx(
            body.friction_resistance + body.form_resistance, rel=1e-15
        )

    def test_effective_power_is_resistance_times_speed(self):
        body = self._case()
        assert body.effective_power == pytest.approx(
            body.total_resistance * body.velocity, rel=1e-15
        )

    def test_flat_plate_defaults_mean_not_measured(self):
        body = bare_hull_resistance(length=10.0, velocity=1.5, wetted_surface=14.879)
        assert body.form_factor == 1.0
        assert body.residual_coefficient == 0.0
        assert body.total_resistance == pytest.approx(body.friction_resistance, rel=1e-15)

    def test_form_factor_below_one_is_refused(self):
        with pytest.raises(ResistanceError, match="form_factor must be >= 1"):
            bare_hull_resistance(length=1.0, velocity=1.0, wetted_surface=1.0, form_factor=0.9)

    def test_negative_residual_is_refused(self):
        with pytest.raises(ResistanceError, match="must be >= 0"):
            bare_hull_resistance(length=1.0, velocity=1.0, wetted_surface=1.0, residual_coefficient=-1e-4)

    def test_resistance_grows_monotonically_with_speed(self):
        speeds = (0.5, 1.0, 1.5, 2.0, 3.0)
        totals = [self._case(velocity=v).total_resistance for v in speeds]
        assert all(later > earlier for earlier, later in zip(totals, totals[1:]))

    def test_doubling_length_at_fixed_speed_lowers_the_coefficient(self):
        # Reynolds effect: the longer body is relatively cheaper per unit area.
        short = self._case(length=5.0)
        long = self._case(length=50.0)
        assert long.friction_coefficient < short.friction_coefficient
        assert long.total_coefficient < short.total_coefficient

    def test_density_scales_force_linearly(self):
        light = self._case(density=999.1)
        heavy = self._case(density=1998.2)
        assert heavy.total_resistance == pytest.approx(2.0 * light.total_resistance, rel=1e-14)

    def test_record_serialises(self):
        body = self._case()
        record = body.to_dict()
        assert json.loads(json.dumps(record))["total_resistance_n"] == pytest.approx(
            body.total_resistance, rel=1e-15
        )
        assert isinstance(body, Resistance)


class TestMeshCoupling:
    """The wetted area that feeds the friction must come off the geometry."""

    @pytest.mark.parametrize(("nx", "nz"), [(41, 13), (81, 21), (161, 41)])
    def test_mesh_wetted_area_converges_to_the_analytic_flank(self, nx, nz):
        parameters = reference_parameters()
        model = generate_wigley(parameters, nx=nx, nz=nz)
        shell = surface_area(model.mesh)
        wetted = shell - _deck_area(model.mesh)
        reference = analytic_wetted_surface_area(
            parameters.length,
            parameters.beam,
            parameters.draft,
            parameters.longitudinal_exponent,
            parameters.vertical_exponent,
        )
        assert abs(wetted - reference) / reference < 2.5e-3, (nx, nz, wetted, reference)

    def test_deck_patch_is_a_large_part_of_the_shell(self):
        # The documented trap: shell area is not wetted area on a closed hull.
        parameters = reference_parameters()
        model = generate_wigley(parameters)
        shell = surface_area(model.mesh)
        deck = _deck_area(model.mesh)
        assert 0.30 < deck / shell < 0.32, deck / shell

    def test_mesh_entry_point_uses_the_shell_area_it_documents(self):
        parameters = reference_parameters()
        model = generate_wigley(parameters)
        body = mesh_bare_hull_resistance(model.mesh, length=parameters.length, velocity=1.5)
        assert body.wetted_surface == pytest.approx(surface_area(model.mesh), rel=1e-15)
        assert body.froude == pytest.approx(
            froude_number(parameters.length, 1.5), rel=1e-15
        )


class TestModelToShip:
    MODEL_LENGTH = 4.36
    SHIP_LENGTH = 174.4
    MODEL_VELOCITY = 2.0
    MODEL_SURFACE = 5.9
    FORM_FACTOR = 1.15

    def _prediction(self, **overrides):
        arguments = dict(
            model_length=self.MODEL_LENGTH,
            ship_length=self.SHIP_LENGTH,
            model_velocity=self.MODEL_VELOCITY,
            model_wetted_surface=self.MODEL_SURFACE,
            ship_wetted_surface=self.MODEL_SURFACE * (self.SHIP_LENGTH / self.MODEL_LENGTH) ** 2,
            form_factor=self.FORM_FACTOR,
            model_total_coefficient=4.2e-3,
        )
        arguments.update(overrides)
        return ship_from_model_test(**arguments)

    def test_velocity_scales_as_the_square_root_of_length(self):
        prediction = self._prediction()
        assert prediction.velocity_scale == pytest.approx(
            math.sqrt(self.SHIP_LENGTH / self.MODEL_LENGTH), rel=1e-15
        )
        assert prediction.ship.velocity == pytest.approx(
            self.MODEL_VELOCITY * prediction.velocity_scale, rel=1e-15
        )
        assert prediction.ship.froude == pytest.approx(prediction.model.froude, rel=1e-13)

    def test_residual_coefficient_is_recovered_from_the_model_total(self):
        prediction = self._prediction()
        expected = 4.2e-3 - self.FORM_FACTOR * prediction.model.friction_coefficient
        assert prediction.residual_coefficient == pytest.approx(expected, rel=1e-13)
        assert prediction.ship.total_coefficient == pytest.approx(
            self.FORM_FACTOR * prediction.ship.friction_coefficient + expected, rel=1e-13
        )

    def test_round_trip_through_the_scaling_recovers_the_model_total(self):
        prediction = self._prediction()
        rebuilt = (
            self.FORM_FACTOR * prediction.model.friction_coefficient
            + prediction.residual_coefficient
        )
        assert rebuilt == pytest.approx(prediction.model_total_coefficient, rel=1e-15)

    def test_resistance_route_and_coefficient_route_agree(self):
        prediction = self._prediction()
        area = prediction.model.dynamic_pressure_area
        resistance_route = self._prediction(model_total_coefficient=None, model_total_resistance=4.2e-3 * area)
        assert resistance_route.residual_coefficient == pytest.approx(
            prediction.residual_coefficient, rel=1e-13
        )

    def test_inconsistent_coefficient_and_resistance_are_refused(self):
        with pytest.raises(ResistanceError, match="disagree"):
            self._prediction(model_total_resistance=1.0)

    def test_explicit_residual_coefficient_supersedes_the_model_total(self):
        explicit = self._prediction(model_total_coefficient=None, residual_coefficient=4.0e-4)
        assert explicit.residual_coefficient == pytest.approx(4.0e-4, rel=1e-15)
        assert explicit.ship.residual_coefficient == pytest.approx(4.0e-4, rel=1e-15)

    def test_allowance_is_applied_to_the_ship_only(self):
        base = self._prediction()
        with_allowance = self._prediction(correlation_allowance=8.0e-4)
        assert with_allowance.model.total_coefficient == pytest.approx(
            base.model.total_coefficient, rel=1e-15
        )
        assert with_allowance.ship.total_coefficient - base.ship.total_coefficient == pytest.approx(
            8.0e-4, rel=1e-10
        )
        assert with_allowance.ship.total_resistance > base.ship.total_resistance

    def test_ship_is_not_a_lambda_cubed_copy_of_the_model(self):
        prediction = self._prediction()
        # CT falls with Re, so R_s / R_m sits below the displacement ratio (4x10^3 here
        # would be the naive geometric answer; the coefficient ratio pulls it down).
        naive = (self.SHIP_LENGTH / self.MODEL_LENGTH) ** 3
        assert prediction.ship_resistance_ratio < naive
        assert prediction.scale_effect > 1.0
        assert prediction.reynolds_scale > 1.0
        assert prediction.reynolds_scale == pytest.approx(
            self.SHIP_LENGTH / self.MODEL_LENGTH * prediction.velocity_scale, rel=1e-13
        )

    def test_no_input_route_is_refused(self):
        with pytest.raises(ResistanceError, match="give one of"):
            ship_from_model_test(
                model_length=4.36,
                ship_length=174.4,
                model_velocity=2.0,
                model_wetted_surface=5.9,
                ship_wetted_surface=5.9 * 1600.0,
                form_factor=1.15,
            )

    def test_superseding_and_model_total_together_is_refused(self):
        with pytest.raises(ResistanceError, match="supersedes"):
            self._prediction(residual_coefficient=4.0e-4)

    def test_record_serialises(self):
        prediction = self._prediction()
        record = prediction.to_dict()
        assert json.loads(json.dumps(record))["ship"]["total_resistance_n"] == pytest.approx(
            prediction.ship.total_resistance, rel=1e-14
        )
        assert isinstance(prediction, ModelShipPrediction)