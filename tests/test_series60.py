"""Tests for the Series 60 hull family."""
import math

import pytest

from tensorcad import mesh as meshlib
from tensorcad import series60
from tensorcad.series60 import (
    Series60Parameters, generate_series60, analytic_series60_hydrostatics,
    analytic_volume, analytic_waterplane_area, analytic_maximum_section_area,
    analytic_vcb_above_keel, analytic_cb, analytic_cp, analytic_cm, analytic_cwp,
    _interpolate_exponents, build_series60_mesh,
    series60_analytic_relative_errors, estimate_series60_discretisation_error,
)


class TestExponentInterpolation:
    def test_parent_form_cb60(self):
        n, m = _interpolate_exponents(0.60)
        assert abs(n - 5.0/3.0) < 1e-10
        assert m > 0

    def test_parent_form_cb80(self):
        n, m = _interpolate_exponents(0.80)
        assert n > 0 and m > 0

    def test_intermediate_cb(self):
        n, m = _interpolate_exponents(0.70)
        assert n > 0 and m > 0

    def test_out_of_range_rejected(self):
        with pytest.raises(ValueError):
            _interpolate_exponents(0.50)
        with pytest.raises(ValueError):
            _interpolate_exponents(0.90)


class TestAnalyticFormulas:
    def test_volume_formula(self):
        vol = analytic_volume(10, 1.5, 0.6, 1.6667, 24.0)
        assert vol > 0
        # Cb = n*m/((n+1)(m+1)) ≈ 0.60
        cb = analytic_cb(1.6667, 24.0)
        assert abs(cb - 0.60) < 0.01

    def test_waterplane_area(self):
        awp = analytic_waterplane_area(10, 1.5, 1.6667)
        assert awp > 0

    def test_vcb_above_keel(self):
        vcb = analytic_vcb_above_keel(0.6, 24.0)
        assert 0 < vcb < 0.6


class TestGeneration:
    def test_generate_cb60(self):
        p = Series60Parameters(cb=0.60)
        model = generate_series60(p)
        assert len(model.mesh.vertices) > 0
        assert len(model.mesh.faces) > 0
        assert model.hydrostatics is not None
        assert model.hydrostatics.displacement_volume > 0

    def test_generate_all_parent_forms(self):
        for cb in [0.60, 0.65, 0.70, 0.75, 0.80]:
            p = Series60Parameters(cb=cb)
            model = generate_series60(p)
            assert model.validation.errors == [], f"Cb={cb} has errors: {model.validation.errors}"

    def test_mesh_is_watertight(self):
        p = Series60Parameters(cb=0.60)
        model = generate_series60(p)
        checks = {c.name: c for c in model.validation.checks}
        assert checks["watertight"].passed

    def test_euler_characteristic(self):
        p = Series60Parameters(cb=0.60)
        model = generate_series60(p)
        checks = {c.name: c for c in model.validation.checks}
        assert checks["euler_characteristic_genus_zero"].passed

    def test_outward_orientation(self):
        p = Series60Parameters(cb=0.60)
        model = generate_series60(p)
        checks = {c.name: c for c in model.validation.checks}
        assert checks["outward_orientation"].passed

    def test_dimensions_match(self):
        p = Series60Parameters(length=12, beam=2.0, draft=0.8, cb=0.65)
        model = generate_series60(p)
        checks = {c.name: c for c in model.validation.checks}
        assert checks["length_matches_request"].passed
        assert checks["beam_matches_request"].passed
        assert checks["draft_matches_request"].passed


class TestAnalyticHydrostatics:
    def test_coefficients_consistent(self):
        p = Series60Parameters(cb=0.60)
        ref = analytic_series60_hydrostatics(p)
        # Cb = Cp * Cm
        assert abs(ref["block_coefficient"] - ref["prismatic_coefficient"] * ref["midship_coefficient"]) < 1e-10
        # Cwp = Cp
        assert abs(ref["waterplane_coefficient"] - ref["prismatic_coefficient"]) < 1e-10

    def test_lcb_at_midship(self):
        p = Series60Parameters(cb=0.70)
        ref = analytic_series60_hydrostatics(p)
        assert abs(ref["longitudinal_centre_buoyancy"] - p.length / 2.0) < 1e-10


class TestDiscretisationError:
    def test_convergence_reported(self):
        p = Series60Parameters(cb=0.60, nx=21, nz=7)
        result = estimate_series60_discretisation_error(p, max_triangles=50000)
        assert "reference" in result
        assert "convergence" in result
        assert len(result["convergence"]) >= 1


class TestVerticalTaper:
    """The section profile is ``zeta**(1/m)``; ``zeta**m`` collapses the hull.

    These two cases are the regression for the bug the geometry audit found
    (mission 5283): the lofted mesh integrated 21 times less volume than the
    module's own closed form for every Cb, and the error did not shrink under
    refinement.  The first test pins the physics, the second reproduces the broken
    exponent so the failure mode stays visible.
    """

    @pytest.mark.parametrize("cb", [0.60, 0.70])
    def test_mesh_volume_matches_closed_form(self, cb):
        p = Series60Parameters(length=10.0, beam=1.5, draft=0.6, nx=161, nz=41, cb=cb)
        model = generate_series60(p)
        reference = analytic_series60_hydrostatics(p)
        assert model.hydrostatics.displacement_volume == pytest.approx(
            reference["displacement_volume"], rel=0.02
        )
        assert model.hydrostatics.block_coefficient == pytest.approx(cb, rel=0.02)

    def test_vertical_centre_of_buoyancy_matches_closed_form(self):
        p = Series60Parameters(length=10.0, beam=1.5, draft=0.6, nx=161, nz=41, cb=0.70)
        model = generate_series60(p)
        reference = analytic_series60_hydrostatics(p)
        above_keel = model.hydrostatics.vertical_centre_buoyancy + p.draft
        assert above_keel == pytest.approx(reference["vertical_centre_buoyancy"], rel=0.02)

    def test_inverted_exponent_loses_most_of_the_volume(self):
        original = series60._vertical_taper
        try:
            def broken(parameters, level):  # zeta**m instead of zeta**(1/m)
                if level >= parameters.nz - 1:
                    return 1.0
                if level <= 0:
                    return 0.0
                return math.pow(level / (parameters.nz - 1), parameters.m)

            series60._vertical_taper = broken
            p = Series60Parameters(length=10.0, beam=1.5, draft=0.6, nx=81, nz=21, cb=0.60)
            mesh, _ = build_series60_mesh(p)
            volume = abs(meshlib.total_signed_volume(mesh))
            reference = analytic_volume(10.0, 1.5, 0.6, *p.exponents)
            assert (volume - reference) / reference < -0.9
        finally:
            series60._vertical_taper = original
