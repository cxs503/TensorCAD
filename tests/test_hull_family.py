"""The unified hull-family entry: ``family="series60"`` through :mod:`tensorcad.hull`.

The Series 60 station law lives in :mod:`tensorcad.series60`; what is tested here is
the *dispatch*: the same :class:`~tensorcad.hull.HullParameters` object can ask for
either family, the delegated surface is bit-for-bit the one the family module
builds, the offsets keep this module's full-girth convention, and the parameters a
family cannot honour are refused instead of silently dropped.
"""
from __future__ import annotations

import math

import pytest

from tensorcad.mesh import is_edge_manifold, is_watertight
from tensorcad.hull import (
    FAMILIES,
    HullParameters,
    analytic_relative_errors,
    estimate_discretisation_error,
    generate_hull,
    resolve_resolution,
)
from tensorcad.series60 import (
    Series60Parameters,
    analytic_series60_hydrostatics,
    generate_series60,
)

PARENTS = (0.60, 0.65, 0.70, 0.75, 0.80)
MAIN = dict(length=10.0, beam=1.5, draft=0.6, water_density=1000.0)


def series60_parameters(**overrides) -> HullParameters:
    base = dict(MAIN, family="series60", block_coefficient=0.70, nx=41, nz=13)
    base.update(overrides)
    return HullParameters(**base)


def direct_parameters(**overrides) -> Series60Parameters:
    base = dict(MAIN, cb=0.70, nx=41, nz=13)
    base.update(overrides)
    return Series60Parameters(**base)


class TestLegacyStaysDefault:
    def test_default_family_is_the_shipped_sine_taper(self):
        parameters = HullParameters()
        assert parameters.family == "sine"
        assert FAMILIES == ("sine", "series60")

    def test_sine_family_regenerates_the_published_volume(self):
        # docs/parametric_hull: 7.446554 m^3 at 41x13 against the analytic
        # 7.512698 m^3. The family dispatch must not touch this path.
        model = generate_hull(HullParameters())
        assert model.hydrostatics is not None
        assert model.hydrostatics.displacement_volume == pytest.approx(7.446554, rel=1e-6)
        assert model.bulb_mesh is None
        assert model.family_report == {}
        assert model.metadata()["family"] == "sine"


class TestSeries60Dispatch:
    def test_mesh_is_identical_to_the_family_module(self):
        unified = generate_hull(series60_parameters(nx=81, nz=21))
        direct = generate_series60(direct_parameters(nx=81, nz=21))
        assert unified.mesh.vertices == direct.mesh.vertices
        assert unified.mesh.faces == direct.mesh.faces
        assert unified.hydrostatics is not None and direct.hydrostatics is not None
        assert unified.hydrostatics.displacement_volume == pytest.approx(
            direct.hydrostatics.displacement_volume, rel=0.0, abs=0.0
        )

    def test_offsets_keep_the_full_girth_convention(self):
        unified = generate_hull(series60_parameters(nx=81, nz=21))
        direct = generate_series60(direct_parameters(nx=81, nz=21))
        widths = [width for _x, width in unified.offsets]
        assert max(widths) == pytest.approx(MAIN["beam"], rel=1e-12)
        for (x_u, width), (x_d, half) in zip(unified.offsets, direct.offsets):
            assert x_u == pytest.approx(x_d, rel=0.0, abs=0.0)
            assert width == pytest.approx(2.0 * half, rel=1e-15)

    def test_five_parent_forms_track_their_closed_forms(self):
        previous = -math.inf
        for cb in PARENTS:
            parameters = series60_parameters(block_coefficient=cb, nx=161, nz=41)
            model = generate_hull(parameters)
            reference = analytic_series60_hydrostatics(
                direct_parameters(cb=cb, nx=161, nz=41)
            )["displacement_volume"]
            measured = model.hydrostatics.displacement_volume
            error = measured / reference - 1.0
            # The loft under-fills the steep section near the keel, at first order
            # in 1/nz; docs/series60.md records -1.0 % to -1.1 % at 161x41.
            assert -0.03 < error < 0.0, f"Cb={cb}: volume error {error:+.4%}"
            assert measured > previous, f"Cb={cb}: displacement must grow with Cb"
            previous = measured

    def test_measured_block_coefficient_matches_the_target(self):
        for cb in PARENTS:
            model = generate_hull(series60_parameters(block_coefficient=cb, nx=161, nz=41))
            measured = model.hydrostatics.block_coefficient
            assert measured == pytest.approx(cb, rel=0.03), f"Cb={cb}"

    def test_refinement_converges_and_keeps_the_family(self):
        coarse = generate_hull(series60_parameters(nx=41, nz=13))
        fine = generate_hull(series60_parameters(nx=161, nz=41))
        reference = analytic_series60_hydrostatics(direct_parameters())["displacement_volume"]
        coarse_error = abs(coarse.hydrostatics.displacement_volume / reference - 1.0)
        fine_error = abs(fine.hydrostatics.displacement_volume / reference - 1.0)
        assert fine_error < coarse_error
        estimate = estimate_discretisation_error(series60_parameters())
        assert estimate["status"] == "measured"
        assert estimate["refined_grid"] == {"nx": 81, "nz": 25}

    def test_relative_error_table_is_signed_and_small(self):
        report = analytic_relative_errors(series60_parameters(nx=161, nz=41))
        assert report["family"] == "series60"
        assert report["relative_error"]["displacement_volume"] < 0.0
        assert report["max_abs_relative_error"] < 0.03
        assert "Todd" in report["reference_basis"]

    def test_resolve_resolution_uses_the_family_closed_form(self):
        plan = resolve_resolution(series60_parameters(nx=41, nz=13), target_volume_error=0.02)
        assert plan.accepted
        assert plan.relative_error < 0.0
        assert plan.parameters.family == "series60"

    def test_metadata_carries_the_family_report(self):
        model = generate_hull(series60_parameters())
        metadata = model.metadata()
        assert metadata["family"] == "series60"
        assert metadata["family_report"]["station_law"]["block_coefficient"] == pytest.approx(0.70)
        assert metadata["family_report"]["exponents"]["n"] > 0.0


class TestTransomAndBulb:
    def test_transom_sets_the_aft_deck_breadth_and_stays_watertight(self):
        parameters = series60_parameters(transom_ratio=0.15, nx=81, nz=21)
        model = generate_hull(parameters)
        assert model.validation.errors == []
        assert model.offsets[0][1] == pytest.approx(0.15 * MAIN["beam"], rel=1e-9)
        assert is_watertight(model.mesh)
        assert is_edge_manifold(model.mesh)
        transom = model.metadata()["family_report"]["transom"]
        assert transom["transom_ratio"] == pytest.approx(0.15)

    def test_bulb_is_a_separate_solid_and_leaves_the_hull_displacement_alone(self):
        plain = generate_hull(series60_parameters(nx=81, nz=21))
        bulbed = generate_hull(
            series60_parameters(
                nx=81, nz=21, bulb_length_ratio=0.03, bulb_area_ratio=0.08
            )
        )
        assert plain.bulb_mesh is None
        assert bulbed.bulb_mesh is not None
        assert bulbed.hydrostatics.displacement_volume == pytest.approx(
            plain.hydrostatics.displacement_volume, rel=0.0, abs=0.0
        )
        bulb = bulbed.metadata()["family_report"]["appendages"]["bulb"]
        assert bulb["volume_m3"] > 0.0
        assert bulb["fp_area_m2"] > 0.0
        assert bulbed.metadata()["appendages"]["bulb_mesh"]["vertex_count"] > 0


class TestRefusalsInsteadOfSilence:
    def test_series60_without_a_target_cb_is_refused(self):
        with pytest.raises(ValueError, match="block_coefficient"):
            series60_parameters(block_coefficient=None).validate()

    def test_series60_knobs_are_refused_on_the_sine_family(self):
        with pytest.raises(ValueError, match="only apply to family='series60'"):
            HullParameters(block_coefficient=0.70).validate()
        with pytest.raises(ValueError, match="only apply to family='series60'"):
            HullParameters(transom_ratio=0.2).validate()

    def test_sine_exponents_are_refused_on_the_series60_family(self):
        with pytest.raises(ValueError, match="sine-family parameters"):
            series60_parameters(fore_exponent=1.5).validate()

    def test_unknown_family_is_named(self):
        with pytest.raises(ValueError, match="family must be one of"):
            HullParameters(family="wigley").validate()

    def test_cb_outside_the_published_ladder_is_refused(self):
        with pytest.raises(ValueError, match="cb must lie in"):
            series60_parameters(block_coefficient=0.95).validate()