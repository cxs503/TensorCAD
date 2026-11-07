"""Tests for the constant-deadrise planing hull family."""
import math

import pytest

from tensorcad import mesh as meshlib
from tensorcad import planing
from tensorcad.planing import (
    PlaningHullParameters,
    analytic_longitudinal_centre_buoyancy,
    analytic_longitudinal_centre_flotation,
    analytic_planing_hydrostatics,
    analytic_vertical_centre_buoyancy,
    analytic_volume,
    analytic_waterplane_area,
    analytic_wetted_area,
    build_planing_mesh,
    chine_half_beam,
    estimate_planing_discretisation_error,
    generate_planing_hull,
    planing_relative_errors,
    section_area,
    simpson_section_integrals,
    waterline_half_breadth,
    wetted_length,
)

#: A loaded deep-V: chine immersed all the way to the transom.
CHINE_WET = PlaningHullParameters(nx=81, nz=17)
#: A light hull: the chine leaves the waterline 72 % aft of the bow.
CHINE_DRY = PlaningHullParameters(draft=0.40, nx=81, nz=17)


def _euler(mesh):
    vertices = len({index for face in mesh.faces for index in face})
    edges = {tuple(sorted((face[i], face[(i + 1) % 3]))) for face in mesh.faces for i in range(3)}
    return vertices - len(edges) + len(mesh.faces)


class TestParameters:
    def test_defaults_are_consistent(self):
        p = PlaningHullParameters()
        assert 0 < p.deadrise_deg < 60
        assert p.chine_height_transom == pytest.approx(1.5 * p.tan_beta)

    def test_chine_regime_flag(self):
        assert CHINE_WET.chine_immersed_at_transom is True
        assert CHINE_DRY.chine_immersed_at_transom is False

    @pytest.mark.parametrize("kwargs", [
        {"deadrise_deg": 0.0},
        {"deadrise_deg": 60.0},
        {"planform_exponent": 0.0},
        {"planform_exponent": 5.0},
        {"length": -1.0},
        {"draft": 0.0},
        {"nx": 2},
        {"nz": 2},
        {"water_density": 10.0},
    ])
    def test_invalid_parameters_rejected(self, kwargs):
        base = PlaningHullParameters(**kwargs)
        errors, _ = base.problems()
        assert errors
        with pytest.raises(ValueError):
            base.validate()


class TestChineGeometry:
    def test_chine_half_beam_endpoints(self):
        p = PlaningHullParameters()
        assert chine_half_beam(p, 0.0) == 0.0
        assert chine_half_beam(p, p.length) == pytest.approx(p.beam / 2)
        assert chine_half_beam(p, 2 * p.length) == pytest.approx(p.beam / 2)

    def test_chine_half_beam_is_monotone(self):
        p = PlaningHullParameters(planform_exponent=0.6)
        values = [chine_half_beam(p, p.length * i / 20) for i in range(21)]
        assert all(a <= b + 1e-12 for a, b in zip(values, values[1:]))

    def test_wetted_length_clipped_to_hull(self):
        assert wetted_length(CHINE_WET) == pytest.approx(CHINE_WET.length)
        x0 = wetted_length(CHINE_DRY)
        assert 0.0 < x0 < CHINE_DRY.length
        # b(x0) tan(beta) == T is the definition of the crossing
        assert chine_half_beam(CHINE_DRY, x0) * CHINE_DRY.tan_beta == pytest.approx(
            CHINE_DRY.draft, rel=1e-9
        )

    def test_waterline_half_breadth_is_capped_by_the_chine(self):
        p = CHINE_DRY
        assert waterline_half_breadth(p, p.length) == pytest.approx(p.draft / p.tan_beta)
        assert waterline_half_breadth(p, 1.0) == pytest.approx(chine_half_beam(p, 1.0))
        assert waterline_half_breadth(CHINE_WET, CHINE_WET.length) == pytest.approx(
            CHINE_WET.beam / 2
        )


class TestSectionArea:
    def test_pure_v_section_is_draft_squared_over_tan(self):
        p = PlaningHullParameters(draft=0.2)  # chine dry everywhere
        assert section_area(p, p.length) == pytest.approx(p.draft**2 / p.tan_beta)

    def test_composite_section(self):
        p = CHINE_WET
        b = chine_half_beam(p, p.length)
        ref = b * b * p.tan_beta + 2.0 * b * (p.draft - b * p.tan_beta)
        assert section_area(p, p.length) == pytest.approx(ref)

    def test_area_is_continuous_across_the_chine_crossing(self):
        p = CHINE_DRY
        x0 = wetted_length(p)
        left = section_area(p, x0 - 1e-6)
        right = section_area(p, x0 + 1e-6)
        assert left == pytest.approx(right, rel=1e-6)


class TestClosedFormAnchors:
    def test_everywhere_immersed_chine_volume(self):
        """Draft above the transom chine height: every section is V plus a wall.

        With ``b = c x**q`` the whole hull is in the composite branch, so the closed forms
        reduce to ``V = 2T c L**(q+1)/(q+1) - tan(beta) c**2 L**(2q+1)/(2q+1)`` and
        ``Awp = 2 c L**(q+1)/(q+1)`` -- i.e. ``Cwp = 1/(q+1)`` on the nose.
        """
        p = CHINE_WET
        assert p.draft > p.chine_height_transom
        assert wetted_length(p) == pytest.approx(p.length)
        c = p.half_beam_transom / p.length**p.planform_exponent
        expected = (2.0 * p.draft * c * p.length ** (p.planform_exponent + 1)
                    / (p.planform_exponent + 1)
                    - p.tan_beta * c * c * p.length ** (2 * p.planform_exponent + 1)
                    / (2 * p.planform_exponent + 1))
        assert analytic_volume(p) == pytest.approx(expected, rel=1e-12)
        assert analytic_waterplane_area(p) == pytest.approx(
            2.0 * c * p.length ** (p.planform_exponent + 1) / (p.planform_exponent + 1),
            rel=1e-12,
        )

    def test_pure_v_region_section_and_centres(self):
        """Aft of the chine crossing the waterline cuts the V: A = T**2/tan(beta)."""
        p = CHINE_DRY
        x0 = wetted_length(p)
        assert section_area(p, p.length) == pytest.approx(p.draft**2 / p.tan_beta)
        assert waterline_half_breadth(p, p.length) == pytest.approx(p.draft / p.tan_beta)
        # the transom is a plain triangle, so its contribution to the wetted area is
        # the same T**2/tan(beta) and the bottom is sec(beta) * Awp
        ref = analytic_planing_hydrostatics(p)
        assert ref["wetted_surface_area"] > ref["waterplane_area"] * math.hypot(1.0, p.tan_beta)
        assert x0 < p.length

    def test_straight_chine_prism_waterplane_coefficient(self):
        """q = 1 gives b = (B/2) x/L, so Cwp = 2/(q+1) / 2 = 1/(q+1) = 1/2 exactly."""
        p = PlaningHullParameters(planform_exponent=1.0, draft=0.6)
        ref = analytic_planing_hydrostatics(p)
        assert ref["waterplane_coefficient"] == pytest.approx(0.5, rel=1e-12)

    def test_straight_chine_prism_volume(self):
        """T B L/2 - tan(beta) B**2 L/12 for an everywhere-immersed chine and q = 1."""
        p = PlaningHullParameters(planform_exponent=1.0, draft=0.6)
        expected = p.draft * p.beam * p.length / 2.0 - p.tan_beta * p.beam**2 * p.length / 12.0
        assert analytic_volume(p) == pytest.approx(expected, rel=1e-12)


class TestAnalyticConsistency:
    def test_simpson_matches_closed_form(self):
        """Two analytic routes: paper integrals against quadrature of the sections."""
        for p in (CHINE_WET, CHINE_DRY):
            ref = analytic_planing_hydrostatics(p)
            quad = simpson_section_integrals(p, intervals=2048)
            assert quad["displacement_volume"] == pytest.approx(ref["displacement_volume"], rel=1e-6)
            assert quad["longitudinal_centre_buoyancy"] == pytest.approx(
                ref["longitudinal_centre_buoyancy"], rel=1e-5
            )
            assert quad["vertical_centre_buoyancy"] == pytest.approx(
                ref["vertical_centre_buoyancy"], rel=1e-5
            )

    def test_centres_are_inside_the_hull(self):
        for p in (CHINE_WET, CHINE_DRY):
            ref = analytic_planing_hydrostatics(p)
            assert 0.0 < ref["longitudinal_centre_buoyancy"] < p.length
            assert 0.0 < ref["longitudinal_centre_flotation"] < p.length
            assert 0.0 < ref["vertical_centre_buoyancy_above_keel"] < p.draft

    def test_coefficients_are_positive_and_ordered(self):
        ref = analytic_planing_hydrostatics(CHINE_WET)
        assert 0.0 < ref["block_coefficient"] < ref["prismatic_coefficient"] <= 1.0000001
        assert 0.0 < ref["waterplane_coefficient"] < 1.0


class TestMesh:
    def test_mesh_is_watertight_and_outward(self):
        mesh = build_planing_mesh(CHINE_WET)
        assert meshlib.is_watertight(mesh) is True
        assert meshlib.is_edge_manifold(mesh) is True
        assert meshlib.zero_area_faces(mesh) == ()
        assert _euler(mesh) == 2
        assert meshlib.total_signed_volume(mesh) > 0.0

    def test_mesh_is_watertight_in_the_chine_dry_regime(self):
        mesh = build_planing_mesh(CHINE_DRY)
        assert meshlib.is_watertight(mesh) is True
        assert _euler(mesh) == 2

    def test_mesh_bounding_box(self):
        p = CHINE_WET
        box = build_planing_mesh(p).bounding_box()
        assert box.minimum[0] == pytest.approx(0.0)
        assert box.maximum[0] == pytest.approx(p.length)
        assert box.maximum[1] == pytest.approx(p.beam / 2)
        assert box.minimum[1] == pytest.approx(-p.beam / 2)
        assert box.minimum[2] == pytest.approx(-p.draft)
        assert box.maximum[2] == pytest.approx(0.0)

    def test_waterline_is_a_mesh_level(self):
        """No face may straddle z = 0, or the waterplane integral is not exact."""
        p = CHINE_WET
        mesh = build_planing_mesh(p)
        above = 0
        for face in mesh.faces:
            zs = [mesh.vertices[i][2] for i in face]
            assert not (min(zs) < -1e-12 and max(zs) > 1e-12)
            if min(zs) >= -1e-12:
                above += 1
        assert above > 0  # the waterline cap

    def test_grid_refinement_of_the_mesh(self):
        coarse = len(build_planing_mesh(PlaningHullParameters(nx=21, nz=9)).faces)
        fine = len(build_planing_mesh(PlaningHullParameters(nx=41, nz=17)).faces)
        assert fine > coarse


class TestMeshAgainstClosedForm:
    @pytest.mark.parametrize("p", [CHINE_WET, CHINE_DRY])
    def test_hydrostatics_match_closed_form(self, p):
        model = generate_planing_hull(p)
        errors = planing_relative_errors(p, model)
        for name, tolerance in (
            ("displacement_volume", 5e-3),
            ("waterplane_area", 5e-3),
            ("wetted_surface_area", 5e-3),
            ("longitudinal_centre_buoyancy", 5e-3),
            ("vertical_centre_buoyancy", 5e-3),
            ("longitudinal_centre_flotation", 5e-3),
        ):
            assert abs(errors[name]) < tolerance, (name, errors[name])

    def test_transom_area_matches_the_analytic_section(self):
        p = CHINE_WET
        model = generate_planing_hull(p)
        assert model.hydrostatics is not None
        assert model.hydrostatics.maximum_section_area == pytest.approx(
            section_area(p, p.length), rel=5e-3
        )

    def test_prismatic_coefficient_matches(self):
        p = CHINE_WET
        model = generate_planing_hull(p)
        assert model.hydrostatics is not None
        ref = analytic_planing_hydrostatics(p)["prismatic_coefficient"]
        assert model.hydrostatics.prismatic_coefficient == pytest.approx(ref, rel=5e-3)

    def test_wetted_area_exceeds_bottom_area(self):
        """The wetted side and the transom must add to the bottom, not replace it."""
        p = CHINE_DRY
        ref = analytic_planing_hydrostatics(p)
        bottom = ref["waterplane_area"] * math.hypot(1.0, p.tan_beta)
        assert ref["wetted_surface_area"] > bottom + ref["transom_section_area"] - 1e-12

    def test_convergence_is_monotone(self):
        study = estimate_planing_discretisation_error(
            PlaningHullParameters(nx=21, nz=9), max_triangles=20_000
        )
        errors = [abs(row["volume_error"]) for row in study["convergence"]]
        assert len(errors) >= 3
        assert all(b < a for a, b in zip(errors, errors[1:]))
        assert errors[-1] < 1e-3

    def test_deterministic(self):
        a = generate_planing_hull(CHINE_WET)
        b = generate_planing_hull(CHINE_WET)
        assert a.hydrostatics is not None and b.hydrostatics is not None
        assert a.hydrostatics.displacement_volume == b.hydrostatics.displacement_volume


class TestModelSurface:
    def test_generate_with_overrides(self):
        model = generate_planing_hull(CHINE_WET, draft=0.40)
        assert model.parameters.draft == 0.40
        assert wetted_length(model.parameters) < model.parameters.length

    def test_metadata_contract(self):
        model = generate_planing_hull(CHINE_WET)
        meta = model.metadata()
        assert meta["model_type"] == "constant_deadrise_planing_hull"
        assert meta["parameters"]["deadrise_deg"] == CHINE_WET.deadrise_deg
        assert "hydrostatics" in meta and "analytic_reference" in meta
        assert meta["generator"]["chine_half_beam"] == "(B/2) (x/L)^q"
        assert meta["engineering_status"].startswith("geometry and hydrostatics only")

    def test_validation_reports_no_errors(self):
        for p in (CHINE_WET, CHINE_DRY):
            model = generate_planing_hull(p)
            assert model.validation.errors == []
            names = {check.name for check in model.validation.checks}
            assert {"watertight", "planing_volume", "planing_wetted_area"} <= names

    def test_validation_catches_a_broken_mesh(self):
        """The check that guards the geometry must fail on a deliberately holed hull."""
        p = CHINE_WET
        mesh = build_planing_mesh(p)
        holed = type(mesh)(vertices=mesh.vertices, faces=mesh.faces[1:],
                           units=mesh.units, coordinate_frame=mesh.coordinate_frame)
        assert meshlib.is_watertight(holed) is False