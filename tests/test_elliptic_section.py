"""The elliptic-section family: a closed form that the mesh integrator has to reproduce.

Every assertion here is a comparison between two independent routes -- an analytic
formula (``sin``/``cos``/one square root) and :func:`tensorcad.righting_arm.righting_arm_curve`
solving the waterline on a mesh at constant displacement -- or between a formula and its
own numerical derivative. Nothing is compared against a remembered number.
"""

from __future__ import annotations

import math

import pytest

from tensorcad.elliptic_section import (
    EllipticSectionError,
    elliptic_gm0,
    elliptic_kn,
    elliptic_prism_mesh,
    elliptic_righting_arm,
    elliptic_righting_arm_curve,
    wall_sided_gm0,
)
from tensorcad.righting_arm import righting_arm_curve

#: (a, b) pairs: an ellipse, the circle limit, a narrow-tall and a wide-flat section.
SECTIONS = ((1.0, 1.0), (1.0, 0.6), (0.5, 1.0), (1.0, 0.35))

HEELS = (0.0, 5.0, 10.0, 20.0, 40.0, 60.0, 80.0)


def _signed_volume(mesh) -> float:
    """Volume enclosed by a closed surface, by the divergence theorem."""
    total = 0.0
    for face in mesh.faces:
        a, b, c = (mesh.vertices[index] for index in face)
        total += (
            a[0] * (b[1] * c[2] - b[2] * c[1])
            - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0])
        ) / 6.0
    return total


class TestClosedForm:
    def test_kn_is_zero_upright_and_b_at_ninety(self):
        for a, b in SECTIONS:
            assert elliptic_kn(a, b, 0.0) == pytest.approx(0.0, abs=1e-15)
            assert elliptic_kn(a, b, 90.0) == pytest.approx(b, rel=1e-14)

    def test_circle_limit_is_r_sin_phi(self):
        for radius in (0.25, 1.0, 3.5):
            for heel in HEELS:
                assert elliptic_kn(radius, radius, heel) == pytest.approx(
                    radius * math.sin(math.radians(heel)), rel=1e-14
                )

    def test_kn_is_antisymmetric_in_heel(self):
        for a, b in SECTIONS:
            for heel in HEELS[1:]:
                assert elliptic_kn(a, b, heel) == pytest.approx(-elliptic_kn(a, b, -heel), rel=1e-14)

    def test_gm0_equals_the_numerical_derivative_of_the_analytic_curve(self):
        for a, b in SECTIONS:
            h = 1e-6
            slope = (
                elliptic_righting_arm(a, b, h, kg=0.4)
                - elliptic_righting_arm(a, b, -h, kg=0.4)
            ) / (2.0 * math.radians(h))
            assert elliptic_gm0(a, b, kg=0.4) == pytest.approx(slope, rel=1e-9)

    def test_gm0_reduces_to_radius_minus_kg_for_the_circle(self):
        for radius in (0.25, 1.0, 3.5):
            assert elliptic_gm0(radius, radius, kg=0.4) == pytest.approx(radius - 0.4, rel=1e-14)

    def test_record_peaks_where_the_curve_peaks(self):
        record = elliptic_righting_arm_curve(1.0, 0.6, kg=0.4)
        assert record["peak_gz"] == max(record["gz"])
        assert record["peak_heel_deg"] in record["heels_deg"]
        # the half-ellipse's reflection identity: KN(phi) + KN(180 - phi) = 2 b sin(phi)
        for heel in (5.0, 20.0, 40.0):
            assert elliptic_kn(1.0, 0.6, heel) + elliptic_kn(1.0, 0.6, 180.0 - heel) == pytest.approx(
                2.0 * 0.6 * math.sin(math.radians(heel)), rel=1e-12
            )


class TestAgainstTheMeshIntegrator:
    @pytest.mark.parametrize("a,b", SECTIONS)
    def test_mesh_gz_reproduces_the_closed_form(self, a, b):
        mesh = elliptic_prism_mesh(a, b, 4.0, ring=384)
        curve = righting_arm_curve(mesh, kg=0.4, heels_deg=HEELS, relative_tolerance=1e-11)
        for point in curve.to_dict()["points"]:
            analytic = elliptic_righting_arm(a, b, point["heel_deg"], kg=0.4)
            assert point["righting_arm"] == pytest.approx(analytic, rel=3e-4, abs=2e-6)

    @pytest.mark.parametrize("a,b", SECTIONS)
    def test_mesh_upright_displacement_is_half_the_prism(self, a, b):
        mesh = elliptic_prism_mesh(a, b, 4.0, ring=384)
        assert _signed_volume(mesh) == pytest.approx(math.pi * a * b * 4.0, rel=2e-4)
        curve = righting_arm_curve(mesh, kg=0.4, heels_deg=(0.0, 30.0))
        assert curve.displacement_volume == pytest.approx(math.pi * a * b * 4.0 / 2.0, rel=2e-4)

    def test_discretisation_error_falls_with_ring_count(self):
        errors = []
        for ring in (48, 96, 192):
            mesh = elliptic_prism_mesh(1.0, 0.6, 4.0, ring=ring)
            curve = righting_arm_curve(mesh, kg=0.4, heels_deg=(20.0, 40.0))
            errors.append(
                max(
                    abs(point["righting_arm"] - elliptic_righting_arm(1.0, 0.6, point["heel_deg"], kg=0.4))
                    for point in curve.to_dict()["points"]
                )
            )
        assert errors[0] > errors[1] > errors[2]

    def test_small_angle_slope_of_the_mesh_curve_matches_gm0(self):
        mesh = elliptic_prism_mesh(1.0, 0.6, 4.0, ring=384)
        curve = righting_arm_curve(mesh, kg=0.4, heels_deg=(0.0, 1.0), relative_tolerance=1e-12)
        gz_one_degree = curve.to_dict()["points"][1]["righting_arm"]
        assert gz_one_degree / math.radians(1.0) == pytest.approx(
            elliptic_gm0(1.0, 0.6, kg=0.4), rel=3e-3
        )


class TestWallSidedCaveat:
    """The reason this family is worth having: ``KB + BM - KG`` is not the answer."""

    @pytest.mark.parametrize("a,b", ((1.0, 0.6), (1.0, 0.35)))
    def test_wall_sided_gm0_disagrees_with_the_exact_slope(self, a, b):
        """A wide, flat ellipse is the case the wall-sided formula gets badly wrong."""
        exact = elliptic_gm0(a, b, kg=0.4)
        wall = wall_sided_gm0(a, b, kg=0.4)
        assert abs(wall - exact) > 0.1 * abs(exact)

    def test_the_disagreement_grows_with_the_flattening(self):
        """Error of ``KB + BM - KG`` against the exact slope, as the section flattens."""
        errors = []
        for b in (2.0, 1.0, 0.6, 0.35):
            exact = elliptic_gm0(1.0, b, kg=0.4)
            errors.append(abs(wall_sided_gm0(1.0, b, kg=0.4) - exact) / abs(exact))
        assert errors == sorted(errors)
        assert errors[-1] > 0.5

    def test_the_two_formulas_meet_in_the_circle_limit(self):
        for radius in (0.25, 1.0):
            # at a = b the coupled term vanishes, but the wall-sided KB + BM term does not:
            # BM = 4R/(3 pi) is an artefact of the wall-sided derivation for a circular section.
            exact = elliptic_gm0(radius, radius, kg=0.4)
            wall = wall_sided_gm0(radius, radius, kg=0.4)
            assert exact == pytest.approx(radius - 0.4, rel=1e-14)
            assert wall == pytest.approx(8.0 * radius / (3.0 * math.pi) - 0.4, rel=1e-14)
            assert exact - wall == pytest.approx((1.0 - 8.0 / (3.0 * math.pi)) * radius, rel=1e-12)


class TestRefusals:
    def test_axes_must_be_positive(self):
        with pytest.raises(EllipticSectionError):
            elliptic_kn(0.0, 1.0, 10.0)
        with pytest.raises(EllipticSectionError):
            elliptic_kn(1.0, -0.5, 10.0)

    def test_heel_must_be_within_a_half_turn(self):
        with pytest.raises(EllipticSectionError):
            elliptic_kn(1.0, 0.6, 180.0)

    def test_ring_count_is_guarded(self):
        with pytest.raises(EllipticSectionError):
            elliptic_prism_mesh(1.0, 0.6, 4.0, ring=4)

    def test_length_must_be_positive(self):
        with pytest.raises(EllipticSectionError):
            elliptic_prism_mesh(1.0, 0.6, 0.0)