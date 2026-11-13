"""Transom stern and bulbous bow: the two ends the Series 60 family did not have.

Every claim here is a number: the transom must hold the integrals it is asked to hold
while actually changing the stern, and the bulb's closed forms are checked against a
second parameterisation, against exact special cases, and against a mesh of the solid
it exports.
"""
from __future__ import annotations

import math

import pytest

from tensorcad import series60
from tensorcad.series60 import (
    Series60Parameters,
    analytic_series60_hydrostatics,
    build_bulb_mesh,
    bulb_centroid_x_ahead,
    bulb_volume,
    generate_series60,
    gauss_legendre_unit,
    series60_appended_hydrostatics,
    series60_bulb,
    transom_area,
)

L, B, T = 10.0, 1.5, 0.6


def parent(cb: float, **overrides) -> Series60Parameters:
    """A published-parent-fidelity hull, at the resolution the audits use."""
    published = {0.60: 0.700, 0.65: 0.727, 0.70: 0.753, 0.75: 0.780, 0.80: 0.807}
    lcb = {0.60: 0.475, 0.65: 0.483, 0.70: 0.490, 0.75: 0.496, 0.80: 0.502}
    base = dict(length=L, beam=B, draft=T, cb=cb, cwp=published[cb], lcb=lcb[cb])
    base.update(overrides)
    return Series60Parameters(**base)


class TestTransomParameterisation:
    def test_zero_transom_is_the_pointed_stern_bit_for_bit(self):
        plain = parent(0.70).shape
        explicit = parent(0.70, transom_ratio=0.0).shape
        assert plain == explicit
        assert plain.transom_ratio == 0.0

    def test_stern_deck_taper_equals_the_transom_ratio(self):
        law = parent(0.70, transom_ratio=0.2).shape
        assert law.deck_taper(0.0) == pytest.approx(0.2, rel=0, abs=1e-15)
        assert law.deck_taper(0.5) == pytest.approx(1.0, rel=0, abs=1e-12)
        assert law.deck_taper(1.0) == 0.0

    @pytest.mark.parametrize("tau", [0.05, 0.10, 0.20, 0.30])
    def test_requested_integrals_still_close(self, tau):
        p = parent(0.70, transom_ratio=tau)
        law = p.shape
        assert law.cwp == pytest.approx(0.753, rel=0, abs=1e-12)
        assert law.lcb_ratio == pytest.approx(0.490, rel=0, abs=1e-12)
        ref = analytic_series60_hydrostatics(p)
        assert ref["block_coefficient"] == pytest.approx(0.70, rel=0, abs=1e-12)
        assert ref["displacement_volume"] == pytest.approx(L * B * T * 0.70, rel=1e-15)

    def test_transom_changes_the_shape_not_just_the_labels(self):
        plain = parent(0.70).shape
        transom = parent(0.70, transom_ratio=0.20).shape
        assert abs(transom.n_aft - plain.n_aft) > 0.1
        assert abs(transom.gain - plain.gain) > 0.01
        assert transom.vcb_ratio != pytest.approx(plain.vcb_ratio, rel=1e-6)

    def test_legacy_single_exponent_holds_cb_with_a_transom(self):
        """Without ``cwp`` the family has one exponent, so it re-solves to keep Cb."""
        p = Series60Parameters(length=L, beam=B, draft=T, cb=0.70, transom_ratio=0.2)
        law = p.shape
        assert law.cwp == pytest.approx(0.70 / law.cm_mid, rel=1e-12)
        assert analytic_series60_hydrostatics(p)["displacement_volume"] == pytest.approx(
            L * B * T * 0.70, rel=1e-12
        )

    def test_transom_is_a_wider_stern(self):
        plain = parent(0.70).shape
        transom = parent(0.70, transom_ratio=0.20).shape
        assert transom.deck_taper(0.05) > plain.deck_taper(0.05)

    def test_fore_half_carries_no_transom(self):
        """The blend is an aft-end story; the fore half keeps the bare taper.

        Added after the mutation audit: a mutant that blended ``xi < 0.6`` survived,
        because the blend value at midship is ``1`` for every tau and the closed forms
        never sample the fore half's deck taper.
        """
        law = parent(0.70, transom_ratio=0.20).shape
        for xi in (0.5, 0.55, 0.7, 0.9):
            base = 1.0 - abs(2.0 * xi - 1.0) ** law.side_exponent(xi)
            assert law.deck_taper(xi) == pytest.approx(base, rel=0, abs=1e-15)

    def test_closed_form_moments_match_the_law_they_describe(self):
        """The solver's own moments, integrated against the law it resolved.

        ``Cwp`` and ``LCB`` are hit by construction -- the solver buys them -- so a
        mis-weighted moment is invisible from the targets alone: it is absorbed by a
        compensating exponent. The only way to see it is to integrate ``deck_taper``
        and ``fullness`` themselves, which is what this does (and what catches the
        transom square-moment mutant that survived the first audit round).
        """
        law = parent(0.70, transom_ratio=0.25).shape
        nodes, weights = gauss_legendre_unit(400)
        cwp = sum(w * law.deck_taper(x) for x, w in zip(nodes, weights))
        assert cwp == pytest.approx(law.cwp, abs=1e-8)  # 400-node rule floor
        volume = sum(w * law.deck_taper(x) * law.fullness(x) for x, w in zip(nodes, weights))
        moment = sum(w * x * law.deck_taper(x) * law.fullness(x)
                     for x, w in zip(nodes, weights))
        assert moment / volume == pytest.approx(law.lcb_ratio, abs=1e-6)

    def test_parameter_range_is_enforced(self):
        with pytest.raises(ValueError):
            Series60Parameters(transom_ratio=0.6).validate()
        with pytest.raises(ValueError):
            Series60Parameters(transom_ratio=-0.1).validate()


class TestTransomGeometry:
    def test_offsets_show_the_transom_breadth(self):
        p = parent(0.70, transom_ratio=0.20, nx=41, nz=13)
        model = generate_series60(p)
        x0, half_breadth = model.offsets[0]
        assert x0 == 0.0
        assert half_breadth == pytest.approx(0.20 * B / 2, rel=1e-12)
        assert model.offsets[-1][1] == 0.0  # the bow is still a point

    def test_mesh_stays_a_closed_manifold(self):
        p = parent(0.70, transom_ratio=0.20, nx=81, nz=21)
        model = generate_series60(p)
        assert model.validation.errors == []
        names = {check.name: check.passed for check in model.validation.checks}
        for required in ("edge_manifold", "vertex_manifold", "euler_characteristic",
                         "outward_orientation", "no_zero_area_faces"):
            assert names.get(required, True), required

    def test_immersed_stern_area_measured_from_the_mesh(self):
        """The flat face the loft added must be the area the closed form claims."""
        p = parent(0.70, transom_ratio=0.20, nx=161, nz=41)
        model = generate_series60(p)
        mesh = model.mesh
        area = 0.0
        for face in mesh.faces:
            points = [mesh.vertices[i] for i in face]
            if all(abs(pt[0]) < 1e-9 for pt in points):
                a, b, c = points
                u = [b[i] - a[i] for i in range(3)]
                v = [c[i] - a[i] for i in range(3)]
                cross = [u[1] * v[2] - u[2] * v[1],
                         u[2] * v[0] - u[0] * v[2],
                         u[0] * v[1] - u[1] * v[0]]
                area += 0.5 * math.dist((0.0, 0.0, 0.0), cross)
        reference = transom_area(p)
        assert reference > 0.0
        assert area == pytest.approx(reference, rel=0.01)

    def test_mesh_displacement_matches_the_closed_form(self):
        p = parent(0.70, transom_ratio=0.20, nx=241, nz=81)
        model = generate_series60(p)
        ref = analytic_series60_hydrostatics(p)
        volume_error = (model.hydrostatics.displacement_volume
                        / ref["displacement_volume"] - 1.0)
        waterplane_error = (model.hydrostatics.waterplane_area
                            / ref["waterplane_area"] - 1.0)
        assert abs(volume_error) < 0.01
        assert abs(waterplane_error) < 0.005


class TestBulbClosedForms:
    def test_absent_by_default(self):
        p = parent(0.70)
        assert series60_bulb(p) is None
        assert build_bulb_mesh(p) is None
        appended = series60_appended_hydrostatics(p, generate_series60(
            parent(0.70, nx=41, nz=13)).hydrostatics)
        assert appended["bulb_displacement_volume_m3"] == 0.0

    def test_partial_specification_is_rejected(self):
        with pytest.raises(ValueError):
            Series60Parameters(bulb_length_ratio=0.02).validate()
        with pytest.raises(ValueError):
            Series60Parameters(bulb_area_ratio=0.08).validate()

    def test_fully_submerged_volume(self):
        a, b = 0.2, 0.15
        assert bulb_volume(a, b, b) == pytest.approx((2.0 / 3.0) * math.pi * a * b * b,
                                                     rel=1e-15)
        assert bulb_volume(a, b, 2.0 * b) == pytest.approx((2.0 / 3.0) * math.pi * a * b * b,
                                                           rel=1e-15)
        # Between one and two radii above the axis the bulb is *still* fully submerged;
        # the closed form would keep integrating past its own pole (it reaches zero at
        # h = 2b), which is how the offset clamp mutant survived the first round.
        assert bulb_volume(a, b, 1.5 * b) == pytest.approx(
            (2.0 / 3.0) * math.pi * a * b * b, rel=1e-15)

    def test_half_and_none_cases(self):
        a, b = 0.2, 0.15
        assert bulb_volume(a, b, 0.0) == pytest.approx(math.pi * a * b * b / 3.0, rel=1e-15)
        assert bulb_volume(a, b, -b) == 0.0
        assert bulb_volume(a, b, -3.0 * b) == 0.0

    @pytest.mark.parametrize("ratio", [-0.9, -0.5, 0.0, 0.5, 0.9, 1.0])
    def test_volume_against_an_independent_parameterisation(self, ratio):
        """Same solid, integrated over horizontal slabs instead of the axial variable.

        The slab area is ``(pi a b / 2)(1 - (z/b)^2)`` -- a quadratic in ``z``, so a
        three-point Gauss rule is exact and this is a genuinely different route to the
        same number.
        """
        a, b = 0.2, 0.15
        h = ratio * b
        nodes, weights = gauss_legendre_unit(3)
        slab = 0.0
        for node, weight in zip(nodes, weights):
            z = -b + (h + b) * node
            slab += weight * (h + b) * 0.5 * math.pi * a * b * (1.0 - (z / b) ** 2)
        assert bulb_volume(a, b, h) == pytest.approx(slab, rel=1e-12, abs=1e-15)

    def test_centroid_collapses_to_the_hemisphere_result(self):
        a, b = 0.2, 0.15
        assert bulb_centroid_x_ahead(a, b, b) == pytest.approx(3.0 * a / 8.0, rel=1e-14)
        assert bulb_centroid_x_ahead(a, b, 0.0) == pytest.approx(3.0 * a / 8.0, rel=1e-14)

    def test_fp_area_is_the_holtrop_ratio(self):
        p = parent(0.70, bulb_length_ratio=0.02, bulb_area_ratio=0.08)
        bulb = series60_bulb(p)
        mid_section = B * T * p.shape.cm_mid
        assert bulb.fp_area == pytest.approx(0.08 * mid_section, rel=1e-15)
        assert bulb.fp_area == pytest.approx(math.pi * bulb.radius ** 2, rel=1e-15)

    def test_piercing_the_waterline_keeps_only_the_submerged_part(self):
        p = parent(0.70, bulb_length_ratio=0.02, bulb_area_ratio=0.08,
                   bulb_centre_height_ratio=1.0)
        bulb = series60_bulb(p)
        assert bulb.pierces_waterline
        assert bulb.volume == pytest.approx(math.pi * bulb.semi_axis * bulb.radius ** 2 / 3.0,
                                            rel=1e-14)
        assert bulb.waterplane_area == pytest.approx(
            math.pi * bulb.semi_axis * bulb.radius / 2.0, rel=1e-14)
        assert bulb.fp_area_submerged == pytest.approx(math.pi * bulb.radius ** 2 / 2.0,
                                                       rel=1e-14)
        model = generate_series60(p)
        crossing = [c for c in model.validation.checks if c.name == "bulb_waterline_crossing"]
        assert crossing and crossing[0].detail  # the warning says what was cut


class TestBulbMesh:
    @staticmethod
    def _measures(mesh) -> tuple[float, float]:
        volume = 0.0
        area = 0.0
        flat = 0.0
        for face in mesh.faces:
            a, b, c = (mesh.vertices[i] for i in face)
            volume += (a[0] * (b[1] * c[2] - b[2] * c[1])
                       - a[1] * (b[0] * c[2] - b[2] * c[0])
                       + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6.0
            u = [b[i] - a[i] for i in range(3)]
            v = [c[i] - a[i] for i in range(3)]
            cross = [u[1] * v[2] - u[2] * v[1],
                     u[2] * v[0] - u[0] * v[2],
                     u[0] * v[1] - u[1] * v[0]]
            face_area = 0.5 * math.dist((0.0, 0.0, 0.0), cross)
            area += face_area
            if all(abs(pt[0] - 10.0) < 1e-9 for pt in (a, b, c)):
                flat += face_area
        return volume, area - flat

    def test_exported_solid_volume_and_lateral_area(self):
        p = parent(0.70, bulb_length_ratio=0.02, bulb_area_ratio=0.08)
        bulb = series60_bulb(p)
        volume, lateral = self._measures(build_bulb_mesh(p, n_theta=33, n_u=17))
        assert volume == pytest.approx(bulb.volume_if_submerged, rel=0.02)
        assert lateral == pytest.approx(bulb.wetted_area, rel=0.02)

    def test_wetted_area_against_a_brute_force_summation(self):
        """A piercing bulb: sum the area element over a dense (u, theta) grid.

        The quadrature's arc measure ``2(pi - acos(h/r))`` is what decides how much of
        each ring is wet, and a fully submerged bulb cannot test it (every ring is then
        a full circle) -- which is how the arc mutant survived the first round.
        """
        p = parent(0.70, bulb_length_ratio=0.02, bulb_area_ratio=0.08,
                   bulb_centre_height_ratio=1.0)
        bulb = series60_bulb(p)
        a, b, h = bulb.semi_axis, bulb.radius, bulb.waterline_offset
        du, dtheta = 1.0 / 4000.0, 2.0 * math.pi / 720.0
        brute = 0.0
        for i in range(4000):
            u = (i + 0.5) * du
            r = b * math.sqrt(max(0.0, 1.0 - u * u))
            if r <= 0.0:
                continue
            ring = math.sqrt((b * b * u) ** 2 + (a * r) ** 2)
            for j in range(720):
                theta = (j + 0.5) * dtheta
                if bulb.centre_z + b * math.sqrt(max(0.0, 1.0 - u * u)) * math.cos(theta) <= T:
                    brute += ring * du * dtheta
        assert bulb.volume == pytest.approx(math.pi * a * b * b / 3.0, rel=1e-14)
        assert bulb.wetted_area == pytest.approx(brute, rel=0.01)

    def test_wetted_area_at_the_axis_is_half_of_the_submerged_area(self):
        p = parent(0.70, bulb_length_ratio=0.02, bulb_area_ratio=0.08,
                   bulb_centre_height_ratio=1.0)
        bulb = series60_bulb(p)
        submerged = series60._bulb_wetted_area(bulb.semi_axis, bulb.radius, bulb.radius)
        assert bulb.wetted_area == pytest.approx(0.5 * submerged, rel=1e-12)

    def test_faceting_error_shrinks_with_resolution(self):
        p = parent(0.70, bulb_length_ratio=0.02, bulb_area_ratio=0.08)
        reference = series60_bulb(p).volume_if_submerged
        coarse, _ = self._measures(build_bulb_mesh(p, n_theta=17, n_u=9))
        fine, _ = self._measures(build_bulb_mesh(p, n_theta=65, n_u=33))
        assert abs(fine / reference - 1.0) < abs(coarse / reference - 1.0)
        assert abs(fine / reference - 1.0) < 0.005

    def test_bulb_mesh_is_homogeneous_and_closed(self):
        p = parent(0.70, bulb_length_ratio=0.02, bulb_area_ratio=0.08)
        model = generate_series60(p)
        assert model.bulb_mesh is not None
        assert model.validation.errors == []
        names = {check.name for check in model.validation.checks}
        assert "bulb_euler_characteristic_genus_zero" in names


class TestAppendedHydrostatics:
    def test_totals_add_up_and_do_not_double_count(self):
        p = parent(0.70, transom_ratio=0.20, bulb_length_ratio=0.02,
                   bulb_area_ratio=0.08, nx=161, nz=41)
        model = generate_series60(p)
        appended = series60_appended_hydrostatics(p, model.hydrostatics)
        bulb = series60_bulb(p)
        assert appended["total_displacement_volume_m3"] == pytest.approx(
            appended["hull_displacement_volume_m3"] + bulb.volume, rel=1e-15)
        hull_lcb = model.hydrostatics.longitudinal_centre_buoyancy
        bulb_x = L + bulb.centroid_x_ahead
        assert hull_lcb < appended["total_lcb_m"] < bulb_x
        assert appended["total_wetted_surface_area_m2"] > model.hydrostatics.wetted_surface_area
        assert appended["total_displacement_mass_kg"] == pytest.approx(
            appended["total_displacement_volume_m3"] * p.water_density, rel=1e-15)

    def test_metadata_declares_two_bodies(self):
        p = parent(0.70, transom_ratio=0.20, bulb_length_ratio=0.02, bulb_area_ratio=0.08)
        metadata = generate_series60(p, nx=41, nz=13).metadata()
        assert metadata["appendages"]["bodies"] == 2
        assert metadata["appendages"]["bulb"]["fp_area_m2"] > 0.0
        assert metadata["transom"]["immersed_area_m2"] > 0.0