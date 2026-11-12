"""Analytic anchors for ``tensorcad.shear_flow`` (thin-walled shear flow and shear centre).

Every number here is recomputed from a closed form written in this file; nothing reads an
artifact, so the gates stay honest after a merge.

======  ============================================================================
A1      Free edges carry no flow and junctions conserve it -- the flow is solenoidal on
        the section graph (open tree).
A2      The flow's resultant equals the applied shear force and has no transverse
        component (open and closed sections); reversing the force reverses the arrows.
A3      Open I-section: the web carries the whole shear force and its flow follows the
        textbook ``V Q(z) / I`` at every station, including the extreme-fibre zero of the
        stress at the free edges and the peak at the neutral axis.
A4      Channel shear centre sits at ``e = 3 b^2 / (6 b + h)`` behind the web centreline
        (classical thin-walled result, re-derived here) -- this is the anchor that a
        section which is not doubly symmetric must hit, so A5 is not vacuous.
A5      Doubly symmetric sections (I, box) put the shear centre on the centroid.
A6      Closed cell: ``CLOSED_INT q / t ds = 0`` for pure bending, and the box webs each
        carry exactly ``V/2``; the open I-section carries none of it in the flanges.
A7      Section properties agree with independent closed forms (thin-walled area, centroid,
        second moment) so the denominator of ``q = V Q / I`` is not taken on trust.
A8      Refusals: multi-cell, doubled walls, disconnected walls, ``I_yz != 0`` and bad
        inputs raise :class:`ShearFlowError` instead of returning a wrong number.
A9      Flow is scale-covariant: ``q`` scales with ``V`` and the residuals stay at
        round-off, which a hidden absolute tolerance somewhere would break.
======  ============================================================================
"""

from __future__ import annotations

import math

import pytest

from tensorcad.shear_flow import (
    AXIS_Y,
    AXIS_Z,
    ShearFlowError,
    ThinWalledSection,
    Wall,
    box_section,
    channel_section,
    i_section,
    max_shear_stress,
    shear_centre,
    shear_flow,
    t_section,
)

#: Residual ceiling for the exact-integration gates: round-off, not percent.
ROUND_OFF = 1e-12


def _i_section_closed_form(depth: float, width: float, t_web: float, t_flange: float):
    """Textbook thin-walled I-section: ``(area, I, Z_max)`` about the centroidal axis."""
    area = 2.0 * width * t_flange + depth * t_web
    inertia = 2.0 * (width * t_flange) * (0.5 * depth) ** 2 + t_web * depth ** 3 / 12.0
    return area, inertia, 0.5 * depth


def _box_closed_form(width: float, depth: float, t_web: float, t_flange: float):
    area = 2.0 * width * t_flange + 2.0 * depth * t_web
    inertia = 2.0 * (width * t_flange) * (0.5 * depth) ** 2 + 2.0 * t_web * depth ** 3 / 12.0
    return area, inertia


# ---------------------------------------------------------------------------
# A7  Section properties against independent closed forms
# ---------------------------------------------------------------------------

class TestSectionProperties:
    def test_i_section_matches_closed_forms(self) -> None:
        depth, width, t_w, t_f = 3.0, 2.0, 0.02, 0.03
        section = i_section(depth, width, t_w, t_f)
        area, inertia, _ = _i_section_closed_form(depth, width, t_w, t_f)
        assert section.is_closed is False
        assert section.area == pytest.approx(area, rel=ROUND_OFF)
        assert section.centroid_y == pytest.approx(0.0, abs=ROUND_OFF)
        assert section.centroid_z == pytest.approx(0.0, abs=ROUND_OFF)
        assert section.second_moment_horizontal == pytest.approx(inertia, rel=ROUND_OFF)

    def test_box_matches_closed_forms(self) -> None:
        width, depth, t_w, t_f = 3.0, 2.0, 0.02, 0.025
        section = box_section(width, depth, t_w, t_f)
        area, inertia = _box_closed_form(width, depth, t_w, t_f)
        assert section.is_closed is True
        assert section.area == pytest.approx(area, rel=ROUND_OFF)
        assert section.second_moment_horizontal == pytest.approx(inertia, rel=ROUND_OFF)
        assert section.perimeter == pytest.approx(2.0 * (width + depth), rel=ROUND_OFF)

    def test_product_of_inertia_vanishes_on_symmetric_sections(self) -> None:
        for section in (i_section(3.0, 2.0, 0.02), box_section(3.0, 2.0, 0.02),
                        channel_section(3.0, 1.0, 0.02), t_section(2.0, 2.0, 0.02)):
            assert abs(section.product_of_inertia) <= 1e-12 * section.area * section.perimeter

    def test_channel_centroid_is_off_the_web(self) -> None:
        # A channel's centroid sits between web and flange tips: A = t(h + 2b),
        # S_y = t b^2 (one flange, since the other sits at y = 0), so c_y = b^2 / (h + 2b).
        depth, width, t = 3.0, 1.0, 0.02
        section = channel_section(depth, width, t)
        assert section.centroid_y == pytest.approx(width ** 2 / (depth + 2.0 * width), rel=1e-12)
        assert section.centroid_z == pytest.approx(0.0, abs=1e-12)


# ---------------------------------------------------------------------------
# A1  Free edges and junctions
# ---------------------------------------------------------------------------

class TestFlowKinematics:
    @pytest.mark.parametrize(
        "factory",
        [
            lambda: i_section(3.0, 2.0, 0.02),
            lambda: channel_section(3.0, 1.0, 0.02),
            lambda: t_section(2.0, 2.0, 0.02),
        ],
    )
    def test_free_edges_carry_no_flow(self, factory) -> None:
        flow = shear_flow(factory(), 1.0e3)
        assert flow.free_edge_residual is not None
        assert flow.free_edge_residual <= 1e-12 * 1.0e3

    @pytest.mark.parametrize(
        "section",
        [
            i_section(3.0, 2.0, 0.02),
            channel_section(3.0, 1.0, 0.02),
            t_section(2.0, 2.0, 0.02),
            box_section(3.0, 2.0, 0.02),
        ],
        ids=["i", "channel", "tee", "box"],
    )
    def test_junctions_conserve_flow(self, section: ThinWalledSection) -> None:
        flow = shear_flow(section, 5.0e3)
        assert flow.node_continuity_residual <= 1e-10 * 5.0e3

    def test_tip_flow_is_exactly_zero(self) -> None:
        flow = shear_flow(i_section(3.0, 2.0, 0.02), 1.0)
        deck = flow.edge("deck")
        assert len(deck) == 2  # the deck wall is split at the web junction
        for half in deck:
            assert half.q_at(1.0) == pytest.approx(0.0, abs=1e-18)
            # The magnitude grows from nothing at the free edge to the junction value; the
            # arrow sense is a sign convention (see the module docstring), so it is asserted
            # through equilibrium in TestEquilibrium instead.
            assert abs(half.q_at(0.0)) > abs(half.q_at(0.5)) > 0.0


# ---------------------------------------------------------------------------
# A2  Equilibrium: the flow is a shear force, not a decoration
# ---------------------------------------------------------------------------

class TestEquilibrium:
    @pytest.mark.parametrize(
        "section",
        [
            i_section(3.0, 2.0, 0.02),
            channel_section(3.0, 1.0, 0.02),
            t_section(2.0, 2.0, 0.02),
            box_section(3.0, 2.0, 0.02),
            box_section(3.0, 2.0, 0.02, 0.05),
        ],
        ids=["i", "channel", "tee", "box", "box-thick-flange"],
    )
    @pytest.mark.parametrize("force", [1.0, -2.5e6])
    def test_resultant_equals_applied_shear(self, section: ThinWalledSection, force: float) -> None:
        flow = shear_flow(section, force)
        assert flow.vertical_resultant == pytest.approx(force, rel=ROUND_OFF)
        assert flow.equilibrium_residual <= ROUND_OFF
        assert flow.transverse_residual <= ROUND_OFF

    def test_horizontal_shear_also_balances(self) -> None:
        section = box_section(3.0, 2.0, 0.02)
        flow = shear_flow(section, 1.0e3, axis=AXIS_Y)
        assert flow.horizontal_resultant == pytest.approx(1.0e3, rel=ROUND_OFF)
        assert flow.transverse_residual <= ROUND_OFF

    def test_reversing_the_force_reverses_every_arrow(self) -> None:
        section = channel_section(3.0, 1.0, 0.02)
        forward = shear_flow(section, 7.0e3)
        backward = shear_flow(section, -7.0e3)
        for first, second in zip(forward.edges, backward.edges):
            for a, b in zip(first.q_values, second.q_values):
                assert a == pytest.approx(-b, rel=1e-12, abs=1e-18)

    def test_zero_force_gives_zero_flow(self) -> None:
        flow = shear_flow(box_section(3.0, 2.0, 0.02), 0.0)
        assert all(value == 0.0 for edge in flow.edges for value in edge.q_values)
        assert flow.circulation == 0.0


# ---------------------------------------------------------------------------
# A3  Open I-section: textbook V Q / I
# ---------------------------------------------------------------------------

class TestOpenISection:
    def test_web_carries_the_whole_shear_force(self) -> None:
        depth, width, t_w, t_f = 3.0, 2.0, 0.02, 0.03
        section = i_section(depth, width, t_w, t_f)
        force = 12.0e3
        flow = shear_flow(section, force)
        web = flow.edge("web")
        assert len(web) == 1
        # Integrate the quadratic flow over the web height with Simpson's rule (exact here).
        edge = web[0]
        half = 0.5 * depth
        integral = (2.0 * half / 6.0) * (
            edge.q_at(0.0) + 4.0 * edge.q_at(0.5) + edge.q_at(1.0)
        )
        assert abs(integral) == pytest.approx(abs(force), rel=1e-12)

    def test_web_flow_matches_closed_form_Q_over_I(self) -> None:
        depth, width, t_w, t_f = 3.0, 2.0, 0.02, 0.03
        section = i_section(depth, width, t_w, t_f)
        _, inertia, _ = _i_section_closed_form(depth, width, t_w, t_f)
        force = 1.0
        flow = shear_flow(section, force)
        edge = flow.edge("web")[0]
        for fraction in (0.0, 0.25, 0.5, 0.75, 1.0):
            z = edge.start[1] + fraction * (edge.end[1] - edge.start[1])
            # Q(z) of the part above the cut: the top flange plus the web strip from z up to
            # the flange. Taking |q| because the sweep direction, not the textbook, sets the
            # sign -- the magnitude is the physical content here.
            q_flange = width * t_f * 0.5 * depth
            q_web = t_w * (0.25 * depth * depth - z * z) / 2.0
            expected = force * abs(q_flange + q_web) / inertia
            assert abs(edge.q_at(fraction)) == pytest.approx(expected, rel=1e-12)

    def test_stress_peaks_at_the_neutral_axis_and_vanishes_at_the_tips(self) -> None:
        section = i_section(3.0, 2.0, 0.02, 0.03)
        peak = max_shear_stress(section, 1.0e3)
        assert peak.label == "web"
        assert peak.fraction == pytest.approx(0.5, abs=1e-12)
        assert peak.point[1] == pytest.approx(0.0, abs=1e-12)
        # tau = q / t with q at the neutral axis from the closed form: flange first moment
        # plus half the web.
        _, inertia, _ = _i_section_closed_form(3.0, 2.0, 0.02, 0.03)
        q_peak = (2.0 * 0.03 * 1.5 + 0.02 * 1.5 ** 2 / 2.0) / inertia * 1.0e3
        assert abs(peak.value) == pytest.approx(q_peak / 0.02, rel=1e-12)

    def test_flange_flow_at_the_junction_matches_the_closed_form(self) -> None:
        width, t_f, depth = 2.0, 0.02, 3.0
        section = i_section(depth, width, 0.02, t_f)
        _, inertia, _ = _i_section_closed_form(depth, width, 0.02, t_f)
        flow = shear_flow(section, 1.0)
        # Each half flange carries q = V (b/2) t_f (h/2) / I where it meets the web, and its
        # flow falls linearly to zero at the free edge.
        expected = (0.5 * width) * t_f * (0.5 * depth) / inertia
        for half in flow.edge("deck"):
            assert abs(half.q_at(0.0)) == pytest.approx(expected, rel=1e-12)
            assert abs(half.q_at(1.0)) == pytest.approx(0.0, abs=1e-18)
            # Monotone towards the free edge (it is linear, so one interior check will do).
            assert abs(half.q_at(0.5)) == pytest.approx(expected / 2.0, rel=1e-12)
        assert flow.transverse_residual <= ROUND_OFF


# ---------------------------------------------------------------------------
# A4/A5  Shear centre: the anchor that must not be satisfied by accident
# ---------------------------------------------------------------------------

class TestShearCentre:
    def test_channel_shear_centre_matches_classical_formula(self) -> None:
        # For a uniform-thickness channel the thin-walled shear centre sits behind the web
        # by e = 3 b^2 / (6 b + h): I = 2 (b t)(h/2)^2 + t h^3 / 12, each flange force is
        # F = INT q ds = t h b^2 / (4 I) * V, and their moment F h / V gives e.
        for depth, width, thickness in ((3.0, 1.0, 0.02), (2.0, 1.5, 0.01), (4.0, 0.8, 0.03)):
            section = channel_section(depth, width, thickness)
            centre = shear_centre(section)
            expected = 3.0 * width ** 2 / (6.0 * width + depth)
            # The web centreline is at y = 0 and the flanges project towards +y, so the
            # shear centre lies on the far side of the web.
            assert centre.y == pytest.approx(-expected, abs=1e-12 * max(1.0, width))
            assert centre.z == pytest.approx(section.centroid_z, abs=1e-12)

    def test_channel_shear_centre_departs_from_the_centroid(self) -> None:
        section = channel_section(3.0, 1.0, 0.02)
        centre = shear_centre(section)
        # The centroid is inside the section; the shear centre is not. That difference is
        # the physical content of A4, so it must not collapse to zero.
        assert abs(centre.y - section.centroid_y) > 0.2

    @pytest.mark.parametrize(
        "section", [i_section(3.0, 2.0, 0.02), box_section(3.0, 2.0, 0.02, 0.05)],
        ids=["i", "box"],
    )
    def test_doubly_symmetric_sections_centre_on_the_centroid(
        self, section: ThinWalledSection
    ) -> None:
        centre = shear_centre(section)
        assert centre.y == pytest.approx(section.centroid_y, abs=1e-15)
        assert centre.z == pytest.approx(section.centroid_z, abs=1e-15)
        assert abs(centre.torque_vertical) <= 1e-15
        assert abs(centre.torque_horizontal) <= 1e-15

    def test_tee_shear_centre_stays_on_the_symmetry_axis(self) -> None:
        section = t_section(2.0, 2.0, 0.02)
        centre = shear_centre(section)
        assert centre.y == pytest.approx(0.0, abs=1e-15)
        # Vertical shear on a tee produces a non-zero horizontal flow resultant? No: the
        # flanges' flows cancel, so the vertical shear centre is a well defined point.
        assert centre.vertical_flow.transverse_residual <= 1e-12


# ---------------------------------------------------------------------------
# A6  Closed cell: circulation, no twist, and where the force goes
# ---------------------------------------------------------------------------

class TestClosedCell:
    def test_no_twist_condition(self) -> None:
        for section in (box_section(3.0, 2.0, 0.02), box_section(2.0, 4.0, 0.03, 0.01)):
            flow = shear_flow(section, 8.0e3)
            assert flow.no_twist_residual is not None
            assert abs(flow.no_twist_residual) <= ROUND_OFF
            assert flow.circulation != 0.0

    def test_open_section_has_no_circulation(self) -> None:
        flow = shear_flow(i_section(3.0, 2.0, 0.02), 1.0)
        assert flow.circulation == 0.0
        assert flow.no_twist_residual is None
        assert flow.free_edge_residual is not None

    def test_box_webs_each_carry_half_the_shear(self) -> None:
        width, depth, t_w, t_f = 3.0, 2.0, 0.02, 0.025
        force = 9.0e3
        flow = shear_flow(box_section(width, depth, t_w, t_f), force)
        for label in ("starboard", "port"):
            edge = flow.edge(label)[0]
            integral = (depth / 6.0) * (
                edge.q_at(0.0) + 4.0 * edge.q_at(0.5) + edge.q_at(1.0)
            )
            assert abs(integral) == pytest.approx(abs(force) / 2.0, rel=1e-12)

    def test_box_flow_peaks_in_the_web_and_crosses_zero_in_the_flange(self) -> None:
        section = box_section(3.0, 2.0, 0.02)
        flow = shear_flow(section, 1.0)
        top = flow.edge("top")[0]
        # Symmetric box: the top flange carries a flow that reverses at mid-width (the
        # symmetry plane), so the flange's flow must change sign.
        assert top.q_at(0.0) * top.q_at(1.0) < 0.0
        peak = max_shear_stress(section, 1.0)
        assert abs(peak.value) >= max(abs(value) for edge in flow.edges for value in edge.stress_values) - 1e-15

    def test_circulation_scales_with_the_force(self) -> None:
        section = box_section(3.0, 2.0, 0.02)
        first = shear_flow(section, 1.0).circulation
        second = shear_flow(section, 250.0).circulation
        assert second == pytest.approx(250.0 * first, rel=1e-12)


# ---------------------------------------------------------------------------
# A9  Scale covariance
# ---------------------------------------------------------------------------

class TestScaleCovariance:
    def test_flow_scales_linearly_with_the_force(self) -> None:
        section = box_section(3.0, 2.0, 0.02)
        base = shear_flow(section, 1.0)
        scaled = shear_flow(section, 1.0e6)
        for first, second in zip(base.edges, scaled.edges):
            for a, b in zip(first.q_values, second.q_values):
                assert b == pytest.approx(1.0e6 * a, rel=1e-12, abs=1e-3)

    def test_residuals_are_scale_free(self) -> None:
        section = channel_section(3.0, 1.0, 0.02)
        for force in (1e-3, 1.0, 1e6):
            flow = shear_flow(section, force)
            assert flow.equilibrium_residual <= ROUND_OFF
            assert flow.node_continuity_residual <= 1e-10 * abs(force)


# ---------------------------------------------------------------------------
# A8  Refusals: a wrong answer is worse than an exception
# ---------------------------------------------------------------------------

class TestRefusals:
    def test_multi_cell_is_rejected(self) -> None:
        # Two adjacent cells share the middle web: 6 vertices, 7 edges, two faces. The
        # flanges are split at the centreline because walls connect at vertices only.
        walls = [
            Wall(((-1.0, -1.0), (0.0, -1.0)), 0.02, "bottom-port"),
            Wall(((0.0, -1.0), (1.0, -1.0)), 0.02, "bottom-starboard"),
            Wall(((1.0, -1.0), (1.0, 1.0)), 0.02, "starboard"),
            Wall(((1.0, 1.0), (0.0, 1.0)), 0.02, "top-starboard"),
            Wall(((0.0, 1.0), (-1.0, 1.0)), 0.02, "top-port"),
            Wall(((-1.0, 1.0), (-1.0, -1.0)), 0.02, "port"),
            Wall(((0.0, -1.0), (0.0, 1.0)), 0.02, "centreline"),
        ]
        with pytest.raises(ShearFlowError, match="Multi-cell"):
            ThinWalledSection(walls)

    def test_doubled_wall_is_rejected(self) -> None:
        walls = [
            Wall(((0.0, 0.0), (1.0, 0.0)), 0.02, "a"),
            Wall(((1.0, 0.0), (0.0, 0.0)), 0.02, "b"),
            Wall(((0.0, 0.0), (0.0, -1.0)), 0.02, "c"),
        ]
        with pytest.raises(ShearFlowError, match="coincide"):
            ThinWalledSection(walls)

    def test_disconnected_walls_are_rejected(self) -> None:
        walls = [
            Wall(((0.0, 0.0), (1.0, 0.0)), 0.02, "a"),
            Wall(((5.0, 0.0), (6.0, 0.0)), 0.02, "b"),
        ]
        with pytest.raises(ShearFlowError, match="one connected"):
            ThinWalledSection(walls)

    def test_product_of_inertia_is_rejected(self) -> None:
        # An L-shaped (angle) section is genuinely uncoupled at the principal axes only.
        walls = [
            Wall(((0.0, 0.0), (2.0, 0.0)), 0.02, "foot"),
            Wall(((0.0, 0.0), (0.0, 2.0)), 0.02, "leg"),
        ]
        section = ThinWalledSection(walls)
        assert abs(section.product_of_inertia) > 0.0
        with pytest.raises(ShearFlowError, match="product of inertia"):
            shear_flow(section, 1.0)

    def test_bad_inputs_are_rejected(self) -> None:
        with pytest.raises(ShearFlowError):
            Wall(((0.0, 0.0),), 0.02)
        with pytest.raises(ShearFlowError):
            Wall(((0.0, 0.0), (1.0, 0.0)), 0.0)
        with pytest.raises(ShearFlowError):
            i_section(-1.0, 2.0, 0.02)
        section = i_section(3.0, 2.0, 0.02)
        with pytest.raises(ShearFlowError):
            shear_flow(section, float("nan"))
        with pytest.raises(ShearFlowError):
            shear_flow(section, 1.0, axis="x")
        with pytest.raises(ShearFlowError):
            shear_flow(section, 1.0, stations=1)
        with pytest.raises(ShearFlowError):
            shear_centre(section, unit_force=0.0)
        with pytest.raises(ShearFlowError):
            shear_flow(section, 1.0).edge("deck")[0].q_at(1.5)

    def test_closed_flag_must_agree_with_the_graph(self) -> None:
        with pytest.raises(ShearFlowError, match="contradicts"):
            ThinWalledSection([Wall(((0.0, 0.0), (1.0, 0.0)), 0.02)], closed=True)


# ---------------------------------------------------------------------------
# Reporting contract
# ---------------------------------------------------------------------------

class TestReporting:
    def test_flow_report_is_json_ready(self) -> None:
        import json

        flow = shear_flow(box_section(3.0, 2.0, 0.02), 1.0e3)
        payload = json.dumps(flow.to_dict())
        restored = json.loads(payload)
        assert len(restored["edges"]) == 4
        assert restored["equilibrium_residual"] <= ROUND_OFF
        assert len(restored["edges"][0]["q"]) == 5

    def test_section_report_round_trips_the_walls(self) -> None:
        import json

        section = channel_section(3.0, 1.0, 0.02)
        payload = json.loads(json.dumps(section.to_dict()))
        assert payload["is_closed"] is False
        assert len(payload["walls"]) == 3
        assert payload["walls"][0]["label"] == "web"

    def test_stations_are_closed_on_both_ends(self) -> None:
        flow = shear_flow(i_section(3.0, 2.0, 0.02), 1.0, stations=9)
        for edge in flow.edges:
            assert edge.s_values[0] == 0.0
            assert edge.s_values[-1] == pytest.approx(edge.length, rel=1e-15)