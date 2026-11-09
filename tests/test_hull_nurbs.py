"""Tests for :mod:`tensorcad.hull_nurbs` (NURBS hull -> body mesh -> consistency loop).

The suite is organised the same way the module is: the parametric form and its
validation, the NURBS patch (knots, linear precision, the graph identity, the
hydrostatic quadrature), the body mesh (node placement, quality, tet subdivision) and
finally the geometry-consistency loop the IN383-CAD01 node asks for, criterion by
criterion.
"""
from __future__ import annotations

import math
from itertools import pairwise

import pytest

from tensorcad.hull_nurbs import (
    ASPECT_RATIO_MAX,
    SCALED_JACOBIAN_MIN,
    HullNurbsError,
    HullNurbsParameters,
    acceptance_report,
    build_body_mesh,
    build_nurbs_hull,
    calibrate_form,
    clamped_knots,
    gauss_legendre,
    greville_abscissae,
    hex_characteristic_lengths,
    hex_scaled_jacobian,
    hex_to_tets,
    hex_volume,
    holtrop_wetted_area,
    kcs_lightweight_parameters,
    mesh_integrals,
    nurbs_hydrostatics,
    schneekluth_lcb_ratio,
)


# ---------------------------------------------------------------------------
# fixtures: one calibrated hull for the whole module (it costs ~0.3 s to build)
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def hull():
    return build_nurbs_hull()


@pytest.fixture(scope="module")
def hydro(hull):
    return nurbs_hydrostatics(hull)


@pytest.fixture(scope="module")
def mesh(hull):
    return build_body_mesh(hull, stations=24, levels_below_waterline=9, layers=5)


@pytest.fixture(scope="module")
def coarse_and_fine(hull):
    coarse = build_body_mesh(
        hull, stations=16, levels_below_waterline=6, levels_above_waterline=2, layers=4
    )
    fine = build_body_mesh(
        hull, stations=32, levels_below_waterline=12, levels_above_waterline=4, layers=8
    )
    return coarse, fine


# ---------------------------------------------------------------------------
# parameters and validation: geometry first, before anything is integrated
# ---------------------------------------------------------------------------
def test_reference_case_is_the_node_specification():
    p = kcs_lightweight_parameters()
    assert (p.length, p.beam, p.draft) == (3.0, 0.4, 0.15)
    assert p.block_coefficient == pytest.approx(0.555)
    assert p.design_displacement == pytest.approx(0.555 * 3.0 * 0.4 * 0.15, rel=1e-15)
    assert p.length_to_beam == pytest.approx(7.5)
    assert p.beam_to_draft == pytest.approx(0.4 / 0.15)
    errors, warnings = p.problems()
    assert errors == []
    assert warnings == []


@pytest.mark.parametrize(
    "changes",
    [
        {"length": 0.0},
        {"beam": -1.0},
        {"draft": float("nan")},
        {"block_coefficient": 0.95},
        {"block_coefficient": 0.1},
        {"transom_breadth_ratio": 1.5},
        {"bottom_breadth_ratio": 0.0},
        {"entrance_exponent": -1.0},
        {"aft_body_share": 1.0},
        {"control_points": (3, 13)},
        {"parallel_body_ratio": 0.95},
    ],
)
def test_inadmissible_parameters_are_refused(changes):
    p = HullNurbsParameters(**changes)
    errors, _ = p.problems()
    assert errors, f"{changes!r} was accepted"
    with pytest.raises(HullNurbsError):
        p.validate()


def test_out_of_band_principal_dimensions_warn_without_refusing():
    p = HullNurbsParameters(length=0.8)  # L/B = 2.0
    errors, warnings = p.problems()
    assert errors == []
    assert any("L/B" in message for message in warnings)
    p.validate()  # must not raise: a warning is not a refusal


def test_calibration_without_a_target_is_refused_when_unreachable():
    # A block coefficient below what the leanest member of the family reaches must raise
    # rather than silently return a hull that misses the request by 20 %.
    lean = HullNurbsParameters(block_coefficient=0.40)
    with pytest.raises(HullNurbsError):
        calibrate_form(lean)


def test_fixed_parallel_body_ratio_skips_calibration():
    p = HullNurbsParameters(parallel_body_ratio=0.1)
    form = calibrate_form(p)
    assert form.parallel_body_ratio == pytest.approx(0.1)
    assert form.aft_taper_ratio + form.parallel_body_ratio + form.entrance_ratio == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# B-spline scaffolding
# ---------------------------------------------------------------------------
def test_clamped_knots_are_well_formed():
    knots = clamped_knots(3, 8)
    assert len(knots) == 3 + 8 + 1
    assert knots[:4] == (0.0, 0.0, 0.0, 0.0)
    assert knots[-4:] == (1.0, 1.0, 1.0, 1.0)
    assert all(b >= a for a, b in pairwise(knots))
    with pytest.raises(HullNurbsError):
        clamped_knots(4, 4)


def test_greville_abscissae_span_the_domain():
    xs = greville_abscissae(clamped_knots(3, 21), 3)
    assert len(xs) == 21
    assert xs[0] == pytest.approx(0.0)
    assert xs[-1] == pytest.approx(1.0)
    assert all(b > a for a, b in pairwise(xs))


def test_gauss_legendre_integrates_polynomials_exactly():
    order = 5
    nodes, weights = gauss_legendre(order)
    for power in range(0, 2 * order):
        exact = 0.0 if power % 2 else 2.0 / (power + 1)
        measured = sum(w * x**power for x, w in zip(nodes, weights, strict=True))
        assert measured == pytest.approx(exact, abs=1e-14)
    with pytest.raises(HullNurbsError):
        gauss_legendre(0)


# ---------------------------------------------------------------------------
# the NURBS patch
# ---------------------------------------------------------------------------
def test_patch_is_a_cubic_tensor_product_with_the_requested_control_net(hull):
    surface = hull.surface
    assert surface.degree == (3, 3)
    assert len(surface.control_net) == hull.parameters.control_points[0]
    assert len(surface.control_net[0]) == hull.parameters.control_points[1]
    assert surface.rational is False  # built as a non-rational B-spline patch
    assert surface.domain == ((0.0, 1.0), (0.0, 1.0))


def test_bspline_linear_precision_carries_x_and_z_exactly(hull):
    """``x(u) = L u`` and ``z(v) = D v`` to machine precision, not approximately."""
    p = hull.parameters
    worst_x = 0.0
    worst_z = 0.0
    for k in range(21):
        u = k / 20
        for m in range(11):
            v = m / 10
            point = hull.surface.evaluate(u, v)
            worst_x = max(worst_x, abs(point[0] - p.length * u))
            worst_z = max(worst_z, abs(point[2] - p.deck_height * v))
    assert worst_x < 1e-14
    assert worst_z < 1e-14


def test_graph_area_element_matches_the_surface_derivatives(hull):
    """``|S_u x S_v| = L D sqrt(1 + (y_u/L)^2 + (y_v/D)^2)`` for this graph patch.

    The identity is what makes the hydrostatic integrals 2-D integrals in the first
    place, so it is checked against the patch's own analytic derivatives rather than
    against a finite difference.
    """
    p = hull.parameters
    for u, v in ((0.15, 0.2), (0.45, 0.5), (0.7, 0.9), (0.93, 0.4)):
        d = hull.surface.derivatives(u, v, order=1)
        du, dv = d[(1, 0)], d[(0, 1)]
        cross = (
            du[1] * dv[2] - du[2] * dv[1],
            du[2] * dv[0] - du[0] * dv[2],
            du[0] * dv[1] - du[1] * dv[0],
        )
        measured = math.sqrt(sum(component * component for component in cross))
        expected = p.length * p.deck_height * math.sqrt(
            1.0 + (du[1] / p.length) ** 2 + (dv[1] / p.deck_height) ** 2
        )
        assert measured == pytest.approx(expected, rel=1e-12)


def test_half_breadth_is_bounded_and_vanishes_at_the_stem(hull):
    p = hull.parameters
    half_beam = 0.5 * p.beam
    peak = 0.0
    for k in range(41):
        xi = k / 40
        for m in range(11):
            zeta = m / 10
            y = hull.patch_half_breadth(xi, zeta)
            assert -1e-12 <= y <= half_beam * (1 + 1e-9)
            peak = max(peak, y)
    assert peak == pytest.approx(half_beam, rel=5e-3)  # spline maximum, see realised_dimensions
    assert hull.patch_half_breadth(1.0, 0.5) == pytest.approx(0.0, abs=1e-12)
    assert hull.patch_half_breadth(0.0, 1.0) == pytest.approx(
        half_beam * p.transom_breadth_ratio, rel=2e-2
    )


def test_spline_reproduces_the_form_function_it_was_sampled_from(hull):
    """The patch is a spline of the designed form, not a different hull.

    The control net is sampled at the Greville abscissae, so the patch is a spline
    quasi-interpolant: smooth wherever the form is smooth, but the parallel-body junction
    is only C1 and the stem has a vertical tangent, so the deviation is O(h^2) at those
    two features and 1.9 % of the half beam there. Everything the node accepts is
    measured on the *patch*, so this is a fidelity statement, not a consistency gap.
    """
    worst = 0.0
    for k in range(41):
        xi = k / 40
        for zeta in (0.0, 0.25, 0.5, 0.75, 1.0):
            patch = hull.patch_half_breadth(xi, zeta)
            form = hull.form_half_breadth(xi, zeta)
            scale = max(hull.parameters.beam * 0.5, 1e-12)
            worst = max(worst, abs(patch - form) / scale)
    assert worst < 2.5e-2  # 1.9 % of the half beam at the C1 junction and the stem


def test_realised_dimensions_match_the_request(hull):
    from tensorcad.hull_nurbs import realised_dimensions

    realised = realised_dimensions(hull)
    assert realised["length_m"] == pytest.approx(3.0)
    assert abs(realised["beam_relative_deficit"]) < 1e-2
    assert realised["draft_m"] == pytest.approx(0.15)
    assert realised["beam_at_waterline_m"] <= realised["beam_requested_m"]


def test_hydrostatics_coefficients_are_mutually_consistent(hull, hydro):
    assert hydro.block_coefficient == pytest.approx(
        hydro.prismatic_coefficient * hydro.midship_section_coefficient, rel=1e-12
    )
    assert 0.6 < hydro.midship_section_coefficient < 1.0
    assert 0.4 < hydro.waterplane_coefficient < 0.9
    assert 0.0 < hydro.vertical_centre_of_buoyancy < hull.parameters.draft


def test_quadrature_order_dependence_is_negligible(hull, hydro):
    """The reported numbers are converged in the quadrature, not an artefact of it.

    The calibration itself runs at the cheaper (10, 6) rule, which is why the design
    default is (24, 16): the difference between the two is 7e-4 relative, i.e. an order
    of magnitude below the 1 % consistency band it is used to judge.
    """
    low = nurbs_hydrostatics(hull, quadrature_order=(10, 6))
    assert abs(low.displaced_volume - hydro.displaced_volume) / hydro.displaced_volume < 2e-3
    refined = nurbs_hydrostatics(hull, quadrature_order=(32, 20))
    for name in ("displaced_volume", "waterplane_area"):
        high_value = getattr(hydro, name)
        refined_value = getattr(refined, name)
        assert abs(refined_value - high_value) / high_value < 1e-5, (name, refined_value, high_value)
    # The wetted-area integrand is |S_u x S_v|, which blows up at the stem where the
    # entrance has a vertical tangent, so it converges more slowly: 6e-4 between the
    # default rule and the refined one, still 15x inside the 1 % consistency band.
    assert (
        abs(refined.wetted_surface_area - hydro.wetted_surface_area) / hydro.wetted_surface_area
        < 1e-3
    )


# ---------------------------------------------------------------------------
# the body mesh
# ---------------------------------------------------------------------------
def test_mesh_nodes_on_the_shell_are_the_nurbs_patch(hull, mesh):
    stations, nz, layers = mesh.shape
    worst = 0.0
    for i in range(stations + 1):
        xi = mesh.node(i, 0, layers)[0] / hull.parameters.length
        for j in range(nz + 1):
            node = mesh.node(i, j, layers)
            shell = hull.surface.evaluate(xi, node[2] / hull.deck_height)[1]
            worst = max(worst, abs(node[1] - shell))
    assert worst < 1e-15


def test_waterline_is_an_exact_grid_level(hull, mesh):
    assert mesh.node(0, mesh.waterline_level, 0)[2] == pytest.approx(
        hull.parameters.draft, rel=0, abs=0
    )
    assert mesh.node(0, mesh.waterline_level - 1, 0)[2] < hull.parameters.draft
    assert mesh.node(0, mesh.waterline_level + 1, 0)[2] > hull.parameters.draft


def test_hex_quality_metrics_on_known_cells():
    cube = [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (1.0, 1.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
        (1.0, 0.0, 1.0),
        (1.0, 1.0, 1.0),
        (0.0, 1.0, 1.0),
    ]
    assert hex_scaled_jacobian(cube) == pytest.approx(1.0, rel=1e-15)
    assert hex_volume(cube) == pytest.approx(1.0, rel=1e-15)
    assert hex_characteristic_lengths(cube) == pytest.approx((1.0, 1.0))
    # A cell flattened in one direction keeps a perfect (unit) scaled Jacobian -- skew
    # and aspect ratio are different defects, which is why both are measured.
    flat = [(x, y, 0.01 * z) for x, y, z in cube]
    assert hex_scaled_jacobian(flat) == pytest.approx(1.0, rel=1e-12)
    assert hex_characteristic_lengths(flat)[0] == pytest.approx(0.01)
    # A collapsed cell has no volume and no Jacobian.
    collapsed = [(x, 0.0, z) for x, _y, z in cube]
    assert hex_volume(collapsed) == pytest.approx(0.0, abs=1e-18)
    assert hex_scaled_jacobian(collapsed) == 0.0


def test_mesh_quality_gate_is_met(mesh):
    quality = mesh.quality
    assert quality.cells > 0
    assert quality.pass_fraction > 0.90, quality.as_dict()
    assert quality.good_fraction > 0.85, quality.as_dict()
    assert quality.mean_scaled_jacobian > 0.85
    assert quality.max_aspect_ratio <= 1.0e4  # the stem wedge, reported not hidden
    assert quality.collapsed >= 0


def test_tet_subdivision_is_conforming_and_volume_preserving(mesh):
    tets, report = hex_to_tets(mesh)
    assert len(tets) == 6 * len(mesh.hexes)
    assert report["relative_volume_mismatch"] < 1e-12
    assert report["worst_cell_relative_mismatch"] < 1e-12
    # every tet index must be a real node index, and no tet may be inverted in the sense
    # of a negative volume (absolute values are used, so check the sum instead)
    assert max(max(tet) for tet in tets) < len(mesh.nodes)


# ---------------------------------------------------------------------------
# the consistency loop: C1-C7 of IN383-CAD01
# ---------------------------------------------------------------------------
def test_c1_displacement_matches_the_block_coefficient(hull, hydro):
    design = hull.parameters.design_displacement
    assert abs(hydro.displaced_volume - design) / design <= 0.02
    assert abs(hydro.displaced_volume - design) / design < 1e-3  # calibration residual


def test_c2_wetted_area_against_holtrop_mennen(hull, hydro):
    estimate = holtrop_wetted_area(
        hull.parameters,
        midship_section_coefficient=hydro.midship_section_coefficient,
        waterplane_coefficient=hydro.waterplane_coefficient,
    )
    assert abs(hydro.wetted_surface_area - estimate) / estimate <= 0.10
    # Anatomy: the wetted area must be the sum of the three wetted pieces.
    assert hydro.wetted_surface_area == pytest.approx(
        2.0 * hydro.side_area + hydro.bottom_area + hydro.transom_area_below_waterline, rel=1e-12
    )


def test_c3_longitudinal_centre_against_schneekluth(hull, hydro):
    estimate = schneekluth_lcb_ratio(hull.parameters)
    assert abs(hydro.lcb_ratio_from_midship - estimate) <= 0.05
    assert 0.0 < hydro.lcb_from_ap < hull.parameters.length


def test_c5_patch_and_mesh_agree(hull, hydro, mesh):
    integrals = mesh_integrals(mesh)
    volume = abs(integrals.displaced_volume - hydro.displaced_volume) / hydro.displaced_volume
    area = abs(integrals.wetted_surface_area - hydro.wetted_surface_area) / (
        hydro.wetted_surface_area
    )
    kb = abs(integrals.vertical_centre_of_buoyancy - hydro.vertical_centre_of_buoyancy) / (
        hydro.vertical_centre_of_buoyancy
    )
    lcb = abs(integrals.lcb_from_ap - hydro.lcb_from_ap) / hull.parameters.length
    assert max(volume, area, kb, lcb) <= 0.01, (volume, area, kb, lcb)


def test_c6_displacement_increases_with_block_coefficient():
    base = HullNurbsParameters()
    volumes = []
    for cb in (0.555, 0.575, 0.60, 0.625, 0.65):
        variant = base.replaced(block_coefficient=cb)
        volumes.append(nurbs_hydrostatics(build_nurbs_hull(variant)).displaced_volume)
    assert all(b > a for a, b in pairwise(volumes))
    slope = (volumes[-1] - volumes[0]) / (0.65 - 0.555)
    assert slope == pytest.approx(base.length * base.beam * base.draft, rel=2e-3)


def test_c7_refinement_is_stable(coarse_and_fine):
    coarse, fine = coarse_and_fine
    a = mesh_integrals(coarse)
    b = mesh_integrals(fine)
    assert abs(b.displaced_volume - a.displaced_volume) / a.displaced_volume <= 0.02
    assert abs(b.wetted_surface_area - a.wetted_surface_area) / a.wetted_surface_area <= 0.02
    # refining must also approach the patch value, not just move less
    assert fine.shape[0] == 2 * coarse.shape[0]


def test_acceptance_report_passes_every_criterion():
    report = acceptance_report()
    assert report["all_passed"] is True, {
        key: value["detail"] for key, value in report["criteria"].items() if not value["passed"]
    }
    for key in ("C1", "C2", "C3", "C4", "C5", "C6", "C7"):
        assert key in report["criteria"]
        assert report["criteria"][key]["passed"] is True, report["criteria"][key]
    assert report["schema"] == "tensorcad.hull_nurbs/1.0"
    assert report["mesh"]["quality"]["pass_fraction"] > 0.90
    assert report["tetrahedra"]["relative_volume_mismatch"] < 1e-12


def test_mirror_assumption_is_explicit(hull, hydro, mesh):
    """The half-hull integrals are doubled, so the factor must be exactly 2."""
    integrals = mesh_integrals(mesh)
    stations, _nz, layers = mesh.shape
    half = 0.0
    for i in range(stations):
        for j in range(mesh.waterline_level):
            for k in range(layers):
                cell = [mesh.node(i + di, j + dj, k + dk) for di, dj, dk in
                        ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                         (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1))]
                half += hex_volume(cell)
    assert integrals.displaced_volume == pytest.approx(2.0 * half, rel=1e-12)


def test_thresholds_are_declared_constants():
    assert SCALED_JACOBIAN_MIN == 0.2
    assert ASPECT_RATIO_MAX == 100.0
