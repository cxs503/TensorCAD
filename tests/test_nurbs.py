"""Acceptance tests for the NURBS parametric-surface service (campaign IN16-K5).

Every assertion is a machine-checkable criterion from the ticket, and the two
negative controls the ticket names are here as explicit rejection tests:

1. a circular arc and a rational cylinder reproduce their analytic length, area and
   volume (``< 1e-6`` relative, in fact round-off);
2. the analytic derivatives agree with central differences of the evaluation, so the
   derivative machinery is verified rather than trusted;
3. the parameter-derivative continuity across interior knots is *measured* from both
   sides, with a cubic ``C2`` control curve proving that a real jump is detected;
4. the tessellation is watertight, its faces all point outward, its normals converge
   to the analytic ones and its area/volume converge at the expected rate;
5. the STL round-trip keeps the shell closed;
6. **negative control 1**: a non-monotone knot vector is rejected;
7. **negative control 2**: a weight of zero (and a negative weight) is rejected;
8. the additive PC6 contract block keeps P1's schema and its parameter parser intact,
   and the L/B/Cb geometry gate rejects a tampered declaration.
"""
from __future__ import annotations

import math
from dataclasses import fields
from itertools import pairwise

import pytest

from tensorcad.export import read_stl, write_stl
from tensorcad.nurbs import (
    NURBS_SCHEMA,
    NURBSCurve,
    NURBSError,
    NURBSSurface,
    circle_curve,
    circular_arc,
    curve_derivative_continuity,
    cylinder_solid_mesh,
    cylinder_surface,
    gauss_legendre,
    mesh_solid_report,
    pc6_compat_block,
    tessellation_report,
    validate_knots,
    validate_weights,
)
from tensorcad.pc6 import PC6Parameters, analytic_form
from tensorcad.pc6_contract import CONTRACT_SCHEMA, _parameters_from_contract

RADIUS = 12.0
HEIGHT = 11.0


def _relative(measured: float, reference: float) -> float:
    return abs(measured - reference) / reference


def _cubic_spline() -> NURBSCurve:
    """A single-span-pair cubic: declared ``C2`` at the interior knot 0.4."""
    return NURBSCurve(
        degree=3,
        knots=(0.0, 0.0, 0.0, 0.0, 0.4, 1.0, 1.0, 1.0, 1.0),
        control_points=((0, 0, 0), (1, 2, 0), (2, 0, 0), (3, 3, 0), (4, 1, 0)),
    )


# ---------------------------------------------------------------------------
# 1. analytic anchors
# ---------------------------------------------------------------------------
def test_circular_arc_is_exact_not_approximated() -> None:
    """A rational quadratic arc lies *on* the circle: the deviation is round-off."""
    radius, sweep = 2.0, 1.0
    arc = circular_arc(radius, sweep)
    deviation = max(
        abs(math.hypot(point[0], point[1]) - radius) for point in arc.polyline(2001)
    )
    assert deviation < 1e-12, deviation
    # The length is r * theta by construction, so the relative error must be round-off
    # rather than the ~1e-3 an inscribed polyline of the same point count would give.
    assert _relative(arc.arc_length(), radius * sweep) < 1e-12


def test_full_circle_length_matches_analytic() -> None:
    circle = circle_curve(RADIUS)
    assert len(circle.control_points) == 7  # three 120 degree arcs: 3 * 2 + 1
    assert _relative(circle.arc_length(), 2.0 * math.pi * RADIUS) < 1e-12
    deviation = max(
        abs(math.hypot(point[0], point[1]) - RADIUS) for point in circle.polyline(2001)
    )
    assert deviation < 1e-12, deviation


def test_cylinder_lateral_area_matches_analytic_below_1e_minus_6() -> None:
    """The acceptance criterion: parametric area vs ``2 pi r h`` under 1e-6."""
    cases = ((RADIUS, HEIGHT), (1.0, 1.0), (0.5, 40.0))
    for radius, height in cases:
        surface = cylinder_surface(radius, height)
        reference = 2.0 * math.pi * radius * height
        assert _relative(surface.area(), reference) < 1e-6, (radius, height)
        # In practice the quadrature is exact to round-off, two orders tighter.
        assert _relative(surface.area(), reference) < 1e-12


def test_cylinder_area_converges_with_quadrature_order() -> None:
    surface = cylinder_surface(RADIUS, HEIGHT)
    rows = surface.area_convergence()["rows"]
    errors = [row["relative_error"] for row in rows]
    assert errors[-1] < 1e-12
    assert all(later <= earlier for earlier, later in pairwise(errors)), errors


def test_gauss_legendre_integrates_polynomials_exactly() -> None:
    """An n-point rule is exact for degree <= 2n - 1: check the identity it relies on."""
    nodes, weights = gauss_legendre(4)
    integral = sum(weight * node**6 for node, weight in zip(nodes, weights, strict=True))
    assert abs(integral - 2.0 / 7.0) < 1e-15
    assert abs(sum(weights) - 2.0) < 1e-14


# ---------------------------------------------------------------------------
# 2. derivatives verified against central differences
# ---------------------------------------------------------------------------
def test_curve_derivatives_match_central_differences() -> None:
    """First derivative on the rational arc, second on a C2 cubic (see the note below).

    The second derivative of the three-arc circle is deliberately *not* asked for: its
    knots carry multiplicity 2 with degree 2, so ``d2C/du2`` is a distribution at the
    joins and the derivative net is not representable there. :func:`derivatives` raises
    instead of returning a plausible-looking number, which
    :func:`test_second_derivative_of_a_c0_join_is_refused` pins down.
    """
    circle = circle_curve(RADIUS)
    target, step = 0.37, 2e-6
    analytic = circle.derivatives(target, 1)
    plus = circle.evaluate(target + step)
    minus = circle.evaluate(target - step)
    first = [(plus[axis] - minus[axis]) / (2 * step) for axis in range(3)]
    for axis in range(3):
        assert abs(analytic[1][axis] - first[axis]) < 1e-6 * max(1.0, abs(first[axis]))

    cubic = _cubic_spline()
    target, step = 0.21, 1e-5
    derivatives = cubic.derivatives(target, 2)
    value = cubic.evaluate(target)
    plus = cubic.evaluate(target + step)
    minus = cubic.evaluate(target - step)
    second = [(plus[axis] - 2 * value[axis] + minus[axis]) / step**2 for axis in range(3)]
    for axis in range(3):
        assert abs(derivatives[1][axis] - (plus[axis] - minus[axis]) / (2 * step)) < 1e-6
        assert abs(derivatives[2][axis] - second[axis]) < 1e-3 * max(1.0, abs(second[axis]))


def test_second_derivative_of_a_c0_join_is_refused() -> None:
    """Documented limitation: the request raises instead of inventing a value."""
    with pytest.raises(NURBSError, match="zero-length knot interval"):
        circle_curve(RADIUS).derivatives(0.37, 2)


def _quadratic_patch() -> NURBSSurface:
    """Single-span degree (2, 2) patch: every partial derivative up to order 2 exists."""
    net = tuple(
        tuple((float(i), float(j), 0.3 * i * i * j + 0.1 * j * j) for j in range(3))
        for i in range(3)
    )
    return NURBSSurface(
        degree=(2, 2),
        knots=((0.0, 0.0, 0.0, 1.0, 1.0, 1.0), (0.0, 0.0, 0.0, 1.0, 1.0, 1.0)),
        control_net=net,
    )


def test_surface_second_derivatives_match_central_differences() -> None:
    surface = _quadratic_patch()
    u, v, step = 0.37, 0.61, 1e-5
    ders = surface.derivatives(u, v, 2)
    centre = surface.evaluate(u, v)
    plus_u, minus_u = surface.evaluate(u + step, v), surface.evaluate(u - step, v)
    plus_v, minus_v = surface.evaluate(u, v + step), surface.evaluate(u, v - step)
    for axis in range(3):
        fd_uu = (plus_u[axis] - 2 * centre[axis] + minus_u[axis]) / step**2
        fd_vv = (plus_v[axis] - 2 * centre[axis] + minus_v[axis]) / step**2
        assert abs(ders[(2, 0)][axis] - fd_uu) < 1e-3 * max(1.0, abs(fd_uu))
        assert abs(ders[(0, 2)][axis] - fd_vv) < 1e-3 * max(1.0, abs(fd_vv))


def test_surface_partials_match_central_differences() -> None:
    surface = cylinder_surface(RADIUS, HEIGHT)
    u, v, step = 0.61, 0.43, 1e-6
    ders = surface.derivatives(u, v, 1)
    plus_u = surface.evaluate(u + step, v)
    minus_u = surface.evaluate(u - step, v)
    plus_v = surface.evaluate(u, v + step)
    minus_v = surface.evaluate(u, v - step)
    for axis in range(3):
        fd_u = (plus_u[axis] - minus_u[axis]) / (2 * step)
        fd_v = (plus_v[axis] - minus_v[axis]) / (2 * step)
        assert abs(ders[(1, 0)][axis] - fd_u) < 1e-5 * max(1.0, abs(fd_u))
        assert abs(ders[(0, 1)][axis] - fd_v) < 1e-5 * max(1.0, abs(fd_v))
    # On the exact cylinder dS/dv must be the unit height vector, to round-off.
    assert abs(ders[(0, 1)][0]) < 1e-12
    assert abs(ders[(0, 1)][1]) < 1e-12
    assert abs(ders[(0, 1)][2] - HEIGHT) < 1e-9


# ---------------------------------------------------------------------------
# 3. parameter-derivative continuity
# ---------------------------------------------------------------------------
def test_circle_joins_are_measured_c1_although_declared_c0() -> None:
    circle = circle_curve(RADIUS)
    records = curve_derivative_continuity(circle, order=2)
    assert len(records) == 2
    for record in records:
        assert record.declared_order == 0  # multiplicity 2 with degree 2
        assert record.max_measurable_order == 1
        assert record.measured_continuous_order == 1
        # |C'(t+0) - C'(t-0)| relative to |C'| ~ 62
        assert record.jumps[1] < 1e-9 * record.scales[1]
        assert record.passed


def test_cubic_control_curve_detects_a_real_derivative_jump() -> None:
    """Positive control: the measurement must *find* a discontinuity, not hide it."""
    cubic = _cubic_spline()
    record = curve_derivative_continuity(cubic, order=3)[0]
    assert record.declared_order == 2
    assert record.jumps[0] == 0.0 and record.jumps[1] == 0.0 and record.jumps[2] == 0.0
    assert record.jumps[3] > 1e-3 * record.scales[3], record.jumps
    assert record.measured_continuous_order == 2


def test_cylinder_surface_continuity_across_u_knots() -> None:
    surface = cylinder_surface(RADIUS, HEIGHT)
    report = surface.continuity(1)
    assert report["passed"]
    assert len(report["u"]) == 2
    for record in report["u"]:
        assert record["jumps"][1] < 1e-9 * record["scales"][1]
    assert report["v"] == []  # the linear direction has no interior knot


# ---------------------------------------------------------------------------
# 4. tessellation: topology, normals, convergence
# ---------------------------------------------------------------------------
def test_lateral_tessellation_welds_its_seam_and_closes_where_it_should() -> None:
    """The revolution seam is welded; only the two rims are open.

    ``2 * u_samples`` boundary edges is the structural statement: the two seam columns
    ``u = 0`` and ``u = 1`` evaluate to the same positions, so :class:`MeshBuilder`
    collapses them into one vertex column and the tube has exactly two open rims. A
    seam that failed to weld would show four.
    """
    surface = cylinder_surface(RADIUS, HEIGHT)
    u_samples = 96
    report = tessellation_report(surface, surface.tessellate(u_samples, 8))
    assert not report["watertight"]  # an open tube, stated rather than implied
    assert report["boundary_edges"] == 2 * u_samples, report["boundary_edges"]
    assert report["edge_manifold"]
    assert report["non_manifold_edges"] == 0
    assert report["components"] == 1
    assert report["inward_faces"] == 0
    assert report["orientation_flips"] == 0
    assert report["max_normal_angle_deg"] < 1.0, report["max_normal_angle_deg"]


def test_tessellated_area_and_normal_error_converge_at_second_order() -> None:
    surface = cylinder_surface(RADIUS, HEIGHT)
    coarse = tessellation_report(surface, surface.tessellate(24, 8))
    fine = tessellation_report(surface, surface.tessellate(48, 8))
    ratio = abs(coarse["area_relative_error"]) / abs(fine["area_relative_error"])
    assert 3.0 < ratio < 5.0, ratio  # O(h^2) -> a factor of ~4 per doubling
    assert abs(fine["area_relative_error"]) < 2e-3
    assert fine["max_normal_angle_deg"] < coarse["max_normal_angle_deg"]
    # The tessellated area is inscribed, so it must stay below the parametric value.
    assert fine["area_relative_error"] < 0.0


def test_closed_solid_is_watertight_and_matches_analytic_volume() -> None:
    volume_reference = math.pi * RADIUS * RADIUS * HEIGHT
    solid = cylinder_solid_mesh(RADIUS, HEIGHT, u_samples=192, v_samples=6)
    report = mesh_solid_report(solid, volume_reference=volume_reference)
    assert report["watertight"]
    assert report["components"] == 1
    assert report["boundary_edges"] == 0
    assert report["orientation_consistent"]
    assert report["signed_volume_m3"] > 0.0
    assert abs(report["volume_relative_error"]) < 1e-3, report["volume_relative_error"]


def test_stl_round_trip_keeps_the_shell_closed(tmp_path) -> None:
    volume_reference = math.pi * RADIUS * RADIUS * HEIGHT
    solid = cylinder_solid_mesh(RADIUS, HEIGHT, u_samples=96, v_samples=4)
    path = write_stl(solid, tmp_path / "cylinder.stl")
    recovered = read_stl(path)
    report = mesh_solid_report(recovered, volume_reference=volume_reference)
    assert report["watertight"]
    assert report["components"] == 1
    # float32 STL versus float64 mesh plus inscribed facets: still within 1e-3.
    assert abs(report["volume_relative_error"]) < 1e-3


# ---------------------------------------------------------------------------
# 5. negative controls
# ---------------------------------------------------------------------------
def test_non_monotone_knot_vector_is_rejected() -> None:
    with pytest.raises(NURBSError, match="non-decreasing"):
        NURBSCurve(
            degree=2,
            knots=(0.0, 0.0, 0.0, 0.6, 0.4, 1.0, 1.0),
            control_points=((0, 0, 0), (1, 1, 0), (2, 0, 0), (3, 1, 0)),
        )
    with pytest.raises(NURBSError, match="non-decreasing"):
        validate_knots((0.0, 0.0, 0.0, 0.5, 0.5, 0.4, 1.0, 1.0, 1.0), 2, 6)


@pytest.mark.parametrize("bad_weight", [0.0, -0.5, -1e-9])
def test_non_positive_weights_are_rejected(bad_weight: float) -> None:
    circle = circle_curve(1.0)
    weights = list(circle.weights)
    weights[1] = bad_weight
    with pytest.raises(NURBSError, match="strictly positive"):
        NURBSCurve(
            degree=circle.degree,
            knots=circle.knots,
            control_points=circle.control_points,
            weights=tuple(weights),
        )
    with pytest.raises(NURBSError, match="strictly positive"):
        validate_weights(weights, len(weights))


def test_parameter_outside_the_domain_is_rejected() -> None:
    circle = circle_curve(RADIUS)
    with pytest.raises(NURBSError, match="outside the domain"):
        circle.evaluate(1.5)
    with pytest.raises(NURBSError, match="outside the domain"):
        cylinder_surface(RADIUS, HEIGHT).evaluate(-0.1, 0.5)


def test_malformed_knot_vectors_are_rejected() -> None:
    with pytest.raises(NURBSError, match="needs"):
        validate_knots((0.0, 0.0, 0.0, 1.0, 1.0), 2, 3)  # one knot short
    with pytest.raises(NURBSError, match="more than"):
        validate_knots((0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0), 2, 7)
    with pytest.raises(NURBSError, match="not clamped"):
        validate_knots((-0.5, 0.0, 0.0, 1.0, 1.0, 1.0), 2, 3)
    with pytest.raises(NURBSError, match="degree"):
        validate_knots((0.0, 0.0, 1.0, 1.0), 0, 3)


def test_control_net_layout_mismatch_is_rejected() -> None:
    circle = circle_curve(1.0)
    with pytest.raises(NURBSError, match="layout"):
        NURBSSurface(
            degree=(2, 1),
            knots=(circle.knots, (0.0, 0.0, 1.0, 1.0)),
            control_net=tuple((point, point) for point in circle.control_points),
            weights=tuple((1.0,) for _ in circle.control_points),
        )


# ---------------------------------------------------------------------------
# 6. PC6 contract compatibility
# ---------------------------------------------------------------------------
def _parameter_table() -> dict[str, dict[str, object]]:
    return {
        field.name: {"value": getattr(PC6Parameters(), field.name)}
        for field in fields(PC6Parameters())
    }


def test_pc6_block_keeps_the_contract_schema_and_parser_intact() -> None:
    surface = cylinder_surface(RADIUS, HEIGHT)
    block = pc6_compat_block(surface, tessellation=surface.tessellate(48, 6))
    assert block["schema"] == CONTRACT_SCHEMA
    assert block["representation"] == NURBS_SCHEMA
    assert block["additive"] is True
    assert block["geometry_gate"]["passed"]

    contract = {"schema": CONTRACT_SCHEMA, "parameters": _parameter_table()}
    before = _parameters_from_contract(contract)
    augmented = dict(contract, nurbs_surface=block)
    after = _parameters_from_contract(augmented)
    assert before == after  # the block is additive, not a schema change


def test_pc6_geometry_gate_rejects_a_tampered_block_coefficient() -> None:
    surface = cylinder_surface(RADIUS, HEIGHT)
    form = analytic_form(PC6Parameters())
    with pytest.raises(NURBSError, match="deviates"):
        pc6_compat_block(
            surface,
            declared={"block_coefficient": form.block_coefficient * 0.9},
        )
    with pytest.raises(NURBSError, match="not part of the gate"):
        pc6_compat_block(surface, declared={"length_overall": 200.0})


def test_pc6_geometry_gate_reports_L_B_and_Cb() -> None:
    surface = cylinder_surface(RADIUS, HEIGHT)
    gate = pc6_compat_block(
        surface, declared={"length_m": 200.0, "beam_m": 24.0}
    )["geometry_gate"]
    assert gate["values"]["length_m"] == pytest.approx(200.0)
    assert gate["values"]["beam_m"] == pytest.approx(24.0)
    assert 0.7 < gate["values"]["block_coefficient"] < 0.85
    assert gate["passed"]
