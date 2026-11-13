"""Tests for tensorcad.kovasznay — exact solution, mesh family, order harness."""
from __future__ import annotations

import json
import math
from dataclasses import FrozenInstanceError

import pytest

from tensorcad.kovasznay import (
    KOVASZNAY_1948,
    KOVASZNAY_ORDER_ACCEPTANCE,
    KovasznayFlow,
    StructuredMesh,
    case_spec,
    cell_averages,
    cell_centre_fields,
    convergence_study,
    fit_order,
    generate_structured_mesh,
    l2_norm,
    lam_for_reynolds,
    mesh_summary,
    mesh_to_json,
    midpoint_consistency_study,
    nodal_fields,
    probe_table,
    truncation_error_study,
    write_case_bundle,
    write_mesh_json,
)

LEVELS = (16, 32, 64)


# ---------------------------------------------------------------------------
# Parameters and the exact solution
# ---------------------------------------------------------------------------
class TestParameters:
    def test_reference_provenance(self):
        assert "Kovasznay" in KOVASZNAY_1948["citation"]
        assert KOVASZNAY_1948["doi"] == "10.1017/S0305004100023999"

    def test_lam_is_negative_root_of_quadratic(self):
        lam = lam_for_reynolds(40.0)
        residual = lam * lam - 40.0 * lam - 4.0 * math.pi ** 2
        assert abs(residual) < 1e-12
        assert lam < 0.0

    def test_lam_re40_value(self):
        # Generated from the closed form, not transcribed: -(sqrt(1600+16 pi^2) - 40)/2
        expected = 0.5 * (40.0 - math.sqrt(40.0 ** 2 + 16.0 * math.pi ** 2))
        assert abs(lam_for_reynolds(40.0) - expected) < 1e-15
        assert abs(lam_for_reynolds(40.0) + 0.963740544) < 5e-10

    def test_lam_approaches_zero_from_below_as_reynolds_grows(self):
        # lam = (Re - sqrt(Re^2+16pi^2))/2 -> -4*pi^2/Re, so it rises towards 0.
        assert lam_for_reynolds(40.0) < lam_for_reynolds(100.0) < 0.0

    def test_lam_rejects_non_positive(self):
        with pytest.raises(ValueError):
            lam_for_reynolds(0.0)

    def test_flow_defaults(self):
        flow = KovasznayFlow()
        assert flow.reynolds_number == 40.0
        assert (flow.x_min, flow.x_max) == (-0.5, 1.5)
        assert (flow.y_min, flow.y_max) == (-0.5, 1.5)
        assert flow.width == 2.0 and flow.height == 2.0

    def test_flow_nu(self):
        assert abs(KovasznayFlow(reynolds_number=40.0).nu - 0.025) < 1e-15

    def test_flow_rejects_bad_inputs(self):
        with pytest.raises(ValueError):
            KovasznayFlow(reynolds_number=-1.0)
        with pytest.raises(ValueError):
            KovasznayFlow(x_min=1.0, x_max=0.0)
        with pytest.raises(ValueError):
            KovasznayFlow(y_min=1.0, y_max=0.0)

    def test_as_dict_round_trips_through_json(self):
        payload = json.loads(json.dumps(KovasznayFlow().as_dict()))
        assert payload["model"] == "kovasznay_1948"
        assert payload["lam"] < 0.0
        assert payload["domain"]["x"] == [-0.5, 1.5]


class TestExactSolution:
    def test_velocity_at_origin_is_stagnation_point(self):
        u, v = KovasznayFlow().velocity(0.0, 0.0)
        assert abs(u) < 1e-15 and abs(v) < 1e-15

    def test_velocity_at_centreline(self):
        u, v = KovasznayFlow().velocity(0.0, 0.5)
        assert abs(u - 2.0) < 1e-14
        assert abs(v) < 1e-15

    def test_pressure_at_x_zero(self):
        assert abs(KovasznayFlow().pressure(0.0, 0.0) + 0.5) < 1e-15

    def test_pressure_offset_applied(self):
        flow = KovasznayFlow(pressure_offset=3.0)
        assert abs(flow.pressure(0.5, 0.25) - (3.0 - 0.5 * math.exp(2 * flow.lam * 0.5))) < 1e-14

    def test_transverse_velocity_vanishes_on_y_walls_and_centreline(self):
        flow = KovasznayFlow()
        for y in (-0.5, 0.0, 0.5):
            for x in (-0.5, 0.0, 0.7, 1.5):
                assert abs(flow.velocity(x, y)[1]) < 1e-14

    def test_continuity_is_zero(self):
        flow = KovasznayFlow()
        for x in (-0.5, -0.1, 0.0, 0.4, 1.5):
            for y in (-0.5, -0.2, 0.0, 0.3, 0.5):
                assert abs(flow.ns_residual(x, y)[2]) < 1e-15

    def test_momentum_residual_is_roundoff(self):
        flow = KovasznayFlow()
        maxima = [0.0, 0.0]
        for x in (-0.5, -0.2, 0.0, 0.35, 1.5):
            for y in (-0.5, -0.1, 0.0, 0.25, 0.5):
                ru, rv, _ = flow.ns_residual(x, y)
                maxima[0] = max(maxima[0], abs(ru))
                maxima[1] = max(maxima[1], abs(rv))
        assert maxima[0] < 1e-12
        assert maxima[1] < 1e-12

    def test_residual_small_at_other_reynolds(self):
        flow = KovasznayFlow(reynolds_number=100.0)
        for x in (-0.5, 0.0, 0.9):
            for y in (-0.4, 0.0, 0.4):
                ru, rv, div = flow.ns_residual(x, y)
                assert max(abs(ru), abs(rv), abs(div)) < 1e-11

    def test_vorticity_routes_agree(self):
        flow = KovasznayFlow()
        for x in (-0.5, 0.0, 0.6):
            for y in (-0.45, -0.1, 0.2, 0.5):
                assert flow.vorticity_mismatch(x, y) < 1e-13

    @pytest.mark.parametrize("reynolds", [10.0, 40.0, 100.0])
    def test_gradient_matches_finite_difference(self, reynolds):
        flow = KovasznayFlow(reynolds_number=reynolds)
        h = 1e-5
        x, y = 0.13, -0.21
        u_x, u_y, v_x, v_y = flow.velocity_gradient(x, y)
        fd_u_x = (flow.velocity(x + h, y)[0] - flow.velocity(x - h, y)[0]) / (2 * h)
        fd_u_y = (flow.velocity(x, y + h)[0] - flow.velocity(x, y - h)[0]) / (2 * h)
        fd_v_x = (flow.velocity(x + h, y)[1] - flow.velocity(x - h, y)[1]) / (2 * h)
        fd_v_y = (flow.velocity(x, y + h)[1] - flow.velocity(x, y - h)[1]) / (2 * h)
        assert abs(u_x - fd_u_x) < 1e-7
        assert abs(u_y - fd_u_y) < 1e-7
        assert abs(v_x - fd_v_x) < 1e-7
        assert abs(v_y - fd_v_y) < 1e-7

    def test_streamfunction_derivatives_reproduce_velocity(self):
        flow = KovasznayFlow()
        h = 1e-6
        for x, y in ((0.0, 0.2), (-0.3, -0.4), (0.9, 0.35)):
            u = (flow.streamfunction(x, y + h) - flow.streamfunction(x, y - h)) / (2 * h)
            v = -(flow.streamfunction(x + h, y) - flow.streamfunction(x - h, y)) / (2 * h)
            ue, ve = flow.velocity(x, y)
            assert abs(u - ue) < 1e-7
            assert abs(v - ve) < 1e-7

    def test_speed(self):
        flow = KovasznayFlow()
        u, v = flow.velocity(0.2, 0.1)
        assert abs(flow.speed(0.2, 0.1) - math.hypot(u, v)) < 1e-15


# ---------------------------------------------------------------------------
# Mesh family
# ---------------------------------------------------------------------------
class TestMesh:
    def test_cell_and_node_counts(self):
        mesh = generate_structured_mesh(KovasznayFlow(), 8)
        assert mesh.nx == 8 and mesh.ny == 8
        assert mesh.n_cells == 64
        assert mesh.n_nodes == 81
        assert len(mesh.quads) == 64

    def test_rectangular_mesh(self):
        mesh = generate_structured_mesh(KovasznayFlow(), 6, 3)
        # hx = 2/6, hy = 2/3 -> the longest/shorter edge ratio is 2.
        assert mesh.summary()["aspect_ratio"] == pytest.approx(2.0)
        assert mesh.n_cells == 18

    def test_exact_area_and_positive_jacobians(self):
        mesh = generate_structured_mesh(KovasznayFlow(), 32)
        summary = mesh.summary()
        assert summary["total_area"] == pytest.approx(4.0, abs=1e-12)
        assert abs(summary["area_error"]) < 1e-12
        assert summary["negative_jacobians"] == 0
        assert mesh_summary(mesh)["n_cells"] == 1024

    def test_node_ordering(self):
        mesh = generate_structured_mesh(KovasznayFlow(), 4)
        assert mesh.node(0, 0) == (-0.5, -0.5)
        assert mesh.node(4, 4) == pytest.approx((1.5, 1.5))
        assert mesh.node_index(0, 0) == 0
        assert mesh.node_index(1, 0) == 1
        assert mesh.node_index(0, 1) == 5

    def test_quad_is_counter_clockwise(self):
        mesh = generate_structured_mesh(KovasznayFlow(), 4)
        assert all(area > 0.0 for area in mesh.areas())

    def test_boundary_sets(self):
        mesh = generate_structured_mesh(KovasznayFlow(), 4)
        sets = mesh.boundary_sets()
        assert len(sets["left"]) == 5 and len(sets["right"]) == 5
        assert len(sets["bottom"]) == 5 and len(sets["top"]) == 5
        assert sets["left"][0] == mesh.node_index(0, 0)
        assert sets["top"][-1] == mesh.node_index(4, 4)

    def test_centroids_are_inside_cells(self):
        mesh = generate_structured_mesh(KovasznayFlow(), 4)
        hx, hy = mesh.hx, mesh.hy
        first = mesh.centroids()[0]
        assert first == pytest.approx((-0.5 + hx / 2, -0.5 + hy / 2))

    def test_rejects_zero_cells(self):
        with pytest.raises(ValueError):
            generate_structured_mesh(KovasznayFlow(), 0)

    def test_fields_have_cell_and_node_length(self):
        mesh = generate_structured_mesh(KovasznayFlow(), 5)
        centre = cell_centre_fields(mesh)
        nodes = nodal_fields(mesh)
        assert len(centre["u"]) == mesh.n_cells
        assert len(nodes["u"]) == mesh.n_nodes
        assert set(centre) == {"u", "v", "p"}


# ---------------------------------------------------------------------------
# Norms and order fitting
# ---------------------------------------------------------------------------
class TestNormsAndFitting:
    def test_l2_norm_unweighted(self):
        assert l2_norm([1.0, 1.0]) == pytest.approx(1.0)
        assert l2_norm([3.0, 4.0]) == pytest.approx(math.sqrt(12.5))

    def test_l2_norm_weighted(self):
        # Weighted norm of a constant is the constant.
        assert l2_norm([2.0, 2.0], [1.0, 3.0]) == pytest.approx(2.0)

    def test_l2_norm_rejects_bad_input(self):
        with pytest.raises(ValueError):
            l2_norm([])
        with pytest.raises(ValueError):
            l2_norm([1.0, 2.0], [1.0])
        with pytest.raises(ValueError):
            l2_norm([1.0, 2.0], [0.0, 0.0])

    def test_fit_order_recovers_exact_power(self):
        hs = [0.1, 0.05, 0.025]
        errors = [3.0 * h ** 2.5 for h in hs]
        fit = fit_order(hs, errors)
        assert fit["order"] == pytest.approx(2.5)
        assert fit["r2"] == pytest.approx(1.0)

    def test_fit_order_recovers_first_order(self):
        hs = [0.2, 0.1, 0.05, 0.025]
        fit = fit_order(hs, [0.5 * h for h in hs])
        assert fit["order"] == pytest.approx(1.0)

    def test_fit_order_rejects_degenerate_input(self):
        with pytest.raises(ValueError):
            fit_order([0.1], [0.01])
        with pytest.raises(ValueError):
            fit_order([0.1, 0.05], [1.0])
        with pytest.raises(ValueError):
            fit_order([0.1, 0.05], [0.0, 1e-3])
        with pytest.raises(ValueError):
            fit_order([0.1, 0.1], [0.01, 0.005])


# ---------------------------------------------------------------------------
# Consistency orders
# ---------------------------------------------------------------------------
class TestConsistency:
    def test_cell_average_close_to_centroid_value(self):
        flow = KovasznayFlow()
        mesh = generate_structured_mesh(flow, 16)
        avg = cell_averages(mesh)
        point = cell_centre_fields(mesh)
        worst = max(abs(a - b) for a, b in zip(avg["u"], point["u"], strict=False))
        assert worst < 5e-2

    def test_midpoint_consistency_is_second_order(self):
        study = midpoint_consistency_study(KovasznayFlow(), LEVELS)
        assert study["order_u"] == pytest.approx(2.0, abs=0.1)
        assert study["order_v"] == pytest.approx(2.0, abs=0.1)
        assert study["order_velocity"] >= 1.8

    def test_midpoint_errors_decrease_monotonically(self):
        study = midpoint_consistency_study(KovasznayFlow(), LEVELS)
        errors = [row["l2_velocity"] for row in study["levels"]]
        assert errors[0] > errors[1] > errors[2]

    def test_truncation_error_is_second_order(self):
        study = truncation_error_study(KovasznayFlow(), LEVELS)
        for key in ("order_momentum_x", "order_momentum_y", "order_continuity"):
            assert study[key] == pytest.approx(2.0, abs=0.15)

    def test_truncation_residual_shrinks_by_four(self):
        study = truncation_error_study(KovasznayFlow(), (32, 64))
        coarse, fine = study["levels"]
        ratio = coarse["l2_momentum_x"] / fine["l2_momentum_x"]
        assert 3.5 < ratio < 4.5


# ---------------------------------------------------------------------------
# Solver-agnostic order harness
# ---------------------------------------------------------------------------
def _injected(flow: KovasznayFlow, order: int, amplitude: float = 0.37):
    def solver(mesh):
        h = mesh.hx
        out = []
        for x, y in mesh.centroids():
            u, v = flow.velocity(x, y)
            out.append((u - amplitude * h ** order * math.sin(3.0 * x) * math.cos(2.0 * y),
                        v + 0.57 * amplitude * h ** order * math.cos(x) * math.sin(4.0 * y)))
        return out

    return solver


class TestConvergenceHarness:
    def test_second_order_solver_is_recovered(self):
        flow = KovasznayFlow()
        study = convergence_study(flow, LEVELS, _injected(flow, 2))
        assert study["order"] == pytest.approx(2.0, abs=0.05)
        assert study["passes_acceptance"] is True
        assert len(study["levels"]) == 3

    def test_first_order_solver_is_rejected(self):
        flow = KovasznayFlow()
        study = convergence_study(flow, LEVELS, _injected(flow, 1))
        assert study["order"] == pytest.approx(1.0, abs=0.05)
        assert study["passes_acceptance"] is False

    def test_third_order_solver_is_recovered(self):
        flow = KovasznayFlow()
        study = convergence_study(flow, LEVELS, _injected(flow, 3))
        assert study["order"] == pytest.approx(3.0, abs=0.05)

    def test_exact_solver_reported_as_roundoff(self):
        flow = KovasznayFlow()
        study = convergence_study(
            flow, LEVELS, lambda mesh: [flow.velocity(x, y) for x, y in mesh.centroids()])
        assert study["order"] is None
        assert "round-off" in study["note"]

    def test_callback_is_required(self):
        with pytest.raises(ValueError):
            convergence_study(KovasznayFlow(), LEVELS)

    def test_wrong_length_is_rejected(self):
        with pytest.raises(ValueError):
            convergence_study(KovasznayFlow(), LEVELS, lambda mesh: [(0.0, 0.0)])

    def test_acceptance_threshold_is_documented(self):
        assert KOVASZNAY_ORDER_ACCEPTANCE["min_acceptable_order"] == 1.8
        assert KOVASZNAY_ORDER_ACCEPTANCE["levels"] == [32, 64, 128]


# ---------------------------------------------------------------------------
# Probe table, case spec and serialisation
# ---------------------------------------------------------------------------
class TestReferenceData:
    def test_probe_table_shape(self):
        rows = probe_table(KovasznayFlow())
        assert len(rows) == 7 * 5
        assert set(rows[0]) >= {"x", "y", "u", "v", "p", "vorticity"}
        assert max(abs(r["ns_residual_u"]) for r in rows) < 1e-12

    def test_probe_table_custom_lattice(self):
        rows = probe_table(KovasznayFlow(), xs=(0.0,), ys=(0.0, -0.5))
        assert len(rows) == 2
        assert rows[0]["u"] == pytest.approx(0.0, abs=1e-15)

    def test_case_spec_contents(self):
        spec = case_spec(KovasznayFlow())
        assert spec["case_code"] == "FVM-KOVASZNAY-ORDER"
        assert spec["acceptance"]["min_acceptable_order"] == 1.8
        assert len(spec["meshes"]) == 3
        assert spec["meshes"][0]["n_cells"] == 32 * 32
        assert set(spec["boundary_conditions"]) == {"left", "right", "bottom", "top"}

    def test_mesh_to_json_is_json_serialisable(self):
        payload = mesh_to_json(generate_structured_mesh(KovasznayFlow(), 4))
        text = json.dumps(payload)
        assert len(payload["nodes"]) == 25
        assert len(payload["quads"]) == 16
        assert len(payload["cell_centre_exact"]["u"]) == 16
        assert len(payload["nodal_exact"]["u"]) == 25
        assert "summary" in json.loads(text)

    def test_mesh_to_json_without_fields(self):
        payload = mesh_to_json(generate_structured_mesh(KovasznayFlow(), 4),
                              include_fields=False)
        assert "cell_centre_exact" not in payload


class TestWriters:
    def test_write_mesh_json(self, tmp_path):
        path = write_mesh_json(generate_structured_mesh(KovasznayFlow(), 4),
                               tmp_path / "mesh.json")
        payload = json.loads(path.read_text())
        assert payload["summary"]["n_cells"] == 16
        assert len(payload["boundary_sets"]["left"]) == 5

    def test_write_case_bundle(self, tmp_path):
        manifest = write_case_bundle(KovasznayFlow(), tmp_path / "bundle", (16, 32))
        assert set(manifest["files"]) == {
            "case_spec.json", "mesh_16.json", "mesh_32.json", "probe_table.json"}
        for name, digest in manifest["files"].items():
            path = tmp_path / "bundle" / name
            assert path.exists()
            assert len(digest) == 64
            json.loads(path.read_text())

    def test_bundle_is_deterministic(self, tmp_path):
        first = write_case_bundle(KovasznayFlow(), tmp_path / "a", (16,))
        second = write_case_bundle(KovasznayFlow(), tmp_path / "b", (16,))
        assert first["files"] == second["files"]


def test_structured_mesh_is_frozen():
    mesh = generate_structured_mesh(KovasznayFlow(), 2)
    assert isinstance(mesh, StructuredMesh)
    with pytest.raises(FrozenInstanceError):
        mesh.nx = 10  # type: ignore[misc]
