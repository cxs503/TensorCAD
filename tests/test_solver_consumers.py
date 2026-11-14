"""Invariants of the solver-consumption layer.

The tests assert *refusals* and *invariants* -- a y-up bundle, a tampered table,
a surface with a hole, an occupancy cell outside the lattice -- and only then the
numbers, because the point of the layer is that a bad handoff cannot reach a
solver.  Solver-dependent assertions are skipped (honestly) when the library is
not importable on this host.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from tensorcad.solver_interface import (
    COORDINATE_CONTRACT,
    HandoffError,
    ellipsoid_mesh,
    grid_quads,
    read_grid_table,
    read_lattice_table,
    signed_volume,
    verify_bundle,
    write_bundle_grid,
    write_bundle_lattice,
    write_bundle_surface,
)
from tensorcad.solver_consumers import (
    bem_added_mass_from_bundle,
    fem_plate_deflection_from_bundle,
    gate,
    lbm_channel_from_bundle,
    watertight_report,
)


def _grid_xs_ys(side: float, n: int) -> tuple[list[float], list[float]]:
    return [side * i / n for i in range(n + 1)], [side * j / n for j in range(n + 1)]


def _channel_bundle(tmp_path: Path, nz: int = 4, nx: int = 8) -> Path:
    cells = [(i, 0, k) for i in range(nx) for k in range(nz)]
    write_bundle_lattice(
        tmp_path,
        shape=(nx, 1, nz),
        origin=(0.0, 0.0, 0.0),
        spacing=(1.0, float(nx), 1.0),
        fluid_cells=cells,
        walls={"normal_axis": "z", "lower_k": 0, "upper_k": nz - 1},
        case="unit_channel",
    )
    return tmp_path


def _sphere_bundle(tmp_path: Path, n_theta: int = 12, n_phi: int = 24) -> Path:
    write_bundle_surface(
        tmp_path, mesh=ellipsoid_mesh((0.5, 0.5, 0.5), n_theta=n_theta, n_phi=n_phi), case="unit_sphere"
    )
    return tmp_path


# ----------------------------------------------------------------- frame contract


def test_surface_bundle_re_reads_and_declares_z_up(tmp_path):
    bundle = _sphere_bundle(tmp_path)
    report = verify_bundle(bundle)
    assert report["ok"] is True
    assert report["route"] == "surface"
    manifest = json.loads((bundle / "handoff.json").read_text())
    assert manifest["frame"]["name"] == COORDINATE_CONTRACT["name"] == "z-up"
    assert manifest["signed_volume_m3"] > 0.0


def test_y_up_manifest_is_refused(tmp_path):
    bundle = _sphere_bundle(tmp_path)
    manifest = json.loads((bundle / "handoff.json").read_text())
    manifest["frame"]["name"] = "y-up"
    (bundle / "handoff.json").write_text(json.dumps(manifest))
    with pytest.raises(HandoffError):
        verify_bundle(bundle)


def test_tampered_surface_table_is_refused(tmp_path):
    """A moved node must fail the digest check, not reach a solver."""
    bundle = _sphere_bundle(tmp_path)
    table = json.loads((bundle / "surface.json").read_text())
    table["vertices"][0][0] += 1e-6
    (bundle / "surface.json").write_text(json.dumps(table))
    with pytest.raises(HandoffError):
        verify_bundle(bundle)
    with pytest.raises(HandoffError):
        bem_added_mass_from_bundle(bundle)  # refused before the solver is imported


def test_surface_hole_is_detected():
    mesh = ellipsoid_mesh((0.5, 0.5, 0.5), n_theta=12, n_phi=24)
    holed = mesh.faces[:-1]
    report = watertight_report(mesh.vertices, holed)
    assert report["closed"] is False and report["bad_edges"] == 3
    full = watertight_report(mesh.vertices, mesh.faces)
    assert full["closed"] and full["orientation_consistent"] and full["outward"]


def test_inward_surface_is_refused(tmp_path):
    mesh = ellipsoid_mesh((0.5, 0.5, 0.5), n_theta=12, n_phi=24)
    flipped = type(mesh)(
        vertices=mesh.vertices,
        faces=tuple((a, c, b) for a, b, c in mesh.faces),
        units=mesh.units,
        coordinate_frame=mesh.coordinate_frame,
    )
    assert signed_volume(flipped.vertices, flipped.faces) < 0.0
    with pytest.raises(HandoffError):
        write_bundle_surface(tmp_path, mesh=flipped, case="inward")


# ---------------------------------------------------------------- grid route


def test_grid_route_round_trip_and_winding(tmp_path):
    xs, ys = _grid_xs_ys(2.0, 4)
    manifest = write_bundle_grid(tmp_path, xs=xs, ys=ys, case="unit_plate")
    report = verify_bundle(tmp_path)
    assert report["route"] == "grid" and report["quads"] == 16
    assert math.isclose(report["surface_area_m2"], 4.0, rel_tol=1e-12)
    table = read_grid_table(tmp_path / "grid.json")
    quad = table["quads"][0]
    a, b, c = (table["vertices"][i] for i in quad[:3])
    cross = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    assert cross > 0.0, "quads must wind counter-clockwise about +z"
    assert len(grid_quads(4, 4)) == 16


def test_grid_tamper_is_refused(tmp_path):
    xs, ys = _grid_xs_ys(2.0, 3)
    write_bundle_grid(tmp_path, xs=xs, ys=ys, case="unit_plate")
    table = json.loads((tmp_path / "grid.json").read_text())
    table["vertices"][-1][1] += 0.5
    (tmp_path / "grid.json").write_text(json.dumps(table))
    with pytest.raises(HandoffError):
        verify_bundle(tmp_path)


def test_grid_rejects_non_increasing_coordinates(tmp_path):
    with pytest.raises(HandoffError):
        write_bundle_grid(tmp_path, xs=[0.0, 0.0, 1.0], ys=[0.0, 1.0], case="bad")


# ------------------------------------------------------------- lattice route


def test_lattice_round_trip(tmp_path):
    _channel_bundle(tmp_path, nz=4, nx=8)
    report = verify_bundle(tmp_path)
    assert report["route"] == "lattice"
    assert report["fluid_cells"] == 32
    assert math.isclose(report["fluid_volume_m3"], 32 * 8.0, rel_tol=1e-12)
    data = read_lattice_table(tmp_path / "lattice.json")
    assert data["walls"]["normal_axis"] == "z"


def test_lattice_rejects_cell_outside_box(tmp_path):
    with pytest.raises(HandoffError):
        write_bundle_lattice(
            tmp_path,
            shape=(2, 1, 2),
            origin=(0.0, 0.0, 0.0),
            spacing=(1.0, 1.0, 1.0),
            fluid_cells=[(2, 0, 0)],
            walls={"normal_axis": "z"},
            case="bad",
        )


def test_lattice_tamper_is_refused(tmp_path):
    _channel_bundle(tmp_path, nz=4, nx=8)
    table = json.loads((tmp_path / "lattice.json").read_text())
    table["fluid_cells"] = table["fluid_cells"][:-1]
    (tmp_path / "lattice.json").write_text(json.dumps(table))
    with pytest.raises(HandoffError):
        verify_bundle(tmp_path)


def test_gate_recomputes_the_error():
    passing = gate("g", 101.0, 100.0, 0.05)
    failing = gate("g", 110.0, 100.0, 0.05)
    assert passing["passed"] and not failing["passed"]
    assert math.isclose(failing["rel_error"], 0.1, rel_tol=1e-12)


# ------------------------------------------------------- solver consumption


def test_fem_plate_consumes_the_grid_file(tmp_path):
    pytest.importorskip("tensorfem.plate")
    xs, ys = _grid_xs_ys(2.0, 8)
    write_bundle_grid(tmp_path, xs=xs, ys=ys, case="plate8")
    result = fem_plate_deflection_from_bundle(
        tmp_path, side=2.0, young=7.0e10, poisson=0.3, thickness=0.05, q0=1.0e4
    )
    assert result.status == "ran"
    rigidity = 7.0e10 * 0.05**3 / (12.0 * (1.0 - 0.3**2))
    reference = 1.0e4 * 2.0**4 / (4.0 * math.pi**4 * rigidity)
    error = abs(result.metrics["centre_deflection_m"] - reference) / reference
    assert error < 0.05, f"plate deflection off by {error:.4%}"
    assert result.metrics["free_dof_residual_norm"] < 1e-6 * result.metrics["load_norm"]


def test_lbm_channel_consumes_the_lattice_file(tmp_path):
    pytest.importorskip("tensorlbm")
    _channel_bundle(tmp_path, nz=8, nx=16)
    result = lbm_channel_from_bundle(tmp_path, tau=0.8, body_force=2.0e-5, max_steps=4000)
    assert result.status == "ran"
    viscosity = (0.8 - 0.5) / 3.0
    reference = 2.0e-5 * 8**2 / (8.0 * viscosity)
    assert abs(result.metrics["u_max_lattice"] - reference) / reference < 0.05
    assert math.isclose(
        result.metrics["curvature_lattice"], 2.0e-5 / viscosity, rel_tol=1e-6
    ), "the steady profile must be parabolic (constant second difference)"
    assert result.metrics["mass_drift"] < 1e-9


def test_bem_sphere_consumes_the_surface_file(tmp_path):
    pytest.importorskip("tensorbem")
    _sphere_bundle(tmp_path, n_theta=16, n_phi=32)
    result = bem_added_mass_from_bundle(tmp_path, density=1000.0, center=(0.0, 0.0, -0.75))
    assert result.status == "ran"
    reference = 2.0 / 3.0 * math.pi * 1000.0 * 0.5**3
    assert abs(result.metrics["heave_kg"] - reference) / reference < 0.03
