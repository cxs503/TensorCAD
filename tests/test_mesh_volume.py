"""Closed-body volume/watertightness checks against analytic solids.

The tests are the same claims the benchmark case makes, at a size a test suite can
afford: the references are closed forms, the box is exact, the icosphere ladder has
the order a piecewise-linear surface gives, and every injected defect is named.
Nothing here compares one mesh estimate with another mesh estimate.
"""
from __future__ import annotations

import importlib.util
import itertools
import math
from pathlib import Path

import numpy as np
import pytest

from tensorcad import mesh_volume as mv
from tensorcad.geometry_audit import audit_watertightness
from tensorcad.mesh_xval import icosphere

_REPO_ROOT = Path(__file__).resolve().parents[1]


# ---------------------------------------------------------------------------
# references
# ---------------------------------------------------------------------------
def test_sphere_reference_has_two_derivations():
    for radius in (0.37, 1.0, 2.5):
        direct = mv.analytic_sphere_volume(radius)
        cylinder = mv.analytic_sphere_volume_via_cylinder(radius)
        assert direct == pytest.approx(cylinder, rel=1e-15)
        # and the area reference is the derivative of the volume in r
        assert mv.analytic_sphere_area(radius) == pytest.approx(4.0 * math.pi * radius**2, rel=1e-15)


def test_icosahedron_closed_form_matches_the_divergence_route():
    """20 skewed triangles: the integrator is checked against a regular-polyhedron formula."""
    mesh = icosphere(1.0, 0)
    reference = mv.analytic_icosahedron_volume(1.0)
    assert mesh.volume() == pytest.approx(reference, rel=1e-12)
    assert mv.volume_by_tet_centroid(mesh) == pytest.approx(reference, rel=1e-12)


# ---------------------------------------------------------------------------
# bodies a triangle mesh can represent exactly
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("dims", [(2.0, 2.0, 2.0), (1.5, 0.75, 3.0), (0.2, 1.7, 3.9)])
def test_box_volume_and_area_are_exact(dims):
    mesh = mv.box_mesh(*dims)
    check = mv.check_closed_body(
        mesh,
        name="box",
        family="box",
        reference_volume=mv.analytic_box_volume(*dims),
        reference_area=mv.analytic_box_area(*dims),
    )
    assert abs(check.volume_error) <= mv.EXACT_TOLERANCE
    assert abs(check.area_error) <= mv.EXACT_TOLERANCE
    assert check.route_difference <= mv.ROUTE_TOLERANCE
    assert check.topologically_closed and check.watertight


def test_box_volume_is_invariant_under_rotation_and_translation():
    dims = (1.5, 0.75, 3.0)
    base = mv.check_closed_body(
        mv.box_mesh(*dims), name="base", family="box",
        reference_volume=mv.analytic_box_volume(*dims),
        reference_area=mv.analytic_box_area(*dims),
    )
    angles = np.radians([31.0, -17.0, 73.0])
    z, y, x = angles
    rx = np.array([[1, 0, 0], [0, math.cos(x), -math.sin(x)], [0, math.sin(x), math.cos(x)]])
    ry = np.array([[math.cos(y), 0, math.sin(y)], [0, 1, 0], [-math.sin(y), 0, math.cos(y)]])
    rz = np.array([[math.cos(z), -math.sin(z), 0], [math.sin(z), math.cos(z), 0], [0, 0, 1]])
    mesh = mv.box_mesh(*dims)
    from tensorcad.mesh_xval import TriMesh

    moved = TriMesh(
        mesh.points @ (rz @ ry @ rx).T + np.array([13.7, -42.1, 7.3]),
        mesh.faces,
        "moved",
    )
    check = mv.check_closed_body(
        moved, name="moved", family="box",
        reference_volume=mv.analytic_box_volume(*dims),
        reference_area=mv.analytic_box_area(*dims),
    )
    assert abs(check.volume_error) <= mv.EXACT_TOLERANCE
    assert abs(check.area_error) <= mv.EXACT_TOLERANCE
    assert check.measured_volume == pytest.approx(base.measured_volume, rel=1e-14)


# ---------------------------------------------------------------------------
# curved body: the deficit is the tessellation, at the mesh's own order
# ---------------------------------------------------------------------------
def test_sphere_ladder_converges_at_second_order_and_crosses_the_acceptance_line():
    rows = mv.convergence_study(1.0, range(0, 6))
    errors = [row.volume_error for row in rows]
    assert all(err < 0.0 for err in errors), "an inscribed polyhedron is below the sphere"
    pairs = list(itertools.pairwise(errors))
    assert all(abs(later) < abs(earlier) for earlier, later in pairs)
    assert mv.observed_order(rows) == pytest.approx(2.0, abs=0.15)
    assert mv.observed_order(rows, metric="area_error") == pytest.approx(2.0, abs=0.15)
    assert abs(rows[-1].volume_error) <= mv.RELATIVE_TOLERANCE
    assert abs(rows[2].volume_error) > mv.RELATIVE_TOLERANCE, "a coarse rung must miss the line"
    assert mv.first_level_within(rows) is not None
    for row in rows:
        assert row.topologically_closed and row.watertight
        assert row.route_difference <= mv.ROUTE_TOLERANCE


def test_sphere_error_is_dimensionless_and_healthily_below_the_line():
    coarse = mv.convergence_study(1.0, (4,))[0]
    scaled = mv.convergence_study(2.5, (4,))[0]
    assert scaled.volume_error == pytest.approx(coarse.volume_error, abs=1e-12)
    assert scaled.measured_volume == pytest.approx(coarse.measured_volume * 2.5**3, rel=1e-13)


# ---------------------------------------------------------------------------
# watertightness: the four defect families, and what volume cannot see
# ---------------------------------------------------------------------------
def test_every_injected_defect_is_named():
    controls = {row["name"]: row for row in mv.defect_detection(icosphere(1.0, 2))}
    assert set(controls) == {"hole", "flipped_facet", "duplicate_facet", "unwelded_soup"}
    for name, row in controls.items():
        assert row["detected"], (name, row)
        assert row["measured"]["watertight"] == 1
        assert row["defects_named"], name


def test_volume_alone_cannot_see_an_unwelded_soup():
    """Volume and watertightness are separate checks; this is why both are reported."""
    clean = icosphere(1.0, 2)
    controls = {row["name"]: row for row in mv.defect_detection(clean)}
    soup = controls["unwelded_soup"]
    assert soup["volume_error"] == pytest.approx(0.0, abs=1e-14)
    assert soup["measured"]["open_edges"] > 0
    holed = controls["hole"]
    assert holed["measured"]["open_edges"] > 0
    assert holed["volume_error"] < 0.0, "removing facets removes volume as well as closing nothing"


def test_audit_rejects_a_body_with_a_hole():
    mesh = icosphere(1.0, 2)
    from tensorcad.geometry import Mesh

    surface = Mesh(tuple(map(tuple, mesh.points.tolist())), tuple(map(tuple, mesh.faces.tolist())))
    assert audit_watertightness(surface, length=2.0).passed
    assert mv.check_closed_body(
        mesh, name="clean", family="sphere",
        reference_volume=mv.analytic_sphere_volume(1.0),
        reference_area=mv.analytic_sphere_area(1.0),
    ).watertight


# ---------------------------------------------------------------------------
# the benchmark case itself
# ---------------------------------------------------------------------------
def test_benchmark_case_passes_its_declared_gates(tmp_path):
    script = _REPO_ROOT / "benchmarks" / "verified" / "mesh_volume_watertight" / "run_benchmark.py"
    spec = importlib.util.spec_from_file_location("mesh_volume_watertight_bench", script)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    result = module.run(tmp_path)
    assert result["passed"], result["failures"]
    assert (tmp_path / "result.json").exists()
    assert (tmp_path / "convergence.csv").exists()
    exemptions = [row for row in result["gates"] if row["case"] == "sphere_r1_sub4"]
    assert any(row["quantity"] == "coarse_must_exceed_tolerance" for row in exemptions)
