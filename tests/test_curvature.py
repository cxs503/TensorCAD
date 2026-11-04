"""Curvature and fairness: the estimator against closed forms, and the audit's claims."""
from __future__ import annotations

import itertools
import math
import subprocess
import sys

import pytest

from tensorcad import mesh as meshlib
from tensorcad.curvature import (
    DEFAULT_CREASE_DEGREES,
    CurvatureError,
    CurvatureSample,
    assess_fairness,
    principal_curvatures,
    vertex_mean_curvature,
)
from tensorcad.geometry import Mesh
from tensorcad.geometry_field import (
    CylinderField,
    EllipsoidField,
    GeometryField,
    PlaneField,
    SphereField,
    to_mesh,
)
from tensorcad.validation import inspect_geometry


def _uv_sphere(radius: float, nu: int, nv: int) -> Mesh:
    """A legal UV sphere: one vertex per pole and a triangle fan onto it.

    ``nu`` is the number of latitude *bands* and ``nv`` the longitude count. The
    poles are single vertices, so no ring of ``nv`` coincident points collapses
    onto itself and no face has zero area -- which matters, because the fairness
    report refuses a degenerate face and a collapsed vertex is *not* a point of
    zero curvature (see ``test_a_collapsed_pole_is_refused_by_the_fairness_report``).
    """
    vertices = [(0.0, 0.0, radius)]
    for i in range(1, nu):
        theta = math.pi * i / nu
        for j in range(nv):
            phi = 2.0 * math.pi * j / nv
            vertices.append(
                (
                    radius * math.sin(theta) * math.cos(phi),
                    radius * math.sin(theta) * math.sin(phi),
                    radius * math.cos(theta),
                )
            )
    south = len(vertices)
    vertices.append((0.0, 0.0, -radius))

    def ring(i: int, j: int) -> int:
        return 1 + (i - 1) * nv + (j % nv)

    faces = []
    for j in range(nv):
        faces.append((0, ring(1, j), ring(1, j + 1)))
        faces.append((south, ring(nu - 1, j + 1), ring(nu - 1, j)))
    for i in range(1, nu - 1):
        for j in range(nv):
            faces.append((ring(i, j), ring(i + 1, j), ring(i, j + 1)))
            faces.append((ring(i, j + 1), ring(i + 1, j), ring(i + 1, j + 1)))
    return Mesh(vertices=vertices, faces=faces)


def _pole_collapsed_sphere(radius: float, nu: int, nv: int) -> Mesh:
    """The historical fixture: ``nu + 1`` rings, so each pole ring is ``nv`` copies.

    Kept only as the negative control for the degeneracy contract; it is not a
    sphere mesh in the sense this package's validator accepts.
    """
    vertices = []
    grid = {}
    for i in range(nu + 1):
        theta = math.pi * i / nu
        for j in range(nv):
            phi = 2.0 * math.pi * j / nv
            vertices.append(
                (
                    radius * math.sin(theta) * math.cos(phi),
                    radius * math.sin(theta) * math.sin(phi),
                    radius * math.cos(theta),
                )
            )
            grid[(i, j)] = len(vertices) - 1
    faces = []
    for i in range(nu):
        for j in range(nv):
            j2 = (j + 1) % nv
            a, b = grid[(i, j)], grid[(i, j2)]
            c, d = grid[(i + 1, j)], grid[(i + 1, j2)]
            faces.append((a, c, b))
            faces.append((b, c, d))
    return Mesh(vertices=vertices, faces=faces)


def _cube(size: float = 1.0) -> Mesh:
    h = 0.5 * size
    vertices = [
        (-h, -h, -h), (h, -h, -h), (h, h, -h), (-h, h, -h),
        (-h, -h, h), (h, -h, h), (h, h, h), (-h, h, h),
    ]
    faces = [
        (0, 2, 1), (0, 3, 2), (4, 5, 6), (4, 6, 7),
        (0, 1, 5), (0, 5, 4), (1, 2, 6), (1, 6, 5),
        (2, 3, 7), (2, 7, 6), (3, 0, 4), (3, 4, 7),
    ]
    return Mesh(vertices=vertices, faces=faces)


class _VanishingGradient(GeometryField):
    """A field whose gradient is zero: the shape operator is undefined on it."""

    kind = "vanishing"

    def signed_distance(self, points):
        return [0.0 for _ in points]

    def gradient(self, points):
        return [(0.0, 0.0, 0.0) for _ in points]


def test_plane_curvature_is_exactly_zero() -> None:
    plane = PlaneField(normal=(0.0, 0.0, 1.0), offset=0.0)
    for sample in principal_curvatures(plane, [(1.0, 2.0, 0.0), (-3.0, 0.5, 0.0)]):
        assert sample.kappa_max == 0.0
        assert sample.kappa_min == 0.0
        assert sample.step_error == 0.0
        assert sample.gradient_norm == pytest.approx(1.0, abs=1e-15)


def test_sphere_curvature_is_one_over_r() -> None:
    radius = 0.5
    sphere = SphereField(center=(0.0, 0.0, 0.0), radius=radius)
    points = [(radius, 0.0, 0.0), (0.0, radius, 0.0), (0.0, 0.0, -radius)]
    for sample in principal_curvatures(sphere, points):
        assert sample.kappa_max == pytest.approx(1.0 / radius, rel=1e-12)
        assert sample.kappa_min == pytest.approx(1.0 / radius, rel=1e-12)
        # the normal points outward, so a convex surface reports positive curvature
        assert sample.kappa_max > 0.0


def test_cylinder_curvature_pair_is_hoop_and_zero() -> None:
    radius = 0.35
    cylinder = CylinderField(radius=radius, height=1.0)
    for sample in principal_curvatures(cylinder, [(radius, 0.0, 0.0), (0.0, radius, 0.2)]):
        assert sample.kappa_max == pytest.approx(1.0 / radius, rel=1e-12)
        assert abs(sample.kappa_min) < 1e-9
        assert sample.mean_curvature == pytest.approx(0.5 / radius, rel=1e-12)
        assert abs(sample.gaussian_curvature) < 1e-9


def test_ellipsoid_poles_match_the_closed_forms() -> None:
    a, b, c = 0.6, 0.4, 0.25
    ellipsoid = EllipsoidField(semi_axes=(a, b, c))
    cases = {
        (a, 0.0, 0.0): (a / b**2, a / c**2),
        (0.0, b, 0.0): (b / c**2, b / a**2),
        (0.0, 0.0, c): (c / b**2, c / a**2),
    }
    for point, expected in cases.items():
        sample = principal_curvatures(ellipsoid, [point])[0]
        measured = sorted((sample.kappa_max, sample.kappa_min))
        assert measured == pytest.approx(sorted(expected), rel=1e-10)


def test_an_implicit_field_is_not_a_signed_distance_and_says_so() -> None:
    """``|grad f| != 1`` on the ellipsoid: the invariance of the estimate is the point."""
    ellipsoid = EllipsoidField(semi_axes=(0.6, 0.4, 0.25))
    sample = principal_curvatures(ellipsoid, [(0.6, 0.0, 0.0)])[0]
    assert sample.gradient_norm != pytest.approx(1.0)  # 1/a = 1.67 on this field
    sphere = SphereField(center=(0.0, 0.0, 0.0), radius=0.5)
    assert principal_curvatures(sphere, [(0.5, 0.0, 0.0)])[0].gradient_norm == pytest.approx(1.0, rel=1e-6)


def test_sample_is_a_frozen_record_with_the_derived_invariants() -> None:
    sphere = SphereField(center=(0.0, 0.0, 0.0), radius=1.0)
    sample = principal_curvatures(sphere, [(1.0, 0.0, 0.0)])[0]
    assert isinstance(sample, CurvatureSample)
    assert sample.mean_curvature == pytest.approx(1.0, rel=1e-12)
    assert sample.gaussian_curvature == pytest.approx(1.0, rel=1e-12)
    payload = sample.to_dict()
    assert payload["schema"] == "tensorcad.curvature.sample/1"
    assert payload["mean_curvature_per_m"] == pytest.approx(1.0, rel=1e-12)
    with pytest.raises(Exception):  # noqa: B017  # any mutation error is the contract
        sample.kappa_max = 2.0  # frozen


def test_refinement_reports_a_step_error_and_can_be_switched_off() -> None:
    sphere = SphereField(center=(0.0, 0.0, 0.0), radius=1.0)
    refined = principal_curvatures(sphere, [(1.0, 0.0, 0.0)])[0]
    assert math.isfinite(refined.step_error)
    assert refined.step_error < 1e-5
    single = principal_curvatures(sphere, [(1.0, 0.0, 0.0)], refine=False)[0]
    assert math.isnan(single.step_error)
    # The unrefined probe is a *second-order* estimator, so it cannot be exact to
    # 1e-12 -- demanding that was the defect, not the estimator. On the sphere field
    # (an exact signed distance) the truncation error is measured at exactly
    # ``0.5 * step**2``, which is what this pins instead of a tolerance nobody can
    # meet. ``test_the_unrefined_probe_obeys_its_measured_truncation_law`` measures
    # the same law over a step ladder and the observed order.
    assert abs(single.kappa_max - 1.0) == pytest.approx(0.5 * single.step**2, rel=1e-3)
    assert refined.kappa_max == pytest.approx(1.0, rel=1e-12)


def test_bad_probe_points_are_rejected() -> None:
    sphere = SphereField(center=(0.0, 0.0, 0.0), radius=1.0)
    with pytest.raises(CurvatureError):
        principal_curvatures(sphere, [(1.0, 0.0)])
    with pytest.raises(CurvatureError):
        principal_curvatures(sphere, [(float("nan"), 0.0, 0.0)])
    for step in (0.0, -1e-3, float("inf")):
        with pytest.raises(CurvatureError):
            principal_curvatures(sphere, [(1.0, 0.0, 0.0)], step=step)


def test_a_vanishing_gradient_is_an_error_not_a_number() -> None:
    with pytest.raises(CurvatureError):
        principal_curvatures(_VanishingGradient(), [(0.0, 0.0, 0.0)])


def test_mesh_estimate_converges_on_a_regular_sphere() -> None:
    radius = 1.0
    errors = []
    for nu in (16, 32):
        values = [abs(value) for value in vertex_mean_curvature(_uv_sphere(radius, nu, 2 * nu))]
        mean = sum(values) / len(values)
        errors.append(abs(mean - 1.0 / radius) / (1.0 / radius))
        assert any(value != 0.0 for value in values)
    assert errors[0] < 0.15
    assert errors[1] < errors[0]


def test_the_unrefined_probe_obeys_its_measured_truncation_law() -> None:
    """Second order in the step, with the measured constant, on a closed form.

    ``kappa(h) - 1/R = 0.5 * h**2`` on the sphere field to better than 1e-3 of the
    law itself, and the observed order over a halving ladder is 2.00. This is the
    measurement the ``refine=False`` path is allowed to be judged by; asserting
    exactness there is a category error, since the difference quotient deletes the
    O(h^2) term only when ``refine=True`` extrapolates it away.
    """
    radius = 1.0
    sphere = SphereField(center=(0.0, 0.0, 0.0), radius=radius)
    probe = [(radius, 0.0, 0.0)]
    errors = []
    for step in (8e-4, 4e-4, 2e-4, 1e-4):
        sample = principal_curvatures(sphere, probe, step=step, refine=False)[0]
        error = abs(sample.kappa_max - 1.0 / radius)
        errors.append(error)
        assert error == pytest.approx(0.5 * step**2, rel=1e-3)
    for coarse, fine in itertools.pairwise(errors):
        order = math.log(coarse / fine) / math.log(2.0)
        assert order == pytest.approx(2.0, abs=1e-3)


def test_richardson_extrapolation_removes_the_second_order_term() -> None:
    """Same ladder, ``refine=True``: the analytic answer to round-off."""
    sphere = SphereField(center=(0.0, 0.0, 0.0), radius=1.0)
    probe = [(1.0, 0.0, 0.0)]
    for step in (8e-4, 4e-4, 2e-4, 1e-4, 5e-5):
        refined = principal_curvatures(sphere, probe, step=step, refine=True)[0]
        assert refined.kappa_max == pytest.approx(1.0, rel=1e-12)
        assert refined.step_error == pytest.approx(1.5 * refined.step**2, rel=1e-3)


def test_the_cotangent_mean_curvature_is_second_order() -> None:
    """The mesh estimator's order, measured on a legal sphere instead of assumed."""
    radius = 1.0
    errors = []
    for level in (16, 32, 64):
        values = [abs(value) for value in vertex_mean_curvature(_uv_sphere(radius, level, 2 * level))]
        errors.append(abs(sum(values) / len(values) - 1.0 / radius))
    order = math.log(errors[-2] / errors[-1]) / math.log(2.0)
    # measured 1.84 (16->32) and 1.93 (32->64): the fan's pole triangles are thin,
    # so the asymptotic rate is approached from below, not exactly hit.
    assert order == pytest.approx(2.0, abs=0.1), order
    assert errors[-1] < 1.0e-3  # 1.8e-4 at level 64


def test_mesh_estimate_conserves_the_integral_on_a_raw_marching_tet_mesh() -> None:
    """Slivers destroy the pointwise value and leave the area-weighted mean alone."""
    radius = 1.0
    field = SphereField(center=(0.0, 0.0, 0.0), radius=radius)
    mesh = to_mesh(field, 16, margin=3).mesh
    values = vertex_mean_curvature(mesh)
    weights = [0.0] * len(mesh.vertices)
    for face in mesh.faces:
        a, b, c = (mesh.vertices[index] for index in face)
        u = [b[axis] - a[axis] for axis in range(3)]
        v = [c[axis] - a[axis] for axis in range(3)]
        cross = (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])
        area = 0.5 * math.sqrt(sum(component * component for component in cross)) / 3.0
        for index in face:
            weights[index] += area
    weighted = sum(value * weight for value, weight in zip(values, weights, strict=False)) / sum(weights)
    assert weighted == pytest.approx(1.0 / radius, rel=2e-3)


def test_mesh_estimate_rejects_an_empty_mesh() -> None:
    with pytest.raises(CurvatureError):
        vertex_mean_curvature(Mesh(vertices=[], faces=[]))


def test_a_sphere_has_no_creases_and_is_watertight() -> None:
    report = assess_fairness(_uv_sphere(1.0, 16, 32))
    assert report.watertight is True
    assert report.edge_manifold is True
    assert report.crease_edges == 0
    # A smooth sphere has no crease edge at any resolution; the largest dihedral in
    # the mesh is the discretisation of the fan, measured at 11.29 deg on 16x32 and
    # 5.63 deg on 32x64. The threshold is therefore the measured value at this
    # resolution, not an aspiration: "no crease" is the ``crease_edges == 0`` claim
    # above, and the angle only has to stay out of the 15 deg crease band.
    assert report.dihedral_degrees["p95"] < 12.0
    assert report.dihedral_degrees["max"] < DEFAULT_CREASE_DEGREES
    assert report.curvature_jump_per_m is None  # no field supplied
    assert report.field_curvature_per_m is None
    assert report.mesh_curvature_per_m["area_weighted_mean"] == pytest.approx(1.0, rel=0.02)


def test_the_dihedral_p95_halves_when_the_sphere_is_refined() -> None:
    """The 16x32 threshold above is a resolution artefact, and this says so."""
    coarse = assess_fairness(_uv_sphere(1.0, 16, 32))
    fine = assess_fairness(_uv_sphere(1.0, 32, 64))
    ratio = coarse.dihedral_degrees["p95"] / fine.dihedral_degrees["p95"]
    assert ratio == pytest.approx(2.0, rel=0.05), ratio  # first order in h, measured 2.006


def test_a_collapsed_pole_is_refused_by_the_fairness_report() -> None:
    """The historical fixture's contract, pinned: a collapsed vertex is not flat.

    A UV sphere built from ``nu + 1`` rings carries ``nv`` coincident vertices at
    each pole, hence ``2 * nv`` zero-area faces. The cotangent Laplacian returns
    *zero* curvature at the collapsed vertices -- a silent "this is flat" -- while
    ``assess_fairness`` refuses the mesh outright. Both halves are pinned here so
    that neither can drift into looking like a sphere.
    """
    mesh = _pole_collapsed_sphere(1.0, 16, 32)
    assert len(meshlib.zero_area_faces(mesh)) == 32
    failed = {check.name for check in inspect_geometry(mesh) if not check.passed}
    assert {"no_degenerate_faces", "watertight", "vertex_manifold"} <= failed
    values = vertex_mean_curvature(mesh)
    assert all(math.isfinite(value) for value in values)
    assert min(values) == pytest.approx(0.0, abs=1e-9)  # collapsed vertex: reported flat
    with pytest.raises(CurvatureError):
        assess_fairness(mesh)


def test_a_cube_has_exactly_twelve_creases() -> None:
    """Twelve 90-degree box edges; the six face diagonals are coplanar and are not creases."""
    report = assess_fairness(_cube(2.0))
    assert report.interior_edges == 18
    assert report.crease_edges == 12
    assert report.dihedral_degrees["max"] == pytest.approx(90.0, abs=1e-6)
    assert report.dihedral_degrees["min"] == pytest.approx(0.0, abs=1e-9)
    expected_crease_fraction = 12.0 * 2.0 / (12.0 * 2.0 + 6.0 * 2.0 * math.sqrt(2.0))
    assert report.crease_length_fraction == pytest.approx(expected_crease_fraction, rel=1e-9)


def test_a_field_turns_the_report_into_a_curvature_jump_measurement() -> None:
    sphere = SphereField(center=(0.0, 0.0, 0.0), radius=1.0)
    mesh = to_mesh(sphere, 16, margin=3).mesh
    report = assess_fairness(mesh, field=sphere)
    assert report.field_curvature_per_m is not None
    assert report.field_curvature_per_m["p50"] == pytest.approx(1.0, rel=0.05)
    # a sphere's true curvature is constant, so the jump across any edge is the
    # discretisation of the probe points only
    assert report.curvature_jump_per_m["p95"] < 0.05


def test_assess_fairness_rejects_an_impossible_crease_threshold() -> None:
    for threshold in (0.0, 180.0, -5.0, float("nan")):
        with pytest.raises(CurvatureError):
            assess_fairness(_cube(), crease_angle_degrees=threshold)


def test_the_report_serialises_with_its_schema() -> None:
    payload = assess_fairness(_cube()).to_dict()
    assert payload["schema"] == "tensorcad.curvature.fairness/1"
    assert payload["crease_angle_degrees"] == 15.0
    assert set(payload["mesh_curvature_per_m"]) >= {"count", "mean", "p50", "p95", "max",
                                                    "area_weighted_mean"}


def test_the_curvature_module_needs_no_third_party_stack() -> None:
    script = (
        "import sys;"
        "sys.modules['numpy'] = None; sys.modules['torch'] = None;"
        "import tensorcad.curvature as curvature;"
        "assert curvature.__all__"
    )
    completed = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
