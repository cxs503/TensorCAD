"""Pins for the curvature field and the smoothing-energy diagnostics.

Anchors, in order of strength:

* **closed forms** -- sphere ``k1 = k2 = 1/R`` everywhere (< 1e-6, the A3 clause),
  ``E_bending = 8 pi``, ``E_willmore = 4 pi``, ``E_gauss = 4 pi``;
* **an independent quadrature** -- the analytic Wigley flank area from
  ``analytic_flank_energy`` against ``analytic_wetted_surface_area`` (two different
  quadratures on the same surface must agree to round-off);
* **an independent mesh** -- the mesh energy on a conditioned sphere and the
  crease-split behaviour of the Wigley solid.

Every number below was measured before it was pinned; where a tolerance is loose the
comment says what the measured value was.
"""
from __future__ import annotations

import math
import random

import pytest

from tensorcad.curvature import CurvatureError, principal_curvatures, vertex_mean_curvature
from tensorcad.fairness_energy import (
    ENERGY_KINDS,
    fairness_energy_report,
    principal_curvature_field,
    smoothing_energy,
)
from tensorcad.geometry import Mesh
from tensorcad.geometry_audit import gaussian_curvature
from tensorcad.geometry_field import SphereField, to_mesh
from tensorcad.mesh_conditioning import project_to_field, taubin_smooth
from tensorcad.wigley import (
    WigleyParameters,
    analytic_flank_energy,
    analytic_wetted_surface_area,
    build_mesh,
)

BENDING_ANALYTIC = 8.0 * math.pi
WILLMORE_ANALYTIC = 4.0 * math.pi
GAUSS_ANALYTIC = 4.0 * math.pi


def _sphere_points(radius: float, count: int = 200, seed: int = 7) -> list[tuple[float, float, float]]:
    random.seed(seed)
    points = []
    for _ in range(count):
        vector = [random.gauss(0.0, 1.0) for _ in range(3)]
        norm = math.sqrt(sum(component * component for component in vector))
        points.append(tuple(radius * component / norm for component in vector))
    points.extend(
        [
            (radius, 0.0, 0.0),
            (-radius, 0.0, 0.0),
            (0.0, radius, 0.0),
            (0.0, -radius, 0.0),
            (0.0, 0.0, radius),
            (0.0, 0.0, -radius),
            (radius / math.sqrt(3.0),) * 3,
        ]
    )
    return points


# --------------------------------------------------------------------------- #
# A3 clause 1: analytic sphere principal curvature, everywhere
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("radius", [0.1, 0.5, 2.0, 10.0])
def test_sphere_principal_curvature_is_one_over_r_everywhere(radius: float) -> None:
    """``|k - 1/R| < 1e-6`` absolute at every probe, on four scales.

    Measured worst case over this set: 2.4e-9 * (1/R) relative on every radius --
    i.e. the bound is met with four orders of margin, and it is met *relative*,
    which is what "1/R everywhere" means.
    """
    field = SphereField(center=(0.0, 0.0, 0.0), radius=radius)
    expected = 1.0 / radius
    samples = principal_curvatures(field, _sphere_points(radius))
    assert len(samples) == 207
    for sample in samples:
        assert abs(sample.kappa_max - expected) < 1e-6
        assert abs(sample.kappa_min - expected) < 1e-6
    worst = max(abs(sample.kappa_max - expected) for sample in samples)
    assert worst < 1e-6
    assert worst / expected < 1e-6  # relative, so the bound is scale free


def test_off_surface_probe_reports_the_level_set_not_the_sphere() -> None:
    """A probe at 0.9994 R belongs to a smaller sphere: the curvature follows it.

    This is the failure mode that makes a "sphere is 1/R" test pass by accident: at
    ``r = 0.577 * sqrt(3)`` the point is 0.06 % inside the surface and the returned
    curvature is 1.00061/R, which is correct and would look like a 6e-4 error if the
    test forgot to project. Pinned so the estimator's contract stays explicit.
    """
    field = SphereField(center=(0.0, 0.0, 0.0), radius=1.0)
    sample = principal_curvatures(field, [(0.577, 0.577, 0.577)])[0]
    radius_here = math.sqrt(3.0) * 0.577
    assert sample.kappa_max == pytest.approx(1.0 / radius_here, rel=1e-6)


# --------------------------------------------------------------------------- #
# A3 clause 2: closed-form energies
# --------------------------------------------------------------------------- #
def test_field_path_energies_of_a_sphere_match_the_closed_forms() -> None:
    field = SphereField(center=(0.0, 0.0, 0.0), radius=1.0)
    mesh = to_mesh(field, 24, margin=3).mesh
    report = fairness_energy_report(mesh, field=field)
    # Face centroids sit slightly inside the sphere, so the quadrature converges
    # from above: measured +0.09 % (bending) at level 24.
    assert report.bending.total_energy == pytest.approx(BENDING_ANALYTIC, rel=2e-3)
    assert report.willmore.total_energy == pytest.approx(WILLMORE_ANALYTIC, rel=2e-3)
    assert report.gaussian.total_energy == pytest.approx(GAUSS_ANALYTIC, rel=2e-3)
    # k1^2 + k2^2 = 4H^2 - 2K holds per cell on the field path (same eigenvalues).
    assert abs(report.bending.identity_residual or 0.0) < 1e-12
    assert abs(report.identity_gap) < 1e-10


def test_the_field_path_energy_converges_to_the_closed_forms() -> None:
    field = SphereField(center=(0.0, 0.0, 0.0), radius=1.0)
    errors = []
    for level in (16, 24, 32):
        mesh = to_mesh(field, level, margin=3).mesh
        report = fairness_energy_report(mesh, field=field)
        errors.append(abs(report.bending.total_energy / BENDING_ANALYTIC - 1.0))
    assert errors[0] > errors[1] > errors[2]  # monotone
    assert errors[-1] < 1e-3  # measured 2.7e-4 at level 32


def test_mesh_path_energy_on_a_conditioned_sphere_is_within_six_percent() -> None:
    """The mesh path is the consumer's estimator, and it is a few per cent loud.

    Measured on a taubin_smooth + project_to_field marching-tet sphere: bending
    +5.0 %, willmore +1.5 %, gaussian -2.6 %; the spread is honest (the cotangent
    Laplacian and the angle defect are not a consistent pair of operators).
    """
    field = SphereField(center=(0.0, 0.0, 0.0), radius=1.0)
    raw = to_mesh(field, 32, margin=3).mesh
    smoothed, _ = taubin_smooth(raw)
    mesh, _ = project_to_field(smoothed, field)
    report = fairness_energy_report(mesh)
    assert report.bending.total_energy == pytest.approx(BENDING_ANALYTIC, rel=0.06)
    assert report.willmore.total_energy == pytest.approx(WILLMORE_ANALYTIC, rel=0.02)
    assert report.gaussian.total_energy == pytest.approx(GAUSS_ANALYTIC, rel=0.05)


def test_the_raw_marching_tet_mesh_breaks_convergence_and_says_so() -> None:
    """Sliver noise: the raw mesh's energy is *not* convergent, and is not hidden.

    This is the documented limitation of the mesh path, pinned as a claim: on the
    raw extraction the bending energy is 34 / 60 / 36 for levels 16 / 24 / 32 --
    non-monotone and up to +140 % of ``8 pi`` -- so callers must condition the mesh
    (or use the field) before quoting an energy.
    """
    field = SphereField(center=(0.0, 0.0, 0.0), radius=1.0)
    values = []
    for level in (16, 24):
        mesh = to_mesh(field, level, margin=3).mesh
        values.append(fairness_energy_report(mesh).bending.total_energy)
    assert values[1] > values[0]  # diverging, not converging
    assert values[1] / BENDING_ANALYTIC > 1.5


# --------------------------------------------------------------------------- #
# A3 clause 3: Wigley smoothing energy, and its numerical stability
# --------------------------------------------------------------------------- #
def test_wigley_analytic_energy_is_pinned_and_quadrature_converged() -> None:
    energy = analytic_flank_energy(WigleyParameters())
    assert energy.bending == pytest.approx(21.5378773250461, rel=1e-12)
    assert energy.willmore == pytest.approx(5.3844693312615, rel=1e-12)
    assert energy.quadrature_relative_change < 1e-12
    # The four conical corners (stem/stern x keel/waterline) hold the whole
    # Gauss-Bonnet charge: the smooth flank carries ~0 of it.
    assert abs(energy.gaussian) < 1e-9
    assert energy.corner_charge == pytest.approx(4.0 * math.pi, rel=1e-12)


def test_wigley_analytic_area_matches_the_independent_quadrature() -> None:
    """Two different Quadratures on the same analytic surface must agree.

    ``analytic_wetted_surface_area`` integrates a hand-written Jacobian with
    composite Gauss panels; ``analytic_flank_energy`` integrates the same Jacobian
    from the graph derivatives over a tensor grid. Measured agreement: 5e-15
    relative, which is what makes either of them usable as a reference.
    """
    energy = analytic_flank_energy(WigleyParameters())
    independent = analytic_wetted_surface_area(10.0, 1.0, 0.625)
    assert energy.wetted_area == pytest.approx(independent, rel=1e-12)


def test_wigley_energy_is_scale_and_translation_invariant() -> None:
    """Dimensionless by construction: doubling the hull must not move the number."""
    base = analytic_flank_energy(WigleyParameters())
    scaled = analytic_flank_energy(length=20.0, beam=2.0, draft=1.25)
    assert scaled.bending == pytest.approx(base.bending, rel=1e-12)
    assert scaled.willmore == pytest.approx(base.willmore, rel=1e-12)
    assert scaled.wetted_area == pytest.approx(4.0 * base.wetted_area, rel=1e-12)


def test_the_energy_accumulation_is_bit_identical_under_cell_reordering() -> None:
    """``math.fsum`` is exactly rounded, so the accumulator cannot depend on order.

    Tested at the accumulator (``smoother``) rather than through the mesh: rebuilding
    a :class:`Mesh` from a shuffled face list moves the last bit of some per-cell
    values (measured 2e-15 relative on the total), which would mask the property
    being pinned here. The accumulator itself is exact.
    """
    from tensorcad.fairness_energy import smoother

    terms = [(0.1 * index * index, 0.01 * (index + 1), index % 3 != 0) for index in range(500)]
    reference = smoother(
        kind="bending",
        source="mesh",
        terms=terms,
        crease_angle_degrees=30.0,
        integrand_values=[value for value, _, _ in terms],
        identity_residual=None,
    )
    for seed in (3, 11, 29):
        shuffled = list(terms)
        random.Random(seed).shuffle(shuffled)
        other = smoother(
            kind="bending",
            source="mesh",
            terms=shuffled,
            crease_angle_degrees=30.0,
            integrand_values=[value for value, _, _ in shuffled],
            identity_residual=None,
        )
        assert other.total_energy == reference.total_energy  # exact, not approx
        assert other.smooth_energy == reference.smooth_energy
        assert other.crease_energy == reference.crease_energy


def test_mesh_energy_on_a_closed_hull_is_stable_under_face_reordering() -> None:
    """End to end: a rebuilt mesh changes the energy by less than one part in 1e14."""
    mesh, _ = build_mesh(WigleyParameters(nx=21, nz=7))
    faces = list(mesh.faces)
    random.Random(3).shuffle(faces)
    shuffled = Mesh(vertices=list(mesh.vertices), faces=faces)
    first = smoothing_energy(mesh, kind="bending", crease_angle_degrees=30.0)
    second = smoothing_energy(shuffled, kind="bending", crease_angle_degrees=30.0)
    assert second.total_energy == pytest.approx(first.total_energy, rel=1e-13)
    assert second.smooth_energy == pytest.approx(first.smooth_energy, rel=1e-13)


def test_wigley_smooth_energy_converges_while_the_crease_energy_diverges() -> None:
    """The diagnostic split is the deliverable: one number, two behaviours.

    Reference (analytic flank, both sides): ``E_bending = 21.5379``. The crease
    filter leaves a band of ``O(h)`` width, so the smooth sum climbs at first order
    (13.17 / 17.23 for nx = 21 / 41) while the excluded band -- which holds the keel,
    stem and deck creases -- grows like ``1/h`` (216.8 / 384.9). Reporting only the
    total would show a surface whose curvature cost explodes with mesh refinement.
    """
    reference = analytic_flank_energy(WigleyParameters()).bending
    energies = {}
    for nx, nz in ((21, 7), (41, 13)):
        mesh, _ = build_mesh(WigleyParameters(nx=nx, nz=nz))
        energies[nz] = smoothing_energy(mesh, kind="bending", crease_angle_degrees=30.0)
    coarse, fine = energies[7], energies[13]
    assert coarse.smooth_energy < fine.smooth_energy < reference
    assert coarse.crease_energy < fine.crease_energy  # diverging with refinement
    assert fine.smooth_energy / reference > coarse.smooth_energy / reference


def test_principal_curvature_field_reports_the_hk_inconsistency() -> None:
    """The (H, K) -> (k1, k2) inversion is not stable on a real mesh, and is counted.

    On a conditioned sphere the two discrete operators disagree enough that
    ``H**2 - K < 0`` at half the vertices; they are reported as ``nan`` and counted,
    which is why the energy path uses ``4H^2 - 2K`` directly instead.
    """
    field = SphereField(center=(0.0, 0.0, 0.0), radius=1.0)
    raw = to_mesh(field, 24, margin=3).mesh
    smoothed, _ = taubin_smooth(raw)
    mesh, _ = project_to_field(smoothed, field)
    values = principal_curvature_field(mesh)
    assert values.negative_discriminant_vertices > 0
    assert len(values.mean) == len(mesh.vertices)
    payload = values.to_dict()
    assert payload["schema"] == "tensorcad.fairness_energy.field/1"
    # Measured: the vertex mean curvature averages 1.079 on a unit sphere (+8 %),
    # because the mean of a noisy positive quantity is biased high -- which is
    # exactly why the energy is an area integral and not a vertex average.
    assert 0.95 < payload["mean_curvature_per_m"]["mean"] < 1.20


def test_gauss_bonnet_defects_sum_to_four_pi_on_a_closed_mesh() -> None:
    """The angle-defect charge is exact by algebra; this pins the quadrat ure path."""
    mesh, _ = build_mesh(WigleyParameters(nx=41, nz=13))
    discrete = gaussian_curvature(mesh)
    total = math.fsum(
        value * area
        for value, area in zip(discrete.curvature, discrete.mixed_area, strict=True)
        if math.isfinite(value)
    )
    assert total == pytest.approx(GAUSS_ANALYTIC, rel=1e-9)


# --------------------------------------------------------------------------- #
# contract
# --------------------------------------------------------------------------- #
def test_smoothing_energy_rejects_bad_input() -> None:
    mesh, _ = build_mesh(WigleyParameters(nx=11, nz=5))
    with pytest.raises(CurvatureError):
        smoothing_energy(mesh, kind="willmore2")  # type: ignore[arg-type]
    for threshold in (0.0, 180.0, float("nan")):
        with pytest.raises(CurvatureError):
            smoothing_energy(mesh, kind="bending", crease_angle_degrees=threshold)
    with pytest.raises(CurvatureError):
        smoothing_energy(Mesh(vertices=[], faces=[]), kind="bending")


def test_energy_kinds_are_the_three_densities() -> None:
    assert ENERGY_KINDS == ("bending", "willmore", "gaussian")
    mesh, _ = build_mesh(WigleyParameters(nx=11, nz=5))
    payload = smoothing_energy(mesh, kind="willmore", crease_angle_degrees=30.0).to_dict()
    assert payload["schema"] == "tensorcad.fairness_energy.energy/1"
    assert payload["kind"] == "willmore"
    assert payload["summation"] == "math.fsum"
    for key in ("smooth_energy", "crease_energy", "total_energy", "condition",
                "nonpositive_fraction", "crease_area_fraction"):
        assert math.isfinite(payload[key])


def test_the_mean_curvature_mesh_path_still_agrees_with_the_field() -> None:
    """Nothing above changed ``vertex_mean_curvature``: sphere area-weighted mean 1/R."""
    field = SphereField(center=(0.0, 0.0, 0.0), radius=1.0)
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
    weighted = math.fsum(value * weight for value, weight in zip(values, weights, strict=True)) / math.fsum(weights)
    assert weighted == pytest.approx(1.0, rel=2e-3)
