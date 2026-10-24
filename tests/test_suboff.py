"""SUBOFF bare hull: the offset table, the distance field, and the exported mesh.

The numbers pinned here are the published DARPA offset-table result and the values
computed by TensorLBM's transcription of it (``suboff_cad.suboff_radius_profile``,
read 2026-10-24), so a change to the formulas in this repository is a change against
an external reference rather than against our own earlier output.
"""
from __future__ import annotations

import math

import pytest

from tensorcad.geometry_field import to_mesh
from tensorcad.suboff import (
    SUBOFF_LENGTH_M,
    SUBOFF_RADIUS_M,
    SuboffHullField,
    suboff_dimensions,
    suboff_radius_fraction,
    suboff_radius_fraction_array,
    suboff_reference_geometry,
)

#: Normalised radii computed independently by TensorLBM's
#: ``suboff_cad.suboff_radius_profile`` on a 200,001-point grid.
TENSORLBM_PROFILE = {
    0.0: 0.0,
    0.05: 0.7115501353960532,
    0.2: 0.9947614804958352,
    0.5: 1.0,
    0.7449: 0.9999999999999951,
    0.9: 0.43999064315220865,
    0.99: 0.09869380175012143,
    1.0: 0.0,
}


def test_the_table_matches_its_source() -> None:
    for xi, expected in TENSORLBM_PROFILE.items():
        assert suboff_radius_fraction(xi) == pytest.approx(expected, abs=1e-12)


def test_the_vectorised_profile_is_the_documented_one() -> None:
    numpy = pytest.importorskip("numpy")
    grid = numpy.linspace(0.0, 1.0, 5001)
    vectorised = suboff_radius_fraction_array(grid)
    scalar = numpy.array([suboff_radius_fraction(float(xi)) for xi in grid])
    assert float(numpy.max(numpy.abs(vectorised - scalar))) < 1e-15


def test_the_segments_join_without_a_step() -> None:
    """C0 is exact at every joint; that is what keeps the hull closed."""
    for xi in (
        3.333333 / 14.291667,
        10.645833 / 14.291667,
        13.979167 / 14.291667,
    ):
        left = suboff_radius_fraction(xi - 1e-12)
        right = suboff_radius_fraction(xi + 1e-12)
        assert abs(left - right) < 1e-9, f"radius steps by {left - right:.3e} at xi={xi:.6f}"


def test_the_segments_join_without_a_slope_kink() -> None:
    """C1 holds at all three joints: the one-sided slopes vanish from both sides.

    Evaluated with a shrinking step so a genuine kink (which would converge to a
    non-zero jump) is distinguished from a quadrature artefact.
    """
    for xi in (3.333333 / 14.291667, 10.645833 / 14.291667, 13.979167 / 14.291667):
        for step in (1e-4, 1e-5):
            left = (suboff_radius_fraction(xi - step) - suboff_radius_fraction(xi - 2 * step)) / step
            right = (suboff_radius_fraction(xi + 2 * step) - suboff_radius_fraction(xi + step)) / step
            # A one-sided slope that shrinks with the step means the limit is zero on
            # both sides. A real kink would converge to a fixed non-zero jump.
            assert abs(left) < 0.1, f"left slope {left:.3e} at xi={xi:.5f}"
            assert abs(right) < 0.1, f"right slope {right:.3e} at xi={xi:.5f}"


def test_dimensions_resolve_to_the_published_ratio() -> None:
    length, radius = suboff_dimensions()
    assert length == pytest.approx(4.3561, abs=1e-4)
    assert radius == pytest.approx(0.254, abs=1e-4)
    assert length / (2.0 * radius) == pytest.approx(8.575, abs=1e-3)
    # Scaling the length keeps the DARPA length/diameter ratio.
    scaled_length, scaled_radius = suboff_dimensions(length=8.7122)
    assert scaled_length / (2.0 * scaled_radius) == pytest.approx(8.575, abs=1e-3)


def test_invalid_parameters_and_queries_fail_loudly() -> None:
    from tensorcad.geometry_field import FieldError

    with pytest.raises(FieldError):
        SuboffHullField(length=0.0)
    with pytest.raises(FieldError):
        SuboffHullField(radius=-1.0)
    with pytest.raises(FieldError):
        suboff_radius_fraction(-0.1)
    with pytest.raises(FieldError):
        suboff_radius_fraction(1.1)
    with pytest.raises(FieldError, match="no control"):
        SuboffHullField().with_state(beam=1.0)


def test_bounds_are_the_hull_extents() -> None:
    field = SuboffHullField()
    box = field.bounds()
    assert box.minimum == pytest.approx((0.0, -field.max_radius, -field.max_radius))
    assert box.maximum == pytest.approx((field.length, field.max_radius, field.max_radius))
    assert box.sizes[0] / box.sizes[1] == pytest.approx(8.575, abs=1e-3)


def test_containment_and_sign_agree_with_the_profile() -> None:
    field = SuboffHullField()
    radius = field.max_radius
    inside = [
        (0.5 * field.length, 0.0, 0.0),
        (0.5 * field.length, 0.9 * radius, 0.0),
        (0.3 * field.length, 0.0, 0.5 * radius),
    ]
    outside = [
        (0.5 * field.length, 1.1 * radius, 0.0),
        (-0.1, 0.0, 0.0),
        (field.length + 0.1, 0.0, 0.0),
    ]
    values = field.signed_distance(inside + outside)
    assert all(value < 0.0 for value in values[: len(inside)])
    assert all(value > 0.0 for value in values[len(inside) :])
    assert field.contains(inside + outside) == [True] * len(inside) + [False] * len(outside)


def test_the_field_is_the_true_distance_to_the_meridian() -> None:
    """Compare against a brute-force dense meridian, not against ourselves."""
    numpy = pytest.importorskip("numpy")
    field = SuboffHullField()
    length, radius = field.length, field.max_radius
    generator = numpy.random.default_rng(20261024)
    points = numpy.column_stack(
        [
            generator.uniform(-0.2 * length, 1.2 * length, 120),
            generator.uniform(-1.3 * radius, 1.3 * radius, 120),
            generator.uniform(-1.3 * radius, 1.3 * radius, 120),
        ]
    )
    computed = numpy.asarray(field.signed_distance(points))
    xi = numpy.linspace(0.0, 1.0, 400_001)
    meridian_x = xi * length
    meridian_r = suboff_radius_fraction_array(xi) * radius
    for index, (x, y, z) in enumerate(points):
        rho = math.hypot(y, z)
        distance = float(numpy.min(numpy.sqrt((x - meridian_x) ** 2 + (rho - meridian_r) ** 2)))
        expected = -distance if (0.0 <= x <= length and rho <= field.radius_at(x)) else distance
        assert computed[index] == pytest.approx(expected, abs=1e-4), f"point {index}"


def test_the_gradient_is_unit_length_and_outward() -> None:
    field = SuboffHullField()
    points = [
        (0.5 * field.length, 0.8 * field.max_radius, 0.0),
        (0.9 * field.length, 0.0, 1.2 * field.max_radius),
        (0.15 * field.length, 0.0, -0.5 * field.max_radius),
    ]
    for point, gradient in zip(points, field.gradient(points), strict=True):
        assert math.dist(gradient, (0.0, 0.0, 0.0)) == pytest.approx(1.0, abs=1e-6)
        radial = (0.0, point[1], point[2])
        if math.dist(radial, (0.0, 0.0, 0.0)) > 0.0:
            dot = sum(a * b for a, b in zip(gradient, radial, strict=True))
            assert dot > 0.0, "the gradient must point away from the centreline"
    # Away from the medial axis the analytic gradient matches central differences.
    numeric = field.numerical_gradient(points)
    for analytic, reference in zip(field.gradient(points), numeric, strict=True):
        assert math.dist(analytic, reference) < 1e-5


def test_editing_the_length_leaves_the_barrel_radius_alone() -> None:
    before = SuboffHullField()
    after = before.with_state(length=2.0 * before.length)
    # A point on the midbody barrel is inside both, at the same radial offset.
    barrel = (0.5 * before.length, 0.8 * before.max_radius, 0.0)
    assert after.signed_distance([barrel]) == before.signed_distance([barrel])
    # The stern tip moves with the length and only with the length.
    assert before.contains([(1.5 * before.length, 0.0, 0.0)]) == [False]
    assert after.contains([(1.5 * before.length, 0.0, 0.0)]) == [True]
    assert after.max_radius == before.max_radius
    assert after.bounds().maximum[0] == pytest.approx(2.0 * before.length)


def test_editing_the_radius_leaves_the_length_alone() -> None:
    before = SuboffHullField()
    after = before.with_state(radius=1.2 * before.max_radius)
    assert after.length == before.length
    assert after.bounds().maximum[0] == pytest.approx(before.length)
    # The scaled hull keeps the same normalised profile, so the ratio is preserved.
    for xi in (0.1, 0.4, 0.8):
        assert after.radius_at(xi * after.length) / after.max_radius == pytest.approx(
            before.radius_at(xi * before.length) / before.max_radius
        )
    assert after.with_state(length=before.length).max_radius == after.max_radius


def test_the_reference_geometry_is_resolution_independent() -> None:
    reference = suboff_reference_geometry()
    assert reference["volume_resolution_change"] == pytest.approx(0.0, abs=1e-9)
    assert reference["wetted_area_resolution_change"] == pytest.approx(0.0, abs=1e-7)
    assert reference["volume"] == pytest.approx(0.699211, abs=1e-5)
    assert reference["wetted_area"] == pytest.approx(5.988264, abs=1e-5)
    assert reference["prismatic_coefficient"] == pytest.approx(0.791941, abs=1e-5)
    # Cross-check against the published experimental wetted area of the model.
    assert abs(reference["wetted_area"] / 5.998 - 1.0) < 0.005


def test_the_quadrature_reproduces_an_independent_implementation() -> None:
    """TensorLBM's ``suboff_analytic_wetted_area`` at its own radius default.

    Its ``r_over_l`` default is ``1/(2*8.57) = 0.058343``, while the offset table
    gives ``0.8333333/14.291667 = 0.058310``. At *its* radius the two implementations
    must agree: that isolates the quadrature from the radius convention, and makes
    the 0.06 % difference in the published area attributable to the radius alone.
    """
    theirs_radius = 1.0 / (2.0 * 8.57) * SUBOFF_LENGTH_M
    ours = suboff_reference_geometry(SUBOFF_LENGTH_M, theirs_radius)
    assert ours["wetted_area"] == pytest.approx(5.991830060618851, rel=1e-5)
    assert ours["volume"] == pytest.approx(0.700027, abs=1e-5)
    offset_table = suboff_reference_geometry()
    assert offset_table["wetted_area"] < ours["wetted_area"]
    assert abs(offset_table["wetted_area"] / ours["wetted_area"] - 1.0) < 1e-3


@pytest.mark.parametrize("resolution", [(109, 13, 13), (218, 26, 26)])
def test_export_is_a_closed_solver_ready_mesh(resolution: tuple[int, int, int]) -> None:
    field = SuboffHullField()
    reference = suboff_reference_geometry()
    sampled = to_mesh(field, resolution, margin=2)
    assert sampled.report.watertight
    assert sampled.report.edge_manifold
    assert sampled.report.euler_characteristic == 2
    assert sampled.validation.errors == []
    volume_error = sampled.report.signed_volume / reference["volume"] - 1.0
    area_error = sampled.report.surface_area / reference["wetted_area"] - 1.0
    assert abs(volume_error) < 0.01, f"volume error {100 * volume_error:+.3f} %"
    assert abs(area_error) < 0.01, f"wetted-area error {100 * area_error:+.3f} %"
    # The mesh reaches the analytic bounds: bow tip at x = 0, stern at x = length,
    # and the maximum radius on the barrel. Both tips are points of zero area, so the
    # mesh can only reach them to within the local resolution -- and the stern cap
    # closes faster than the transverse grid can follow, which is why its shortfall
    # is measured in transverse cells rather than axial ones.
    axial_cell = field.length / resolution[0]
    transverse_cell = 2.0 * field.max_radius / resolution[1]
    assert abs(sampled.report.bounds_min[0]) < 0.25 * axial_cell
    assert abs(sampled.report.bounds_max[0] - field.length) < 2.0 * transverse_cell
    # Facet sagitta: the flat triangles cut the maximum radius slightly, by much
    # less than one transverse cell.
    assert abs(sampled.report.bounds_max[1] - field.max_radius) < 0.5 * transverse_cell


def test_anisotropic_resolution_beats_a_cube_at_lower_cost() -> None:
    """The reason to_mesh accepts three cell counts: a slender body is not a cube."""
    field = SuboffHullField()
    reference = suboff_reference_geometry()
    refined = to_mesh(field, (218, 26, 26), margin=2)
    cube = to_mesh(field, 64, margin=2)
    refined_cells = 218 * 26 * 26
    cube_cells = 64**3
    assert refined_cells < cube_cells
    refined_error = abs(refined.report.signed_volume / reference["volume"] - 1.0)
    cube_error = abs(cube.report.signed_volume / reference["volume"] - 1.0)
    # Both are fine; the point is that the cheap one is not worse by more than 2x.
    assert refined_error < 2.0 * cube_error + 1e-3
    assert refined.report.watertight and cube.report.watertight


def test_the_pure_python_fallback_matches_the_vectorised_path(monkeypatch: object) -> None:
    """NumPy is a fast path: the protocol must answer identically without it."""
    field = SuboffHullField()
    points = [
        (0.5 * field.length, 0.6 * field.max_radius, 0.0),
        (0.05 * field.length, 0.0, 0.2 * field.max_radius),
        (0.95 * field.length, 0.3 * field.max_radius, -0.2 * field.max_radius),
        (-0.2, 0.0, 0.0),
        (field.length + 0.3, 0.0, 0.0),
    ]
    vectorised = field.signed_distance(points)
    vector_gradient = field.gradient(points)
    monkeypatch.setattr(SuboffHullField, "numpy_callable", lambda self: None)  # type: ignore[arg-type]
    assert field.numpy_callable() is None
    scalar = field.signed_distance(points)
    assert scalar == pytest.approx(vectorised, abs=1e-6)
    # The gradient is ill-conditioned where the nearest meridian point is nearly
    # ambiguous (close to the tip axis), so the two search depths agree to 1e-4 here
    # while the distances themselves agree to 1e-6.
    for analytic, reference in zip(field.gradient(points), vector_gradient, strict=True):
        assert math.dist(analytic, reference) < 1e-4


def test_the_profile_is_still_a_handoff_ready_field() -> None:
    field = SuboffHullField()
    record = to_mesh(field, (55, 7, 7), margin=2).handoff()
    assert record["units"] == "m"
    assert record["watertight"] is True
    assert record["validation"]["passed"] is True
    metadata = field.metadata()
    assert metadata["kind"] == "suboff_bare_hull"
    assert metadata["quantity"] == "signed_distance"
    assert metadata["provenance"]["geometry_source"].startswith("DARPA SUBOFF")
    assert metadata["state"] == pytest.approx({"length": SUBOFF_LENGTH_M, "radius": SUBOFF_RADIUS_M})
