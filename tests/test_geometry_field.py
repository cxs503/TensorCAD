"""The continuous-geometry protocol: exactness, edits, and the export bridge.

Every test here is a *deterministic* claim about an analytic primitive, because
that is what makes it a reference implementation: if the sphere's field is not
exactly ``|x - c| - r``, it cannot judge a learned field.
"""
from __future__ import annotations

import math
import pathlib

import pytest

from tensorcad.geometry_field import (
    IMPLICIT_FUNCTION,
    SIGNED_DISTANCE,
    CylinderField,
    EllipsoidField,
    FieldError,
    GeometryField,
    PlaneField,
    SphereField,
    surface_point,
    to_mesh,
)

PROBES = [
    (0.0, 0.0, 0.0),
    (0.3, 0.0, 0.0),
    (0.0, 0.4, -0.2),
    (-0.25, 0.1, 0.35),
    (0.9, -0.7, 0.5),
]


def distance_to_sphere(
    point: tuple[float, float, float], center: tuple[float, float, float], r: float
) -> float:
    return math.dist(point, center) - r


def test_sphere_field_is_the_exact_signed_distance() -> None:
    field = SphereField(radius=0.5, center=(0.1, -0.2, 0.05))
    values = field.signed_distance(PROBES)
    for point, value in zip(PROBES, values, strict=True):
        assert value == pytest.approx(distance_to_sphere(point, field.center, 0.5), abs=1e-12)
    assert field.quantity == SIGNED_DISTANCE
    assert field.contains([(0.1, -0.2, 0.05)])


def test_plane_field_normalises_its_normal_and_reports_the_offset() -> None:
    field = PlaneField(normal=(0.0, 0.0, 2.0), offset=0.25)
    assert field.normal == (0.0, 0.0, 1.0)
    values = field.signed_distance([(0.0, 0.0, 0.75), (0.0, 0.0, -0.25), (5.0, 5.0, 0.25)])
    assert values == pytest.approx([0.5, -0.5, 0.0], abs=1e-12)
    assert field.bounds() is None, "a plane is unbounded; that has to be visible"


def test_capped_cylinder_matches_its_closed_form_distance() -> None:
    field = CylinderField(radius=0.3, height=1.0)
    # Beside the barrel: radial distance only.
    assert field.signed_distance([(0.5, 0.0, 0.0)])[0] == pytest.approx(0.2)
    # Above the lid: axial distance only.
    assert field.signed_distance([(0.0, 0.0, 0.8)])[0] == pytest.approx(0.3)
    # Past a corner: distance to the rim circle.
    corner = math.hypot(0.4 - 0.3, 0.6 - 0.5)
    assert field.signed_distance([(0.4, 0.0, 0.6)])[0] == pytest.approx(corner)
    # On the axis inside: the nearest exit is the barrel (0.3), not the lid (0.5).
    assert field.signed_distance([(0.0, 0.0, 0.0)])[0] == pytest.approx(-0.3)


def test_ellipsoid_is_an_implicit_function_not_a_distance() -> None:
    field = EllipsoidField(semi_axes=(0.5, 0.3, 0.2))
    on_surface = field.signed_distance([(0.5, 0.0, 0.0), (0.0, 0.3, 0.0), (0.0, 0.0, 0.2)])
    assert on_surface == pytest.approx([0.0, 0.0, 0.0], abs=1e-12)
    assert field.quantity == IMPLICIT_FUNCTION, "the magnitude is not a distance and must not claim to be"
    centre = field.signed_distance([(0.0, 0.0, 0.0)])[0]
    assert centre == -1.0


@pytest.mark.parametrize(
    "field",
    [
        SphereField(radius=0.5, center=(0.05, -0.1, 0.2)),
        CylinderField(radius=0.3, height=1.0, axis=(1.0, 1.0, 0.5)),
        EllipsoidField(semi_axes=(0.5, 0.3, 0.2), center=(0.0, 0.1, 0.0)),
        PlaneField(normal=(1.0, -2.0, 0.5), offset=0.1),
    ],
)
def test_analytic_gradient_matches_central_differences(field: object) -> None:
    points = [(0.31, 0.17, -0.09), (-0.12, 0.44, 0.05), (0.05, 0.02, 0.61)]
    analytic = field.gradient(points)  # type: ignore[attr-defined]
    numeric = field.numerical_gradient(points)  # type: ignore[attr-defined]
    worst = max(
        abs(a - b)
        for pair_a, pair_b in zip(analytic, numeric, strict=True)
        for a, b in zip(pair_a, pair_b, strict=True)
    )
    assert worst < 1e-6, f"gradient and finite differences disagree by {worst:.2e}"


@pytest.mark.parametrize(
    "field",
    [
        SphereField(radius=0.5),
        CylinderField(radius=0.3, height=1.0),
        EllipsoidField(semi_axes=(0.5, 0.3, 0.2)),
        PlaneField(normal=(0.0, 0.0, 1.0), offset=0.0),
    ],
)
def test_batch_queries_equal_point_queries(field: object) -> None:
    batch = field.signed_distance(PROBES)  # type: ignore[attr-defined]
    single = [field.signed_distance([point])[0] for point in PROBES]  # type: ignore[attr-defined]
    assert batch == single, "the batch contract must not change the answer"
    assert field.signed_distance([]) == []  # type: ignore[attr-defined]


@pytest.mark.parametrize(
    "field",
    [
        SphereField(radius=0.5),
        CylinderField(radius=0.3, height=1.0),
        EllipsoidField(semi_axes=(0.5, 0.3, 0.2)),
        PlaneField(normal=(0.0, 0.0, 1.0), offset=0.0),
    ],
)
def test_numpy_path_agrees_with_the_pure_python_path(field: object) -> None:
    numpy = pytest.importorskip("numpy")
    evaluate = field.numpy_callable()  # type: ignore[attr-defined]
    assert evaluate is not None
    pure = field.signed_distance(PROBES)  # type: ignore[attr-defined]
    vectorised = [float(value) for value in evaluate(numpy.asarray(PROBES, dtype=numpy.float64))]
    assert vectorised == pytest.approx(pure, abs=1e-12)


def test_contains_follows_the_declared_sign_convention() -> None:
    field = SphereField(radius=0.5)
    assert field.contains([(0.0, 0.0, 0.0)]) == [True]
    assert field.contains([(0.6, 0.0, 0.0)]) == [False]
    assert field.contains([(0.5, 0.0, 0.0)]) == [False], "the surface itself is not inside"


# ---------------------------- the protocol itself --------------------------- #


class _BareSphere(GeometryField):
    """A field written against the protocol only: no NumPy, no shortcuts."""

    kind = "bare_sphere"

    def __init__(self, radius: float = 0.5) -> None:
        self.radius = radius
        self.units = "m"

    def signed_distance(self, points: object) -> list[float]:
        return [math.dist(point, (0.0, 0.0, 0.0)) - self.radius for point in points]  # type: ignore[union-attr]

    def gradient(self, points: object) -> list[tuple[float, float, float]]:
        return self.numerical_gradient(points)  # type: ignore[arg-type]

    def bounds(self) -> object:
        from tensorcad.geometry import BoundingBox

        return BoundingBox((-self.radius,) * 3, (self.radius,) * 3)

    def state(self) -> dict[str, float]:
        return {"radius": self.radius}

    def with_state(self, **values: float) -> _BareSphere:
        unknown = set(values) - {"radius"}
        if unknown:
            raise FieldError(f"{self.kind} has no control(s) {sorted(unknown)}")
        return _BareSphere(float(values.get("radius", self.radius)))

    def parameters(self) -> dict[str, object]:
        return {"radius": self.radius}


class _ComplementBall(GeometryField):
    """Inside is where the field is *positive*: the other sign convention."""

    kind = "complement_ball"
    quantity = IMPLICIT_FUNCTION
    inside_sign = "positive"

    def signed_distance(self, points: object) -> list[float]:
        return [0.5 - math.dist(point, (0.0, 0.0, 0.0)) for point in points]  # type: ignore[union-attr]

    def gradient(self, points: object) -> list[tuple[float, float, float]]:
        return [(-x, -y, -z) for x, y, z in self.numerical_gradient(points)]  # type: ignore[arg-type]


def test_the_base_protocol_refuses_to_guess() -> None:
    field = GeometryField()
    assert field.bounds() is None
    assert field.state() == {}
    assert field.parameters() == {}
    assert field.numpy_callable() is None
    with pytest.raises(NotImplementedError):
        field.signed_distance([(0.0, 0.0, 0.0)])
    with pytest.raises(NotImplementedError):
        field.gradient([(0.0, 0.0, 0.0)])
    with pytest.raises(NotImplementedError):
        field.with_state(radius=1.0)


def test_a_field_without_a_vectorised_path_still_exports() -> None:
    """NumPy is a fast path, not a requirement: this field samples in pure Python."""
    field = _BareSphere(0.5)
    assert field.numpy_callable() is None
    assert field.metadata()["bounds"]["size"] == [1.0, 1.0, 1.0]
    sampled = to_mesh(field, 24)
    exact = 4.0 / 3.0 * math.pi * 0.125
    assert sampled.report.watertight
    assert sampled.report.signed_volume == pytest.approx(exact, rel=1e-2)


def test_the_inside_sign_convention_is_honoured() -> None:
    complement = _ComplementBall()
    assert complement.contains([(0.0, 0.0, 0.0)]) == [True]
    assert complement.contains([(0.9, 0.0, 0.0)]) == [False]
    assert complement.metadata()["inside_sign"] == "positive"


def test_metadata_is_json_ready() -> None:
    import json

    payload = json.dumps(SphereField(radius=0.5).metadata())
    assert "signed_distance" in payload


# --------------------------------------------------------------------------- #
# edits: the acceptance criterion is that the target moves and nothing else does
# --------------------------------------------------------------------------- #


def test_editing_the_sphere_radius_does_not_move_the_centre() -> None:
    before = SphereField(radius=0.5, center=(0.2, 0.0, 0.0))
    after = before.with_state(radius=0.75)
    assert after.center == before.center
    assert after.radius == 0.75
    probe = (0.9, 0.0, 0.0)
    assert before.signed_distance([probe])[0] - after.signed_distance([probe])[0] == pytest.approx(0.25)
    assert after.bounds().maximum[0] - before.bounds().maximum[0] == pytest.approx(0.25)


def test_editing_the_cylinder_height_does_not_change_the_radius() -> None:
    before = CylinderField(radius=0.3, height=1.0)
    after = before.with_state(height=1.6)
    assert after.radius == before.radius
    radial_probe = (0.8, 0.0, 0.0)  # beside the barrel, axially inside both
    assert after.signed_distance([radial_probe]) == before.signed_distance([radial_probe])
    cap_probe = (0.0, 0.0, 0.6)  # above the short lid, under the tall one
    assert before.signed_distance([cap_probe])[0] > 0.0
    assert after.signed_distance([cap_probe])[0] < 0.0


def test_editing_one_ellipsoid_axis_leaves_the_other_two_alone() -> None:
    before = EllipsoidField(semi_axes=(0.5, 0.3, 0.2))
    after = before.with_state(semi_y=0.45)
    assert after.semi_axes == (0.5, 0.45, 0.2)
    assert after.state() == {"semi_x": 0.5, "semi_y": 0.45, "semi_z": 0.2}
    assert before.semi_axes == (0.5, 0.3, 0.2), "with_state must not mutate the original"


def test_editing_the_plane_offset_keeps_the_normal() -> None:
    before = PlaneField(normal=(0.0, 0.0, 2.0), offset=0.0)
    after = before.with_state(offset=0.4)
    assert after.normal == before.normal
    assert after.signed_distance([(0.0, 0.0, 0.0)])[0] == pytest.approx(-0.4)


@pytest.mark.parametrize(
    "field",
    [SphereField(), CylinderField(), EllipsoidField(), PlaneField()],
)
def test_an_unknown_control_is_rejected_instead_of_ignored(field: object) -> None:
    with pytest.raises(FieldError, match="no control"):
        field.with_state(thickness=1.0)  # type: ignore[attr-defined]


@pytest.mark.parametrize(
    "factory",
    [
        lambda: SphereField(radius=0.0),
        lambda: SphereField(radius=float("nan")),
        lambda: CylinderField(radius=0.3, height=-1.0),
        lambda: CylinderField(axis=(0.0, 0.0, 0.0)),
        lambda: EllipsoidField(semi_axes=(0.5, 0.0, 0.2)),
        lambda: PlaneField(normal=(0.0, 0.0, 0.0)),
    ],
)
def test_degenerate_parameters_fail_loudly(factory: object) -> None:
    with pytest.raises(FieldError):
        factory()  # type: ignore[operator]


def test_malformed_query_points_are_rejected() -> None:
    field = SphereField()
    with pytest.raises(FieldError, match="expected 3"):
        field.signed_distance([(0.0, 0.0)])  # type: ignore[list-item]
    with pytest.raises(FieldError, match="not finite"):
        field.signed_distance([(0.0, float("inf"), 0.0)])


def test_metadata_records_what_a_consumer_needs() -> None:
    field = SphereField(radius=0.5, center=(0.0, 0.0, 0.0))
    metadata = field.metadata()
    assert metadata["kind"] == "sphere"
    assert metadata["quantity"] == SIGNED_DISTANCE
    assert metadata["inside_sign"] == "negative"
    assert metadata["units"] == "m"
    assert metadata["coordinate_frame"] == "object_local"
    assert metadata["state"] == {"radius": 0.5}
    assert metadata["bounds"]["size"] == [1.0, 1.0, 1.0]
    assert metadata["provenance"]["representation"] == "analytic"


# --------------------------------------------------------------------------- #
# the export bridge: field -> closed mesh -> the same acceptance pipeline
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    ("field", "exact_volume"),
    [
        (SphereField(radius=0.5), 4.0 / 3.0 * math.pi * 0.125),
        (CylinderField(radius=0.3, height=1.0), math.pi * 0.09),
        (EllipsoidField(semi_axes=(0.5, 0.3, 0.2)), 4.0 / 3.0 * math.pi * 0.5 * 0.3 * 0.2),
    ],
)
def test_analytic_fields_export_to_a_closed_validated_mesh(field: object, exact_volume: float) -> None:
    sampled = to_mesh(field, 32)  # type: ignore[arg-type]
    assert sampled.report.watertight
    assert sampled.report.edge_manifold
    assert sampled.report.euler_characteristic == 2
    assert sampled.report.signed_volume == pytest.approx(exact_volume, rel=5e-3)
    assert sampled.validation.errors == []
    assert {check.name for check in sampled.quality} == {
        "mesh_min_angle",
        "mesh_aspect_ratio",
        "mesh_sliver_share",
    }
    assert sampled.mesh.units == field.units  # type: ignore[attr-defined]


def test_the_handoff_record_is_self_describing_and_reproducible() -> None:
    import json

    sampled = to_mesh(SphereField(radius=0.5), 24)
    record = sampled.handoff()
    assert record["schema"] == "tensorcad.geometry.handoff/1"
    assert record["units"] == "m"
    assert record["watertight"] is True
    assert record["triangles"] == len(sampled.mesh.faces)
    assert record["validation"]["passed"] is True
    assert set(record["quality"]) == {"mesh_min_angle", "mesh_aspect_ratio", "mesh_sliver_share"}
    json.dumps(record)  # must be serialisable for a manifest
    assert sampled.mesh_digest() == to_mesh(SphereField(radius=0.5), 24).mesh_digest()
    other = to_mesh(SphereField(radius=0.55), 24)
    assert other.mesh_digest() != sampled.mesh_digest()


def test_an_unbounded_field_cannot_be_sampled() -> None:
    with pytest.raises(FieldError, match="unbounded"):
        to_mesh(PlaneField(), 8)


def test_surface_point_lands_on_the_surface() -> None:
    for field, radius in [(SphereField(radius=0.5), 0.5), (EllipsoidField(), None)]:
        point = surface_point(field, (1.0, 0.4, 0.2))
        assert abs(field.signed_distance([point])[0]) < 1e-6
        if radius is not None:
            assert math.dist(point, (0.0, 0.0, 0.0)) == pytest.approx(radius, abs=1e-6)


def test_mesh_quality_reports_slivers_instead_of_hiding_them() -> None:
    """A sliver is a real property of marching tetrahedra, so it must be reported."""
    sampled = to_mesh(SphereField(radius=0.5), 24)
    by_name = {check.name: check for check in sampled.quality}
    share = by_name["mesh_sliver_share"]
    assert "triangles" in share.detail
    assert share.severity == "warning", "a valid solid can still be a poor discretisation"


def _degenerate_faces(mesh: object) -> int:
    """Faces with exactly zero area -- the signature of a grid-vertex coincidence."""
    count = 0
    for face in mesh.faces:  # type: ignore[attr-defined]
        a, b, c = (mesh.vertices[index] for index in face)  # type: ignore[attr-defined]
        ux, uy, uz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
        vx, vy, vz = c[0] - a[0], c[1] - a[1], c[2] - a[2]
        area2 = math.sqrt(
            (uy * vz - uz * vy) ** 2 + (uz * vx - ux * vz) ** 2 + (ux * vy - uy * vx) ** 2
        )
        if area2 < 1e-12:
            count += 1
    return count


@pytest.mark.parametrize(
    "field",
    [
        SphereField(radius=0.5),
        CylinderField(radius=0.3, height=1.0),
        EllipsoidField(semi_axes=(0.5, 0.3, 0.2)),
    ],
)
def test_sampling_does_not_put_the_surface_on_grid_vertices(field: object) -> None:
    """r = 0.3 in a +-0.3 box is exactly 20 cells per axis, and that used to hurt.

    Marching tetrahedra used to emit 8,108 exactly-zero-area faces for the
    cylinder at 40 cells per axis, because the level set passed through grid
    vertices. ``to_mesh`` now offsets its sampling box by an irrational fraction
    of a cell, so the count must be zero -- a solver that sees zero-area elements
    is entitled to reject the mesh.
    """
    sampled = to_mesh(field, 40)  # type: ignore[arg-type]
    assert _degenerate_faces(sampled.mesh) == 0
    assert sampled.report.watertight
    assert sampled.report.euler_characteristic == 2


# --------------------------------------------------------------------------- #
# the learned representation, queried through the very same protocol
# --------------------------------------------------------------------------- #

BUNDLE = (
    pathlib.Path(__file__).resolve().parents[1]
    / "docs/experiments/2026-10-24-neural-geometry-object/bundle"
)


def test_the_learned_revision_speaks_the_same_protocol() -> None:
    pytest.importorskip("torch")
    from tensorcad.geometry_field import TRUNCATED_SIGNED_DISTANCE
    from tensorcad.neural.objects import field_from_object

    field = field_from_object(BUNDLE)
    assert field.kind == "neural_geometry_object"
    assert field.quantity == TRUNCATED_SIGNED_DISTANCE, (
        "outside the truncation band the magnitude is not a distance and must say so"
    )
    assert field.state() == {"radius": 0.5}
    assert field.units == "m"
    assert field.bounds().sizes == pytest.approx((1.16, 1.16, 1.16))
    assert field.parameters()["weights_sha256"].startswith("fefbd7be")
    assert field.provenance()["representation"] == "neural_field"

    # The sphere is inside and the far corner is outside -- through the base API.
    assert field.contains([(0.0, 0.0, 0.0)]) == [True]
    assert field.contains([(1.4, 1.4, 1.4)]) == [False]
    assert field.signed_distance([]) == []

    # Editing the conditioning value needs no retraining. The +x axis is dented,
    # so probe a direction the edit did not touch.
    probe = (0.0, 0.52, 0.0)
    assert field.signed_distance([probe])[0] > 0.0, "0.52 is outside the fitted 0.5 sphere"
    wider = field.with_state(radius=0.55)
    assert wider.state()["radius"] == 0.55
    assert wider.signed_distance([probe])[0] < 0.0, "the conditioned edit really moved the surface"
    assert field.state()["radius"] == 0.5, "with_state must not mutate the original"
    with pytest.raises(FieldError, match="no control"):
        field.with_state(depth=0.01)


def test_the_learned_field_gradient_is_autograd_and_points_outward() -> None:
    """The learned field's gradient is taken with autograd, not a finite difference.

    A float32 network makes central differences unreliable below ~1e-4 of the
    object size; the adapter therefore overrides ``gradient`` with autograd and
    keeps the protocol's numerical version only as a fallback. The two must agree
    once the step is coarse enough to see past the rounding.
    """
    pytest.importorskip("torch")
    from tensorcad.neural.objects import field_from_object

    field = field_from_object(BUNDLE)
    points = [(0.35, 0.0, 0.0), (0.0, 0.35, 0.0), (0.0, 0.0, 0.35)]
    gradients = field.gradient(points)
    numeric = field.numerical_gradient(points, step=1e-4 * max(field.bounds().sizes))
    for point, gradient, reference in zip(points, gradients, numeric, strict=True):
        norm = math.dist(gradient, (0.0, 0.0, 0.0))
        dot = sum(a * b for a, b in zip(gradient, point, strict=True))
        assert dot / (norm * math.dist(point, (0.0, 0.0, 0.0))) > 0.8, (
            "the gradient must point away from the solid"
        )
        assert math.dist(gradient, reference) < 1e-3, "autograd and central differences disagree"
