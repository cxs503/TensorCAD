"""Surface--surface intersection: analytic anchors, one per surface pair.

No number in this file is typed in from a report. Every reference is computed here, by
a route that does not share code with the module under test:

* the **meridian** and **truncated-arc** cases have elementary closed forms;
* the **ellipse** and the **Viviani curve** need a complete elliptic integral, evaluated
  by the arithmetic-geometric mean (a different algorithm from marching squares) and,
  for the Viviani curve, cross-checked against a composite midpoint quadrature of the
  arc-length integrand;
* the **circle** cases are ``2 pi r`` with ``r`` fixed by plane/cylinder geometry.

The acceptance criterion the module is judged on is ``GAP-CAD-002``'s: arc length within
3% of the analytic value, and every point on both surfaces. The measured errors are
five orders of magnitude inside that, which is why the assertions below use ``1e-4``.
"""
from __future__ import annotations

import math

import pytest

from tensorcad.geometry_field import SphereField
from tensorcad.nurbs import cylinder_surface, sphere_surface
from tensorcad.surface_intersection import (
    CylinderPatch,
    IntersectionCurve,
    IntersectionResult,
    PlanePatch,
    SpherePatch,
    SurfaceIntersectionError,
    cylinder_implicit,
    plane_implicit,
    sphere_implicit,
    surface_surface_intersection,
)


# --------------------------------------------------------------------------- #
# independent references
# --------------------------------------------------------------------------- #
def complete_elliptic_e(m: float) -> float:
    """``E(m) = integral_0^{pi/2} sqrt(1 - m sin^2 t) dt`` by the AGM of Gauss.

    Four lines of AGM, deliberately not the quadrature the module's caller would use:
    the point of a second route is that it fails differently.
    """

    a, b, c = 1.0, math.sqrt(1.0 - m), math.sqrt(m)
    total = 0.5 * c * c
    power = 1.0
    for _ in range(40):
        a_next = 0.5 * (a + b)
        b_next = math.sqrt(a * b)
        c = 0.5 * (a - b)
        power *= 2.0
        total += 0.5 * power * c * c
        a, b = a_next, b_next
        if abs(c) < 1e-18:
            break
    return (math.pi / (2.0 * a)) * (1.0 - total)


def midpoint_arc_length(integrand, count: int = 200_000) -> float:
    """Composite midpoint rule on ``[0, 2 pi]``: an O(h**2) oracle for a smooth period."""

    step = 2.0 * math.pi / count
    return step * sum(integrand((index + 0.5) * step) for index in range(count))


def relative_error(measured: float, reference: float) -> float:
    return abs(measured - reference) / abs(reference)


# --------------------------------------------------------------------------- #
# surface pairs with elementary answers
# --------------------------------------------------------------------------- #
def test_plane_through_the_cylinder_axis_gives_two_meridian_lines() -> None:
    curve = surface_surface_intersection(
        CylinderPatch(radius=1.0, height=4.0),
        plane_implicit((0.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
        samples=(96, 12),
    )
    assert len(curve.curves) == 2
    assert all(not branch.closed for branch in curve.curves)
    # Two meridian lines, each exactly the patch height: nothing to discretise.
    for branch in curve.curves:
        assert branch.length == pytest.approx(4.0, rel=1e-12)
    assert curve.total_length == pytest.approx(8.0, rel=1e-12)


def test_crossing_that_lands_exactly_on_a_grid_node_is_still_found() -> None:
    """``u = 0`` and ``u = pi`` are grid nodes for an even ``u`` sample count.

    The zero of ``g`` then sits on a node, where a sign-change test sees nothing on one
    side. The crossing has to be recovered from the edge that leaves the node, otherwise
    one of the two meridians silently disappears.
    """

    even = surface_surface_intersection(
        CylinderPatch(radius=1.0, height=4.0),
        plane_implicit((0.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
        samples=(96, 12),
    )
    odd = surface_surface_intersection(
        CylinderPatch(radius=1.0, height=4.0),
        plane_implicit((0.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
        samples=(95, 13),
    )
    assert len(even.curves) == len(odd.curves) == 2
    assert even.total_length == pytest.approx(odd.total_length, rel=1e-12)


def test_inclined_plane_cuts_an_ellipse_whose_perimeter_matches_the_elliptic_integral() -> None:
    radius, alpha = 1.0, math.radians(30.0)
    normal = (0.0, math.sin(alpha), math.cos(alpha))
    result = surface_surface_intersection(
        CylinderPatch(radius=radius, height=4.0),
        plane_implicit((0.0, 0.0, 0.0), normal),
        samples=(128, 24),
    )
    assert len(result.curves) == 1
    assert result.curves[0].closed, "the ellipse is a closed loop across the u seam"
    # Semi-axes R and R / cos(alpha); perimeter 4 a E(1 - (b/a)**2).
    semi_major = radius / math.cos(alpha)
    reference = 4.0 * semi_major * complete_elliptic_e(1.0 - math.cos(alpha) ** 2)
    assert relative_error(result.total_length, reference) < 1e-4


def test_sphere_patch_cut_by_a_plane_gives_a_circle() -> None:
    radius, height = 1.0, 0.5
    result = surface_surface_intersection(
        SpherePatch(radius=radius),
        plane_implicit((0.0, 0.0, height), (0.0, 0.0, 1.0)),
        samples=(96, 48),
    )
    assert len(result.curves) == 1 and result.curves[0].closed
    reference = 2.0 * math.pi * math.sqrt(radius * radius - height * height)
    assert relative_error(result.total_length, reference) < 1e-4


def test_cylinder_cut_by_a_concentric_sphere_gives_two_circles() -> None:
    cylinder_radius, sphere_radius = 0.5, 1.0
    result = surface_surface_intersection(
        CylinderPatch(radius=cylinder_radius, height=3.0),
        sphere_implicit(sphere_radius),
        samples=(96, 24),
    )
    assert len(result.closed_curves) == 2
    assert relative_error(result.total_length, 4.0 * math.pi * cylinder_radius) < 1e-4


def test_perpendicular_equal_cylinders_trace_the_viviani_curve() -> None:
    """Two closed branches, each of length ``4 sqrt(2) E(1/2) R``.

    This is the case ``GAP-CAD-002`` names, and it is the hard one: on the parametric
    cylinder ``z = +/- R sin u``, so the two branches cross at ``(u, v) = (0, 0)`` where
    ``grad g`` vanishes. The level set is an X in parameter space and the marching
    squares cell containing it is the four-crossing saddle case.
    """

    radius = 1.0
    result = surface_surface_intersection(
        CylinderPatch(radius=radius, height=3.0),
        cylinder_implicit(radius, axis="x"),
        samples=(128, 32),
    )
    reference_branch = 4.0 * math.sqrt(2.0) * complete_elliptic_e(0.5) * radius
    quadrature = midpoint_arc_length(
        lambda t: radius * math.sqrt(1.0 + math.cos(t) ** 2)
    )
    assert reference_branch == pytest.approx(quadrature, rel=1e-9), "the two oracles agree"

    assert len(result.curves) == 2
    assert all(branch.closed for branch in result.curves)
    for branch in result.curves:
        assert relative_error(branch.length, reference_branch) < 1e-4
    assert relative_error(result.total_length, 2.0 * reference_branch) < 1e-4


def test_every_point_lies_on_both_surfaces() -> None:
    radius = 1.0
    result = surface_surface_intersection(
        CylinderPatch(radius=radius, height=3.0),
        cylinder_implicit(radius, axis="x"),
        samples=(128, 32),
    )
    scale = 2.0 * radius  # L of the problem
    assert result.max_residual <= 1e-9 * scale
    for branch in result.curves:
        for x, y, z in branch.points:
            # On the parametric cylinder and on the implicit one, checked independently.
            assert abs(math.hypot(x, y) - radius) <= 1e-9 * scale
            assert abs(math.hypot(y, z) - radius) <= 1e-9 * scale
    assert max(branch.closure_defect for branch in result.closed_curves) <= 1e-9 * scale


# --------------------------------------------------------------------------- #
# NURBS patches, seams and boundaries
# --------------------------------------------------------------------------- #
def test_nurbs_cylinder_cut_by_a_plane_closes_across_its_duplicated_seam() -> None:
    radius = 1.0
    surface = cylinder_surface(radius, 2.0, base_z=-1.0)
    result = surface_surface_intersection(
        surface, plane_implicit((0.0, 0.0, 0.0), (0.0, 0.0, 1.0)), samples=(96, 16)
    )
    assert len(result.curves) == 1
    assert result.curves[0].closed, "a clamped patch closes by position, not by edge key"
    assert relative_error(result.total_length, 2.0 * math.pi * radius) < 1e-4


def test_nurbs_sphere_cut_by_a_plane_gives_a_circle() -> None:
    height = 0.3
    surface = sphere_surface(1.0)
    result = surface_surface_intersection(
        surface, plane_implicit((0.0, 0.0, height), (0.0, 0.0, 1.0)), samples=(96, 48)
    )
    reference = 2.0 * math.pi * math.sqrt(1.0 - height * height)
    assert relative_error(result.total_length, reference) < 1e-4


def test_a_curve_leaving_the_patch_is_truncated_and_reported_open() -> None:
    radius = 0.75
    patch = PlanePatch(origin=(0.0, 0.0, 0.0), du=(3.0, 0.0, 0.0), dv=(0.0, 2.0, 0.0))
    result = surface_surface_intersection(patch, cylinder_implicit(radius), samples=(64, 64))
    assert len(result.curves) == 1
    branch = result.curves[0]
    assert not branch.closed
    # Only the first-quadrant quarter of the circle is inside x, y in [0, 3] x [0, 2].
    assert relative_error(branch.length, math.pi * radius / 2.0) < 1e-4


def test_geometry_field_can_be_the_implicit_side() -> None:
    """The package's batch ``signed_distance`` API is adapted, not demanded of the caller."""

    result = surface_surface_intersection(
        SpherePatch(radius=1.0),
        SphereField(radius=1.0, center=(0.0, 0.0, 1.0)),
        samples=(96, 48),
    )
    assert relative_error(result.total_length, 2.0 * math.pi * math.sqrt(0.75)) < 1e-4


# --------------------------------------------------------------------------- #
# convergence and honesty about what is not found
# --------------------------------------------------------------------------- #
def test_length_converges_second_order_and_refinement_pays() -> None:
    reference = 2.0 * math.pi * math.sqrt(1.0 - 0.25)
    errors = []
    for samples in (48, 96, 192):
        result = surface_surface_intersection(
            SpherePatch(radius=1.0),
            plane_implicit((0.0, 0.0, 0.5), (0.0, 0.0, 1.0)),
            samples=(samples, samples // 2),
            refine_levels=2,
        )
        errors.append(relative_error(result.total_length, reference))
    assert errors[1] < errors[0] and errors[2] < errors[1], errors
    order = math.log2(errors[0] / errors[2])
    assert order > 1.8, f"chord error should halve fourfold, measured order {order:.2f}"


def test_refinement_reduces_the_chord_underestimate() -> None:
    """A chord sum is always short; the refined length has to be closer, not just longer."""

    reference = 2.0 * math.pi * math.sqrt(0.75)
    coarse = surface_surface_intersection(
        SpherePatch(radius=1.0),
        plane_implicit((0.0, 0.0, 0.5), (0.0, 0.0, 1.0)),
        samples=(48, 24),
        refine_levels=0,
    )
    refined = surface_surface_intersection(
        SpherePatch(radius=1.0),
        plane_implicit((0.0, 0.0, 0.5), (0.0, 0.0, 1.0)),
        samples=(48, 24),
        refine_levels=2,
    )
    assert coarse.total_length < refined.total_length
    assert coarse.total_chord_length == pytest.approx(coarse.total_length, rel=1e-12)
    assert relative_error(refined.total_length, reference) < relative_error(
        coarse.total_length, reference
    )


def test_separated_surfaces_intersect_nowhere() -> None:
    result = surface_surface_intersection(
        SpherePatch(radius=1.0), plane_implicit((0.0, 0.0, 2.0), (0.0, 0.0, 1.0))
    )
    assert result.curves == ()
    assert result.total_length == 0.0
    assert result.max_residual == 0.0


def test_tangential_contact_is_not_reported_as_a_curve() -> None:
    """Documented limitation, asserted so that it cannot regress into a fake answer.

    A plane touching a sphere meets it in a *point*. There is no sign change, so the
    march cannot see it; what comes back is a rounding-scale trace around the pole, and
    the test pins that trace to ``1e-6 L`` rather than pretending it is a curve.
    """

    result = surface_surface_intersection(
        SpherePatch(radius=1.0),
        plane_implicit((0.0, 0.0, 1.0), (0.0, 0.0, 1.0)),
        samples=(64, 32),
    )
    assert result.total_length < 1e-6


def test_coincident_surfaces_raise_instead_of_returning_a_cloud() -> None:
    with pytest.raises(SurfaceIntersectionError, match="coincide"):
        surface_surface_intersection(
            SpherePatch(radius=1.0), SphereField(radius=1.0, center=(0.0, 0.0, 0.0))
        )


# --------------------------------------------------------------------------- #
# contract and inputs
# --------------------------------------------------------------------------- #
def test_reported_parameters_stay_inside_the_domain() -> None:
    surface = CylinderPatch(radius=1.0, height=4.0)
    result = surface_surface_intersection(
        surface, plane_implicit((0.0, 0.0, 0.0), (0.0, 0.0, 1.0)), samples=(32, 8)
    )
    (u_lo, u_hi), (v_lo, v_hi) = surface.domain
    for branch in result.curves:
        for u, v in branch.params:
            assert u_lo <= u <= u_hi and v_lo <= v <= v_hi


def test_result_and_curve_invariants() -> None:
    result = surface_surface_intersection(
        CylinderPatch(radius=1.0, height=3.0), sphere_implicit(1.0), samples=(48, 12)
    )
    assert isinstance(result, IntersectionResult)
    assert all(isinstance(branch, IntersectionCurve) for branch in result.curves)
    assert result.point_count == sum(len(branch.points) for branch in result.curves)
    assert result.samples == (48, 12)
    assert result.refine_levels == 2
    assert len(result.closed_curves) + len(result.open_curves) == len(result.curves)


def test_curve_rejects_mismatched_points_and_params() -> None:
    with pytest.raises(SurfaceIntersectionError):
        IntersectionCurve(
            points=((0.0, 0.0, 0.0),),
            params=((0.0, 0.0), (1.0, 1.0)),
            closed=False,
            length=0.0,
            chord_length=0.0,
            max_residual=0.0,
        )


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"samples": (1, 4)}, "at least"),
        ({"samples": (4,)}, "pair"),
        ({"refine_levels": 4}, "refine_levels"),
        ({"tolerance": 0.0}, "tolerance"),
        ({"tolerance": float("nan")}, "tolerance"),
    ],
)
def test_bad_sampling_arguments_are_rejected(kwargs, message) -> None:
    with pytest.raises(SurfaceIntersectionError, match=message):
        surface_surface_intersection(
            CylinderPatch(), plane_implicit((0.0, 0.0, 0.0)), **kwargs
        )


@pytest.mark.parametrize(
    "factory, message",
    [
        (lambda: CylinderPatch(radius=0.0), "radius"),
        (lambda: CylinderPatch(height=-1.0), "height"),
        (lambda: SpherePatch(radius=0.0), "radius"),
        (lambda: PlanePatch(du=(0.0, 0.0, 0.0)), "degenerate"),
        (lambda: sphere_implicit(0.0), "radius"),
        (lambda: cylinder_implicit(-1.0), "radius"),
        (lambda: cylinder_implicit(1.0, axis="w"), "axis"),
        (lambda: plane_implicit((0.0, 0.0, 0.0), (0.0, 0.0, 0.0)), "non-zero"),
    ],
)
def test_degenerate_geometry_is_rejected(factory, message) -> None:
    with pytest.raises(SurfaceIntersectionError, match=message):
        factory()


def test_periodic_flags_match_the_parameterisation() -> None:
    assert CylinderPatch().periodic == (True, False)
    assert SpherePatch().periodic == (True, False)
    assert PlanePatch().periodic == (False, False)


def test_implicit_must_be_usable() -> None:
    with pytest.raises(SurfaceIntersectionError, match="callable"):
        surface_surface_intersection(CylinderPatch(), object())


def test_surface_must_expose_a_domain_and_an_evaluate() -> None:
    class Bare:
        def evaluate(self, u, v):  # pragma: no cover - never reached
            return (u, v, 0.0)

    with pytest.raises(SurfaceIntersectionError, match="domain"):
        surface_surface_intersection(Bare(), plane_implicit((0.0, 0.0, 0.0)))

    class Degenerate:
        domain = ((1.0, 1.0), (0.0, 1.0))

        def evaluate(self, u, v):  # pragma: no cover - never reached
            return (u, v, 0.0)

    with pytest.raises(SurfaceIntersectionError, match="empty"):
        surface_surface_intersection(Degenerate(), plane_implicit((0.0, 0.0, 0.0)))

    with pytest.raises(SurfaceIntersectionError, match="evaluate"):
        surface_surface_intersection(object(), plane_implicit((0.0, 0.0, 0.0)))