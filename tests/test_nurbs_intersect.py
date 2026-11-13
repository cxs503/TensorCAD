"""NURBS x NURBS intersection: analytic anchors, one per pair, plus the honest failures.

No reference number in this file is typed in from a report; every one is computed here by
a route that shares no code with the module under test:

* the **circle** radii come from the two-sphere / plane-sphere geometry itself
  (``r = sqrt(R**2 - h**2)``, ``r = sqrt(R1**2 - x**2)``), so the check is not "does the
  march agree with the march" but "does the march land on the analytic circle";
* the **ellipse** perimeter needs a complete elliptic integral ``E(m)``, evaluated by the
  arithmetic-geometric mean -- a different algorithm from marching squares;
* the **trim areas** are compared with exact parameter-space integrals (a latitude band on
  a lat/long sphere, a disc on a plane patch whose parameters are orthonormal);
* the **cross-check** runs :func:`tensorcad.surface_intersection.surface_surface_intersection`
  -- a parametric-x-implicit marching-squares implementation with its own grid, its own
  chaining and no 4-D Newton at all -- on the same pair and compares lengths.

The acceptance criterion the module is judged on is ``GAP-CAD-A1-NURBS-INTERSECT``'s:
points on **both** surfaces, a circle radius within ``1e-8 R`` of the analytic value, and
the cross-check agreeing. The measured radius error is around ``1e-15 R``, six orders
inside that, which is why the assertions below can be that tight.
"""
from __future__ import annotations

import math

import pytest

from tensorcad.nurbs import cylinder_surface, sphere_surface
from tensorcad.nurbs_intersect import (
    DEFAULT_SAMPLES,
    NURBSIntersectionError,
    nurbs_surface_intersection,
)
from tensorcad.surface_intersection import (
    CylinderPatch,
    PlanePatch,
    SpherePatch,
    plane_implicit,
    sphere_implicit,
    surface_surface_intersection,
)


# --------------------------------------------------------------------------- #
# independent references
# --------------------------------------------------------------------------- #
def complete_elliptic_e(m: float) -> float:
    """``E(m)`` for ``0 <= m < 1`` by the arithmetic-geometric mean.

    ``K = pi / (2 a_N)`` and ``E = K (1 - sum_{n>=0} 2**(n-1) c_n**2)``; the AGM converges
    quadratically, so twenty terms is round-off. ``E(0) = pi/2`` and ``E(0.5) = 1.3506438``
    are pinned in the test that uses it, so the reference itself is checked.
    """

    a, b = 1.0, math.sqrt(1.0 - m)
    c = math.sqrt(m)
    total = 0.5 * c * c
    weight = 1.0
    for _ in range(20):
        a_next = 0.5 * (a + b)
        b_next = math.sqrt(a * b)
        c = 0.5 * (a - b)
        a, b = a_next, b_next
        total += weight * c * c
        weight *= 2.0
    return (math.pi / (2.0 * a)) * (1.0 - total)


def elliptic_e_quadrature(m: float, points: int = 200001) -> float:
    """``E(m) = int_0^{pi/2} sqrt(1 - m sin^2 t) dt`` by the midpoint rule."""

    total = 0.0
    half_pi = 0.5 * math.pi
    step = half_pi / points
    for i in range(points):
        t = (i + 0.5) * step
        total += math.sqrt(1.0 - m * math.sin(t) ** 2)
    return total * step


def sphere_sphere_circle(r1: float, r2: float, d: float) -> tuple[float, float]:
    """``(centre offset along the axis of centres, circle radius)`` for two spheres."""

    x = (d * d + r1 * r1 - r2 * r2) / (2.0 * d)
    return x, math.sqrt(max(0.0, r1 * r1 - x * x))


def radial_errors(points, centre_xy=(0.0, 0.0)) -> tuple[float, float]:
    """Radius of each point about the ``z`` axis, and the spread of ``z``."""

    radii = [math.hypot(p[0] - centre_xy[0], p[1] - centre_xy[1]) for p in points]
    zs = [p[2] for p in points]
    return max(radii), max(zs) - min(zs)


# --------------------------------------------------------------------------- #
# sphere x plane: the canonical circle, radius to 1e-8 R
# --------------------------------------------------------------------------- #
def test_sphere_plane_circle_radius_to_one_part_in_1e8():
    radius, height = 2.0, 0.5
    sphere = sphere_surface(radius)
    plane = PlanePatch(origin=(-3.0, -3.0, height), du=(6.0, 0.0, 0.0), dv=(0.0, 6.0, 0.0))

    result = nurbs_surface_intersection(sphere, plane)

    assert len(result.curves) == 1
    curve = result.curves[0]
    assert curve.closed
    assert result.diagnostics["dropped_crossings"] == 0.0

    reference = math.sqrt(radius * radius - height * height)
    measured, z_spread = radial_errors(curve.points)
    # The acceptance bound of the anchor is 1e-8 R; the measured error is ~1e-15 R, and
    # asserting the measured value (not just the bound) is what makes a regression visible.
    assert abs(measured - reference) <= 1e-8 * radius
    assert abs(measured - reference) <= 1e-12 * radius
    assert z_spread <= 1e-12 * radius

    # Every vertex is on *both* surfaces, which is the other half of the criterion.
    for x, y, z in curve.points:
        assert abs(math.sqrt(x * x + y * y + z * z) - radius) <= 1e-12 * radius
        assert abs(z - height) <= 1e-12 * radius
    assert curve.max_gap <= 1e-13 * radius

    reference_length = 2.0 * math.pi * reference
    assert abs(curve.length - reference_length) / reference_length < 1e-4


def test_length_converges_at_second_order():
    """The chord argument predicts ``O(h**2)``; measure it instead of asserting it."""

    radius = 1.0
    plane = PlanePatch(origin=(-2.0, -2.0, 0.25), du=(4.0, 0.0, 0.0), dv=(0.0, 4.0, 0.0))
    reference = 2.0 * math.pi * math.sqrt(radius * radius - 0.0625)

    errors = []
    for cells in (32, 64, 128):
        result = nurbs_surface_intersection(
            sphere_surface(radius), plane, samples=(cells, max(8, cells // 3)), refine_levels=1
        )
        errors.append(abs(result.total_length - reference) / reference)

    orders = [
        math.log(errors[i] / errors[i + 1]) / math.log(2.0)
        for i in range(len(errors) - 1)
    ]
    assert errors[-1] < errors[0]
    for order in orders:
        assert order > 1.5


# --------------------------------------------------------------------------- #
# sphere x sphere: two analytic cases
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "r1, r2, separation",
    [
        (1.0, 1.0, 1.2),   # equal radii
        (1.0, 1.2, 1.0),   # unequal, centre of the small one inside the large
    ],
)
def test_sphere_sphere_circle_radius(r1, r2, separation):
    offset, reference_radius = sphere_sphere_circle(r1, r2, separation)
    first = sphere_surface(r1)
    second = SpherePatch(radius=r2, centre=(0.0, 0.0, separation))

    result = nurbs_surface_intersection(first, second)

    assert len(result.curves) == 1
    curve = result.curves[0]
    assert curve.closed
    measured, z_spread = radial_errors(curve.points)
    assert abs(measured - reference_radius) <= 1e-8 * r1
    assert abs(measured - reference_radius) <= 1e-12 * r1
    assert abs(z_spread) <= 1e-12 * r1
    # The circle also has to *sit* at the analytic offset along the axis of centres.
    assert abs(sum(p[2] for p in curve.points) / len(curve.points) - offset) <= 1e-8 * r1

    reference_length = 2.0 * math.pi * reference_radius
    assert abs(curve.length - reference_length) / reference_length < 1e-4


# --------------------------------------------------------------------------- #
# the periodic-parameter path: a cylinder table, u wrapped
# --------------------------------------------------------------------------- #
def test_periodic_cylinder_closes_through_the_loop_path():
    cylinder = CylinderPatch(radius=1.0, height=2.0, centre=(0.0, 0.0, 0.0))
    plane = PlanePatch(origin=(-2.0, -2.0, 0.3), du=(4.0, 0.0, 0.0), dv=(0.0, 4.0, 0.0))

    result = nurbs_surface_intersection(cylinder, plane)

    assert len(result.curves) == 1
    curve = result.curves[0]
    assert curve.closed
    # A wrapped angle parameter closes through the cell topology, not on position.
    assert curve.closed_by == "loop"
    reference = 2.0 * math.pi
    assert abs(curve.length - reference) / reference < 1e-4


def test_clamped_nurbs_cylinder_closes_on_the_seam():
    cylinder = cylinder_surface(1.0, 4.0, base_z=-2.0)
    plane = PlanePatch(origin=(-2.0, -2.0, 0.3), du=(4.0, 0.0, 0.0), dv=(0.0, 4.0, 0.0))

    result = nurbs_surface_intersection(cylinder, plane)

    assert len(result.curves) == 1
    assert result.curves[0].closed
    # ``cylinder_surface`` is clamped, so the seam is two different parameters that
    # happen to coincide in space: closure has to be detected on position.
    assert result.curves[0].closed_by == "seam"


# --------------------------------------------------------------------------- #
# NURBS x NURBS with a genuinely tilted second patch: the ellipse
# --------------------------------------------------------------------------- #
def test_nurbs_cylinder_cut_by_inclined_plane_is_an_ellipse():
    tilt = math.radians(30.0)
    normal = (0.0, -math.sin(tilt), math.cos(tilt))
    # Two in-plane axes, built by hand so the reference integral is the only machinery.
    along = _unit(_cross(normal, (1.0, 0.0, 0.0)))
    across = _cross(normal, along)
    origin = tuple(-3.0 * along[i] - 3.0 * across[i] for i in range(3))
    plane = PlanePatch(origin=origin, du=tuple(6.0 * c for c in along), dv=tuple(6.0 * c for c in across))

    cylinder = cylinder_surface(1.0, 4.0, base_z=-2.0)
    result = nurbs_surface_intersection(cylinder, plane)

    # The reference integral is itself cross-checked: the AGM value against a midpoint
    # quadrature of the defining integral, which shares no step with it.
    assert abs(complete_elliptic_e(0.0) - math.pi / 2.0) < 1e-14
    assert abs(complete_elliptic_e(0.5) - elliptic_e_quadrature(0.5)) < 1e-9

    semi_minor, semi_major = 1.0, 1.0 / math.cos(tilt)
    reference = 4.0 * semi_major * complete_elliptic_e(1.0 - (semi_minor / semi_major) ** 2)
    assert len(result.curves) == 1
    assert result.curves[0].closed
    assert abs(result.total_length - reference) / reference < 1e-4


def _unit(v):
    length = math.sqrt(sum(c * c for c in v))
    return tuple(c / length for c in v)


def _cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


# --------------------------------------------------------------------------- #
# honest failure modes
# --------------------------------------------------------------------------- #
def test_tangential_contact_is_reported_not_invented():
    """A plane grazing a sphere meets it in a point, and the module says so."""

    sphere = sphere_surface(1.0)
    plane = PlanePatch(origin=(-2.0, -2.0, 1.0), du=(4.0, 0.0, 0.0), dv=(0.0, 4.0, 0.0))

    result = nurbs_surface_intersection(sphere, plane)

    assert result.curves == ()
    assert result.total_length == 0.0
    # A round-off cloud of crossings does appear; it is collapsed and counted rather than
    # returned as a zero-length "curve".
    assert result.diagnostics["collapsed_branches"] >= 1.0
    assert result.tangential_candidate


def test_disjoint_patches_return_nothing_and_are_not_flagged_tangential():
    sphere = sphere_surface(1.0)
    plane = PlanePatch(origin=(-2.0, -2.0, 5.0), du=(4.0, 0.0, 0.0), dv=(0.0, 4.0, 0.0))

    result = nurbs_surface_intersection(sphere, plane)

    assert result.curves == ()
    assert not result.tangential_candidate


def test_coincident_patches_raise_instead_of_returning_a_cloud():
    with pytest.raises(NURBSIntersectionError, match="coincide"):
        nurbs_surface_intersection(sphere_surface(1.0), sphere_surface(1.0))


def test_curve_truncated_by_a_small_patch_is_reported_open():
    sphere = sphere_surface(1.0)
    plane = PlanePatch(origin=(-0.5, -2.0, 0.4), du=(1.0, 0.0, 0.0), dv=(0.0, 4.0, 0.0))

    result = nurbs_surface_intersection(sphere, plane)

    assert len(result.curves) == 1
    curve = result.curves[0]
    assert not curve.closed
    assert curve.closed_by == ""
    # Both ends are still on the sphere; they are cut off by the *patch* boundary.
    for point in (curve.points[0], curve.points[-1]):
        assert abs(math.sqrt(sum(c * c for c in point)) - 1.0) <= 1e-9
        assert -0.5 - 1e-9 <= point[0] <= 0.5 + 1e-9


# --------------------------------------------------------------------------- #
# trim domains
# --------------------------------------------------------------------------- #
def test_trim_domain_on_a_plane_patch_is_the_analytic_disc():
    radius, height = 0.5, 0.2
    sphere = sphere_surface(radius)
    # Unit parameter axes on a unit patch: one parameter unit is one length unit, so the
    # parameter area *is* the surface area and ``pi r**2`` is a legitimate reference.
    plane = PlanePatch(origin=(-0.5, -0.5, height), du=(1.0, 0.0, 0.0), dv=(0.0, 1.0, 0.0))

    result = nurbs_surface_intersection(sphere, plane)

    reference_radius = math.sqrt(radius * radius - height * height)
    reference_area = math.pi * reference_radius ** 2

    area = result.trim_b.surface_area((24, 24))
    assert abs(area - reference_area) / reference_area < 2e-2
    assert abs(result.trim_b.area_uv((24, 24)) - reference_area) / reference_area < 2e-2

    # The disc's centre is kept, a corner of the patch is not.
    assert result.trim_b.contains(0.5, 0.5)
    assert not result.trim_b.contains(0.02, 0.02)

    # The loop polygon itself is a genuine closed polygon in this parameter space, and its
    # shoelace area is the inscribed-polygon area: within ``O(h**2)`` of the disc.
    assert len(result.trim_b.loops) == 1
    loop = result.trim_b.loops[0]
    assert loop.closed and loop.area_reliable and not loop.seam_spanning
    assert abs(loop.area_uv - reference_area) / reference_area < 5e-3
    assert loop.perimeter_uv < 2.0 * math.pi * reference_radius * 1.01


def test_trim_domain_on_the_sphere_reports_a_latitude_band_honestly():
    radius, height = 2.0, 0.5
    result = nurbs_surface_intersection(
        sphere_surface(radius),
        PlanePatch(origin=(-3.0, -3.0, height), du=(6.0, 0.0, 0.0), dv=(0.0, 6.0, 0.0)),
    )

    # The kept region is the part of the sphere below the plane: on this lat/long patch
    # that is a latitude band starting at the south pole, ``v in [0, v*]`` with
    # ``2 sin(pi (v-1)/2) = height``, and the parameter area is ``v*`` exactly.
    v_star = 1.0 + 2.0 * math.asin(height / radius) / math.pi
    assert abs(result.trim_a.area_uv((16, 32)) - v_star) / v_star < 5e-2

    # Surface area of that band: the sphere minus the cap, ``4 pi R^2 - 2 pi R (R - h)``.
    reference = 4.0 * math.pi * radius ** 2 - 2.0 * math.pi * radius * (radius - height)
    assert abs(result.trim_a.surface_area((16, 32)) - reference) / reference < 2e-2

    # A latitude line encloses no parameter area, so the shoelace number is flagged
    # unreliable rather than printed as if it meant something.
    loop = result.trim_a.loops[0]
    assert loop.seam_spanning
    assert not loop.area_reliable


# --------------------------------------------------------------------------- #
# cross-check against the independent implicit-field marcher
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "make_pair, implicit_factory",
    [
        (lambda: (sphere_surface(2.0),
                  PlanePatch(origin=(-3.0, -3.0, 0.5), du=(6.0, 0.0, 0.0), dv=(0.0, 6.0, 0.0))),
         lambda: plane_implicit((0.0, 0.0, 0.5), (0.0, 0.0, 1.0))),
        (lambda: (sphere_surface(2.0), SpherePatch(radius=1.2, centre=(0.0, 0.0, 1.0))),
         lambda: sphere_implicit(1.2, (0.0, 0.0, 1.0))),
        (lambda: (cylinder_surface(1.0, 4.0, base_z=-2.0),
                  PlanePatch(origin=(-2.0, -2.0, 0.3), du=(4.0, 0.0, 0.0), dv=(0.0, 4.0, 0.0))),
         lambda: plane_implicit((0.0, 0.0, 0.3), (0.0, 0.0, 1.0))),
    ],
)
def test_cross_check_against_the_implicit_marcher(make_pair, implicit_factory):
    """Two disjoint implementations, same pair, same length.

    ``surface_surface_intersection`` marches the zero set of ``f(S(u, v))`` on its own grid
    and never forms a 4-D Newton system, so agreement here is evidence about the *method*
    and not just about self-consistency.
    """

    surface, partner = make_pair()
    ours = nurbs_surface_intersection(surface, partner)
    theirs = surface_surface_intersection(
        surface, implicit_factory(), samples=DEFAULT_SAMPLES, refine_levels=2
    )

    assert ours.curves and theirs.curves
    assert abs(ours.total_length - theirs.total_length) / theirs.total_length < 1e-4


# --------------------------------------------------------------------------- #
# contract
# --------------------------------------------------------------------------- #
def test_samples_and_refinement_are_validated():
    sphere = sphere_surface(1.0)
    plane = PlanePatch(origin=(-2.0, -2.0, 0.2), du=(4.0, 0.0, 0.0), dv=(0.0, 4.0, 0.0))
    with pytest.raises(NURBSIntersectionError, match="samples"):
        nurbs_surface_intersection(sphere, plane, samples=(2, 2))
    with pytest.raises(NURBSIntersectionError, match="refine_levels"):
        nurbs_surface_intersection(sphere, plane, refine_levels=9)
    with pytest.raises(NURBSIntersectionError, match="closure_tolerance"):
        nurbs_surface_intersection(sphere, plane, closure_tolerance=0.0)


def test_surface_without_a_domain_is_rejected():
    class Nothing:
        def evaluate(self, u, v):
            return (u, v, 0.0)

    with pytest.raises(NURBSIntersectionError, match="domain"):
        nurbs_surface_intersection(Nothing(), sphere_surface(1.0))