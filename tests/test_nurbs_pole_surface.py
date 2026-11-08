"""Closed NURBS patch with singular pole rings: the sphere anchor.

Why this file exists
--------------------
``docs/nurbs-surface.md`` anchors an *open* tube and a closed cylinder whose caps are
triangle fans added by hand. A closed patch that pinches **by itself** -- a control net
whose boundary row is one repeated point, so the tessellator has to weld a whole row
into a single vertex and drop the degenerate quads -- had no committed anchor at all.
That is exactly the case where the parametrisation is singular, and it is where three
separate things can silently go wrong:

* the pole ring welds only ``to the last bit`` (the evaluation at ``v = v_min`` differs
  by 1 ulp from sample to sample), so a shell that looks closed is not;
* a quad at the pole collapses to one triangle, so a tessellator that assumes two
  faces per quad mis-indexes every following quad;
* the *order* of the poles decides the sign of the surface normal, so a closed shell
  can be uniformly inward and still be watertight.

Every number below is recomputed here; ``docs/nurbs_pole_surface.md`` reports the same
quantities from ``examples/nurbs_sphere/run_example.py``.
"""
from __future__ import annotations

import math

import pytest

from tensorcad import mesh as meshlib
from tensorcad.export import read_stl, write_stl
from tensorcad.geometry import MeshBuilder
from tensorcad.nurbs import (
    NURBSError,
    NURBSSurface,
    SQRT_HALF,
    circle_curve,
    cylinder_surface,
    mesh_solid_report,
    sphere_surface,
    tessellation_report,
)

RADII = (1.0, 12.0)
DENSITIES = ((24, 8), (48, 16), (96, 32))


def sphere_volume(radius: float) -> float:
    return 4.0 / 3.0 * math.pi * radius**3


def sphere_area(radius: float) -> float:
    return 4.0 * math.pi * radius**2


# ---------------------------------------------------------------------------
# 1. the patch is the sphere, not a fit to it
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("radius", RADII)
def test_every_evaluation_lies_on_the_sphere(radius: float) -> None:
    surface = sphere_surface(radius)
    worst = 0.0
    for i in range(9):
        for j in range(9):
            point = surface.evaluate(i / 8.0, 2.0 * j / 8.0)
            worst = max(worst, abs(math.sqrt(sum(c * c for c in point)) - radius))
    # Round-off only: the surface is the rational quotient of two exact arcs.
    assert worst < 1e-13, f"|S| deviates from R by {worst:.3e}"


def test_the_control_net_carries_the_pole_rows() -> None:
    surface = sphere_surface(12.0)
    low = surface.singular_boundary_point("v")
    high = surface.singular_boundary_point("v", upper=True)
    assert low == (0.0, 0.0, -12.0), "the south pole is one repeated control point"
    assert high == (0.0, 0.0, 12.0)
    # The u direction is a genuine closed circle, not a pinched row.
    assert surface.singular_boundary_point("u") is None
    assert surface.singular_boundary_point("u", upper=True) is None


def test_an_open_cylinder_has_no_singular_rows() -> None:
    cylinder = cylinder_surface(12.0, 11.0)
    for axis in ("u", "v"):
        assert cylinder.singular_boundary_point(axis) is None
        assert cylinder.singular_boundary_point(axis, upper=True) is None


def test_singular_boundary_point_rejects_a_bad_axis() -> None:
    with pytest.raises(NURBSError, match="axis must be"):
        sphere_surface(1.0).singular_boundary_point("w")


def test_zero_and_negative_radius_are_rejected() -> None:
    for bad in (0.0, -3.0):
        with pytest.raises(NURBSError, match="radius must be positive"):
            sphere_surface(bad)


# ---------------------------------------------------------------------------
# 2. parametric area: 4 pi r^2 to round-off, by quadrature per knot span
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("radius", RADII)
def test_parametric_area_matches_the_closed_form(radius: float) -> None:
    surface = sphere_surface(radius)
    relative = surface.area() / sphere_area(radius) - 1.0
    assert abs(relative) < 1e-12, f"area error {relative:.3e}"


def test_area_quadrature_converges_to_the_round_off_floor() -> None:
    surface = sphere_surface(1.0)
    exact = sphere_area(1.0)
    errors = {
        order: abs(surface.area(order) / exact - 1.0) for order in (4, 8, 12, 16, 24)
    }
    assert errors[4] > 1e-6, "a 4-point rule must still show a real quadrature error"
    assert errors[8] < 1e-7
    assert errors[16] < 1e-13, "16 points per span is already at round-off"
    assert errors[24] < 1e-13


# ---------------------------------------------------------------------------
# 3. tessellation: the pole rows weld, so the shell is closed
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("u_samples,v_samples", DENSITIES)
def test_the_pole_ring_welds_to_a_single_vertex(u_samples: int, v_samples: int) -> None:
    surface = sphere_surface(12.0)
    tess = surface.tessellate(u_samples, v_samples)
    poles = [
        point
        for point in tess.mesh.vertices
        if abs(abs(point[2]) - 12.0) < 1e-12 and math.hypot(point[0], point[1]) < 1e-12
    ]
    assert len(poles) == 2, f"expected one vertex per pole, got {len(poles)}: {poles}"
    # Closed in u (the seam welds) and pinched at both poles:
    # (v_samples - 1) rings of u_samples distinct points, plus two poles.
    assert len(tess.mesh.vertices) == u_samples * (v_samples - 1) + 2
    assert tess.degenerate_skipped == 2 * u_samples, "one collapsed quad per pole quad"


@pytest.mark.parametrize("u_samples,v_samples", DENSITIES)
def test_the_tessellated_sphere_is_a_closed_outward_solid(
    u_samples: int, v_samples: int
) -> None:
    surface = sphere_surface(12.0)
    tess = surface.tessellate(u_samples, v_samples)
    report = mesh_solid_report(tess.mesh, volume_reference=sphere_volume(12.0))
    assert report["watertight"] is True
    assert report["boundary_edges"] == 0
    assert report["non_manifold_edges"] == 0
    assert report["edge_manifold"] is True
    assert report["components"] == 1
    assert report["signed_volume_m3"] > 0.0, "a pole-to-pole meridian must be outward"
    assert report["orientation_consistent"] is True
    assert meshlib.zero_area_faces(tess.mesh) == (), (
        "a zero-area face would hide a collapsed quad"
    )
    # Inscribed facets: the volume approaches 4/3 pi r^3 from below.
    assert -0.06 < report["volume_relative_error"] < 0.0


@pytest.mark.parametrize("u_samples,v_samples", DENSITIES)
def test_faces_per_quad_stays_aligned_with_the_face_list(
    u_samples: int, v_samples: int
) -> None:
    """The pole quads emit one triangle, so the quad->face map must not assume two."""
    surface = sphere_surface(12.0)
    tess = surface.tessellate(u_samples, v_samples)
    flat = [index for quad in tess.faces_per_quad for index in quad]
    assert len(flat) == len(tess.mesh.faces), (
        "the quad->face map and the face list disagree: "
        f"{len(flat)} indices for {len(tess.mesh.faces)} faces"
    )
    assert max(flat) == len(tess.mesh.faces) - 1
    for i in range(u_samples):
        assert len(tess.faces_of_quad((i, 0))) == 1, "pole quad collapses to one triangle"
        assert len(tess.faces_of_quad((i, v_samples - 1))) == 1
        for j in range(1, v_samples - 1):
            assert len(tess.faces_of_quad((i, j))) == 2


@pytest.mark.parametrize("u_samples,v_samples", DENSITIES)
def test_no_face_points_inward_and_normals_converge(
    u_samples: int, v_samples: int
) -> None:
    surface = sphere_surface(12.0)
    tess = surface.tessellate(u_samples, v_samples)
    report = tessellation_report(surface, tess)
    assert report["inward_faces"] == 0
    assert report["max_normal_angle_deg"] < 2.0
    assert report["watertight"] is True
    assert report["area_relative_error"] < 0.0, "inscribed facets underestimate the area"


def test_the_row_major_flip_counter_is_not_a_winding_defect() -> None:
    """``orientation_flips`` compares consecutive faces *in walk order*.

    The walk is row-major over the quad grid, so on a patch closed in ``u`` it jumps
    straight from the south-pole row to the north-pole row and reports one "flip" per
    u step. The topology-aware test is that every directed edge appears exactly once,
    i.e. the shell is consistently oriented.
    """
    surface = sphere_surface(12.0)
    tess = surface.tessellate(24, 8)
    report = tessellation_report(surface, tess)
    assert report["orientation_flips"] == 24 - 1, "one artifact flip per u step"
    directed: dict[tuple[int, int], int] = {}
    for face in tess.mesh.faces:
        for k in range(3):
            edge = (face[k], face[(k + 1) % 3])
            directed[edge] = directed.get(edge, 0) + 1
    for (a, b), count in directed.items():
        assert count == 1, f"edge {a}->{b} used {count} times"
        assert directed.get((b, a), 0) == 1, f"edge {a}->{b} has no opposite twin"


# ---------------------------------------------------------------------------
# 4. volume convergence and STL round-trip
# ---------------------------------------------------------------------------
def test_volume_converges_as_o_h_squared() -> None:
    surface = sphere_surface(12.0)
    exact = sphere_volume(12.0)
    errors = []
    for u_samples, v_samples in ((24, 8), (48, 16), (96, 32)):
        tess = surface.tessellate(u_samples, v_samples)
        errors.append(abs(meshlib.total_signed_volume(tess.mesh) / exact - 1.0))
    assert errors[1] / errors[0] == pytest.approx(0.25, rel=0.05), errors
    assert errors[2] / errors[1] == pytest.approx(0.25, rel=0.05), errors
    assert errors[0] > 1e-2, "the coarsest grid must still show a discretisation error"


def test_a_refined_tessellation_meets_the_volume_threshold() -> None:
    surface = sphere_surface(12.0)
    tess = surface.tessellate(192, 64)
    report = mesh_solid_report(tess.mesh, volume_reference=sphere_volume(12.0))
    assert report["watertight"] is True
    assert abs(report["volume_relative_error"]) < 1e-3, report["volume_relative_error"]


def test_stl_round_trip_keeps_the_shell_closed(tmp_path) -> None:
    surface = sphere_surface(12.0)
    tess = surface.tessellate(48, 16)
    path = tmp_path / "sphere.stl"
    write_stl(tess.mesh, str(path))
    recovered = read_stl(str(path))
    assert len(recovered.faces) == len(tess.mesh.faces)
    assert meshlib.is_watertight(recovered) is True
    exact = sphere_volume(12.0)
    # The STL text carries a finite number of digits, so the volume matches to that
    # precision and no closer.
    assert meshlib.total_signed_volume(recovered) / exact - 1.0 == pytest.approx(
        meshlib.total_signed_volume(tess.mesh) / exact - 1.0, abs=1e-6
    )


# ---------------------------------------------------------------------------
# 5. negative control: the pole order decides the sign, and says so
# ---------------------------------------------------------------------------
def _sphere_with_reversed_meridian(radius: float) -> NURBSSurface:
    """The same sphere with the meridian running north -> south."""
    circle = circle_curve(1.0, arcs=3)
    meridian = ((0.0, 1.0), (1.0, 1.0), (1.0, 0.0), (1.0, -1.0), (0.0, -1.0))
    weights = (1.0, SQRT_HALF, 1.0, SQRT_HALF, 1.0)
    net = tuple(
        tuple((radius * r * p[0], radius * r * p[1], radius * z) for r, z in meridian)
        for p in circle.control_points
    )
    combined = tuple(
        tuple(circle.weights[i] * weights[j] for j in range(len(meridian)))
        for i in range(len(circle.control_points))
    )
    return NURBSSurface(
        degree=(2, 2),
        knots=(circle.knots, (0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 2.0)),
        control_net=net,
        weights=combined,
    )


def test_the_reversed_meridian_is_watertight_but_inward() -> None:
    """Reversing the poles keeps the geometry and flips every normal: a real trap."""
    surface = _sphere_with_reversed_meridian(12.0)
    worst = max(
        abs(math.sqrt(sum(c * c for c in surface.evaluate(i / 8.0, 2.0 * j / 8.0))) - 12.0)
        for i in range(9)
        for j in range(9)
    )
    assert worst < 1e-12, "same sphere geometrically"
    assert surface.area() / sphere_area(12.0) - 1.0 == pytest.approx(0.0, abs=1e-12)
    tess = surface.tessellate(24, 8)
    report = mesh_solid_report(tess.mesh, volume_reference=sphere_volume(12.0))
    assert report["watertight"] is True
    assert report["signed_volume_m3"] < 0.0
    assert report["orientation_consistent"] is False
    assert report["orientation_flipped_faces"] == len(tess.mesh.faces)


# ---------------------------------------------------------------------------
# 6. the builder contract the pole quads rely on
# ---------------------------------------------------------------------------
def test_a_pole_quad_reports_the_face_it_actually_emitted() -> None:
    builder = MeshBuilder()
    pole = builder.add_vertex((0.0, 0.0, 1.0))
    left = builder.add_vertex((1.0, 0.0, 0.0))
    right = builder.add_vertex((0.0, 1.0, 0.0))
    # A pole quad has two coincident corners: only ``(a, c, d)`` survives.
    assert builder.add_quad(pole, pole, right, left) == 1
    assert len(builder.build().faces) == 1
    assert builder.degenerate_skipped == 1