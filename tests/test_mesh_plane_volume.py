"""The plane-closure surface integral: volume below a waterplane on an *open* shell.

An open hull sheet has no divergence-theorem volume, which is why the DTMB 5415 golden
case publishes its displacement from the station-area route and withholds the surface
integral. :func:`tensorcad.mesh.volume_below_plane` supplies the missing closure
analytically (waterplane disc plus vertical openings contribute nothing) and is tested
here against geometry whose integral is known in closed form -- so that when it is used on
a real shell, a disagreement is a statement about the shell, not about the integrator.
"""

from __future__ import annotations

import math

import pytest

from tensorcad import mesh as meshlib
from tensorcad.geometry import Mesh


def open_box(length: float, beam: float, height: float) -> Mesh:
    """Five faces of a rectangular box: everything but the deck at ``z = height``."""
    vertices = (
        (0.0, 0.0, 0.0),
        (length, 0.0, 0.0),
        (length, beam, 0.0),
        (0.0, beam, 0.0),
        (0.0, 0.0, height),
        (length, 0.0, height),
        (length, beam, height),
        (0.0, beam, height),
    )
    faces = (
        (0, 3, 2), (0, 2, 1),          # bottom, normal -z
        (0, 1, 5), (0, 5, 4),          # y = 0 wall
        (1, 2, 6), (1, 6, 5),          # x = L wall
        (2, 3, 7), (2, 7, 6),          # y = B wall
        (3, 0, 4), (3, 4, 7),          # x = 0 wall
    )
    return Mesh(vertices=vertices, faces=faces)


def flipped(mesh: Mesh) -> Mesh:
    return Mesh(
        vertices=mesh.vertices,
        faces=tuple(tuple(reversed(face)) for face in mesh.faces),
        units=mesh.units,
        coordinate_frame=mesh.coordinate_frame,
    )


def uv_sphere(radius: float = 1.0, n_theta: int = 48, n_phi: int = 24) -> Mesh:
    vertices = []
    for i in range(n_phi + 1):
        phi = math.pi * i / n_phi
        for j in range(n_theta):
            theta = 2.0 * math.pi * j / n_theta
            vertices.append(
                (
                    radius * math.sin(phi) * math.cos(theta),
                    radius * math.sin(phi) * math.sin(theta),
                    radius * math.cos(phi),
                )
            )
    faces = []
    for i in range(n_phi):
        for j in range(n_theta):
            a = i * n_theta + j % n_theta
            b = i * n_theta + (j + 1) % n_theta
            c = (i + 1) * n_theta + (j + 1) % n_theta
            d = (i + 1) * n_theta + j % n_theta
            faces.extend(((a, c, b), (a, d, c)))
    return meshlib.orient_outward(Mesh(vertices=tuple(vertices), faces=tuple(faces))).mesh


def test_open_box_matches_the_closed_form_exactly() -> None:
    """Box without a deck: the waterplane closure is the deck, so the answer is exact."""
    mesh = open_box(2.0, 3.0, 4.0)
    assert not meshlib.is_watertight(mesh)
    result = meshlib.volume_below_plane(mesh, 4.0)
    assert result.volume == pytest.approx(2.0 * 3.0 * 4.0, rel=1e-15)
    assert result.waterplane_area == pytest.approx(2.0 * 3.0, rel=1e-15)
    assert result.submerged_surface_area == pytest.approx(2.0 * 3.0 + 2 * (2.0 + 3.0) * 4.0, rel=1e-15)
    assert result.orientation == 1
    assert result.clipped_triangle_count == 0


def test_open_box_below_the_deck_clips_triangles_exactly() -> None:
    """A plane through the walls: 8 wall triangles are cut, the integrals stay exact."""
    mesh = open_box(2.0, 3.0, 4.0)
    result = meshlib.volume_below_plane(mesh, 2.0)
    assert result.clipped_triangle_count == 8
    assert result.volume == pytest.approx(2.0 * 3.0 * 2.0, rel=1e-15)
    assert result.waterplane_area == pytest.approx(2.0 * 3.0, rel=1e-15)


def test_winding_does_not_change_the_answer() -> None:
    """The sign of the shell's winding is resolved, not assumed."""
    mesh = open_box(2.0, 3.0, 4.0)
    inward = flipped(mesh)
    outward = meshlib.volume_below_plane(mesh, 4.0)
    inward_result = meshlib.volume_below_plane(inward, 4.0)
    assert outward.orientation == 1
    assert inward_result.orientation == -1
    assert inward_result.volume == pytest.approx(outward.volume, rel=1e-15)
    assert inward_result.waterplane_area == pytest.approx(outward.waterplane_area, rel=1e-15)


def test_closed_shell_has_no_waterplane_theorem_applies() -> None:
    """On a closed mesh the closure contributes nothing: volume only, waterplane area 0."""
    mesh = open_box(2.0, 3.0, 4.0)
    closed = Mesh(
        vertices=mesh.vertices,
        faces=mesh.faces + ((4, 5, 6), (4, 6, 7)),
    )
    assert meshlib.is_watertight(closed)
    result = meshlib.volume_below_plane(closed, 4.0)
    assert result.volume == pytest.approx(24.0, rel=1e-15)
    assert result.waterplane_area == pytest.approx(0.0, abs=1e-12)


def test_hemisphere_matches_the_polyhedral_limit() -> None:
    """Lower half of a unit sphere: V -> 2/3 pi, A -> pi as the tessellation refines."""
    coarse = meshlib.volume_below_plane(uv_sphere(1.0, 24, 12), 0.0)
    fine = meshlib.volume_below_plane(uv_sphere(1.0, 96, 48), 0.0)
    exact_volume = 2.0 / 3.0 * math.pi
    exact_area = math.pi
    assert fine.volume > coarse.volume  # an inscribed polyhedron converges from below
    assert abs(fine.volume / exact_volume - 1.0) < 5e-3
    assert abs(fine.waterplane_area / exact_area - 1.0) < 1e-2
    assert abs(fine.volume / exact_volume - 1.0) < abs(coarse.volume / exact_volume - 1.0)


def test_shell_entirely_above_the_plane_has_no_column() -> None:
    mesh = Mesh(vertices=((0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (0.0, 1.0, 2.0)), faces=((0, 1, 2),))
    result = meshlib.volume_below_plane(mesh, 0.0)
    assert result.volume == 0.0
    assert result.waterplane_area == 0.0
    assert result.clipped_triangle_count == 0


def test_non_finite_plane_is_refused() -> None:
    with pytest.raises(ValueError):
        meshlib.volume_below_plane(open_box(1.0, 1.0, 1.0), float("nan"))