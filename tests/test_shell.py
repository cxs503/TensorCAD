"""Shelling, held to two exact references and to its own refusals.

The module claims three things that can be checked without trusting it:

1. **An eroded box.** ``shell(box, t)`` must reproduce the analytic ``(L - 2 t)^3`` inner
   volume. This is the ``GAP-CAD-011`` acceptance criterion, run over the whole ``t/L`` sweep
   the criterion names (0.01 to 0.1).
2. **An independent erosion.** For a convex pentagon prism the eroded solid is also the
   intersection of the offset half-planes, which ``tests/shell_oracle.py`` computes by 2-D
   Sutherland-Hodgman clipping and multiplies by the shortened length. The box case cannot
   distinguish the vertex solve from "shrink along the normals" (three orthogonal planes make
   every weighting equal); the pentagon can, because its facets are not mutually orthogonal.
3. **A measured wall, not an assumed one.** ``min_wall_thickness`` is exercised on a
   *deliberately under-offset* surface: it must report the smaller distance it is given, which
   is what makes the same number meaningful when the module reports ``t`` for a correct shell.

The refusals are part of the contract, so they are tested as first-class behaviour: a plate
thinner than twice the thickness, a collapsing taper, an open surface, an inward-oriented
surface and a non-positive thickness all have to fail loudly and say why. Nothing here checks
that a *bad* mesh is repaired -- the module refuses rather than repairs, and the tests assert
that refusal.
"""

from __future__ import annotations

import itertools
import json
import math

import pytest
from shell_oracle import eroded_prism_volume, eroded_section, section_area

from tensorcad import mesh as meshlib
from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.section_gz import prism_mesh
from tensorcad.shell import (
    ShellError,
    hollow_body,
    min_wall_thickness,
    offset_inward,
    shell,
)

RECTANGLE = ((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))
PENTAGON = ((0.0, 0.0), (1.0, 0.0), (1.3, 0.6), (0.6, 1.2), (-0.2, 0.7))
THICKNESSES = tuple(0.01 * step for step in range(1, 11))


def box(length: float = 1.0) -> Mesh:
    """Unit box as a convex prism: section in (y, z), extruded along +x."""
    return prism_mesh(RECTANGLE, length)


def plate(thickness: float = 0.1) -> Mesh:
    return prism_mesh(((0.0, 0.0), (1.0, 0.0), (1.0, thickness), (0.0, thickness)), 1.0)


def taper() -> Mesh:
    """Quadrilateral prism whose wall thins from 0.5 m to 0.1 m along +y."""
    return prism_mesh(((0.0, 0.0), (1.0, 0.0), (1.0, 0.1), (0.0, 0.5)), 1.0)


def uv_sphere(radius: float, meridians: int, parallels: int) -> Mesh:
    """Latitude/longitude sphere: poles fanned, body quads split into two triangles each."""
    builder = MeshBuilder()
    north = builder.add_vertex((0.0, 0.0, radius))
    south = builder.add_vertex((0.0, 0.0, -radius))
    rings = []
    for ring in range(1, parallels):
        phi = math.pi * ring / parallels
        z = radius * math.cos(phi)
        rho = radius * math.sin(phi)
        rings.append(
            [
                builder.add_vertex(
                    (
                        rho * math.cos(2.0 * math.pi * column / meridians),
                        rho * math.sin(2.0 * math.pi * column / meridians),
                        z,
                    )
                )
                for column in range(meridians)
            ]
        )
    for column in range(meridians):
        builder.add_triangle(north, rings[0][column], rings[0][(column + 1) % meridians])
    for upper, lower in itertools.pairwise(rings):
        for column in range(meridians):
            nxt = (column + 1) % meridians
            builder.add_quad(upper[column], lower[column], lower[nxt], upper[nxt])
    for column in range(meridians):
        builder.add_triangle(south, rings[-1][(column + 1) % meridians], rings[-1][column])
    return builder.build()


class TestExactReferences:
    @pytest.mark.parametrize("thickness", THICKNESSES)
    def test_box_inner_volume_is_the_analytic_eroded_box(self, thickness: float) -> None:
        """The GAP-CAD-011 sweep: t/L from 0.01 to 0.1, analytic inner volume."""
        length = 1.0
        result = shell(box(length), thickness)
        expected_inner = (length - 2.0 * thickness) ** 3
        relative = abs(result.inner_volume - expected_inner) / expected_inner
        assert relative < 1e-12, f"inner volume off by {relative:.3e}"
        expected_shell = length**3 - expected_inner
        assert abs(result.shell_volume - expected_shell) / expected_shell < 1e-12
        assert abs(result.outer_volume - 1.0) < 1e-15

    def test_box_volume_error_is_round_off_not_just_inside_the_two_percent(self) -> None:
        """The 2 % acceptance threshold is loose; the achieved error is nine orders tighter."""
        worst = 0.0
        for thickness in THICKNESSES:
            result = shell(box(), thickness)
            expected = (1.0 - 2.0 * thickness) ** 3
            worst = max(worst, abs(result.inner_volume - expected) / expected)
        assert worst < 1e-12
        assert worst < 0.02 / 1e6

    @pytest.mark.parametrize("thickness", [0.02, 0.05, 0.10, 0.15])
    def test_pentagon_prism_matches_independent_half_plane_clipping(self, thickness: float) -> None:
        """A facet set that is *not* mutually orthogonal, against a 2-D clipping oracle."""
        length = 1.0
        result = shell(prism_mesh(PENTAGON, length), thickness)
        expected = eroded_prism_volume(PENTAGON, length, thickness)
        relative = abs(result.inner_volume - expected) / expected
        assert relative < 1e-12, f"inner volume off by {relative:.3e}"
        # The oracle is not a copy of the answer: the erosion really removes material.
        assert section_area(PENTAGON) > section_area(eroded_section(PENTAGON, thickness))

    def test_every_offset_plane_is_satisfied_at_a_pentagon_vertex(self) -> None:
        result = shell(prism_mesh(PENTAGON, 1.0), 0.1)
        assert result.max_plane_residual < 1e-12

    def test_sphere_shell_converges_to_the_analytic_spherical_shell(self) -> None:
        """Faceting is the only error left on a curved body, so it must fall at second order."""
        radius, thickness = 1.0, 0.2
        analytic = 4.0 / 3.0 * math.pi * (radius**3 - (radius - thickness) ** 3)
        errors = []
        for meridians, parallels in ((12, 6), (24, 12), (48, 24)):
            result = shell(uv_sphere(radius, meridians, parallels), thickness)
            errors.append(abs(result.shell_volume - analytic) / analytic)
        assert errors[0] > errors[1] > errors[2]
        assert errors[2] < 0.01
        for coarse, fine in itertools.pairwise(errors):
            assert coarse / fine > 3.0, f"observed order below 1.6: {errors}"

    def test_inner_volume_is_strictly_decreasing_in_the_thickness(self) -> None:
        volumes = [shell(box(), t).inner_volume for t in THICKNESSES]
        assert all(later < earlier for earlier, later in itertools.pairwise(volumes))


class TestTheWallIsMeasuredNotAssumed:
    @pytest.mark.parametrize("thickness", THICKNESSES)
    def test_measured_wall_equals_the_request_on_the_box(self, thickness: float) -> None:
        result = shell(box(), thickness)
        assert abs(result.min_wall_thickness - thickness) <= 1e-12 * thickness
        assert abs(result.wall_shortfall) <= 1e-12 * thickness

    @pytest.mark.parametrize("factor", [0.2, 0.5, 0.7, 1.0])
    def test_measurement_detects_an_under_offset_surface(self, factor: float) -> None:
        """Feed the measurement a surface offset by ``factor * t``: it must report that distance."""
        target = 0.1
        under = offset_inward(box(), factor * target)
        measured = min_wall_thickness(box(), under)
        assert abs(measured - factor * target) <= 1e-12 * target

    def test_measurement_on_a_sphere_is_the_request_to_round_off(self) -> None:
        result = shell(uv_sphere(1.0, 24, 12), 0.1)
        assert abs(result.min_wall_thickness - 0.1) < 1e-9


class TestRefusals:
    def test_plate_beyond_half_thickness_is_refused(self) -> None:
        with pytest.raises(ShellError, match="inverts the solid"):
            shell(plate(0.1), 0.08)

    def test_plate_exactly_at_half_thickness_is_refused(self) -> None:
        """The offset degenerates at exactly w/2; it must fail rather than emit slivers."""
        with pytest.raises(ShellError):
            shell(plate(0.1), 0.05)

    def test_thin_plate_inside_half_thickness_shells_with_the_analytic_volume(self) -> None:
        thickness = 0.04
        result = shell(plate(0.1), thickness)
        expected = (1.0 - 2.0 * thickness) ** 2 * (0.1 - 2.0 * thickness)
        assert abs(result.inner_volume - expected) / expected < 1e-12

    def test_collapsing_taper_is_refused(self) -> None:
        with pytest.raises(ShellError):
            shell(taper(), 0.1)

    def test_open_surface_is_refused(self) -> None:
        opened = Mesh(vertices=box().vertices, faces=box().faces[1:])
        with pytest.raises(ShellError, match="not closed"):
            shell(opened, 0.05)

    def test_inward_oriented_surface_is_refused_and_orient_outward_cures_it(self) -> None:
        flipped = Mesh(
            vertices=box().vertices,
            faces=tuple((a, c, b) for a, b, c in box().faces),
        )
        with pytest.raises(ShellError, match="not outward-oriented"):
            shell(flipped, 0.05)
        cured = meshlib.orient_outward(flipped).mesh
        result = shell(cured, 0.05)
        assert abs(result.inner_volume - 0.9**3) / 0.9**3 < 1e-12

    @pytest.mark.parametrize("thickness", [0.0, -0.1, math.nan, math.inf])
    def test_non_positive_or_non_finite_thickness_is_refused(self, thickness: float) -> None:
        with pytest.raises(ShellError, match="positive and finite"):
            shell(box(), thickness)

    def test_empty_mesh_is_refused(self) -> None:
        with pytest.raises(ShellError, match="empty"):
            shell(Mesh(vertices=(), faces=()), 0.05)


class TestTheExport:
    def test_hollow_body_is_watertight_and_reports_both_volumes(self) -> None:
        result = shell(box(), 0.1)
        body = hollow_body(result)
        assert meshlib.is_watertight(body.mesh)
        assert meshlib.is_edge_manifold(body.mesh)
        assert len(body.mesh.faces) == 2 * len(result.outer.faces)
        assert body.material_volume == pytest.approx(result.outer_volume - result.inner_volume)
        assert body.cavity_volume == pytest.approx(result.inner_volume)
        # The reversed inner component makes the surface bound the material, so its integral is
        # the wall volume -- not the sum a naive "both outward" export would produce.
        assert (
            abs(meshlib.total_signed_volume(body.mesh) - body.material_volume) < 1e-12 * body.material_volume
        )

    def test_topology_is_preserved_index_for_index(self) -> None:
        outer = box()
        result = shell(outer, 0.1)
        assert result.inner.faces == outer.faces
        assert len(result.inner.vertices) == len(outer.vertices)
        assert result.inner.vertices != outer.vertices

    def test_the_input_mesh_is_not_mutated(self) -> None:
        outer = box()
        before = (outer.vertices, outer.faces)
        shell(outer, 0.1)
        assert (outer.vertices, outer.faces) == before

    def test_result_dict_is_json_serialisable(self) -> None:
        payload = shell(box(), 0.1).as_dict()
        restored = json.loads(json.dumps(payload))
        assert restored["thickness"] == 0.1
        assert restored["outer_faces"] == len(box().faces)
        assert abs(restored["shell_volume"] - (1.0 - 0.8**3)) < 1e-12
