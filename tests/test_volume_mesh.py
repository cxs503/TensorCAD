"""Tests for the body-fitted volume mesher (``tensorcad.volume_mesh``)."""
from __future__ import annotations

import math

import pytest

from tensorcad.hull import HullParameters, analytic_volume, build_mesh
from tensorcad.mesh import total_signed_volume
from tensorcad.volume_mesh import (
    DISTORTION_THRESHOLD,
    SectionGrid,
    VolumeMeshError,
    box_section_grid,
    build_volume_mesh,
    tet_scaled_jacobian,
    volume_mesh_from_hull,
)

LENGTH, BEAM, DRAFT = 10.0, 2.0, 0.8
ANALYTIC = analytic_volume(LENGTH, BEAM, DRAFT)


def hull_parameters(nx: int = 41, nz: int = 13) -> HullParameters:
    return HullParameters(length=LENGTH, beam=BEAM, draft=DRAFT, nx=nx, nz=nz)


# --------------------------------------------------------------------------- the metric


def test_regular_tetrahedron_scores_its_known_value():
    """The scaled Jacobian of a regular tetrahedron is 1/sqrt(2) at every corner."""
    a = 2.0
    vertices = (
        (0.0, 0.0, 0.0),
        (a, 0.0, 0.0),
        (a / 2.0, a * math.sqrt(3.0) / 2.0, 0.0),
        (a / 2.0, a * math.sqrt(3.0) / 6.0, a * math.sqrt(2.0 / 3.0)),
    )
    assert tet_scaled_jacobian(vertices, (0, 1, 2, 3)) == pytest.approx(1.0 / math.sqrt(2.0))
    # a flat tet has no volume: the worst corner is 0 whatever the winding
    flat = ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (1.0, 1.0, 0.0))
    assert tet_scaled_jacobian(flat, (0, 1, 2, 3)) == pytest.approx(0.0)
    # inverting the tet reports the inverted volume, not a smaller good one
    assert tet_scaled_jacobian(vertices, (0, 2, 1, 3)) < 0.0


def test_unit_cube_kuhn_tets_are_not_better_than_regular():
    cube = build_volume_mesh(
        box_section_grid(length=1.0, beam=1.0, draft=1.0, nx=2, nz=2),
        first_layer=0.25,
        layers=2,
        growth=1.0,
    )
    values = sorted(
        tet_scaled_jacobian(cube.vertices, tet) for tet in cube.tets
    )
    # every Kuhn tet is positive and none beats a regular tetrahedron
    assert values[0] > 0.0
    assert values[-1] <= 1.0 / math.sqrt(2.0) + 1e-12
    # a 1x1x1 cell split into 2x2x2 breadth cells of side 0.5: the cube's own Kuhn tets
    assert values[0] == pytest.approx(1.0 / math.sqrt(6.0) / 2.0 * 2.0, rel=0.75)


# --------------------------------------------------------------------------- exact case


def test_box_volume_first_layer_and_growth_are_exact():
    mesh = build_volume_mesh(
        box_section_grid(length=LENGTH, beam=BEAM, draft=DRAFT, nx=11, nz=7),
        first_layer=0.01,
        growth=1.2,
        layers=4,
    )
    assert mesh.volume() == pytest.approx(LENGTH * BEAM * DRAFT, rel=1e-14)
    assert mesh.tet_volume() == pytest.approx(mesh.volume(), rel=1e-14)
    column = mesh.columns[0]
    assert column[0] == pytest.approx(0.01, rel=0, abs=1e-15)
    assert column[1] / column[0] == pytest.approx(1.2, rel=1e-12)
    assert sum(column) == pytest.approx(BEAM / 2, rel=1e-14)
    # a box has orthogonal cells, so nothing is distorted however thin the layers are
    quality = mesh.quality()
    assert quality["hex"]["all"]["min"] == pytest.approx(1.0, abs=1e-12)
    assert quality["hex"]["all"]["below_threshold"] == 0
    assert quality["collapsed_cells"] == 0


def test_box_symmetry_places_the_centreline_node_exactly():
    mesh = build_volume_mesh(
        box_section_grid(length=1.0, beam=2.0, draft=1.0, nx=3, nz=3),
        first_layer=0.05,
        growth=1.3,
        layers=3,
    )
    half = mesh.shape[2] // 2
    for i in range(mesh.shape[0]):
        for j in range(mesh.shape[1]):
            assert mesh.vertices[mesh.node_id(i, j, half)][1] == 0.0
            for k in range(half):
                left = mesh.vertices[mesh.node_id(i, j, k)][1]
                right = mesh.vertices[mesh.node_id(i, j, 2 * half - k)][1]
                assert left == pytest.approx(-right, rel=1e-15, abs=1e-18)


# --------------------------------------------------------------------------- hull case


def test_volume_mesh_reproduces_the_hull_surface_volume():
    """The volume the cells fill must be the volume the surface encloses, not a resampling."""
    parameters = hull_parameters()
    mesh = volume_mesh_from_hull(parameters, first_layer=0.005, layers=4)
    surface, _ = build_mesh(parameters)
    assert mesh.volume() == pytest.approx(total_signed_volume(surface), rel=1e-13)
    # and the distance to the analytic volume is the surface discretisation error, which the
    # repository's own audit records as -0.88 % for this grid
    assert mesh.volume() / ANALYTIC - 1.0 == pytest.approx(-0.0088, abs=2e-4)


def test_volume_error_falls_with_refinement():
    errors = []
    for nx, nz in ((41, 13), (81, 21)):
        mesh = volume_mesh_from_hull(hull_parameters(nx, nz), first_layer=0.005, layers=4)
        errors.append(abs(mesh.volume() / ANALYTIC - 1.0))
    assert errors[1] < errors[0]
    assert errors[1] < 0.005


def test_the_split_is_a_closed_conforming_tiling():
    mesh = volume_mesh_from_hull(hull_parameters(), first_layer=0.005, layers=4)
    counts = mesh.face_multiplicities()
    assert max(counts.values()) <= 2
    assert mesh.is_conforming()
    assert len(mesh.boundary_faces()) == len(counts) - sum(
        1 for value in counts.values() if value == 2
    )
    # the positive-volume subset alone is also closed: the zero-volume tets that come from
    # collapsed wedges can be dropped without opening the tiling
    positive = [tet for tet in mesh.tets if tet_scaled_jacobian(mesh.vertices, tet) > 0.0]
    assert 0 < len(positive) < len(mesh.tets)
    subset: dict[tuple[int, int, int], int] = {}
    for tet in positive:
        for a in range(4):
            for b in range(a + 1, 4):
                for c in range(b + 1, 4):
                    face = tuple(sorted((tet[a], tet[b], tet[c])))
                    subset[face] = subset.get(face, 0) + 1
    assert max(subset.values()) <= 2


def test_collapsed_columns_are_counted_not_hidden():
    mesh = volume_mesh_from_hull(hull_parameters(), first_layer=0.005, layers=4)
    summary = mesh.summary()
    assert summary["collapsed_columns"] > 0
    assert summary["quality"]["collapsed_cells"] > 0
    assert summary["quality"]["collapsed_cells"] < summary["counts"]["hexes"]
    # the keel line of this hull is at zero breadth: nothing may be reported there as flat
    assert summary["quality"]["hex"]["core"]["min"] > DISTORTION_THRESHOLD


def test_wall_normal_clearance_matches_the_request():
    parameters = hull_parameters(81, 21)
    for first_layer in (0.002, 0.01):
        mesh = volume_mesh_from_hull(parameters, first_layer=first_layer, layers=3, growth=1.4)
        i, j, layers = (parameters.nx - 1) // 2, (parameters.nz - 1) // 2, 3
        wall = mesh.vertices[mesh.node_id(i, j, 2 * layers)]
        inner = mesh.vertices[mesh.node_id(i, j, 2 * layers - 1)]
        a = mesh.vertices[mesh.node_id(i + 1, j, 2 * layers)]
        b = mesh.vertices[mesh.node_id(i - 1, j, 2 * layers)]
        c = mesh.vertices[mesh.node_id(i, min(parameters.nz - 1, j + 1), 2 * layers)]
        d = mesh.vertices[mesh.node_id(i, max(0, j - 1), 2 * layers)]
        t1 = tuple(a[m] - b[m] for m in range(3))
        t2 = tuple(c[m] - d[m] for m in range(3))
        normal = (
            t1[1] * t2[2] - t1[2] * t2[1],
            t1[2] * t2[0] - t1[0] * t2[2],
            t1[0] * t2[1] - t1[1] * t2[0],
        )
        length = math.sqrt(sum(component * component for component in normal))
        normal = tuple(component / length for component in normal)
        gap = tuple(inner[m] - wall[m] for m in range(3))
        clearance = abs(sum(gap[m] * normal[m] for m in range(3)))
        assert clearance == pytest.approx(first_layer, rel=1e-9)


# --------------------------------------------------------------------------- output


def test_vtk_round_trip(tmp_path):
    mesh = build_volume_mesh(
        box_section_grid(length=2.0, beam=1.0, draft=1.0, nx=3, nz=3),
        first_layer=0.02,
        layers=2,
    )
    path = mesh.to_vtk(tmp_path / "box.vtk")
    lines = path.read_text().splitlines()
    assert lines[4] == f"POINTS {len(mesh.vertices)} double"
    n_points = len(mesh.vertices)
    n_cells = int(lines[5 + n_points].split()[1])
    assert n_cells == len(mesh.hexes)
    assert lines[6 + n_points].split()[0] == "8"

    points = [tuple(float(v) for v in line.split()) for line in lines[5 : 5 + n_points]]
    total = 0.0
    for line in lines[6 + n_points : 6 + n_points + n_cells]:
        ids = [int(v) for v in line.split()[1:]]
        p0, p1, p2, p3, p4, p5, p6, p7 = (points[i] for i in ids)
        for index in ((0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6), (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)):
            a, b, c, d = (points[ids[k]] for k in index)
            ab = [b[k] - a[k] for k in range(3)]
            ac = [c[k] - a[k] for k in range(3)]
            ad = [d[k] - a[k] for k in range(3)]
            total += (
                ab[0] * (ac[1] * ad[2] - ac[2] * ad[1])
                - ab[1] * (ac[0] * ad[2] - ac[2] * ad[0])
                + ab[2] * (ac[0] * ad[1] - ac[1] * ad[0])
            ) / 6.0
    assert total > 0.0
    assert total == pytest.approx(mesh.volume(), rel=1e-12)


def test_vtk_refuses_an_unknown_cell_type(tmp_path):
    mesh = build_volume_mesh(
        box_section_grid(length=1.0, beam=1.0, draft=1.0, nx=2, nz=2), first_layer=0.1, layers=1
    )
    with pytest.raises(VolumeMeshError):
        mesh.to_vtk(tmp_path / "x.vtk", cells="wedge")


# --------------------------------------------------------------------------- refusals


@pytest.mark.parametrize(
    "first_layer, growth, layers",
    [(-0.1, 1.2, 4), (0.0, 1.2, 4), (0.1, 0.9, 4), (0.1, 1.2, 0), (float("nan"), 1.2, 4)],
)
def test_layer_spec_is_validated(first_layer, growth, layers):
    with pytest.raises(VolumeMeshError):
        build_volume_mesh(
            box_section_grid(length=1.0, beam=1.0, draft=1.0),
            first_layer=first_layer,
            growth=growth,
            layers=layers,
        )


@pytest.mark.parametrize(
    "xs, zs, breadth",
    [
        ((0.0, 1.0), (0.0, 1.0), lambda i, j: math.nan),
        ((0.0, 1.0), (0.0, 1.0), lambda i, j: -1.0),
        ((0.0, 0.0), (0.0, 1.0), lambda i, j: 1.0),
        ((1.0, 0.0), (0.0, 1.0), lambda i, j: 1.0),
        ((0.0,), (0.0, 1.0), lambda i, j: 1.0),
        ((0.0, 1.0), (1.0,), lambda i, j: 1.0),
        ((0.0, float("inf")), (0.0, 1.0), lambda i, j: 1.0),
    ],
)
def test_section_grid_is_validated(xs, zs, breadth):
    with pytest.raises(VolumeMeshError):
        SectionGrid(xs=xs, zs=zs, breadth=breadth)


def test_node_and_cell_accessors_refuse_out_of_range():
    mesh = build_volume_mesh(
        box_section_grid(length=1.0, beam=1.0, draft=1.0, nx=3, nz=3), first_layer=0.05, layers=2
    )
    with pytest.raises(VolumeMeshError):
        mesh.node_id(3, 0, 0)
    with pytest.raises(VolumeMeshError):
        mesh.cell_nodes(2, 0, 0)


def test_zero_breadth_columns_build_wedges_but_do_not_invert():
    grid = SectionGrid(
        xs=(0.0, 1.0, 2.0), zs=(0.0, 1.0), breadth=lambda i, j: 0.0 if i in (0, 2) else 1.0
    )
    mesh = build_volume_mesh(grid, first_layer=0.05, layers=2)
    assert mesh.conforming if hasattr(mesh, "conforming") else mesh.is_conforming()
    assert mesh.counts()["zero_volume_tets"] > 0
    assert mesh.volume() > 0.0


def test_summary_is_json_safe():
    import json

    mesh = volume_mesh_from_hull(hull_parameters(), first_layer=0.005, layers=2)
    summary = mesh.summary()
    json.dumps(summary)
    assert summary["conforming"] is True
    assert summary["volume_split_relative_difference"] < 1e-12
    assert summary["layers_per_side"] == 2
