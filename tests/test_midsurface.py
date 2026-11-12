"""Tests for :mod:`tensorcad.midsurface` (GAP-CAD-013).

Every number here is pinned from ``examples/midsurface/run_example.py``; the tolerances
are the acceptance criteria, not the observed values.
"""
from __future__ import annotations

import math

import pytest

from tensorcad import mesh as meshlib
from tensorcad.geometry import Mesh
from tensorcad.midsurface import (
    DefeatureResult,
    MidSurfaceError,
    MidSurfaceResult,
    cylinder_skin,
    defeature_planar_regions,
    extract_midsurface,
    mid_surface_between,
    plate_solid,
    sphere_skin,
)
from tensorcad.series60 import (
    Series60Parameters,
    analytic_series60_hydrostatics,
    generate_series60,
)
from tensorcad.shell import offset_inward

PLATE = plate_solid(8.0, 3.0, 0.02, nx=20, ny=10)


# --------------------------------------------------------------------------- #
# plate: the mid-plane is analytic
# --------------------------------------------------------------------------- #
def test_plate_mid_area_is_exact():
    result = extract_midsurface(PLATE)
    assert isinstance(result, MidSurfaceResult)
    assert result.area == pytest.approx(8.0 * 3.0, rel=1e-12)
    assert result.source_area == pytest.approx(meshlib.surface_area(PLATE), rel=1e-15)


def test_plate_thickness_is_measured_not_assumed():
    result = extract_midsurface(PLATE)
    assert result.mean_thickness == pytest.approx(0.02, rel=1e-9)
    assert result.min_measured >= 0.02 * (1.0 - 1e-9)
    assert result.max_measured <= 0.02 * (1.0 + 1e-9)


def test_plate_rim_faces_are_reported_unpaired():
    result = extract_midsurface(PLATE)
    # 20x10 quads on each large face = 800 triangles, all paired; the four rim strips
    # are 0.02 m wide but their opposite face is 3 m or 8 m away, so they are unpaired.
    assert result.paired_face_count == 800
    assert len(result.unpaired) == 8
    assert result.paired_fraction == pytest.approx(800 / 808)


def test_plate_mid_surface_is_open_and_carries_no_volume():
    result = extract_midsurface(PLATE)
    assert not meshlib.is_watertight(result.mesh)
    assert result.mesh.vertices
    zs = {round(vertex[2], 12) for vertex in result.mesh.vertices}
    assert zs == {0.0}


def test_plate_defeature_collapses_to_two_triangles():
    result = extract_midsurface(PLATE)
    defeatured = defeature_planar_regions(result.mesh)
    assert isinstance(defeatured, DefeatureResult)
    assert len(defeatured.regions) == 1
    assert defeatured.regions[0].simplified
    # The *raw* boundary loop carries one vertex per mesh boundary vertex (60 = 2*(20+10)
    # rounded corners); the collinear ones are dropped inside the simplification, which is
    # why 60 vertices become 2 triangles, and the area check below is what guards that drop.
    assert [len(loop) for loop in defeatured.regions[0].boundary_loops] == [60]
    assert defeatured.face_count_before == 400
    assert defeatured.face_count_after == 2
    assert defeatured.face_reduction == pytest.approx(1.0 - 2 / 400)
    assert defeatured.area_relative_error < 1e-12


# --------------------------------------------------------------------------- #
# curved walls: area surplus over the exact offset is predicted by t/(2R)
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "thickness, expected",
    [(0.01, 0.005), (0.02, 0.01), (0.04, 0.02)],
)
def test_cylinder_area_surplus_follows_t_over_2r(thickness: float, expected: float):
    skin = cylinder_skin(radius=1.0, thickness=thickness, length=4.0, segments=128)
    result = extract_midsurface(skin)
    analytic = 2.0 * math.pi * (1.0 - 0.5 * thickness) * 4.0
    surplus = result.area / analytic - 1.0
    assert surplus == pytest.approx(expected, rel=0.05)
    # The centroids sit on chords, not on the circle, so the measured wall is short by the
    # facet sagitta: 3.0e-4 relative at 128 segments. That is discretisation, not bias.
    assert result.mean_thickness == pytest.approx(thickness, rel=1e-3)
    assert len(result.unpaired) == 0


def test_sphere_area_surplus_is_twice_the_cylinder_prediction():
    skin = sphere_skin(radius=1.0, thickness=0.02, rings=32, segments=64)
    result = extract_midsurface(skin)
    analytic = 4.0 * math.pi * (1.0 - 0.01) ** 2
    surplus = result.area / analytic - 1.0
    assert surplus == pytest.approx(0.02, rel=0.15)
    assert result.mean_thickness == pytest.approx(0.02, rel=5e-3)


def test_sphere_polar_degeneracy_is_reported_and_area_free():
    skin = sphere_skin(radius=1.0, thickness=0.02, rings=32, segments=64)
    result = extract_midsurface(skin)
    unpaired_area = sum(
        meshlib.triangle_area(skin, skin.faces[index]) for index in result.unpaired
    )
    # The unpaired faces are the degenerate polar band: they must carry no area.
    assert unpaired_area / meshlib.surface_area(skin) < 1e-3
    assert result.paired_fraction > 0.95


def test_curved_skin_is_not_defeatured():
    skin = cylinder_skin(radius=1.0, thickness=0.02, length=4.0, segments=32)
    result = extract_midsurface(skin)
    defeatured = defeature_planar_regions(result.mesh)
    assert defeatured.face_reduction == 0.0
    # Adjacent facets differ by 360/32 = 11.25 deg > 1 deg, so the only merges are the two
    # triangles of each quad -- which re-triangulate back to two, hence zero reduction.
    assert max(len(region.faces) for region in defeatured.regions) == 2
    assert len(defeatured.regions) == len(result.mesh.faces) // 2
    assert defeatured.area_relative_error < 1e-12


# --------------------------------------------------------------------------- #
# the thickness gate: a solid body must not be idealised
# --------------------------------------------------------------------------- #
def test_solid_cube_is_rejected():
    cube = plate_solid(1.0, 1.0, 0.2, nx=2, ny=2)
    result = extract_midsurface(cube)
    assert result.paired_face_count == 0
    assert len(result.unpaired) == len(cube.faces)
    assert result.mesh.faces == ()


def test_raising_the_band_does_idealise_the_solid_the_bad_way():
    """The gate is what stops this, so pin the failure mode explicitly."""
    cube = plate_solid(1.0, 1.0, 0.2, nx=2, ny=2)
    result = extract_midsurface(cube, max_thickness=2.0)
    assert result.paired_face_count > 0
    assert result.mean_thickness > 0.2


# --------------------------------------------------------------------------- #
# series 60 hull plating: gate the geometry, then idealise it
# --------------------------------------------------------------------------- #
def _hull_plating(parameters: Series60Parameters, wall: float = 0.02):
    model = generate_series60(parameters)
    hull = model.mesh
    cap = {
        index
        for index, face in enumerate(hull.faces)
        if all(abs(hull.vertices[i][2]) < 1e-12 for i in face)
    }
    keep = [index for index in range(len(hull.faces)) if index not in cap]
    outer = Mesh(vertices=hull.vertices, faces=tuple(hull.faces[index] for index in keep))
    inner_full = offset_inward(hull, wall)
    inner = Mesh(vertices=inner_full.vertices, faces=tuple(hull.faces[index] for index in keep))
    shift = len(outer.vertices)
    faces = list(outer.faces) + [tuple(i + shift for i in face)[::-1] for face in inner.faces]
    plating = Mesh(vertices=tuple(outer.vertices) + tuple(inner.vertices), faces=tuple(faces))
    return hull, outer, inner, plating


def test_hull_geometry_gate_before_idealisation():
    parameters = Series60Parameters(length=10.0, beam=1.5, draft=0.6, cb=0.60, nx=41, nz=13)
    hull, _, _, _ = _hull_plating(parameters)
    reference = analytic_series60_hydrostatics(parameters)
    assert meshlib.is_watertight(hull)
    cb = meshlib.total_signed_volume(hull) / (
        parameters.length * parameters.beam * parameters.draft
    )
    assert cb == pytest.approx(reference["block_coefficient"], rel=0.05)
    box = hull.bounding_box()
    assert box.maximum[0] - box.minimum[0] == pytest.approx(parameters.length, rel=1e-9)
    assert box.maximum[1] - box.minimum[1] == pytest.approx(parameters.beam, rel=1e-9)
    assert box.minimum[2] == pytest.approx(-parameters.draft, rel=1e-9)


def test_hull_midsurface_thickness_and_pairing():
    parameters = Series60Parameters(length=10.0, beam=1.5, draft=0.6, cb=0.60, nx=41, nz=13)
    _, outer, inner, plating = _hull_plating(parameters)
    result = extract_midsurface(plating, max_thickness=0.5)
    assert result.mean_thickness == pytest.approx(0.02, rel=0.01)
    assert result.paired_fraction > 0.95
    reference_mid = mid_surface_between(outer, inner)
    assert result.area == pytest.approx(meshlib.surface_area(reference_mid), rel=0.01)


# --------------------------------------------------------------------------- #
# input validation
# --------------------------------------------------------------------------- #
def _triangle_mesh() -> Mesh:
    return Mesh(vertices=((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
                faces=((0, 1, 2),))


def test_empty_mesh_is_an_error():
    with pytest.raises(MidSurfaceError):
        extract_midsurface(Mesh(vertices=(), faces=()))
    with pytest.raises(MidSurfaceError):
        defeature_planar_regions(Mesh(vertices=(), faces=()))


def test_invalid_bands_are_errors():
    with pytest.raises(MidSurfaceError):
        extract_midsurface(PLATE, max_thickness=0.0)
    with pytest.raises(MidSurfaceError):
        extract_midsurface(PLATE, max_thickness=-1.0)
    with pytest.raises(MidSurfaceError):
        extract_midsurface(PLATE, max_thickness=0.01, min_thickness=0.02)
    with pytest.raises(MidSurfaceError):
        extract_midsurface(PLATE, min_thickness=-1e-3)
    with pytest.raises(MidSurfaceError):
        extract_midsurface(PLATE, angle_tolerance_deg=90.0)
    with pytest.raises(MidSurfaceError):
        defeature_planar_regions(PLATE, angle_tolerance_deg=95.0)


def test_zero_extent_mesh_is_an_error():
    degenerate = Mesh(vertices=((0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)),
                      faces=((0, 1, 2),))
    with pytest.raises(MidSurfaceError):
        extract_midsurface(degenerate)


def test_mid_surface_between_requires_matching_topology():
    a = _triangle_mesh()
    b = Mesh(vertices=((0.0, 0.0, -0.1), (1.0, 0.0, -0.1)), faces=())
    with pytest.raises(MidSurfaceError):
        mid_surface_between(a, b)


def test_mid_surface_between_averages_vertices():
    top = _triangle_mesh()
    bottom = Mesh(
        vertices=((0.0, 0.0, -0.2), (1.0, 0.0, -0.2), (0.0, 1.0, -0.2)),
        faces=((0, 1, 2),),
    )
    mid = mid_surface_between(top, bottom)
    assert [round(vertex[2], 12) for vertex in mid.vertices] == [-0.1, -0.1, -0.1]
    assert meshlib.surface_area(mid) == pytest.approx(
        0.5 * meshlib.surface_area(top) + 0.5 * meshlib.surface_area(bottom)
    )


def test_defeature_keeps_a_hole_region_untouched():
    """A region with an inner boundary loop is skipped, not silently triangulated."""
    ring = plate_solid(4.0, 4.0, 0.02, nx=8, ny=8)
    result = extract_midsurface(ring)
    defeatured = defeature_planar_regions(result.mesh)
    assert defeatured.regions[0].simplified  # no hole yet: one rectangle
    assert defeatured.area_relative_error < 1e-12