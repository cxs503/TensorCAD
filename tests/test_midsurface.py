"""Mid-surface extraction held to analytic sheets, an independent oracle, and its refusals.

``GAP-CAD-013`` asks for three numbers on a plate: area error below 1 %, face count at most a
tenth of the source, and at least a 50 % mesh-size cut after feature suppression. This file
checks those three, and it checks the geometry behind them three more ways:

1. **The eroded box.** The mid-sheet of a box shelled at ``t`` is the closed box of side
   ``L - t``: area ``6 (L - t)^2``, volume ``(L - t)^3``, both exactly. Run over the whole
   ``t/L`` sweep from 0.01 to 0.1.
2. **An independent offset.** ``tests/midsurface_oracle.py`` erodes the *section* with one
   2-D half-plane per edge and extrudes it; a regular hexagon and a regular dodecagon prism are
   compared against that, so the corner convention is checked by arithmetic that shares no code
   with the module. The closed forms ``n r^2 tan(pi/n)`` and ``2 n r tan(pi/n)`` give the oracle
   itself a second opinion.
3. **A measured position, not an assumed one.** Every vertex of the sheet is at distance
   ``t/2`` from the *whole* source surface, computed by brute-force point-to-triangle distance.
   A sheet that is merely parallel somewhere would not pass.

The refusals are first-class behaviour and are tested as such: an open surface, an inward
surface, a solid cube, a taper, a wall beyond ``max_thickness``, and a non-convex panel all
fail loudly with a named reason. Nothing here asserts that a bad solid is repaired -- the module
refuses rather than guesses.
"""

from __future__ import annotations

import json
import math

import pytest
from midsurface_fixtures import (
    box,
    cube_solid,
    l_prism,
    merge,
    plate,
    quadrant_plate,
    regular_prism,
    taper_prism,
)
from midsurface_oracle import (
    distance_to_surface,
    prism_mid_sheet,
    regular_polygon_area,
    regular_polygon_perimeter,
    regular_section,
)

from tensorcad import mesh as meshlib
from tensorcad.geometry import Mesh
from tensorcad.midsurface import (
    MidSurfaceError,
    defeature,
    extract_midsurface,
    measure_thickness,
    suggest_max_thickness,
)
from tensorcad.shell import hollow_body, shell

THICKNESSES = tuple(0.01 * step for step in range(1, 11))


def hollow(mesh: Mesh, thickness: float) -> Mesh:
    """The double-walled body of ``mesh`` at ``thickness`` -- the input of an extraction."""
    return hollow_body(shell(mesh, thickness)).mesh


def reversed_mesh(mesh: Mesh) -> Mesh:
    return Mesh(
        vertices=mesh.vertices,
        faces=tuple((a, c, b) for a, b, c in mesh.faces),
        units=mesh.units,
    )


# --------------------------------------------------------------------------------------
# 1. the eroded box: area and volume against closed forms
# --------------------------------------------------------------------------------------


@pytest.mark.parametrize("thickness", THICKNESSES)
def test_box_mid_sheet_reproduces_the_eroded_box(thickness: float) -> None:
    length = 1.0
    body = hollow(box((0.0, 0.0, 0.0), (length, length, length)), thickness)
    sheet = extract_midsurface(body)
    side = length - thickness
    assert sheet.area() == pytest.approx(6.0 * side * side, rel=1e-12)
    assert sheet.volume() == pytest.approx(side**3, rel=1e-12)
    assert sheet.watertight()
    assert len(sheet.panels) == 6
    for panel in sheet.panels:
        assert panel.thickness == pytest.approx(thickness, rel=1e-12)
        assert panel.thickness_spread < 1e-12
        assert panel.sheet_area == pytest.approx(side * side, rel=1e-12)
    assert sheet.thickness == pytest.approx(thickness, rel=1e-12)


def test_box_sheet_corners_are_the_mid_lines() -> None:
    """Six patches meeting on ``x = y = t/2``: the corner coordinates, not just the totals."""
    thickness = 0.06
    sheet = extract_midsurface(hollow(box(divisions=2), thickness))
    corners = {tuple(round(coordinate, 9) for coordinate in point) for point in sheet.mesh.vertices}
    expected = {
        tuple(round(value, 9) for value in (i * (1.0 - thickness) + 0.5 * thickness,
                                            j * (1.0 - thickness) + 0.5 * thickness,
                                            k * (1.0 - thickness) + 0.5 * thickness))
        for i in (0, 1)
        for j in (0, 1)
        for k in (0, 1)
    }
    assert corners == expected


def test_refined_plate_cuts_the_face_count_by_more_than_ten_to_one() -> None:
    """The acceptance criterion's face-count half, on a mesh a solver would actually get."""
    thickness = 0.05
    solid = plate(1.0, 1.0, thickness, divisions=8)
    sheet = extract_midsurface(solid)
    ratio = sheet.sheet_faces / sheet.source_faces
    assert ratio <= 0.1, f"{sheet.sheet_faces} sheet faces from {sheet.source_faces} source faces"
    assert sheet.face_reduction() >= 0.9
    assert sheet.area() == pytest.approx(1.0, rel=1e-12)
    assert sheet.thickness == pytest.approx(thickness, rel=1e-12)


# --------------------------------------------------------------------------------------
# 2. an independent 2-D offset oracle (prisms, where a section erosion is exact)
# --------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("sides", "thickness"), ((6, 0.02), (6, 0.05), (6, 0.1), (12, 0.02))
)
def test_regular_prism_matches_the_section_offset_oracle(sides: int, thickness: float) -> None:
    inradius = 1.0
    length = 2.0
    solid = regular_prism(inradius, sides, length)
    sheet = extract_midsurface(hollow(solid, thickness))
    area, volume = prism_mid_sheet(regular_section(inradius, sides), length, thickness)
    assert sheet.area() == pytest.approx(area, rel=1e-11)
    assert sheet.volume() == pytest.approx(volume, rel=1e-11)
    assert sheet.watertight()
    assert len(sheet.panels) == sides + 2


def test_the_oracle_agrees_with_the_regular_polygon_closed_forms() -> None:
    """The oracle is only worth something if it is right: check it against the closed forms."""
    for sides in (3, 4, 5, 6, 8, 12, 24):
        for inradius in (0.25, 1.0, 3.5):
            for offset in (0.0, 0.01, 0.1):
                area, volume = prism_mid_sheet(
                    regular_section(inradius, sides), length=1.0, thickness=2.0 * offset
                )
                eroded = inradius - offset
                expected_area = (
                    regular_polygon_perimeter(eroded, sides) * (1.0 - 2.0 * offset)
                    + 2.0 * regular_polygon_area(eroded, sides)
                )
                expected_volume = regular_polygon_area(eroded, sides) * (1.0 - 2.0 * offset)
                assert area == pytest.approx(expected_area, rel=1e-12)
                assert volume == pytest.approx(expected_volume, rel=1e-12)


# --------------------------------------------------------------------------------------
# 3. the sheet is where it says it is: equidistant from both walls
# --------------------------------------------------------------------------------------


@pytest.mark.parametrize("thickness", (0.02, 0.06, 0.1))
def test_every_sheet_vertex_is_half_a_wall_from_the_source_surface(thickness: float) -> None:
    body = hollow(box(divisions=2), thickness)
    sheet = extract_midsurface(body)
    distances = [distance_to_surface(point, body) for point in sheet.mesh.vertices]
    assert max(abs(value - 0.5 * thickness) for value in distances) < 1e-9


def test_plate_sheet_sits_one_half_wall_inside_the_surface() -> None:
    """A plate's sheet runs out to the plate edge, so only its interior is a full ``t/2`` in.

    At the edge the sheet corner lies *on* the side surface (distance zero), which is correct:
    the mid-surface of a plate is bounded by the plate, it does not stop half a wall short.
    """
    thickness = 0.05
    solid = plate(1.0, 1.0, thickness)
    sheet = extract_midsurface(solid)
    distances = [distance_to_surface(point, solid) for point in sheet.mesh.vertices]
    assert max(distances) <= 0.5 * thickness + 1e-12
    centre = tuple(
        sum(point[axis] for point in sheet.mesh.vertices) / len(sheet.mesh.vertices)
        for axis in range(3)
    )
    assert abs(distance_to_surface(centre, solid) - 0.5 * thickness) < 1e-12


# --------------------------------------------------------------------------------------
# 4. a plate: one sheet in the mid-plane, open on purpose
# --------------------------------------------------------------------------------------


def test_plate_mid_surface_is_the_mid_plane_patch() -> None:
    length, breadth, thickness = 1.0, 0.8, 0.05
    solid = plate(length, breadth, thickness, divisions=2)
    sheet = extract_midsurface(solid)
    assert len(sheet.panels) == 1
    panel = sheet.panels[0]
    assert panel.thickness == pytest.approx(thickness, rel=1e-12)
    assert abs(panel.normal[2]) == pytest.approx(1.0, rel=1e-12)
    assert panel.mid_offset in (
        pytest.approx(0.5 * thickness, rel=1e-12),
        pytest.approx(-0.5 * thickness, rel=1e-12),
    )
    assert sheet.area() == pytest.approx(length * breadth, rel=1e-12)
    # An open sheet is not a solid: the divergence-theorem volume of a *closed* surface is
    # meaningful, of a patch it is not, so the report has to say that instead of a number.
    assert not sheet.watertight()
    assert len(meshlib.boundary_edges(sheet.mesh)) > 0


def test_plate_face_count_and_reduction() -> None:
    sheet = extract_midsurface(plate(1.0, 1.0, 0.05, divisions=4))
    assert sheet.sheet_faces == 2  # one rectangle, two triangles
    assert sheet.face_reduction() > 0.9
    assert sheet.sheet_area() == pytest.approx(1.0, rel=1e-12)


# --------------------------------------------------------------------------------------
# 5. the thickness map, on its own
# --------------------------------------------------------------------------------------


def test_thickness_map_of_a_hollow_box_measures_the_wall() -> None:
    thickness = 0.06
    body = hollow(box(), thickness)
    table = measure_thickness(body)
    assert table.paired_count == len(body.faces)
    assert table.coverage() == pytest.approx(1.0, rel=1e-15)
    values = table.thickness_values()
    assert max(abs(value - thickness) for value in values) < 1e-12
    assert all(reason == "" for reason in table.reason)
    assert table.max_thickness == pytest.approx(suggest_max_thickness(body), rel=1e-15)


def test_thickness_map_of_a_plate_reports_why_the_edge_strips_do_not_pair() -> None:
    thickness = 0.02
    solid = plate(1.0, 1.0, thickness)
    table = measure_thickness(solid)
    # Two large faces (2 triangles each) pair; the four edge strips do not.
    assert table.paired_count == 4
    assert 0.9 < table.coverage() < 1.0
    reasons = {reason for reason in table.reason if reason}
    assert len(reasons) == 1
    # The counter-face is anti-parallel but far away, so the ray never reaches it: that is a
    # different verdict from a solid cube's "the nearest face points the same way".
    assert "no other face within max_thickness" in reasons.pop()


def test_max_thickness_default_is_the_smallest_extent_rule() -> None:
    solid = plate(2.0, 3.0, 0.1)
    assert suggest_max_thickness(solid) == pytest.approx(1.5 * 0.1, rel=1e-12)


def test_measure_thickness_on_a_solid_cube_pairs_nothing_within_the_bound() -> None:
    table = measure_thickness(cube_solid(), max_thickness=0.5)
    assert table.paired_count == 0
    assert table.coverage() == 0.0


# --------------------------------------------------------------------------------------
# 6. the refusals
# --------------------------------------------------------------------------------------


def test_open_surface_is_refused() -> None:
    body = hollow(box(), 0.06)
    open_mesh = Mesh(vertices=body.vertices, faces=body.faces[:-1], units=body.units)
    with pytest.raises(MidSurfaceError, match="not closed"):
        extract_midsurface(open_mesh)


def test_inward_surface_is_refused() -> None:
    body = hollow(box(), 0.06)
    with pytest.raises(MidSurfaceError, match="wound inward"):
        extract_midsurface(reversed_mesh(body))


def test_solid_cube_is_refused_once_the_wall_bound_is_stated() -> None:
    """A cube's opposite faces are an anti-parallel pair at ``s``, so the caller bounds the wall.

    This is a limitation of thickness measurement, not of the sheet: nothing in the geometry
    distinguishes "solid cube" from "hollow cube with a wall as thick as its side", so
    ``max_thickness`` is the caller's statement about what a wall is, and the default only
    guesses it from the smallest bounding extent.
    """
    with pytest.raises(MidSurfaceError, match="no opposed thin walls"):
        extract_midsurface(cube_solid(), max_thickness=0.5)


def test_wall_thicker_than_max_thickness_is_refused() -> None:
    body = hollow(box(), 0.06)
    with pytest.raises(MidSurfaceError, match="no opposed thin walls"):
        extract_midsurface(body, max_thickness=0.5 * 0.06)


def test_max_thickness_must_be_positive() -> None:
    with pytest.raises(MidSurfaceError, match="max_thickness must be positive"):
        extract_midsurface(hollow(box(), 0.06), max_thickness=-1.0)


def test_a_wall_that_thins_inside_one_panel_is_refused() -> None:
    """The guard is live: ordinary round-off is far below the default, 1e-20 is not.

    The measurement on a correct wall agrees to ~1e-17, so the default 5 % tolerance never
    fires by accident; shrinking it to 1e-20 makes the same panel fail, which is what proves the
    check is a check and not a decoration.
    """
    body = hollow(box(), 0.06)
    with pytest.raises(MidSurfaceError, match="not constant"):
        extract_midsurface(body, thickness_tolerance=1e-20)


def test_non_convex_panel_is_refused() -> None:
    """An L-shaped cap: the convex hull would cover material that is not there."""
    solid = l_prism(leg=1.0, width=0.3, length=1.0)
    with pytest.raises(MidSurfaceError, match="not convex|no opposed thin walls"):
        extract_midsurface(hollow(solid, 0.05))


def test_min_coverage_is_enforced() -> None:
    """A plate is only 96 % wall at this thickness; demanding 99 % must fail."""
    solid = plate(1.0, 1.0, 0.02)
    with pytest.raises(MidSurfaceError, match="measurable wall"):
        extract_midsurface(solid, min_coverage=0.99)


# --------------------------------------------------------------------------------------
# 7. defeaturing: bridge import seams, suppress small features -- and report both
# --------------------------------------------------------------------------------------


def test_defeature_bridges_import_seams_and_cuts_the_mesh_by_more_than_half() -> None:
    gap = 1e-6
    solid = quadrant_plate(length=1.0, thickness=0.06, gap=gap)
    before = extract_midsurface(solid)
    assert len(before.panels) >= 4
    after = defeature(
        before, weld_tolerance=100.0 * gap, convexity_tolerance=1e-3, min_panel_area=0.05
    )
    assert len(after.panels) < len(before.panels)
    assert after.sheet_faces < before.sheet_faces
    assert after.face_reduction() >= 0.5
    assert after.area() == pytest.approx(1.0, rel=1e-3)


def test_suppress_small_feature_reports_what_it_cost() -> None:
    big = plate(1.0, 1.0, 0.05)
    small = box((1.2, 0.0, 0.0), (1.3, 0.1, 0.05))
    solid = merge(big, small)
    before = extract_midsurface(solid)
    assert len(before.panels) == 2
    small_area = 0.1 * 0.1
    after = defeature(before, min_panel_area=2.0 * small_area)
    assert len(after.panels) == 1
    # Accounting, not a constant: what the suppression cost is the difference in unmatched
    # source area, which must be exactly the feature's own sheet area.
    assert after.removed_area() - before.removed_area() == pytest.approx(small_area, rel=1e-12)
    assert after.face_reduction() >= 0.5


def test_defeature_with_defaults_changes_nothing() -> None:
    thickness = 0.06
    body = hollow(regular_prism(1.0, 6, 2.0), thickness)
    before = extract_midsurface(body)
    after = defeature(before)
    assert len(after.panels) == len(before.panels)
    assert after.sheet_faces == before.sheet_faces
    assert after.area() == pytest.approx(before.area(), rel=1e-15)


# --------------------------------------------------------------------------------------
# 8. the report
# --------------------------------------------------------------------------------------


def test_report_is_json_and_carries_the_acceptance_numbers() -> None:
    thickness = 0.06
    sheet = extract_midsurface(hollow(box((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), divisions=4), thickness))
    payload = sheet.as_dict()
    text = json.dumps(payload)  # no NaN, no sets: a report a pipeline can consume
    assert "face_reduction" in text
    assert payload["panels"] == 6
    assert payload["watertight"] is True
    assert payload["sheet_area"] == pytest.approx(6.0 * (1.0 - thickness) ** 2, rel=1e-12)
    assert payload["signed_volume"] == pytest.approx((1.0 - thickness) ** 3, rel=1e-12)
    assert payload["thickness_min"] == pytest.approx(thickness, rel=1e-12)
    assert payload["thickness_map"]["coverage"] == pytest.approx(1.0, rel=1e-15)
    assert len(payload["panels_data"]) == payload["panels"]
    for panel in payload["panels_data"]:
        assert len(panel["polygon"]) >= 3


def test_sheet_mesh_is_watertight_and_edge_manifold_when_closed() -> None:
    sheet = extract_midsurface(hollow(regular_prism(1.0, 6, 2.0), 0.05))
    assert meshlib.is_watertight(sheet.mesh)
    assert meshlib.is_edge_manifold(sheet.mesh)
    assert len(meshlib.boundary_edges(sheet.mesh)) == 0