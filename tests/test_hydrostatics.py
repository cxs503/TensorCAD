"""Tests for mesh-integrated hydrostatics."""
from __future__ import annotations

import itertools
import math

import pytest

from tensorcad import HullParameters, compute_hydrostatics, generate_hull
from tensorcad import mesh as meshlib
from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.hydrostatics import max_section_area, section_area_at, stations


def box(minimum: tuple[float, float, float], maximum: tuple[float, float, float]) -> Mesh:
    """Closed, outward-oriented axis-aligned box built from 12 raw triangles."""
    x0, y0, z0 = minimum
    x1, y1, z1 = maximum
    builder = MeshBuilder()
    corners = [
        (x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
        (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1),
    ]
    ids = [builder.add_vertex(point) for point in corners]
    for face in (
        (0, 3, 2), (0, 2, 1), (4, 5, 6), (4, 6, 7), (0, 1, 5), (0, 5, 4),
        (2, 3, 7), (2, 7, 6), (1, 2, 6), (1, 6, 5), (0, 4, 7), (0, 7, 3),
    ):
        builder.add_triangle(*[ids[index] for index in face])
    return meshlib.orient_outward(builder.build()).mesh


@pytest.fixture()
def unit_box() -> Mesh:
    return box((0.0, -0.5, 0.0), (1.0, 0.5, 1.0))


def test_unit_box_integrates_exactly(unit_box: Mesh) -> None:
    hydro = compute_hydrostatics(unit_box, waterline_z=1.0)
    assert hydro.displacement_volume == pytest.approx(1.0)
    assert hydro.wetted_surface_area == pytest.approx(5.0), "bottom plus four sides, no deck"
    assert hydro.waterplane_area == pytest.approx(1.0), "from the vertical flux of the wetted faces"
    assert hydro.longitudinal_centre_buoyancy == pytest.approx(0.5)
    assert hydro.longitudinal_centre_flotation == pytest.approx(0.5)
    assert hydro.vertical_centre_buoyancy == pytest.approx(0.5)
    assert hydro.block_coefficient == pytest.approx(1.0)
    assert hydro.prismatic_coefficient == pytest.approx(1.0)
    assert hydro.maximum_section_area == pytest.approx(1.0)
    assert hydro.displacement_mass == pytest.approx(1025.0)
    assert hydro.basis == "closed surface integrals"


def test_water_density_is_an_explicit_input(unit_box: Mesh) -> None:
    fresh = compute_hydrostatics(unit_box, 1.0, water_density=1000.0)
    assert fresh.displacement_mass == pytest.approx(1000.0)
    with pytest.raises(ValueError, match="water_density"):
        compute_hydrostatics(unit_box, 1.0, water_density=0.0)
    with pytest.raises(ValueError, match="water_density"):
        compute_hydrostatics(unit_box, 1.0, water_density=float("nan"))


def test_open_surface_is_flagged_rather_than_silently_wrong(unit_box: Mesh) -> None:
    opened = Mesh(
        vertices=unit_box.vertices,
        faces=tuple(face for face in unit_box.faces if meshlib.face_centroid(unit_box, face)[2] < 0.99),
    )
    assert not meshlib.is_watertight(opened)
    hydro = compute_hydrostatics(opened, waterline_z=1.0)
    assert "open surface" in hydro.basis
    # The flux-based waterplane integral survives an open deck; the volume does not.
    assert hydro.waterplane_area == pytest.approx(1.0)
    # What is missing is exactly the flux of the removed cap: 1 - 2*(1/6) = 2/3.
    assert hydro.displacement_volume == pytest.approx(2.0 / 3.0)
    assert hydro.block_coefficient == pytest.approx(2.0 / 3.0)


def test_zero_volume_short_circuits_to_nan_instead_of_dividing_by_zero() -> None:
    flat = box((0.0, 0.0, 0.0), (1.0, 1.0, 0.0))
    hydro = compute_hydrostatics(flat, waterline_z=0.0)
    assert hydro.displacement_volume == 0.0
    assert math.isnan(hydro.longitudinal_centre_buoyancy)
    assert math.isnan(hydro.vertical_centre_buoyancy)
    assert math.isnan(hydro.block_coefficient)
    assert math.isnan(hydro.prismatic_coefficient)
    assert hydro.wetted_surface_area == 0.0
    assert hydro.waterplane_area == 0.0
    payload = hydro.to_dict()
    # NaN survives into the API rather than being silently zeroed; export turns it to null.
    nan = payload["longitudinal_centre_buoyancy_m"]
    assert nan != nan


def test_section_area_functions(unit_box: Mesh) -> None:
    assert stations(unit_box) == (0.0, 1.0)
    assert section_area_at(unit_box, 0.0) == pytest.approx(1.0)
    assert max_section_area(unit_box) == pytest.approx(1.0)
    assert math.isnan(section_area_at(unit_box, 99.0))


def test_section_area_documents_the_symmetry_assumption() -> None:
    """A section is treated as symmetric about y=0, so an offset box doubles."""
    offset = box((0.0, 0.0, 0.0), (1.0, 1.0, 1.0))
    # Real area is 1.0; 2 * integral of max|y| over z counts the mirror image.
    assert section_area_at(offset, 0.0) == pytest.approx(2.0)


def test_hull_hydrostatics_are_reported_in_the_metadata() -> None:
    model = generate_hull(HullParameters(nx=11, nz=5))
    data = model.metadata()["hydrostatics"]
    assert data["basis"] == "closed surface integrals"
    assert data["water_density_kg_per_m3"] == 1025.0
    assert data["displacement_volume_m3"] == pytest.approx(model.hydrostatics.displacement_volume)
    assert all(math.isfinite(value) for value in data.values() if isinstance(value, float))
    assert data["block_coefficient"] < 0.5 < data["prismatic_coefficient"]


def brute_force_section_area(mesh: Mesh, x: float, tolerance: float = 1e-9) -> float:
    """The pre-optimisation definition, kept here as the reference result."""
    profile: dict[float, float] = {}
    for px, py, pz in mesh.vertices:
        if abs(px - x) <= tolerance:
            previous = profile.get(pz)
            breadth = abs(py)
            if previous is None or breadth > previous:
                profile[pz] = breadth
    levels = sorted(profile)
    if len(levels) < 2:
        return float("nan")
    area = 0.0
    for lower, upper in itertools.pairwise(levels):
        area += 0.5 * (profile[lower] + profile[upper]) * (upper - lower)
    return 2.0 * area


def test_section_areas_match_the_brute_force_definition_on_every_station() -> None:
    model = generate_hull(HullParameters(nx=41, nz=13))
    for station in stations(model.mesh):
        assert section_area_at(model.mesh, station) == pytest.approx(
            brute_force_section_area(model.mesh, station)
        )
    assert max_section_area(model.mesh) == pytest.approx(
        max(brute_force_section_area(model.mesh, station) for station in stations(model.mesh))
    )
