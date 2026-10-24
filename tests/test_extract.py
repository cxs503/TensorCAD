"""Tests for turning an implicit field into a closed, solver-ready mesh."""
from __future__ import annotations

import math

import numpy as np
import pytest

from tensorcad import mesh as meshlib
from tensorcad.neural.extract import ExtractionError, extract_surface, sample_grid, summarise

pytest.importorskip("numpy")


def sphere_field(points: np.ndarray) -> np.ndarray:
    return np.linalg.norm(points, axis=1) - 0.5


def ellipsoid_field(points: np.ndarray) -> np.ndarray:
    return np.sqrt((points[:, 0] / 0.6) ** 2 + (points[:, 1] / 0.35) ** 2 + (points[:, 2] / 0.45) ** 2) - 1.0


@pytest.mark.parametrize("resolution", [24, 48])
def test_sphere_extraction_is_closed_and_oriented(resolution: int) -> None:
    mesh = extract_surface(sphere_field, (-0.6,) * 3, (0.6,) * 3, resolution)
    report = summarise(mesh, resolution)
    exact = 4.0 / 3.0 * math.pi * 0.5**3
    assert report.watertight
    assert report.edge_manifold
    assert report.euler_characteristic == 2
    assert report.signed_volume > 0.0
    assert abs(report.signed_volume / exact - 1.0) < 0.01, "volume must match the analytic sphere"
    assert report.smallest_triangle_area > 0.0

    vertices = np.array(mesh.vertices)
    faces = np.array(mesh.faces)
    normals = np.cross(
        vertices[faces[:, 1]] - vertices[faces[:, 0]], vertices[faces[:, 2]] - vertices[faces[:, 0]]
    )
    centroids = vertices[faces].mean(axis=1)
    outward = np.einsum("ij,ij->i", normals, centroids)
    assert np.all(outward > 0.0), "every triangle must face away from the enclosed volume"


def test_refining_the_grid_reduces_the_volume_error() -> None:
    exact = 4.0 / 3.0 * math.pi * 0.5**3
    coarse = summarise(extract_surface(sphere_field, (-0.6,) * 3, (0.6,) * 3, 24), 24)
    fine = summarise(extract_surface(sphere_field, (-0.6,) * 3, (0.6,) * 3, 48), 48)
    assert abs(fine.signed_volume - exact) < abs(coarse.signed_volume - exact)


def test_ellipsoid_extraction_matches_its_analytic_volume() -> None:
    bounds = ((-0.7, -0.45, -0.55), (0.7, 0.45, 0.55))
    mesh = extract_surface(ellipsoid_field, bounds[0], bounds[1], 48)
    report = summarise(mesh, 48)
    exact = 4.0 / 3.0 * math.pi * 0.6 * 0.35 * 0.45
    assert report.watertight and report.euler_characteristic == 2
    assert abs(report.signed_volume / exact - 1.0) < 0.01


def test_units_and_frame_are_carried_into_the_mesh() -> None:
    mesh = extract_surface(
        sphere_field, (-0.6,) * 3, (0.6,) * 3, 24, units="m", coordinate_frame="hull_local"
    )
    assert mesh.units == "m"
    assert mesh.coordinate_frame == "hull_local"


def test_summarise_rejects_a_unit_mismatch() -> None:
    mesh = extract_surface(sphere_field, (-0.6,) * 3, (0.6,) * 3, 16)
    with pytest.raises(ExtractionError, match="units"):
        summarise(mesh, 16, expected_units="mm")


def test_sampling_rejects_bad_inputs() -> None:
    with pytest.raises(ExtractionError, match="positive extent"):
        sample_grid(sphere_field, (0.0, 0.0, 0.0), (1.0, 0.0, 1.0), 16)
    with pytest.raises(ExtractionError, match="at least 2 cells"):
        sample_grid(sphere_field, (-1.0,) * 3, (1.0,) * 3, 1)
    with pytest.raises(ExtractionError, match="at least 1 cell"):
        sample_grid(sphere_field, (-1.0,) * 3, (1.0,) * 3, 8, margin=0)


def test_a_field_without_a_surface_is_reported_precisely() -> None:
    def always_outside(points: np.ndarray) -> np.ndarray:
        return np.ones(len(points))

    def always_inside(points: np.ndarray) -> np.ndarray:
        return -np.ones(len(points))

    with pytest.raises(ExtractionError, match="no inside"):
        extract_surface(always_outside, (-1.0,) * 3, (1.0,) * 3, 16)
    with pytest.raises(ExtractionError, match="inside everywhere"):
        extract_surface(always_inside, (-1.0,) * 3, (1.0,) * 3, 16)


def test_non_finite_fields_are_refused() -> None:
    def broken(points: np.ndarray) -> np.ndarray:
        return np.full(len(points), np.nan)

    with pytest.raises(ExtractionError, match="non-finite"):
        extract_surface(broken, (-1.0,) * 3, (1.0,) * 3, 16)


def test_a_field_returning_the_wrong_number_of_values_is_refused() -> None:
    def wrong(points: np.ndarray) -> np.ndarray:
        return np.zeros(len(points) // 2)

    with pytest.raises(ExtractionError, match="returned"):
        sample_grid(wrong, (-1.0,) * 3, (1.0,) * 3, 8)


def test_padding_keeps_a_surface_away_from_the_grid_boundary() -> None:
    """A sphere touching the bound must still close, because of the padding cells."""
    mesh = extract_surface(sphere_field, (-0.5,) * 3, (0.5,) * 3, 32, margin=2)
    assert meshlib.is_watertight(mesh)
    assert meshlib.total_signed_volume(mesh) > 0.0
