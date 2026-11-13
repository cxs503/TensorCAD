"""Tests for the geometry contract validator (tensorcad.geom_contract)."""
from __future__ import annotations

import math
import pytest

from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.geom_contract import (
    ContractViolation,
    ContractReport,
    VALID_UNITS,
    MIN_DIAGONAL_M,
    MAX_DIAGONAL_M,
    validate_contract,
    _check_coordinate_system,
    _check_units,
    _check_scale,
    _check_orientation,
)
from tensorcad import mesh as meshlib


def _make_cube(size: float = 1.0) -> Mesh:
    """Build a unit cube mesh with outward normals."""
    builder = MeshBuilder(units="m", coordinate_frame="right-handed")
    s = size
    # 6 faces, each as 2 triangles
    # Front (z = s/2)
    v = [
        builder.add_vertex((-s/2, -s/2, s/2)),
        builder.add_vertex((s/2, -s/2, s/2)),
        builder.add_vertex((s/2, s/2, s/2)),
        builder.add_vertex((-s/2, s/2, s/2)),
    ]
    builder.add_triangle(v[0], v[1], v[2])
    builder.add_triangle(v[0], v[2], v[3])
    # Back (z = -s/2)
    v = [
        builder.add_vertex((-s/2, -s/2, -s/2)),
        builder.add_vertex((-s/2, s/2, -s/2)),
        builder.add_vertex((s/2, s/2, -s/2)),
        builder.add_vertex((s/2, -s/2, -s/2)),
    ]
    builder.add_triangle(v[0], v[1], v[2])
    builder.add_triangle(v[0], v[2], v[3])
    # Right (x = s/2)
    v = [
        builder.add_vertex((s/2, -s/2, -s/2)),
        builder.add_vertex((s/2, s/2, -s/2)),
        builder.add_vertex((s/2, s/2, s/2)),
        builder.add_vertex((s/2, -s/2, s/2)),
    ]
    builder.add_triangle(v[0], v[1], v[2])
    builder.add_triangle(v[0], v[2], v[3])
    # Left (x = -s/2)
    v = [
        builder.add_vertex((-s/2, -s/2, -s/2)),
        builder.add_vertex((-s/2, -s/2, s/2)),
        builder.add_vertex((-s/2, s/2, s/2)),
        builder.add_vertex((-s/2, s/2, -s/2)),
    ]
    builder.add_triangle(v[0], v[1], v[2])
    builder.add_triangle(v[0], v[2], v[3])
    # Top (y = s/2)
    v = [
        builder.add_vertex((-s/2, s/2, -s/2)),
        builder.add_vertex((-s/2, s/2, s/2)),
        builder.add_vertex((s/2, s/2, s/2)),
        builder.add_vertex((s/2, s/2, -s/2)),
    ]
    builder.add_triangle(v[0], v[1], v[2])
    builder.add_triangle(v[0], v[2], v[3])
    # Bottom (y = -s/2)
    v = [
        builder.add_vertex((-s/2, -s/2, -s/2)),
        builder.add_vertex((s/2, -s/2, -s/2)),
        builder.add_vertex((s/2, -s/2, s/2)),
        builder.add_vertex((-s/2, -s/2, s/2)),
    ]
    builder.add_triangle(v[0], v[1], v[2])
    builder.add_triangle(v[0], v[2], v[3])
    return builder.build()


class TestCheckCoordinateSystem:
    def test_right_handed_passes(self):
        mesh = _make_cube()
        violations = _check_coordinate_system(mesh)
        assert not violations

    def test_left_handed_fails(self):
        mesh = Mesh(
            vertices=_make_cube().vertices,
            faces=_make_cube().faces,
            units="m",
            coordinate_frame="left-handed",
        )
        violations = _check_coordinate_system(mesh)
        assert len(violations) == 1
        assert violations[0].severity == "error"
        assert "left-handed" in violations[0].message

    def test_unspecified_warns(self):
        mesh = Mesh(
            vertices=_make_cube().vertices,
            faces=_make_cube().faces,
            units="m",
            coordinate_frame="some unknown frame",
        )
        violations = _check_coordinate_system(mesh)
        assert len(violations) == 1
        assert violations[0].severity == "warning"


class TestCheckUnits:
    @pytest.mark.parametrize("unit", VALID_UNITS)
    def test_valid_units_pass(self, unit):
        mesh = Mesh(vertices=_make_cube().vertices, faces=_make_cube().faces, units=unit)
        violations = _check_units(mesh)
        assert not violations

    def test_invalid_unit_fails(self):
        mesh = Mesh(vertices=_make_cube().vertices, faces=_make_cube().faces, units="inch")
        violations = _check_units(mesh)
        assert len(violations) == 1
        assert violations[0].severity == "error"
        assert "inch" in violations[0].message


class TestCheckScale:
    def test_normal_scale_passes(self):
        mesh = _make_cube(size=1.0)  # diagonal = sqrt(3) ≈ 1.73 m
        violations = _check_scale(mesh)
        assert not violations

    def test_tiny_model_fails(self):
        mesh = _make_cube(size=0.001)  # diagonal ≈ 0.0017 m < MIN_DIAGONAL_M
        violations = _check_scale(mesh)
        assert len(violations) == 1
        assert violations[0].severity == "error"
        assert "below minimum" in violations[0].message

    def test_huge_model_fails(self):
        mesh = _make_cube(size=1000.0)  # diagonal ≈ 1732 m > MAX_DIAGONAL_M
        violations = _check_scale(mesh)
        assert len(violations) == 1
        assert violations[0].severity == "error"
        assert "exceeds maximum" in violations[0].message

    def test_empty_mesh_fails(self):
        mesh = Mesh(vertices=(), faces=(), units="m")
        violations = _check_scale(mesh)
        assert len(violations) == 1
        assert "empty mesh" in violations[0].message


class TestCheckOrientation:
    def test_outward_normals_pass(self):
        mesh = _make_cube()
        violations = _check_orientation(mesh)
        assert not violations

    def test_inward_normals_fail(self):
        cube = _make_cube()
        flipped = tuple((f[2], f[1], f[0]) for f in cube.faces)
        mesh = Mesh(vertices=cube.vertices, faces=flipped, units="m")
        violations = _check_orientation(mesh)
        assert len(violations) == 1
        assert violations[0].severity == "error"
        assert "non-positive" in violations[0].message

    def test_no_faces_fails(self):
        mesh = Mesh(vertices=((0, 0, 0), (1, 0, 0), (0, 1, 0)), faces=(), units="m")
        violations = _check_orientation(mesh)
        assert len(violations) == 1
        assert "no faces" in violations[0].message


class TestValidateContract:
    def test_valid_cube_passes(self):
        mesh = _make_cube()
        report = validate_contract(mesh, source="test.stl")
        assert report.passed
        assert not report.violations

    def test_inverted_cube_fails(self):
        cube = _make_cube()
        flipped = tuple((f[2], f[1], f[0]) for f in cube.faces)
        mesh = Mesh(vertices=cube.vertices, faces=flipped, units="m", coordinate_frame="left-handed")
        report = validate_contract(mesh, source="test.stl")
        assert not report.passed
        errors = [v for v in report.violations if v.severity == "error"]
        assert len(errors) >= 1  # at least orientation error

    def test_report_contains_mesh_info(self):
        mesh = _make_cube()
        report = validate_contract(mesh, source="test.stl")
        assert "vertex_count" in report.mesh_info
        assert "triangle_count" in report.mesh_info
        assert "watertight" in report.mesh_info
        assert "diagonal_m" in report.mesh_info
        assert "signed_volume" in report.mesh_info

    def test_report_to_dict(self):
        mesh = _make_cube()
        report = validate_contract(mesh, source="test.stl")
        d = report.to_dict()
        assert d["source"] == "test.stl"
        assert d["passed"] is True
        assert isinstance(d["violations"], list)
        assert isinstance(d["mesh_info"], dict)


class TestContractViolation:
    def test_to_dict(self):
        v = ContractViolation(
            convention="units",
            severity="error",
            message="bad unit",
            detail="expected m",
        )
        d = v.to_dict()
        assert d["convention"] == "units"
        assert d["severity"] == "error"
        assert d["message"] == "bad unit"
        assert d["detail"] == "expected m"

    def test_to_dict_no_detail(self):
        v = ContractViolation(convention="scale", severity="error", message="too small")
        d = v.to_dict()
        assert "detail" not in d
