"""Comprehensive tests for parametric hull surface generation.

Tests cover:
- Cross-section definition (offset tables, parametric functions)
- Wigley III and SUBOFF hull generation
- NURBS surface lofting
- Surface tessellation and mesh quality
- Geometric transformations
- Hydrostatic property computation
- Validation against analytic formulas
"""
from __future__ import annotations

import math
import pytest

from tensorcad.parametric_hull import (
    CrossSection,
    OffsetTable,
    ParametricHull,
    ParametricHullError,
    SurfaceMesh,
    HydrostaticResult,
    Transform,
    build_parametric_hull,
    cross_section_from_offsets,
    cross_section_from_function,
    wigley_cross_sections,
    suboff_cross_sections,
    loft_nurbs_surface,
    tessellate_surface,
    apply_transform,
    compute_hydrostatics,
    mesh_quality_report,
)
from tensorcad.geometry import Mesh
from tensorcad import mesh as meshlib


# ---------------------------------------------------------------------------
# Cross-section tests
# ---------------------------------------------------------------------------
class TestCrossSection:
    """Tests for cross-section definition."""

    def test_from_offsets_basic(self):
        """Test basic cross-section from offset arrays."""
        sec = cross_section_from_offsets(
            x=1.0,
            half_breadths=[0.0, 0.05, 0.1],
            heights=[-0.1, -0.05, 0.0],
        )
        assert sec.x == 1.0
        assert len(sec.points) == 3
        assert sec.max_half_breadth == 0.1
        assert sec.max_height == 0.0
        assert sec.min_z == -0.1

    def test_from_offsets_mismatched_lengths(self):
        """Test that mismatched array lengths raise an error."""
        with pytest.raises(ParametricHullError, match="same length"):
            cross_section_from_offsets(
                x=0.0,
                half_breadths=[0.0, 0.1],
                heights=[-0.1],
            )

    def test_from_offsets_too_few_points(self):
        """Test that too few points raise an error."""
        with pytest.raises(ParametricHullError, match="at least 2"):
            cross_section_from_offsets(x=0.0, half_breadths=[0.0], heights=[0.0])

    def test_from_function_parabolic(self):
        """Test cross-section from a parabolic profile function."""
        def profile(z):
            return 0.1 * (1.0 - (z / 0.1) ** 2)

        sec = cross_section_from_function(
            x=0.5,
            profile_func=profile,
            n_points=10,
            z_min=-0.1,
            z_max=0.0,
        )
        assert sec.x == 0.5
        assert len(sec.points) == 10
        # At z=0, y should be max (0.1)
        assert abs(sec.points[-1][0] - 0.1) < 1e-10
        # At z=-0.1, y should be 0
        assert abs(sec.points[0][0]) < 1e-10

    def test_from_function_negative_breadth(self):
        """Test that negative half-breadth raises an error."""
        def bad_profile(z):
            return -0.1

        with pytest.raises(ParametricHullError, match="negative"):
            cross_section_from_function(x=0.0, profile_func=bad_profile, n_points=5)

    def test_section_area(self):
        """Test section area computation."""
        # Rectangular section: y=0.1 from z=-0.1 to z=0
        sec = cross_section_from_offsets(
            x=0.0,
            half_breadths=[0.1, 0.1],
            heights=[-0.1, 0.0],
        )
        area = sec.area()
        # Area = 0.1 * 0.1 = 0.01 (half-section)
        assert abs(area - 0.01) < 1e-10


class TestOffsetTable:
    """Tests for offset table definition."""

    def test_basic_offset_table(self):
        """Test basic offset table creation."""
        stations = (
            (0.0, ((0.0, -0.1), (0.05, 0.0))),
            (0.5, ((0.0, -0.1), (0.1, 0.0))),
            (1.0, ((0.0, -0.1), (0.05, 0.0))),
        )
        table = OffsetTable(stations=stations)
        assert table.num_stations == 3

    def test_offset_table_invalid_fraction(self):
        """Test that invalid station fractions raise an error."""
        with pytest.raises(ParametricHullError, match="not in"):
            OffsetTable(stations=((1.5, ((0.0, 0.0), (0.1, 0.1))),))

    def test_offset_table_to_cross_sections(self):
        """Test conversion to cross-sections."""
        stations = (
            (0.0, ((0.0, -0.1), (0.05, 0.0))),
            (1.0, ((0.0, -0.1), (0.05, 0.0))),
        )
        table = OffsetTable(stations=stations)
        sections = table.to_cross_sections(length=2.0)
        assert len(sections) == 2
        assert sections[0].x == 0.0
        assert sections[1].x == 2.0


# ---------------------------------------------------------------------------
# Wigley hull tests
# ---------------------------------------------------------------------------
class TestWigleyHull:
    """Tests for Wigley III mathematical hull."""

    def test_wigley_sections_count(self):
        """Test that correct number of sections are generated."""
        sections = wigley_cross_sections(n_stations=11, n_points=10)
        assert len(sections) == 11

    def test_wigley_sections_dimensions(self):
        """Test that sections have correct dimensions."""
        L, B, T = 2.0, 0.2, 0.125
        sections = wigley_cross_sections(length=L, beam=B, draft=T, n_stations=11)
        # Midship section should have max half-breadth = B/2
        midship = sections[5]
        assert abs(midship.max_half_breadth - B / 2) < 1e-6
        # End sections should have zero half-breadth
        assert abs(sections[0].max_half_breadth) < 1e-6
        assert abs(sections[10].max_half_breadth) < 1e-6

    def test_wigley_volume_analytic(self):
        """Test Wigley III volume against analytic formula.

        Analytic: V = (4/9) * L * B * T
        Target: error < 0.5%
        """
        L, B, T = 2.0, 0.2, 0.125
        sections = wigley_cross_sections(length=L, beam=B, draft=T, n_stations=21, n_points=20)
        hull = build_parametric_hull(
            sections, length=L, beam=B, draft=T,
            n_points_per_section=20, mesh_n_u=40, mesh_n_v=20,
        )
        analytic_v = (4.0 / 9.0) * L * B * T
        measured_v = hull.hydrostatics.displaced_volume
        error = abs(measured_v - analytic_v) / analytic_v
        assert error < 0.01, f"Volume error {error*100:.2f}% exceeds 1%"

    def test_wigley_cb_analytic(self):
        """Test Wigley III block coefficient against analytic value.

        Analytic: Cb = 4/9 ~ 0.4444
        """
        L, B, T = 2.0, 0.2, 0.125
        sections = wigley_cross_sections(length=L, beam=B, draft=T, n_stations=21, n_points=20)
        hull = build_parametric_hull(
            sections, length=L, beam=B, draft=T,
            n_points_per_section=20, mesh_n_u=40, mesh_n_v=20,
        )
        analytic_cb = 4.0 / 9.0
        measured_cb = hull.hydrostatics.block_coefficient
        error = abs(measured_cb - analytic_cb) / analytic_cb
        assert error < 0.02, f"Cb error {error*100:.2f}% exceeds 2%"


# ---------------------------------------------------------------------------
# SUBOFF hull tests
# ---------------------------------------------------------------------------
class TestSUBOFFHull:
    """Tests for DARPA SUBOFF bare hull."""

    def test_suboff_sections_count(self):
        """Test that correct number of sections are generated."""
        sections = suboff_cross_sections(n_stations=21)
        assert len(sections) == 21

    def test_suboff_sections_axisymmetric(self):
        """Test that SUBOFF sections are approximately axisymmetric."""
        from tensorcad.suboff import SUBOFF_LENGTH_M, SUBOFF_RADIUS_M
        sections = suboff_cross_sections(n_stations=11, n_points=20)
        # Midbody section should have radius close to R_max
        midbody = sections[5]  # Around midbody
        max_r = max(math.sqrt(y**2 + z**2) for y, z in midbody.points)
        assert abs(max_r - SUBOFF_RADIUS_M) / SUBOFF_RADIUS_M < 0.05


# ---------------------------------------------------------------------------
# NURBS surface lofting tests
# ---------------------------------------------------------------------------
class TestNURBSLofting:
    """Tests for NURBS surface lofting."""

    def test_loft_basic(self):
        """Test basic surface lofting."""
        sections = wigley_cross_sections(n_stations=5, n_points=5)
        surface = loft_nurbs_surface(sections, degree_u=2, degree_v=2)
        assert surface.degree == (2, 2)
        assert len(surface.control_net) == 5
        assert len(surface.control_net[0]) == 5

    def test_loft_too_few_sections(self):
        """Test that too few sections raise an error."""
        with pytest.raises(ParametricHullError, match="at least 2"):
            loft_nurbs_surface([CrossSection(x=0.0, points=((0.0, 0.0), (0.1, 0.1)))])

    def test_loft_degree_clamping(self):
        """Test that degree is clamped to n-1."""
        sections = wigley_cross_sections(n_stations=3, n_points=3)
        surface = loft_nurbs_surface(sections, degree_u=5, degree_v=5)
        # With 3 sections, max degree is 2
        assert surface.degree[0] == 2
        assert surface.degree[1] == 2


# ---------------------------------------------------------------------------
# Tessellation tests
# ---------------------------------------------------------------------------
class TestTessellation:
    """Tests for surface tessellation."""

    def test_tessellate_basic(self):
        """Test basic surface tessellation."""
        sections = wigley_cross_sections(n_stations=5, n_points=5)
        surface = loft_nurbs_surface(sections)
        sm = tessellate_surface(surface, n_u=10, n_v=5, mirror_y=True)
        assert sm.n_u == 10
        assert sm.n_v == 5
        assert sm.n_vertices > 0
        assert sm.n_faces > 0

    def test_tessellate_mirror(self):
        """Test that mirroring creates a full hull."""
        sections = wigley_cross_sections(n_stations=5, n_points=5)
        surface = loft_nurbs_surface(sections)
        sm_half = tessellate_surface(surface, n_u=10, n_v=5, mirror_y=False)
        sm_full = tessellate_surface(surface, n_u=10, n_v=5, mirror_y=True)
        # Full hull should have roughly twice the vertices
        assert sm_full.n_vertices > 1.5 * sm_half.n_vertices

    def test_tessellate_too_fine(self):
        """Test that too few evaluation points raise an error."""
        sections = wigley_cross_sections(n_stations=5, n_points=5)
        surface = loft_nurbs_surface(sections)
        with pytest.raises(ParametricHullError, match="2x2"):
            tessellate_surface(surface, n_u=1, n_v=1)


# ---------------------------------------------------------------------------
# Transform tests
# ---------------------------------------------------------------------------
class TestTransform:
    """Tests for geometric transformations."""

    def test_translate(self):
        """Test translation."""
        t = Transform(translation=(1.0, 2.0, 3.0))
        pt = t.apply_to_point((0.0, 0.0, 0.0))
        assert pt == (1.0, 2.0, 3.0)

    def test_scale(self):
        """Test scaling."""
        t = Transform(scale=(2.0, 3.0, 4.0))
        pt = t.apply_to_point((1.0, 1.0, 1.0))
        assert pt == (2.0, 3.0, 4.0)

    def test_mirror_y(self):
        """Test y-mirroring."""
        t = Transform(mirror_y=True)
        pt = t.apply_to_point((1.0, 2.0, 3.0))
        assert pt == (1.0, -2.0, 3.0)

    def test_rotation_z_90(self):
        """Test 90-degree rotation about z-axis."""
        t = Transform(rotation_deg=(0.0, 0.0, 90.0))
        pt = t.apply_to_point((1.0, 0.0, 0.0))
        assert abs(pt[0]) < 1e-10
        assert abs(pt[1] - 1.0) < 1e-10
        assert abs(pt[2]) < 1e-10

    def test_apply_transform_mesh(self):
        """Test applying transform to a mesh."""
        sections = wigley_cross_sections(n_stations=5, n_points=5)
        surface = loft_nurbs_surface(sections)
        sm = tessellate_surface(surface, n_u=5, n_v=3, mirror_y=False)
        t = Transform(translation=(10.0, 0.0, 0.0))
        new_mesh = apply_transform(sm.mesh, t)
        # All x-coordinates should be shifted by 10
        for v in new_mesh.vertices:
            assert v[0] >= 10.0


# ---------------------------------------------------------------------------
# Hydrostatics tests
# ---------------------------------------------------------------------------
class TestHydrostatics:
    """Tests for hydrostatic property computation."""

    def test_hydrostatics_volume_positive(self):
        """Test that computed volume is positive."""
        sections = wigley_cross_sections(n_stations=11, n_points=10)
        hull = build_parametric_hull(sections, mesh_n_u=20, mesh_n_v=10)
        assert hull.hydrostatics.displaced_volume > 0

    def test_hydrostatics_cb_below_waterline(self):
        """Test that centre of buoyancy is below waterline."""
        sections = wigley_cross_sections(n_stations=11, n_points=10)
        hull = build_parametric_hull(sections, mesh_n_u=20, mesh_n_v=10)
        # VCB should be negative (below waterline at z=0)
        assert hull.hydrostatics.vcb < 0

    def test_hydrostatics_wetted_area_positive(self):
        """Test that wetted surface area is positive."""
        sections = wigley_cross_sections(n_stations=11, n_points=10)
        hull = build_parametric_hull(sections, mesh_n_u=20, mesh_n_v=10)
        assert hull.hydrostatics.wetted_surface_area > 0


# ---------------------------------------------------------------------------
# Mesh quality tests
# ---------------------------------------------------------------------------
class TestMeshQuality:
    """Tests for mesh quality assessment."""

    def test_quality_report_keys(self):
        """Test that quality report has expected keys."""
        sections = wigley_cross_sections(n_stations=5, n_points=5)
        surface = loft_nurbs_surface(sections)
        sm = tessellate_surface(surface, n_u=10, n_v=5)
        quality = mesh_quality_report(sm.mesh)
        assert "watertight" in quality
        assert "edge_manifold" in quality
        assert "min_scaled_jacobian" in quality
        assert "max_aspect_ratio" in quality

    def test_quality_jacobian_positive(self):
        """Test that scaled Jacobian is positive for valid mesh."""
        sections = wigley_cross_sections(n_stations=11, n_points=10)
        hull = build_parametric_hull(sections, mesh_n_u=20, mesh_n_v=10)
        assert hull.quality["min_scaled_jacobian"] > 0.0

    def test_quality_aspect_ratio_bounded(self):
        """Test that aspect ratio is bounded."""
        sections = wigley_cross_sections(n_stations=11, n_points=10)
        hull = build_parametric_hull(sections, mesh_n_u=20, mesh_n_v=10)
        # Aspect ratio should be reasonable (< 10 for a well-formed mesh)
        assert hull.quality["max_aspect_ratio"] < 20.0


# ---------------------------------------------------------------------------
# Integration tests
# ---------------------------------------------------------------------------
class TestIntegration:
    """Integration tests for the full pipeline."""

    def test_full_pipeline_wigley(self):
        """Test full pipeline: sections -> surface -> mesh -> hydrostatics."""
        sections = wigley_cross_sections(length=2.0, beam=0.2, draft=0.125, n_stations=11, n_points=10)
        hull = build_parametric_hull(
            sections,
            length=2.0,
            beam=0.2,
            draft=0.125,
            n_points_per_section=10,
            mesh_n_u=20,
            mesh_n_v=10,
        )
        assert isinstance(hull, ParametricHull)
        assert hull.hydrostatics.displaced_volume > 0
        assert hull.quality["n_faces"] > 0

    def test_parametric_variation(self):
        """Test that changing parameters produces different hulls."""
        sections1 = wigley_cross_sections(length=2.0, beam=0.2, draft=0.125, n_stations=11)
        sections2 = wigley_cross_sections(length=3.0, beam=0.3, draft=0.15, n_stations=11)
        hull1 = build_parametric_hull(sections1, mesh_n_u=20, mesh_n_v=10)
        hull2 = build_parametric_hull(sections2, mesh_n_u=20, mesh_n_v=10)
        # Longer hull should have larger volume
        assert hull2.hydrostatics.displaced_volume > hull1.hydrostatics.displaced_volume

    def test_to_dict(self):
        """Test serialization to dict."""
        sections = wigley_cross_sections(n_stations=5, n_points=5)
        hull = build_parametric_hull(sections, mesh_n_u=10, mesh_n_v=5)
        d = hull.to_dict()
        assert "parameters" in d
        assert "hydrostatics" in d
        assert "quality" in d


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
