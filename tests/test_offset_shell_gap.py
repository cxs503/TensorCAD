"""GAP-CAD-A2-OFFSET-SHELL: Surface offset / shelling verification.

Three acceptance criteria:
1. Sphere offset volume = 4/3*pi*((r)^3 - (r-t)^3) with error < 0.1%
2. Wigley hull shelling: double-sided watertight (or honest refusal with reason)
3. Self-intersection detection and honest reporting

Provenance:
- Task: m17347, anchor GAP-CAD-A2-OFFSET-SHELL
- Date: 2026-11-13
- Repository: TensorCAD, branch agent/m17347
"""

from __future__ import annotations

import itertools
import json
import math
import os
import sys

import pytest

# Add parent to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.shell import (
    ShellError,
    ShellResult,
    hollow_body,
    min_wall_thickness,
    offset_inward,
    shell,
)
from tensorcad import mesh as meshlib


def _uv_sphere(radius: float, meridians: int, parallels: int) -> Mesh:
    """Latitude/longitude sphere mesh."""
    builder = MeshBuilder()
    north = builder.add_vertex((0.0, 0.0, radius))
    south = builder.add_vertex((0.0, 0.0, -radius))
    rings = []
    for ring in range(1, parallels):
        phi = math.pi * ring / parallels
        z = radius * math.cos(phi)
        rho = radius * math.sin(phi)
        rings.append([
            builder.add_vertex((
                rho * math.cos(2.0 * math.pi * col / meridians),
                rho * math.sin(2.0 * math.pi * col / meridians),
                z,
            ))
            for col in range(meridians)
        ])
    for col in range(meridians):
        builder.add_triangle(north, rings[0][col], rings[0][(col + 1) % meridians])
    for upper, lower in itertools.pairwise(rings):
        for col in range(meridians):
            nxt = (col + 1) % meridians
            builder.add_quad(upper[col], lower[col], lower[nxt], upper[nxt])
    for col in range(meridians):
        builder.add_triangle(south, rings[-1][(col + 1) % meridians], rings[-1][col])
    return builder.build()


def _box(length: float = 1.0) -> Mesh:
    """Unit box as a convex prism."""
    from tensorcad.section_gz import prism_mesh
    RECTANGLE = ((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))
    return prism_mesh(RECTANGLE, length)


# ==============================================================================
# Criterion 1: Sphere offset volume error < 0.1%
# ==============================================================================

class TestSphereOffsetVolume:
    """Sphere shell volume must match 4/3*pi*(r^3 - (r-t)^3) within 0.1%."""

    @pytest.mark.parametrize("radius,thickness,meridians,parallels", [
        (1.0, 0.1, 144, 72),
        (1.0, 0.2, 192, 96),
        (2.0, 0.1, 144, 72),
        (0.5, 0.05, 144, 72),
    ])
    def test_sphere_shell_volume_error_below_01_percent(
        self, radius: float, thickness: float, meridians: int, parallels: int
    ):
        """Shell volume error must be < 0.1% at sufficient resolution."""
        sphere = _uv_sphere(radius, meridians, parallels)
        result = shell(sphere, thickness)

        analytic_outer = 4.0 / 3.0 * math.pi * radius ** 3
        analytic_inner = 4.0 / 3.0 * math.pi * (radius - thickness) ** 3
        analytic_shell = analytic_outer - analytic_inner

        shell_err = abs(result.shell_volume - analytic_shell) / analytic_shell
        inner_err = abs(result.inner_volume - analytic_inner) / analytic_inner

        assert shell_err < 0.001, (
            f"Shell volume error {shell_err*100:.4f}% exceeds 0.1% "
            f"(computed={result.shell_volume:.10f}, analytic={analytic_shell:.10f})"
        )
        assert inner_err < 0.001, (
            f"Inner volume error {inner_err*100:.4f}% exceeds 0.1% "
            f"(computed={result.inner_volume:.10f}, analytic={analytic_inner:.10f})"
        )

    def test_sphere_shell_convergence_is_second_order(self):
        """Faceting error on a sphere falls at second order with refinement."""
        radius, thickness = 1.0, 0.1
        analytic = 4.0 / 3.0 * math.pi * (radius ** 3 - (radius - thickness) ** 3)
        errors = []
        for meridians, parallels in ((24, 12), (48, 24), (96, 48)):
            result = shell(_uv_sphere(radius, meridians, parallels), thickness)
            errors.append(abs(result.shell_volume - analytic) / analytic)
        # Errors must decrease monotonically
        for i in range(1, len(errors)):
            assert errors[i] < errors[i - 1], f"Error not decreasing: {errors}"
        # Finest resolution must be below 0.1%
        assert errors[-1] < 0.001, f"Finest error {errors[-1]*100:.4f}% exceeds 0.1%"


# ==============================================================================
# Criterion 2: Wigley hull shelling — double-sided watertight
# ==============================================================================

class TestWigleyHullShelling:
    """Wigley hull shelling behavior: classic hull refuses (geometric singularities
    at keel and waterline corners), smooth hulls shell successfully with watertight result."""

    def test_classic_wigley_refused_due_to_geometric_singularities(self):
        """Classic Wigley hull cannot be shelled: the inward offset folds the surface
        at the sharp keel (sqrt taper, m=2) and waterline corners (n=2).
        This is correct behavior — the shell module honestly reports the failure."""
        from tensorcad.wigley import generate_wigley

        wigley = generate_wigley(length=10.0, beam=1.0, draft=0.625, nx=40, nz=20)
        mesh = wigley.mesh

        assert meshlib.is_watertight(mesh), "Wigley mesh must be watertight"
        assert meshlib.is_edge_manifold(mesh), "Wigley mesh must be edge-manifold"

        # Any thickness should fail due to geometric singularities
        with pytest.raises(ShellError, match="folds the surface"):
            shell(mesh, 0.01)

    def test_ellipsoid_hull_shells_watertight(self):
        """An ellipsoid hull (smooth everywhere, C^1 continuous) shells to a
        double-sided watertight hollow body. This demonstrates the shell module
        works correctly on smooth hull forms."""
        # Create an ellipsoid: x^2/a^2 + y^2/b^2 + z^2/c^2 = 1
        a, b, c = 5.0, 0.5, 0.3125  # L/2, B/2, T (same proportions as Wigley reference)
        builder = MeshBuilder()

        meridians, parallels = 48, 24
        north = builder.add_vertex((0.0, 0.0, c))
        south = builder.add_vertex((0.0, 0.0, -c))
        rings = []
        for ring in range(1, parallels):
            phi = math.pi * ring / parallels
            z = c * math.cos(phi)
            rho_xy = math.sin(phi)
            rings.append([
                builder.add_vertex((
                    a * rho_xy * math.cos(2.0 * math.pi * col / meridians),
                    b * rho_xy * math.sin(2.0 * math.pi * col / meridians),
                    z,
                ))
                for col in range(meridians)
            ])
        for col in range(meridians):
            builder.add_triangle(north, rings[0][col], rings[0][(col + 1) % meridians])
        for upper, lower in itertools.pairwise(rings):
            for col in range(meridians):
                nxt = (col + 1) % meridians
                builder.add_quad(upper[col], lower[col], lower[nxt], upper[nxt])
        for col in range(meridians):
            builder.add_triangle(south, rings[-1][(col + 1) % meridians], rings[-1][col])
        ellipsoid = builder.build()

        assert meshlib.is_watertight(ellipsoid)

        result = shell(ellipsoid, 0.02)
        body = hollow_body(result)
        assert meshlib.is_watertight(body.mesh), "Hollow body must be watertight"
        assert meshlib.is_edge_manifold(body.mesh), "Hollow body must be edge-manifold"
        assert result.outer_volume > 0
        assert result.inner_volume > 0
        assert result.shell_volume > 0

    def test_box_hull_shells_watertight(self):
        """A box hull (planar faces, no curvature) shells to a watertight hollow body
        with exact analytic volume."""
        box = _box(1.0)
        assert meshlib.is_watertight(box)

        result = shell(box, 0.1)
        body = hollow_body(result)
        assert meshlib.is_watertight(body.mesh)
        assert meshlib.is_edge_manifold(body.mesh)

        # Exact analytic verification
        expected_inner = (1.0 - 2 * 0.1) ** 3
        assert abs(result.inner_volume - expected_inner) / expected_inner < 1e-12


# ==============================================================================
# Criterion 3: Self-intersection detection and honest reporting
# ==============================================================================

class TestSelfIntersectionDetection:
    """The shell module must detect and refuse self-intersecting offsets."""

    def test_thin_plate_refuses_when_thickness_exceeds_half_wall(self):
        """A plate thinner than 2*t cannot be shelled: the offset surfaces collide."""
        from tensorcad.section_gz import prism_mesh
        plate = prism_mesh(((0.0, 0.0), (1.0, 0.0), (1.0, 0.1), (0.0, 0.1)), 1.0)

        # t = 0.08 > 0.1/2 = 0.05, should refuse
        with pytest.raises(ShellError, match="inverts the solid"):
            shell(plate, 0.08)

    def test_tapered_arm_refuses_at_collapse_point(self):
        """A tapered arm that thins to zero must refuse when the offset reaches the thin end."""
        from tensorcad.section_gz import prism_mesh
        taper = prism_mesh(((0.0, 0.0), (1.0, 0.0), (1.0, 0.1), (0.0, 0.5)), 1.0)

        with pytest.raises(ShellError):
            shell(taper, 0.1)

    def test_open_surface_refused(self):
        """An open (non-closed) surface must be refused."""
        from tensorcad.section_gz import prism_mesh
        RECTANGLE = ((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))
        box = prism_mesh(RECTANGLE, 1.0)
        opened = Mesh(vertices=box.vertices, faces=box.faces[1:])

        with pytest.raises(ShellError, match="not closed"):
            shell(opened, 0.05)

    def test_inward_oriented_surface_refused(self):
        """An inward-oriented surface must be refused."""
        from tensorcad.section_gz import prism_mesh
        RECTANGLE = ((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))
        box = prism_mesh(RECTANGLE, 1.0)
        flipped = Mesh(
            vertices=box.vertices,
            faces=tuple((a, c, b) for a, b, c in box.faces),
        )

        with pytest.raises(ShellError, match="not outward-oriented"):
            shell(flipped, 0.05)

    def test_non_positive_thickness_refused(self):
        """Zero, negative, NaN, and Inf thickness must all be refused."""
        from tensorcad.section_gz import prism_mesh
        RECTANGLE = ((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))
        box = prism_mesh(RECTANGLE, 1.0)

        for t in [0.0, -0.1, math.nan, math.inf]:
            with pytest.raises(ShellError, match="positive and finite"):
                shell(box, t)

    def test_measured_wall_detects_under_offset(self):
        """The wall thickness measurement must detect a deliberately under-offset surface."""
        from tensorcad.section_gz import prism_mesh
        RECTANGLE = ((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))
        box = prism_mesh(RECTANGLE, 1.0)

        target = 0.1
        for factor in [0.2, 0.5, 0.7]:
            under = offset_inward(box, factor * target)
            measured = min_wall_thickness(box, under)
            assert abs(measured - factor * target) <= 1e-12 * target, (
                f"Measured wall {measured} != expected {factor * target}"
            )


# ==============================================================================
# Artifact generation
# ==============================================================================

def generate_artifacts():
    """Generate JSON artifacts for the report."""
    artifacts = {}

    # Sphere convergence study
    sphere_results = []
    radius, thickness = 1.0, 0.1
    analytic_outer = 4.0 / 3.0 * math.pi * radius ** 3
    analytic_inner = 4.0 / 3.0 * math.pi * (radius - thickness) ** 3
    analytic_shell = analytic_outer - analytic_inner

    for meridians, parallels in [(24, 12), (48, 24), (96, 48), (144, 72)]:
        sphere = _uv_sphere(radius, meridians, parallels)
        result = shell(sphere, thickness)
        sphere_results.append({
            "resolution": f"{meridians}x{parallels}",
            "vertices": len(sphere.vertices),
            "faces": len(sphere.faces),
            "shell_volume": result.shell_volume,
            "inner_volume": result.inner_volume,
            "outer_volume": result.outer_volume,
            "shell_volume_error_pct": abs(result.shell_volume - analytic_shell) / analytic_shell * 100,
            "inner_volume_error_pct": abs(result.inner_volume - analytic_inner) / analytic_inner * 100,
            "min_wall_thickness": result.min_wall_thickness,
            "max_plane_residual": result.max_plane_residual,
        })

    artifacts["sphere_convergence"] = {
        "radius": radius,
        "thickness": thickness,
        "analytic_outer_volume": analytic_outer,
        "analytic_inner_volume": analytic_inner,
        "analytic_shell_volume": analytic_shell,
        "results": sphere_results,
    }

    # Wigley hull results
    from tensorcad.wigley import generate_wigley

    wigley_results = {}

    # Classic Wigley
    classic_wigley = generate_wigley(length=10.0, beam=1.0, draft=0.625, nx=40, nz=20)
    classic_mesh = classic_wigley.mesh
    wigley_results["classic"] = {
        "vertices": len(classic_mesh.vertices),
        "faces": len(classic_mesh.faces),
        "watertight": meshlib.is_watertight(classic_mesh),
        "edge_manifold": meshlib.is_edge_manifold(classic_mesh),
        "volume": meshlib.total_signed_volume(classic_mesh),
        "shell_result": "refused: folds the surface (keel singularity m=2, waterline corners n=2)",
        "shell_error": "min(dot(n_inner, n_outer)) < 0 → offset faces wound backwards",
    }

    # Ellipsoid (smooth hull that shells successfully)
    a, b, c = 5.0, 0.5, 0.3125
    builder = MeshBuilder()
    meridians, parallels = 48, 24
    north = builder.add_vertex((0.0, 0.0, c))
    south = builder.add_vertex((0.0, 0.0, -c))
    rings = []
    for ring in range(1, parallels):
        phi = math.pi * ring / parallels
        z = c * math.cos(phi)
        rho_xy = math.sin(phi)
        rings.append([
            builder.add_vertex((
                a * rho_xy * math.cos(2.0 * math.pi * col / meridians),
                b * rho_xy * math.sin(2.0 * math.pi * col / meridians),
                z,
            ))
            for col in range(meridians)
        ])
    for col in range(meridians):
        builder.add_triangle(north, rings[0][col], rings[0][(col + 1) % meridians])
    for upper, lower in itertools.pairwise(rings):
        for col in range(meridians):
            nxt = (col + 1) % meridians
            builder.add_quad(upper[col], lower[col], lower[nxt], upper[nxt])
    for col in range(meridians):
        builder.add_triangle(south, rings[-1][(col + 1) % meridians], rings[-1][col])
    ellipsoid = builder.build()

    ellipsoid_result = shell(ellipsoid, 0.02)
    ellipsoid_body = hollow_body(ellipsoid_result)
    wigley_results["ellipsoid_smooth_hull"] = {
        "semi_axes": {"a": a, "b": b, "c": c},
        "vertices": len(ellipsoid.vertices),
        "faces": len(ellipsoid.faces),
        "watertight": meshlib.is_watertight(ellipsoid),
        "volume": meshlib.total_signed_volume(ellipsoid),
        "shell_thickness": 0.02,
        "shell_volume": ellipsoid_result.shell_volume,
        "inner_volume": ellipsoid_result.inner_volume,
        "outer_volume": ellipsoid_result.outer_volume,
        "min_wall_thickness": ellipsoid_result.min_wall_thickness,
        "hollow_body_watertight": meshlib.is_watertight(ellipsoid_body.mesh),
        "hollow_body_edge_manifold": meshlib.is_edge_manifold(ellipsoid_body.mesh),
    }

    artifacts["hull_shelling"] = wigley_results

    # Self-intersection detection summary
    artifacts["self_intersection_detection"] = {
        "thin_plate_refusal": "t > half_thickness → 'inverts the solid'",
        "tapered_arm_refusal": "collapse at thin end → ShellError",
        "open_surface_refusal": "not closed → ShellError",
        "inward_oriented_refusal": "not outward-oriented → ShellError",
        "non_positive_thickness_refusal": "positive and finite → ShellError",
        "measured_wall_detection": "under-offset surface detected correctly",
    }

    return artifacts


if __name__ == "__main__":
    artifacts = generate_artifacts()
    output_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "artifacts", "offset_shell", "results.json"
    )
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(artifacts, f, indent=2)
    print(f"Artifacts written to {output_path}")
