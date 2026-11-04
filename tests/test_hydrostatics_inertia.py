"""Verification of the waterplane second moments on the ``Hydrostatics`` record.

Three independent references, because each one closes a different hole:

* A **rectangular barge** has closed-form waterplane inertias (``L*B**3/12`` and
  ``B*L**3/12``). The barge mesh is built here vertex by vertex, so the answer is
  checked against arithmetic rather than against another mesh, and it is checked on
  three tessellations so a lucky cancellation cannot pass.
* A **circular waterplane** (a vertical cylinder) has ``I_T = I_L = pi*R**4/4``; the
  inscribed polygon used to mesh it converges to that from below, which tests the
  transverse and longitudinal paths with a shape whose surface has no flat bottom
  and no station grid.
* The **analytic hull** of :mod:`tensorcad.hull` reports the same quantities from
  the closed forms in :mod:`tensorcad.hydrostatics_analytic`, so the mesh integrator
  can be measured at the resolutions users actually run.

Everything is standard library plus the package itself.
"""
from __future__ import annotations

import itertools
import math

import pytest

from tensorcad import HullParameters, generate_hull
from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.hull import analytic_hydrostatics
from tensorcad.hydrostatics import compute_hydrostatics
from tensorcad.hydrostatics_analytic import (
    tanh_sinh_quadrature,
    waterplane_inertia_longitudinal,
    waterplane_inertia_transverse,
)


def box_mesh(
    length: float, beam: float, depth: float, nx: int = 1, ny: int = 1, nz: int = 1
) -> Mesh:
    """Outward-oriented barge: ``x`` aft->fore, ``y`` across, ``z`` keel (-depth) -> 0.

    The waterline is the top face at ``z = 0``, so the submerged body is the whole
    box and every quantity has a closed form. Vertices are shared by construction
    through :class:`MeshBuilder`, so the surface is watertight without any welding.
    """
    builder = MeshBuilder(units="m", coordinate_frame="box test frame")

    def corner(i: int, j: int, k: int) -> int:
        """Vertex index of a grid corner (``MeshBuilder`` deduplicates by position)."""
        return builder.add_vertex(
            (length * i / nx, -0.5 * beam + beam * j / ny, -depth + depth * k / nz)
        )

    def quad(a: int, b: int, c: int, d: int) -> None:
        builder.add_quad(a, b, c, d)

    for i in range(nx):
        for j in range(ny):
            # bottom (outward -z) and top (outward +z)
            quad(
                corner(i, j, 0), corner(i, j + 1, 0), corner(i + 1, j + 1, 0), corner(i + 1, j, 0)
            )
            quad(
                corner(i, j, nz),
                corner(i + 1, j, nz),
                corner(i + 1, j + 1, nz),
                corner(i, j + 1, nz),
            )
    for j in range(ny):
        for k in range(nz):
            quad(corner(0, j, k), corner(0, j, k + 1), corner(0, j + 1, k + 1), corner(0, j + 1, k))
            quad(
                corner(nx, j, k),
                corner(nx, j + 1, k),
                corner(nx, j + 1, k + 1),
                corner(nx, j, k + 1),
            )
    for i in range(nx):
        for k in range(nz):
            quad(corner(i, 0, k), corner(i + 1, 0, k), corner(i + 1, 0, k + 1), corner(i, 0, k + 1))
            quad(
                corner(i, ny, k),
                corner(i, ny, k + 1),
                corner(i + 1, ny, k + 1),
                corner(i + 1, ny, k),
            )
    return builder.build()


def cylinder_mesh(radius: float, depth: float, segments: int) -> Mesh:
    """Vertical cylinder, closed, with the waterline as its top face at ``z = 0``."""
    builder = MeshBuilder(units="m", coordinate_frame="cylinder test frame")
    angles = [2.0 * math.pi * i / segments for i in range(segments)]

    def rim(theta: float, z: float) -> int:
        return builder.add_vertex(
            (radius * math.cos(theta), radius * math.sin(theta), z)
        )

    bottom = builder.add_vertex((0.0, 0.0, -depth))
    top = builder.add_vertex((0.0, 0.0, 0.0))
    for index in range(segments):
        a = angles[index]
        b = angles[(index + 1) % segments]
        builder.add_quad(rim(a, -depth), rim(b, -depth), rim(b, 0.0), rim(a, 0.0))
        builder.add_triangle(bottom, rim(b, -depth), rim(a, -depth))
        builder.add_triangle(top, rim(a, 0.0), rim(b, 0.0))
    return builder.build()


# --------------------------------------------------------------------------- barge


@pytest.mark.parametrize("tessellation", [(1, 1, 1), (4, 3, 2), (20, 6, 4)])
def test_rectangular_barge_is_exact(tessellation: tuple[int, int, int]) -> None:
    """A box waterplane is a rectangle: the record's inertias have exact closed forms."""
    length, beam, depth = 20.0, 6.0, 2.0
    mesh = box_mesh(length, beam, depth, *tessellation)
    hydro = compute_hydrostatics(mesh, 0.0, water_density=1000.0)

    assert hydro.basis == "closed surface integrals"
    assert hydro.displacement_volume == pytest.approx(length * beam * depth, rel=1e-12)
    assert hydro.waterplane_area == pytest.approx(length * beam, rel=1e-12)
    assert hydro.transverse_waterplane_inertia == pytest.approx(length * beam**3 / 12.0, rel=1e-12)
    assert hydro.longitudinal_waterplane_inertia == pytest.approx(beam * length**3 / 12.0, rel=1e-12)
    assert hydro.keel_z == pytest.approx(-depth, abs=1e-12)


def test_barge_inertia_does_not_depend_on_tessellation() -> None:
    """Refining a flat waterplane must not move the integrals at all."""
    coarse = compute_hydrostatics(box_mesh(12.0, 4.0, 1.5, 1, 1, 1), 0.0)
    fine = compute_hydrostatics(box_mesh(12.0, 4.0, 1.5, 6, 4, 3), 0.0)
    assert coarse.transverse_waterplane_inertia == pytest.approx(
        fine.transverse_waterplane_inertia, rel=1e-12
    )
    assert coarse.longitudinal_waterplane_inertia == pytest.approx(
        fine.longitudinal_waterplane_inertia, rel=1e-12
    )


# ------------------------------------------------------------------------ cylinder


def test_circular_waterplane_converges_to_the_closed_form() -> None:
    """An inscribed polygon approaches ``pi*R**4/4`` from below, and stays below."""
    radius, depth = 1.5, 2.0
    exact = math.pi * radius**4 / 4.0
    errors = []
    for segments in (24, 48, 96, 192):
        hydro = compute_hydrostatics(cylinder_mesh(radius, depth, segments), 0.0)
        transverse = hydro.transverse_waterplane_inertia
        longitudinal = hydro.longitudinal_waterplane_inertia
        # A waterplane that is symmetric in both axes must give the same two
        # moments, and an inscribed polygon can only be smaller than the circle.
        assert transverse == pytest.approx(longitudinal, rel=1e-9)
        assert transverse < exact
        errors.append(exact - transverse)
    # Inscribed regular polygon: the area deficit is O(1/n**2), so each doubling of
    # the segment count should cut the error by about four.
    ratios = [coarse / fine for coarse, fine in itertools.pairwise(errors)]
    assert all(3.5 < ratio < 4.5 for ratio in ratios), ratios
    assert errors[-1] / exact < 5e-4


# -------------------------------------------------------------------- analytic hull


def test_analytic_closed_forms_match_their_own_quadrature() -> None:
    """Transverse: gamma identity vs quadrature of ``int (2/3)*b**3 dx``.

    The closed form and the quadrature are derived separately (a gamma function and
    a double-exponential rule), so agreement is a real check rather than a
    restatement of the same algebra.
    """
    length, beam, exponent = 10.0, 2.0, 0.7
    closed = waterplane_inertia_transverse(length, beam, breadth_exponent=exponent)

    def integrand(xi: float) -> float:
        return (2.0 / 3.0) * (0.5 * beam) ** 3 * math.sin(math.pi * xi) ** (3.0 * exponent)

    numerical = length * tanh_sinh_quadrature(integrand, 0.0, 1.0, levels=7)
    assert closed == pytest.approx(numerical, rel=1e-12)


def test_longitudinal_inertia_quadrature_reports_its_own_error() -> None:
    """The reference carries its convergence evidence instead of asking to be believed."""
    value, error = waterplane_inertia_longitudinal(10.0, 2.0)
    assert value > 0.0
    assert error < 1e-10 * value
    # Scaling laws: I_L goes with beam and with length cubed.
    doubled = waterplane_inertia_longitudinal(20.0, 2.0)[0] / 8.0
    assert doubled == pytest.approx(value, rel=1e-9)


@pytest.mark.parametrize("parameters", [(10.0, 2.0, 0.8), (60.0, 14.0, 6.0), (4.3561, 1.0, 0.4)])
def test_hull_mesh_inertia_matches_the_analytic_surface(
    parameters: tuple[float, float, float]
) -> None:
    """The sampled hull against the closed forms of the surface it samples."""
    length, beam, draft = parameters
    reference = analytic_hydrostatics(length, beam, draft)

    default = generate_hull(HullParameters(length=length, beam=beam, draft=draft))
    measured = default.hydrostatics
    assert measured.transverse_waterplane_inertia == pytest.approx(
        reference.transverse_waterplane_inertia, rel=0.03
    )
    assert measured.longitudinal_waterplane_inertia == pytest.approx(
        reference.longitudinal_waterplane_inertia, rel=0.03
    )

    fine = generate_hull(HullParameters(length=length, beam=beam, draft=draft, nx=161, nz=41))
    assert fine.hydrostatics.transverse_waterplane_inertia == pytest.approx(
        reference.transverse_waterplane_inertia, rel=0.005
    )
    assert fine.hydrostatics.longitudinal_waterplane_inertia == pytest.approx(
        reference.longitudinal_waterplane_inertia, rel=0.005
    )
    # Refining has to help, not just land inside the band by luck.
    assert abs(
        fine.hydrostatics.transverse_waterplane_inertia
        - reference.transverse_waterplane_inertia
    ) < abs(measured.transverse_waterplane_inertia - reference.transverse_waterplane_inertia)


def test_hydrostatics_payload_carries_the_inertias() -> None:
    record = generate_hull().hydrostatics.to_dict()
    assert record["waterplane_inertia_transverse_m4"] > 0.0
    assert record["waterplane_inertia_longitudinal_m4"] > record["waterplane_inertia_transverse_m4"]
    assert record["keel_z_m"] == pytest.approx(-0.8, abs=1e-9)


# ------------------------------------------------------------------------ rejection


def test_a_body_without_a_waterplane_reports_nan() -> None:
    """A closed surface entirely below the waterline has no free surface to float on."""
    mesh = box_mesh(4.0, 2.0, 1.0, 2, 2, 2)
    submerged = compute_hydrostatics(mesh, waterline_z=2.0)
    assert math.isnan(submerged.transverse_waterplane_inertia)
    assert math.isnan(submerged.longitudinal_waterplane_inertia)


def test_open_surface_inertia_is_not_claimed() -> None:
    """An unclosed hull shell gives no trustworthy moments; the value is nan."""
    mesh = generate_hull().mesh
    open_faces = tuple(face for face in mesh.faces if max(mesh.vertices[v][2] for v in face) < -0.2)
    shell = Mesh(vertices=mesh.vertices, faces=open_faces)
    hydro = compute_hydrostatics(shell, 0.0)
    assert hydro.basis == "open surface (volume not trustworthy)"
    assert math.isnan(hydro.transverse_waterplane_inertia)
    assert math.isnan(hydro.longitudinal_waterplane_inertia)
