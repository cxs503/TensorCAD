"""The section grid must agree with the hull's own ring points, level for level.

Why this file exists
--------------------
``SectionGrid.from_hull`` is the bridge between the lofted hull surface and the
body-fitted volume mesh every downstream solver consumes (``GAP-CAD-B3-SOLVER-INTERFACE``).
Its ``zs[j]``/``breadth(i, j)`` pairing is not observable in the displacement volume:
mirroring a body in ``z`` leaves the volume exactly unchanged, so a flipped pairing
passes every volume gate.  The vertical centre of buoyancy (VCB) is what sees it, and
the hull family carries an analytic VCB (``-2/5 * draft``), so the pairing is gated
against it here.

The defect this file guards: ``breadth(i, j)`` used to read ring point ``j``, whose
height is ``-draft * j / (nz - 1)``, while ``zs[j]`` was ``-draft * (1 - j / (nz - 1))``.
Ring index 0 is the design waterline, so the grid claimed ``b = B/2`` at the keel and a
closed waterline: the hull turned upside down, VCB off by 50.7 % of the analytic value,
while the volume stayed at 7.446554 m3.
"""

from __future__ import annotations

import math

import pytest

from tensorcad.hull import (
    HullParameters,
    _section_point,
    analytic_vertical_centre_buoyancy,
)
from tensorcad.volume_mesh import SectionGrid, box_section_grid, volume_mesh_from_hull

LENGTH, BEAM, DRAFT = 10.0, 2.0, 0.8


def hull_parameters(nx: int = 41, nz: int = 13) -> HullParameters:
    return HullParameters(length=LENGTH, beam=BEAM, draft=DRAFT, nx=nx, nz=nz)


def grid_vcb(grid: SectionGrid) -> tuple[float, float]:
    """Volume and VCB of the solid a section grid describes, by bilinear quadrature."""
    xs, zs = grid.xs, grid.zs
    volume = 0.0
    moment = 0.0
    for i in range(len(xs) - 1):
        dx = xs[i + 1] - xs[i]
        for j in range(len(zs) - 1):
            dz = zs[j + 1] - zs[j]
            b = sum(grid.breadth(i + a, j + c) for a in (0, 1) for c in (0, 1)) / 4.0
            volume += dx * dz * 2.0 * b
            moment += dx * dz * 2.0 * b * 0.5 * (zs[j] + zs[j + 1])
    return volume, moment / volume


def test_every_level_carries_the_ring_point_at_that_height():
    """``zs[j]`` must be the height of the ring point the breadth was read from."""
    parameters = hull_parameters()
    grid = SectionGrid.from_hull(parameters)
    for i in range(parameters.nx):
        for j in range(parameters.nz):
            ring = _section_point(parameters, i, parameters.nz - 1 - j)
            assert grid.zs[j] == ring[2]
            assert grid.breadth(i, j) == pytest.approx(abs(ring[1]), rel=0.0, abs=1e-15)


def test_waterline_is_wide_and_the_keel_closes():
    parameters = hull_parameters()
    grid = SectionGrid.from_hull(parameters)
    middle = (parameters.nx - 1) // 2
    assert grid.zs[0] == pytest.approx(-DRAFT, abs=1e-12)
    assert grid.zs[-1] == 0.0
    assert grid.breadth(middle, 0) == pytest.approx(0.0, abs=1e-15)  # keel
    assert grid.breadth(middle, parameters.nz - 1) == pytest.approx(BEAM / 2.0, rel=1e-12)


def test_vcb_matches_the_analytic_value():
    """The gate the volume gates cannot provide: vertical centre of buoyancy."""
    grid = SectionGrid.from_hull(hull_parameters())
    _, vcb = grid_vcb(grid)
    reference = analytic_vertical_centre_buoyancy(DRAFT)
    assert reference < 0.0
    assert abs(vcb / reference - 1.0) < 0.03, (vcb, reference)


def test_a_mirrored_pairing_is_caught_by_vcb_and_not_by_volume():
    """Negative control: the pre-fix pairing keeps the volume and breaks the VCB."""
    parameters = hull_parameters()
    correct = SectionGrid.from_hull(parameters)
    mirrored = SectionGrid(
        xs=correct.xs,
        zs=correct.zs,
        breadth=lambda i, j: abs(_section_point(parameters, i, j)[1]),
    )
    volume_correct, vcb_correct = grid_vcb(correct)
    volume_mirrored, vcb_mirrored = grid_vcb(mirrored)
    assert abs(volume_mirrored / volume_correct - 1.0) < 1e-12  # blind to the flip
    reference = analytic_vertical_centre_buoyancy(DRAFT)
    assert abs(vcb_mirrored / reference - 1.0) > 0.15  # and this is what sees it
    assert abs(vcb_correct / reference - 1.0) < 0.03


def test_volume_mesh_centroid_sits_below_the_waterline():
    """The FVM/FEM body mesh inherits the grid: its centroid must be below the waterline."""
    mesh = volume_mesh_from_hull(hull_parameters(), first_layer=0.005, layers=4)
    volume = 0.0
    moment = 0.0
    for tet in mesh.tets:
        points = [mesh.vertices[index] for index in tet]
        a, b, c, d = points
        u = tuple(b[k] - a[k] for k in range(3))
        v = tuple(c[k] - a[k] for k in range(3))
        w = tuple(d[k] - a[k] for k in range(3))
        det = (
            u[0] * (v[1] * w[2] - v[2] * w[1])
            - u[1] * (v[0] * w[2] - v[2] * w[0])
            + u[2] * (v[0] * w[1] - v[1] * w[0])
        )
        cell = abs(det) / 6.0
        if cell == 0.0:
            continue
        volume += cell
        moment += cell * sum(point[2] for point in points) / 4.0
    centroid = moment / volume
    reference = analytic_vertical_centre_buoyancy(DRAFT)
    assert reference < centroid < 0.0
    assert math.isclose(centroid, reference, rel_tol=0.05), (centroid, reference)


def test_the_box_reference_case_is_unchanged():
    """A box's breadth is constant, so the pairing has nothing to flip."""
    grid = box_section_grid(length=1.0, beam=1.0, draft=1.0, nx=3, nz=3)
    for j in range(grid.nz):
        assert grid.breadth(1, j) == pytest.approx(0.5)