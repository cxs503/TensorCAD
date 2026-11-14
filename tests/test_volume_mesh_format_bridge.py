"""Bridge from a body-fitted hull grid (tensorcad.volume_mesh) to the volume formats.

Why this file exists
====================
``meshio_exchange`` could send a *star-tessellated surface* to SU2/CGNS, but star
tessellation is refused for the package's own hull generator (the rings collapse to zero
breadth at the stem and post, so no single interior apex sees the whole shell).  The
documented remedy -- "use :mod:`tensorcad.volume_mesh` and feed the tetrahedra instead" --
had no implementation: nothing turned a :class:`VolumeMesh` into a
:class:`TetrahedralMesh`.  That is the gap these tests pin down.

The invariants asserted here are identities, not smoke tests
============================================================
* **Split conservation**: the six Kuhn tetrahedra of a hexahedron tile it without overlap,
  so ``sum(tet volumes) == volume_mesh.volume()`` to round-off (1e-12 relative).  Nothing
  is allowed to be lost between the hex grid and the tetrahedra.
* **Analytic agreement**: for a box the grid volume must equal ``L * B * H`` exactly, and
  for a ``SectionGrid`` it must equal ``2 * integral(b dx dz)`` (trapezoid, exact because
  both the cells and the integration are linear).  These are closed-form anchors, not
  self-comparison.
* **Orientation**: every kept tetrahedron has strictly positive signed volume, and the
  zero-volume ones are *counted* in ``degenerate_cells`` (collapsed keel/stem columns).
* **Round trip**: SU2 and CGNS preserve the volume through meshio's own writer/reader.
"""
from __future__ import annotations

import dataclasses
import math

import pytest

from tensorcad import volume_mesh as vmlib
from tensorcad.meshio_exchange import (
    FORMATS_VOLUME,
    MeshioFormatError,
    TetrahedralMesh,
    round_trip_meshio,
    volume_mesh_to_tetrahedral,
)

pytest.importorskip("meshio")

REL = 1e-12


# ---------------------------------------------------------------------------------------
# Fixtures: one analytic box, one real hull body
# ---------------------------------------------------------------------------------------

def box_grid(nx: int = 4, nz: int = 3, b: float = 0.5) -> vmlib.SectionGrid:
    """A rectangular box of half breadth ``b``: volume is exactly ``2b * L * H``."""
    xs = tuple(0.0 + i * (2.0 / (nx - 1)) for i in range(nx))
    zs = tuple(0.0 + j * (0.6 / (nz - 1)) for j in range(nz))
    return vmlib.SectionGrid(xs=xs, zs=zs, breadth=lambda i, j: b)


def hull_grid() -> vmlib.SectionGrid:
    """A real hull body: Wigley-type half breadth over a station/level grid.

    Analytic half breadth ``b(x, z) = (B/2)(1 - (2x/L)^2)(1 - (z/T)^2)`` vanishes at both
    ends and at the keel, so this fixture exercises the collapsed columns that make the
    surface-tessellation path fail on real hulls.
    """
    length, beam, draft = 10.0, 1.0, 0.625
    xs = tuple(-length / 2 + i * (length / 24) for i in range(25))
    zs = tuple(j * (draft / 8) for j in range(9))

    def breadth(i: int, j: int) -> float:
        x, z = xs[i], zs[j]
        value = (beam / 2) * (1 - (2 * x / length) ** 2) * (1 - (z / draft) ** 2)
        return max(0.0, value)

    return vmlib.SectionGrid(xs=xs, zs=zs, breadth=breadth)


def trapezoid_volume(grid: vmlib.SectionGrid) -> float:
    """``2 * integral(b)`` over the (x, z) box, exact for the piecewise-linear cells."""
    total = 0.0
    for i in range(grid.nx - 1):
        dx = grid.xs[i + 1] - grid.xs[i]
        for j in range(grid.nz - 1):
            dz = grid.zs[j + 1] - grid.zs[j]
            corners = (
                grid.breadth(i, j),
                grid.breadth(i + 1, j),
                grid.breadth(i, j + 1),
                grid.breadth(i + 1, j + 1),
            )
            total += dx * dz * sum(corners) / 4.0
    return 2.0 * total


@pytest.fixture(scope="module")
def box_mesh() -> vmlib.VolumeMesh:
    return vmlib.build_volume_mesh(box_grid(), first_layer=0.02, growth=1.2, layers=3)


@pytest.fixture(scope="module")
def hull_mesh() -> vmlib.VolumeMesh:
    return vmlib.build_volume_mesh(hull_grid(), first_layer=0.004, growth=1.25, layers=4)


# ---------------------------------------------------------------------------------------
# The split itself
# ---------------------------------------------------------------------------------------

def test_box_split_conserves_volume_and_counts() -> None:
    """6 tetrahedra per hexahedron, and their volumes add up to the box volume."""
    vm = vmlib.build_volume_mesh(box_grid(), first_layer=0.02, growth=1.2, layers=3)
    assert len(vm.tets) == 6 * len(vm.hexes)
    tet = volume_mesh_to_tetrahedral(vm)
    assert tet.degenerate_cells == 0
    assert len(tet.tets) == 6 * len(vm.hexes)
    assert tet.volume() == pytest.approx(vm.volume(), rel=REL)


def test_box_volume_matches_analytic() -> None:
    """A box has a closed-form volume: 2b * L * H, checked through the tetrahedra."""
    vm = vmlib.build_volume_mesh(box_grid(nx=5, nz=4, b=0.5), first_layer=0.02, layers=3)
    tet = volume_mesh_to_tetrahedral(vm)
    analytic = 2 * 0.5 * 2.0 * 0.6
    assert tet.volume() == pytest.approx(analytic, rel=REL)
    assert tet.volume() == pytest.approx(vm.volume(), rel=REL)


def test_grid_volume_matches_trapezoid_integral() -> None:
    """Anchor against an independent quadrature of the same breadth field."""
    vm = vmlib.build_volume_mesh(hull_grid(), first_layer=0.004, growth=1.25, layers=4)
    tet = volume_mesh_to_tetrahedral(vm)
    assert tet.volume() == pytest.approx(trapezoid_volume(hull_grid()), rel=REL)
    assert tet.volume() == pytest.approx(vm.volume(), rel=REL)


def test_hull_split_is_certified_and_zero_cells_counted() -> None:
    """Real hull: collapsed columns produce zero-volume cells, which are counted."""
    vm = vmlib.build_volume_mesh(hull_grid(), first_layer=0.004, growth=1.25, layers=4)
    tet = volume_mesh_to_tetrahedral(vm)
    assert tet.degenerate_cells > 0, "collapsed keel/stem columns must be detected"
    assert tet.degenerate_cells + len(tet.tets) == len(vm.tets)
    assert tet.volume() == pytest.approx(vm.volume(), rel=REL)


def test_all_kept_tetrahedra_strictly_positive() -> None:
    """Orientation is an invariant: no kept cell may be flat or inverted."""
    vm = vmlib.build_volume_mesh(hull_grid(), first_layer=0.004, growth=1.25, layers=4)
    tet = volume_mesh_to_tetrahedral(vm)
    from tensorcad.meshio_exchange import tetrahedral_volume

    for index, cell in enumerate(tet.tets):
        assert tetrahedral_volume(tet.vertices, cell) > 0.0, f"cell {index} is not positive"
    assert tet.signed_volume() > 0.0


# ---------------------------------------------------------------------------------------
# Failure modes: refused, loudly, never silently
# ---------------------------------------------------------------------------------------

def test_unsplit_grid_is_refused_with_the_remedy() -> None:
    vm = vmlib.build_volume_mesh(box_grid(), first_layer=0.02, layers=2, split_tets=False)
    assert vm.tets == ()
    with pytest.raises(MeshioFormatError, match="split_tets=True"):
        volume_mesh_to_tetrahedral(vm)


def test_inverted_cell_is_refused_not_reordered() -> None:
    """A deliberately flipped tet must raise: reordering here would hide a builder bug."""
    vm = vmlib.build_volume_mesh(box_grid(), first_layer=0.02, layers=2)
    broken = list(vm.tets)
    a, b, c, d = broken[3]
    broken[3] = (a, c, b, d)  # swap two nodes -> negative signed volume
    tampered = dataclasses.replace(vm, tets=tuple(broken))
    with pytest.raises(MeshioFormatError, match="inverted"):
        volume_mesh_to_tetrahedral(tampered)


def test_surface_mesh_is_rejected_by_type() -> None:
    with pytest.raises(TypeError, match="tetrahedralise_surface"):
        volume_mesh_to_tetrahedral(object())  # type: ignore[arg-type]


def test_tolerance_is_scale_aware() -> None:
    """The default tolerance scales with the body, so a metre-sized hull still passes."""
    vm = vmlib.build_volume_mesh(hull_grid(), first_layer=0.004, layers=3)
    tet = volume_mesh_to_tetrahedral(vm)
    assert math.isfinite(tet.volume())
    # An absurdly tight explicit tolerance must not turn round-off into a false alarm.
    assert volume_mesh_to_tetrahedral(vm, tolerance=1e-12 * vm.volume()).volume() == pytest.approx(
        vm.volume(), rel=REL
    )


# ---------------------------------------------------------------------------------------
# Round trip through the volume formats
# ---------------------------------------------------------------------------------------

@pytest.mark.parametrize("fmt", FORMATS_VOLUME)
def test_round_trip_preserves_volume(fmt: str, hull_mesh: vmlib.VolumeMesh) -> None:
    tet = volume_mesh_to_tetrahedral(hull_mesh)
    back = round_trip_meshio(tet, file_format=fmt)
    assert isinstance(back, TetrahedralMesh)
    assert len(back.tets) == len(tet.tets)
    assert back.volume() == pytest.approx(tet.volume(), rel=1e-9)
    assert back.signed_volume() > 0.0


@pytest.mark.parametrize("fmt", FORMATS_VOLUME)
def test_round_trip_preserves_vertices_and_cells(fmt: str, box_mesh: vmlib.VolumeMesh) -> None:
    tet = volume_mesh_to_tetrahedral(box_mesh)
    back = round_trip_meshio(tet, file_format=fmt)
    assert len(back.vertices) == len(tet.vertices)
    for original, restored in zip(tet.vertices, back.vertices, strict=True):
        for a, b in zip(original, restored, strict=True):
            assert abs(a - b) < 1e-12
    assert set(back.tets) == set(tet.tets) or back.tets == tet.tets


def test_bridge_scales_with_a_finer_grid() -> None:
    """Refinement must not change the volume: same body, more cells."""
    coarse = vmlib.build_volume_mesh(box_grid(nx=4, nz=3), first_layer=0.02, layers=2)
    fine = vmlib.build_volume_mesh(box_grid(nx=9, nz=7), first_layer=0.01, layers=2)
    a = volume_mesh_to_tetrahedral(coarse)
    b = volume_mesh_to_tetrahedral(fine)
    assert len(b.tets) > len(a.tets)
    assert b.volume() == pytest.approx(a.volume(), rel=1e-9)
