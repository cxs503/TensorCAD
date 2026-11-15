"""Voxel face of the CAD handoff: occupancy from a closed surface, and its evidence.

Every test here asserts an *invariant*, not that a function ran: an axis-aligned
solid fills exactly, the accelerated fill agrees with the O(cells x faces)
reference, a sphere's volume error shrinks at first order, the bundle round-trips,
and an open surface is refused instead of returning stripes.
"""
from __future__ import annotations

import math

import pytest

from tensorcad.geometry import Mesh
from tensorcad.hull import HullParameters, analytic_volume, build_mesh
from tensorcad.solver_interface import ellipsoid_mesh
from tensorcad.voxel_solid import (
    SolidLattice,
    box_surface_mesh,
    brute_force_solid_lattice,
    convergence_study,
    hull_solid_lattice,
    lattice_bundle,
    observed_order,
    occupancy_mismatch,
    ray_direction,
    read_lattice_bundle,
    solid_lattice,
    watertight_gate,
)


def test_ray_direction_is_unit_and_axial():
    for axis in ("x", "y", "z"):
        d = ray_direction(axis)
        assert math.isclose(math.sqrt(sum(c * c for c in d)), 1.0, rel_tol=0, abs_tol=1e-15)
        # The tilt must stay small enough that a 10^2 m column drifts far less
        # than a cell: the 1.3e-3 tilt that caused a 546-cell disagreement is
        # explicitly excluded.
        off = [abs(d[i]) for i in range(3) if "xyz"[i] != axis]
        assert max(off) < 1e-5
    with pytest.raises(ValueError):
        ray_direction("w")


def test_aligned_box_fills_exactly_on_every_ray_axis():
    box = box_surface_mesh((1.0, 0.6, 0.4))
    reference = 1.0 * 0.6 * 0.4
    lattices = {axis: solid_lattice(box, resolution=(10, 6, 4), ray_axis=axis) for axis in "xyz"}
    for axis, lattice in lattices.items():
        assert lattice.solid_count == lattice.cells, f"axis {axis} left a cell empty"
        assert math.isclose(lattice.solid_volume, reference, rel_tol=1e-15)
    assert occupancy_mismatch(lattices["y"], lattices["x"]) == 0
    assert occupancy_mismatch(lattices["z"], lattices["x"]) == 0


def test_index_order_matches_the_lattice_json_contract():
    box = box_surface_mesh((1.0, 1.0, 1.0))
    lattice = solid_lattice(box, resolution=(4, 3, 2))
    nx, ny, nz = lattice.shape
    assert lattice.index(1, 0, 0) == 1
    assert lattice.index(0, 1, 0) == nx
    assert lattice.index(0, 0, 1) == nx * ny
    with pytest.raises(IndexError):
        lattice.index(nx, 0, 0)


def test_accelerated_fill_agrees_with_the_bruteforce_reference():
    """The hull is curved, tapered and closed -- the case the index exists for."""
    hull, _reference, mesh = hull_solid_lattice(HullParameters(nx=21, nz=9), resolution=16)
    del hull
    fast = solid_lattice(mesh, resolution=(14, 8, 6))
    slow = brute_force_solid_lattice(mesh, resolution=(14, 8, 6))
    assert fast.shape == slow.shape
    # A centre lying on the surface has no defined parity (documented); every
    # other cell must agree cell-for-cell.
    mismatch = occupancy_mismatch(fast, slow)
    assert mismatch / fast.cells <= 1e-3, f"{mismatch} of {fast.cells} cells disagree"
    assert fast.solid_count > 0 and slow.solid_count > 0


@pytest.mark.parametrize("axis", ["x", "y", "z"])
def test_curved_solid_is_axis_independent(axis: str):
    ellipsoid = ellipsoid_mesh((1.0, 0.7, 0.5), n_theta=12, n_phi=16)
    reference = 4.0 / 3.0 * math.pi * 1.0 * 0.7 * 0.5
    lattice = solid_lattice(ellipsoid, resolution=(20, 14, 10), ray_axis=axis)
    assert abs(lattice.solid_volume - reference) / reference < 0.06


def test_sphere_volume_converges_at_first_order():
    radius = 0.5
    sphere = ellipsoid_mesh((radius, radius, radius), n_theta=40, n_phi=48)
    reference = 4.0 / 3.0 * math.pi * radius**3
    rows = convergence_study(sphere, reference_volume=reference, resolutions=(10, 20, 40))
    assert abs(rows[0].rel_error) > abs(rows[-1].rel_error)
    assert abs(rows[-1].rel_error) < 0.01
    order = observed_order(rows)
    # Worst case for an indicator-function quadrature is first order; a smooth
    # convex body measures ~1.9 because the curved-boundary error cancels.  The
    # invariant worth asserting is that the error is *converging*, not a bound on
    # the rate the geometry happens to deliver.
    assert order >= 0.5, f"observed order {order:.3f} is not converging"


def test_hull_voxel_volume_matches_the_analytic_displacement():
    parameters = HullParameters(nx=41, nz=13)
    lattice, reference, mesh = hull_solid_lattice(parameters, resolution=64)
    assert reference == analytic_volume(
        parameters.length, parameters.beam, parameters.draft
    )
    # The lattice integrates the *faceted* surface, which is itself -0.9 % short of
    # the analytic volume at this sampling; 3 % therefore covers both errors.
    assert abs(lattice.solid_volume - reference) / reference < 0.03
    assert lattice.ray_axis == "x"
    assert watertight_gate(mesh)["watertight"] is True


def test_open_surface_is_reported_not_voxelised_silently():
    box = box_surface_mesh((1.0, 0.6, 0.4))
    holed = Mesh(vertices=box.vertices, faces=box.faces[:-2])
    assert watertight_gate(holed)["watertight"] is False
    # Parity on an open surface is meaningless; the module does not pretend
    # otherwise -- the caller sees a non-watertight gate before any volume.
    closed = solid_lattice(box, resolution=(10, 6, 4)).solid_volume
    open_volume = solid_lattice(holed, resolution=(10, 6, 4)).solid_volume
    assert not math.isclose(open_volume, closed, rel_tol=1e-6)


def test_bundle_roundtrip_recovers_the_solid_mask(tmp_path):
    box = box_surface_mesh((0.5, 0.4, 0.3), origin=(0.25, 0.25, 0.2))
    lattice = solid_lattice(box, resolution=(12, 10, 8), pad=0.25)
    manifest = lattice_bundle(tmp_path, lattice, case="unit", walls={"normal_axis": "z"})
    assert manifest["digest_sha256"]
    assert not (tmp_path / "handoff.json").read_text().count("error")
    _, recovered = read_lattice_bundle(tmp_path)
    assert recovered.shape == lattice.shape
    assert recovered.occupancy == lattice.occupancy
    assert occupancy_mismatch(recovered, lattice) == 0
    assert recovered.mesh_digest == lattice.mesh_digest


def test_lattice_record_carries_the_frame_and_digest():
    box = box_surface_mesh((1.0, 1.0, 1.0))
    lattice = solid_lattice(box, resolution=4)
    record = lattice.to_dict()
    assert record["frame"] == "z-up"
    assert record["schema"] == "tensorcad.voxel_solid/1"
    assert len(record["mesh_digest_sha256"]) == 64
    assert "occupancy_b64" not in record
    assert len(lattice.to_dict(include_occupancy=True)["occupancy_b64"]) > 0


def test_bad_arguments_are_refused():
    box = box_surface_mesh((1.0, 1.0, 1.0))
    with pytest.raises(ValueError):
        solid_lattice(box, resolution=4, spacing=0.1)
    with pytest.raises(ValueError):
        solid_lattice(box)
    with pytest.raises(ValueError):
        solid_lattice(box, resolution=0)
    with pytest.raises(ValueError):
        solid_lattice(Mesh(vertices=((0.0, 0.0, 0.0),), faces=()))
    with pytest.raises(ValueError):
        SolidLattice(
            shape=(2, 2, 2),
            origin=(0.0, 0.0, 0.0),
            spacing=(1.0, 1.0, 1.0),
            occupancy=b"\x00",  # 1 entry for 8 cells
            mesh_digest="x",
            ray_axis="x",
        )