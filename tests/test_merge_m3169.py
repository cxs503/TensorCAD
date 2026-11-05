"""PR #40's union: the righting-arm family and the transom/hydrostatics family together.

Both sides of the merge are alive on the merged tree, and — more interesting — they
compose. These three tests pin that composition, because a merge that only *compiles*
would let the next change break one side silently:

* the two modules must integrate the same geometry to the same displacement and the same
  buoyancy centre when the hull is upright (``compute_hydrostatics`` from ``main`` versus
  ``heeled_displacement`` from the righting-arm module);
* the heeled solver must accept ``main``'s transom mesh and hold the displacement it
  solved at every heel, with ``kg`` entering exactly as ``-sin(phi)``;
* the public surface must carry both families at once.

The known cross-module inconsistency on a *fine* hull (the small-angle slope misses the
``I_T/nabla`` term, exactly ``BM``) is deliberately **not** asserted away here: it is
measured and reported in ``artifacts/m3923/merge_verification.json`` and
``docs/m3923_pr40_merge_cure.md``.
"""
from __future__ import annotations

import math

import pytest

import tensorcad
from tensorcad import (
    HullParameters,
    TransomHullParameters,
    compute_hydrostatics,
    generate_hull,
    generate_transom_hull,
    righting_arm,
    righting_arm_curve,
)
from tensorcad.righting_arm import heeled_displacement

HEELS = (0.0, 5.0, 10.0)


def transom_mesh(tau: float = 0.2):  # type: ignore[no-untyped-def]
    return generate_transom_hull(TransomHullParameters(transom_fraction=tau)).mesh


def test_the_merged_public_surface_carries_both_families() -> None:
    names = set(tensorcad.__all__)
    # m3169's family
    assert {"righting_arm", "righting_arm_curve", "heeled_displacement",
            "wall_sided_righting_arm"} <= names
    # main's families
    assert {"generate_transom_hull", "analytic_transom_hydrostatics", "build_transom_mesh",
            "compute_hydrostatics", "transom_fraction_for_area_ratio"} <= names
    assert {"transverse_waterplane_inertia", "longitudinal_waterplane_inertia"}.isdisjoint(names)


@pytest.mark.parametrize("mesh_factory", [lambda: generate_hull(HullParameters()).mesh,
                                          transom_mesh])
def test_the_two_modules_agree_on_the_upright_displacement(mesh_factory) -> None:  # type: ignore[no-untyped-def]
    """Two independent integrations of the same geometry, one geometry."""
    mesh = mesh_factory()
    record = compute_hydrostatics(mesh, 0.0)
    state = heeled_displacement(mesh, heel_deg=0.0, waterline_z=0.0)
    assert state.volume == pytest.approx(record.displacement_volume, rel=1e-9)
    assert state.centre_of_buoyancy[2] == pytest.approx(record.vertical_centre_buoyancy,
                                                        abs=1e-9)


def test_righting_arm_on_a_transom_hull_holds_displacement_and_the_kg_lever() -> None:
    """``main``'s transom mesh, solved by m3169's heeled integrator."""
    mesh = transom_mesh(0.2)
    hydro = compute_hydrostatics(mesh, 0.0)
    curve = righting_arm_curve(mesh, kg=hydro.vertical_centre_buoyancy + 0.30, heels_deg=HEELS)
    for point in curve.points:
        assert abs(point.volume_residual) < 1e-9
        assert math.isfinite(point.righting_arm)
    # GZ(0) is the upright buoyancy centre's own y offset, which on this transom mesh
    # is 4.0e-5 m rather than 0: a y-symmetry property of the tessellation, not of the
    # merge. Pinned at 0.1 mm so a real asymmetry would fail.
    assert abs(curve.points[0].righting_arm) < 1e-4

    base = righting_arm(mesh, heel_deg=5.0, kg=0.0).righting_arm
    for kg in (0.3, 0.6):
        predicted = base - kg * math.sin(math.radians(5.0))
        assert righting_arm(mesh, heel_deg=5.0, kg=kg).righting_arm == pytest.approx(
            predicted, abs=1e-12
        )