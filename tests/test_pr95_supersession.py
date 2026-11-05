"""PR #95 (``agent/m3923``, originally PR #40) cured: one integrator, its evidence kept.

The branch shipped ``tensorcad/righting_arm.py``, a complete heeled-displacement
integrator with a public surface (``RightingArmCurve``, ``righting_arm_curve``,
``wall_sided_righting_arm``) that collides name-for-name with ``tensorcad/gz.py``, which
landed on ``main`` first (PR #39 plus the m4029 union cure). Merging the branch as it
stood would have left **two parallel integrators for the same righting arm**, with the
package-level names silently rebound to whichever import ran last -- the same reason
PR #40 was closed instead of merged.

So the branch is cured by keeping the *engine* that is on ``main`` and porting the
branch's verification claims onto it, which is what these tests pin:

* there is exactly one integrator, and the public names belong to it (no duplicate
  module left behind to shadow it later);
* the branch's **sphere reference**: a floating sphere's transverse metacentre is its
  centre, so ``GZ(phi) = (R - KG) sin(phi)`` at every heel, and the meshed sphere
  converges to it as the sampling is refined (the branch measured worst < 1 mm);
* the branch's **constant-displacement claim** on both generated hull families
  (parametric and transom stern): every point of the curve is solved back onto the
  upright displacement, so the arm is a righting arm and not a curve of varying load.
"""

from __future__ import annotations

import importlib.util
import itertools
import math

import pytest

import tensorcad
from tensorcad.geometry import Mesh
from tensorcad.geometry_field import SphereField, to_mesh
from tensorcad.gz import heeled_hydrostatics, righting_arm_curve, with_freeboard
from tensorcad.hull import HullParameters, generate_hull
from tensorcad.transom import TransomHullParameters, generate_transom_hull

SPHERE_RADIUS = 1.0
SPHERE_KG = 0.5  # above the keel, i.e. below the section centre: GZ > 0 in every heel
SPHERE_HEELS = (0.0, 10.0, 20.0, 30.0, 45.0, 60.0)


def sphere_mesh(resolution: int) -> Mesh:
    """The branch's sphere built the way the surviving engine consumes it: a closed mesh."""
    sampled = to_mesh(SphereField(radius=SPHERE_RADIUS, center=(0.0, 0.0, 0.0)), resolution)
    return sampled.mesh


def test_only_the_surviving_integrator_exists() -> None:
    """The branch's parallel module is gone and the package names point at ``gz``."""
    assert importlib.util.find_spec("tensorcad.righting_arm") is None, (
        "tensorcad.righting_arm must not come back: it shadows tensorcad.gz's names"
    )
    for name in ("RightingArmCurve", "righting_arm_curve", "wall_sided_righting_arm"):
        obj = getattr(tensorcad, name)
        assert getattr(obj, "__module__", "tensorcad.gz") == "tensorcad.gz", (name, obj)


def test_sphere_reproduces_its_exact_righting_arm() -> None:
    """``GZ(phi) = (R - KG) sin(phi)`` -- exact for a sphere, meshing error only."""
    worst = {}
    for resolution in (24, 32):
        curve = righting_arm_curve(sphere_mesh(resolution), SPHERE_KG, SPHERE_HEELS)
        worst[resolution] = max(
            abs(gz - (SPHERE_RADIUS - SPHERE_KG) * math.sin(math.radians(heel)))
            for heel, gz in zip(curve.angles_deg, curve.gz, strict=True)
        )
        # The branch's own promise: worst |error| below 1 mm at resolution 32.
        assert worst[resolution] < 1e-3, (resolution, worst[resolution])
        assert curve.max_volume_error < 1e-12, curve.max_volume_error
    assert worst[32] < worst[24], worst  # refinement has to help, not wobble


def test_sphere_slope_at_zero_is_the_metacentric_height() -> None:
    """``dGZ/dphi`` at upright is ``R - KG``: the metacentre sits at the centre."""
    curve = righting_arm_curve(sphere_mesh(32), SPHERE_KG, (0.0, 5.0))
    slope = curve.gz[1] / math.sin(math.radians(5.0))
    assert slope == pytest.approx(SPHERE_RADIUS - SPHERE_KG, rel=5e-3), slope


@pytest.mark.parametrize(
    ("label", "mesh"),
    [
        ("parametric hull", generate_hull(HullParameters(length=10.0, beam=2.0, draft=0.8)).mesh),
        ("transom stern hull", generate_transom_hull(TransomHullParameters()).mesh),
    ],
)
def test_generated_hulls_hold_one_displacement_along_the_curve(label: str, mesh: Mesh) -> None:
    """The branch's cross-family claim: the arm is measured at constant displacement."""
    banded = with_freeboard(mesh, 0.25)  # a hull modelled to its waterline cannot be heeled
    upright = heeled_hydrostatics(banded, 0.0).displacement_volume
    plain = heeled_hydrostatics(mesh, 0.0).displacement_volume
    assert upright == pytest.approx(plain, rel=1e-12), (label, upright, plain)

    curve = righting_arm_curve(banded, kg=0.5, angles_deg=(0.0, 10.0, 20.0, 30.0))
    assert curve.max_volume_error < 1e-9, (label, curve.max_volume_error)
    assert abs(curve.gz[0]) < 1e-4, (label, curve.gz[0])  # upright: KN(0) is zero
    assert all(later > earlier for earlier, later in itertools.pairwise(curve.gz)), (
        label,
        curve.gz,
    )
