"""Pins for the cure of PR #44: the union of ``main`` and ``agent/m3347``.

The merge of ``main`` (waterplane second moments, asymmetric taper, transom and Wigley
families) with ``agent/m3347`` (finite-heel righting arms) touched a shared export list
and a shared CHANGELOG section. These tests exist so the resolution cannot rot back into
a silent loss of one side:

* every symbol either side added is exported exactly once;
* the branch's route still holds the displacement fixed and still puts ``GZ(0) == 0``;
* the branch's upright ``GM`` is the arithmetic ``initial_stability`` reports, i.e. the
  two integrations were not accidentally mixed;
* ``main``'s waterplane-inertia and transom routes still answer their own references.

The full numerical evidence (181 heel angles against an independent 2-D section cut, a
grid ladder, and the two hull families) lives in ``tools/verify_merge_m3347.py``; these
are the cheap invariants that would break first if the merge were redone differently.
"""

from __future__ import annotations

import math

import pytest

import tensorcad as tc
from tensorcad import transom
from tensorcad.geometry import MeshBuilder
from tensorcad.hull import HullParameters, generate_hull
from tensorcad.hydrostatics import compute_hydrostatics
from tensorcad.stability import initial_stability, righting_arm_curve

# The symbols each side introduced, spelled out so a dropped export fails here.
BRANCH_SYMBOLS = (
    "GZCurve",
    "GZPoint",
    "SCHEMA_GZ",
    "DEFAULT_HEEL_ANGLES",
    "righting_arm_curve",
)
MAIN_SYMBOLS = (
    "MAX_TRANSOM_FRACTION",
    "TransomHullParameters",
    "TransomModel",
    "analytic_transom_hydrostatics",
    "build_transom_mesh",
    "estimate_transom_discretisation_error",
    "generate_transom_hull",
    "parent_length",
    "sinc_power_coefficients",
    "transom_area_ratio",
    "transom_fraction_for_area_ratio",
    "truncated_taper_integral",
    "truncated_taper_moment",
    "waterplane_inertia_longitudinal",
    "waterplane_inertia_transverse",
)


def box_mesh(length: float, beam: float, draft: float, freeboard: float):
    """Closed box: keel on ``z = -draft``, deck on ``z = +freeboard``, waterline ``z = 0``."""
    builder = MeshBuilder()
    corner = {}
    for x in (0.0, length):
        for y in (-0.5 * beam, 0.5 * beam):
            for z in (-draft, freeboard):
                corner[(x, y, z)] = builder.add_vertex((x, y, z))

    def quad(p0, p1, p2, p3) -> None:
        a, b, c, d = (corner[p] for p in (p0, p1, p2, p3))
        builder.add_triangle(a, b, c)
        builder.add_triangle(a, c, d)

    lo, hi = -0.5 * beam, 0.5 * beam
    bottom, top = -draft, freeboard
    quad((0.0, lo, bottom), (0.0, hi, bottom), (length, hi, bottom), (length, lo, bottom))
    quad((0.0, lo, top), (length, lo, top), (length, hi, top), (0.0, hi, top))
    quad((0.0, lo, bottom), (length, lo, bottom), (length, lo, top), (0.0, lo, top))
    quad((0.0, hi, bottom), (0.0, hi, top), (length, hi, top), (length, hi, bottom))
    quad((0.0, lo, bottom), (0.0, lo, top), (0.0, hi, top), (0.0, hi, bottom))
    quad((length, lo, bottom), (length, hi, bottom), (length, hi, top), (length, lo, top))
    return builder.build()


@pytest.fixture(scope="module")
def box_curve():
    return righting_arm_curve(
        box_mesh(10.0, 2.0, 0.5, 1.5), kg=0.45, waterline_z=0.0, angles_deg=range(0, 181, 5)
    )


def test_both_sides_of_the_merge_are_exported_once() -> None:
    exported = list(tc.__all__)
    for name in BRANCH_SYMBOLS + MAIN_SYMBOLS:
        assert name in exported, f"{name} was lost by the merge"
        assert exported.count(name) == 1, f"{name} is exported more than once"
        assert hasattr(tc, name), f"{name} is declared but not importable"


def test_the_righting_arm_route_holds_the_displacement_fixed(box_curve) -> None:
    volumes = [point.displacement_volume for point in box_curve.points]
    assert max(abs(volume - 10.0) for volume in volumes) < 1e-9
    assert box_curve.points[0].gz == pytest.approx(0.0, abs=1e-12)


def test_the_righting_arm_recovers_the_wall_sided_formula_inside_its_window(box_curve) -> None:
    # tan(phi) < 2T/B = 0.5 keeps both waterline ends on the vertical sides.
    bm = 2.0**2 / (12.0 * 0.5)
    bg = 0.45 - 0.5 * 0.5
    gz_15 = next(point.gz for point in box_curve.points if point.heel_deg == 15)
    phi = math.radians(15.0)
    expected = math.sin(phi) * (bm + 0.5 * bm * math.tan(phi) ** 2 - bg)
    assert gz_15 == pytest.approx(expected, rel=1e-12)


def test_the_two_upright_integrations_still_agree() -> None:
    mesh = generate_hull(HullParameters(length=10.0, beam=2.0, draft=0.8)).mesh
    upright = initial_stability(mesh, kg=0.6)
    curve = righting_arm_curve(mesh, kg=0.6, waterline_z=0.0, freeboard=0.4, angles_deg=range(0, 91, 10))
    assert curve.initial_gm_transverse == pytest.approx(upright.gm_transverse, rel=1e-12)


def test_the_waterplane_inertia_route_still_answers_its_reference() -> None:
    barge = compute_hydrostatics(box_mesh(10.0, 2.0, 0.8, 0.0))
    assert barge.transverse_waterplane_inertia == pytest.approx(10.0 * 2.0**3 / 12.0, rel=1e-12)
    assert barge.longitudinal_waterplane_inertia == pytest.approx(2.0 * 10.0**3 / 12.0, rel=1e-12)
    assert barge.keel_z == pytest.approx(-0.8, abs=1e-12)


def test_the_transom_family_still_answers_its_reference() -> None:
    params = transom.TransomHullParameters()
    model = transom.generate_transom_hull(params)
    mesh = transom.compute_hydrostatics(model.mesh)
    reference = transom.analytic_transom_hydrostatics(params)
    assert mesh.block_coefficient == pytest.approx(reference["block_coefficient"], rel=1e-2)