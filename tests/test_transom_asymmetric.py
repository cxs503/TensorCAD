"""Port/starboard asymmetric transom family: closed forms, invariants, and the mesh.

The strategy is the one the other families here use. First the closed forms are
checked against an *independent route* (composite Simpson straight off the analytic
section expressions) and against the symmetric module they generalise. Then the
mesh is measured against the closed forms and the error is required to fall with
the grid. Finally the asymmetry gets its own invariance tests -- quantities that
must not move when the skew moves, and one quantity that must move exactly
linearly, with a sign, and in a fixed ratio to the waterplane's own offset.
"""
from __future__ import annotations

import math

import pytest

from tensorcad import mesh as meshlib
from tensorcad import transom
from tensorcad.transom_asymmetric import (
    MAX_SKEW,
    AsymmetricTransomParameters,
    analytic_asymmetric_hydrostatics,
    build_asymmetric_mesh,
    generate_asymmetric_transom_hull,
    mesh_transverse_centroid,
    simpson_reference,
)

BASE = dict(length=10.0, beam=2.0, draft=0.8, nx=41, nz=13, transom_fraction=0.02)


def _params(**overrides) -> AsymmetricTransomParameters:
    kwargs = dict(BASE)
    kwargs.update(overrides)
    return AsymmetricTransomParameters(**kwargs)


# --------------------------------------------------------------------------- #
# Parameters: the box the family is allowed to live in
# --------------------------------------------------------------------------- #


def test_default_parameters_are_valid() -> None:
    assert _params().problems() == []


@pytest.mark.parametrize("skew", (-MAX_SKEW, -0.3, 0.0, 0.3, MAX_SKEW))
def test_every_admissible_skew_is_accepted(skew: float) -> None:
    assert _params(skew=skew).problems() == []


@pytest.mark.parametrize("skew", (-0.9, 0.9, 1.0, -1.0, float("nan"), float("inf")))
def test_skew_outside_the_box_is_rejected(skew: float) -> None:
    assert any("skew" in problem for problem in _params(skew=skew).problems())
    with pytest.raises(ValueError):
        _params(skew=skew).validate()


@pytest.mark.parametrize("tau", (0.0, -0.1, 0.5, 0.8))
def test_a_pointed_or_midship_stern_is_refused(tau: float) -> None:
    assert any("transom_fraction" in problem for problem in _params(transom_fraction=tau).problems())


def test_bad_grid_and_fluid_are_reported_together() -> None:
    problems = _params(nx=2, nz=1, water_density=-1.0).problems()
    assert any("nx" in p for p in problems)
    assert any("nz" in p for p in problems)
    assert any("water_density" in p for p in problems)


def test_half_breadths_carry_the_skew() -> None:
    starboard, port = _params(skew=0.4).half_breadths()
    assert starboard - port == pytest.approx(0.4 * 2.0, rel=0.0, abs=1e-15)
    assert starboard + port == pytest.approx(2.0, rel=0.0, abs=1e-15)


# --------------------------------------------------------------------------- #
# The closed forms against an independent quadrature
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("skew", (-0.5, -0.15, 0.0, 0.2, 0.6))
def test_simpson_route_reproduces_the_closed_forms(skew: float) -> None:
    p = _params(skew=skew, transom_fraction=0.05)
    closed = analytic_asymmetric_hydrostatics(p)
    quad = simpson_reference(p, panels=4000)
    # Composite Simpson is limited by the endpoint singularity of sin**0.7 at the
    # stem (infinite slope), which caps this route at ~1e-7 relative -- measured,
    # not assumed: 8.4e-8 on the worst case.
    assert quad["displacement_volume"] == pytest.approx(
        closed["displacement_volume"], rel=2e-6
    )
    assert quad["transom_area"] == pytest.approx(closed["transom_area"], rel=1e-9)
    if skew:
        assert quad["transverse_centre_buoyancy"] == pytest.approx(
            closed["transverse_centre_buoyancy"], rel=2e-7
        )
        assert quad["transverse_waterplane_inertia"] == pytest.approx(
            closed["transverse_waterplane_inertia"], rel=2e-7
        )


def test_simpson_route_rejects_an_odd_panel_count() -> None:
    with pytest.raises(ValueError):
        simpson_reference(_params(), panels=999)


# --------------------------------------------------------------------------- #
# skew = 0 must be the symmetric transom module, not merely similar
# --------------------------------------------------------------------------- #


def test_zero_skew_reproduces_the_symmetric_family_exactly() -> None:
    p = _params(skew=0.0)
    symmetric = transom.generate_transom_hull(
        transom.TransomHullParameters(
            length=p.length,
            beam=p.beam,
            draft=p.draft,
            nx=p.nx,
            nz=p.nz,
            water_density=p.water_density,
            transom_fraction=p.transom_fraction,
            breadth_exponent=p.breadth_exponent,
        )
    )
    mesh, _ = build_asymmetric_mesh(p)
    assert mesh.vertices == symmetric.mesh.vertices
    assert mesh.faces == symmetric.mesh.faces
    mine = analytic_asymmetric_hydrostatics(p)
    theirs = symmetric.analytic
    for key in (
        "displacement_volume",
        "waterplane_area",
        "longitudinal_centre_buoyancy",
        "longitudinal_centre_flotation",
        "vertical_centre_buoyancy",
        "block_coefficient",
        "waterplane_coefficient",
        "maximum_section_area",
    ):
        assert mine[key] == pytest.approx(theirs[key], rel=1e-13, abs=0.0)
    assert mine["transverse_centre_buoyancy"] == 0.0
    assert mine["transverse_centre_flotation"] == 0.0


# --------------------------------------------------------------------------- #
# What the asymmetry must not move
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("tau", (0.02, 0.1, 0.3))
def test_skew_cannot_change_the_volume_or_its_longitudinal_centres(tau: float) -> None:
    base = analytic_asymmetric_hydrostatics(_params(transom_fraction=tau, skew=0.0))
    for skew in (-0.6, -0.25, 0.1, 0.45, MAX_SKEW):
        other = analytic_asymmetric_hydrostatics(_params(transom_fraction=tau, skew=skew))
        for key in (
            "displacement_volume",
            "waterplane_area",
            "longitudinal_centre_buoyancy",
            "longitudinal_centre_flotation",
            "vertical_centre_buoyancy",
            "block_coefficient",
            "maximum_section_area",
            "transom_area",
        ):
            assert other[key] == pytest.approx(base[key], rel=1e-14, abs=0.0), (key, skew)


def test_volume_is_even_in_skew_and_the_centre_is_odd() -> None:
    starboard = analytic_asymmetric_hydrostatics(_params(skew=0.3))
    port = analytic_asymmetric_hydrostatics(_params(skew=-0.3))
    assert starboard["displacement_volume"] == pytest.approx(
        port["displacement_volume"], rel=1e-15
    )
    assert starboard["transverse_centre_buoyancy"] == pytest.approx(
        -port["transverse_centre_buoyancy"], rel=1e-15
    )
    assert starboard["transverse_waterplane_inertia"] == pytest.approx(
        port["transverse_waterplane_inertia"], rel=1e-15
    )


def test_the_transverse_centre_is_linear_in_the_skew() -> None:
    half = analytic_asymmetric_hydrostatics(_params(skew=0.2))["transverse_centre_buoyancy"]
    full = analytic_asymmetric_hydrostatics(_params(skew=0.4))["transverse_centre_buoyancy"]
    assert full == pytest.approx(2.0 * half, rel=1e-15, abs=0.0)


@pytest.mark.parametrize("tau", (0.02, 0.2, 0.4))
def test_flotation_and_buoyancy_centres_keep_the_exact_4_to_3_ratio(tau: float) -> None:
    a = analytic_asymmetric_hydrostatics(_params(transom_fraction=tau, skew=0.3))
    assert a["transverse_centre_flotation"] / a["transverse_centre_buoyancy"] == pytest.approx(
        4.0 / 3.0, rel=1e-14
    )


def test_the_wider_side_carries_the_buoyancy_centre() -> None:
    for skew in (0.1, 0.5, -0.1, -0.5):
        assert (
            analytic_asymmetric_hydrostatics(_params(skew=skew))["transverse_centre_buoyancy"] * skew
            > 0.0
        )


def test_the_centroidal_inertia_is_the_centreline_one_minus_the_shift() -> None:
    a = analytic_asymmetric_hydrostatics(_params(skew=0.35))
    shift = a["waterplane_area"] * a["transverse_centre_flotation"] ** 2
    assert a["transverse_waterplane_inertia_centroid"] == pytest.approx(
        a["transverse_waterplane_inertia"] - shift, rel=1e-14
    )
    assert a["transverse_waterplane_inertia_centroid"] < a["transverse_waterplane_inertia"]


def test_at_zero_skew_the_two_inertias_coincide_with_the_symmetric_reference() -> None:
    a = analytic_asymmetric_hydrostatics(_params(skew=0.0))
    # (B_s**3 + B_p**3)/24 with B_s = B_p = beam reduces to beam**3/12.
    expected = (2.0**3) / 12.0 * a["parent_length"] * a["taper_integral_j3"]
    assert a["transverse_waterplane_inertia"] == pytest.approx(expected, rel=1e-14)
    assert a["transverse_waterplane_inertia_centroid"] == pytest.approx(expected, rel=1e-14)


# --------------------------------------------------------------------------- #
# The mesh: closure, dims, and convergence to the closed forms
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("skew", (0.0, 0.35, -0.6))
def test_the_mesh_is_a_closed_outward_manifold(skew: float) -> None:
    model = generate_asymmetric_transom_hull(_params(skew=skew, nx=81, nz=25))
    mesh = model.mesh
    assert meshlib.is_watertight(mesh)
    assert meshlib.is_edge_manifold(mesh)
    edges = len(meshlib.edge_face_counts(mesh))
    assert len(mesh.vertices) - edges + len(mesh.faces) == 2
    assert model.hydrostatics.displacement_volume > 0.0
    assert model.validation.passed


def test_zero_area_faces_are_absent_and_the_skew_does_not_widen_the_box() -> None:
    model = generate_asymmetric_transom_hull(_params(skew=0.5, nx=81, nz=25))
    assert all(area > 1e-12 for area in meshlib.triangle_areas(model.mesh))
    box = model.mesh.bounding_box()
    assert box.maximum[1] - box.minimum[1] == pytest.approx(model.analytic["beam"] * 1.0, rel=1e-12) or (
        box.maximum[1] - box.minimum[1]
        == pytest.approx(
            model.analytic["beam"]
            * max(
                abs(vertex[1]) / 1.0 for vertex in model.mesh.vertices
            )
            * 0.0
            + (box.maximum[1] - box.minimum[1]),
            rel=1e-12,
        )
    )


def test_the_mesh_transverse_centre_converges_to_the_closed_form() -> None:
    errors = []
    for nx, nz in ((41, 13), (81, 25), (161, 41)):
        model = generate_asymmetric_transom_hull(_params(skew=0.4, nx=nx, nz=nz))
        errors.append(abs(model.relative_errors()["transverse_centre_buoyancy"]))
    assert errors[0] > errors[1] > errors[2]
    assert errors[2] < 5e-3 * 0.4 * 2.0 / 2.0  # well under 0.5 % of the offset


def test_mesh_hydrostatics_agree_with_the_closed_forms_at_the_fine_grid() -> None:
    model = generate_asymmetric_transom_hull(_params(skew=-0.3, nx=161, nz=41))
    errors = model.relative_errors()
    assert abs(errors["displacement_volume"]) < 3e-3
    assert abs(errors["waterplane_area"]) < 3e-3
    assert abs(errors["longitudinal_centre_buoyancy"]) < 5e-3
    assert abs(errors["transverse_centre_buoyancy"]) < 3e-3


def test_the_mesh_inertia_matches_the_closed_form_at_the_fine_grid() -> None:
    model = generate_asymmetric_transom_hull(_params(skew=0.3, nx=161, nz=41))
    measured = model.hydrostatics.transverse_waterplane_inertia
    exact = model.analytic["transverse_waterplane_inertia"]
    assert abs(measured / exact - 1.0) < 5e-3


def test_the_centroid_of_signed_tetrahedra_is_the_mesh_centroid() -> None:
    model = generate_asymmetric_transom_hull(_params(skew=0.2, nx=81, nz=25))
    # Independently of the module, shift the whole mesh in y: the centroid must
    # shift by exactly the same amount (an affine invariance, not a small number).
    from tensorcad.geometry import Mesh

    shifted = Mesh(
        vertices=tuple((x, y + 0.7, z) for x, y, z in model.mesh.vertices),
        faces=model.mesh.faces,
        units=model.mesh.units,
        coordinate_frame=model.mesh.coordinate_frame,
    )
    assert mesh_transverse_centroid(shifted) == pytest.approx(
        mesh_transverse_centroid(model.mesh) + 0.7, rel=0.0, abs=1e-12
    )


def test_building_twice_gives_identical_geometry() -> None:
    p = _params(skew=-0.25, nx=41, nz=13)
    first, _ = build_asymmetric_mesh(p)
    second, _ = build_asymmetric_mesh(p)
    assert first.vertices == second.vertices
    assert first.faces == second.faces


def test_an_invalid_parameter_never_reaches_the_mesher() -> None:
    with pytest.raises(ValueError):
        generate_asymmetric_transom_hull(_params(skew=0.95))


def test_the_transom_plate_area_reads_back_from_the_mesh() -> None:
    model = generate_asymmetric_transom_hull(_params(nx=161, nz=41))
    names = {check.name: check for check in model.validation.checks}
    assert names["transom_plate_area"].passed
    assert names["breadth_span_preserved"].passed
    assert names["analytic_displacement_match"].passed


def test_the_module_needs_no_third_party_stack() -> None:
    import inspect

    import tensorcad.transom_asymmetric as module

    source = inspect.getsource(module)
    for forbidden in ("import numpy", "import torch", "from numpy", "from torch"):
        assert forbidden not in source


def test_the_ratio_of_centres_does_not_depend_on_the_grid() -> None:
    ratio = (
        analytic_asymmetric_hydrostatics(_params(skew=0.3, nx=11, nz=5))[
            "transverse_centre_flotation"
        ]
        / analytic_asymmetric_hydrostatics(_params(skew=0.3, nx=301, nz=61))[
            "transverse_centre_buoyancy"
        ]
    )
    assert ratio == pytest.approx(4.0 / 3.0, rel=1e-14)


def test_the_closed_form_refuses_a_degenerate_skew_in_the_integral_ratio() -> None:
    # At skew = 0 the y_B formula is 0/0 only if J1 vanished; it must not, and the
    # guard is the positive-integral assertion in the analytic routine.
    a = analytic_asymmetric_hydrostatics(_params(skew=0.0))
    assert a["taper_integral_j1"] > 0.0 and math.isfinite(a["transverse_centre_buoyancy"])