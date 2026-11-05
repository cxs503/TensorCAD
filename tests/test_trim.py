"""Trim and the centimetre constants: closed forms, two routes, and the solves."""
from __future__ import annotations

import math

import pytest

from tensorcad import HullParameters, analytic_volume, generate_hull, generate_wigley
from tensorcad.bonjean import section_properties
from tensorcad.geometry import Mesh
from tensorcad.trim import (
    SEAWATER_DENSITY,
    TrimError,
    _SectionCutter,
    _stations,
    box_analytic_particulars,
    box_trim_tangent,
    equilibrium_trim,
    hydrostatic_particulars,
    hydrostatic_table,
    small_angle_trim,
    trimmed_hydrostatics,
    waterline_for_displacement,
)

BOX = (20.0, 4.0, 5.0)


def box_mesh(length: float, beam: float, draft: float, freeboard: float = 2.0) -> Mesh:
    """A watertight rectangular hull; the one geometry in the package with closed forms."""
    xs = [0.0, length]
    ys = [-0.5 * beam, 0.5 * beam]
    zs = [-draft, freeboard]
    vertices: list[tuple[float, float, float]] = []
    index: dict[tuple[float, float, float], int] = {}
    faces: list[tuple[int, int, int]] = []

    def vid(point):
        if point not in index:
            index[point] = len(vertices)
            vertices.append(point)
        return index[point]

    def quad(a, b, c, d):
        faces.append((vid(a), vid(b), vid(c)))
        faces.append((vid(a), vid(c), vid(d)))

    quad((xs[0], ys[0], zs[0]), (xs[0], ys[1], zs[0]), (xs[1], ys[1], zs[0]), (xs[1], ys[0], zs[0]))
    quad((xs[0], ys[0], zs[1]), (xs[1], ys[0], zs[1]), (xs[1], ys[1], zs[1]), (xs[0], ys[1], zs[1]))
    quad((xs[0], ys[0], zs[0]), (xs[1], ys[0], zs[0]), (xs[1], ys[0], zs[1]), (xs[0], ys[0], zs[1]))
    quad((xs[1], ys[1], zs[0]), (xs[0], ys[1], zs[0]), (xs[0], ys[1], zs[1]), (xs[1], ys[1], zs[1]))
    quad((xs[0], ys[1], zs[0]), (xs[0], ys[0], zs[0]), (xs[0], ys[0], zs[1]), (xs[0], ys[1], zs[1]))
    quad((xs[1], ys[0], zs[0]), (xs[1], ys[1], zs[0]), (xs[1], ys[1], zs[1]), (xs[1], ys[0], zs[1]))
    return Mesh(vertices=tuple(vertices), faces=tuple(faces))


@pytest.fixture(scope="module")
def box() -> Mesh:
    return box_mesh(*BOX)


@pytest.fixture(scope="module")
def hull():
    return generate_hull(HullParameters(length=10.0, beam=2.0, draft=0.8))


@pytest.fixture(scope="module")
def box_particulars(box):
    return hydrostatic_particulars(box, 0.0, station_count=41)


@pytest.mark.parametrize(
    "field",
    [
        "volume",
        "lcb",
        "lcf",
        "kb",
        "waterplane_area",
        "transverse_second_moment",
        "longitudinal_second_moment",
        "tpc",
        "mct",
    ],
)
def test_box_particulars_reproduce_the_closed_forms(box_particulars, field):
    """A box is the one hull whose particulars are elementary integrals, so they are exact."""
    exact = box_analytic_particulars(*BOX)
    measured = getattr(box_particulars, field)
    reference = getattr(exact, field)
    assert measured == pytest.approx(reference, rel=1e-12), field


def test_the_centimetre_constants_have_the_units_the_drawing_office_uses(box_particulars):
    exact = box_analytic_particulars(*BOX)
    assert box_particulars.tpc == pytest.approx(1025.0 * BOX[0] * BOX[1] / 1e5, rel=1e-12)
    drawing_office_mct = 1.025 * BOX[1] * BOX[0] ** 3 / 12.0 / (100.0 * BOX[0])
    assert box_particulars.mct == pytest.approx(drawing_office_mct, rel=1e-12)
    assert box_particulars.mct == pytest.approx(exact.mct, rel=1e-12)


@pytest.mark.parametrize("weight,shift", [(24000.0, 4.0), (48000.0, 8.0)])
def test_box_trim_matches_the_wallsided_closed_form_at_finite_angle(box, weight, shift):
    """No linearisation in the solve: 8 degrees of trim, agreement to round-off."""
    solution = equilibrium_trim(box, weight=weight, shift=shift, lcb_tolerance=1e-10)
    closed = math.degrees(math.atan(box_trim_tangent(BOX[0], BOX[1], weight, shift)))
    assert solution.trim_angle_deg == pytest.approx(closed, rel=1e-8)
    assert abs(solution.residual_m) < 1e-8


def test_the_two_centimetre_routes_agree(box_particulars):
    """``w d / MCT`` centimetres and ``atan(w d / (rho I_L))`` are the same answer."""
    _, _, change_cm = small_angle_trim(box_particulars, 48000.0, 8.0)
    loading_route = 48000.0 / 1000.0 * 8.0 / box_particulars.mct
    assert change_cm == pytest.approx(loading_route, rel=1e-12)


def test_trim_answers_the_moment_not_the_weights(box):
    """Halving the weight and doubling the shift is the same moment, hence the same trim."""
    first = equilibrium_trim(box, weight=24000.0, shift=8.0, lcb_tolerance=1e-9)
    second = equilibrium_trim(box, weight=48000.0, shift=4.0, lcb_tolerance=1e-9)
    assert first.trim_angle_deg == pytest.approx(second.trim_angle_deg, rel=1e-9)


def test_a_weight_moved_aft_trims_the_other_way(box):
    solution = equilibrium_trim(box, weight=24000.0, shift=-4.0, lcb_tolerance=1e-9)
    assert solution.trim_angle_deg < 0.0
    assert solution.lcb < solution.lcb_0


def test_the_indexed_cutter_matches_the_reference_route(hull):
    """The fast path and :func:`tensorcad.bonjean.section_properties` are independent."""
    stations = _stations(hull.mesh, 21)
    cutter = _SectionCutter(hull.mesh, stations)
    for index, station in enumerate(stations):
        for waterline in (-0.7, -0.4, -0.16, 0.0):
            reference = section_properties(hull.mesh, station, waterline).area
            assert cutter.area(index, waterline) == pytest.approx(reference, abs=1e-12)


def test_the_waterline_solve_reproduces_the_volume(hull):
    upright = hydrostatic_particulars(hull.mesh, 0.0, station_count=41)
    height = waterline_for_displacement(
        hull.mesh, trim_angle=0.0, displaced_volume=0.8 * upright.volume, station_count=41
    )
    measured = hydrostatic_particulars(hull.mesh, height, station_count=41)
    assert measured.volume == pytest.approx(0.8 * upright.volume, rel=1e-6)


def test_a_ring_waterline_is_nudged_and_says_so(hull):
    """The generator puts rings on even levels; a degenerate plane cut is nudged, not hidden."""
    row = hydrostatic_particulars(hull.mesh, -0.25 * hull.parameters.draft, station_count=21)
    assert row.waterline_nudged is True
    assert row.waterplane_area > 0.0
    neighbour = hydrostatic_particulars(hull.mesh, -0.25 * hull.parameters.draft + 1e-3, station_count=21)
    assert row.waterplane_area == pytest.approx(neighbour.waterplane_area, rel=1e-3)


def test_hull_part_load_volume_against_the_analytic_scaling(hull):
    z = -0.2 * hull.parameters.draft
    row = hydrostatic_particulars(hull.mesh, z, station_count=41)
    reference = (
        analytic_volume(hull.parameters.length, hull.parameters.beam, hull.parameters.draft)
        * (row.draft / hull.parameters.draft) ** 1.5
    )
    assert row.volume == pytest.approx(reference, rel=1.5e-2)


def test_tpc_is_the_derivative_of_displacement(hull):
    """TPC is defined as dV/dz, so the two routes must meet -- they use different cuts."""
    z = -0.2 * hull.parameters.draft
    row = hydrostatic_particulars(hull.mesh, z, station_count=41)
    step = 1e-3
    up = hydrostatic_particulars(hull.mesh, z + step, station_count=41).volume
    down = hydrostatic_particulars(hull.mesh, z - step, station_count=41).volume
    derivative = 0.5 * (up - down) / step * SEAWATER_DENSITY / 1e5
    assert row.tpc == pytest.approx(derivative, rel=5e-3)


def test_hull_trim_agrees_with_the_small_angle_answer(hull):
    z = -0.2 * hull.parameters.draft
    solution = equilibrium_trim(hull.mesh, weight=300.0, shift=1.0, waterline_z=z, station_count=21)
    assert solution.trim_angle_deg == pytest.approx(solution.small_angle_deg, rel=5e-3)
    small_angle_cm = 100.0 * hull.parameters.length * solution.trim_angle_rad
    assert solution.change_of_trim_cm == pytest.approx(small_angle_cm, rel=6e-3)


def test_wigley_trim_agrees_with_the_small_angle_answer():
    wigley = generate_wigley()
    z = -0.2 * wigley.parameters.draft
    solution = equilibrium_trim(wigley.mesh, weight=300.0, shift=0.5, waterline_z=z, station_count=21)
    assert solution.trim_angle_deg == pytest.approx(solution.small_angle_deg, rel=5e-3)


def test_the_hydrostatic_table_is_monotone_in_draft(hull):
    drafts = [-0.6, -0.5, -0.4, -0.3, -0.2, -0.1]  # shallow first: draft from keel rises
    rows = hydrostatic_table(hull.mesh, drafts, station_count=21)
    assert [row.draft for row in rows] == pytest.approx([0.2, 0.3, 0.4, 0.5, 0.6, 0.7])
    tpc = [row.tpc for row in rows]
    mct = [row.mct for row in rows]
    assert all(tpc[index] < tpc[index + 1] for index in range(len(tpc) - 1))
    assert all(mct[index] < mct[index + 1] for index in range(len(mct) - 1))


def test_trimmed_hydrostatics_holds_the_displacement(hull):
    upright = hydrostatic_particulars(hull.mesh, -0.16, station_count=21)
    trimmed = trimmed_hydrostatics(
        hull.mesh, trim_angle=math.radians(0.5), preserved=upright.volume, station_count=21
    )
    assert trimmed.volume == pytest.approx(upright.volume, rel=1e-6)
    assert trimmed.lcb > upright.lcb


def test_records_serialise_with_their_schema(box_particulars):
    payload = box_particulars.to_dict()
    assert payload["schema"] == "tensorcad.hydrostatic_particulars/1"
    assert payload["tpc_t_per_cm"] == pytest.approx(0.82, rel=1e-12)
    assert payload["mct_t_m_per_cm"] == pytest.approx(box_particulars.mct, rel=1e-8)


@pytest.mark.parametrize(
    "call",
    [
        lambda mesh: hydrostatic_particulars(mesh, -10.0),
        lambda mesh: hydrostatic_particulars(mesh, 0.0, water_density=-1.0),
        lambda mesh: hydrostatic_particulars(mesh, 0.0, station_count=2),
        lambda mesh: hydrostatic_table(mesh, ()),
        lambda mesh: equilibrium_trim(mesh, weight=-1.0, shift=1.0),
        lambda mesh: waterline_for_displacement(mesh, trim_angle=0.0, displaced_volume=-1.0),
    ],
)
def test_impossible_questions_are_refused(box, call):
    with pytest.raises(TrimError):
        call(box)


def test_a_hull_capped_at_its_waterline_cannot_be_trimmed_at_full_displacement(hull):
    """Honest limit: with no freeboard there is no waterline to trim to, and we say so."""
    with pytest.raises(TrimError):
        equilibrium_trim(hull.mesh, weight=600.0, shift=1.0, waterline_z=0.0, station_count=21)


def test_the_module_needs_no_third_party_stack():
    import subprocess
    import sys

    script = (
        "import sys, tensorcad.trim;"
        "loaded = {name.split('.')[0] for name in sys.modules};"
        "bad = loaded & {'numpy', 'scipy', 'torch'};"
        "assert not bad, sorted(bad)"
    )
    completed = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
