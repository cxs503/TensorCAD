"""Verification debt closed: Bonjean closed forms for the asymmetric hull family.

The asymmetric generator (0.17.0) gave the parametric hull a taper exponent per
half, and ``analytic_hydrostatics`` followed. The Bonjean module did not: its
closed form still hard-wired ``BREADTH_EXPONENT``, so the one tool that answers
"how much area is under this waterline at this station" could not be checked
against a hull whose two halves differ. These tests close that gap and pin the
numbers with tolerances rather than with prose.

Every expected value here is either a closed form or a previously measured
quantity; nothing is read back from the code under test.
"""
from __future__ import annotations

import math

import pytest

from tensorcad import HullParameters, analytic_volume, bonjean_table, generate_hull
from tensorcad.bonjean import (
    analytic_section_area,
    analytic_section_centroid_z,
    section_properties,
)
from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.hull import (
    BREADTH_EXPONENT,
    analytic_longitudinal_centre,
    analytic_vertical_centre_buoyancy,
)
from tensorcad.mesh import orient_outward

LENGTH, BEAM, DRAFT = 10.0, 2.0, 0.8
FORE, AFT = 0.5, 0.9


def _asymmetric(nx: int = 41, nz: int = 13):
    parameters = HullParameters(
        length=LENGTH,
        beam=BEAM,
        draft=DRAFT,
        nx=nx,
        nz=nz,
        fore_exponent=FORE,
        aft_exponent=AFT,
    )
    return generate_hull(parameters)


def _box(minimum, maximum) -> Mesh:
    x0, y0, z0 = minimum
    x1, y1, z1 = maximum
    builder = MeshBuilder()
    corners = [
        (x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
        (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1),
    ]
    ids = [builder.add_vertex(point) for point in corners]
    for face in (
        (0, 3, 2), (0, 2, 1), (4, 5, 6), (4, 6, 7), (0, 1, 5), (0, 5, 4),
        (2, 3, 7), (2, 7, 6), (1, 2, 6), (1, 6, 5), (0, 4, 7), (0, 7, 3),
    ):
        builder.add_triangle(*[ids[index] for index in face])
    return orient_outward(builder.build()).mesh


# --------------------------------------------------------------------------
# Back-compatibility: the symmetric closed form must not move a digit
# --------------------------------------------------------------------------


@pytest.mark.parametrize("station", [1.0, 2.5, 5.0, 7.5, 9.0])
@pytest.mark.parametrize("waterline", [-0.8, -0.6, -0.2, 0.0])
def test_default_exponents_reproduce_the_shipped_expression_bit_for_bit(
    station: float, waterline: float
) -> None:
    """The 0.18 closed form, recomputed here from the literal expression."""
    zeta = waterline + DRAFT
    shipped = (2.0 / 3.0) * BEAM * math.pow(
        math.sin(math.pi * station / LENGTH), BREADTH_EXPONENT
    ) * math.pow(zeta, 1.5) / math.sqrt(DRAFT)
    measured = analytic_section_area(LENGTH, BEAM, DRAFT, station, waterline)
    assert measured == shipped


def test_default_exponents_equal_the_explicit_symmetric_pair() -> None:
    for station in (0.5, 3.0, 5.0, 8.0):
        assert analytic_section_area(
            LENGTH, BEAM, DRAFT, station, -0.3
        ) == analytic_section_area(
            LENGTH,
            BEAM,
            DRAFT,
            station,
            -0.3,
            fore_exponent=BREADTH_EXPONENT,
            aft_exponent=BREADTH_EXPONENT,
        )


# --------------------------------------------------------------------------
# The split itself
# --------------------------------------------------------------------------


def test_forward_station_uses_the_fore_exponent_and_aft_the_aft_one() -> None:
    """Ablaft midship only the aft exponent may appear, forward only the fore."""
    waterline = -0.2
    zeta = waterline + DRAFT
    for station in (1.0, 3.0, 4.0):
        expected = (2.0 / 3.0) * BEAM * math.sin(math.pi * station / LENGTH) ** AFT * (
            zeta**1.5
        ) / math.sqrt(DRAFT)
        assert analytic_section_area(
            LENGTH, BEAM, DRAFT, station, waterline, fore_exponent=FORE, aft_exponent=AFT
        ) == pytest.approx(expected, rel=1e-14)
    for station in (6.0, 7.0, 9.0):
        expected = (2.0 / 3.0) * BEAM * math.sin(math.pi * station / LENGTH) ** FORE * (
            zeta**1.5
        ) / math.sqrt(DRAFT)
        assert analytic_section_area(
            LENGTH, BEAM, DRAFT, station, waterline, fore_exponent=FORE, aft_exponent=AFT
        ) == pytest.approx(expected, rel=1e-14)


def test_midship_is_a_fixed_point_of_the_split() -> None:
    """``sin(pi/2) = 1`` for every exponent, so midship cannot see the split."""
    midship = analytic_section_area(LENGTH, BEAM, DRAFT, LENGTH / 2, -0.2)
    split = analytic_section_area(
        LENGTH, BEAM, DRAFT, LENGTH / 2, -0.2, fore_exponent=FORE, aft_exponent=AFT
    )
    assert midship == split
    zeta = -0.2 + DRAFT
    assert midship == pytest.approx((2.0 / 3.0) * BEAM * zeta**1.5 / math.sqrt(DRAFT), rel=1e-15)
    assert midship == pytest.approx(0.692820323027551, rel=1e-15)


def test_the_split_is_monotone_in_each_exponent() -> None:
    station, waterline = 2.0, -0.2
    areas = [
        analytic_section_area(
            LENGTH, BEAM, DRAFT, station, waterline, aft_exponent=exponent
        )
        for exponent in (0.4, 0.7, 1.0, 1.4)
    ]
    assert areas == sorted(areas, reverse=True)  # sin < 1 => larger exponent is thinner
    forward = [
        analytic_section_area(
            LENGTH, BEAM, DRAFT, LENGTH - station, waterline, fore_exponent=exponent
        )
        for exponent in (0.4, 0.7, 1.0, 1.4)
    ]
    assert forward == sorted(forward, reverse=True)


def test_swapping_the_exponents_mirrors_the_closed_form_about_midship() -> None:
    for station in (1.0, 2.5, 7.5, 9.0):
        mirrored = LENGTH - station
        left = analytic_section_area(
            LENGTH, BEAM, DRAFT, station, -0.2, fore_exponent=FORE, aft_exponent=AFT
        )
        right = analytic_section_area(
            LENGTH, BEAM, DRAFT, mirrored, -0.2, fore_exponent=AFT, aft_exponent=FORE
        )
        assert left == pytest.approx(right, rel=1e-14)


def test_analytic_section_area_integrates_to_the_analytic_volume_of_the_split() -> None:
    """The station-wise closed form and the whole-hull closed form are one model."""
    stations = 4001
    spacing = LENGTH / (stations - 1)
    values = [
        analytic_section_area(
            LENGTH,
            BEAM,
            DRAFT,
            spacing * index,
            0.0,
            fore_exponent=FORE,
            aft_exponent=AFT,
        )
        for index in range(stations)
    ]
    total = values[0] + values[-1]
    for index in range(1, stations - 1):
        total += (4.0 if index % 2 == 1 else 2.0) * values[index]
    simpson = total * spacing / 3.0
    reference = analytic_volume(LENGTH, BEAM, DRAFT, fore_exponent=FORE, aft_exponent=AFT)
    assert simpson == pytest.approx(reference, rel=2e-6)


# --------------------------------------------------------------------------
# The section centroid, exact for the family and independent of the exponent
# --------------------------------------------------------------------------


def test_section_centroid_formula_reaches_the_analytic_vcb_at_the_waterline() -> None:
    assert analytic_section_centroid_z(DRAFT, 0.0) == pytest.approx(
        analytic_vertical_centre_buoyancy(DRAFT), rel=1e-15
    )
    assert analytic_section_centroid_z(DRAFT, 0.0) == pytest.approx(-0.32, rel=1e-15)


def test_section_centroid_is_three_fifths_of_the_draft_above_the_keel() -> None:
    for waterline in (-0.79, -0.6, -0.2, 0.0):
        zeta = waterline + DRAFT
        assert analytic_section_centroid_z(DRAFT, waterline) == pytest.approx(
            -DRAFT + 0.6 * zeta, rel=1e-15
        )


def test_section_centroid_does_not_see_the_station_or_the_exponent() -> None:
    """The vertical profile is station-independent, so this is a cut error test."""
    model = _asymmetric()
    values = [
        section_properties(model.mesh, station, 0.0).centroid_z
        for station in (2.0, 3.5, 5.0, 6.5, 8.0)
    ]
    exact = analytic_section_centroid_z(DRAFT, 0.0)
    spread = max(values) - min(values)
    assert spread < 1e-9  # measured 4.2e-12 m: the cut is exact at every station
    for value in values:
        assert value == pytest.approx(exact, rel=0.01)


def test_section_centroid_mesh_error_is_below_one_percent_and_converges() -> None:
    exact = analytic_section_centroid_z(DRAFT, 0.0)
    errors = []
    for nz in (13, 21, 41):
        model = _asymmetric(nx=81, nz=nz)
        measured = section_properties(model.mesh, LENGTH / 2, 0.0).centroid_z
        errors.append(abs(measured - exact) / abs(exact))
    assert errors[0] < 0.01
    assert errors[1] < errors[0]
    assert errors[2] < errors[1]
    assert errors[2] < 0.002


def test_section_centroid_on_a_barge_is_the_mid_draft() -> None:
    mesh = _box((0.0, -1.0, -1.0), (4.0, 1.0, 0.0))
    properties = section_properties(mesh, 2.0, -0.25)
    assert properties.centroid_z == pytest.approx(-0.625, rel=1e-15)


# --------------------------------------------------------------------------
# Bonjean curves of the asymmetric hull against both closed forms
# --------------------------------------------------------------------------


def test_asymmetric_bonjean_areas_track_the_split_closed_form() -> None:
    """Worst error over 5 stations x 4 waterlines, and it converges vertically."""
    errors = []
    for nz in (13, 21):
        model = _asymmetric(nz=nz)
        worst = 0.0
        for station in (LENGTH * 0.1, LENGTH * 0.3, LENGTH / 2, LENGTH * 0.7, LENGTH * 0.9):
            for waterline in (-0.6, -0.4, -0.1, 0.0):
                measured = section_properties(model.mesh, station, waterline).area
                expected = analytic_section_area(
                    LENGTH,
                    BEAM,
                    DRAFT,
                    station,
                    waterline,
                    fore_exponent=FORE,
                    aft_exponent=AFT,
                )
                worst = min(worst, (measured - expected) / expected)
        errors.append(worst)
    assert errors[0] == pytest.approx(-0.05308, abs=3e-4)  # -5.308 % at nz = 13
    assert errors[1] > errors[0]  # monotone towards zero from below
    assert errors[1] > -0.03


def test_asymmetric_mesh_is_actually_asymmetric() -> None:
    """Guards the split against silently collapsing back to the symmetric hull."""
    model = _asymmetric()
    aft = section_properties(model.mesh, LENGTH * 0.25, -0.2).area
    fore = section_properties(model.mesh, LENGTH * 0.75, -0.2).area
    assert fore > aft * 1.05  # exponent 0.9 thins the afterbody, 0.5 fills the forebody
    symmetric = generate_hull(HullParameters(length=LENGTH, beam=BEAM, draft=DRAFT))
    half = section_properties(symmetric.mesh, LENGTH * 0.25, -0.2).area
    assert aft < half * 0.95


def test_bonjean_volume_matches_the_split_analytic_volume() -> None:
    model = _asymmetric()
    table = bonjean_table(model.mesh, station_count=41, waterline_count=13)
    reference = analytic_volume(LENGTH, BEAM, DRAFT, fore_exponent=FORE, aft_exponent=AFT)
    relative = (table.integrated_volume(0.0) - reference) / reference
    assert relative == pytest.approx(-0.0079, abs=1e-3)
    assert abs(relative) < 0.01


def test_bonjean_longitudinal_centre_matches_the_split_analytic_lcb() -> None:
    """The independent measurement: no surface integral enters this number."""
    model = _asymmetric()
    table = bonjean_table(model.mesh, station_count=41, waterline_count=13)
    measured = table.longitudinal_centre(0.0)
    reference = analytic_longitudinal_centre(LENGTH, fore_exponent=FORE, aft_exponent=AFT)
    assert measured > LENGTH / 2 + 0.2  # the thin afterbody pulls the centre forward
    assert (measured - reference) / reference == pytest.approx(-7.0e-4, abs=3e-4)


def test_bonjean_longitudinal_centre_of_the_symmetric_hull_is_midship() -> None:
    model = generate_hull(HullParameters(length=LENGTH, beam=BEAM, draft=DRAFT))
    table = bonjean_table(model.mesh, station_count=41, waterline_count=13)
    assert table.longitudinal_centre(0.0) == pytest.approx(LENGTH / 2, rel=1e-4)


def test_longitudinal_centre_identity_survives_a_pure_aft_taper() -> None:
    """``LCB(p, q) + LCB(q, p) = L`` exactly; a joint test of both halves."""
    for fore, aft in ((0.5, 0.9), (0.4, 1.2), (1.0, 0.6)):
        total = analytic_longitudinal_centre(
            LENGTH, fore_exponent=fore, aft_exponent=aft
        ) + analytic_longitudinal_centre(LENGTH, fore_exponent=aft, aft_exponent=fore)
        assert total == pytest.approx(LENGTH, rel=1e-14)


def test_bonjean_first_moment_of_the_symmetric_hull_is_half_the_volume() -> None:
    model = generate_hull(HullParameters(length=LENGTH, beam=BEAM, draft=DRAFT))
    table = bonjean_table(model.mesh, station_count=41, waterline_count=13)
    volume = table.integrated_volume(0.0)
    moment = table.integrated_first_moment(0.0)
    assert moment / volume == pytest.approx(LENGTH / 2, rel=1e-4)
    assert moment == pytest.approx(0.5 * LENGTH * volume, rel=1e-4)


# --------------------------------------------------------------------------
# Input validation
# --------------------------------------------------------------------------


@pytest.mark.parametrize("exponent", [0.0, -0.5, float("nan"), float("inf")])
def test_bad_exponents_are_rejected(exponent: float) -> None:
    with pytest.raises(ValueError):
        analytic_section_area(
            LENGTH, BEAM, DRAFT, 2.0, -0.2, fore_exponent=exponent
        )
    with pytest.raises(ValueError):
        analytic_section_area(
            LENGTH, BEAM, DRAFT, 2.0, -0.2, aft_exponent=exponent
        )


@pytest.mark.parametrize("station", [-0.1, 10.1, float("nan")])
def test_bad_stations_are_rejected(station: float) -> None:
    with pytest.raises(ValueError):
        analytic_section_area(LENGTH, BEAM, DRAFT, station, -0.2)


@pytest.mark.parametrize("draft", [0.0, -1.0, float("nan")])
def test_bad_draft_is_rejected_by_the_centroid(draft: float) -> None:
    with pytest.raises(ValueError):
        analytic_section_centroid_z(draft, -0.2)


def test_centroid_below_the_keel_is_rejected() -> None:
    with pytest.raises(ValueError):
        analytic_section_centroid_z(DRAFT, -DRAFT)
    assert analytic_section_centroid_z(DRAFT, -DRAFT + 1e-12) == pytest.approx(
        -DRAFT + 0.6e-12, abs=1e-15
    )
