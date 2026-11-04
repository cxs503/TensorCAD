"""Verification of the Bonjean curves against references the mesh never saw.

Three kinds of reference, in increasing order of independence:

1. **Closed form on a barge.** A wall-sided box has a rectangular section, so
   area, perimeter and centroid are exact rationals and the section integral can
   be compared to the box volume without a quadrature error at all.
2. **Closed form on the hull family.** ``(2/3) beam sin(pi x/L)^0.7 zeta^(3/2) /
   sqrt(draft)`` is an integral of the analytic surface; the mesh only samples
   it, so the residual is pure discretisation and it must shrink at second order
   as the vertical resolution grows.
3. **An identity across two independent integrations.** For a loft that is linear
   in ``x`` between stations -- which the generator's loft is, exactly -- the
   volume between two stations is ``(A_i + A_{i+1}) / 2 * dx``. So the trapezoid
   integral of the *section areas* must equal the *surface-flux* displacement
   from :func:`tensorcad.compute_hydrostatics`, which sums tetrahedra and never
   builds a section. It does, to round-off. That turns the remaining error into a
   statement about surface sampling only, which is what the closed form then
   measures.
"""
from __future__ import annotations

import itertools
import json
import math

import pytest

from tensorcad import (
    HullParameters,
    analytic_volume,
    bonjean_curve,
    bonjean_table,
    compute_hydrostatics,
    generate_hull,
    section_loops,
    section_properties,
)
from tensorcad import mesh as meshlib
from tensorcad.bonjean import analytic_section_area, default_waterlines
from tensorcad.geometry import Mesh, MeshBuilder

LENGTH, BEAM, DRAFT = 10.0, 2.0, 0.8


def box(minimum: tuple[float, float, float], maximum: tuple[float, float, float]) -> Mesh:
    """Closed, outward-oriented axis-aligned box, as the hydrostatics tests build it."""
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
    return meshlib.orient_outward(builder.build()).mesh


@pytest.fixture(scope="module")
def barge() -> Mesh:
    return box((0.0, -BEAM / 2, -DRAFT), (LENGTH, BEAM / 2, 0.0))


@pytest.fixture(scope="module")
def hull_mesh() -> Mesh:
    return generate_hull(HullParameters(length=LENGTH, beam=BEAM, draft=DRAFT, nx=41, nz=13)).mesh


# --------------------------------------------------------------------------- #
# 1. Exact geometry: the barge
# --------------------------------------------------------------------------- #
def test_barge_section_is_one_closed_loop(barge: Mesh) -> None:
    loops = section_loops(barge, 5.0)
    assert len(loops) == 1
    assert loops[0][0] == loops[0][-1] or len(loops[0]) >= 4


@pytest.mark.parametrize("waterline", [-0.8, -0.6, -0.4, -0.000001, 0.0])
def test_barge_section_area_is_exact(barge: Mesh, waterline: float) -> None:
    """``A = beam * (waterline - keel)`` with nothing to discretise."""
    expected = BEAM * (waterline + DRAFT)
    measured = section_properties(barge, 5.0, waterline).area
    assert measured == pytest.approx(expected, rel=1e-12, abs=1e-15)


def test_barge_section_centroid_and_perimeter_are_exact(barge: Mesh) -> None:
    properties = section_properties(barge, 5.0)
    assert properties.centroid_y == pytest.approx(0.0, abs=1e-12)
    assert properties.centroid_z == pytest.approx(-DRAFT / 2, rel=1e-12)
    # bottom + two sides up to the waterline: 2 + 2 * 0.8 (the top edge is not
    # part of the wetted outline, but the cut polygon does close along it).
    assert properties.perimeter == pytest.approx(2 * BEAM + 2 * DRAFT, rel=1e-12)


def test_barge_bonjean_is_linear_in_the_waterline(barge: Mesh) -> None:
    curve = bonjean_curve(barge, 5.0, count=9)
    for waterline, area in zip(curve.waterlines, curve.areas, strict=False):
        assert area == pytest.approx(BEAM * (waterline + DRAFT), rel=1e-12, abs=1e-15)
    assert curve.area_at(-0.4) == pytest.approx(BEAM * 0.4, rel=1e-12)


def test_barge_bonjean_recovers_the_box_volume(barge: Mesh) -> None:
    table = bonjean_table(barge, station_count=11, waterline_count=5)
    measured = table.integrated_volume(0.0)
    assert measured == pytest.approx(LENGTH * BEAM * DRAFT, rel=1e-12)


# --------------------------------------------------------------------------- #
# 2. Closed form: the parametric hull family
# --------------------------------------------------------------------------- #
def test_analytic_section_area_matches_the_midship_reference() -> None:
    """At midship and the design waterline the closed form is the old reference."""
    from tensorcad.hull import analytic_maximum_section_area

    measured = analytic_section_area(LENGTH, BEAM, DRAFT, LENGTH / 2, 0.0)
    assert measured == pytest.approx(analytic_maximum_section_area(BEAM, DRAFT), rel=1e-12)


def test_analytic_section_area_integrates_to_the_analytic_volume() -> None:
    """The station integral of the closed form is the closed-form volume.

    Integrated with Gauss-Legendre at 400 nodes, i.e. a quadrature the closed form
    has never been fitted to: agreement at 1e-9 relative is the internal
    consistency of the analytic pair, and it is what licenses using the station
    values as the reference for the mesh below.
    """
    nodes, weights = _gauss_legendre(400)
    total = 0.0
    for node, weight in zip(nodes, weights, strict=False):
        station = 0.5 * LENGTH * (node + 1.0)
        total += weight * analytic_section_area(LENGTH, BEAM, DRAFT, station, 0.0)
    total *= 0.5 * LENGTH
    assert total == pytest.approx(analytic_volume(LENGTH, BEAM, DRAFT), rel=1e-9)


def _gauss_legendre(count: int) -> tuple[list[float], list[float]]:
    """Newton on the Legendre polynomial: nodes and weights on ``[-1, 1]``."""
    nodes: list[float] = []
    weights: list[float] = []
    for index in range(1, count + 1):
        x = math.cos(math.pi * (index - 0.25) / (count + 0.5))
        for _ in range(100):
            p0, p1 = 1.0, x
            for order in range(2, count + 1):
                p0, p1 = p1, ((2 * order - 1) * x * p1 - (order - 1) * p0) / order
            derivative = count * (x * p1 - p0) / (x * x - 1.0)
            step = p1 / derivative
            x -= step
            if abs(step) < 1e-15:
                break
        nodes.append(x)
        weights.append(2.0 / ((1.0 - x * x) * derivative * derivative))
    return nodes, weights


@pytest.mark.parametrize("station", [1.0, 2.5, 5.0, 7.5, 9.0])
@pytest.mark.parametrize("waterline", [-0.6, -0.4, 0.0])
def test_hull_section_area_tracks_the_closed_form(
    hull_mesh: Mesh, station: float, waterline: float
) -> None:
    """The default 41x13 grid sits within 6 % everywhere at the worst waterline.

    The worst case is the shallow waterline, where only three of the thirteen
    levels are immersed and the piecewise-linear skin is a coarse stand-in for
    ``sqrt(z)``; :func:`test_hull_section_area_converges_vertically` is the proof
    that this is discretisation and not a wrong formula.
    """
    measured = section_properties(hull_mesh, station, waterline).area
    reference = analytic_section_area(LENGTH, BEAM, DRAFT, station, waterline)
    assert measured == pytest.approx(reference, rel=0.06)
    # And the sign is systematic: a linear skin under a concave profile always
    # under-estimates, never over-estimates.
    assert measured <= reference


def test_hull_section_area_converges_vertically() -> None:
    """Second-order convergence in ``nz`` at the shallow waterline."""
    errors = []
    for nz in (13, 21, 41, 81):
        mesh = generate_hull(HullParameters(length=LENGTH, beam=BEAM, draft=DRAFT, nx=21, nz=nz)).mesh
        measured = section_properties(mesh, LENGTH / 2, -0.6).area
        reference = analytic_section_area(LENGTH, BEAM, DRAFT, LENGTH / 2, -0.6)
        errors.append(abs(measured / reference - 1.0))
    assert 0.05 < errors[0] < 0.06
    assert errors[-1] < 0.004
    # Halving the vertical spacing cuts the error by at least 1.8x over the range
    # measured (2.09, 2.75, 2.77 from 13->21->41->81 levels): monotone, and second
    # order once the shallow waterline stops being resolved by three segments.
    for coarse, fine in itertools.pairwise(errors):
        assert coarse > 1.8 * fine


# --------------------------------------------------------------------------- #
# 3. The identity across two independent integrations
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("nx", [11, 21, 41])
def test_section_integral_equals_the_surface_flux(nx: int) -> None:
    """Trapezoid over the station areas == the tetrahedron sum, to round-off.

    The loft is linear in ``x`` between stations, so the exact volume of the
    piecewise-linear body is the trapezoid rule on its exact station areas. Any
    disagreement beyond round-off would mean the section cut and the flux are not
    measuring the same body.
    """
    parameters = HullParameters(length=LENGTH, beam=BEAM, draft=DRAFT, nx=nx, nz=13)
    mesh = generate_hull(parameters).mesh
    stations = tuple(LENGTH * index / (nx - 1) for index in range(nx))
    table = bonjean_table(mesh, stations=stations, waterline_count=3)

    by_sections = table.integrated_volume(0.0, method="trapezoid")
    by_flux = compute_hydrostatics(mesh, 0.0).displacement_volume
    assert by_sections == pytest.approx(by_flux, rel=1e-12)


def test_simpson_beats_the_trapezoid_on_the_analytic_volume(hull_mesh: Mesh) -> None:
    """Against the closed-form volume, Simpson is the better rule -- as it should be.

    The station areas the mesh reports are a slightly low *sample* of the analytic
    area (piecewise-linear skin under a concave profile), so the trapezoid rule
    inherits that bias while Simpson, weighting the interior more heavily, lands
    closer to the analytic volume. Both are recorded so the claim is not a single
    number chosen after the fact.
    """
    stations = tuple(LENGTH * index / 40 for index in range(41))
    table = bonjean_table(hull_mesh, stations=stations, waterline_count=3)
    analytic = analytic_volume(LENGTH, BEAM, DRAFT)
    simpson_error = abs(table.integrated_volume(0.0, method="simpson") / analytic - 1.0)
    trapezoid_error = abs(table.integrated_volume(0.0, method="trapezoid") / analytic - 1.0)
    assert simpson_error < trapezoid_error
    assert simpson_error < 0.008


# --------------------------------------------------------------------------- #
# 4. Contract: curves, tables, edges
# --------------------------------------------------------------------------- #
def test_default_waterlines_span_the_section(hull_mesh: Mesh) -> None:
    properties = section_properties(hull_mesh, LENGTH / 2)
    waterlines = default_waterlines(properties, count=5)
    assert waterlines[0] == pytest.approx(properties.minimum_z)
    assert waterlines[-1] == pytest.approx(properties.maximum_z)


def test_area_is_monotonic_in_the_waterline(hull_mesh: Mesh) -> None:
    curve = bonjean_curve(hull_mesh, LENGTH / 2, count=13)
    assert all(b >= a - 1e-15 for a, b in itertools.pairwise(curve.areas))
    assert curve.areas[-1] == pytest.approx(
        section_properties(hull_mesh, LENGTH / 2).area, rel=1e-12
    )


def test_waterline_above_the_section_is_the_whole_section(hull_mesh: Mesh) -> None:
    whole = section_properties(hull_mesh, LENGTH / 2)
    above = section_properties(hull_mesh, LENGTH / 2, whole.maximum_z + 1.0)
    assert above.area == pytest.approx(whole.area, rel=1e-12)


def test_waterline_below_the_keel_is_empty(hull_mesh: Mesh) -> None:
    properties = section_properties(hull_mesh, LENGTH / 2, -10.0)
    assert properties.area == 0.0
    assert properties.loops >= 1
    assert math.isnan(properties.centroid_z)


def test_station_outside_the_mesh_has_no_section(hull_mesh: Mesh) -> None:
    assert section_loops(hull_mesh, LENGTH + 1.0) == ()
    properties = section_properties(hull_mesh, LENGTH + 1.0)
    assert properties.area == 0.0
    assert properties.loops == 0


def test_station_must_be_finite(hull_mesh: Mesh) -> None:
    with pytest.raises(ValueError):
        section_loops(hull_mesh, float("nan"))


def test_analytic_section_area_rejects_bad_input() -> None:
    with pytest.raises(ValueError):
        analytic_section_area(0.0, BEAM, DRAFT, 0.0, 0.0)
    with pytest.raises(ValueError):
        analytic_section_area(LENGTH, BEAM, DRAFT, -1.0, 0.0)
    assert analytic_section_area(LENGTH, BEAM, DRAFT, 5.0, -DRAFT) == 0.0


def test_table_and_curves_serialise(hull_mesh: Mesh) -> None:
    table = bonjean_table(hull_mesh, station_count=5, waterline_count=4)
    payload = table.to_dict()
    assert json.loads(json.dumps(payload)) == payload
    assert len(payload["stations_m"]) == 5
    assert len(payload["waterlines_m"]) == 4
    assert all(len(row) == 4 for row in payload["areas_m2"])
    # Bilinear reads agree with the station curves they came from. The end
    # stations of a hull close on the centreline, so they carry no waterlines.
    for station, curve in zip(table.stations, table.curves, strict=False):
        if not curve.waterlines:
            continue
        assert table.area_at(station, curve.waterlines[-1]) == pytest.approx(
            curve.areas[-1], rel=1e-12
        )
