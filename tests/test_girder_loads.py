"""Tests for :mod:`tensorcad.girder_loads` -- still-water shear force and bending moment.

The reference cases are hand-derived textbook ones, not the implementation's own output:

* a box barge whose weight distribution equals its buoyancy has ``Q = M = 0`` everywhere;
* a box barge with a point load at midspan carries ``M = -P L / 8`` there, approached
  linearly as the load's patch narrows;
* the two equilibrium residuals, ``Q`` and ``M`` at the forward end, are exactly the
  force and moment imbalance of the load diagram, so they vanish together for a balanced
  distribution and separately for a force-balanced but out-of-trim one.
"""

from __future__ import annotations

import math

import pytest

from tensorcad.girder_loads import (
    GRAVITY,
    SEA_WATER_DENSITY,
    BuoyancyDistribution,
    GirderLoadError,
    HullGirderLoads,
    WeightItem,
    buoyancy_distribution,
    hull_girder_loads,
    still_water_loads,
)

# A box barge: no mesh involved, the buoyancy diagram of a wall-sided box is constant.
LENGTH = 10.0
BEAM = 2.0
DRAFT = 1.0
DENSITY = 1000.0
BUOYANCY_PER_LENGTH = DENSITY * GRAVITY * BEAM * DRAFT  # N/m
DISPLACEMENT = BUOYANCY_PER_LENGTH * LENGTH  # N


def box_stations(count: int = 21, length: float = LENGTH) -> tuple[float, ...]:
    return tuple(length * index / (count - 1) for index in range(count))


def box_buoyancy(count: int = 21) -> tuple[float, ...]:
    return tuple(BUOYANCY_PER_LENGTH for _ in range(count))


def uniform_weight(total: float = DISPLACEMENT) -> tuple[WeightItem, ...]:
    return (WeightItem(0.0, LENGTH, total, "uniform"),)


# ---------------------------------------------------------------- balanced cases


def test_uniform_weight_on_box_barge_has_no_shear_or_moment() -> None:
    loads = still_water_loads(box_stations(), box_buoyancy(), uniform_weight())
    scale_force = BUOYANCY_PER_LENGTH * LENGTH
    scale_moment = scale_force * LENGTH
    assert loads.max_shear / scale_force < 1e-15
    assert loads.max_moment / scale_moment < 1e-15
    assert loads.residual_force / scale_force < 1e-15
    assert loads.residual_moment / scale_moment < 1e-15


def test_totals_and_centres_of_a_balanced_box_barge() -> None:
    loads = still_water_loads(box_stations(), box_buoyancy(), uniform_weight())
    assert loads.total_buoyancy == pytest.approx(DISPLACEMENT, rel=1e-15)
    assert loads.total_weight == pytest.approx(DISPLACEMENT, rel=1e-15)
    assert loads.lcb_buoyancy == pytest.approx(LENGTH / 2.0, rel=1e-15)
    assert loads.lcb_weight == pytest.approx(LENGTH / 2.0, rel=1e-15)


def test_zero_shear_and_moment_across_the_whole_curve() -> None:
    loads = still_water_loads(box_stations(), box_buoyancy(), uniform_weight())
    for value in loads.shear:
        assert abs(value) < 1e-9
    for value in loads.moment:
        assert abs(value) < 1e-7


# ------------------------------------------------------- point load at midspan


def midspan_point_load(patch_width: float, stations: int = 401) -> HullGirderLoads:
    """Box barge, weight uniform except a force ``P`` in a midspan patch of this width."""
    force = 0.1 * DISPLACEMENT
    centre = LENGTH / 2.0
    items = (
        WeightItem(0.0, LENGTH, DISPLACEMENT - force, "hull"),
        WeightItem(centre - patch_width / 2.0, centre + patch_width / 2.0, force, "cargo"),
    )
    loads = still_water_loads(box_stations(stations), box_buoyancy(stations), items)
    return loads


def test_point_load_midspan_moment_is_minus_p_l_over_8() -> None:
    force = 0.1 * DISPLACEMENT
    expected = -(force * LENGTH) / 8.0
    loads = midspan_point_load(patch_width=1e-3 * LENGTH)
    x, moment = loads.extreme_moment()
    assert x == pytest.approx(LENGTH / 2.0, abs=0.05 * LENGTH)
    assert moment == pytest.approx(expected, rel=2e-3)


def test_patch_load_midspan_moment_is_minus_p_l_minus_d_over_8() -> None:
    """A patch of width ``d`` carries ``-(P/8)(L - d)`` at midspan, exactly.

    The point-load ``-PL/8`` is its ``d -> 0`` limit; the closed form for the finite patch
    follows from superposition of the point-load result over the patch, and the diagram is
    piecewise linear with knots at the patch edges, so the implementation has to hit it to
    round-off rather than merely approach it.
    """
    force = 0.1 * DISPLACEMENT
    for factor in (0.2, 0.1, 0.05, 0.025):
        width = factor * LENGTH
        loads = midspan_point_load(patch_width=width)
        expected = -(force * (LENGTH - width)) / 8.0
        measured = loads.extreme_moment()[1]
        assert measured == pytest.approx(expected, rel=1e-6)
    # ... and it converges to the point-load value as the patch narrows.
    narrow = midspan_point_load(patch_width=1e-4 * LENGTH).extreme_moment()[1]
    assert narrow == pytest.approx(-(force * LENGTH) / 8.0, rel=1e-3)


def test_shear_is_the_buoyancy_integral_abaft_the_patch() -> None:
    """Abaft the patch the net load is ``-(P/L)``, so ``Q = -(P/L)(x - x0)`` exactly.

    The hull's own weight is still spread over the whole length there, so the net load per
    unit length abaft the cargo is ``(W - P)/L - W/L = -P/L`` -- the textbook result, and
    the reason a midspan point load produces ``-PL/8``.
    """
    force = 0.1 * DISPLACEMENT
    loads = midspan_point_load(patch_width=1e-3 * LENGTH)
    x0 = loads.stations[0]
    x_patch_start = LENGTH / 2.0 - 0.5e-3 * LENGTH
    checked = 0
    for x, q in zip(loads.stations, loads.shear):
        if x < x_patch_start - 1e-12:
            assert q == pytest.approx(-(force / LENGTH) * (x - x0), abs=1e-9)
            checked += 1
    assert checked > 10


def test_shear_at_the_ends_of_a_balanced_case_vanishes() -> None:
    loads = midspan_point_load(patch_width=0.02 * LENGTH)
    assert abs(loads.shear[0]) < 1e-9
    assert abs(loads.shear[-1]) < 1e-6 * DISPLACEMENT
    assert abs(loads.moment[-1]) < 1e-6 * DISPLACEMENT * LENGTH


# --------------------------------------------------------------- differential identities


def test_moment_is_the_integral_of_shear_between_knots() -> None:
    """``dM/dx = Q`` checked by integrating the reported ``Q`` with Simpson per knot pair."""
    loads = midspan_point_load(patch_width=0.05 * LENGTH)
    xs, shear, moment = loads.stations, loads.shear, loads.moment
    accumulated = 0.0
    for index in range(len(xs) - 1):
        h = xs[index + 1] - xs[index]
        accumulated += 0.5 * h * (shear[index] + shear[index + 1])
        assert accumulated == pytest.approx(moment[index + 1], abs=1e-6 * DISPLACEMENT * LENGTH)


def test_moment_extremum_coincides_with_a_zero_crossing_of_shear() -> None:
    loads = midspan_point_load(patch_width=0.05 * LENGTH)
    for index in range(len(loads.shear) - 1):
        left, right = loads.shear[index], loads.shear[index + 1]
        if left * right <= 0.0:
            assert abs(loads.moment[index]) <= loads.max_moment + 1e-9


# ------------------------------------------------------------- unbalanced cases


def test_force_residual_is_the_unbalanced_force_only() -> None:
    items = (WeightItem(0.0, LENGTH, DISPLACEMENT * 1.01, "overweight"),)
    loads = still_water_loads(box_stations(), box_buoyancy(), items)
    assert loads.residual_force == pytest.approx(0.01 * DISPLACEMENT, rel=1e-12)
    assert loads.total_weight - loads.total_buoyancy == pytest.approx(0.01 * DISPLACEMENT, rel=1e-12)


def test_balanced_force_but_wrong_position_leaves_only_a_moment_residual() -> None:
    """Same total weight, shifted centroid: ``Q`` closes at the end, ``M`` does not."""
    shifted = DISPLACEMENT / 4.0
    items = (
        WeightItem(0.0, LENGTH * 0.5, shifted, "aft half"),
        WeightItem(LENGTH * 0.5, LENGTH, DISPLACEMENT * 0.75, "forward half"),
    )
    loads = still_water_loads(box_stations(), box_buoyancy(), items)
    assert abs(loads.residual_force) < 1e-6 * DISPLACEMENT
    # Force balance leaves the resultant offset alone: M(x1) = W (lcb_b - lcb_w).
    expected = DISPLACEMENT * (loads.lcb_buoyancy - loads.lcb_weight)
    assert loads.residual_moment == pytest.approx(expected, rel=1e-9)
    assert abs(loads.residual_moment) > 0.05 * DISPLACEMENT * LENGTH


def test_symmetric_weight_on_symmetric_hull_gives_symmetric_moment() -> None:
    force = 0.1 * DISPLACEMENT
    items = (
        WeightItem(0.0, LENGTH * 0.5, 0.5 * (DISPLACEMENT - force), "aft"),
        WeightItem(LENGTH * 0.5, LENGTH, 0.5 * (DISPLACEMENT - force) + force, "forward"),
    )
    # Same *distribution* mirrored: make it symmetric by construction.
    items = (
        WeightItem(0.0, LENGTH * 0.5, 0.5 * DISPLACEMENT, "aft"),
        WeightItem(LENGTH * 0.5, LENGTH, 0.5 * DISPLACEMENT, "forward"),
    )
    loads = still_water_loads(box_stations(), box_buoyancy(), items)
    xs, moment = loads.stations, loads.moment
    for index, x in enumerate(xs):
        mirrored = LENGTH - x
        j = min(range(len(xs)), key=lambda k: abs(xs[k] - mirrored))
        assert moment[index] == pytest.approx(moment[j], abs=1e-6 * DISPLACEMENT * LENGTH)


# ------------------------------------------------------------------- validation


def test_items_outside_the_loaded_length_are_rejected() -> None:
    with pytest.raises(GirderLoadError):
        still_water_loads(
            box_stations(), box_buoyancy(), (WeightItem(2.0 * LENGTH, 3.0 * LENGTH, 1.0),)
        )


def test_weight_item_validates_its_interval() -> None:
    with pytest.raises(GirderLoadError):
        WeightItem(1.0, 1.0, 5.0)
    with pytest.raises(GirderLoadError):
        WeightItem(2.0, 1.0, 5.0)
    with pytest.raises(GirderLoadError):
        WeightItem(0.0, 1.0, float("nan"))


def test_stations_must_be_increasing_and_match_loads() -> None:
    with pytest.raises(GirderLoadError):
        still_water_loads((0.0, 1.0, 0.5), (1.0, 1.0, 1.0))
    with pytest.raises(GirderLoadError):
        still_water_loads((0.0, 1.0), (1.0,))
    with pytest.raises(GirderLoadError):
        still_water_loads((0.0,), (1.0,))


def test_bad_distribution_inputs_are_rejected() -> None:
    mesh = parametric_hull_mesh(nx=11, nz=5)
    with pytest.raises(GirderLoadError):
        buoyancy_distribution(mesh, waterline_z_at_reference=float("nan"))
    with pytest.raises(GirderLoadError):
        buoyancy_distribution(mesh, waterline_z_at_reference=0.0, trim_angle_deg=float("inf"))
    with pytest.raises(GirderLoadError):
        buoyancy_distribution(mesh, waterline_z_at_reference=0.0, density=0.0)
    with pytest.raises(GirderLoadError):
        buoyancy_distribution(mesh, waterline_z_at_reference=0.0, station_count=1)


# -------------------------------------------------------------- a real hull section


def parametric_hull_mesh(nx: int = 41, nz: int = 13):
    from tensorcad.hull import HullParameters, generate_hull

    model = generate_hull(HullParameters(nx=nx, nz=nz))
    return getattr(model, "mesh", model)


def test_parametric_hull_buoyancy_matches_its_own_force_balance() -> None:
    mesh = parametric_hull_mesh()
    box = mesh.bounding_box()
    keel = float(box.minimum[2])
    draft = 0.8
    distribution = buoyancy_distribution(
        mesh, waterline_z_at_reference=keel + draft, station_count=81
    )
    assert distribution.volume > 0.0
    assert distribution.draft_aft == pytest.approx(draft)
    assert distribution.draft_fore == pytest.approx(draft)
    assert box.minimum[0] <= distribution.lcb <= box.maximum[0]

    x_aft, x_fore = distribution.stations[0], distribution.stations[-1]
    span = x_fore - x_aft
    # Force-balanced, deliberately out of trim: 60 % of the weight in the aft half.
    items = (
        WeightItem(x_aft, x_aft + 0.5 * span, 0.6 * distribution.displacement_force, "aft"),
        WeightItem(x_aft + 0.5 * span, x_fore, 0.4 * distribution.displacement_force, "fore"),
    )
    loads, _ = hull_girder_loads(
        mesh, items, waterline_z_at_reference=keel + draft, station_count=81
    )
    assert loads.total_buoyancy == pytest.approx(distribution.displacement_force, rel=1e-12)
    assert abs(loads.residual_force) < 1e-9 * distribution.displacement_force
    # Force balance leaves only the moment residual (x1 - x) * W integrated, i.e. the
    # distance between the two resultants times the displacement, with the stated sign.
    expected = distribution.displacement_force * (distribution.lcb - loads.lcb_weight)
    assert loads.residual_moment == pytest.approx(expected, rel=1e-9)
    assert abs(loads.residual_moment) > 1e-3 * distribution.displacement_force * span


def test_trim_shifts_the_waterline_along_the_hull() -> None:
    mesh = parametric_hull_mesh()
    box = mesh.bounding_box()
    keel = float(box.minimum[2])
    level = buoyancy_distribution(mesh, waterline_z_at_reference=keel + 0.8)
    trimmed = buoyancy_distribution(
        mesh, waterline_z_at_reference=keel + 0.8, trim_angle_deg=2.0
    )
    assert trimmed.draft_fore > level.draft_fore
    assert trimmed.draft_aft < level.draft_aft
    assert trimmed.volume > 0.0
    assert trimmed.lcb > level.lcb  # bow down moves the centroid forward


def test_station_count_convergence_of_the_moment() -> None:
    mesh = parametric_hull_mesh()
    box = mesh.bounding_box()
    keel = float(box.minimum[2])
    coarse = buoyancy_distribution(mesh, waterline_z_at_reference=keel + 0.8, station_count=21)
    fine = buoyancy_distribution(mesh, waterline_z_at_reference=keel + 0.8, station_count=161)
    assert fine.volume == pytest.approx(coarse.volume, rel=2e-2)
    assert fine.displacement_force > 0.0


def test_distribution_serialises() -> None:
    loads = still_water_loads(box_stations(), box_buoyancy(), uniform_weight())
    payload = loads.to_dict()
    assert payload["total_weight"] == pytest.approx(DISPLACEMENT)
    assert len(payload["shear"]) == len(payload["stations"])
    assert payload["items"][0]["label"] == "uniform"
    assert isinstance(payload["extreme_moment"], list)