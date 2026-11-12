"""Anchors for the multi-cell shear flow, checked against the merged single-cell solver.

The references are, in order of independence:

1. **Reduction to the single-cell module**: a two-cell section whose internal web goes to
   zero thickness must reproduce ``tensorcad.shear_flow.shear_flow`` on the surviving walls
   (that solver is the validated one), and the same for a one-cell section and an open tree.
2. **Reduction by symmetry**: a two-cell box with equal cells is symmetric about its centre,
   so the shear centre must land on the centroid to round-off and the transverse resultant
   must vanish -- numbers neither equation set was told.
3. **The classical channel result** ``e = 3 b^2 / (6 b + h)``, which the graph solver
   reproduces without ever being given a circulation formula.
4. **Stationarity**: adding a unit circulation to one cell must leave the complementary
   energy unchanged to first order and raise it by ``delta^2 * SUM L/t`` to second order --
   the no-twist condition *is* the minimum-energy condition, measured, not asserted.
"""
from __future__ import annotations

import math

import pytest
from shear_flow_reduction import (
    canonical_span,
    difference_on_points,
    difference_on_spans,
    flow_scale,
    points_along,
)

from tensorcad.shear_flow import (
    AXIS_Y,
    ShearFlowError,
    ThinWalledSection,
    Wall,
    channel_section,
    i_section,
    shear_centre,
    shear_flow,
)
from tensorcad.shear_flow_multicell import (
    MultiCellSection,
    multi_cell_box_section,
    multi_cell_shear_centre,
    multi_cell_shear_flow,
)

ROUND_OFF = 1e-10


def energy(section: MultiCellSection, flow, *, circulation: float = 0.0, cell: int = 0) -> float:
    """Complementary energy ``SUM INT q^2 / t ds``, optionally with a circulation added."""
    senses = dict(section.cells[cell])
    total = 0.0
    steps = 400
    for index, edge in enumerate(flow.edge_flows):
        length = edge.length
        thickness = edge.thickness
        bump = circulation * senses.get(index, 0)
        for step in range(steps):
            s = length * (step + 0.5) / steps
            c0, c1, c2 = edge.flow_coefficients
            q = c0 + c1 * s + c2 * s * s + bump
            total += length / steps * q * q / thickness
    return total


def shifted(section, dy: float, dz: float):
    """The same walls moved by ``(dy, dz)``: the comparison bodies are the same shape."""
    return [Wall(tuple((y + dy, z + dz) for y, z in wall.points), wall.thickness, wall.label)
            for wall in section.walls]


def close(measured: float, expected: float) -> float:
    if expected == 0.0:
        return abs(measured)
    return abs(measured - expected) / expected


def wall_at(flow, midpoint):
    """The wall whose midpoint is ``midpoint`` (unique: two walls cannot share one)."""
    for edge in flow.edge_flows:
        mid = tuple(0.5 * (a + b) for a, b in zip(edge.start, edge.end))
        if abs(mid[0] - midpoint[0]) < 1e-9 and abs(mid[1] - midpoint[1]) < 1e-9:
            return edge
    raise AssertionError(f"no wall with midpoint {midpoint}")


def _abs_fields(flow):
    """``|q|`` per wall, as a callable of the fraction along it (module and graph alike)."""
    for edge in flow.edge_flows:
        yield lambda fraction, edge=edge: abs(edge.q_at(fraction))


def integral(field, steps: int = 200) -> float:
    """Simpson-free midpoint rule on a polygon-free wall: exact enough for a quadratic."""
    return sum(field((step + 0.5) / steps) for step in range(steps)) / steps


def partner(section: MultiCellSection, flow, start, end):
    """The multi-cell wall whose endpoints are ``{start, end}`` (either direction)."""
    for edge in flow.edge_flows:
        if {tuple(round(v, 9) for v in edge.start), tuple(round(v, 9) for v in edge.end)} == {
            tuple(round(v, 9) for v in start),
            tuple(round(v, 9) for v in end),
        }:
            return edge
    raise AssertionError(f"no wall between {start} and {end}")


def flow_at_point(flow, point) -> float:
    """``q`` at a physical point of the section, on whichever wall contains it."""
    py, pz = point
    for edge in flow.edge_flows:
        (y0, z0), (y1, z1) = edge.start, edge.end
        length = edge.length
        dy, dz = (y1 - y0) / length, (z1 - z0) / length
        along = (py - y0) * dy + (pz - z0) * dz
        if not -1e-9 <= along <= length + 1e-9:
            continue
        if abs((py - y0) * dz - (pz - z0) * dy) > 1e-9:
            continue
        fraction = min(max(along / length, 0.0), 1.0)
        return abs(edge.q_at(fraction))
    raise AssertionError(f"no wall passes through {point}")


def same_physical_point(edge, fraction: float, start, end) -> float:
    """``|q|`` at the point ``start + fraction (end - start)``, whatever the arc sense."""
    forward = tuple(round(v, 9) for v in edge.start) == tuple(round(v, 9) for v in start)
    return abs(edge.q_at(fraction if forward else 1.0 - fraction))


class TestTopology:
    def test_bounded_faces_are_found_for_a_two_cell_box(self):
        section = multi_cell_box_section([2.0, 2.0], 3.0, 0.02)
        assert section.cell_count == 2
        assert section.cycle_rank == 2
        assert section.cell_area(0) == pytest.approx(6.0, rel=0, abs=1e-12)
        assert section.cell_area(1) == pytest.approx(6.0, rel=0, abs=1e-12)
        assert not section.is_open

    def test_three_cells_and_an_asymmetric_row(self):
        section = multi_cell_box_section([1.0, 2.0, 3.0], 2.0, 0.01)
        assert section.cell_count == 3
        assert [section.cell_area(i) for i in range(3)] == pytest.approx([2.0, 4.0, 6.0], abs=1e-12)

    def test_an_open_tree_has_no_cell(self):
        section = MultiCellSection(i_section(depth=3.0, flange_width=2.0, web_thickness=0.02,
                                            flange_thickness=0.03).walls)
        assert section.cell_count == 0
        assert section.is_open

    def test_a_crossing_wall_is_refused(self):
        walls = [
            Wall(((0.0, 0.0), (4.0, 0.0)), 0.02, "a"),
            Wall(((0.0, 2.0), (4.0, 2.0)), 0.02, "b"),
            Wall(((0.0, 1.0), (0.0, 3.0)), 0.02, "c"),
        ]
        with pytest.raises(ShearFlowError):
            MultiCellSection(walls)

    def test_a_doubled_wall_is_refused(self):
        walls = [
            Wall(((0.0, 0.0), (2.0, 0.0)), 0.02, "a"),
            Wall(((2.0, 0.0), (0.0, 0.0)), 0.02, "b"),
        ]
        with pytest.raises(ShearFlowError):
            MultiCellSection(walls)


def reference_box() -> ThinWalledSection:
    """The merged solver on the same outer walls the two-cell builder lays down.

    Labelled as in the smoke comparison: ``bottom``/``side_s``/``top``/``side_p`` run
    counter-clockwise, so the module's closed-loop plan and this module's edge order differ
    only by the global arc sense (which both resolve from equilibrium).
    """
    return ThinWalledSection(
        [
            Wall(((0.0, 0.0), (4.0, 0.0)), 0.02, "bottom"),
            Wall(((4.0, 0.0), (4.0, 3.0)), 0.02, "side_s"),
            Wall(((4.0, 3.0), (0.0, 3.0)), 0.02, "top"),
            Wall(((0.0, 3.0), (0.0, 0.0)), 0.02, "side_p"),
        ]
    )


class TestReductions:
    """The graph solver reduced to the merged one, compared as a *field*.

    These three comparisons were ``xfail`` items with the note "the pointwise reduction is
    not yet wired up correctly in this test": the two solvers agreed at round-off while the
    *comparison* did not. ``tests/shear_flow_reduction.py`` records the diagnosis (a wall's
    arc sense is not shared between the modules; a relative tolerance is meaningless at a
    zero crossing; a spliced flange is not an unspliced one). Everything here is now
    compared on physical points or on canonical wall spans, and the tolerance is a fraction
    of the section's own flow scale -- so a wall carrying zero flow is not compared against
    a tolerance proportional to zero.
    """

    TOLERANCE = 1e-9  # of the largest |q| anywhere in the section

    def test_two_cells_with_a_vanishing_web_reproduce_the_single_cell_box(self):
        """The limit of a vanishing internal web is the merged single-cell box, at O(t_web).

        The limit is not round-off and cannot be: a wall whose thickness vanishes still
        carries a finite *stress*, so it carries a vanishing *flow* (``q_web = tau t_web``),
        and the perturbation it exerts on the rest of the section is first order in
        ``t_web``. The assertion is the convergence *rate*: refine the web by 100x and the
        difference falls by ~100x. Measured over seven decades; the ratio
        ``difference / t_web`` holds two significant digits.
        """
        module = shear_flow(reference_box(), 8000.0)
        scale = flow_scale(module)
        side = {edge.label: edge for edge in module.edges if edge.label.startswith("side")}
        spans = [canonical_span(side[label]) for label in ("side_s", "side_p")]
        previous = 0.0
        ratios: list[float] = []
        for inner in (1e-3, 1e-5, 1e-7, 1e-9):
            section = multi_cell_box_section([2.0, 2.0], 3.0, 0.02, inner_thickness=inner)
            flow = multi_cell_shear_flow(section, 8000.0)
            worst = difference_on_spans(flow, module, spans)
            if previous:
                ratios.append(previous / worst)
            previous = worst
        assert previous < 1e-6 * scale, previous / scale
        # Two refinements of two decades each: a first-order limit gives 100, not 1 or 1e4.
        assert min(ratios) > 50.0, ratios
        assert max(ratios) < 200.0, ratios

    def test_one_cell_matches_the_merged_box_pointwise(self):
        walls = multi_cell_box_section([4.0], 3.0, 0.02).walls
        section = MultiCellSection(walls)
        assert section.cell_count == 1
        flow = multi_cell_shear_flow(section, 8000.0)
        module = shear_flow(reference_box(), 8000.0)
        scale = flow_scale(module)
        # Wall by wall, on the spans both bodies have...
        assert difference_on_spans(flow, module) < TestReductions.TOLERANCE * scale
        # ...and at physical points, which is where the flange midpoints live: q is zero
        # there (symmetry), and that zero is what a relative tolerance called a 100 % error.
        points = points_along(reference_box().walls)
        assert difference_on_points(flow, module, points) < TestReductions.TOLERANCE * scale

    def test_an_open_tree_matches_the_merged_open_section(self):
        section = i_section(depth=3.0, flange_width=2.0, web_thickness=0.02,
                            flange_thickness=0.03)
        graph = MultiCellSection(section.walls)
        assert graph.cell_count == 0
        flow = multi_cell_shear_flow(graph, 12000.0)
        module = shear_flow(section, 12000.0)
        assert flow.resultant[1] == pytest.approx(12000.0, rel=ROUND_OFF)
        assert flow.resultant[0] == pytest.approx(0.0, abs=1e-9)
        assert len(flow.edge_flows) == len(module.edges)
        scale = flow_scale(module)
        # The merged module sweeps this tree from its free edges inwards, so half of its
        # walls run the other way; on physical points that is not visible.
        points = points_along(section.walls, count=11)
        assert difference_on_points(flow, module, points) < TestReductions.TOLERANCE * scale

    def test_the_channel_shear_centre_still_lands_on_the_classical_formula(self):
        channel = channel_section(depth=3.0, flange_width=1.0, thickness=0.02)
        graph = MultiCellSection(channel.walls)
        centre = multi_cell_shear_centre(graph)
        formula = 3 * 1.0**2 / (6 * 1.0 + 3.0)
        merged = shear_centre(channel)
        assert abs(abs(centre.y) - formula) / formula < ROUND_OFF
        assert centre.y == pytest.approx(merged.y, rel=ROUND_OFF)
        assert abs(centre.z - merged.z) < 1e-9


class TestTwoCellPhysics:
    def test_the_coupled_cells_satisfy_every_condition_at_round_off(self):
        section = multi_cell_box_section([2.0, 2.0], 3.0, 0.02)
        flow = multi_cell_shear_flow(section, 8000.0)
        assert flow.equilibrium_residual < ROUND_OFF
        assert flow.transverse_residual < ROUND_OFF
        assert flow.continuity_residual < 1e-9 * 8000.0
        assert flow.twist_residual < ROUND_OFF
        assert all(junction.net_outflow == pytest.approx(0.0, abs=1e-9)
                   for junction in flow.junctions)

    def test_a_symmetric_two_cell_box_has_its_shear_centre_on_the_centroid(self):
        section = multi_cell_box_section([2.0, 2.0], 3.0, 0.02)
        centre = multi_cell_shear_centre(section)
        assert centre.y == pytest.approx(section.centroid_y, abs=1e-12)
        assert centre.z == pytest.approx(section.centroid_z, abs=1e-12)

    def test_the_shared_web_carries_flow_so_the_coupling_is_actually_exercised(self):
        """A solver that quietly solved each cell apart would leave this wall at zero."""
        section = multi_cell_box_section([2.0, 2.0], 3.0, 0.02)
        flow = multi_cell_shear_flow(section, 8000.0)
        web = flow.edge("box.web1")[0]
        peak = max(abs(web.q_at(f)) for f in (0.0, 0.25, 0.5, 0.75, 1.0))
        assert peak > 100.0, peak
        # The two cells are mirror images in y = 2, so the flange flows must mirror too:
        # equal magnitude at mirrored points, with the arc sense flipping. (The shared web
        # is NOT odd about mid-depth -- a section symmetric in y still lets the middle wall
        # carry flow -- which is exactly the trap a "shear centre = centroid, so every wall
        # is symmetric" shortcut would fall into.)
        top0 = flow.edge("box.top0")[0]
        top1 = flow.edge("box.top1")[0]
        for fraction in (0.0, 0.25, 0.5, 0.75, 1.0):
            point0 = (2.0 * fraction, 3.0)
            point1 = (4.0 - point0[0], 3.0)
            assert flow_at_point(flow, point0) == pytest.approx(
                flow_at_point(flow, point1), rel=1e-9
            )
        assert top0.q_at(0.0) != pytest.approx(top1.q_at(0.0), rel=1e-9)

    def test_unequal_cells_move_the_shear_centre_off_the_axis(self):
        section = multi_cell_box_section([1.0, 3.0], 2.0, 0.02)
        centre = multi_cell_shear_centre(section)
        assert abs(centre.y - section.centroid_y) > 1e-3
        flow = multi_cell_shear_flow(section, 5000.0)
        assert flow.continuity_residual < 1e-9 * 5000.0
        assert flow.twist_residual < ROUND_OFF

    def test_the_merged_module_refuses_what_this_module_solves(self):
        section = multi_cell_box_section([2.0, 2.0], 3.0, 0.02)
        with pytest.raises(ShearFlowError):
            ThinWalledSection(section.walls)
        assert multi_cell_shear_flow(section, 1.0).cell_count == 2

    def test_a_dead_end_branch_carries_flow_and_stops_at_its_free_edge(self):
        walls = multi_cell_box_section([2.0, 2.0], 3.0, 0.02).walls + (
            Wall(((2.0, 3.0), (2.0, 4.0)), 0.02, "stub"),
        )
        section = MultiCellSection(walls)
        assert section.cell_count == 2
        flow = multi_cell_shear_flow(section, 4000.0)
        stub = flow.edge("stub")[0]
        assert flow.continuity_residual < 1e-9 * 4000.0
        assert max(abs(stub.q_at(f)) for f in (0.0, 0.5, 1.0)) > 1.0


class TestEnergyAndRefusals:
    def test_the_no_twist_condition_is_the_minimum_energy_condition(self):
        section = multi_cell_box_section([2.0, 3.0], 2.5, 0.02)
        flow = multi_cell_shear_flow(section, 6000.0)
        base = energy(section, flow)
        length_over_t = sum(
            flow.edge_flows[index].length / flow.edge_flows[index].thickness
            for index, _ in section.cells[0]
        )
        delta = 50.0
        bumped = energy(section, flow, circulation=delta)
        # First order vanishes (the twist integral is zero), second order is delta^2 SUM L/t.
        assert abs(bumped - base - delta * delta * length_over_t) / base < 1e-9
        assert bumped > base

    def test_a_product_of_inertia_is_refused(self):
        angle = math.radians(20.0)
        c, s = math.cos(angle), math.sin(angle)

        def rotate(point):
            y, z = point
            return (c * y - s * z, s * y + c * z)

        base = i_section(depth=3.0, flange_width=2.0, web_thickness=0.03,
                         flange_thickness=0.05)
        walls = [Wall(tuple(rotate(p) for p in wall.points), wall.thickness, wall.label)
                 for wall in base.walls]
        section = MultiCellSection(walls)
        assert abs(section.product_of_inertia) > 0.0
        with pytest.raises(ShearFlowError, match="product of inertia"):
            multi_cell_shear_flow(section, 1000.0)

    def test_a_disconnected_or_zero_force_input_is_refused(self):
        section = multi_cell_box_section([2.0, 2.0], 3.0, 0.02)
        with pytest.raises(ShearFlowError):
            multi_cell_shear_flow(section, 0.0)
        with pytest.raises(ShearFlowError):
            multi_cell_shear_flow(section, 1000.0, axis="x")
        with pytest.raises(ShearFlowError):
            multi_cell_shear_flow(section, 1000.0, stations=1)

    def test_the_horizontal_shear_is_solved_too(self):
        section = multi_cell_box_section([2.0, 1.0], 3.0, 0.02)
        flow = multi_cell_shear_flow(section, 2000.0, axis=AXIS_Y)
        assert flow.resultant[0] == pytest.approx(2000.0, rel=ROUND_OFF)
        assert flow.twist_residual < ROUND_OFF
        assert flow.axis == AXIS_Y
        assert flow.cell_count == 2
