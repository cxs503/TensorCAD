"""Thin-walled shear flow on a **multi-cell** section, and the shear centre that follows.

``tensorcad.shear_flow`` solves two graph topologies: one open tree, or one closed loop. Both
topologies a real hull girder section has -- a double bottom, a double hull, a wing-tank box --
are *neither*: they are two or more closed cells sharing walls, so the merged module refuses
them by design (its docstring says so: "the two-box section of a bulk carrier is *not* this").
This module answers that case, with the same thin-walled idealisation so the two cannot
disagree about what the geometry is.

The extra physics, and nothing else
-----------------------------------
A single closed cell leaves one unknown (the circulation ``q0``) fixed by Bredt's no-twist
condition ``CLOSED_INT q/t ds = 0``. With ``m`` cells there are ``m`` unknowns -- one
circulation per cell -- and the *same* condition applies to each cell separately, because pure
bending must produce **no twist in any cell**:

    ``CLOSED_INT_cell_c q / t ds = 0``   for every cell ``c``.

A shared wall appears in both of its cells' integrals with opposite sense, and that is the
whole coupling: it is why a two-box section is not "two single cells solved apart".

Method: a graph, not a topology switch
--------------------------------------
The section is an undirected graph of straight walls in the ``(y, z)`` plane (``y``
starboard, ``z`` up, metres -- the frame the rest of the package uses). The flow on each wall
is written

    ``q_e(s) = c_e + (V / I) t_e [ (v_a - v_c) s + (v_b - v_a) s^2 / (2 L_e) ]``

which is the integrated form of the local axial-equilibrium relation
``dq/ds = (V / I) t (v - v_c)`` for a shear force ``V`` about the axis whose coordinate is
``v`` (``v = z`` for the vertical shear ``V_z``, ``v = y`` for ``V_y``). The unknown is one
constant per wall, ``c_e = q_e(0)``, so there are ``E`` unknowns and the equations are

* **continuity at every vertex but one** (``V - 1`` equations): the net outflow of a vertex
  vanishes. This also *implies* ``q = 0`` at the free edge of a dead-end branch -- no special
  case, exactly the way :mod:`tensorcad.girder` treats a level cut;
* **no twist in every cell** (``m`` equations).

``V - 1 + m = E`` by Euler's formula for a connected planar graph whose bounded faces are the
cells, so the system is square -- and the cells are *found*, not supplied: a half-edge walk
around the planar embedding recovers the bounded faces, each as a signed set of walls. A
dead-end branch (a longitudinal stiffener, say) belongs to no face; its constant is fixed by
continuity alone, which is right, because it is statically determinate.

Because the walls are straight, every integral above is a polynomial in the arc length and is
evaluated in closed form, so the residuals are round-off rather than a mesh error: there is no
discretisation here to converge.

What this module is **not**
---------------------------
* Not a **product of inertia** solver. ``I_yz != 0`` couples the two directions through the
  principal-axis rotation; such a section is refused (rotate it to its principal axes first,
  which is what the refusal says). Same rule as the single-cell module.
* Not **shear lag / effective breadth**, not **plastic** shear, not a **warping** or torsional
  solution, and not a strength, buckling or fatigue check -- ``q/t`` is the elastic
  thin-walled shear stress for the shear force the caller supplies.
* Not a **variable-thickness wall**: one thickness per wall, as in the single-cell module.
* Not a **curved** wall model: walls are polylines, and the chord error of a curved wall is
  the caller's, exactly as in :mod:`tensorcad.shear_flow`.

Sign convention
---------------
``q`` is positive from a wall's ``start`` to its ``end``. The equations fix the field up to
one global sign, which corresponds to nothing physical (it is the choice of which side of a
cut the flow is measured on); the conventional one -- the one whose resultant reproduces the
applied shear force -- is imposed last and asserted by the tests, the same way
:mod:`tensorcad.shear_flow` fixes it.

Only the Python standard library is required.
"""
from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

from .shear_flow import (
    AXIS_Y,
    AXIS_Z,
    DEFAULT_MERGE_TOLERANCE,
    DEFAULT_STATIONS,
    ShearFlowError,
    Wall,
)

Point = tuple[float, float]

__all__ = [
    "CellEdgeFlow",
    "CellReport",
    "JunctionReport",
    "MultiCellFlow",
    "MultiCellSection",
    "MultiCellShearCentre",
    "multi_cell_box_section",
    "multi_cell_shear_centre",
    "multi_cell_shear_flow",
]

#: Relative tolerance below which a section is treated as having ``I_yz = 0``.
PRODUCT_TOLERANCE = 1e-12

#: Three-point Gauss-Legendre on ``[-1, 1]``: exact for the cubic integrands below.
_GAUSS3 = (
    (-math.sqrt(0.6), 5.0 / 9.0),
    (0.0, 8.0 / 9.0),
    (math.sqrt(0.6), 5.0 / 9.0),
)


@dataclass(frozen=True)
class _Segment:
    a: int
    b: int
    thickness: float
    label: str


class MultiCellSection:
    """A connected wall graph of straight polyline walls in the ``(y, z)`` section plane.

    Coincident vertices are merged, a segment covered twice is refused, and the bounded faces
    of the planar embedding are recovered as the cells, so a caller never declares the
    topology by hand (and cannot declare it wrongly).
    """

    def __init__(
        self,
        walls: Iterable[Wall],
        *,
        merge_tolerance: float = DEFAULT_MERGE_TOLERANCE,
    ) -> None:
        wall_list = list(walls)
        if not wall_list:
            raise ShearFlowError("a section needs at least one wall")
        if not math.isfinite(merge_tolerance) or merge_tolerance <= 0.0:
            raise ShearFlowError("merge_tolerance must be finite and positive")

        nodes: list[Point] = []
        lookup: dict[tuple[int, int], int] = {}

        def node_id(point: Point) -> int:
            key = (
                int(round(point[0] / merge_tolerance)),
                int(round(point[1] / merge_tolerance)),
            )
            found = lookup.get(key)
            if found is not None:
                return found
            lookup[key] = len(nodes)
            nodes.append(point)
            return len(nodes) - 1

        segments: list[_Segment] = []
        seen: dict[tuple[int, int], str] = {}
        for wall in wall_list:
            ids = [node_id(point) for point in wall.points]
            for first, second in zip(ids, ids[1:]):
                if first == second:
                    raise ShearFlowError(
                        f"wall {wall.label!r} has a zero-length segment (vertices merged); "
                        "drop the duplicated point or use a finer merge_tolerance"
                    )
                key = (min(first, second), max(first, second))
                if key in seen:
                    raise ShearFlowError(
                        f"walls {seen[key]!r} and {wall.label!r} coincide on one segment; "
                        "a doubled wall would be counted twice"
                    )
                seen[key] = wall.label
                segments.append(
                    _Segment(a=first, b=second, thickness=wall.thickness, label=wall.label)
                )

        if len(nodes) < 2:
            raise ShearFlowError("the section needs at least two distinct vertices")

        self.walls = tuple(wall_list)
        self._nodes = tuple(nodes)
        self._segments = tuple(segments)
        self.merge_tolerance = float(merge_tolerance)

        adjacency: list[list[tuple[int, int]]] = [[] for _ in range(len(nodes))]
        for index, segment in enumerate(segments):
            adjacency[segment.a].append((segment.b, index))
            adjacency[segment.b].append((segment.a, index))
        self._adjacency = tuple(tuple(row) for row in adjacency)
        if any(len(row) == 0 for row in adjacency):
            raise ShearFlowError("the section has an isolated vertex")

        # One connected component, or the cell walk and the continuity system are both wrong.
        reachable = [False] * len(nodes)
        reachable[0] = True
        stack = [0]
        while stack:
            current = stack.pop()
            for neighbour, _ in adjacency[current]:
                if not reachable[neighbour]:
                    reachable[neighbour] = True
                    stack.append(neighbour)
        if not all(reachable):
            raise ShearFlowError(
                "the walls do not form one connected section; a disconnected section is "
                "several shear-flow problems, not one"
            )

        self._lengths = tuple(self._length(segment) for segment in segments)
        self._cells = self._find_cells()
        cycle_rank = len(segments) - len(nodes) + 1
        if cycle_rank != len(self._cells):
            raise ShearFlowError(
                f"the wall graph has cycle rank {cycle_rank} but {len(self._cells)} bounded "
                "faces; the walls are not a planar embedding (crossing walls are not a "
                "thin-walled section)"
            )
        self.cycle_rank = cycle_rank
        self._compute_properties()

    # -- geometry -----------------------------------------------------------
    def _length(self, segment: _Segment) -> float:
        (y0, z0), (y1, z1) = self._nodes[segment.a], self._nodes[segment.b]
        return math.hypot(y1 - y0, z1 - z0)

    @property
    def is_open(self) -> bool:
        """``True`` when the graph encloses no cell at all (a tree of walls)."""
        return not self._cells

    @property
    def cell_count(self) -> int:
        return len(self._cells)

    @property
    def node_count(self) -> int:
        return len(self._nodes)

    @property
    def edge_count(self) -> int:
        return len(self._segments)

    @property
    def cells(self) -> tuple[tuple[tuple[int, int], ...], ...]:
        """Each cell as ``(segment index, sense)`` pairs, ``sense = +1`` for ``a -> b``."""
        return self._cells

    def _find_cells(self) -> tuple[tuple[tuple[int, int], ...], ...]:
        """Recover the bounded faces of the planar embedding by a half-edge walk.

        Around every vertex the incident half-edges are sorted by direction angle. Leaving a
        vertex on the half-edge that is the next one *clockwise* from the edge just arrived on
        keeps the enclosed area on the left, so the walk closes on the faces of the embedding;
        bounded faces come out with positive signed area and the outer face negative, which is
        how they are told apart without a bounding box.
        """
        nodes = self._nodes
        segments = self._segments
        order: list[list[int]] = []
        for node, row in enumerate(self._adjacency):
            y0, z0 = nodes[node]
            decorated = []
            for neighbour, index in row:
                y1, z1 = nodes[neighbour]
                decorated.append((math.atan2(z1 - z0, y1 - y0), index))
            decorated.sort()
            order.append([index for _, index in decorated])
        position = {}
        for node, row in enumerate(order):
            for slot, index in enumerate(row):
                position[(node, index)] = slot

        visited: set[tuple[int, int]] = set()
        faces: list[tuple[tuple[int, int], ...]] = []
        for index, segment in enumerate(segments):
            for direction in (1, -1):
                start = (index, direction)
                if start in visited:
                    continue
                loop: list[tuple[int, int]] = []
                key = start
                guard = 0
                while key not in visited:
                    visited.add(key)
                    edge_index, sense = key
                    walked = segments[edge_index]
                    here = walked.a if sense > 0 else walked.b
                    loop.append((edge_index, sense))
                    there = walked.b if here == walked.a else walked.a
                    row = order[there]
                    slot = position[(there, edge_index)]
                    following = row[(slot - 1) % len(row)]
                    nxt = segments[following]
                    key = (following, 1 if nxt.a == there else -1)
                    guard += 1
                    if guard > 4 * len(segments) + 4:
                        raise ShearFlowError("the face walk did not close")
                faces.append(tuple(loop))

        cells: list[tuple[tuple[int, int], ...]] = []
        for loop in faces:
            if self._signed_area(loop) <= 0.0:
                continue
            # Net multiplicity: a dead-end branch is walked out and back and cancels.
            net: dict[int, int] = {}
            for edge_index, sense in loop:
                net[edge_index] = net.get(edge_index, 0) + sense
            kept = tuple((index, sign) for index, sign in sorted(net.items()) if sign != 0)
            if kept:
                cells.append(kept)
        return tuple(cells)

    def _signed_area(self, loop: Sequence[tuple[int, int]]) -> float:
        total = 0.0
        for edge_index, sense in loop:
            segment = self._segments[edge_index]
            first = segment.a if sense > 0 else segment.b
            second = segment.b if sense > 0 else segment.a
            y0, z0 = self._nodes[first]
            y1, z1 = self._nodes[second]
            total += y0 * z1 - y1 * z0
        return 0.5 * total

    def cell_area(self, cell: int) -> float:
        """Enclosed area of a cell (m^2): positive for every cell this module keeps."""
        return self._signed_area(self._cells[cell])

    # -- section properties -------------------------------------------------
    def _compute_properties(self) -> None:
        area = 0.0
        sum_y = 0.0
        sum_z = 0.0
        for segment, length in zip(self._segments, self._lengths):
            (y0, z0), (y1, z1) = self._nodes[segment.a], self._nodes[segment.b]
            strip = length * segment.thickness
            area += strip
            sum_y += strip * 0.5 * (y0 + y1)
            sum_z += strip * 0.5 * (z0 + z1)
        if area <= 0.0:
            raise ShearFlowError("the section has no wall area")
        c_y = sum_y / area
        c_z = sum_z / area

        second_z = 0.0  # INT (z - c_z)^2 dA: the bending-about-y inertia
        second_y = 0.0  # INT (y - c_y)^2 dA: the bending-about-z inertia
        product = 0.0
        for segment, length in zip(self._segments, self._lengths):
            (y0, z0), (y1, z1) = self._nodes[segment.a], self._nodes[segment.b]
            weight = length * segment.thickness
            ay, by = y0 - c_y, y1 - c_y
            az, bz = z0 - c_z, z1 - c_z
            second_z += weight * (az * az + az * bz + bz * bz) / 3.0
            second_y += weight * (ay * ay + ay * by + by * by) / 3.0
            product += weight * (ay * az + 0.5 * (ay * bz + by * az) + by * bz) / 3.0

        self.area = area
        self.centroid_y = c_y
        self.centroid_z = c_z
        self.second_moment_horizontal = second_z
        self.second_moment_vertical = second_y
        self.product_of_inertia = product
        self.perimeter = sum(self._lengths)

    def require_uncoupled(self) -> None:
        """Refuse a section whose product of inertia is not zero (principal-axis rotation)."""
        scale = math.sqrt(self.second_moment_horizontal * self.second_moment_vertical)
        if scale <= 0.0:
            raise ShearFlowError("the section has no second moment about one of the axes")
        if abs(self.product_of_inertia) > PRODUCT_TOLERANCE * scale:
            raise ShearFlowError(
                "the section has a non-zero product of inertia (I_yz = "
                f"{self.product_of_inertia:.6g}); vertical and horizontal shear flow are "
                "coupled through the principal-axis rotation. Rotate the walls to the "
                "principal axes and solve there"
            )

    def to_dict(self) -> dict[str, object]:
        return {
            "cell_count": self.cell_count,
            "cycle_rank": self.cycle_rank,
            "vertex_count": self.node_count,
            "edge_count": self.edge_count,
            "area": self.area,
            "centroid": [self.centroid_y, self.centroid_z],
            "second_moment_horizontal": self.second_moment_horizontal,
            "second_moment_vertical": self.second_moment_vertical,
            "product_of_inertia": self.product_of_inertia,
            "perimeter": self.perimeter,
            "cell_areas": [self.cell_area(index) for index in range(self.cell_count)],
        }


def _solve_linear(matrix: Sequence[Sequence[float]], rhs: Sequence[float]) -> list[float]:
    """Gaussian elimination with partial pivoting (standard library only)."""
    size = len(rhs)
    augmented = [list(row) + [float(value)] for row, value in zip(matrix, rhs)]
    scale = max((abs(value) for row in matrix for value in row), default=0.0)
    if scale == 0.0:
        raise ShearFlowError("the shear-flow system is empty")
    for column in range(size):
        pivot = max(range(column, size), key=lambda row: abs(augmented[row][column]))
        if abs(augmented[pivot][column]) <= 1e-14 * scale:
            raise ShearFlowError(
                "the shear-flow system is singular: the cells do not fix one circulation "
                "each (check that no cell is a duplicate of another)"
            )
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        lead = augmented[column][column]
        for row in range(column + 1, size):
            factor = augmented[row][column] / lead
            if factor == 0.0:
                continue
            for since in range(column, size + 1):
                augmented[row][since] -= factor * augmented[column][since]
    solution = [0.0] * size
    for row in range(size - 1, -1, -1):
        total = augmented[row][size]
        for column in range(row + 1, size):
            total -= augmented[row][column] * solution[column]
        solution[row] = total / augmented[row][row]
    return solution


@dataclass(frozen=True)
class CellEdgeFlow:
    """The flow on one wall, ``q(s) = c0 + c1 s + c2 s^2`` from ``start`` to ``end``."""

    label: str
    thickness: float
    start: Point
    end: Point
    length: float
    flow_coefficients: tuple[float, float, float]
    s_values: tuple[float, ...]
    q_values: tuple[float, ...]
    stress_values: tuple[float, ...]

    def q_at(self, fraction: float) -> float:
        if not 0.0 <= fraction <= 1.0:
            raise ShearFlowError("fraction must lie in [0, 1]")
        c0, c1, c2 = self.flow_coefficients
        s = fraction * self.length
        return c0 + c1 * s + c2 * s * s

    def stress_at(self, fraction: float) -> float:
        return self.q_at(fraction) / self.thickness

    def peak_stress(self) -> tuple[float, float]:
        """``(stress, fraction)`` at the stationary point of the quadratic, else the ends."""
        c0, c1, c2 = self.flow_coefficients
        candidates = [(self.stress_at(0.0), 0.0), (self.stress_at(1.0), 1.0)]
        if c2 != 0.0:
            stationary = -c1 / (2.0 * c2)
            if 0.0 < stationary < self.length:
                q = c0 + c1 * stationary + c2 * stationary * stationary
                candidates.append((q / self.thickness, stationary / self.length))
        return max(candidates, key=lambda item: abs(item[0]))

    def to_dict(self) -> dict[str, object]:
        return {
            "label": self.label,
            "thickness": self.thickness,
            "start": list(self.start),
            "end": list(self.end),
            "length": self.length,
            "flow_coefficients": list(self.flow_coefficients),
            "s_values": list(self.s_values),
            "q_values": list(self.q_values),
            "stress_values": list(self.stress_values),
        }


@dataclass(frozen=True)
class CellReport:
    """One cell's twist condition and the area it encloses."""

    index: int
    area: float
    edges: tuple[int, ...]
    twist_integral: float
    twist_residual: float

    def to_dict(self) -> dict[str, object]:
        return {
            "index": self.index,
            "area": self.area,
            "edges": list(self.edges),
            "twist_integral": self.twist_integral,
            "twist_residual": self.twist_residual,
        }


@dataclass(frozen=True)
class JunctionReport:
    """Net outflow of a vertex: zero when the flow is divergence-free on the graph."""

    point: Point
    degree: int
    net_outflow: float

    def to_dict(self) -> dict[str, object]:
        return {
            "point": list(self.point),
            "degree": self.degree,
            "net_outflow": self.net_outflow,
        }


@dataclass(frozen=True)
class MultiCellFlow:
    """The solved flow field of a multi-cell section for one shear force."""

    axis: str
    shear_force: float
    second_moment: float
    cell_count: int
    edge_flows: tuple[CellEdgeFlow, ...]
    cells: tuple[CellReport, ...]
    junctions: tuple[JunctionReport, ...]
    resultant: tuple[float, float]
    torque_about_centroid: float
    equilibrium_residual: float
    transverse_residual: float
    continuity_residual: float
    twist_residual: float
    global_sign: int

    def edge(self, label: str) -> tuple[CellEdgeFlow, ...]:
        """Every wall carrying ``label`` (a wall split at a junction keeps its label)."""
        found = tuple(flow for flow in self.edge_flows if flow.label == label)
        if not found:
            raise ShearFlowError(f"no wall labelled {label!r} in this section")
        return found

    def to_dict(self) -> dict[str, object]:
        return {
            "axis": self.axis,
            "shear_force": self.shear_force,
            "second_moment": self.second_moment,
            "cell_count": self.cell_count,
            "edges": [flow.to_dict() for flow in self.edge_flows],
            "cells": [cell.to_dict() for cell in self.cells],
            "junctions": [junction.to_dict() for junction in self.junctions],
            "resultant": list(self.resultant),
            "torque_about_centroid": self.torque_about_centroid,
            "equilibrium_residual": self.equilibrium_residual,
            "transverse_residual": self.transverse_residual,
            "continuity_residual": self.continuity_residual,
            "twist_residual": self.twist_residual,
            "global_sign": self.global_sign,
        }


@dataclass(frozen=True)
class MultiCellShearCentre:
    """The point a shear force must pass through to twist no cell."""

    y: float
    z: float
    torque_vertical: float
    torque_horizontal: float
    vertical_flow: MultiCellFlow
    horizontal_flow: MultiCellFlow

    def to_dict(self) -> dict[str, object]:
        return {
            "y": self.y,
            "z": self.z,
            "torque_vertical": self.torque_vertical,
            "torque_horizontal": self.torque_horizontal,
        }


def _centroid(section: MultiCellSection, axis: str) -> float:
    return section.centroid_z if axis == AXIS_Z else section.centroid_y


def _coordinate(section: MultiCellSection, node: int, axis: str) -> float:
    y, z = section._nodes[node]
    return z if axis == AXIS_Z else y


def _solve_field(
    section: MultiCellSection, axis: str, shear_force: float, stations: int
) -> MultiCellFlow:
    section.require_uncoupled()
    if not math.isfinite(shear_force):
        raise ShearFlowError("shear_force must be finite")
    if shear_force == 0.0:
        raise ShearFlowError("a zero shear force has no shear flow to report")

    second_moment = (
        section.second_moment_horizontal if axis == AXIS_Z else section.second_moment_vertical
    )
    if second_moment <= 0.0:
        raise ShearFlowError("the section has no second moment about the requested axis")
    scale = shear_force / second_moment
    centroid = _centroid(section, axis)

    segments = section._segments
    lengths = section._lengths
    edge_count = len(segments)

    increments: list[float] = []
    base_integrals: list[float] = []
    for segment, length in zip(segments, lengths):
        v0 = _coordinate(section, segment.a, axis)
        v1 = _coordinate(section, segment.b, axis)
        thickness = segment.thickness
        increments.append(scale * thickness * length * (0.5 * (v0 + v1) - centroid))
        # INT_0^L (L - s)(v(s) - v_c) ds = L^2 (2 v0 + v1 - 3 v_c) / 6
        base_integrals.append(
            scale * thickness * length * length * (2.0 * v0 + v1 - 3.0 * centroid) / 6.0
        )

    rows: list[list[float]] = []
    rhs: list[float] = []

    # Continuity at every vertex but one: the dropped row is minus the sum of the others
    # (the total outflow of a closed graph vanishes), so it carries no information.
    for node in range(section.node_count - 1):
        row = [0.0] * edge_count
        value = 0.0
        for _, index in section._adjacency[node]:
            segment = segments[index]
            if segment.a == node:
                row[index] += 1.0
            else:
                row[index] -= 1.0
                value += increments[index]
        rows.append(row)
        rhs.append(value)

    # No twist in each cell: SUM (q / t) ds = 0. The t in q cancels in the base part.
    for cell in section._cells:
        row = [0.0] * edge_count
        value = 0.0
        for index, sense in cell:
            row[index] += sense * lengths[index] / segments[index].thickness
            value -= sense * base_integrals[index] / segments[index].thickness
        rows.append(row)
        rhs.append(value)

    if len(rows) != edge_count:
        raise ShearFlowError(
            f"the shear-flow system has {len(rows)} equations for {edge_count} unknowns; "
            "the section is degenerate (a zero-length wall or a doubled cell)"
        )

    constants = _solve_linear(rows, rhs)

    def coefficients(sign: int) -> list[tuple[float, float, float]]:
        built: list[tuple[float, float, float]] = []
        for segment, length, constant in zip(segments, lengths, constants):
            v0 = _coordinate(section, segment.a, axis)
            v1 = _coordinate(section, segment.b, axis)
            thickness = segment.thickness
            built.append(
                (
                    sign * constant,
                    sign * scale * thickness * (v0 - centroid),
                    sign * scale * thickness * (v1 - v0) / (2.0 * length),
                )
            )
        return built

    flow = _assemble(section, axis, shear_force, second_moment, coefficients(1), stations)
    measured = flow.resultant[1] if axis == AXIS_Z else flow.resultant[0]
    if measured * shear_force < 0.0:
        flow = _assemble(
            section, axis, shear_force, second_moment, coefficients(-1), stations, sign=-1
        )
    return flow


def _assemble(
    section: MultiCellSection,
    axis: str,
    shear_force: float,
    second_moment: float,
    coefficients: Sequence[tuple[float, float, float]],
    stations: int,
    *,
    sign: int = 1,
) -> MultiCellFlow:
    segments = section._segments
    lengths = section._lengths
    nodes = section._nodes

    edge_flows: list[CellEdgeFlow] = []
    for segment, length, poly in zip(segments, lengths, coefficients):
        s_values = tuple(length * index / (stations - 1) for index in range(stations))
        q_values = tuple(poly[0] + poly[1] * s + poly[2] * s * s for s in s_values)
        edge_flows.append(
            CellEdgeFlow(
                label=segment.label,
                thickness=segment.thickness,
                start=nodes[segment.a],
                end=nodes[segment.b],
                length=length,
                flow_coefficients=(poly[0], poly[1], poly[2]),
                s_values=s_values,
                q_values=q_values,
                stress_values=tuple(q / segment.thickness for q in q_values),
            )
        )

    def integrate(index: int, weight_of) -> float:
        """Three-point Gauss on a cubic integrand (exact for a cubic)."""
        segment = segments[index]
        length = lengths[index]
        poly = coefficients[index]
        y0, z0 = nodes[segment.a]
        y1, z1 = nodes[segment.b]
        dy = (y1 - y0) / length
        dz = (z1 - z0) / length
        total = 0.0
        for node, weight in _GAUSS3:
            s = 0.5 * length * (1.0 + node)
            q = poly[0] + poly[1] * s + poly[2] * s * s
            y = y0 + dy * s
            z = z0 + dz * s
            total += 0.5 * length * weight * q * weight_of(y, z, dy, dz)
        return total

    resultant_y = sum(
        integrate(index, lambda y, z, dy, dz: dy) for index in range(len(segments))
    )
    resultant_z = sum(
        integrate(index, lambda y, z, dy, dz: dz) for index in range(len(segments))
    )
    c_y, c_z = section.centroid_y, section.centroid_z
    torque = sum(
        integrate(index, lambda y, z, dy, dz: (y - c_y) * dz - (z - c_z) * dy)
        for index in range(len(segments))
    )

    junctions: list[JunctionReport] = []
    worst_continuity = 0.0
    for node, row in enumerate(section._adjacency):
        net = 0.0
        for _, index in row:
            segment = segments[index]
            poly = coefficients[index]
            if segment.a == node:
                net += poly[0]
            else:
                length = lengths[index]
                net -= poly[0] + poly[1] * length + poly[2] * length * length
        worst_continuity = max(worst_continuity, abs(net))
        junctions.append(JunctionReport(point=nodes[node], degree=len(row), net_outflow=net))

    cells: list[CellReport] = []
    worst_twist = 0.0
    length_over_t = sum(
        length / segments[index].thickness for index, length in enumerate(lengths)
    )
    peak = max((abs(value) for poly in coefficients for value in poly), default=0.0)
    for index, cell in enumerate(section._cells):
        integral = 0.0
        for edge_index, sense in cell:
            segment = segments[edge_index]
            length = lengths[edge_index]
            poly = coefficients[edge_index]
            total = (
                poly[0] * length
                + poly[1] * length * length / 2.0
                + poly[2] * length**3 / 3.0
            )
            integral += sense * total / segment.thickness
        denominator = peak * length_over_t
        residual = integral / denominator if denominator else 0.0
        worst_twist = max(worst_twist, abs(residual))
        cells.append(
            CellReport(
                index=index,
                area=section.cell_area(index),
                edges=tuple(edge_index for edge_index, _ in cell),
                twist_integral=integral,
                twist_residual=residual,
            )
        )

    magnitude = abs(shear_force)
    measured = resultant_z if axis == AXIS_Z else resultant_y
    other = resultant_y if axis == AXIS_Z else resultant_z
    return MultiCellFlow(
        axis=axis,
        shear_force=shear_force,
        second_moment=second_moment,
        cell_count=section.cell_count,
        edge_flows=tuple(edge_flows),
        cells=tuple(cells),
        junctions=tuple(junctions),
        resultant=(resultant_y, resultant_z),
        torque_about_centroid=torque,
        equilibrium_residual=abs(measured - shear_force) / magnitude,
        transverse_residual=abs(other) / magnitude,
        continuity_residual=worst_continuity,
        twist_residual=worst_twist,
        global_sign=sign,
    )


def multi_cell_shear_flow(
    section: MultiCellSection,
    shear_force: float,
    *,
    axis: str = AXIS_Z,
    stations: int = DEFAULT_STATIONS,
) -> MultiCellFlow:
    """Solve the thin-walled shear flow on a section with any number of closed cells.

    ``shear_force`` is in newtons and acts along ``axis`` (``"z"`` vertical by default);
    ``stations`` is how many points per wall the report carries (the flow is quadratic
    everywhere, so five is already the whole field).
    """
    if axis not in (AXIS_Y, AXIS_Z):
        raise ShearFlowError(f"axis must be {AXIS_Z!r} or {AXIS_Y!r}, got {axis!r}")
    if stations < 2:
        raise ShearFlowError("stations must be at least 2 so that both wall ends appear")
    return _solve_field(section, axis, shear_force, stations)


def multi_cell_shear_centre(
    section: MultiCellSection, *, unit_force: float = 1.0
) -> MultiCellShearCentre:
    """Locate the shear centre of a multi-cell section from the flow fields themselves.

    Both directions are solved, because the shear centre is a point: the vertical shear's
    resultant acts at ``y_s = c_y + T_z / V_z``, the horizontal shear's at
    ``z_s = c_z - T_y / V_y``. Nothing is assumed about where it is.
    """
    if not math.isfinite(unit_force) or unit_force == 0.0:
        raise ShearFlowError("unit_force must be finite and non-zero")
    vertical = multi_cell_shear_flow(section, unit_force, axis=AXIS_Z)
    horizontal = multi_cell_shear_flow(section, unit_force, axis=AXIS_Y)
    return MultiCellShearCentre(
        y=section.centroid_y + vertical.torque_about_centroid / unit_force,
        z=section.centroid_z - horizontal.torque_about_centroid / unit_force,
        torque_vertical=vertical.torque_about_centroid,
        torque_horizontal=horizontal.torque_about_centroid,
        vertical_flow=vertical,
        horizontal_flow=horizontal,
    )


def multi_cell_box_section(
    cell_widths: Sequence[float],
    height: float,
    thickness: float,
    *,
    inner_thickness: float | None = None,
    label: str = "box",
) -> MultiCellSection:
    """A row of ``len(cell_widths)`` rectangular cells sharing their internal webs.

    ``cell_widths`` are the clear widths in order (from ``y = 0`` to starboard); ``height`` is
    the overall depth; ``thickness`` is used for the outer walls and ``inner_thickness``
    (default: the same) for the internal webs. This is the geometry a double bottom is
    idealised as before the floors are cut in, and the geometry the anchors are written for.
    """
    widths = [float(value) for value in cell_widths]
    if not widths:
        raise ShearFlowError("a multi-cell box needs at least one cell")
    for width in widths:
        if not math.isfinite(width) or width <= 0.0:
            raise ShearFlowError(f"cell widths must be finite and positive, got {width!r}")
    if not math.isfinite(height) or height <= 0.0:
        raise ShearFlowError(f"height must be finite and positive, got {height!r}")
    if not math.isfinite(thickness) or thickness <= 0.0:
        raise ShearFlowError(f"thickness must be finite and positive, got {thickness!r}")
    inner = thickness if inner_thickness is None else float(inner_thickness)
    if not math.isfinite(inner) or inner <= 0.0:
        raise ShearFlowError(f"inner_thickness must be finite and positive, got {inner!r}")

    stations = [0.0]
    for width in widths:
        stations.append(stations[-1] + width)
    # The flanges are built one span per cell: a wall must be split where another wall meets
    # it, otherwise the meeting point is an interior point of a segment rather than a vertex
    # and the walls are not connected in the graph. (Same rule `tensorcad.shear_flow`'s own
    # section builders follow, and the trap a caller notices only when the graph is walked.)
    walls: list[Wall] = []
    for index in range(len(widths)):
        y0, y1 = stations[index], stations[index + 1]
        walls.append(Wall(((y0, height), (y1, height)), thickness, f"{label}.top{index}"))
        walls.append(Wall(((y0, 0.0), (y1, 0.0)), thickness, f"{label}.bottom{index}"))
    walls.append(Wall(((stations[-1], 0.0), (stations[-1], height)), thickness, f"{label}.side_s"))
    walls.append(Wall(((stations[0], height), (stations[0], 0.0)), thickness, f"{label}.side_p"))
    for index, station in enumerate(stations[1:-1], start=1):
        walls.append(Wall(((station, 0.0), (station, height)), inner, f"{label}.web{index}"))
    return MultiCellSection(walls)
